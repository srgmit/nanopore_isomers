function [tf,tknock_timeseries,num_dangling_bonds_timeseries,num_dangling_bonds_CH_timeseries,num_dangling_bonds_CH2_timeseries,num_AC,num_ZZ,num_UA,num_5R]=kmc_isomers_edge_diffusion(isoindex,poresize,basedir)
% This function performs a single kinetic Monte Carlo simulation for the
% formation of a nanopore in a 2D hexagonal lattice and returns the final
% configuration of the nanopore, in the presence of edge diffusion effects

% Inputs: isoindex - index of the isomer, poresize - size of the nanopore
% in terms of the number of atoms removed, basedir - directory in which to
% store the generated files

% Outputs: tf - formation time of nanopore, tknock_timeseries - time
% required for atom to be knocked out at each timestep,
% num_dangling_bonds_timeseries - number of dangling bonds as a function of
% time, num_dangling_bonds_CH_timeseries - number of doubly bonded carbons, 
% num_dangling_bonds_CH2_timeseries - number of singly bonded carbons, 
% num_AC - number of armchair edge atoms, num_ZZ - number of zigzag edge atoms,
% num_UA - number of unassigned edge atoms (in this study, etching and diffusion
% barriers of UA = that of AC),num_5R - number of 5-membered rings 

global basedirectory N  T R NA kb kb_eV count_etches count_diffusions nuperp nuperp_etch area_site lattice_const num_dt isomer_index pore_size;

% Initialize time
t0=0;

% Set isomer index, pore size, and base directory global variables
isomer_index = isoindex;
pore_size = poresize;
basedirectory = basedir;

% Set grid size
N = 18;

% Number of KMC steps required to form pore of given size
num_dt = pore_size-1;

% Counters for etching and diffusion events of various types
% Indices (1 through 5) - SC,AC,ZZ,UA,FC
count_etches = zeros(1,5);
count_diffusions = zeros(5,5,3);

% Values of physical constants
R=8.314;        % Gas constant in J/mol-K
NA=6.022e23;    % Avogadro's number in mol^-1
kb=R/NA;        % Boltzmann constant in SI units
kb_eV=8.6173324e-5;% Boltzmann constant in eV/K

% Set attempt frequency for A and B sublattices for etching and edge
% diffusion events
nuperp = 1e13*[1,1];
nuperp_etch = 1e13*[1,1];

% Set temperature of system in Kelvin
T = 500+273.15;

% Set lattice size of graphene
lattice_const = 1.422591*sqrt(3);
area_site = 0.50*sqrt(3)*(lattice_const^2)*1e-20;

% Initialize 2D material layer with a monoatom vacancy
[~,C,num_atoms]=init(N);

% Do the KMC algorithm
tic
[tf,tknock_timeseries,num_dangling_bonds_timeseries,num_dangling_bonds_CH_timeseries,num_dangling_bonds_CH2_timeseries,num_AC,num_ZZ,num_UA,num_5R]=kmc(C,num_atoms,t0,N);
toc

end

function [t,tknock_timeseries,num_dangling_bonds_timeseries,num_dangling_bonds_CH_timeseries,num_dangling_bonds_CH2_timeseries,num_AC,num_ZZ,num_UA,num_5R]=kmc(C0,num_atoms,t0,N)
% This function performs the KMC algorithm and returns final configuration
% of the system

% Inputs: initial surface configuration (C0), number of atoms on A and B sublattices (num_atoms), initial time (t0), size of
% grid (N)

% Outputs: tf - formation time of nanopore, tknock_timeseries - time
% required for atom to be knocked out at each timestep,
% num_dangling_bonds_timeseries - number of dangling bonds as a function of
% time, num_dangling_bonds_CH_timeseries - number of doubly bonded carbons, 
% num_dangling_bonds_CH2_timeseries - number of singly bonded carbons, 
% num_AC - number of armchair edge atoms, num_ZZ - number of zigzag edge atoms,
% num_UA - number of unassigned edge atoms (in this study, etching and diffusion
% barriers of UA atoms = that of AC atoms),num_5R - number of 5-membered rings 

global num_dt time current_pore_size superbasin_sites;

t=t0; % initial time
count=0; % counter for number of time steps
C=C0; % store initial state of system in variable C

% Variable to store lattice sites which are located in a superbasin state
superbasin_sites = [];

% Determine the list of sites which can be etched
sites_list=sites(C,N);
r = [sites_list(:,7); sites_list(:,8)]; % vector of all rates; 7th column is etching, 8th column is edge diffusion
num=length(r);

% Initialize variables to store timeseries
num_dangling_bonds_timeseries =  zeros(num_dt+1,1);
num_dangling_bonds_CH_timeseries =  zeros(num_dt+1,1);
num_dangling_bonds_CH2_timeseries =  zeros(num_dt+1,1);
tknock_timeseries = zeros(num_dt+1,1);

current_pore_size = 1;

max_steps=2000;
dt=0;
while (num>0 &&  current_pore_size<=num_dt && count<=max_steps) 

    % Generate two uniformly-distributed random numbers u1 and u2
    random = rand(1,2);
    u1=random(1);
    u2=random(2);
      
    % Compute effective rate vector
    num=length(r);
    rcum = cumsum(r);
    rcum=rcum/rcum(num); % normalize with total rate
    
    % Store number of dangling bonds vs time
    [data1,data2,data3] = count_dangling_bonds(sites_list);
    num_dangling_bonds_timeseries(count+1) = data1;
    num_dangling_bonds_CH_timeseries(count+1) = data2;
    num_dangling_bonds_CH2_timeseries(count+1) = data3;
    tknock_timeseries(count+1) = 1/sum(r);

    % Find which of the sites is going to have a reaction, based on
    % Gillespie's algorithm
    site_sel=0;
    for i=1:num       
        if (i==1 && u1<=rcum(i))
            site_sel=i;
        elseif (i>1 && u1>rcum(i-1) && u1<=rcum(i))
            site_sel=i;
        end     
    end
    
    % Increment the system time based on a Poisson process    
    rnet=sum(r);
    if (rnet>0)
        dt = log(1/u2)/rnet;
        t=t+dt;  
        time = t;
        count=count+1;
    end
    
    if (site_sel>0)
        
        if (site_sel<=num/2)
            % the case for vacancy addition, i.e., atom etching (mode=0)
            mode = sites_list(site_sel,6);           
        elseif (site_sel>num/2)
            % the case for site diffusion
            mode = 2;
            site_sel = site_sel-(num/2);
        end       
        
        % Find (r,c) position of the lattice site selected to be etched
        pos_sel_r = sites_list(site_sel,4); 
        pos_sel_c = sites_list(site_sel,5);

        % Find sublattice at which etching is going to occurr
        sublattice_sel = sites_list(site_sel,1);

        site_changed = [sublattice_sel,pos_sel_r,pos_sel_c,mode];
        
        % Update list of available sites
        [sites_list,C]=sites_change(C,sublattice_sel,[pos_sel_r,pos_sel_c],sites_list,mode,count);      
        
        if (mode==0)  % case of a vacancy being added
           current_pore_size = current_pore_size+1;
           % C matrix will be changed as part of sites_change function            

        else  % case of atom/vacancy diffusing from source to destination lattice site

            % current_pore_size need not change
            % C matrix will be changed as part of sites_change function
        end           
        
        % Compute the next step's rates
        r = [sites_list(:,7); sites_list(:,8)];        
    end
    
    if (mod(count,100)==0)
        disp(['count = ', num2str(count)]);
        disp(['t = ',num2str(t),' s']);
        disp(['Current pore size = ', num2str(current_pore_size)]);
    end
end

% Extract and save the number of dangling bonds, and the knocking time at
% each time step
[data1,data2,data3] = count_dangling_bonds(sites_list);
num_dangling_bonds_timeseries(end+1) = data1;
num_dangling_bonds_CH_timeseries(end+1) = data2;
num_dangling_bonds_CH2_timeseries(end+1) = data3;
tknock_timeseries(end+1) = 1/sum(r);

% Save XYZ file and antimolecule adjacency matrix at the last timestep
[num_AC,num_ZZ,num_UA,num_5R]=visualize(C,N,count+1);       
end

function [num_dangling_bonds,num_dangling_bonds_CH,num_dangling_bonds_CH2]=count_dangling_bonds(sites_list)
% This function determines the number of total, singly-bonded, and doubly-bonded
% dangling bonds in the system

% Input: list of atomic sites in the system
% Output: total number of dangling bonds (num_dangling_bonds),
% doubly-bonded dangling bond (num_dangling_bonds_CH), 
% singly-bonded dangling bond (num_dangling_bonds_CH2)

% Initialize the count for all dangling bonds to zero
num_dangling_bonds = 0;
num_dangling_bonds_CH = 0;
num_dangling_bonds_CH2 = 0;

% Loop through all the atomic sites
for i=1:size(sites_list,1)
    if (sites_list(i,6)==1)  % If the site is a vacancy
        if (sites_list(i,2)>0) % If number of nearest neighboring vacancies is greater than zero
             num_dangling_bonds = num_dangling_bonds + (3-sites_list(i,2));
        end
    else                     % If the site has an atom
        if (sites_list(i,2)==2)     % 2 vacancy neighbors, so a 'CH2' type dangling bond
            num_dangling_bonds_CH2 = num_dangling_bonds_CH2 + 1; 
        elseif (sites_list(i,2)==1) % 1 vacancy neighbor,  so a 'CH'  type dangling bond
            num_dangling_bonds_CH = num_dangling_bonds_CH + 1;
        end
    end
end
end

function y=atomtype(i,j,sublattice,C)
% This function classifies a given lattice site to be an armchair, zigzag,
% singly-bonded, lone, or fully-bonded carbon atom

% Inputs: location (i,j) of lattice site, sublattice, state of the system (C)
% Output: atomtype (AC, ZZ, SB, LC, FC)

% Find the list of nearest neighbor lattice sites to the given lattice site
[nearest,~]=neighbor_sites(i,j,sublattice);

% Initialize number of vacancy neighbors, number of filled neighbors, and
% number of filled neighbors of filled neighbors, all to be zero
num_vacancy_neigh = 0;
num_atomic_neigh  = 0;
num_atomic_neigh_atomic_neigh = zeros(2,1);

% Loop through all nearest neighbor lattice sites
for i=1:size(nearest,1)
    if (C(nearest(i,1),nearest(i,2),nearest(i,3))==0) % nearest site is not a vacancy, so it is filled with an atom
        num_atomic_neigh = num_atomic_neigh+1;
        how_many_vacancy_neigh_atomic_neigh = count_vacancy_neighbors(nearest(i,1),nearest(i,2),nearest(i,3),C); % this function counts the number of vacancy neighbors
        num_atomic_neigh_atomic_neigh(num_atomic_neigh) = 3-(how_many_vacancy_neigh_atomic_neigh(1));
    else
        num_vacancy_neigh = num_vacancy_neigh+1;
    end
end

if (num_vacancy_neigh==0) % 0 filled vacancy neighbors, i.e. 3 filled atom neighbors
    y = 'FC';  % fully coordinated
elseif (num_vacancy_neigh==1) % 1 filled vacancy neighbors, i.e. 2 filled atomic neighbors - armchair or zigzag
    if ((num_atomic_neigh_atomic_neigh(1)==2 && num_atomic_neigh_atomic_neigh(2)==3) || (num_atomic_neigh_atomic_neigh(1)==3 && num_atomic_neigh_atomic_neigh(2)==2))     %  armchair
        y = 'AC';
    elseif (num_atomic_neigh_atomic_neigh(1)==3 && num_atomic_neigh_atomic_neigh(2)==3)  % zigzag
        y = 'ZZ';
    else % unassigned edge, treat similar to armchair
        y = 'UA';
    end
elseif (num_vacancy_neigh==2) % 1 filled atomic neighbor - singly-coordinated atom
    if ((num_atomic_neigh_atomic_neigh(1)==3 && num_atomic_neigh_atomic_neigh(2)==0) || (num_atomic_neigh_atomic_neigh(1)==0 && num_atomic_neigh_atomic_neigh(2)==3))  % actual KL
        y = 'KL';
    else
        y = 'LC';  % lone C-C cluster
    end
elseif (num_vacancy_neigh==3) % 0 filled atomic neighbors, is possible sometimes (e.g., zigzag edge etching)
    y = 'LC'; % lone carbon
end
end

function y=barrier(varargin)
% This function returns the etching barrier for a given atomtype, or the
% edge diffusion barrier for diffusion from a given lattice site to another
% lattice site
% Input: atomtype(s) - one (for etching) or two (for edge diffusion) in
% number, which neighboring (1st, 2nd, or 3rd) site for edge diffusion,
% optional flags denoting whether the initial and final lattice sites are
% superbasin states
% Output: etching or edge diffusion barrier (in eV)

global superbasin_sites;

if (nargin<=2)
    atom_type_init = varargin{1};
        switch atom_type_init
            case 'LC'
                y = 0;
            case 'AC'
                y = 2.28986495;
            case 'ZZ'
                y = 2.3010539487;
            case 'UA'
                y = 2.28986495;
            case 'KL'
                y = 1.02986589;
            case 'FC'
                y = inf;
        end

    if (size(superbasin_sites,1)>=50)
        is_superbasin_site = varargin{2};
        if (~is_superbasin_site)
            y=inf;
        end
    end

else
    atom_type_init  = varargin{1};
    atom_type_final = varargin{2};
    which_neighbor  = varargin{3};
   
    switch atom_type_init
        case 'LC'
            atom_type_init_id = 1;
        case 'KL'
            atom_type_init_id = 2;
        case 'UA'
            atom_type_init_id = 3;
        case 'AC'
            atom_type_init_id = 3;
        case 'ZZ'
            atom_type_init_id = 4;
        case 'FC'
            atom_type_init_id = 5;
    end
    
    switch atom_type_final
        case 'LC'
            atom_type_final_id = 1;
        case 'KL'
            atom_type_final_id = 2;            
        case 'UA'
            atom_type_final_id = 3;
        case 'AC'
            atom_type_final_id = 3;
        case 'ZZ'
            atom_type_final_id = 4;
        case 'FC'
            atom_type_final_id = 5;
    end
    
    barrier_matrix = zeros(5,5,3);

    barrier_matrix(:,:,1) =  [ 0         0         0       0       0
                               8.8029    0         0	   0       inf
                               7.6386    inf       2.691   1.372   inf
                               9.6079    inf       7.409   1.211   inf
                              16.9716    inf       inf     inf     inf];

    barrier_matrix(:,:,2) = [  0          0        0        0       0     
                              8.8029      0	       2.689	2.100	1.947
                              7.6386      4.794	   4.853	5.357	2.325
                              9.6079      3.463	   3.580	6.592	2.906   
                             16.9716      inf	   5.580	2.906	3.479];

    barrier_matrix(:,:,3) = [  0          0        0        0       0
                               8.8029     0	       2.510	0.867	4.566
                               7.6386     4.615	   6.885	6.479	inf
                               9.6079     inf	   4.666	2.856	inf
                              16.9716     4.819	   inf      inf	    5.336];

    y  = barrier_matrix(atom_type_init_id,atom_type_final_id,which_neighbor);

    if (size(superbasin_sites,1)>=50)
            is_source_superbasin_site = varargin{4};
            is_destination_superbasin_site = varargin{5};
        if (~is_source_superbasin_site || is_destination_superbasin_site)
            y=inf;
        end
    end
end

end

function y=atomtype2num(atomtype)
switch atomtype
    case 'LC'
        y = 1;
    case 'KL'
        y = 2;
    case 'UA'
        y = 3;
    case 'AC'
        y = 3;
    case 'ZZ'
        y = 4;
    case 'FC'
        y = 5;
end
end

function sites_list=sites(C,N)
% Compute the positions of all possible relavent interactions
% The columns are 1. sublattice, 2. number of vacancy neighbors, 
% 3. number of vacancy next neighbors, 4. row-position (y),
% 5. column-position (x), 6. vacancy flag (1-occupied by vacancy,
% 0-occupied by atom), 7. etching rate, 8. total edge diffusion rate

global T kb_eV nuperp nuperp_etch;

% Initialize monoatom vacancy in lattice, if need other shape, change
% init() function
[initial_vacancy,~,~]=init(N);

[num,~]=size(initial_vacancy);
sites_list=zeros(num,7);
sites_list(:,[1,4,5])=initial_vacancy(:,[1,2,3]); % store sublattice, row location, column location in 1st, 4th, and 5th columns of sites_list

for i=1:num
    % This has been made a loop so that if one wants to initialize with a
    % defect other than a monoatom vacancy, the code can be easily modified simply by changing `initial_vacancy`
    sublattice = initial_vacancy(i,1);
    pos = initial_vacancy(i,2:3);

    [nearest,next_nearest]=neighbor_sites(pos(1),pos(2),sublattice);
    [n_near, ~] = size(nearest);
    [n_nextnear, ~] = size(next_nearest);

    % Update vacancy flag (i.e., 6th column) to be 1
    [~,locb] =   ismember([sublattice,pos(1),pos(2)],sites_list(:,[1,4,5]),'rows');
    sites_list(locb,6) = 1;

    % Loop through nearest neighbor lattice sites    
    for k=1:n_near
        posr = nearest(k,1);
        posc = nearest(k,2);
        neigh_spec = nearest(k,3);
        
        [lia,locb] =   ismember([neigh_spec,posr,posc],sites_list(:,[1,4,5]),'rows');        
        if (~any(lia)) % If lattice site is not in sites_list, add it there
            sites_list = [sites_list; neigh_spec,1,0,posr,posc,ismember([neigh_spec,posr,posc],initial_vacancy,'rows'),0];
        else % If lattice site is already there, update number of vacancy neighbors
            sites_list(locb,2) = sites_list(locb,2)+1;
        end
    end
    
    for k=1:n_nextnear
        posr = next_nearest(k,1);
        posc = next_nearest(k,2);
        neigh_spec = next_nearest(k,3);
            
        [lia,locb] =   ismember([neigh_spec,posr,posc],sites_list(:,[1,4,5]),'rows');
        if (~any(lia)) % If lattice site is not in sites_list, add it there
            sites_list = [sites_list; neigh_spec,0,1,posr,posc,ismember([neigh_spec,posr,posc],initial_vacancy,'rows'),0];
        else % If lattice site is already there, update number of vacancy next neighbors
            sites_list(locb,3) = sites_list(locb,3)+1;
        end
    end
    
end   

% Loop through all lattice sites in the sites list
for i=1:size(sites_list,1)
    sublattice = sites_list(i,1);
    num_vacancy_neigh = sites_list(i,2); 
    num_next_neigh = sites_list(i,3);
    vacancy_flag = sites_list(i,6);
    posr = sites_list(i,4);
    posc = sites_list(i,5);
    
    if (vacancy_flag==0)  % site doesn't have vacancy
        if (num_vacancy_neigh>0) % site has at least one vacancy neighbor
            Ea = barrier(atomtype(posr,posc,sublattice,C));
            sites_list(i,7)=nuperp_etch(sublattice)*exp(-Ea/(kb_eV*T));  % Rate constant for etching of atom
        else
            sites_list(i,7)=0; % site has no vacancy neighbor, so cannot be etched
        end
        sites_list(i,8)=0; % site has no vacancy, so diffusion rate is zero
        
    else  % site has a vacancy
        sites_list(i,7) = 0; % vacancies cannot be removed, but can diffuse around    
        
        [num,possible_sites]=diffusion_sites(posr,posc,sublattice,C);
        r = zeros(num,1);
        for j=1:num
            C_try = C;
            C_try(posr,               posc,               sublattice)             = 0;  % remove atom from source site
            C_try(possible_sites(j,1),possible_sites(j,2),possible_sites(j,3)) = possible_sites(j,3);  % add atom to destination site
                       
            if (C(possible_sites(j,1),possible_sites(j,2),possible_sites(j,3)==0))   % if diffusion site is empty

                % Remember that the barriers are listed from the initial
                % and final perspectives of a filled atom, and hence the
                % entries for source and destination are interchanged
                % However, the initial type of the filled atom will still be at
                % the initial state vector of the system
                Ea = barrier(atomtype(possible_sites(j,1),possible_sites(j,2),possible_sites(j,3),C), atomtype(posr,posc,sublattice,C_try), possible_sites(j,4));
                r(j) = nuperp(sublattice)*exp(-Ea/(kb_eV*T));  % Rate constant for diffusion
            else
                r(j)  = 0; % if adjacent diffusion site is filled 
            end
        end               
        rtot = sum(r);
        sites_list(i,8)=rtot;       
    end   
end
end

function y=ismembercus(varargin)
    A=varargin{1};
    B=varargin{2};
    if (nargin>2)
        C=varargin{3};       
    end
    if (isempty(B))
        y = 0;
    else
        y = ismember(A,B,C);
    end
end

function [updated_sites_list,C_new]=sites_change(C,sublattice,pos,sites_list,mode,count)
% This function updates the sites_list for the simulation and returns an
% updated list of potential sites for etching and edge diffusion

% Inputs: state of the system (C), properties of currently etched site
% (sublattice, pos), current sites list (sites_list)
% Outputs: updated sites list and list of modified lattice sites

global T kb_eV nuperp count_etches count_diffusions nuperp_etch superbasin_sites;

% Create variables for new sites list and new state of the system
updated_sites_list = sites_list;
C_new = C;

% Initialize list of indices and counter for modified sites
modified_sites_index = zeros(15,1);
modified_sites_count = 0;

if (mode==2)    % edge diffusion
    [lia,~] =   ismember([sublattice,pos(1),pos(2)],sites_list(:,[1,4,5]),'rows');
    % Perform housekeeping check
    if (~any(lia))
        error('Selected site not in candidate list');
    end
    
    % Generate random diffusion site starting at candidate site
    [num,possible_sites]=diffusion_sites(pos(1),pos(2),sublattice,C);
    % Perform another housekeeping check
    if (num==0)
        disp(['Count = ', num2str(count)]);
        selected_site = [sublattice, pos(1), pos(2)];
        sites_list
        tot_rate=sum(sites_list(:,7))
        error('No sites possible for edge diffusion')
    end   
    
    r = zeros(num,1);
    Ea_list = zeros(num,1);
    sites_info = cell(num,5);
    for j=1:num
        
        C_try = C;
        C_try(pos(1), pos(2), sublattice) = 0;  % remove vacancy from source site
        C_try(possible_sites(j,1),possible_sites(j,2),possible_sites(j,3)) = possible_sites(j,3);  % add vacancy to destination site
        
        possible_site_neighbors = count_vacancy_neighbors(possible_sites(j,1),possible_sites(j,2),possible_sites(j,3),C_try);
        possible_site_num_neigh = possible_site_neighbors(1);
        
        if (C(possible_sites(j,1),possible_sites(j,2),possible_sites(j,3))==0) % &&  possible_site_num_neigh>0)   % if diffusion site is empty and it has atleast one neighbor
            
            % Remember that the barriers are listed from the initial
            % and final perspectives of a filled atom, and hence the
            % entries for source and destination are interchanged
            % However, the initial type of the filled atom will still be at
            % the initial state vector of the system
                        
            filled_site_type_initial = atomtype(possible_sites(j,1),possible_sites(j,2),possible_sites(j,3),C);
            filled_site_type_final = atomtype(pos(1),pos(2),sublattice,C_try);
            is_source_superbasin_site = any(ismembercus([possible_sites(j,1),possible_sites(j,2),possible_sites(j,3)],superbasin_sites,'rows'));
            is_destination_superbasin_site = any(ismembercus([pos(1),pos(2),sublattice],superbasin_sites,'rows'));
            Ea = barrier(filled_site_type_initial, filled_site_type_final , possible_sites(j,4),is_source_superbasin_site,is_destination_superbasin_site);
            r(j) = nuperp(sublattice)*exp(-Ea/(kb_eV*T));  % Rate constant for diffusion
            Ea_list(j) = Ea;
            sites_info{j,1} = filled_site_type_initial;
            sites_info{j,2} = filled_site_type_final;
            sites_info{j,3} = possible_sites(j,4);
            sites_info{j,4} = Ea;
            sites_info{j,5} = r(j);
            sites_info{j,6} = possible_sites(j,1);
            sites_info{j,7} = possible_sites(j,2);
            sites_info{j,8} = possible_sites(j,3);
        else
            r(j)  = 0; % if adjacent diffusion site is filled
        end
    end
    rcum = cumsum(r);
    rcum=rcum/rcum(num); % normalize with total rate
    
    % find which of the sites is going to be the destination diffusion site
    random = rand(2,1);
    u1 = random(1);
    u2 = random(2);
    site_sel=0;
    for j=1:num       
        if (j==1 && u1<=rcum(j))
            site_sel=j;
        elseif (j>1 && u1>rcum(j-1) && u1<=rcum(j))
            site_sel=j;
        end     
    end
    
    rnet=sum(r);
    if (rnet>0)
        dt = log(1/u2)/rnet;
    else
        sites_info
        vacancy_source_site=[pos(1),pos(2),sublattice]
    end
    selected_site = possible_sites(site_sel,:);
    Ea_selected_site = sites_info{site_sel,4};
    if (Ea_selected_site == 0 || Ea_selected_site == 1.986)
        superbasin_sites = [superbasin_sites;
                            pos(1), pos(2), sublattice;
                            selected_site(1), selected_site(2), selected_site(3)                       
                            ];
    else
        superbasin_sites = [];
    end

    type_selected_site_init = atomtype2num(sites_info{site_sel,1});
    type_selected_site_final = atomtype2num(sites_info{site_sel,2});
    type_diffusion_event = sites_info{site_sel,3};
    count_diffusions(type_selected_site_init,type_selected_site_final,type_diffusion_event) =  count_diffusions(type_selected_site_init,type_selected_site_final,type_diffusion_event)+1;
        
    % Check whether generated site is suitable
    C_try=C;
    if (C_try(pos(1),pos(2),sublattice)==0)
        mode
        [sublattice,pos]        
        error('source site cannot be empty');
    end
    if (C_try(selected_site(1),selected_site(2),selected_site(3))==selected_site(3))
        error('destination site cannot be filled');
    end
    C_try(pos(1),pos(2),sublattice) = 0;  % remove atom from source site
    C_try(selected_site(1),selected_site(2),selected_site(3))=selected_site(3);       % add atom to destination site
    
    % compute number of neighbors of the source and destination sites
    updated_num_neighbors_source = count_vacancy_neighbors(pos(1),pos(2),sublattice,C_try); % of the source site
    updated_num_neighbors_destination = count_vacancy_neighbors(selected_site(1),selected_site(2),selected_site(3),C_try); % of the destination site
    
    % at the end of this step, no vacancy should exist on its own (have zero neighbors), 
    % i.e. the destination site (which is now a vacancy) should not have
    % zero neighbors, i.e. updated_num_neighbors_destination(1) > 0
    
    % further, the source vacancy site should have less than 3 neighbors, 
    % i.e. the destination atom site should not have zero neighbors   
    allow_diffusion  = 0;
    if (1==1) %updated_num_neighbors_destination(1) > 0 && updated_num_neighbors_source(1)<3)
            allow_diffusion=1;
    end
      
    if (allow_diffusion==1) %updated_num_neighbors(1)>=sites_list(locb,2) && updated_num_neighbors(2)>=sites_list(locb,3)) % check if number of neighors and next-neighbors have increased
              
        % Accept the selected diffusion site
        % Update the concentration matrix
        C_new=C_try;

        added = [];
        % Increase the number of neighbors for the neighbors of the
        % destination site, except if it is the source site
        % Increase the number of next-nearest neighbors for the next-nearest neighbors of the destination site
        [nearest,next_nearest]=neighbor_sites(selected_site(1),selected_site(2),selected_site(3));
        [n_near, ~] = size(nearest);
        [n_nextnear, ~] = size(next_nearest);
        for k=1:n_near
            posr = nearest(k,1);
            posc = nearest(k,2);
            neigh_spec = nearest(k,3);
            [lia,locb] =   ismember([neigh_spec,posr,posc],updated_sites_list(:,[1,4,5]),'rows');
            
            if (any(lia))% && (posr~=pos(1) && posc~=pos(2)))
                updated_sites_list(locb,2) = updated_sites_list(locb,2)+1;
                modified_sites_count = modified_sites_count + 1;
                modified_sites_index(modified_sites_count) = locb;
            else
               updated_sites_list(end+1,:) = [neigh_spec, 1 , 0, posr, posc, C_new(posr,posc,neigh_spec)==neigh_spec  , 0 , 0]; 
               modified_sites_count = modified_sites_count + 1;
               modified_sites_index(modified_sites_count) = size(updated_sites_list,1);
               added = [added;neigh_spec,posr,posc,modified_sites_index(modified_sites_count)];               
            end
        end
        for k=1:n_nextnear
            posr = next_nearest(k,1);
            posc = next_nearest(k,2);
            neigh_spec = next_nearest(k,3);
            [lia,locb] =   ismember([neigh_spec,posr,posc],updated_sites_list(:,[1,4,5]),'rows');
            
            if (any(lia))% && (posr~=pos(1) && posc~=pos(2)))
                updated_sites_list(locb,3) = updated_sites_list(locb,3)+1;
                modified_sites_count = modified_sites_count + 1;
                modified_sites_index(modified_sites_count) = locb;
            else
               updated_sites_list(end+1,:) = [neigh_spec, 0, 1, posr, posc, C_new(posr,posc,neigh_spec)==neigh_spec  , 0 , 0]; 
               modified_sites_count = modified_sites_count + 1;
               modified_sites_index(modified_sites_count) = size(updated_sites_list,1);
               added = [added;neigh_spec,posr,posc,modified_sites_index(modified_sites_count)];                
            end
        end
        
        % Reduce the number of neighbors for the neighbors of the source
        % site, except if it is the destination site
        % Reduce the number of next-nearest neighbors for the next-nearest
        % neighbors of the source site, except for the destination site
        [nearest,next_nearest]=neighbor_sites(pos(1),pos(2),sublattice);
        [n_near, ~] = size(nearest);
        [n_nextnear, ~] = size(next_nearest);
        for k=1:n_near
            posr = nearest(k,1); 
            posc = nearest(k,2);
            neigh_spec = nearest(k,3);
            [lia,locb] =   ismember([neigh_spec,posr,posc],updated_sites_list(:,[1,4,5]),'rows');
            
            if (any(lia))
                updated_sites_list(locb,2) = updated_sites_list(locb,2)-1;
                modified_sites_count = modified_sites_count + 1;
                modified_sites_index(modified_sites_count) = locb;
            end
        end
        for k=1:n_nextnear
            posr = next_nearest(k,1);
            posc = next_nearest(k,2);
            neigh_spec = next_nearest(k,3);
            [lia,locb] =   ismember([neigh_spec,posr,posc],updated_sites_list(:,[1,4,5]),'rows');
            
            if (any(lia))
                found_neighbor = 1;
                updated_sites_list(locb,3) = updated_sites_list(locb,3)-1;
                modified_sites_count = modified_sites_count + 1;
                modified_sites_index(modified_sites_count) = locb;
            end
        end
        
        % correct the number of neighbors and next-nearest neighbors of the
        % added sites
        for i=1:size(added,1)
            updated_num_neighbors_added = count_vacancy_neighbors(added(i,2),added(i,3),added(i,1),C_new); % of the added site
            updated_sites_list(added(i,4),2:3) = updated_num_neighbors_added;
        end
        
        % Update the number of neighbors and next-nearest neighbors of new site
        % The columns are 1. species, 2. num of neighbors, 3. num of next
        % neighbors, 4. row-position (y), 5. column-position (x), 6. type of site (0-unoccupied,1-occupied).       
        [~,locb_source] =   ismember([sublattice,pos(1),pos(2)],sites_list(:,[1,4,5]),'rows');
        updated_sites_list(locb_source,1:6) = [sublattice, updated_num_neighbors_source(1), updated_num_neighbors_source(2), pos(1), pos(2),0];
        modified_sites_count = modified_sites_count + 1;
        modified_sites_index(modified_sites_count) = locb_source;
        
        [lia_destination,locb_destination] =   ismember([selected_site(3),selected_site(1),selected_site(2)],sites_list(:,[1,4,5]),'rows');
        if (any(lia_destination))
           updated_sites_list(locb_destination,1:6) = [selected_site(3), updated_num_neighbors_destination(1), updated_num_neighbors_destination(2), selected_site(1), selected_site(2),1];
        else
           updated_sites_list(end+1,1:6) =            [selected_site(3), updated_num_neighbors_destination(1), updated_num_neighbors_destination(2), selected_site(1), selected_site(2),1];
           locb_destination = size(updated_sites_list,1);
        end
        modified_sites_count = modified_sites_count + 1;
        modified_sites_index(modified_sites_count) = locb_destination;

    else
        % Reject the selected diffusion site as it is not leading to
        % betterment of the pore
        updated_sites_list=sites_list;
        C_new=C;
        modified_sites_count = 0;
    end
    
% We do not allow for desorption of vacancy (i.e., addition of an atom from the vapor), but we can surely have addition of vacancy
elseif (mode == 0)  % i.e. addition of vacancy
    
    superbasin_sites = [];

    site_type = atomtype2num(atomtype(pos(1),pos(2),sublattice,C));
    count_etches(site_type) = count_etches(site_type)+1;
    
    % Update the C matrix
    C_new = C;
    C_new(pos(1),pos(2),sublattice) = sublattice;
    
    % get nearest and next nearest neighbors of the new vacancy
    [nearest,next_nearest]=neighbor_sites(pos(1),pos(2),sublattice);
    
    % update the number of neighbors of nearest sites  
    added = [];
    [n_near, ~] = size(nearest);
    [n_nextnear, ~] = size(next_nearest);
    for k=1:n_near
        posr = nearest(k,1);
        posc = nearest(k,2);
        neigh_spec = nearest(k,3);
        [lia,locb] =   ismember([neigh_spec,posr,posc],updated_sites_list(:,[1,4,5]),'rows');
        
        if (any(lia))% && (posr~=pos(1) && posc~=pos(2)))
            updated_sites_list(locb,2) = updated_sites_list(locb,2)+1;
            modified_sites_count = modified_sites_count + 1;
            modified_sites_index(modified_sites_count) = locb;
        else
            updated_sites_list(end+1,:) = [neigh_spec, 1 , 0, posr, posc, C_new(posr,posc,neigh_spec)==neigh_spec  , 0 , 0];
            modified_sites_count = modified_sites_count + 1;
            modified_sites_index(modified_sites_count) = size(updated_sites_list,1);
            added = [added;neigh_spec,posr,posc,modified_sites_index(modified_sites_count)];
        end
    end
       
    % update the number of next-nearest neighbors of next-nearest sites
    for k=1:n_nextnear
        posr = next_nearest(k,1);
        posc = next_nearest(k,2);
        neigh_spec = next_nearest(k,3);
        [lia,locb] =   ismember([neigh_spec,posr,posc],updated_sites_list(:,[1,4,5]),'rows');
        
        if (any(lia))% && (posr~=pos(1) && posc~=pos(2)))
            updated_sites_list(locb,3) = updated_sites_list(locb,3)+1;
            modified_sites_count = modified_sites_count + 1;
            modified_sites_index(modified_sites_count) = locb;
        else
            updated_sites_list(end+1,:) = [neigh_spec, 0, 1, posr, posc, C_new(posr,posc,neigh_spec)==neigh_spec  , 0 , 0];
            modified_sites_count = modified_sites_count + 1;
            modified_sites_index(modified_sites_count) = size(updated_sites_list,1);
            added = [added;neigh_spec,posr,posc,modified_sites_index(modified_sites_count)];
        end
    end
    
    % correct the number of neighbors and next-nearest neighbors of the
    % added sites
    for i=1:size(added,1)
        updated_num_neighbors_added = count_vacancy_neighbors(added(i,2),added(i,3),added(i,1),C_new); % of the added site
        updated_sites_list(added(i,4),2:3) = updated_num_neighbors_added;
    end
        
    % update the number of neighbors and next-neighbors of the new vacancy
    % site
    updated_num_neighbors_source = count_vacancy_neighbors(pos(1),pos(2),sublattice,C_new); % of the source site 
    [~,locb_source] =   ismember([sublattice,pos(1),pos(2)],sites_list(:,[1,4,5]),'rows');
    updated_sites_list(locb_source,1:6) = [sublattice, updated_num_neighbors_source(1), updated_num_neighbors_source(2), pos(1), pos(2),1];
    modified_sites_count = modified_sites_count + 1;
    modified_sites_index(modified_sites_count) = locb_source;
end

% Making the program check every site, instead of just the updated ones
modified_sites_count = size(updated_sites_list,1);
modified_sites_index = 1:modified_sites_count;

for i=1:modified_sites_count
    index = modified_sites_index(i);
    curr_species = updated_sites_list(index,1);
    num_neigh = updated_sites_list(index,2);
    num_next_neigh = updated_sites_list(index,3);
    curr_posr = updated_sites_list(index,4);
    curr_posc = updated_sites_list(index,5);
    site_status = updated_sites_list(index,6);
    
    if (site_status==0)  % empty site
        if (num_neigh>0)
            is_superbasin_site = any(ismembercus([curr_posr,curr_posc,curr_species],superbasin_sites,'rows'));
            Ea =  barrier(atomtype(curr_posr, curr_posc, curr_species, C_new),is_superbasin_site);
            updated_sites_list(index,7)=nuperp_etch(curr_species)*exp(-Ea/(kb_eV*T));  % Rate constant for addition
        else
            updated_sites_list(index,7)=0;
        end
        updated_sites_list(index,8)=0;
        
    else  % filled sites cannot be removed, but can diffuse around    
        updated_sites_list(index,7) = 0;
        
        [num,possible_sites]=diffusion_sites(curr_posr, curr_posc,curr_species,C_new);
        r = zeros(num,1);
        for j=1:num
            
            C_try = C_new;
            C_try(curr_posr, curr_posc, curr_species) = 0;  % remove atom from source site
            C_try(possible_sites(j,1),possible_sites(j,2),possible_sites(j,3)) = possible_sites(j,3);  % add atom to destination site
            
            possible_site_neighbors = count_vacancy_neighbors(possible_sites(j,1),possible_sites(j,2),possible_sites(j,3),C_try);
            possible_site_num_neigh = possible_site_neighbors(1);           

            if (C(possible_sites(j,1),possible_sites(j,2),possible_sites(j,3))==0 &&  possible_site_num_neigh>0)   % if diffusion site is empty and it has atleast one neighbor
                           
                % Remember that the barriers are listed from the initial
                % and final perspectives of a filled atom, and hence the
                % entries for source and destination are interchanged
                destination_site = atomtype(curr_posr,curr_posc,curr_species,C_try);
                source_site = atomtype(possible_sites(j,1),possible_sites(j,2),possible_sites(j,3),C_new);
                is_destination_superbasin_site = any(ismembercus([curr_posr,curr_posc,curr_species],superbasin_sites,'rows'));
                is_source_superbasin_site = any(ismembercus([possible_sites(j,1),possible_sites(j,2),possible_sites(j,3)],superbasin_sites,'rows'));
                Ea = barrier(source_site,destination_site , possible_sites(j,4), is_source_superbasin_site, is_destination_superbasin_site); 
                r(j) = nuperp(curr_species)*exp(-Ea/(kb_eV*T));  % Rate constant for diffusion
           
            else
                r(j)  = 0; % if adjacent diffusion site is filled 
            end
        end
        rtot = sum(r);
        updated_sites_list(index,8)=rtot;       
    end       
end
for m=1:size(updated_sites_list,1)
    [num, possible_sites] = diffusion_sites(updated_sites_list(m,4), updated_sites_list(m,5), updated_sites_list(m,1),C_new);
    num_sites(m,1) = num;
end
end

function [num,possible_sites]=diffusion_sites(i,j,species,C)
[nearest,next_nearest,third_nearest]=neighbor_sites_beyond(i,j,species);
[n_near, ~] = size(nearest);
[n_nextnear, ~] = size(next_nearest);
[n_next_to_nextnear,~] = size(third_nearest);

possible_sites = zeros(n_near+n_nextnear+n_next_to_nextnear,4);
num=0;

for k=1:n_near
    posr = nearest(k,1);
    posc = nearest(k,2);
    neigh_spec = nearest(k,3);
         
    if (C(posr,posc,neigh_spec)==0)
        C_try = C;
        C_try(i,j,species) = 0;  % remove atom from source site
        C_try(posr,posc,neigh_spec)=neigh_spec;       % add atom to destination site
        
        % compute number of neighbors of the source and destination sites
        updated_num_neighbors_source = count_vacancy_neighbors(i,j,species,C_try); % of the source site
        updated_num_neighbors_destination = count_vacancy_neighbors(posr,posc,neigh_spec,C_try); % of the destination site
        
        allow_diffusion  = 1;
        if (count_pores(C_try)==1) %updated_num_neighbors_destination(1) > 0 && updated_num_neighbors_source(1)<3)
            allow_diffusion=1;
        end
        
        if (allow_diffusion==1)
            num=num+1;
            possible_sites(num,:) = [posr,posc,neigh_spec,1];
        end
    end
end

for k=1:n_nextnear
    posr = next_nearest(k,1);
    posc = next_nearest(k,2);
    next_neigh_spec = next_nearest(k,3);
    
    if (C(posr,posc,next_neigh_spec)==0)
        C_try = C;
        C_try(i,j,species) = 0;  % remove atom from source site
        C_try(posr,posc,next_neigh_spec)=next_neigh_spec;       % add atom to destination site
               
        % compute number of neighbors of the source and destination sites
        updated_num_neighbors_source = count_vacancy_neighbors(i,j,species,C_try); % of the source site
        updated_num_neighbors_destination = count_vacancy_neighbors(posr,posc,next_neigh_spec,C_try); % of the destination site
        
        allow_diffusion  = 1;
        if (count_pores(C_try)==1) %updated_num_neighbors_destination(1) > 0 && updated_num_neighbors_source(1)<3)
            allow_diffusion=1;
        end
        
        if (allow_diffusion==1)
            num=num+1;
            possible_sites(num,:) = [posr,posc,next_neigh_spec,2];
        end
    end
end

for k=1:n_next_to_nextnear
    posr = third_nearest(k,1);
    posc = third_nearest(k,2);
    next_to_nextneigh_spec = third_nearest(k,3);
    
    if (C(posr,posc,next_to_nextneigh_spec)==0)
        C_try = C;
        C_try(i,j,species) = 0;  % remove atom from source site
        C_try(posr,posc,next_to_nextneigh_spec)=next_to_nextneigh_spec;       % add atom to destination site
        
        % compute number of neighbors of the source and destination sites
        updated_num_neighbors_source = count_vacancy_neighbors(i,j,species,C_try); % of the source site
        updated_num_neighbors_destination = count_vacancy_neighbors(posr,posc,next_to_nextneigh_spec,C_try); % of the destination site
        
        allow_diffusion  = 1;
        if (count_pores(C_try)==1) %updated_num_neighbors_destination(1) > 0 && updated_num_neighbors_source(1)<3)
            allow_diffusion=1;
        end
        
        if (allow_diffusion==1)
            num=num+1;
            possible_sites(num,:) = [posr,posc,next_to_nextneigh_spec,3];
        end
    end
end
end

function y=pbc(x)
% This function implements periodic boundary conditions for lattice
% indices

% Input: row/column position before wrapping
% Output: row/column position after periodic boundary wrapping

global N;
if (x==0 || x==N)
    y=N;
else
    y=mod(x,N);
end
end

function [nearest,next_nearest]=neighbor_sites(i,j,sublattice)
% This function returns the nearest and next-nearest lattice sites for a
% given site (i, j, sublattice)

% Inputs: row (i), column (j), and sublattice (1 or 2)
% Outputs: list of nearest neighbor and next-nearest neighbor lattice sites
% in the format (row, column, sublattice)

if (sublattice==1)
    
    nearest = [pbc(i-1),    pbc(j-1+(~mod(i,2))),        2;
               pbc(i-1),    pbc(j+(~mod(i,2))),          2;
               i,           j,                           2];
    
    next_nearest = [pbc(i-1),    pbc(j-1+(~mod(i,2))),        sublattice;
        pbc(i-1),                pbc(j+(~mod(i,2))),          sublattice;
        i,                       pbc(j-1),                    sublattice;
        i,                       pbc(j+1),                    sublattice;
        pbc(i+1),                pbc(j-1+(~mod(i,2))),        sublattice;
        pbc(i+1),                pbc(j+(~mod(i,2))),          sublattice];
    
elseif (sublattice==2)

    nearest = [i,       j,                       1;
        pbc(i+1),       pbc(j+(~mod(i,2))),      1;
        pbc(i+1),       pbc(j-1+(~mod(i,2))),    1];
    
    next_nearest = [pbc(i-1),       pbc(j-1+(~mod(i,2))),        sublattice;
        pbc(i-1),                   pbc(j+(~mod(i,2))),          sublattice;
        i,                          pbc(j-1),                    sublattice;
        i,                          pbc(j+1),                    sublattice;
        pbc(i+1),                   pbc(j-1+(~mod(i,2))),        sublattice;
        pbc(i+1),                   pbc(j+(~mod(i,2))),          sublattice];
end
end

function [nearest,next_nearest,third_nearest]=neighbor_sites_beyond(i,j,sublattice)
% This function returns the nearest, next-nearest, and next-to-next-nearest
% lattice sites for a given site (i, j, sublattice)

% Inputs: row (i), column (j), and sublattice (1 or 2)
% Outputs: list of nearest neighbor, next-nearest neighbor, and 
% next-to-next-nearest neighbor lattice sites % in the format (row, column,
% sublattice)

if (sublattice==1)
    
    nearest = [pbc(i-1),    pbc(j-1+(~mod(i,2))),        2;
               pbc(i-1),    pbc(j+(~mod(i,2))),          2;
               i,           j,                           2];
    
    next_nearest = [pbc(i-1),    pbc(j-1+(~mod(i,2))),        sublattice;
        pbc(i-1),                pbc(j+(~mod(i,2))),          sublattice;
        i,                       pbc(j-1),                    sublattice;
        i,                       pbc(j+1),                    sublattice;
        pbc(i+1),                pbc(j-1+(~mod(i,2))),        sublattice;
        pbc(i+1),                pbc(j+(~mod(i,2))),          sublattice];
    
    third_nearest = [pbc(i-2)    ,   pbc(j-1)              ,  2
        pbc(i-2)    ,   pbc(j)                ,  2;
        pbc(i-2)    ,   pbc(j+1)              ,  2;
        pbc(i-1)    ,   pbc(j+1+(~mod(i,2)))  ,  2;
        pbc(i)      ,   pbc(j+1)              ,  2;
        pbc(i+1)    ,   pbc(j+(~mod(i,2)))    ,  2;
        pbc(i+1)    ,   pbc(j-1+(~mod(i,2)))  ,  2;
        pbc(i)      ,   pbc(j-1)              ,  2;
        pbc(i-1)    ,   pbc(j-2+(~mod(i,2)))  ,  2];
    
elseif (sublattice==2)

    nearest = [i,       j,                       1;
        pbc(i+1),       pbc(j+(~mod(i,2))),      1;
        pbc(i+1),       pbc(j-1+(~mod(i,2))),    1];
    
    next_nearest = [pbc(i-1),       pbc(j-1+(~mod(i,2))),        sublattice;
        pbc(i-1),                   pbc(j+(~mod(i,2))),          sublattice;
        i,                          pbc(j-1),                    sublattice;
        i,                          pbc(j+1),                    sublattice;
        pbc(i+1),                   pbc(j-1+(~mod(i,2))),        sublattice;
        pbc(i+1),                   pbc(j+(~mod(i,2))),          sublattice];
    
    third_nearest  = [pbc(i-1),  pbc(j-1+(~mod(i,2))),  1;
        pbc(i-1),  pbc(j+(~mod(i,2))),    1;
        pbc(i),    pbc(j+1)            ,  1;
        pbc(i+1),  pbc(j+1+(~mod(i,2))),  1;
        pbc(i+2),  pbc(j+1)           ,   1;
        pbc(i+2),  pbc(j)           ,     1;
        pbc(i+2),  pbc(j-1)           ,   1;
        pbc(i+1),  pbc(j-2+(~mod(i,2))),  1;
        pbc(i),    pbc(j-1)            ,  1];
end
end

function y=count_vacancy_neighbors(i,j,sublattice,C)
% This function counts the number of vacancy neighbors to a given lattice
% site

% Input: row (i), column (j), sublattice (1 or 2), state of the system (C)
% Output: number of vacancy neighbors in the format (nearest, next-nearest)

[nearest,next_nearest]=neighbor_sites(i,j,sublattice);
[n_near, ~] = size(nearest);
[n_nextnear, ~] = size(next_nearest);

y=[0;0];
for k=1:n_near
    posr = nearest(k,1);
    posc = nearest(k,2);
    neigh_sublattice = nearest(k,3);
    if (C(posr,posc,neigh_sublattice)==neigh_sublattice) % site contains vacancy
        y(1)=y(1)+1;
    end
end
for k=1:n_nextnear
    posr = next_nearest(k,1);
    posc = next_nearest(k,2);
    next_neigh_sublattice = next_nearest(k,3);
    if (C(posr,posc,next_neigh_sublattice)==next_neigh_sublattice) % site contains vacancy
        y(2)=y(2)+1;
    end
end
end

function [y,complete_members,incomplete_members]=is_complete_hexagon(i,j,C)
% This function determines if the (i,j)th hexagon in the graphene lattice
% has all 6 possible atomic sites filled with carbon atoms. This function
% is useful to determine the rim atoms of a nanopore, because the rim atoms
% are all those atoms which belong to an "incomplete" hexagon

% Inputs: row (i), column (j), state of the system (C)
% Output: flag (y), list of filled and vacant members of the hexagon

row_mod_two = ~mod(i,2);

count_filled=0;
incomplete_members=[];
complete_members  =[];

% We will sequentially go through all 6 possible lattice sites in the
% (i,j)th hexagon and determine the filled and unfilled lattice sites

if (C(i,j,1)==0)   % 0 means vacancy not added so site is filled
    count_filled=count_filled+1;
    complete_members = [complete_members; i  j 1];
else
    incomplete_members = [incomplete_members; i  j 1];
end
if (C(i,j,2)==0)   % 0 means vacancy not added so site is filled
    count_filled=count_filled+1;
    complete_members = [complete_members; i  j 2];
else
    incomplete_members = [incomplete_members; i  j 2];
end
if (C(pbc(i-1),pbc(j+row_mod_two),2)==0)   % 0 means vacancy not added so site is filled
    count_filled=count_filled+1;
    complete_members = [complete_members; pbc(i-1),pbc(j+row_mod_two),2];
else
    incomplete_members = [incomplete_members; pbc(i-1),pbc(j+row_mod_two),2];
end
if (C(pbc(i+1),pbc(j+row_mod_two),1)==0)   % 0 means vacancy not added so site is filled
    count_filled=count_filled+1;
    complete_members = [complete_members; pbc(i+1),pbc(j+row_mod_two),1]; 
else
    incomplete_members = [incomplete_members; pbc(i+1),pbc(j+row_mod_two),1]; 
end
if (C(i,pbc(j+1),1)==0)   % 0 means vacancy not added so site is filled
    count_filled=count_filled+1;
    complete_members = [complete_members; i,pbc(j+1),1];     
else
    incomplete_members = [incomplete_members; i,pbc(j+1),1];     
end
if (C(i,pbc(j+1),2)==0)   % 0 means vacancy not added so site is filled
    count_filled=count_filled+1;
    complete_members = [complete_members; i,pbc(j+1),2];         
else
    incomplete_members = [incomplete_members; i,pbc(j+1),2];         
end

y=0;  % assuming hexagon is NOT complete
if (count_filled==6)
    y=1;  % hexagon is complete
end
end

function y=hexagon_neighbors(i,j,sublattice)

hex_to_the_right = [   i,       j,                (sublattice==1)*2 + (sublattice==2)*1
                      pbc(i-1), pbc(j+~mod(i,2)),  2
                       i,       pbc(j+1),          1
                       i,       pbc(j+1),          2
                      pbc(i+1), pbc(j+~mod(i,2)),  1];
                  
hex_to_the_left = [   i,        j,                 (sublattice==1)*2 + (sublattice==2)*1
                      pbc(i-1), pbc(j-1+~mod(i,2)), 2
                      i,        pbc(j-1),           1
                      i,        pbc(j-1),           2
                      pbc(i+1), pbc(j-1+~mod(i,2)), 1];

hex_to_the_top = [    pbc(i-1), pbc(j-1+~mod(i,2)), 2
                      pbc(i-1), pbc(j-1+~mod(i,2)), 1
                      pbc(i-2), j,                  2
                      pbc(i-1), pbc(j+~mod(i,2)),   1
                      pbc(i-1), pbc(j+~mod(i,2)),   2];

hex_to_the_bottom = [ pbc(i+1), pbc(j+~mod(i,2)),   1
                      pbc(i+1), pbc(j+~mod(i,2)),   2
                      pbc(i+2), j,                  1
                      pbc(i+1), pbc(j-1+~mod(i,2)), 2
                      pbc(i+1), pbc(j-1+~mod(i,2)), 1];
if (sublattice==1)
    y = [hex_to_the_right;
         hex_to_the_left;
         hex_to_the_top];
elseif (sublattice==2)
    y = [hex_to_the_right;
         hex_to_the_left;
         hex_to_the_bottom];
end
y = unique(y,'rows');
end
     
function [initial_vacancy,C0,num_atoms]=init(N)
% This function initializes the state of the system

% Input: size of the lattice (N). Number of lattice sites is: N*N*2, the
% factor of 2 due to there being two sublattices in graphene
% Output: initial state of the system (initial_vacancy and C0), number of
% vacancies in each sublattice

% Initialize system with monoatom vacancy, can suitably tweak as required
C0 = zeros(N,N,3);
center=round(N/2);

initial_vacancy = [      1       center  center];
C0(center,center,1) = 1;
num_atoms = [1,0];
end

function [num_pores,size_each_pore,node_details]=count_pores(C)
% To decide if a set of vacancies correspond to a single pore, we look at
% the adjacency matrix, wherein two vacancies are neighbors if they belong
% to the same hexagon in graphene (note that the conventional definition of
% neighboring vacancy sites being connected by a "bond" is not used here) 
N = size(C,1);

num_removed_atoms = sum(sum(sum(C>0)));
node_details=zeros(num_removed_atoms,4);
count_removed_atoms=0;
for i=1:N
    for j=1:N       
        for k=1:2
            if (C(i,j,k)>0)
                count_removed_atoms = count_removed_atoms+1;
                C(i,j,k) = count_removed_atoms;
                
                node_details(count_removed_atoms,1)=i;
                node_details(count_removed_atoms,2)=j;
                node_details(count_removed_atoms,3)=k;
            end
        end
    end
end

adjmat_removed = zeros(count_removed_atoms, count_removed_atoms);
for i=1:N
    for j=1:N
        for k=1:2
            curr_index = C(i,j,k);
            if (curr_index>0)  % current atom is removed atom
                
                hex_neigh = hexagon_neighbors(i,j,k);
                n_hex_neigh = size(hex_neigh,1);
                
                for l=1:n_hex_neigh
                    posr = hex_neigh(l,1);
                    posc = hex_neigh(l,2);
                    neigh_spec = hex_neigh(l,3);
                    
                    neigh_index = C(posr,posc,neigh_spec);
                    if (neigh_index>0)  % atom in the same hexagon is removed atom
                        adjmat_removed(curr_index,neigh_index)=1;
                        adjmat_removed(neigh_index,curr_index)=1;
                    end
                end
            end
        end
    end
end
[num_pores, which_component_node_belongs_to] = graphconncomp(sparse(adjmat_removed),'Directed',false, 'Weak', true);

% C is a vector indicating to which component each node belongs
node_details(:,4) = which_component_node_belongs_to';
size_each_pore = zeros(num_pores,1);
for i=1:num_pores
    size_each_pore(i) = sum(which_component_node_belongs_to==i);
end
end
 
function [num_AC_atoms,num_ZZ_atoms,num_UA_atoms,num_5_membered_ring]=visualize(C,N,count)
% This function saves the current state of the system, given the surface
% state matrix, as an XYZ file and a directed adjacency matrix for the
% antimolecule of the formed nanopore

% Inputs: Surface state matrix (C), Size of grid (N), and current timestep (count)
% Outputs: number of armchair, zigzag, and unassigned atoms, and 5-membered
% rings

global lattice_const isomer_index pore_size basedirectory;

% Define box sizes in the x and y directions
box_x = N*2*(lattice_const/sqrt(3))*cos(pi/6);
box_y = N*(lattice_const/sqrt(3))*(2+2*cos(pi/3))/2;

if (count<pore_size && count>=3)
    
else
    fid_xyz = fopen(      [basedirectory,'/pore',num2str(pore_size),'/pore',            num2str(isomer_index)  ,'.xyz'],'wt');
    fid_txt = fopen(      [basedirectory,'/pore',num2str(pore_size),'/pore',            num2str(isomer_index)  ,'.txt'],'wt');
    fid_adjmat_removed =  [basedirectory,'/pore',num2str(pore_size),'/adjmat_antimol_', num2str(isomer_index)  ,'.txt'];
end

% Create variables for number of total atoms, number of remaining atoms,
% and number of removed atoms
[m,n]=size(C(:,:,1));
num_total_atoms = m*n*2;
num_remaining_atoms = sum(sum(C(:,:,1)==0))+sum(sum(C(:,:,2)==0));
num_removed_atoms = m*n*2-num_remaining_atoms;

rim_atoms=[];
% rim_atoms contains four columns: (row, col, sublattice, edge_atom_type)
% fourth column in rim_atoms is the type of edge site (0- not an edge site,
% 1- AC site, 2-ZZ site, 3-UA site)
num_5_membered_ring=0;
for i=1:N
    for j=1:N        
            [y,complete_members,~]=is_complete_hexagon(i,j,C);
            if (y==0) % if hexagon is incomplete
                rim_atoms = [rim_atoms; complete_members];
                
                if (size(complete_members,1)==5)
                    num_5_membered_ring = num_5_membered_ring + 1;
                end
            end
   end
end

% remove duplicate entries from the rim atoms
rim_atoms = unique(rim_atoms,'rows');
num_rim_atoms = size(rim_atoms,1);

% add extra column at end to indicate kind of edge 
rim_atoms(:,end+1) = zeros(num_rim_atoms,1);

% count the number of AC, ZZ, and UA edge atoms
num_AC_atoms=0;
num_ZZ_atoms=0;
num_UA_atoms=0;

for i=1:num_rim_atoms
    % rim atom is itself a filled site, so no vacancy, i.e. C(row_pos,col_pos,species)=0
    row_pos = rim_atoms(i,1);
    col_pos = rim_atoms(i,2);
    species = rim_atoms(i,3);
    
    [nearest,~]=neighbor_sites(row_pos,col_pos,species);    
    num_neigh_neigh = 0;
    for j=1:size(nearest,1)
        if (C(nearest(j,1), nearest(j,2), nearest(j,3)) == 0) % the neighboring site is occupied, so count neighbors of neighbors to assign ZZ, AC, etc.
            neigh_neigh = count_vacancy_neighbors(nearest(j,1), nearest(j,2), nearest(j,3),C);
            num_neigh_neigh = num_neigh_neigh + neigh_neigh(1);         
        end
    end
    
    y = count_vacancy_neighbors(row_pos, col_pos, species,C);
    num_neighbors = y(1);
    
    % num_neighbors can be (a) 0, i.e. fully filled coordination sphere
    % (cannot be rim atom), (b) 1, i.e. 1 neighbor empty, (c) 2, i.e. 2
    % neighbors empty, (d) 3, i.e. no neighbors (cannot be rim atom)
    % 1 and 2 are plausible options
    
    if (num_neighbors > 0) % this means atleast one neighbor is a vacancy so it is an edge atom (not all rim atoms are edge atoms)               
        switch num_neigh_neigh
            case 0
                % this is a zigzag (ZZ) site
                num_ZZ_atoms = num_ZZ_atoms+1;
                rim_atoms(i,4) = 2;               
                
            case 1
                % this is an armchair (AC) site
                num_AC_atoms = num_AC_atoms+1;
                rim_atoms(i,4) = 1;
                
            case 2
                % this is an unassigned (UA) site
                num_UA_atoms = num_UA_atoms+1;
                rim_atoms(i,4) = 3;
        end
    end
end
num_edge_atoms = num_AC_atoms + num_ZZ_atoms + num_UA_atoms;

fprintf(fid_xyz,[num2str(num_total_atoms+num_rim_atoms+num_edge_atoms),'\n # of carbons removed: ', num2str(num_removed_atoms)  ...
            ,'; # of rim atoms: ', num2str(num_rim_atoms),  '; # of edge atoms: ', num2str(num_edge_atoms),  ' \n']);      
fprintf(fid_txt,'%d\n',N);

count_removed_atoms=0;
count_rim_atoms=0;
% Cycle through all locations on the 2D lattice
for i=1:N
    for j=1:N       
        if (C(i,j,1)==0)
            % atom is present, i.e. vacancy not added
            coord_sublattice_1 = [ ~mod(i,2)*(lattice_const/2) + (j-1)*lattice_const, (i-1)*sqrt(3)*lattice_const/2  ];
            fprintf(fid_xyz,['C ', num2str(coord_sublattice_1(1)), '  ', num2str(coord_sublattice_1(2)), '   0.000\n']);           
            
            % check if rim atom, and assign index
            [lia,~] = ismember([i,j,1],rim_atoms(:,1:3),'rows');
            if (any(lia))
               count_rim_atoms = count_rim_atoms+1;
               C(i,j,1) = num_removed_atoms + count_rim_atoms;
            end
        else
            % print removed atoms with different name
            coord_sublattice_1 = [ ~mod(i,2)*(lattice_const/2) + (j-1)*lattice_const, (i-1)*sqrt(3)*lattice_const/2  ];           
            fprintf(fid_xyz,['Re ', num2str(coord_sublattice_1(1)), '  ', num2str(coord_sublattice_1(2)), '   0.000\n']);
          
            % assign index to removed atoms
            fprintf(fid_txt,'%d %d %d\n',1,i,j);
            count_removed_atoms = count_removed_atoms+1;
            C(i,j,1) = count_removed_atoms;
        end
        
        if (C(i,j,2)==0)
            % atom is present, i.e. vacancy not added
            coord_sublattice_2 = [ ~mod(i,2)*(lattice_const/2) + (j-1)*lattice_const, (i-1)*sqrt(3)*lattice_const/2 + (lattice_const/sqrt(3)) ];
            fprintf(fid_xyz,['C ', num2str(coord_sublattice_2(1)), '  ', num2str(coord_sublattice_2(2)), '   0.000\n']);
            
            % check if rim atom, and assign index
            [lia,~] = ismember([i,j,2],rim_atoms(:,1:3),'rows');
            if (any(lia))
               count_rim_atoms = count_rim_atoms+1;
               C(i,j,2) = num_removed_atoms + count_rim_atoms;
            end
        else
            % print removed atoms with different name
            coord_sublattice_2 = [ ~mod(i,2)*(lattice_const/2) + (j-1)*lattice_const, (i-1)*sqrt(3)*lattice_const/2 + (lattice_const/sqrt(3)) ];
            fprintf(fid_xyz,['Re ', num2str(coord_sublattice_2(1)), '  ', num2str(coord_sublattice_2(2)), '   0.000\n']);
            
            % assign index to removed atoms
            fprintf(fid_txt,'%d %d %d\n',2,i,j);
            count_removed_atoms = count_removed_atoms+1;
            C(i,j,2) = count_removed_atoms;
        end
    end
end
rim_atoms_sublattice_1 = rim_atoms(rim_atoms(:,3)==1,:);
rim_atoms_sublattice_2 = rim_atoms(rim_atoms(:,3)==2,:);

% print rim atoms with different element name
coord_rim_sublattice_1 = [ ~mod(rim_atoms_sublattice_1(:,1),2)*(lattice_const/2) + (rim_atoms_sublattice_1(:,2)-1)*lattice_const, (rim_atoms_sublattice_1(:,1)-1)*sqrt(3)*lattice_const/2  ];
for i=1:size(coord_rim_sublattice_1,1)
    fprintf(fid_xyz,['Ri ', num2str(coord_rim_sublattice_1(i,1)), '  ', num2str(coord_rim_sublattice_1(i,2)), '   0.000\n']);

    label = {'AC','ZZ','UA'};
    which_edge_atom = rim_atoms_sublattice_1(i,4);
    if (which_edge_atom)
        fprintf(fid_xyz,[label{which_edge_atom},' ', num2str(coord_rim_sublattice_1(i,1)), '  ', num2str(coord_rim_sublattice_1(i,2)), '   0.000\n']);   
    end
end

coord_rim_sublattice_2 = [ ~mod(rim_atoms_sublattice_2(:,1),2)*(lattice_const/2) + (rim_atoms_sublattice_2(:,2)-1)*lattice_const, (rim_atoms_sublattice_2(:,1)-1)*sqrt(3)*lattice_const/2 + (lattice_const/sqrt(3)) ];
for i=1:size(coord_rim_sublattice_2,1)
    fprintf(fid_xyz,['Ri ', num2str(coord_rim_sublattice_2(i,1)), '  ', num2str(coord_rim_sublattice_2(i,2)), '   0.000\n']);

    label = {'AC','ZZ','UA'};
    which_edge_atom = rim_atoms_sublattice_2(i,4);
    if (which_edge_atom)
        fprintf(fid_xyz,[label{which_edge_atom},' ', num2str(coord_rim_sublattice_2(i,1)), '  ', num2str(coord_rim_sublattice_2(i,2)), '   0.000\n']);   
    end
end

% make adjacency matrix out of removed atoms, including bond orientations
adjmat_removed = zeros(num_removed_atoms, num_removed_atoms);
for i=1:size(C,1)
    for j=1:size(C,2)
        curr_index_1 = C(i,j,1);
        if (curr_index_1>0)  % current atom is removed atom or rim atom
            [nearest,~]=neighbor_sites(i,j,1);
            [n_near, ~] = size(nearest);
            coord_A = [ ~mod(i,2)*(lattice_const/2) + (j-1)*lattice_const, (i-1)*sqrt(3)*lattice_const/2  ];

            for k=1:n_near
                posr = nearest(k,1);
                posc = nearest(k,2);
                neigh_spec = nearest(k,3);
                coord_B = [ ~mod(posr,2)*(lattice_const/2) + (posc-1)*lattice_const, (posr-1)*sqrt(3)*lattice_const/2 + (lattice_const/sqrt(3)) ];

                direction = coord_B-coord_A; % find direction of bond
                
                if (direction(1)>box_x/2.0) % right
                    direction(1) = -(box_x-direction(1));
                end
                if (direction(1)<-box_x/2.0) % left
                    direction(1) = (box_x+direction(1));                   
                end
                if (direction(2)>box_y/2.0) % top
                    direction(2) = -(box_y-direction(2));
                end
                if (direction(2)<-box_y/2.0) % down
                    direction(2) = (box_y+direction(2));
                end                  
                
                slope_direction = atan(abs(direction(2)/direction(1)))*180/pi;
                
                if (slope_direction > 85 && slope_direction < 95)
                    if (direction(2)>0)
                        adj_direction = 1; % top
                    else
                        adj_direction = 2; % bottom
                    end
                elseif (slope_direction>25 && slope_direction<35)
                    if (direction(2)>0) % top
                        if (direction(1)<0)
                            adj_direction = 3; % top left
                        else
                            adj_direction = 4; % top right
                        end
                    else                % bottom
                        if (direction(1)<0)
                            adj_direction = 5; % bottom left
                        else
                            adj_direction = 6; % bottom right
                        end
                    end
                end
                                   
                neigh_index = C(posr,posc,neigh_spec); 
                if (neigh_index>0) 
                   if (curr_index_1 <= num_removed_atoms && neigh_index <= num_removed_atoms)  % both current and neighboring atoms are removed atom only
                        adjmat_removed(curr_index_1,neigh_index)=adj_direction;
                   end
                end
            end
        end
        
        curr_index_2 = C(i,j,2);
        if (curr_index_2>0) 
            [nearest,~]=neighbor_sites(i,j,2);
            [n_near, ~] = size(nearest);
            coord_B = [ ~mod(i,2)*(lattice_const/2) + (j-1)*lattice_const, (i-1)*sqrt(3)*lattice_const/2 + (lattice_const/sqrt(3)) ];

            for k=1:n_near
                posr = nearest(k,1);
                posc = nearest(k,2);
                neigh_spec = nearest(k,3);
                coord_A = [ ~mod(posr,2)*(lattice_const/2) + (posc-1)*lattice_const, (posr-1)*sqrt(3)*lattice_const/2  ];

                direction = coord_A-coord_B;
                
                if (direction(1)>box_x/2.0) % right
                    direction(1) = -(box_x-direction(1));
                end
                if (direction(1)<-box_x/2.0) % left
                    direction(1) = (box_x+direction(1));                   
                end
                if (direction(2)>box_y/2.0) % top
                    direction(2) = -(box_y-direction(2));
                end
                if (direction(2)<-box_y/2.0) % down
                    direction(2) = (box_y+direction(2));
                end                 
                                                 
                slope_direction = atan(abs(direction(2)/direction(1)))*180/pi;
                
                if (slope_direction > 85 && slope_direction < 95)
                    if (direction(2)>0)
                        adj_direction = 1; % top
                    else
                        adj_direction = 2; % bottom
                    end
                elseif (slope_direction>25 && slope_direction<35)
                    
                    if (direction(2)>0) % top
                        if (direction(1)<0)
                            adj_direction = 3; % top left
                        else
                            adj_direction = 4; % top right
                        end
                    else                % bottom
                        if (direction(1)<0)
                            adj_direction = 5; % bottom left
                        else
                            adj_direction = 6; % bottom right
                        end
                    end
                end
                                
                neigh_index = C(posr,posc,neigh_spec); 
                if (neigh_index>0)   % neighboring atom is removed or rim atom
                    if (curr_index_2 <= num_removed_atoms && neigh_index <= num_removed_atoms)  % both current and neighboring atoms are removed atom only
                        adjmat_removed(curr_index_2,neigh_index)=adj_direction;
                    end
                end
            end
        end
        
    end
end
dlmwrite(fid_adjmat_removed, adjmat_removed);
fclose('all');
end
