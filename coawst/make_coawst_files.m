%% create COAWST set up 
% for Hawke Bay New Zealand
% Ted Conroy

%% Specify files to create (logicals 0 or 1)
make_tide_file = 0; % make a tidal forcing file from txpo atlas 9.2
make_init = 1  ;    % make an initial conditions netcdf file
make_bc_file = 0;   % make open boundary forcing file
make_forcing = 0;   % make atmospheric forcing files

%% model information

% TIME specifics
start_time = datenum(2017,3,1,0,0,0); % start time of run
end_time = datenum(2017,5,31);
num_months = 1;
roms_time_reference = datenum(2010,1,1); % will be seconds since this time
% set reference time
ref_time = (roms_time_reference-datenum('01-Jan-1900'))*86400;

% FILE DIRECTORIES
grid_file = 'sept_2020_v2.nc';
bc_file_path ='/nesi/nobackup/mocean02574/NZB_N50/'; % path to directory of boundary forcing files
input_file_base = '/nesi/nobackup/mocean02574/NZB_N50/nz5km_his_';


% using nested grids? if yes list the grid files in structure nest
nesting = 1; % 0 or 1
nest(1).grid = 'l1_south.nc';

% create list of input files based on date, here for monthly output files
input_files(1).name = '/nesi/nobackup/mocean02574/NZB_N50/nz5km_his_201703.nc';
%input_files(2).name = '/nesi/nobackup/mocean02574/NZB_N50/nz5km_his_201704.nc';
%input_files(3).name = '/nesi/nobackup/mocean02574/NZB_N50/nz5km_his_201705.nc';

% automate this
%for n=1:num_months
 %   input_file(n).name = [input_file_base 
  %      end
        
% OUTPUT FILE NAMES
tide_ncfile=[ datestr(start_time,'ddmmmyyyy') '.nc']; % tide forcing file
init_file_name = 'ocean_ini.nc'; % initial conditions file
bc_file_name = 'ocean_bry_'; % open boundary base
frc_file_name = 'ocean_frc_'; % forcing file base

% add (recursive) paths to m-files needed
addpath(genpath('/scale_wlg_persistent/filesets/home/tc196/COAWST'))
% add roms netcdf tools, which are different than native matlab!
addpath(genpath('/scale_wlg_persistent/filesets/home/tc196/COAWST/Tools/mfiles/rutgers/netcdf/'))
cd nctoolbox-1.1.3
setup_nctoolbox
cd ../


%% -----------------------------------------------------------------------------------
%% -----------------------------------------------------------------------------------
%% END USER DEFS ---------------------------------------------------------------------
%% -----------------------------------------------------------------------------------
%% -----------------------------------------------------------------------------------

%% Print ROMS time information
disp(['start day to enter in ocean.in file: '])
days(start_time-roms_time_reference)
disp(['from the reference time: ' datestr(roms_time_reference)])
disp(['forcing files will end on ' datestr(end_time)])

%% TIDAL FORCING: TXPO tidal model 
% txpo atlas v9.2, has regional domains nested in global
% downloaded txpo2roms from austides, had to ask for txpo netcdf file

if make_tide_file
ROMSnames = {'MM' 'MF' 'Q1' 'O1' 'P1' 'K1' 'N2' 'M2' 'S2' 'K2' 'MN4' 'M4' 'MS4' '2N2' 'S1'};
TPXO2ROMS_v5pt1(start_time,ROMSnames,grid_file,tide_ncfile,days(end_time - start_time))
disp(['created tidal forcing file from TXPO: ' tide_ncfile])
end

%% create initial condition file from existing .his file
% In cppdefs.h you should have #undef ana_initial
% interpolates onto grid, removes NaNs, ...              

if make_init 
    disp('-------------------------------------------------------------------')
    disp (['creating initial conditions file from first input file: ' input_files(1).name ])
    disp('finding desired time to interpolate from')
    disp('check output below to make sure time selected matches start time')
    make_init_file (grid_file,input_files(1).name,init_file_name,start_time,ref_time)
    if nesting % need to make init files for the children
        for n=1:length(nest)
            make_init_file(nest(n).grid,input_files(1).name,['child' num2str(n) init_file_name],start_time,ref_time)
        end
    end
    disp('done initializing')
    disp('-------------------------------------------------------------------')
end

%% boundary forcing

if make_bc_file
    make_bry_files(bc_file_name,grid_file,input_files)
    disp('-------------------------------------------------------------------')
    disp (['created boundary conditions file: ' bc_file_name ])
    disp('-------------------------------------------------------------------')
end

%% Forcing files
% currently sets these = 0
% sustr svstr bhflux

if make_forcing
    disp('-------------------------------------------------------------------')
    make_frc_file(input_files,frc_file_name,grid_file,1,'sustr','svstr','shflux');
    disp('made frc file(s) ------------------------------------------------------')
end


%% nesting files

%% atmospheric forcing
% Hau-moana will be done ~september
% constant forcing for now, or reanalysis 

%% waves

% roms2swan

%% sediment
%% other files
%add_mask add_sponge add_drag
