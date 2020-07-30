%% creating a COAWST set up 
% for hawke bay , new zealand
% ted conroy 2020

%% FILE CREATION (logicals)

make_tide_file = 0; % make a tidal forcing file from txpo atlas 9.2
make_init = 1; % make an initial conditions netcdf file
make_bc_file = 1; % make open boundary forcing file
make_rivers = 0; % create rirver forcing files
make_forcing = 0; % make forcing file (e.g. atmospheric)

% nesting
% sediment
% waves
% clm file

%% TIME specifics

% For ROMS .in file use TIME_REF=-2; means in MJD
% more info about time in roms

start_time = datenum(2017,2,1); % start time of run
mjd_init = datenum(1968,5,23);
tval_for_infile = start_time - mjd_init; % enter this value in the .in file
disp(['enter ' num2str(tval_for_infile) ' for TIDE_START value'])
days_to_simulate=365;  % length of model run in days

%% FILE DIRECTORIES

grid_dir = '/home/tc196/setup/files/';
grid_file_name = 'roms_grid_001_1km.nc';
grid_mat_file = '/home/tc196/setup/files/roms_grid_001_1km.mat'; % output from Gridbuilder
grid_file  = [grid_dir '/' grid_file_name]; % full path of grid file
data_file ='/nesi/nobackup/mocean02574/NZB_N50/nz5km_his_201701.nc'; % file to initialize from (automate this based on input date)
bc_file_path ='/nesi/nobackup/mocean02574/NZB_N50/'; % path to directory of boundary forcing files

%% OUTPUT FILE NAMES

tide_ncfile=['C:\Users\tc196\Dropbox\research\hawkes_bay\model\coawst\input_files\tide\' datestr(start_time,'ddmmmyyyy') '.nc']; % tide forcing file
init_file_name = 'ocean_ini.nc'; % initial conditions file
bc_file_name = 'ocean_bry.nc'; % open boundary file

%% add (recursive) paths to m-files needed

addpath(genpath('scale_wlg_persistent/filesets/home/tc196/'))
cd nctoolbox-1.1.3
setup_nctoolbox
cd ../
%-------------------------------------------------------------------------------------------------
%% end user defs ---------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------

%% grid, bathy, mask
% i used GridBuilder (google: Austides)
% can import coastline and bathymetry into GridBuilder

%% tide from TXPO tidal model 
% txpo atlas v9.2, has regional domains nested in global
% downloaded txpo2roms from austides, had to ask for txpo netcdf file

if make_tide_file
ROMSnames={'MM' 'MF' 'Q1' 'O1' 'P1' 'K1' 'N2' 'M2' 'S2' 'K2' 'MN4' 'M4' 'MS4' '2N2' 'S1'};
TPXO2ROMS_v5pt1(start_time,ROMSnames,grid_file,tide_ncfile,days_to_simulate) 
disp(['created tidal forcing file:' tide_ncfile])
end

%% create initial condition file from existing .his file
% In cppdefs.h you should have #undef ana_initial
% interpolates onto grid, removes NaNs, ...

if make_init
make_init_file (grid_file,data_file,init_file_name,start_time,grid_mat_file)
disp (['created initial conditions file:' init_file_name ])
end

%% boundary forcing
% first find the files to use (can add this to function below)
% allfiles = dir([bc_file_path '\*his']); % currently using history files
% will be like ^ that on nesi, not here though

if make_bc_file
  % boundary forcing file(s)
  % takes awhile, worth looking into speed up..
    input_f = data_file; % temporary
    make_bry_files(bc_file_name,grid_file,input_f,grid_mat_file)
    
   % nudging time scale file
end

%% rivers

if make_rivers
make_riv_files
end

%% nesting files

%% atmospheric forcing
% Hau-moana will be done ~september
% constant forcing for now, or reanalysis 
if make_forcing
    create_roms_forcings
end

%% waves

% roms2swan

%% sediment
%% other files
%add_mask add_sponge add_drag
