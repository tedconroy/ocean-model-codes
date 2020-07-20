%% creating a COAWST set up for hawke bay 
% ted conroy 2020

%% specify what to do (logicals)
make_tide_file = 0; % make a tidal forcing file from txpo atlas 9.2
make_init = 1;
make_bc_files = 0;
make_rivers = 0;
% nesting
% sediment
% waves
% clm file

%% specify time interval

% For ROMS .in file use TIME_REF=-2; means in MJD
start_time = datenum(2017,2,1); % start time of run
mjd_init = datenum(1968,5,23);
tval_for_infile = start_time - mjd_init; % enter this value in the .in file
disp(['enter ' num2str(tval_for_infile) ' for TIDE_START value'])
days_to_simulate=365;  % length of model run in days

%% file directories
grid_dir = 'C:\Users\tc196\Dropbox\research\hawkes_bay\model\coawst\input_files\grid\';
grid_file_name = 'roms_grid_002_500m.nc';
grid_mat_file = 'C:\Users\tc196\Dropbox\research\hawkes_bay\model\coawst\input_files\grid\roms_grid_002_500m.mat'; % output from Gridbuilder
grid_file  = [grid_dir '\' grid_file_name]; % full path
data_file ='http://thredds.moanaproject.org:8080/thredds/dodsC/moana/ocean/NZB/v1.9/raw_3D/nz5km_his_201702.nc'; % file to initialize from
bc_file_path ='http://thredds.moanaproject.org:8080/thredds/dodsC/moana/ocean/NZB/v1.9/raw_3D/'; % path to directory of boundary forcing files

%% output file names
tide_ncfile=['C:\Users\tc196\Dropbox\research\hawkes_bay\model\coawst\input_files\tide\' datestr(start_time,'ddmmmyyyy') '.nc'];
init_file_name = 'initial_file.nc';

%% add (recursive) paths to m-files 
addpath(genpath('C:\Users\tc196\Dropbox\research\hawkes_bay\model\coawst\'));
addpath(genpath('/Users/ted/Dropbox/research/hawkes_bay/model/coawst/'))
addpath(genpath('C:\Users\tc196\Dropbox\research\hawkes_bay\model\coawst\input_files\tide'))

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
% In cppdefs.h you should have 
% #undef ana_initial
% interpolates nans

if make_init
make_init_file(grid_file,data_file,init_file_name,start_time,grid_mat_file)
disp(['created initial conditions file:' init_file_name ])
end

%% boundary forcing

% if make_bc_files
% interp_boundary extract_bry, 
% 
% 
% B = obc_roms2roms(ncfile,D,R,VarList,Tindex,boundary, ...
%                     Hmethod,offset,RemoveNaN);
% 
% B = obc_roms2roms(ncfile,D,R,VarList,Tindex,boundary,varargin)
% 
% make_bry_files

%% rivers
%create_roms_river

%% nesting files

%% atmospheric forcing
% Hau-moana will be done ~september
% constant forcing for now, or reanalysis 

%% waves

%% sediment
%% other files
%add_mask add_sponge add_drag
