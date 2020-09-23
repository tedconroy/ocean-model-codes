% make refinement grids for coawst roms model
% grids will be referenced based on levels of nesting
%
%
% the points for the refinement grid need to be exact rho points from the
% coarse grid
%
% might also look into making or obtaining a very high resolution coastline
% data
%
%
% make grids then interp bathy and mask? landsea.m, editmask.m

addpath(genpath( 'C:\Users\tc196\Dropbox\research\hawkes_bay\model\coawst'))

%% specify information about grids

%  coarsest grid file
c_g = 'sept_2020_v2.nc';

lon_rho = ncread(c_g,'lon_rho');
lat_rho = ncread(c_g,'lat_rho');
dep = ncread(c_g,'h');
mask_rho = logical(ncread(c_g,'mask_rho'));
dep(~mask_rho) = NaN;

%% make structure (ng) for new grid(s) to be created below
% going with other users convention, grids go from L0 --> LN refinement
% levels

%% L1: 

% Napier Area
ng(1).name = 'l1_south.nc'; % name of file to be created
% specify coordinate bounds
ng(1).Istr = 12; % lower left i coord
ng(1).Iend = 63; % upper right i coord
ng(1).Jstr = 51; % lower left j coord
ng(1).Jend = 105; % upper right j coord
ng(1).ref = 3; % resolution refinement factor

% Northern bay
%ng(1).name = 'l1_north.nc'; % name of file to be created
% specify coordinate bounds
%ng(1).Istr = 73; % lower left i coord
%ng(1).Iend = 158; % upper right i coord
%ng(1).Jstr = 183; % lower left j coord
%ng(1).Jend = 189; % upper right j coord
%ng(1).ref = 5; % resolution refinement factor

%% L2:
%% L3: 

% initiate name list for contact file
Gnames{1} = c_g; % add coarse grid file

% loop to create new grids
for n = 1:length(ng)
    F = coarse2fine(c_g, ng(n).name, ng(n).ref, ng(n).Istr, ng(n).Iend, ng(n).Jstr, ng(n).Jend);
    Gnames{n+1} = ng(n).name; % store file name for contact file generation
end


%% ----------
%% bathy, mask

% load compiled bathy data
data = load('C:\Users\tc196\Dropbox\research\hawkes_bay\data\bathymetry\merged_xyz_4roms.mat');
i = isinf(data.x); data.x(i) = []; data.y(i)=[]; data.z(i)=[];

% -----------
% coarse grid
% -----------

% already edited bathy in gridbuilder

% griddata v4 doesnt work - convert to nztm? doesnt work. tries to make a huge
% matrix
%[data.xx,data.yy]=pyproj_transform(data.x,data.y,4326,2193);
%[x,y] = pyproj_transform(lon_rho,lat_rho,4326,2193);
%nh = griddata(data.x,data.y,data.z,lon_rho,lat_rho,'natural'); 
% set min depth
MIN_dep = 4;
%nh(find(nh < MIN_dep)) = MIN_dep;
% any smoothing?
% write out new depth
%ncwrite(c_g,'h',nh)
%disp(['wrote new depth to ' c_g])

% -----------
% nested grid
% -----------

L1 = ng(1).name;
lon = ncread(L1,'lon_rho'); lat = ncread(L1,'lat_rho');
h = ncread(L1,'h');
nh = griddata(data.x,data.y,data.z,lon,lat,'natural'); 

% set min depth
nh(find(nh < MIN_dep)) = MIN_dep;

% write out new depth
ncwrite(L1,'h',nh)
disp(['wrote new depth to ' L1])


% change mask
f = landsea(L1,'f'); % looks ok, not great
editmask(L1)

% need to update vertical grid information
list = {'Tcline' 'theta_b' 'theta_s' 'Vstretching' 'Vtransform' 's_rho' 's_w' 'Cs_w' 'Cs_r' 'hc'}
for var = list
    v = ncread(c_g,char(var));
    ncwrite(L1,char(var),v)
end

%%  make contact files
% there are 2*(ngrids-1) contact regions
[S,G] = contact(Gnames, 'contact_file.nc',0,0,1);
disp(['NCONTACT = ' num2str(2*(length(Gnames)-1))])


%% plot
close all
pcolor(lon_rho,lat_rho,dep);shading flat;
caxis([0 100]); colorbar; colormap(cmocean('deep'))
latratio(-40)
hold on
nmask = logical( ncread(L1,'mask_rho'));
nh(~nmask)=NaN;
pcolor(lon,lat,nh);shading flat

