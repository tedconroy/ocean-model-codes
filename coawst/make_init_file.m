function make_init_file(grid_file,data_file,init_file_name,start_time,grid_mat_file)
% 
%
% interpolates existing_output_file onto grid_file
%
% todo: 
%  - add sediment (e.g. from create_roms_init), use this to create file create_roms_netcdf_init_mw
%  - add option to use same grid (no interpolation needed)
%
% tc2020

%% Grid info of child
% a lot of this might not be necessary...

% structure with grid info % edit these to not manually enter 
S.ncname = init_file_name; % name of file
load(grid_mat_file)
if SG.grid.coord == 'spherical'; S.spherical = 1; else; S.spherical = 0; end
S.Vtransform=SG.Z.ROMS.Vtransform; %vertical transformation equation
S.N = SG.Z.ROMS.N; %number of vertical levels
S.NT = 2; % number tracers
S.Lm =258; S.Mm=518; 
S.Vstretching =4;  %vertical stretching function
S.theta_s     =8;      %surface control parameter
S.theta_b     =4;      %bottom  control parameter
S.Tcline      =20;       %surface/bottom stretching width

h=ncread(grid_file,'h');
hmin=min(h(:));
if (S.Vtransform==1)
    hmin=min(h(:));
    hc=min(max(hmin,0),S.Tcline);
elseif (S.Vtransform==2)
    hc=S.Tcline;
end
S.hc          =hc;           %stretching width used in ROMS
[LP,MP]=size(h);
L  = LP-1;
M  = MP-1;
xi_psi  = L;
xi_rho  = LP;
xi_u    = L;
xi_v    = LP;
eta_psi = M;
eta_rho = MP;
eta_u   = MP;
eta_v   = M;

Gout=get_roms_grid(grid_file,S); % grid structure

%% make initial (template) file 
[error] = c_initial(S) 

%% more grid stuff

%  Set attributes for "ocean_time".
avalue='seconds since 0001-01-01 00:00:00';
[~]=nc_attadd(init_file_name,'units',avalue,'ocean_time');
  
avalue='360.0 days in every year';
[~]=nc_attadd(init_file_name,'calendar',avalue,'ocean_time');

%  Set grid variables.
V=nc_vnames(grid_file);
nvars=length(V.Variables);

%  Horizontal grid variables. Read in for input GRID NetCDF file.

if (S.spherical),
  S.lon_rho = nc_read(grid_file, 'lon_rho');
  S.lat_rho = nc_read(grid_file, 'lat_rho');
  S.lon_u   = nc_read(grid_file, 'lon_u');
  S.lat_u   = nc_read(grid_file, 'lat_u');  
  S.lon_v   = nc_read(grid_file, 'lon_v');
  S.lat_v   = nc_read(grid_file, 'lat_v');
else  
  S.x_rho   = nc_read(grid_file, 'x_rho');
  S.y_rho   = nc_read(grid_file, 'y_rho');  
  S.x_u     = nc_read(grid_file, 'x_u');
  S.y_u     = nc_read(grid_file, 'y_u');  
  S.x_v     = nc_read(grid_file, 'x_v');
  S.y_v     = nc_read(grid_file, 'y_v');  
end  

%  Read in Land/Sea mask, if appropriate.
for n=1:nvars,
  name=char(V.Variables(n).Name);
  switch (name),
    case 'mask_rho'
      S.mask_rho = nc_read(grid_file, 'mask_rho');
    case 'mask_u'
      S.mask_u   = nc_read(grid_file, 'mask_u');
    case 'mask_v'
      S.mask_v   = nc_read(grid_file, 'mask_v');
  end,
end,

%  Bathymetry.
S.h = nc_read(grid_file, 'h');
%  Set vertical grid variables.
[S.s_rho, S.Cs_r]=stretching(S.Vstretching, ...
                             S.theta_s, S.theta_b, S.hc, S.N,           ...
			     0, 1);
[S.s_w,   S.Cs_w]=stretching(S.Vstretching, ...
                             S.theta_s, S.theta_b, S.hc, S.N,           ...
			     1, 1);
%  Set zero initial conditions.
Lr = S.Lm+2;   Lu = Lr-1;   Lv = Lr;
Mr = S.Mm+2;   Mu = Mr;     Mv = Mr-1;
S.zeta = zeros([Lr Mr]);
S.ubar = zeros([Lu Mu]);
S.vbar = zeros([Lv Mv]);
S.u    = zeros([Lu Mu S.N]);
S.v    = zeros([Lv Mv S.N]);
S.temp = zeros([Lr Mr S.N]);
S.salt = zeros([Lr Mr S.N]);

%  If Land/Sea masking arrays are not found, initialize them to unity.
if (~isfield(S, 'mask_rho')),  S.mask_rho = ones([Lr Mr]);  end,
if (~isfield(S, 'mask_u'  )),  S.mask_u   = ones([Lu Mu]);  end,
if (~isfield(S, 'mask_v'  )),  S.mask_v   = ones([Lv Mv]);  end,

%  Write out grid variables.
			 
[~]=nc_write(init_file_name,   'spherical',   S.spherical);
[~]=nc_write(init_file_name,   'Vtransform',  S.Vtransform);
[~]=nc_write(init_file_name,   'Vstretching', S.Vstretching);
[~]=nc_write(init_file_name,   'theta_s',     S.theta_s);
[~]=nc_write(init_file_name,   'theta_b',     S.theta_b);
[~]=nc_write(init_file_name,   'Tcline',      S.Tcline);
[~]=nc_write(init_file_name,   'hc',          S.hc);
[~]=nc_write(init_file_name,   's_rho',       S.s_rho);
[~]=nc_write(init_file_name,   's_w',         S.s_w);
[~]=nc_write(init_file_name,   'Cs_r',        S.Cs_r);
[~]=nc_write(init_file_name,   'Cs_w',        S.Cs_w);
[~]=nc_write(init_file_name,   'h',           S.h);
if (S.spherical),
  [~]=nc_write(init_file_name, 'lon_rho',     S.lon_rho);
  [~]=nc_write(init_file_name, 'lat_rho',     S.lat_rho);
  [~]=nc_write(init_file_name, 'lon_u',       S.lon_u);
  [~]=nc_write(init_file_name, 'lat_u',       S.lat_u);
  [~]=nc_write(init_file_name, 'lon_v',       S.lon_v);
  [~]=nc_write(init_file_name, 'lat_v',       S.lat_v);
else
  [~]=nc_write(init_file_name, 'x_rho',       S.x_rho);
  [~]=nc_write(init_file_name, 'y_rho',       S.y_rho);
  [~]=nc_write(init_file_name, 'x_u',         S.x_u);
  [~]=nc_write(init_file_name, 'y_u',         S.y_u);
  [~]=nc_write(init_file_name, 'x_v',         S.x_v);
  [~]=nc_write(init_file_name, 'y_v',         S.y_v);
end

%  Compute depths at horizontal and vertical RHO-points.
igrid = 1;
[z_r] = set_depth(S.Vtransform, S.Vstretching,                          ...
                  S.theta_s, S.theta_b, S.hc, S.N,                      ...
                  igrid, S.h, S.zeta);

%% interpolate onto initial file
% variables needed: zeta, salt, temp, u, v, ubar, vbar
% make sure that using stretching.m and set_depth.m from rutgers/utility
% folder and NOT mtools folder, which doesn't support vstretching=5

% find nearest time index
t = datenum(0,0,0,0,0,ncread(data_file,'ocean_time'))+datenum(2010,1,1); % moana project date format
ind  = dsearchn(t,start_time); % find nearest time match
disp(['selecting nearest time to start_time (' datestr(start_time) ') = ' datestr(t(ind))])

% convert time to moana project time format
init_time  = etime(datevec(start_time),datevec(datenum(2010,1,1)));

d = get_roms_grid(data_file); % donor grid info

% list of variable names to interpolate from data_file (some are Moana Project specific)
field_names = {'zeta' 'ubar_eastward' 'vbar_northward' 'u_eastward' 'v_northward' ...
    'salt' 'temp'};
% change variable names in init_file back to ROMS standard (same order)
field_names_update = {'zeta' 'ubar' 'vbar' 'u' 'v' 'salt' 'temp'};

% loop through variables, interpolate, write out
for n=1:length(field_names)
v = roms2roms (data_file,d,Gout,field_names{n},ind,1,'natural',5,'true');
ncwrite(init_file_name,field_names_update{n},v); % write out var
disp(['writing var to init file: ' field_names_update{n}])
end

% write out static variables
ncwrite(init_file_name,'theta_s',Gout.theta_s);
ncwrite(init_file_name,'theta_b',Gout.theta_b);
ncwrite(init_file_name,'Tcline',Gout.Tcline);
ncwrite(init_file_name,'Cs_r',Gout.Cs_r);
ncwrite(init_file_name,'Cs_w',Gout.Cs_w);
ncwrite(init_file_name,'sc_w',Gout.s_w);
ncwrite(init_file_name,'sc_r',Gout.s_rho);
ncwrite(init_file_name,'hc',Gout.hc);
ncwrite(init_file_name,'Vtransform',Gout.Vtransform);
ncwrite(init_file_name,'Vstretching',Gout.Vstretching);
ncwrite(init_file_name,'spherical',Gout.spherical);
ncwrite(init_file_name,'ocean_time',init_time);

disp(['created ', init_file_name])



