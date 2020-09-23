function make_init_file(grid_file,data_file,init_file_name,start_time,grid_mat_file)
% interpolates existing_output_file onto grid_file using roms tools
%
% input
% grid_file: full path to grid nc file
% data_file: full path to file to initialize from
% init_file_name: initial file name to be created
% grid_mat_file: sourced from gridbuilder
%
%
%
% todo: 
%  - remove hard coded stuff
%  - add sediment (e.g. from create_roms_init), use this to create file create_roms_netcdf_init_mw
%  - add option to use same grid (no interpolation needed)
%
% tc2020


%% Grid info of child
% structure with grid info % edit these to not manually enter
S.ncname = init_file_name; % name of file

% grid info
[Lr,Mr] = size(ncread(grid_file,'h'));
Lu = Lr-1;
Lv = Lr;
Mu = Mr;
Mv = Mr-1;
% create structure S with grid info
S.spherical = ncread(grid_file,'spherical');
S.Vtransform = ncread(grid_file,'Vtransform'); %vertical transformation equation
S.Lm          = Lr-2;       % number of interior RHO-points, X-direction
S.Mm          = Mr-2;       % number of interior RHO-points, Y-direction
S.N = length(ncread(grid_file,'s_rho')); % number of vertical levels at rho
S.NT          = 2;          % total number of tracers
S.Vstretching = ncread(grid_file,'Vstretching');  %vertical stretching function
S.theta_s     = ncread(grid_file,'theta_s');      %surface control parameter
S.theta_b     = ncread(grid_file,'theta_b');      %bottom  control parameter
S.Tcline      = ncread(grid_file,'Tcline');       %surface/bottom stretching width
S.hc          = S.Tcline;   % S-coordinate stretching width

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

%% make initial file
[error] = c_initial(S)

%  Set attributes for "ocean_time".
avalue='seconds since 2010-01-01 00:00:00';
[~]=nc_attadd(init_file_name,'units',avalue,'ocean_time');

avalue='time units';
[~]=nc_attadd(init_file_name,'calendar',avalue,'ocean_time');

%  Set grid variables.
V=nc_vnames(grid_file);
nvars=length(V.Variables);

%  Horizontal grid variables. Read in for input GRID NetCDF file.

if (S.spherical),
    S.lon_rho = ncread(grid_file, 'lon_rho');
    S.lat_rho = ncread(grid_file, 'lat_rho');
    S.lon_u   = ncread(grid_file, 'lon_u');
    S.lat_u   = ncread(grid_file, 'lat_u');
    S.lon_v   = ncread(grid_file, 'lon_v');
    S.lat_v   = ncread(grid_file, 'lat_v');
else
    S.x_rho   = ncread(grid_file, 'x_rho');
    S.y_rho   = ncread(grid_file, 'y_rho');
    S.x_u     = ncread(grid_file, 'x_u');
    S.y_u     = ncread(grid_file, 'y_u');
    S.x_v     = ncread(grid_file, 'x_v');
    S.y_v     = ncread(grid_file, 'y_v');
end

%  Read in Land/Sea mask, if appropriate.
for n=1:nvars,
    name=char(V.Variables(n).Name);
    switch (name),
        case 'mask_rho'
            S.mask_rho = ncread(grid_file, 'mask_rho');
        case 'mask_u'
            S.mask_u   = ncread(grid_file, 'mask_u');
        case 'mask_v'
            S.mask_v   = ncread(grid_file, 'mask_v');
    end,
end,

%  Bathymetry.
S.h = ncread(grid_file, 'h');
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
ncwrite(init_file_name,   'spherical',   S.spherical);
ncwrite(init_file_name,   'Vtransform',  S.Vtransform);
ncwrite(init_file_name,   'Vstretching', S.Vstretching);
ncwrite(init_file_name,   'theta_s',     S.theta_s);
ncwrite(init_file_name,   'theta_b',     S.theta_b);
ncwrite(init_file_name,   'Tcline',      S.Tcline);
ncwrite(init_file_name,   'hc',          S.hc);
ncwrite(init_file_name,   's_rho',       S.s_rho);
ncwrite(init_file_name,   's_w',         S.s_w);
ncwrite(init_file_name,   'Cs_r',        S.Cs_r);
ncwrite(init_file_name,   'Cs_w',        S.Cs_w);
ncwrite(init_file_name,   'h',           S.h);
if (S.spherical),
    ncwrite(init_file_name, 'lon_rho',     S.lon_rho);
    ncwrite(init_file_name, 'lat_rho',     S.lat_rho);
    ncwrite(init_file_name, 'lon_u',       S.lon_u);
    ncwrite(init_file_name, 'lat_u',       S.lat_u);
    ncwrite(init_file_name, 'lon_v',       S.lon_v);
    ncwrite(init_file_name, 'lat_v',       S.lat_v);
else
    ncwrite(init_file_name, 'x_rho',       S.x_rho);
    ncwrite(init_file_name, 'y_rho',       S.y_rho);
    ncwrite(init_file_name, 'x_u',         S.x_u);
    ncwrite(init_file_name, 'y_u',         S.y_u);
    ncwrite(init_file_name, 'x_v',         S.x_v);
    ncwrite(init_file_name, 'y_v',         S.y_v);
end

%  Compute depths at horizontal and vertical RHO-points.
igrid = 1;
[z_r] = set_depth(S.Vtransform, S.Vstretching,                          ...
    S.theta_s, S.theta_b, S.hc, S.N,                      ...
    igrid, S.h, S.zeta);

%% interpolate onto initial file
% make sure that using stretching.m and set_depth.m from rutgers/utility
% folder and NOT mtools folder, which doesn't support vstretching=5

% find nearest time index
t = datenum(0,0,0,0,0,ncread(data_file,'ocean_time'))+datenum(2010,1,1); % moana project date format
ind  = dsearchn(t,start_time) ; % find nearest time match
disp(['selecting nearest time to start_time (' datestr(start_time) ') as ' datestr(t(ind))])

% convert time to moana project time format
init_time  = etime(datevec(start_time),datevec(datenum(2010,1,1)));

d = get_roms_grid(data_file); % donor grid info

% list of variable names to interpolate from data_file
field_names = {'zeta' 'ubar' 'vbar' 'u' 'v' ...
    'salt' 'temp'};

% loop through variables, interpolate
for var=field_names
    field = char(var);
    % call roms2roms interpolation
    I.(field) = roms2roms (data_file,d,Gout,field,ind,1,'natural',5,'true');
end


%  Rotate interpolated 3D velocity at RHO-points to TRUE North and East.
%  Need to interpolate Parent grid rotation angle to Target grid.
irotate = 1;               % rotate for (XI,ETA) to (lon,lat)
method = 'linear';             % linear interpolation
offset = 10;                   % number of extra points for sampling
RemoveNaN = true;              % remove NaN with nearest-neighbor
Rvector = true;                % interpolate vectors to RHO-points

Gout = get_roms_grid(grid_file, S);

%  If vector rotation is required in the parent grid, interpolate
%  rotation angle (parent to target) and add it to target grid
%  structure.
Gout.parent_angle = roms2roms(data_file, d, Gout, 'angle', [], Rvector, method, offset, RemoveNaN);
[Urho,Vrho] = rotate_vec(I.u, I.v, Gout.parent_angle, irotate);
%  Rotate resulting 3D velocity (RHO-points) to target grid angle and
%  average to staggered C-grid locations.
[I.u,I.v] = roms_vectors(Urho, Vrho, Gout.angle, Gout.mask_u, Gout.mask_v);
%  Compute barotropic velocities by vertically integrating (u,v).
[I.ubar,I.vbar] = uv_barotropic(I.u, I.v, Gout.Hz);

% write variables
for var = field_names
    field = char(var);
    [err.(field)] = nc_write(init_file_name, field, I.(field), 1);
    disp(['writing var to init file: ' var])
end

% write out static variables
ncwrite(init_file_name,'theta_s',Gout.theta_s);
ncwrite(init_file_name,'theta_b',Gout.theta_b);
ncwrite(init_file_name,'Tcline',Gout.Tcline);
ncwrite(init_file_name,'Cs_r',Gout.Cs_r);
ncwrite(init_file_name,'Cs_w',Gout.Cs_w);
ncwrite(init_file_name,'s_w',Gout.s_w);
ncwrite(init_file_name,'s_rho',Gout.s_rho);
ncwrite(init_file_name,'hc',Gout.hc);
ncwrite(init_file_name,'Vtransform',Gout.Vtransform);
ncwrite(init_file_name,'Vstretching',Gout.Vstretching);
ncwrite(init_file_name,'spherical',Gout.spherical);
ncwrite(init_file_name,'ocean_time',init_time);

disp(['created ', init_file_name])
end


