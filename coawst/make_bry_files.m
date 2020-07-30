function make_bry_files(BRYname,grid_file,input_files,grid_mat_file)
% make_bry_files ( BRYname , grid_file , data_file , input_files )
% make boundary forcing files (one way nesting) from other roms domain
%
% input:
% BRYname = name of boundary condition file to be created
%
% input_files = output of dir command
%
% todo: 
% remove hard coded stuff
% automate which boundary segments to do.
% consistent time in outer file
% 
% notes:
% - could also use interp_boundary.m?
% - obc_roms2roms interpolates u,v to rho, had to average back to u,v
% points
% 
%
% tc2020


%% grid stuff
% hard coded

[Lr,Mr] = size(nc_read(grid_file,'h'));
Lu = Lr-1;
Lv = Lr;
Mu = Mr;
Mv = Mr-1;
S.ncname      = BRYname;    % output NetCDF file
load(grid_mat_file)
if SG.grid.coord == 'spherical'; S.spherical = 1; else; S.spherical = 0; end
S.Vtransform=SG.Z.ROMS.Vtransform; %vertical transformation equation
S.Lm          = Lr-2;       % number of interior RHO-points, X-direction
S.Mm          = Mr-2;       % number of interior RHO-points, Y-direction
S.N           = 15;         % number of vertical levels at RHO-points
S.NT          = 2;          % total number of tracers
S.Vstretching =4;  %vertical stretching function
S.theta_s     =8;      %surface control parameter
S.theta_b     =4;      %bottom  control parameter
S.Tcline      =20;       %surface/bottom stretching width
S.hc          = S.Tcline;   % S-coordinate stretching width

%  Set switches for boundary segments to process.
% *automate this based on mask
OBC.west  = false;          % process western  boundary segment
OBC.east  = true;           % process eastern  boundary segment
OBC.south = true;           % process southern boundary segment
OBC.north = true;           % process northern boundary segment
S.boundary(1) = OBC.west;
S.boundary(2) = OBC.east;
S.boundary(3) = OBC.south;
S.boundary(4) = OBC.north;

%  Get parent and target grids structures
%%f1 = ([input_files(1).folder '/' input_files(1).name]); % get first file for grid info
f1 = input_files; % grab the first file to use for grid info etc

d = get_roms_grid(f1); % donor/parent grid info
T = get_roms_grid(grid_file,S); % child/receiver/target grid structure 

%  interpolate rotation angle (parent to target) and add it to target grid
T.parent_angle = roms2roms(f1, d, T, 'angle', [], true, 'linear', 5, true); % linear best option?

%%  Set variables to process.
VarGrd = {'spherical',                                                ...
          'Vtransform', 'Vstretching',                                ...
          'theta_s', 'theta_b', 'Tcline', 'hc',                       ...
          's_rho', 'Cs_r', 's_w', 'Cs_w'};
if (S.spherical)
  if (OBC.west)
    VarGrd = [VarGrd, 'lon_rho_west',  'lat_rho_west',                ...
                      'lon_u_west',    'lat_u_west',                  ...
                      'lon_v_west',    'lat_v_west'];
  end
  if (OBC.east)
    VarGrd = [VarGrd, 'lon_rho_east',  'lat_rho_east',                ...
                      'lon_u_east',    'lat_u_east',                  ...
                      'lon_v_east',    'lat_v_east'];
  end
  if (OBC.south)
    VarGrd = [VarGrd, 'lon_rho_south', 'lat_rho_south',               ...
                      'lon_u_south',   'lat_u_south',                 ...
                      'lon_v_south',   'lat_v_south'];
  end
  if (OBC.north)
    VarGrd = [VarGrd, 'lon_rho_north', 'lat_rho_north',               ...
                      'lon_u_north',   'lat_u_north',                 ...
                      'lon_v_north',   'lat_v_north'];
  end
else
  if (OBC.west)
    VarGrd = [VarGrd, 'x_rho_west',  'y_rho_west',                    ...
                      'x_u_west',    'y_u_west',                      ...
                      'x_v_west',    'y_v_west'];
  end
  if (OBC.east)
    VarGrd = [VarGrd, 'x_rho_east',  'y_rho_east',                    ...
                      'x_u_east',    'y_u_east',                      ...
                      'x_v_east',    'y_v_east'];
  end
  if (OBC.south)
    VarGrd = [VarGrd, 'x_rho_south', 'y_rho_south',                   ...
                      'x_u_south',   'y_u_south',                     ...
                      'x_v_south',   'y_v_south'];
  end
  if (OBC.north)
    VarGrd = [VarGrd, 'x_rho_north', 'y_rho_north',                   ...
                      'x_u_north',   'y_u_north',                     ...
                      'x_v_north',   'y_v_north'];
  end
end

%  ROMS state variables to process. (moana project specific)
VarBry  = {'zeta', 'u', 'v', 'temp', 'salt', 'ubar', 'vbar'};
VarList = [VarBry];

% var names to write (roms general)
VarBry_write  = {'zeta', 'u', 'v', 'temp', 'salt', 'ubar', 'vbar'};
VarList_write = [VarBry];
                      
%%  Create boundary condition Netcdf file.

  [status]=c_boundary(S); % create template file

%  Set attributes for "bry_time".
  avalue='seconds since 2010-01-01 00:00:00';
  [status]=nc_attadd(BRYname,'units',avalue,'bry_time');
  avalue='gregorian';
  [status]=nc_attadd(BRYname,'calendar',avalue,'bry_time');

%  Set global attribute.
  avalue='for hawke bay model, grid L2';
  [status]=nc_attadd(BRYname,'title',avalue);
  avalue='MOANA HINDCAST MODEL, grid L1';
  [status]=nc_attadd(BRYname,'source',avalue);
  [status]=nc_attadd(BRYname,'grd_file',grid_file);
  
  % reference time
  ref_time = (datenum('01-Jan-2010')-datenum('01-Jan-1900'))*86400;

%  Write out grid data.
  for var = VarGrd
    field = char(var);
    [err.(field)] = nc_write(BRYname, field, T.(field));
  end                     
 
%%  Interpolate boundary conditions to regional grid
% use input_files (change later)

% loop through monthly files
for f=1% :length(input_files)
    curfile = input_files; % (f)
    for Rec = 1:length(ncread(input_files,'ocean_time'))
        timecur = ncread(curfile,'ocean_time',[Rec],[1]);
        
        %  Interpolate boundary conditions
        method = 'linear';
        offset = 5;
        RemoveNaN = true;
        Rvector = true;         % interpolate vectors to RHO-points
        B = obc_roms2roms(curfile, d, T, VarBry, Rec, OBC, method, offset, RemoveNaN);
        
        %  Set Target Grid vertical level thicknesses, Hz.
        % use time dependepent values
        
        zeta = roms2roms(curfile, d, T, 'zeta', Rec, Rvector,             ...
            method, offset, RemoveNaN);
        
        N = S.N;
        igrid = 5;
        [z_w] = set_depth(T.Vtransform, T.Vstretching,                    ...
            T.theta_s, T.theta_b, T.hc, N,                  ...
            igrid, T.h, zeta, 0);
        Hz = z_w(:,:,2:N+1) -z_w(:,:,1:N);
        
        % move velocity from rho points to u,v grid. They were interpolated for
        % rotation.
        %  Vector components that require rotation were interpolated at
        %  RHO-points.  They are interpolated over two points adjacent to
        %  the boundary edge to allow averaging to the appropriate C-grid
        %  location.
        % renaming same variable, bad practice
        %
        % B.u_north = (B.u_north(1:end-1,:)+B.u_north(2:end,:))./2;
        % B.u_south = (B.u_south(1:end-1,:)+B.u_south(2:end,:))./2;
        % B.u_east = B.u_east;
        % B.v_north = B.v_north;
        % B.v_south = B.v_south;
        % B.v_east = (B.v_east(1:end-1,:)+B.v_east(2:end,:))./2;
        % B.ubar_north = (B.ubar_north(1:end-1,:)+B.ubar_north(2:end,:))./2;
        % B.ubar_south = (B.ubar_south(1:end-1,:)+B.ubar_south(2:end,:))./2;
        % B.ubar_east = B.ubar_east;
        % B.vbar_north = B.vbar_north;
        % B.vbar_south = B.vbar_south;
        % B.vbar_northward_east = B.vbar_east';
        % B.vbar_east = (B.vbar_east(1:end-1,:)+B.vbar_east(2:end,:))./2;
        
        %  Set boundary conditions time (seconds since 2000-01-01 00:00:00).
        B.bry_time = timecur ; %- ref_time; % using same time format
        
        %  Write out boundary conditions.
        %BryRec = BryRec+1;
        varlist = fieldnames(B)';
        varlist(contains(varlist,'ward'))=[]; % get rid of eastward, northward varibales..
        for var = varlist
            field = char(var);
            [err.(field)] = nc_write(BRYname, field, B.(field), Rec);
        end
        
        %  Process next boundary record. If processing OpenDAP files, force Java
        %  garbage collection.
        [~,url,~] = nc_interface(curfile);
        if (url),
            java.lang.System.gc
        end
    end
end
end

