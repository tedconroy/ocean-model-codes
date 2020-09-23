function make_bry_files(BRY_base,grid_file,input_files,ref_time)
% make_bry_files ( BRYname , grid_file , data_file , input_files )
% make boundary forcing files (one way nesting) from other roms domain
% will make length(input_files) number of bry files
%
% input:
% BRYname = name of boundary condition file to be created
% grid file netcdf
% input files: files to create forcing from
% roms reference time to use
%
% example:
% make_bry_files(bc_file_name,grid_file,input_files,ref_time)
%
% todo: 
% automate which boundary segments to do.
% remove any hard coded stuff (NT,..)
% 
%
% tc2020



%% loop through monthly files, make a bry file for each one
for f=1:length(input_files)
    
    % file name to create
    BRYname = [BRY_base num2str(f) '.nc'];
    
    % grid info
    [Lr,Mr] = size(ncread(grid_file,'h'));
    Lu = Lr-1;
    Lv = Lr;
    Mu = Mr;
    Mv = Mr-1;
    % create structure S with grid info
    S.ncname      = BRYname;    % output NetCDF file
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
    S.hc          = S.Tcline;   % S-coordinate stretching width, ok for vtransform = 2
    
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
    f1 = input_files(1).name; % grab the first file to use for grid info etc
    
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
    
    %  ROMS state variables to read from input files
    VarBry  = {'zeta', 'u', 'v', 'temp', 'salt', 'ubar', 'vbar'};
    VarList = VarBry;
    
    %%  Create boundary condition Netcdf file.
    
    [status]=c_boundary(S); % create file
    
    %  Set attributes for "bry_time".
    avalue='seconds since 2010-01-01 00:00:00';
    [status]=nc_attadd(BRYname,'units',avalue,'bry_time');
    avalue='time units';
    [status]=nc_attadd(BRYname,'calendar',avalue,'bry_time');
    
    %  Set global attribute.
    avalue='for hawke bay model';
    [status]=nc_attadd(BRYname,'title',avalue);
    avalue='sourced from Moana hindcast mdoel';
    [status]=nc_attadd(BRYname,'source',avalue);
    [status]=nc_attadd(BRYname,'grd_file',grid_file);
    
    %  Write out grid data.
    for var = VarGrd
        field = char(var);
        ncwrite(BRYname, field, T.(field));
    end
    
    %%  Interpolate boundary conditions to regional grid into multiple bry files
    
    curfile = input_files(f).name;
    cur_time_vec = ncread(curfile,'ocean_time');
    
    for Rec = 1:length(cur_time_vec)
        
        timecur = cur_time_vec(Rec);
        
        %  Interpolate boundary conditions
        method = 'linear';
        offset = 5;
        RemoveNaN = true;
        % added a 'rotated vectors' print statement to obc_roms2roms to
        % make sure it is rotating to receiver grid
        B = obc_roms2roms(curfile, d, T, VarBry, Rec, OBC, method, offset, RemoveNaN);
                
        %  Set boundary conditions time
        B.bry_time = timecur; % existing time format is correct
        
        %  Write out boundary conditions.
        for var = fieldnames(B)'
            field = char(var);
            [err.(field)] = nc_write(BRYname, field, B.(field), Rec);
        end
        
        % if processing OpenDAP files, force Java garbage collection.
        [~,url,~] = nc_interface(curfile);
        if (url)
            java.lang.System.gc
        end
        
    end
    
    % clear structures, keep variables needed in loop for next file
    disp(['CREATED: ' BRYname])
    clearvars -except BRY_base grid_file input_files grid_mat_file f
    
end
end

