% generate a COOS Bay FVCOM setup using Geoff's codes, SMS, and the FVCOM toolboxes
% ted conroy, 2018 (U. of Oregon)

%***************************************************************************
clear all; close all;
addpath ~/Dropbox/CoosBay/model/fvcom_prepro
addpath ~/Dropbox/CoosBay/model
%--------------------------------------------------------------------------
% output file names 
%--------------------------------------------------------------------------
modelid = 'coos_f18_tc';
modeldir = 'Users/ted/Dropbox/CoosBay/model/high/runfiles/run_fall2018_tc/'; % ~ doesnt work for creating ncfiles
if ~exist(modeldir); mkdir(modeldir); end

do_elj = 1; % do elevation forcing-- or just copy from old run
do_rivs = 1;
do_IC = 0;  % initialize t/s field
do_wind = 1;  % do wind forcing-- or just copy from old run

probefile = [modeldir, modelid '_probes.nml'];
obcfile = [modeldir, modelid '_obc.dat'];
corfile = [modeldir, modelid '_cor.dat'];
depfile = [modeldir, modelid '_dep.dat'];
spgfile = [modeldir, modelid '_spg.dat'];
grdfile = [modeldir, modelid '_grd.dat'];
elvfile = [modeldir, modelid '_elj.nc'];
rivfile_base = [modeldir, modelid '_riv.nc'];
wndfile = [modeldir, modelid '_wnd.nc'];
z0bfile = [modeldir, modelid '_z0.nc'];
inifile = [modeldir, modelid '_ini.nc'];

%--------------------------------------------------------------------------
% read in the SMS mesh and add lat/lon, metrics, etc. 
%--------------------------------------------------------------------------
mesh_dir = '/Users/ted/Dropbox/CoosBay/model/high/gridfiles/2018_march/';
mesh_file = [mesh_dir,'CE_GRID_V3.2dm']; 
Mobj = read_sms_mesh('2dm',mesh_file); % mesh resolution, etc. all set in SMS (see SMS how to guide and _xy.dat)

% if x,y switched in SMS, switch back
% x = Mobj.y; y = Mobj.x;
%Mobj.x = x; Mobj.y = y;
% x,y switched in SMS, so also switch tri
%Mobj.tri = Mobj.tri(:,[1 3 2]); % had to do for v.3.1.6

% fix cell geometry near open boundary 
node_ids=[1:52];
[adjx,adjy]=fix_inside_boundary(Mobj.x,Mobj.y,node_ids); 
Mobj.x=adjx; Mobj.y=adjy;
% mesh is in state plane meters, convert to lat/lon using sp_proj (zone 3602)
[Mobj.lon,Mobj.lat] = sp_proj('3602','inverse',Mobj.x,Mobj.y,'m');  
Mobj.have_lonlat = 1;
% calculate the Corolis
Mobj = add_coriolis(Mobj); 
% add metrics
Mobj = setup_metrics(Mobj);

%--------------------------------------------------------------------
% Bathymetry: interpolate onto grid and smooth 
%--------------------------------------------------------------------
% best to make non-overlapping data set then interpolate
interp_bathy = 1; % note: changing wet points to positive values when interpolating
bathymetry_dir1 = '~/Dropbox/CoosBay/model/bathymetry/'
Mobj.have_bath = 1;

if(interp_bathy)
    
% first interpolate noaa_dem as it covers entire area + lower res
load([bathymetry_dir1 'noaa_dem.mat']);
F=scatteredInterpolant(lon(:),lat(:),-z(:)); [Mobj.h]=F(Mobj.lon,Mobj.lat);
clearvars lon lat z F

 % interp port orford dem in its bounds
load([bathymetry_dir1 'portorford_dem.mat']);
xx=[min(lon(:)) max(lon(:)) max(lon(:)) min(lon(:)) min(lon(:))];
yy=[min(lat(:)) min(lat(:)) max(lat(:)) max(lat(:)) min(lat(:))];
in = inpolygon(Mobj.lon,Mobj.lat,xx,yy);
F=scatteredInterpolant(lon(:),lat(:),-z(:));
[Mobj.h(in)]=F(Mobj.lon(in),Mobj.lat(in));
clearvars lon lat z in xx yy F

% interp gridded 2meter resolution coos bay bathy data
load([bathymetry_dir1 'coosbay2m.mat']);

% interp in bounds of data
xx=[min(lon(:)) max(lon(:)) max(lon(:)) min(lon(:)) min(lon(:))];
yy=[min(lat(:)) min(lat(:)) max(lat(:)) max(lat(:)) min(lat(:))];
in = inpolygon(Mobj.lon,Mobj.lat,xx,yy);
F=scatteredInterpolant(lon(:),lat(:),-z(:));
[Mobj.h(in)] = F(Mobj.lon(in), Mobj.lat(in));
clearvars lon lat z in xx yy F         
% convert from MLLW to MSL
Mobj.h=Mobj.h+1.244;     

% edit bathymetry - use other script
change_mobj_bathy; close all;

%--------------------------------------------------------------------
% Open boundaries -- deal with these for OBC and rivers
%--------------------------------------------------------------------
% if grid changes need to change river nodes
open_boundary_bc = 41; 
river_bc = 42;
load('~/Dropbox/model/forcing_files/river/river_names.mat'); riv(14).name='talbot'; riv(15).name='johnb';
Mobj.nRivers=15;
Mobj.obc_nodes = zeros(Mobj.nRivers+1,500); % set up matrix for all different bndy's
Mobj.nmark = zeros(Mobj.nVerts,1);
[e,te,e2t,bnd] = connectivity([Mobj.x,Mobj.y],Mobj.tri);
Mobj.nmark(bnd) = 1;
nn=1:52;  % for this grid, open boundary nodes are index 1-52
Mobj.nObcNodes(1) = length(nn);
Mobj.obc_type(1) = 1;
Mobj.obc_name(1) = {'open boundary node'};
Mobj.obc_nodes(1,1:length(nn)) = nn;
Mobj.nmark(nn) = open_boundary_bc;

inr=[102863 103064 99889 102006 90663 67697 34542 56693 68072 81914 79284 77087 78089 59880 58212]; % river input nodes (in order)

for i=1:Mobj.nRivers
  nn = length(inr(i));
  Mobj.obc_nodes(i+1,1:nn) = inr(i); % first row is OBC, then rivers
  Mobj.riv_nodes(i, 1:nn) = Mobj.obc_nodes(i+1,1:nn);
  Mobj.nRivNodes(i) = nn;
  Mobj.nObs = 1;
  Mobj.obc_type(i+1) = 2; % changed from 1 to 2
  Mobj.obc_name(i+1) = {'river open boundary node'};
  Mobj.nmark(inr(i)) = river_bc;
  Mobj.nObcNodes(i+1) = nn;
end;

%--------------------------------------------------------------------
% Add sponge layer
%--------------------------------------------------------------------
sponge_rad =5000; % sponge layer radius (m)
sponge_coef = 1e-4; %sponge coefficient (inverse timescale in 1/seconds)
Mobj.nSponge = 1;
% dump sponge layer file
ob=1:52; % number of obc nodes that want sponge
Nlist = Mobj.obc_nodes(ob,1:Mobj.nObcNodes(ob));
Mobj = add_sponge_nodes_list(Mobj,Nlist,'SpongeLayer',sponge_rad,sponge_coef);
write_FVCOM_sponge(Mobj,spgfile);

%--------------------------------------------------------------------
% Create bottom roughness file, OR specify constant value in namelist
%--------------------------------------------------------------------
% make spatially constant z0 -- can just specify constant in namelist...
%my_mode = bitor ( nc_clobber_mode, nc_64bit_offset_mode );
%disp(['*** creating ',z0bfile,' ***']);
%nc_create_empty ( z0bfile, my_mode );
%nc_padheader (z0bfile, 10000 );
%nc_attput(z0bfile, nc_global, 'type', 'Example scheme for bottom roughness');
%nc_add_dimension(z0bfile, 'nele', Mobj.nElems);
%varstruct.Name = 'z0b';
%     varstruct.Dimension = {'nele'};
%        long_name = ['bottom roughness'];
%        units = 'number elements';
%        varstruct.Attribute = struct('Name', ...
%            {'long_name','units'},'Value',{long_name,units});
%     nc_addvar(z0bfile, varstruct);
%hc = nodes2elems(Mobj.h,Mobj);  %bathy on the elements
%nc_varput(z0bfile, 'z0b', .001*ones(Mobj.nElems,1)); % here specify z0

%------------ D.Ralston formulation from 2017 paper, based on potential bedformheight
% Create spatially varying bed roughness, values on elements. (need to make sure these are in meters)
% find average depth from 3 nodes to center of cell, make sure to use current h
%H=[Mobj.h(Mobj.tri(:,1)) Mobj.h(Mobj.tri(:,2)) Mobj.h(Mobj.tri(:,3))]; H=mean(H,2);
%Ah=30; % fudge factor
%z0=H/(Ah*6*30); % following Ralston 2017
% set min and max value
%z0=max(z0,0.0005);
%z0=min(z0,0.004);
%write_FVCOM_z0(z0,'/Users/ted/Dropbox/model/high/runfiles/run_2_2018/coos_z0_ah30.nc','coosz0');
%[lonc latc]=sp_proj('3602','inverse',Mobj.xc,Mobj.yc,'m');
%tri = delaunayTriangulation(lonc,latc);
%K = nearestNeighbor(tri,[Mobj.lon(:) Mobj.lat(:)]);
%plot_field(Mobj,z0(K))

%------------ EFE file Nov2017 - based on measured grain sizes
%load z0_fromgrainsize.mat 
% % set min bottom roughness (for values above MSL)
%i=find(z0<0.05); sum(i)
%z0(i)=0.0005;
%write_FVCOM_z0(z0,'z0_coos.nc','z0 2017 obs');
            
%--------------------------------------------------------------------
% dump the Coriolis, mesh, dep,  open boundary node list 
%--------------------------------------------------------------------
write_FVCOM_cor(Mobj,corfile) ;                 
write_FVCOM_grid(Mobj,grdfile);
write_FVCOM_bath(Mobj,depfile);               
write_FVCOM_obc(Mobj,obcfile);

%--------------------------------------------------------------------
% boundary elevation forcing 
%--------------------------------------------------------------------
% use TXPO modeled elevations, adds charleston subtidal water level
if(do_elj)  
%  txpo modeled tide--- has 13 contstituents & is evaluated at each boundary node
load('~/Dropbox/CoosBay/model/forcing_files/txpo_tide_1min.mat')
time=time-datenum(0,0,0,8,0,0); MJD=MJD-datenum(0,0,0,8,0,0);
Mobj.surfaceElevation=[];
%for i=1:52;  Mobj.surfaceElevation(i,:)=hts(i).ts; end  
% number of obc nodes = nmark == 41
nn = length(find(Mobj.nmark==41));
% add subtidal elevation from charleston to txpo
wl=load('~/Dropbox/CoosBay/model/forcing_files/charleston_tide.mat','height');wl=wl.height;
wltime=load('~/Dropbox/CoosBay/model/forcing_files/charleston_tide.mat','time');wltime=wltime.time;
wlfilt=godinfilt(wl(:)); wlfilt_txpo=interp1(wltime,wlfilt,time);
ht35=hts(35).ts;  clear hts
for i=1:length(wlfilt_txpo)
txpo_sub(i)=ht35(i) + (wlfilt_txpo(i)); % using mid-boundary txpo (think it was causing problems using values for each point)
end 
Mobj.surfaceElevation = repmat(txpo_sub',1,nn)'; 
write_FVCOM_elevtide(Mobj, MJD, elvfile, modelid)

disp(['*** creating ',elvfile,' ***']);
else
    disp(['skipped writing ',elvfile,' *** copy from elsewhere!']);
end
%--------------------------------------------------------------------
% river forcing & flux estimates
%--------------------------------------------------------------------
if(do_rivs)
load('~/Dropbox/CoosBay/Obs_Data/Rivers/CWA_excel_discharge_data/cooswa_river_data.mat')
time=river_data(5).time;

% scaling ungauged discharges:
% linear relationship winchester/marlow creek (Marlow is smallest of gauged
% creeks during 2014, and most similar to winchester qr when it was gauged.)
win=0.56.*river_data(3).qr +0.13; % linear fit
win(2923:4018)=river_data(7).qr(1:1096);  % use real winchester data when available 
win_marlow_relarea=river_data(3).qr.*(24.68/15.56); % relative area
win=min(win,win_marlow_relarea); % this more accurately captures summer qr (only changes during this time); (now using both linear fit and relative area)

% notes: watershed relative areas don't always work too well (as opposed to fits) (e.g. winchester,marlow
% significantly overestimates flux in winchester - winchester/pony
% significantly overestimates flux in pony)
pony=min(win.*0.13 + 0.1,win.*(17/24.68)); % same thing, using both linear fit or ratio
% now can either scale by winchester or pony. these areas are in km^2
ar.win=24.68; ar.elliot=7.99; ar.north=33; ar.haynes=28; ar.kentuck=34; ar.willanch=20; ar.catching=40; ar.isthmus=22; ar.coalbank=10; ar. pony=17; ar.joekney=3;

% make Flux matrix ---- use watershed areas to estimate fluxes for ungauged
clear flux; flux=zeros(length(river_data(1).time),15);
flux(:,1) = river_data(1).qr + river_data(3).qr + river_data(6).qr;   % north fork coos river= WF Mill + EF Mill + Marlow
flux(:,2) =river_data(5).qr; % south fork coos river
flux(:,9) = win;            % Winchester Creek  % 9
flux(:,5) = win.*0.39;      % Coalbank Slough % 5
flux(:,4) = win.*0.84;      % Isthmus Slough  % 4
flux(:,3) = win.*1.5;       % Catching Slough  %3
flux(:,11) =pony.*1.96;     % North Slough  % 11
flux(:,10) = win.*1.07;     % Palouse (Haynes in model) % 10
flux(:,8) = win.*0.29;      % Elliott Creek % 8
flux(:,7) = win.*0.1;       % Joe Ney Creek  % 7
flux(:,12) = win.*1.29;     % Kentuck Inlet 12
flux(:,13) = win.*0.76;     % Willanch Creek 13
flux(:,6) = pony; % Pony Slough   % 6
flux(:,14)= win.*(9.83/ar.win).*0.5; % john b and talbot creek, have total area, so split in two
flux(:,15)= win.*(9.83/ar.win).*0.5; % john b and talbot creek, have total area, so split in two

% make flux, salt, temp be Ntimes x Nrivers matrix
%load('~/Dropbox/model/forcing_files/Q_TEMP_interp_CWA.mat')
for n=1:15
    nanx = isnan(flux(:,n));
    t= 1:numel(flux(:,n));
    flux(nanx,n) = interp1(t(~nanx), flux(~nanx,n), t(nanx));
end
salt = zeros(size(flux)); temp=zeros(size(flux)); % not using temp in model
MJD=time-datenum('Nov 17, 1858,00:00');
nTimes = numel(time);

% temp
%for i=1:4; temp(:,i) = TEMP(:,1);end % temps for coos river, using also for catching, isthmus
%load('~/Documents/coos_bay/obs_data/south slough 2007-2017.mat','winchester')
 %K=dsearchn(winchester.time,TD);
%for i=5:13;  t=nanfastsmooth(winchester.temp,10000); temp(:,i)=t(K);end %dif temp signal from winchester compared with larger (not as cold winter, and warmer summer), using smoothed 
% summer signal from winchester for all smaller creeks (?)

for i=1:15
rname=char(riv(i).name);
    rivfile=[rivfile_base(1:end-3) num2str(i) '.nc'];
    write_FVCOM_river_das(Mobj,rivfile,rname,i,1,MJD,flux(:,i),temp(:,i),salt(:,i),'Coos Freshwater Fluxes','source: CWA + estimated')
end

% report river nodes
fprintf('these rivers have been added at the following nodes: \n');
for i=1:Mobj.nRivers
  disp( char(riv(i).name) );
  fprintf('%f\n',Mobj.riv_nodes(i,1:Mobj.nRivNodes(1)))
end;
else
    disp(['skipped writing river file(s) *** copy from elsewhere!']);
end

%--------------------------------------------------------------------
% wind forcing 
%--------------------------------------------------------------------
if(do_wind)
% --------------------- original forcing --------------
% use spatially constant, but time varying based on Port Orford met station, but adjusted for Coos Bay
%  -- very well correlated with SS met station, but twice the magnitude
% wnddir = '~/Dropbox/Ted/Obs_Data/Met_data/Port_Orford_Buoy/';
% load([wnddir,'buoy_portorford_2012_2014'],'td','u','v');
% subs = 4; 
% time_wnd = td(1:subs:end) - datenum('Nov 17, 1858,00:00'); % originally every hour, make every 4th hour
% u10 = 0.48.*u(1:subs:end); v10 = 0.48.*v(1:subs:end);
% clear u v td
% write_FVCOM_wind_ts_speed_das(Mobj, wndfile, time_wnd, u10, v10)

% ---------- use charleston met station instead of Port Orford ------------
% (puts value on every grid point)
%subs = 2; % make 1/2 hourly
%load('~/Dropbox/CoosBay/model/forcing_files/charleston_met_data.mat')
%chmet.time=chmet.time(1:subs:end); chmet.u=chmet.u(1:subs:end); chmet.v=chmet.v(1:subs:end);
%ind1=dsearchn(chmet.time,datenum(2013,12,30)); ind2=dsearchn(chmet.time,datenum(2015,1,1)); % do only 2014
%chmet.time=chmet.time(ind1:ind2); chmet.u=chmet.u(ind1:ind2); chmet.v=chmet.v(ind1:ind2); 
%chmet.time=chmet.time-datenum('Nov 17, 1858,00:00');
%write_FVCOM_wind(Mobj, wndfile, chmet.time,chmet.u,chmet.v) 
%clear chmet

%-------------- use northbend airport met data ---------------------------
addpath ~/Dropbox/CoosBay/model/high/wind % thanks to efe and dkr
%load('~/Dropbox/CoosBay/model/forcing_files/northbend_airport_MET')
%cut1=dsearchn(airport.time,datenum(2013,9,1)); cut2=dsearchn(airport.time,datenum(2015,1,1));
load('~/Dropbox/CoosBay/model/high/wind/nbwind')
wspd=fillmissing(wspd,'linear');
wdir=fillmissing(wdir,'nearest');
u = wspd.*sin(pi*wdir/180);
v = wspd.*cos(pi*wdir/180);
write_windnc_uniform(dn,u,v,'/Users/ted/Dropbox/CoosBay/model/high/runfiles/coos_wind_2013_2015.nc','speed') 
else
    disp(['skipped writing ',wndfile,' *** copy from elsewhere!']);
end
                      
 %--------------------------------------------------------------------
% surface heat flux
%--------------------------------------------------------------------
% spatially constant, time varying forcing from ncep/reanalysisII
load('~/Dropbox/model/forcing_files/ncep_heat_flux/ncep_surface.mat')
data.nshf=ncep.swhf; data.lon=Mobj.lon; data.lat=Mobj.lat; data.x=Mobj.x; data.y=Mobj.y; 
data.rhum=repmat(ncep.rhum,1,92004); data.air=repmat(ncep.air,1,92004); data.slp=repmat(ncep.slp,1,92004); data.Times=ncep.time; data.time=ncep.time- datenum('Nov 17, 1858,00:00');
addpath('forcing_files/ncep_heat_flux/')
write_heat(Mobj,'coos_heat',data,'','3.1.6')
clear ncep data
%--------------------------------------------------------------------
% IC file
%--------------------------------------------------------------------
% pick timespan to search for CTD data and time of IC; 
if(do_IC)
 timespan=[datenum(2012,11,2) datenum(2012,11,5)];
 timeIC = datenum(2012,10,10,0,0,0);
 polyFile = '~/Dropbox/Ted/Model_Maker/CBAYpolys.mat';
 MobjFile = '~/Dropbox/Ted/Model_Maker/coosbay_SMS_Mobj_riveredit.mat';
 % run IC code
 Z_CBAY_ini_polys(timespan, timeIC, inifile, polyFile, MobjFile)
 % now do rivers
 rivsFile = '~/Dropbox/Ted/Model_Maker/riverpolys_CBAY.mat';
 RivTempFile='~/Dropbox/Ted/Obs_Data/Rivers/Q_TEMP_interp_CWA.mat';
 % run rivs IC code
 Z_CBAY_ini_rivs(inifile, rivsFile, MobjFile, RivTempFile)
end
% plot initial Salinity
if(do_IC)
 ss = nc_varget(inifile,'ssl');
 plot_field(Mobj,ss(1,:),'title','S initial condition');
 hold on; grid on;
end

% report estimate minimum time step
Mobj = estimate_ts(Mobj);
fprintf('estimated minimum time step in seconds: %f\n',min(Mobj.ts));
% plot everything
plot_field(Mobj,Mobj.h,'title','domain','withextra',true,'showgrid',true);
dump_dtascii([modeldir 'input/check'],Mobj);

%--------------------------------------------------------------------
% Boundary conditions from Regional ROMS (P. MacCready)
%--------------------------------------------------------------------
cd ~/Dropbox/CoosBay/model/forcing_files/LiveOcean_ROMS/2014
get_bc; cd ~/Dropbox/CoosBay/model/high
bc_dir='~/Dropbox/CoosBay/model/forcing_files/roms_bc/';
sigmafile='/Users/ted/Dropbox/CoosBay/model/high/runfiles/run_8_27_2017/input/sigma.dat';
Mobj.read_obc_nodes={1:52};
Mobj = read_sigma(Mobj, sigmafile)
bc_mjd_time=Mobj.bc_time-datenum('Nov 17, 1858,00:00'); % need to make time mjd
modeldir = '/Users/ted/Dropbox/CoosBay/model/high/runfiles/run_3_2018/coos'; % needs full dir where other run files are located
write_FVCOM_tsobc(modeldir,bc_mjd_time,20,Mobj.bc_temp,Mobj.bc_salt);
% give salt value
date_2_begin='Jan 1, 2014,00:00';
disp(['surface salinity from bc is ' num2str(Mobj.bc_salt(1,1,dsearchn(Mobj.bc_time,datenum(date_2_begin))))])
% advective transports at open boundary --haven't tried yet
%write_FVCOM_meanflow(Mobj, 'meanflow_obc_test',meanflow_data);

% 2017 B.C. from LiveOcean
addpath('~/Dropbox/CoosBay/model/forcing_files/LiveOcean_ROMS/2017'); load LiveOcean_2017
bc_dir='~/Dropbox/CoosBay/model/forcing_files/roms_bc/';
sigmafile='/Users/ted/Dropbox/CoosBay/model/high/runfiles/run_8_27_2017/input/sigma.dat';
Mobj.read_obc_nodes={1:52}; Mobj = read_sigma(Mobj, sigmafile); 
bc_mjd_time=time-datenum('Nov 17, 1858,00:00'); % need to make time mjd
Mobj.bc_salt=reshape(repmat(points(3).salt(:,3,end),1,52,20),[52,20,length(time)]);; % node,depth,time
Mobj.bc_temp=reshape(repmat(points(3).temp(:,3,end),1,52,20),[52,20,length(time)]);; % node,depth,time
modeldir = '/Users/ted/Dropbox/CoosBay/model/high/runfiles/run_fall2018_tc/coos'; % needs full dir where other run files are located
write_FVCOM_tsobc(modeldir,bc_mjd_time,20,Mobj.bc_temp,Mobj.bc_salt);
% advective transports at open boundary --haven't tried yet
%write_FVCOM_meanflow(Mobj, 'meanflow_obc_test',meanflow_data);

%% surface heat flux (not sucessfully used)
data.time=Mobj.bc_time; data.nshf=nshf; data.Times=[]; data.slp=[]; data.rhum=[]; data.air=[];
addpath('~/Dropbox/model/forcing_files/ncep_heat_flux/')
write_heat(Mobj,'coos_heat',data,'','3.1.6')
