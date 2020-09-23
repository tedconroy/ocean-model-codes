function make_frc_file(input_files,frc_file,grid_file,blank,varargin)
% specify what fields to pull from Moana history file. Then write these out
% in the frc_file. ROMS will interpolate in house.
%
% after the input fields specified, need to input a string of the variables
% to write
%
% also has the option to write = 0 to all terms using blank = 1 
% note: setting shflux = 0 is bad idea
%
% From Moana project use the following variables: wind surface stress, net
% heat atm flux, set some of the others = 0
%
% Accepts any combinations of desired parameter(s):
% (You must use these specific names on the input line)
% Uwind: surface u-wind component (m/s)
% Vwind: surface v-wind component (m/s)
% Pair: surface air pressure (mbar)
% Tair: surface air temperature (Celsius)
% Qair: surface air relative humidity (%)
% rain: rain fall rate (kg/m2/s)
% swrad: solar shortwave radiation (W/m2)
% lwrad: solar longwave radiation (W/m2)
% sustr: surface u-stress (N/m2)
% svstr: surface v-stress (N/m2)


% Loop through # of files to be created
for n=1:length(input_files)
   
cur_file = input_files(n).name;
fn = [frc_file num2str(n) '.nc']; % file name to be created in this loop
time = ncread(cur_file,'ocean_time');
lon = ncread(grid_file,'lon_rho');
lat = ncread(grid_file,'lat_rho');

% create file
nc=netcdf.create(fn,'clobber');
if isempty(nc), return, end

disp(' ## Defining Global Attributes...')
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'history', ['Created by ' mfilename ' on ' datestr(now)]);
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'type', 'forcing file');

% Dimensions:
disp(' ## Defining Dimensions...')
[ix,iy]=size(lon);
t=length(time);
lon_dimID = netcdf.defDim(nc,'xrho',ix);
lat_dimID = netcdf.defDim(nc,'yrho',iy);
t_dimID = netcdf.defDim(nc,'time',t);

% Variables and attributes:
disp(' ## Defining Variables, and Attributes...')
tID = netcdf.defVar(nc,'time','double',t_dimID);
netcdf.putAtt(nc,tID,'long_name','forcing time');
netcdf.putAtt(nc,tID,'units','seconds since 2010-01-01 00:00:00');
netcdf.putAtt(nc,tID,'field','time, scalar, series');

lonID = netcdf.defVar(nc,'lon','double',[lon_dimID lat_dimID]);
netcdf.putAtt(nc,lonID,'long_name','longitude');
netcdf.putAtt(nc,lonID,'units','degrees_east');
netcdf.putAtt(nc,lonID,'field','lon_rho, scalar');

latID = netcdf.defVar(nc,'lat','double',[lon_dimID lat_dimID]);
netcdf.putAtt(nc,latID,'long_name','latitude');
netcdf.putAtt(nc,latID,'units','degrees_north');
netcdf.putAtt(nc,latID,'field','lat_rho, scalar');

%% put variable information into NC file:
% uses varagin strings to determine which vars to use

% if u,v wind defined
if sum(strcmpi(varargin,'Uwind'))>0
  wt_dimID = netcdf.defDim(nc,'wind_time',t);
  wtID = netcdf.defVar(nc,'wind_time','double',wt_dimID);
  netcdf.putAtt(nc,wtID,'long_name','wind_time');
  netcdf.putAtt(nc,wtID,'units','days');
  netcdf.putAtt(nc,wtID,'field','Uwind_time, scalar, series');
netcdf.putAtt(nc,lonID,'long_name','longitude');

  UwindID = netcdf.defVar(nc,'Uwind','double',[lon_dimID lat_dimID wt_dimID]);
  netcdf.putAtt(nc,UwindID,'long_name','surface u-wind component');
  netcdf.putAtt(nc,UwindID,'units','meter second-1');
  netcdf.putAtt(nc,UwindID,'field','Uwind, scalar, series');
  netcdf.putAtt(nc,UwindID,'coordinates','lon lat');
  netcdf.putAtt(nc,UwindID,'time','wind_time');
end
if sum(strcmpi(varargin,'Vwind'))>0
  if ~sum(strcmpi(varargin,'Uwind'))>0
    wt_dimID = netcdf.defDim(nc,'wind_time',t);
    wtID = netcdf.defVar(nc,'wind_time','double',wt_dimID);
    netcdf.putAtt(nc,wtID,'long_name','wind_time');
    netcdf.putAtt(nc,wtID,'units','days');
    netcdf.putAtt(nc,wtID,'field','Vwind_time, scalar, series');
  end
  VwindID = netcdf.defVar(nc,'Vwind','double',[lon_dimID lat_dimID wt_dimID]);
  netcdf.putAtt(nc,VwindID,'long_name','surface v-wind component');
  netcdf.putAtt(nc,VwindID,'units','meter second-1');
  netcdf.putAtt(nc,VwindID,'field','Vwind, scalar, series');
  netcdf.putAtt(nc,VwindID,'coordinates','lon lat');
  netcdf.putAtt(nc,VwindID,'time','wind_time');
end

% air
if sum(strcmpi(varargin,'Pair'))>0
  Pat_dimID = netcdf.defDim(nc,'Pair_time',t);
  PatID = netcdf.defVar(nc,'Pair_time','double',Pat_dimID);
  netcdf.putAtt(nc,PatID,'long_name','Pair_time');
  netcdf.putAtt(nc,PatID,'units','days');
  netcdf.putAtt(nc,PatID,'field','Pair_time, scalar, series');

  PairID = netcdf.defVar(nc,'Pair','double',[lon_dimID lat_dimID Pat_dimID]);
  netcdf.putAtt(nc,PairID,'long_name','surface air pressure');
  netcdf.putAtt(nc,PairID,'units','millibar');
  netcdf.putAtt(nc,PairID,'field','Pair, scalar, series');
  netcdf.putAtt(nc,PairID,'coordinates','lon lat');
  netcdf.putAtt(nc,PairID,'time','Pair_time');
end
if sum(strcmpi(varargin,'Tair'))>0
  Tat_dimID = netcdf.defDim(nc,'Tair_time',t);
  TatID = netcdf.defVar(nc,'Tair_time','double',Tat_dimID);
  netcdf.putAtt(nc,TatID,'long_name','Tair_time');
  netcdf.putAtt(nc,TatID,'units','days');
  netcdf.putAtt(nc,TatID,'field','Tair_time, scalar, series');

  TairID = netcdf.defVar(nc,'Tair','double',[lon_dimID lat_dimID Tat_dimID]);
  netcdf.putAtt(nc,TairID,'long_name','surface air temperature');
  netcdf.putAtt(nc,TairID,'units','Celsius');
  netcdf.putAtt(nc,TairID,'field','Tair, scalar, series');
  netcdf.putAtt(nc,TairID,'coordinates','lon lat');
  netcdf.putAtt(nc,TairID,'time','Tair_time');
end
if sum(strcmpi(varargin,'Qair'))>0
  Qat_dimID = netcdf.defDim(nc,'Qair_time',t);
  QatID = netcdf.defVar(nc,'Qair_time','double',Qat_dimID);
  netcdf.putAtt(nc,QatID,'long_name','Qair_time');
  netcdf.putAtt(nc,QatID,'units','days');
  netcdf.putAtt(nc,QatID,'field','Qair_time, scalar, series');
  QairID = netcdf.defVar(nc,'Qair','double',[lon_dimID lat_dimID Qat_dimID]);
  netcdf.putAtt(nc,QairID,'long_name','surface air relative humidity');
  netcdf.putAtt(nc,QairID,'units','percentage');
  netcdf.putAtt(nc,QairID,'field','Qair, scalar, series');
  netcdf.putAtt(nc,QairID,'coordinates','lon lat');
  netcdf.putAtt(nc,QairID,'time','Qair_time');
end

% rain
if sum(strcmpi(varargin,'rain'))>0
  rt_dimID = netcdf.defDim(nc,'rain_time',t);
  rtID = netcdf.defVar(nc,'rain_time','double',rt_dimID);
  netcdf.putAtt(nc,rtID,'long_name','rain_time');
  netcdf.putAtt(nc,rtID,'units','days');
  netcdf.putAtt(nc,rtID,'field','rain_time, scalar, series');

  rainID = netcdf.defVar(nc,'rain','double',[lon_dimID lat_dimID rt_dimID]);
  netcdf.putAtt(nc,rainID,'long_name','rain fall rate');
  netcdf.putAtt(nc,rainID,'units','kilogram meter-2 second-1');
  netcdf.putAtt(nc,rainID,'field','rain, scalar, series');
  netcdf.putAtt(nc,rainID,'coordinates','lon lat');
  netcdf.putAtt(nc,rainID,'time','rain_time');
end

%  fluxes
if sum(strcmpi(varargin,'swrad'))>0
  swrt_dimID = netcdf.defDim(nc,'swrad_time',t);
  swrtID = netcdf.defVar(nc,'swrad_time','double',swrt_dimID);
  netcdf.putAtt(nc,swrtID,'long_name','swrad_time');
  netcdf.putAtt(nc,swrtID,'units','days');
  netcdf.putAtt(nc,swrtID,'field','swrad_time, scalar, series');
  swradID = netcdf.defVar(nc,'swrad','double',[lon_dimID lat_dimID swrt_dimID]);
  netcdf.putAtt(nc,swradID,'long_name','solar shortwave radiation');
  netcdf.putAtt(nc,swradID,'units','Watts meter-2');
  netcdf.putAtt(nc,swradID,'positive_value','downward flux, heating');
  netcdf.putAtt(nc,swradID,'negative_value','upward flux, cooling');
  netcdf.putAtt(nc,swradID,'field','swrad, scalar, series');
  netcdf.putAtt(nc,swradID,'coordinates','lon lat');
  netcdf.putAtt(nc,swradID,'time','swrad_time');
end
if sum(strcmpi(varargin,'lwrad'))>0
  lwrt_dimID = netcdf.defDim(nc,'lwrad_time',t);
  lwrtID = netcdf.defVar(nc,'lwrad_time','double',lwrt_dimID);
  netcdf.putAtt(nc,lwrtID,'long_name','lwrad_time');
  netcdf.putAtt(nc,lwrtID,'units','days');
  netcdf.putAtt(nc,lwrtID,'field','lwrad_time, scalar, series');
  lwradID = netcdf.defVar(nc,'lwrad','double',[lon_dimID lat_dimID lwrt_dimID]);
  netcdf.putAtt(nc,lwradID,'long_name','solar longwave radiation');
  netcdf.putAtt(nc,lwradID,'units','Watts meter-2');
  netcdf.putAtt(nc,lwradID,'positive_value','downward flux, heating');
  netcdf.putAtt(nc,lwradID,'negative_value','upward flux, cooling');
  netcdf.putAtt(nc,lwradID,'field','lwrad, scalar, series');
  netcdf.putAtt(nc,lwradID,'coordinates','lon lat');
  netcdf.putAtt(nc,lwradID,'time','lwrad_time');
end
if sum(strcmpi(varargin,'shflux'))>0
  lwrt_dimID = netcdf.defDim(nc,'shflux_time',t);
  lwrtID = netcdf.defVar(nc,'shflux_time','double',lwrt_dimID);
  netcdf.putAtt(nc,lwrtID,'long_name','shflux_time');
  netcdf.putAtt(nc,lwrtID,'units','seconds since 2010-01-01 00:00:00');
  netcdf.putAtt(nc,lwrtID,'field','shflux_time, scalar, series');
  lwradID = netcdf.defVar(nc,'shflux','double',[lon_dimID lat_dimID lwrt_dimID]);
  netcdf.putAtt(nc,lwradID,'long_name','net surface heat flux');
  netcdf.putAtt(nc,lwradID,'units','Watts meter-2');
  netcdf.putAtt(nc,lwradID,'positive_value','downward flux, heating');
  netcdf.putAtt(nc,lwradID,'negative_value','upward flux, cooling');
  netcdf.putAtt(nc,lwradID,'field','shflux, scalar, series');
  netcdf.putAtt(nc,lwradID,'coordinates','lon lat');
  netcdf.putAtt(nc,lwradID,'time','shflux_time');
end
if sum(strcmpi(varargin,'bhflux'))>0
  lwrt_dimID = netcdf.defDim(nc,'bhflux_time',t);
  lwrtID = netcdf.defVar(nc,'bhflux_time','double',lwrt_dimID);
  netcdf.putAtt(nc,lwrtID,'long_name','bhflux_time');
  netcdf.putAtt(nc,lwrtID,'units','seconds since 2010-01-01 00:00:00');
  netcdf.putAtt(nc,lwrtID,'field','bhflux_time, scalar, series');
  lwradID = netcdf.defVar(nc,'bhflux','double',[lon_dimID lat_dimID lwrt_dimID]);
  netcdf.putAtt(nc,lwradID,'long_name','net bottom heat flux');
  netcdf.putAtt(nc,lwradID,'units','Watts meter-2');
  netcdf.putAtt(nc,lwradID,'positive_value','downward flux, heating');
  netcdf.putAtt(nc,lwradID,'negative_value','upward flux, cooling');
  netcdf.putAtt(nc,lwradID,'field','bhflux, scalar, series');
  netcdf.putAtt(nc,lwradID,'coordinates','lon lat');
  netcdf.putAtt(nc,lwradID,'time','bhflux_time');
end
if sum(strcmpi(varargin,'ssflux'))>0
  lwrt_dimID = netcdf.defDim(nc,'ssflux_time',t);
  lwrtID = netcdf.defVar(nc,'ssflux_time','double',lwrt_dimID);
  netcdf.putAtt(nc,lwrtID,'long_name','ssflux_time');
  netcdf.putAtt(nc,lwrtID,'units','seconds since 2010-01-01 00:00:00');
  netcdf.putAtt(nc,lwrtID,'field','ssflux_time, scalar, series');
  lwradID = netcdf.defVar(nc,'ssflux','double',[lon_dimID lat_dimID lwrt_dimID]);
  netcdf.putAtt(nc,lwradID,'long_name','net surface salt flux');
  netcdf.putAtt(nc,lwradID,'units','cm / day');
  netcdf.putAtt(nc,lwradID,'positive_value','flux');
  netcdf.putAtt(nc,lwradID,'negative_value','flux');
  netcdf.putAtt(nc,lwradID,'field','ssflux, scalar, series');
  netcdf.putAtt(nc,lwradID,'coordinates','lon lat');
  netcdf.putAtt(nc,lwradID,'time','ssflux_time');
end
if sum(strcmpi(varargin,'bwflux'))>0
  lwrt_dimID = netcdf.defDim(nc,'bwflux_time',t);
  lwrtID = netcdf.defVar(nc,'bwflux_time','double',lwrt_dimID);
  netcdf.putAtt(nc,lwrtID,'long_name','bwflux_time');
  netcdf.putAtt(nc,lwrtID,'units','seconds since 2010-01-01 00:00:00');
  netcdf.putAtt(nc,lwrtID,'field','bwflux_time, scalar, series');
  lwradID = netcdf.defVar(nc,'bwflux','double',[lon_dimID lat_dimID lwrt_dimID]);
  netcdf.putAtt(nc,lwradID,'long_name','net bottom freshwater flux');
  netcdf.putAtt(nc,lwradID,'units','cm / day');
  netcdf.putAtt(nc,lwradID,'positive_value','flux');
  netcdf.putAtt(nc,lwradID,'negative_value','flux');
  netcdf.putAtt(nc,lwradID,'field','bwflux, scalar, series');
  netcdf.putAtt(nc,lwradID,'coordinates','lon lat');
  netcdf.putAtt(nc,lwradID,'time','bwflux_time');
end

% surface stresses
if sum(strcmpi(varargin,'sustr'))>0
  su_dimID = netcdf.defDim(nc,'sustr_time',t);
  suID = netcdf.defVar(nc,'sustr_time','double',su_dimID);
  netcdf.putAtt(nc,suID,'long_name','sustr_time');
  netcdf.putAtt(nc,suID,'units','seconds since 2010-01-01 00:00:00');
  netcdf.putAtt(nc,suID,'field','sustr_time, scalar, series');
  sustrID = netcdf.defVar(nc,'sustr','double',[lon_dimID lat_dimID su_dimID]);
  netcdf.putAtt(nc,sustrID,'long_name','surface u-stress');
  netcdf.putAtt(nc,sustrID,'units','Newtons meter-2');
  netcdf.putAtt(nc,sustrID,'field','sustr, scalar, series');
  netcdf.putAtt(nc,sustrID,'coordinates','lon lat');
  netcdf.putAtt(nc,sustrID,'time','sustr_time');
end
if sum(strcmpi(varargin,'svstr'))>0
  sv_dimID = netcdf.defDim(nc,'svstr_time',t);
  svID = netcdf.defVar(nc,'svstr_time','double',sv_dimID);
  netcdf.putAtt(nc,svID,'long_name','svstr_time');
  netcdf.putAtt(nc,svID,'units','seconds since 2010-01-01 00:00:00');
  netcdf.putAtt(nc,svID,'field','svstr_time, scalar, series');
  svstrID = netcdf.defVar(nc,'svstr','double',[lon_dimID lat_dimID sv_dimID]);
  netcdf.putAtt(nc,svstrID,'long_name','surface v-stress');
  netcdf.putAtt(nc,svstrID,'units','Newtons meter-2');
  netcdf.putAtt(nc,svstrID,'field','svstr, scalar, series');
  netcdf.putAtt(nc,svstrID,'coordinates','lon lat');
  netcdf.putAtt(nc,svstrID,'time','svstr_time');
end

netcdf.close(nc)

%% open the file for writing
nc=netcdf.open(fn,'NC_WRITE');

% write time
ID=netcdf.inqVarID(nc,'time');
netcdf.putVar(nc,ID,time);

% write lon and lat
ID=netcdf.inqVarID(nc,'lon');
netcdf.putVar(nc,ID,lon);
ID=netcdf.inqVarID(nc,'lat');
netcdf.putVar(nc,ID,lat);

if sum(strcmpi(varargin,'Uwind'))>0
  ID=netcdf.inqVarID(nc,'wind_time');
  netcdf.putVar(nc,ID,time);
  if blank
      Uwind = zeros(size(lon));
  else
  Uwind=ncread(cur_file,'Uwind');
  end
  ID=netcdf.inqVarID(nc,'Uwind');
  netcdf.putVar(nc,ID,Uwind);
end
if sum(strcmpi(varargin,'Vwind'))>0
  ID=netcdf.inqVarID(nc,'wind_time');
  netcdf.putVar(nc,ID,time);
  if blank
      Vwind = zeros(size(lon));
  else
      Vwind = ncread(curfile,'Vwind');
  end
  ID=netcdf.inqVarID(nc,'Vwind');
  netcdf.putVar(nc,ID,Vwind);
end
if sum(strcmpi(varargin,'Pair'))>0
  ID=netcdf.inqVarID(nc,'Pair_time');
  netcdf.putVar(nc,ID,time);
  Pair=evalin('base','Pair');
  ID=netcdf.inqVarID(nc,'Pair');
  netcdf.putVar(nc,ID,Pair);
end
if sum(strcmpi(varargin,'Tair'))>0
  ID=netcdf.inqVarID(nc,'Tair_time');
  netcdf.putVar(nc,ID,time);
  Tair=evalin('base','Tair');
  ID=netcdf.inqVarID(nc,'Tair');
  netcdf.putVar(nc,ID,Tair);
end
if sum(strcmpi(varargin,'Qair'))>0
  ID=netcdf.inqVarID(nc,'Qair_time');
  netcdf.putVar(nc,ID,time);
  Qair=evalin('base','Qair');
  ID=netcdf.inqVarID(nc,'Qair');
  netcdf.putVar(nc,ID,Qair);
end
if sum(strcmpi(varargin,'rain'))>0
  ID=netcdf.inqVarID(nc,'rain_time');
  netcdf.putVar(nc,ID,time);
  rain=evalin('base','rain');
  ID=netcdf.inqVarID(nc,'rain');
  netcdf.putVar(nc,ID,rain);
end
if sum(strcmpi(varargin,'swrad'))>0
  ID=netcdf.inqVarID(nc,'swrad_time');
  netcdf.putVar(nc,ID,time);
  swrad=evalin('base','swrad');
  ID=netcdf.inqVarID(nc,'swrad');
  netcdf.putVar(nc,ID,swrad);
end
if sum(strcmpi(varargin,'lwrad'))>0
  ID=netcdf.inqVarID(nc,'lwrad_time');
  netcdf.putVar(nc,ID,time);
  lwrad=evalin('base','lwrad');
  ID=netcdf.inqVarID(nc,'lwrad');
  netcdf.putVar(nc,ID,lwrad);
end
if sum(strcmpi(varargin,'shflux'))>0
  ID=netcdf.inqVarID(nc,'shflux_time');
  netcdf.putVar(nc,ID,time);
  shflux = ncread(fn,'shflux');
  end
  ID=netcdf.inqVarID(nc,'shflux');
  netcdf.putVar(nc,ID,shflux);
end
if sum(strcmpi(varargin,'bhflux'))>0
  ID=netcdf.inqVarID(nc,'bhflux_time');
  netcdf.putVar(nc,ID,time);
  if blank
      tmp = ncread(fn,'bhflux');
      shflux = zeros(size(tmp));
  else
      shflux = ncread(cur_file,'bhflux');
  end
  ID=netcdf.inqVarID(nc,'bhflux');
  netcdf.putVar(nc,ID,shflux);
end
if sum(strcmpi(varargin,'ssflux'))>0
  ID=netcdf.inqVarID(nc,'ssflux_time');
  netcdf.putVar(nc,ID,time);
  if blank
      tmp = ncread(fn,'ssflux');
      ssflux = zeros(size(tmp));
  else
      ssflux = ncread(cur_file,'ssflux');
  end
  ID=netcdf.inqVarID(nc,'ssflux');
  netcdf.putVar(nc,ID,ssflux);
end
if sum(strcmpi(varargin,'bwflux'))>0
  ID=netcdf.inqVarID(nc,'bwflux_time');
  netcdf.putVar(nc,ID,time);
  if blank
      tmp = ncread(fn,'bwflux');
      bwflux = zeros(size(tmp));
  else
      bwflux = ncread(cur_file,'bwflux');
  end
  ID=netcdf.inqVarID(nc,'bwflux');
  netcdf.putVar(nc,ID,bwflux);
end
if sum(strcmpi(varargin,'sustr'))>0
  ID=netcdf.inqVarID(nc,'sustr_time');
  netcdf.putVar(nc,ID,time);
  if blank
      tmp = ncread(fn,'sustr');
      sustr = zeros(size(tmp));
  else
      sustr = ncread(cur_file,'sustr');
  end
  ID=netcdf.inqVarID(nc,'sustr');
  netcdf.putVar(nc,ID,sustr);
end
if sum(strcmpi(varargin,'svstr'))>0
  ID=netcdf.inqVarID(nc,'svstr_time');
  netcdf.putVar(nc,ID,time);
  if blank 
      tmp = ncread(fn,'svstr');
      svstr = zeros(size(tmp));
  else
      svstr = ncread(cur_file,'svstr');
  end
  ID=netcdf.inqVarID(nc,'svstr');
  netcdf.putVar(nc,ID,svstr);
end
netcdf.close(nc)
disp(['done ' num2str(n)])

end
