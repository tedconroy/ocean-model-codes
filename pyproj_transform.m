function [xout,yout] = pyproj_transform(x,y, incord,outcord)
% use python pyproj library
% [xout,yout] = pyproj_trans(x,y, incord,outcord)
% in/out cord is integer epsg number
%
% to use, first need to install pyproj module on command line. for windows:
% pip install pyproj
% then in matlab:
% py.importlib.import_module('pyproj')
%
% example:
% [xout,yout] = pyproj_trans(lon,lat,3994,2193);
%
% some common epsg codes:
% nztm 2193
% wgs84 4326
% niwa wgs84 3994

inProj   = py.pyproj.Proj(pyargs('init',['epsg:'  int2str((incord))]));
outProj  = py.pyproj.Proj(pyargs('init',['epsg:'  int2str((outcord))]));

a= py.pyproj.transform(inProj, outProj, x(:)', y(:)');

% restore origonal shape and type
s=size(x);
xout= reshape(double(a{1}),s);
yout = reshape(double(a{2}),s);
    
end
  
