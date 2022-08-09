function xfilt = godinfilt(xin, skip_ind,dt_spacing)
%-------------------------------------------------
% function that runs filter:
%     'godin': running averages of 24,24,25 hrs applied successively
%
% x_filtered = godinfilt(x_in,skip_ind,dt_spacing);
%         x_in = time series to be filtered
%          
% e.g.   x_filtered = godinfilt(x_in,0,1);
%   or if hourly: x_filtered=godinfilt(x_in);
%
% make skip_ind == 1 to not filter at all, default is to filter 
%
% dt_spacing - number is fractional denominator (1/dt_spacing) - '1' default, hourly input vector
%             '2' half hourly '6' 10 minute '12' 5 minute
%
% DAS, Oct. 2010, tc added dif. time vector inputs
%-------------------------------------------------
nanend = 0; % to make ends NaN's

if ~iscolumn(xin) % fix size
    xin=xin';
end

if(nargin<2);
    skip_ind = 0;
end

    
if(skip_ind~=1)
 [nrows,ncols] = size(xin);
 xnew = zeros(nrows,ncols);
 
 % deal with dif. time inputs, need to change filter by some factor
if(nargin>2);
 % build 24 hr filter
 filter24 = ones(24*dt_spacing,1); 
 filter24 = filter24 ./ sum(filter24);
 
 % build 25 hr filter
 filter25 = ones(25*dt_spacing,1); 
 filter25 = filter25 ./ sum(filter25);
 
 % covolve filters together, works because conv is associative
 temp_filter = conv(filter24, filter24);
 filter = conv(temp_filter, filter25);
 
 a = round(length(filter)/2);
 n = length(filter); 
% disp(['input vector ' num2str((1/dt_spacing)*60) ' minute intervals'])
    
else
 % build 24 hr filter
 filter24 = ones(24,1); 
 filter24 = filter24 ./ sum(filter24);
 
 % build 25 hr filter
 filter25 = ones(25,1); 
 filter25 = filter25 ./ sum(filter25);
 
 % covolve filters together, works because conv is associative
 temp_filter = conv(filter24, filter24);
 filter = conv(temp_filter, filter25);
 
 a = round(length(filter)/2);
 n = length(filter);
 
end

 for i=1:ncols
     temp = conv(xin(:,i), filter);
     xnew(:,i) = temp(a:a+nrows-1);
         if nanend  %make ends NaN's
               xnew(1:n,i) = nan*ones(n,1);
               xnew(nrows-n+1:nrows,i) = nan*ones(n,1);
         end
 end
 
 xfilt = xnew;
elseif(skip_ind==1)
 xfilt = xin;   
end %end skip_ind
 
%----------------------------------------------------------
%----------------------------------------------------------
