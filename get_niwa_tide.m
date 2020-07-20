function [time,water_level]=get_niwa_tide(long,lat,datum,output_interval,start_date,number_days,api_key)
% Get tidal elevation from niwa tide model api
% https://developer.niwa.co.nz/docs/tide-api/1/routes/data.csv/get
% tc april2020
% *all input parameters should be strings*
%
% input parameters:
% lon/lat in NZGD coordinates.
% lon range:[-29 to -53] lat range: [160 to 180 and -175 to -180]
% datum: 'MSL' mean sea level or 'LAT' lowest astronomical tide
% output_interval (minutes) range [10-1440]
% start_date: format 'yyyy-mm-dd'
% number days: minimum 1
% api_key: create account on website, get api key
% (https://developer.niwa.co.nz/get-started)
%
% example:
% [time,water_level]=get_niwa_tide('177','-39','MSL','10','2019-01-01','31','your api key here');

% api only allows 31 day chunks, if want more than that need to combine 31
% day parts
if str2num(number_days) > 31
    loops=ceil(str2num(number_days)/31);
    day_end=datenum(start_date,29)+datenum(0,0,str2num(number_days));
    number_days='31';
else
    loops=1;
end

for n=1:loops
    output=table2array(webread(['https://api.niwa.co.nz/tides/data.csv?lat=' lat '&long=' long '&numberOfDays=' number_days ...
        '&startDate=' start_date '&datum=' datum '&interval=' output_interval '&apikey=' api_key]));
    tmp=char(output(9:end,1));
    
    e=str2double(string(char((output(9:end,2))))); % need to get better at this...
    t=datenum(str2num(tmp(:,1:4)),str2num(tmp(:,6:7)),str2num(tmp(:,9:10)),str2num(tmp(:,12:13)) ...
        ,str2num(tmp(:,15:16)),str2num(tmp(:,18:19))); %#ok<*ST2NM>
    
    if n==1
        time=t;
        water_level=e;
    else
        % store full time series
        time=[time;t]; %#ok<AGROW>
        water_level=[water_level;e]; %#ok<AGROW>
    end
    
    % if >31 days reset start_time number_days for next loop
    start_date=datestr(time(end)+datenum(0,0,1),29);
    if loops>1 && n==loops
    number_days=string(day(day_end - datenum(start_date,29)));
    end
    
end
end