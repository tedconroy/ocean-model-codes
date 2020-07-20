function [flux]=calcFluxSect(fn,xysec,dn12,strin,nstep,sedname)
%   flux=calcFluxSect(fn,xysec,dn12,strin,nstep,sedname)
% 
% calculate flux through a section 
% based on code from Charlie Flagg (transport_calc_6.m in skagit/fvcom)
% written by David Ralston (WHOI), modified Ted Conroy
%
%   fn=input netcdf file (or object from mDataset)
%
%   xysec = section location
%       --> to define start and end of section (assumes straight line)
%      xysec = [x0 y0 x1 y1];
%       	-> can have xysec(nsect,4) for multiple sections
%       --> to define a path of points for the section
%      xysec{1:nsec}[x0 y0; x1 y1; x2 y2; ....]; 
%           -> can have xysec{nsect} cells for multiple sections
%
%   dn12 = beginning and end time of flux calc
%
%   strin = string w/ what to calc 
%      'v'=volume
%      's'=salinity
%      't'=temperature
%      'f'=freshwater (Smax-s)/Smax
%      'd'=sediment (for dirt)
%      'p'=plot x-sect locations
%      'x'=returns integrated & 2d structure of flux
%      
%
%   nstep = time steps to calculate [default=1]
%   
%   sedname = cells w/ names of sediment types [optional]


if nargin<5
    nstep=1;
end

if strfind(strin,'d') & nargin < 6;
    error('define sednames in input');
end

dn00=datenum(1858,11,17,0,0,0);  % start of modified julian day
nc=ncgeodataset(fn);
xn=double(nc.data('x')');  % 0.3 s to do this snippet
yn=double(nc.data('y')'); 
nv=double(nc.data('nv')); 
siglev=double(nc.data('siglev')); 
nbe=double(nc.data('nbe')); 
hn=double(nc.data('h')'); 
dn=double(nc.data('time'))+dn00; 

M=length(xn);
N=length(nv);
KB = size(siglev,1);
KBM1 = KB - 1; 
xc=mean(xn(nv));
yc=mean(yn(nv));

% This portion determines the section and the indices for the cells and TCE segments    
if iscell(xysec);
    nsect = length(xysec);
    for n=1:nsect
        xsect{n} = xysec{n}(:,1);
        ysect{n} = xysec{n}(:,2);
    end    
else
    nsect = size(xysec,1);
    for n=1:nsect
        xsect{n} = xysec(n,[1 3]);
        ysect{n} = xysec(n,[2 4]);
    end
end    

if strfind(strin,'p');
    figure(1); clf;
    plot(xn,yn,'.'); axis equal; hold on;    
end

for n=1:nsect        
    xs0=xsect{n}(:); ys0=ysect{n}(:);
    ns0 = findNodes(xc,yc,xs0,ys0);
    ietmp=[]; dastmp=[];
    for i=2:length(ns0);
        nnc=shortestPath(xc,yc,nbe,ns0(i-1),ns0(i));  
        nnc=nnc(~ismember(nnc,ietmp));
        ietmp=[ietmp nnc];
        % project waypoints onto lines
        ltmp= createLine([xc(ns0(i-1)) yc(ns0(i-1))],[xc(ns0(i)) yc(ns0(i))]);
        dtmp=sqrt(ltmp(3)^2+ltmp(4)^2);
        dstmp = linePosition([xc(nnc)' yc(nnc)'],ltmp)'.*dtmp;    % dist along line (linePosition=normalized)    
        if i==2; dastmp=dstmp; else dastmp=[dastmp dstmp+max(dastmp)]; end
    end
    ie{n}=double(ietmp);
    ncells(n)=length(ie{n});
    xas{n}=xc(ie{n});
    yas{n}=yc(ie{n});
    das{n}=dastmp;
        
    % Now find the information for the 2*ncells TCE segments
    ibn1 = nv(find(nbe(:,ie{n}(1)) ~= 0),ie{n}(1));   % Boundary nodes for the first boundary cell
    ibnn = nv(find(nbe(:,ie{n}(end)) ~= 0),ie{n}(end));   % Boundary nodes for the ending boundary cell
    xm{n}(1) = mean(xn(ibn1)); ym{n}(1) = mean(yn(ibn1));                 % mid-point of first boundary cell edge
    hnm{n}(1) = mean(hn(ibn1));                                      % Undisturbed water depth at mid-point
    xm{n}(ncells(n)+1) = mean(xn(ibnn)); ym{n}(ncells(n)+1) = mean(yn(ibnn));   % mid-point of last boundary cell edge
    hnm{n}(ncells(n)+1) = mean(hn(ibnn));                               % Undisturbed water depth at mid-point
    for i = 2:ncells(n)
        icns = nv(find([1:3]~=find(nbe(:,ie{n}(i-1))==ie{n}(i))),ie{n}(i-1));    % common nodes
        xm{n}(i) = mean(xn(icns)); ym{n}(i) = mean(yn(icns));                   % mid-point of the common cell edge
        hnm{n}(i) = mean(hn(icns));                                        % Undisturbed water depth at mid-point
    end
    hc{n} = mean(hn(nv(:,ie{n})));                                     % Undisturbed water depths at cell centers

    if strfind(strin,'p');
       plot(xc(ie{n}),yc(ie{n}),'r.-');
       text(mean(double(xc(ie{n}(1)))),mean(double(yc(ie{n}(1)))),num2str(n),'color','r','fontsize',14);
    end
end % setup section loop


% Now calculate the transports through the TCE segments for each time step
if isempty(dn12);
    nt1=1;
    nt2=length(dn);
else
    [jnk,nt1]=min(abs(dn-dn12(1)));
    [jnk,nt2]=min(abs(dn-dn12(2)));
end
flux.dn=dn(nt1:nstep:nt2);
if strfind(strin,'x');
    flux.x = xas;
    flux.y = yas;
    flux.d = das;
end
idt=0;

for it = nt1:nstep:nt2
    idt = idt+1;
    
    % get vars
    u=sqz(nc.data('u',[it 1 1],[it KBM1 N]));
    v=sqz(nc.data('v',[it 1 1],[it KBM1 N]));
    el=nc.data('zeta',[it 1],[it M]);
    %wet=double(nc.data('wet_nodes')); % 0=dry 1=wet
    %dry=find(wet==0);el_wnan=el; el_wnan(dry)=NaN; 
    
    if ~isempty(strfind(strin,'s')) || ~isempty(strfind(strin,'f'));
        sal=sqz(nc.data('salinity',[it 1 1],[it KBM1 M]));
    end
    if ~isempty(strfind(strin,'t'));
        tem=sqz(nc.data('temp',[it 1 1],[it KBM1 M]));
    end
    if ~isempty(strfind(strin,'d'));
        nsed=length(sedname);
        for m=1:nsed
            sed{m}=sqz(nc.data(sedname{m},[it 1 1],[it KBM1 M]));
        end
    end
        
    for n=1:nsect
        elc = mean(el(nv(:,ie{n})));        % water elevations at the cell centroids
        %eln = nanmean(el_wnan(nv(:,ie{n}))); % same but with nans where dry
        zz=calcCellZ(21,hc{n},elc); %  depth of each layer (size xsection)
        %zcord=calcCellZ(20,hc{n},elc); %  depth of each layer (size xsection) (changed from 21)
        
        if ~isempty(strfind(strin,'s')) || ~isempty(strfind(strin,'f'));
            clear salc sal_t fwf_t
            for i = 1:ncells(n)
                salc(:,i) = mean(sal(:,nv(:,ie{n}(i)))')';
            end
        end
        if ~isempty(strfind(strin,'t'));
            clear temc tem_t
            for i = 1:ncells(n)
                temc(:,i) = mean(tem(:,nv(:,ie{n}(i)))')';
            end
        end
        if ~isempty(strfind(strin,'d'));
            clear sedc sed_t
            for m=1:nsed
                for i = 1:ncells(n)
                    sedc{m}(:,i) = mean(sed{m}(:,nv(:,ie{n}(i)))')';
                end
            end
        end
        clear sigc vol_t
        for i = 1:ncells(n)
            sigc(:,i) = mean(siglev(:,nv(:,ie{n}(i)))')';
        end

        if ~isempty(strfind(strin,'f'));
            smax=max(sal(:));
        	fwfc=(smax-salc)./smax;
        end
        
        % Determine transports
         Area(n,idt) = 0;
        for i = 1:ncells(n)                      
            iseg = 1+2*(i-1);
            dx = xc(ie{n}(i))-xm{n}(i);  dy = yc(ie{n}(i))-ym{n}(i);
            TD = hc{n}(i) + elc(i);  % mean cell depth 
            area = TD*sqrt(dx^2 + dy^2);
            % for each cell want area of YZ face --- use 'Area'(integrated) to check 
            for dep=1:KBM1 % area for element in cross section direction
               A(n,idt,i,dep) = 0; % preallocate
               dvertical=abs(zz(dep+1,i)-zz(dep,i));
               dx = xc(ie{n}(i))-xm{n}(i);  dy = yc(ie{n}(i))-ym{n}(i);
               ar=dvertical*sqrt(dx^2 + dy^2);
               A(n,idt,i,dep) = A(n,idt,i,dep) + ar;
               dx = abs(xm{n}(i+1) - xc(ie{n}(i)));  dy = abs(ym{n}(i+1) - yc(ie{n}(i))); 
               ar=dvertical*sqrt(dx^2 + dy^2);
               A(n,idt,i,dep) = A(n,idt,i,dep) + ar;
               flux.area{n}(idt,i,dep)=A(n,idt,i,dep);
            end
            iseg = 1+2*(i-1);
            dx = xc(ie{n}(i))-xm{n}(i);  dy = yc(ie{n}(i))-ym{n}(i);
            TD = hc{n}(i) + elc(i);  % mean cell depth 
            area = TD*sqrt(dx^2 + dy^2);
            
            vol_t(:,iseg) = (TD*diff(-sigc(:,i))).*(u(:,ie{n}(i))*dy - v(:,ie{n}(i))*dx);
            if ~isempty(strfind(strin,'s'));
               sal_t(:,iseg) = (TD*diff(-sigc(:,i))).*(u(:,ie{n}(i)).*salc(:,i)*dy - v(:,ie{n}(i)).*salc(:,i)*dx); % psu m^3/s
            end
            if ~isempty(strfind(strin,'f'));
               fwf_t(:,iseg) = (TD*diff(-sigc(:,i))).*(u(:,ie{n}(i)).*fwfc(:,i)*dy - v(:,ie{n}(i)).*fwfc(:,i)*dx); % m^3/s
            end
            if ~isempty(strfind(strin,'t'));
               tem_t(:,iseg) = (TD*diff(-sigc(:,i))).*(u(:,ie{n}(i)).*temc(:,i)*dy - v(:,ie{n}(i)).*temc(:,i)*dx); % Deg C m^3/s
            end
            if ~isempty(strfind(strin,'d'));
                for m=1:nsed
                    sed_t{m}(:,iseg) = (TD*diff(-sigc(:,i))).*(u(:,ie{n}(i)).*sedc{m}(:,i)*dy - v(:,ie{n}(i)).*sedc{m}(:,i)*dx); % kg/s
                end
            end
             Area(n,idt) = Area(n,idt) + area;

            iseg = 2*i;
            dx = xm{n}(i+1) - xc(ie{n}(i));  dy = ym{n}(i+1) - yc(ie{n}(i));
            area = TD*sqrt(dx^2 + dy^2);
            vol_t(:,iseg) = (TD*diff(-sigc(:,i))).*(u(:,ie{n}(i))*dy - v(:,ie{n}(i))*dx);
            if ~isempty(strfind(strin,'s'));
               sal_t(:,iseg) = (TD*diff(-sigc(:,i))).*(u(:,ie{n}(i)).*salc(:,i)*dy - v(:,ie{n}(i)).*salc(:,i)*dx); % psu m^3/s
            end
            if ~isempty(strfind(strin,'f'));
               fwf_t(:,iseg) = (TD*diff(-sigc(:,i))).*(u(:,ie{n}(i)).*fwfc(:,i)*dy - v(:,ie{n}(i)).*fwfc(:,i)*dx); % m^3/s
            end
            if ~isempty(strfind(strin,'t'));
               tem_t(:,iseg) = (TD*diff(-sigc(:,i))).*(u(:,ie{n}(i)).*temc(:,i)*dy - v(:,ie{n}(i)).*temc(:,i)*dx); % Deg C m^3/s
            end
            if ~isempty(strfind(strin,'d'));
                for m=1:nsed
                    sed_t{m}(:,iseg) = (TD*diff(-sigc(:,i))).*(u(:,ie{n}(i)).*sedc{m}(:,i)*dy - v(:,ie{n}(i)).*sedc{m}(:,i)*dx); % kg/s
                end
            end
             Area(n,idt) = Area(n,idt) + area;
        end

        % Sum up transports for this time step
        if ~isempty(strfind(strin,'x'));
            flux.vol(n,idt)= nansum(vol_t(:)); % integrated
            tmp=sum(vol_t); % 2x cells
            flux.vol_x{n}(idt,:) = tmp(1:2:end-1)+tmp(2:2:end); % 1d
            flux.vol_d{n}(idt,:,:)=vol_t(:,1:2:end-1)+vol_t(:,2:2:end); %2d
        else
            flux.vol(n,idt)= sum(vol_t(:));
        end
        if ~isempty(strfind(strin,'s'));
            if ~isempty(strfind(strin,'x'));
                flux.sal(n,idt)= nansum(sal_t(:)); % integrated
                tmp=nansum(sal_t);
                flux.sal_x{n}(idt,:) = tmp(1:2:end-1)+tmp(2:2:end); % 1d
                flux.sal_d{n}(idt,:,:)=sal_t(:,1:2:end-1)+sal_t(:,2:2:end); % 2d
            else
                flux.sal(n,idt)= sum(sal_t(:));
            end
        end
        if ~isempty(strfind(strin,'f'));
            if ~isempty(strfind(strin,'x'));
                flux.fw(n,idt)= nansum(fwf_t(:)); % integrated
                tmp=nansum(fwf_t);
                flux.fw_x{n}(idt,:) = tmp(1:2:end-1)+tmp(2:2:end); % 1d
                flux.fw_d{n}(idt,:,:)=fwf_t(:,1:2:end-1)+fwf_t(:,2:2:end); % 2d
            else
                flux.fw(n,idt)= sum(fwf_t(:));
            end
        end
        if ~isempty(strfind(strin,'t'));
            if ~isempty(strfind(strin,'x'));
                flux.tem(n,idt)= nansum(tem_t(:)); % integrated 
                tmp=nansum(tem_t);
                flux.tem_x{n}(idt,:) = tmp(1:2:end-1)+tmp(2:2:end); % 1d
                flux.tem_d{n}(idt,:,:)=tem_t(:,1:2:end-1)+tem_t(:,2:2:end); % 2d
            else
                flux.tem(n,idt)= sum(tem_t(:));
            end
        end       
        if ~isempty(strfind(strin,'d'));
            if ~isempty(strfind(strin,'x'));
                for m=1:nsed
                    flux.sed(n,idt,m)= nansum(sed_t{m}(:));
                    tmp=nansum(sed_t{m});
                    flux.sed_x{n}(idt,:,m) = tmp(1:2:end-1)+tmp(2:2:end);
                    flux.sed_d{n}(idt,:,:,m) = sed_t{m}(:,1:2:end-1)+sed_t{m}(:,2:2:end);
                end
            else
                for m=1:nsed
                    flux.sed(n,idt,m)= sum(sed_t{m}(:));
                end
            end
        end       
    flux.area_integr{n}=Area(n,idt); % integrated area
    end  % loop over sections
end
close(nc);
return
