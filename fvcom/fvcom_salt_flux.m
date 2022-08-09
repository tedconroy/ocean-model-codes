% calculate salt fluxes through cross sections, calculate TEF
% ted conroy

% file directory
fdir=('/projects/oceans/tconroy/FVCOM/FVCOM3.2.1/run/w14/out_w_mixing/out/');  addpath(fdir);
savename='w14_w_mixing';
fnames=dir([fdir, '/*.nc']);

% what is time output? (for filtering)
dt=4; % e.g.:  1=hourly, 4=15minute, 3=20minute

% remove restart and avg files
rem=0; for n=1:length(fnames); if strfind(fnames(n).name,'restart') 
rem=[rem; n]; end ; end; rem=rem(2:end);

rem=0; for n=1:length(fnames); if strfind(fnames(n).name,'avg')
rem=[rem; n]; end ; end; rem=rem(2:end);

% set up stuff
if length(rem)>0;fnames=fnames(1:rem(1)-1);end; clear n;
cd /projects/oceans/tconroy/FVCOM/FVCOM3.2.1/Post-processing/matlab_functions/nctoolbox-1.1.3/; setup_nctoolbox;
addpath(genpath('/projects/oceans/tconroy/FVCOM/FVCOM3.2.1/Post-processing/matlab_functions'))
cd /projects/oceans/tconroy/FVCOM/FVCOM3.2.1/Post-processing
load('/projects/oceans/tconroy/FVCOM/FVCOM3.2.1/Post-processing/input_matfiles/mobj_&_transects.mat','Mobj')
load('/projects/oceans/tconroy/FVCOM/FVCOM3.2.1/Post-processing/input_matfiles/cross_sections_4fluxes.mat')
[lonc latc]=sp_proj('3602','inverse',Mobj.xc,Mobj.yc,'m');
tri = delaunayTriangulation(lonc,latc);
K = nearestNeighbor(tri,[Mobj.lon(:) Mobj.lat(:)]);
h=ncread([fdir, fnames(1).name],'h');
clear cross_section lxsect; lxsect=xsect;

%% loop through cross sections
for n=1:length(lxsect) 
    for i=1:length(fnames) % get fluxes at each output
        fname=([fdir, fnames(i).name]);
        tcur=datenum(nc_varget([fdir, fnames(i).name],'time'))+datenum('1858-11-17 00:00:00'); 
        if n==1; if i==1; time=tcur; end; time=[time;tcur]; end
        [cross_section.fluxes(i).f]=fvcom_calcfluxsect(fname,[lxsect(n).x0 lxsect(n).y0 lxsect(n).x1 lxsect(n).y1],[tcur tcur],'v s x f',1);
    end
    cross_section.vol=0; % make timeseries
    
    for i=1:length(cross_section.fluxes)
        % integrals over cross section
        cross_section.vol=[cross_section.vol; cross_section.fluxes(i).f.vol ];
        area_x(i,:,:)= sqz(cross_section.fluxes(i).f.area{:});
        % total area of cross section
        cross_section.area(i)=squeeze(nansum(nansum(cross_section.fluxes(i).f.area{:}))); % only wet area
        cross_section.area_integr_w_dry(i)=cell2mat(cross_section.fluxes(i).f.area_integr);  % includes dry+wet area
    end
        tef(n).vol=cross_section.vol(2:end);     
        tef(n).fw(i) = cross_section.fluxes(i).f.fw;
    
    %% fluxes, total exchange flow
    nbins=1000; % salt bins
    range=[0 36]; % salinity range
    bins=linspace(range(1),range(2),nbins);
    binsize=bins(end)-bins(end-1);
    
    % get output
    for f=1:length(cross_section.fluxes) % each output time
        salt=cross_section.fluxes(f).f.sal_d{:}./cross_section.fluxes(f).f.vol_d{:};
        salt_flux_timeseries(f,:,:)=sqz(cross_section.fluxes(f).f.sal_d{:});
        salt_timeseries(f,:,:)=sqz(salt);
        vol_timeseries(f,:,:)=sqz(cross_section.fluxes(f).f.vol_d{:});
        vel_timeseries(f,:,:)=sqz(cross_section.fluxes(f).f.vol_d{:})./sqz(cross_section.fluxes(f).f.area{:})';
        area_timeseries(f,:,:)=sqz(cross_section.fluxes(f).f.area{:})';
        
        % each salt bin
        for i=2:length(bins) 
            j=find( salt >= bins(i)-binsize & salt < bins(i));
            cross_section.q(f).salt(i).s=salt(j); % bin salt in salt class
            cross_section.q(f).vol(i).v=cross_section.fluxes(f).f.vol_d{:}(j);  % bin volume flux in current salt class
            cross_section.q(f).vol_sum(i).v=sum(cross_section.fluxes(f).f.vol_d{:}(j));
            tef(n).dq(i,f)=cross_section.q(f).vol_sum(i).v;
        end
        
    end
    
    % Tidally filter volume fluxes in each salt bin
    for s=1:length(bins)
        tofilt=tef(n).dq(s,:);
        tef(n).filt_dq(s,:)=godinfilt(tofilt(:),0,dt);
    end
   
    % Get Eulerian component (t,z,y), and 'eulerian' decomposition components
    for z=1:length(salt_timeseries(1,:,1)) % depth
        for y=1:length(salt_timeseries(1,1,:)) % cross estuary
            salt_filt(:,z,y)=godinfilt(salt_timeseries(:,z,y),0,dt); % s1
            vel_filt(:,z,y)=godinfilt(vel_timeseries(:,z,y),0,dt); % u1
            vol_filt(:,z,y)=godinfilt(vol_timeseries(:,z,y),0,dt); % u1
            area_filt(:,z,y)=godinfilt(area_timeseries(:,z,y),0,dt); % filtered dA 
            udA_f(:,z,y)=godinfilt(sqz(vel_timeseries(:,z,y) .* area_timeseries(:,z,y)),0,dt);
            sdA_f(:,z,y)=godinfilt(sqz(salt_timeseries(:,z,y) .* area_timeseries(:,z,y)),0,dt);     
        end
    end    
    
    % salt flux terms
    tef(n).total_salt_flux = sqz(nansum(nansum(salt_flux_timeseries,3),2)); %
    tef(n).a0=godinfilt(sqz(nansum(nansum(area_timeseries,3),2)),0,dt);% a0
    tef(n).u0 = godinfilt(sqz(nansum(nansum(vol_timeseries,3),2)),0,dt)./tef(n).a0;
    tef(n).s0 = godinfilt(sqz(nansum(nansum(salt_timeseries .* area_timeseries,3),2)),0,dt) ./ godinfilt(sqz(nansum(nansum(area_timeseries,3),2)),0,dt);
    u1 = udA_f ./ area_filt - tef(n).u0;
    s1 = sdA_f ./ area_filt - tef(n).s0;
    u2 = vel_timeseries - tef(n).u0 - u1;
    s2 = salt_timeseries - tef(n).s0 - s1;
    tef(n).Fr = tef(n).u0 .* tef(n).s0 .* tef(n).a0;
    tef(n).Fe = sqz(nansum(nansum(u1.*s1.*area_filt,3),2));
    tef(n).Ft = godinfilt(sqz(nansum(nansum(u2.*s2.*area_timeseries,3),2)),0,dt);
    
    % Integrate
    for f=1:length(cross_section.fluxes) % each output time
            s_filt=sqz(salt_filt(f,:,:));
            ve_filt=sqz(vel_filt(f,:,:)); 
            vo_filt_tmp=sqz(vol_filt(f,:,:));
            io_ind=find(vo_filt_tmp>0); io_ind2=find(vo_filt_tmp<0);
            vol_filtc=vo_filt_tmp(:); salt_filtc=s_filt(:);
            tef(n).q_class_eul_in(f)=sqz(nansum(vol_filtc(io_ind)));
            tef(n).q_class_eul_out(f)=sqz(nansum(vol_filtc(io_ind2)));
            tef(n).f_class_eul_in(f)=sqz(nansum(vol_filtc(io_ind).*salt_filtc(io_ind)));
            tef(n).f_class_eul_out(f)=sqz(nansum(vol_filtc(io_ind2).*salt_filtc(io_ind2)));
            a_filt=sqz(area_filt(f,:,:));
        for i=1:length(bins) % salt bins
            jj=find( s_filt >= bins(i)-binsize & s_filt < bins(i)); % find data to bin
            if isempty(jj);
                eul.vol(i,f)=0;
            else
                eul.vol(i,f)=nansum(a_filt(jj).*ve_filt(jj));
            end
        end
    end
    % Eulerian DqDs
    tef(n).eul_dqds=eul.vol; 
    
    % total exchange flow
    for f=1:length(cross_section.fluxes) % nfiles
        % New method, use Q(s) to find integration points
        tef(n).filt_Q(:,f)=cumsum(tef(n).filt_dq(:,f),'reverse'); % Q(s)
        [loc,pks]=peakseek(tef(n).filt_Q(:,f),1000000000,1); % peakseek 10x faster than findpeaks!
        [mloc,mpks]=peakseek(-1.*tef(n).filt_Q(:,f),1000000000,1);  mpks=mpks*-1; 
        % find number of inflow/outflow layers
        si=sign(tef(n).filt_Q(:,f)); is0si=find(si==0); if is0si(end)-is0si(end-1)==1; start=is0si(end); else; disp('something went wrong: change code'); end
        dsi=diff(si); nlay=length(find(abs(dsi) > 0)); % delta function at each sign change
        if loc > mloc; normal=1;else; normal=0; end % location of maxima vs minima to determine in/out going from highest salinity
        tef(n).nlayers(f)=nlay;
        
        % do binning, here assume 2 layer normal (inflow higher salinity)
        sp=loc:length(bins);
        sm=1:loc-1;
        tef(n).q_in(f)=sum(tef(n).filt_dq(sp,f));
        tef(n).q_out(f)=sum(tef(n).filt_dq(sm,f));
        tef(n).f_in(f)=sum(tef(n).filt_dq(sp,f).*(bins(sp))');
        tef(n).f_out(f)=sum(tef(n).filt_dq(sm,f).*(bins(sm))');
        tef(n).s_in(f)=tef(n).f_in(f)/tef(n).q_in(f);
        tef(n).s_out(f)=tef(n).f_out(f)/tef(n).q_out(f);
        
        % eulerian component
        tef(n).filt_eul_Q(:,f)=cumsum(tef(n).eul_dqds(:,f),'reverse'); % Eul Q(s)
        [loc_eul,pks]=peakseek(tef(n).filt_eul_Q(:,f),1000000000,1);
        sfp=loc_eul:length(bins);
        sfm=1:loc_eul-1;
        tef(n).q_eul_in(f)=sum(eul.vol(sfp,f));
        tef(n).q_eul_out(f)=sum(eul.vol(sfm,f));
        tef(n).f_eul_in(f)=sum(tef(n).eul_dqds(sfp,f).*(bins(sfp))');
        tef(n).f_eul_out(f)=sum(tef(n).eul_dqds(sfm,f).*(bins(sfm))');
        tef(n).s_eul_in(f)=tef(n).f_eul_in(f)/tef(n).q_eul_in(f);
        tef(n).s_eul_out(f)=tef(n).f_eul_out(f)/tef(n).q_eul_out(f);
        
        % old method, use dQdS, integrate +/-
        clear sm sp sfm sfp nlay loce mloce loc mloc
        sm=find(tef(n).filt_dq(:,f) < 0 );
        sp=find(tef(n).filt_dq(:,f) > 0 );
        sfm=find(eul.vol(:,f) < 0 );
        sfp=find(eul.vol(:,f) > 0 );
        tef_q(n).q_eul_in(f)=sum(eul.vol(sfp,f));
        tef_q(n).q_eul_out(f)=sum(eul.vol(sfm,f));
        tef_q(n).q_in(f)=sum(tef(n).filt_dq(sp,f));
        tef_q(n).q_out(f)=sum(tef(n).filt_dq(sm,f));
        tef_q(n).f_in(f)=sum(tef(n).filt_dq(sp,f).*(bins(sp))');
        tef_q(n).f_out(f)=sum(tef(n).filt_dq(sm,f).*(bins(sm))');
        tef_q(n).s_in(f)=tef_q(n).f_in(f)/tef_q(n).q_in(f);
        tef_q(n).s_out(f)=tef_q(n).f_out(f)/tef_q(n).q_out(f);
    end

    clear salt_timeseries vol_timeseries area_timeseries salt vel_timeseries eul area_x U S j el loc U VOL S salt_filt vel_filt area_filt vol_filt
    clear cross_section st ut u2 u1 udA_f sdA_f s1 s2 salt_flux_timeseries
    disp(['cross section' num2str(n) ' out of ' num2str(length(lxsect)) ' done']);
    % save every cross section
    %save(['./matfile_output/saltflux_x_' num2str(n) '_' savename '.mat']','-v7.3')
            
end % loop for each cross section
            
save(['./matfile_output/saltflux' savename '.mat']','-v7.3')
