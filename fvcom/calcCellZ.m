function zz = calcCellZ(nk,hin,etain)
%    zz=calcCellZ(nk,hin,etain)
% calc fvcom z levels from # sigma levels, h, and eta
%
% can give this fxn actual sigma layers in nk
% e.g. zz=calcCellZ(21,h,zeta);


    if length(nk)==1;
        dsi = 1/nk ;
        si = (-dsi/2:-dsi:-1+dsi/2)' ;           
    else
        si=nk;
    end
    
    etain=etain(:)';   
    hin=hin(:)';   
    if length(etain)==1; 
        etain=etain*ones(size(hin));
    end
    htmp=etain+hin; % water depth
    
    zz = si*htmp + ones(size(si))*etain; % cell elevation (m)
