function nn = shortestPath(xx,yy,nextnc,n1,n2,sf)

% shortest path b/t 2 cells or nodes
%   nn = shortestPath(xx,yy,nextnc,n1,n2);
% for nodes or cells:
%  nodes --> xx=xn,yy=yn,nextnc=nbsn;
%  cells --> xx=xc,yy=yc,nextnc=nbe;
%    can choose 2 grid points (n1,n2)
%    or a bunch of grid points in n1

% add a safety factor so get an answer w/ highly regular grids
% (simple_mesh_maker or burchest)
if nargin<6
    sf=1;
end

% reshape nextnc if need be
if size(nextnc,2)>size(nextnc,1);
    nextnc=nextnc';
end

if length(n1)>1
    n2=n1(2:end);
    n1=n1(1:end-1);
end
nna=[];
for m=1:length(n1)
    % draw a line b/t two points, and follow path that is closest to the line
    L = createLine([xx(n1(m)) yy(n1(m))],[xx(n2(m)) yy(n2(m))]); 

    % move from pt 1 to pt 2, staying near line
    nn = [n1(m)]; n=1;
    dnL = distancePointLine([xx(:) yy(:)],L);
    while nn(n) ~= n2(m)
        n=n+1;
        % adjacent cells
        ii = nextnc(nn(n-1),:);
        ii = ii(~isnan(ii) & ii > 0); % boundary or non-existent cells/nodes      
        % not already in path
        ii = setdiff(ii,nn);    
        % closer to end of path than previously
        dtmpn =sqrt((xx(nn(n-1)) - xx(n2(m))).^2 +(yy(nn(n-1)) - yy(n2(m))).^2);
        dtmp =[];
        for i=1:length(ii);
            dtmp(i) = sqrt((xx(ii(i)) - xx(n2(m))).^2 +(yy(ii(i)) - yy(n2(m))).^2);
        end
        jj=dtmp < sf*dtmpn*ones(size(ii));  % safety factor for regular grids
        if sum(jj)<1
%             disp('warning: no shortestPath');
%             nn=[];
%             break
            disp('warning: had to make path longer');
            [jnk,jj]=min(dtmp);
        end
        ii = ii(jj);    

        % closest to straight line path
        [a,i] = min(dnL(ii));
        nn = [nn ii(i)];
    end
    nna=[nna nn];
end

nn=unique_no_sort(nna);