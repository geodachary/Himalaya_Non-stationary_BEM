function [el,nd] = make_tri_mesh(contrkm,intv)



%% identify the points of the shallowest and deepest contours
% nodet: the shallowest points
% nodeb: the deepest points
ind = (contrkm(:,3) == contrkm(1,3));% hardwired since first element is shallowest
nodet0 = contrkm(ind,:);
ind = (contrkm(:,3) == contrkm(end,3));
nodeb0 = contrkm(ind,:);

% here I remove edge points to make sure that griddata works properly
nodet0(end,:) = [];
nodet0(1,:) = [];
nodeb0(end,:) = [];
nodeb0(1,:) = [];


%% set the intervals of shallowest and deepest points as specified 
N = round(curvlength(nodet0)./intv);
nodet = curvspace(nodet0,N);
N = round(curvlength(nodeb0)./intv);
nodeb = curvspace(nodeb0,N);


%% mesh generation
[nd0,el] = meshfrac2(nodet,nodeb,intv);
% nd0  are coordinates of nodes on the mesh
% el   are the elements, the values specify the x,y,z indices of the nodes

%% shift the height of the nodes based on the contour info
% the mesh previously made is a straight-line interpolation from the top
% nodes to bottom nodes
[xi,yi,zi] = griddata(contrkm(:,1),contrkm(:,2),contrkm(:,3),...
    nd0(:,1),nd0(:,2),'linear');
ind = find(nd0(:,3) ~= nd0(1,3) & nd0(:,3) ~= nd0(end,3));
nd = nd0;
nd(ind,:) = [xi(ind),yi(ind),zi(ind)];

