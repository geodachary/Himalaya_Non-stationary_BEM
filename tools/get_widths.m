function [U_xyz, width] = get_widths(ypts,depths_U,depths_L, centroids, nd)

mesh_xyz = [nd;centroids];

%This function computes the width of the locked zone at node positions
%(ypts). For each (ypts,depths_U) node, find the nearest point on the lower
%boundary and compute the distance. 


pts_U = [ypts depths_U];
pts_L = [ypts depths_L];
   
 
%spline fit to point, cscvn is spline interpolation
%c = fnplt(cscvn(pts_U'),'r',2);
c = fnplt(pchip(pts_U(:,1)',pts_U(:,2)'),'r',2);

for k=1:3
    %add more points
    cnew = zeros(2,2*size(c,2)-1);
    cnew(:,1:2:end) = c;
    cnew(1,2:2:end) = (c(1,1:end-1)+c(1,2:end))/2;
    cnew(2,2:2:end) = (c(2,1:end-1)+c(2,2:end))/2;
    c = cnew;
end

%toss out repeats
[y,ind] = unique(c(1,:));
c = c(:,ind);
c_U = c';



%spline fit to points
%c = fnplt(cscvn(pts_L'),'r',2);
c = fnplt(pchip(pts_L(:,1)',pts_L(:,2)'),'r',2);

for k=1:3
    %add more points
    cnew = zeros(2,2*size(c,2)-1);
    cnew(:,1:2:end) = c;
    cnew(1,2:2:end) = (c(1,1:end-1)+c(1,2:end))/2;
    cnew(2,2:2:end) = (c(2,1:end-1)+c(2,2:end))/2;
    c = cnew;
end

%toss out repeats
[y,ind] = unique(c(1,:));
c = c(:,ind);
c_L = c';



%% find x-coordinate of L boundary points
c_L_x = griddata(mesh_xyz(:,2),-mesh_xyz(:,3),mesh_xyz(:,1),c_L(:,1),c_L(:,2));
c_L_x(isnan(c_L_x)) = griddata(mesh_xyz(:,2),-mesh_xyz(:,3),mesh_xyz(:,1),c_L(isnan(c_L_x),1),c_L(isnan(c_L_x),2),'nearest');


%% find x-coordinate of U boundary points
c_U_x = griddata(mesh_xyz(:,2),-mesh_xyz(:,3),mesh_xyz(:,1),c_U(:,1),c_U(:,2));
c_U_x(isnan(c_U_x)) = griddata(mesh_xyz(:,2),-mesh_xyz(:,3),mesh_xyz(:,1),c_U(isnan(c_U_x),1),c_U(isnan(c_U_x),2),'nearest');

%% find width and distances of zones
%width of transition zone
width = sqrt( (c_U(:,1)-c_L(:,1)).^2 + (c_U_x-c_L_x).^2 + (c_U(:,2)-c_L(:,2)).^2);



%(x,y,z) poistion of upper locking depth points
U_xyz = [c_U_x c_U];