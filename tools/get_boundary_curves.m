function [c_U,c_L]= get_boundary_curves(ypts,depths_U,depths_L,centroids,centroids_rot)



pts_U = [ypts depths_U];
pts_L = [ypts depths_L];
   
 
%spline fit to points
c = fnplt(cscvn(pts_U'),'r',2);

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
c_U = c;



%spline fit to points
c = fnplt(cscvn(pts_L'),'r',2);

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
c_L = c;


%interpolate to get x-positions of boundary nodes
cx_L = griddata(-centroids_rot(:,3),centroids_rot(:,2)-min(centroids_rot(:,2)),centroids_rot(:,1),c_L(2,:),c_L(1,:));
cx_U = griddata(-centroids_rot(:,3),centroids_rot(:,2)-min(centroids_rot(:,2)),centroids_rot(:,1),c_U(2,:),c_U(1,:));

%coordinates in rotated coordinates
c_L = [cx_L' c_L(1,:)' c_L(2,:)'];
c_U = [cx_U' c_U(1,:)' c_U(2,:)'];





