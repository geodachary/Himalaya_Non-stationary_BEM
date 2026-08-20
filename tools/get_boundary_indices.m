function i_lock= get_boundary_indices(ypts,depths_U,depths_L,centroids)



pts_U = [ypts depths_U];
pts_L = [ypts depths_L];
   
 
%spline fit to points
c = fnplt(cscvn(pts_U'),'r',2);
%c = fnplt(pchip(pts_U(:,1)',pts_U(:,2)'),'r',2);

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
%c = fnplt(pchip(pts_L(:,1)',pts_L(:,2)'),'r',2);

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



%map centroids to nearst boundary locations
N = length(c_L(1,:));
cy_L_ind = interp1(c_L(1,:),1:N,centroids(:,2),'nearest');
N = length(c_U(1,:));
cy_U_ind = interp1(c_U(1,:),1:N,centroids(:,2),'nearest');

%find centroids between boundaries
i_lock = -centroids(:,3) < c_L(2,cy_L_ind)' & -centroids(:,3) > c_U(2,cy_U_ind)';




