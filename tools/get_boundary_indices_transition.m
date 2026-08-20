function [i_lock,i_transition] = get_boundary_indices_transition(ypts,depths_U,depths_L,depths_T,centroids)



pts_U = [ypts depths_U];
pts_L = [ypts depths_L];
pts_T = [ypts depths_T];
   
 
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



%spline fit to points
c = fnplt(cscvn(pts_T'),'r',2);

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
c_T = c;


%map centroids to nearst boundary locations
N = length(c_L(1,:));
cy_L_ind = interp1(c_L(1,:),1:N,centroids(:,2),'nearest');
N = length(c_U(1,:));
cy_U_ind = interp1(c_U(1,:),1:N,centroids(:,2),'nearest');
N = length(c_T(1,:));
cy_T_ind = interp1(c_T(1,:),1:N,centroids(:,2),'nearest');

%find centroids between boundaries
i_lock = -centroids(:,3) < c_L(2,cy_L_ind)' & -centroids(:,3) > c_U(2,cy_U_ind)';

i_transition = -centroids(:,3) < c_T(2,cy_L_ind)' & -centroids(:,3) > c_U(2,cy_U_ind)';



