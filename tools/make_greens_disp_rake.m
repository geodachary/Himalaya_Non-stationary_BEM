function [Grake_east,Grake_north,Grake_up] = make_greens_disp_rake(el,nd,rake,xysites)


%INPUTS:
%   el,nd = Nx3 triangular element parameterization for mesh with N elements
%   xysites = Mx2 vector of M surface position coords in xy coordinate system
%
%OUTPUTS:
%   Disp.ss = 3MxN matrix of surface displacements due to strike-slip on triangular elements
%   Disp.ds = 3MxN matrix of surface displacements due to dip-slip on triangular elements
%

nu = 0.25;
staxy = [xysites'; zeros(1, size(xysites,1))];

ss = cos(rake*pi/180);
ds = sin(rake*pi/180);
    
for j=1:size(el,1)


temp1{1}=nd;
temp2{1}=el(j,:);    

[U, D, S] = tridisloc3d(staxy, temp1{1}', temp2{1}', [ss(j) ds(j) 0]', 1, .25);
%[Uds2, D, S] = tridisloc3d(staxy, temp1{1}', temp2{1}', [ss(j) 0 0]', 1, .25);  %negative sign because of chosen sense of slip convention 

    
 
Grake_east(:,j) = U(1,:)';
Grake_north(:,j) = U(2,:)';
Grake_up(:,j) = U(3,:)';


       
if mod(j,50)==0; 
        disp(['Displacements: Finished node number ' num2str(j)])
end


end %j

%Disp.ds = Gds;
%Disp.ss = Gss;
