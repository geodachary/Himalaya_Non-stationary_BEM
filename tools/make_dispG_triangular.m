function [Geast,Gnorth,Gup] = make_dispG_triangular(tri,nd,staxy,rake,rel)


%initialize to zero
Geast = zeros(size(staxy,2),size(tri,1));
Gnorth = zeros(size(staxy,2),size(tri,1));
Gup = zeros(size(staxy,2),size(tri,1));


%strike-slip and dip-slip components
ss = cos(rake*pi/180);
ds = sin(rake*pi/180);


for j=1:size(tri,1)



    %Compute displacements for triangular sources
    [U, D, S] = tridisloc3d(staxy, nd', tri(j,:)', [-ss(j) -ds(j) 0]', 1, .25);  %negative indicates normal slip and right-lateral -- this is correct for backslip and the definition of rake
    

    %rel is a vector of station indices for computing relative
    %displacements (can be empty)
    if ~isempty(rel)
        Geast(:,j) =  (U(1,:) - mean(U(1,rel)))';
        Gnorth(:,j) = (U(2,:) - mean(U(2,rel)))';
        Gup(:,j) = (U(3,:) - mean(U(3,rel)))';
    else
        Geast(:,j) =  U(1,:)';
        Gnorth(:,j) = U(2,:)';
        Gup(:,j) = U(3,:)';
    end

    disp(['Completed ' num2str(j/size(tri,1)*100) '% of displacement Greens function caclulations.'] )


end %j



