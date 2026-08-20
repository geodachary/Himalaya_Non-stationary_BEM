function [Gexx,Gexy,Geyy] = make_strainG_triangular(xy_obs,tri,nd,rake)


% Compute surface strain rates due to slip on triangular disloction with
% given rake
% xy_obs = Mx2 matrix of observation coordinates (x,y) (km)
% tri = matrix of triangle indices
% nd = matrix of triangle nodes (x,y coordinates) (km)
% rake = vector of rake angles for each slip patch


%get patch geometry
patch_stuff = make_triangular_patch_stuff(tri,nd);


centroids = patch_stuff.centroids_faces;
NormalFaces = patch_stuff.normal_faces;
StrikeVecFaces = patch_stuff.strikevec_faces;


 
ss = cos(rake*pi/180);
ds = sin(rake*pi/180);

%convert from km to m
nd = 1000*nd;
xy_obs = 1000*xy_obs;



for j=1:size(tri,1)

    
    t=tri(j,:);    

    %call triangular dislocation code
    [d,Strain] = CalcTriStrains_O(xy_obs(:,1), xy_obs(:,2), xy_obs(:,3), nd(t,1), nd(t,2), nd(t,3), 0.25, ss(j), 0, -ds(j));    
   
 
   %save strain components in matrices
    Gexx(:,j) = Strain.xx;  
    Gexy(:,j) = Strain.xy;  
    Geyy(:,j) = Strain.yy;  
   


end %j



