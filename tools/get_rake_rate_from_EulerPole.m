

%% compute strike and dip vectors for each triangular patch


patch_stuff_llh = make_triangular_patch_stuff(el,nd_ll); %this function returns a structure with triangle properties

%rename some geometric variables
centers_llh_small = patch_stuff_llh.centroids_faces;  %centroids of triangles
strikevecs_small = patch_stuff_llh.strikevec_faces; %strike vectorss
dipvecs_small = patch_stuff_llh.dipvec_faces;  %dip vectors
dip_angle = patch_stuff_llh.dip_faces;  %dip angle (degrees)

%plate_vel computes horizontal velocities (Veast, Vnorth) for given Euler
%pole (latp, lonp, omega)
vel_plate = plate_vel(centers_llh_small(:,2), centers_llh_small(:,1), latp, lonp, omega)';
%vel_plate = - vel_plate;  %switch sign so that it is motion of over-riding plate relative to subducting
vel_plate = fliplr(vel_plate(:,1:2));  %flip to (ve,vn)


%compute rake 
norm_vel_plate = sqrt(vel_plate(:,1).^2+vel_plate(:,2).^2);
rake = acos(dot(strikevecs_small(:,1:2),vel_plate,2)./norm_vel_plate)*180/pi;
rake_sign = -sign(dot(dipvecs_small(:,1:2),vel_plate,2));
rake_patches = rake.*rake_sign;

%compute slip rate on patches
rate_patches = norm_vel_plate./cos(dip_angle*pi/180);


