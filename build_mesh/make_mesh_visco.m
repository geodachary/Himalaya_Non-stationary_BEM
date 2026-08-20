%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script builes a trianglur mesh of a surface defined
% by depth contours based on code of Yo Fukushima, 9 Jul 2010

%origin to convert from llh to local Caresian coordinates
origin = [origin_lon, origin_lat];  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath ../tools/
addpath ../ne_10m_coastline/


% rakes_rates = load(conv_rate_file);
rakes_rates = rates_interp;


figure;scatter(contrkm(:,1),contrkm(:,2),10,contrkm(:,3),'filled'); colorbar; 
xlabel('Easting'); ylabel('Northing'); daspect([1,1,1]);
hold on
bbox = [73.04 26.6; 95.18 38]; plot_coast_xy(bbox,origin,'k')



figure;
h = trisurf(el,nd(:,1),nd(:,2),nd(:,3)); colorbar; view(2); daspect([1,1,1]);
xlabel('Easting'); ylabel('Northing'); daspect([1,1,1]);
hold on
bbox = [73.04 26.6; 95.18 38]; plot_coast_xy(bbox,origin,'k')

% %generate triangular mesh geometry parameters
patch_stuff=make_triangular_patch_stuff(el,nd);

rakes = rakes_interp;
rates = rates_interp;


% plot the rakes
figure;
h = trisurf(el,nd(:,1),nd(:,2),nd(:,3),rakes); colorbar; view(2); daspect([1,1,1]);
xlabel('Easting'); ylabel('Northing'); daspect([1,1,1]);
hold on
bbox = [73.04 26.6; 95.18 38]; plot_coast_xy(bbox,origin,'k')
title('convergence rake')

% plot the rates
figure;
h = trisurf(el,nd(:,1),nd(:,2),nd(:,3),rates); colorbar; view(2); daspect([1,1,1]);
xlabel('Easting'); ylabel('Northing'); daspect([1,1,1]);
hold on
bbox = [73.04 26.6; 95.18 38]; plot_coast_xy(bbox,origin,'k')
title('plate convergence rate')

% plot the vectors of strike and dips
figure;
h = trimesh(el,nd(:,1),nd(:,2),nd(:,3)); view(2); daspect([1,1,1]);
hold on
quiver3(centroids(:,1),centroids(:,2),centroids(:,3),patch_stuff.strikevec_faces(:,1),patch_stuff.strikevec_faces(:,2),patch_stuff.strikevec_faces(:,3))
quiver3(centroids(:,1),centroids(:,2),centroids(:,3),patch_stuff.dipvec_faces(:,1),patch_stuff.dipvec_faces(:,2),patch_stuff.dipvec_faces(:,3))
title('patch strike and dip vectors')
