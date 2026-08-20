%clear all
%close all

% load setup.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%  BEGIN SETUP   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%number of boundary nodes along strike 
num_Dpts = 20;

%upper locked boundary depths (km)
locked_depths_U = 0*ones(num_Dpts,1);

%lower locked boundary depths (km)
%locked_depths_L = 10*ones(num_Dpts,1);
locked_depths_L = [15*ones(1,8) 15*ones(1,4) 15*ones(1,4) 15*ones(1,4)]';

%accumulated stress on ring around asperities (Pa)
ring_tau = 0.15e6*ones(num_Dpts,1); %multiply it with zero if stationary model is preferred


    
%load inversion results?
load_inversion = false;

%folder name for inversion results (if load_inversion is true)
folder_name = 'Block_I_outputs_Bline';


D = 3;  %process-zone width (km)


%time step for propagation
dT = 100;


%if inverting velocities, provide long-term velocity field in a file with
%format:
%lon lat Ve(mm/yr) Vn(mm/yr) (in same order as data file)  -- Here Ve and Vn are
%long-term velocities (not GPS-derived)
longterm_velocity_file = 'Himalaya_block_velocities.txt';


%elastic modulus
mu = 3e10;


%specify start and end points of slip rate profiles across interface
%each row is a profile:  start x, start y, end x, end y
profiles = [-775 550 -630 685; -520 250 -350 450; -235 75 -60 265; 400 -120 475 125; 810 -122 890 130; 1070 -5 1060 225];

%%%%  END SETUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if load_inversion

    eval(['load ./' folder_name '/M_locked_depths_U.txt'])
    locked_depths_U =  M_locked_depths_U(end,:)';
    
    eval(['load ./' folder_name '/M_locked_depths_L.txt'])
    locked_depths_L =  M_locked_depths_L(end,:)';
    
    eval(['load ./' folder_name '/M_ring_tau.txt'])
    ring_tau =  M_ring_tau(end,:)';

end

%% load long-term velocities
if invert_vel
    lt_vels = load(longterm_velocity_file);
    lt_vels = lt_vels(:,3:4);  %columns are Ve, Vn
end


% ---- Software configuration
% Make hmmvp available
addpath hmmvp0.16/
addpath tools/
addpath hmat/

[gamb.hmat.id nnz] = hm_mvp('init',gamb.hmat.savefn,4);


% Generate triangular patch data and extract face centroids and strike angles
patch_stuff = make_triangular_patch_stuff(el,nd);
centroids = patch_stuff.centroids_faces;
strikes = patch_stuff.strike_faces;

%rotate and shift interface so that average strike is aligned with y axis
%this is needed for computing boundary positions
strikes = mod(strikes, 360);    % strikes in degrees (any range)

theta = deg2rad(strikes);
str = atan2(mean(sin(theta)), mean(cos(theta)));   % radians, robust mean angle

R = [cos(str) -sin(str); sin(str) cos(str)];
centroids_rot = (R*centroids(:,1:2)')';
centroids_rot(:,3) = centroids(:,3);
centroids_rot(:,2) = centroids_rot(:,2) - min(centroids_rot(:,2));


%node points along strike
Dpts = linspace(0,max(centroids_rot(:,2)),num_Dpts)';

i_locked = get_boundary_indices(Dpts,locked_depths_U,locked_depths_L,centroids_rot);
i_ring = get_boundary_indices(Dpts,locked_depths_U-D,locked_depths_L+D,centroids_rot);

%include ring of accumulated stress
ring_stress = interp1(Dpts,ring_tau,centroids_rot(:,2));
ring_tau = zeros(size(el,1),1);
ring_tau(i_ring) = ring_stress(i_ring);

%slip rate in mm/yr to m/yr
srate = rates/1000;
% %adjust for stiffness
nel = size(el,1);
% loading_rate = -hm_mvp('mvp',gamb.hmat.id,srate);
% [max_depth,i] = max(abs(centroids(i_ring)));
% lr = loading_rate(i_ring);
% stiffness = loading_rate/lr(i);

stiffness = self_rate;


% Change the Convergence Criteria
num = 1:size(el,1);
maxit = 200;
tol = 1e-4;

%need to scale stresses
%GFs computed with km, convert to meters
%GFs computed for mu=1
scale = mu*10^-3;
 

asp_index = i_locked;

bslip = zeros(size(el,1),1);
bslip(asp_index) = srate(asp_index);  %mm/yr, backslip rate
dtau_rate = scale*hm_mvp('mvp',gamb.hmat.id,bslip);
dtau = dtau_rate*dT; %accumulated stress (stressing rate times T1)

% include ring of accumulated stress with stiffness
% ring_tau = zeros(size(el,1),1);
% ring_tau(i_ring) = ring_stress(i_ring);
ring_tau = ring_tau.*stiffness(:);


%include ring stress
dtau = -dtau + ring_tau;



%increment of displacement due to far field load and co+after stress
rs = num(~asp_index);
s = gmres(@(x)mvp(x,gamb.hmat.id,rs,rs,~asp_index),-dtau(~asp_index)/scale,[],tol,maxit);
inc_slip = zeros(size(el,1),1);
inc_slip(~asp_index)=s; %increment of forward slip
U = srate*dT + inc_slip;  %srate*T1 to convert to displacement over time T1
U(asp_index)=0;  %locked asperity
creep_rate = U/dT;



%% Plot creep rate

figure
h=trisurf(el,nd(:,1),nd(:,2),nd(:,3),creep_rate,'edgecolor','none'); colorbar;  daspect([1,1,1])
colormap(jet)
colorbar
axis equal
axis tight      
% clim([0 1.25*max(srate)])
set(gca,'fontsize',15)
title('Creep Rate (m/yr)')
view(2)
hold on
bbox = [73.04 26.6; 95.18 38]; plot_coast_xy(bbox,origin,'k')
ylim([-600 600])
plot([profiles(:,1) profiles(:,3)]', [profiles(:,2) profiles(:,4)]','linewidth',2)

%% Plot Slip Deficit Rate
figure
h=trisurf(el,nd(:,1),nd(:,2),nd(:,3),srate-creep_rate,'edgecolor','none'); colorbar;  daspect([1,1,1])
colormap(jet)
colorbar
axis equal
axis tight      
set(gca,'fontsize',15)
title('Slip Deficit Rate (m/yr)')
view(2)
hold on
bbox = [73.04 26.6; 95.18 38]; plot_coast_xy(bbox,origin,'k')
ylim([-200 1000])
   
ring_stress = ring_tau;
ring_stress(asp_index) = 0;
ring_stress(~i_ring) = 0;

% Plot Ring Stress
figure
h=trisurf(el,nd(:,1),nd(:,2),nd(:,3),ring_stress,'edgecolor','none'); colorbar;  daspect([1,1,1])
colormap(jet)
colorbar
axis equal
axis tight      
set(gca,'fontsize',15)
title('Ring Stress')
view(2)
hold on
bbox = [73.04 26.6; 95.18 38]; plot_coast_xy(bbox,origin,'k')
ylim([-200 950])
   

%plot baselines (if computed)
if ~invert_vel

    dhat_base = Gbase*(srate-creep_rate)*1000; %convert rates to mm/yr

    figure; 
    subplot(121)
    hold on;  
    axis equal
    for k=1:size(BaseEnds,1)
        cline([BaseEnds(k,1) BaseEnds(k,3)],[BaseEnds(k,2) BaseEnds(k,4)],[Vbase(k) Vbase(k)])
    end
    title('Observed Baseline length changes rates')
    bbox = [73.04 26.6; 95.18 38]; plot_coast_xy(bbox,origin,'k')
    load cmap
    colormap(cmap)
    colorbar
    caxis([-max(abs(Vbase)) max(abs(Vbase))])

    subplot(122)
    hold on;  
    axis equal
    for k=1:size(BaseEnds,1)
        cline([BaseEnds(k,1) BaseEnds(k,3)],[BaseEnds(k,2) BaseEnds(k,4)],[dhat_base(k) dhat_base(k)])
    end
    title('Model Baseline length changes rates')
    bbox = [73.04 26.6; 95.18 38]; plot_coast_xy(bbox,origin,'k')
    load cmap
    colormap(cmap)
    colorbar
    caxis([-max(abs(Vbase)) max(abs(Vbase))])


end

%% plot velocities
if invert_vel

    figure   
    hold on;  
    scale=2;
    quiver(xysites(:,1),xysites(:,2),scale*Ve,scale*Vn,0,'b')
    axis equal
    title('Velocities')
    
    %velocities due to backslip 
    dhat_e = Ge*(srate-creep_rate)*1000; %convert to mm/yr
    dhat_n = Gn*(srate-creep_rate)*1000; %convert to mm/yr

%add long-term velocities
    if invert_vel
        dhat_e = dhat_e + lt_vels(:,1);
        dhat_n = dhat_n + lt_vels(:,2);
    end

    quiver(xysites(:,1),xysites(:,2),scale*dhat_e,scale*dhat_n,0,'r')
    bbox = [73.04 26.6; 95.18 38]; plot_coast_xy(bbox,origin,'k')

    legend('observed','model')


end

%% plot slip rate profiles

for k=1:size(profiles,1)

    v1 = [profiles(k,1:2) 0];
    v2 = [profiles(k,3:4) 0];
    c = [centroids(:,1:2) zeros(size(centroids,1),1)];
    d = point_to_line(c, v1, v2);

    ind = d<=intv;
    dist = sqrt((centroids(ind,1) - v1(1)).^2 + (centroids(ind,2) - v1(2)).^2);
    rate = creep_rate(ind);
    [y,i] = sort(dist);

    figure
    plot(dist(i),medfilt1(rate(i),10),'linewidth',2)
    title(['slip rate profile, x=' num2str(v1(1)) ', y=' num2str(v1(2))])
    xlabel('horizontal distance from end of profile (km)')
    ylabel('slip rate, m/yr')

end







