close all;

% % Select result type
% result_type = 'visco';   % 'elastic' or 'visco'
% 
% % % Define the main folder containing MCMC results
% if strcmpi(result_type, 'elastic')
%     results_folder = 'Elastic_results';
% else
%     results_folder = 'Viso_results';
% end
% 
% %folder name with output files
% results_name = 'visco_block_II_1000eta';
% folder_name = fullfile('MCMC_results', results_folder, results_name);



folder_name = 'test_outputs_noringtau';

%discard burn-in samples
discard = 10;

eval(['load ./' folder_name '/M_locked_depths_U.txt'])
eval(['load ./' folder_name '/M_locked_depths_L.txt'])
eval(['load ./' folder_name '/M_ring_tau.txt'])
eval(['load ./' folder_name '/dhat.txt'])
eval(['load ./' folder_name '/logrho.txt'])
eval(['load ./' folder_name '/locked_index.txt'])
eval(['load ./' folder_name '/creep_rates.txt']) 


figure
plot(logrho)

% %toss out burn-in samples
M_locked_depths_U(1:discard,:) = [];
M_locked_depths_L(1:discard,:) = [];
M_ring_tau(1:discard,:) = [];
dhat(1:discard,:) = [];
locked_index(1:discard,:) = [];
creep_rates(1:discard,:) = [];


% % Optionally select the sample range for subsequent processing and plotting.
% M_locked_depths_U = M_locked_depths_U(1:9200,:); 
% M_locked_depths_L = M_locked_depths_L(1:9200,:); 
% M_ring_tau = M_ring_tau(1:9200,:); 
% dhat = dhat(1:9200,:); 
% locked_index = locked_index(1:9200,:); 
% creep_rates = creep_rates(1:9200,:); 



%number of boundary nodes along strike 
num_Dpts = 20;

D = 3;  %process-zone width (km)

%elastic modulus
mu = 3e10;

%time step for propagation 
dT = 10;

% Generate triangular patch data and extract face centroids and strike angles
patch_stuff = make_triangular_patch_stuff(el,nd);
centroids = patch_stuff.centroids_faces;
strikes = patch_stuff.strike_faces;

%construct boundary curves

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


k=1;
[c_U,c_L] = get_boundary_curves(Dpts,M_locked_depths_U(k,:)',M_locked_depths_L(k,:)',centroids,centroids_rot);
boundary_curves_U_x = zeros(size(c_U,1),size(dhat,1));
boundary_curves_U_y = zeros(size(c_U,1),size(dhat,1));
boundary_curves_U_z = zeros(size(c_U,1),size(dhat,1));
boundary_curves_L_x = zeros(size(c_U,1),size(dhat,1));
boundary_curves_L_y = zeros(size(c_U,1),size(dhat,1));
boundary_curves_L_z = zeros(size(c_U,1),size(dhat,1));


% %% Construct the processing zone boundaries
all_ring_tau = zeros(length(el),size(M_locked_depths_U, 1));

for k = 1:size(M_locked_depths_U, 1)

    [c_U, c_L] = get_boundary_curves(Dpts, M_locked_depths_U(k,:)', M_locked_depths_L(k,:)', centroids, centroids_rot);

    c_U(:,2) = c_U(:,2) + min(centroids_rot(:,2));
    c_L(:,2) = c_L(:,2) + min(centroids_rot(:,2));

    % Shift and rotate boundary nodes to true position
    str = -mean(strikes) * pi / 180; % Negative sign to rotate back
    R = [cos(str) -sin(str); sin(str) cos(str)];
    c_U_rot = (R * c_U(:,1:2)')';
    c_L_rot = (R * c_L(:,1:2)')';

    % Dynamically adjust array sizes to match c_U and c_L sizes for this iteration
    boundary_curves_U_x(1:size(c_U, 1), k) = c_U_rot(:,1);
    boundary_curves_U_y(1:size(c_U, 1), k) = c_U_rot(:,2);
    boundary_curves_U_z(1:size(c_U, 1), k) = c_U(:,3);

    boundary_curves_L_x(1:size(c_L, 1), k) = c_L_rot(:,1);
    boundary_curves_L_y(1:size(c_L, 1), k) = c_L_rot(:,2);
    boundary_curves_L_z(1:size(c_L, 1), k) = c_L(:,3);


    %construct ring stress
    %node points along strike
    Dpts = linspace(0,max(centroids_rot(:,2)),num_Dpts)';

    i_locked = get_boundary_indices(Dpts,M_locked_depths_U(k,:)',M_locked_depths_L(k,:)',centroids_rot);

    i_ring = get_boundary_indices(Dpts,M_locked_depths_U(k,:)'-D,M_locked_depths_L(k,:)'+D,centroids_rot);

    ring_stress = interp1(Dpts, M_ring_tau(k,:)',centroids_rot(:,2));


    %include ring of accumulated stress
    ring_tau = zeros(size(el,1),1);
    ring_tau(i_ring) = ring_stress(i_ring);

    ring_tau(i_locked) = 0;
    all_ring_tau(:,k) = ring_tau;


end


%compute stressing rates
for k = 1:size(M_locked_depths_U, 1)
    scale = mu*10^-3;
    stressing_rates(:,k) = hm_mvp('mvp',gamb.hmat.id,creep_rates(k,:)');
end

%specify start and end points of slip rate profiles across interface
%each row is a profile:  start x, start y, end x, end y
profiles = [-775 550 -630 685; -520 250 -350 450; -235 75 -60 265; 400 -120 475 125; 810 -122 890 130; 1070 -5 1060 225];


% Plot code for the Locking Boundaries

figure
% trimesh(el,nd(:,1),nd(:,2), 'color',[0.8 0.8 0.8])
axis equal
axis tight  
set(gca,'fontsize',15)
title('Locking Boundaries')
hold on
ylim([200 1200])

plot(boundary_curves_L_x,boundary_curves_L_y,'Color',[0.75 0 0 0.005], 'LineWidth',1)
hold on
plot(boundary_curves_U_x,boundary_curves_U_y,'Color',[0 0 .5 0.005], 'LineWidth',1)



% Plot width of locked zone
str = mean(strikes)*pi/180;
R = [cos(str) -sin(str); sin(str) cos(str)];
nd_rot = (R*nd(:,1:2)')';
nd_rot(:,3) = nd(:,3);
nd_rot(:,2) = nd_rot(:,2)-min(nd_rot(:,2));
for k=1:size(M_locked_depths_U,1)
    [U_xyz{k}, widths(:,k)] = get_widths(Dpts,M_locked_depths_U(k,:)',M_locked_depths_L(k,:)',centroids_rot, nd_rot);
end
figure
plot(U_xyz{end}(:,2),widths,'b','Color',[0 0 .5 0.01])
xlabel('distance along the fault (km)')
ylabel('width of locked zone (km)')
set(gca,'fontsize',15)
title('Width of Locked Zone')


% Plot code for the Probability of Locking
figure
h=trisurf(el,nd(:,1),nd(:,2),nd(:,3),mean(locked_index),'edgecolor','none'); colorbar;  daspect([1,1,1])
colormap(jet)
colorbar
axis equal
axis tight      
set(gca,'fontsize',15)
title('Probability of Locking')
view(2)
hold on
bbox = [73.04 26.6; 95.18 38]; plot_coast_xy(bbox,origin,'k')
ylim([-200 900])
   

% Plot code for the Mean Ring Stress
ringtau_actual = mean(all_ring_tau,2)/dT;

% Compute mean ring stress in kPa/year
ringtau_actual_kpa = ringtau_actual/1000;

figure
h=trisurf(el,nd(:,1),nd(:,2),nd(:,3),ringtau_actual_kpa,'edgecolor','none'); colorbar;  daspect([1,1,1])
colormap(jet)
colorbar
axis equal
axis tight      
set(gca,'fontsize',15)
title('Mean Ring Stress in kPa/year')
view(2)
% caxis([0 max(ringtau_actual_kpa)/10])
hold on
bbox = [73.04 26.6; 95.18 38]; plot_coast_xy(bbox,origin,'k')
ylim([-200 950])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot mean creep rate %%%%%%%%%%%%%%%%%%%%%

figure
h=trisurf(el,nd(:,1),nd(:,2),nd(:,3),mean(creep_rates),'edgecolor','none'); colorbar;  daspect([1,1,1])
colormap(jet)
colorbar
axis equal
axis tight      
clim([0 max(mean(creep_rates)/2)])
set(gca,'fontsize',15)
title('Creep Rate (mm/yr)')
view(2)
hold on
bbox = [73.04 26.6; 95.18 38]; plot_coast_xy(bbox,origin,'k')
ylim([-200 900])
plot([profiles(:,1) profiles(:,3)]', [profiles(:,2) profiles(:,4)]','linewidth',2)


%% % Plot slip rate profiles based on all rows of creep_rates, creating a new figure for each profile
for k = 1:size(profiles, 1)
    figure
    hold on
    
    for row = 1:size(creep_rates, 1)
        v1 = [profiles(k, 1:2) 0];
        v2 = [profiles(k, 3:4) 0];
        c = [centroids(:, 1:2) zeros(size(centroids, 1), 1)];
        d = point_to_line(c, v1, v2);

        ind = d <= intv;
        dist = sqrt((centroids(ind, 1) - v1(1)).^2 + (centroids(ind, 2) - v1(2)).^2);
        rate = creep_rates(row, ind);
        [y, i] = sort(dist);

        plot(dist(i), medfilt1(rate(i), 10), 'linewidth', 2)
    end

    title(['Creep Rate Profile, x=' num2str(v1(1)) ', y=' num2str(v1(2))])
    xlabel('Horizontal Distance from End of Profile (km)')
    ylabel('Creep Rate (mm/yr)')
    hold off
end


%% %%%%%%%%%%%% For Observed and Modelled mean baseline plot of all data%%%%%%%%%%%%%%%%%%%

% plot baselines (if computed)
if ~invert_vel
     dhat_base = mean(dhat,1); % Convert rates to mm/yr
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
    caxis([-max(abs(Vbase))/2 max(abs(Vbase))/2])

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
    caxis([-max(abs(Vbase))/2 max(abs(Vbase))/2])





%% %%%%%%%%%%%% For Reduced Chi-squared plot of all data in histogram%%%%%%%%%%%%%%%%%%%%%%%%%%
    dd = Vbase./Sigbase;
    mean_dhat = mean(dhat, 1);
    weighted_residuals = dd' - mean_dhat./Sigbase';
    
    p = num_Dpts * 2; % Number of parameters
    N = length(weighted_residuals);
    chi_square = sum(weighted_residuals.^2);
    reduced_chi_square = chi_square / (N - p);
    
    % Plot histogram with normal fit
    figure;
    histogram(weighted_residuals, 'Normalization', 'pdf', 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'k');
    hold on;
    pd = fitdist(weighted_residuals', 'Normal');  
    x_values = linspace(min(weighted_residuals), max(weighted_residuals), 100);
    y_values = pdf(pd, x_values);
    plot(x_values, y_values, 'k-', 'LineWidth', 1.5);
    title(sprintf('Reduced \\chi^2 = %.2f', reduced_chi_square), 'FontSize', 14);
    xlabel('Residuals (mm/yr)', 'FontSize', 14);
    ylabel('Probability Density', 'FontSize', 14);
    hold off;

end


%% %%%%%%%%%%%% For Observed and Modelled velocities plot of all data%%%%%%%%%%%%%%%%%%%

if invert_vel
    % Number of data points
    N_e = length(dhat_e); % 215
    N_n = length(dhat_n); % 215
    
    % Split dhat into two parts for modelled data
    dhat_e_model = dhat(:, 1:N_e);
    dhat_n_model = dhat(:, N_e+1:N_e+N_n);
    
    % Plotting observed and modelled velocities
    figure;
    hold on;
    scale = 2;
    
    % Plot observed velocities
    quiver(xysites(:,1), xysites(:,2), scale*Ve, scale*Vn, 0, 'b');
    
    % Ensure axis equal
    axis equal;
    title('Velocities');
    
    % Plot all modelled velocities
    for i = 1:size(dhat_e_model, 1)
        quiver(xysites(:,1), xysites(:,2), scale*dhat_e_model(i, :)', scale*dhat_n_model(i, :)', 0, 'r');
    end

    % Define and plot the bounding box
    bbox = [73.04 26.6; 95.18 38]; % Example bounding box
    plot_coast_xy(bbox, origin, 'k');
    
    % Add legend
    legend('observed', 'model');


    %%%%%%%%%%%%%% Observed and Modelled mean velocitiy plot of the data %%%%%%%%%%%%%%%%%%%%%%%
    
    % Number of data points
    N_e = length(dhat_e); % 215
    N_n = length(dhat_n); % 215
    
    % Split dhat into two parts for modelled data
    dhat_e_model = dhat(:, 1:N_e);
    dhat_n_model = dhat(:, N_e+1:N_e+N_n);
    
    
    % Calculate mean of columns for dhat_e_model and dhat_n_model
    mean_dhat_e_model = mean(dhat_e_model, 1); % Mean across columns
    mean_dhat_n_model = mean(dhat_n_model, 1); % Mean across columns
    
    % Plotting observed and modelled velocities
    figure;
    hold on;
    scale = 2;
    quiver(xysites(:,1), xysites(:,2), scale*Ve, scale*Vn, 0, 'b');
    axis equal;
    title('Mean Velocities');
    
    % Use the mean modelled data for quiver
    quiver(xysites(:,1), xysites(:,2), scale*mean_dhat_e_model', scale*mean_dhat_n_model', 0, 'r');
    bbox = [73.04 26.6; 95.18 38]; % Example bounding box
    plot_coast_xy(bbox, origin, 'k');
    
    legend('observed', 'model');



    %% %%%%%%%%%%%% For Reduced Chi-squared plot of all data in histogram%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Number of data points
    N_e = length(dhat_e);
    N_n = length(dhat_n);
    
    % Split dhat into two parts for modeled data
    dhat_e = mean(dhat(:, 1:N_e))';
    dhat_n = mean(dhat(:, N_e+1:N_e+N_n))';
    dd = [Ve./Sige; Vn./Sign];
    
    ddhat = [dhat_e./Sige;dhat_n./Sign];  %weighted model vector
    
    residuals = dd - ddhat; %weighted residuals
    
    % Calculate the combined reduced chi-square value
    p = num_Dpts*2; % Number of parameters in the model
    N = length(dd); % Total number of residuals
    chi_square = sum((residuals).^2);
    reduced_chi_square = chi_square / (N - p);
    
    
    
    % Plotting the histogram of residuals
    figure;
    histogram(residuals, 'Normalization', 'pdf', 'FaceColor', [0.7 0.7 0.7]);
    hold on;
    
    % Fit a normal distribution to the residuals
    pd = fitdist(residuals, 'Normal');
    x_values = -8:0.1:8; % Adjusted range to better fit the distribution
    y_values = pdf(pd, x_values);
    plot(x_values, y_values, 'k-', 'LineWidth', 1.5);
    
    % Adding the reduced chi-square value to the plot
    title(sprintf('reduced \\chi^2 %.2f', reduced_chi_square), 'FontSize', 14);
    xlabel('Residuals (mm/yr)', 'FontSize', 14);
    ylabel('Counts', 'FontSize', 14);
    hold off;

end

%%

plot_moment_rates = true;

if plot_moment_rates
    moment_rates_plot;
else

end