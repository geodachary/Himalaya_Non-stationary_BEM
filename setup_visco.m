clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT 
addpath Green_function_visco/
addpath visco_mesh/
addpath hmat/

% Block model: use 'himalaya_visco' (Lindsey et al., 2018) or 'visco_block_I', 'visco_block_II', or 'visco_block_III' (Panda & Lindsey, 2024).
block_id = 'himalaya_visco';

% Load the mesh from visco3d
mesh_file = fullfile('visco_mesh', sprintf('%s_mesh_inversion.mat', block_id));
load(mesh_file);

% Set the folder where Green's functions are stored
green_function_folder = 'Green_function_visco';

% Generate the name of the Green's function file or folder
% Viscosity models:
% 10eta       : Mantle = 10^19 Pa·s, Lower crust = 10^19 Pa·s
% 100eta      : Mantle = 10^19 Pa·s, Lower crust = 10^20 Pa·s
% 1000eta     : Mantle = 10^19 Pa·s, Lower crust = 10^21 Pa·s
% 100-1000eta : Mantle = 10^20 Pa·s, Lower crust = 10^21 Pa·s
gf_name = sprintf('%s_interseismic_vel_10eta', block_id);


% Create the full path to the desired folder
gf_path = fullfile(green_function_folder, gf_name);

load([gf_name, '.mat']);

Ge = Ge_inter;
Gn = Gn_inter;
Gu = Gu_inter;

%invert velocities or baselines?
%true for velocities, false for baselines
invert_vel = false; 

%compute stress H-matrix?
% false uses existing h-matrix
compute_hmat = false; 

%compute displacement Greens functions?
%false uses existing displacement GFs save in 
%matlab file with name disp_filename
compute_disp = true;
% disp_filename = 'test_disp';


%name of file containing data
%format of columns: lon, lat, Ve(mm/yr), Vn(mm/yr), Sige, Sign
data_filename = 'observed_vel_subset.txt';


%path to build mesh script
mesh_path = './build_mesh/make_mesh_visco.m';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath tools
addpath ne_10m_coastline/


if any(strcmpi(block_id, {'visco_block_I'}))
    hsub = 'hmat_I';
elseif any(strcmpi(block_id, {'visco_block_II'}))
    hsub = 'hmat_II';
elseif any(strcmpi(block_id, {'visco_block_III'}))
    hsub = 'hmat_III';
elseif any(strcmpi(block_id, {'himalaya_visco'}))
    hsub = 'hmat_ES';
else
    error('Unknown block name: %s', block_id);
end

addpath(fullfile('hmat', hsub));


%build mesh
run(mesh_path)


%% load and plot data


%load GPS data
data = load(data_filename);
xysites = llh2local(data(:,1:2)',origin)';
Ve = data(:,3);
Vn = data(:,4);

if size(data,2) == 6
    Sige = data(:,5);
    Sign = data(:,6);

    % Apply the condition to each element of Sige and Sign
    % Multiply by 2 only if the value is less than 1
    Sige(Sige < 1) = Sige(Sige < 1) .* 2;
    Sign(Sign < 1) = Sign(Sign < 1) .* 2;


else  %if no uncertainties, assign value of 1
    Sige = ones(size(Ve));
    Sign = ones(size(Vn));
end


%compute baselines if requested
if ~invert_vel

    %find baselines that cross elements breaking the ground surface (need
    %to integrate strain rates across these baselines)
    nodes_z = [nd(el(:,1),3) nd(el(:,2),3) nd(el(:,3),3)];  %depth of nodes in all triangles
    nodes_x = [nd(el(:,1),1) nd(el(:,2),1) nd(el(:,3),1)];  %x position of nodes in all triangles
    nodes_y = [nd(el(:,1),2) nd(el(:,2),2) nd(el(:,3),2)];  %y position of nodes in all triangles
    zero_nodes = nodes_z==0;
    surf_break = sum(zero_nodes,2)==2; %find triangles with two nodes at surface

    %make patch endpoints for elements that break surface
    PatchEnds = nan(size(el,1),4);
    for j=1:size(PatchEnds,1)
        if surf_break(j)
            PatchEnds(j,[1 3]) = nodes_x(j,zero_nodes(j,:));
            PatchEnds(j,[2 4]) = nodes_y(j,zero_nodes(j,:));
        end
    end
        
    [Vbase,Sigbase,BaseEnds,L,baselines,Vec_unit,crossind] = make_baseline_rate_changes(xysites,Ve,Vn,Sige,Sign,PatchEnds);

end


%plot velocities
figure; hold on;  
quiver(xysites(:,1),xysites(:,2),Ve,Vn)
axis equal
title('Observed velocities')
bbox = [73.04 26.6; 95.18 38]; plot_coast_xy(bbox,origin,'k')


%plot baselines (if computed)
if ~invert_vel
    figure; hold on;  
    axis equal
    for k=1:size(BaseEnds,1)
        cline([BaseEnds(k,1) BaseEnds(k,3)],[BaseEnds(k,2) BaseEnds(k,4)],[Vbase(k)/L(k) Vbase(k)/L(k)])
    end
    title('Baseline elongation rates')
    bbox = [73.04 26.6; 95.18 38]; plot_coast_xy(bbox,origin,'k')
    load cmap
    colormap(cmap)
    colorbar
    caxis([-max(abs(Vbase./L)) max(abs(Vbase./L))])
end


%% compute displacement GFs or load file

if compute_disp  

    if invert_vel
        %velocity GFs from visco 3d
        Ge = Ge_inter;
        Gn = Gn_inter;
        Gu = Gu_inter;

        disp(['loaded the displacement Greens function calculations from Visco3d.'] )

    else
        % % load the velocity GFs function from visco 3d
        Ge = Ge_inter;
        Gn = Gn_inter;
        Gu = Gu_inter;

        %convert velocity GFs to baseline GFs
        Gbase = Get_Gs_baselines_visco(Ge, Gn, xysites, baselines, Vec_unit, crossind, L, el, nd, rakes);

    end

else
    % load precomputed displacement file
    load(disp_filename)

end




%% Hmatrix setup
setup_Hmat   %this is a script


%for computing stiffness
[gamb.hmat.id nnz] = hm_mvp('init',gamb.hmat.savefn,4);

for k=1:size(el,1)
    
    r = zeros(size(el,1),1); r(k) = 1;    
    sr = -hm_mvp('mvp',gamb.hmat.id,r);
    self_rate(k) = sr(k);

    if mod(k,100)==0
    disp(['Completed ' num2str(k/size(el,1)*100) '% of stiffness calculations.'])
    end
    
end


