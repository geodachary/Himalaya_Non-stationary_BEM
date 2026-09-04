%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT 

% Block model: use 'lindsey_2018' (Lindsey et al., 2018) or 'Block_I', 'Block_II', or 'Block_III' (Panda & Lindsey, 2024).
block_id = 'lindsey_2018';

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
disp_filename = 'test_disp';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add folder path for input data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath tools
addpath vel_field_lindsey/
addpath hmat/
addpath ne_10m_coastline/

%name of file containing data
%format of columns: lon, lat, Ve(mm/yr), Vn(mm/yr), Sige, Sign

data_filename = fullfile('vel_field_lindsey', 'observed_vel_subset.txt'); % Data from Lindsey et al. (2018)


%path to build mesh script
mesh_path = './build_mesh/make_mesh_elastic.m';





%build mesh
run(mesh_path)


%% load and plot data

%load GPS data
data = load(data_filename);
xysites = llh2local(data(:,1:2)',origin)';
Ve = data(:,3);
Vn = data(:,4);


%% Add floor mat on the gps velocities data
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


%% If you want to try with more Uncertainties value by Hardwired assumptions
% Sige = data(:,5);
% Sign = data(:,6);
% 
% Sige = Sige*4;
% Sign = Sign*4;


%% compute baselines if requested
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
        cline([BaseEnds(k,1) BaseEnds(k,3)],[BaseEnds(k,2) BaseEnds(k,4)],[Vbase(k) Vbase(k)])
    end
    title('Baseline Baseline length changes rates')
    bbox = [73.04 26.6; 95.18 38]; plot_coast_xy(bbox,origin,'k')
    load cmap
    colormap(cmap)
    colorbar
    caxis([-max(abs(Vbase)) max(abs(Vbase))])
end


%% compute displacement GFs or load file

if compute_disp  

    if invert_vel
        %velocity GFs
        [Ge,Gn,Gu] = make_dispG_triangular(el,nd,[xysites zeros(size(xysites,1),1)]',rakes,[]);
    else
        %baseline GFs
        Gbase = Get_Gs_baselines(el,nd,xysites,rakes,baselines,Vec_unit,crossind,L);
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
self_rate = self_rate/max(self_rate);




