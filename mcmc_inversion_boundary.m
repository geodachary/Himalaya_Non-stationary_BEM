% addpath ./vel_field_nicholas/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%  BEGIN SETUP   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%number of boundary nodes along strike 
num_Dpts = 20;

%starting value -- upper locked boundary
locked_depths_U = 0*ones(num_Dpts,1);

%MCMC step size
stepsize_locked_depths_U = 6*ones(num_Dpts,1);

%starting value -- lower locked boundary
locked_depths_L = [5*ones(1,8) 12*ones(1,4) 12*ones(1,4) 20*ones(1,4)]';

%MCMC step size
stepsize_locked_depths_L = 6*ones(num_Dpts,1);

%starting value -- accumulated stress on ring around asperities (Pa)
ring_tau = 0*[0.15e6*ones(1,8) 0.25e6*ones(1,4) 0.05e6*ones(1,4) 0*0.05e6*ones(1,4)]';  %ringth width (km)

%MCMC step size
stepsize_ring_tau = 0.5e6*ones(num_Dpts,1);
% stepsize_ring_tau = 0e6*ones(num_Dpts,1); %if stationary mode

%continue sampling, or start new?
%setting to false overwrites output files
%setting to true appends existing files
continuing = false;

%time step for propgation 
dT = 10;

D = 3;  %process-zone width (km)

%elastic modulus
mu = 3e10;


%if inverting velocities, provide long-term velocity field in a file with
%format:
%lon lat Ve(mm/yr) Vn(mm/yr) (in same order as data file)  -- Here Ve and Vn are
%long-term velocities (not GPS-derived)
longterm_velocity_file = 'Himalaya_block_velocities.txt';


% folder for saving output files
folder_name = fullfile('MCMC_results', 'test_outputs_noringtau');

% Create the folder if it does not exist
if ~exist(folder_name, 'dir')
    mkdir(folder_name);
end


%%%%  END SETUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ---- Software configuration
% Make hmmvp available
addpath hmmvp0.16
addpath tools

[gamb.hmat.id nnz] = hm_mvp('init',gamb.hmat.savefn,4);



%load long-term velocities
if invert_vel
    lt_vels = load(longterm_velocity_file);
    lt_vels = lt_vels(:,3:4);  %columns are Ve, Vn
end

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

%slip rate in mm/yr
srate = rates/1000;
%adjust for stiffness
nel = size(el,1);
stiffness = self_rate;



%Greens function matrices  
if invert_vel

    G = [Ge;Gn];
    GG = [Ge./repmat(Sige,1,size(Ge,2)); Gn./repmat(Sign,1,size(Gn,2)) ];
    dd = [Ve./Sige; Vn./Sign];

else

    G = Gbase;
    GG = Gbase./repmat(Sigbase,1,size(Gbase,2));
    dd = Vbase./Sigbase;

end

% Change the Convergence Criteria
num = 1:size(el,1);
maxit = 200;
tol = 1e-4;

%need to scale stresses
%GFs computed with km, convert to meters
%GFs computed for mu=1
scale = mu*10^-3;



if continuing
eval(['load ./' folder_name '/M_locked_depths_U.txt'])
locked_depths_U =  M_locked_depths_U(end,:)';

eval(['load ./' folder_name '/M_locked_depths_L.txt'])
locked_depths_L =  M_locked_depths_L(end,:)';

eval(['load ./' folder_name '/M_ring_tau.txt'])
ring_tau =  M_ring_tau(end,:)';

end



%vector of unknown values
%starting values
X = [locked_depths_U; locked_depths_L; ring_tau];
stepsize = [stepsize_locked_depths_U; stepsize_locked_depths_L; stepsize_ring_tau];
%%%%%%


%forward model calculation
i_locked = get_boundary_indices(Dpts,locked_depths_U,locked_depths_L,centroids_rot);
i_ring = get_boundary_indices(Dpts,locked_depths_U-D,locked_depths_L+D,centroids_rot);

%include ring of accumulated stress
ring_stress = interp1(Dpts,ring_tau,centroids_rot(:,2));
ring_tau = zeros(size(el,1),1); %name reused: from here ring_tau is per-element, not per-node
ring_tau(i_ring) = ring_stress(i_ring);


asp_index = i_locked;
bslip = zeros(size(el,1),1);
bslip(asp_index) = srate(asp_index);  %mm/yr, backslip rate
dtau_rate = scale*hm_mvp('mvp',gamb.hmat.id,bslip);
dtau = dtau_rate*dT; %accumulated stress (stressing rate times T1)

%include ring stress
dtau = -dtau + ring_tau;

%compute creep rate 
rs = num(~asp_index);
s = gmres(@(x)mvp(x,gamb.hmat.id,rs,rs,~asp_index),-dtau(~asp_index)/scale,[],tol,maxit);
inc_slip = zeros(size(el,1),1);
inc_slip(~asp_index)=s; %increment of forward slip
U = srate*dT + inc_slip;  %srate*T1 to convert to displacement over time T1
U(asp_index)=0;  %locked asperity
creep_rate = U/dT*1000;  %convert from m/yr to mm/yr;


dhat = G*(rates-creep_rate);  %unweighted    %slip deficit drives the surface signal       
ddhat = GG*(rates-creep_rate);  %weighted by uncertainties

%add long-term velocities
if invert_vel
    dhat = dhat + [lt_vels(:,1);lt_vels(:,2)];    
    ddhat = ddhat + [lt_vels(:,1)./Sige; lt_vels(:,2)./Sign ];
end


%Gaussian log likelihood; constants dropped since only ratios matter
logrho = -.5*(dd-ddhat)'*(dd-ddhat);



%keep the last accepted state so a rejected step can be rolled back
Xprev = X;
logrhoprev = logrho;
dhatprev = dhat;
i_locked_prev = i_locked;
locked_depths_U_prev = locked_depths_U;
locked_depths_L_prev = locked_depths_L;
ring_tau_prev = ring_tau;
creep_rate_prev = creep_rate;


if ~continuing
%open with 'w' first to truncate any previous run
fid = fopen(['./' folder_name '/M_locked_depths_U.txt'],'w'); fclose(fid);
fid = fopen(['./' folder_name '/M_locked_depths_L.txt'],'w'); fclose(fid);
fid = fopen(['./' folder_name '/M_ring_tau.txt'],'w'); fclose(fid);
fid = fopen(['./' folder_name '/logrho.txt'],'w'); fclose(fid);
fid = fopen(['./' folder_name '/dhat.txt'],'w'); fclose(fid);
fid = fopen(['./' folder_name '/locked_index.txt'],'w'); fclose(fid);
fid = fopen(['./' folder_name '/creep_rates.txt'],'w'); fclose(fid);
end

%kept open for the whole run to avoid reopening every sample
fidM_U = fopen(['./' folder_name '/M_locked_depths_U.txt'],'a'); 
fidM_L = fopen(['./' folder_name '/M_locked_depths_L.txt'],'a'); 
fidM_tau = fopen(['./' folder_name '/M_ring_tau.txt'],'a'); 
fidrho = fopen(['./' folder_name '/logrho.txt'],'a'); 
fiddhat = fopen(['./' folder_name '/dhat.txt'],'a');
fid_locked = fopen(['./' folder_name '/locked_index.txt'],'a');
fid_creeping = fopen(['./' folder_name '/creep_rates.txt'],'a');


numparams = length(X);
numaccept = 0;

rand('state', sum(100*clock)); %different seed per run so parallel chains are independent
opts =  optimset('display','off');

% Increase the Number of Iterations

Niterations = 10*10^5;

for iter=1:Niterations

    %step through all params and vary one at a time
    count = mod(iter,numparams)+1; 
    r=(-1)^round(rand(1))*rand(1); %symmetric proposal, so the Metropolis ratio needs no correction term
    r=r*stepsize(count); 
    X(count)=X(count)+r;


    %%%%%%


    locked_depths_U = X(1:num_Dpts);
    locked_depths_L = X(1+num_Dpts:2*num_Dpts);
    ring_taus = X(1+2*num_Dpts:3*num_Dpts);
    
   
  %apply bounds  
  if sum(locked_depths_U>locked_depths_L | locked_depths_U<0 | ring_taus<0 | ring_taus>5e6 | locked_depths_L>40  )==0   %upper depth < lower depth
          
            
        %forward model calculation
        i_locked = get_boundary_indices(Dpts,locked_depths_U,locked_depths_L,centroids_rot);
        i_ring = get_boundary_indices(Dpts,locked_depths_U-D,locked_depths_L+D,centroids_rot);
        
        %include ring of accumulated stress
        ring_stress = interp1(Dpts,ring_taus,centroids_rot(:,2)); %node values to element values
        ring_tau = zeros(size(el,1),1);
        ring_tau(i_ring) = ring_stress(i_ring);
        
        
        asp_index = i_locked;
        bslip = zeros(size(el,1),1);
        bslip(asp_index) = srate(asp_index);  %mm/yr, backslip rate
        dtau_rate = scale*hm_mvp('mvp',gamb.hmat.id,bslip);
        dtau = dtau_rate*dT; %accumulated stress (stressing rate times T1)
        
        %include ring stress
        dtau = -dtau + ring_tau;
        
        %compute creep rate
        rs = num(~asp_index);
        s = gmres(@(x)mvp(x,gamb.hmat.id,rs,rs,~asp_index),-dtau(~asp_index)/scale,[],tol,maxit,[],[],inc_slip(~asp_index));
        inc_slip = zeros(size(el,1),1);
        inc_slip(~asp_index)=s; %increment of forward slip
        U = srate*dT + inc_slip;  %srate*T1 to convert to displacement over time T1
        U(asp_index)=0;  %locked asperity
        creep_rate = U/dT*1000;  %convert from m/yr to mm/yr
        
        
        dhat = G*(rates-creep_rate);  
        ddhat = GG*(rates-creep_rate);  


        if invert_vel
            dhat = dhat + [lt_vels(:,1);lt_vels(:,2)];    
            ddhat = ddhat + [lt_vels(:,1)./Sige; lt_vels(:,2)./Sign ];
        end


        logrho2 = -.5*(dd-ddhat)'*(dd-ddhat);
 
        %compared in log space to avoid underflow
        accept=metropolis_log(1,logrho,1,logrho2);

    else

        accept = 0;

    end
    

    if accept==1

       logrho=logrho2;
      

       Xprev=X;
       dhatprev = dhat;
       logrhoprev = logrho;
       
        i_locked_prev = i_locked;
        locked_depths_U_prev = locked_depths_U;
        locked_depths_L_prev = locked_depths_L;
        ring_tau_prev = ring_taus;
    
        creep_rate_prev = creep_rate;


       numaccept = numaccept+1;

    else
        %revert everything, otherwise the rejected model leaks into the next warm start
        X=Xprev;
        dhat=dhatprev;
        logrho=logrhoprev;
        i_locked = i_locked_prev;
        locked_depths_U = locked_depths_U_prev;
        locked_depths_L = locked_depths_L_prev;
        ring_taus = ring_tau_prev;
        creep_rate = creep_rate_prev;


    end

    
    if count==1 %thins the chain to one record per full parameter sweep


        fprintf(fidM_U,'\n',' ');fprintf(fidM_U,'%6.8f\t',locked_depths_U');
        fprintf(fidM_L,'\n',' ');fprintf(fidM_L,'%6.8f\t',locked_depths_L');

        fprintf(fidrho,'\n',' ');fprintf(fidrho,'%6.8f\t',logrho);
        fprintf(fiddhat,'\n',' ');fprintf(fiddhat,'%6.8f\t',dhat');
        
        fprintf(fidM_tau,'\n',' ');fprintf(fidM_tau,'%6.8f\t',ring_taus');
        fprintf(fid_locked,'\n',' ');fprintf(fid_locked,'%6.8f\t',i_locked');

        fprintf(fid_creeping,'\n',' ');fprintf(fid_creeping,'%6.8f\t',creep_rate');


        disp(['Completed sample number ' num2str(iter)])
        disp(['Acceptance rate: ' num2str(numaccept/iter*100)]) %watch this to retune stepsize

                
        disp(['Log likelihood: ' num2str(logrho)])
                
    
    end
    
end



