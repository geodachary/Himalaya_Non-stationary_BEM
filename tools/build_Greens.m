%This computes the interface mesh and displacement and stress Greens functions. This is slow 
%and needs to be run only after changes to the mesh or locations of
%observation coordinates.

%Need to load observation coordinates 
%xysites = Mx2 vector of [x,y] coordinates of M GPS sites
%load interseismic_1998

%load slipvectors  %slipvectors from block model (will assume backslip in this direction)
%get unit vector in slip direction
%n = sqrt(Vp.^2 + Vn.^2);
%vs = Vp./n;  %strike component
%vd = Vn./n;  %dipcomponent


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  FAST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%First build the triangular mesh for the subduction interface


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  SLOW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute stress and displacement greens functions




% ---- Software configuration
% Make hmmvp available
addpath hmmvp0.16;


model_name = 'Hik_intv3';

global gamb;
gamb.use = 1;
% Optional H-matrix compression of the BEM matrix:
% Use H-matrix compression as implemented in hmmvp?
gamb.hmat.use = 1;
% Element-wise relative error of approximation. Should be ~ RelTol in the
% ODE integration.
gamb.hmat.rerr = 1e-3;
% Directory to which to save H-mat for future use:
dr = '.';


% ODE software parameters:
% File to which to save ODE solution
ne = size(el,1);
fn_dec = sprintf('ne%dhmat%drerr%3.2f',...
       ne,gamb.hmat.use,log10(gamb.hmat.rerr));
gamb.ode.savefn = [dr '/ode_' fn_dec model_name];

% Save at every k'th ODE time step
gamb.ode.saveEvery = 10; 

% H-matrix file name. 
gamb.hmat.savefn = sprintf('%s/Hmat_%s.dat',dr,fn_dec);

  
%make stress-slip matrix, stressG, such that stress = stressG*slip;
CalcG_rake('Make',el,nd,rakes,gamb.hmat.rerr,gamb.hmat.savefn);






%save outptuts to mat file
save setup_mesh