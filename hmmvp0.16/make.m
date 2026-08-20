% Configure options
isunix = 1;
use_omp = 0;

% Not sure I have all the flags right. The UNIX part works on my Ubuntu
% 64-bit 4-CPU system with use_omp = 0 or 1. At least one 32-bit Windows user
% has built the code successfully using this script.
p = '.';
flags = '-O -DHMMVP_MEXSVD';
if(use_omp) flags = [flags ' -DHMMVP_OMP ']; end
if(isunix)
  if(use_omp)
    flags = [flags [' CXXFLAGS="\$CXXFLAGS -fopenmp" ',...
		    'LDFLAGS="\$LDFLAGS -fopenmp"']];
  end
  flags = [flags ' CXXLIBS="\$CXXLIBS -lmwblas"'];
  mc = sprintf('mex -f %s/mexopts.sh -outdir %s -largeArrayDims %s',p,p,flags);
else
  % Definitely not sure about the Windows stuff.
  use_omp = 0;
  flags = [flags ' CXXLIBS="$CXXLIBS -lmwblas'];
  if(~isempty(dir(sprintf('%s/mexopts.bat',p))))
    mc = sprintf('mex -f %s/mexopts.bat -outdir %s %s',p,p,flags);
  else
    mc = sprintf('mex -outdir %s %s',p,flags);
  end
end

%changed by KMJ, July 21, 2022, because mexopts.sh file no longer supported
%eval(sprintf('%s %s/hm_mvp.cpp',mc,p));
%eval(sprintf('%s %s/hmcc_hd_mex.cpp %s/Hd.cpp %s/HmmvpMexUtil.cpp',mc,p,p,p));
mex  -outdir . -largeArrayDims -O -DHMMVP_MEXSVD ./hm_mvp.cpp -lmwblas
mex  -outdir . -largeArrayDims -O -DHMMVP_MEXSVD ./hmcc_hd_mex.cpp ./Hd.cpp ./HmmvpMexUtil.cpp -lmwblas
 
 