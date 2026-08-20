function [bs p q] = hmcc_hd(D,R,varargin)
% Spatially decompose domain and range point clouds in preparation to form an
% H-matrix.
%
% [bs p q] = hmcc_hd(D,R,[split_method],[eta])
%   D is a 3xnD matrix of points and similarly for R. D corresponds to the
% influence points: in other words, the columns of the matrix. R corresponds
% to the influenced points: in other words, the rows. In many cases, D = R,
% in which case pass the same matrix twice.
%   bs,p,q are inputs to hm('Compress'). These are the H-matrix blocks and
% the permutations associated with the range and domain spaces.
%   The optional parameter split_method is either
%     0 [default]: Split normal to the major axis of a point cloud. This is a
%        good value unless your problem is a rectangle or cube.
%     1: Split aligned with Cartesian axes.
%   The optional parameter eta > 0 controls admissibility of blocks. It
% defaults to 3, which seems like a good value.

% AMB ambrad@cs.stanford.edu
  
  CheckX(D,'D');
  CheckX(R,'R');
  [bs p q] = hmcc_hd_mex(D,R,varargin{:});

function CheckX(X,name)
  if (size(X,1) ~= 3)
    error(sprintf('%s is not 3xN',name));
  end
  [foo p] = sort(X(1,:));
  if (any(sum(abs(diff(X(:,p),1,2))) == 0))
    error(sprintf('%s contains duplicate points.',name));
  end
  