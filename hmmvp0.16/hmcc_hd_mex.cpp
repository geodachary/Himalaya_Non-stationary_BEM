#include <iostream>
#include "mex.h"
#include "Hd.hpp"
#include "HmmvpMexUtil.hpp"

#if 0
// I discovered that Matlab struct arrays take ~15x as much memory as one would
// think they should in this situation. Hence I can't use them here for large
// problems.
void HdToMex(const Hm::Hd& hd, mxArray*& mbs, mxArray*& mp, mxArray*& mq)
{
  const char* flds[] = {"r0", "m", "c0", "n"};
  mbs = mxCreateStructMatrix(1, hd.NbrBlocks(), 4, flds);
  int k = 0;
  for (Hm::Hd::iterator it = hd.begin(), end = hd.end();
       it != end; ++it, ++k) {
    SetField(mbs, k, flds[0], it->r0);
    SetField(mbs, k, flds[1], it->m);
    SetField(mbs, k, flds[2], it->c0);
    SetField(mbs, k, flds[3], it->n);
  }
  Vecui p, q;
  hd.Permutations(p, q);
  mp = VectorToMex(p);
  mq = VectorToMex(q);
}
#else
// Each column has [r0 c0 nr nc].
void HdToMex(const Hm::Hd& hd, mxArray*& mbs, mxArray*& mp, mxArray*& mq)
{
  mbs = mxCreateDoubleMatrix(4, hd.NbrBlocks(), mxREAL);
  double* bs = mxGetPr(mbs);
  int k = 0;
  for (Hm::Hd::iterator it = hd.begin(), end = hd.end();
       it != end; ++it, k += 4) {
    bs[k  ] = it->r0;
    bs[k+1] = it->c0;
    bs[k+2] = it->m;
    bs[k+3] = it->n;
  }
  Vecui p, q;
  hd.Permutations(p, q);
  mp = VectorToMex(p);
  mq = VectorToMex(q);
}
#endif

// Implements
//     [bs p q] = hmcc_hd(D,R,[split_method],[eta])
// D corresponds to the *column*, not row, coordinates. I might switch this in
// the future; the original logic is that the matrix maps slip in space D to
// stress in R: ie, D -> R.
void mexFunction(int nlhs, mxArray **plhs, int nrhs, const mxArray **prhs)
{
  if (nrhs < 1 || nlhs != 3)
    mexErrMsgTxt("Calling sequence is [bs p q] = "
		 "hmcc_hd(D,R,[split_method],[eta]).");

  Matd D, R;
  MexToMatrix(prhs[0], D);
  if (D.Size(1) != 3) mexErrMsgTxt("D is 3xN.");
  MexToMatrix(prhs[1], R);
  if (R.Size(1) != 3) mexErrMsgTxt("R is 3xN.");
 
  Hm::HdCc::SplitMethod sm = Hm::HdCc::splitMajorAxis;
  if (nrhs > 2) {
    int smi = (int)mxGetScalar(prhs[2]);
    if (smi != 0 && smi != 1)
      mexErrMsgTxt("split_method must be 0 (major axis) or 1 (axis-aligned).");
    if (smi == 1) sm = Hm::HdCc::splitAxisAligned;
  }

  double eta = 3.0;
  if (nrhs > 3) {
    eta = mxGetScalar(prhs[3]);
    if (eta <= 0.0)
      mexErrMsgTxt("eta must be > 0.");
  }

  Hm::HdCcNonsym hd(D, R);
  hd.Opt_eta(eta);
  hd.Opt_split_method(sm);
  if (sm == Hm::HdCc::splitAxisAligned)
    hd.Opt_distance_method(Hm::HdCc::distAxisAligned);
  else
    hd.Opt_distance_method(Hm::HdCc::distCentroid);
  hd.Compute();
  HdToMex(hd, plhs[0], plhs[1], plhs[2]);
}
