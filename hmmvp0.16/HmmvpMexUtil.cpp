#include <algorithm>
#include "mex.h"
#include "Hd.hpp"
#include "HmmvpMexUtil.hpp"
using namespace std;

template<> mxArray* MatrixToMex<double>(const Matrix<double>& m)
  throw (OutOfMemoryException)
{
  int nr = m.Size(1), nc = m.Size(2);
  mxArray* ma = mxCreateDoubleMatrix(nr, nc, mxREAL);
  if (!ma) throw OutOfMemoryException("MatrixToMex");
  memcpy(mxGetPr(ma), m.GetPtr(), nr*nc*sizeof(double));
  return ma;
}

template<> void MexToMatrix<double>(const mxArray* ma, Matrix<double>& m)
  throw (OutOfMemoryException)
{
  int nr = mxGetM(ma), nc = mxGetN(ma);
  m.Resize(nr,nc);
  memcpy(m.GetPtr(), mxGetPr(ma), nr*nc*sizeof(double));
}

// [U S V] = svd(A,'econ')
void MexSvd(const Matd& A, Matd* U, Matd* s, Matd* V)
{
  if (A.Size() == 0) {
    if (U) U->Resize(0);
    if (s) s->Resize(0);
    if (V) V->Resize(0);
    return;
  }
  int m = A.Size(1), n = A.Size(2), k = min(m, n);
  mxArray *mplhs[3], *mprhs[2];
  mprhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
  memcpy(mxGetPr(mprhs[0]), A.GetPtr(), m*n*sizeof(double));
  mprhs[1] = mxCreateString("econ");
  mexCallMATLAB(3, mplhs, 2, mprhs, "svd");
  mxDestroyArray(mprhs[0]);
  mxDestroyArray(mprhs[1]);
  mxArray *mU = mplhs[0], *mS = mplhs[1], *mV = mplhs[2];
  if (U) {
    U->Resize(m, k);
    memcpy(U->GetPtr(), mxGetPr(mU), m*k*sizeof(double));
  }
  mxDestroyArray(mU);
  if (s) {
    s->Resize(k);
    double* ps = s->GetPtr();
    double* pmS = mxGetPr(mS);
    for (int i = 0; i < k; i++) ps[i] = pmS[i*(k + 1)];
  }
  mxDestroyArray(mS);
  if (V) {
    V->Resize(n, k);
    memcpy(V->GetPtr(), mxGetPr(mV), n*k*sizeof(double));
  }
  mxDestroyArray(mV);
}
