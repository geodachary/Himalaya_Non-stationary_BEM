#include <algorithm>
#include "HmmvpUtil.hpp"

namespace Hm {

#ifdef HMMVP_LAPACK
extern "C" void dgesvd_
(char*, char*, int*, int*, double*, int*, double*, double*, int*,
 double*, int*, double*, int*, int*);

void LapackSvd(const Matd& A, Matd& U, Matd& s, Matd& Vt)
  throw (LapackException)
{
  if (A.Size() == 0)
    throw LapackException("LapackSvd: A.Size() == 0");
  int m, lda, n, ldu, ldvt, info, lwork;
  m = lda = ldu = A.Size(1);
  n = A.Size(2);
  ldvt = std::min(m, n);
  U.Resize(m, ldvt);
  Matd work(1), A1(A);
  Vt.Resize(ldvt, n);
  s.Resize(ldvt);
  char jobu = 'S', jobvt = 'S';
  lwork = -1; // workspace query
  dgesvd_
    (&jobu, &jobvt, &m, &n, A1.GetPtr(), &lda, s.GetPtr(), U.GetPtr(), &ldu,
     Vt.GetPtr(), &ldvt, work.GetPtr(), &lwork, &info);
  lwork = (int)work(1);
  work.Resize(lwork);
  dgesvd_
    (&jobu, &jobvt, &m, &n, A1.GetPtr(), &lda, s.GetPtr(), U.GetPtr(), &ldu,
     Vt.GetPtr(), &ldvt, work.GetPtr(), &lwork, &info);
}
#endif

}
