#include <string>
#include "CFHmat.h"
#include "Hmat.hpp"
using namespace Hm;

// Since C and Fortran can't handle C++ objects, we need to hide a global one in
// this file.
static Hmat* hm = 0;

int cfhmat_init_(const char* hmat_fn, fint* ncol, fint* max_nthreads,
		 fint len_hmat_fn)
{
  cfhmat_cleanup_();
  try {
    string my_hmat_fn(hmat_fn);
    my_hmat_fn.resize(len_hmat_fn);
    hm = NewHmat(my_hmat_fn, *ncol, *max_nthreads);
  } catch (const Exception& e) {
    printf("cfhmat_init_: %s\n", e.GetMsg().c_str());
    return -1;
  }
  return 0;
}

int cfhmat_mvp_(const Real* x, Real* y, fint* ncol)
{
  if (!hm) return -1;
  try {
    hm->Mvp(x, y, *ncol);
  } catch (const Exception& e) {
    printf("cfhmat_mvp_: %s\n", e.GetMsg().c_str());
    return -1;
  }
  return 0;
}

void cfhmat_cleanup_()
{
  if (hm) delete hm;
  hm = 0;
}
