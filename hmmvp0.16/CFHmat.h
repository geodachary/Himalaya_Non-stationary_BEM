/* Fortran interface to some of Hmat.hpp.  */

#ifndef INCLUDE_CFHMAT
#define INCLUDE_CFHMAT

typedef float Real;
typedef long long int fint;

extern "C" {
  /* For gfortran, we use the following Fortran -> C conversion rules:
     1. All lowercase function name.
     2. Underscore at end of function name.
     3. If a char* is an argument, then there is an extra argument at the end
        that gives the string length.
     4. All arguments are pointers.  */
  int cfhmat_init_(const char* hmat_fn, fint* ncol, fint* max_nthreads,
		   /* Do not pass this next argument in your Fortran file; it
		      is automatically added by gfortran.  */
		   fint len_hmat_fn);
  int cfhmat_mvp_(const Real* x, Real* y, fint* ncol);
  void cfhmat_cleanup_();
}

#endif
