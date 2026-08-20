! An example of how to use Hmat.hpp in Fortran code by using the CFHmat.h
! interface.

! On my machine, I build this file as follows:
!   g++ -c -DHMMVP_OMP -fopenmp CFHmat.cpp ;
!   gfortran -O -fdefault-integer-8 ex_mvp.f90 CFHmat.o \
!     ~/dl/lapack-3.2.1/libblas.a -fopenmp -lstdc++ -o ex_mvp ;

program main
  real*4 x(235420), y(117710)
  integer ncol, nthreads

  ncol = 1
  nthreads = 4
  x(:) = 1.0e0
  call CFHmat_Init('/scratch/ambrad/sd4v0.8/bem_id1tol-5.7nx298nz395Bfro4.1e+10sym_comp11.dat', ncol, nthreads)
  do i = 1, 20
     call CFHmat_Mvp(x, y, ncol)
  end do
  write(*,*) y(1), y(100), y(117710)
  call CFHmat_Cleanup()
end program main
