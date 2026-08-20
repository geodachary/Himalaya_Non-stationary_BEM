// mpic++ -O3 mpi_mvp_speed.cpp ~/dl/lapack-3.2.1/libblas.a -lgfortran -o mpi_mvp_speed

// Hmmm: MPI_Init, etc., doesn't seem to play nicely with my OpenMP code:
// mpic++ -O -DHMMVP_OMP mpi_mvp_speed.cpp ~/dl/lapack-3.2.1/libblas.a -lgfortran -fopenmp -o m,.acdflmnorzpi_mvp_speed

// mpirun -np 16 mpi_mvp_speed

#include <iostream>
#include <math.h>
#include "MpiHmat.hpp"
#include <sys/time.h>
using namespace std;

typedef Matrix<float> Matf;

static double difftime(struct timeval& t1, struct timeval& t2)
{
  static const double us = 1.0e6;
  return (t2.tv_sec*us + t2.tv_usec - t1.tv_sec*us - t1.tv_usec)/us;
}

double relerr(const Matf& y1, const Matf& y2)
{
  assert(y1.Size() == y2.Size());
  double num = 0.0, den = 0.0;
  for (uint i = 1; i <= y1.Size(); i++) {
    double d = y1(i) - y2(i);
    num += d*d;
    den += y1(i)*y1(i);
  }
  return sqrt(num/den);
}

int main(int argc, char** argv)
{
  const char* fn =
    "/scratch/ambrad/evolve/bem_id1tol-2.0nx1337nz2478sym_comp11.dat";
  struct timeval t1, t2;
  const int nrhs = 2, nrepeat = 20;

  MPI_Init(&argc, &argv);
  int nproc, pid;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  //printf("nproc = %d  pid = %d\n", nproc, pid);

  Matf x, y, y1;
  Hm::Blint m, n;
  bool do_perms = false;

  MPI_Barrier(MPI_COMM_WORLD);
  Hm::MpiHmat<float>* hm;
  try {
    hm = new Hm::MpiHmat<float>(fn, nrhs);
  } catch (const Exception& e) {
    cout << e.GetMsg() << endl;
    exit(-1);
  }
  m = hm->GetM();
  n = hm->GetN();
  y.Resize(m,nrhs);
  if (pid == 0 && x.Size() == 0) {
    x.Resize(n,nrhs);
    for (uint j = 1; j <= nrhs; j++)
      for (uint i = 1; i <= n; i++)
	x(i,j) = drand48() - 0.5;
  }
  if (pid != 0) x.Resize(n,nrhs);
  MPI_Bcast(x.GetPtr(), nrhs*hm->GetN(), MPI_FLOAT, 0, MPI_COMM_WORLD);    

  if (!do_perms) hm->TurnOffPermute();
  hm->Mvp(x.GetPtr(), y.GetPtr(), nrhs);
  MPI_Barrier(MPI_COMM_WORLD);
  if (pid == 0) gettimeofday(&t1, 0);
  for (int j = 0; j < nrepeat; j++)
    hm->Mvp(x.GetPtr(), y.GetPtr(), nrhs);
  MPI_Barrier(MPI_COMM_WORLD);
  if (pid == 0) {
    gettimeofday(&t2, 0);
    double et = difftime(t1, t2)/nrepeat;
    printf("mpi(%d): et = %f\n", nproc, et);
  }

  do_perms = !do_perms;
  if (do_perms) hm->TurnOnPermute();
  else hm->TurnOffPermute();
  hm->Mvp(x.GetPtr(), y.GetPtr(), nrhs);
  MPI_Barrier(MPI_COMM_WORLD);
  if (pid == 0) gettimeofday(&t1, 0);
  for (int j = 0; j < nrepeat; j++)
    hm->Mvp(x.GetPtr(), y.GetPtr(), nrhs);
  MPI_Barrier(MPI_COMM_WORLD);
  if (pid == 0) {
    gettimeofday(&t2, 0);
    double et = difftime(t1, t2)/nrepeat;
    printf("mpi(%d): et = %f with do_perms = %d\n", nproc, et, do_perms);
  }

  delete hm;

  MPI_Finalize();
}
