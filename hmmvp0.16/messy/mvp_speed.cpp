// g++ -O3 -DHMMVP_OMP mvp_speed.cpp ~/dl/lapack-3.2.1/libblas.a -lgfortran -fopenmp -o mvp_speed

#include <iostream>
#include <math.h>
#include "Hmat.hpp"
#include <sys/time.h>
using namespace std;

static double difftime(struct timeval& t1, struct timeval& t2)
{
  static const double us = 1.0e6;
  return (t2.tv_sec*us + t2.tv_usec - t1.tv_sec*us - t1.tv_usec)/us;
}

double relerr(const Matd& y1, const Matd& y2)
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
  const int nrhs = 2, max_nthreads = 16, nrepeat = 20;

  Hm::Hmat* hm;
  try {
    hm = Hm::NewHmat(fn, nrhs, max_nthreads);
  } catch (const Exception& e) {
    cout << e.GetMsg() << endl;
    exit(-1);
  }

  printf("Bfro = %e\n", sqrt(hm->NormFrobenius2()));

  const Hm::Blint n = hm->GetN(), m = hm->GetM();
  Matd x(n,nrhs);
  for (uint j = 1; j <= nrhs; j++)
    for (uint i = 1; i <= n; i++)
      x(i,j) = drand48() - 0.5;

  Matd y(m,nrhs), y1;
  double et1, re;

  if (true) {
    printf("perms off\n");
    hm->TurnOffPermute();
  }

  // New implementation.
  int nthrs[8] = {1, 4, 8, 10, 12, 14, 15, 16};
  for (int i = 0; i < 8; i++) {
    int nthr = nthrs[i];
    nthr = hm->SetThreads(nthr);
    hm->Mvp(x.GetPtr(), y.GetPtr(), nrhs);
    gettimeofday(&t1, 0);
    for (int j = 0; j < nrepeat; j++)
      hm->Mvp(x.GetPtr(), y.GetPtr(), nrhs);
    gettimeofday(&t2, 0);
    double et = difftime(t1, t2);
    if (i == 0) {
      et1 = et;
      y1 = y;
    }
    printf("%d %f %f %e\n", nthr, et/nrepeat, et1/et, relerr(y1, y));
  }

  delete hm;
}
