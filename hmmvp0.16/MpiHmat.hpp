// Andrew M. Bradley ambrad@cs.stanford.edu

/* Routines to parallelize Hmat::Mvp(). We assume that MPI_Init has already been
   called.
     Basic usage:

       MPI_Init(&argc, &argv);
       Hm::MpiHmat<float>* hm;
       try {
         hm = new Hm::MpiHmat<float>(hmat_fn, nrhs);
       } catch (const Exception& e) {
         cout << e.GetMsg() << endl;
         exit(-1);
       }
       hm->Mvp(x, y, nrhs);
       delete hm;
       MPI_Finalize();

     Compile with mpic++. Run with a command of the sort (in Unix):
       mpirun -np 4 ./a.out
*/


/* TODO: In Mvp(x, y, ncol), one could imagine having a class for x and y that
   puts different pieces on each proc automatically. Then all of x and y
   wouldn't end up on every proc.  */

#ifndef INCLUDE_MPI_HMAT
#define INCLUDE_MPI_HMAT

#include <mpi.h>
#include "Hmat.hpp"

namespace Hm {
  using namespace std;

  // Corresponds to Hm::Blint. MPI_BLINT is not a type, but rather a macro, so a
  // typedef doesn't work here.
#define MPI_BLINT MPI_LONG
  
  template<typename Real> class MpiHmat {
  public:
    MpiHmat(const string& filename, uint ncol = 1) throw (FileException);
    ~MpiHmat();

    // Number of rows.
    Blint GetM() const { return hm->GetM(); }
    // Number of columns.
    Blint GetN() const { return hm->GetN(); }
    // Number of nonzeros in the H-matrix approximation.
    Blint GetNnz() const { return nnz; }

    // y = B*x, where x has ncol columns. y is allocated by the caller.
    void Mvp(const Real *x, Real *y, uint ncol) const
      throw (Exception);

    void TurnOnPermute()  { do_perms = true;  }
    void TurnOffPermute() { do_perms = false; }

  private:
    Hmat* hm;
    mutable vector<Real> workr;
    bool am_root;
    MPI_Datatype mpi_datatype;
    Blint nnz;
    bool do_perms;
  };

}

// As templates can't be compiled separately:
#include "MpiHmat.cpp"

#endif
