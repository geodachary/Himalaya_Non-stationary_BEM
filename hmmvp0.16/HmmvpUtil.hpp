#ifndef INCLUDE_HMMVPUTIL
#define INCLUDE_HMMVPUTIL

#include <vector>
#include <string>
#include <sstream>

class Exception {
public:
  Exception(const std::string& msg) { _msg = msg; }
  virtual ~Exception() {}
  const std::string& GetMsg() const { return _msg; }
private:
  std::string _msg;
};

class FileException : public Exception {
public:
  FileException(const std::string& msg) : Exception(msg) {}
};

class OutOfMemoryException : public Exception {
public:
  OutOfMemoryException(const std::string& msg = "oom") : Exception(msg) {}
};

#include "Matrix.hpp"
#include "Defs.hpp"

namespace Hm {

template<typename T> void read(T* v, size_t sz, FILE* fid)
  throw (FileException)
{
  size_t nread;
  if (!fid) throw FileException("read: fid is null");
  if ((nread = fread(v, sizeof(T), sz, fid)) != sz) {
    std::stringstream s;
    s << "read: nread = " << nread << " sz = " << sz;
    throw FileException(s.str());
  }
}

template<typename T> void write(T* v, size_t sz, FILE* fid)
  throw (FileException)
{
  size_t nwrite;
  if (!fid) throw FileException("write: fid is null");
  if ((nwrite = fwrite(v, sizeof(T), sz, fid)) != sz) {
    std::stringstream s;
    s << "write: nwrite = " << nwrite << " sz = " << sz;
    throw FileException(s.str());
  }
}

template<typename T1, typename T2> void MatrixToVector
(const Matrix<T1>& A, std::vector<T2>& v, T2 a = 0)
{
  int n = A.Size();
  v.reserve(n);
  for (int i = 1; i <= n; i++)
    v.push_back(static_cast<T1>(A(i)) + a);
}

template<typename T1, typename T2> void VectorToMatrix
(const std::vector<T1>& v, Matrix<T2>& A, T2 a = 0)
  throw (OutOfMemoryException)
{
  int n = v.size();
  A.Resize(n);
  for (int i = 0; i < n; i++)
    A(i+1) = static_cast<T1>(v[i] + a);
}

#ifdef HMMVP_LAPACK
class LapackException : public Exception {
public:
  LapackException(const std::string& msg) : Exception(msg) {}
};

// Corresponds to Matlab's [U S V] = svd(A,'econ').
void LapackSvd(const Matd& A, Matd& U, Matd& s, Matd& Vt)
  throw (LapackException);
#endif

// LAPACK/BLAS interface
typedef long long blas_int;

extern "C" void sgemm_
 (char*, char*, blas_int*, blas_int*, blas_int*, float*, const float*,
  blas_int*, const float*, blas_int*, float*, float*, blas_int*);
extern "C" void sgemv_
 (char*, blas_int*, blas_int*, float*, const float*, blas_int*, const float*,
  blas_int*, float*, float*, blas_int*);
extern "C" float sdot_
 (blas_int*, const float*, blas_int*, const float*, blas_int*);
extern "C" float saxpy_
 (blas_int*, float*, const float*, blas_int*, float*, blas_int*);
extern "C" void dgemm_
 (char*, char*, blas_int*, blas_int*, blas_int*, double*, const double*,
  blas_int*, const double*, blas_int*, double*, double*, blas_int*);
extern "C" void dgemv_
 (char*, blas_int*, blas_int*, double*, const double*, blas_int*, const double*,
  blas_int*, double*, double*, blas_int*);
extern "C" double ddot_
 (blas_int*, const double*, blas_int*, const double*, blas_int*);
extern "C" double daxpy_
 (blas_int*, double*, const double*, blas_int*, double*, blas_int*);

template<typename T> void gemm
 (char transa, char transb, blas_int m, blas_int n, blas_int k,
  T alpha, const T* A, blas_int lda, const T* B, blas_int ldb, T beta,
  T* C, blas_int ldc);
template<typename T> void gemv
 (char trans, blas_int m, blas_int n, T alpha, const T* a, blas_int lda,
  const T* x, blas_int incx, T beta, T* y, blas_int incy);
template<typename T> T dot
 (blas_int n, const T* x, blas_int incx, const T* y, blas_int incy);
template<typename T> T axpy
 (blas_int n, T a, const T* x, blas_int incx, T* y, blas_int incy);

template<> inline void gemm<double>
 (char transa, char transb, blas_int m, blas_int n, blas_int k,
  double alpha, const double* A, blas_int lda,
  const double* B, blas_int ldb, double beta,
  double* C, blas_int ldc)
{
  dgemm_(&transa, &transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta,
	 C, &ldc);
}

template<> inline void gemm<float>
 (char transa, char transb, blas_int m, blas_int n, blas_int k,
  float alpha, const float* A, blas_int lda,
  const float* B, blas_int ldb, float beta,
  float* C, blas_int ldc)
{
  sgemm_(&transa, &transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta,
	 C, &ldc);
}

template<> inline void gemv<double>
(char trans, blas_int m, blas_int n, double alpha,
 const double* a, blas_int lda, const double* x, blas_int incx, double beta,
 double* y, blas_int incy)
{
  dgemv_(&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

template<> inline void gemv<float>
(char trans, blas_int m, blas_int n, float alpha,
 const float* a, blas_int lda, const float* x, blas_int incx, float beta,
 float* y, blas_int incy)
{
  sgemv_(&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

template<> inline double dot<double>
(blas_int n, const double* x, blas_int incx, const double* y, blas_int incy)
{
  return ddot_(&n, x, &incx, y, &incy);
}

template<> inline float dot<float>
(blas_int n, const float* x, blas_int incx, const float* y, blas_int incy)
{
  return sdot_(&n, x, &incx, y, &incy);
}

template<> inline double axpy<double>
(blas_int n, double a, const double* x, blas_int incx, double* y, blas_int incy)
{
  return daxpy_(&n, &a, x, &incx, y, &incy);
}

template<> inline float axpy<float>
(blas_int n, float a, const float* x, blas_int incx, float* y, blas_int incy)
{
  return saxpy_(&n, &a, x, &incx, y, &incy);
}

}

#endif
