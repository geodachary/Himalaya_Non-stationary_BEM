#ifndef INCLUDE_MATRIX
#define INCLUDE_MATRIX

#include <iostream>
#include <stdlib.h>
#include <string.h>
//#include <assert.h>

template<typename T,int idxStart=1> class Matrix {
private:
  T* _data;
  int _m,_n;
  bool _amResponsible; // to delete _data

public:
  Matrix()               : _m(0),_n(0),_amResponsible(false) { Resize(0,0); }
  explicit Matrix(int m) : _m(0),_n(0),_amResponsible(false) { Resize(m,1); }
  Matrix(int m,int n)    : _m(0),_n(0),_amResponsible(false) { Resize(m,n); }
  ~Matrix()           { if(_amResponsible && _data) free(_data); }

  Matrix(const Matrix<T>& A) throw (OutOfMemoryException)
    : _m(0),_n(0),_amResponsible(false)
  {
    Resize(A._m,A._n);
    memcpy(_data,A._data,_m*_n*sizeof(T));
  }

  template<typename T1> Matrix(const Matrix<T1>& A)
    throw (OutOfMemoryException)
    : _m(0),_n(0),_amResponsible(false)
  {
    Resize(A.Size(1),A.Size(2));
    int sz = Size();
    for(int i=0; i<sz; ++i) _data[i] = (T)A(i+1);
  }

  Matrix<T>& operator=(const Matrix<T>& A) throw (OutOfMemoryException)
  {
    if(this != &A) {
      Resize(A._m,A._n);
      memcpy(_data,A._data,_m*_n*sizeof(T)); }
    return *this;
  }
  
  template<typename T1> Matrix<T>& operator=(const Matrix<T1>& A)
    throw (OutOfMemoryException)
  {
    if(this != (const Matrix<T>*) &A) {
      Resize(A.Size(1),A.Size(2));
      int sz = Size();
      for(int i=0; i<sz; ++i) _data[i] = (T)A(i+1); }
    return *this;
  }

  // Get size of row (dim=1), col (dim=2), or row*col (dim=0)
  int Size(int dim=0) const
  {
    switch(dim) {
    case 1: return _m;
    case 2: return _n;
    default: return _m*_n; }
  }

  Matrix<T>& Resize(int m) throw (OutOfMemoryException) { return Resize(m,1); }
  Matrix<T>& Reshape(int m) throw (OutOfMemoryException) { return Reshape(m,1); }

  // Resize matrix. Memory now handled by Matrix.
  Matrix<T>& Resize(int m,int n) throw (OutOfMemoryException)
  {
    if(_m*_n != m*n) {
      if(_amResponsible) {
	if(_data) free(_data);
	if(m==0 || n==0) _amResponsible = false;
	else _data = (T*)malloc(m*n*sizeof(T));
	if(!_data) throw OutOfMemoryException(); }
      else if(m>0 && n>0) {
	_data = (T*)malloc(m*n*sizeof(T));
	if(!_data) throw OutOfMemoryException();
	_amResponsible = true; }}
    _m = m; _n = n;
    return *this;
  }
  
  // Same as Resize, but copies old data over. The extra is set to zero.
  Matrix<T>& Reshape(int m,int n) throw (OutOfMemoryException)
  {
    if(_m*_n != m*n) {
      if(_m == 0 || _n == 0) {
	if(_data) free(_data);
	_amResponsible = false; }
      else {
	if(_amResponsible) {
	  _data = (T*)realloc(_data,m*n*sizeof(T));
	  if(!_data) throw OutOfMemoryException(); }
	else if(m>0 && n>0) {
	  T* data = (T*)malloc(m*n*sizeof(T));
	  if(!data) throw OutOfMemoryException();
	  if(_m>0 && _n>0) memcpy(data,_data,_m*_n*sizeof(T));
	  _data = data; }
	// Zero the remaining memory
	if(m*n>_m*_n) memset(_data+_m*_n,0,(m*n-_m*_n)*sizeof(T)); }}
    _m = m; _n = n;
    return *this;
  }
  
  Matrix<T>& SetPtr(int m,T* data) { return SetPtr(m,1,data); }
  
  Matrix<T>& SetPtr(int m,int n,T* data)
  {
    if(_amResponsible) {
      if(_data) free(_data);
      _amResponsible = false; }
    _data = data;
    _m = m; _n = n;
    return *this;
  }

  T* GetPtr() { return _data; }
  const T* GetPtr() const { return _data; }

  T& operator()(int i)
  {
    //assert(i >= idxStart && i < _m*_n+idxStart);
    return _data[i-idxStart];
  }

  const T& operator()(int i) const
  {
    //assert(i >= idxStart && i < _m*_n+idxStart);
    return _data[i-idxStart];
  }

  T& operator()(int i,int j)
  {
    //assert(i>=idxStart && j>=idxStart && i<_m+idxStart && j<_n+idxStart);
    return _data[(j-idxStart)*_m+i-idxStart];
  }
  
  const T& operator()(int i,int j) const
  {
    //assert(i>=idxStart && j>=idxStart && i<_m+idxStart && j<_n+idxStart);
    return _data[(j-idxStart)*_m+i-idxStart];
  }

  Matrix<T>& Bzero() { memset(_data,0,_m*_n*sizeof(T)); return *this; }

  Matrix<T>& Set(T v)
  {
    for(int i=0; i<_m*_n; ++i) _data[i] = v;
    return *this;
  }

  Matrix<T>& Zero() { return Set((T)0); }

  Matrix<T>& operator*=(T s)
  {
    for(int i=0; i<_m*_n; ++i) _data[i] *= s;
    return *this;
  }

  Matrix<T>& operator+=(T s)
  {
    for(int i=0; i<_m*_n; ++i) _data[i] += s;
    return *this;
  }
};

template<typename T1,typename T2>
Matrix<T2>& Transpose(const Matrix<T1>& A,Matrix<T2>& At)
{
  At.Resize(A.Size(2),A.Size(1));
  for(int i = 1; i <= A.Size(1); ++i)
    for(int j = 1; j <= A.Size(2); ++j)
      At(j,i) = A(i,j);
  return At;
}

template<typename T>
void mvp(const T* A, int m, int n, const T* x, T* y)
{
  int os = 0;
  memset(y, 0, m*sizeof(double));
  for (int j = 0; j < n; j++) {
    double xj = x[j];
    for (int i = 0; i < m; i++)
      y[i] += A[os + i]*xj;
    os += m;
  }
}

template<typename T>
Matrix<T>& mvp(const Matrix<T>& A, const Matrix<T>& x, Matrix<T>& y)
{
  int m = A.Size(1), n = A.size(2);
  y.Resize(m);
  mvp(A.GetPtr(), m, n, x.GetPtr(), y.GetPtr());
  return y;
}

template<typename T>
Matrix<T>& mmp(const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C)
{
  int m = A.Size(1), n = A.Size(2), Bm = B.Size(1), Bn = B.Size(2);
  int Bos = 0, Cos = 0;
  C.Resize(m, Bn);
  for (int k = 0; k < Bn; k++) {
    mvp(A.GetPtr(), m, n, B.GetPtr() + Bos, C.GetPtr() + Cos);
    Bos += Bm;
    Cos += m;
  }
  return C;
}

template<typename T>
std::ostream& operator<<(std::ostream& os,const Matrix<T>& A)
{
  int m,n;
  m = A.Size(1);
  n = A.Size(2);
  os<<"[";
  if(n == 1)
    for(int r = 1; r <= m; ++r) os<<A(r)<<" ";
  else
    for(int r = 1; r <= m; ++r) {
      if(r>1) os<<"\n";
      for(int c = 1; c < n; ++c) os<<A(r,c)<<" ";
      os<<A(r,n)<<";"; }
  os<<"]";
  return os;
}

#endif
