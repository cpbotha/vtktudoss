// -*- C++ -*-

/*! 
  \file ext/vector.h
  \brief Functions for vector.
*/

#if !defined(__ext_vector_h__)
#define __ext_vector_h__

#include <vector>
#include <iterator>
#include <numeric>
#include <iostream>

#include <cmath>
#include <cassert>

namespace std {

//----------------------------------------------------------------------------
//! \defgroup VectorAssignmentScalar Vector Assignment Operators with a Scalar Operand.
//@{

//! Vector-scalar addition.
template<typename _T>
inline
vector<_T>&
operator+=(vector<_T>& x, const _T& value) {
  for (typename vector<_T>::iterator i = x.begin(); i != x.end(); ++i) {
    *i += value;
  }
  return x;
}

//! Vector-scalar subtraction.
template<typename _T>
inline
vector<_T>&
operator-=(vector<_T>& x, const _T& value) {
  for (typename vector<_T>::iterator i = x.begin(); i != x.end(); ++i) {
    *i -= value;
  }
  return x;
}

//! Vector-scalar multiplication.
template<typename _T>
inline
vector<_T>&
operator*=(vector<_T>& x, const _T& value) {
  for (typename vector<_T>::iterator i = x.begin(); i != x.end(); ++i) {
    *i *= value;
  }
  return x;
}

//! Vector-scalar division.
template<typename _T>
inline
vector<_T>&
operator/=(vector<_T>& x, const _T& value) {
#ifdef DEBUG_stlib
  assert(value != 0);
#endif
  for (typename vector<_T>::iterator i = x.begin(); i != x.end(); ++i) {
    *i /= value;
  }
  return x;
}

//! Vector-scalar modulus.
template<typename _T>
inline
vector<_T>&
operator%=(vector<_T>& x, const _T& value) {
#ifdef DEBUG_stlib
  assert(value != 0);
#endif
  for (typename vector<_T>::iterator i = x.begin(); i != x.end(); ++i) {
    *i %= value;
  }
  return x;
}

//! Left shift.
template<typename _T>
inline
vector<_T>&
operator<<=(vector<_T>& x, const int offset) {
  for (typename vector<_T>::iterator i = x.begin(); i != x.end(); ++i) {
    *i <<= offset;
  }
  return x;
}

//! Right shift.
template<typename _T>
inline
vector<_T>&
operator>>=(vector<_T>& x, const int offset) {
  for (typename vector<_T>::iterator i = x.begin(); i != x.end(); ++i) {
    *i >>= offset;
  }
  return x;
}

//@}
//----------------------------------------------------------------------------
//! \defgroup VectorAssignmentVector Vector Assignment Operators with a Vector Operand.
//@{


//! Vector-vector addition.
template<typename _T1, typename _T2>
inline
vector<_T1>&
operator+=(vector<_T1>& x, const vector<_T2>& y) {
#ifdef DEBUG_stlib
  assert(x.size() == y.size());
#endif
  for (size_t n = 0; n != x.size(); ++n) {
    x[n] += y[n];
  }
  return x;
}

//! Vector-vector subtraction.
template<typename _T1, typename _T2>
inline
vector<_T1>&
operator-=(vector<_T1>& x, const vector<_T2>& y) {
#ifdef DEBUG_stlib
  assert(x.size() == y.size());
#endif
  for (size_t n = 0; n != x.size(); ++n) {
    x[n] -= y[n];
  }
  return x;
}

//! Vector-vector multiplication.
template<typename _T1, typename _T2>
inline
vector<_T1>&
operator*=(vector<_T1>& x, const vector<_T2>& y) {
#ifdef DEBUG_stlib
  assert(x.size() == y.size());
#endif
  for (size_t n = 0; n != x.size(); ++n) {
    x[n] *= y[n];
  }
  return x;
}

//! Vector-vector division.
template<typename _T1, typename _T2>
inline
vector<_T1>&
operator/=(vector<_T1>& x, const vector<_T2>& y) {
#ifdef DEBUG_stlib
  assert(x.size() == y.size());
#endif
  for (size_t n = 0; n != x.size(); ++n) {
#ifdef DEBUG_stlib
    assert(y[n] != 0);
#endif
    x[n] /= y[n];
  }
  return x;
}

//! Vector-vector modulus.
template<typename _T1, typename _T2>
inline
vector<_T1>&
operator%=(vector<_T1>& x, const vector<_T2>& y) {
#ifdef DEBUG_stlib
  assert(x.size() == y.size());
#endif
  for (size_t n = 0; n != x.size(); ++n) {
#ifdef DEBUG_stlib
    assert(y[n] != 0);
#endif
    x[n] %= y[n];
  }
  return x;
}

//! Vector-vector left shift.
template<typename _T1, typename _T2>
inline
vector<_T1>&
operator<<=(vector<_T1>& x, const vector<_T2>& y) {
#ifdef DEBUG_stlib
  assert(x.size() == y.size());
#endif
  for (size_t n = 0; n != x.size(); ++n) {
    x[n] <<= y[n];
  }
  return x;
}

//! Vector-vector right shift.
template<typename _T1, typename _T2>
inline
vector<_T1>&
operator>>=(vector<_T1>& x, const vector<_T2>& y) {
#ifdef DEBUG_stlib
  assert(x.size() == y.size());
#endif
  for (size_t n = 0; n != x.size(); ++n) {
    x[n] >>= y[n];
  }
  return x;
}

//@}
//----------------------------------------------------------------------------
//! \defgroup VectorFile Vector File I/O
//@{

//! Write the size and then the space-separated elements.
/*! 
  Format:
  x.size()
  x[0] x[1] x[2] ...
*/
template<typename _T>
inline
ostream&
operator<<(ostream& out, const vector<_T>& x) {
  out << x.size() << '\n';
  copy(x.begin(), x.end(), ostream_iterator<_T>(out, " "));
  return out;
}
    
//! Read the size and then the elements.
/*!
  The vector will be resized.
*/
template<typename _T>
inline
istream&
operator>>(istream& in, vector<_T>& x) {
  size_t size;
  in >> size;
  x.resize(size);
  for (size_t n = 0; n != x.size(); ++n) {
    in >> x[n];
  }
  return in;
}

//@}
//----------------------------------------------------------------------------
//! \defgroup VectorMathematical Vector Mathematical Functions
//@{

//! Return the sum of the components.
template<typename _T>
inline
_T
sum(const vector<_T>& x) {
  return accumulate(x.begin(), x.end(), _T(0));
}
    
//! Return the product of the components.
template<typename _T>
inline
_T
product(const vector<_T>& x) {
  return accumulate(x.begin(), x.end(), _T(1), multiplies<_T>());
}

//! Return the minimum component.  Use < for comparison.
template<typename _T>
inline
_T
min(const vector<_T>& x) {
#ifdef DEBUG_stlib
  assert(x.size() != 0);
#endif
  return *min_element(x.begin(), x.end());
}
    
//! Return the maximum component.  Use > for comparison.
template<typename _T>
inline
_T
max(const vector<_T>& x) {
#ifdef DEBUG_stlib
  assert(x.size() != 0);
#endif
  return *max_element(x.begin(), x.end());
}

//! Return the dot product of the two vectors.
template<typename _T1, typename _T2>
inline
_T1
dot(const vector<_T1>& x, const vector<_T2>& y) {
  return inner_product(x.begin(), x.end(), y.begin(), _T1(0));
}

//@}
//----------------------------------------------------------------------------
//! \defgroup VectorApply Apply the Standard Math Functions.
//@{

//! Apply the absolute value (\f$|x|\f$) to each element.
template<typename _T>
inline
void
applyAbs(vector<_T>* x) {
  for (typename vector<_T>::iterator i = x->begin(); i != x->end(); ++i) {
    *i = abs(*i);
  }
}

//! Apply the inverse cosine (\f$ \cos^{-1}(x) \f$) to each element.
template<typename _T>
inline
void
applyAcos(vector<_T>* x) {
  for (typename vector<_T>::iterator i = x->begin(); i != x->end(); ++i) {
    *i = acos(*i);
  }
}

//! Apply the inverse sine (\f$ \sin^{-1}(x) \f$) to each element.
template<typename _T>
inline
void
applyAsin(vector<_T>* x) {
  for (typename vector<_T>::iterator i = x->begin(); i != x->end(); ++i) {
    *i = asin(*i);
  }
}

//! Apply the inverse tangent (\f$ \tan^{-1}(x) \f$) to each element.
template<typename _T>
inline
void
applyAtan(vector<_T>* x) {
  for (typename vector<_T>::iterator i = x->begin(); i != x->end(); ++i) {
    *i = atan(*i);
  }
}

//! Apply the ceiling function (\f$ \lceil x \rceil \f$) to each element.
template<typename _T>
inline
void
applyCeil(vector<_T>* x) {
  for (typename vector<_T>::iterator i = x->begin(); i != x->end(); ++i) {
    *i = ceil(*i);
  }
}

//! Apply the cosine (\f$ \cos(x) \f$) to each element.
template<typename _T>
inline
void
applyCos(vector<_T>* x) {
  for (typename vector<_T>::iterator i = x->begin(); i != x->end(); ++i) {
    *i = cos(*i);
  }
}

//! Apply the hyperbolic cosine (\f$ \cosh(x) \f$) to each element.
template<typename _T>
inline
void
applyCosh(vector<_T>* x) {
  for (typename vector<_T>::iterator i = x->begin(); i != x->end(); ++i) {
    *i = cosh(*i);
  }
}

//! Apply the exponential function (\f$ \mathrm{e}^x \f$) to each element.
template<typename _T>
inline
void
applyExp(vector<_T>* x) {
  for (typename vector<_T>::iterator i = x->begin(); i != x->end(); ++i) {
    *i = exp(*i);
  }
}

//! Apply the floor function (\f$ \lfloor x \rfloor \f$) to each element.
template<typename _T>
inline
void
applyFloor(vector<_T>* x) {
  for (typename vector<_T>::iterator i = x->begin(); i != x->end(); ++i) {
    *i = floor(*i);
  }
}

//! Apply the natural logarithm (\f$ \ln(x) \f$) to each element.
template<typename _T>
inline
void
applyLog(vector<_T>* x) {
  for (typename vector<_T>::iterator i = x->begin(); i != x->end(); ++i) {
    *i = log(*i);
  }
}

//! Apply the logarithm base 10 (\f$ \log_{10}(x) \f$) to each element.
template<typename _T>
inline
void
applyLog10(vector<_T>* x) {
  for (typename vector<_T>::iterator i = x->begin(); i != x->end(); ++i) {
    *i = log10(*i);
  }
}

//! Apply the sine (\f$ \sin(x) \f$) to each element.
template<typename _T>
inline
void
applySin(vector<_T>* x) {
  for (typename vector<_T>::iterator i = x->begin(); i != x->end(); ++i) {
    *i = sin(*i);
  }
}

//! Apply the hyperbolic sine (\f$ \sinh(x) \f$) to each element.
template<typename _T>
inline
void
applySinh(vector<_T>* x) {
  for (typename vector<_T>::iterator i = x->begin(); i != x->end(); ++i) {
    *i = sinh(*i);
  }
}

//! Apply the square root (\f$ \sqrt{x} \f$) to each element.
template<typename _T>
inline
void
applySqrt(vector<_T>* x) {
  for (typename vector<_T>::iterator i = x->begin(); i != x->end(); ++i) {
    *i = sqrt(*i);
  }
}

//! Apply the tangent (\f$ \tan(x) \f$) to each element.
template<typename _T>
inline
void
applyTan(vector<_T>* x) {
  for (typename vector<_T>::iterator i = x->begin(); i != x->end(); ++i) {
    *i = tan(*i);
  }
}

//! Apply the hyperbolic tangent (\f$ \tanh(x) \f$) to each element.
template<typename _T>
inline
void
applyTanh(vector<_T>* x) {
  for (typename vector<_T>::iterator i = x->begin(); i != x->end(); ++i) {
    *i = tanh(*i);
  }
}

//@}

}

#endif
