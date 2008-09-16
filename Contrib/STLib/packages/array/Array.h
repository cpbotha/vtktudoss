// -*- C++ -*-

/*! 
  \file array/Array.h
  \brief Fixed size array that derives from std::tr1::array.
*/

#if !defined(__array_Array_h__)
#define __array_Array_h__

#include "defs.h"

#include <boost/call_traits.hpp>
#include <boost/static_assert.hpp>

#include <tr1/array>

#include <numeric>
#include <iostream>

#include <cmath>
#include <cassert>

#if defined(DEBUG_array) && !defined(DEBUG_array_Array)
#define DEBUG_array_Array
#endif

BEGIN_NAMESPACE_ARRAY

//! Fixed size array that derives from std::tr1::array.
/*!
  <b>Constructors etc.</b>

  The default constructor makes an array with uninitialized data.
  \verbatim
  array::Array<int, 3> a; \endverbatim

  The copy constructor creates an array with the same element values as the 
  argument.
  \verbatim
  array::Array<int, 3> a;
  ...
  array::Array<int, 3> b(a);
  array::Array<int, 3> c = a; \endverbatim

  There is also a copy constructor for arrays with different value types.
  (However the size must be the same.)
  \verbatim
  array::Array<int, 3> a;
  ...
  array::Array<double, 3> b(a); \endverbatim

  The assignment operators copy the element values.
  \verbatim
  array::Array<int, 3> a;
  ...
  array::Array<int, 3> b;
  b = a;
  array::Array<double, 3> c;
  c = a; \endverbatim

  There are constructors that take the element values as arguments.
  \verbatim
  array::Array<int, 1> a(2);
  array::Array<int, 2> a(2, 3);
  array::Array<int, 3> a(2, 3, 5);
  array::Array<int, 4> a(2, 3, 5, 7); \endverbatim

  You can also construct an array an initialize each element to a specified 
  value.
  \verbatim
  array::Array<int, 3> a(2); \endverbatim

  There is a copy constructor and assignment operator for copying 
  a std::tr1::array .
  \verbatim
  std::tr1::array<int, 3> a;
  array::Array<int, 3> b(a);
  array::Array<double, 3> c;
  c = a; \endverbatim

  Finally, you can construct an array and copy the elements from a C array.
  \verbatim
  int a[3];
  array::Array<int, 3> b(a); \endverbatim

  <b>Member Functions</b>

  Array inherits the following member functions from std::tr1::array .
  - assign()
  - swap()
  - begin()
  - end()
  - rbegin()
  - rend()
  - size()
  - max_size()
  - empty()
  - at()
  - front()
  - back()
  - data()

  <b>Free Operators and Functions</b>

  Array supports the usual mathematical operations.
  These functions are grouped into the following categories.
  - \ref ArrayUnaryOperators
  - \ref ArrayAssignmentScalar
  - \ref ArrayAssignmentArray
  - \ref ArrayBinary
  - \ref ArrayComparison
  - \ref ArrayFile
  - \ref ArrayMathematical
  - \ref ArrayElement
  - \ref ArrayApply
  - \ref ArrayStandard
*/
template<typename _T, std::size_t _N>
class Array : public std::tr1::array<_T, _N> {
  //
  // Types.
  //
private:

  //! The base class.
  typedef std::tr1::array<_T, _N> Base;

public:

  //! The value type.
  typedef typename Base::value_type value_type;
  //! Reference to the value type.
  typedef typename Base::reference reference;
  //! Constant reference to the value type.
  typedef typename Base::const_reference const_reference;
  //! Iterator in the container.
  typedef typename Base::iterator iterator;
  //! Constant iterator in the container.
  typedef typename Base::const_iterator const_iterator;
  //! The size type.
  typedef typename Base::size_type size_type;
  //! The pointer difference type.
  typedef typename Base::difference_type difference_type;
  //! Reverse iterator.
  typedef typename Base::reverse_iterator reverse_iterator;
  //! Constant reverse iterator.
  typedef typename Base::const_reverse_iterator const_reverse_iterator;

  //! Type for efficiently passing the value type as a parameter.
  typedef typename boost::call_traits<value_type>::param_type param_type;

  //
  // Using.
  //
public:

  // CONTINUE: assign is broken in GCC 4.0.
  //using Base::assign;

  //! Assign to each element the specified value.
  void 
  assign(param_type x) {
    std::fill(begin(), end(), x);
  }

  using Base::swap;
  using Base::begin;
  using Base::end;
  using Base::rbegin;
  using Base::rend;
  using Base::size;
  using Base::max_size;
  using Base::empty;
  using Base::at;
  using Base::front;
  using Base::back;


  // CONTINUE: assign is broken in GCC 4.0.
  //using Base::data;

  //! Return a pointer to the data.
  _T* 
  data() {
    return begin();
  }

  //! Return a const pointer to the data.
  const _T* 
  data() const {
    return begin();
  }

  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  //@{
public:

  //! Default constructor. Uninitialized data.
  Array() {
  }

  //! Copy constructor for std::tr1::array.
  Array(const std::tr1::array<_T, _N>& x) :
    Base(x) {
  }

  //! Copy constructor for a different value type.
  template<typename _T2>
  Array(const Array<_T2, _N>& x) :
    Base() {
    *this = x;
  }

  //! Copy constructor.
  Array(const Array& other) :
    Base(other) {
  }

  //! Assignment operator for other array types.
  template<typename _T2>
  Array&
  operator=(const std::tr1::array<_T2, _N>& x) {
    std::copy(x.begin(), x.end(), begin());
    return *this;
  }

  //! Assignment operator.
  Array&
  operator=(const Array& other) {
    if (this != &other) {
      Base::operator=(other);
    }
    return *this;
  }

  //! Construct and initialize all of the elements.
  Array(param_type x) :
    Base() {
    std::fill(begin(), end(), x);
  }

  //! Construct and initialize the first two elements.
  Array(param_type x0, param_type x1) :
    Base() {
    (*this)[0] = x0;
    (*this)[1] = x1;
  }

  //! Construct and initialize the first three elements.
  Array(param_type x0, param_type x1, param_type x2) :
    Base() {
    (*this)[0] = x0;
    (*this)[1] = x1;
    (*this)[2] = x2;
  }

  //! Construct and initialize the first four elements.
  Array(param_type x0, param_type x1, param_type x2, param_type x3) :
    Base() {
    (*this)[0] = x0;
    (*this)[1] = x1;
    (*this)[2] = x2;
    (*this)[3] = x3;
  }

  //! Construct and initialize the elements from the array.
  Array(const_iterator i) :
    Base() {
    std::copy(i, i + size(), begin());
  }

  //! Destructor.
  ~Array() {
  }

  //@}
};

//----------------------------------------------------------------------------
//! \defgroup ArrayUnaryOperators Array Unary Operators
//@{

//! Unary positive operator.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
Array<_T, _N>&
operator+(Array<_T, _N>& x) {
  return x;
}
    
//! Unary negate operator.
/*! \relates Array */
template<typename _T, std::size_t _N>
Array<_T, _N>
operator-(const Array<_T, _N>& x) {
  Array<_T, _N> a(x);
  for (std::size_t n = 0; n != _N; ++n) {
    a[n] = - a[n];
  }
  return a;
}

//@}
//----------------------------------------------------------------------------
//! \defgroup ArrayAssignmentScalar Array Assignment Operators with a Scalar Operand.
//@{

//! Array-scalar addition.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
Array<_T, _N>&
operator+=(Array<_T, _N>& x,
	   typename boost::call_traits<_T>::param_type value) {
  for (std::size_t n = 0; n != _N; ++n) {
    x[n] += value;
  }
  return x;
}

//! Array-scalar subtraction.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
Array<_T, _N>&
operator-=(Array<_T, _N>& x, 
	   typename boost::call_traits<_T>::param_type value) {
  for (std::size_t n = 0; n != _N; ++n) {
    x[n] -= value;
  }
  return x;
}

//! Array-scalar multiplication.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
Array<_T, _N>&
operator*=(Array<_T, _N>& x, 
	   typename boost::call_traits<_T>::param_type value) {
  for (std::size_t n = 0; n != _N; ++n) {
    x[n] *= value;
  }
  return x;
}

//! Array-scalar division.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
Array<_T, _N>&
operator/=(Array<_T, _N>& x, 
	   typename boost::call_traits<_T>::param_type value) {
#ifdef DEBUG_array_Array
  assert(value != 0);
#endif
  for (std::size_t n = 0; n != _N; ++n) {
    x[n] /= value;
  }
  return x;
}

//! Array-scalar modulus.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
Array<_T, _N>&
operator%=(Array<_T, _N>& x, 
	   typename boost::call_traits<_T>::param_type value) {
#ifdef DEBUG_array_Array
  assert(value != 0);
#endif
  for (std::size_t n = 0; n != _N; ++n) {
    x[n] %= value;
  }
  return x;
}

//! Left shift.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
Array<_T, _N>&
operator<<=(Array<_T, _N>& x, 
	    const int offset) {
  for (std::size_t n = 0; n != _N; ++n) {
    x[n] <<= offset;
  }
  return x;
}

//! Right shift.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
Array<_T, _N>&
operator>>=(Array<_T, _N>& x, 
	    const int offset) {
  for (std::size_t n = 0; n != _N; ++n) {
    x[n] >>= offset;
  }
  return x;
}

//@}
//----------------------------------------------------------------------------
//! \defgroup ArrayAssignmentArray Array Assignment Operators with an Array Operand.
//@{


//! Array-array addition.
/*! \relates Array */
template<typename _T1, typename _T2, std::size_t _N>
inline
Array<_T1, _N>&
operator+=(Array<_T1, _N>& x, const Array<_T2, _N>& y) {
  for (std::size_t n = 0; n != _N; ++n) {
    x[n] += y[n];
  }
  return x;
}

//! Array-array subtraction.
/*! \relates Array */
template<typename _T1, typename _T2, std::size_t _N>
inline
Array<_T1, _N>&
operator-=(Array<_T1, _N>& x, const Array<_T2, _N>& y) {
  for (std::size_t n = 0; n != _N; ++n) {
    x[n] -= y[n];
  }
  return x;
}

//! Array-array multiplication.
/*! \relates Array */
template<typename _T1, typename _T2, std::size_t _N>
inline
Array<_T1, _N>&
operator*=(Array<_T1, _N>& x, const Array<_T2, _N>& y) {
  for (std::size_t n = 0; n != _N; ++n) {
    x[n] *= y[n];
  }
  return x;
}

//! Array-array division.
/*! \relates Array */
template<typename _T1, typename _T2, std::size_t _N>
inline
Array<_T1, _N>&
operator/=(Array<_T1, _N>& x, const Array<_T2, _N>& y) {
  for (std::size_t n = 0; n != _N; ++n) {
#ifdef DEBUG_array_Array
    assert(y[n] != 0);
#endif
    x[n] /= y[n];
  }
  return x;
}

//! Array-array modulus.
/*! \relates Array */
template<typename _T1, typename _T2, std::size_t _N>
inline
Array<_T1, _N>&
operator%=(Array<_T1, _N>& x, const Array<_T2, _N>& y) {
  for (std::size_t n = 0; n != _N; ++n) {
#ifdef DEBUG_array_Array
    assert(y[n] != 0);
#endif
    x[n] %= y[n];
  }
  return x;
}

//! Array-array left shift.
/*! \relates Array */
template<typename _T1, typename _T2, std::size_t _N>
inline
Array<_T1, _N>&
operator<<=(Array<_T1, _N>& x, const Array<_T2, _N>& y) {
  for (std::size_t n = 0; n != _N; ++n) {
    x[n] <<= y[n];
  }
  return x;
}

//! Array-array right shift.
/*! \relates Array */
template<typename _T1, typename _T2, std::size_t _N>
inline
Array<_T1, _N>&
operator>>=(Array<_T1, _N>& x, const Array<_T2, _N>& y) {
  for (std::size_t n = 0; n != _N; ++n) {
    x[n] >>= y[n];
  }
  return x;
}

//@}
//----------------------------------------------------------------------------
//! \defgroup ArrayBinary Array Binary Operators
//@{

//! Array-scalar addition.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
Array<_T, _N>
operator+(const Array<_T, _N>& x,
	  typename boost::call_traits<_T>::param_type value) {
  Array<_T, _N> result(x);
  return result += value;
}

//! Scalar-Array addition.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
Array<_T, _N>
operator+(typename boost::call_traits<_T>::param_type value,
	  const Array<_T, _N>& x) {
  return x + value;
}

//! Array-array addition.
/*! \relates Array */
template<typename _T1, typename _T2, std::size_t _N>
inline
Array<_T1, _N>
operator+(const Array<_T1, _N>& x, const Array<_T2, _N>& y) {
  Array<_T1, _N> result(x);
  return result += y;
}

//! Array-scalar subtraction.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
Array<_T, _N>
operator-(const Array<_T, _N>& x, 
	  typename boost::call_traits<_T>::param_type value) {
  Array<_T, _N> result(x);
  return result -= value;
}

//! Scalar-array subtraction.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
Array<_T, _N>
operator-(typename boost::call_traits<_T>::param_type value,
	  const Array<_T, _N>& x) {
  return -(x - value);
}

//! Array-array subtraction.
/*! \relates Array */
template<typename _T1, typename _T2, std::size_t _N>
inline
Array<_T1, _N>
operator-(const Array<_T1, _N>& x, const Array<_T2, _N>& y) {
  Array<_T1, _N> result(x);
  return result -= y;
}

//! Array-scalar multiplication.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
Array<_T, _N>
operator*(const Array<_T, _N>& x, 
	  typename boost::call_traits<_T>::param_type value) {
  Array<_T, _N> result(x);
  return result *= value;
}

//! Scalar-array multiplication.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
Array<_T, _N>
operator*(typename boost::call_traits<_T>::param_type value,
	  const Array<_T, _N>& x) {
  return x * value;
}

//! Array-array multiplication.
/*! \relates Array */
template<typename _T1, typename _T2, std::size_t _N>
inline
Array<_T1, _N>
operator*(const Array<_T1, _N>& x, const Array<_T2, _N>& y) {
  Array<_T1, _N> result(x);
  return result *= y;
}

//! Array-scalar division.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
Array<_T, _N>
operator/(const Array<_T, _N>& x, 
	  typename boost::call_traits<_T>::param_type value) {
  Array<_T, _N> result(x);
  return result /= value;
}

//! Scalar-array division.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
Array<_T, _N>
operator/(typename boost::call_traits<_T>::param_type value,
	  const Array<_T, _N>& x) {
  Array<_T, _N> result(value);
  return result /= x;
}

//! Array-array division.
/*! \relates Array */
template<typename _T1, typename _T2, std::size_t _N>
inline
Array<_T1, _N>
operator/(const Array<_T1, _N>& x, const Array<_T2, _N>& y) {
  Array<_T1, _N> result(x);
  return result /= y;
}

//! Array-scalar modulus.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
Array<_T, _N>
operator%(const Array<_T, _N>& x, 
	  typename boost::call_traits<_T>::param_type value) {
  Array<_T, _N> result(x);
  return result %= value;
}

//! Scalar-array modulus.
/*! \relates Array */
template<typename _T, std::size_t _N>
Array<_T, _N>
operator%(typename boost::call_traits<_T>::param_type value, 
	  const Array<_T, _N>& x) {
  Array<_T, _N> result(value);
  return result %= x;
}

//! Array-array modulus.
/*! \relates Array */
template<typename _T1, typename _T2, std::size_t _N>
inline
Array<_T1, _N>
operator%(const Array<_T1, _N>& x, const Array<_T2, _N>& y) {
  Array<_T1, _N> result(x);
  return result %= y;
}

//@}
//----------------------------------------------------------------------------
//! \defgroup ArrayComparison Array Comparison
//@{

//! Return true if the arrays are equal.
/*! \relates Array */
template<typename _T1, typename _T2, std::size_t _N>
inline
bool
operator==(const Array<_T1, _N>& x, const Array<_T2, _N>& y) {
  return std::equal(x.begin(), x.end(), y.begin());
}

// CONTINUE: Do I need the rest?
    
//@}
//----------------------------------------------------------------------------
//! \defgroup ArrayFile Array File I/O
//@{

//! Write an array as space-separated numbers.
/*! \relates Array */
template<typename _T, std::size_t _N>
std::ostream&
operator<<(std::ostream& out, const Array<_T, _N>& x);
    
//! Read white space-separated numbers into an array.
/*! \relates Array */
template<typename _T, std::size_t _N>
std::istream&
operator>>(std::istream& in, Array<_T, _N>& x);

//! Write the array elements in binary format.
/*! \relates Array */
template<typename _T, std::size_t _N>
void
writeElementsBinary(std::ostream& out, const Array<_T, _N>& x);
    
//! Read the array elements in binary format.
/*! \relates Array */
template<typename _T, std::size_t _N>
void
readElementsBinary(std::istream& in, Array<_T, _N>& x);

//@}
//----------------------------------------------------------------------------
//! \defgroup ArrayMathematical Array Mathematical Functions
//@{

//! Return the sum of the components.
/*! \relates Array */
template<typename _T, std::size_t _N>
_T
sum(const Array<_T, _N>& x);
    
//! Return the product of the components.
/*! \relates Array */
template<typename _T, std::size_t _N>
_T
product(const Array<_T, _N>& x);

//! Return the minimum component.  Use < for comparison.
/*! \relates Array */
template<typename _T, std::size_t _N>
_T
min(const Array<_T, _N>& x);
    
//! Return the maximum component.  Use > for comparison.
/*! \relates Array */
template<typename _T, std::size_t _N>
_T
max(const Array<_T, _N>& x);

//! Return the dot product of the two vectors.
/*! \relates Array */
template<typename _T1, typename _T2, std::size_t _N>
_T1
dot(const Array<_T1, _N>& x, const Array<_T2, _N>& y);

//@}
//----------------------------------------------------------------------------
//! \defgroup ArrayElement Array Element-wise Functions
//@{

//! Return an Array that is element-wise the maximum of the two.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
Array<_T, _N>
computeMaximum(const Array<_T, _N>& x, const Array<_T, _N>& y) {
  Array<_T, _N> z;
  for (int n = 0; n != _N; ++n) {
    z[n] = std::max(x[n], y[n]);
  }
  return z;
}


//! Return an Array that is element-wise the minimum of the two.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
Array<_T, _N>
computeMinimum(const Array<_T, _N>& x, const Array<_T, _N>& y) {
  Array<_T, _N> z;
  for (int n = 0; n != _N; ++n) {
    z[n] = std::min(x[n], y[n]);
  }
  return z;
}





//@}
//----------------------------------------------------------------------------
//! \defgroup ArrayApply Array Apply the Standard Math Functions
//@{

//! Apply the absolute value (\f$|x|\f$) to each array element.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
void
applyAbs(Array<_T, _N>* x) {
  for (typename Array<_T, _N>::iterator i = x->begin(); i != x->end(); 
       ++i) {
    *i = std::abs(*i);
  }
}

//! Apply the inverse cosine (\f$ \cos^{-1}(x) \f$) to each array element.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
void
applyAcos(Array<_T, _N>* x) {
  for (typename Array<_T, _N>::iterator i = x->begin(); i != x->end(); 
       ++i) {
    *i = std::acos(*i);
  }
}

//! Apply the inverse sine (\f$ \sin^{-1}(x) \f$) to each array element.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
void
applyAsin(Array<_T, _N>* x) {
  for (typename Array<_T, _N>::iterator i = x->begin(); i != x->end(); 
       ++i) {
    *i = std::asin(*i);
  }
}

//! Apply the inverse tangent (\f$ \tan^{-1}(x) \f$) to each array element.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
void
applyAtan(Array<_T, _N>* x) {
  for (typename Array<_T, _N>::iterator i = x->begin(); i != x->end(); 
       ++i) {
    *i = std::atan(*i);
  }
}

//! Apply the ceiling function (\f$ \lceil x \rceil \f$) to each array element.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
void
applyCeil(Array<_T, _N>* x) {
  for (typename Array<_T, _N>::iterator i = x->begin(); i != x->end(); 
       ++i) {
    *i = std::ceil(*i);
  }
}

//! Apply the cosine (\f$ \cos(x) \f$) to each array element.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
void
applyCos(Array<_T, _N>* x) {
  for (typename Array<_T, _N>::iterator i = x->begin(); i != x->end(); 
       ++i) {
    *i = std::cos(*i);
  }
}

//! Apply the hyperbolic cosine (\f$ \cosh(x) \f$) to each array element.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
void
applyCosh(Array<_T, _N>* x) {
  for (typename Array<_T, _N>::iterator i = x->begin(); i != x->end(); 
       ++i) {
    *i = std::cosh(*i);
  }
}

//! Apply the exponential function (\f$ \mathrm{e}^x \f$) to each array element.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
void
applyExp(Array<_T, _N>* x) {
  for (typename Array<_T, _N>::iterator i = x->begin(); i != x->end(); 
       ++i) {
    *i = std::exp(*i);
  }
}

//! Apply the floor function (\f$ \lfloor x \rfloor \f$) to each array element.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
void
applyFloor(Array<_T, _N>* x) {
  for (typename Array<_T, _N>::iterator i = x->begin(); i != x->end(); 
       ++i) {
    *i = std::floor(*i);
  }
}

//! Apply the natural logarithm (\f$ \ln(x) \f$) to each array element.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
void
applyLog(Array<_T, _N>* x) {
  for (typename Array<_T, _N>::iterator i = x->begin(); i != x->end(); 
       ++i) {
    *i = std::log(*i);
  }
}

//! Apply the logarithm base 10 (\f$ \log_{10}(x) \f$) to each array element.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
void
applyLog10(Array<_T, _N>* x) {
  for (typename Array<_T, _N>::iterator i = x->begin(); i != x->end(); 
       ++i) {
    *i = std::log10(*i);
  }
}

//! Apply the sine (\f$ \sin(x) \f$) to each array element.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
void
applySin(Array<_T, _N>* x) {
  for (typename Array<_T, _N>::iterator i = x->begin(); i != x->end(); 
       ++i) {
    *i = std::sin(*i);
  }
}

//! Apply the hyperbolic sine (\f$ \sinh(x) \f$) to each array element.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
void
applySinh(Array<_T, _N>* x) {
  for (typename Array<_T, _N>::iterator i = x->begin(); i != x->end(); 
       ++i) {
    *i = std::sinh(*i);
  }
}

//! Apply the square root (\f$ \sqrt{x} \f$) to each array element.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
void
applySqrt(Array<_T, _N>* x) {
  for (typename Array<_T, _N>::iterator i = x->begin(); i != x->end(); 
       ++i) {
    *i = std::sqrt(*i);
  }
}

//! Apply the tangent (\f$ \tan(x) \f$) to each array element.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
void
applyTan(Array<_T, _N>* x) {
  for (typename Array<_T, _N>::iterator i = x->begin(); i != x->end(); 
       ++i) {
    *i = std::tan(*i);
  }
}

//! Apply the hyperbolic tangent (\f$ \tanh(x) \f$) to each array element.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
void
applyTanh(Array<_T, _N>* x) {
  for (typename Array<_T, _N>::iterator i = x->begin(); i != x->end(); 
       ++i) {
    *i = std::tanh(*i);
  }
}



//@}
//----------------------------------------------------------------------------
//! \defgroup ArrayStandard Array Standard Math Functions
//@{

//! Return the absolute value (\f$|x|\f$) applied to each array element.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
Array<_T, _N>
abs(Array<_T, _N> x) {
  applyAbs(&x);
  return x;
}

//! Return the inverse cosine (\f$ \cos^{-1}(x) \f$) appiled to each array element.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
Array<_T, _N>
acos(Array<_T, _N> x) {
  applyAcos(&x);
  return x;
}

//! Return the inverse sine (\f$ \sin^{-1}(x) \f$) applied to each array element.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
Array<_T, _N>
asin(Array<_T, _N> x) {
  applyAsin(&x);
  return x;
}

//! Return the inverse tangent (\f$ \tan^{-1}(x) \f$) applied to each array element.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
Array<_T, _N>
atan(Array<_T, _N> x) {
  applyAtan(&x);
  return x;
}

//! Return the ceiling function (\f$ \lceil x \rceil \f$) applied to each array element.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
Array<_T, _N>
ceil(Array<_T, _N> x) {
  applyCeil(&x);
  return x;
}

//! Return the cosine (\f$ \cos(x) \f$) applied to each array element.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
Array<_T, _N>
cos(Array<_T, _N> x) {
  applyCos(&x);
  return x;
}

//! Return the hyperbolic cosine (\f$ \cosh(x) \f$) applied to each array element.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
Array<_T, _N>
cosh(Array<_T, _N> x) {
  applyCosh(&x);
  return x;
}

//! Return the exponential function (\f$ \mathrm{e}^x \f$) applied to each array element.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
Array<_T, _N>
exp(Array<_T, _N> x) {
  applyExp(&x);
  return x;
}

//! Return the floor function (\f$ \lfloor x \rfloor \f$) applied to each array element.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
Array<_T, _N>
floor(Array<_T, _N> x) {
  applyFloor(&x);
  return x;
}

//! Return the natural logarithm (\f$ \ln(x) \f$) applied to each array element.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
Array<_T, _N>
log(Array<_T, _N> x) {
  applyLog(&x);
  return x;
}

//! Return the logarithm base 10 (\f$ \log_{10}(x) \f$) applied to each array element.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
Array<_T, _N>
log10(Array<_T, _N> x) {
  applyLog10(&x);
  return x;
}

//! Return the sine (\f$ \sin(x) \f$) applied to each array element.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
Array<_T, _N>
sin(Array<_T, _N> x) {
  applySin(&x);
  return x;
}

//! Return the hyperbolic sine (\f$ \sinh(x) \f$) applied to each array element.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
Array<_T, _N>
sinh(Array<_T, _N> x) {
  applySinh(&x);
  return x;
}

//! Return the square root (\f$ \sqrt{x} \f$) applied to each array element.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
Array<_T, _N>
sqrt(Array<_T, _N> x) {
  applySqrt(&x);
  return x;
}

//! Return the tangent (\f$ \tan(x) \f$) applied to each array element.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
Array<_T, _N>
tan(Array<_T, _N> x) {
  applyTan(&x);
  return x;
}

//! Return the hyperbolic tangent (\f$ \tanh(x) \f$) applied to each array element.
/*! \relates Array */
template<typename _T, std::size_t _N>
inline
Array<_T, _N>
tanh(Array<_T, _N> x) {
  applyTanh(&x);
  return x;
}

//@}

END_NAMESPACE_ARRAY

#define __array_Array_ipp__
#include "array.ipp"
#undef __array_Array_ipp__

#endif
