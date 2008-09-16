// -*- C++ -*-

#if !defined(__array_h__)
#define __array_h__

#include "array/Array.h"
#include "array/MultiArray.h"

BEGIN_NAMESPACE_ARRAY

/*! 
  \file array.h
  \brief Includes the array classes.
*/

//=============================================================================
//=============================================================================
/*!
\mainpage Arrays

<!--------------------------------------------------------------------------->
\section array_introduction Introduction

This package provides classes for multidimensional arrays as well as a 
fixed-size %array. All of the classes and functions are defined in the
array namespace. To use this package you can either include the header
for the class you are using
\verbatim
#include "array/MultiArray.h" \endverbatim
or you can include the convenience header.
\verbatim
#include "array.h" \endverbatim
(Here I assume that <tt>stlib/packages</tt> is in your include path, i.e.
<tt>g++ -Istlib/packages ...</tt>.)

<!-- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -->
\subsection array_classes Classes

Array is a fixed-size %array that derives from std::tr1::array . For information
on "Technical Report 1" (TR1) check out 
<a href="http://www.aristeia.com/EC3E/TR1_info_frames.html">Scott Meyer's page</a>
or the 
<a href="http://en.wikipedia.org/wiki/Technical_Report_1">Wikipedia page</a>.
There is also the book "The C++ Standard Library Extentsions" by Pete Becker.
Array just adds mathematical operations to std::tr1::array .
I use it internally to describe the %array extents and as an 
index for the multidimensional arrays.

The Array class has two template parameters. Below is its declaration.
\verbatim
template<typename _T, std::size_t _N> class Array; \endverbatim
<tt>_T</tt> is the value type and <tt>_N</tt> is the size of the %array.
Below we construct an %array of double precision numbers of size 10.
\verbatim
array::Array<double, 10> a; \endverbatim

There are five multidimensional %array classes:
- MultiArray allocates its memory and has contiguous storage.
- MultiArrayRef references memory and has contiguous storage.
- MultiArrayConstRef is the constant version of MultiArrayRef.
- MultiArrayView is a view of multi-dimensional %array data.
- MultiArrayConstView is the constant version of MultiArrayView.

Each of these classes is templated on the value type and the dimension
(rank). For example MultiArray has the following declaration.
\verbatim
template<typename _T, std::size_t _Dimension>
class MultiArray; \endverbatim
The value type can be any class or built-in type. Below we construct a 
3-dimensional %array of integers.
\verbatim
array::MultiArray<int, 3> a; \endverbatim
Each of the multidimensional %array classes have different constructors.
Consult the class documentation for details.

<!-- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -->
\subsection arrayArrayExamples Array Examples

We will introduce the %array classes with a few simple examples.
We can use the Array class to store a Cartesian point in 3-D.
\verbatim
typedef array::Array<double, 3> Point; \endverbatim
We can construct a point an then set its values.
\verbatim
Point x;
x[0] = 2.0;
x[1] = 3.0;
x[2] = 5.0; \endverbatim
Alternatively, we can set the initial values with the constructor.
\verbatim
Point x(2.0, 3.0, 5.0); \endverbatim

The Array class supports the usual mathematical operations. Let 
<tt>x</tt>, <tt>y</tt>, and <tt>z</tt> be Array variables of the same
size. We can perform %array-scalar arithmetic. For example, we can add
a quantity to each element.
\verbatim
x += 1;
y = x + 3; \endverbatim

We can perform array-array arithmetic in which the %array is treated 
as a vector and the operation is performed component-wise.
\verbatim
x *= 2;
z = x - y; \endverbatim

We can compute the sum or product of the %array elements, the minimum or
maximum element, or the dot product of two arrays.
\verbatim
double sumOfElements = sum(x);
double productOfElements = product(x);
double minimum = min(x);
double maximum = max(x);
double dotProduct = dot(x, y); \endverbatim
(Note that because of the C++ function lookup rules, you don't need to specify
the array namespace to use these functions. However, explicit namespace 
qualification is safer.)

You can use the standard mathematical functions with an Array as an argument.
This returns an Array with the functions applied to each element.
\verbatim
y = abs(x);
z = exp(x); \endverbatim

Alternatively, you can apply a function to the %array elements.
\verbatim
applyAbs(&x);
applyExp(&x); \endverbatim

Equality comparison returns a boolean.
\verbatim
if (x == y) {
...
} \endverbatim

The output operator writes the elements of the %array separated by spaces.
The input operator reads white-space separated elements.
Below we read two arrays and write the sum.
\verbatim
std::cin >> x >> y;
std::cout << x + y << '\n'; \endverbatim




<!-- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -->
\subsection arrayMultidimensionalArrayExamples Multidimensional Array Examples

We can construct a multidimensional %array by specifying the index extents.
We use an Array of <tt>std::size_t</tt> for the extents.
Below we define some types and construct a 3x4 2-D %array.
\verbatim
typedef array::MultiArray<double, 2> MultiArray;
typedef MultiArray::SizeList SizeList;
MultiArray a(SizeList(3, 4)); \endverbatim
The analagous C multidimensional %array would be declared as:
\verbatim
double ca[3][4]; \endverbatim

We can treat 2-D %array as simply a container of 12 elements. The 
multidimensional %array classes support the standard STL-style interface.
For those of you who are familiar with the C++ Standard Template Library (STL),
each of the multidimensional arrays classes fulfill the requirements of
a random access container. For the rest of you, I suggest you read
"Generic Programming and the STL" by Matthew H. Austern.
\verbatim
assert(! a.empty());
assert(a.size() == 12);
assert(a.max_size() == 12);
std::fill(a.begin(), a.end(), 0);
for (std::size_t i = 0; i != a.size(); ++i) {
  a[i] = i;
} \endverbatim
Note that <tt>operator[]</tt> performs container indexing and <em>not</em>
multidimensional %array indexing.

The %array classes use <tt>operator()</tt> to perform multidimensional
indexing. The argument is an Array of integers. The index range for the
%array we have declared is [0..2]x[0..3]. Below we define the 
multidimensional index type and iterate over the elements and assign each
the sum of the indices.
\verbatim
typedef MultiArray::IndexList IndexList;
IndexList i;
for (i[0] = 0; i[0] != 3; ++i[0]) {
  for (i[1] = 0; i[1] != 4; ++i[1]) {
    a(i) = sum(i);
  }
} \endverbatim

The %array we declared has indices that are zero-offset. That is, the lower
bounds for the index ranges are zero. We can also declare an %array with
different index bases. The following %array has index ranges [1..3]x[1..4].
\verbatim
MultiArray b(SizeList(3, 4), IndexList(1, 1)); \endverbatim
The <tt>extents()</tt> and <tt>bases()</tt> accessors give the %array extents
and the index bases.
\verbatim
const SizeList extents = b.extents();
const IndexList bases = b.bases();
const IndexList upper = bases + extents;
IndexList i;
for (i[0] = bases[0]; i[0] != upper[0]; ++i[0]) {
  for (i[1] = bases[0]; i[1] != upper[1]; ++i[1]) {
    b(i) = sum(i);
  }
} \endverbatim




<!--------------------------------------------------------------------------->
\section arrayTypes Multidimensional Array Types

Each of the multidimensional %array classes is an STL-compliant random access
containers. The following types are related to this functionality. First
the types defined in both the mutable and constant %array classes.


- \c value_type is the element type of the %array.
- \c const_pointer is a pointer to a constant %array element.
- \c const_iterator is an iterator on constant elements in the %array.
- \c const_reverse_iterator is a reverse iterator on constant elements in the 
%array.
- \c const_reference is a reference to a constant %array element.
- \c size_type is the size type.
- \c difference_type is the pointer difference type.

The mutable %array classes define mutable versions of the pointers, iterators,
and references.

- \c pointer is a pointer to an %array element.
- \c iterator is an iterator on elements in the %array.
- \c reverse_iterator is a reverse iterator on elements in the %array.
- \c reference is a reference to an %array element.

The remainder of the types start with a capital letter. They support
%array indexing and working with different views of an %array.

- \c Parameter is the parameter type. This is used for passing the value
type as an argument. If the value type is a built-in type then it is
\c value_type, otherwise it is <tt>const value_type&</tt>.
- \c Index is a single %array index, which is a signed integer.
- \c IndexList is an Array of indices: <tt>Array<int, N></tt>. This type is 
used for multidimensional %array indexing.
- \c SizeList is an Array of sizes: <tt>Array<std:size_t, N></tt>. This type is 
used to describe the index extents.
- \c Storage is a class that specifies the storage order.
- \c Range is a class that represents the index range for the %array.
- \c ConstView is the class for a constant view of the %array:
<tt>MultiArrayConstView<value_type, N></tt>.
- \c View is the class for a mutable view of the %array:
<tt>MultiArrayView<value_type, N></tt>.




<!--------------------------------------------------------------------------->
\section arrayContainer The Multidimensional Array as a Random Access Container

The multidimensional %array classes provide a number of member functions for
using the %array as an STL random access container.

- \c empty() returns \c true if the %array has zero elements.
- \c size() returns the number of elements.
- \c max_size() returns the number of elements as well. (Dynamically-sized
containers like \c std::vector return the maximum number of elements that 
the container could hold, which is determined by the integer precision. 
For statically-sized containers, the \c max_size() is the same as the 
\c size().
- \c operator[]() returns the specified element in the %array.
Below we sum the elements of an %array and check it against the \c sum() 
function.
\verbatim
value_type s = 0;
for (size_type i = 0; i != a.size(); ++i) {
  s += a[i];
}
assert(s == sum(a)); \endverbatim
Note again that \c operator[] performs container indexing and not 
%array indexing. In the following example we create a 1-D %array with index 
range [-5..5] but use container indexing to initialize the elements.
\verbatim
typedef array::MultiArray<double, 1> MultiArray;
typedef MultiArray::SizeList SizeList;
typedef MultiArray::IndexList IndexList;
typedef MultiArray::size_type size_type;
MultiArray a(SizeList(11), IndexList(-5));
for (size_type i = 0; i != a.size(); ++i) {
  a[i] = i;
} \endverbatim
- \c begin() returns a random access iterator to the first element.
- \c end() returns a random access iterator to one past the last element.
Below we copy the elements of an %array \c a to a buffer.
\verbatim
std::copy(a.begin(), a.end(), buffer); \endverbatim
- \c rbegin() returns a random access reverse iterator to the last element.
- \c rend()returns a random access reverse iterator to one past the first 
element.
Below we copy the elements from \c a to \c b in reverse order.
\verbatim
assert(a.size() == b.size());
std::copy(a.rbegin(), a.rend(), b.begin()); \endverbatim
- \c fill() fills the %array with the specified value. The following two
lines are equivalent.
\verbatim
std::fill(a.begin(), a.end(), 1);
a.fill(1); \endverbatim


<!--------------------------------------------------------------------------->
\section arrayIndexing Indexing operations.

The %array classes offer only one way to perform multidimensional %array 
indexing: use \c operator() with an \c IndexList as an argument. This helps
keep the interface simple and makes it unlikely that one will confuse %array
indexing and container indexing. Below we create a 3-D %array of integers
with index range [-5..5]x[-5..5]x[-5..5] and set the element values to the 
product of the indices.
\verbatim
typedef array::MultiArray<int, 3> MultiArray;
typedef MultiArray::SizeList SizeList;
typedef MultiArray::IndexList IndexList;
MultiArray a(SizeList(11, 11, 11), IndexList(-5, -5, -5));
const IndexList lower = a.bases();
const IndexList upper = a.bases() + a.extents();
IndexList i;
for (i[0] = lower[0]; i[0] != upper[0]; ++i[0]) {
  for (i[1] = lower[1]; i[1] != upper[1]; ++i[1]) {
    for (i[2] = lower[2]; i[2] != upper[2]; ++i[2]) {
      a(i) = product(i);
    }
  }
} \endverbatim

<!--------------------------------------------------------------------------->
\section arrayIndexRange Index Ranges and Their Iterators

The IndexRange class, as the name suggests, is used to represent index
ranges. We can describe a continuous index range by its extents and bases. 
Consider an index range in 2-D. Below we construct the range [0..2]x[0..3]
\verbatim
typedef array::IndexRange<2> Range;
typedef Range::SizeList SizeList;
Range range(SizeList(3, 4)); \endverbatim
If the index range is not zero-offset, we can specify the index bases.
Below we construct the range [1..3]x[1..4].
\verbatim
typedef Range::IndexList IndexList;
Range range(SizeList(3, 4), IndexList(1, 1)); \endverbatim

We can specify a non-continuous index range by specifying steps to take
in each dimension. Below we construct the index range
[1, 3, 5]x[1, 4, 7, 10].
\verbatim
Range range(SizeList(3, 4), IndexList(1, 1), IndexList(2, 3)); \endverbatim

The \c range() member function returns the index range for each of the
the multidimensional %array classes. Below we store the range for an %array.
We check that the extents, bases, and steps are correct using the IndexRange
accessor member functions.
\verbatim
typedef MultiArray<double, 2> MultiArray;
typedef MultiArray::SizeList SizeList;
typedef MultiArray::IndexList IndexList;
typedef MultiArray::Range Range;
MultiArray a(SizeList(10, 20));
Range range = a.range();
assert(range.extents() == a.extents());
assert(range.bases() == a.bases());
assert(range.steps() == IndexList(1, 1)); \endverbatim

IndexRange is not useful by itself. Rather, we use it to construct index
range iterators. IndexRangeIterator is a random access iterator over an 
index range. Below is a function that sets each element of an %array to 
the sum of its indices.
\verbatim
template<std::size_t _N>
void
setToSum(array::MultiArray<int, _N>* a) {
  typedef array::IndexRangeIterator<_N> Iterator;
  const Iterator end = Iterator::end(a->range());
  for (Iterator i = Iterator::begin(a->range()); i != end; ++i) {
    (*a)(*i) = sum(*i);
  }
} \endverbatim
We construct iterators to the beginning and end of the index range with 
the static member functions \c begin() and \c end(). Dereferencing the 
iterator yields an index.

Next we consider a more sophisticated example. Given a N-D arrays \c a 
and \c b, we set the boundary elements of \c b to the same value 
as in \c a, and set the interior elements 
to the average of the adjacent neighbors in \c a. In 1-D,
<tt>b</tt><sub>i</sub> = 
(<tt>a</tt><sub>i-1</sub> + <tt>a</tt><sub>i+1</sub>) / 2
for the interior elements. In 2-D we have
<tt>b</tt><sub>i,j</sub> = 
(<tt>a</tt><sub>i-1,j</sub> + <tt>a</tt><sub>i+1,j</sub>
+ <tt>a</tt><sub>i,j-1</sub> + <tt>a</tt><sub>i,j+1</sub>) / 4.
\verbatim
template<typename _T, std::size_t _N>
void
laplacianAverage(const array::MultiArray<T, _N>& a, array::MultiArray<_T, _N>* b) {
  assert(a.range() == b->range());
  typedef array::IndexRangeIterator<_N> Iterator;
  typedef typename Iterator::IndexList IndexList;

  // Get the boundary values by copying the entire array.
  *b = a;

  // Skip if there are no interior points.
  if (min(a.extents()) <= 2) {
    return;
  }

  // The range for the interior elements.
  const array::IndexRange<_N> range(a.extents() - 2, a.bases() + 1);
  // Compute the interior values.
  _T s;
  IndexList index;
  const Iterator end = Iterator::end(range);
  for (Iterator i = Iterator::begin(range); i != end; ++i) {
    s = 0;
    for (std::size_t n = 0; n != _N; ++n) {
      index = *i;
      index[n] -= 1;
      s += a(index);
      index[n] += 2;
      s += a(index);
    }
    (*b)(*i) = s / _T(2 * _N);
  }
} \endverbatim


<!--------------------------------------------------------------------------->
\section arrayReferences Multidimensional Array References

The MultiArrayRef and MultiArrayConstRef classes are multidimensional arrays
that reference externally allocated data. Their constructors differ from 
MultiArray in that the first argument is a pointer to the data. Below is 
a function that receives C arrays for the data, extents and bases and 
constructs a 3-D MultiArrayRef.
\verbatim
void
foo(double* data, const std::size_t[3] extents, const int[3] bases) {
  typedef MultiArrayRef<double, 3> MultiArrayRef;
  typedef MultiArrayRef::SizeList SizeList;
  typedef MultiArrayRef::IndexList IndexList;
  MultiArrayRef a(data, SizeList(extents), IndexList(bases));
  ...
} \endverbatim
For MultiArrayConstRef the constructors take const pointers.

<!--------------------------------------------------------------------------->
\section arrayViews Multidimensional Array Views

MultiArrayView and MultiArrayConstView are %array classes that can be used
to create views of existing arrays. The easiest way to construct them is to
use the \c view() member function with a specified index range. Below we 
create a view of the interior elements of an %array. 
\verbatim
typedef array::MultiArray<double, 2> MultiArray;
typedef MultiArray::SizeList SizeList;
typedef MultiArray::IndexList IndexList;
typedef MultiArray::Range Range;
typedef MultiArray::View View;
MultiArray a(SizeList(12, 12), IndexList(-1, -1));
View interior(Range(SizeList(10, 10))); \endverbatim

When using the \c view() function, the index bases of the %array view are 
the same as the index bases of the specified range. Below is a function 
that returns a view of the interior elements of an %array.
\verbatim
template<typename _T, std::size_t _N>
array::MultiArrayView<_T, _N>
interior(array::MultiArray<_T, _N>& a) {
  typedef typename array::MultiArray<_T, _N>::Range Range;
  assert(min(a.extents()) > 2);
  return a.view(Range(a.extents() - 2, a.bases() + 1));
} \endverbatim

The %array view classes are fully fledged multidimensional arrays, just like
MultiArray, MultiArrayRef, and MultiArrayConstRef. In particular, they have 
iterators. Below we use the \c interior() function be defined above and 
set the interior elements of an %array to zero.
\verbatim
MultiArray a(SizeList(12, 12), IndexList(-1, -1));
MultiArrayView<double, 2> = interior(a);
typedef MultiArrayView<double, 2>::iterator iterator;
const iterator end = interior.end();
for (iterator i = interior.begin(); i != end; ++i) {
  *i = 0;
} \endverbatim
These iterators are standard random access iterators. In addition, you can 
get the index list that corresponds to the iterator position with the 
\c indexList member function. Below we verify that this works correctly 
for the \c interior %array.
\verbatim
for (iterator i = interior.begin(); i != end; ++i) {
  assert(*i == interior(i->indexList));
} \endverbatim

Note that the iterators for MultiArrayView and MultiArrayConstView are less
efficient than those for the other multidimensional %array types. This is
because the data for the view classes is not contiguous in memory. The 
view %array classes use ViewIterator, while the rest of the classes use
raw pointers for iterators.

<!--------------------------------------------------------------------------->
\section arrayInheritance Inheritance

MultiArrayBase is a base class for the multidimensional arrays. It stores
the index extents and bases as well as the number of elements. 

MultiArrayConstView derives from MultiArrayBase. It adds accessors 
(like \c begin() and \c end()) that make
it a constant, random access container. It also adds constant %array indexing
with \c operator().

MultiArrayView derives from MultiArrayConstView using virtual inheritance.
It adds the mutable interface.

MultiArrayConstRef derives from MultiArrayConstView using virtual inheritance.
It provides more efficient constant iterators and adds constant container 
indexing.

MultiArrayRef derives from both MultiArrayConstRef and MultiArrayView.
It provides more 
efficient mutable iterators and adds mutable container indexing.

MultiArray derives from MultiArrayRef. Other than the constructors
(and related functions), it hase the same interface as its base. It adds 
memory allocation.

When writing functions which take these arrays as arguments, it is best to
use the least derived class which supplies the needed interface.
Consider the following function.
\verbatim
template<typename _T, std::size_t _N>
void
foo(array::MultiArray<_T, _N>& a); \endverbatim
Since the function only uses the constant interface, we can make the function
more general by accepting either a MultiArrayConstRef or a MultiArrayConstView.
If the function uses iterators on the %array, then use the former because
arrays with contiguous data have more efficient iterators.
\verbatim
template<typename _T, std::size_t _N>
void
foo(array::MultiArrayConstRef<_T, _N>& a); \endverbatim
Otherwise use the latter.

Consider a function that modifies an array.
\verbatim
template<typename _T, std::size_t _N>
void
bar(array::MultiArray<_T, _N>* a); \endverbatim
We can make the function more general by accepting either a MultiArrayRef
or a MultiArrayView. Again the choice is usually based on whether we use 
the %array iterators. One needs to use a MultiArray as the argument type
only when changing the size of the %array.


<!--------------------------------------------------------------------------->
\section array_other Other Multidimensional Array Libraries.

Why, oh why, did I write a multidimensional %array package? Hasn't someone 
already done that? Yes, there are a couple of good libraries. I wrote this 
package because I wasn't entirely satisfied with other libraries.
Also, I was worried about depending on external libraries.
In writing this package I have stolen many ideas from the following.

<a href="http://www.oonumerics.org/blitz/">Blitz++</a> is a full-featured
library and has good documentation. It is a powerful library that uses 
template meta-programming to generate efficient code.
One thing that I don't like about Blitz++
is that it is not entirely STL-compliant and breaks some usual
C++ conventions. For instance, the copy 
constructor does not copy an %array, but references the data. (This is a feature
to keep naive programmers from shooting themselves in the foot.) Also, I've 
had some compilation problems with the template expressions.

Multidimensional arrays are provided in the multi_array package of the
<a href="http://www.boost.org/">Boost</a> C++ libraries. It is a nice design
and has some cool tricks. However, in using this library I found it hard to
write dimension-independent code, i.e. functions and classes where the 
dimension is a template parameter.

*/

END_NAMESPACE_ARRAY

#endif
