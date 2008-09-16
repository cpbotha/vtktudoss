// -*- C++ -*-

/*! 
  \file numerical/random/discreteFinite.h
  \brief Includes the discrete, finite random number generator classes.
*/

#if !defined(__numerical_random_discreteFinite_h__)
#define __numerical_random_discreteFinite_h__

#include "discreteFinite/DiscreteFiniteGeneratorBinarySearch.h"
#include "discreteFinite/DiscreteFiniteGeneratorBinned.h"
#include "discreteFinite/DiscreteFiniteGeneratorBinsSplitting.h"
#include "discreteFinite/DiscreteFiniteGeneratorBinsSplittingStacking.h"
#include "discreteFinite/DiscreteFiniteGeneratorCdfInversionUsingPartialPmfSums.h"
#include "discreteFinite/DiscreteFiniteGeneratorCdfInversionUsingPartialRecursiveCdf.h"
#include "discreteFinite/DiscreteFiniteGeneratorLinearSearch.h"
#include "discreteFinite/DiscreteFiniteGeneratorLinearSearchInteger.h"
#include "discreteFinite/DiscreteFiniteGeneratorRejectionBinsSplitting.h"

BEGIN_NAMESPACE_NUMERICAL

/*!
  \page numerical_random_discreteFinite General Discrete Random Number Generators

  <!----------------------------------------------------------------------->
  \section numerical_random_discreteFinite_introduction Introduction

  Many common discrete distributions have probablity masses that can be written
  in the form of a function. For example 
  \f$\mathrm{pmf}_{\mu}(n) = \mathrm{e}^{-\mu} \mu^n / n!\f$ for
  the Poisson distribution with mean \f$\mu\f$. Generators for such a
  distribution usually take advantage of the structure of the 
  probablity mass function (PMF).
  For the general discrete, finite distribution, the PMF
  does not have any special form. (Note the "finite" qualification. 
  One can work with the infinite case on paper, but It doesn't
  make sense to design algorithms for them. Representing the PMF would require
  an infinite amount of storage since the PMF has no special structure.
  The "finite" qualification will hence be taken for granted.)

  The static case is well studied. Of course one can iterate over the PMF and
  perform the inversion. This linear search has linear computational 
  complexity. A better approach is to store the cumulative mass function (CMF)
  and do a binary search. This improves the complexity to \f$\log_2(N)\f$
  for \e N events. The best method is Walker's algorithm, which uses lookup
  tables to achieve constant complexity.

  This package was designed for the dynamic case; the PMF changes after
  generating each deviate. Both a binary search on the CMF and 
  Walker's algorithm require at least \f$\mathcal{O}(N)\f$
  operations to initialize their data structures. This makes them ill-suited
  to the dynamic case. This package has algorithms that can generate deviates
  and modify probabilities efficiently.

  Note that if most of the probability masses change after
  generating each deviate, one cannot do better than the simple linear 
  search. For this dense update case, even setting the new probability 
  masses has linear complexity. Thus the overall algorithm for generating 
  deviates cannot have better than linear complexity. Most of the algorithms in 
  this package are intended for the sparse update case: after generating
  a deviate, a small number of probablity masses are changed.

  Note that in dealing with the dynamic case, it is necessary to work
  with scaled probablities. That is, the "probabilities" are not
  required to sum to unity. (With regular probabilities, changing one
  of them necessitates changing all of them.)  Then the probability
  mass is the scaled probability divided by the sum of the scaled
  probabilities.  In the following we will refer to scaled
  probabilities simply as probabilities.

  <!----------------------------------------------------------------------->
  \section numerical_random_discreteFinite_valueSearch Value Search

  Most of the algorithms for generating discrete deviates use searching.
  We will consider a few search algorithms to set the stage for our subsequent
  work. First we present algorithms for a value search. Given an array \em a of
  \em N elements (depicted below), find an element with a specified value. 
  The simplest algorithm for this is a linear search. Loop over the array 
  until we find an element with the value we are looking for. The linear 
  search has computational complexity O(<em>N</em>). 

  \image html random/discreteFinite/ElementSearchArray.png "An array of elements."

  If the elements can be ordered by value, then we can construct
  an array of pointers to the elements and sort these pointers by
  the element values. This is depicted below. At the bottom is the array of 
  elements; at the top is the array of pointers. The gray letters show
  the pointee values. (Alternatively one could use a single array of pairs
  of elements and indices.) By performing a binary search on the array of 
  sorted pointers, one can find a specified element in O(log <em>N</em>)
  operations.

  \image html random/discreteFinite/ElementSearchBinary.png "An array of pointers sorted by element value."

  For a third alternative, we could store pointers to the elements in
  a hash table. Below is a depiction of a hash table. In this case we
  show hashing with chaining. Each item in the hash table is a
  container of pointers to the original elements. To find a value we use 
  the hash function to convert that value to a bin index. Then we search
  the container in that bin.  If we determine a
  suitable hash function and ensure that the hash table has no more
  than unit load, then we can find elements in constant time. 

  \image html random/discreteFinite/ElementSearchHashing.png "A hash table containing pointers to the array elements."

  So far we have introduced three standard algorithms; we have methods 
  with linear, logarithmic and constant complexities.
  Another family of algorithms for finding an element with a specified value
  is multi-dimensional
  searching. Consider a two-dimensional search. We construct a 2-D array
  of pointers to the elements. The 2-D array has dimensions
  O(\f$\sqrt{N}\f$) by O(\f$\sqrt{N}\f$). The pointers are partially
  sorted; pointers in a given row have element values less than the 
  pointers in subsequent rows. We will also need an array of length 
  O(\f$\sqrt{N}\f$) which stores the minimum element value for the 
  pointers in each row. Then we can find an element with two linear
  searches. First we search in the 1-D array of minimum row values 
  to determine the appropriate row. Then we search over the element values in
  that row. Thus the complexity of a 2-D search is O(2 \f$\sqrt{N}\f$).
  In the figure below we show the array of minimum row values and
  the 2-D array of pointers to the elements.

  \image html random/discreteFinite/ElementSearch2D.png "The data structure for a 2-D search."

  This method can easily be generalized to higher dimensions. For an 
  M-D search we need an M-D array of pointers whose length in each dimension
  is O(\f$N^{1/M}\f$). We also need a heirarchical set of arrays each 
  of length O(\f$N^{1/M}\f$) for storing minimum element values in rows.
  The pointers in the M-D array are partially sorted, analogous to the 2-D
  case. We can find an element with M linear searches over arrays
  of length O(\f$N^{1/M}\f$). Thus an M-D search has complexity
  O(\f$M N^{1/M}\f$).

  Now we have have search methods with complexities
  O(\em N), O(\f$M N^{1/M}\f$), O(log \em N), and O(1). Note that for the case
  \em M = log \em N, the multidimensional search becomes equivalent to the
  binary search (though it has a much more complicated data structure).
  \f[
  (\log N) N^{1/ \log N} =
  \log N (2^{\log N})^{1/\log N} =
  2 \log N =
  \mathcal{O}(\log N)
  \f]
  
  <!----------------------------------------------------------------------->
  \section numerical_random_discreteFinite_accumulation Accumulation Search

  The previous section served to remind us of some standard algorithms for 
  value searching. In generating discrete deviates we will do something 
  a little different, but related. Here we will perform accumulation searches
  on an array of elements \f$\{a_n\}\f$. Assume that the element values 
  are numeric and non-negative.
  Instead of searching for the element with a particular value, we determine
  the smallest index \em n such that
  \f[
  \sum_{i=0}^n a_i > r
  \f]
  for a specified value \em r. 
  We will assume that the sum of the elements of \em a is nonzero and that
  the value \em r is positive and less than that sum.

  Accumulation searches are closely related to value searches.
  Let \f$\{c_n\}\f$ be the partial sums of 
  \em a, that is \f$c_n = \sum_{i=0}^n a_i\f$. Because the elements of 
  \em a are non-negative, the sequence \f$\{c_n\}\f$ is non-decreasing.
  In an accumulation search with the value \em r, we are looking for the 
  semi-open interval \f$(c_{n-1}, c_n]\f$ which contains the value.
  (For convenience we can define \f$c_{-1} = 0\f$.) This shows the
  relationship between value searches and accumulation searches.
  The latter is a "value" search where the "value" is the interval
  that contains \em r. We will see that algorithms for accumulation searches
  are very similar to those for value searches and that these 
  algorithms often use partial sums.

  For a linear search we can first compute the array of partial sums and then
  loop over this array to find the first element that exceeds \em r. 
  Alternatively, we could track the partial sum as we loop over the 
  array \em a. The former approach is better when we will perform several 
  searches. The computational complexity of the linear accumulation search is
  O(\em N).

  Note that the array of partial sums is sorted. If we construct this 
  array, at cost O(\em N), then we can perform a binary search to find
  the first partial sum that exceeds \em r. The binary search has 
  computational complexity O(log \em N).

  We can also apply multi-dimensional techniques to accumulation searches.
  Consider the 2-D case. Again we construct a 2-D array of pointers to 
  the elements. However, there is no need to partially sort the pointers
  by pointee value. The 1-D array holds the sum of the pointee values
  in each row. (By contrast, it holds the minimum pointee value when
  performing a value search.) To perform an accumulation search on
  the \em a, we first perform an accumulation search on the 1-D array
  of row sums. This indicates the appropriate row for the second search.
  This second accumulation search in the row of pointers determines 
  the correct element. Since there are two linear searches on arrays 
  of length \f$\mathcal{O}(\sqrt{N})\f$, the total computational complexity
  is \f$\mathcal{O}(2 \sqrt{N})\f$. For an M-D search, the 
  computational complexity is \f$\mathcal{O}(M N^{1/M})\f$.
  
  Finally, we can apply hashing techniques to perform an accumulation 
  search. One method is table lookup. We can make a table of pointers
  into the partial sums array. 
  Let \em T be the size of the table and \em S be the sum of the elements.
  The formula \f$ \lfloor T r / S \rfloor \f$, converts \em r to an 
  index in the lookup table. The n<sup>th<\sup> element in the table
  points to the first partial sum that exceeds \f$S T / n\f$. The 
  lookup table gives up a good place to start a linear search on the 
  partial sums. For some distributions of element values, this 
  approach works very well, even attaining constant expected complexity.
  However in other cases, particularly when the element values differ 
  greatly, the expected complexity is higher.
 
  Now we are ready to tackle algorithms for generating discrete deviates.
  Here the array values describe a probability mass function. The 
  n<sup>th<\sup> element is the scaled probability of the n<sup>th<\sup>
  event. Many methods for generating discrete deviates use an accumulation 
  search. Note that if \em r is a random number between 0 and
  \em S then an accumulation search will generate a discrete deviate.
  

  <!----------------------------------------------------------------------->
  \section numerical_random_discreteFinite_linear Linear Search

  The simplest method of computing a discrete deviate is CMF inversion:
  - Compute a continuous uniform deviate \em r in the range (0..1).
  - Scale \r by the sum of the probabilities.
  - The deviate is the smallest \em n such that \f$\mathrm{cmf}(n) > r\f$.
  .
  There are various methods for CMF inversion. They differ in whether they
  store the CMF, how they order the events, and how they search for \em n.

  The simplest algorithm for CMF inversion is a linear search on the PMF.
  A PMF array with 16 events is depicted below.

  \image html random/discreteFinite/Pmf.png "A probablity mass function."

  Below is one way of implementing the linear search.
  \code
  template<typename RandomAccessConstIterator, typename T>
  int
  linearSearchChopDownUnguarded(RandomAccessConstIterator begin,
                                RandomAccessConstIterator end, T r) {
    RandomAccessConstIterator i = begin;
    for ( ; i != end && (r -= *i) > 0; ++i) {
    }
    return i - begin - (i == end);
  } \endcode
  \c begin and \c end are iterators to the beginning and end of the PMF 
  array. \c r is the scaled uniform deviate. This is called a chop-down 
  search because we chop-down the value of \c r until it is non-positive.
  It is a guarded search because we check that we do not go past the end 
  of the array. Note that the final line handles the special case that 
  round-off errors make us reach the end of the array. (Using an if statement
  to do this would be a little more expensive.) The function returns the 
  discrete deviate.

  There are many ways to implement a linear search. The differ in performance,
  but they all have  (surprise) linear computational complexity. On the
  other hand, modifying a probability has constant complexity; we simply change
  an array value. The linear search method is suitable for small problems.

  <!----------------------------------------------------------------------->
  \section numerical_random_discreteFinite_linearSorting Linear Search With Sorting

  Sorting the events in the order of descending probability may improve the 
  performance of the linear search. In order to use the sorted PMF array, one 
  needs two additional arrays of integers. The first stores the index of
  the event in the original PMF array. The index array is useful when generating
  the deviate. We can efficiently go from an element in the sorted PMF array
  to an event index. The second array stores the rank of the elements in 
  the original PMF array. This is useful in modifying probabilities. Here one 
  needs to access event probabilities by their index. More concretely,
  <tt>sortedPmf[rank[i]]</tt> is the same as <tt>pmf[i]</tt>.
  The sorted PMF array along with the index and rank arrays are depicted below.

  \image html random/discreteFinite/PmfSorted.png "The probability mass function sorted in descending order."

  Note that as the event probabilities change, one needs to re-sort the PMF 
  array to maintain it in approximately sorted order. 
  Sorting the events may improve performance if the probabilities differ by 
  a large amount. Otherwise, it may just add overhead. See the performance 
  results below for timings with various probability distributions.

  <!----------------------------------------------------------------------->
  \section numerical_random_discreteFinite_binary Binary Search

  Another method for CMF inversion is to store the CMF in an array and 
  perform a binary search to generate the deviate. The CMF array is depicted
  below.

  \image html random/discreteFinite/Cmf.png "Cumulative mass function."
  
  Generating a deviate has logarithmic computational complexity, which is 
  pretty good. However, modifying a probablity has linear complexity.
  After modifying an event's probability, the CMF must be recomputed starting
  at that event. The binary search method is suitable for small problems.

  Sorting the event probabilities is applicable to the binary search method.
  For each event one accumulates the probabilities of the influencing events.
  (Event \e a influences event \e b if the occurence of the former changes
  the subsequent probability for the latter.) Here one sorts the events in
  ascending order of accumulated influencing probability. The idea is to 
  minimize the portion of the CMF one needs to rebuild after modifying 
  a probability.

  <!----------------------------------------------------------------------->
  \section numerical_random_discreteFinite_double Double Linear Search

  One can speed up the linear search method by storing an array with 
  partial PMF sums. For instance one could have an array of length 
  \e N / 2 where each element holds the sum of two probabilities. Specifically,
  element \e i holds the sum of probabilities 2 \e i and 2 \e i + 1.
  One first performs a linear search on the short array. If the first search
  returns \e n, then the deviate is either 2 \e n or 2 \e n + 1. Examining 
  those elements in the PMF array determines which one. The cost of searching
  has roughly been cut in half. 
  
  If one stores the sum of three probablities in each element of the 
  additional array, then the cost of searching is 
  \f$\mathcal{O}(N / 3 + 3)\f$. There is a linear 
  search on an array of length \e N / 3, followed by a linear search on 
  three elements in the PMF array. Choosing the additional array to have 
  size \f$\sqrt{N}\f$ yields the best complexity, namely 
  \f$\mathcal{O}(\sqrt{N})\f$. This partial PMF sum array is depicted 
  below.
  
  \image html random/discreteFinite/PartialPmfSums.png "Partial PMF sums."
  
  The double linear search does have great complexity for generating deviates,
  but it has constant complexity for modifying probabilities. To change a 
  probability, one sets an element in the PMF array and then uses the 
  difference between the old and new values to update the appropriate element
  in the partial PMF sums array. Because of its simple design, the double 
  linear search has good performance for a wide range of problem sizes.
  
  "Tripling", "Quadrupling", etc., will also work. In general, by using \e p
  arrays of sizes \e N, \f$N^{(p-1)/p}\f$, \f$\ldots\f$, \f$N^{1/p}\f$, 
  generating a deviate has complexity \f$\mathcal{O}(p N^{1/p})\f$ and
  modifying a probability has complexity \f$\mathcal{O}(p)\f$. The performance 
  of these higher order methods depends on the problem size. Doubling often
  yields the most bang for the buck.

  <!----------------------------------------------------------------------->
  \section numerical_random_discreteFinite_partial Partial Recursive CMF

  One can build a partial recursive CMF that enables generating deviates
  and modifying probabilities in \f$\mathcal{O}(\log_2 N)\f$ time. The 
  process of building this data structure is shown below.

  \image html random/discreteFinite/PartialRecursiveCmf.png "The process of building the partial recursive CMF."
  
  One starts with the PMF. To every second element, add the previous element.
  Then to every fourth element \e i, add the element at position \e i - 2.
  Then to every eighth element \e i, add the element at position \e i - 4.
  After \f$\log_2 N\f$ steps, the final element holds the sum of the 
  probabilities.

  One can generate a deviate in \f$\mathcal{O}(\log_2 N) + 1\f$ steps.
  Modifying a probablity necessitates updating at most 
  \f$\mathcal{O}(\log_2 N) + 1\f$ elements of the partial recursive CMF array.
  See the source code for details.

  There are re-orderings of the above partial recursive CMF that also work.
  Below we show an alternative. In the order above, searching progresses 
  back-to-front and probability modifications proceed front-to-back.
  Vice-versa for the ordering below.

  \image html random/discreteFinite/PartialRecursiveCmfAlternate.png "Another partial recursive CMF."

  <!----------------------------------------------------------------------->
  \section numerical_random_discreteFinite_rejection Rejection
  
  Draw a rectangular around the PMF array. This is illustrated below.
  The rejection method for generating deviates is:
  - Randomly pick a point in the bounding box. 
  - If you hit one of the event boxes, its index is the discrete deviate.
  - Otherwise, pick again.

  \image html random/discreteFinite/Rejection.png "The rejection method."
  
  To pick a random point, you could use two random numbers, one for each 
  coordinate. However, a better method is to split a single random integer.
  For the case above, we could use the first four bits to pick a bin and use
  the remaining bits to determine the height.
  
  The efficiency of the rejection method is determined by the area of the 
  event boxes divided by the area of the bounding box. If the event 
  probabilities are similar, the efficiency will be high. If the probabilities
  differ by large amounts, the efficiency will be low. 
  
  <!----------------------------------------------------------------------->
  \section numerical_random_discreteFinite_binning Rejection with Binning.

  The rejection method is interesting, but it is not a useful
  technique by itself. Now we'll combine it with binning to obtain on
  optimal method for generating discrete deviates and modifying
  probabilities.  We distribute the event probabilities accross a
  number of bins in order to maximize the efficiency of the rejection
  method. Below we depict packing the probabilities into 32 bins. Each
  event is given one or more bins.  For this case, we can use the
  first five bits of a random integer to pick a bin and rest to
  compute a height. This gives us our random point for the rejection
  method.

  \image html random/discreteFinite/RejectionBinsSplitting.png "The rejection method with bins and splitting."
  
  As with many sophisticated methods, the devil is in the details. What is 
  an appropriate number of bins? How do you pack the bins in order to minimize
  the height of the bounding box? As the probabilities change, the efficiency
  may degrade. When should should you re-pack the bins? For answers to these
  questions, consult the class documentation and source code. Skipping to 
  the punchline: You can guarantee a high efficiency, and for an efficiency 
  \e E, the expected computational complexity of generating a deviate is 
  \f$\mathcal{O}(1/E)\f$. The data structure can be designed so that modifying
  an event probability involves updating 
  a single bin. Thus we have constant complexity for generating deviates
  and modifying probabilities.
  
  <!----------------------------------------------------------------------->
  \section numerical_random_discreteFinite_classes Classes

  This package provides the the following functors for computing 
  discrete, finite random deviates.
  - DiscreteFiniteGeneratorBinarySearch<false,Generator,T>
  - DiscreteFiniteGeneratorBinarySearch<true,Generator,T>
  - DiscreteFiniteGeneratorBinned
  - DiscreteFiniteGeneratorBinned
  - DiscreteFiniteGeneratorBinsSplitting<false,UseImmediateUpdate,BinConstants,Generator>
  - DiscreteFiniteGeneratorBinsSplitting<true,UseImmediateUpdate,BinConstants,Generator>
  - DiscreteFiniteGeneratorBinsSplittingStacking<false,UseImmediateUpdate,BinConstants,Generator>
  - DiscreteFiniteGeneratorBinsSplittingStacking<true,UseImmediateUpdate,BinConstants,Generator>
  - DiscreteFiniteGeneratorCdfInversionUsingPartialPmfSums
  - DiscreteFiniteGeneratorCdfInversionUsingPartialRecursiveCdf
  - DiscreteFiniteGeneratorLinearSearch
  - DiscreteFiniteGeneratorLinearSearchInteger
  - DiscreteFiniteGeneratorRejectionBinsSplitting<false,ExactlyBalance,BinConstants,Generator>
  - DiscreteFiniteGeneratorRejectionBinsSplitting<true,ExactlyBalance,BinConstants,Generator>

  I have implemented each of the discrete, finite deviate generators as 
  an <em>adaptable generator</em>, a functor that takes no arguments.
  (See \ref numerical_random_austern1999 "Generic Programming and the STL".)
  The classes are templated on the floating point number type
  (\c double by default) and the discrete, uniform generator type
  (DiscreteUniformGeneratorMt19937 by default).
  Below are a few ways of constructing a discrete, finite generator.
  \code
  // Use the default number type (double) and the default discrete, uniform generator.
  typedef numerical::DiscreteFiniteGeneratorBinarySearch<> Generator;
  Generator::DiscreteUniformGenerator uniform;
  Generator generator(&uniform); \endcode
  \code
  // Use single precision numbers and the default uniform deviate generator.
  typedef numerical::DiscreteFiniteGeneratorBinary<float> Generator;
  Generator::DiscreteUniformGenerator uniform;
  Generator generator(&uniform); \endcode
  \code
  // Use double precision numbers and the Mersenne twister generator from the GSL.
  typedef numerical::DiscreteFiniteGeneratorBinarySearch<double, numerical::DiscreteUniformGeneratorMt19937Gsl> Generator;
  Generator::DiscreteUniformGenerator uniform;
  Generator generator(&uniform); \endcode


  Each class defines the following types.
  - \c Number is the floating point number type.
  - \c DiscreteUniformGenerator is the discrete uniform generator type.
  - \c argument_type is \c void.
  - \c result_type is \c int.
  .
  The generators may be seeded with the seed() member function.  You can
  access the discrete uniform generator with getDiscreteUniformGenerator() .

  Each generator has the \c operator()() member function which generates 
  a deviate.
  \code
  int deviate = generator(); \endcode

  The tables below give the execution times (in nanoseconds) 
  for calling the discrete, finite deviate generators.

  Unit probabilities.
  \htmlinclude randomDiscreteFiniteStaticUnit.txt

  Uniformly distributed probabilities.
  \htmlinclude randomDiscreteFiniteStaticUniform.txt

  Geometric series with a factor of 10.
  \htmlinclude randomDiscreteFiniteStaticGeometric1.txt

  Geometric series with a factor of 100.
  \htmlinclude randomDiscreteFiniteStaticGeometric2.txt

  Geometric series with a factor of 1000.
  \htmlinclude randomDiscreteFiniteStaticGeometric3.txt

  Next we measure the execution time for generating a deviate and 
  trivially modifying that probability.  That is, we call setPmf()
  for the drawn deviate, but do not change the value of the probability.

  Unit probabilities.
  \htmlinclude randomDiscreteFiniteDynamicTrivialUnit.txt

  Uniformly distributed probabilities.
  \htmlinclude randomDiscreteFiniteDynamicTrivialUniform.txt
*/

/*
3 * x + 5 % 16

0  7
1  10
2  3
3  14
4  15
5  2
6  11
7  6
8  8*
9  13
10 12
11 9
12 0
13 5
14 4
15 1

0  h
1  k
2  d
3  o
4  p
5  c
6  l
7  g
8  i
9  n
10 m
11 j
12 a
13 f
14 e
15 b

*/

END_NAMESPACE_NUMERICAL

#endif
