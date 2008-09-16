// -*- C++ -*-

/*! 
  \file hj/hj.h
  \brief N-D Hamilton-Jacobi interface functions.
*/

//----------------------------------------------------------------------------

#if !defined(__hj_hj_h__)
#define __hj_hj_h__

#include "DiffSchemeAdj.h"
#include "DiffSchemeAdjDiag.h"

#include "DistanceAdj1st.h"
#include "DistanceAdj2nd.h"
#include "DistanceAdjDiag1st.h"
#include "DistanceAdjDiag2nd.h"

#include "EikonalAdj1st.h"
#include "EikonalAdj2nd.h"
#include "EikonalAdjDiag1st.h"
#include "EikonalAdjDiag2nd.h"

#include "GridBFS.h"
#include "GridFM_BH.h"
#include "GridFM_BHDK.h"
#include "GridMCC.h"
#include "GridMCC_CA.h"
#include "GridSort.h"

#include "../geom/kernel/Point.h"
#include "../geom/grid/RegularGrid.h"

BEGIN_NAMESPACE_HJ

/*! \defgroup utility Utility Functions
  All of the functions take the solution array as an argument.
*/
// @{

//! Initialize the array.
/*!
  Set all the values in the array to infinity.
*/
template<int N, typename T, bool A>
inline
void
initialize(ads::Array<N,T,A>& array) {
  array = std::numeric_limits<T>::max();
}

//! Compute the unsigned distance from known values.
/*!
  \param array is the solution array.  All unknown values should be 
  infinity.  The distance will be computed from known values.
  \param dx is the grid spacing.
  \param maximumDistance: The distance will be computed up to maximumDistance.
  If \c maximumDistance is zero, the distance will be computed for the 
  entire grid.  The default value of \c maximumDistance is zero.
*/
template<int N, typename T, bool A>
void
computeUnsignedDistance(ads::Array<N,T,A>& array, T dx, T maximumDistance = 0);

//! Compute the signed distance from known values.
/*!
  \param array is the solution array.  All unknown values should be 
  infinity.  The distance will be computed from known values.
  \param dx is the grid spacing.
  \param maximumDistance: The distance will be computed up to +-maximumDistance.
  If \c maximumDistance is zero, the distance will be computed for the 
  entire grid.  The default value of \c maximumDistance is zero.
*/
template<int N, typename T, bool A>
void
computeSignedDistance(ads::Array<N,T,A>& array, T dx, T maximumDistance = 0);


//! Flood fill the unsigned distance.
/*!
  All unknown distances and distances greater than \c maximumDistance are
  set to \c fillValue.  The default value of \c fillValue is 
  \c maximumDistance.
*/
template<int N, typename T, bool A>
void
floodFillUnsignedDistance(ads::Array<N,T,A>& array,
			  T maximumDistance, T fillValue = 0);

//! Flood fill the signed distance.
/*!
  All unknown distances and distances greater than \c maximumDistance are
  set to \c +-fillValue.  By default, \c fillValue is set to
  \c maximumDistance.  If there are no known distances, set all the grid points
  to +fillValue.
*/
template<int N, typename T, bool A>
void
floodFillSignedDistance(ads::Array<N,T,A>& array,
			T maximumDistance, T fillValue = 0);


//! Convert a level set to the signed distance function from the iso-surface.
/*!
  \param array is the field array.  All the values must be finite.
  \param dx is the grid spacing.
  \param maximumDistance is the maximum distance to compute the signed
  distance.  If \c maximumDistance is zero the distance will be computed
  for the entire grid.  The default value of \c maximumDistance is zero.
  \param isoValue is the value of the iso-surface.  The default value
  is zero.
  \param fillValue: If the distance is not computed for the entire 
  grid, the point far away from the surface will be flood filled with
  \c +-fillValue.  By default, \c fillValue is set to \c maximumDistance.
*/
template<int N, typename T, bool A>
void
convertLevelSetToSignedDistance(ads::Array<N,T,A>& array, T dx,
				T isoValue = 0,
				T maximumDistance = 0, 
				T fillValue = 0);


// CONTINUE: Add an interface that is efficient for many fields.
//! Constant advection of the fields.
/*!
  \param field The field on with to adject values.

  \param grid is the computational grid which holds the grid extents and the 
  Cartesian domain.

  \param distance is the distance array.

  \param maximumDistance  The signed distance must have been computed up to 
  at least \c maximumDistance.  The boundary condition will be set for grid 
  points with distances in the range \f$(-maximumDistance .. 0)\f$.

  \param defaultValue  In the field, grid points with distances in the 
  range \f$(- \infty .. -maximumDistance]\f$ will be set to \c defaultValue.
*/
template<typename T, typename F, bool A1, bool A2>
void
advectConstantIntoNegativeDistance(ads::Array<3,F,A1>& field,
				   const geom::RegularGrid<3,T>& grid,
				   const ads::Array<3,T,A2>& distance, 
				   T maximumDistance,
				   F defaultValue);

// @}















/*! \defgroup utility_pointer Utility Functions: Pointer Interface
  These functions take a pointer as an argument.
*/
// @{


//! Initialize the array.
/*!
  Wrapper for
  initialize(ads::Array<2,T,A>& array)
*/
template<typename T>
inline
void
initialize(const int extentX, const int extentY, T* data) {
  initialize(ads::Array<2,T,false>(extentX, extentY, data));
}


//! Initialize the array.
/*!
  Wrapper for
  initialize(ads::Array<3,T,A>& array)
*/
template<typename T>
inline
void
initialize(const int extentX, const int extentY, const int extentZ, 
	   T* data) {
  initialize(ads::Array<3,T,false>(extentX, extentY, extentZ, data));
}




//! Compute the unsigned distance from known values.
/*!
  Wrapper for
  computeUnsignedDistance(ads::Array<2,T,A>& array,T dx,T maximumDistance)
*/
template<typename T>
inline
void
computeUnsignedDistance(const int extentX, const int extentY,
			T* data, const T dx, const T maximumDistance = 0) {
  computeUnsignedDistance(ads::Array<2,T,false>(extentX, extentY, data),
			  dx, maximumDistance);
}


//! Compute the unsigned distance from known values.
/*!
  Wrapper for
  computeUnsignedDistance(ads::Array<3,T,A>& array,T dx,T maximumDistance)
*/
template<typename T>
inline
void
computeUnsignedDistance(const int extentX, const int extentY, const int extentZ, 
			T* data, const T dx, const T maximumDistance = 0) {
  computeUnsignedDistance(ads::Array<3,T,false>
			  (extentX, extentY, extentZ, data),
			  dx, maximumDistance);
}



//! Compute the signed distance from known values.
/*!
  Wrapper for
  computeSignedDistance(ads::Array<2,T,A>& array,T dx,T maximumDistance)
*/
template<typename T>
inline
void
computeSignedDistance(const int extentX, const int extentY,
		      T* data, const T dx, const T maximumDistance = 0) {
  computeSignedDistance(ads::Array<2,T,false>(extentX, extentY, data),
			dx, maximumDistance);
}


//! Compute the signed distance from known values.
/*!
  Wrapper for
  computeSignedDistance(ads::Array<3,T,A>& array,T dx,T maximumDistance)
*/
template<typename T>
inline
void
computeSignedDistance(const int extentX, const int extentY, const int extentZ, 
		      T* data, const T dx, const T maximumDistance = 0) {
  computeSignedDistance(ads::Array<3,T,false>(extentX, extentY, extentZ, data),
			dx, maximumDistance);
}



//! Flood fill the unsigned distance.
/*!
  Wrapper for
  floodFillUnsignedDistance(ads::Array<2,T,A>& array,T maximumDistance,T fillValue)
*/
template<typename T>
inline
void
floodFillUnsignedDistance(const int extentX, const int extentY, T* data,
			  const T maximumDistance, const T fillValue = 0) {
  floodFillUnsignedDistance(ads::Array<2,T,false>
			    (extentX, extentY, data),
			    maximumDistance, fillValue);
}


//! Flood fill the unsigned distance.
/*!
  Wrapper for
  floodFillUnsignedDistance(ads::Array<3,T,A>& array,T maximumDistance,T fillValue)
*/
template<typename T>
inline
void
floodFillUnsignedDistance(const int extentX, const int extentY, 
			  const int extentZ, T* data,
			  const T maximumDistance, const T fillValue = 0) {
  floodFillUnsignedDistance(ads::Array<3,T,false>
			    (extentX, extentY, extentZ, data),
			    maximumDistance, fillValue);
}



//! Flood fill the signed distance.
/*!
  Wrapper for
  floodFillSignedDistance(ads::Array<2,T,A>& array,T maximumDistance,T fillValue)
*/
template<typename T>
inline
void
floodFillSignedDistance(const int extentX, const int extentY, T* data,
			const T maximumDistance, const T fillValue = 0) {
  floodFillSignedDistance(ads::Array<2,T,false>
			  (extentX, extentY, data),
			  maximumDistance, fillValue);
}


//! Flood fill the signed distance.
/*!
  Wrapper for
  floodFillSignedDistance(ads::Array<3,T,A>& array,T maximumDistance,T fillValue)
*/
template<typename T>
inline
void
floodFillSignedDistance(const int extentX, const int extentY, 
			const int extentZ, T* data,
			const T maximumDistance, const T fillValue = 0) {
  floodFillSignedDistance(ads::Array<3,T,false>
			  (extentX, extentY, extentZ, data),
			  maximumDistance, fillValue);
}



//! Convert a level set to the signed distance function from the iso-surface.
/*!
  Wrapper for
  convertLevelSetToSignedDistance(ads::Array<2,T,A>& array,T dx,T isoValue,T maximumDistance,T fillValue)
*/
template<typename T>
inline
void
convertLevelSetToSignedDistance(const int extentX, const int extentY, 
				T* data,
				const T dx,
				const T isoValue = 0,
				const T maximumDistance = 0, 
				const T fillValue = 0) {
  convertLevelSetToSignedDistance(ads::Array<2,T,false>
				  (extentX, extentY, data),
				  dx, isoValue, maximumDistance, fillValue);
}


//! Convert a level set to the signed distance function from the iso-surface.
/*!
  Wrapper for
  convertLevelSetToSignedDistance(ads::Array<3,T,A>& array,T dx,T isoValue,T maximumDistance,T fillValue)
*/
template<typename T>
inline
void
convertLevelSetToSignedDistance(const int extentX, const int extentY, 
				const int extentZ, T* data,
				const T dx,
				const T isoValue = 0,
				const T maximumDistance = 0, 
				const T fillValue = 0) {
  convertLevelSetToSignedDistance(ads::Array<3,T,false>
				  (extentX, extentY, extentZ, data),
				  dx, isoValue, maximumDistance, fillValue);
}



// @}

END_NAMESPACE_HJ

#define __hj_hj_ipp__
#include "hj.ipp"
#undef __hj_hj_ipp__

#endif
