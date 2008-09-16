// -*- C++ -*-

#if !defined(__hj_debug_h__)
#define __hj_debug_h__

#include "defs.h"

#include "../ads/array/FixedArray.h"

BEGIN_NAMESPACE_HJ

namespace debug {
  
  //
  // 2-D.
  //

  int _adj2[6][3] = 
    {
      {  1,  0 },
      { -1,  0 },
      {  0,  1 },
      {  0, -1 }
    };

  int _diag2[12][3] = 
    {
      {  1,  1 },
      {  1, -1 },
      { -1,  1 },
      { -1, -1 }
    };

  //! Return true if the difference is in an adjacent direction.
  inline
  bool
  is_adjacent( const ads::FixedArray<2,int>& di )
  { 
    for ( int i = 0; i != 4; ++i ) {
      if ( di[0] == _adj2[i][0] && di[1] == _adj2[i][1] ) {
	return true; 
      }
    }
    return false;
  }

  //! Return true if the difference is in a diagonal direction.
  inline
  bool
  is_diagonal( const ads::FixedArray<2,int>& di )
  { 
    for ( int i = 0; i != 4; ++i ) {
      if ( di[0] == _diag2[i][0] && di[1] == _diag2[i][1] ) {
	return true;
      }
    }
    return false;
  }

  //
  // 3-D.
  //

  int _adj3[6][3] = 
    {
      {  1,  0,  0 },
      { -1,  0,  0 },
      {  0,  1,  0 },
      {  0, -1,  0 },
      {  0,  0,  1 },
      {  0,  0, -1 }
    };

  int _diag3[12][3] = 
    {
      {  1,  1,  0 },
      {  1, -1,  0 },
      { -1,  1,  0 },
      { -1, -1,  0 },
      {  0,  1,  1 },
      {  0,  1, -1 },
      {  0, -1,  1 },
      {  0, -1, -1 },
      {  1,  0,  1 },
      { -1,  0,  1 },
      {  1,  0, -1 },
      { -1,  0, -1 }
    };

  //! Return true if the difference is in an adjacent direction.
  inline
  bool
  is_adjacent( const ads::FixedArray<3,int>& di )
  { 
    for ( int i = 0; i != 6; ++i ) {
      if ( di[0] == _adj3[i][0] && di[1] == _adj3[i][1] && 
	   di[2] == _adj3[i][2] ) {
	return true; 
      }
    }
    return false;
  }

  //! Return true if the difference is in a diagonal direction.
  inline
  bool
  is_diagonal( const ads::FixedArray<3,int>& di )
  { 
    for ( int i = 0; i != 12; ++i ) {
      if ( di[0] == _diag3[i][0] && di[1] == _diag3[i][1] && 
	   di[2] == _diag3[i][2] ) {
	return true; 
      }
    }
    return false;
  }

} // namespace debug

END_NAMESPACE_HJ

#endif
