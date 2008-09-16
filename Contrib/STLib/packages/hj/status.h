// -*- C++ -*-

#if !defined(__hj_status_h__)
#define __hj_status_h__

#include "defs.h"

BEGIN_NAMESPACE_HJ

//! The status of a grid point.
typedef char Status;
//! Status for a grid point with a known solution.
const Status KNOWN = 0;
//! Status for a grid point that has been labeled, but has not yet been determined to be correct.
const Status LABELED = 1;
//! Status for a grid point that has not been labeled.
const Status UNLABELED = 2;
//! Status for a point that is outside the grid and should not be used.
const Status VOID = 3;
//! Status for a grid point that is set in the initial condition.
const Status INITIAL = 4;

// The method below is cleaner.  However, it is both slower and uses more 
// memory as the Status type will be an integer.
//enum Status { KNOWN, LABELED, UNLABELED, VOID, INITIAL };

// CONTINUE
//! An array of status variables.
//typedef ads::Array<3,Status> StatusArray;
//! Status array reference.
//typedef ads::ArrayRef<3,Status> StatusArrayRef;

END_NAMESPACE_HJ

#endif
