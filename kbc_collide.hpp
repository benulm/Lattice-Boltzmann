#ifndef STANDARD_COLLIDE_
#define STANDARD_COLLIDE_
#include <iomanip>  

#define FABIAN_


#ifdef STELLIS_
#include "no_implementation_trick.hpp"
#endif // STELLIS_

#ifdef FABIAN_
#include "implementation_trick.hpp"
#endif // FABIAN_

#ifdef COMPARE_
#include "compare.hpp"
#endif // COMPARE_

#endif // STANDARD_COLLIDE_