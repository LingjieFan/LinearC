#ifndef _MSC_COMPLEX_H_
#define _MSC_COMPLEX_H_

#include <complex.h>

#ifdef _MSC_VER

#define Img _Cbuild(0.0,1.0)

extern _Dcomplex _Caddcc(_Dcomplex x, _Dcomplex y);

extern _Dcomplex _Csubcc(_Dcomplex x, _Dcomplex y);

extern _Dcomplex _Cdivcc(_Dcomplex x, _Dcomplex y);

extern _Dcomplex _Cdivrc(double x, _Dcomplex y);

#endif /*_MSC_VER*/

#endif /*_MSC_COMPLEX_H_*/