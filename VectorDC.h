#ifndef _VECTORDC_H_
#define _VECTORDC_H_

#include "Vector.h"

extern VectorDC *VectorDC_New(int size);

extern VectorDC *VectorDC_Del(VectorDC *pIVector);

#ifdef _MSC_VER
extern VectorDC *VectorDC_Wrap(_Dcomplex *work, int size, int inc);
#else
extern VectorDC *VectorDC_Wrap(double complex *work, int size, int inc);
#endif /*_MSC_VER*/

extern VectorDC *VectorDC_UnWrap(VectorDC *pIVector);

extern void VectorDC_Show(VectorDC *pIVector);

#ifdef _MSC_VER
extern VectorDC *VectorDC_Full(VectorDC *pIVector, _Dcomplex num);
#else
extern VectorDC *VectorDC_Full(VectorDC *pIVector, double complex num);
#endif /*_MSC_VER*/

extern VectorDC *VectorDC_Full2(VectorDC *pIVector, double num);

#ifdef _MSC_VER
extern VectorDC *VectorDC_Linspace(VectorDC *pIVector, _Dcomplex start, _Dcomplex end);
#else
extern VectorDC *VectorDC_Linspace(VectorDC *pIVector, double complex start, double complex end);
#endif /*_MSC_VER*/

extern VectorDC *VectorDC_Linspace1(VectorDC *pIVector, double start, double end);

#ifdef _MSC_VER
extern VectorDC *VectorDC_Linspace2(VectorDC *pIVector, double start, _Dcomplex end);
#else
extern VectorDC *VectorDC_Linspace2(VectorDC *pIVector, double start, double complex end);
#endif /*_MSC_VER*/

#ifdef _MSC_VER
extern VectorDC *VectorDC_Linspace3(VectorDC *pIVector, _Dcomplex start, double end);
#else
extern VectorDC *VectorDC_Linspace3(VectorDC *pIVector, double complex start, double end);
#endif /*_MSC_VER*/

extern VectorDC *VectorDC_Copy(VectorDC *pIVector1, VectorDC *pIVector2);

extern VectorD *VectorDC_Copy2(VectorDC *pIVector1, VectorD *pIVector2);

#ifdef _MSC_VER
extern _Dcomplex VectorDC_Max(VectorDC *pVector);
#else
extern double complex VectorDC_Max(VectorDC *pVector);
#endif /*_MSC_VER*/

extern double VectorDC_Max2(VectorDC *pVector);

#ifdef _MSC_VER
extern _Dcomplex VectorDC_Min(VectorDC *pVector);
#else
extern double complex VectorDC_Min(VectorDC *pVector);
#endif /*_MSC_VER*/

extern double VectorDC_Min2(VectorDC *pVector);

extern VectorDC *VectorDC_AddVectorD(VectorDC *pIVector1, VectorD *pIVector2, VectorDC *pOVector);

extern VectorD *VectorDC_AddVectorD2(VectorDC *pIVector1, VectorD *pIVector2, VectorD *pOVector);

extern VectorDC *VectorDC_AddVectorDC(VectorDC *pIVector1, VectorDC *pIVector2, VectorDC *pOVector);

extern VectorD *VectorDC_AddVectorDC2(VectorDC *pIVector1, VectorDC *pIVector2, VectorD *pOVector);

extern VectorDC *VectorDC_SubVectorD(VectorDC *pIVector1, VectorD *pIVector2, VectorDC *pOVector);

extern VectorD *VectorDC_SubVectorD2(VectorDC *pIVector1, VectorD *pIVector2, VectorD *pOVector);

extern VectorDC *VectorDC_SubVectorDC(VectorDC *pIVector1, VectorDC *pIVector2, VectorDC *pOVector);

extern VectorD *VectorDC_SubVectorDC2(VectorDC *pIVector1, VectorDC *pIVector2, VectorD *pOVector);

extern VectorDC *VectorDC_MulNumD(VectorDC *pIVector, double num, VectorDC *pOVector);

extern VectorD *VectorDC_MulNumD2(VectorDC *pIVector, double num, VectorD *pOVector);

#ifdef _MSC_VER
extern VectorDC *VectorDC_MulNumDC(VectorDC *pIVector, _Dcomplex num, VectorDC *pOVector);
#else
extern VectorDC *VectorDC_MulNumDC(VectorDC *pIVector, double complex num, VectorDC *pOVector);
#endif /*_MSC_VER*/

#ifdef _MSC_VER
extern VectorD *VectorDC_MulNumDC2(VectorDC *pIVector, _Dcomplex num, VectorD *pOVector);
#else
extern VectorD *VectorDC_MulNumDC2(VectorDC *pIVector, double complex num, VectorD *pOVector);
#endif /*_MSC_VER*/

#ifdef _MSC_VER
extern _Dcomplex VectorDC_MulVectorD(VectorDC *pIVector1, VectorD *pIVector2);
#else
extern double complex VectorDC_MulVectorD(VectorDC *pIVector1, VectorD *pIVector2);
#endif /*_MSC_VER*/

extern double VectorDC_MulVectorD2(VectorDC *pIVector1, VectorD *pIVector2);

#ifdef _MSC_VER
extern _Dcomplex VectorDC_MulVectorDC(VectorDC *pIVector1, VectorDC *pIVector2);
#else
extern double complex VectorDC_MulVectorDC(VectorDC *pIVector1, VectorDC *pIVector2);
#endif /*_MSC_VER*/

extern double VectorDC_MulVectorDC2(VectorDC *pIVector1, VectorDC *pIVector2);

extern VectorDC *VectorDC_EMulVectorD(VectorDC *pIVector1, VectorD *pIVector2, VectorDC *pOVector);

extern VectorD *VectorDC_EMulVectorD2(VectorDC *pIVector1, VectorD *pIVector2, VectorD *pOVector);

extern VectorDC *VectorDC_EMulVectorDC(VectorDC *pIVector1, VectorDC *pIVector2, VectorDC *pOVector);

extern VectorD *VectorDC_EMulVectorDC2(VectorDC *pIVector1, VectorDC *pIVector2, VectorD *pOVector);

#ifdef _MSC_VER
extern _Dcomplex VectorDC_CMulVectorD(VectorDC *pIVector1, VectorD *pIVector2);
#else
extern double complex VectorDC_CMulVectorD(VectorDC *pIVector1, VectorD *pIVector2);
#endif /*_MSC_VER*/

extern double VectorDC_CMulVectorD2(VectorDC *pIVector1, VectorD *pIVector2);

#ifdef _MSC_VER
extern _Dcomplex VectorDC_CMulVectorDC(VectorDC *pIVector1, VectorDC *pIVector2);
#else
extern double complex VectorDC_CMulVectorDC(VectorDC *pIVector1, VectorDC *pIVector2);
#endif /*_MSC_VER*/

extern double VectorDC_CMulVectorDC2(VectorDC *pIVector1, VectorDC *pIVector2);

extern double VectorDC_Norm(VectorDC *pIVector);

#ifdef _MSC_VER
extern _Dcomplex VectorDC_Norm2(VectorDC *pIVector);
#else
extern double complex VectorDC_Norm2(VectorDC *pIVector);
#endif /*_MSC_VER*/

extern double VectorDC_NormSq(VectorDC *pIVector);

#ifdef _MSC_VER
extern _Dcomplex VectorDC_NormSq2(VectorDC *pIVector);
#else
extern double complex VectorDC_NormSq2(VectorDC *pIVector);
#endif /*_MSC_VER*/

#endif /*_VECTORD_H_*/
