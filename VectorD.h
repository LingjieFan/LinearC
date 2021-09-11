#ifndef _VECTORD_H_
#define _VECTORD_H_

#include "Vector.h"

extern VectorD *VectorD_New(int size);

extern VectorD *VectorD_Del(VectorD *pIVector);

extern VectorD *VectorD_Wrap(double *work, int size, int inc);

extern VectorD *VectorD_UnWrap(VectorD *pIVector);

extern void VectorD_Show(VectorD *pIVector);

extern VectorD *VectorD_Full(VectorD *pIVector, double num);

#ifdef _MSC_VER
extern VectorD *VectorD_Full2(VectorD *pIVector, _Dcomplex num);
#else
extern VectorD *VectorD_Full2(VectorD *pIVector, double complex num);
#endif /*_MSC_VER*/

extern VectorD *VectorD_Linspace(VectorD *pIVector, double start, double end);

#ifdef _MSC_VER
extern VectorD *VectorD_Linspace2(VectorD *pIVector, _Dcomplex start, double end);
#else
extern VectorD *VectorD_Linspace2(VectorD *pIVector, double complex start, double end);
#endif /*_MSC_VER*/

#ifdef _MSC_VER
extern VectorD *VectorD_Linspace3(VectorD *pIVector, double start, _Dcomplex end);
#else
extern VectorD *VectorD_Linspace3(VectorD *pIVector, double start, double complex end);
#endif /*_MSC_VER*/

#ifdef _MSC_VER
extern VectorD *VectorD_Linspace4(VectorD *pIVector, _Dcomplex start, _Dcomplex end);
#else
extern VectorD *VectorD_Linspace4(VectorD *pIVector, double complex start, double complex end);
#endif /*_MSC_VER*/

extern VectorD *VectorD_Copy(VectorD *pIVector, VectorD *pOVector);

extern VectorDC *VectorD_Copy2(VectorD *pIVector, VectorDC *pOVector);

extern double VectorD_Max(VectorD *pVector);

#ifdef _MSC_VER
extern _Dcomplex VectorD_Max2(VectorD *pVector);
#else
extern double complex VectorD_Max2(VectorD *pVector);
#endif /*_MSC_VER*/

extern double VectorD_Min(VectorD *pVector);

#ifdef _MSC_VER
extern _Dcomplex VectorD_Min2(VectorD *pVector);
#else
extern double complex VectorD_Min2(VectorD *pVector);
#endif /*_MSC_VER*/

extern VectorD *VectorD_AddVectorD(VectorD *pIVector1, VectorD *pIVector2, VectorD *pOVector);

extern VectorDC *VectorD_AddVectorD2(VectorD *pIVector1, VectorD *pIVector2, VectorDC *pOVector);

extern VectorDC *VectorD_AddVectorDC(VectorD *pIVector1, VectorDC *pIVector2, VectorDC *pOVector);

extern VectorD *VectorD_AddVectorDC2(VectorD *pIVector1, VectorDC *pIVector2, VectorD *pOVector);

extern VectorD *VectorD_SubVectorD(VectorD *pIVector1, VectorD *pIVector2, VectorD *pOVector);

extern VectorDC *VectorD_SubVectorD2(VectorD *pIVector1, VectorD *pIVector2, VectorDC *pOVector);

extern VectorDC *VectorD_SubVectorDC(VectorD *pIVector1, VectorDC *pIVector2, VectorDC *pOVector);

extern VectorD *VectorD_SubVectorDC2(VectorD *pIVector1, VectorDC *pIVector2, VectorD *pOVector);

extern VectorD *VectorD_MulNumD(VectorD *pIVector, double num, VectorD *pOVector);

extern VectorDC *VectorD_MulNumD2(VectorD *pIVector, double num, VectorDC *pOVector);

#ifdef _MSC_VER
extern VectorDC *VectorD_MulNumDC(VectorD *pIVector, _Dcomplex num, VectorDC *pOVector);
#else
extern VectorDC *VectorD_MulNumDC(VectorD *pIVector, double complex num, VectorDC *pOVector);
#endif /*_MSC_VER*/

#ifdef _MSC_VER
extern VectorD *VectorD_MulNumDC2(VectorD *pIVector, _Dcomplex num, VectorD *pOVector);
#else
extern VectorD *VectorD_MulNumDC2(VectorD *pIVector, double complex num, VectorD *pOVector);
#endif /*_MSC_VER*/

extern double VectorD_MulVectorD(VectorD *pIVector1, VectorD *pIVector2);

#ifdef _MSC_VER
extern _Dcomplex VectorD_MulVectorDC(VectorD *pIVector1, VectorDC *pIVector2);
#else
extern double complex VectorD_MulVectorDC(VectorD *pIVector1, VectorDC *pIVector2);
#endif /*_MSC_VER*/

extern double VectorD_MulVectorDC2(VectorD *pIVector1, VectorDC *pIVector2);

extern VectorD *VectorD_EMulVectorD(VectorD *pIVector1, VectorD *pIVector2, VectorD *pOVector);

extern VectorDC *VectorD_EMulVectorD2(VectorD *pIVector1, VectorD *pIVector2, VectorDC *pOVector);

extern VectorDC *VectorD_EMulVectorDC(VectorD *pIVector1, VectorDC *pIVector2, VectorDC *pOVector);

extern VectorD *VectorD_EMulVectorDC2(VectorD *pIVector1, VectorDC *pIVector2, VectorD *pOVector);

extern double VectorD_CMulVectorD(VectorD *pIVector1, VectorD *pIVector2);

#ifdef _MSC_VER
extern _Dcomplex VectorD_CMulVectorD2(VectorD *pIVector1, VectorD *pIVector2);
#else
extern double complex VectorD_CMulVectorD2(VectorD *pIVector1, VectorD *pIVector2);
#endif /*_MSC_VER*/

#ifdef _MSC_VER
extern _Dcomplex VectorD_CMulVectorDC(VectorD *pIVector1, VectorDC *pIVector2);
#else
extern double complex VectorD_CMulVectorDC(VectorD *pIVector1, VectorDC *pIVector2);
#endif /*_MSC_VER*/

extern double VectorD_CMulVectorDC2(VectorD *pIVector1, VectorDC *pIVector2);

extern double VectorD_Norm(VectorD *pIVector);

#ifdef _MSC_VER
extern _Dcomplex VectorD_Norm2(VectorD *pIVector);
#else
extern double complex VectorD_Norm2(VectorD *pIVector);
#endif /*_MSC_VER*/

extern double VectorD_NormSq(VectorD *pIVector);

#ifdef _MSC_VER
extern _Dcomplex VectorD_NormSq2(VectorD *pIVector);
#else
extern double complex VectorD_NormSq2(VectorD *pIVector);
#endif /*_MSC_VER*/

#endif /*_VECTORD_H_*/