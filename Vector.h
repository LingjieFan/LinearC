#ifndef _VECTOR_H_
#define _VECTOR_H_

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>

#include <complex.h>
#include "MSC_Complex.h"

#include "Blas.h"
#include "Lapack.h"

#include "Object.h"
#include "TensorClass.h"
#include "Num.h"
#include "VectorClass.h"
#include "VectorD.h"
#include "VectorDC.h"

extern Vector *Vector_Del(Vector *pIVector);

extern Vector *Vector_UnWrap(Vector *pIVector);

extern void Vector_Show(Vector *pIVector);

extern Vector *Vector_Full(Vector *pIVector, Num *pNum);

extern Vector *Vector_Linspace(Vector *pIVector, Num *start, Num *end);

extern Vector *Vector_Copy(Vector *pIVector, Vector *pOVector);

extern Num *Vector_Max(Vector *pVector, Num *pONum);

extern Num *Vector_Min(Vector *pVector, Num *pONum);

extern Vector *Vector_AddVector(Vector *pVector1, Vector *pVector2, Vector *pOVector);

extern Vector *Vector_SubVector(Vector *pVector1, Vector *pVector2, Vector *pOVector);

extern Vector *Vector_MulNum(Vector *pVector1, Num *pNum, Vector *pOVector);

extern Num *Vector_MulVector(Vector *pVector1, Vector *pVector2, Num *pONum);

extern Vector *Vector_EMulVector(Vector *pVector1, Vector *pVector2, Vector *pOVector);

extern Num *Vector_CMulVector(Vector *pVector1, Vector *pVector2, Num *pONum);

extern Num *Vector_Norm(Vector *pIVector, Num *pONum);

extern Num *Vector_NormSq(Vector *pIVector, Num *pONum);

#endif /*_Vector_H_*/