#ifndef _VECTORCLASS_H_
#define _VECTORCLASS_H_

#include "Type.h"

typedef struct _Vector Vector;
typedef struct _VectorD VectorD;
typedef struct _VectorDC VectorDC;

struct _Vector
{
    Tensor parent;
    Vector *(*Full)(Vector *pIVector, Num *pNum);
    Vector *(*Linspace)(Vector *pIVector, Num *start, Num *end);
    Vector *(*Copy)(Vector *pIVector, Vector *pOVector);
    Num *(*Max)(Vector *pVector, Num *pONum);
    Num *(*Min)(Vector *pVector, Num *pONum);
    Vector *(*AddVector)(Vector *pVector1, Vector *pVector2, Vector *pOVector);
    Vector *(*SubVector)(Vector *pVector1, Vector *pVector2, Vector *pOVector);
    Vector *(*MulNum)(Vector *pIVector, Num *pNum, Vector *pOVector);
    Vector *(*EMulVector)(Vector *pVector1, Vector *pVector2, Vector *pOVector);
    Num *(*MulVector)(Vector *pVector1, Vector *pVector2, Num *pONum);
    Num *(*CMulVector)(Vector *pVector1, Vector *pVector2, Num *pONum);
    Num *(*Norm)(Vector *pIVector, Num *pONum);
    Num *(*NormSq)(Vector *pIVector, Num *pONum);
};

struct _VectorD
{
    Vector parent;
    double *vector;
    int inc;
    int size;
};

struct _VectorDC
{
    Vector parent;
    #ifdef _MSC_VER
    _Dcomplex *vector;
    #else
    double complex *vector;
    #endif /*_MSC_VER*/
    int inc;
    int size;
};

#endif /*_VECTORCLASS_H_*/