#ifndef _MATRIXCLASS_H_
#define _MATRIXCLASS_H_

#include "Type.h"

typedef struct _Matrix Matrix;
typedef struct _MatrixD MatrixD;
typedef struct _MatrixDC MatrixDC;

struct _Matrix
{
    Tensor parent;
    Matrix *(*Del)(Matrix *pIMatrix);
    Matrix *(*UnWrap)(Matrix *pIMatrix);
    void (*Show)(Matrix *pIMatrix);
    Matrix *(*Set)(Matrix *pIMatrix, Num *diag, Num *offDiag);
    Matrix *(*Copy)(Matrix *pIMatrix, Matrix *pOMatrix);
    Matrix *(*T)(Matrix *pIMatrix, Matrix *pOMatrix);
    Num *(*Det)(Matrix *pIMatrix, Num *pNum);
    Matrix *(*Inv)(Matrix *pIMatrix, Matrix *pOMatrix);
    Vector *(*MulVector)(Matrix *pIMatrix, Vector *pIVector, Vector *pOVector);
    Vector *(*TMulVector)(Matrix *pIMatrix, Vector *pIVector, Vector *pOVector);
    Vector *(*CMulVector)(Matrix *pIMatrix, Vector *pIVector, Vector *pOVector);
    Matrix *(*TMulSelf)(Matrix *pIMatrix, Matrix *pOMatrix);
    Matrix *(*O2Mul)(Matrix *pIMatrix1, Matrix *pIMatrix2, Matrix *pOMatrix);
    Matrix *(*MulMatrix)(Matrix *pIMatrix1, Matrix *pIMatrix2, Matrix *pOMatrix);
    Num *(*MaxDiagElem)(Matrix *pIMatrix, Num *pONum);
    Vector *(*LinearSolver)(Matrix *pILUMatrix, Vector *pIVector, Vector *pOVector);
    void (*EigenEquation)(Matrix *pIMatrix, Vector *pOVector, Matrix *pOMatrix);
};

struct _MatrixD
{
    Matrix parent;
    double *matrix;
    int ld;
    int size[2];
};

struct _MatrixDC
{
    Matrix parent;
    #ifdef _MSC_VER
    _Dcomplex *matrix;
    #else
    double complex *matrix;
    #endif /*_MSC_VER*/
    int ld;
    int size[2];
};

#endif /*_MATRIXCLASS_H_*/