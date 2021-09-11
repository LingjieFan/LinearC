#ifndef _MATRIXD_H_
#define _MATRIXD_H_

#include "Matrix.h"

extern MatrixD *MatrixD_New(int row, int col);

extern MatrixD *MatrixD_Del(MatrixD *pIMatrix);

extern MatrixD *MatrixD_Wrap(double *work, int row, int col, int ld);

extern MatrixD *MatrixD_UnWrap(MatrixD *pIMatrix);

extern void MatrixD_Show(MatrixD *pIMatrix);

extern MatrixD *MatrixD_Set(MatrixD *pIMatrix, double diag, double offDiag);

#ifdef _MSC_VER
MatrixD *MatrixD_Set2(MatrixD *pIMatrix, double diag, _Dcomplex offDiag);
#else
MatrixD *MatrixD_Set2(MatrixD *pIMatrix, double diag, double complex offDiag);
#endif /*_MSC_VER*/

#ifdef _MSC_VER
MatrixD *MatrixD_Set3(MatrixD *pIMatrix, _Dcomplex diag, double offDiag);
#else
MatrixD *MatrixD_Set3(MatrixD *pIMatrix, double complex diag, double offDiag);
#endif /*_MSC_VER*/

#ifdef _MSC_VER
MatrixD *MatrixD_Set4(MatrixD *pIMatrix, _Dcomplex diag, _Dcomplex offDiag);
#else
MatrixD *MatrixD_Set4(MatrixD *pIMatrix, double complex diag, double complex offDiag);
#endif /*_MSC_VER*/

extern MatrixD *MatrixD_Copy(MatrixD *pIMatrix, MatrixD *pOMatrix);

extern MatrixDC *MatrixD_Copy2(MatrixD *pIMatrix, MatrixDC *pOMatrix);

extern MatrixD *MatrixD_T(MatrixD *pIMatrix, MatrixD *pOMatrix);

extern MatrixDC *MatrixD_T2(MatrixD *pIMatrix, MatrixDC *pOMatrix);

extern double MatrixD_Det(MatrixD *pIMatrix);

#ifdef _MSC_VER
extern _Dcomplex MatrixD_Det2(MatrixD *pIMatrix);
#else
extern double complex MatrixD_Det2(MatrixD *pIMatrix);
#endif /*_MSC_VER*/

extern MatrixD *MatrixD_Inv(MatrixD *pIMatrix, MatrixD *pOMatrix);

extern MatrixDC *MatrixD_Inv2(MatrixD *pIMatrix, MatrixDC *pOMatrix);

extern VectorD *MatrixD_MulVectorD(MatrixD *pIMatrix, VectorD *pIVector, VectorD *pOVector);

extern VectorDC *MatrixD_MulVectorD2(MatrixD *pIMatrix, VectorD *pIVector, VectorDC *pOVector);

extern VectorDC *MatrixD_MulVectorDC(MatrixD *pIMatrix, VectorDC *pIVector, VectorDC *pOVector);

extern VectorD *MatrixD_MulVectorDC2(MatrixD *pIMatrix, VectorDC *pIVector, VectorD *pOVector);

extern VectorD *MatrixD_TMulVectorD(MatrixD *pIMatrix, VectorD *pIVector, VectorD *pOVector);

extern VectorDC *MatrixD_TMulVectorD2(MatrixD *pIMatrix, VectorD *pIVector, VectorDC *pOVector);

extern VectorDC *MatrixD_TMulVectorDC(MatrixD *pIMatrix, VectorDC *pIVector, VectorDC *pOVector);

extern VectorD *MatrixD_TMulVectorDC2(MatrixD *pIMatrix, VectorDC *pIVector, VectorD *pOVector);

extern VectorD *MatrixD_CMulVectorD(MatrixD *pIMatrix, VectorD *pIVector, VectorD *pOVector);

extern VectorDC *MatrixD_CMulVectorD2(MatrixD *pIMatrix, VectorD *pIVector, VectorDC *pOVector);

extern VectorDC *MatrixD_CMulVectorDC(MatrixD *pIMatrix, VectorDC *pIVector, VectorDC *pOVector);

extern VectorD *MatrixD_CMulVectorDC2(MatrixD *pIMatrix, VectorDC *pIVector, VectorD *pOVector);

extern MatrixD *MatrixD_TMulSelf(MatrixD *pIMatrix, MatrixD *pOMatrix);

extern MatrixDC *MatrixD_TMulSelf2(MatrixD *pIMatrix, MatrixDC *pOMatrix);

extern MatrixD *MatrixD_O2Mul(MatrixD *pIMatrix1, MatrixD *pIMatrix2, MatrixD *pOMatrix);

extern MatrixDC *MatrixD_O2Mul2(MatrixD *pIMatrix1, MatrixD *pIMatrix2, MatrixDC *pOMatrix);

extern MatrixDC *MatrixD_O2Mul3(MatrixD *pIMatrix1, MatrixDC *pIMatrix2, MatrixDC *pOMatrix);

extern MatrixD *MatrixD_O2Mul4(MatrixD *pIMatrix1, MatrixDC *pIMatrix2, MatrixD *pOMatrix);

extern MatrixD *MatrixD_MulMatrixD(MatrixD *pIMatrix1, MatrixD *pIMatrix2, MatrixD *pOMatrix);

extern MatrixDC *MatrixD_MulMatrixD2(MatrixD *pIMatrix1, MatrixD *pIMatrix2, MatrixDC *pOMatrix);

extern MatrixDC *MatrixD_MulMatrixDC(MatrixD *pIMatrix1, MatrixDC *pIMatrix2, MatrixDC *pOMatrix);

extern MatrixD *MatrixD_MulMatrixDC2(MatrixD *pIMatrix1, MatrixDC *pIMatrix2, MatrixD *pOMatrix);

extern double MatrixD_MaxDiagElem(MatrixD *pIMatrix);

#ifdef _MSC_VER
extern _Dcomplex MatrixD_MaxDiagElem2(MatrixD *pIMatrix);
#else
extern double complex MatrixD_MaxDiagElem2(MatrixD *pIMatrix);
#endif /*_MSC_VER*/

extern VectorD *MatrixD_LinearSolver(MatrixD *pIMatrix, VectorD *pIVector, VectorD *pOVector);

extern VectorDC *MatrixD_LinearSolver2(MatrixD *pIMatrix, VectorD *pIVector, VectorDC *pOVector);

extern VectorDC *MatrixD_LinearSolver3(MatrixD *pIMatrix, VectorDC *pIVector, VectorDC *pOVector);

extern VectorD *MatrixD_LinearSolver4(MatrixD *pIMatrix, VectorDC *pIVector, VectorD *pOVector);

void MatrixD_EigenEquation(MatrixD *pIMatrix, VectorD *pOVector, MatrixD *pOMatrix);

void MatrixD_EigenEquation2(MatrixD *pIMatrix, VectorD *pOVector, MatrixDC *pOMatrix);

void MatrixD_EigenEquation3(MatrixD *pIMatrix, VectorDC *pOVector, MatrixD *pOMatrix);

void MatrixD_EigenEquation4(MatrixD *pIMatrix, VectorDC *pOVector, MatrixDC *pOMatrix);

#endif /*_MATRIXD_H_*/
