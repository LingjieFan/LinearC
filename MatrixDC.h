#ifndef _MATRIXDC_H_
#define _MATRIXDC_H_

#include "Matrix.h"

extern MatrixDC *MatrixDC_New(int row, int col);

extern MatrixDC *MatrixDC_Del(MatrixDC *pIMatrix);

#ifdef _MSC_VER
extern MatrixDC *MatrixDC_Wrap(_Dcomplex *work, int row, int col, int ld);
#else
extern MatrixDC *MatrixDC_Wrap(double complex *work, int row, int col, int ld);
#endif /*_MSC_VER*/

extern MatrixDC *MatrixDC_UnWrap(MatrixDC *pIMatrix);

extern void MatrixDC_Show(MatrixDC *pIMatrix);

#ifdef _MSC_VER
extern MatrixDC *MatrixDC_Set(MatrixDC *pIMatrix, _Dcomplex diag, _Dcomplex offDiag);
#else
extern MatrixDC *MatrixDC_Set(MatrixDC *pIMatrix, double complex diag, double complex offDiag);
#endif /*_MSC_VER*/

extern MatrixDC *MatrixDC_Set2(MatrixDC *pIMatrix, double diag, double offDiag);

#ifdef _MSC_VER
extern MatrixDC *MatrixDC_Set3(MatrixDC *pIMatrix, _Dcomplex diag, double offDiag);
#else
extern MatrixDC *MatrixDC_Set3(MatrixDC *pIMatrix, double complex diag, double offDiag);
#endif /*_MSC_VER*/

#ifdef _MSC_VER
extern MatrixDC *MatrixDC_Set4(MatrixDC *pIMatrix, double diag, _Dcomplex offDiag);
#else
extern MatrixDC *MatrixDC_Set4(MatrixDC *pIMatrix, double diag, double complex offDiag);
#endif /*_MSC_VER*/

extern MatrixDC *MatrixDC_Copy(MatrixDC *pIMatrix, MatrixDC *pOMatrix);

extern MatrixD *MatrixDC_Copy2(MatrixDC *pIMatrix, MatrixD *pOMatrix);

extern MatrixDC *MatrixDC_T(MatrixDC *pIMatrix, MatrixDC *pOMatrix);

extern MatrixD *MatrixDC_T2(MatrixDC *pIMatrix, MatrixD *pOMatrix);

#ifdef _MSC_VER
extern _Dcomplex MatrixDC_Det(MatrixDC *pIMatrix);
#else
extern double complex MatrixDC_Det(MatrixDC *pIMatrix);
#endif /*_MSC_VER*/

extern double MatrixDC_Det2(MatrixDC *pIMatrix);

extern MatrixDC *MatrixDC_Inv(MatrixDC *pIMatrix, MatrixDC *pOMatrix);

extern MatrixD *MatrixDC_Inv2(MatrixDC *pIMatrix, MatrixD *pOMatrix);

extern VectorDC *MatrixDC_MulVectorD(MatrixDC *pIMatrix, VectorD *pIVector, VectorDC *pOVector);

extern VectorD *MatrixDC_MulVectorD2(MatrixDC *pIMatrix, VectorD *pIVector, VectorD *pOVector);

extern VectorDC *MatrixDC_MulVectorDC(MatrixDC *pIMatrix, VectorDC *pIVector, VectorDC *pOVector);

extern VectorD *MatrixDC_MulVectorDC2(MatrixDC *pIMatrix, VectorDC *pIVector, VectorD *pOVector);

extern VectorDC *MatrixDC_TMulVectorD(MatrixDC *pIMatrix, VectorD *pIVector, VectorDC *pOVector);

extern VectorD *MatrixDC_TMulVectorD2(MatrixDC *pIMatrix, VectorD *pIVector, VectorD *pOVector);

extern VectorDC *MatrixDC_TMulVectorDC(MatrixDC *pIMatrix, VectorDC *pIVector, VectorDC *pOVector);

extern VectorD *MatrixDC_TMulVectorDC2(MatrixDC *pIMatrix, VectorDC *pIVector, VectorD *pOVector);

extern VectorDC *MatrixDC_CMulVectorD(MatrixDC *pIMatrix, VectorD *pIVector, VectorDC *pOVector);

extern VectorD *MatrixDC_CMulVectorD2(MatrixDC *pIMatrix, VectorD *pIVector, VectorD *pOVector);

extern VectorDC *MatrixDC_CMulVectorDC(MatrixDC *pIMatrix, VectorDC *pIVector, VectorDC *pOVector);

extern VectorD *MatrixDC_CMulVectorDC2(MatrixDC *pIMatrix, VectorDC *pIVector, VectorD *pOVector);

extern MatrixDC *MatrixDC_TMulSelf(MatrixDC *pIMatrix, MatrixDC *pOMatrix);

extern MatrixD *MatrixDC_TMulSelf2(MatrixDC *pIMatrix, MatrixD *pOMatrix);

extern MatrixDC *MatrixDC_O2Mul(MatrixDC *pIMatrix1, MatrixDC *pIMatrix2, MatrixDC *pOMatrix);

extern MatrixDC *MatrixDC_O2Mul2(MatrixDC *pIMatrix1, MatrixD *pIMatrix2, MatrixDC *pOMatrix);

extern MatrixD *MatrixDC_O2Mul3(MatrixDC *pIMatrix1, MatrixDC *pIMatrix2, MatrixD *pOMatrix);

extern MatrixD *MatrixDC_O2Mul4(MatrixDC *pIMatrix1, MatrixD *pIMatrix2, MatrixD *pOMatrix);

extern MatrixDC *MatrixDC_MulMatrixD(MatrixDC *pIMatrix1, MatrixD *pIMatrix2, MatrixDC *pOMatrix);

extern MatrixD *MatrixDC_MulMatrixD2(MatrixDC *pIMatrix1, MatrixD *pIMatrix2, MatrixD *pOMatrix);

extern MatrixDC *MatrixDC_MulMatrixDC(MatrixDC *pIMatrix1, MatrixDC *pIMatrix2, MatrixDC *pOMatrix);

extern MatrixD *MatrixDC_MulMatrixDC2(MatrixDC *pIMatrix1, MatrixDC *pIMatrix2, MatrixD *pOMatrix);

#ifdef _MSC_VER
extern _Dcomplex MatrixDC_MaxDiagElem(MatrixDC *pIMatrix);
#else
extern double complex MatrixDC_MaxDiagElem(MatrixDC *pIMatrix);
#endif /*_MSC_VER*/

extern double MatrixDC_MaxDiagElem2(MatrixDC *pIMatrix);

extern VectorDC *MatrixDC_LinearSolver(MatrixDC *pIMatrix, VectorDC *pIVector, VectorDC *pOVector);

extern VectorDC *MatrixDC_LinearSolver2(MatrixDC *pIMatrix, VectorD *pIVector, VectorDC *pOVector);

extern VectorD *MatrixDC_LinearSolver3(MatrixDC *pIMatrix, VectorD *pIVector, VectorD *pOVector);

extern VectorD *MatrixDC_LinearSolver4(MatrixDC *pIMatrix, VectorDC *pIVector, VectorD *pOVector);

extern void MatrixDC_EigenEquation(MatrixDC *pIMatrix, VectorDC *pOVector, MatrixDC *pOMatrix);

extern void MatrixDC_EigenEquation2(MatrixDC *pIMatrix, VectorD *pOVector, MatrixD *pOMatrix);

extern void MatrixDC_EigenEquation3(MatrixDC *pIMatrix, VectorDC *pOVector, MatrixD *pOMatrix);

extern void MatrixDC_EigenEquation4(MatrixDC *pIMatrix, VectorD *pOVector, MatrixDC *pOMatrix);

#endif /*_MSC_VER*/
