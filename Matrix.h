#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>

#include <complex.h>
#include "MSC_Complex.h"

#include "Blas.h"
#include "Lapack.h"

#include "Num.h"
#include "Vector.h"
#include "MatrixType.h"
#include "MatrixD.h"
#include "MatrixDC.h"

Matrix *Matrix_Del(Matrix *pIMatrix);

Matrix *Matrix_UnWrap(Matrix *pIMatrix);

void Matrix_Show(Matrix *pIMatrix);

Matrix *Matrix_Set(Matrix *pIMatrix, Num *diag, Num *offDiag);

Matrix *Matrix_Copy(Matrix *pIMatrix, Matrix *pOMatrix);

Matrix *Matrix_T(Matrix *pIMatrix, Matrix *pOMatrix);

Num *Matrix_Det(Matrix *pIMatrix, Num *pNum);

Matrix *Matrix_Inv(Matrix *pIMatrix, Matrix *pOMatrix);

Vector *Matrix_MulVector(Matrix *pIMatrix, Vector *pIVector, Vector *pOVector);

Vector *Matrix_TMulVector(Matrix *pIMatrix, Vector *pIVector, Vector *pOVector);

Vector *Matrix_CMulVector(Matrix *pIMatrix, Vector *pIVector, Vector *pOVector);

Matrix *Matrix_TMulSelf(Matrix *pIMatrix, Matrix *pOMatrix);

Matrix *Matrix_O2Mul(Matrix *pIMatrix1, Matrix *pIMatrix2, Matrix *pOMatrix);

Matrix *Matrix_MulMatrix(Matrix *pIMatrix1, Matrix *pIMatrix2, Matrix *pOMatrix);

Num *Matrix_MaxDiagElem(Matrix *pIMatrix, Num *pONum);

Vector *Matrix_LinearSolver(Matrix *pILUMatrix, Vector *pIVector, Vector *pOVector);

void Matrix_EigenEquation(Matrix *pIMatrix, Vector *pOVector, Matrix *pOMatrix);

#endif /*_MATRIX_H_*/