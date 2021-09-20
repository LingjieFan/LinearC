#include "Matrix.h"

Matrix *Matrix_Del(Matrix *pIMatrix)
{
    if(pIMatrix == NULL)
    {
        return NULL;
    }

    return (Matrix *)(*((Object *)pIMatrix)->Del)((Object *)pIMatrix);
}

Object *_MatrixD_Del(Object *pIMatrix)
{
    return (Object *)MatrixD_Del((MatrixD *)pIMatrix);
}

Object *_MatrixDC_Del(Object *pIMatrix)
{
    return (Object *)MatrixDC_Del((MatrixDC *)pIMatrix);
}

Matrix *Matrix_UnWrap(Matrix *pIMatrix)
{
    if(pIMatrix == NULL)
    {
        return NULL;
    }

    return (Matrix *)(*((Object *)pIMatrix)->Del)((Object *)pIMatrix);
}

Object *_MatrixD_UnWrap(Object *pIMatrix)
{
    return (Object *)MatrixD_UnWrap((MatrixD *)pIMatrix);
}

Object *_MatrixDC_UnWrap(Object *pIMatrix)
{
    return (Object *)MatrixDC_UnWrap((MatrixDC *)pIMatrix);
}

void Matrix_Show(Matrix *pIMatrix)
{
    (*((Tensor *)pIMatrix)->Show)((Tensor *)pIMatrix);
}

void _MatrixD_Show(Tensor *pIMatrix)
{
    MatrixD_Show((MatrixD *)pIMatrix);
}

void _MatrixDC_Show(Tensor *pIMatrix)
{
    MatrixDC_Show((MatrixDC *)pIMatrix);
}

Matrix *Matrix_Set(Matrix *pIMatrix, Num *diag, Num *offDiag)
{
    return (*pIMatrix->Set)(pIMatrix, diag, offDiag);
}

Matrix *_MatrixD_Set(Matrix *pIMatrix, Num *diag, Num *offDiag)
{
    if(instanceof(diag, NUMD))
    {
        if(instanceof(offDiag, NUMD))
        {
            MatrixD_Set((MatrixD *)pIMatrix, ((NumD *)diag)->num, ((NumD *)offDiag)->num);
        }
        else if(instanceof(offDiag, NUMDC))
        {
            MatrixD_Set2((MatrixD *)pIMatrix, ((NumD *)diag)->num, ((NumDC *)offDiag)->num);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else if(instanceof(offDiag, NUMDC))
    {
        if(instanceof(offDiag, NUMD))
        {
            MatrixD_Set3((MatrixD *)pIMatrix, ((NumDC *)diag)->num, ((NumD *)offDiag)->num);
        }
        else if(instanceof(offDiag, NUMDC))
        {
            MatrixD_Set4((MatrixD *)pIMatrix, ((NumDC *)diag)->num, ((NumDC *)offDiag)->num);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }

    return pIMatrix;
}

Matrix *_MatrixDC_Set(Matrix *pIMatrix, Num *diag, Num *offDiag)
{
    if(instanceof(diag, NUMD))
    {
        if(instanceof(offDiag, NUMD))
        {
            MatrixDC_Set2((MatrixDC *)pIMatrix, ((NumD *)diag)->num, ((NumD *)offDiag)->num);
        }
        else if(instanceof(offDiag, NUMDC))
        {
            MatrixDC_Set4((MatrixDC *)pIMatrix, ((NumD *)diag)->num, ((NumDC *)offDiag)->num);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else if(instanceof(diag, NUMDC))
    {
        if(instanceof(offDiag, NUMD))
        {
            MatrixDC_Set3((MatrixDC *)pIMatrix, ((NumDC *)diag)->num, ((NumD *)offDiag)->num);
        }
        else if(instanceof(offDiag, NUMDC))
        {
            MatrixDC_Set((MatrixDC *)pIMatrix, ((NumDC *)diag)->num, ((NumDC *)offDiag)->num);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pIMatrix;
}

Matrix *Matrix_Copy(Matrix *pIMatrix, Matrix *pOMatrix)
{
    return (*pIMatrix->Copy)(pIMatrix, pOMatrix);
}

Matrix *_MatrixD_Copy(Matrix *pIMatrix, Matrix *pOMatrix)
{
    if(instanceof(pOMatrix, MATRIXD))
    {
        MatrixD_Copy((MatrixD *)pIMatrix, (MatrixD *)pOMatrix);
    }
    else if(instanceof(pOMatrix, MATRIXDC))
    {
        MatrixD_Copy2((MatrixD *)pIMatrix, (MatrixDC *)pOMatrix);
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pOMatrix;
}

Matrix *_MatrixDC_Copy(Matrix *pIMatrix, Matrix *pOMatrix)
{
    if(instanceof(pOMatrix, MATRIXD))
    {
        MatrixDC_Copy2((MatrixDC *)pIMatrix, (MatrixD *)pOMatrix);
    }
    else if(instanceof(pOMatrix, MATRIXDC))
    {
        MatrixDC_Copy((MatrixDC *)pIMatrix, (MatrixDC *)pOMatrix);
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pOMatrix;
}

Matrix *Matrix_T(Matrix *pIMatrix, Matrix *pOMatrix)
{
    return (*pIMatrix->T)(pIMatrix, pOMatrix);
}

Matrix *_MatrixD_T(Matrix *pIMatrix, Matrix *pOMatrix)
{
    if(instanceof(pOMatrix, MATRIXD))
    {
        MatrixD_T((MatrixD *)pIMatrix, (MatrixD *)pOMatrix);
    }
    else if(instanceof(pOMatrix, MATRIXDC))
    {
        MatrixD_T2((MatrixD *)pIMatrix, (MatrixDC *)pOMatrix);
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pOMatrix;
}

Matrix *_MatrixDC_T(Matrix *pIMatrix, Matrix *pOMatrix)
{
    if(instanceof(pOMatrix, MATRIXD))
    {
        MatrixDC_T2((MatrixDC *)pIMatrix, (MatrixD *)pOMatrix);
    }
    else if(instanceof(pOMatrix, MATRIXDC))
    {
        MatrixDC_T((MatrixDC *)pIMatrix, (MatrixDC *)pOMatrix);
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pOMatrix;
}

Num *Matrix_Det(Matrix *pIMatrix, Num *pNum)
{
    return (*pIMatrix->Det)(pIMatrix, pNum);
}

Num *_MatrixD_Det(Matrix *pIMatrix, Num *pNum)
{
    if(instanceof(pNum, NUMD))
    {
        ((NumD *)pNum)->num = MatrixD_Det((MatrixD *)pIMatrix);
    }
    else if(instanceof(pNum, NUMDC))
    {
        ((NumDC *)pNum)->num = MatrixD_Det2((MatrixD *)pIMatrix);
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pNum;
}

Num *_MatrixDC_Det(Matrix *pIMatrix, Num *pNum)
{
    if(instanceof(pNum, NUMD))
    {
        ((NumD *)pNum)->num = MatrixDC_Det2((MatrixDC *)pIMatrix);
    }
    else if(instanceof(pNum, NUMDC))
    {
        ((NumDC *)pNum)->num = MatrixDC_Det((MatrixDC *)pIMatrix);
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pNum;
}

Matrix *Matrix_Inv(Matrix *pIMatrix, Matrix *pOMatrix)
{
    return (*pIMatrix->Inv)(pIMatrix, pOMatrix);
}

Matrix *_MatrixD_Inv(Matrix *pIMatrix, Matrix *pOMatrix)
{
    if(instanceof(pOMatrix, MATRIXD))
    {
        MatrixD_Inv((MatrixD *)pIMatrix, (MatrixD *)pOMatrix);
    }
    else if(instanceof(pOMatrix, MATRIXDC))
    {
        MatrixD_Inv2((MatrixD *)pIMatrix, (MatrixDC *)pOMatrix);
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pOMatrix;
}

Matrix *_MatrixDC_Inv(Matrix *pIMatrix, Matrix *pOMatrix)
{
    if(instanceof(pOMatrix, MATRIXD))
    {
        MatrixDC_Inv2((MatrixDC *)pIMatrix, (MatrixD *)pOMatrix);
    }
    else if(instanceof(pOMatrix, MATRIXDC))
    {
        MatrixDC_Inv((MatrixDC *)pIMatrix, (MatrixDC *)pOMatrix);
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pOMatrix;
}

Vector *Matrix_MulVector(Matrix *pIMatrix, Vector *pIVector, Vector *pOVector)
{
    return (*pIMatrix->MulVector)(pIMatrix, pIVector, pOVector);
}

Vector *_MatrixD_MulVector(Matrix *pIMatrix, Vector *pIVector, Vector *pOVector)
{
    if(instanceof(pOVector, VECTORD))
    {
        if(instanceof(pIVector, VECTORD))
        {
            MatrixD_MulVectorD((MatrixD *)pIMatrix, (VectorD *)pIVector, (VectorD *)pOVector);
        }
        else if(instanceof(pIVector, VECTORDC))
        {
            MatrixD_MulVectorDC2((MatrixD *)pIMatrix, (VectorDC *)pIVector, (VectorD *)pOVector);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else if(instanceof(pOVector, VECTORDC))
    {
        if(instanceof(pIVector, VECTORD))
        {
            MatrixD_MulVectorD2((MatrixD *)pIMatrix, (VectorD *)pIVector, (VectorDC *)pOVector);
        }
        else if(instanceof(pIVector, VECTORDC))
        {
            MatrixD_MulVectorDC((MatrixD *)pIMatrix, (VectorDC *)pIVector, (VectorDC *)pOVector);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pOVector;
}

Vector *_MatrixDC_MulVector(Matrix *pIMatrix, Vector *pIVector, Vector *pOVector)
{
    if(instanceof(pOVector, VECTORD))
    {
        if(instanceof(pIVector, VECTORD))
        {
            MatrixDC_MulVectorD2((MatrixDC *)pIMatrix, (VectorD *)pIVector, (VectorD *)pOVector);
        }
        else if(instanceof(pIVector, VECTORDC))
        {
            MatrixDC_MulVectorDC2((MatrixDC *)pIMatrix, (VectorDC *)pIVector, (VectorD *)pOVector);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else if(instanceof(pOVector, VECTORDC))
    {
        if(instanceof(pIVector, VECTORD))
        {
            MatrixDC_MulVectorD((MatrixDC *)pIMatrix, (VectorD *)pIVector, (VectorDC *)pOVector);
        }
        else if(instanceof(pIVector, VECTORDC))
        {
            MatrixDC_MulVectorDC((MatrixDC *)pIMatrix, (VectorDC *)pIVector, (VectorDC *)pOVector);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pOVector;
}

Vector *Matrix_TMulVector(Matrix *pIMatrix, Vector *pIVector, Vector *pOVector)
{
    return (*pIMatrix->MulVector)(pIMatrix, pIVector, pOVector);
}

Vector *_MatrixD_TMulVector(Matrix *pIMatrix, Vector *pIVector, Vector *pOVector)
{
    if(instanceof(pOVector, VECTORD))
    {
        if(instanceof(pIVector, VECTORD))
        {
            MatrixD_TMulVectorD((MatrixD *)pIMatrix, (VectorD *)pIVector, (VectorD *)pOVector);
        }
        else if(instanceof(pIVector, VECTORDC))
        {
            MatrixD_TMulVectorDC2((MatrixD *)pIMatrix, (VectorDC *)pIVector, (VectorD *)pOVector);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else if(instanceof(pOVector, VECTORDC))
    {
        if(instanceof(pIVector, VECTORD))
        {
            MatrixD_TMulVectorD2((MatrixD *)pIMatrix, (VectorD *)pIVector, (VectorDC *)pOVector);
        }
        else if(instanceof(pIVector, VECTORDC))
        {
            MatrixD_TMulVectorDC((MatrixD *)pIMatrix, (VectorDC *)pIVector, (VectorDC *)pOVector);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pOVector;
}

Vector *_MatrixDC_TMulVector(Matrix *pIMatrix, Vector *pIVector, Vector *pOVector)
{
    if(instanceof(pOVector, VECTORD))
    {
        if(instanceof(pIVector, VECTORD))
        {
            MatrixDC_TMulVectorD2((MatrixDC *)pIMatrix, (VectorD *)pIVector, (VectorD *)pOVector);
        }
        else if(instanceof(pIVector, VECTORDC))
        {
            MatrixDC_TMulVectorDC2((MatrixDC *)pIMatrix, (VectorDC *)pIVector, (VectorD *)pOVector);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else if(instanceof(pOVector, VECTORDC))
    {
        if(instanceof(pIVector, VECTORD))
        {
            MatrixDC_TMulVectorD((MatrixDC *)pIMatrix, (VectorD *)pIVector, (VectorDC *)pOVector);
        }
        else if(instanceof(pIVector,VECTORDC))
        {
            MatrixDC_TMulVectorDC((MatrixDC *)pIMatrix, (VectorDC *)pIVector, (VectorDC *)pOVector);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pOVector;
}

Vector *Matrix_CMulVector(Matrix *pIMatrix, Vector *pIVector, Vector *pOVector)
{
    return (*pIMatrix->CMulVector)(pIMatrix, pIVector, pOVector);
}

Vector *_MatrixD_CMulVector(Matrix *pIMatrix, Vector *pIVector, Vector *pOVector)
{
    if(instanceof(pOVector, VECTORD))
    {
        if(instanceof(pIVector, VECTORD))
        {
            MatrixD_CMulVectorD((MatrixD *)pIMatrix, (VectorD *)pIVector, (VectorD *)pOVector);
        }
        else if(instanceof(pIVector, VECTORDC))
        {
            MatrixD_CMulVectorDC2((MatrixD *)pIMatrix, (VectorDC *)pIVector, (VectorD *)pOVector);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else if(instanceof(pOVector, VECTORDC))
    {
        if(instanceof(pIVector, VECTORD))
        {
            MatrixD_CMulVectorD2((MatrixD *)pIMatrix, (VectorD *)pIVector, (VectorDC *)pOVector);
        }
        else if(instanceof(pIVector, VECTORDC))
        {
            MatrixD_CMulVectorDC((MatrixD *)pIMatrix, (VectorDC *)pIVector, (VectorDC *)pOVector);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pOVector;
}

Vector *_MatrixDC_CMulVector(Matrix *pIMatrix, Vector *pIVector, Vector *pOVector)
{
    if(instanceof(pOVector, VECTORD))
    {
        if(instanceof(pIVector, VECTORD))
        {
            MatrixDC_CMulVectorD2((MatrixDC *)pIMatrix, (VectorD *)pIVector, (VectorD *)pOVector);
        }
        else if(instanceof(pIVector, VECTORDC))
        {
            MatrixDC_CMulVectorDC2((MatrixDC *)pIMatrix, (VectorDC *)pIVector, (VectorD *)pOVector);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else if(instanceof(pOVector, VECTORDC))
    {
        if(instanceof(pIVector, VECTORD))
        {
            MatrixDC_CMulVectorD((MatrixDC *)pIMatrix, (VectorD *)pIVector, (VectorDC *)pOVector);
        }
        else if(instanceof(pIVector, VECTORDC))
        {
            MatrixDC_CMulVectorDC((MatrixDC *)pIMatrix, (VectorDC *)pIVector, (VectorDC *)pOVector);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pOVector;
}

Matrix *Matrix_TMulSelf(Matrix *pIMatrix, Matrix *pOMatrix)
{
    return (*pIMatrix->TMulSelf)(pIMatrix, pOMatrix);
}

Matrix *_MatrixD_TMulSelf(Matrix *pIMatrix, Matrix *pOMatrix)
{
    if(instanceof(pOMatrix, MATRIXD))
    {
        MatrixD_TMulSelf((MatrixD *)pIMatrix, (MatrixD *)pOMatrix);
    }
    else if(instanceof(pOMatrix, MATRIXDC))
    {
        MatrixD_TMulSelf2((MatrixD *)pIMatrix, (MatrixDC *)pOMatrix);
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pOMatrix;
}

Matrix *_MatrixDC_TMulSelf(Matrix *pIMatrix, Matrix *pOMatrix)
{
    if(instanceof(pOMatrix, MATRIXD))
    {
        MatrixDC_TMulSelf2((MatrixDC *)pIMatrix, (MatrixD *)pOMatrix);
    }
    else if(instanceof(pOMatrix, MATRIXDC))
    {
        MatrixDC_TMulSelf((MatrixDC *)pIMatrix, (MatrixDC *)pOMatrix);
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pOMatrix;
}

Matrix *Matrix_O2Mul(Matrix *pIMatrix1, Matrix *pIMatrix2, Matrix *pOMatrix)
{
    return (*pIMatrix1->O2Mul)(pIMatrix1,pIMatrix2,pOMatrix);
}

Matrix *_MatrixD_O2Mul(Matrix *pIMatrix1, Matrix *pIMatrix2, Matrix *pOMatrix)
{
    if(instanceof(pOMatrix, MATRIXD))
    {
        if(instanceof(pIMatrix2, MATRIXD))
        {
            MatrixD_O2Mul((MatrixD *)pIMatrix1, (MatrixD *)pIMatrix2, (MatrixD *)pOMatrix);
        }
        else if(instanceof(pIMatrix2, MATRIXDC))
        {
            MatrixD_O2Mul4((MatrixD *)pIMatrix1, (MatrixDC *)pIMatrix2, (MatrixD *)pOMatrix);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else if(instanceof(pOMatrix, MATRIXDC))
    {
        if(instanceof(pIMatrix2, MATRIXD))
        {
            MatrixD_O2Mul2((MatrixD *)pIMatrix1, (MatrixD *)pIMatrix2, (MatrixDC *)pOMatrix);
        }
        else if(instanceof(pIMatrix2, MATRIXDC))
        {
            MatrixD_O2Mul3((MatrixD *)pIMatrix1, (MatrixDC *)pIMatrix2, (MatrixDC *)pOMatrix);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pOMatrix;
}

Matrix *_MatrixDC_O2Mul(Matrix *pIMatrix1, Matrix *pIMatrix2, Matrix *pOMatrix)
{
    if(instanceof(pOMatrix, MATRIXD))
    {
        if(instanceof(pIMatrix1, MATRIXD))
        {
            MatrixDC_O2Mul4((MatrixDC *)pIMatrix1, (MatrixD *)pIMatrix2, (MatrixD *)pOMatrix);
        }
        else if(instanceof(pIMatrix1, MATRIXDC))
        {
            MatrixDC_O2Mul3((MatrixDC *)pIMatrix1, (MatrixDC *)pIMatrix2, (MatrixD *)pOMatrix);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else if(instanceof(pOMatrix, MATRIXDC))
    {
        if(instanceof(pIMatrix1, MATRIXD))
        {
            MatrixDC_O2Mul2((MatrixDC *)pIMatrix1, (MatrixD *)pIMatrix2, (MatrixDC *)pOMatrix);
        }
        else if(instanceof(pIMatrix1, MATRIXDC))
        {
            MatrixDC_O2Mul((MatrixDC *)pIMatrix1, (MatrixDC *)pIMatrix2, (MatrixDC *)pOMatrix);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pOMatrix;
}

Matrix *Matrix_MulMatrix(Matrix *pIMatrix1, Matrix *pIMatrix2, Matrix *pOMatrix)
{
    return (*pIMatrix1->MulMatrix)(pIMatrix1, pIMatrix2, pOMatrix);
}

Matrix *_MatrixD_MulMatrix(Matrix *pIMatrix1, Matrix *pIMatrix2, Matrix *pOMatrix)
{
    if(instanceof(pOMatrix, MATRIXD))
    {
        if(instanceof(pIMatrix2, MATRIXD))
        {
            MatrixD_MulMatrixD((MatrixD *)pIMatrix1, (MatrixD *)pIMatrix2, (MatrixD *)pOMatrix);
        }
        else if(instanceof(pIMatrix2, MATRIXDC))
        {
            MatrixD_MulMatrixDC2((MatrixD *)pIMatrix1, (MatrixDC *)pIMatrix2, (MatrixD *)pOMatrix);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else if(instanceof(pOMatrix, MATRIXDC))
    {
        if(instanceof(pIMatrix2, MATRIXD))
        {
            MatrixD_MulMatrixD2((MatrixD *)pIMatrix1, (MatrixD *)pIMatrix2, (MatrixDC *)pOMatrix);
        }
        else if(instanceof(pIMatrix2, MATRIXDC))
        {
            MatrixD_MulMatrixDC((MatrixD *)pIMatrix1, (MatrixDC *)pIMatrix2, (MatrixDC *)pOMatrix);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pOMatrix;
}

Matrix *_MatrixDC_MulMatrix(Matrix *pIMatrix1, Matrix *pIMatrix2, Matrix *pOMatrix)
{
    if(instanceof(pOMatrix, MATRIXD))
    {
        if(instanceof(pIMatrix2, MATRIXD))
        {
            MatrixDC_MulMatrixD2((MatrixDC *)pIMatrix1, (MatrixD *)pIMatrix2, (MatrixD *)pOMatrix);
        }
        else if(instanceof(pIMatrix2, MATRIXDC))
        {
            MatrixDC_MulMatrixDC2((MatrixDC *)pIMatrix1, (MatrixDC *)pIMatrix2, (MatrixD *)pOMatrix);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else if(instanceof(pOMatrix, MATRIXDC))
    {
        if(instanceof(pIMatrix2, MATRIXD))
        {
            MatrixDC_MulMatrixD((MatrixDC *)pIMatrix1, (MatrixD *)pIMatrix2, (MatrixDC *)pOMatrix);
        }
        else if(instanceof(pIMatrix2, MATRIXDC))
        {
            MatrixDC_MulMatrixDC((MatrixDC *)pIMatrix1, (MatrixDC *)pIMatrix2, (MatrixDC *)pOMatrix);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pOMatrix;
}

Num *Matrix_MaxDiagElem(Matrix *pIMatrix, Num *pONum)
{
    return (*pIMatrix->MaxDiagElem)(pIMatrix, pONum);
}

Num *_MatrixD_MaxDiagElem(Matrix *pIMatrix, Num *pONum)
{
    if(instanceof(pONum, NUMD))
    {
        ((NumD *)pONum)->num = MatrixD_MaxDiagElem((MatrixD *)pIMatrix);
    }
    else if(instanceof(pONum, NUMDC))
    {
        ((NumDC *)pONum)->num = MatrixD_MaxDiagElem2((MatrixD *)pIMatrix);
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pONum;
}

Num *_MatrixDC_MaxDiagElem(Matrix *pIMatrix, Num *pONum)
{
    if(instanceof(pONum, NUMD))
    {
        ((NumD *)pONum)->num = MatrixDC_MaxDiagElem2((MatrixDC *)pIMatrix);
    }
    else if(instanceof(pONum, NUMDC))
    {
        ((NumDC *)pONum)->num = MatrixDC_MaxDiagElem((MatrixDC *)pIMatrix);
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pONum;
}

Vector *Matrix_LinearSolver(Matrix *pIMatrix, Vector *pIVector, Vector *pOVector)
{
    return (*pIMatrix->LinearSolver)(pIMatrix, pIVector, pOVector);
}

Vector *_MatrixD_LinearSolver(Matrix *pIMatrix, Vector *pIVector, Vector *pOVector)
{
    if(instanceof(pOVector, VECTORD))
    {
        if(instanceof(pIVector, VECTORD))
        {
            MatrixD_LinearSolver((MatrixD *)pIMatrix, (VectorD *)pIVector, (VectorD *)pOVector);
        }
        else if(instanceof(pIVector, VECTORD))
        {
            MatrixD_LinearSolver4((MatrixD *)pIMatrix, (VectorDC *)pIVector, (VectorD *)pOVector);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else if(instanceof(pOVector, VECTORDC))
    {
        if(instanceof(pIVector, VECTORD))
        {
            MatrixD_LinearSolver2((MatrixD *)pIMatrix, (VectorD *)pIVector, (VectorDC *)pOVector);
        }
        else if(instanceof(pIVector, VECTORD))
        {
            MatrixD_LinearSolver3((MatrixD *)pIMatrix, (VectorDC *)pIVector, (VectorDC *)pOVector);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pOVector;
}

Vector *_MatrixDC_LinearSolver(Matrix *pIMatrix, Vector *pIVector, Vector *pOVector)
{
    if(instanceof(pOVector, VECTORD))
    {
        if(instanceof(pIVector, VECTORD))
        {
            MatrixDC_LinearSolver3((MatrixDC *)pIMatrix, (VectorD *)pIVector, (VectorD *)pOVector);
        }
        else if(instanceof(pIVector, VECTORDC))
        {
            MatrixDC_LinearSolver4((MatrixDC *)pIMatrix, (VectorDC *)pIVector, (VectorD *)pOVector);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else if(instanceof(pOVector, VECTORDC))
    {
        if(instanceof(pIVector, VECTORD))
        {
            MatrixDC_LinearSolver2((MatrixDC *)pIMatrix, (VectorD *)pIVector, (VectorDC *)pOVector);
        }
        else if(instanceof(pIVector, VECTORDC))
        {
            MatrixDC_LinearSolver((MatrixDC *)pIMatrix, (VectorDC *)pIVector, (VectorDC *)pOVector);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pOVector;
}

void Matrix_EigenEquation(Matrix *pIMatrix, Vector *pOVector, Matrix *pOMatrix)
{
    (*pIMatrix->EigenEquation)(pIMatrix, pOVector, pOMatrix);
}

void _MatrixD_EigenEquation(Matrix *pIMatrix, Vector *pOVector, Matrix *pOMatrix)
{
    if(instanceof(pOMatrix, MATRIXD))
    {
        if(instanceof(pOVector, VECTORD))
        {
            MatrixD_EigenEquation((MatrixD *)pIMatrix, (VectorD *)pOVector, (MatrixD *)pOMatrix);
        }
        else if(instanceof(pOVector, VECTORDC))
        {
            MatrixD_EigenEquation2((MatrixD *)pIMatrix, (VectorD *)pOVector, (MatrixDC *)pOMatrix);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else if(instanceof(pOMatrix, MATRIXD))
    {
        if(instanceof(pOVector, VECTORD))
        {
            MatrixD_EigenEquation3((MatrixD *)pIMatrix, (VectorDC *)pOVector, (MatrixD *)pOMatrix);
        }
        else if(instanceof(pOVector, VECTORDC))
        {
            MatrixD_EigenEquation4((MatrixD *)pIMatrix, (VectorDC *)pOVector, (MatrixDC *)pOMatrix);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }
}

void _MatrixDC_EigenEquation(Matrix *pIMatrix, Vector *pOVector, Matrix *pOMatrix)
{
    if(instanceof(pOMatrix, MATRIXD))
    {
        if(instanceof(pOVector, VECTORD))
        {
            MatrixDC_EigenEquation2((MatrixDC *)pIMatrix, (VectorD *)pOVector, (MatrixD *)pOMatrix);
        }
        else if(instanceof(pOVector, VECTORDC))
        {
            MatrixDC_EigenEquation3((MatrixDC *)pIMatrix, (VectorDC *)pOVector, (MatrixD *)pOMatrix);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else if(instanceof(pOMatrix, MATRIXDC))
    {
        if(instanceof(pOVector, VECTORD))
        {
            MatrixDC_EigenEquation4((MatrixDC *)pIMatrix, (VectorD *)pOVector, (MatrixDC *)pOMatrix);
        }
        else if(instanceof(pOVector, VECTORDC))
        {
            MatrixDC_EigenEquation((MatrixDC *)pIMatrix, (VectorDC *)pOVector, (MatrixDC *)pOMatrix);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }
}
