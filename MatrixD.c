#include "Matrix.h"

extern Matrix *_MatrixD_Del(Matrix *pIMatrix);
extern Matrix *_MatrixD_UnWrap(Matrix *pIMatrix);
extern void _MatrixD_Show(Matrix *pIMatrix);
extern Matrix *_MatrixD_Set(Matrix *pIMatrix, Num *diag, Num *offDiag);
extern Matrix *_MatrixD_Copy(Matrix *pIMatrix, Matrix *pOMatrix);
extern Matrix *_MatrixD_T(Matrix *pIMatrix, Matrix *pOMatrix);
extern Num *_MatrixD_Det(Matrix *pIMatrix, Num *pNum);
extern Matrix *_MatrixD_Inv(Matrix *pIMatrix, Matrix *pOMatrix);
extern Vector *_MatrixD_MulVector(Matrix *pIMatrix, Vector *pIVector, Vector *pOVector);
extern Vector *_MatrixD_TMulVector(Matrix *pIMatrix, Vector *pIVector, Vector *pOVector);
extern Vector *_MatrixD_CMulVector(Matrix *pIMatrix, Vector *pIVector, Vector *pOVector);
extern Matrix *_MatrixD_TMulSelf(Matrix *pIMatrix, Matrix *pOMatrix);
extern Matrix *_MatrixD_O2Mul(Matrix *pIMatrix1, Matrix *pIMatrix2, Matrix *pOMatrix);
extern Matrix *_MatrixD_MulMatrix(Matrix *pIMatrix1, Matrix *pIMatrix2, Matrix *pOMatrix);
extern Num *_MatrixD_MaxDiagElem(Matrix *pIMatrix, Num *pONum);
extern Vector *_MatrixD_LinearSolver(Matrix *pILUMatrix, Vector *pIVector, Vector *pOVector);
extern void _MatrixD_EigenEquation(Matrix *pIMatrix, Vector *pOVector, Matrix *pOMatrix);

MatrixD *MatrixD_New(int row, int col)
{
    MatrixD *pNewMatrixD;

    pNewMatrixD = (MatrixD *)malloc(sizeof(MatrixD));
    if(pNewMatrixD == NULL)
    {
        fprintf(stderr, "Error: Can't create a new matrixD. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, __LINE__);
        return NULL;
    }
    pNewMatrixD->matrix = (double *)malloc(sizeof(double)*col*row);
    if(pNewMatrixD->matrix == NULL)
    {
        free(pNewMatrixD);
        fprintf(stderr, "Error: Can't create a new matrixD. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, __LINE__);
        return NULL;
    }
    pNewMatrixD->size[0] = row;
    pNewMatrixD->size[1] = col;
    pNewMatrixD->ld = col;
    pNewMatrixD->parent.type = MATRIXD;
    pNewMatrixD->parent.Del = _MatrixD_Del;
    pNewMatrixD->parent.UnWrap = _MatrixD_UnWrap;
    pNewMatrixD->parent.Show = _MatrixD_Show;
    pNewMatrixD->parent.Set = _MatrixD_Set;
    pNewMatrixD->parent.Copy = _MatrixD_Copy;
    pNewMatrixD->parent.T = _MatrixD_T;
    pNewMatrixD->parent.Det = _MatrixD_Det;
    pNewMatrixD->parent.Inv = _MatrixD_Inv;
    pNewMatrixD->parent.MulVector = _MatrixD_MulVector;
    pNewMatrixD->parent.TMulVector = _MatrixD_TMulVector;
    pNewMatrixD->parent.CMulVector = _MatrixD_CMulVector;
    pNewMatrixD->parent.TMulSelf = _MatrixD_TMulSelf;
    pNewMatrixD->parent.O2Mul = _MatrixD_O2Mul;
    pNewMatrixD->parent.MulMatrix = _MatrixD_MulMatrix;
    pNewMatrixD->parent.MaxDiagElem = _MatrixD_MaxDiagElem;
    pNewMatrixD->parent.LinearSolver = _MatrixD_LinearSolver;
    pNewMatrixD->parent.EigenEquation = _MatrixD_EigenEquation;

    return pNewMatrixD;
}

MatrixD *MatrixD_Del(MatrixD *pIMatrix)
{
    if(pIMatrix != NULL)
    {
        free(pIMatrix->matrix);
        free(pIMatrix);
    }

    return NULL;
}

MatrixD *MatrixD_Wrap(double *work, int row, int col, int ld)
{
    MatrixD *pNewMatrixD;

    pNewMatrixD = (MatrixD *)malloc(sizeof(MatrixD));
    if(pNewMatrixD == NULL)
    {
        fprintf(stderr, "Error: Can't create a new matrixD. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, __LINE__);
        return NULL;
    }
    pNewMatrixD->matrix = work;
    if(pNewMatrixD->matrix == NULL)
    {
        free(pNewMatrixD);
        fprintf(stderr, "Error: Can't create a new matrixD. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, __LINE__);
        return NULL;
    }
    pNewMatrixD->size[0] = row;
    pNewMatrixD->size[1] = col;
    pNewMatrixD->ld = ld;
    pNewMatrixD->parent.type = MATRIXD;
    pNewMatrixD->parent.Del = _MatrixD_Del;
    pNewMatrixD->parent.UnWrap = _MatrixD_UnWrap;
    pNewMatrixD->parent.Show = _MatrixD_Show;
    pNewMatrixD->parent.Set = _MatrixD_Set;
    pNewMatrixD->parent.Copy = _MatrixD_Copy;
    pNewMatrixD->parent.T = _MatrixD_T;
    pNewMatrixD->parent.Det = _MatrixD_Det;
    pNewMatrixD->parent.Inv = _MatrixD_Inv;
    pNewMatrixD->parent.MulVector = _MatrixD_MulVector;
    pNewMatrixD->parent.TMulVector = _MatrixD_TMulVector;
    pNewMatrixD->parent.CMulVector = _MatrixD_CMulVector;
    pNewMatrixD->parent.TMulSelf = _MatrixD_TMulSelf;
    pNewMatrixD->parent.O2Mul = _MatrixD_O2Mul;
    pNewMatrixD->parent.MulMatrix = _MatrixD_MulMatrix;
    pNewMatrixD->parent.MaxDiagElem = _MatrixD_MaxDiagElem;
    pNewMatrixD->parent.LinearSolver = _MatrixD_LinearSolver;
    pNewMatrixD->parent.EigenEquation = _MatrixD_EigenEquation;

    return pNewMatrixD;
}

MatrixD *MatrixD_UnWrap(MatrixD *pIMatrix)
{
    if(pIMatrix != NULL)
    {
        free(pIMatrix);
    }

    return NULL;
}

void MatrixD_Show(MatrixD *pIMatrix)
{
    register int i, j, row, col, ld;

    row = pIMatrix->size[0];
    col = pIMatrix->size[1];
    ld = pIMatrix->ld;
    printf("MatrixD_Show:\n");
    for(i=0;i<row;i++)
    {
        printf(" [");
        for(j=0;j<col;j++)
        {
            printf("%lf, ", *(pIMatrix->matrix+i*ld+j));
        }
        printf("]\n");
    }
    printf(" row:%d \n col:%d \n ld:%d \n\n", row, col, ld);
}

MatrixD *MatrixD_Set(MatrixD *pIMatrix, double diag, double offDiag)
{
    register int i, j, row, col, ld;
    double *matrix;

    row = pIMatrix->size[0];
    col = pIMatrix->size[1];
    ld = pIMatrix->ld;
    matrix = pIMatrix->matrix;

    for(i=0;i<row;i++)
    {
        for(j=0;j<col;j++)
        {
            if(i==j)
            {
                *(matrix+i*col+j) = diag;
            }
            else
            {
                *(matrix+i*col+j) = offDiag;
            }
        }
    }

    return pIMatrix;
}

#ifdef _MSC_VER
MatrixD *MatrixD_Set2(MatrixD *pIMatrix, double diag, _Dcomplex offDiag)
#else
MatrixD *MatrixD_Set2(MatrixD *pIMatrix, double diag, double complex offDiag)
#endif /*_MSC_VER*/
{
    register int i, j, row, col, ld;
    double *matrix;
    double tmp;

    row = pIMatrix->size[0];
    col = pIMatrix->size[1];
    ld = pIMatrix->ld;
    matrix = pIMatrix->matrix;
    #ifdef _MSC_VER
    tmp = offDiag._Val[0];
    #else
    tmp = creal(offDiag);
    #endif /*_MSC_VER*/

    for(i=0;i<row;i++)
    {
        for(j=0;j<col;j++)
        {
            if(i==j)
            {
                *(matrix+i*col+j) = diag;
            }
            else
            {
                *(matrix+i*col+j) = tmp;
            }
        }
    }

    return pIMatrix;
}

#ifdef _MSC_VER
MatrixD *MatrixD_Set3(MatrixD *pIMatrix, _Dcomplex diag, double offDiag)
#else
MatrixD *MatrixD_Set3(MatrixD *pIMatrix, double complex diag, double offDiag)
#endif /*_MSC_VER*/
{
    register int i, j, row, col, ld;
    double *matrix;
    double tmp;

    row = pIMatrix->size[0];
    col = pIMatrix->size[1];
    ld = pIMatrix->ld;
    matrix = pIMatrix->matrix;
    #ifdef _MSC_VER
    tmp = diag._Val[0];
    #else
    tmp = creal(diag);
    #endif /*_MSC_VER*/

    for(i=0;i<row;i++)
    {
        for(j=0;j<col;j++)
        {
            if(i==j)
            {
                *(matrix+i*col+j) = tmp;
            }
            else
            {
                *(matrix+i*col+j) = offDiag;
            }
        }
    }

    return pIMatrix;
}

#ifdef _MSC_VER
MatrixD *MatrixD_Set4(MatrixD *pIMatrix, _Dcomplex diag, _Dcomplex offDiag)
#else
MatrixD *MatrixD_Set4(MatrixD *pIMatrix, double complex diag, double complex offDiag)
#endif /*_MSC_VER*/
{
    register int i, j, row, col, ld;
    double *matrix;
    double tmp1, tmp2;

    row = pIMatrix->size[0];
    col = pIMatrix->size[1];
    ld = pIMatrix->ld;
    matrix = pIMatrix->matrix;
    #ifdef _MSC_VER
    tmp1 = diag._Val[0];
    #else
    tmp1 = creal(diag);
    #endif /*_MSC_VER*/
    #ifdef _MSC_VER
    tmp2 = offDiag._Val[0];
    #else
    tmp2 = creal(offDiag);
    #endif /*_MSC_VER*/

    for(i=0;i<row;i++)
    {
        for(j=0;j<col;j++)
        {
            if(i==j)
            {
                *(matrix+i*col+j) = tmp1;
            }
            else
            {
                *(matrix+i*col+j) = tmp2;
            }
        }
    }

    return pIMatrix;
}

MatrixD *MatrixD_Copy(MatrixD *pIMatrix, MatrixD *pOMatrix)
{
    register int i, j, ld1, ld2, row, col;
    register double *iMatrix, *oMatrix, *ioMatrix, *iiMatrix;

    if(pIMatrix->size[0] != pOMatrix->size[0] || pIMatrix->size[1] != pOMatrix->size[1])
    {
        fprintf(stderr, "Error: pIMatrix and pOMatrix have incorrect size. File: %s Func: %s Line: %d\n", \
            __FILE__, __FUNCTION__, __LINE__);
        return pOMatrix;
    }

    ld1 = pIMatrix->ld;
    ld2 = pOMatrix->ld;
    row = pIMatrix->size[0];
    col = pIMatrix->size[1];
    iMatrix = pIMatrix->matrix;
    oMatrix = pOMatrix->matrix;

    for(i=0;i<row;i++)
    {
        iiMatrix = iMatrix + i * ld1;
        ioMatrix = oMatrix + i * ld2;
        for(j=0;j<col;j++)
        {
            *(ioMatrix+j) = *(iiMatrix+j);
        }
    }

    return pOMatrix;
}

MatrixDC *MatrixD_Copy2(MatrixD *pIMatrix, MatrixDC *pOMatrix)
{
    register int i, j, ld1, ld2, row, col;
    register double *iMatrix, *iiMatrix;
    #ifdef _MSC_VER
    register _Dcomplex *oMatrix, *ioMatrix;
    #else
    register double complex *oMatrix, *ioMatrix;
    #endif /*_MSC_VER*/

    if(pIMatrix->size[0] != pOMatrix->size[0] || pIMatrix->size[1] != pOMatrix->size[1])
    {
        fprintf(stderr, "Error: pIMatrix and pOMatrix have incorrect size. File: %s Func: %s Line: %d\n", \
            __FILE__, __FUNCTION__, __LINE__);
        return pOMatrix;
    }

    ld1 = pIMatrix->ld;
    ld2 = pOMatrix->ld;
    row = pIMatrix->size[0];
    col = pIMatrix->size[1];
    iMatrix = pIMatrix->matrix;
    oMatrix = pOMatrix->matrix;

    for(i=0;i<row;i++)
    {
        iiMatrix = iMatrix + i * ld1;
        ioMatrix = oMatrix + i * ld2;
        for(j=0;j<col;j++)
        {
            #ifdef _MSC_VER
            (ioMatrix+j)->_Val[0] = *(iiMatrix+j);
            (ioMatrix+j)->_Val[1] = 0;
            #else
            *(ioMatrix+j) = *(iiMatrix+j);
            #endif /*_MSC_VER*/
        }
    }

    return pOMatrix;
}

MatrixD *MatrixD_T(MatrixD *pIMatrix, MatrixD *pOMatrix)
{
    register int i, j, row, col, ld1, ld2;
    register double *iMatrix, *oMatrix, *joMatrix;

    row = pIMatrix->size[0];
    col = pIMatrix->size[1];
    ld1 = pIMatrix->ld;
    ld2 = pOMatrix->ld;
    iMatrix = pIMatrix->matrix;
    oMatrix = pOMatrix->matrix;

    if(row != pOMatrix->size[1] || col != pOMatrix->size[0])
    {
        fprintf(stderr, "Error: pIMatrix and pOMatrix have incorrect size. File: %s Func: %s Line: %d\n", \
            __FILE__, __FUNCTION__, __LINE__);
        return pOMatrix;
    }

    for(j=col; j-- > 0;)
    {
        joMatrix = oMatrix + j*ld2;
        for(i=row; i-- >0;)
        {
            *(joMatrix+i) = *(iMatrix+i*ld1+j);
        }
    }

    return pOMatrix;
}

MatrixDC *MatrixD_T2(MatrixD *pIMatrix, MatrixDC *pOMatrix)
{
    register int i, j, row, col, ld1, ld2;
    register double *iMatrix;
    #ifdef _MSC_VER
    register _Dcomplex *oMatrix, *joMatrix;
    #else
    register double complex *oMatrix, *joMatrix;
    #endif /*_MSC_VER*/

    row = pIMatrix->size[0];
    col = pIMatrix->size[1];
    ld1 = pIMatrix->ld;
    ld2 = pOMatrix->ld;
    iMatrix = pIMatrix->matrix;
    oMatrix = pOMatrix->matrix;

    if(row != pOMatrix->size[1] || col != pOMatrix->size[0])
    {
        fprintf(stderr, "Error: pIMatrix and pOMatrix have incorrect size. File: %s Func: %s Line: %d\n", \
            __FILE__, __FUNCTION__, __LINE__);
        return pOMatrix;
    }

    for(j=col; j-- > 0;)
    {
        joMatrix = oMatrix + j*ld2;
        for(i=row; i-- >0;)
        {
            #ifdef _MSC_VER
            (joMatrix+i)->_Val[0] = *(iMatrix+i*ld1+j);
            (joMatrix+i)->_Val[1] = 0;
            #else
            *(joMatrix+i) = *(iMatrix+i*ld1+j);
            #endif /*_MSC_VER*/
        }
    }

    return pOMatrix;
}

double MatrixD_Det(MatrixD *pIMatrix)
{
    int row, col, info;
    int *ipiv;
    register int i;
    double det;
    MatrixD *pTmpMatrix;

    row = pIMatrix->size[0];
    col = pIMatrix->size[1];
    ipiv = (int *)malloc(sizeof(int) * ((row > col) ? col : row));
    det = 1.0;
    pTmpMatrix = MatrixD_New(row,col);
    pTmpMatrix = MatrixD_Copy(pIMatrix, pTmpMatrix);

    dgetrf_(&row, &col, pTmpMatrix->matrix, &row, ipiv, &info);

    for(i=0;i<row;i++)
    {
        det *= *(pTmpMatrix->matrix+i*col+i);
    }

    MatrixD_Del(pTmpMatrix);
    free(ipiv);

    return det;
}

#ifdef _MSC_VER
_Dcomplex MatrixD_Det2(MatrixD *pIMatrix)
#else
double complex MatrixD_Det2(MatrixD *pIMatrix)
#endif /*_MSC_VER*/
{
    int row, col, info;
    int *ipiv;
    register int i;
    #ifdef _MSC_VER
    _Dcomplex det;
    #else
    double complex det;
    #endif /*_MSC_VER*/
    MatrixD *pTmpMatrix;

    row = pIMatrix->size[0];
    col = pIMatrix->size[1];
    ipiv = (int *)malloc(sizeof(int) * ((row > col) ? col : row));
    #ifdef _MSC_VER
    det._Val[0] = 1;
    det._Val[1] = 0;
    #else
    det = 1.0;
    #endif /*_MSC_VER*/
    pTmpMatrix = MatrixD_New(row,col);
    pTmpMatrix = MatrixD_Copy(pIMatrix, pTmpMatrix);

    dgetrf_(&row, &col, pTmpMatrix->matrix, &row, ipiv, &info);

    for(i=0;i<row;i++)
    {
        #ifdef _MSC_VER
        det._Val[0] *= *(pTmpMatrix->matrix+i*col+i);
        #else
        det *= *(pTmpMatrix->matrix+i*col+i);
        #endif /*_MSC_VER*/
    }

    MatrixD_Del(pTmpMatrix);
    free(ipiv);

    return det;
}

MatrixD *MatrixD_Inv(MatrixD *pIMatrix, MatrixD *pOMatrix)
{
    int row, col, info, ld2;
    int *ipiv;
    double *work;

    row = pIMatrix->size[0];
    col = pIMatrix->size[1];
    ld2 = pOMatrix->ld;
    ipiv = (int *)malloc(sizeof(int) * ((row > col) ? col : row));
    work = (double *)malloc(sizeof(double) * col * row);

    pOMatrix = MatrixD_T(pIMatrix, pOMatrix);
    dgetrf_(&row, &col, pOMatrix->matrix, &ld2, ipiv, &info);
    dgetri_(&row, pOMatrix->matrix, &ld2, ipiv, work, &row, &info);

    free(ipiv);
    free(work);
    return pOMatrix;
}

MatrixDC *MatrixD_Inv2(MatrixD *pIMatrix, MatrixDC *pOMatrix)
{
    int row, col, info, ld2;
    int *ipiv;
    #ifdef _MSC_VER
    _Dcomplex *work;
    #else
    double complex *work;
    #endif /*_MSC_VER*/

    row = pIMatrix->size[0];
    col = pIMatrix->size[1];
    ld2 = pOMatrix->ld;
    ipiv = (int *)malloc(sizeof(int) * ((row > col) ? col : row));
    #ifdef _MSC_VER
    work = (_Dcomplex *)malloc(sizeof(_Dcomplex) * col *row);
    #else
    work = (double complex *)malloc(sizeof(double complex) * col * row);
    #endif /*_MSC_VER*/

    pOMatrix = MatrixD_T2(pIMatrix, pOMatrix);
    zgetrf_(&row, &col, pOMatrix->matrix, &ld2, ipiv, &info);
    zgetri_(&row, pOMatrix->matrix, &ld2, ipiv, work, &row, &info);

    free(ipiv);
    free(work);
    return pOMatrix;
}

VectorD *MatrixD_MulVectorD(MatrixD *pIMatrix, VectorD *pIVector, VectorD *pOVector)
{
    register int i, l, iSize, oSize, ld, inc1, inc2;
    register double *iVector, *oVector, *iMatrix, *iiMatrix;

    iSize = pIMatrix->size[1];
    oSize = pIMatrix->size[0];
    ld = pIMatrix->ld;
    inc1 = pIVector->inc;
    inc2 = pOVector->inc;
    iVector = pIVector->vector;
    iMatrix = pIMatrix->matrix;
    oVector = pOVector->vector;

    if(iSize != pIVector->size || oSize != pOVector->size)
    {
        fprintf(stderr, "Error: pIMatrix and pIVector or pOVector have incorrect size. File: %s Func: %s Line: %d\n", \
            __FILE__, __FUNCTION__, __LINE__);
        return pOVector;
    }

    for(i=0; i<oSize; ++i)
    {
        *(oVector+i*inc2) = 0.0;
        iiMatrix = iMatrix + i*ld;
        for(l=0; l<iSize; ++l)
        {
            *(oVector+i*inc2) += *(iiMatrix+l) * *(iVector+l*inc1);
        }
    }
    return pOVector;
}

VectorDC *MatrixD_MulVectorD2(MatrixD *pIMatrix, VectorD *pIVector, VectorDC *pOVector)
{
    register int i, l, iSize, oSize, ld, inc1, inc2;
    register double *iVector, *iMatrix, *iiMatrix;
    #ifdef _MSC_VER
    register _Dcomplex *oVector;
    #else
    register double complex *oVector;
    #endif /*_MSC_VER*/

    iSize = pIMatrix->size[1];
    oSize = pIMatrix->size[0];
    ld = pIMatrix->ld;
    inc1 = pIVector->inc;
    inc2 = pOVector->inc;
    iVector = pIVector->vector;
    iMatrix = pIMatrix->matrix;
    oVector = pOVector->vector;

    if(iSize != pIVector->size || oSize != pOVector->size)
    {
        fprintf(stderr, "Error: pIMatrix and pIVector or pOVector have incorrect size. File: %s Func: %s Line: %d\n", \
            __FILE__, __FUNCTION__, __LINE__);
        return pOVector;
    }

    for(i=0; i<oSize; ++i)
    {
        #ifdef _MSC_VER
        (oVector+i*inc2)->_Val[0] = 0.0;
        (oVector+i*inc2)->_Val[1] = 0.0;
        #else
        *(oVector+i*inc2) = 0.0;
        #endif /*_MSC_VER*/

        iiMatrix = iMatrix + i*ld;
        for(l=0; l<iSize; ++l)
        {
            #ifdef _MSC_VER
            (oVector+i*inc2)->_Val[0] += *(iiMatrix+l) * *(iVector+l*inc1);
            #else
            *(oVector+i*inc2) += *(iiMatrix+l) * *(iVector+l*inc1);
            #endif /*_MSC_VER*/
        }
    }

    return pOVector;
}

VectorDC *MatrixD_MulVectorDC(MatrixD *pIMatrix, VectorDC *pIVector, VectorDC *pOVector)
{
    register int i, l, iSize, oSize, ld, inc1, inc2;
    register double *iMatrix, *iiMatrix;
    #ifdef _MSC_VER
    register _Dcomplex *iVector, *oVector;
    #else
    register double complex *iVector, *oVector;
    #endif /*_MSC_VER*/

    iSize = pIMatrix->size[1];
    oSize = pIMatrix->size[0];
    ld = pIMatrix->ld;
    inc1 = pIVector->inc;
    inc2 = pOVector->inc;
    iVector = pIVector->vector;
    iMatrix = pIMatrix->matrix;
    oVector = pOVector->vector;

    if(iSize != pIVector->size || oSize != pOVector->size)
    {
        fprintf(stderr, "Error: pIMatrix and pIVector or pOVector have incorrect size. File: %s Func: %s Line: %d\n", \
            __FILE__, __FUNCTION__, __LINE__);
        return pOVector;
    }

    for(i=0; i<oSize; ++i)
    {
        #ifdef _MSC_VER
        (oVector+i*inc2)->_Val[0] = 0.0;
        (oVector+i*inc2)->_Val[1] = 0.0;
        #else
        *(oVector+i*inc2) = 0.0;
        #endif /*_MSC_VER*/

        iiMatrix = iMatrix + i*ld;
        for(l=0; l<iSize; ++l)
        {
            #ifdef _MSC_VER
            (oVector+i*inc2)->_Val[0] += *(iiMatrix+l) * (iVector+l*inc1)->_Val[0];
            (oVector+i*inc2)->_Val[1] += *(iiMatrix+l) * (iVector+l*inc1)->_Val[1];
            #else
            *(oVector+i*inc2) += *(iiMatrix+l) * *(iVector+l*inc1);
            #endif /*_MSC_VER*/
        }
    }

    return pOVector;
}

VectorD *MatrixD_MulVectorDC2(MatrixD *pIMatrix, VectorDC *pIVector, VectorD *pOVector)
{
    register int i, l, iSize, oSize, ld, inc1, inc2;
    register double *iMatrix, *iiMatrix, *oVector;
    #ifdef _MSC_VER
    register _Dcomplex *iVector;
    #else
    register double complex *iVector;
    #endif /*_MSC_VER*/

    iSize = pIMatrix->size[1];
    oSize = pIMatrix->size[0];
    ld = pIMatrix->ld;
    inc1 = pIVector->inc;
    inc2 = pOVector->inc;
    iVector = pIVector->vector;
    iMatrix = pIMatrix->matrix;
    oVector = pOVector->vector;

    if(iSize != pIVector->size || oSize != pOVector->size)
    {
        fprintf(stderr, "Error: pIMatrix and pIVector or pOVector have incorrect size. File: %s Func: %s Line: %d\n", \
            __FILE__, __FUNCTION__, __LINE__);
        return pOVector;
    }

    for(i=0; i<oSize; ++i)
    {
        *(oVector+i*inc2) = 0.0;

        iiMatrix = iMatrix + i*ld;
        for(l=0; l<iSize; ++l)
        {
            #ifdef _MSC_VER
            *(oVector+i*inc2) += *(iiMatrix+l) * (iVector+l*inc1)->_Val[0];
            #else
            *(oVector+i*inc2) += *(iiMatrix+l) * creal(*(iVector+l*inc1));
            #endif /*_MSC_VER*/
        }
    }

    return pOVector;
}

VectorD *MatrixD_TMulVectorD(MatrixD *pIMatrix, VectorD *pIVector, VectorD *pOVector)
{
    register int i, j, row, col, ld, inc1, inc2;
    register double *iMatrix, *iVector, *oVector, *iiMatrix;
    register double tmp;

    row = pIMatrix->size[0];
    col = pIMatrix->size[1];
    ld = pIMatrix->ld;
    inc1 = pIVector->inc;
    inc2 = pOVector->inc;

    iMatrix = pIMatrix->matrix;
    iVector = pIVector->vector;
    oVector = pOVector->vector;

    if(row != pIVector->size || col != pOVector->size)
    {
        fprintf(stderr, "Error: pIMatrix and pIVector or pOVector have incorrect size. File: %s Func: %s Line: %d\n", \
            __FILE__, __FUNCTION__, __LINE__);
        return pOVector;
    }

    for(j=col; j-- >0;)
    {
        *(oVector+j*inc2) = 0;
    }

    for(i=row; i-- >0;)
    {
        iiMatrix = iMatrix+i*ld;
        tmp = *(iVector+i*inc1);
        for(j=col; j-- >0;)
        {
            *(oVector+j*inc2) +=  *(iiMatrix+j) * tmp;
        }
    }

    return pOVector;
}

VectorDC *MatrixD_TMulVectorD2(MatrixD *pIMatrix, VectorD *pIVector, VectorDC *pOVector)
{
    register int i, j, row, col, ld, inc1, inc2;
    register double *iMatrix, *iVector, *iiMatrix;
    register double tmp;
    #ifdef _MSC_VER
    register _Dcomplex *oVector;
    #else
    register double complex *oVector;
    #endif /*_MSC_VER*/

    row = pIMatrix->size[0];
    col = pIMatrix->size[1];
    ld = pIMatrix->ld;
    inc1 = pIVector->inc;
    inc2 = pOVector->inc;

    iMatrix = pIMatrix->matrix;
    iVector = pIVector->vector;
    oVector = pOVector->vector;

    if(row != pIVector->size || col != pOVector->size)
    {
        fprintf(stderr, "Error: pIMatrix and pIVector or pOVector have incorrect size. File: %s Func: %s Line: %d\n", \
            __FILE__, __FUNCTION__, __LINE__);
        return pOVector;
    }

    for(j=col; j-- >0;)
    {
        #ifdef _MSC_VER
        (oVector+j*inc2)->_Val[0] = 0;
        (oVector+j*inc2)->_Val[1] = 0;
        #else
        *(oVector+j*inc2) = 0;
        #endif /*_MSC_VER*/
    }

    for(i=row; i-- >0;)
    {
        iiMatrix = iMatrix+i*ld;
        tmp = *(iVector+i*inc1);
        for(j=col; j-- >0;)
        {
            #ifdef _MSC_VER
            (oVector+j*inc2)->_Val[0] +=  *(iiMatrix+j) * tmp;
            #else
            *(oVector+j*inc2) += *(iiMatrix+j) * tmp;
            #endif /*_MSC_VER*/
        }
    }

    return pOVector;
}

VectorDC *MatrixD_TMulVectorDC(MatrixD *pIMatrix, VectorDC *pIVector, VectorDC *pOVector)
{
    register int i, j, row, col, ld, inc1, inc2;
    register double *iMatrix, *iiMatrix;
    #ifdef _MSC_VER
    register _Dcomplex *iVector, *oVector;
    register _Dcomplex tmp;
    #else
    register double complex *iVector, *oVector;
    register double complex tmp;
    #endif /*_MSC_VER*/

    row = pIMatrix->size[0];
    col = pIMatrix->size[1];
    ld = pIMatrix->ld;
    inc1 = pIVector->inc;
    inc2 = pOVector->inc;

    iMatrix = pIMatrix->matrix;
    iVector = pIVector->vector;
    oVector = pOVector->vector;

    if(row != pIVector->size || col != pOVector->size)
    {
        fprintf(stderr, "Error: pIMatrix and pIVector or pOVector have incorrect size. File: %s Func: %s Line: %d\n", \
            __FILE__, __FUNCTION__, __LINE__);
        return pOVector;
    }

    for(j=col; j-- >0;)
    {
        #ifdef _MSC_VER
        (oVector+j*inc2)->_Val[0] = 0;
        (oVector+j*inc2)->_Val[1] = 0;
        #else
        *(oVector+j*inc2) = 0;
        #endif /*_MSC_VER*/
    }

    for(i=row; i-- >0;)
    {
        iiMatrix = iMatrix+i*ld;
        tmp = *(iVector+i*inc1);
        for(j=col; j-- >0;)
        {
            #ifdef _MSC_VER
            (oVector+j*inc2)->_Val[0] += *(iiMatrix+j) * tmp._Val[0];
            (oVector+j*inc2)->_Val[1] += *(iiMatrix+j) * tmp._Val[1];
            #else
            *(oVector+j*inc2) += *(iiMatrix+j) * tmp;
            #endif /*_MSC_VER*/
        }
    }

    return pOVector;
}

VectorD *MatrixD_TMulVectorDC2(MatrixD *pIMatrix, VectorDC *pIVector, VectorD *pOVector)
{
    register int i, j, row, col, ld, inc1, inc2;
    register double *iMatrix, *iiMatrix;
    double *oVector;
    #ifdef _MSC_VER
    register _Dcomplex *iVector;
    register _Dcomplex tmp;
    #else
    register double complex *iVector;
    register double complex tmp;
    #endif /*_MSC_VER*/

    row = pIMatrix->size[0];
    col = pIMatrix->size[1];
    ld = pIMatrix->ld;
    inc1 = pIVector->inc;
    inc2 = pOVector->inc;

    iMatrix = pIMatrix->matrix;
    iVector = pIVector->vector;
    oVector = pOVector->vector;

    if(row != pIVector->size || col != pOVector->size)
    {
        fprintf(stderr, "Error: pIMatrix and pIVector or pOVector have incorrect size. File: %s Func: %s Line: %d\n", \
            __FILE__, __FUNCTION__, __LINE__);
        return pOVector;
    }

    for(j=col; j-- >0;)
    {
        *(oVector+j*inc2) = 0;
    }

    for(i=row; i-- >0;)
    {
        iiMatrix = iMatrix+i*ld;
        tmp = *(iVector+i*inc1);
        for(j=col; j-- >0;)
        {
            #ifdef _MSC_VER
            *(oVector+j*inc2) += *(iiMatrix+j) * tmp._Val[0];
            #else
            *(oVector+j*inc2) += *(iiMatrix+j) * creal(tmp);
            #endif /*_MSC_VER*/
        }
    }

    return pOVector;
}

VectorD *MatrixD_CMulVectorD(MatrixD *pIMatrix, VectorD *pIVector, VectorD *pOVector)
{
    register int i, j, row, col, ld, inc1, inc2;
    register double *iMatrix, *iVector, *oVector, *iiMatrix;
    register double tmp;

    row = pIMatrix->size[0];
    col = pIMatrix->size[1];
    ld = pIMatrix->ld;
    inc1 = pIVector->inc;
    inc2 = pOVector->inc;
    iMatrix = pIMatrix->matrix;
    iVector = pIVector->vector;
    oVector = pOVector->vector;

    if(row != pIVector->size || col != pOVector->size)
    {
        fprintf(stderr, "Error: pIMatrix and pIVector or pOVector have incorrect size. File: %s Func: %s Line: %d\n", \
            __FILE__, __FUNCTION__, __LINE__);
        return pOVector;
    }

    for(j=col; j-- >0;)
    {
        *(oVector+j*inc2) = 0;
    }

    for(i=row; i-- >0;)
    {
        iiMatrix = iMatrix+i*ld;
        tmp = *(iVector+i*inc1);
        for(j=col; j-- >0;)
        {
            *(oVector+j) +=  *(iiMatrix+j) * tmp;
        }
    }

    return pOVector;
}

VectorDC *MatrixD_CMulVectorD2(MatrixD *pIMatrix, VectorD *pIVector, VectorDC *pOVector)
{
    register int i, j, row, col, ld, inc1, inc2;
    register double *iMatrix, *iVector, *iiMatrix;
    register double tmp;
    #ifdef _MSC_VER
    register _Dcomplex *oVector;
    #else
    register double complex *oVector;
    #endif /*_MSC_VER*/

    row = pIMatrix->size[0];
    col = pIMatrix->size[1];
    ld = pIMatrix->ld;
    inc1 = pIVector->inc;
    inc2 = pOVector->inc;
    iMatrix = pIMatrix->matrix;
    iVector = pIVector->vector;
    oVector = pOVector->vector;

    if(row != pIVector->size || col != pOVector->size)
    {
        fprintf(stderr, "Error: pIMatrix and pIVector or pOVector have incorrect size. File: %s Func: %s Line: %d\n", \
            __FILE__, __FUNCTION__, __LINE__);
        return pOVector;
    }

    for(j=col; j-- >0;)
    {
        #ifdef _MSC_VER
        (oVector+j*inc2)->_Val[0] = 0;
        (oVector+j*inc2)->_Val[1] = 0;
        #else
        *(oVector+j*inc2) = 0;
        #endif /*_MSC_VER*/
    }

    for(i=row; i-- >0;)
    {
        iiMatrix = iMatrix+i*ld;
        tmp = *(iVector+i*inc1);
        for(j=col; j-- >0;)
        {
            #ifdef _MSC_VER
            (oVector+j)->_Val[0] += *(iiMatrix+j) * tmp;
            #else
            *(oVector+j) +=  *(iiMatrix+j) * tmp;
            #endif /*_MSC_VER*/
        }
    }

    return pOVector;
}

VectorDC *MatrixD_CMulVectorDC(MatrixD *pIMatrix, VectorDC *pIVector, VectorDC *pOVector)
{
    register int i, j, row, col, ld, inc1, inc2;
    register double *iMatrix, *iiMatrix;
    #ifdef _MSC_VER
    register _Dcomplex tmp;
    register _Dcomplex *iVector, *oVector;
    #else
    register double complex tmp;
    register double complex *iVector, *oVector;
    #endif /*_MSC_VER*/

    row = pIMatrix->size[0];
    col = pIMatrix->size[1];
    ld = pIMatrix->ld;
    inc1 = pIVector->inc;
    inc2 = pOVector->inc;
    iMatrix = pIMatrix->matrix;
    iVector = pIVector->vector;
    oVector = pOVector->vector;

    if(row != pIVector->size || col != pOVector->size)
    {
        fprintf(stderr, "Error: pIMatrix and pIVector or pOVector have incorrect size. File: %s Func: %s Line: %d\n", \
            __FILE__, __FUNCTION__, __LINE__);
        return pOVector;
    }

    for(j=col; j-- >0;)
    {
        #ifdef _MSC_VER
        (oVector+j*inc2)->_Val[0] = 0;
        (oVector+j*inc2)->_Val[1] = 0;
        #else
        *(oVector+j*inc2) = 0;
        #endif /*_MSC_VER*/
    }

    for(i=row; i-- >0;)
    {
        iiMatrix = iMatrix+i*ld;
        tmp = *(iVector+i*inc1);
        for(j=col; j-- >0;)
        {
            #ifdef _MSC_VER
            (oVector+j)->_Val[0] += *(iiMatrix+j) * tmp._Val[0];
            (oVector+j)->_Val[1] += *(iiMatrix+j) * tmp._Val[1];
            #else
            *(oVector+j) +=  *(iiMatrix+j) * tmp;
            #endif /*_MSC_VER*/
        }
    }

    return pOVector;
}

VectorD *MatrixD_CMulVectorDC2(MatrixD *pIMatrix, VectorDC *pIVector, VectorD *pOVector)
{
    register int i, j, row, col, ld, inc1, inc2;
    register double *iMatrix, *iiMatrix;
    #ifdef _MSC_VER
    register _Dcomplex tmp;
    register _Dcomplex *iVector;
    #else
    register double complex tmp;
    register double complex *iVector;
    #endif /*_MSC_VER*/
    register double *oVector;

    row = pIMatrix->size[0];
    col = pIMatrix->size[1];
    ld = pIMatrix->ld;
    inc1 = pIVector->inc;
    inc2 = pOVector->inc;
    iMatrix = pIMatrix->matrix;
    iVector = pIVector->vector;
    oVector = pOVector->vector;

    if(row != pIVector->size || col != pOVector->size)
    {
        fprintf(stderr, "Error: pIMatrix and pIVector or pOVector have incorrect size. File: %s Func: %s Line: %d\n", \
            __FILE__, __FUNCTION__, __LINE__);
        return pOVector;
    }

    for(j=col; j-- >0;)
    {
        *(oVector+j*inc2) = 0;
    }

    for(i=row; i-- >0;)
    {
        iiMatrix = iMatrix+i*ld;
        tmp = *(iVector+i*inc1);
        for(j=col; j-- >0;)
        {
            #ifdef _MSC_VER
            *(oVector+j) += *(iiMatrix+j) * tmp._Val[0];
            #else
            *(oVector+j) += *(iiMatrix+j) * creal(tmp);
            #endif /*_MSC_VER*/
        }
    }

    return pOVector;
}

MatrixD *MatrixD_TMulSelf(MatrixD *pIMatrix, MatrixD *pOMatrix)
{
    register int i, j, l, row, col, ld1, ld2;
    register double *iMatrix, *oMatrix, *liMatrix, *ioMatrix;
    register double tmp;

    row = pIMatrix->size[0];
    col = pIMatrix->size[1];
    ld1 = pIMatrix->ld;
    ld2 = pOMatrix->ld;
    iMatrix = pIMatrix->matrix;
    oMatrix = pOMatrix->matrix;

    if(pOMatrix->size[0]!=pOMatrix->size[1] || col!=pOMatrix->size[1])
    {
        fprintf(stderr, "Error: pI/OMatrix isn't a square matrix or has incorrect size. File: %s Func: %s Line: %d\n", \
            __FILE__, __FUNCTION__, __LINE__);
        return pOMatrix;
    }

    for(i = col; i-- > 0;)
    {
        for(j = col; j-- > 0;)
        {
            *(oMatrix+i*ld2+j) = 0;
        }
    }

    for(l=row; l-- > 0;)
    {
        liMatrix = iMatrix + l*ld1;
        for(i=col; i-- > 0;)
        {
            ioMatrix = oMatrix + i*ld2;
            tmp = *(liMatrix+i);
            for(j=i+1; j-- > 0;)
            {
                *(ioMatrix+j) += tmp * *(liMatrix+j);
            }
        }
    }

    for(i = col;i-- > 0;)
    {
        for(j=i+1;j<col;++j)
        {
            *(oMatrix+i*ld2+j) = *(oMatrix+j*ld2+i);
        }
    }
    return pOMatrix;
}

MatrixDC *MatrixD_TMulSelf2(MatrixD *pIMatrix, MatrixDC *pOMatrix)
{
    register int i, j, l, row, col, ld1, ld2;
    register double *iMatrix, *liMatrix;
    register double tmp;
    #ifdef _MSC_VER
    register _Dcomplex *oMatrix, *ioMatrix;
    #else
    register double complex *oMatrix, *ioMatrix;
    #endif /*_MSC_VER*/

    row = pIMatrix->size[0];
    col = pIMatrix->size[1];
    ld1 = pIMatrix->ld;
    ld2 = pOMatrix->ld;
    iMatrix = pIMatrix->matrix;
    oMatrix = pOMatrix->matrix;

    if(pOMatrix->size[0]!=pOMatrix->size[1] || col!=pOMatrix->size[1])
    {
        fprintf(stderr, "Error: pOMatrix isn't a square matrix or has incorrect size. File: %s Func: %s Line: %d\n", \
            __FILE__, __FUNCTION__, __LINE__);
        return pOMatrix;
    }

    for(i = col; i-- > 0;)
    {
        for(j = col; j-- > 0;)
        {
            #ifdef _MSC_VER
            (oMatrix+i*ld2+j)->_Val[0] = 0;
            (oMatrix+i*ld2+j)->_Val[1] = 0;
            #else
            *(oMatrix+i*ld2+j) = 0;
            #endif /*_MSC_VER*/
        }
    }

    for(l=row; l-- > 0;)
    {
        liMatrix = iMatrix + l*ld1;
        for(i=col; i-- > 0;)
        {
            ioMatrix = oMatrix + i*ld2;
            tmp = *(liMatrix+i);
            for(j=i+1; j-- > 0;)
            {
                #ifdef _MSC_VER
                (ioMatrix+j)->_Val[0] += tmp * *(liMatrix+j);
                #else
                *(ioMatrix+j) += tmp * *(liMatrix+j);
                #endif /*_MSC_VER*/
            }
        }
    }

    for(i = col;i-- > 0;)
    {
        for(j=i+1;j<col;++j)
        {
            *(oMatrix+i*ld2+j) = *(oMatrix+j*ld2+i);
        }
    }

    return pOMatrix;
}

MatrixD *MatrixD_O2Mul(MatrixD *pIMatrix1, MatrixD *pIMatrix2, MatrixD *pOMatrix)
{
    int ld1, ld2, ld3;
    register double *iMatrix1, *iMatrix2, *oMatrix;
    double m11, m12, m13, m14, m21, m22, m23, m24;

    ld1 = pIMatrix1->ld;
    ld2 = pIMatrix2->ld;
    ld3 = pOMatrix->ld;
    iMatrix1 = pIMatrix1->matrix;
    iMatrix2 = pIMatrix2->matrix;
    oMatrix = pOMatrix->matrix;

    m11 = *(iMatrix1);
    m12 = *(iMatrix1+1);
    m13 = *(iMatrix1+ld1);
    m14 = *(iMatrix1+ld1+1);
    m21 = *(iMatrix2);
    m22 = *(iMatrix2+1);
    m23 = *(iMatrix2+ld2);
    m24 = *(iMatrix2+ld2+1);

    *(oMatrix) = m11 * m21 + m12 * m23;
    *(oMatrix+1) = m11 * m22 + m12 * m24;
    *(oMatrix+ld3) = m13 * m21 + m14 * m23;
    *(oMatrix+ld3+1) = m13 * m22 + m14 * m24;

    return pOMatrix;
}

MatrixDC *MatrixD_O2Mul2(MatrixD *pIMatrix1, MatrixD *pIMatrix2, MatrixDC *pOMatrix)
{
    int ld1, ld2, ld3;
    register double *iMatrix1, *iMatrix2;
    double m11, m12, m13, m14, m21, m22, m23, m24;
    #ifdef _MSC_VER
    register _Dcomplex *oMatrix;
    #else
    register double complex *oMatrix;
    #endif /*_MSC_VER*/

    ld1 = pIMatrix1->ld;
    ld2 = pIMatrix2->ld;
    ld3 = pOMatrix->ld;
    iMatrix1 = pIMatrix1->matrix;
    iMatrix2 = pIMatrix2->matrix;
    oMatrix = pOMatrix->matrix;

    m11 = *(iMatrix1);
    m12 = *(iMatrix1+1);
    m13 = *(iMatrix1+ld1);
    m14 = *(iMatrix1+ld1+1);
    m21 = *(iMatrix2);
    m22 = *(iMatrix2+1);
    m23 = *(iMatrix2+ld2);
    m24 = *(iMatrix2+ld2+1);

    #ifdef _MSC_VER
    (oMatrix)->_Val[0] = m11 * m21 + m12 * m23;
    (oMatrix)->_Val[1] = 0;
    (oMatrix+1)->_Val[0] = m11 * m22 + m12 * m24;
    (oMatrix+1)->_Val[1] = 0;
    (oMatrix+ld3)->_Val[0] = m13 * m21 + m14 * m23;
    (oMatrix+ld3)->_Val[1] = 0;
    (oMatrix+ld3+1)->_Val[0] = m13 * m22 + m14 * m24;
    (oMatrix+ld3+1)->_Val[1] = 0;
    #else
    *(oMatrix) = m11 * m21 + m12 * m23;
    *(oMatrix+1) = m11 * m22 + m12 * m24;
    *(oMatrix+ld3) = m13 * m21 + m14 * m23;
    *(oMatrix+ld3+1) = m13 * m22 + m14 * m24;
    #endif /*_MSC_VER*/

    return pOMatrix;
}

MatrixDC *MatrixD_O2Mul3(MatrixD *pIMatrix1, MatrixDC *pIMatrix2, MatrixDC *pOMatrix)
{
    int ld1, ld2, ld3;
    register double *iMatrix1;
    double m11, m12, m13, m14;
    #ifdef _MSC_VER
    register _Dcomplex *iMatrix2, *oMatrix;
    _Dcomplex m21, m22, m23, m24;
    #else
    register double complex *iMatrix2, *oMatrix;
    double complex m21, m22, m23, m24;
    #endif /*_MSC_VER*/

    ld1 = pIMatrix1->ld;
    ld2 = pIMatrix2->ld;
    ld3 = pOMatrix->ld;
    iMatrix1 = pIMatrix1->matrix;
    iMatrix2 = pIMatrix2->matrix;
    oMatrix = pOMatrix->matrix;

    m11 = *(iMatrix1);
    m12 = *(iMatrix1+1);
    m13 = *(iMatrix1+ld1);
    m14 = *(iMatrix1+ld1+1);
    m21 = *(iMatrix2);
    m22 = *(iMatrix2+1);
    m23 = *(iMatrix2+ld2);
    m24 = *(iMatrix2+ld2+1);

    #ifdef _MSC_VER
    (oMatrix)->_Val[0] = m11 * m21._Val[0] + m12 * m23._Val[0];
    (oMatrix)->_Val[1] = m11 * m21._Val[1] + m12 * m23._Val[1];
    (oMatrix+1)->_Val[0] = m11 * m22._Val[0] + m12 * m24._Val[0];
    (oMatrix+1)->_Val[1] = m11 * m22._Val[1] + m12 * m24._Val[1];
    (oMatrix+ld3)->_Val[0] = m13 * m21._Val[0] + m14 * m23._Val[0];
    (oMatrix+ld3)->_Val[1] = m13 * m21._Val[1] + m14 * m23._Val[1];
    (oMatrix+ld3+1)->_Val[0] = m13 * m22._Val[0] + m14 * m24._Val[0];
    (oMatrix+ld3+1)->_Val[1] = m13 * m22._Val[1] + m14 * m24._Val[1];
    #else
    *(oMatrix) = m11 * m21 + m12 * m23;
    *(oMatrix+1) = m11 * m22 + m12 * m24;
    *(oMatrix+ld3) = m13 * m21 + m14 * m23;
    *(oMatrix+ld3+1) = m13 * m22 + m14 * m24;
    #endif /*_MSC_VER*/

    return pOMatrix;
}

MatrixD *MatrixD_O2Mul4(MatrixD *pIMatrix1, MatrixDC *pIMatrix2, MatrixD *pOMatrix)
{
    int ld1, ld2, ld3;
    register double *iMatrix1, *oMatrix;
    double m11, m12, m13, m14;
    #ifdef _MSC_VER
    register _Dcomplex *iMatrix2;
    _Dcomplex m21, m22, m23, m24;
    #else
    register double complex *iMatrix2;
    double complex m21, m22, m23, m24;
    #endif /*_MSC_VER*/

    ld1 = pIMatrix1->ld;
    ld2 = pIMatrix2->ld;
    ld3 = pOMatrix->ld;
    iMatrix1 = pIMatrix1->matrix;
    iMatrix2 = pIMatrix2->matrix;
    oMatrix = pOMatrix->matrix;

    m11 = *(iMatrix1);
    m12 = *(iMatrix1+1);
    m13 = *(iMatrix1+ld1);
    m14 = *(iMatrix1+ld1+1);
    m21 = *(iMatrix2);
    m22 = *(iMatrix2+1);
    m23 = *(iMatrix2+ld2);
    m24 = *(iMatrix2+ld2+1);

    #ifdef _MSC_VER
    *(oMatrix) = m11 * m21._Val[0] + m12 * m23._Val[0];
    *(oMatrix+1) = m11 * m22._Val[0] + m12 * m24._Val[0];
    *(oMatrix+ld3) = m13 * m21._Val[0] + m14 * m23._Val[0];
    *(oMatrix+ld3+1) = m13 * m22._Val[0] + m14 * m24._Val[0];
    #else
    *(oMatrix) = m11 * creal(m21) + m12 * creal(m23);
    *(oMatrix+1) = m11 * creal(m22) + m12 * creal(m24);
    *(oMatrix+ld3) = m13 * creal(m21) + m14 * creal(m23);
    *(oMatrix+ld3+1) = m13 * creal(m22) + m14 * creal(m24);
    #endif /*_MSC_VER*/

    return pOMatrix;
}

MatrixD *MatrixD_MulMatrixD(MatrixD *pIMatrix1, MatrixD *pIMatrix2, MatrixD *pOMatrix)
{
    double alpha, beta;
    int m, n, k, lda, ldb, ldc;

    alpha = 1;
    beta = 0;
    lda = pIMatrix1->ld;
    ldb = pIMatrix2->ld;
    ldc = pOMatrix->ld;
    m = pIMatrix1->size[1];
    n = pIMatrix2->size[0];
    k = pIMatrix1->size[0];

    dgemm_("T","T",&m,&n,&k,&alpha,pIMatrix1->matrix,&lda,pIMatrix2->matrix,&ldb,&beta,pOMatrix->matrix,&ldc);

    return pOMatrix;
}

MatrixDC *MatrixD_MulMatrixD2(MatrixD *pIMatrix1, MatrixD *pIMatrix2, MatrixDC *pOMatrix)
{
    double alpha, beta;
    int m, n, k, lda, ldb, ldc;
    MatrixD *pTmpMatrix;

    alpha = 1;
    beta = 0;
    lda = pIMatrix1->ld;
    ldb = pIMatrix2->ld;
    m = pIMatrix1->size[1];
    n = pIMatrix2->size[0];
    k = pIMatrix1->size[0];
    pTmpMatrix = MatrixD_New(pOMatrix->size[0],pOMatrix->size[1]);
    ldc = pTmpMatrix->ld;

    dgemm_("T","T",&m,&n,&k,&alpha,pIMatrix1->matrix,&lda,pIMatrix2->matrix,&ldb,&beta,pTmpMatrix->matrix,&ldc);
    pOMatrix = MatrixD_Copy2(pTmpMatrix, pOMatrix);
    pTmpMatrix = MatrixD_Del(pTmpMatrix);

    return pOMatrix;
}

MatrixDC *MatrixD_MulMatrixDC(MatrixD *pIMatrix1, MatrixDC *pIMatrix2, MatrixDC *pOMatrix)
{
    #ifdef _MSC_VER
    _Dcomplex alpha, beta;
    #else
    double complex alpha, beta;
    #endif /*_MSC_VER*/
    int m, n, k, lda, ldb, ldc;
    MatrixDC *pTmpMatrix;

    #ifdef _MSC_VER
    alpha = _Cbuild(1,0);
    beta = _Cbuild(0,0);
    #else
    alpha = 1;
    beta = 0;
    #endif /*_MSC_VER*/
    ldb = pIMatrix2->ld;
    ldc = pOMatrix->ld;
    m = pIMatrix1->size[1];
    n = pIMatrix2->size[0];
    k = pIMatrix1->size[0];
    pTmpMatrix = MatrixDC_New(k,m);
    pTmpMatrix = MatrixD_Copy2(pIMatrix1, pTmpMatrix);
    lda = pTmpMatrix->ld;

    zgemm_("T","T",&m,&n,&k,&alpha,pTmpMatrix->matrix,&lda,pIMatrix2->matrix,&ldb,&beta,pOMatrix->matrix,&ldc);
    pTmpMatrix = MatrixDC_Del(pTmpMatrix);

    return pOMatrix;
}

MatrixD *MatrixD_MulMatrixDC2(MatrixD *pIMatrix1, MatrixDC *pIMatrix2, MatrixD *pOMatrix)
{
    double alpha, beta;
    int m, n, k, lda, ldb, ldc;
    MatrixD *pTmpMatrix;

    alpha = 1;
    beta = 0;
    lda = pIMatrix1->ld;
    ldc = pOMatrix->ld;
    m = pIMatrix1->size[1];
    n = pIMatrix2->size[0];
    k = pIMatrix1->size[0];
    pTmpMatrix = MatrixD_New(n,k);
    MatrixDC_Copy2(pIMatrix2, pTmpMatrix);
    ldb = pTmpMatrix->ld;

    dgemm_("T","T",&m,&n,&k,&alpha,pIMatrix1->matrix,&lda,pTmpMatrix->matrix,&ldb,&beta,pOMatrix->matrix,&ldc);
    pTmpMatrix = MatrixD_Del(pTmpMatrix);

    return pOMatrix;
}

double MatrixD_MaxDiagElem(MatrixD *pIMatrix)
{
    register int i, size, ld;
    register double *iMatrix, max;

    size = pIMatrix->size[0];
    ld = pIMatrix->ld;
    iMatrix = pIMatrix->matrix;
    max = -DBL_MAX;

    if(size != pIMatrix->size[1])
    {
        fprintf(stderr, "Error: pIMatrix isn't a square matrix. File: %s Func: %s Line: %d\n", \
            __FILE__, __FUNCTION__, __LINE__);
        return max;
    }

    for(i=0;i<size;++i)
    {
        max = (max > *(iMatrix+i*ld+i)) ? max : *(iMatrix+i*ld+i);
    }

    return max;
}

#ifdef _MSC_VER
_Dcomplex MatrixD_MaxDiagElem2(MatrixD *pIMatrix)
#else
double complex MatrixD_MaxDiagElem2(MatrixD *pIMatrix)
#endif /*_MSC_VER*/
{
    register int i, size, ld;
    register double *iMatrix, max;
    #ifdef _MSC_VER
    _Dcomplex max2;
    #else
    double complex max2;
    #endif /*_MSC_VER*/

    size = pIMatrix->size[0];
    ld = pIMatrix->ld;
    iMatrix = pIMatrix->matrix;
    max = -DBL_MAX;
    #ifdef _MSC_VER
    max2._Val[0] = max;
    max2._Val[1] = 0;
    #else
    max2 = max;
    #endif /*_MSC_VER*/

    if(size != pIMatrix->size[1])
    {
        fprintf(stderr, "Error: pIMatrix isn't a square matrix. File: %s Func: %s Line: %d\n", \
            __FILE__, __FUNCTION__, __LINE__);
        return max2;
    }

    for(i=0;i<size;++i)
    {
        max = (max > *(iMatrix+i*ld+i)) ? max : *(iMatrix+i*ld+i);
    }

    #ifdef _MSC_VER
    max2._Val[0] = max;
    max2._Val[1] = 0;
    #else
    max2 = max;
    #endif /*_MSC_VER*/

    return max2;
}

VectorD *MatrixD_LinearSolver(MatrixD *pIMatrix, VectorD *pIVector, VectorD *pOVector)
{
    int row, col, info, nrhs;
    int *ipiv;
    MatrixD *pTmpMatrix;

    nrhs = 1;
    row = pIMatrix->size[0];
    col = pIMatrix->size[1];
    ipiv = (int *)malloc(sizeof(int) * ((row > col) ? col : row));
    pTmpMatrix = MatrixD_New(col, row);

    pTmpMatrix = MatrixD_T(pIMatrix, pTmpMatrix);
    pOVector = VectorD_Copy(pIVector, pOVector);
    dgetrf_(&row, &col, pTmpMatrix->matrix, &row, ipiv, &info);
    dgetrs_("N", &row, &nrhs, pTmpMatrix->matrix, &row, ipiv, pOVector->vector, &row, &info);

    free(ipiv);
    pTmpMatrix = MatrixD_Del(pTmpMatrix);

    return pOVector;
}

VectorDC *MatrixD_LinearSolver2(MatrixD *pIMatrix, VectorD *pIVector, VectorDC *pOVector)
{
    int row, col, info, nrhs;
    int *ipiv;
    VectorD *pTmpVector;
    MatrixD *pTmpMatrix;

    nrhs = 1;
    row = pIMatrix->size[0];
    col = pIMatrix->size[1];
    ipiv = (int *)malloc(sizeof(int) * ((row > col) ? col : row));
    pTmpVector = VectorD_New(row);
    pTmpMatrix = MatrixD_New(col, row);

    pTmpMatrix = MatrixD_T(pIMatrix, pTmpMatrix);
    pTmpVector = VectorD_Copy(pIVector, pTmpVector);
    dgetrf_(&row, &col, pTmpMatrix->matrix, &row, ipiv, &info);
    dgetrs_("N", &row, &nrhs, pTmpMatrix->matrix, &row, ipiv, pTmpVector->vector, &row, &info);
    pOVector = VectorD_Copy2(pTmpVector, pOVector);

    free(ipiv);
    pTmpVector = VectorD_Del(pTmpVector);
    pTmpMatrix = MatrixD_Del(pTmpMatrix);

    return pOVector;
}

VectorDC *MatrixD_LinearSolver3(MatrixD *pIMatrix, VectorDC *pIVector, VectorDC *pOVector)
{
    int row, col, info, nrhs;
    int *ipiv;
    MatrixDC *pTmpMatrix;

    nrhs = 1;
    row = pIMatrix->size[0];
    col = pIMatrix->size[1];
    ipiv = (int *)malloc(sizeof(int) * ((row > col) ? col : row));
    pTmpMatrix = MatrixDC_New(col, row);

    MatrixD_T2(pIMatrix, pTmpMatrix);
    VectorDC_Copy(pIVector, pOVector);
    zgetrf_(&row, &col, pTmpMatrix->matrix, &row, ipiv, &info);
    zgetrs_("N", &row, &nrhs, pTmpMatrix->matrix, &row, ipiv, pOVector->vector, &row, &info);

    free(ipiv);
    pTmpMatrix = MatrixDC_Del(pTmpMatrix);

    return pOVector;
}

VectorD *MatrixD_LinearSolver4(MatrixD *pIMatrix, VectorDC *pIVector, VectorD *pOVector)
{
    int row, col, info, nrhs;
    int *ipiv;
    MatrixDC *pTmpMatrix;
    VectorDC *pTmpVector;

    nrhs = 1;
    row = pIMatrix->size[0];
    col = pIMatrix->size[1];
    ipiv = (int *)malloc(sizeof(int) * ((row > col) ? col : row));
    pTmpMatrix = MatrixDC_New(col, row);
    pTmpVector = VectorDC_New(row);

    MatrixD_T2(pIMatrix, pTmpMatrix);
    VectorDC_Copy(pIVector, pTmpVector);
    zgetrf_(&row, &col, pTmpMatrix->matrix, &row, ipiv, &info);
    zgetrs_("N", &row, &nrhs, pTmpMatrix->matrix, &row, ipiv, pTmpVector->vector, &row, &info);
    VectorDC_Copy2(pTmpVector, pOVector);

    free(ipiv);
    MatrixDC_Del(pTmpMatrix);
    VectorDC_Del(pTmpVector);

    return pOVector;
}

void MatrixD_EigenEquation(MatrixD *pIMatrix, VectorD *pOVector, MatrixD *pOMatrix)
{
    int n, ldvl, ldvr, lwork, info;
    double *rwork, *work;
    double wkopt;
    MatrixD *pTmpMatrix,*pTmpMatrix2;

    n = pIMatrix->size[0];
    ldvl = n;
    ldvr = pOMatrix->ld;
    lwork = -1;
    rwork = (double *)malloc(sizeof(double) * 2 * n);
    pTmpMatrix = MatrixD_New(n,n);
    pTmpMatrix2 = MatrixD_New(n,n);

    pTmpMatrix = MatrixD_T(pIMatrix, pTmpMatrix);
    dgeev_("N","V",&n,pTmpMatrix->matrix,&n,pOVector->vector,NULL,&ldvl,pTmpMatrix2->matrix,&ldvr, &wkopt, &lwork, rwork, &info);
    lwork = (int) wkopt;
    work = (double *)malloc(sizeof(double) * lwork);
    dgeev_("N","V",&n,pTmpMatrix->matrix,&n,pOVector->vector,NULL,&ldvl,pTmpMatrix2->matrix,&ldvr, work, &lwork, rwork, &info);
    pOMatrix = MatrixD_T(pTmpMatrix2,pOMatrix);

    pTmpMatrix = MatrixD_Del(pTmpMatrix);
    pTmpMatrix2 = MatrixD_Del(pTmpMatrix2);
    free(work);
    free(rwork);
}

void MatrixD_EigenEquation2(MatrixD *pIMatrix, VectorD *pOVector, MatrixDC *pOMatrix)
{
    int n, ldvl, ldvr, lwork, info;
    double *rwork, *work;
    double wkopt;
    MatrixD *pTmpMatrix,*pTmpMatrix2;

    n = pIMatrix->size[0];
    ldvl = n;
    ldvr = pOMatrix->ld;
    lwork = -1;
    rwork = (double *)malloc(sizeof(double) * 2 * n);
    pTmpMatrix = MatrixD_New(n,n);
    pTmpMatrix2 = MatrixD_New(n,n);

    pTmpMatrix = MatrixD_T(pIMatrix, pTmpMatrix);
    dgeev_("N","V",&n,pTmpMatrix->matrix,&n,pOVector->vector,NULL,&ldvl,pTmpMatrix2->matrix,&ldvr, &wkopt, &lwork, rwork, &info);
    lwork = (int) wkopt;
    work = (double *)malloc(sizeof(double) * lwork);
    dgeev_("N","V",&n,pTmpMatrix->matrix,&n,pOVector->vector,NULL,&ldvl,pTmpMatrix2->matrix,&ldvr, work, &lwork, rwork, &info);
    pOMatrix = MatrixD_T2(pTmpMatrix2,pOMatrix);

    pTmpMatrix = MatrixD_Del(pTmpMatrix);
    pTmpMatrix2 = MatrixD_Del(pTmpMatrix2);
    free(work);
    free(rwork);
}

void MatrixD_EigenEquation3(MatrixD *pIMatrix, VectorDC *pOVector, MatrixD *pOMatrix)
{
    int n, ldvl, ldvr, lwork, info;
    double *rwork, *work;
    double wkopt;
    VectorD *pTmpVector;
    MatrixD *pTmpMatrix,*pTmpMatrix2;

    n = pIMatrix->size[0];
    ldvl = n;
    ldvr = pOMatrix->ld;
    lwork = -1;
    rwork = (double *)malloc(sizeof(double) * 2 * n);
    pTmpVector = VectorD_New(n);
    pTmpMatrix = MatrixD_New(n,n);
    pTmpMatrix2 = MatrixD_New(n,n);

    pTmpMatrix = MatrixD_T(pIMatrix, pTmpMatrix);
    dgeev_("N","V",&n,pTmpMatrix->matrix,&n,pTmpVector->vector,NULL,&ldvl,pTmpMatrix2->matrix,&ldvr, &wkopt, &lwork, rwork, &info);
    lwork = (int) wkopt;
    work = (double *)malloc(sizeof(double) * lwork);
    dgeev_("N","V",&n,pTmpMatrix->matrix,&n,pTmpVector->vector,NULL,&ldvl,pTmpMatrix2->matrix,&ldvr, work, &lwork, rwork, &info);
    pOMatrix = MatrixD_T(pTmpMatrix2,pOMatrix);
    pOVector = VectorD_Copy2(pTmpVector, pOVector);

    pTmpVector = VectorD_Del(pTmpVector);
    pTmpMatrix = MatrixD_Del(pTmpMatrix);
    pTmpMatrix2 = MatrixD_Del(pTmpMatrix2);
    free(work);
    free(rwork);
}

void MatrixD_EigenEquation4(MatrixD *pIMatrix, VectorDC *pOVector, MatrixDC *pOMatrix)
{
    int n, ldvl, ldvr, lwork, info;
    double *rwork, *work;
    double wkopt;
    VectorD *pTmpVector;
    MatrixD *pTmpMatrix,*pTmpMatrix2;

    n = pIMatrix->size[0];
    ldvl = n;
    ldvr = pOMatrix->ld;
    lwork = -1;
    rwork = (double *)malloc(sizeof(double) * 2 * n);
    pTmpVector = VectorD_New(n);
    pTmpMatrix = MatrixD_New(n,n);
    pTmpMatrix2 = MatrixD_New(n,n);

    pTmpMatrix = MatrixD_T(pIMatrix, pTmpMatrix);
    dgeev_("N","V",&n,pTmpMatrix->matrix,&n,pTmpVector->vector,NULL,&ldvl,pTmpMatrix2->matrix,&ldvr, &wkopt, &lwork, rwork, &info);
    lwork = (int) wkopt;
    work = (double *)malloc(sizeof(double) * lwork);
    dgeev_("N","V",&n,pTmpMatrix->matrix,&n,pTmpVector->vector,NULL,&ldvl,pTmpMatrix2->matrix,&ldvr, work, &lwork, rwork, &info);
    pOMatrix = MatrixD_T2(pTmpMatrix2,pOMatrix);
    pOVector = VectorD_Copy2(pTmpVector, pOVector);

    pTmpVector = VectorD_Del(pTmpVector);
    pTmpMatrix = MatrixD_Del(pTmpMatrix);
    pTmpMatrix2 = MatrixD_Del(pTmpMatrix2);
    free(work);
    free(rwork);
}
