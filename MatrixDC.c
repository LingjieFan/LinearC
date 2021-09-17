#include "Matrix.h"

extern Matrix *_MatrixDC_Del(Matrix *pIMatrix);
extern Matrix *_MatrixDC_UnWrap(Matrix *pIMatrix);
extern void _MatrixDC_Show(Matrix *pIMatrix);
extern Matrix *_MatrixDC_Set(Matrix *pIMatrix, Num *diag, Num *offDiag);
extern Matrix *_MatrixDC_Copy(Matrix *pIMatrix, Matrix *pOMatrix);
extern Matrix *_MatrixDC_T(Matrix *pIMatrix, Matrix *pOMatrix);
extern Num *_MatrixDC_Det(Matrix *pIMatrix, Num *pNum);
extern Matrix *_MatrixDC_Inv(Matrix *pIMatrix, Matrix *pOMatrix);
extern Vector *_MatrixDC_MulVector(Matrix *pIMatrix, Vector *pIVector, Vector *pOVector);
extern Vector *_MatrixDC_TMulVector(Matrix *pIMatrix, Vector *pIVector, Vector *pOVector);
extern Vector *_MatrixDC_CMulVector(Matrix *pIMatrix, Vector *pIVector, Vector *pOVector);
extern Matrix *_MatrixDC_TMulSelf(Matrix *pIMatrix, Matrix *pOMatrix);
extern Matrix *_MatrixDC_O2Mul(Matrix *pIMatrix1, Matrix *pIMatrix2, Matrix *pOMatrix);
extern Matrix *_MatrixDC_MulMatrix(Matrix *pIMatrix1, Matrix *pIMatrix2, Matrix *pOMatrix);
extern Num *_MatrixDC_MaxDiagElem(Matrix *pIMatrix, Num *pONum);
extern Vector *_MatrixDC_LinearSolver(Matrix *pILUMatrix, Vector *pIVector, Vector *pOVector);
extern void _MatrixDC_EigenEquation(Matrix *pIMatrix, Vector *pOVector, Matrix *pOMatrix);

MatrixDC *MatrixDC_New(int row, int col)
{
    MatrixDC *pNewMatrixDC;

    pNewMatrixDC = (MatrixDC *)malloc(sizeof(MatrixDC));
    if(pNewMatrixDC == NULL)
    {
        fprintf(stderr, "Error: Can't create a new matrixDC. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, __LINE__);
        return NULL;
    }
    #ifdef _MSC_VER
    pNewMatrixDC->matrix = (_Dcomplex *)malloc(sizeof(_Dcomplex)*col*row);
    #else
    pNewMatrixDC->matrix = (double complex *)malloc(sizeof(double complex)*col*row);
    #endif /*_MSC_VER*/
    if(pNewMatrixDC->matrix == NULL)
    {
        free(pNewMatrixDC);
        fprintf(stderr, "Error: Can't create a new matrixDC. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, __LINE__);
        return NULL;
    }
    pNewMatrixDC->size[0] = row;
    pNewMatrixDC->size[1] = col;
    pNewMatrixDC->ld = col;
    ((Object *)pNewMatrixDC)->type = MATRIXDC;
    pNewMatrixDC->parent.Del = _MatrixDC_Del;
    pNewMatrixDC->parent.UnWrap = _MatrixDC_UnWrap;
    pNewMatrixDC->parent.Show = _MatrixDC_Show;
    pNewMatrixDC->parent.Set = _MatrixDC_Set;
    pNewMatrixDC->parent.Copy = _MatrixDC_Copy;
    pNewMatrixDC->parent.T = _MatrixDC_T;
    pNewMatrixDC->parent.Det = _MatrixDC_Det;
    pNewMatrixDC->parent.Inv = _MatrixDC_Inv;
    pNewMatrixDC->parent.MulVector = _MatrixDC_MulVector;
    pNewMatrixDC->parent.TMulVector = _MatrixDC_TMulVector;
    pNewMatrixDC->parent.CMulVector = _MatrixDC_CMulVector;
    pNewMatrixDC->parent.TMulSelf = _MatrixDC_TMulSelf;
    pNewMatrixDC->parent.O2Mul = _MatrixDC_O2Mul;
    pNewMatrixDC->parent.MulMatrix = _MatrixDC_MulMatrix;
    pNewMatrixDC->parent.MaxDiagElem = _MatrixDC_MaxDiagElem;
    pNewMatrixDC->parent.LinearSolver = _MatrixDC_LinearSolver;
    pNewMatrixDC->parent.EigenEquation = _MatrixDC_EigenEquation;

    return pNewMatrixDC;
}

MatrixDC *MatrixDC_Del(MatrixDC *pIMatrix)
{
    if (pIMatrix != NULL)
    {
        free(pIMatrix->matrix);
        free(pIMatrix);
    }

    return NULL;
}

#ifdef _MSC_VER
MatrixDC *MatrixDC_Wrap(_Dcomplex *work, int row, int col, int ld)
#else
MatrixDC *MatrixDC_Wrap(double complex *work, int row, int col, int ld)
#endif /*_MSC_VER*/
{
    MatrixDC *pNewMatrixDC;

    pNewMatrixDC = (MatrixDC *)malloc(sizeof(MatrixDC));
    if(pNewMatrixDC == NULL)
    {
        fprintf(stderr, "Error: Can't create a new matrixDC. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, __LINE__);
        return NULL;
    }
    pNewMatrixDC->matrix = work;
    if(pNewMatrixDC->matrix == NULL)
    {
        free(pNewMatrixDC);
        fprintf(stderr, "Error: Can't create a new matrixDC. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, __LINE__);
        return NULL;
    }
    pNewMatrixDC->size[0] = row;
    pNewMatrixDC->size[1] = col;
    pNewMatrixDC->ld = ld;
    ((Object *)pNewMatrixDC)->type = MATRIXDC;
    pNewMatrixDC->parent.Del = _MatrixDC_Del;
    pNewMatrixDC->parent.UnWrap = _MatrixDC_UnWrap;
    pNewMatrixDC->parent.Show = _MatrixDC_Show;
    pNewMatrixDC->parent.Set = _MatrixDC_Set;
    pNewMatrixDC->parent.Copy = _MatrixDC_Copy;
    pNewMatrixDC->parent.T = _MatrixDC_T;
    pNewMatrixDC->parent.Det = _MatrixDC_Det;
    pNewMatrixDC->parent.Inv = _MatrixDC_Inv;
    pNewMatrixDC->parent.MulVector = _MatrixDC_MulVector;
    pNewMatrixDC->parent.TMulVector = _MatrixDC_TMulVector;
    pNewMatrixDC->parent.CMulVector = _MatrixDC_CMulVector;
    pNewMatrixDC->parent.TMulSelf = _MatrixDC_TMulSelf;
    pNewMatrixDC->parent.O2Mul = _MatrixDC_O2Mul;
    pNewMatrixDC->parent.MulMatrix = _MatrixDC_MulMatrix;
    pNewMatrixDC->parent.MaxDiagElem = _MatrixDC_MaxDiagElem;
    pNewMatrixDC->parent.LinearSolver = _MatrixDC_LinearSolver;
    pNewMatrixDC->parent.EigenEquation = _MatrixDC_EigenEquation;

    return pNewMatrixDC;
}

MatrixDC *MatrixDC_UnWrap(MatrixDC *pIMatrix)
{
    if (pIMatrix != NULL)
    {
        free(pIMatrix);
    }

    return NULL;
}

void MatrixDC_Show(MatrixDC *pIMatrix)
{
    register int i, j, row, col, ld;

    row = pIMatrix->size[0];
    col = pIMatrix->size[1];
    ld = pIMatrix->ld;
    printf("MatrixDC_Show:\n");
    for(i=0;i<row;i++)
    {
        printf(" [");
        for(j=0;j<col;j++)
        {
            printf("%f+%f j, ", creal(*(pIMatrix->matrix+i*ld+j)),cimag(*(pIMatrix->matrix+i*ld+j)));
        }
        printf("]\n");
    }
    printf(" row:%d \n col:%d \n ld:%d \n\n", row, col, ld);
}

#ifdef _MSC_VER
MatrixDC *MatrixDC_Set(MatrixDC *pIMatrix, _Dcomplex diag, _Dcomplex offDiag)
#else
MatrixDC *MatrixDC_Set(MatrixDC *pIMatrix, double complex diag, double complex offDiag)
#endif /*_MSC_VER*/
{
    register int i, j, row, col, ld;
    #ifdef _MSC_VER
    register _Dcomplex *matrix;
    #else
    register double complex *matrix;
    #endif /*_MSC_VER*/

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

MatrixDC *MatrixDC_Set2(MatrixDC *pIMatrix, double diag, double offDiag)
{
    register int i, j, row, col, ld;
    #ifdef _MSC_VER
    register _Dcomplex *matrix;
    #else
    register double complex *matrix;
    #endif /*_MSC_VER*/

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
                #ifdef _MSC_VER
                (matrix+i*col+j)->_Val[0] = diag;
                (matrix+i*col+j)->_Val[1] = 0;
                #else
                *(matrix+i*col+j) = diag;
                #endif /*_MSC_VER*/
            }
            else
            {
                #ifdef _MSC_VER
                (matrix+i*col+j)->_Val[0] = offDiag;
                (matrix+i*col+j)->_Val[1] = 0;
                #else
                *(matrix+i*col+j) = offDiag;
                #endif /*_MSC_VER*/
            }
        }
    }

    return pIMatrix;
}

#ifdef _MSC_VER
MatrixDC *MatrixDC_Set3(MatrixDC *pIMatrix, _Dcomplex diag, double offDiag)
#else
MatrixDC *MatrixDC_Set3(MatrixDC *pIMatrix, double complex diag, double offDiag)
#endif /*_MSC_VER*/
{
    register int i, j, row, col, ld;
    #ifdef _MSC_VER
    register _Dcomplex *matrix;
    #else
    register double complex *matrix;
    #endif /*_MSC_VER*/

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
                #ifdef _MSC_VER
                (matrix+i*col+j)->_Val[0] = offDiag;
                (matrix+i*col+j)->_Val[1] = 0;
                #else
                *(matrix+i*col+j) = offDiag;
                #endif /*_MSC_VER*/
            }
        }
    }

    return pIMatrix;
}

#ifdef _MSC_VER
MatrixDC *MatrixDC_Set4(MatrixDC *pIMatrix, double diag, _Dcomplex offDiag)
#else
MatrixDC *MatrixDC_Set4(MatrixDC *pIMatrix, double diag, double complex offDiag)
#endif /*_MSC_VER*/
{
    register int i, j, row, col, ld;
    #ifdef _MSC_VER
    register _Dcomplex *matrix;
    #else
    register double complex *matrix;
    #endif /*_MSC_VER*/

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
                #ifdef _MSC_VER
                (matrix+i*col+j)->_Val[0] = diag;
                (matrix+i*col+j)->_Val[1] = 0;
                #else
                *(matrix+i*col+j) = diag;
                #endif /*_MSC_VER*/
            }
            else
            {
                *(matrix+i*col+j) = offDiag;
            }
        }
    }

    return pIMatrix;
}

MatrixDC *MatrixDC_Copy(MatrixDC *pIMatrix, MatrixDC *pOMatrix)
{
    register int i, j, ld1, ld2, row, col;
    #ifdef _MSC_VER
    register _Dcomplex *iMatrix, *oMatrix, *ioMatrix, *iiMatrix;
    #else
    register double complex *iMatrix, *oMatrix, *ioMatrix, *iiMatrix;
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
            *(ioMatrix+j) = *(iiMatrix+j);
        }
    }

    return pOMatrix;
}

MatrixD *MatrixDC_Copy2(MatrixDC *pIMatrix, MatrixD *pOMatrix)
{
    register int i, j, ld1, ld2, row, col;
    #ifdef _MSC_VER
    register _Dcomplex *iMatrix, *iiMatrix;
    #else
    register double complex *iMatrix, *iiMatrix;
    #endif /*_MSC_VER*/
    double *oMatrix, *ioMatrix;

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
            *(ioMatrix+j) = (iiMatrix+j)->_Val[0];
            #else
            *(ioMatrix+j) = creal(*(iiMatrix+j));
            #endif /*_MSC_VER*/
        }
    }

    return pOMatrix;
}

MatrixDC *MatrixDC_T(MatrixDC *pIMatrix, MatrixDC *pOMatrix)
{
    register int i, j, row, col, ld1, ld2;
    #ifdef _MSC_VER
    register _Dcomplex *iMatrix, *oMatrix, *joMatrix;
    #else
    register double complex *iMatrix, *oMatrix, *joMatrix;
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
            *(joMatrix+i) = *(iMatrix+i*ld1+j);
        }
    }

    return pOMatrix;
}

MatrixD *MatrixDC_T2(MatrixDC *pIMatrix, MatrixD *pOMatrix)
{
    register int i, j, row, col, ld1, ld2;
    #ifdef _MSC_VER
    register _Dcomplex *iMatrix;
    #else
    register double complex *iMatrix;
    #endif /*_MSC_VER*/
    register double *oMatrix, *joMatrix;

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
            *(joMatrix+i) = (iMatrix+i*ld1+j)->_Val[0];
            #else
            *(joMatrix+i) = creal(*(iMatrix+i*ld1+j));
            #endif /*_MSC_VER*/
        }
    }

    return pOMatrix;
}

#ifdef _MSC_VER
_Dcomplex MatrixDC_Det(MatrixDC *pIMatrix)
#else
double complex MatrixDC_Det(MatrixDC *pIMatrix)
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
    MatrixDC *pTmpMatrix;

    row = pIMatrix->size[0];
    col = pIMatrix->size[1];
    ipiv = (int *)malloc(sizeof(int) * ((row > col) ? col : row));
    #ifdef _MSC_VER
    det = _Cbuild(1.0,0.0);
    #else
    det = 1.0;
    #endif /*_MSC_VER*/
    pTmpMatrix = MatrixDC_New(row,col);
    pTmpMatrix = MatrixDC_Copy(pIMatrix, pTmpMatrix);

    zgetrf_(&row, &col, pTmpMatrix->matrix, &row, ipiv, &info);

    for(i=0;i<row;i++)
    {
        #ifdef _MSC_VER
        det = _Cmulcc(det, *(pTmpMatrix->matrix+i*col+i));
        #else
        det *= *(pTmpMatrix->matrix+i*col+i);
        #endif /*_MSC_VER*/
    }

    MatrixDC_Del(pTmpMatrix);
    free(ipiv);
    return det;
}

double MatrixDC_Det2(MatrixDC *pIMatrix)
{
    int row, col, info;
    int *ipiv;
    register int i;
    #ifdef _MSC_VER
    _Dcomplex det;
    #else
    double complex det;
    #endif /*_MSC_VER*/
    double tmp;
    MatrixDC *pTmpMatrix;

    row = pIMatrix->size[0];
    col = pIMatrix->size[1];
    ipiv = (int *)malloc(sizeof(int) * ((row > col) ? col : row));
    #ifdef _MSC_VER
    det = _Cbuild(1.0,0.0);
    #else
    det = 1.0;
    #endif /*_MSC_VER*/
    pTmpMatrix = MatrixDC_New(row,col);
    pTmpMatrix = MatrixDC_Copy(pIMatrix, pTmpMatrix);

    zgetrf_(&row, &col, pTmpMatrix->matrix, &row, ipiv, &info);

    for(i=0;i<row;i++)
    {
        #ifdef _MSC_VER
        det = _Cmulcc(det, *(pTmpMatrix->matrix+i*col+i));
        #else
        det *= *(pTmpMatrix->matrix+i*col+i);
        #endif /*_MSC_VER*/
    }

    tmp = creal(det);
    MatrixDC_Del(pTmpMatrix);
    free(ipiv);
    return tmp;
}

MatrixDC *MatrixDC_Inv(MatrixDC *pIMatrix, MatrixDC *pOMatrix)
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

    MatrixDC_T(pIMatrix, pOMatrix);
    zgetrf_(&row, &col, pOMatrix->matrix, &ld2, ipiv, &info);
    zgetri_(&row, pOMatrix->matrix, &ld2, ipiv, work, &row, &info);

    free(ipiv);
    free(work);
    return pOMatrix;
}

MatrixD *MatrixDC_Inv2(MatrixDC *pIMatrix, MatrixD *pOMatrix)
{
    int row, col, info, ld2;
    int *ipiv;
    #ifdef _MSC_VER
    _Dcomplex *work;
    #else
    double complex *work;
    #endif /*_MSC_VER*/
    MatrixDC *pTmpMatrix;

    row = pIMatrix->size[0];
    col = pIMatrix->size[1];
    ld2 = pOMatrix->ld;
    ipiv = (int *)malloc(sizeof(int) * ((row > col) ? col : row));
    #ifdef _MSC_VER
    work = (_Dcomplex *)malloc(sizeof(_Dcomplex) * col *row);
    #else
    work = (double complex *)malloc(sizeof(double complex) * col * row);
    #endif /*_MSC_VER*/

    pTmpMatrix = MatrixDC_New(row, ld2);
    MatrixDC_T(pIMatrix, pTmpMatrix);
    zgetrf_(&row, &col, pTmpMatrix->matrix, &ld2, ipiv, &info);
    zgetri_(&row, pTmpMatrix->matrix, &ld2, ipiv, work, &row, &info);
    MatrixDC_Copy2(pTmpMatrix, pOMatrix);

    free(ipiv);
    free(work);
    return pOMatrix;
}

VectorDC *MatrixDC_MulVectorD(MatrixDC *pIMatrix, VectorD *pIVector, VectorDC *pOVector)
{
    register int i, l, iSize, oSize, ld, inc1, inc2;
    register double *iVector;
    #ifdef _MSC_VER
    register _Dcomplex *oVector, *iMatrix, *iiMatrix;
    #else
    register double complex *oVector, *iMatrix, *iiMatrix;
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
        *(oVector+i*inc2) = _Cbuild(0.0,0.0);
        #else
        *(oVector+i*inc2) = 0.0;
        #endif /*_MSC_VER*/
        iiMatrix = iMatrix + i*ld;
        for(l=0; l<iSize; ++l)
        {
            #ifdef _MSC_VER
            *(oVector+i*inc2) = _Caddcc(*(oVector+i*inc2), _Cmulcr(*(iiMatrix+l), *(iVector+l*inc1)));
            #else
            *(oVector+i*inc2) += *(iiMatrix+l) * *(iVector+l*inc1);
            #endif/*_MSC_VER*/
        }
    }
    return pOVector;
}

VectorD *MatrixDC_MulVectorD2(MatrixDC *pIMatrix, VectorD *pIVector, VectorD *pOVector)
{
    register int i, l, iSize, oSize, ld, inc1, inc2;
    register double *iVector, *oVector;
    #ifdef _MSC_VER
    register _Dcomplex *iMatrix, *iiMatrix;
    #else
    register double complex *iMatrix, *iiMatrix;
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
            *(oVector+i*inc2) += (iiMatrix+l)->_Val[0] * *(iVector+l*inc1);
            #else
            *(oVector+i*inc2) += creal(*(iiMatrix+l)) * *(iVector+l*inc1);
            #endif/*_MSC_VER*/
        }
    }
    return pOVector;
}

VectorDC *MatrixDC_MulVectorDC(MatrixDC *pIMatrix, VectorDC *pIVector, VectorDC *pOVector)
{
    register int i, l, iSize, oSize, ld, inc1, inc2;
    #ifdef _MSC_VER
    register _Dcomplex *iVector, *oVector, *iMatrix, *iiMatrix;
    #else
    register double complex *iVector, *oVector, *iMatrix, *iiMatrix;
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
        *(oVector+i*inc2) = _Cbuild(0.0,0.0);
        #else
        *(oVector+i*inc2) = 0.0;
        #endif /*_MSC_VER*/
        iiMatrix = iMatrix + i*ld;
        for(l=0; l<iSize; ++l)
        {
            #ifdef _MSC_VER
            *(oVector+i*inc2) = _Caddcc(*(oVector+i*inc2), _Cmulcc(*(iiMatrix+l), *(iVector+l*inc1)));
            #else
            *(oVector+i*inc2) += *(iiMatrix+l) * *(iVector+l*inc1);
            #endif/*_MSC_VER*/
        }
    }
    return pOVector;
}

VectorD *MatrixDC_MulVectorDC2(MatrixDC *pIMatrix, VectorDC *pIVector, VectorD *pOVector)
{
    register int i, l, iSize, oSize, ld, inc1, inc2;
    #ifdef _MSC_VER
    register _Dcomplex *iVector, *iMatrix, *iiMatrix;
    #else
    register double complex *iVector, *iMatrix, *iiMatrix;
    #endif /*_MSC_VER*/
    register double *oVector;

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
            *(oVector+i*inc2) +=  ((iiMatrix+l)->_Val[0] * (iVector+l*inc1)->_Val[0] - (iiMatrix+l)->_Val[1] * (iVector+l*inc1)->_Val[1]);
            #else
            *(oVector+i*inc2) += creal(*(iiMatrix+l) * *(iVector+l*inc1));
            #endif/*_MSC_VER*/
        }
    }
    return pOVector;
}

VectorDC *MatrixDC_TMulVectorD(MatrixDC *pIMatrix, VectorD *pIVector, VectorDC *pOVector)
{
    register int i, j, row, col, ld, inc1, inc2;
    register double tmp;
    register double *iVector;
    #ifdef _MSC_VER
    register _Dcomplex *iMatrix, *oVector, *iiMatrix;
    #else
    register double complex *iMatrix, *oVector, *iiMatrix;
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
        *(oVector+j*inc2) = _Cbuild(0.0,0.0);
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
            *(oVector+j) = _Caddcc(*(oVector+j),_Cmulcr(*(iiMatrix+j),tmp));
            #else
            *(oVector+j) +=  *(iiMatrix+j) * tmp;
            #endif /*_MSC_VER*/
        }
    }

    return pOVector;
}

VectorD *MatrixDC_TMulVectorD2(MatrixDC *pIMatrix, VectorD *pIVector, VectorD *pOVector)
{
    register int i, j, row, col, ld, inc1, inc2;
    register double tmp;
    register double *iVector, *oVector;
    #ifdef _MSC_VER
    register _Dcomplex *iMatrix, *iiMatrix;
    #else
    register double complex *iMatrix, *iiMatrix;
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
            *(oVector+j) += (iiMatrix+j)->_Val[0] * tmp;
            #else
            *(oVector+j) +=  creal(*(iiMatrix+j)) * tmp;
            #endif /*_MSC_VER*/
        }
    }

    return pOVector;
}

VectorDC *MatrixDC_TMulVectorDC(MatrixDC *pIMatrix, VectorDC *pIVector, VectorDC *pOVector)
{
    register int i, j, row, col, ld, inc1, inc2;
    #ifdef _MSC_VER
    register _Dcomplex *iMatrix, *iVector, *oVector, *iiMatrix;
    register _Dcomplex tmp;
    #else
    register double complex *iMatrix, *iVector, *oVector, *iiMatrix;
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
        *(oVector+j*inc2) = _Cbuild(0.0,0.0);
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
            *(oVector+j) = _Caddcc(*(oVector+j),_Cmulcc(*(iiMatrix+j),tmp));
            #else
            *(oVector+j) +=  *(iiMatrix+j) * tmp;
            #endif /*_MSC_VER*/
        }
    }

    return pOVector;
}

VectorD *MatrixDC_TMulVectorDC2(MatrixDC *pIMatrix, VectorDC *pIVector, VectorD *pOVector)
{
    register int i, j, row, col, ld, inc1, inc2;
    #ifdef _MSC_VER
    register _Dcomplex *iMatrix, *iVector, *iiMatrix;
    register _Dcomplex tmp;
    #else
    register double complex *iMatrix, *iVector, *iiMatrix;
    register double complex tmp;
    #endif /*_MSC_VER*/
    register double  *oVector;

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
            *(oVector+j) += ((iiMatrix+j)->_Val[0] * tmp._Val[0] - (iiMatrix+j)->_Val[1] * tmp._Val[1]);
            #else
            *(oVector+j) +=  creal(*(iiMatrix+j) * tmp);
            #endif /*_MSC_VER*/
        }
    }

    return pOVector;
}

VectorDC *MatrixDC_CMulVectorD(MatrixDC *pIMatrix, VectorD *pIVector, VectorDC *pOVector)
{
    register int i, j, row, col, ld, inc1, inc2;
    register double tmp;
    register double *iVector;
    #ifdef _MSC_VER
    register _Dcomplex *iMatrix, *oVector, *iiMatrix;
    #else
    register double complex *iMatrix, *oVector, *iiMatrix;
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
        *(oVector+j*inc2) = _Cbuild(0.0,0.0);
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
            *(oVector+j*inc2) = _Caddcc(*(oVector+j*inc2),_Cmulcr(conj(*(iiMatrix+j)), tmp));
            #else
            *(oVector+j*inc2) +=  conj(*(iiMatrix+j)) * tmp;
            #endif/*_MSC_VER*/
        }
    }

    return pOVector;
}

VectorD *MatrixDC_CMulVectorD2(MatrixDC *pIMatrix, VectorD *pIVector, VectorD *pOVector)
{
    register int i, j, row, col, ld, inc1, inc2;
    register double tmp;
    register double *iVector, *oVector;
    #ifdef _MSC_VER
    register _Dcomplex *iMatrix, *iiMatrix;
    #else
    register double complex *iMatrix, *iiMatrix;
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
            *(oVector+j*inc2) += (iiMatrix+j)->_Val[0] * tmp;
            #else
            *(oVector+j*inc2) +=  creal(*(iiMatrix+j)) * tmp;
            #endif/*_MSC_VER*/
        }
    }

    return pOVector;
}

VectorDC *MatrixDC_CMulVectorDC(MatrixDC *pIMatrix, VectorDC *pIVector, VectorDC *pOVector)
{
    register int i, j, row, col, ld, inc1, inc2;
    #ifdef _MSC_VER
    register _Dcomplex *iMatrix, *iVector, *oVector, *iiMatrix;
    register _Dcomplex tmp;
    #else
    register double complex *iMatrix, *iVector, *oVector, *iiMatrix;
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
        *(oVector+j*inc2) = _Cbuild(0.0,0.0);
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
            *(oVector+j) = _Caddcc(*(oVector+j*inc2),_Cmulcc(conj(*(iiMatrix+j)), tmp));
            #else
            *(oVector+j) +=  conj(*(iiMatrix+j)) * tmp;
            #endif/*_MSC_VER*/
        }
    }

    return pOVector;
}

VectorD *MatrixDC_CMulVectorDC2(MatrixDC *pIMatrix, VectorDC *pIVector, VectorD *pOVector)
{
    register int i, j, row, col, ld, inc1, inc2;
    #ifdef _MSC_VER
    register _Dcomplex *iMatrix, *iVector, *iiMatrix;
    register _Dcomplex tmp;
    #else
    register double complex *iMatrix, *iVector, *iiMatrix;
    register double complex tmp;
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
            *(oVector+j*inc2) += ((iiMatrix+j)->_Val[0] * tmp._Val[0] + (iiMatrix+j)->_Val[1] * tmp._Val[1]);
            #else
            *(oVector+j*inc2) += creal(conj(*(iiMatrix+j)) * tmp);
            #endif/*_MSC_VER*/
        }
    }

    return pOVector;
}

MatrixDC *MatrixDC_TMulSelf(MatrixDC *pIMatrix, MatrixDC *pOMatrix)
{
    register int i, j, l, row, col, ld1, ld2;
    #ifdef _MSC_VER
    register _Dcomplex *iMatrix, *oMatrix, *liMatrix, *ioMatrix;
    register _Dcomplex tmp;
    #else
    register double complex *iMatrix, *oMatrix, *liMatrix, *ioMatrix;
    register double complex tmp;
    #endif /*_MSC_VER*/

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
                (ioMatrix+j)->_Val[0] += (tmp._Val[0] * (liMatrix+j)->_Val[0] - tmp._Val[1] * (liMatrix+j)->_Val[1]);
                (ioMatrix+j)->_Val[1] += tmp._Val[0] * (liMatrix+j)->_Val[1] + tmp._Val[1] * (liMatrix+j)->_Val[0];
                #else
                *(ioMatrix+j) += tmp * *(liMatrix+j);
                #endif/*_MSC_VER*/
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

MatrixD *MatrixDC_TMulSelf2(MatrixDC *pIMatrix, MatrixD *pOMatrix)
{
    register int i, j, l, row, col, ld1, ld2;
    #ifdef _MSC_VER
    register _Dcomplex *iMatrix, *liMatrix;
    register _Dcomplex tmp;
    #else
    register double complex *iMatrix, *liMatrix;
    register double complex tmp;
    #endif /*_MSC_VER*/
    register double *oMatrix, *ioMatrix;

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
                #ifdef _MSC_VER
                *(ioMatrix+j) += (tmp._Val[0] * (liMatrix+j)->_Val[0] - tmp._Val[1] * (liMatrix+j)->_Val[1]);
                #else
                *(ioMatrix+j) += creal(tmp * *(liMatrix+j));
                #endif/*_MSC_VER*/
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

MatrixDC *MatrixDC_O2Mul(MatrixDC *pIMatrix1, MatrixDC *pIMatrix2, MatrixDC *pOMatrix)
{
    int ld1, ld2, ld3;
    #ifdef _MSC_VER
    register _Dcomplex *iMatrix1, *iMatrix2, *oMatrix;
    _Dcomplex m11, m12, m13, m14, m21, m22, m23, m24;
    #else
    register double complex *iMatrix1, *iMatrix2, *oMatrix;
    register double complex m11, m12, m13, m14, m21, m22, m23, m24;
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
    *(oMatrix) = _Caddcc(_Cmulcc(m11, m21), _Cmulcc(m12, m23));
    *(oMatrix+1) = _Caddcc(_Cmulcc(m11, m22), _Cmulcc(m12, m24));
    *(oMatrix+ld3) = _Caddcc(_Cmulcc(m13, m21), _Cmulcc(m14, m23));
    *(oMatrix+ld3+1) = _Caddcc(_Cmulcc(m13, m22), _Cmulcc(m14, m24));
    #else
    *(oMatrix) = m11 * m21 + m12 * m23;
    *(oMatrix+1) = m11 * m22 + m12 * m24;
    *(oMatrix+ld3) = m13 * m21 + m14 * m23;
    *(oMatrix+ld3+1) = m13 * m22 + m14 * m24;
    #endif /*_MSC_VER*/

    return pOMatrix;
}

MatrixDC *MatrixDC_O2Mul2(MatrixDC *pIMatrix1, MatrixD *pIMatrix2, MatrixDC *pOMatrix)
{
    int ld1, ld2, ld3;
    register double *iMatrix2;
    register double  m21, m22, m23, m24;
    #ifdef _MSC_VER
    register _Dcomplex *iMatrix1, *oMatrix;
    _Dcomplex m11, m12, m13, m14;
    #else
    register double complex *iMatrix1, *oMatrix;
    register double complex m11, m12, m13, m14;
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
    *(oMatrix) = _Caddcc(_Cmulcr(m11, m21), _Cmulcr(m12, m23));
    *(oMatrix+1) = _Caddcc(_Cmulcr(m11, m22), _Cmulcr(m12, m24));
    *(oMatrix+ld3) = _Caddcc(_Cmulcr(m13, m21), _Cmulcr(m14, m23));
    *(oMatrix+ld3+1) = _Caddcc(_Cmulcr(m13, m22), _Cmulcr(m14, m24));
    #else
    *(oMatrix) = m11 * m21 + m12 * m23;
    *(oMatrix+1) = m11 * m22 + m12 * m24;
    *(oMatrix+ld3) = m13 * m21 + m14 * m23;
    *(oMatrix+ld3+1) = m13 * m22 + m14 * m24;
    #endif /*_MSC_VER*/

    return pOMatrix;
}

MatrixD *MatrixDC_O2Mul3(MatrixDC *pIMatrix1, MatrixDC *pIMatrix2, MatrixD *pOMatrix)
{
    int ld1, ld2, ld3;
    register double *oMatrix;
    #ifdef _MSC_VER
    register _Dcomplex *iMatrix1, *iMatrix2;
    _Dcomplex m11, m12, m13, m14, m21, m22, m23, m24;
    #else
    register double complex *iMatrix1, *iMatrix2;
    register double complex m11, m12, m13, m14, m21, m22, m23, m24;
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
    *(oMatrix) = m11._Val[0] * m21._Val[0] - m11._Val[1] + m21._Val[1] + \
        m12._Val[0] * m23._Val[0] - m12._Val[1] * m23._Val[1];
    *(oMatrix+1) = m11._Val[0] * m22._Val[0] - m11._Val[1] * m22._Val[1] + \
        m12._Val[0] * m24._Val[0] - m12._Val[1] * m24._Val[1];
    *(oMatrix+ld3) = m13._Val[0] * m21._Val[0] - m13._Val[1] * m21._Val[1] + \
        m14._Val[0] * m23._Val[0] - m14._Val[1] * m23._Val[1];
    *(oMatrix+ld3+1) = m13._Val[0] * m22._Val[0] - m13._Val[1] * m22._Val[1] + \
        m14._Val[0] * m24._Val[0] - m14._Val[1] * m24._Val[1];
    #else
    *(oMatrix) = creal(m11 * m21 + m12 * m23);
    *(oMatrix+1) = creal(m11 * m22 + m12 * m24);
    *(oMatrix+ld3) = creal(m13 * m21 + m14 * m23);
    *(oMatrix+ld3+1) = creal(m13 * m22 + m14 * m24);
    #endif /*_MSC_VER*/

    return pOMatrix;
}

MatrixD *MatrixDC_O2Mul4(MatrixDC *pIMatrix1, MatrixD *pIMatrix2, MatrixD *pOMatrix)
{
    int ld1, ld2, ld3;
    register double *iMatrix2, *oMatrix;
    register double  m21, m22, m23, m24;
    #ifdef _MSC_VER
    register _Dcomplex *iMatrix1;
    _Dcomplex m11, m12, m13, m14;
    #else
    register double complex *iMatrix1;
    register double complex m11, m12, m13, m14;
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
    *(oMatrix) = m11._Val[0] * m21 + m12._Val[0] * m23;
    *(oMatrix+1) = m11._Val[0] * m22 + m12._Val[0] * m24;
    *(oMatrix+ld3) = m13._Val[0] * m21 + m14._Val[0] * m23;
    *(oMatrix+ld3+1) = m13._Val[0] * m22 + m14._Val[0] * m24;
    #else
    *(oMatrix) = creal(m11) * m21 + creal(m12) * m23;
    *(oMatrix+1) = creal(m11) * m22 + creal(m12) * m24;
    *(oMatrix+ld3) = creal(m13) * m21 + creal(m14) * m23;
    *(oMatrix+ld3+1) = creal(m13) * m22 + creal(m14) * m24;
    #endif /*_MSC_VER*/

    return pOMatrix;
}

MatrixDC *MatrixDC_MulMatrixD(MatrixDC *pIMatrix1, MatrixD *pIMatrix2, MatrixDC *pOMatrix)
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
    lda = pIMatrix1->ld;
    ldc = pOMatrix->ld;
    m = pIMatrix1->size[1];
    n = pIMatrix2->size[0];
    k = pIMatrix1->size[0];
    pTmpMatrix = MatrixDC_New(n,k);
    pTmpMatrix = MatrixD_Copy2(pIMatrix2, pTmpMatrix);
    ldb = pIMatrix2->ld;

    zgemm_("T","T",&m,&n,&k,&alpha,pIMatrix1->matrix,&lda,pTmpMatrix->matrix,&ldb,&beta,pOMatrix->matrix,&ldc);

    pTmpMatrix = MatrixDC_Del(pTmpMatrix);
    return pOMatrix;
}

MatrixD *MatrixDC_MulMatrixD2(MatrixDC *pIMatrix1, MatrixD *pIMatrix2, MatrixD *pOMatrix)
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
    pTmpMatrix = MatrixD_New(k,m);
    MatrixDC_Copy2(pIMatrix1, pTmpMatrix);
    ldb = pIMatrix2->ld;

    dgemm_("T","T",&m,&n,&k,&alpha,pTmpMatrix->matrix,&lda,pIMatrix2->matrix,&ldb,&beta,pOMatrix->matrix,&ldc);

    pTmpMatrix = MatrixD_Del(pTmpMatrix);
    return pOMatrix;
}

MatrixDC *MatrixDC_MulMatrixDC(MatrixDC *pIMatrix1, MatrixDC *pIMatrix2, MatrixDC *pOMatrix)
{
    #ifdef _MSC_VER
    _Dcomplex alpha, beta;
    #else
    double complex alpha, beta;
    #endif /*_MSC_VER*/
    int m, n, k, lda, ldb, ldc;

    #ifdef _MSC_VER
    alpha = _Cbuild(1,0);
    beta = _Cbuild(0,0);
    #else
    alpha = 1;
    beta = 0;
    #endif /*_MSC_VER*/
    lda = pIMatrix1->ld;
    ldb = pIMatrix2->ld;
    ldc = pOMatrix->ld;
    m = pIMatrix1->size[1];
    n = pIMatrix2->size[0];
    k = pIMatrix1->size[0];

    zgemm_("T","T",&m,&n,&k,&alpha,pIMatrix1->matrix,&lda,pIMatrix2->matrix,&ldb,&beta,pOMatrix->matrix,&ldc);

    return pOMatrix;
}

MatrixD *MatrixDC_MulMatrixDC2(MatrixDC *pIMatrix1, MatrixDC *pIMatrix2, MatrixD *pOMatrix)
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
    lda = pIMatrix1->ld;
    ldb = pIMatrix2->ld;
    ldc = pOMatrix->ld;
    m = pIMatrix1->size[1];
    n = pIMatrix2->size[0];
    k = pIMatrix1->size[0];
    pTmpMatrix = MatrixDC_New(pOMatrix->size[0],pOMatrix->size[1]);

    zgemm_("T","T",&m,&n,&k,&alpha,pIMatrix1->matrix,&lda,pIMatrix2->matrix,&ldb,&beta,pTmpMatrix->matrix,&ldc);
    MatrixDC_Copy2(pTmpMatrix, pOMatrix);

    return pOMatrix;
}

#ifdef _MSC_VER
_Dcomplex MatrixDC_MaxDiagElem(MatrixDC *pIMatrix)
#else
double complex MatrixDC_MaxDiagElem(MatrixDC *pIMatrix)
#endif /*_MSC_VER*/
{
    register int i, size, ld;
    register double max;
    #ifdef _MSC_VER
    register _Dcomplex *iMatrix;
    register _Dcomplex max_value;
    #else
    register double complex *iMatrix;
    register double complex max_value;
    #endif /*_MSC_VER*/

    size = pIMatrix->size[0];
    ld = pIMatrix->ld;
    iMatrix = pIMatrix->matrix;
    max = -DBL_MAX;
    #ifdef _MSC_VER
    max_value = _Cbuild(0.0,0.0);
    #else
    max_value = 0;
    #endif /*_MSC_VER*/

    if(size != pIMatrix->size[1])
    {
        fprintf(stderr, "Error: pIMatrix isn't a square matrix. File: %s Func: %s Line: %d\n", \
            __FILE__, __FUNCTION__, __LINE__);
        return max_value;
    }

    for(i=0;i<size;++i)
    {
        if(max < cabs(*(iMatrix+i*ld+i)))
        {
            max = cabs(*(iMatrix+i*ld+i));
            max_value = *(iMatrix+i*ld+i);
        }
    }

    return max_value;
}

double MatrixDC_MaxDiagElem2(MatrixDC *pIMatrix)
{
    register int i, size, ld;
    register double max;
    #ifdef _MSC_VER
    register _Dcomplex *iMatrix;
    register _Dcomplex max_value;
    #else
    register double complex *iMatrix;
    register double complex max_value;
    #endif /*_MSC_VER*/

    size = pIMatrix->size[0];
    ld = pIMatrix->ld;
    iMatrix = pIMatrix->matrix;
    max = -DBL_MAX;
    #ifdef _MSC_VER
    max_value = _Cbuild(0.0,0.0);
    #else
    max_value = 0;
    #endif /*_MSC_VER*/

    if(size != pIMatrix->size[1])
    {
        fprintf(stderr, "Error: pIMatrix isn't a square matrix. File: %s Func: %s Line: %d\n", \
            __FILE__, __FUNCTION__, __LINE__);
        return max;
    }

    for(i=0;i<size;++i)
    {
        if(max < cabs(*(iMatrix+i*ld+i)))
        {
            max = cabs(*(iMatrix+i*ld+i));
            max_value = *(iMatrix+i*ld+i);
        }
    }

    #ifdef _MSC_VER
    max = max_value._Val[0];
    #else
    max = creal(max_value);
    #endif /*_MSC_VER*/

    return max;
}

VectorDC *MatrixDC_LinearSolver(MatrixDC *pIMatrix, VectorDC *pIVector, VectorDC *pOVector)
{
    int row, col, info, nrhs;
    int *ipiv;
    MatrixDC *pTmpMatrix;

    nrhs = 1;
    row = pIMatrix->size[0];
    col = pIMatrix->size[1];
    ipiv = (int *)malloc(sizeof(int) * ((row > col) ? col : row));
    pTmpMatrix = MatrixDC_New(col, row);

    pTmpMatrix = MatrixDC_T(pIMatrix, pTmpMatrix);
    pOVector = VectorDC_Copy(pIVector, pOVector);
    zgetrf_(&row, &col, pTmpMatrix->matrix, &row, ipiv, &info);
    zgetrs_("N", &row, &nrhs, pTmpMatrix->matrix, &row, ipiv, pOVector->vector, &row, &info);

    free(ipiv);
    pTmpMatrix = MatrixDC_Del(pTmpMatrix);
    return pOVector;
}

VectorDC *MatrixDC_LinearSolver2(MatrixDC *pIMatrix, VectorD *pIVector, VectorDC *pOVector)
{
    int row, col, info, nrhs;
    int *ipiv;
    MatrixDC *pTmpMatrix;

    nrhs = 1;
    row = pIMatrix->size[0];
    col = pIMatrix->size[1];
    ipiv = (int *)malloc(sizeof(int) * ((row > col) ? col : row));
    pTmpMatrix = MatrixDC_New(col, row);

    pTmpMatrix = MatrixDC_T(pIMatrix, pTmpMatrix);
    pOVector = VectorD_Copy2(pIVector, pOVector);
    zgetrf_(&row, &col, pTmpMatrix->matrix, &row, ipiv, &info);
    zgetrs_("N", &row, &nrhs, pTmpMatrix->matrix, &row, ipiv, pOVector->vector, &row, &info);

    free(ipiv);
    pTmpMatrix = MatrixDC_Del(pTmpMatrix);
    return pOVector;
}

VectorD *MatrixDC_LinearSolver3(MatrixDC *pIMatrix, VectorD *pIVector, VectorD *pOVector)
{
    int row, col, info, nrhs;
    int *ipiv;
    VectorDC *pTmpVector;
    MatrixDC *pTmpMatrix;

    nrhs = 1;
    row = pIMatrix->size[0];
    col = pIMatrix->size[1];
    ipiv = (int *)malloc(sizeof(int) * ((row > col) ? col : row));
    pTmpVector = VectorDC_New(row);
    pTmpMatrix = MatrixDC_New(col, row);

    pTmpMatrix = MatrixDC_T(pIMatrix, pTmpMatrix);
    pTmpVector = VectorD_Copy2(pIVector, pTmpVector);
    zgetrf_(&row, &col, pTmpMatrix->matrix, &row, ipiv, &info);
    zgetrs_("N", &row, &nrhs, pTmpMatrix->matrix, &row, ipiv, pTmpVector->vector, &row, &info);
    VectorDC_Copy2(pTmpVector, pOVector);

    free(ipiv);
    pTmpMatrix = MatrixDC_Del(pTmpMatrix);
    return pOVector;
}

VectorD *MatrixDC_LinearSolver4(MatrixDC *pIMatrix, VectorDC *pIVector, VectorD *pOVector)
{
    int row, col, info, nrhs;
    int *ipiv;
    VectorDC *pTmpVector;
    MatrixDC *pTmpMatrix;

    nrhs = 1;
    row = pIMatrix->size[0];
    col = pIMatrix->size[1];
    ipiv = (int *)malloc(sizeof(int) * ((row > col) ? col : row));
    pTmpVector = VectorDC_New(row);
    pTmpMatrix = MatrixDC_New(col, row);

    pTmpMatrix = MatrixDC_T(pIMatrix, pTmpMatrix);
    pTmpVector = VectorDC_Copy(pIVector, pTmpVector);
    zgetrf_(&row, &col, pTmpMatrix->matrix, &row, ipiv, &info);
    zgetrs_("N", &row, &nrhs, pTmpMatrix->matrix, &row, ipiv, pTmpVector->vector, &row, &info);
    VectorDC_Copy2(pTmpVector, pOVector);

    free(ipiv);
    pTmpMatrix = MatrixDC_Del(pTmpMatrix);
    return pOVector;
}

void MatrixDC_EigenEquation(MatrixDC *pIMatrix, VectorDC *pOVector, MatrixDC *pOMatrix)
{
    int n, ldvl, ldvr, lwork, info;
    double *rwork;
    #ifdef _MSC_VER
    _Dcomplex wkopt;
    _Dcomplex *work;
    #else
    double complex wkopt;
    double complex *work;
    #endif /*_MSC_VER*/
    MatrixDC *pTmpMatrix, *pTmpMatrix2;

    n = pIMatrix->size[0];
    ldvl = n;
    ldvr = pOMatrix->ld;
    lwork = -1;
    rwork = (double *)malloc(sizeof(double) * 2 * n);
    pTmpMatrix = MatrixDC_New(n,n);
    pTmpMatrix2 = MatrixDC_New(n,n);

    pTmpMatrix = MatrixDC_T(pIMatrix, pTmpMatrix);
    zgeev_("N","V",&n,pTmpMatrix->matrix,&n,pOVector->vector,NULL,&ldvl,pTmpMatrix2->matrix,&ldvr, &wkopt, &lwork, rwork, &info);
    lwork = (int) creal(wkopt);
    #ifdef _MSC_VER
    work = (_Dcomplex *)malloc(sizeof(_Dcomplex) * lwork);
    #else
    work = (double complex *)malloc(sizeof(double complex) * lwork);
    #endif /*_MSC_VER*/
    zgeev_("N","V",&n,pTmpMatrix->matrix,&n,pOVector->vector,NULL,&ldvl,pTmpMatrix2->matrix,&ldvr, work, &lwork, rwork, &info);
    pOMatrix = MatrixDC_T(pTmpMatrix2,pOMatrix);

    pTmpMatrix = MatrixDC_Del(pTmpMatrix);
    pTmpMatrix2 = MatrixDC_Del(pTmpMatrix2);
    free(work);
    free(rwork);
}

void MatrixDC_EigenEquation2(MatrixDC *pIMatrix, VectorD *pOVector, MatrixD *pOMatrix)
{
    int n, ldvl, ldvr, lwork, info;
    double *rwork;
    #ifdef _MSC_VER
    _Dcomplex wkopt;
    _Dcomplex *work;
    #else
    double complex wkopt;
    double complex *work;
    #endif /*_MSC_VER*/
    MatrixDC *pTmpMatrix, *pTmpMatrix2;
    VectorDC *pTmpVector;

    n = pIMatrix->size[0];
    ldvl = n;
    ldvr = pOMatrix->ld;
    lwork = -1;
    rwork = (double *)malloc(sizeof(double) * 2 * n);
    pTmpMatrix = MatrixDC_New(n,n);
    pTmpMatrix2 = MatrixDC_New(n,n);
    pTmpVector = VectorDC_New(n);

    pTmpMatrix = MatrixDC_T(pIMatrix, pTmpMatrix);
    zgeev_("N","V",&n,pTmpMatrix->matrix,&n,pTmpVector->vector,NULL,&ldvl,pTmpMatrix2->matrix,&ldvr, &wkopt, &lwork, rwork, &info);
    lwork = (int) creal(wkopt);
    #ifdef _MSC_VER
    work = (_Dcomplex *)malloc(sizeof(_Dcomplex) * lwork);
    #else
    work = (double complex *)malloc(sizeof(double complex) * lwork);
    #endif /*_MSC_VER*/
    zgeev_("N","V",&n,pTmpMatrix->matrix,&n,pTmpVector->vector,NULL,&ldvl,pTmpMatrix2->matrix,&ldvr, work, &lwork, rwork, &info);
    pOMatrix = MatrixDC_T2(pTmpMatrix2,pOMatrix);
    pOVector = VectorDC_Copy2(pTmpVector, pOVector);

    pTmpMatrix = MatrixDC_Del(pTmpMatrix);
    pTmpMatrix2 = MatrixDC_Del(pTmpMatrix2);
    pTmpVector = VectorDC_Del(pTmpVector);
    free(work);
    free(rwork);
}

void MatrixDC_EigenEquation3(MatrixDC *pIMatrix, VectorDC *pOVector, MatrixD *pOMatrix)
{
    int n, ldvl, ldvr, lwork, info;
    double *rwork;
    #ifdef _MSC_VER
    _Dcomplex wkopt;
    _Dcomplex *work;
    #else
    double complex wkopt;
    double complex *work;
    #endif /*_MSC_VER*/
    MatrixDC *pTmpMatrix, *pTmpMatrix2;
    VectorDC *pTmpVector;

    n = pIMatrix->size[0];
    ldvl = n;
    ldvr = pOMatrix->ld;
    lwork = -1;
    rwork = (double *)malloc(sizeof(double) * 2 * n);
    pTmpMatrix = MatrixDC_New(n,n);
    pTmpMatrix2 = MatrixDC_New(n,n);
    pTmpVector = VectorDC_New(n);

    pTmpMatrix = MatrixDC_T(pIMatrix, pTmpMatrix);
    zgeev_("N","V",&n,pTmpMatrix->matrix,&n,pTmpVector->vector,NULL,&ldvl,pTmpMatrix2->matrix,&ldvr, &wkopt, &lwork, rwork, &info);
    lwork = (int) creal(wkopt);
    #ifdef _MSC_VER
    work = (_Dcomplex *)malloc(sizeof(_Dcomplex) * lwork);
    #else
    work = (double complex *)malloc(sizeof(double complex) * lwork);
    #endif /*_MSC_VER*/
    zgeev_("N","V",&n,pTmpMatrix->matrix,&n,pTmpVector->vector,NULL,&ldvl,pTmpMatrix2->matrix,&ldvr, work, &lwork, rwork, &info);
    pOMatrix = MatrixDC_T2(pTmpMatrix2,pOMatrix);
    pOVector = VectorDC_Copy(pTmpVector, pOVector);

    pTmpMatrix = MatrixDC_Del(pTmpMatrix);
    pTmpMatrix2 = MatrixDC_Del(pTmpMatrix2);
    pTmpVector = VectorDC_Del(pTmpVector);
    free(work);
    free(rwork);
}

void MatrixDC_EigenEquation4(MatrixDC *pIMatrix, VectorD *pOVector, MatrixDC *pOMatrix)
{
    int n, ldvl, ldvr, lwork, info;
    double *rwork;
    #ifdef _MSC_VER
    _Dcomplex wkopt;
    _Dcomplex *work;
    #else
    double complex wkopt;
    double complex *work;
    #endif /*_MSC_VER*/
    MatrixDC *pTmpMatrix, *pTmpMatrix2;
    VectorDC *pTmpVector;

    n = pIMatrix->size[0];
    ldvl = n;
    ldvr = pOMatrix->ld;
    lwork = -1;
    rwork = (double *)malloc(sizeof(double) * 2 * n);
    pTmpMatrix = MatrixDC_New(n,n);
    pTmpMatrix2 = MatrixDC_New(n,n);
    pTmpVector = VectorDC_New(n);

    pTmpMatrix = MatrixDC_T(pIMatrix, pTmpMatrix);
    zgeev_("N","V",&n,pTmpMatrix->matrix,&n,pTmpVector->vector,NULL,&ldvl,pTmpMatrix2->matrix,&ldvr, &wkopt, &lwork, rwork, &info);
    lwork = (int) creal(wkopt);
    #ifdef _MSC_VER
    work = (_Dcomplex *)malloc(sizeof(_Dcomplex) * lwork);
    #else
    work = (double complex *)malloc(sizeof(double complex) * lwork);
    #endif /*_MSC_VER*/
    zgeev_("N","V",&n,pTmpMatrix->matrix,&n,pTmpVector->vector,NULL,&ldvl,pTmpMatrix2->matrix,&ldvr, work, &lwork, rwork, &info);
    pOMatrix = MatrixDC_T(pTmpMatrix2,pOMatrix);
    pOVector = VectorDC_Copy2(pTmpVector, pOVector);

    pTmpMatrix = MatrixDC_Del(pTmpMatrix);
    pTmpMatrix2 = MatrixDC_Del(pTmpMatrix2);
    pTmpVector = VectorDC_Del(pTmpVector);
    free(work);
    free(rwork);
}