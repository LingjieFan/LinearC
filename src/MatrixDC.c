#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include "Blas.h"
#include "Lapack.h"
#include "IObject.h"
#include "Class.h"
#include "Object.h"
#include "VectorDC.h"
#include "MatrixDC.h"

Class MatrixDC_Class = {OBJECT};

struct _MatrixDC
{
    Class* class;
    IObject* iObject;
    double complex* matrix;
    int row;
    int col;
    int ld;
};

void* _MatrixDC_New(void* implementor);
void* _MatrixDC_Del(void* implementor);
int _MatrixDC_Copy(void* implementor, void* object);
int _MatrixDC_Equal(void* implementor, void* object);

MatrixDC* MatrixDC_New(int row, int col)
{
    MatrixDC* matrix;
    
    matrix = (MatrixDC*)malloc(sizeof(MatrixDC)+sizeof(IObject)+sizeof(double complex)*col*row);
    matrix->class = MATRIXDC;
    matrix->iObject = (IObject*)(matrix+1);
    matrix->iObject->implementor = matrix;
    matrix->iObject->New = _MatrixDC_New;
    matrix->iObject->Del = _MatrixDC_Del;
    matrix->iObject->Copy = _MatrixDC_Copy;
    matrix->iObject->Equal = _MatrixDC_Equal;
    matrix->row = row;
    matrix->col = col;
    matrix->ld = col;
    matrix->matrix = (double complex*)(matrix->iObject+1);    
    
    return matrix;
}

MatrixDC* MatrixDC_Del(MatrixDC* this)
{
    free(this);
    
    return NULL;
}

MatrixDC* MatrixDC_Wrap(double complex* work, int row, int col, int ld)
{
    MatrixDC* matrix;
    
    matrix = (MatrixDC*)malloc(sizeof(MatrixDC)+sizeof(IObject));
    matrix->class = MATRIXDC;
    matrix->iObject = (IObject*)(matrix+1);
    matrix->iObject->implementor = matrix->matrix;
    matrix->iObject->New = _MatrixDC_New;
    matrix->iObject->Del = _MatrixDC_Del;
    matrix->iObject->Copy = _MatrixDC_Copy;
    matrix->iObject->Equal = _MatrixDC_Equal;
    matrix->row = row;
    matrix->col = col;
    matrix->ld = ld;
    matrix->matrix = work;    
    
    return matrix;    
}

MatrixDC* MatrixDC_ReWrap(MatrixDC* this, double complex* work, int row, int col, int ld)
{
    this->matrix = work;
    this->row = row;
    this->col = col;
    this->ld = ld;

    return this;
}

MatrixDC* MatrixDC_UnWrap(MatrixDC* this)
{
    free(this->matrix);
    free(this);
    
    return NULL;
}

IObject* MatrixDC_GetIObject(MatrixDC* this)
{
    return this->iObject;
}

int MatrixDC_GetRow(MatrixDC* this)
{    
    return this->row;
}

int MatrixDC_GetCol(MatrixDC* this)
{
    return this->col;
}

int MatrixDC_GetLd(MatrixDC* this)
{
    return this->ld;
}

double complex* MatrixDC_GetMatrix(MatrixDC* this)
{
    return this->matrix;
}

double complex MatrixDC_GetElement(MatrixDC* this, int row_index, int col_index)
{
    return *(this->matrix+row_index*this->ld+col_index);
}

MatrixDC* MatrixDC_SetElement(MatrixDC* this, int row_index, int col_index, double complex value)
{
    *(this->matrix+row_index*this->ld+col_index) = value;
    
    return this;
}

void MatrixDC_Show(MatrixDC* this)
{
    register int i, j, row, col, ld;
    register double complex* matrix;

    row = this->row;
    col = this->col;
    ld = this->ld;
    matrix = this->matrix;
    printf("MatrixDC_Show:\n");
    for(i=0;i<row;i++)
    {
        printf(" [");
        for(j=0;j<col;j++)
        {
            printf("%g+%g j, ", creal(*(matrix+i*ld+j)), cimag(*(matrix+i*ld+j)));
        }
        printf("]\n");
    }
    printf(" row:%d \n col:%d \n ld:%d \n\n", row, col, ld);
}

MatrixDC* MatrixDC_Set(MatrixDC* this, double complex diag, double complex off_diag)
{
    int row, col;
    register int i, j, ld;
    register double complex* matrix;

    row = this->row;
    col = this->col;
    ld = this->ld;
    matrix = this->matrix;
    for(i=0;i<row;i++)
    {
        for(j=0;j<col;j++)
        {
            if(i==j)
            {
                *(matrix+i*ld+j) = diag;
            }
            else
            {
                *(matrix+i*ld+j) = off_diag;
            }
        }
    }
    
    return this;
}

MatrixDC* MatrixDC_Copy(MatrixDC* this, MatrixDC* out_matrix)
{
    int row, col;
    register int i, j, ld1, ld2;
    register double complex* matrix1,* matrix2;
    
    row = this->row;
    col = this->col;
    ld1 = this->ld;
    ld2 = out_matrix->ld;
    matrix1 = this->matrix;
    matrix2 = out_matrix->matrix;
    for(i=0;i<row;i++)
    {
        for(j=0;j<col;j++)
        {
            *(matrix2+i*ld2+j) = *(matrix1+i*ld1+j);
        }
    }
    
    return out_matrix;
}

MatrixDC* MatrixDC_T(MatrixDC* this, MatrixDC* out_matrix)
{
    int row, col;
    register int i, j, ld1, ld2;
    register double complex tmp;
    register double complex* matrix1,* matrix2;
    
    matrix1 = this->matrix;
    matrix2 = out_matrix->matrix;
    row = this->row;
    col = this->col;
    ld1 = this->ld;
    ld2 = out_matrix->ld;
    if(row == col)
    {
        for(j=col;j-->0;)
        {
            for(i=j+1;i-->0;)
            {
                tmp = *(matrix1+j*ld1+i);
                *(matrix2+j*ld2+i) = *(matrix1+i*ld1+j);
                *(matrix2+i*ld2+j) = tmp;
            }
        }
    }
    else
    {
        for(j=col;j-->0;)
        {
            for(i=row;i-->0;)
            {
                *(matrix2+j*ld2+i) = *(matrix1+i*ld1+j);
            }
        }
    }
    
    return out_matrix;
}

MatrixDC* MatrixDC_CT(MatrixDC* this, MatrixDC* out_matrix)
{
    int row, col;
    register int i, j, ld1, ld2;
    register double complex tmp;
    register double complex* matrix1,* matrix2;
    
    matrix1 = this->matrix;
    matrix2 = out_matrix->matrix;
    row = this->row;
    col = this->col;
    ld1 = this->ld;
    ld2 = out_matrix->ld;
    if(row == col)
    {
        for(j=col;j-->0;)
        {
            for(i=j+1;i-->0;)
            {
                tmp = *(matrix1+j*ld1+i);
                *(matrix2+j*ld2+i) = conj(*(matrix1+i*ld1+j));
                *(matrix2+i*ld2+j) = conj(tmp);
            }
        }
    }
    else
    {
        for(j=col;j-->0;)
        {
            for(i=row;i-->0;)
            {
                *(matrix2+j*ld2+i) = conj(*(matrix1+i*ld1+j));
            }
        }
    }
    
    return out_matrix;    
}

double complex MatrixDC_Det(MatrixDC* this)
{
    int row, col, info;
    int *ipiv;
    register int i;
    double complex det;
    MatrixDC* tmp_matrix;
    double complex* matrix;
     
    row = this->row;
    col = this->col;
    ipiv = (int*)malloc(sizeof(int)*((row>col)?col:row));
    det = 1.0;
    tmp_matrix = MatrixDC_New(row, col);
    MatrixDC_Copy(this, tmp_matrix);
    matrix = tmp_matrix->matrix;
    zgetrf_(&row, &col, matrix, &row, ipiv, &info);
    for(i=0;i<row;i++)
    {
        det *= *(matrix+i*col+i);
    }
    MatrixDC_Del(tmp_matrix);
    free(ipiv);
    
    return det;
}

MatrixDC* MatrixDC_Inv(MatrixDC* this, MatrixDC* out_matrix)
{
    int row, col, info, ld;
    int* ipiv;
    double complex* work,* matrix;
    
    row = this->row;
    col = this->col;
    ld = out_matrix->ld;
    ipiv = (int*)malloc(sizeof(int)*((row>col)?col:row));
    work = (double complex*)malloc(sizeof(double complex)*col*row);
    matrix = out_matrix->matrix;

    MatrixDC_Copy(this, out_matrix);
    zgetrf_(&row, &col, matrix, &ld, ipiv, &info);
    zgetri_(&row, matrix, &ld, ipiv, work, &row, &info);
    free(ipiv);
    free(work);
    
    return out_matrix;
}

MatrixDC* MatrixDC_PInv(MatrixDC* this, MatrixDC* out_matrix)
{
    register int i, j;
    int row, col, ld1, ld2, number_of_singular, lwork, info;
    double rwork_tmp;
    double* s,* rwork;
    double complex* u,* vt,* work,* matrix;
    double complex work_size;
    MatrixDC* tmp_matrix;

    row = this->row;
    col = this->col;
    ld1 = this->ld;
    ld2 = out_matrix->ld;
    work = &work_size;
    matrix = out_matrix->matrix;
    lwork = -1;
    rwork = &rwork_tmp;
    info = 0;
    number_of_singular = row<col ? row : col;
    
    s = (double*)calloc(col, sizeof(double));
    u = (double complex*)malloc(row*row * sizeof(double complex));
    vt = (double complex*)malloc(col*col* sizeof(double complex));

    MatrixDC_T(this, out_matrix);
    zgesvd_("A", "A", &row, &col, matrix, &row, s, u, &row, vt, &col, work, &lwork, rwork, &info);
    lwork = work_size;
    work = (double complex*)malloc(lwork*sizeof(double complex));
    rwork = (double*)malloc(lwork*sizeof(double));
    zgesvd_("A", "A", &row, &col, matrix, &row, s, u, &row, vt, &col, work, &lwork, rwork, &info);
    free(work);
    free(rwork);
    tmp_matrix = MatrixDC_Wrap(u, row, row, row);
    MatrixDC_CT(tmp_matrix, tmp_matrix);
    MatrixDC_ReWrap(tmp_matrix, vt, col, col, col);
    MatrixDC_CT(tmp_matrix, tmp_matrix);
    for(i=0;i<col;i++)
    {
        register double tmp;

        if(*(s+i)!=0)
        {
            tmp = 1 / *(s+i);
        }
        else
        {
            tmp = 0;
        }
        for(j=0;j<row;j++)
        {
            *(matrix+i*ld2+j) = tmp * *(u+i*row+j);
        }
    }
    MatrixDC_MulMatrixDC(tmp_matrix, out_matrix, out_matrix);
    free(s);
    free(u);
    free(vt);

    return out_matrix;
}

MatrixDC* MatrixDC_AddMatrixDC(MatrixDC* this, MatrixDC* in_matrix, MatrixDC* out_matrix)
{
    register int row, col, i, j, ld1, ld2, ld3;
    register double complex* matrix1,* matrix2,* matrix3;

    row = this->row;
    col = this->col;
    ld1 = this->ld;
    ld2 = in_matrix->ld;
    ld3 = out_matrix->ld;
    matrix1 = this->matrix;
    matrix2 = in_matrix->matrix;
    matrix3 = out_matrix->matrix;
    for(i=0;i<row;i++)
    {
        for(j=0;j<col;j++)
        {
            *(matrix3+i*ld3+j) = *(matrix1+i*ld1+j) + *(matrix2+i*ld2+j);
        }
    }

    return out_matrix;
}

MatrixDC* MatrixDC_SubMatrixDC(MatrixDC* this, MatrixDC* in_matrix, MatrixDC* out_matrix)
{
    register int row, col, i, j, ld1, ld2, ld3;
    register double complex* matrix1,* matrix2,* matrix3;

    row = this->row;
    col = this->col;
    ld1 = this->ld;
    ld2 = in_matrix->ld;
    ld3 = out_matrix->ld;
    matrix1 = this->matrix;
    matrix2 = in_matrix->matrix;
    matrix3 = out_matrix->matrix;
    for(i=0;i<row;i++)
    {
        for(j=0;j<col;j++)
        {
            *(matrix3+i*ld3+j) = *(matrix1+i*ld1+j) - *(matrix2+i*ld2+j);
        }
    }

    return out_matrix;
}

MatrixDC* MatrixDC_MulNumDC(MatrixDC* this, double complex num, MatrixDC* out_matrix)
{
    register int i, j, ld1, ld2, row, col;
    register double complex* matrix1,* matrix2;

    row = out_matrix->row;
    col = out_matrix->col;
    ld1 = this->ld;
    ld2 = out_matrix->ld;
    matrix1 = this->matrix;
    matrix2 = out_matrix->matrix;
    for(i=0;i<row;i++)
    {
        for(j=0;j<col;j++)
        {
            *(matrix2+i*ld2+j) = num * *(matrix1+i*ld1+j);
        }
    }
    
    return out_matrix;
}

VectorDC* MatrixDC_MulVectorDC(MatrixDC* this, VectorDC* in_vector, VectorDC* out_vector)
{
    int in_size, out_size;
    register int i, l, ld, inc1, inc2;
    register double complex result;
    register double complex* matrix,* vector1,* vector2;
    VectorDC* tmp_vector;
    
    out_size = this->row;
    in_size = this->col;
    ld = this->ld;
    tmp_vector = VectorDC_New(in_size);
    inc1 = VectorDC_GetInc(tmp_vector);
    inc2 = VectorDC_GetInc(out_vector);
    matrix = this->matrix;
    VectorDC_Copy(in_vector, tmp_vector);
    vector1 = VectorDC_GetVector(tmp_vector);
    vector2 = VectorDC_GetVector(out_vector);
    
    result = 0;
    
    for(i=0;i<out_size;i++)
    {
        result = 0;
        for(l=0;l<in_size;l++)
        {
            result += *(matrix+i*ld+l) * *(vector1+l*inc1);
        }
        *(vector2+inc2*i) = result;
    }
    VectorDC_Del(tmp_vector);
    
    return out_vector;
}

VectorDC* MatrixDC_TMulVectorDC(MatrixDC* this, VectorDC* in_vector, VectorDC* out_vector)
{
    int row, col;
    register int i, j, ld, inc1, inc2;
    register double complex tmp;
    register double complex* matrix,* vector1,* vector2;
    VectorDC* tmp_vector;

    row = this->row;
    col = this->col;
    ld = this->ld;
    tmp_vector = VectorDC_New(row);
    inc1 = VectorDC_GetInc(tmp_vector);
    inc2 = VectorDC_GetInc(out_vector);
    matrix = this->matrix;
    vector1 = VectorDC_GetVector(tmp_vector);
    vector2 = VectorDC_GetVector(out_vector);
    VectorDC_Copy(in_vector, tmp_vector);
    
    for(j=col;j-->0;)
    {
        *(vector2+j*inc2) = 0;
    }
    
    for(i=row;i-->0;)
    {
        tmp = *(vector1+i*inc1);
        for(j=col;j-->0;)
        {
            *(vector2+j*inc2) += *(matrix+i*ld+j) * tmp;
        }
    }
    VectorDC_Del(tmp_vector);
    
    return out_vector;
}

VectorDC* MatrixDC_CMulVectorDC(MatrixDC* this, VectorDC* in_vector, VectorDC* out_vector)
{
    int row, col;
    register int i, j, ld, inc1, inc2;
    register double complex tmp;
    register double complex* matrix,* vector1,* vector2;
    VectorDC* tmp_vector;
    
    row = this->row;
    col = this->col;
    ld = this->ld;
    tmp_vector = VectorDC_New(row);
    inc1 = VectorDC_GetInc(tmp_vector);
    inc2 = VectorDC_GetInc(out_vector);
    matrix = this->matrix;    
    vector1 = VectorDC_GetVector(tmp_vector);    
    vector2 = VectorDC_GetVector(out_vector);
    VectorDC_Copy(in_vector, tmp_vector);

    for(j=col;j-->0;)
    {
        *(vector2+inc2*j) = 0;
    }
    
    for(i=row;i-->0;)
    {
        tmp = *(vector1+i*inc1);
        for(j=col;j-->0;)
        {
            *(vector2+inc2*j) += conj(*(vector1+i*ld+j))*tmp;
        }
    }
    VectorDC_Del(tmp_vector);
    
    return out_vector;
}   

MatrixDC* MatrixDC_TMulSelf(MatrixDC* this, MatrixDC* out_matrix)
{
    int row, col;
    register int i, j, l, ld1, ld2;
    double complex tmp;
    register double complex* matrix1,* matrix2;

    row = this->row;
    col = this->col;
    ld1 = this->ld;
    ld2 = out_matrix->ld;
    matrix1 = this->matrix;
    matrix2 = out_matrix->matrix;
    
    for(i=col;i-->0;)
    {
        for(j=col;j-->0;)
        {
            *(matrix2+i*ld2+j) = 0;
        }
    }
    
    for(l=row;l-->0;)
    {
        for(i=col;i-->0;)
        {
            tmp = *(matrix1+l*ld1+i);
            for(j=i+1;j-->0;)
            {
                *(matrix2+i*ld2+j) += tmp * *(matrix1+l*ld1+j);
            }
        }
    }

    for(i=col;i-->0;)
    {
        for(j=i+1;j<col;j++)
        {
            *(matrix2+i*ld2+j) = *(matrix2+j*ld2+i);
        }
    }
    
    return out_matrix;
}

MatrixDC* MatrixDC_O2Mul(MatrixDC* this, MatrixDC* in_matrix, MatrixDC* out_matrix)
{
    register int ld1, ld2, ld3;
    register double complex m11, m12, m13, m14, m21, m22, m23, m24;
    register double complex* matrix1,* matrix2,* matrix3;
    
    ld1 = this->ld;
    ld2 = in_matrix->ld;
    ld3 = out_matrix->ld;
    matrix1 = this->matrix;
    matrix2 = in_matrix->matrix;
    matrix3 = out_matrix->matrix;
    m11 = *(matrix1);
    m12 = *(matrix1+1);
    m13 = *(matrix1+ld1);
    m14 = *(matrix1+ld1+1);
    m21 = *(matrix2);
    m22 = *(matrix2+1);
    m23 = *(matrix2+ld2);
    m24 = *(matrix2+ld2+1);
    
    *(matrix3) = m11 * m21 + m12 * m23;
    *(matrix3+1) = m11 * m22 + m12 * m24;
    *(matrix3+ld3) = m13 * m21 + m14 * m23;
    *(matrix3+ld3+1) = m13 * m22 + m14 * m24;
    
    return out_matrix;
}

MatrixDC* MatrixDC_DiagLeftMul(MatrixDC* this, VectorDC* diag_matrix, MatrixDC* out_matrix)
{
    int row, col;
    register int i, j, ld1, ld2, inc;
    register double complex tmp;
    register double complex* matrix1,* matrix2,* diag;

    row = this->row;
    col = this->col;
    ld1 = this->ld;
    ld2 = this->ld;
    matrix1 = this->matrix;
    matrix2 = out_matrix->matrix;
    inc = VectorDC_GetInc(diag_matrix);
    diag = VectorDC_GetVector(diag_matrix);
    for(i=0;i<row;i++)
    {
        tmp = *(diag+inc*i);
        for(j=0;j<col;j++)
        {
            *(matrix2+i*ld2+j) = *(matrix1+i*ld1+j) * tmp;
        }
    }

    return out_matrix;
}

MatrixDC* MatrixDC_DiagRightMul(MatrixDC* this, VectorDC* diag_matrix, MatrixDC* out_matrix)
{
    int row, col;
    register int i, j, ld1, ld2, inc;
    register double complex* matrix1,* matrix2,* diag;

    row = this->row;
    col = this->col;
    ld1 = this->ld;
    ld2 = this->ld;
    matrix1 = this->matrix;
    matrix2 = out_matrix->matrix;
    inc = VectorDC_GetInc(diag_matrix);
    diag = VectorDC_GetVector(diag_matrix);
    for(i=0;i<row;i++)
    {
        for(j=0;j<col;j++)
        {
            
            *(matrix2+i*ld2+j) = *(matrix1+i*ld1+j) * *(diag+inc*j);
        }
    }

    return out_matrix;
}

MatrixDC* MatrixDC_MulMatrixDC(MatrixDC* this, MatrixDC* in_matrix, MatrixDC* out_matrix)
{
    register int i, j;
    double complex alpha, beta;
    int m, n, k, lda, ldb, ldc;
    register int ld1, ld2;
    double complex* matrix1,* matrix2,* matrix3,* tmp_matrix;
    
    alpha = 1;
    beta = 0;
    lda = this->ld;
    ldb = in_matrix->ld;
    ldc = out_matrix->col;
    ld1 = out_matrix->ld;
    ld2 = ldc;
    m = this->row;
    n = in_matrix->col;
    k = in_matrix->row;
    matrix1 = this->matrix;
    matrix2 = in_matrix->matrix;
    matrix3 = out_matrix->matrix;
    tmp_matrix = (double complex*)malloc(sizeof(double complex)*out_matrix->row*out_matrix->col);
    zgemm_("T","T",&m,&n,&k,&alpha,matrix1,&lda,matrix2,&ldb,&beta,tmp_matrix,&ldc);
    for(i=out_matrix->row;i-->0;)
    {
        for(j=out_matrix->col;j-->0;)
        {
            *(matrix3+i*ld1+j)=tmp_matrix[j*ld2+i];
        }
    }
    free(tmp_matrix);
    
    return out_matrix;
}

VectorDC* MatrixDC_LinearSolver(MatrixDC* this, VectorDC* in_vector, VectorDC* out_vector)
{
    int row, col, info, nrhs;
    int* ipiv;
    double complex* matrix,* vector;
    MatrixDC* tmp_matrix;
    
    nrhs = 1;
    row = this->row;
    col = this->col;
    ipiv = (int*)malloc(sizeof(int)*((row>col)?col:row));
    tmp_matrix = MatrixDC_New(col, row);
    matrix = tmp_matrix->matrix;
    vector = VectorDC_GetVector(out_vector);
    
    MatrixDC_T(this, tmp_matrix);
    VectorDC_Copy(in_vector, out_vector);
    zgetrf_(&row, &col, matrix, &row, ipiv, &info);
    zgetrs_("N", &row, &nrhs, matrix, &row, ipiv, vector, &row, &info);
    free(ipiv);
    MatrixDC_Del(tmp_matrix);
    
    return out_vector;
}

void MatrixDC_EigenEquation(MatrixDC* this, VectorDC* out_vector, MatrixDC* out_matrix)
{
    int n, ldvl, ldvr, lwork, info;
    double* rwork;
    double complex wkopt;
    double complex* work,* matrix1,* matrix2,* vector;
    MatrixDC* tmp_matrix1;
    
    n = this->row;
    ldvl = n;
    ldvr = out_matrix->ld;
    lwork = -1;
    rwork = (double*)malloc(sizeof(double)*(n<<1));
    tmp_matrix1 = MatrixDC_New(n,n);
    matrix1 = tmp_matrix1->matrix;
    matrix2 = out_matrix->matrix;
    vector = VectorDC_GetVector(out_vector);
    MatrixDC_T(this, tmp_matrix1);
    zgeev_("N","V",&n,matrix1,&n,vector,NULL,&ldvl,matrix2,&ldvr,&wkopt,&lwork,rwork,&info);
    lwork = (int)creal(wkopt);
    work = (double complex*)malloc(sizeof(double complex)*lwork);
    zgeev_("N","V",&n,matrix1,&n,vector,NULL,&ldvl,matrix2,&ldvr, work, &lwork, rwork, &info);
    MatrixDC_T(out_matrix, out_matrix);
    MatrixDC_Del(tmp_matrix1);
    free(work);
    free(rwork);
}

int MatrixDC_Equal(MatrixDC* this, MatrixDC* matrix)
{
    int row, col;
    register int i, j, ld1, ld2;
    register double complex* matrix1,* matrix2;
    
    row = this->row;
    col = this->col;
    ld1 = this->ld;
    ld2 = matrix->ld;
    matrix1 = this->matrix;
    matrix2 = matrix->matrix;
    
    if((row == matrix->row)&&(col == matrix->col))
    {
        for(i=0;i<row;i++)
        {
            for(j=0;j<col;j++)
            {
                if(*(matrix1+i*ld1+j)!= *(matrix2+i*ld2+j))
                {
                    return 0;
                }
            }
        }
        
        return 1;
    }
    
    return 0;        
}

void* _MatrixDC_New(void* implementor)
{
    MatrixDC* matrix;
    
    matrix = (MatrixDC*)implementor;
    
    return MatrixDC_New(matrix->row, matrix->col);
}

void* _MatrixDC_Del(void* implementor)
{
    MatrixDC* matrix;
    
    matrix = (MatrixDC*)implementor;
    
    return MatrixDC_Del(matrix);
}

int _MatrixDC_Copy(void* implementor, void* object)
{
    MatrixDC* matrix1,* matrix2;
    
    matrix1 = (MatrixDC*)implementor;
    matrix2 = (MatrixDC*)object;
    MatrixDC_Copy(matrix2, matrix1);
    
    return 1;
}

int _MatrixDC_Equal(void* implementor, void* object)
{
    MatrixDC* matrix1,* matrix2;
   
    matrix1 = (MatrixDC*)implementor;
    matrix2 = (MatrixDC*)object;
   
    return MatrixDC_Equal(matrix1, matrix2);
}
