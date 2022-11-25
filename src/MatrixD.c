#include <stdio.h>
#include <stdlib.h>
#include "Blas.h"
#include "Lapack.h"
#include "IObject.h"
#include "Class.h"
#include "Object.h"
#include "VectorD.h"
#include "MatrixD.h"

Class MatrixD_Class = {OBJECT};

struct _MatrixD
{
    Class* class;
    IObject* iObject;
    double* matrix;
    int row;
    int col;
    int ld;
};

void* _MatrixD_New(void* implementor);
void* _MatrixD_Del(void* implementor);
int _MatrixD_Copy(void* implementor, void* object);
int _MatrixD_Equal(void* implementor, void* object);

MatrixD* MatrixD_New(int row, int col)
{
    MatrixD* matrix;
    
    matrix = (MatrixD*)malloc(sizeof(MatrixD)+sizeof(IObject)+sizeof(double)*row*col);
    matrix->class = MATRIXD;   
    matrix->iObject = (IObject*)(matrix+1);
    matrix->iObject->implementor = matrix;
    matrix->iObject->New = _MatrixD_New;
    matrix->iObject->Del = _MatrixD_Del;
    matrix->iObject->Copy = _MatrixD_Copy;
    matrix->iObject->Equal = _MatrixD_Equal;
    matrix->row = row;
    matrix->col = col;
    matrix->ld = col;
    matrix->matrix = (double*)(matrix->iObject+1); 
    
    return matrix;
}

MatrixD* MatrixD_Del(MatrixD* this)
{    
    free(this);
    
    return NULL;
}

MatrixD* MatrixD_Wrap(double* work, int row, int col, int ld)
{
    MatrixD* matrix;
    
    matrix = (MatrixD*)malloc(sizeof(MatrixD)+sizeof(IObject));
    matrix->class = MATRIXD;
    matrix->iObject = (IObject*)(matrix+1);
    matrix->iObject->implementor = matrix;
    matrix->iObject->New = _MatrixD_New;
    matrix->iObject->Del = _MatrixD_Del;
    matrix->iObject->Copy = _MatrixD_Copy;
    matrix->iObject->Equal = _MatrixD_Equal;
    matrix->row = row;
    matrix->col = col;
    matrix->ld = ld;
    matrix->matrix = work;    
    
    return matrix;
}

MatrixD* MatrixD_UnWrap(MatrixD* this)
{  
    free(this);
    
    return NULL;
}

IObject* MatrixD_GetIObject(MatrixD* this)
{
    return this->iObject;
}

int MatrixD_GetRow(MatrixD* this)
{   
    return this->row; 
}

int MatrixD_GetCol(MatrixD* this)
{    
    return this->col; 
}

int MatrixD_GetLd(MatrixD* this)
{    
    return this->ld; 
}

double* MatrixD_GetMatrix(MatrixD* this)
{
    return this->matrix;
}

double MatrixD_GetElement(MatrixD* this, int row_index, int col_index)
{    
    return *(this->matrix+row_index*this->ld+col_index);
}

MatrixD* MatrixD_SetElement(MatrixD* this, int row_index, int col_index, double value)
{
    *(this->matrix+row_index*this->ld+col_index) = value;
 
    return this;
}

void MatrixD_Show(MatrixD* this)
{
    register int i, j, row, col, ld;
    register double* matrix;

    row = this->row;
    col = this->col;
    ld = this->ld;
    matrix = this->matrix;
    printf("MatrixD_Show:\n");
    for(i=0;i<row;i++)
    {
        printf(" [");
        for(j=0;j<col;j++)
        {
            printf("%lf, ", *(matrix+i*ld+col));
        }
        printf("]\n");
    }
    printf(" row:%d \n col:%d \n ld:%d \n\n", row, col, ld);    
}

MatrixD* MatrixD_Set(MatrixD* this, double diag, double off_diag)
{
    int row, col, ld;
    register int i, j;
    register double* matrix;
    
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

MatrixD* MatrixD_Copy(MatrixD* this, MatrixD* out_matrix)
{
    register int i, j, row, col, ld1, ld2;
    double* matrix1,* matrix2;

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

MatrixD* MatrixD_T(MatrixD* this, MatrixD* out_matrix)
{
    register int i, j, row, col, ld1, ld2;
    double* matrix1,* matrix2;
    
    row = this->row;
    col = out_matrix->col;
    ld1 = this->ld;
    ld2 = out_matrix->ld;
    matrix1 = this->matrix;
    matrix2 = out_matrix->matrix;
    for(i=0;i<row;i++)
    {
        for(j=0; j<col;j++)
        {
            *(matrix2+j*ld2+i) = *(matrix1+i*ld1+j);
        }
    }
    
    return out_matrix;
}

double MatrixD_Det(MatrixD* this)
{
    int row, col, info;
    int* ipiv;
    register int i;
    double* matrix;
    register double det;
    MatrixD* tmp_matrix;
    
    row = this->row;
    col = this->ld;
    ipiv = (int*)malloc(sizeof(int)*((row > col)? col : row));
    det = 1.0;
    tmp_matrix = MatrixD_New(row, col);
    tmp_matrix = MatrixD_Copy(this, tmp_matrix);
    matrix = tmp_matrix->matrix;
    dgetrf_(&row, &col, matrix, &row, ipiv, &info);
    for(i=0;i<row;i++)
    {
        det *= *(matrix+i*row+i);
    }
    MatrixD_Del(tmp_matrix);
    free(ipiv);
    
    return det;
}

MatrixD* MatrixD_Inv(MatrixD* this, MatrixD* out_matrix)
{
    int row, col, info, ld;
    int* ipiv;
    double* work,* matrix;
    
    row = this->row;
    col = this->col;
    ld = this->ld;
    ipiv = (int*)malloc(sizeof(int)*((row > col) ? col : row));
    work = (double*)malloc(sizeof(double)*col*row);
    matrix = out_matrix->matrix;
    
    MatrixD_T(this, out_matrix);
    dgetrf_(&row, &col, matrix, &ld, ipiv, &info);
    dgetri_(&row, matrix, &ld, ipiv, work, &row, &info);
    free(ipiv);
    free(work);
    
    return out_matrix;
}

VectorD* MatrixD_MulVectorD(MatrixD* this, VectorD* in_vector, VectorD* out_vector)
{
    register int i, l, ld, inc1, inc2, in_size, out_size;
    register double result;
    register double* matrix,* vector1,* vector2;
    VectorD* tmp_vector;
    
    in_size = this->row;
    out_size = this->col;
    ld = this->ld;
    tmp_vector = VectorD_New(in_size);
    inc1 = VectorD_GetInc(tmp_vector);
    inc2 = VectorD_GetInc(out_vector);
    matrix = this->matrix;    
    VectorD_Copy(in_vector, tmp_vector);
    vector1 = VectorD_GetVector(tmp_vector);
    vector2 = VectorD_GetVector(out_vector);
    result = 0;
    
    for(i=0;i<out_size;i++)
    {
        result = 0;
        for(l=0; l<in_size; l++)
        {
            result += *(matrix+i*ld+l) * *(vector1+l*inc1);
        }
        *(vector2+i*inc2) = result;
    }
    VectorD_Del(tmp_vector);
    
    return out_vector;
}

VectorD* MatrixD_TMulVectorD(MatrixD* this, VectorD* in_vector, VectorD* out_vector)
{
    register int i, j, row, col, ld, inc1, inc2;
    register double tmp;
    register double* matrix,* vector1,* vector2;
    VectorD* tmp_vector;
    
    row = this->row;
    col = this->col;
    ld = this->ld;
    matrix = this->matrix;
    tmp_vector = VectorD_New(row);
    inc1 = VectorD_GetInc(in_vector);
    inc2 = VectorD_GetInc(tmp_vector);
    vector1 = VectorD_GetVector(in_vector);
    vector2 = VectorD_GetVector(tmp_vector);
    
    for(j=col;j-->0;)
    {
        *(vector2+inc2*j) = 0;
    }
    
    for(i=row;i-->0;)
    {
        tmp = *(vector1+inc1*i);
        for(j=col;j-->0;)
        {
            *(vector2+j*inc2) += *(matrix+i*ld+j) * tmp;
        }
    }
    VectorD_Copy(tmp_vector, out_vector);
    VectorD_Del(tmp_vector);
    
    return out_vector;
}

MatrixD* MatrixD_TMulSelf(MatrixD* this, MatrixD* out_matrix)
{
    register int i, j, l, row, col, ld1, ld2;
    register double tmp;
    double* matrix1,* matrix2;
    
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

MatrixD* MatrixD_O2Mul(MatrixD* this, MatrixD* in_matrix, MatrixD* out_matrix)
{
    register double m11, m12, m13, m14, m21, m22, m23, m24;
    int ld1, ld2, ld3;
    double* matrix1,* matrix2,* matrix3;

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

MatrixD* MatrixD_MulMatrixD(MatrixD* this, MatrixD* in_matrix, MatrixD* out_matrix)
{
    double alpha, beta;
    int m, n, k, lda, ldb, ldc;
    double* matrix1,* matrix2,* matrix3;

    alpha = 1;
    beta = 0;    
    lda = this->ld;
    ldb = in_matrix->ld;
    ldc = out_matrix->ld;    
    m = this->col;
    n = in_matrix->row;
    k = this->row;
    matrix1 = this->matrix;
    matrix2 = in_matrix->matrix;
    matrix3 = out_matrix->matrix;
    dgemm_("T","T",&m,&n,&k,&alpha,matrix1,&lda,matrix2,&ldb,&beta,matrix3,&ldc);

    return out_matrix;
}

double MatrixD_MaxDiagElem(MatrixD* this)
{
    register int i, size, ld;
    register double max;
    double* matrix;
    
    size = this->row;
    ld = this->ld;
    matrix = this->matrix;
    max = *(matrix);
    for(i=1;i<size;++i)
    {
        max = (max > *(matrix+i*ld+i)) ? max : *(matrix+i*ld+i);
    }

    return max;
}

VectorD* MatrixD_LinearSolver(MatrixD* this, VectorD* in_vector, VectorD* out_vector)
{
    int row, col, info, nrhs;
    int* ipiv;
    MatrixD* tmp_matrix;
    double* matrix,* vector;
    
    nrhs = 1;
    row = this->row;
    col = this->col;
    ipiv = (int*)malloc(sizeof(int)*((row>col)?col:row));
    tmp_matrix = MatrixD_New(col, row);
    matrix = tmp_matrix->matrix;
    vector = VectorD_GetVector(out_vector);
    
    MatrixD_T(this, tmp_matrix);
    VectorD_Copy(in_vector, out_vector);
    dgetrf_(&row, &col, matrix, &row, ipiv, &info);
    dgetrs_("N", &row, &nrhs, matrix, &row, ipiv, vector, &row, &info);
    free(ipiv);
    MatrixD_Del(tmp_matrix);
    
    return out_vector;
}

void MatrixD_EigenEquation(MatrixD* this, VectorD* out_vector, MatrixD* out_matrix)
{
    int n, ldvl, ldvr, lwork, info;
    double* rwork,* work,* matrix1,* matrix2,* vector;
    double wkopt;
    MatrixD* tmp_matrix1,* tmp_matrix2;
    
    n = this->row;
    ldvl = n;
    ldvr = this->ld;
    lwork = -1;
    rwork = (double*)malloc(sizeof(double)*(n<<1));
    tmp_matrix1 = MatrixD_New(n,n);
    tmp_matrix2 = MatrixD_New(n,n);
    MatrixD_T(this, tmp_matrix1);
    matrix1 = tmp_matrix1->matrix;
    matrix2 = tmp_matrix2->matrix;
    vector = VectorD_GetVector(out_vector);
    dgeev_("N","V",&n,matrix1,&n,vector,NULL,&ldvl,matrix2,&ldvr, &wkopt, &lwork, rwork, &info);
    lwork = (int) wkopt;
    work = (double*)malloc(sizeof(double)*lwork);
    dgeev_("N","V",&n,matrix1,&n,vector,NULL,&ldvl,matrix2,&ldvr, work, &lwork, rwork, &info);
    MatrixD_Del(tmp_matrix1);
    MatrixD_Del(tmp_matrix2);
    free(work);
    free(rwork);    
}

void* _MatrixD_New(void* implementor)
{
    MatrixD* matrix;
    
    matrix = (MatrixD*)implementor;
    
    return MatrixD_New(matrix->row, matrix->col);
}

void* _MatrixD_Del(void* implementor)
{
    MatrixD* matrix;

    matrix = (MatrixD*)implementor;
    
    return MatrixD_Del(matrix);
}

int _MatrixD_Copy(void* implementor, void* object)
{
    MatrixD* matrix1,* matrix2;
    
    matrix1 = (MatrixD*)implementor;
    matrix2 = (MatrixD*)object;
    MatrixD_Copy(object, implementor);
    
    return 1;
}

int _MatrixD_Equal(void* implementor, void* object)
{
    register int i, j, ld1, ld2;
    int row, col;
    MatrixD* matrix_self,* matrix_cmp;
    double* matrix1,* matrix2;
   
    matrix_self = (MatrixD*)implementor;
    matrix_cmp = (MatrixD*)object;
    row = matrix_self->row;
    col = matrix_self->col;
    ld1 = matrix_self->ld;
    ld2 = matrix_cmp->ld;
    matrix1 = matrix_self->matrix;
    matrix2 = matrix_cmp->matrix;
    
    if((row == matrix_cmp->row)&&(col == matrix_cmp->col))
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
