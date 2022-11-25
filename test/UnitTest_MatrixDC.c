#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "Class.h"
#include "IObject.h"
#include "IUnitTest.h"
#include "Lapack.h"
#include "MatrixDC.h"
#include "VectorDC.h"
#include "UnitTest_MatrixDC.h"

#define FABS(x) ((x)>0 ? (x) : -(x))

Class UnitTest_MatrixDC_Class = {CLASS};

struct _UnitTest_MatrixDC
{
    Class* class;
    IUnitTest* iUnitTest;
};

struct _MatrixDC
{
    Class* class;
    IObject* iObject;
    double complex* matrix;
    int row;
    int col;
    int ld;
};

void _UnitTest_MatrixDC_TestAll(void* implementor);

static VectorDC* _Stub_MatrixDC_MulVectorDC(MatrixDC* this, VectorDC* in_vector, VectorDC* out_vector)
{
    int in_size, out_size;
    register int i, l, ld, inc1, inc2;
    register double complex result;
    register double complex* matrix,* vector1,* vector2;
    
    out_size = this->row;
    in_size = this->col;
    ld = this->ld;
    inc1 = VectorDC_GetInc(in_vector);
    inc2 = VectorDC_GetInc(out_vector);
    matrix = this->matrix;
    vector1 = VectorDC_GetVector(in_vector);
    vector2 = VectorDC_GetVector(out_vector);
    result = 0;
    
    for(i=0;i<out_size;i++)
    {
        result = 0;
        for(l=0;l<in_size;l++)
        {
            result += *(matrix+i*ld+l) * *(vector1+l*inc1);
            printf("%g,%g\n", result);
        }
        *(vector2+inc2*i) = result;
    }
    
    return out_vector;
}

static MatrixDC* _Stub_MatrixDC_PInv(MatrixDC* this, MatrixDC* out_matrix)
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
    //MatrixDC_Show(out_matrix);
    zgesvd_("A", "A", &row, &col, matrix, &row, s, u, &row, vt, &col, work, &lwork, rwork, &info);
    lwork = (int)work_size;
    //printf("%d\n",lwork);
    //printf("%d\n",info);
    work = (double complex*)malloc(lwork*sizeof(double complex));
    //printf("%p\n",work);

    rwork = (double*)malloc(lwork*sizeof(double));
    zgesvd_("A", "A", &row, &col, matrix, &row, s, u, &row, vt, &col, work, &lwork, rwork, &info);
    free(work);
    free(rwork);
    tmp_matrix = MatrixDC_Wrap(u, row, row, row);
    MatrixDC_CT(tmp_matrix, tmp_matrix);
    //MatrixDC_Show(tmp_matrix);
    MatrixDC_ReWrap(tmp_matrix, vt, col, col, col);
    MatrixDC_CT(tmp_matrix, tmp_matrix);
    //MatrixDC_Show(tmp_matrix);
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
            //printf("OK!\n");
        }
        //printf("%g\n",*(s+i));
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

static VectorDC* _Stub_MatrixDC_LinearSolver(MatrixDC* this, VectorDC* in_vector, VectorDC* out_vector)
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
    zgesv_(&row, &nrhs, matrix, &row, ipiv, vector, &col, &info);
    free(ipiv);
    MatrixDC_Del(tmp_matrix);
    
    return out_vector;
}

static void Test_MatrixDC_Big()
{
    MatrixDC* matrix,* tmp_matrix;

    matrix = MatrixDC_New(3368, 3368);
    tmp_matrix = MatrixDC_New(3368, 3368);
    MatrixDC_MulMatrixDC(matrix, tmp_matrix, matrix);
    MatrixDC_Del(matrix);
    MatrixDC_Del(tmp_matrix);
    printf("Test_MatrixDC_Big Success.\n");
}

static void Test_MatrixDC_IObject_New_Copy()
{
    MatrixDC* matrix,* tmp_matrix;
    IObject* iObject;
    
    matrix = MatrixDC_New(2,2);
    iObject = MatrixDC_GetIObject(matrix);
    tmp_matrix = (MatrixDC*)IObject_New(iObject);  
    iObject->implementor = tmp_matrix;   
    MatrixDC_Set(matrix, 1, 1);
    IObject_Copy(iObject, matrix);   
    iObject->implementor = matrix;
    assert(MatrixDC_GetElement(tmp_matrix, 0, 0) == 1);
    assert(MatrixDC_GetElement(tmp_matrix, 0, 1) == 1);
    assert(MatrixDC_GetElement(tmp_matrix, 1, 0) == 1);
    assert(MatrixDC_GetElement(tmp_matrix, 1, 1) == 1);
    MatrixDC_Del(matrix);
    MatrixDC_Del(tmp_matrix);
    printf("Test_MatrixDC_IObject Success.\n");
}

static void Test_MatrixDC_Copy()
{
    MatrixDC* matrix,* tmp_matrix;
    
    matrix = MatrixDC_New(2,2);
    tmp_matrix = MatrixDC_New(2,2);
    MatrixDC_Set(matrix, 1, 1);
    MatrixDC_Copy(matrix, tmp_matrix);
    assert(MatrixDC_GetElement(tmp_matrix, 0, 0) == 1);
    assert(MatrixDC_GetElement(tmp_matrix, 0, 1) == 1);
    assert(MatrixDC_GetElement(tmp_matrix, 1, 0) == 1);
    assert(MatrixDC_GetElement(tmp_matrix, 1, 1) == 1);    
    MatrixDC_Del(matrix);
    MatrixDC_Del(tmp_matrix);
    printf("Test_MatrixDC_Copy Success. \n");
}

static void Test_MatrixDC_Set()
{
    MatrixDC* matrix;
    
    matrix = MatrixDC_New(2,2);
    MatrixDC_Set(matrix, 1, 1);
    assert(MatrixDC_GetElement(matrix, 0, 0) == 1);
    assert(MatrixDC_GetElement(matrix, 0, 1) == 1);
    assert(MatrixDC_GetElement(matrix, 1, 0) == 1);
    assert(MatrixDC_GetElement(matrix, 1, 1) == 1);
    MatrixDC_Del(matrix);
    printf("Test_MatrixDC_Set Success. \n");
}

static void Test_MatrixDC_Inv()
{
    MatrixDC* matrix;

    matrix = MatrixDC_New(3,3);
    MatrixDC_Set(matrix, 1, 1);
    MatrixDC_SetElement(matrix, 0, 2, 0);
    MatrixDC_SetElement(matrix, 1, 1, 0);
    MatrixDC_SetElement(matrix, 2, 0, 0);
    MatrixDC_Inv(matrix, matrix);
    assert(FABS(creal(MatrixDC_GetElement(matrix, 0, 0))-0.5)<1E-6);
    assert(FABS(creal(MatrixDC_GetElement(matrix, 0, 1))-0.5)<1E-6);
    assert(FABS(creal(MatrixDC_GetElement(matrix, 0, 2))+0.5)<1E-6);
    assert(FABS(creal(MatrixDC_GetElement(matrix, 1, 0))-0.5)<1E-6);
    assert(FABS(creal(MatrixDC_GetElement(matrix, 1, 1))+0.5)<1E-6);
    assert(FABS(creal(MatrixDC_GetElement(matrix, 1, 2))-0.5)<1E-6);
    assert(FABS(creal(MatrixDC_GetElement(matrix, 2, 0))+0.5)<1E-6);
    assert(FABS(creal(MatrixDC_GetElement(matrix, 2, 1))-0.5)<1E-6);
    assert(FABS(creal(MatrixDC_GetElement(matrix, 2, 2))-0.5)<1E-6);
    MatrixDC_Del(matrix);
    printf("Test_MatrixDC_Inv Success.\n");
}

static void Test_MatrixDC_MulMatrixDC()
{
    MatrixDC* matrix1,* matrix2;

    matrix1 = MatrixDC_New(3,3);
    matrix2 = MatrixDC_New(3,3);
    MatrixDC_Set(matrix1, 2, 1);
    MatrixDC_Set(matrix2, 0, 1);
    MatrixDC_MulMatrixDC(matrix1, matrix2, matrix1);
    assert(FABS(creal(MatrixDC_GetElement(matrix1, 0, 0))-2)<1E-6);
    assert(FABS(creal(MatrixDC_GetElement(matrix1, 0, 1))-3)<1E-6);
    assert(FABS(creal(MatrixDC_GetElement(matrix1, 0, 2))-3)<1E-6);
    assert(FABS(creal(MatrixDC_GetElement(matrix1, 1, 0))-3)<1E-6);
    assert(FABS(creal(MatrixDC_GetElement(matrix1, 1, 1))-2)<1E-6);
    assert(FABS(creal(MatrixDC_GetElement(matrix1, 1, 2))-3)<1E-6);
    assert(FABS(creal(MatrixDC_GetElement(matrix1, 2, 0))-3)<1E-6);
    assert(FABS(creal(MatrixDC_GetElement(matrix1, 2, 1))-3)<1E-6);
    assert(FABS(creal(MatrixDC_GetElement(matrix1, 2, 2))-2)<1E-6);
    MatrixDC_Del(matrix1);
    MatrixDC_Del(matrix2);
    printf("Test_MatrixDC_MulMatrixDC Success.\n");
}

static void Test_MatrixDC_MulMatrixDC_NonSquare()
{
    MatrixDC* matrix1,* matrix2,* matrix3;

    matrix1 = MatrixDC_New(2,3);
    matrix2 = MatrixDC_New(3,2);
    matrix3 = MatrixDC_New(2,2);
    MatrixDC_Set(matrix1, 1, 2);
    MatrixDC_SetElement(matrix1, 0, 2, 3);
    MatrixDC_Show(matrix1);
    MatrixDC_Set(matrix2, 2, 3);
    MatrixDC_Show(matrix2);
    MatrixDC_SetElement(matrix2, 0, 1, 1);
    MatrixDC_MulMatrixDC(matrix1, matrix2, matrix3);
    MatrixDC_Show(matrix3);
    assert(FABS(creal(MatrixDC_GetElement(matrix3, 0, 0))-17)<1E-6);
    assert(FABS(creal(MatrixDC_GetElement(matrix3, 0, 1))-14)<1E-6);
    assert(FABS(creal(MatrixDC_GetElement(matrix3, 1, 0))-13)<1E-6);
    assert(FABS(creal(MatrixDC_GetElement(matrix3, 1, 1))-10)<1E-6);
    MatrixDC_Del(matrix1);
    MatrixDC_Del(matrix2);
    printf("Test_MatrixDC_MulMatrixDC_NonSquare Success.\n");    
}

static void Test_MatrixDC_T_Overlapped()
{
    MatrixDC* matrix;

    matrix = MatrixDC_New(3,3);
    MatrixDC_Set(matrix, 0, 1);
    MatrixDC_SetElement(matrix, 0, 2, 3);
    MatrixDC_T(matrix, matrix);
    assert(MatrixDC_GetElement(matrix, 2, 0) == 3);
    MatrixDC_Del(matrix);
    printf("Test MatrixDC_T_Overlapped Success.\n");
}

static void Test_MatrixDC_T_NonOverlapped()
{
    MatrixDC* matrix,* out_matrix;

    matrix = MatrixDC_New(3,3);
    out_matrix = MatrixDC_New(3,3);
    MatrixDC_Set(matrix, 0, 1);
    MatrixDC_SetElement(matrix, 0, 2, 3);
    MatrixDC_T(matrix, out_matrix);
    assert(MatrixDC_GetElement(out_matrix, 2, 0) == 3);
    MatrixDC_Del(out_matrix);
    MatrixDC_Del(matrix);
    printf("Test MatrixDC_T_NonOverlapped Success.\n");    
}

static void Test_MatrixDC_EigenEquation()
{
    MatrixDC* matrix;
    VectorDC* vector;
    
    matrix = MatrixDC_New(3,3);
    vector = VectorDC_New(3);
    MatrixDC_Set(matrix, 0, 1);
    MatrixDC_EigenEquation(matrix, vector, matrix);
    assert(FABS(creal(VectorDC_GetElement(vector, 0))+1)<1E-6);
    assert(FABS(creal(VectorDC_GetElement(vector, 1))-2)<1E-6);
    assert(FABS(creal(VectorDC_GetElement(vector, 2))+1)<1E-6);
    assert(FABS(creal(MatrixDC_GetElement(matrix, 0, 0)+MatrixDC_GetElement(matrix, 1, 0)+MatrixDC_GetElement(matrix, 2, 0)))<1E-6);
    assert(FABS(creal(MatrixDC_GetElement(matrix, 0, 1) - MatrixDC_GetElement(matrix, 1, 1)))<1E-6);
    assert(FABS(creal(MatrixDC_GetElement(matrix, 1, 1) - MatrixDC_GetElement(matrix, 2, 1)))<1E-6);
    assert(FABS(creal(MatrixDC_GetElement(matrix, 0, 2)+MatrixDC_GetElement(matrix, 1, 2)+MatrixDC_GetElement(matrix, 2, 2)))<1E-6);
    printf("Test_MatrixDC_EigenEquation Success.\n");
}

static void Test_MatrixDC_DiagLeftMul()
{
    MatrixDC* matrix;
    VectorDC* diag_matrix;

    matrix = MatrixDC_New(3,3);
    diag_matrix = VectorDC_New(3);
    MatrixDC_Set(matrix, 1, 2);
    MatrixDC_SetElement(matrix, 0, 2, 3);
    VectorDC_Full(diag_matrix, 3);
    MatrixDC_DiagLeftMul(matrix, diag_matrix, matrix);
    assert(FABS(creal(MatrixDC_GetElement(matrix, 0, 0))-3)<1E-6);
    assert(FABS(creal(MatrixDC_GetElement(matrix, 0, 1))-6)<1E-6);
    assert(FABS(creal(MatrixDC_GetElement(matrix, 0, 2))-9)<1E-6);
    assert(FABS(creal(MatrixDC_GetElement(matrix, 1, 0))-6)<1E-6);
    assert(FABS(creal(MatrixDC_GetElement(matrix, 1, 1))-3)<1E-6);
    assert(FABS(creal(MatrixDC_GetElement(matrix, 1, 2))-6)<1E-6);
    assert(FABS(creal(MatrixDC_GetElement(matrix, 2, 0))-6)<1E-6);
    assert(FABS(creal(MatrixDC_GetElement(matrix, 2, 1))-6)<1E-6);
    assert(FABS(creal(MatrixDC_GetElement(matrix, 2, 2))-3)<1E-6);
    MatrixDC_Del(matrix);
    VectorDC_Del(diag_matrix);
    printf("Test_MatrixDC_DiagLeftMul Success.\n");
}

static void Test_MatrixDC_DiagRightMul()
{
    MatrixDC* matrix;
    VectorDC* diag_matrix;

    matrix = MatrixDC_New(3,3);
    diag_matrix = VectorDC_New(3);
    MatrixDC_Set(matrix, 1, 2);
    MatrixDC_SetElement(matrix, 0, 2, 3);
    VectorDC_SetElement(diag_matrix, 0, 1);
    VectorDC_SetElement(diag_matrix, 1, 2);
    VectorDC_SetElement(diag_matrix, 2, 3);
    MatrixDC_DiagRightMul(matrix, diag_matrix, matrix);
    assert(FABS(creal(MatrixDC_GetElement(matrix, 0, 0))-1)<1E-6);
    assert(FABS(creal(MatrixDC_GetElement(matrix, 0, 1))-4)<1E-6);
    assert(FABS(creal(MatrixDC_GetElement(matrix, 0, 2))-9)<1E-6);
    assert(FABS(creal(MatrixDC_GetElement(matrix, 1, 0))-2)<1E-6);
    assert(FABS(creal(MatrixDC_GetElement(matrix, 1, 1))-2)<1E-6);
    assert(FABS(creal(MatrixDC_GetElement(matrix, 1, 2))-6)<1E-6);
    assert(FABS(creal(MatrixDC_GetElement(matrix, 2, 0))-2)<1E-6);
    assert(FABS(creal(MatrixDC_GetElement(matrix, 2, 1))-4)<1E-6);
    assert(FABS(creal(MatrixDC_GetElement(matrix, 2, 2))-3)<1E-6);
    MatrixDC_Del(matrix);
    VectorDC_Del(diag_matrix);
    printf("Test_MatrixDC_DiagRightMul Success.\n");    
}

static void Test_MatrixDC_LinearSolver()
{
    MatrixDC* matrix,* matrix2;
    VectorDC* vector;

    matrix = MatrixDC_New(2,2);
    matrix2 = MatrixDC_New(2,2);
    vector = VectorDC_New(2);
    MatrixDC_SetElement(matrix, 0, 0, -0.0783603-0.385762*I);
    MatrixDC_SetElement(matrix, 0, 1, -1.1737+0.238416*I);
    MatrixDC_SetElement(matrix, 1, 0, -0.408132+0.0503494*I);
    MatrixDC_SetElement(matrix, 1, 1, -0.148936-1.20727*I);
    VectorDC_Full(vector, 0);
    VectorDC_SetElement(vector, 0, 1);
    VectorDC_SetElement(vector, 1, 0);
    MatrixDC_LinearSolver(matrix, vector, vector);
    printf("Vector\n");
    VectorDC_Show(vector);
    MatrixDC_MulVectorDC(matrix, vector, vector);
    VectorDC_Show(vector);
    MatrixDC_Del(matrix);
    VectorDC_Del(vector);
    printf("Test_MatrixDC_LinearSolver Success.\n");
}

void Test_MatrixDC_PInv()
{
    MatrixDC* matrix1,* matrix2;

    matrix1 = MatrixDC_New(2,2);
    matrix2 = MatrixDC_New(2,2);
    MatrixDC_SetElement(matrix1, 0, 0, 1);
    MatrixDC_SetElement(matrix1, 0, 1, -1);
    MatrixDC_SetElement(matrix1, 1, 0, 1);
    MatrixDC_SetElement(matrix1, 1, 1, 1);
    MatrixDC_PInv(matrix1, matrix2);
    MatrixDC_Show(matrix2);
    printf("Test_MatrixDC_PInv Success.\n");
}

UnitTest_MatrixDC* UnitTest_MatrixDC_New()
{
    UnitTest_MatrixDC* unit_test;
    
    unit_test = (UnitTest_MatrixDC*)malloc(sizeof(UnitTest_MatrixDC));
    unit_test->class = UNITTEST_MATRIXDC;
    unit_test->iUnitTest = (IUnitTest*)malloc(sizeof(IUnitTest));
    unit_test->iUnitTest->implementor = unit_test;
    unit_test->iUnitTest->TestAll = _UnitTest_MatrixDC_TestAll;
    
    return unit_test;
}

UnitTest_MatrixDC* UnitTest_MatrixDC_Del(UnitTest_MatrixDC* this)
{
    free(this->iUnitTest);
    free(this);
    
    return NULL;
}

IUnitTest* UnitTest_MatrixDC_GetIUnitTest(UnitTest_MatrixDC* this)
{
    return this->iUnitTest;
}

void UnitTest_MatrixDC_TestAll(UnitTest_MatrixDC* this)
{
    Test_MatrixDC_Big();
    Test_MatrixDC_Set();
    Test_MatrixDC_Copy();
    Test_MatrixDC_IObject_New_Copy();
    Test_MatrixDC_Inv();
    Test_MatrixDC_MulMatrixDC();
    Test_MatrixDC_MulMatrixDC_NonSquare();
    Test_MatrixDC_EigenEquation();
    Test_MatrixDC_T_Overlapped();
    Test_MatrixDC_T_NonOverlapped();
    Test_MatrixDC_DiagLeftMul();
    Test_MatrixDC_DiagRightMul();
    Test_MatrixDC_LinearSolver();
    Test_MatrixDC_PInv();
}

void _UnitTest_MatrixDC_TestAll(void* implementor)
{
    UnitTest_MatrixDC* unit_test;
    
    unit_test = (UnitTest_MatrixDC*)implementor;
    UnitTest_MatrixDC_TestAll(unit_test);
}
