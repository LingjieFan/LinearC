#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "Class.h"
#include "IUnitTest.h"
#include "VectorD.h"
#include "MatrixD.h"
#include "UnitTest_MatrixD.h"

#define FABS(x) (x)>0?(x):-(x)

Class UnitTest_MatrixD_Class = {CLASS};

struct _UnitTest_MatrixD
{
    Class* class;
    IUnitTest* iUnitTest;
};

void _UnitTest_MatrixD_TestAll(void* implementor);

void Test_MatirxD_LinearSolver()
{
    MatrixD* matrix;
    VectorD* in_vector,* out_vector;
    
    matrix = MatrixD_New(2,2);
    in_vector = VectorD_New(2);
    out_vector = VectorD_New(2);
    MatrixD_Set(matrix, 1, 1);
    MatrixD_SetElement(matrix, 1, 1, -1);
    VectorD_SetElement(in_vector, 0, 2);
    VectorD_SetElement(in_vector, 1, 4);
    MatrixD_LinearSolver(matrix, in_vector, out_vector);
    assert(VectorD_GetElement(out_vector, 0) == 3.0);
    assert(VectorD_GetElement(out_vector, 1) == -1.0);
    MatrixD_Del(matrix);
    VectorD_Del(in_vector);
    VectorD_Del(out_vector);
    printf("Test_MatirxD_LinearSolver Success.\n");
}

void Test_MatrixD_TMulSelf()
{
    MatrixD* matrix1,* matrix2;
    
    matrix1 = MatrixD_New(2,2);
    matrix2 = MatrixD_New(2,2);
    MatrixD_SetElement(matrix1, 0, 0, 1);
    MatrixD_SetElement(matrix1, 0, 1, 2);
    MatrixD_SetElement(matrix1, 1, 0, 3);
    MatrixD_SetElement(matrix1, 1, 1, 4);
    MatrixD_TMulSelf(matrix1, matrix2);
    assert(MatrixD_GetElement(matrix2, 0, 0) == 10);
    assert(MatrixD_GetElement(matrix2, 0, 1) == 14);
    assert(MatrixD_GetElement(matrix2, 1, 0) == 14);
    assert(MatrixD_GetElement(matrix2, 1, 1) == 20);
    MatrixD_Del(matrix1);
    MatrixD_Del(matrix2);
    printf("Test_MatrixD_TMulSelf Success.\n");
}

void Test_MatrixD_TMulVectorD()
{
    MatrixD* matrix;
    VectorD* vector1,* vector2;
    
    matrix = MatrixD_New(2,2);
    vector1 = VectorD_New(2);
    vector2 = VectorD_New(2);
    MatrixD_SetElement(matrix, 0, 0, 1);
    MatrixD_SetElement(matrix, 0, 1, 2);
    MatrixD_SetElement(matrix, 1, 0, 3);
    MatrixD_SetElement(matrix, 1, 1, 4);    
    VectorD_SetElement(vector1, 0, 3);
    VectorD_SetElement(vector1, 1, 4);
    MatrixD_TMulVectorD(matrix, vector1, vector2);
    assert(FABS(VectorD_GetElement(vector2, 0)-15)<1E-6);
    assert(FABS(VectorD_GetElement(vector2, 1)-22)<1E-6);
    MatrixD_Del(matrix);
    VectorD_Del(vector1);
    VectorD_Del(vector2);
    printf("Test_MatrixD_TMulVectorD Success.\n");
}

void Test_MatrixD_MaxDiagElem()
{
    double result;
    MatrixD* matrix;
    
    matrix = MatrixD_New(2,2);
    MatrixD_SetElement(matrix, 0, 0, 1);
    MatrixD_SetElement(matrix, 0, 1, 2);
    MatrixD_SetElement(matrix, 1, 0, 3);
    MatrixD_SetElement(matrix, 1, 1, 4);  
    result = MatrixD_MaxDiagElem(matrix);
    assert(result==4);
    MatrixD_Del(matrix);
    printf("Test_MatrixD_MaxDiagElem Success.\n");
}

UnitTest_MatrixD* UnitTest_MatrixD_New()
{
    UnitTest_MatrixD* unit_test;
    
    unit_test = (UnitTest_MatrixD*)malloc(sizeof(UnitTest_MatrixD)+sizeof(IUnitTest));
    unit_test->class = UNITTEST_MATRIXD;
    unit_test->iUnitTest = (IUnitTest*)(unit_test+1);
    unit_test->iUnitTest->implementor = unit_test;
    unit_test->iUnitTest->TestAll = _UnitTest_MatrixD_TestAll;
    
    return unit_test;
}

UnitTest_MatrixD* UnitTest_MatrixD_Del(UnitTest_MatrixD* this)
{
    free(this);
    
    return NULL;
}

IUnitTest* UnitTest_MatrixD_GetIUnitTest(UnitTest_MatrixD* this)
{
    return this->iUnitTest;
}

void UnitTest_MatrixD_TestAll(UnitTest_MatrixD* this)
{
    Test_MatirxD_LinearSolver();
    Test_MatrixD_TMulSelf();
    Test_MatrixD_TMulVectorD();
    Test_MatrixD_MaxDiagElem();
}

void _UnitTest_MatrixD_TestAll(void* implementor)
{
    UnitTest_MatrixD* unit_test;
    
    unit_test = (UnitTest_MatrixD*)implementor;
    UnitTest_MatrixD_TestAll(unit_test);
}
