#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "Class.h"
#include "IUnitTest.h"
#include "VectorD.h"
#include "UnitTest_VectorD.h"

Class UnitTest_VectorD_Class = {CLASS};

struct _UnitTest_VectorD
{
    Class* class;
    IUnitTest* iUnitTest;
};

void _UnitTest_VectorD_TestAll(void* implementor);

void Test_VectorD_MulVectorD()
{
    double result;
    VectorD* in_vector1,* in_vector2;
    
    in_vector1 = VectorD_New(3);
    in_vector2 = VectorD_New(3);
    in_vector1 = VectorD_Linspace(in_vector1, 1, 3);
    in_vector2 = VectorD_Full(in_vector2, 3);
    result = VectorD_MulVectorD(in_vector1, in_vector2);
    assert(result == 18.0);
    VectorD_Del(in_vector1);
    VectorD_Del(in_vector2);
    printf("Test_VectorD_MulVectorD Success.\n");
}

void Test_VectorD_GetSize()
{
    int size;
    VectorD* vector;
    
    vector = VectorD_New(12);
    size = VectorD_GetSize(vector);
    assert(size == 12);
    VectorD_Del(vector);
    printf("Test_VectorD_GetSize Success.\n");
}

void Test_VectorD_Copy()
{
    VectorD* vector1,* vector2;
    
    vector1 = VectorD_New(2);
    vector2 = VectorD_New(2);
    VectorD_Linspace(vector1, 0, 1);
    VectorD_Full(vector2, 11);
    VectorD_Copy(vector1, vector2);
    assert(VectorD_GetElement(vector2, 0) == 0);
    assert(VectorD_GetElement(vector2, 1) == 1);
    VectorD_Del(vector1);
    VectorD_Del(vector2);
    printf("Test_VectorD_Copy Success.\n");
}

void Test_VectorD_SubVectorD()
{
    VectorD* vector1,* vector2,* vector3;
    
    vector1 = VectorD_New(2);
    vector2 = VectorD_New(2);
    vector3 = VectorD_New(2);
    VectorD_Linspace(vector1, 0, 1);
    VectorD_Full(vector2, 11);
    VectorD_SubVectorD(vector1, vector2, vector3);
    assert(VectorD_GetElement(vector3, 0) == -11);
    assert(VectorD_GetElement(vector3, 1) == -10);
    VectorD_Del(vector1);
    VectorD_Del(vector2);
    VectorD_Del(vector3);
    printf("Test_VectorD_SubVectorD Success.\n");
}

void Test_VectorD_AddVectorD()
{
    VectorD* vector1,* vector2,* vector3;
    
    vector1 = VectorD_New(2);
    vector2 = VectorD_New(2);
    vector3 = VectorD_New(2);
    VectorD_Linspace(vector1, 0, 1);
    VectorD_Full(vector2, 11);
    VectorD_AddVectorD(vector1, vector2, vector3);
    assert(VectorD_GetElement(vector3, 0) == 11);
    assert(VectorD_GetElement(vector3, 1) == 12);
    VectorD_Del(vector1);
    VectorD_Del(vector2);
    VectorD_Del(vector3);
    printf("Test_VectorD_SubVectorD Success.\n");
}

void Test_VectorD_NormSq()
{
    double result;
    VectorD* vector1;
    
    vector1 = VectorD_New(2);
    VectorD_Linspace(vector1, 3,4);
    result = VectorD_NormSq(vector1);
    assert(result == 25);
    VectorD_Del(vector1);
    printf("Test_VectorD_NormSq Success.\n");
}

void Test_VectorD_GetInc()
{
    int inc;
    VectorD* vector;
    
    vector = VectorD_New(2);
    inc = VectorD_GetInc(vector);
    assert(inc == 1);
    VectorD_Del(vector);
    printf("Test_VectorD_GetInc Success.\n");
}

void Test_VectorD_Norm()
{
    double norm;
    VectorD* vector;
    
    vector = VectorD_New(2);
    VectorD_Linspace(vector, 3, 4);
    norm = VectorD_Norm(vector);
    assert(norm == 5);
    VectorD_Del(vector);
    printf("Test_VectorD_Norm Success.\n");
}

void Test_VectorD_MulNumD()
{
    VectorD* vector1,* vector2;
    
    vector1 = VectorD_New(2);
    vector2 = VectorD_New(2);
    VectorD_Linspace(vector1, 3, 4);
    VectorD_MulNumD(vector1, 3, vector2);
    assert(VectorD_GetElement(vector2, 0) == 9);
    assert(VectorD_GetElement(vector2, 1) == 12);
    VectorD_Del(vector1);
    VectorD_Del(vector2);
    printf("Test_VectorD_MulNumD Success.\n");
}

UnitTest_VectorD* UnitTest_VectorD_New()
{
    UnitTest_VectorD* unit_test;
    
    unit_test = (UnitTest_VectorD*)malloc(sizeof(UnitTest_VectorD)+sizeof(IUnitTest));
    unit_test->class = UNITTEST_VECTORD;
    unit_test->iUnitTest = (IUnitTest*)(unit_test+1);
    unit_test->iUnitTest->implementor = unit_test;
    unit_test->iUnitTest->TestAll = _UnitTest_VectorD_TestAll;
    
    return unit_test;
}

UnitTest_VectorD* UnitTest_VectorD_Del(UnitTest_VectorD* this)
{
    free(this);
    
    return NULL;
}

IUnitTest* UnitTest_VectorD_GetIUnitTest(UnitTest_VectorD* this)
{
    return this->iUnitTest;
}

void UnitTest_VectorD_TestAll(UnitTest_VectorD* this)
{
    Test_VectorD_MulVectorD();
    Test_VectorD_Copy();
    Test_VectorD_GetSize();
    Test_VectorD_NormSq();
    Test_VectorD_SubVectorD();
    Test_VectorD_AddVectorD();
    Test_VectorD_GetInc();
    Test_VectorD_MulNumD();
    Test_VectorD_Norm();
}

void _UnitTest_VectorD_TestAll(void* implementor)
{
    UnitTest_VectorD* unit_test;
    
    unit_test = (UnitTest_VectorD*)implementor;
    UnitTest_VectorD_TestAll(unit_test);
}
