#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "Class.h"
#include "IUnitTest.h"
#include "IObject.h"
#include "VectorDC.h"
#include "UnitTest_VectorDC.h"

#define FABS(x) (x)>0?(x):-(x)

Class UnitTest_VectorDC_Class = {CLASS};

struct _UnitTest_VectorDC
{
    Class* class;
    IUnitTest* iUnitTest;
};

struct _VectorDC
{
    Class* class;
    IObject* iObject;
    double complex* vector;
    int size;
    int inc;
};

void _UnitTest_VectorDC_TestAll(void* implementor);

void Test_VectorDC_New()
{
    VectorDC* vector;

    vector  = VectorDC_New(3);
    assert(vector->class == VECTORDC);
    assert(vector->size == 3);
    assert(vector->inc == 1);
    VectorDC_Del(vector);
    printf("Test_VectorDC_New Success.\n");
}

void Test_VectorDC_IObject_New()
{
    VectorDC* vector1,* vector2;
    IObject* iObject;

    vector1 = VectorDC_New(3);
    iObject = VectorDC_GetIObject(vector1);
    vector2 = (VectorDC*)IObject_New(iObject);
    VectorDC_Del(vector1);
    VectorDC_Del(vector2);
    printf("Test_VectorDC_IObject_New Success.\n");
}

void Test_VectorDC_CMulVectorDC()
{
    VectorDC* vector1,* vector2;
    double complex result;

    vector1 = VectorDC_New(2);
    vector2 = VectorDC_New(2);

    VectorDC_Full(vector1, 1+1*I);
    VectorDC_Full(vector2, 1+1*I);
    result = VectorDC_CMulVectorDC(vector1, vector2);
    assert(FABS(creal(result)-4)<1E-6);
    printf("Test_VectorDC_CMulVectorDC Success.\n");
}

UnitTest_VectorDC* UnitTest_VectorDC_New()
{
    UnitTest_VectorDC* unit_test;
    
    unit_test = (UnitTest_VectorDC*)malloc(sizeof(UnitTest_VectorDC)+sizeof(IUnitTest));
    unit_test->class = UNITTEST_VECTORDC;
    unit_test->iUnitTest = (IUnitTest*)(unit_test+1);
    unit_test->iUnitTest->implementor = unit_test;
    unit_test->iUnitTest->TestAll = _UnitTest_VectorDC_TestAll;
    
    return unit_test;
}

UnitTest_VectorDC* UnitTest_VectorDC_Del(UnitTest_VectorDC* this)
{
    free(this);
    
    return NULL;
}

IUnitTest* UnitTest_VectorDC_GetIUnitTest(UnitTest_VectorDC* this)
{
    return this->iUnitTest;
}

void UnitTest_VectorDC_TestAll(UnitTest_VectorDC* this)
{
    Test_VectorDC_New();
    Test_VectorDC_IObject_New();
    Test_VectorDC_CMulVectorDC();
}

void _UnitTest_VectorDC_TestAll(void* implementor)
{
    UnitTest_VectorDC* unit_test;
    
    unit_test = (UnitTest_VectorDC*)implementor;
    UnitTest_VectorDC_TestAll(unit_test);
}

