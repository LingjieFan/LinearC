#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "Class.h"
#include "IUnitTest.h"
#include "UnitTest_VectorD.h"
#include "UnitTest_VectorDC.h"
#include "UnitTest_MatrixD.h"
#include "UnitTest_MatrixDC.h"
#include "UnitTest_LinearC.h"

Class UnitTest_LinearC_Class = {CLASS};

struct _UnitTest_LinearC
{
    Class* class;
    IUnitTest* iUnitTest;
    UnitTest_VectorD* unitTest_vectorD;
    UnitTest_VectorDC* unitTest_vectorDC;
    UnitTest_MatrixD* unitTest_matrixD;
    UnitTest_MatrixDC* unitTest_matrixDC;
};

void _UnitTest_LinearC_TestAll(void* implementor);

UnitTest_LinearC* UnitTest_LinearC_New()
{
    UnitTest_LinearC* unit_test;
    
    unit_test = (UnitTest_LinearC*)malloc(sizeof(UnitTest_LinearC)+sizeof(IUnitTest));
    unit_test->class = UNITTEST_LINEARC;
    unit_test->iUnitTest = (IUnitTest*)(unit_test+1);
    unit_test->iUnitTest->implementor = unit_test;
    unit_test->iUnitTest->TestAll = _UnitTest_LinearC_TestAll;
    unit_test->unitTest_vectorD = UnitTest_VectorD_New();
    unit_test->unitTest_vectorDC = UnitTest_VectorDC_New();
    unit_test->unitTest_matrixD = UnitTest_MatrixD_New();
    unit_test->unitTest_matrixDC = UnitTest_MatrixDC_New();
    
    return unit_test;
}

UnitTest_LinearC* UnitTest_LinearC_Del(UnitTest_LinearC* this)
{
    UnitTest_VectorD_Del(this->unitTest_vectorD);
    UnitTest_VectorDC_Del(this->unitTest_vectorDC);
    UnitTest_MatrixD_Del(this->unitTest_matrixD);
    UnitTest_MatrixDC_Del(this->unitTest_matrixDC);
    free(this);
    
    return NULL;
}

IUnitTest* UnitTest_LinearC_GetIUnitTest(UnitTest_LinearC* this)
{
    return this->iUnitTest;
}

void UnitTest_LinearC_TestAll(UnitTest_LinearC* this)
{
    UnitTest_VectorD_TestAll(this->unitTest_vectorD);
    UnitTest_VectorDC_TestAll(this->unitTest_vectorDC);
    UnitTest_MatrixD_TestAll(this->unitTest_matrixD);
    UnitTest_MatrixDC_TestAll(this->unitTest_matrixDC);
}

void _UnitTest_LinearC_TestAll(void* implementor)
{
    UnitTest_LinearC* unit_test;
    
    unit_test = (UnitTest_LinearC*)implementor;
    UnitTest_LinearC_TestAll(unit_test);
}
