#ifndef _UNITTEST_VECTORD_H_
#define _UNITTEST_VECTORD_H_

#define UNITTEST_VECTORD &UnitTest_VectorD_Class

typedef struct _Class Class;
typedef struct _IUnitTest IUnitTest;
typedef struct _UnitTest_VectorD UnitTest_VectorD;

extern Class UnitTest_VectorD_Class;

extern UnitTest_VectorD* UnitTest_VectorD_New();

extern UnitTest_VectorD* UnitTest_VectorD_Del(UnitTest_VectorD* this);

extern IUnitTest* UnitTest_VectorD_GetIUnitTest(UnitTest_VectorD* this);

extern void UnitTest_VectorD_TestAll(UnitTest_VectorD* this);

#endif /*_UNITTEST_VECTORD_H_*/