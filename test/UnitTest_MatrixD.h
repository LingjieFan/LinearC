#ifndef _UNITTEST_MATRIXD_H_
#define _UNITTEST_MATRIXD_H_

#define UNITTEST_MATRIXD &UnitTest_MatrixD_Class

typedef struct _Class Class;
typedef struct _IUnitTest IUnitTest;
typedef struct _UnitTest_MatrixD UnitTest_MatrixD;

extern Class UnitTest_MatrixD_Class;

extern UnitTest_MatrixD* UnitTest_MatrixD_New();

extern UnitTest_MatrixD* UnitTest_MatrixD_Del(UnitTest_MatrixD* this);

extern IUnitTest* UnitTest_MatrixD_GetIUnitTest(UnitTest_MatrixD* this);

extern void UnitTest_MatrixD_TestAll(UnitTest_MatrixD* this);

#endif /*_UNITTEST_MATRIXD_H_*/