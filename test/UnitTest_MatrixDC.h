#ifndef _UNITTEST_MATRIXDC_H_
#define _UNITTEST_MATRIXDC_H_

#define UNITTEST_MATRIXDC &UnitTest_MatrixDC_Class

typedef struct _Class Class;
typedef struct _IUnitTest IUnitTest;
typedef struct _UnitTest_MatrixDC UnitTest_MatrixDC;

extern Class UnitTest_MatrixDC_Class;

extern UnitTest_MatrixDC* UnitTest_MatrixDC_New();

extern UnitTest_MatrixDC* UnitTest_MatrixDC_Del(UnitTest_MatrixDC* this);

extern IUnitTest* UnitTest_MatrixDC_GetIUnitTest(UnitTest_MatrixDC* this);

extern void UnitTest_MatrixDC_TestAll(UnitTest_MatrixDC* this);

#endif /*_UNITTEST_MATRIXDC_H_*/