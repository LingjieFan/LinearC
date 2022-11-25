#ifndef _UNITTEST_VECTORDC_H_
#define _UNITTEST_VECTORDC_H_

#define UNITTEST_VECTORDC &UnitTest_VectorDC_Class

typedef struct _Class Class;
typedef struct _IUnitTest IUnitTest;
typedef struct _UnitTest_VectorDC UnitTest_VectorDC;

extern Class UnitTest_VectorDC_Class;

extern UnitTest_VectorDC* UnitTest_VectorDC_New();

extern UnitTest_VectorDC* UnitTest_VectorDC_Del(UnitTest_VectorDC* this);

extern IUnitTest* UnitTest_VectorDC_GetIUnitTest(UnitTest_VectorDC* this);

extern void UnitTest_VectorDC_TestAll(UnitTest_VectorDC* this);

#endif /*_UNITTEST_VECTORDC_H_*/
