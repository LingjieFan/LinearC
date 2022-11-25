#ifndef _UNITTEST_LINEARC_H_
#define _UNITTEST_LINEARC_H_

#define UNITTEST_LINEARC &UnitTest_LinearC_Class

typedef struct _Class Class;
typedef struct _IUnitTest IUnitTest;
typedef struct _UnitTest_LinearC UnitTest_LinearC;

extern Class UnitTest_LinearC_Class;

extern UnitTest_LinearC* UnitTest_LinearC_New();

extern UnitTest_LinearC* UnitTest_LinearC_Del(UnitTest_LinearC* this);

extern IUnitTest* UnitTest_LinearC_GetIUnitTest(UnitTest_LinearC* this);

extern void UnitTest_LinearC_TestAll(UnitTest_LinearC* this);

#endif /*_UNITTEST_LINEARC_H_*/