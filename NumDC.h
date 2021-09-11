#ifndef _NUMDC_H_
#define _NUMDC_H_

#include "Num.h"

#ifdef _MSC_VER
extern NumDC *NumDC_New(_Dcomplex num);
#else
extern NumDC *NumDC_New(double complex num);
#endif /*_MSC_VER*/

extern NumDC *NumDC_Del(NumDC *pNumDC);

extern NumDC *NumDC_AddNumD(NumDC *pNumDC1, NumD *pNumD2, NumDC *pONumDC);

extern NumD *NumDC_AddNumD2(NumDC *pNumDC1, NumD *pNumD2, NumD *pONumD);

extern NumDC *NumDC_AddNumDC(NumDC *pNumDC1, NumDC *pNumDC2, NumDC *pONumDC);

extern NumD *NumDC_AddNumDC2(NumDC *pNumDC1, NumDC *pNumDC2, NumD *pONumD);

extern NumDC *NumDC_SubNumD(NumDC *pNumDC1, NumD *pNumD2, NumDC *pONumDC);

extern NumD *NumDC_SubNumD2(NumDC *pNumDC1, NumD *pNumD2, NumD *pONumD);

extern NumDC *NumDC_SubNumDC(NumDC *pNumDC1, NumDC *pNumDC2, NumDC *pONumDC);

extern NumD *NumDC_SubNumDC2(NumDC *pNumDC1, NumDC *pNumDC2, NumD *pONumD);

extern NumDC *NumDC_MulNumD(NumDC *pNumDC1, NumD *pNumD2, NumDC *pONumDC);

extern NumD *NumDC_MulNumD2(NumDC *pNumDC1, NumD *pNumD2, NumD *pONumD);

extern NumDC *NumDC_MulNumDC(NumDC *pNumDC1, NumDC *pNumDC2, NumDC *pONumDC);

extern NumD *NumDC_MulNumDC2(NumDC *pNumDC1, NumDC *pNumDC2, NumD *pONumD);

extern NumDC *NumDC_DivNumD(NumDC *pNumDC1, NumD *pNumDC2, NumDC *pONumDC);

extern NumD *NumDC_DivNumD2(NumDC *pNumDC1, NumD *pNumDC2, NumD *pONumD);

extern NumDC *NumDC_DivNumDC(NumDC *pNumDC1, NumDC *pNumDC2, NumDC *pONumDC);

extern NumD *NumDC_DivNumDC2(NumDC *pNumDC1, NumDC *pNumDC2, NumD *pONumD);

#endif /*_NUMDC_H_*/
