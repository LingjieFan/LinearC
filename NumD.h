#ifndef _NUMD_H_
#define _NUMD_H_

#include "Num.h"

extern NumD *NumD_New(double num);

extern NumD *NumD_Del(NumD *pNumD);

extern NumD *NumD_AddNumD(NumD *pNumD1, NumD *pNumD2, NumD *pONumD);

extern NumDC *NumD_AddNumD2(NumD *pNumD1, NumD *pNumD2, NumDC *pONumDC);

extern NumDC *NumD_AddNumDC(NumD *pNumD1, NumDC *pNumDC2, NumDC *pONumDC);

extern NumD *NumD_AddNumDC2(NumD *pNumD1, NumDC *pNumDC2, NumD *pONumD);

extern NumD *NumD_SubNumD(NumD *pNumD1, NumD *pNumD2, NumD *pONumD);

extern NumDC *NumD_SubNumD2(NumD *pNumD1, NumD *pNumD2, NumDC *pONumDC);

extern NumDC *NumD_SubNumDC(NumD *pNumD1, NumDC *pNumDC2, NumDC *pONumDC);

extern NumD *NumD_SubNumDC2(NumD *pNumD1, NumDC *pNumDC2, NumD *pONumD);

extern NumD *NumD_MulNumD(NumD *pNumD1, NumD *pNumD2, NumD *pONumD);

extern NumDC *NumD_MulNumD2(NumD *pNumD1, NumD *pNumD2, NumDC *pONumDC);

extern NumDC *NumD_MulNumDC(NumD *pNumD1, NumDC *pNumDC2, NumDC *pONumDC);

extern NumD *NumD_MulNumDC2(NumD *pNumD1, NumDC *pNumDC2, NumD *pONumD);

extern NumD *NumD_DivNumD(NumD *pNumD1, NumD *pNumD2, NumD *pONumD);

extern NumDC *NumD_DivNumD2(NumD *pNumD1, NumD *pNumD2, NumDC *pONumDC);

extern NumDC *NumD_DivNumDC(NumD *pNumD1, NumDC *pNumDC2, NumDC *pONumDC);

extern NumD *NumD_DivNumDC2(NumD *pNumD1, NumDC *pNumDC2, NumD *pONumD);

#endif /*_NUMD_H_*/