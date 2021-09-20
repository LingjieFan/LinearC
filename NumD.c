#include "Num.h"

extern Object *_NumD_Del(Object *pNum);
extern void _NumD_Show(Tensor *pTensor);
extern Num *_NumD_AddNum(Num *pNum1, Num *pNum2, Num *pONum);
extern Num *_NumD_SubNum(Num *pNum1, Num *pNum2, Num *pONum);
extern Num *_NumD_MulNum(Num *pNum1, Num *pNum2, Num *pONum);
extern Num *_NumD_DivNum(Num *pNum1, Num *pNum2, Num *pONum);

NumD *NumD_New(double num)
{
    NumD *pNewNumD;

    pNewNumD = (NumD *)malloc(sizeof(NumD));
    if(pNewNumD == NULL)
    {
        fprintf(stderr, "Error: Can't create a new numD. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
        return NULL;
    }
    pNewNumD->num = num;
    ((Object *)pNewNumD)->type = NUMD;
    ((Object *)pNewNumD)->Del = _NumD_Del;
    ((Object *)pNewNumD)->UnWrap = _NumD_Del;
    ((Tensor *)pNewNumD)->Show = _NumD_Show;
    pNewNumD->parent.AddNum = _NumD_AddNum;
    pNewNumD->parent.SubNum = _NumD_SubNum;
    pNewNumD->parent.MulNum = _NumD_MulNum;
    pNewNumD->parent.DivNum = _NumD_DivNum;

    return pNewNumD;
}

NumD *NumD_Del(NumD *pNumD)
{
    if(pNumD != NULL)
    {
        free(pNumD);
    }

    return NULL;
}

void NumD_Show(NumD *pNumD)
{
    printf("NumD_Show: \n num:%lf\n\n", pNumD->num);
}

NumD *NumD_AddNumD(NumD *pNumD1, NumD *pNumD2, NumD *pONumD)
{
    pONumD->num = pNumD1->num + pNumD2->num;
    return pONumD;
}

NumDC *NumD_AddNumD2(NumD *pNumD1, NumD *pNumD2, NumDC *pONumDC)
{
    #ifdef _MSC_VER
    pONumDC->num._Val[0] = pNumD1->num + pNumD1->num;
    pONumDC->num._Val[1] = 0;
    #else
    pONumDC->num = pNumD1->num + pNumD2->num;
    #endif /*_MSC_VER*/

    return pONumDC;
}

NumDC *NumD_AddNumDC(NumD *pNumD1, NumDC *pNumDC2, NumDC *pONumDC)
{
    #ifdef _MSC_VER
    pONumDC->num._Val[0] = pNumD1->num + pNumDC2->num._Val[0];
    pONumDC->num._Val[1] = pNumDC2->num._Val[1];
    #else
    pONumDC->num = pNumD1->num + pNumDC2->num;
    #endif /*_MSC_VER*/

    return pONumDC;
}

NumD *NumD_AddNumDC2(NumD *pNumD1, NumDC *pNumDC2, NumD *pONumD)
{
    #ifdef _MSC_VER
    pONumD->num = pNumD1->num + pNumDC2->num._Val[0];
    #else
    pONumD->num = pNumD1->num + creal(pNumDC2->num);
    #endif /*_MSC_VER*/

    return pONumD;
}

NumD *NumD_SubNumD(NumD *pNumD1, NumD *pNumD2, NumD *pONumD)
{
    pONumD->num = pNumD1->num - pNumD2->num;
    return pONumD;
}

NumDC *NumD_SubNumD2(NumD *pNumD1, NumD *pNumD2, NumDC *pONumDC)
{
    #ifdef _MSC_VER
    pONumDC->num._Val[0] = pNumD1->num - pNumD1->num;
    pONumDC->num._Val[1] = 0;
    #else
    pONumDC->num = pNumD1->num - pNumD2->num;
    #endif /*_MSC_VER*/

    return pONumDC;
}

NumDC *NumD_SubNumDC(NumD *pNumD1, NumDC *pNumDC2, NumDC *pONumDC)
{
    #ifdef _MSC_VER
    pONumDC->num._Val[0] = pNumD1->num - pNumDC2->num._Val[0];
    pONumDC->num._Val[1] = -pNumDC2->num._Val[1];
    #else
    pONumDC->num = pNumD1->num - pNumDC2->num;
    #endif /*_MSC_VER*/

    return pONumDC;
}

NumD *NumD_SubNumDC2(NumD *pNumD1, NumDC *pNumDC2, NumD *pONumD)
{
    #ifdef _MSC_VER
    pONumD->num = pNumD1->num - pNumDC2->num._Val[0];
    #else
    pONumD->num = pNumD1->num - creal(pNumDC2->num);
    #endif /*_MSC_VER*/

    return pONumD;
}

NumD *NumD_MulNumD(NumD *pNumD1, NumD *pNumD2, NumD *pONumD)
{
    pONumD->num = pNumD1->num * pNumD2->num;
    return pONumD;
}

NumDC *NumD_MulNumD2(NumD *pNumD1, NumD *pNumD2, NumDC *pONumDC)
{
    #ifdef _MSC_VER
    pONumDC->num._Val[0] = pNumD1->num * pNumD1->num;
    pONumDC->num._Val[1] = 0;
    #else
    pONumDC->num = pNumD1->num * pNumD2->num;
    #endif /*_MSC_VER*/

    return pONumDC;
}

NumDC *NumD_MulNumDC(NumD *pNumD1, NumDC *pNumDC2, NumDC *pONumDC)
{
    #ifdef _MSC_VER
    pONumDC->num._Val[0] = pNumD1->num * pNumDC2->num._Val[0];
    pONumDC->num._Val[1] = pNumD1->num * pNumDC2->num._Val[1];
    #else
    pONumDC->num = pNumD1->num * pNumDC2->num;
    #endif /*_MSC_VER*/

    return pONumDC;
}

NumD *NumD_MulNumDC2(NumD *pNumD1, NumDC *pNumDC2, NumD *pONumD)
{
    #ifdef _MSC_VER
    pONumD->num = pNumD1->num * pNumDC2->num._Val[0];
    #else
    pONumD->num = pNumD1->num * creal(pNumDC2->num);
    #endif /*_MSC_VER*/

    return pONumD;
}

NumD *NumD_DivNumD(NumD *pNumD1, NumD *pNumD2, NumD *pONumD)
{
    pONumD->num = pNumD1->num / pNumD2->num;
    return pONumD;
}

NumDC *NumD_DivNumD2(NumD *pNumD1, NumD *pNumD2, NumDC *pONumDC)
{
    #ifdef _MSC_VER
    pONumDC->num._Val[0] = pNumD1->num / pNumD1->num;
    pONumDC->num._Val[1] = 0;
    #else
    pONumDC->num = pNumD1->num / pNumD2->num;
    #endif /*_MSC_VER*/

    return pONumDC;
}

NumDC *NumD_DivNumDC(NumD *pNumD1, NumDC *pNumDC2, NumDC *pONumDC)
{
    #ifdef _MSC_VER
    pONumDC->num = _Cdivrc(pNumD1->num, pNumDC2->num);
    #else
    pONumDC->num = pNumD1->num / pNumDC2->num;
    #endif /*_MSC_VER*/

    return pONumDC;
}

NumD *NumD_DivNumDC2(NumD *pNumD1, NumDC *pNumDC2, NumD *pONumD)
{
    #ifdef _MSC_VER
    pONumD->num = pNumD1->num * pNumDC2->num._Val[0] / (pNumDC2->num._Val[0] * pNumDC2->num._Val[0] + \
        pNumDC2->num._Val[1] * pNumDC2->num._Val[1]);
    #else
    pONumD->num = creal(pNumD1->num / pNumDC2->num);
    #endif /*_MSC_VER*/

    return pONumD;
}
