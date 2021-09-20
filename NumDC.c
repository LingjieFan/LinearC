#include "Num.h"

extern Object *_NumDC_Del(Object *pObject);
extern void _NumDC_Show(Tensor *pTensor);
extern Num *_NumDC_AddNum(Num *pNum1, Num *pNum2, Num *pONum);
extern Num *_NumDC_SubNum(Num *pNum1, Num *pNum2, Num *pONum);
extern Num *_NumDC_MulNum(Num *pNum1, Num *pNum2, Num *pONum);
extern Num *_NumDC_DivNum(Num *pNum1, Num *pNum2, Num *pONum);

#ifdef _MSC_VER
NumDC *NumDC_New(_Dcomplex num)
#else
NumDC *NumDC_New(double complex num)
#endif /*_MSC_VER*/
{
    NumDC *pNewNumDC;

    pNewNumDC = (NumDC *)malloc(sizeof(NumDC));
    if(pNewNumDC == NULL)
    {
        fprintf(stderr, "Error: Can't create a new numDC. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
        return NULL;
    }
    pNewNumDC->num = num;
    ((Object *)pNewNumDC)->type = NUMDC;
    ((Object *)pNewNumDC)->Del = _NumDC_Del;
    ((Object *)pNewNumDC)->UnWrap = _NumDC_Del;
    ((Tensor *)pNewNumDC)->Show = _NumDC_Show;
    pNewNumDC->parent.AddNum = _NumDC_AddNum;
    pNewNumDC->parent.SubNum = _NumDC_SubNum;
    pNewNumDC->parent.MulNum = _NumDC_MulNum;
    pNewNumDC->parent.DivNum = _NumDC_DivNum;

    return pNewNumDC;
}

NumDC *NumDC_Del(NumDC *pNumDC)
{
    if(pNumDC != NULL)
    {
        free(pNumDC);
    }

    return NULL;
}

void NumDC_Show(NumDC *pNumDC)
{
    printf("NumDC_Show: \n num:%lf+%lf*j\n\n", creal(pNumDC->num),cimag(pNumDC->num));
}

NumDC *NumDC_AddNumD(NumDC *pNumDC1, NumD *pNumD2, NumDC *pONumDC)
{
    #ifdef _MSC_VER
    pONumDC->num._Val[0] = pNumDC1->num._Val[0] + pNumD2->num;
    pONumDC->num._Val[1] = pNumDC1->num._Val[1];
    #else
    pONumDC->num = pNumDC1->num + pNumD2->num;
    #endif /*_MSC_VER*/

    return pONumDC;
}

NumD *NumDC_AddNumD2(NumDC *pNumDC1, NumD *pNumD2, NumD *pONumD)
{
    #ifdef _MSC_VER
    pONumD->num = pNumDC1->num._Val[0] + pNumD2->num;
    #else
    pONumD->num = creal(pNumDC1->num) + pNumD2->num;
    #endif /*_MSC_VER*/

    return pONumD;
}

NumDC *NumDC_AddNumDC(NumDC *pNumDC1, NumDC *pNumDC2, NumDC *pONumDC)
{
    #ifdef _MSC_VER
    pONumDC->num = _Caddcc(pNumDC1->num, pNumDC2->num);
    #else
    pONumDC->num = pNumDC1->num + pNumDC2->num;
    #endif /*_MSC_VER*/

    return pONumDC;
}

NumD *NumDC_AddNumDC2(NumDC *pNumDC1, NumDC *pNumDC2, NumD *pONumD)
{
    #ifdef _MSC_VER
    pONumD->num = pNumDC1->num._Val[0] + pNumDC2->num._Val[0];
    #else
    pONumD->num = creal(pNumDC1->num) + creal(pNumDC2->num);
    #endif /*_MSC_VER*/

    return pONumD;
}

NumDC *NumDC_SubNumD(NumDC *pNumDC1, NumD *pNumD2, NumDC *pONumDC)
{
    #ifdef _MSC_VER
    pONumDC->num._Val[0] = pNumDC1->num._Val[0] - pNumD2->num;
    pONumDC->num._Val[1] = pNumDC1->num._Val[1];
    #else
    pONumDC->num = pNumDC1->num - pNumD2->num;
    #endif /*_MSC_VER*/

    return pONumDC;
}

NumD *NumDC_SubNumD2(NumDC *pNumDC1, NumD *pNumD2, NumD *pONumD)
{
    #ifdef _MSC_VER
    pONumD->num = pNumDC1->num._Val[0] - pNumD2->num;
    #else
    pONumD->num = creal(pNumDC1->num) - pNumD2->num;
    #endif /*_MSC_VER*/

    return pONumD;
}

NumDC *NumDC_SubNumDC(NumDC *pNumDC1, NumDC *pNumDC2, NumDC *pONumDC)
{
    #ifdef _MSC_VER
    pONumDC->num = _Csubcc(pNumDC1->num, pNumDC2->num);
    #else
    pONumDC->num = pNumDC1->num - pNumDC2->num;
    #endif /*_MSC_VER*/

    return pONumDC;
}

NumD *NumDC_SubNumDC2(NumDC *pNumDC1, NumDC *pNumDC2, NumD *pONumD)
{
    #ifdef _MSC_VER
    pONumD->num = pNumDC1->num._Val[0] + pNumDC2->num._Val[0];
    #else
    pONumD->num = creal(pNumDC1->num) + creal(pNumDC2->num);
    #endif /*_MSC_VER*/

    return pONumD;
}

NumDC *NumDC_MulNumD(NumDC *pNumDC1, NumD *pNumD2, NumDC *pONumDC)
{
    #ifdef _MSC_VER
    pONumDC->num = _Cmulcr(pNumDC1->num,pNumD2->num);
    #else
    pONumDC->num = pNumDC1->num * pNumD2->num;
    #endif /*_MSC_VER*/

    return pONumDC;
}

NumD *NumDC_MulNumD2(NumDC *pNumDC1, NumD *pNumD2, NumD *pONumD)
{
    #ifdef _MSC_VER
    pONumD->num = pNumDC1->num._Val[0] + pNumD2->num;
    #else
    pONumD->num = creal(pNumDC1->num) + pNumD2->num;
    #endif /*_MSC_VER*/

    return pONumD;
}

NumDC *NumDC_MulNumDC(NumDC *pNumDC1, NumDC *pNumDC2, NumDC *pONumDC)
{
    #ifdef _MSC_VER
    pONumDC->num = _Cmulcc(pNumDC1->num, pNumDC2->num);
    #else
    pONumDC->num = pNumDC1->num * pNumDC2->num;
    #endif /*_MSC_VER*/

    return pONumDC;
}

NumD *NumDC_MulNumDC2(NumDC *pNumDC1, NumDC *pNumDC2, NumD *pONumD)
{
    #ifdef _MSC_VER
    pONumD->num = pNumDC1->num._Val[0] + pNumDC2->num._Val[0];
    #else
    pONumD->num = creal(pNumDC1->num) + creal(pNumDC2->num);
    #endif /*_MSC_VER*/

    return pONumD;
}

NumDC *NumDC_DivNumD(NumDC *pNumDC1, NumD *pNumD2, NumDC *pONumDC)
{
    #ifdef _MSC_VER
    pONumDC->num = _Cmulcr(pNumDC1->num, 1/pNumD2->num);
    #else
    pONumDC->num = pNumDC1->num / pNumD2->num;
    #endif /*_MSC_VER*/

    return pONumDC;
}

NumD *NumDC_DivNumD2(NumDC *pNumDC1, NumD *pNumD2, NumD *pONumD)
{
    #ifdef _MSC_VER
    pONumD->num = pNumDC1->num._Val[0] / pNumD2->num;
    #else
    pONumD->num = creal(pNumDC1->num / pNumD2->num);
    #endif /*_MSC_VER*/

    return pONumD;
}

NumDC *NumDC_DivNumDC(NumDC *pNumDC1, NumDC *pNumDC2, NumDC *pONumDC)
{
    #ifdef _MSC_VER
    pONumDC->num = _Cdivcc(pNumDC1->num, pNumDC2->num);
    #else
    pONumDC->num = pNumDC1->num / pNumDC2->num;
    #endif /*_MSC_VER*/

    return pONumDC;
}

NumD *NumDC_DivNumDC2(NumDC *pNumDC1, NumDC *pNumDC2, NumD *pONumD)
{
    #ifdef _MSC_VER
    pONumD->num = (pNumDC1->num._Val[0] * pNumDC2->num._Val[0] + pNumDC1->num._Val[1] * pNumDC2->num._Val[1]) / \
        (pNumDC1->num._Val[0] * pNumDC1->num._Val[0] + pNumDC1->num._Val[1] * pNumDC1->num._Val[1]);
    #else
    pONumD->num = creal(pNumDC1->num / pNumDC2->num);
    #endif /*_MSC_VER*/

    return pONumD;
}
