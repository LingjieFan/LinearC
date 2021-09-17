#include "Num.h"

Num *Num_Del(Num *pNum)
{
    if(pNum != NULL)
    {
        return (*pNum->Del)(pNum);
    }
}

Num *_NumD_Del(Num *pNum)
{
    return (Num *)NumD_Del((NumD *)pNum);
}

Num *_NumDC_Del(Num *pNum)
{
    return (Num *)NumDC_Del((NumDC *)pNum);
}

Num *Num_AddNum(Num *pNum1, Num *pNum2, Num *pONum)
{
    return (*pNum1->AddNum)(pNum1, pNum2, pONum);
}

Num *_NumD_AddNum(Num *pNum1, Num *pNum2, Num *pONum)
{
    if(instanceof(pNum2,NUMD))
    {
        if(instanceof(pONum,NUMD))
        {
            NumD_AddNumD((NumD *)pNum1, (NumD *)pNum2, (NumD *)pONum);
        }
        else if(instanceof(pONum, NUMDC))
        {
            NumD_AddNumD2((NumD *)pNum1, (NumD *)pNum2, (NumDC *)pONum);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else if(instanceof(pNum2, NUMDC))
    {
        if(instanceof(pONum, NUMD))
        {
            NumD_AddNumDC2((NumD *)pNum1, (NumDC *)pNum2, (NumD *)pONum);
        }
        else if(instanceof(pONum, NUMDC))
        {
            NumD_AddNumDC((NumD *)pNum1, (NumDC *)pNum2,  (NumDC *)pONum);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pONum;
}

Num *_NumDC_AddNum(Num *pNum1, Num *pNum2, Num *pONum)
{
    if(instanceof(pNum2,NUMD))
    {
        if(instanceof(pONum, NUMD))
        {
            NumDC_AddNumD2((NumDC *)pNum1, (NumD *)pNum2, (NumD *)pONum);
        }
        else if(instanceof(pONum, NUMDC))
        {
            NumDC_AddNumD((NumDC *)pNum1, (NumD *)pNum2, (NumDC *)pONum);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else if(instanceof(pNum2, NUMDC))
    {
        if(instanceof(pONum, NUMD))
        {
            NumDC_AddNumDC2((NumDC *)pNum1, (NumDC *)pNum2, (NumD *)pONum);
        }
        else if(instanceof(pONum, NUMDC))
        {
            NumDC_AddNumDC((NumDC *)pNum1, (NumDC *)pNum2, (NumDC *)pONum);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pONum;
}

Num *Num_SubNum(Num *pNum1, Num *pNum2, Num *pONum)
{
    return (*pNum1->SubNum)(pNum1, pNum2, pONum);
}

Num *_NumD_SubNum(Num *pNum1, Num *pNum2, Num *pONum)
{
    if(instanceof(pNum2,NUMD))
    {
        if(instanceof(pONum, NUMD))
        {
            NumD_SubNumD((NumD *)pNum1, (NumD *)pNum2, (NumD *)pONum);
        }
        else if(instanceof(pONum, NUMDC))
        {
            NumD_SubNumD2((NumD *)pNum1, (NumD *)pNum2, (NumDC *)pONum);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else if(instanceof(pNum2, NUMDC))
    {
        if(instanceof(pONum, NUMD))
        {
            NumD_SubNumDC2((NumD *)pNum1, (NumDC *)pNum2, (NumD *)pONum);
        }
        else if(instanceof(pONum, NUMDC))
        {
            NumD_SubNumDC((NumD *)pNum1, (NumDC *)pNum2,  (NumDC *)pONum);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pONum;
}

Num *_NumDC_SubNum(Num *pNum1, Num *pNum2, Num *pONum)
{
    if(instanceof(pNum2,NUMD))
    {
        if(instanceof(pONum, NUMD))
        {
            NumDC_SubNumD2((NumDC *)pNum1, (NumD *)pNum2, (NumD *)pONum);
        }
        else if(instanceof(pONum, NUMDC))
        {
            NumDC_SubNumD((NumDC *)pNum1, (NumD *)pNum2, (NumDC *)pONum);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else if(instanceof(pNum2, NUMDC))
    {
        if(instanceof(pONum, NUMD))
        {
            NumDC_SubNumDC2((NumDC *)pNum1, (NumDC *)pNum2, (NumD *)pONum);
        }
        else if(instanceof(pONum, NUMDC))
        {
            NumDC_SubNumDC((NumDC *)pNum1, (NumDC *)pNum2, (NumDC *)pONum);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pONum;
}

Num *Num_MulNum(Num *pNum1, Num *pNum2, Num *pONum)
{
    return (*pNum1->MulNum)(pNum1, pNum2, pONum);
}

Num *_NumD_MulNum(Num *pNum1, Num *pNum2, Num *pONum)
{
    if(instanceof(pNum2,NUMD))
    {
        if(instanceof(pONum, NUMD))
        {
            NumD_MulNumD((NumD *)pNum1, (NumD *)pNum2, (NumD *)pONum);
        }
        else if(instanceof(pONum, NUMDC))
        {
            NumD_MulNumD2((NumD *)pNum1, (NumD *)pNum2, (NumDC *)pONum);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else if(instanceof(pNum2, NUMDC))
    {
        if(instanceof(pONum, NUMD))
        {
            NumD_MulNumDC2((NumD *)pNum1, (NumDC *)pNum2, (NumD *)pONum);
        }
        else if(instanceof(pONum, NUMDC))
        {
            NumD_MulNumDC((NumD *)pNum1, (NumDC *)pNum2,  (NumDC *)pONum);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pONum;
}

Num *_NumDC_MulNum(Num *pNum1, Num *pNum2, Num *pONum)
{
    if(instanceof(pNum2,NUMD))
    {
        if(instanceof(pONum, NUMD))
        {
            NumDC_MulNumD2((NumDC *)pNum1, (NumD *)pNum2, (NumD *)pONum);
        }
        else if(instanceof(pONum, NUMDC))
        {
            NumDC_MulNumD((NumDC *)pNum1, (NumD *)pNum2, (NumDC *)pONum);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else if(instanceof(pNum2, NUMDC))
    {
        if(instanceof(pONum, NUMD))
        {
            NumDC_MulNumDC2((NumDC *)pNum1, (NumDC *)pNum2, (NumD *)pONum);
        }
        else if(instanceof(pONum, NUMDC))
        {
            NumDC_MulNumDC((NumDC *)pNum1, (NumDC *)pNum2, (NumDC *)pONum);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pONum;
}

Num *Num_DivNum(Num *pNum1, Num *pNum2, Num *pONum)
{
    return (*pNum1->DivNum)(pNum1, pNum2, pONum);
}

Num *_NumD_DivNum(Num *pNum1, Num *pNum2, Num *pONum)
{
    if(instanceof(pNum2,NUMD))
    {
        if(instanceof(pONum, NUMD))
        {
            NumD_DivNumD((NumD *)pNum1, (NumD *)pNum2, (NumD *)pONum);
        }
        else if(instanceof(pONum, NUMDC))
        {
            NumD_DivNumD2((NumD *)pNum1, (NumD *)pNum2, (NumDC *)pONum);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else if(instanceof(pNum2, NUMDC))
    {
        if(instanceof(pONum, NUMD))
        {
            NumD_DivNumDC2((NumD *)pNum1, (NumDC *)pNum2, (NumD *)pONum);
        }
        else if(instanceof(pONum, NUMDC))
        {
            NumD_DivNumDC((NumD *)pNum1, (NumDC *)pNum2,  (NumDC *)pONum);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pONum;
}

Num *_NumDC_DivNum(Num *pNum1, Num *pNum2, Num *pONum)
{
    if(instanceof(pNum2,NUMD))
    {
        if(instanceof(pONum, NUMD))
        {
            NumDC_DivNumD2((NumDC *)pNum1, (NumD *)pNum2, (NumD *)pONum);
        }
        else if(instanceof(pONum, NUMDC))
        {
            NumDC_DivNumD((NumDC *)pNum1, (NumD *)pNum2, (NumDC *)pONum);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else if(instanceof(pNum2, NUMDC))
    {
        if(instanceof(pONum, NUMD))
        {
            NumDC_DivNumDC2((NumDC *)pNum1, (NumDC *)pNum2, (NumD *)pONum);
        }
        else if(instanceof(pONum, NUMDC))
        {
            NumDC_DivNumDC((NumDC *)pNum1, (NumDC *)pNum2, (NumDC *)pONum);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pONum;
}
