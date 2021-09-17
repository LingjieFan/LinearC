#include "Vector.h"

Vector *Vector_Del(Vector *pIVector)
{
    if(pIVector != NULL)
    {
        return (*pIVector->Del)(pIVector);
    }

    return NULL;
}

Vector *_VectorD_Del(Vector *pIVector)
{
    return (Vector *)VectorD_Del((VectorD *)pIVector);
}

Vector *_VectorDC_Del(Vector *pIVector)
{
    return (Vector *)VectorDC_Del((VectorDC *)pIVector);
}

Vector *Vector_UnWrap(Vector *pIVector)
{
    if(pIVector != NULL)
    {
        return (*pIVector->UnWrap)(pIVector);
    }

    return NULL;
}

Vector *_VectorD_UnWrap(Vector *pIVector)
{
    return (Vector *)VectorD_UnWrap((VectorD *)pIVector);
}

Vector *_VectorDC_UnWrap(Vector *pIVector)
{
    return (Vector *)VectorD_UnWrap((VectorD *)pIVector);
}

void Vector_Show(Vector *pIVector)
{
    (*pIVector->Show)(pIVector);
}

void _VectorD_Show(Vector *pIVector)
{
    VectorD_Show((VectorD *)pIVector);
}

void _VectorDC_Show(Vector *pIVector)
{
    VectorDC_Show((VectorDC *)pIVector);
}

Vector *Vector_Full(Vector *pIVector, Num *pNum)
{
    return (*pIVector->Full)(pIVector, pNum);
}

Vector *_VectorD_Full(Vector *pIVector, Num *pNum)
{
    if(instanceof(pNum, NUMD))
    {
        VectorD_Full((VectorD *)pIVector, ((NumD *)pNum)->num);
    }
    else if(instanceof(pNum, NUMDC))
    {
        VectorD_Full2((VectorD *)pIVector, ((NumDC *)pNum)->num);
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pIVector;
}

Vector *_VectorDC_Full(Vector *pVector, Num *pNum)
{
    if(instanceof(pNum,NUMD))
    {
        VectorDC_Full2((VectorDC *)pVector, ((NumD *)pNum)->num);
    }
    else if(instanceof(pNum, NUMDC))
    {
        VectorDC_Full((VectorDC *)pVector, ((NumDC *)pNum)->num);
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pVector;
}

Vector *Vector_Linspace(Vector *pIVector, Num *start, Num *end)
{
    return (*pIVector->Linspace)(pIVector, start, end);
}

Vector *_VectorD_Linspace(Vector *pIVector, Num *pStart, Num *pEnd)
{
    if(instanceof(pStart, NUMD))
    {
        if(instanceof(pEnd, NUMD))
        {
            VectorD_Linspace((VectorD *)pIVector, ((NumD *)pStart)->num, ((NumD *)pEnd)->num);
        }
        else if(instanceof(pEnd, NUMDC))
        {
            VectorD_Linspace3((VectorD *)pIVector, ((NumD *)pStart)->num, ((NumDC *)pEnd)->num);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else if(instanceof(pStart, NUMDC))
    {
        if(instanceof(pEnd, NUMD))
        {
            VectorD_Linspace2((VectorD *)pIVector, ((NumDC *)pStart)->num, ((NumD *)pEnd)->num);
        }
        else if(instanceof(pEnd, NUMDC))
        {
            VectorD_Linspace4((VectorD *)pIVector, ((NumDC *)pStart)->num, ((NumDC *)pEnd)->num);
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

    return pIVector;
}

Vector *_VectorDC_Linspace(Vector *pIVector, Num *pStart, Num *pEnd)
{
    if(instanceof(pEnd, NUMD))
    {
        if(instanceof(pStart, NUMD))
        {
            VectorDC_Linspace1((VectorDC *)pIVector, ((NumD *)pStart)->num, ((NumD *)pEnd)->num);
        }
        else if(instanceof(pStart, NUMDC))
        {
            VectorDC_Linspace3((VectorDC *)pIVector, ((NumDC *)pStart)->num, ((NumD *)pEnd)->num);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else if(instanceof(pEnd, NUMDC))
    {
        if(instanceof(pStart, NUMD))
        {
            VectorDC_Linspace2((VectorDC *)pIVector, ((NumD *)pStart)->num, ((NumDC *)pEnd)->num);
        }
        else if(instanceof(pStart, NUMDC))
        {
            VectorDC_Linspace((VectorDC *)pIVector, ((NumDC *)pStart)->num, ((NumDC *)pEnd)->num);
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

    return pIVector;
}

Vector *Vector_Copy(Vector *pIVector, Vector *pOVector)
{
    return (*pIVector->Copy)(pIVector, pOVector);
}

Vector *_VectorD_Copy(Vector *pIVector, Vector *pOVector)
{
    if(instanceof(pOVector, VECTORD))
    {
        VectorD_Copy((VectorD *)pIVector, (VectorD *)pOVector);
    }
    else if(instanceof(pOVector, VECTORDC))
    {
        VectorD_Copy2((VectorD *)pIVector, (VectorDC *)pOVector);
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pOVector;
}

Vector *_VectorDC_Copy(Vector *pIVector, Vector *pOVector)
{
    if(instanceof(pOVector, VECTORD))
    {
        VectorDC_Copy2((VectorDC *)pIVector, (VectorD *)pOVector);
    }
    else if(instanceof(pOVector, VECTORDC))
    {
        VectorDC_Copy((VectorDC *)pIVector, (VectorDC *)pOVector);
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pOVector;
}

Num *Vector_Max(Vector *pVector, Num *pONum)
{
    return (*pVector->Max)(pVector, pONum);
}

Num *_VectorD_Max(Vector *pVector, Num *pONum)
{
    if(instanceof(pONum, NUMD))
    {
        ((NumD *)pONum)->num = VectorD_Max((VectorD *)pVector);
    }
    else if(instanceof(pONum, NUMDC))
    {
        ((NumDC *)pONum)->num = VectorD_Max2((VectorD *)pVector);
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pONum;
}

Num *_VectorDC_Max(Vector *pVector, Num *pONum)
{
    if(instanceof(pONum, NUMD))
    {
        ((NumD *)pONum)->num = VectorDC_Max2((VectorDC *)pVector);
    }
    else if(instanceof(pONum, NUMDC))
    {
        ((NumDC *)pONum)->num = VectorDC_Max((VectorDC *)pVector);
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pONum;
}

Num *Vector_Min(Vector *pVector, Num *pONum)
{
    return (*pVector->Min)(pVector, pONum);
}

Num *_VectorD_Min(Vector *pVector, Num *pONum)
{
    if(instanceof(pONum, NUMD))
    {
        ((NumD *)pONum)->num = VectorD_Min((VectorD *)pVector);
    }
    else if(instanceof(pONum, NUMDC))
    {
        ((NumDC *)pONum)->num = VectorD_Min2((VectorD *)pVector);
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pONum;
}

Num *_VectorDC_Min(Vector *pVector, Num *pONum)
{
    if(instanceof(pONum, NUMD))
    {
        ((NumD *)pONum)->num = VectorDC_Min2((VectorDC *)pVector);
    }
    else if(instanceof(pONum, NUMDC))
    {
        ((NumDC *)pONum)->num = VectorDC_Min((VectorDC *)pVector);
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pONum;
}

Vector *Vector_AddVector(Vector *pVector1, Vector *pVector2, Vector *pOVector)
{
    return (*pVector1->AddVector)(pVector1, pVector2, pOVector);
}

Vector *_VectorD_AddVector(Vector *pVector1, Vector *pVector2, Vector *pOVector)
{
    if(instanceof(pOVector, VECTORD))
    {
        if(instanceof(pVector2, VECTORD))
        {
            VectorD_AddVectorD((VectorD *)pVector1, (VectorD *)pVector2, (VectorD *)pOVector);
        }
        else if(instanceof(pVector2, VECTORDC))
        {
            VectorD_AddVectorDC2((VectorD *)pVector1, (VectorDC *)pVector2, (VectorD *)pOVector);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else if(instanceof(pOVector, VECTORDC))
    {
        if(instanceof(pVector2, VECTORD))
        {
            VectorD_AddVectorD2((VectorD *)pVector1, (VectorD *)pVector2, (VectorDC *)pOVector);
        }
        else if(instanceof(pVector2, VECTORDC))
        {
            VectorD_AddVectorDC((VectorD *)pVector1, (VectorDC *)pVector2, (VectorDC *)pOVector);
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

    return pOVector;
}

Vector *_VectorDC_AddVector(Vector *pVector1, Vector *pVector2, Vector *pOVector)
{
    if(instanceof(pOVector, VECTORD))
    {
        if(instanceof(pVector2, VECTORD))
        {
            VectorDC_AddVectorD2((VectorDC *)pVector1, (VectorD *)pVector2, (VectorD *)pOVector);
        }
        else if(instanceof(pVector2, VECTORDC))
        {
            VectorDC_AddVectorDC2((VectorDC *)pVector1, (VectorDC *)pVector2, (VectorD *)pOVector);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else if(instanceof(pOVector, VECTORDC))
    {
        if(instanceof(pVector2, VECTORD))
        {
            VectorDC_AddVectorD((VectorDC *)pVector1, (VectorD *)pVector2, (VectorDC *)pOVector);
        }
        else if(instanceof(pVector2, VECTORDC))
        {
            VectorDC_AddVectorDC((VectorDC *)pVector1, (VectorDC *)pVector2, (VectorDC *)pOVector);
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

    return pOVector;
}

Vector *Vector_SubVector(Vector *pVector1, Vector *pVector2, Vector *pOVector)
{
    return (*pVector1->SubVector)(pVector1, pVector2, pOVector);
}

Vector *_VectorD_SubVector(Vector *pVector1, Vector *pVector2, Vector *pOVector)
{
    if(instanceof(pOVector, VECTORD))
    {
        if(instanceof(pVector2, VECTORD))
        {
            VectorD_SubVectorD((VectorD *)pVector1, (VectorD *)pVector2, (VectorD *)pOVector);
        }
        else if(instanceof(pVector2, VECTORDC))
        {
            VectorD_SubVectorDC2((VectorD *)pVector1, (VectorDC *)pVector2, (VectorD *)pOVector);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else if(instanceof(pOVector, VECTORDC))
    {
        if(instanceof(pVector2, VECTORD))
        {
            VectorD_SubVectorD2((VectorD *)pVector1, (VectorD *)pVector2, (VectorDC *)pOVector);
        }
        else if(instanceof(pVector2, VECTORDC))
        {
            VectorD_SubVectorDC((VectorD *)pVector1, (VectorDC *)pVector2, (VectorDC *)pOVector);
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

    return pOVector;
}

Vector *_VectorDC_SubVector(Vector *pVector1, Vector *pVector2, Vector *pOVector)
{
    if(instanceof(pOVector, VECTORD))
    {
        if(instanceof(pVector2, VECTORD))
        {
            VectorDC_SubVectorD2((VectorDC *)pVector1, (VectorD *)pVector2, (VectorD *)pOVector);
        }
        else if(instanceof(pVector2, VECTORDC))
        {
            VectorDC_SubVectorDC2((VectorDC *)pVector1, (VectorDC *)pVector2, (VectorD *)pOVector);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else if(instanceof(pOVector, VECTORDC))
    {
        if(instanceof(pVector2, VECTORD))
        {
            VectorDC_SubVectorD((VectorDC *)pVector1, (VectorD *)pVector2, (VectorDC *)pOVector);
        }
        else if(instanceof(pVector2, VECTORDC))
        {
            VectorDC_SubVectorDC((VectorDC *)pVector1, (VectorDC *)pVector2, (VectorDC *)pOVector);
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

    return pOVector;
}

Vector *Vector_MulNum(Vector *pIVector, Num *pNum, Vector *pOVector)
{
    return (*pIVector->MulNum)(pIVector, pNum, pOVector);
}

Vector *_VectorD_MulNum(Vector *pIVector, Num *pNum, Vector *pOVector)
{
    if(instanceof(pOVector, VECTORD))
    {
        if(instanceof(pNum, NUMD))
        {
            VectorD_MulNumD((VectorD *)pIVector, ((NumD *)pNum)->num, (VectorD *)pOVector);
        }
        else if(instanceof(pNum, NUMDC))
        {
            VectorD_MulNumDC2((VectorD *)pIVector, ((NumDC *)pNum)->num, (VectorD *)pOVector);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else if(instanceof(pOVector, VECTORDC))
    {
        if(instanceof(pNum, NUMD))
        {
            VectorD_MulNumD2((VectorD *)pIVector, ((NumD *)pNum)->num, (VectorDC *)pOVector);
        }
        else if(instanceof(pNum, NUMDC))
        {
            VectorD_MulNumDC((VectorD *)pIVector, ((NumDC *)pNum)->num, (VectorDC *)pOVector);
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

    return pOVector;
}

Vector *_VectorDC_MulNum(Vector *pIVector, Num *pNum, Vector *pOVector)
{
    if(instanceof(pOVector, VECTORD))
    {
        if(instanceof(pNum, NUMD))
        {
            VectorDC_MulNumD2((VectorDC *)pIVector, ((NumD *)pNum)->num, (VectorD *)pOVector);
        }
        else if(instanceof(pNum, NUMDC))
        {
            VectorDC_MulNumDC2((VectorDC *)pIVector, ((NumDC *)pNum)->num, (VectorD *)pOVector);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else if(instanceof(pOVector, VECTORDC))
    {
        if(instanceof(pNum, NUMD))
        {
            VectorDC_MulNumD((VectorDC *)pIVector, ((NumD *)pNum)->num, (VectorDC *)pOVector);
        }
        else if(instanceof(pNum, NUMDC))
        {
            VectorDC_MulNumDC((VectorDC *)pIVector, ((NumDC *)pNum)->num, (VectorDC *)pOVector);
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

    return pOVector;
}

Num *Vector_MulVector(Vector *pVector1, Vector *pVector2, Num *pONum)
{
    return (*pVector1->MulVector)(pVector1, pVector2, pONum);
}

Num *_VectorD_MulVector(Vector *pVector1, Vector *pVector2, Num *pONum)
{
    if(instanceof(pONum, NUMD))
    {
        if(instanceof(pVector2, VECTORD))
        {
            ((NumD *)pONum)->num = VectorD_MulVectorD((VectorD *)pVector1, (VectorD *)pVector2);
        }
        else if(instanceof(pVector2, VECTORDC))
        {
            ((NumD *)pONum)->num = VectorD_MulVectorDC2((VectorD *)pVector1, (VectorDC *)pVector2);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else if(instanceof(pONum, NUMDC))
    {
        if(instanceof(pVector2, VECTORD))
        {
            #ifdef _MSC_VER
            ((NumDC *)pONum)->num._Val[0] = VectorD_MulVectorD((VectorD *)pVector1, (VectorD *)pVector2);
            ((NumDC *)pONum)->num._Val[1] = 0;
            #else
            ((NumDC *)pONum)->num = VectorD_MulVectorD((VectorD *)pVector1, (VectorD *)pVector2);
            #endif /*_MSC_VER*/
        }
        else if(instanceof(pVector2, VECTORDC))
        {
            ((NumDC *)pONum)->num = VectorD_MulVectorDC((VectorD *)pVector1, (VectorDC *)pVector2);
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

Num *_VectorDC_MulVector(Vector *pVector1, Vector *pVector2, Num *pONum)
{
    if(instanceof(pONum, NUMD))
    {
        if(instanceof(pVector2, VECTORD))
        {
            ((NumD *)pONum)->num = VectorDC_MulVectorD2((VectorDC *)pVector1, (VectorD *)pVector2);
        }
        else if(instanceof(pVector2, VECTORDC))
        {
            ((NumD *)pONum)->num = VectorDC_MulVectorDC2((VectorDC *)pVector1, (VectorDC *)pVector2);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else if(instanceof(pONum, NUMDC))
    {
        if(instanceof(pVector2, VECTORD))
        {
            #ifdef _MSC_VER
            ((NumDC *)pONum)->num = VectorDC_MulVectorD((VectorDC *)pVector1, (VectorD *)pVector2);
            #else
            ((NumDC *)pONum)->num = VectorDC_MulVectorD((VectorDC *)pVector1, (VectorD *)pVector2);
            #endif /*_MSC_VER*/
        }
        else if(instanceof(pVector2, VECTORDC))
        {
            ((NumDC *)pONum)->num = VectorDC_MulVectorDC((VectorDC *)pVector1, (VectorDC *)pVector2);
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

Vector *Vector_EMulVector(Vector *pVector1, Vector *pVector2, Vector *pOVector)
{
    return (*pVector1->EMulVector)(pVector1, pVector2, pOVector);
}

Vector *_VectorD_EMulVector(Vector *pIVector1, Vector *pIVector2, Vector *pOVector)
{
    if(instanceof(pOVector, VECTORD))
    {
        if(instanceof(pIVector2, VECTORD))
        {
            VectorD_EMulVectorD((VectorD *)pIVector1, (VectorD *)pIVector2, (VectorD *)pOVector);
        }
        else if(instanceof(pIVector2, VECTORDC))
        {
            VectorD_EMulVectorDC2((VectorD *)pIVector1, (VectorDC *)pIVector2, (VectorD *)pOVector);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else if(instanceof(pOVector, VECTORDC))
    {
        if(instanceof(pIVector2, VECTORD))
        {
            VectorD_EMulVectorD2((VectorD *)pIVector1, (VectorD *)pIVector2, (VectorDC *)pOVector);
        }
        else if(instanceof(pIVector2, VECTORDC))
        {
            VectorD_EMulVectorDC((VectorD *)pIVector1, (VectorDC *)pIVector2, (VectorDC *)pOVector);
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

    return pOVector;
}

Vector *_VectorDC_EMulVector(Vector *pIVector1, Vector *pIVector2, Vector *pOVector)
{
    if(instanceof(pOVector, VECTORD))
    {
        if(instanceof(pIVector2, VECTORD))
        {
            VectorDC_EMulVectorD2((VectorDC *)pIVector1, (VectorD *)pIVector2, (VectorD *)pOVector);
        }
        else if(instanceof(pIVector2, VECTORDC))
        {
            VectorDC_EMulVectorDC2((VectorDC *)pIVector1, (VectorDC *)pIVector2, (VectorD *)pOVector);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else if(instanceof(pOVector, VECTORDC))
    {
        if(instanceof(pIVector2, VECTORD))
        {
            VectorDC_EMulVectorD((VectorDC *)pIVector1, (VectorD *)pIVector2, (VectorDC *)pOVector);
        }
        else if(instanceof(pIVector2, VECTORDC))
        {
            VectorDC_EMulVectorDC((VectorDC *)pIVector1, (VectorDC *)pIVector2, (VectorDC *)pOVector);
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

    return pOVector;
}

Num *Vector_CMulVector(Vector *pVector1, Vector *pVector2, Num *pONum)
{
    return (*pVector1->CMulVector)(pVector1, pVector2, pONum);
}

Num *_VectorD_CMulVector(Vector *pVector1, Vector *pVector2, Num *pONum)
{
    if(instanceof(pONum, NUMD))
    {
        if(instanceof(pVector2, VECTORD))
        {
            ((NumD *)pONum)->num=VectorD_CMulVectorD((VectorD *)pVector1, (VectorD *)pVector2);
        }
        else if(instanceof(pVector2, VECTORDC))
        {
            ((NumD *)pONum)->num=VectorD_CMulVectorDC2((VectorD *)pVector1, (VectorDC *)pVector2);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else if(instanceof(pONum, NUMDC))
    {
        if(instanceof(pVector2, VECTORD))
        {
            ((NumDC *)pONum)->num=VectorD_CMulVectorD2((VectorD *)pVector1, (VectorD *)pVector2);
        }
        else if(instanceof(pVector2, VECTORDC))
        {
            ((NumDC *)pONum)->num=VectorD_CMulVectorDC((VectorD *)pVector1, (VectorDC *)pVector2);
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

Num *_VectorDC_CMulVector(Vector *pVector1, Vector *pVector2, Num *pONum)
{
    if(instanceof(pONum, NUMD))
    {
        if(instanceof(pVector2, VECTORD))
        {
            ((NumD *)pONum)->num=VectorDC_CMulVectorD2((VectorDC *)pVector1, (VectorD *)pVector2);
        }
        else if(instanceof(pVector2, VECTORDC))
        {
            ((NumD *)pONum)->num=VectorDC_CMulVectorDC2((VectorDC *)pVector1, (VectorDC *)pVector2);
        }
        else
        {
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
        }
    }
    else if(instanceof(pONum, NUMDC))
    {
        if(instanceof(pVector2, VECTORD))
        {
            ((NumDC *)pONum)->num=VectorDC_CMulVectorD((VectorDC *)pVector1, (VectorD *)pVector2);
        }
        else if(instanceof(pVector2, VECTORDC))
        {
            ((NumDC *)pONum)->num=VectorDC_CMulVectorDC((VectorDC *)pVector1, (VectorDC *)pVector2);
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

Num *Vector_Norm(Vector *pIVector, Num *pONum)
{
    return (*pIVector->Norm)(pIVector, pONum);
}

Num *_VectorD_Norm(Vector *pIVector, Num *pONum)
{
    if(instanceof(pONum, NUMD))
    {
        ((NumD *)pONum)->num = VectorD_Norm((VectorD *)pIVector);
    }
    else if(instanceof(pONum, NUMDC))
    {
        ((NumDC *)pONum)->num = VectorD_Norm2((VectorD *)pIVector);
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pONum;
}

Num *_VectorDC_Norm(Vector *pIVector, Num *pONum)
{
    if(instanceof(pONum, NUMD))
    {
        ((NumD *)pONum)->num = VectorDC_Norm((VectorDC *)pIVector);
    }
    else if(instanceof(pONum, NUMDC))
    {
        ((NumDC *)pONum)->num = VectorDC_Norm2((VectorDC *)pIVector);
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pONum;
}

Num *Vector_NormSq(Vector *pIVector, Num *pONum)
{
    return (*pIVector->NormSq)(pIVector, pONum);
}

Num *_VectorD_NormSq(Vector *pIVector, Num *pONum)
{
    if(instanceof(pONum, NUMD))
    {
        ((NumD *)pONum)->num = VectorD_NormSq((VectorD *)pIVector);
    }
    else if(instanceof(pONum, NUMDC))
    {
        ((NumDC *)pONum)->num = VectorD_NormSq2((VectorD *)pIVector);
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pONum;
}

Num *_VectorDC_NormSq(Vector *pIVector, Num *pONum)
{
    if(instanceof(pONum, NUMD))
    {
        ((NumD *)pONum)->num = VectorDC_NormSq((VectorDC *)pIVector);
    }
    else if(instanceof(pONum, NUMDC))
    {
        ((NumDC *)pONum)->num = VectorDC_NormSq2((VectorDC *)pIVector);
    }
    else
    {
        fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
    }

    return pONum;
}
