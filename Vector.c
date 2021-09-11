#include "Vector.h"

Vector *Vector_Del(Vector *pIVector)
{
    return (*pIVector->Del)(pIVector);
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
    return (*pIVector->UnWrap)(pIVector);
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
    switch(pNum->type)
    {
        case NUMD:
            VectorD_Full((VectorD *)pIVector, ((NumD *)pNum)->num);
            break;
        case NUMDC:
            VectorD_Full2((VectorD *)pIVector, ((NumDC *)pNum)->num);
            break;
        default:
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
    }

    return pIVector;
}

Vector *_VectorDC_Full(Vector *pVector, Num *pNum)
{
    switch(pNum->type)
    {
        case NUMD:
            VectorDC_Full2((VectorDC *)pVector, ((NumD *)pNum)->num);
            break;
        case NUMDC:
            VectorDC_Full((VectorDC *)pVector, ((NumDC *)pNum)->num);
            break;
        default:
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
    switch(pStart->type)
    {
        case NUMD:
            switch(pEnd->type)
            {
                case NUMD:
                    VectorD_Linspace((VectorD *)pIVector, ((NumD *)pStart)->num, ((NumD *)pEnd)->num);
                    break;
                case NUMDC:
                    VectorD_Linspace3((VectorD *)pIVector, ((NumD *)pStart)->num, ((NumDC *)pEnd)->num);
                    break;
                default:
                    fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                        __LINE__);
            }
            break;
        case NUMDC:
            switch(pEnd->type)
            {
                case NUMD:
                    VectorD_Linspace2((VectorD *)pIVector, ((NumDC *)pStart)->num, ((NumD *)pEnd)->num);
                    break;
                case NUMDC:
                    VectorD_Linspace4((VectorD *)pIVector, ((NumDC *)pStart)->num, ((NumDC *)pEnd)->num);
                    break;
                default:
                    fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                        __LINE__);
            }
            break;
        default:
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
    }

    return pIVector;
}

Vector *_VectorDC_Linspace(Vector *pIVector, Num *start, Num *end)
{
    switch(end->type)
    {
        case NUMD:
            switch(start->type)
            {
                case NUMD:
                    VectorDC_Linspace1((VectorDC *)pIVector, ((NumD *)start)->num, ((NumD *)end)->num);
                    break;
                case NUMDC:
                    VectorDC_Linspace3((VectorDC *)pIVector, ((NumDC *)start)->num, ((NumD *)end)->num);
                    break;
                default:
                    fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                        __LINE__);
            }
            break;
        case NUMDC:
            switch(start->type)
            {
                case NUMD:
                    VectorDC_Linspace2((VectorDC *)pIVector, ((NumD *)start)->num, ((NumDC *)end)->num);
                    break;
                case NUMDC:
                    VectorDC_Linspace((VectorDC *)pIVector, ((NumDC *)start)->num, ((NumDC *)end)->num);
                    break;
                default:
                    fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                        __LINE__);
            }
            break;
        default:
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
    switch(pOVector->type)
    {
        case VECTORD:
            VectorD_Copy((VectorD *)pIVector, (VectorD *)pOVector);
            break;
        case VECTORDC:
            VectorD_Copy2((VectorD *)pIVector, (VectorDC *)pOVector);
            break;
        default:
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
    }

    return pOVector;
}

Vector *_VectorDC_Copy(Vector *pIVector, Vector *pOVector)
{
    switch(pOVector->type)
    {
        case VECTORD:
            VectorDC_Copy2((VectorDC *)pIVector, (VectorD *)pOVector);
            break;
        case VECTORDC:
            VectorDC_Copy((VectorDC *)pIVector, (VectorDC *)pOVector);
            break;
        default:
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
    switch(pONum->type)
    {
        case NUMD:
            ((NumD *)pONum)->num = VectorD_Max((VectorD *)pVector);
            break;
        case NUMDC:
            ((NumDC *)pONum)->num = VectorD_Max2((VectorD *)pVector);
            break;
        default:
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
    }

    return pONum;
}

Num *_VectorDC_Max(Vector *pVector, Num *pONum)
{
    switch(pONum->type)
    {
        case NUMD:
            ((NumD *)pONum)->num = VectorDC_Max2((VectorDC *)pVector);
            break;
        case NUMDC:
            ((NumDC *)pONum)->num = VectorDC_Max((VectorDC *)pVector);
            break;
        default:
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
    switch(pONum->type)
    {
        case NUMD:
            ((NumD *)pONum)->num = VectorD_Min((VectorD *)pVector);
            break;
        case NUMDC:
            ((NumDC *)pONum)->num = VectorD_Min2((VectorD *)pVector);
            break;
        default:
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
    }

    return pONum;
}

Num *_VectorDC_Min(Vector *pVector, Num *pONum)
{
    switch(pONum->type)
    {
        case NUMD:
            ((NumD *)pONum)->num = VectorDC_Min2((VectorDC *)pVector);
            break;
        case NUMDC:
            ((NumDC *)pONum)->num = VectorDC_Min((VectorDC *)pVector);
            break;
        default:
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
    switch(pOVector->type)
    {
        case VECTORD:
            switch(pVector2->type)
            {
                case VECTORD:
                    VectorD_AddVectorD((VectorD *)pVector1, (VectorD *)pVector2, (VectorD *)pOVector);
                    break;
                case VECTORDC:
                    VectorD_AddVectorDC2((VectorD *)pVector1, (VectorDC *)pVector2, (VectorD *)pOVector);
                    break;
                default:
                    fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                        __LINE__);
            }
            break;
        case VECTORDC:
            switch(pVector2->type)
            {
                case VECTORD:
                    VectorD_AddVectorD2((VectorD *)pVector1, (VectorD *)pVector2, (VectorDC *)pOVector);
                    break;
                case VECTORDC:
                    VectorD_AddVectorDC((VectorD *)pVector1, (VectorDC *)pVector2, (VectorDC *)pOVector);
                    break;
                default:
                    fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                        __LINE__);
            }
            break;
        default:
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
    }

    return pOVector;
}

Vector *_VectorDC_AddVector(Vector *pVector1, Vector *pVector2, Vector *pOVector)
{
    switch(pOVector->type)
    {
        case VECTORD:
            switch(pVector2->type)
            {
                case VECTORD:
                    VectorDC_AddVectorD2((VectorDC *)pVector1, (VectorD *)pVector2, (VectorD *)pOVector);
                    break;
                case VECTORDC:
                    VectorDC_AddVectorDC2((VectorDC *)pVector1, (VectorDC *)pVector2, (VectorD *)pOVector);
                    break;
                default:
                    fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                        __LINE__);
            }
            break;
        case VECTORDC:
            switch(pVector2->type)
            {
                case VECTORD:
                    VectorDC_AddVectorD((VectorDC *)pVector1, (VectorD *)pVector2, (VectorDC *)pOVector);
                    break;
                case VECTORDC:
                    VectorDC_AddVectorDC((VectorDC *)pVector1, (VectorDC *)pVector2, (VectorDC *)pOVector);
                    break;
                default:
                    fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                        __LINE__);
            }
            break;
        default:
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
    switch(pOVector->type)
    {
        case VECTORD:
            switch(pVector2->type)
            {
                case VECTORD:
                    VectorD_SubVectorD((VectorD *)pVector1, (VectorD *)pVector2, (VectorD *)pOVector);
                    break;
                case VECTORDC:
                    VectorD_SubVectorDC2((VectorD *)pVector1, (VectorDC *)pVector2, (VectorD *)pOVector);
                    break;
                default:
                    fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                        __LINE__);
            }
            break;
        case VECTORDC:
            switch(pVector2->type)
            {
                case VECTORD:
                    VectorD_SubVectorD2((VectorD *)pVector1, (VectorD *)pVector2, (VectorDC *)pOVector);
                    break;
                case VECTORDC:
                    VectorD_SubVectorDC((VectorD *)pVector1, (VectorDC *)pVector2, (VectorDC *)pOVector);
                    break;
                default:
                    fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                        __LINE__);
            }
            break;
        default:
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
    }

    return pOVector;
}

Vector *_VectorDC_SubVector(Vector *pVector1, Vector *pVector2, Vector *pOVector)
{
    switch(pOVector->type)
    {
        case VECTORD:
            switch(pVector2->type)
            {
                case VECTORD:
                    VectorDC_SubVectorD2((VectorDC *)pVector1, (VectorD *)pVector2, (VectorD *)pOVector);
                    break;
                case VECTORDC:
                    VectorDC_SubVectorDC2((VectorDC *)pVector1, (VectorDC *)pVector2, (VectorD *)pOVector);
                    break;
                default:
                    fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                        __LINE__);
            }
            break;
        case VECTORDC:
            switch(pVector2->type)
            {
                case VECTORD:
                    VectorDC_SubVectorD((VectorDC *)pVector1, (VectorD *)pVector2, (VectorDC *)pOVector);
                    break;
                case VECTORDC:
                    VectorDC_SubVectorDC((VectorDC *)pVector1, (VectorDC *)pVector2, (VectorDC *)pOVector);
                    break;
                default:
                    fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                        __LINE__);
            }
            break;
        default:
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
    switch(pOVector->type)
    {
        case VECTORD:
            switch(pNum->type)
            {
                case NUMD:
                    VectorD_MulNumD((VectorD *)pIVector, ((NumD *)pNum)->num, (VectorD *)pOVector);
                    break;
                case NUMDC:
                    VectorD_MulNumDC2((VectorD *)pIVector, ((NumDC *)pNum)->num, (VectorD *)pOVector);
                    break;
                default:
                    fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                        __LINE__);
            }
            break;
        case VECTORDC:
            switch(pNum->type)
            {
                case NUMD:
                    VectorD_MulNumD2((VectorD *)pIVector, ((NumD *)pNum)->num, (VectorDC *)pOVector);
                    break;
                case NUMDC:
                    VectorD_MulNumDC((VectorD *)pIVector, ((NumDC *)pNum)->num, (VectorDC *)pOVector);
                    break;
                default:
                    fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                        __LINE__);
            }
            break;
        default:
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
    }

    return pOVector;
}

Vector *_VectorDC_MulNum(Vector *pIVector, Num *pNum, Vector *pOVector)
{
    switch(pOVector->type)
    {
        case VECTORD:
            switch(pNum->type)
            {
                case NUMD:
                    VectorDC_MulNumD2((VectorDC *)pIVector, ((NumD *)pNum)->num, (VectorD *)pOVector);
                    break;
                case NUMDC:
                    VectorDC_MulNumDC2((VectorDC *)pIVector, ((NumDC *)pNum)->num, (VectorD *)pOVector);
                    break;
                default:
                    fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                        __LINE__);
            }
            break;
        case VECTORDC:
            switch(pNum->type)
            {
                case NUMD:
                    VectorDC_MulNumD((VectorDC *)pIVector, ((NumD *)pNum)->num, (VectorDC *)pOVector);
                    break;
                case NUMDC:
                    VectorDC_MulNumDC((VectorDC *)pIVector, ((NumDC *)pNum)->num, (VectorDC *)pOVector);
                    break;
                default:
                    fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                        __LINE__);
            }
            break;
        default:
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
    switch(pONum->type)
    {
        case NUMD:
            switch(pVector2->type)
            {
                case VECTORD:
                    ((NumD *)pONum)->num = VectorD_MulVectorD((VectorD *)pVector1, (VectorD *)pVector2);
                    break;
                case VECTORDC:
                    ((NumD *)pONum)->num = VectorD_MulVectorDC2((VectorD *)pVector1, (VectorDC *)pVector2);
                    break;
                default:
                    fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                        __LINE__);
            }
            break;
        case NUMDC:
            switch(pVector2->type)
            {
                case VECTORD:
                    #ifdef _MSC_VER
                    ((NumDC *)pONum)->num._Val[0] = VectorD_MulVectorD((VectorD *)pVector1, (VectorD *)pVector2);
                    ((NumDC *)pONum)->num._Val[1] = 0;
                    #else
                    ((NumDC *)pONum)->num = VectorD_MulVectorD((VectorD *)pVector1, (VectorD *)pVector2);
                    #endif /*_MSC_VER*/
                    break;
                case VECTORDC:
                    ((NumDC *)pONum)->num = VectorD_MulVectorDC((VectorD *)pVector1, (VectorDC *)pVector2);
                    break;
                default:
                    fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                        __LINE__);
            }
            break;
        default:
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
    }

    return pONum;
}

Num *_VectorDC_MulVector(Vector *pVector1, Vector *pVector2, Num *pONum)
{
    switch(pONum->type)
    {
        case NUMD:
            switch(pVector2->type)
            {
                case VECTORD:
                    ((NumD *)pONum)->num = VectorDC_MulVectorD2((VectorDC *)pVector1, (VectorD *)pVector2);
                    break;
                case VECTORDC:
                    ((NumD *)pONum)->num = VectorDC_MulVectorDC2((VectorDC *)pVector1, (VectorDC *)pVector2);
                    break;
                default:
                    fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                        __LINE__);
            }
            break;
        case NUMDC:
            switch(pVector2->type)
            {
                case VECTORD:
                    #ifdef _MSC_VER
                    ((NumDC *)pONum)->num = VectorDC_MulVectorD((VectorDC *)pVector1, (VectorD *)pVector2);
                    #else
                    ((NumDC *)pONum)->num = VectorDC_MulVectorD((VectorDC *)pVector1, (VectorD *)pVector2);
                    #endif /*_MSC_VER*/
                    break;
                case VECTORDC:
                    ((NumDC *)pONum)->num = VectorDC_MulVectorDC((VectorDC *)pVector1, (VectorDC *)pVector2);
                    break;
                default:
                    fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                        __LINE__);
            }
            break;
        default:
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
    switch(pOVector->type)
    {
        case VECTORD:
            switch(pIVector2->type)
            {
                case VECTORD:
                    VectorD_EMulVectorD((VectorD *)pIVector1, (VectorD *)pIVector2, (VectorD *)pOVector);
                    break;
                case VECTORDC:
                    VectorD_EMulVectorDC2((VectorD *)pIVector1, (VectorDC *)pIVector2, (VectorD *)pOVector);
                    break;
                default:
                    fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                        __LINE__);
            }
            break;
        case VECTORDC:
            switch(pIVector2->type)
            {
                case VECTORD:
                    VectorD_EMulVectorD2((VectorD *)pIVector1, (VectorD *)pIVector2, (VectorDC *)pOVector);
                    break;
                case VECTORDC:
                    VectorD_EMulVectorDC((VectorD *)pIVector1, (VectorDC *)pIVector2, (VectorDC *)pOVector);
                    break;
                default:
                    fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                        __LINE__);
            }
            break;
        default:
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
    }

    return pOVector;
}

Vector *_VectorDC_EMulVector(Vector *pIVector1, Vector *pIVector2, Vector *pOVector)
{
    switch(pOVector->type)
    {
        case VECTORD:
            switch(pIVector2->type)
            {
                case VECTORD:
                    VectorDC_EMulVectorD2((VectorDC *)pIVector1, (VectorD *)pIVector2, (VectorD *)pOVector);
                    break;
                case VECTORDC:
                    VectorDC_EMulVectorDC2((VectorDC *)pIVector1, (VectorDC *)pIVector2, (VectorD *)pOVector);
                    break;
                default:
                    fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                        __LINE__);
            }
            break;
        case VECTORDC:
            switch(pIVector2->type)
            {
                case VECTORD:
                    VectorDC_EMulVectorD((VectorDC *)pIVector1, (VectorD *)pIVector2, (VectorDC *)pOVector);
                    break;
                case VECTORDC:
                    VectorDC_EMulVectorDC((VectorDC *)pIVector1, (VectorDC *)pIVector2, (VectorDC *)pOVector);
                    break;
                default:
                    fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                        __LINE__);
            }
            break;
        default:
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
    switch(pONum->type)
    {
        case NUMD:
            switch(pVector2->type)
            {
                case VECTORD:
                    ((NumD *)pONum)->num=VectorD_CMulVectorD((VectorD *)pVector1, (VectorD *)pVector2);
                    break;
                case VECTORDC:
                    ((NumD *)pONum)->num=VectorD_CMulVectorDC2((VectorD *)pVector1, (VectorDC *)pVector2);
                    break;
                default:
                    fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                        __LINE__);
            }
            break;
        case NUMDC:
            switch(pVector2->type)
            {
                case VECTORD:
                    ((NumDC *)pONum)->num=VectorD_CMulVectorD2((VectorD *)pVector1, (VectorD *)pVector2);
                    break;
                case VECTORDC:
                    ((NumDC *)pONum)->num=VectorD_CMulVectorDC((VectorD *)pVector1, (VectorDC *)pVector2);
                    break;
                default:
                    fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                        __LINE__);
            }
            break;
        default:
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
    }

    return pONum;
}

Num *_VectorDC_CMulVector(Vector *pVector1, Vector *pVector2, Num *pONum)
{
    switch(pONum->type)
    {
        case NUMD:
            switch(pVector2->type)
            {
                case VECTORD:
                    ((NumD *)pONum)->num=VectorDC_CMulVectorD2((VectorDC *)pVector1, (VectorD *)pVector2);
                    break;
                case VECTORDC:
                    ((NumD *)pONum)->num=VectorDC_CMulVectorDC2((VectorDC *)pVector1, (VectorDC *)pVector2);
                    break;
                default:
                    fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                        __LINE__);
            }
            break;
        case NUMDC:
            switch(pVector2->type)
            {
                case VECTORD:
                    ((NumDC *)pONum)->num=VectorDC_CMulVectorD((VectorDC *)pVector1, (VectorD *)pVector2);
                    break;
                case VECTORDC:
                    ((NumDC *)pONum)->num=VectorDC_CMulVectorDC((VectorDC *)pVector1, (VectorDC *)pVector2);
                    break;
                default:
                    fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                        __LINE__);
            }
            break;
        default:
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
    switch(pONum->type)
    {
        case NUMD:
            ((NumD *)pONum)->num = VectorD_Norm((VectorD *)pIVector);
            break;
        case NUMDC:
            ((NumDC *)pONum)->num = VectorD_Norm2((VectorD *)pIVector);
            break;
        default:
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
    }

    return pONum;
}

Num *_VectorDC_Norm(Vector *pIVector, Num *pONum)
{
    switch(pONum->type)
    {
        case NUMD:
            ((NumD *)pONum)->num = VectorDC_Norm((VectorDC *)pIVector);
            break;
        case NUMDC:
            ((NumDC *)pONum)->num = VectorDC_Norm2((VectorDC *)pIVector);
            break;
        default:
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
    switch(pONum->type)
    {
        case NUMD:
            ((NumD *)pONum)->num = VectorD_NormSq((VectorD *)pIVector);
            break;
        case NUMDC:
            ((NumDC *)pONum)->num = VectorD_NormSq2((VectorD *)pIVector);
            break;
        default:
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
    }

    return pONum;
}

Num *_VectorDC_NormSq(Vector *pIVector, Num *pONum)
{
    switch(pONum->type)
    {
        case NUMD:
            ((NumD *)pONum)->num = VectorDC_NormSq((VectorDC *)pIVector);
            break;
        case NUMDC:
            ((NumDC *)pONum)->num = VectorDC_NormSq2((VectorDC *)pIVector);
            break;
        default:
            fprintf(stderr, "Error: There is no such type. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
                __LINE__);
    }

    return pONum;
}
