#include "Vector.h"

extern Vector *_VectorD_Del(Vector *pIVector);
extern Vector *_VectorD_UnWrap(Vector *pIVector);
extern void _VectorD_Show(Vector *pIVector);
extern Vector *_VectorD_Full(Vector *pIVector, Num *pNum);
extern Vector *_VectorD_Linspace(Vector *pIVector, Num *start, Num *end);
extern Vector *_VectorD_Copy(Vector *pIVector, Vector *pOVector);
extern Num *_VectorD_Max(Vector *pVector, Num *pONum);
extern Num *_VectorD_Min(Vector *pVector, Num *pONum);
extern Vector *_VectorD_AddVector(Vector *pVector1, Vector *pVector2, Vector *pOVector);
extern Vector *_VectorD_SubVector(Vector *pVector1, Vector *pVector2, Vector *pOVector);
extern Vector *_VectorD_MulNum(Vector *pIVector, Num *pNum, Vector *pOVector);
extern Num *_VectorD_MulVector(Vector *pVector1, Vector *pVector2, Num *pONum);
extern Vector *_VectorD_EMulVector(Vector *pVector1, Vector *pVector2, Vector *pOVector);
extern Num *_VectorD_CMulVector(Vector *pVector1, Vector *pVector2, Num *pONum);
extern Num *_VectorD_Norm(Vector *pIVector, Num *pONum);
extern Num *_VectorD_NormSq(Vector *pIVector, Num *pONum);

VectorD *VectorD_New(int size)
{
    VectorD *pNewVectorD;

    pNewVectorD = (VectorD *)malloc(sizeof(VectorD));
    if(pNewVectorD == NULL)
    {
        fprintf(stderr, "Error: Can't create a new vectorD. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
        return NULL;
    }
    pNewVectorD->vector = (double *)malloc(sizeof(double)*size);
    if(pNewVectorD->vector == NULL)
    {
        free(pNewVectorD);
        fprintf(stderr, "Error: Can't create a new vectorD. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
        return NULL;
    }
    pNewVectorD->size = size;
    pNewVectorD->inc = 1;
    ((Object *)pNewVectorD)->type = VECTORD;
    pNewVectorD->parent.Del = _VectorD_Del;
    pNewVectorD->parent.UnWrap = _VectorD_UnWrap;
    pNewVectorD->parent.Show = _VectorD_Show;
    pNewVectorD->parent.Full = _VectorD_Full;
    pNewVectorD->parent.Linspace = _VectorD_Linspace;
    pNewVectorD->parent.Copy = _VectorD_Copy;
    pNewVectorD->parent.Max = _VectorD_Max;
    pNewVectorD->parent.Min = _VectorD_Min;
    pNewVectorD->parent.AddVector = _VectorD_AddVector;
    pNewVectorD->parent.SubVector = _VectorD_SubVector;
    pNewVectorD->parent.MulNum = _VectorD_MulNum;
    pNewVectorD->parent.EMulVector = _VectorD_EMulVector;
    pNewVectorD->parent.MulVector = _VectorD_MulVector;
    pNewVectorD->parent.CMulVector = _VectorD_CMulVector;
    pNewVectorD->parent.Norm = _VectorD_Norm;
    pNewVectorD->parent.NormSq = _VectorD_NormSq;

    return pNewVectorD;
}

VectorD *VectorD_Del(VectorD *pIVector)
{
    if (pIVector != NULL)
    {
        free(pIVector->vector);
        free(pIVector);
    }

    return NULL;
}

VectorD *VectorD_Wrap(double *work, int size, int inc)
{
    if(work == NULL)
    {
        return NULL;
    }

    VectorD *pNewVectorD;

    pNewVectorD = (VectorD *)malloc(sizeof(VectorD));
    if(pNewVectorD == NULL)
    {
        fprintf(stderr, "Error: Can't create a new vectorD. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
        return NULL;
    }
    pNewVectorD->vector = work;
    pNewVectorD->size = size;
    pNewVectorD->inc = inc;
    ((Object *)pNewVectorD)->type = VECTORD;
    pNewVectorD->parent.Del = _VectorD_Del;
    pNewVectorD->parent.UnWrap = _VectorD_UnWrap;
    pNewVectorD->parent.Show = _VectorD_Show;
    pNewVectorD->parent.Full = _VectorD_Full;
    pNewVectorD->parent.Linspace = _VectorD_Linspace;
    pNewVectorD->parent.Copy = _VectorD_Copy;
    pNewVectorD->parent.Max = _VectorD_Max;
    pNewVectorD->parent.Min = _VectorD_Min;
    pNewVectorD->parent.AddVector = _VectorD_AddVector;
    pNewVectorD->parent.SubVector = _VectorD_SubVector;
    pNewVectorD->parent.MulNum = _VectorD_MulNum;
    pNewVectorD->parent.MulVector = _VectorD_MulVector;
    pNewVectorD->parent.EMulVector = _VectorD_EMulVector;
    pNewVectorD->parent.CMulVector = _VectorD_CMulVector;
    pNewVectorD->parent.Norm = _VectorD_Norm;
    pNewVectorD->parent.NormSq = _VectorD_NormSq;

    return pNewVectorD;
}

VectorD *VectorD_UnWrap(VectorD *pIVector)
{
    if (pIVector != NULL)
    {
        free(pIVector);
    }

    return NULL;
}

void VectorD_Show(VectorD *pIVector)
{
    register int i, size, inc;

    size = pIVector->size;
    inc = pIVector->inc;
    printf("VectorD_Show: \n vector:[");
    for(i=0;i<size;i++)
    {
        printf("%f, ", *(pIVector->vector+i*inc));
    }
    printf("] \n size:%d \n inc:%d \n\n", pIVector->size, pIVector->inc);
}

VectorD *VectorD_Full(VectorD *pIVector, double num)
{
    register int i;
    int size;
    register double *vector;

    size = pIVector->size;
    vector = pIVector->vector;

    for(i=0;i<size;i++)
    {
        *(vector+i) = num;
    }

    return pIVector;
}

#ifdef _MSC_VER
VectorD *VectorD_Full2(VectorD *pIVector, _Dcomplex num)
#else
VectorD *VectorD_Full2(VectorD *pIVector, double complex num)
#endif /*_MSC_VER*/
{
    register int i;
    int size;
    register double *vector;

    size = pIVector->size;
    vector = pIVector->vector;

    for(i=0;i<size;i++)
    {
        #ifdef _MSC_VER
        *(vector+i) = num._Val[0];
        #else
        *(vector+i) = creal(num);
        #endif /*_MSC_VER*/
    }

    return pIVector;
}

VectorD *VectorD_Linspace(VectorD *pIVector, double start, double end)
{
    register int i;
    int size;
    double step;
    register double *vector;

    size = pIVector->size;
    vector = pIVector->vector;
    step = (end - start) / ((double)size-1);

    for(i=0;i<size;i++)
    {
        *(vector+i) = start + step*i;
    }

    return pIVector;
}

#ifdef _MSC_VER
VectorD *VectorD_Linspace2(VectorD *pIVector, _Dcomplex start, double end)
#else
VectorD *VectorD_Linspace2(VectorD *pIVector, double complex start, double end)
#endif /*_MSC_VER*/
{
    register int i;
    int size;
    double step;
    register double *vector;

    size = pIVector->size;
    vector = pIVector->vector;
    #ifdef _MSC_VER
    step = (end - start._Val[0]) / ((double)size-1);
    #else
    step = (end - creal(start)) / (size-1);
    #endif /*_MSC_VER*/

    for(i=0;i<size;i++)
    {
        #ifdef _MSC_VER
        *(vector+i) = start._Val[0] + step*i;
        #else
        *(vector+i) = creal(start) + step*i;
        #endif
    }

    return pIVector;
}

#ifdef _MSC_VER
VectorD *VectorD_Linspace3(VectorD *pIVector, double start, _Dcomplex end)
#else
VectorD *VectorD_Linspace3(VectorD *pIVector, double start, double complex end)
#endif /*_MSC_VER*/
{
    register int i;
    int size;
    double step;
    register double *vector;

    size = pIVector->size;
    vector = pIVector->vector;
    #ifdef _MSC_VER
    step = (end._Val[0] - start) / ((double)size-1);
    #else
    step = (creal(end) - start) / (size-1);
    #endif /*_MSC_VER*/

    for(i=0;i<size;i++)
    {
        *(vector+i) = start + step*i;
    }

    return pIVector;
}

#ifdef _MSC_VER
VectorD *VectorD_Linspace4(VectorD *pIVector, _Dcomplex start, _Dcomplex end)
#else
VectorD *VectorD_Linspace4(VectorD *pIVector, double complex start, double complex end)
#endif /*_MSC_VER*/
{
    register int i;
    int size;
    double step;
    register double *vector;

    size = pIVector->size;
    vector = pIVector->vector;
    #ifdef _MSC_VER
    step = (end._Val[0] - start._Val[0]) / ((double)size-1);
    #else
    step = (creal(end) - creal(start)) / (size-1);
    #endif /*_MSC_VER*/

    for(i=0;i<size;i++)
    {
        #ifdef _MSC_VER
        *(vector+i) = start._Val[0] + step*i;
        #else
        *(vector+i) = creal(start) + step*i;
        #endif
    }

    return pIVector;
}

VectorD *VectorD_Copy(VectorD *pIVector, VectorD *pOVector)
{
    int inc1, inc2, size;

    inc1 = pIVector->inc;
    inc2 = pOVector->inc;
    size = (pIVector->size>pOVector->size) ? pOVector->size : pIVector->size;

    dcopy_(&size, pIVector->vector, &inc1, pOVector->vector, &inc2);

    return pOVector;
}

VectorDC *VectorD_Copy2(VectorD *pIVector, VectorDC *pOVector)
{
    register int i;
    int inc1, inc2, size;

    inc1 = pIVector->inc;
    inc2 = pOVector->inc;
    size = (pIVector->size>pOVector->size) ? pOVector->size : pIVector->size;

    for(i=0;i<size;i++)
    {
        #ifdef _MSC_VER
        (pOVector->vector+i*inc2)->_Val[0] = *(pIVector->vector+i*inc1);
        (pOVector->vector+i*inc2)->_Val[1] = 0;
        #else
        *(pOVector->vector+i*inc2) = *(pIVector->vector+i*inc1);
        #endif /*_MSC_VER*/
    }

    return pOVector;
}

double VectorD_Max(VectorD *pVector)
{
    register int i, size;
    register double *vector;
    register double max;

    max = DBL_MIN;
    size = pVector->size;
    vector = pVector->vector;

    for(i=0;i<size;i++)
    {
        max = (*(vector+i) > max) ? *(vector+i) : max;
    }

    return max;
}

#ifdef _MSC_VER
_Dcomplex VectorD_Max2(VectorD *pVector)
#else
double complex VectorD_Max2(VectorD *pVector)
#endif /*_MSC_VER*/
{
    register int i, size;
    register double *vector;
    register double max;
    #ifdef _MSC_VER
    _Dcomplex max2;
    #else
    double complex max2;
    #endif /*_MSC_VER*/

    max = DBL_MIN;
    size = pVector->size;
    vector = pVector->vector;

    for(i=0;i<size;i++)
    {
        max = (*(vector+i) > max) ? *(vector+i) : max;
    }
    #ifdef _MSC_VER
    max2._Val[0] = max;
    max2._Val[1] = 0;
    #else
    max2 = max;
    #endif /*MSC_VER*/

    return max2;
}

double VectorD_Min(VectorD *pVector)
{
    register int i, size;
    register double *vector;
    register double min;

    min = DBL_MAX;
    size = pVector->size;
    vector = pVector->vector;

    for(i=0;i<size;i++)
    {
        min = (*(vector+i) < min) ? *(vector+i) : min;
    }

    return min;
}

#ifdef _MSC_VER
_Dcomplex VectorD_Min2(VectorD *pVector)
#else
double complex VectorD_Min2(VectorD *pVector)
#endif /*_MSC_VER*/
{
    register int i, size;
    register double *vector;
    register double min;
    #ifdef _MSC_VER
    _Dcomplex min2;
    #else
    double complex min2;
    #endif /*_MSC_VER*/

    min = DBL_MAX;
    size = pVector->size;
    vector = pVector->vector;

    for(i=0;i<size;i++)
    {
        min = (*(vector+i) < min) ? *(vector+i) : min;
    }
    #ifdef _MSC_VER
    min2._Val[0] = min;
    min2._Val[1] = 0;
    #else
    min2 = min;
    #endif /*MSC_VER*/

    return min2;
}

VectorD *VectorD_AddVectorD(VectorD *pIVector1, VectorD *pIVector2, VectorD *pOVector)
{
    register int i, size, inc1, inc2, inc3;
    register double *iVector1, *iVector2, *oVector;

    size = pIVector1->size;
    inc1 = pIVector1->inc;
    inc2 = pIVector2->inc;
    inc3 = pOVector->inc;
    iVector1 = pIVector1->vector;
    iVector2 = pIVector2->vector;
    oVector = pOVector->vector;

    if(size != pIVector2->size || pOVector->size != size)
    {
        fprintf(stderr, "Error: pIVector1 and pIVector2 have incorrect size. File: %s Func: %s Line: %d\n", \
            __FILE__, __FUNCTION__, __LINE__);
        return pOVector;
    }

    for(i=0;i<size;++i)
    {
        *(oVector+i*inc3) = *(iVector1+i*inc1) + *(iVector2+i*inc2);
    }

    return pOVector;
}

VectorDC *VectorD_AddVectorD2(VectorD *pIVector1, VectorD *pIVector2, VectorDC *pOVector)
{
    register int i, size, inc1, inc2, inc3;
    register double *iVector1, *iVector2;
    #ifdef _MSC_VER
    register _Dcomplex *oVector;
    #else
    register double complex *oVector;
    #endif /*_MSC_VER*/

    size = pIVector1->size;
    inc1 = pIVector1->inc;
    inc2 = pIVector2->inc;
    inc3 = pOVector->inc;
    iVector1 = pIVector1->vector;
    iVector2 = pIVector2->vector;
    oVector = pOVector->vector;

    if(size != pIVector2->size || pOVector->size != size)
    {
        fprintf(stderr, "Error: pIVector1 and pIVector2 have incorrect size. File: %s Func: %s Line: %d\n", \
            __FILE__, __FUNCTION__, __LINE__);
        return pOVector;
    }

    for(i=0;i<size;++i)
    {
        #ifdef _MSC_VER
        (oVector+i*inc3)->_Val[0] = *(iVector1+i*inc1) + *(iVector2+i*inc2);
        (oVector+i*inc3)->_Val[1] = 0;
        #else
        *(oVector+i*inc3) = *(iVector1+i*inc1) + *(iVector2+i*inc2);
        #endif /*_MSC_VER*/
    }

    return pOVector;
}

VectorDC *VectorD_AddVectorDC(VectorD *pIVector1, VectorDC *pIVector2, VectorDC *pOVector)
{
    register int i, size, inc1, inc2, inc3;
    register double *iVector1;
    #ifdef _MSC_VER
    register _Dcomplex *iVector2, *oVector;
    #else
    register double complex *iVector2, *oVector;
    #endif /*_MSC_VER*/

    size = pIVector1->size;
    inc1 = pIVector1->inc;
    inc2 = pIVector2->inc;
    inc3 = pOVector->inc;
    iVector1 = pIVector1->vector;
    iVector2 = pIVector2->vector;
    oVector = pOVector->vector;

    if(size != pIVector2->size || pOVector->size != size)
    {
        fprintf(stderr, "Error: pIVector1 and pIVector2 have incorrect size. File: %s Func: %s Line: %d\n", \
            __FILE__, __FUNCTION__, __LINE__);
        return pOVector;
    }

    for(i=0;i<size;++i)
    {
        #ifdef _MSC_VER
        (oVector+i*inc3)->_Val[0] = *(iVector1+i*inc1) + (iVector2+i*inc2)->_Val[0];
        (oVector+i*inc3)->_Val[1] = (iVector2+i*inc2)->_Val[1];
        #else
        *(oVector+i*inc3) = *(iVector1+i*inc1) + *(iVector2+i*inc2);
        #endif /*_MSC_VER*/
    }

    return pOVector;
}

VectorD *VectorD_AddVectorDC2(VectorD *pIVector1, VectorDC *pIVector2, VectorD *pOVector)
{
    register int i, size, inc1, inc2, inc3;
    register double *iVector1, *oVector;
    #ifdef _MSC_VER
    register _Dcomplex *iVector2;
    #else
    register double complex *iVector2;
    #endif /*_MSC_VER*/

    size = pIVector1->size;
    inc1 = pIVector1->inc;
    inc2 = pIVector2->inc;
    inc3 = pOVector->inc;
    iVector1 = pIVector1->vector;
    iVector2 = pIVector2->vector;
    oVector = pOVector->vector;

    if(size != pIVector2->size || pOVector->size != size)
    {
        fprintf(stderr, "Error: pIVector1 and pIVector2 have incorrect size. File: %s Func: %s Line: %d\n", \
            __FILE__, __FUNCTION__, __LINE__);
        return pOVector;
    }

    for(i=0;i<size;++i)
    {
        #ifdef _MSC_VER
        *(oVector+i*inc3) = *(iVector1+i*inc1) + (iVector2+i*inc2)->_Val[0];
        #else
        *(oVector+i*inc3) = *(iVector1+i*inc1) + creal(*(iVector2+i*inc2));
        #endif /*_MSC_VER*/
    }

    return pOVector;
}

VectorD *VectorD_SubVectorD(VectorD *pIVector1, VectorD *pIVector2, VectorD *pOVector)
{
    register int i, size, inc1, inc2, inc3;
    register double *iVector1, *iVector2, *oVector;

    size = pIVector1->size;
    inc1 = pIVector1->inc;
    inc2 = pIVector2->inc;
    inc3 = pOVector->inc;
    iVector1 = pIVector1->vector;
    iVector2 = pIVector2->vector;
    oVector = pOVector->vector;

    if(size != pIVector2->size || pOVector->size != size)
    {
        fprintf(stderr, "Error: pIVector1 and pIVector2 have incorrect size. File: %s Func: %s Line: %d\n", \
            __FILE__, __FUNCTION__, __LINE__);

        return pOVector;
    }

    for(i=0;i<size;++i)
    {
        *(oVector+i*inc3) = *(iVector1+i*inc1) - *(iVector2+i*inc2);
    }

    return pOVector;
}

VectorDC *VectorD_SubVectorD2(VectorD *pIVector1, VectorD *pIVector2, VectorDC *pOVector)
{
    register int i, size, inc1, inc2, inc3;
    register double *iVector1, *iVector2;
    #ifdef _MSC_VER
    register _Dcomplex *oVector;
    #else
    register double complex *oVector;
    #endif /*_MSC_VER*/

    size = pIVector1->size;
    inc1 = pIVector1->inc;
    inc2 = pIVector2->inc;
    inc3 = pOVector->inc;
    iVector1 = pIVector1->vector;
    iVector2 = pIVector2->vector;
    oVector = pOVector->vector;

    if(size != pIVector2->size || pOVector->size != size)
    {
        fprintf(stderr, "Error: pIVector1 and pIVector2 have incorrect size. File: %s Func: %s Line: %d\n", \
            __FILE__, __FUNCTION__, __LINE__);

        return pOVector;
    }

    for(i=0;i<size;++i)
    {
        #ifdef _MSC_VER
        (oVector+i*inc3)->_Val[0] = *(iVector1+i*inc1) - *(iVector2+i*inc2);
        (oVector+i*inc3)->_Val[1] = 0;
        #else
        *(oVector+i*inc3) = *(iVector1+i*inc1) - *(iVector2+i*inc2);
        #endif /*_MSC_VER*/
    }

    return pOVector;
}

VectorDC *VectorD_SubVectorDC(VectorD *pIVector1, VectorDC *pIVector2, VectorDC *pOVector)
{
    register int i, size, inc1, inc2, inc3;
    register double *iVector1;
    #ifdef _MSC_VER
    register _Dcomplex *iVector2, *oVector;
    #else
    register double complex *iVector2, *oVector;
    #endif /*_MSC_VER*/

    size = pIVector1->size;
    inc1 = pIVector1->inc;
    inc2 = pIVector2->inc;
    inc3 = pOVector->inc;
    iVector1 = pIVector1->vector;
    iVector2 = pIVector2->vector;
    oVector = pOVector->vector;

    if(size != pIVector2->size || pOVector->size != size)
    {
        fprintf(stderr, "Error: pIVector1 and pIVector2 have incorrect size. File: %s Func: %s Line: %d\n", \
            __FILE__, __FUNCTION__, __LINE__);

        return pOVector;
    }

    for(i=0;i<size;++i)
    {
        #ifdef _MSC_VER
        (oVector+i*inc3)->_Val[0] = *(iVector1+i*inc1) - (iVector2+i*inc2)->_Val[0];
        (oVector+i*inc3)->_Val[1] = -(iVector2+i*inc2)->_Val[1];
        #else
        *(oVector+i*inc3) = *(iVector1+i*inc1) - *(iVector2+i*inc2);
        #endif /*_MSC_VER*/
    }

    return pOVector;
}

VectorD *VectorD_SubVectorDC2(VectorD *pIVector1, VectorDC *pIVector2, VectorD *pOVector)
{
    register int i, size, inc1, inc2, inc3;
    register double *iVector1;
    #ifdef _MSC_VER
    register _Dcomplex *iVector2;
    #else
    register double complex *iVector2;
    #endif /*_MSC_VER*/
    register double *oVector;

    size = pIVector1->size;
    inc1 = pIVector1->inc;
    inc2 = pIVector2->inc;
    inc3 = pOVector->inc;
    iVector1 = pIVector1->vector;
    iVector2 = pIVector2->vector;
    oVector = pOVector->vector;

    if(size != pIVector2->size || pOVector->size != size)
    {
        fprintf(stderr, "Error: pIVector1 and pIVector2 have incorrect size. File: %s Func: %s Line: %d\n", \
            __FILE__, __FUNCTION__, __LINE__);

        return pOVector;
    }

    for(i=0;i<size;++i)
    {
        #ifdef _MSC_VER
        *(oVector+i*inc3) = *(iVector1+i*inc1) - (iVector2+i*inc2)->_Val[0];
        #else
        *(oVector+i*inc3) = *(iVector1+i*inc1) - creal(*(iVector2+i*inc2));
        #endif /*_MSC_VER*/
    }

    return pOVector;
}

VectorD *VectorD_MulNumD(VectorD *pIVector, double num, VectorD *pOVector)
{
    int size, inc1, inc2;

    size = pIVector->size;
    inc1 = pIVector->inc;
    inc2 = pOVector->inc;

    dcopy_(&size, pIVector->vector, &inc1, pOVector->vector, &inc2);
    dscal_(&size, &num, pOVector->vector, &inc2);

    return pOVector;
}

VectorDC *VectorD_MulNumD2(VectorD *pIVector, double num, VectorDC *pOVector)
{
    register int i;
    int size, inc1, inc2;
    double inv_num;

    size = pIVector->size;
    inc1 = pIVector->inc;
    inc2 = pOVector->inc;
    inv_num = 1 / num;

    dscal_(&size, &num, pIVector->vector, &inc1);
    for(i=0;i<size;i++)
    {
        #ifdef _MSC_VER
        (pOVector->vector+i*inc2)->_Val[0] = *(pIVector->vector+i*inc1);
        (pOVector->vector+i*inc2)->_Val[1] = 0;
        #else
        *(pOVector->vector+i*inc2) = *(pIVector->vector+i*inc1);
        #endif /*_MSC_VER*/
    }
    dscal_(&size, &inv_num, pIVector->vector, &inc1);

    return pOVector;
}

#ifdef _MSC_VER
VectorDC *VectorD_MulNumDC(VectorD *pIVector, _Dcomplex num, VectorDC *pOVector)
#else
VectorDC *VectorD_MulNumDC(VectorD *pIVector, double complex num, VectorDC *pOVector)
#endif /*_MSC_VER*/
{
    register int i;
    int size, inc1, inc2;

    size = pIVector->size;
    inc1 = pIVector->inc;
    inc2 = pOVector->inc;

    for(i=0;i<size;i++)
    {
        #ifdef _MSC_VER
        (pOVector->vector+i*inc2)->_Val[0] = *(pIVector->vector+i*inc1);
        (pOVector->vector+i*inc2)->_Val[1] = 0;
        #else
        *(pOVector->vector+i*inc2) = *(pIVector->vector+i*inc1);
        #endif /*_MSC_VER*/
    }
    zscal_(&size, &num, pOVector->vector, &inc2);

    return pOVector;
}

#ifdef _MSC_VER
VectorD *VectorD_MulNumDC2(VectorD *pIVector, _Dcomplex num, VectorD *pOVector)
#else
VectorD *VectorD_MulNumDC2(VectorD *pIVector, double complex num, VectorD *pOVector)
#endif /*_MSC_VER*/
{
    register int i;
    int size, inc1, inc2;
    double tmp;

    size = pIVector->size;
    inc1 = pIVector->inc;
    inc2 = pOVector->inc;
    #ifdef _MSC_VER
    tmp = num._Val[0];
    #else
    tmp = creal(num);
    #endif /*_MSC_VER*/

    for(i=0;i<size;i++)
    {
        *(pOVector->vector+i*inc2) = *(pIVector->vector+i*inc1);
    }
    dscal_(&size, &tmp, pOVector->vector, &inc2);

    return pOVector;
}

double VectorD_MulVectorD(VectorD *pIVector1, VectorD *pIVector2)
{
    int size, inc1, inc2;
    register double result;

    size = pIVector1->size;
    inc1 = pIVector1->inc;
    inc2 = pIVector2->inc;

    result = ddot_(&size, pIVector1->vector, &inc1, pIVector2->vector, &inc2);

    return result;
}

#ifdef _MSC_VER
_Dcomplex VectorD_MulVectorDC(VectorD *pIVector1, VectorDC *pIVector2)
#else
double complex VectorD_MulVectorDC(VectorD *pIVector1, VectorDC *pIVector2)
#endif /*_MSC_VER*/
{
    register int i;
    int size, inc1, inc2;
    #ifdef _MSC_VER
    register _Dcomplex result;
    #else
    register double complex result;
    #endif /*_MSC_VER*/
    register double *vector1;
    #ifdef _MSC_VER
    register _Dcomplex *vector2;
    #else
    register double complex *vector2;
    #endif /*_MSC_VER*/

    size = pIVector1->size;
    inc1 = pIVector1->inc;
    inc2 = pIVector2->inc;
    #ifdef _MSC_VER
    result._Val[0] = 0;
    result._Val[1] = 0;
    #else
    result = 0;
    #endif /*_MSC_VER*/
    vector1 = pIVector1->vector;
    vector2 = pIVector2->vector;

    for(i=0;i<size;i++)
    {
        #ifdef _MSC_VER
        result._Val[0] += *(vector1+i*inc1) * (vector2+i*inc2)->_Val[0];
        result._Val[1] += *(vector1+i*inc1) * (vector2+i*inc2)->_Val[1];
        #else
        result += *(vector1+i*inc1) * *(vector2+i*inc2);
        #endif /*_MSC_VER*/
    }

    return result;
}

double VectorD_MulVectorDC2(VectorD *pIVector1, VectorDC *pIVector2)
{
    register int i;
    int size, inc1, inc2;
    register double result;
    register double *vector1;
    #ifdef _MSC_VER
    register _Dcomplex *vector2;
    #else
    register double complex *vector2;
    #endif /*_MSC_VER*/

    size = pIVector1->size;
    inc1 = pIVector1->inc;
    inc2 = pIVector2->inc;
    result = 0;
    vector1 = pIVector1->vector;
    vector2 = pIVector2->vector;

    for(i=0;i<size;i++)
    {
        #ifdef _MSC_VER
        result += *(vector1+i*inc1) * (vector2+i*inc2)->_Val[0];
        #else
        result += *(vector1+i*inc1) * creal(*(vector2+i*inc2));
        #endif /*_MSC_VER*/
    }

    return result;
}

VectorD *VectorD_EMulVectorD(VectorD *pIVector1, VectorD *pIVector2, VectorD *pOVector)
{
    int size, inc1, inc2, inc3;
    register int i;

    size = pIVector1->size;
    inc1 = pIVector1->inc;
    inc2 = pIVector2->inc;
    inc3 = pOVector->inc;

    for(i=0;i<size;i++)
    {
        *(pOVector->vector+i*inc3) = *(pIVector1->vector+i*inc1) * *(pIVector2->vector+i*inc2);
    }

    return pOVector;
}

VectorDC *VectorD_EMulVectorD2(VectorD *pIVector1, VectorD *pIVector2, VectorDC *pOVector)
{
    int size, inc1, inc2, inc3;
    register int i;

    size = pIVector1->size;
    inc1 = pIVector1->inc;
    inc2 = pIVector2->inc;
    inc3 = pOVector->inc;

    for(i=0;i<size;i++)
    {
        #ifdef _MSC_VER
        (pOVector->vector+i*inc3)->_Val[0] = *(pIVector1->vector+i*inc1) * *(pIVector2->vector+i*inc2);
        (pOVector->vector+i*inc3)->_Val[1] = 0;
        #else
        *(pOVector->vector+i*inc3) = *(pIVector1->vector+i*inc1) * *(pIVector2->vector+i*inc2);
        #endif /*_MSC_VER*/
    }

    return pOVector;
}

VectorDC *VectorD_EMulVectorDC(VectorD *pIVector1, VectorDC *pIVector2, VectorDC *pOVector)
{
    int size, inc1, inc2, inc3;
    register int i;

    size = pIVector1->size;
    inc1 = pIVector1->inc;
    inc2 = pIVector2->inc;
    inc3 = pOVector->inc;

    for(i=0;i<size;i++)
    {
        #ifdef _MSC_VER
        (pOVector->vector+i*inc3)->_Val[0] = *(pIVector1->vector+i*inc1) * (pIVector2->vector+i*inc2)->_Val[0];
        (pOVector->vector+i*inc3)->_Val[1] = *(pIVector1->vector+i*inc1) * (pIVector2->vector+i*inc2)->_Val[1];
        #else
        *(pOVector->vector+i*inc3) = *(pIVector1->vector+i*inc1) * *(pIVector2->vector+i*inc2);
        #endif /*_MSC_VER*/
    }

    return pOVector;
}

VectorD *VectorD_EMulVectorDC2(VectorD *pIVector1, VectorDC *pIVector2, VectorD *pOVector)
{
    int size, inc1, inc2, inc3;
    register int i;

    size = pIVector1->size;
    inc1 = pIVector1->inc;
    inc2 = pIVector2->inc;
    inc3 = pOVector->inc;

    for(i=0;i<size;i++)
    {
        #ifdef _MSC_VER
        *(pOVector->vector+i*inc3) = *(pIVector1->vector+i*inc1) * (pIVector2->vector+i*inc2)->_Val[0];
        #else
        *(pOVector->vector+i*inc3) = *(pIVector1->vector+i*inc1) * creal(*(pIVector2->vector+i*inc2));
        #endif /*_MSC_VER*/
    }

    return pOVector;
}

double VectorD_CMulVectorD(VectorD *pIVector1, VectorD *pIVector2)
{
    register int i;
    int size, inc1, inc2;
    register double result;
    register double *vector1;
    register double *vector2;

    size = pIVector1->size;
    inc1 = pIVector1->inc;
    inc2 = pIVector2->inc;
    result = 0;
    vector1 = pIVector1->vector;
    vector2 = pIVector2->vector;

    for(i=0;i<size;i++)
    {
        result += *(vector1+i*inc1) * *(vector2+i*inc2);
    }

    return result;
}

#ifdef _MSC_VER
_Dcomplex VectorD_CMulVectorD2(VectorD *pIVector1, VectorD *pIVector2)
#else
double complex VectorD_CMulVectorD2(VectorD *pIVector1, VectorD *pIVector2)
#endif /*_MSC_VER*/
{
    register int i;
    int size, inc1, inc2;
    #ifdef _MSC_VER
    register _Dcomplex result;
    #else
    register double complex result;
    #endif /*_MSC_VER*/
    register double *vector1;
    register double *vector2;

    size = pIVector1->size;
    inc1 = pIVector1->inc;
    inc2 = pIVector2->inc;
    #ifdef _MSC_VER
    result._Val[0] = 0;
    result._Val[1] = 0;
    #else
    result = 0;
    #endif /*_MSC_VER*/
    vector1 = pIVector1->vector;
    vector2 = pIVector2->vector;

    for(i=0;i<size;i++)
    {
        #ifdef _MSC_VER
        result._Val[0] += *(vector1+i*inc1) * *(vector2+i*inc2);
        #else
        result += *(vector1+i*inc1) * *(vector2+i*inc2);
        #endif /*_MSC_VER*/
    }

    return result;
}

#ifdef _MSC_VER
_Dcomplex VectorD_CMulVectorDC(VectorD *pIVector1, VectorDC *pIVector2)
#else
double complex VectorD_CMulVectorDC(VectorD *pIVector1, VectorDC *pIVector2)
#endif /*_MSC_VER*/
{
    register int i;
    int size, inc1, inc2;
    #ifdef _MSC_VER
    register _Dcomplex result;
    #else
    register double complex result;
    #endif /*_MSC_VER*/
    register double *vector1;
    #ifdef _MSC_VER
    register _Dcomplex *vector2;
    #else
    register double complex *vector2;
    #endif /*_MSC_VER*/

    size = pIVector1->size;
    inc1 = pIVector1->inc;
    inc2 = pIVector2->inc;
    #ifdef _MSC_VER
    result._Val[0] = 0;
    result._Val[1] = 0;
    #else
    result = 0;
    #endif /*_MSC_VER*/
    vector1 = pIVector1->vector;
    vector2 = pIVector2->vector;

    for(i=0;i<size;i++)
    {
        #ifdef _MSC_VER
        result._Val[0] += *(vector1+i*inc1) * (vector2+i*inc2)->_Val[0];
        result._Val[1] += *(vector1+i*inc1) * (vector2+i*inc2)->_Val[1];
        #else
        result += *(vector1+i*inc1) * *(vector2+i*inc2);
        #endif /*_MSC_VER*/
    }

    return result;
}

double VectorD_CMulVectorDC2(VectorD *pIVector1, VectorDC *pIVector2)
{
    register int i;
    int size, inc1, inc2;
    register double result;
    register double *vector1;
    #ifdef _MSC_VER
    register _Dcomplex *vector2;
    #else
    register double complex *vector2;
    #endif /*_MSC_VER*/

    size = pIVector1->size;
    inc1 = pIVector1->inc;
    inc2 = pIVector2->inc;
    result = 0;
    vector1 = pIVector1->vector;
    vector2 = pIVector2->vector;

    for(i=0;i<size;i++)
    {
        #ifdef _MSC_VER
        result += *(vector1+i*inc1) * (vector2+i*inc2)->_Val[0];
        #else
        result += *(vector1+i*inc1) * creal(*(vector2+i*inc2));
        #endif /*_MSC_VER*/
    }

    return result;
}

double VectorD_Norm(VectorD *pIVector)
{
    int size, inc;
    double result;

    inc = pIVector->inc;
    size = pIVector->size;

    result = dnrm2_(&size, pIVector->vector, &inc);

    return result;
}

#ifdef _MSC_VER
_Dcomplex VectorD_Norm2(VectorD *pIVector)
#else
double complex VectorD_Norm2(VectorD *pIVector)
#endif /*_MSC_VER*/
{
    int size, inc;
    #ifdef _MSC_VER
    _Dcomplex result;
    #else
    double complex result;
    #endif /*_MSC_VER*/

    inc = pIVector->inc;
    size = pIVector->size;

    #ifdef _MSC_VER
    result._Val[0] = dnrm2_(&size, pIVector->vector, &inc);
    result._Val[1] = 0;
    #else
    result = dnrm2_(&size, pIVector->vector, &inc);
    #endif /*_MSC_VER*/

    return result;
}

double VectorD_NormSq(VectorD *pIVector)
{
    register int i, size, inc;
    register double normSq;
    register double *iVector;

    size = pIVector->size;
    inc = pIVector->inc;
    iVector = pIVector->vector;
    normSq = 0.0;

    for(i=size; i-- > 0;)
    {
        normSq += *(iVector+i*inc) * *(iVector+i*inc);
    }

    return normSq;
}

#ifdef _MSC_VER
_Dcomplex VectorD_NormSq2(VectorD *pIVector)
#else
double complex VectorD_NormSq2(VectorD *pIVector)
#endif /*_MSC_VER*/
{
    register int i, size, inc;
    register double *iVector;
    #ifdef _MSC_VER
    register _Dcomplex normSq;
    #else
    register double complex normSq;
    #endif /*_MSC_VER*/

    size = pIVector->size;
    inc = pIVector->inc;
    iVector = pIVector->vector;
    #ifdef _MSC_VER
    normSq._Val[0] = 0;
    normSq._Val[1] = 0;
    #else
    normSq = 0.0;
    #endif /*_MSC_VER*/

    for(i=size; i-- > 0;)
    {
        #ifdef _MSC_VER
        normSq._Val[0] +=  *(iVector+i*inc) * *(iVector+i*inc);
        #else
        normSq += *(iVector+i*inc) * *(iVector+i*inc);
        #endif /*_MSC_VER*/
    }

    return normSq;
}
