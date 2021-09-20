#include "Vector.h"

extern Object *_VectorDC_Del(Object *pIVector);
extern Object *_VectorDC_UnWrap(Object *pIVector);
extern void _VectorDC_Show(Tensor *pIVector);
extern Vector *_VectorDC_Full(Vector *pIVector, Num *pNum);
extern Vector *_VectorDC_Linspace(Vector *pIVector, Num *start, Num *end);
extern Vector *_VectorDC_Copy(Vector *pIVector, Vector *pOVector);
extern Num *_VectorDC_Max(Vector *pVector, Num *pONum);
extern Num *_VectorDC_Min(Vector *pVector, Num *pONum);
extern Vector *_VectorDC_AddVector(Vector *pVector1, Vector *pVector2, Vector *pOVector);
extern Vector *_VectorDC_SubVector(Vector *pVector1, Vector *pVector2, Vector *pOVector);
extern Vector *_VectorDC_MulNum(Vector *pIVector, Num *pNum, Vector *pOVector);
extern Num *_VectorDC_MulVector(Vector *pVector1, Vector *pVector2, Num *pONum);
extern Vector *_VectorDC_EMulVector(Vector *pVector1, Vector *pVector2, Vector *pOVector);
extern Num *_VectorDC_CMulVector(Vector *pVector1, Vector *pVector2, Num *pONum);
extern Num *_VectorDC_Norm(Vector *pIVector, Num *pONum);
extern Num *_VectorDC_NormSq(Vector *pIVector, Num *pONum);

VectorDC *VectorDC_New(int size)
{
    VectorDC *pNewVectorDC;

    pNewVectorDC = (VectorDC *)malloc(sizeof(VectorDC));
    if(pNewVectorDC == NULL)
    {
        fprintf(stderr, "Error: Can't create a new vectorDC. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
        return NULL;
    }
    #ifdef _MSC_VER
    pNewVectorDC->vector = (_Dcomplex *)malloc(sizeof(_Dcomplex)*size);
    #else
    pNewVectorDC->vector = (double complex *)malloc(sizeof(double complex)*size);
    #endif /*_MSC_VER*/
    if(pNewVectorDC->vector == NULL)
    {
        free(pNewVectorDC);
        fprintf(stderr, "Error: Can't create a new vectorDC. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
        return NULL;
    }
    pNewVectorDC->size = size;
    pNewVectorDC->inc = 1;
    ((Object *)pNewVectorDC)->type = VECTORDC;
    ((Object *)pNewVectorDC)->Del = _VectorDC_Del;
    ((Object *)pNewVectorDC)->UnWrap = _VectorDC_UnWrap;
    ((Tensor *)pNewVectorDC)->Show = _VectorDC_Show;
    pNewVectorDC->parent.Full = _VectorDC_Full;
    pNewVectorDC->parent.Linspace = _VectorDC_Linspace;
    pNewVectorDC->parent.Copy = _VectorDC_Copy;
    pNewVectorDC->parent.Max = _VectorDC_Max;
    pNewVectorDC->parent.Min = _VectorDC_Min;
    pNewVectorDC->parent.AddVector = _VectorDC_AddVector;
    pNewVectorDC->parent.SubVector = _VectorDC_SubVector;
    pNewVectorDC->parent.MulNum = _VectorDC_MulNum;
    pNewVectorDC->parent.EMulVector = _VectorDC_EMulVector;
    pNewVectorDC->parent.MulVector = _VectorDC_MulVector;
    pNewVectorDC->parent.CMulVector = _VectorDC_CMulVector;
    pNewVectorDC->parent.Norm = _VectorDC_Norm;
    pNewVectorDC->parent.NormSq = _VectorDC_NormSq;

    return pNewVectorDC;
}

VectorDC *VectorDC_Del(VectorDC *pIVector)
{
    if(pIVector != NULL)
    {
        free(pIVector->vector);
        free(pIVector);
    }

    return NULL;
}

#ifdef _MSC_VER
VectorDC *VectorDC_Wrap(_Dcomplex *work, int size, int inc)
#else
VectorDC *VectorDC_Wrap(double complex *work, int size, int inc)
#endif /*_MSC_VER*/
{
    if(work == NULL)
    {
        return NULL;
    }

    VectorDC *pNewVectorDC;

    pNewVectorDC = (VectorDC *)malloc(sizeof(VectorDC));
    if(pNewVectorDC == NULL)
    {
        fprintf(stderr, "Error: Can't create a new vectorDC. File: %s Func: %s Line: %d\n", __FILE__, __FUNCTION__, \
            __LINE__);
        return NULL;
    }
    pNewVectorDC->vector = work;
    pNewVectorDC->size = size;
    pNewVectorDC->inc = inc;
    ((Object *)pNewVectorDC)->type = VECTORDC;
    ((Object *)pNewVectorDC)->Del = _VectorDC_Del;
    ((Object *)pNewVectorDC)->UnWrap = _VectorDC_UnWrap;
    ((Tensor *)pNewVectorDC)->Show = _VectorDC_Show;
    pNewVectorDC->parent.Full = _VectorDC_Full;
    pNewVectorDC->parent.Linspace = _VectorDC_Linspace;
    pNewVectorDC->parent.Copy = _VectorDC_Copy;
    pNewVectorDC->parent.Max = _VectorDC_Max;
    pNewVectorDC->parent.Min = _VectorDC_Min;
    pNewVectorDC->parent.AddVector = _VectorDC_AddVector;
    pNewVectorDC->parent.SubVector = _VectorDC_SubVector;
    pNewVectorDC->parent.MulNum = _VectorDC_MulNum;
    pNewVectorDC->parent.MulVector = _VectorDC_MulVector;
    pNewVectorDC->parent.EMulVector = _VectorDC_EMulVector;
    pNewVectorDC->parent.CMulVector = _VectorDC_CMulVector;
    pNewVectorDC->parent.Norm = _VectorDC_Norm;
    pNewVectorDC->parent.NormSq = _VectorDC_NormSq;

    return pNewVectorDC;
}

VectorDC *VectorDC_UnWrap(VectorDC *pIVector)
{
    if (pIVector != NULL)
    {
        free(pIVector);
    }

    return NULL;
}

void VectorDC_Show(VectorDC *pIVector)
{
    register int i, size, inc;

    size = pIVector->size;
    inc = pIVector->inc;
    printf("VectorDC_Show: \n vector:[");
    for(i=0;i<size;i++)
    {
        printf("%f+%f*j, ", creal(*(pIVector->vector+i*inc)), cimag(*(pIVector->vector+i*inc)));
    }
    printf("] \n size:%d \n inc:%d \n\n", pIVector->size, pIVector->inc);
}

#ifdef _MSC_VER
VectorDC *VectorDC_Full(VectorDC *pIVector, _Dcomplex num)
#else
VectorDC *VectorDC_Full(VectorDC *pIVector, double complex num)
#endif /*_MSC_VER*/
{
    register int i;
    int size;
    #ifdef _MSC_VER
    register _Dcomplex *vector;
    #else
    register double complex *vector;
    #endif /*_MSC_VER*/

    size = pIVector->size;
    vector = pIVector->vector;

    for(i=0;i<size;i++)
    {
        *(vector+i) = num;
    }

    return pIVector;
}

VectorDC *VectorDC_Full2(VectorDC *pIVector, double num)
{
    register int i;
    int size;
    #ifdef _MSC_VER
    register _Dcomplex *vector;
    #else
    register double complex *vector;
    #endif /*_MSC_VER*/

    size = pIVector->size;
    vector = pIVector->vector;

    for(i=0;i<size;i++)
    {
        #ifdef _MSC_VER
        (vector+i)->_Val[0] = num;
        (vector+i)->_Val[0] = 0;
        #else
        *(vector+i) = num;
        #endif /*_MSC_VER*/
    }

    return pIVector;
}

#ifdef _MSC_VER
VectorDC *VectorDC_Linspace(VectorDC *pIVector, _Dcomplex start, _Dcomplex end)
#else
VectorDC *VectorDC_Linspace(VectorDC *pIVector, double complex start, double complex end)
#endif /*_MSC_VER*/
{
    register int i;
    int size;
    #ifdef _MSC_VER
    register _Dcomplex step;
    register _Dcomplex *vector;
    #else
    register double complex step;
    register double complex *vector;
    #endif /*_MSC_VER*/

    size = pIVector->size;
    vector = pIVector->vector;
    #ifdef _MSC_VER
    step._Val[0] = (end._Val[0] - start._Val[0]) / ((double)size-1.);
    step._Val[1] = (end._Val[1] - start._Val[1]) / ((double)size-1.);
    #else
    step = (end-start) / (size-1);
    #endif /*_MSC_VER*/

    for(i=0;i<size;i++)
    {
        #ifdef _MSC_VER
        (vector+i)->_Val[0] = start._Val[0] + step._Val[0]*i;
        (vector+i)->_Val[1] = start._Val[1] + step._Val[1]*i;
        #else
        *(vector+i) = start + step*i;
        #endif /*_MSC_VER*/
    }

    return pIVector;
}

VectorDC *VectorDC_Linspace1(VectorDC *pIVector, double start, double end)
{
    register int i;
    int size;
    double step;
    #ifdef _MSC_VER
    register _Dcomplex *vector;
    #else
    register double complex *vector;
    #endif /*_MSC_VER*/

    size = pIVector->size;
    vector = pIVector->vector;
    step = (end - start) / ((double)size-1.);

    for(i=0;i<size;i++)
    {
        #ifdef _MSC_VER
        (vector+i)->_Val[0] = start + step*i;
        (vector+i)->_Val[1] = 0;
        #else
        *(vector+i) = start + step*i;
        #endif /*_MSC_VER*/
    }

    return pIVector;
}

#ifdef _MSC_VER
VectorDC *VectorDC_Linspace2(VectorDC *pIVector, double start, _Dcomplex end)
#else
VectorDC *VectorDC_Linspace2(VectorDC *pIVector, double start, double complex end)
#endif /*_MSC_VER*/
{
    register int i;
    int size;
    #ifdef _MSC_VER
    register _Dcomplex step;
    register _Dcomplex *vector;
    #else
    register double complex step;
    register double complex *vector;
    #endif /*_MSC_VER*/

    size = pIVector->size;
    vector = pIVector->vector;
    #ifdef _MSC_VER
    step._Val[0] = (end._Val[0] - start) / ((double)size-1.);
    step._Val[1] = (end._Val[1] - 0) / ((double)size-1.);
    #else
    step = (end - start) / (size-1);
    #endif /*_MSC_VER*/

    for(i=0;i<size;i++)
    {
        #ifdef _MSC_VER
        (vector+i)->_Val[0] = start + step._Val[0]*i;
        (vector+i)->_Val[1] = 0 + step._Val[1]*i;
        #else
        *(vector+i) = start + step*i;
        #endif /*_MSC_VER*/
    }

    return pIVector;
}

#ifdef _MSC_VER
VectorDC *VectorDC_Linspace3(VectorDC *pIVector, _Dcomplex start, double end)
#else
VectorDC *VectorDC_Linspace3(VectorDC *pIVector, double complex start, double end)
#endif /*_MSC_VER*/
{
    register int i;
    int size;
    #ifdef _MSC_VER
    register _Dcomplex step;
    register _Dcomplex *vector;
    #else
    register double complex step;
    register double complex *vector;
    #endif /*_MSC_VER*/

    size = pIVector->size;
    vector = pIVector->vector;
    #ifdef _MSC_VER
    step._Val[0] = (end - start._Val[0]) / ((double)size-1);
    step._Val[1] = (0 - start._Val[1]) / ((double)size-1);
    #else
    step = (end - start) / (size-1);
    #endif /*_MSC_VER*/

    for(i=0;i<size;i++)
    {
        #ifdef _MSC_VER
        (vector+i)->_Val[0] = start._Val[0] + step._Val[0]*i;
        (vector+i)->_Val[1] = start._Val[1] + step._Val[1]*i;
        #else
        *(vector+i) = start + step*i;
        #endif /*_MSC_VER*/
    }

    return pIVector;
}

VectorDC *VectorDC_Copy(VectorDC *pIVector, VectorDC *pOVector)
{
    int inc1, inc2, size;

    inc1 = pIVector->inc;
    inc2 = pOVector->inc;
    size = (pIVector->size>pOVector->size) ? pOVector->size : pIVector->size;

    zcopy_(&size, pIVector->vector, &inc1, pOVector->vector, &inc2);

    return pOVector;
}

VectorD *VectorDC_Copy2(VectorDC *pIVector, VectorD *pOVector)
{
    register int i;
    int inc1, inc2, size;

    inc1 = pIVector->inc;
    inc2 = pOVector->inc;
    size = (pIVector->size>pOVector->size) ? pOVector->size : pIVector->size;

    for(i=0;i<size;i++)
    {
        #ifdef _MSC_VER
        *(pOVector->vector+i*inc2) = (pIVector->vector+i*inc1)->_Val[0];
        #else
        *(pOVector->vector+i*inc2) = creal(*(pIVector->vector+i*inc1));
        #endif /*_MSC_VER*/
    }

    return pOVector;
}

#ifdef _MSC_VER
_Dcomplex VectorDC_Max(VectorDC *pVector)
#else
double complex VectorDC_Max(VectorDC *pVector)
#endif /*_MSC_VER*/
{
    register int i, size;
    #ifdef _MSC_VER
    register _Dcomplex *vector;
    register _Dcomplex max_value;
    #else
    register double complex *vector;
    register double complex max_value;
    #endif /*_MSC_VER*/
    register double max, tmp;

    max = DBL_MIN;
    #ifdef _MSC_VER
    max_value = _Cbuild(0, 0);
    #else
    max_value = 0;
    #endif /*_MSC_VER*/
    size = pVector->size;
    vector = pVector->vector;

    for(i=0;i<size;i++)
    {
        tmp = cabs(*(vector+i));
        max = (tmp > max) ? tmp : max;
        max_value = *(vector+i);
    }

    return max_value;
}

double VectorDC_Max2(VectorDC *pVector)
{
    register int i, size;
    #ifdef _MSC_VER
    register _Dcomplex *vector;
    register _Dcomplex max_value;
    #else
    register double complex *vector;
    register double complex max_value;
    #endif /*_MSC_VER*/
    register double max, tmp;

    max = DBL_MIN;
    #ifdef _MSC_VER
    max_value = _Cbuild(0, 0);
    #else
    max_value = 0;
    #endif
    size = pVector->size;
    vector = pVector->vector;

    for(i=0;i<size;i++)
    {
        tmp = cabs(*(vector+i));
        max = (tmp > max) ? tmp : max;
        max_value = *(vector+i);
    }
    #ifdef _MSC_VER
    max = max_value._Val[0];
    #else
    max = creal(max_value);
    #endif /*_MSC_VER*/

    return max;
}

#ifdef _MSC_VER
_Dcomplex VectorDC_Min(VectorDC *pVector)
#else
double complex VectorDC_Min(VectorDC *pVector)
#endif /*_MSC_VER*/
{
    register int i, size;
    #ifdef _MSC_VER
    register _Dcomplex *vector;
    register _Dcomplex min_value;
    #else
    register double complex *vector;
    register double complex min_value;
    #endif /*_MSC_VER*/
    register double min, tmp;

    min = DBL_MAX;
    #ifdef _MSC_VER
    min_value = _Cbuild(0, 0);
    #else
    min_value = 0;
    #endif /*_MSC_VER*/
    size = pVector->size;
    vector = pVector->vector;

    for(i=0;i<size;i++)
    {
        tmp = cabs(*(vector+i));
        min = (tmp < min) ? tmp : min;
        min_value = *(vector+i);
    }

    return min_value;
}

double VectorDC_Min2(VectorDC *pVector)
{
    register int i, size;
    #ifdef _MSC_VER
    register _Dcomplex *vector;
    register _Dcomplex min_value;
    #else
    register double complex *vector;
    register double complex min_value;
    #endif /*_MSC_VER*/
    register double min, tmp;

    min = DBL_MAX;
    #ifdef _MSC_VER
    min_value = _Cbuild(0,0);
    #else
    min_value = 0;
    #endif /*_MSC_VER*/
    size = pVector->size;
    vector = pVector->vector;

    for(i=0;i<size;i++)
    {
        tmp = cabs(*(vector+i));
        min = (tmp < min) ? tmp : min;
        min_value = *(vector+i);
    }
    #ifdef MSC_VER
    min = min_value._Val[0];
    #else
    min = creal(min_value);
    #endif /*_MSC_VER*/
    
    return min;
}

VectorDC *VectorDC_AddVectorD(VectorDC *pIVector1, VectorD *pIVector2, VectorDC *pOVector)
{
    register int i, size, inc1, inc2, inc3;
    #ifdef _MSC_VER
    register _Dcomplex *iVector1, *oVector;
    #else
    register double complex *iVector1, *oVector;
    #endif /*_MSC_VER*/
    register double *iVector2;

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
        (oVector+i*inc3)->_Val[0] = (iVector1+i*inc1)->_Val[0] + *(iVector2+i*inc2);
        (oVector+i*inc3)->_Val[1] = (iVector1+i*inc1)->_Val[1];
        #else
        *(oVector+i*inc3) = *(iVector1+i*inc1) + *(iVector2+i*inc2);
        #endif /*_MSC_VER*/
    }

    return pOVector;
}

VectorD *VectorDC_AddVectorD2(VectorDC *pIVector1, VectorD *pIVector2, VectorD *pOVector)
{
    register int i, size, inc1, inc2, inc3;
    #ifdef _MSC_VER
    register _Dcomplex *iVector1;
    #else
    register double complex *iVector1;
    #endif /*_MSC_VER*/
    register double *iVector2, *oVector;

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
        *(oVector+i*inc3) = (iVector1+i*inc1)->_Val[0] + *(iVector2+i*inc2);
        #else
        *(oVector+i*inc3) = creal(*(iVector1+i*inc1)) + *(iVector2+i*inc2);
        #endif /*_MSC_VER*/
    }

    return pOVector;
}

VectorDC *VectorDC_AddVectorDC(VectorDC *pIVector1, VectorDC *pIVector2, VectorDC *pOVector)
{
    register int i, size, inc1, inc2, inc3;
    #ifdef _MSC_VER
    register _Dcomplex *iVector1, *iVector2, *oVector;
    #else
    register double complex *iVector1, *iVector2, *oVector;
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
        *(oVector+i*inc3) = _Caddcc(*(iVector1+i*inc1), *(iVector2+i*inc2));
        #else
        *(oVector+i*inc3) = *(iVector1+i*inc1) + *(iVector2+i*inc2);
        #endif /*_MSC_VER*/
    }

    return pOVector;
}

VectorD *VectorDC_AddVectorDC2(VectorDC *pIVector1, VectorDC *pIVector2, VectorD *pOVector)
{
    register int i, size, inc1, inc2, inc3;
    #ifdef _MSC_VER
    register _Dcomplex *iVector1, *iVector2;
    #else
    register double complex *iVector1, *iVector2;
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
        *(oVector+i*inc3) = (iVector1+i*inc1)->_Val[0] + (iVector2+i*inc2)->_Val[0];
        #else
        *(oVector+i*inc3) = creal(*(iVector1+i*inc1) + *(iVector2+i*inc2));
        #endif /*_MSC_VER*/
    }

    return pOVector;
}

VectorDC *VectorDC_SubVectorD(VectorDC *pIVector1, VectorD *pIVector2, VectorDC *pOVector)
{
    register int i, size, inc1, inc2, inc3;
    double *iVector2;
    #ifdef _MSC_VER
    register _Dcomplex *iVector1, *oVector;
    #else
    register double complex *iVector1, *oVector;
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
        (oVector+i*inc3)->_Val[0] = (iVector1+i*inc1)->_Val[0] - *(iVector2+i*inc2);
        (oVector+i*inc3)->_Val[1] = (iVector1+i*inc1)->_Val[1];
        #else
        *(oVector+i*inc3) = *(iVector1+i*inc1) - *(iVector2+i*inc2);
        #endif /*_MSC_VER*/
    }

    return pOVector;
}

VectorD *VectorDC_SubVectorD2(VectorDC *pIVector1, VectorD *pIVector2, VectorD *pOVector)
{
    register int i, size, inc1, inc2, inc3;
    #ifdef _MSC_VER
    register _Dcomplex *iVector1;
    #else
    register double complex *iVector1;
    #endif /*_MSC_VER*/
    register double *iVector2, *oVector;

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
        *(oVector+i*inc3) = (iVector1+i*inc1)->_Val[0] - *(iVector2+i*inc2);
        #else
        *(oVector+i*inc3) = creal(*(iVector1+i*inc1)) - *(iVector2+i*inc2);
        #endif /*_MSC_VER*/
    }

    return pOVector;
}

VectorDC *VectorDC_SubVectorDC(VectorDC *pIVector1, VectorDC *pIVector2, VectorDC *pOVector)
{
    register int i, size, inc1, inc2, inc3;
    #ifdef _MSC_VER
    register _Dcomplex *iVector1, *iVector2, *oVector;
    #else
    register double complex *iVector1, *iVector2, *oVector;
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
        *(oVector+i*inc3) = _Csubcc(*(iVector1+i*inc1), *(iVector2+i*inc2));
        #else
        *(oVector+i*inc3) = *(iVector1+i*inc1) - *(iVector2+i*inc2);
        #endif /*_MSC_VER*/
    }

    return pOVector;
}

VectorD *VectorDC_SubVectorDC2(VectorDC *pIVector1, VectorDC *pIVector2, VectorD *pOVector)
{
    register int i, size, inc1, inc2, inc3;
    #ifdef _MSC_VER
    register _Dcomplex *iVector1, *iVector2;
    #else
    register double complex *iVector1, *iVector2;
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
        *(oVector+i*inc3) = (iVector1+i*inc1)->_Val[0] - (iVector2+i*inc2)->_Val[0];
        #else
        *(oVector+i*inc3) = creal(*(iVector1+i*inc1) - *(iVector2+i*inc2));
        #endif /*_MSC_VER*/
    }

    return pOVector;
}

VectorDC *VectorDC_MulNumD(VectorDC *pIVector, double num, VectorDC *pOVector)
{
    register int i;
    int size, inc1, inc2;

    size = pIVector->size;
    inc1 = pIVector->inc;
    inc2 = pOVector->inc;

    for(i=0;i<size;i++)
    {
        #ifdef _MSC_VER
        (pOVector->vector+i*inc2)->_Val[0] = (pIVector->vector+i*inc1)->_Val[0]*num;
        (pOVector->vector+i*inc2)->_Val[1] = (pIVector->vector+i*inc1)->_Val[1]*num;
        #else
        *(pOVector->vector+i*inc2) = *(pIVector->vector+i*inc1) * num;
        #endif /*_MSC_VER*/
    }

    return pOVector;
}

VectorD *VectorDC_MulNumD2(VectorDC *pIVector, double num, VectorD *pOVector)
{
    register int i;
    int size, inc1, inc2;

    size = pIVector->size;
    inc1 = pIVector->inc;
    inc2 = pOVector->inc;

    for(i=0;i<size;i++)
    {
        #ifdef _MSC_VER
        *(pOVector->vector+i*inc2) = (pIVector->vector+i*inc1)->_Val[0]*num;
        #else
        *(pOVector->vector+i*inc2) = creal(*(pIVector->vector+i*inc1)) * num;
        #endif /*_MSC_VER*/
    }

    return pOVector;
}

#ifdef _MSC_VER
VectorDC *VectorDC_MulNumDC(VectorDC *pIVector, _Dcomplex num, VectorDC *pOVector)
#else
VectorDC *VectorDC_MulNumDC(VectorDC *pIVector, double complex num, VectorDC *pOVector)
#endif /*_MSC_VER*/
{
    int size, inc1, inc2;

    size = pIVector->size;
    inc1 = pIVector->inc;
    inc2 = pOVector->inc;

    zcopy_(&size, pIVector->vector, &inc1, pOVector->vector, &inc2);
    zscal_(&size, &num, pOVector->vector, &inc2);

    return pOVector;
}

#ifdef _MSC_VER
VectorD *VectorDC_MulNumDC2(VectorDC *pIVector, _Dcomplex num, VectorD *pOVector)
#else
VectorD *VectorDC_MulNumDC2(VectorDC *pIVector, double complex num, VectorD *pOVector)
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
        *(pOVector->vector+i*inc2) = (pIVector->vector+i*inc1)->_Val[0]*num._Val[0] - \
            (pIVector->vector+i*inc1)->_Val[1]*num._Val[1];
        #else
        *(pOVector->vector+i*inc2) = creal(*(pIVector->vector+i*inc1) * num);
        #endif /*_MSC_VER*/
    }

    return pOVector;
}

#ifdef _MSC_VER
_Dcomplex VectorDC_MulVectorD(VectorDC *pIVector1, VectorD *pIVector2)
#else
double complex VectorDC_MulVectorD(VectorDC *pIVector1, VectorD *pIVector2)
#endif /*_MSC_VER*/
{
    register int i;
    int size, inc1, inc2;
    #ifdef _MSC_VER
    register _Dcomplex result;
    #else
    register double complex result;
    #endif /*_MSC_VER*/
    #ifdef _MSC_VER
    register _Dcomplex *vector1;
    #else
    register double complex *vector1;
    #endif /*_MSC_VER*/
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
        result._Val[0] += (vector1+i*inc1)->_Val[0] * *(vector2+i*inc2);
        result._Val[1] += (vector1+i*inc1)->_Val[1] * *(vector2+i*inc2);
        #else
        result += *(vector1+i*inc1) * *(vector2+i*inc2);
        #endif /*_MSC_VER*/
    }

    return result;
}

double VectorDC_MulVectorD2(VectorDC *pIVector1, VectorD *pIVector2)
{
    register int i;
    int size, inc1, inc2;
    register double result;
    #ifdef _MSC_VER
    register _Dcomplex *vector1;
    #else
    register double complex *vector1;
    #endif /*_MSC_VER*/
    register double *vector2;

    size = pIVector1->size;
    inc1 = pIVector1->inc;
    inc2 = pIVector2->inc;
    result = 0;
    vector1 = pIVector1->vector;
    vector2 = pIVector2->vector;

    for(i=0;i<size;i++)
    {
        #ifdef _MSC_VER
        result += (vector1+i*inc1)->_Val[0] * *(vector2+i*inc2);
        #else
        result += creal(*(vector1+i*inc1)) * *(vector2+i*inc2);
        #endif /*_MSC_VER*/
    }

    return result;
}

#ifdef _MSC_VER
_Dcomplex VectorDC_MulVectorDC(VectorDC *pIVector1, VectorDC *pIVector2)
#else
double complex VectorDC_MulVectorDC(VectorDC *pIVector1, VectorDC *pIVector2)
#endif /*_MSC_VER*/
{
    int size, inc1, inc2;
    #ifdef _MSC_VER
    _Dcomplex result;
    #else
    double complex result;
    #endif /*_MSC_VER*/

    size = pIVector1->size;
    inc1 = pIVector1->inc;
    inc2 = pIVector2->inc;

    result = zdotu_(&size, pIVector1->vector, &inc1, pIVector2->vector, &inc2);

    return result;
}

double VectorDC_MulVectorDC2(VectorDC *pIVector1, VectorDC *pIVector2)
{
    register int i;
    int size, inc1, inc2;
    register double result;
    #ifdef _MSC_VER
    register _Dcomplex *vector1, *vector2;
    #else
    register double complex *vector1, *vector2;
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
        result += ((vector1+i*inc1)->_Val[0] * (vector2+i*inc2)->_Val[0] - \
            (vector1+i*inc1)->_Val[1] * (vector2+i*inc2)->_Val[1]);
        #else
        result += creal(*(vector1+i*inc1) * *(vector2+i*inc2));
        #endif /*_MSC_VER*/
    }

    return result;
}

VectorDC *VectorDC_EMulVectorD(VectorDC *pIVector1, VectorD *pIVector2, VectorDC *pOVector)
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
        (pOVector->vector+i*inc3)->_Val[0] = (pIVector1->vector+i*inc1)->_Val[0] * *(pIVector2->vector+i*inc2);
        (pOVector->vector+i*inc3)->_Val[1] = (pIVector1->vector+i*inc1)->_Val[1] * *(pIVector2->vector+i*inc2);
        #else
        *(pOVector->vector+i*inc3) = *(pIVector1->vector+i*inc1) * *(pIVector2->vector+i*inc2);
        #endif /*_MSC_VER*/
    }

    return pOVector;
}

VectorD *VectorDC_EMulVectorD2(VectorDC *pIVector1, VectorD *pIVector2, VectorD *pOVector)
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
        *(pOVector->vector+i*inc3) = (pIVector1->vector+i*inc1)->_Val[0] * *(pIVector2->vector+i*inc2);
        #else
        *(pOVector->vector+i*inc3) = creal(*(pIVector1->vector+i*inc1)) * *(pIVector2->vector+i*inc2);
        #endif /*_MSC_VER*/
    }

    return pOVector;
}

VectorDC *VectorDC_EMulVectorDC(VectorDC *pIVector1, VectorDC *pIVector2, VectorDC *pOVector)
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
        *(pOVector->vector+i*inc3) = _Cmulcc(*(pIVector1->vector+i*inc1), *(pIVector2->vector+i*inc2));
        #else
        *(pOVector->vector+i*inc3) = *(pIVector1->vector+i*inc1) * *(pIVector2->vector+i*inc2);
        #endif /*_MSC_VER*/
    }

    return pOVector;
}

VectorD *VectorDC_EMulVectorDC2(VectorDC *pIVector1, VectorDC *pIVector2, VectorD *pOVector)
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
        *(pOVector->vector+i*inc3) = (pIVector1->vector+i*inc1)->_Val[0] * (pIVector2->vector+i*inc2)->_Val[0] - \
            (pIVector1->vector+i*inc1)->_Val[1] * (pIVector2->vector+i*inc2)->_Val[1];
        #else
        *(pOVector->vector+i*inc3) = creal(*(pIVector1->vector+i*inc1) * *(pIVector2->vector+i*inc2));
        #endif /*_MSC_VER*/
    }

    return pOVector;
}

#ifdef _MSC_VER
_Dcomplex VectorDC_CMulVectorD(VectorDC *pIVector1, VectorD *pIVector2)
#else
double complex VectorDC_CMulVectorD(VectorDC *pIVector1, VectorD *pIVector2)
#endif /*_MSC_VER*/
{
    register int i;
    int size, inc1, inc2;
    #ifdef _MSC_VER
    register _Dcomplex result;
    #else
    register double complex result;
    #endif /*_MSC_VER*/
    #ifdef _MSC_VER
    register _Dcomplex *vector1;
    #else
    register double complex *vector1;
    #endif /*_MSC_VER*/
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
        result._Val[0] += (vector1+i*inc1)->_Val[0] * *(vector2+i*inc2);
        result._Val[1] -= (vector1+i*inc1)->_Val[1] * *(vector2+i*inc2);
        #else
        result += conj(*(vector1+i*inc1)) * *(vector2+i*inc2);
        #endif /*_MSC_VER*/
    }

    return result;
}

double VectorDC_CMulVectorD2(VectorDC *pIVector1, VectorD *pIVector2)
{
    register int i;
    int size, inc1, inc2;
    register double result;
    #ifdef _MSC_VER
    register _Dcomplex *vector1;
    #else
    register double complex *vector1;
    #endif /*_MSC_VER*/
    register double *vector2;

    size = pIVector1->size;
    inc1 = pIVector1->inc;
    inc2 = pIVector2->inc;
    result = 0;
    vector1 = pIVector1->vector;
    vector2 = pIVector2->vector;

    for(i=0;i<size;i++)
    {
        #ifdef _MSC_VER
        result += (vector1+i*inc1)->_Val[0] * *(vector2+i*inc2);
        #else
        result += creal(*(vector1+i*inc1)) * *(vector2+i*inc2);
        #endif /*_MSC_VER*/
    }

    return result;
}

#ifdef _MSC_VER
_Dcomplex VectorDC_CMulVectorDC(VectorDC *pIVector1, VectorDC *pIVector2)
#else
double complex VectorDC_CMulVectorDC(VectorDC *pIVector1, VectorDC *pIVector2)
#endif /*_MSC_VER*/
{
    int size, inc1, inc2;
    #ifdef _MSC_VER
    _Dcomplex result;
    #else
    double complex result;
    #endif /*_MSC_VER*/

    size = pIVector1->size;
    inc1 = pIVector1->inc;
    inc2 = pIVector2->inc;

    result = zdotc_(&size, pIVector1->vector, &inc1, pIVector2->vector, &inc2);

    return result;
}

double VectorDC_CMulVectorDC2(VectorDC *pIVector1, VectorDC *pIVector2)
{
    register int i;
    int size, inc1, inc2;
    register double result;
    #ifdef _MSC_VER
    register _Dcomplex *vector1, *vector2;
    #else
    register double complex *vector1, *vector2;
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
        result += ((vector1+i*inc1)->_Val[0] * (vector2+i*inc2)->_Val[0] + \
            (vector1+i*inc1)->_Val[1] * (vector2+i*inc2)->_Val[1]);
        #else
        result += creal(conj(*(vector1+i*inc1)) * *(vector2+i*inc2));
        #endif /*_MSC_VER*/
    }

    return result;
}

double VectorDC_Norm(VectorDC *pIVector)
{
    int size, inc;
    double result;
    register int i;
    #ifdef _MSC_VER
    register _Dcomplex *vector;
    #else
    register double complex *vector;
    #endif /*_MSC_VER*/

    inc = pIVector->inc;
    size = pIVector->size;
    result = 0;
    vector = pIVector->vector;

    for(i=0;i<size;i++)
    {
        result += cabs(*(vector+i*inc));
    }

    return result;
}

#ifdef _MSC_VER
_Dcomplex VectorDC_Norm2(VectorDC *pIVector)
#else
double complex VectorDC_Norm2(VectorDC *pIVector)
#endif /*_MSC_VER*/
{
    int size, inc;
    #ifdef _MSC_VER
    _Dcomplex result;
    #else
    double complex result;
    #endif /*_MSC_VER*/
    register int i;
    #ifdef _MSC_VER
    register _Dcomplex *vector;
    #else
    register double complex *vector;
    #endif /*_MSC_VER*/

    inc = pIVector->inc;
    size = pIVector->size;
    #ifdef _MSC_VER
    result._Val[0] = 0;
    result._Val[1] = 0;
    #else
    result = 0;
    #endif /*_MSC_VER*/
    vector = pIVector->vector;

    for(i=0;i<size;i++)
    {
        #ifdef _MSC_VER
        result._Val[0] += cabs(*(vector+i*inc));
        #else
        result += cabs(*(vector+i*inc));
        #endif /*_MSC_VER*/
    }

    return result;
}

double VectorDC_NormSq(VectorDC *pIVector)
{
    register int i, size, inc;
    register double normSq;
    #ifdef _MSC_VER
    register _Dcomplex *iVector;
    #else
    register double complex *iVector;
    #endif /*_MSC_VER*/

    size = pIVector->size;
    inc = pIVector->inc;
    iVector = pIVector->vector;
    normSq = 0.0;

    for(i=size; i-- > 0;)
    {
        normSq += cabs(*(iVector+i*inc));
    }

    return normSq;
}

#ifdef _MSC_VER
_Dcomplex VectorDC_NormSq2(VectorDC *pIVector)
#else
double complex VectorDC_NormSq2(VectorDC *pIVector)
#endif /*_MSC_VER*/
{
    register int i, size, inc;
    #ifdef _MSC_VER
    register _Dcomplex normSq;
    register _Dcomplex *iVector;
    #else
    register double complex normSq;
    register double complex *iVector;
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
        normSq._Val[0] +=  cabs(*(iVector+i*inc));
        #else
        normSq += *(iVector+i*inc) * *(iVector+i*inc);
        #endif /*_MSC_VER*/
    }

    return normSq;
}
