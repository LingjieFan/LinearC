#include <stdio.h>
#include <stdlib.h>
#include "Blas.h"
#include "IObject.h"
#include "Class.h"
#include "Object.h"
#include "VectorD.h"

Class VectorD_Class = {OBJECT};

struct _VectorD
{
    Class* class;
    IObject* iObject;
    int size;
    int inc;
    double* vector;
};

void* _VectorD_New(void* implementor);
void* _VectorD_Del(void* implementor);
int _VectorD_Copy(void* implementor, void* object);
int _VectorD_Equal(void* implementor, void* object);

VectorD* VectorD_New(int size)
{
    VectorD* vector;
    
    vector = (VectorD*)malloc(sizeof(VectorD)+sizeof(IObject)+sizeof(double)*size);
    vector->class = VECTORD;
    vector->iObject = (IObject*)(vector+1);
    vector->iObject->implementor = vector;
    vector->iObject->New = _VectorD_New;
    vector->iObject->Del = _VectorD_Del;
    vector->iObject->Copy = _VectorD_Copy;
    vector->iObject->Equal = _VectorD_Equal;
    vector->size = size;
    vector->inc = 1;
    vector->vector = (double*)(vector->iObject+1); 
    
    return vector;
}

VectorD* VectorD_Del(VectorD* this)
{
    free(this);
    
    return NULL;
}

VectorD* VectorD_Wrap(double* work, int size, int inc)
{
    VectorD* vector;
    
    vector = (VectorD*)malloc(sizeof(VectorD)+sizeof(IObject));
    vector->class = VECTORD;
    vector->iObject = (IObject*)(vector+1);
    vector->iObject->implementor = vector;
    vector->iObject->New = _VectorD_New;
    vector->iObject->Del = _VectorD_Del;
    vector->iObject->Copy = _VectorD_Copy;
    vector->iObject->Equal = _VectorD_Equal;
    vector->size = size;
    vector->inc = inc;
    vector->vector = work;
    
    return vector;
}

VectorD* VectorD_ReWrap(VectorD* this, double* work, int size, int inc)
{
    this->vector = work;
    this->size = size;
    this->inc = inc;

    return this;
}

VectorD* VectorD_UnWrap(VectorD* this)
{
    free(this);
    
    return NULL;
}

IObject* VectorD_GetIObject(VectorD* this)
{
    return this->iObject;
}

int VectorD_GetInc(VectorD* this)
{
    return this->inc;
}

int VectorD_GetSize(VectorD* this)
{
    return this->size;
}

double* VectorD_GetVector(VectorD* this)
{
    return this->vector;
}

double VectorD_GetElement(VectorD* this, int index)
{
    return *(this->vector+index*this->inc);
}

VectorD* VectorD_SetElement(VectorD* this, int index, double value)
{
    *(this->vector+index*this->inc) = value;
    
    return this;
}

void VectorD_Show(VectorD* this)
{
    int size;
    register int i, inc;
    register double* vector;
    
    size = this->size;
    inc = this->inc;
    vector = this->vector;
    printf("VectorD_Show: \n vector:[");
    for(i=0;i<size;i++)
    {
        printf("%g, ", *(vector+inc*i));
    }
    printf("] \n size:%d \n inc:%d \n\n", size, inc);
}

VectorD* VectorD_Full(VectorD* this, double num)
{
    register int i, inc;
    register double* vector;
    int size;
    
    size = this->size;
    inc = this->inc;
    vector = this->vector;
    for(i=0;i<size;i++)
    {
        *(vector+inc*i) = num;
    }
    
    return this;
}

VectorD* VectorD_Linspace(VectorD* this, double start, double end)
{
    register int i, inc;
    int size;
    register double step;
    register double* vector;
    
    size = this->size;
    inc = this->inc;
    step = (end - start) / ((double)size-1);
    vector = this->vector;
    for(i=0;i<size;i++)
    {
       *(vector+i*inc) = start;
       start += step;
    }
    
    return this;
}

VectorD* VectorD_Copy(VectorD* this, VectorD* out_vector)
{
    int inc1, inc2, size;
    double* vector1,* vector2;
    
    inc1 = this->inc;
    inc2 = out_vector->inc;
    size = this->size;
    vector1 = this->vector;
    vector2 = out_vector->vector;
    dcopy_(&size, vector1, &inc1, vector2, &inc2);
    
    return out_vector;
}

double VectorD_Max(VectorD* this)
{
    register int i, size, inc;
    register double max;
    register double* vector;
    
    vector = this->vector;
    max = *(vector);
    size = this->size;
    inc = this->inc;
    for(i=1;i<size;i++)
    {
        max = (*(vector+i*inc) > max) ? *(vector+i*inc) : max;
    }
    
    return max;
}

double VectorD_Min(VectorD* this)
{
    register int i, size, inc;
    register double min;
    register double* vector;
    
    vector = this->vector;
    min = *(vector);
    size = this->size;
    inc = this->inc;
    for(i=1;i<size;i++)
    {
        min = (*(vector+i*inc) < min) ? *(vector+i*inc) : min;
    }
    
    return min;
}

VectorD* VectorD_AddVectorD(VectorD* this, VectorD* in_vector, VectorD* out_vector)
{
    register int i, size, inc1, inc2, inc3;
    register double* vector1,* vector2,* vector3;
    
    size = this->size;
    inc1 = this->inc;
    inc2 = this->inc;
    inc3 = this->inc;
    vector1 = this->vector;
    vector2 = in_vector->vector;
    vector3 = out_vector->vector;
    for(i=0;i<size;i++)
    {
        *(vector3+i*inc3) = *(vector1+i*inc1)+*(vector2+i*inc2);
    }
    
    return out_vector;
}    

VectorD* VectorD_SubVectorD(VectorD* this, VectorD* in_vector, VectorD* out_vector)
{
    register int i, size, inc1, inc2, inc3;
    register double* vector1,* vector2,* vector3;
    
    size = this->size;
    inc1 = this->inc;
    inc2 = this->inc;
    inc3 = this->inc;
    vector1 = this->vector;
    vector2 = in_vector->vector;
    vector3 = out_vector->vector;
    for(i=0;i<size;i++)
    {
        *(vector3+i*inc3) = *(vector1+i*inc1)-*(vector2+i*inc2);
    }
    
    return out_vector;
}    

VectorD* VectorD_MulNumD(VectorD* this, double num, VectorD* out_vector)
{
    int size, inc1, inc2;
    double* vector1,* vector2;
    
    size = this->size;
    inc1 = this->inc;
    inc2 = out_vector->inc;
    vector1 = this->vector;
    vector2 = out_vector->vector;
    dcopy_(&size, vector1, &inc1, vector2, &inc2);
    dscal_(&size, &num, vector2, &inc2);
    
    return out_vector;
}

double VectorD_MulVectorD(VectorD* this, VectorD* in_vector)
{
    int size, inc1, inc2;
    double result;
    double* vector1,* vector2;
    
    size = this->size;
    inc1 = this->inc;
    inc2 = in_vector->inc;
    vector1 = this->vector;
    vector2 = in_vector->vector;
    result = ddot_(&size, vector1, &inc1, vector2, &inc2);
    
    return result;
}

VectorD* VectorD_EMulVectorD(VectorD* this, VectorD* in_vector, VectorD* out_vector)
{
    int size;
    register int i, inc1, inc2, inc3;
    register double* vector1,* vector2,* vector3;
    
    size = this->size;
    inc1 = this->inc;
    inc2 = in_vector->inc;
    inc3 = out_vector->inc;
    vector1 = this->vector;
    vector2 = in_vector->vector;
    vector3 = out_vector->vector;
    for(i=0;i<size;i++)
    {
        *(vector3+i*inc3) = *(vector1+i*inc1) * *(vector2+i*inc2);
    }
    
    return out_vector;
}

double VectorD_Norm(VectorD* this)
{
    int size, inc;
    double result;
    double* vector;
    
    inc = this->inc;
    size = this->size;
    vector = this->vector;
    result = dnrm2_(&size, vector, &inc);
    
    return result;
}

double VectorD_NormSq(VectorD* this)
{
    register int i, size, inc;
    register double normSq;
    register double* vector;
    
    size = this->size;
    inc = this->inc;
    vector = this->vector;
    normSq = 0.0;
    for(i=size;i-->0;)
    {
        normSq += *(vector+i*inc) * *(vector+i*inc);
    }
    
    return normSq;
}

void* _VectorD_New(void* implementor)
{
    VectorD* vector;
    
    vector = (VectorD*)implementor;
    
    return VectorD_New(vector->size);
}

void* _VectorD_Del(void* implementor)
{
    VectorD* vector;
    
    vector = (VectorD*)implementor;
    
    return VectorD_Del(vector);
}

int _VectorD_Copy(void* implementor, void* object)
{
    VectorD* vector1,* vector2;
    
    vector1 = (VectorD*)implementor;
    vector2 = (VectorD*)object;
    VectorD_Copy(vector2, vector1);
    
    return 1;
}

int _VectorD_Equal(void* implementor, void* object)
{
    register int i, inc1, inc2, size;
    VectorD* vector_this,* vector_cmp;
    register double* vector1,* vector2;
    
    vector_this = (VectorD*)implementor;
    vector_cmp = (VectorD*)object;
    size = vector_this->size;
    inc1 = vector_this->inc;
    inc2 = vector_cmp->inc;
    vector1 = vector_this->vector;
    vector2 = vector_cmp->vector;
    if(size == vector_cmp->size)
    {
        for(i=0;i<size;i++)
        {
            if(*(vector1+i*inc1)!=*(vector2+i*inc2))
            {
                return 0;
            }
        }
        
        return 1;
    }
    
    return 0;
}
