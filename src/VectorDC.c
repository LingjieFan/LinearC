#include <stdio.h>
#include <stdlib.h>
#include "Blas.h"
#include "IObject.h"
#include "Class.h"
#include "Object.h"
#include "MatrixDC.h"
#include "VectorDC.h"

Class VectorDC_Class = {OBJECT};

struct _VectorDC
{
    Class* class;
    IObject* iObject;
    double complex* vector;
    int size;
    int inc;
};

static void* _VectorDC_New(void* implementor);
static void* _VectorDC_Del(void* implementor);
static int _VectorDC_Copy(void* implementor, void* object);
static int _VectorDC_Equal(void* implementor, void* object);

VectorDC* VectorDC_New(int size)
{
    VectorDC* vector;
    
    vector = (VectorDC*)malloc(sizeof(VectorDC)+sizeof(IObject)+sizeof(double complex)*size);
    vector->class = VECTORDC;
    vector->iObject = (IObject*)(vector+1);
    vector->iObject->implementor = vector;
    vector->iObject->New = _VectorDC_New;
    vector->iObject->Del = _VectorDC_Del;
    vector->iObject->Copy = _VectorDC_Copy;
    vector->iObject->Equal = _VectorDC_Equal;
    vector->size = size;
    vector->inc = 1;
    vector->vector = (double complex*)(vector->iObject+1);
    
    return vector;
}

VectorDC* VectorDC_Del(VectorDC* this)
{
    free(this);
    
    return NULL;
}

VectorDC* VectorDC_Wrap(double complex* work, int size, int inc)
{
    VectorDC* vector;
    
    vector = (VectorDC*)malloc(sizeof(VectorDC)+sizeof(IObject));
    vector->class = VECTORDC;
    vector->iObject = (IObject*)(vector+1);
    vector->iObject->implementor = vector;
    vector->iObject->New = _VectorDC_New;
    vector->iObject->Del = _VectorDC_Del;
    vector->iObject->Copy = _VectorDC_Copy;
    vector->iObject->Equal = _VectorDC_Equal;
    vector->size = size;
    vector->inc = inc;
    vector->vector = work;
    
    return vector;
}

VectorDC* VectorDC_ReWrap(VectorDC* this, double complex* work, int size, int inc)
{
    this->vector = work;
    this->size = size;
    this->inc = inc;

    return this;
}

VectorDC* VectorDC_UnWrap(VectorDC* this)
{
    free(this);
    
    return NULL;
}

IObject* VectorDC_GetIObject(VectorDC* this)
{
    return this->iObject;
}

int VectorDC_GetInc(VectorDC* this)
{
    return this->inc;
}

int VectorDC_GetSize(VectorDC* this)
{
    return this->size;
}

double complex* VectorDC_GetVector(VectorDC* this)
{
    return this->vector;
}

double complex VectorDC_GetElement(VectorDC* this, int index)
{
    return *(this->vector+index*this->inc);
}

VectorDC* VectorDC_SetElement(VectorDC* this, int index, double complex value)
{
    *(this->vector+index*this->inc) = value;
    
    return this;
}

void VectorDC_Show(VectorDC* this)
{
    int size;
    register int i, inc;
    register double complex* vector;
    
    size = this->size;
    inc = this->inc;
    vector = this->vector;
    printf("VectorDC_Show: \n vector:[");
    for(i=0;i<size;i++)
    {
        printf("%g+%g j, ", creal(*(vector+inc*i)), cimag(*(vector+inc*i)));
    }
    printf("] \n size:%d \n inc:%d \n\n", size, inc);
}

VectorDC* VectorDC_Full(VectorDC* this, double complex num)
{
    register int i, inc;
    register double complex* vector;
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

VectorDC* VectorDC_Linspace(VectorDC* this, double complex start, double complex end)
{
    register int i, inc;
    int size;
    register double complex step;
    register double complex* vector;
    
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

VectorDC* VectorDC_Copy(VectorDC* this, VectorDC* out_vector)
{
    int inc1, inc2, size;
    double complex* vector1,* vector2;
    
    inc1 = this->inc;
    inc2 = out_vector->inc;
    size = this->size;
    vector1 = this->vector;
    vector2 = out_vector->vector;
    zcopy_(&size, vector1, &inc1, vector2, &inc2);
    
    return out_vector;
}

VectorDC* VectorDC_AddVectorDC(VectorDC* this, VectorDC* in_vector, VectorDC* out_vector)
{
    register int i, size, inc1, inc2, inc3;
    register double complex* vector1,* vector2,* vector3;
    
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

VectorDC* VectorDC_SubVectorDC(VectorDC* this, VectorDC* in_vector, VectorDC* out_vector)
{
    register int i, size, inc1, inc2, inc3;
    register double complex* vector1,* vector2,* vector3;
    
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

VectorDC* VectorDC_MulNumDC(VectorDC* this, double complex num, VectorDC* out_vector)
{
    int size, inc1, inc2;
    double complex* vector1,* vector2;
    
    size = this->size;
    inc1 = this->inc;
    inc2 = out_vector->inc;
    vector1 = this->vector;
    vector2 = out_vector->vector;
    zcopy_(&size, vector1, &inc1, vector2, &inc2);
    zscal_(&size, &num, vector2, &inc2);
    
    return out_vector;
}

double complex VectorDC_MulVectorDC(VectorDC* this, VectorDC* in_vector)
{
    int size, inc1, inc2;
    double complex result;
    double complex* vector1,* vector2;
    
    size = this->size;
    inc1 = this->inc;
    inc2 = this->inc;
    vector1 = this->vector;
    vector2 = in_vector->vector;
    result = zdotu_(&size, vector1, &inc1, vector2, &inc2);
    
    return result;
}

double complex VectorDC_CMulVectorDC(VectorDC* this, VectorDC* in_vector)
{
    int size, inc1, inc2;
    double complex result;
    double complex* vector1,* vector2;
    
    size = this->size;
    inc1 = this->inc;
    inc2 = this->inc;
    vector1 = this->vector;
    vector2 = in_vector->vector;
    result = zdotc_(&size, vector1, &inc1, vector2, &inc2);
    
    return result;
}

VectorDC* VectorDC_EMulVectorDC(VectorDC* this, VectorDC* in_vector, VectorDC* out_vector)
{
    int size;
    register int i, inc1, inc2, inc3;
    register double complex* vector1,* vector2,* vector3;
    
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

MatrixDC* VectorDC_TMulVectorDC(VectorDC* this, VectorDC* in_vector, MatrixDC* out_matrix)
{
    int size1, size2;
    register int i, j, inc1, inc2, ld;
    register double complex* vector1,* vector2,* matrix;

    size1 = this->size;
    size2 = in_vector->size;
    vector1 = this->vector;
    vector2 = in_vector->vector;
    matrix = MatrixDC_GetMatrix(out_matrix);
    inc1 = this->inc;
    inc2 = in_vector->inc;
    ld = MatrixDC_GetLd(out_matrix);

    for(i=size1;i-->0;)
    {
        register double tmp;

        tmp = *(vector1+i*inc1);
        for(j=size2;j-->0;)
        {
            *(matrix+i*ld+j) = tmp * *(vector2+i*inc2);
        }
    }

    return out_matrix;
}

static void* _VectorDC_New(void* implementor)
{
    VectorDC* vector;
    
    vector = (VectorDC*)implementor;
    
    return VectorDC_New(vector->size);
}

static void* _VectorDC_Del(void* implementor)
{
    VectorDC* vector;
    
    vector = (VectorDC*)implementor;
    
    return VectorDC_Del(vector);
}

static int _VectorDC_Copy(void* implementor, void* object)
{
    VectorDC* vector1,* vector2;
    
    vector1 = (VectorDC*)implementor;
    vector2 = (VectorDC*)object;
    VectorDC_Copy(vector2, vector1);
    
    return 1;
}

static int _VectorDC_Equal(void* implementor, void* object)
{
    register int i, inc1, inc2, size;
    VectorDC* vector_this,* vector_cmp;
    register double complex* vector1,* vector2;
    
    vector_this = (VectorDC*)implementor;
    vector_cmp = (VectorDC*)object;
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

