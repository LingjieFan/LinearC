#ifndef _TENSORCLASS_H_
#define _TENSORCLASS_H_

#include "ObjectClass.h"

typedef struct _Tensor Tensor;

struct _Tensor
{
    Object parent;
    void (*Show)(Tensor *pTensor);
};

#endif /*_TENSOR_H_*/