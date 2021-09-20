#ifndef _TENSOR_H_
#define _TENSOR_H_

#include <stdio.h>
#include "TensorClass.h"

extern Tensor *Tensor_Del(Tensor *pTensor);

extern Tensor *Tensor_UnWrap(Tensor *pTensor);

extern void Tensor_Show(Tensor *pTensor);

#endif /*_TENSOR_H_*/
