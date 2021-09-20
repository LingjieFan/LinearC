#include "Tensor.h"

Tensor *Tensor_Del(Tensor *pTensor)
{
    if(pTensor == NULL)
    {
        return NULL;
    }

    return (*((Object *)pTensor)->Del)((Object *)pTensor);
}

Tensor *Tensor_UnWrap(Tensor *pTensor)
{
    if(pTensor == NULL)
    {
        return NULL;
    }

    return (*((Object *)pTensor)->UnWrap)((Object *)pTensor);
}

void *Tensor_Show(Tensor *pTensor)
{
    return (*pTensor->Show)(pTensor);
}