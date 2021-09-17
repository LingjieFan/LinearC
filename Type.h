#ifndef _TYPE_H_
#define _TYPE_H_

#include "stdlib.h"

typedef struct _Type Type;

struct _Type
{
    const Type *super;
};

extern const Type *NUMD, *NUMDC;

extern const Type *VECTORD, *VECTORDC;

extern const Type *MATRIXD, *MATRIXDC;

#endif /*_TYPE_H_*/
