#ifndef _OBJECTCLASS_H_
#define _OBJECTCLASS_H_

#include "Type.h"

typedef struct _Object Object;

struct _Object
{
    const Type *type;
};

#endif /*OBJECTCLASS_H_*/