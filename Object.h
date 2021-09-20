#ifndef _OBJECT_H_
#define _OBJECT_H_

#include <stddef.h>
#include "ObjectClass.h"

extern int instanceof(void *object, const Type *type);

extern Object *Object_Del(Object *pObject);

extern Object *Object_UnWrap(Object *pObject);

#endif /*_OBJECT_H_*/