#include "Object.h"

int instanceof(void *object, const Type *type)
{
    const Object *_object;
    const Type *_type;

    _type = type;
    _object = object;

    if(_object == NULL)
    {
        return 0;
    }

    while(_type != NULL)
    {
        if(_object->type == type)
        {
            return 1;
        }
        _type = _type->super;
    }

    return 0;
}

Object *Object_Del(Object *pObject)
{
    if(pObject == NULL)
    {
        return NULL;
    }
    return (*pObject->Del)(pObject);
}

Object *Object_UnWrap(Object *pObject)
{
    if(pObject == NULL)
    {
        return NULL;
    }
    return (*pObject->UnWrap)(pObject);
}
