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
        _type = type->super;
    }

    return 0;
}
