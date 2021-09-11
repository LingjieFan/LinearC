#include "MSC_Complex.h"

#ifdef _MSC_VER

inline _Dcomplex _Caddcc(_Dcomplex x, _Dcomplex y)
{
    _Dcomplex result;
    result._Val[0] = x._Val[0] + y._Val[0];
    result._Val[1] = x._Val[1] + y._Val[1];

    return result;
}

inline _Dcomplex _Csubcc(_Dcomplex x, _Dcomplex y)
{
    _Dcomplex result;
    result._Val[0] = x._Val[0] - y._Val[0];
    result._Val[1] = x._Val[1] - y._Val[1];

    return result;
}

inline _Dcomplex _Cdivcc(_Dcomplex x, _Dcomplex y)
{
    _Dcomplex result, tmp;
    double inv;

    inv = 1.0 / (y._Val[0] * y._Val[0] + y._Val[1] * y._Val[1]);
    tmp._Val[0] = y._Val[0] * inv;
    tmp._Val[1] = -y._Val[1] * inv;
    result = _Cmulcc(x, tmp);

    return result;
}

inline _Dcomplex _Cdivrc(double x, _Dcomplex y)
{
    _Dcomplex result, tmp;
    double inv;

    inv = 1.0 / (y._Val[0] * y._Val[0] + y._Val[1] * y._Val[1]);
    tmp._Val[0] = y._Val[0] * inv;
    tmp._Val[1] = -y._Val[1] * inv;
    result = _Cmulcr(tmp, x);

    return result;
}

#endif /*_MSC_VER*/