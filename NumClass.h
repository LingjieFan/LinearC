#ifndef _NUMCLASS_H_
#define _NUMCLASS_H_

typedef struct _Num Num;
typedef struct _NumD NumD;
typedef struct _NumDC NumDC;

struct _Num
{
    Tensor parent;
    Num *(*Del)(Num *pNum);
    Num *(*AddNum)(Num *pNum1, Num *pNum2, Num *pONum);
    Num *(*SubNum)(Num *pNum1, Num *pNum2, Num *pONum);
    Num *(*MulNum)(Num *pNum1, Num *pNum2, Num *pONum);
    Num *(*DivNum)(Num *pNum1, Num *pNum2, Num *pONum);
};

struct _NumD
{
    Num parent;
    double num;
};

struct _NumDC
{
    Num parent;
    #ifdef _MSC_VER
    _Dcomplex num;
    #else
    double complex num;
    #endif /*_MSC_VER*/
};

#endif /*_NUMCLASS_H_*/