#ifndef _NUM_H_
#define _NUM_H_

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>

#include <complex.h>
#include "MSC_Complex.h"

#include "Object.h"
#include "TensorClass.h"
#include "NumClass.h"
#include "NumD.h"
#include "NumDC.h"

extern Num *Num_Del(Num *pNum);

extern Num *Num_AddNum(Num *pNum1, Num *pNum2, Num *pONum);

extern Num *Num_SubNum(Num *pNum1, Num *pNum2, Num *pONum);

extern Num *Num_MulNum(Num *pNum1, Num *pNum2, Num *pONum);

extern Num *Num_DivNum(Num *pNum1, Num *pNum2, Num *pONum);

#endif /*_NUM_H_*/
