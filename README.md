# LinearC

## Description
C Math Library

## Software Architecture
---------------------------

LinearC

---------------------------

Matrix: MatrixD MatrixDC...

---------------------------

Vector: VectorD VectorDC...

---------------------------

Num: NumD NumDC...

---------------------------

## Ing.......

1. This project is buiding.
2. Num, Vector Matrix LinearC are available.
3. It could run on Linux and Windows.

## TODO

1. test it on Windows

2. complete some unreasonable cases for Matrix; those of Num, Vector have been done.

3. optimize the code for speed.

## How to use?

1. #include "Num.h" or #include "Vector.h" or #include "Matrix.h"

2. or directly #include "LinearC.h"

## Example

#### test_Num.c

```C
#include "Num.h"
#include <stdio.h>

void main()
{
    Num *num1, *num2, *num3;

    num1 = (Num *)NumD_New(2);
    num2 = (Num *)NumDC_New(1+1*I);
    num3 = (Num *)NumDC_New(2);

    num3 = Num_AddNum(num1, num2, num3);
    printf("%lf,%lf\n",creal(((NumDC *)num3)->num),cimag(((NumDC *)num3)->num));

    num3 = Num_SubNum(num1, num2, num3);
    printf("%lf,%lf\n",creal(((NumDC *)num3)->num),cimag(((NumDC *)num3)->num));

    num3 = Num_MulNum(num1, num2, num3);
    printf("%lf,%lf\n",creal(((NumDC *)num3)->num),cimag(((NumDC *)num3)->num));

    num3 = Num_DivNum(num1, num2, num3);
    printf("%lf,%lf\n",creal(((NumDC *)num3)->num),cimag(((NumDC *)num3)->num));

    num1 = Num_Del(num1);
    num2 = Num_Del(num2);
    num3 = Num_Del(num3);
}
```

#### test_Vector.c

```C
#include <stdio.h>
#include "Num.h"
#include "Vector.h"

void main()
{
    Num *num1, *num2;
    Vector *vector1, *vector2, *vector3;

    num1 = (Num *)NumD_New(1);
    num2 = (Num *)NumDC_New(2+3*I);
    vector1 = (Vector *)VectorD_New(2);
    vector1 = Vector_Full(vector1, num1);
    Vector_Show(vector1);
    vector2 = (Vector *)VectorDC_New(2);
    vector2 = Vector_Linspace(vector2, num1, num2);
    Vector_Show(vector2);
    vector3 = (Vector *)VectorDC_New(2);

    num1 = Vector_Max(vector1, num1);
    printf("%lf",((NumD *)num1)->num);

    num1 = Vector_Min(vector1, num1);
    printf("%lf", ((NumD *)num1)->num);

    vector3 = Vector_AddVector(vector1, vector2, vector3);
    Vector_Show(vector3);

    vector3 = Vector_SubVector(vector1, vector2, vector3);
    Vector_Show(vector3);

    vector3 = Vector_MulNum(vector1, num1,vector3);
    Vector_Show(vector3);

    Num_Del(num1);
    Num_Del(num2);
    Vector_Del(vector1);
    Vector_Del(vector2);
    Vector_Del(vector3);
}
```

#### test_Matrix.c

```C
#include "Matrix.h"

void main()
{
    Num *diag, *offDiag1, *offDiag2;
    Matrix *matrix1, *matrix2, *matrix3;

    offDiag1 = (Num *)NumD_New(2);
    offDiag2 = (Num *)NumDC_New(2+1*I);
    diag = (Num *)NumD_New(1);

    matrix1 = (Matrix *)MatrixD_New(2,2);
    matrix2 = (Matrix *)MatrixDC_New(2,2);
    matrix3 = (Matrix *)MatrixDC_New(2,2);

    matrix1 = Matrix_Set(matrix1, diag, offDiag1);
    Matrix_Show(matrix1);
    matrix2 = Matrix_Set(matrix2, diag, offDiag2);
    Matrix_Show(matrix2);

    matrix3 = Matrix_MulMatrix(matrix1, matrix2, matrix3);
    Matrix_Show(matrix3);

    Num_Del(diag);
    Num_Del(offDiag1);
    Num_Del(offDiag2);
    Matrix_Del(matrix1);
    Matrix_Del(matrix2);
    Matrix_Del(matrix3);
}
```

## Contact Us

The author's email: ljfan20@fudan.edu.cn
