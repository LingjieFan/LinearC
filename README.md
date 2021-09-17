# LinearC

## Makefile
```Makefile
gcc -I ../LinearC main.c ../LinearC/LinearC.c -o test -lblas -llapack -lm
```

## How to use?

1. Use Num, Vector, Matrix separately
```C
#include "Num.h"
```
or 
```C
#include "Vector.h"
```
or 
```C
#include "Matrix.h"
```

2. Use LinearC directly
```C
#include "LinearC.h"
```

## Abstract VS Strong type
There are two kind of objects in LinearC:

Abstract objects:
Num, Vector, Matrix

Strongly typed objects:
NumD, NumDC, VectorD, VectorDC, MatrixD, MatrixDC

### 1.Abstract objects
Writing codes with abstract object makes your codes extensible, 
but may cost a little more time for running program.

### 2.Strongly typed objects
Writing codes with strongly typed object makes your codes faster, 
but non-extensible.

### 3.Comparison
codes with abstract objects
```C
#include "LinearC.h"

void main()
{
    Num *num1;
    Vector *vector1, *vector2, *vector3;
    
    num1 = (Num *)NumD_New(5);
    vector1 = (Vector *)VectorD_New(2);
    vector2 = (Vector *)VectorDC_New(2);
    vector3 = (Vector *)VectorDC_New(2);
    
    vector1 = Vector_Full(vector1, num1);
    vector2 = Vector_Full(vector2, num2);
    vector3 = Vector_AddVector(vector1,vector2,vector3);
    Vector_Show(vector3);
    
    Num_Del(num1);
    Vector_Del(vector1);
    Vector_Del(vector2);
    Vector_Del(vector3);
}
```

codes with strongly typed objects
```C
#include "LinearC.h"

void main()
{
    VectorD *vector1;
    VectorDC *vector2, *vector3;
    
    vector1 = VectorD_New(2);
    vector2 = VectorDC_New(2);
    vector3 = VectorDC_New(2);
    
    vector1 = VectorD_Full(vector1, 5);
    vector2 = VectorDC_Full2(vector2, 5);
    vector3 = VectorD_AddVectorDC(vector1,vector2,vector3);
    VectorDC_Show(vector3);
    
    VectorD_Del(vector1);
    VectorDC_Del(vector2);
    VectorDC_Del(vector3);
}
```

## Wrap and UnWrap
You can use Wrap and Unwrap to create Vector, Matrix from double *pointer or double complex *pointer.

Some flexible operations could be done with LinearC
### 1. Work space
```C
void func(..., ..., double *work)
{
    VectorD *vector1, *vector2;
    
    vector1=VectorD_Wrap(work,2,1);
    vector2=VectorD_Wrap(work+2,2,1);
    
    .....
    
    VectorD_UnWrap(vector1);
    VectorD_UnWrap(vector2);
}
```
### 2. Block matrix/vector
```C
void func()
{
    MatrixD *matrix, *block_matrix;
    
    matrix = MatrixD_New(100,100);
    block_matrix = MatrixD_Wrap(matrix->matrix, 50, 50, 50);
    ...
    
    MatrixD_Del(matrix);
    MatrixD_UnWrap(block_matrix);
}
```

## Example

#### 1.test_Num.c

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

#### 2.test_Vector.c

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

#### 3.test_Matrix.c

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

## Contact us

The author's email: ljfan20@fudan.edu.cn
