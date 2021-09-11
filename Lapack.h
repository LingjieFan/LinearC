#ifndef _LAPACK_H_
#define _LAPACK_H_

#ifdef _MSC_VER
#ifdef _WIN64
#pragma comment(lib, "./Win64/liblapack.lib")
#else
#pragma comment(lib, "./Win32/liblapack.lib")
#endif /*_WIN64*/
#endif /*_MSC_VER*/

extern void dgetrf_(int *row, int *col, double *matrixA, int *lda, int *ipiv, int *info);

extern void dgetri_(int *sizeA, double *matrixA, int *lda, int *ipiv, double *work, \
    int *lwork, int *info);

extern void dgetrs_(char *trans, int *sizeA, int *sizeB, double *matrixA, int *lda, int *ipiv, double *matrixB, \
    int *ldb, int *info);

extern void dgeev_(char *jobvl, char *jobvr, int *sizeA, double *matrixA, int *lda, double *w, \
    double *matrixVL, int *ldvl, double *matrixVR, int *ldvr, double *work, int *lwork, \
    double *rwork, int *info);

#ifdef _MSC_VER

extern void zgetrf_(int *row, int *col, _Dcomplex *matrixA, int *lda, int *ipiv, int *info);

extern void zgetri_(int *sizeA, _Dcomplex *matrixA, int *lda, int *ipiv, _Dcomplex *work, \
    int *lwork, int *info);

extern void zgeev_(char *jobvl, char *jobvr, int *sizeA, _Dcomplex *matrixA, int *lda, _Dcomplex *w, \
    _Dcomplex *matrixVL, int *ldvl, _Dcomplex *matrixVR, int *ldvr, _Dcomplex *work, int *lwork, \
    double *rwork, int *info);

extern void zgetrs_(char *trans, int *sizeA, int *sizeB, _Dcomplex *matrixA, int *lda, int *ipiv, _Dcomplex *matrixB, \
    int *ldb, int *info);

#else

extern void zgetrf_(int *row, int *col, double complex *matrixA, int *lda, int *ipiv, int *info);

extern void zgetri_(int *sizeA, double complex *matrixA, int *lda, int *ipiv, double complex *work, \
    int *lwork, int *info);

extern void zgeev_(char *jobvl, char *jobvr, int *sizeA, double complex *matrixA, int *lda, double complex *w, \
    double complex *matrixVL, int *ldvl, double complex *matrixVR, int *ldvr, double complex *work, int *lwork, \
    double *rwork, int *info);

extern void zgetrs_(char *trans, int *sizeA, int *sizeB, double complex *matrixA, int *lda, int *ipiv, double complex *matrixB, \
    int *ldb, int *info);

#endif /*_MSC_VER*/

#endif /*_LAPACK_H_*/
