#include "dmrg_lapack.h"
#include "dmrg_types.h"

void Xlacpy_(const char*                 uplo,
             const IntegerType*          m,
             const IntegerType*          n,
             const std::complex<double>* A,
             const IntegerType*          lda,
             std::complex<double>*       B,
             const IntegerType*          ldb)
{
	zlacpy_(uplo, m, n, A, lda, B, ldb);
}

void Xgemm_(const char*                 transA,
            const char*                 transB,
            const IntegerType*          m,
            const IntegerType*          n,
            const IntegerType*          k,
            const std::complex<double>* alpha,
            const std::complex<double>* A,
            const IntegerType*          lda,
            const std::complex<double>* B,
            const IntegerType*          ldb,
            const std::complex<double>* beta,
            std::complex<double>*       C,
            const IntegerType*          ldc)
{
	zgemm_(transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}

void Xcopy_(const IntegerType*          n,
            const std::complex<double>* X,
            const IntegerType*          incX,
            std::complex<double>*       Y,
            const IntegerType*          incY)
{
	zcopy_(n, X, incX, Y, incY);
}

void Xlacpy_(const char*        uplo,
             const IntegerType* m,
             const IntegerType* n,
             const double*      A,
             const IntegerType* lda,
             double*            B,
             const IntegerType* ldb)
{
	dlacpy_(uplo, m, n, A, lda, B, ldb);
}

void Xgemm_(const char*        transA,
            const char*        transB,
            const IntegerType* m,
            const IntegerType* n,
            const IntegerType* k,
            const double*      alpha,
            const double*      A,
            const IntegerType* lda,
            const double*      B,
            const IntegerType* ldb,
            const double*      beta,
            double*            C,
            const IntegerType* ldc)
{
	dgemm_(transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}

void Xcopy_(const IntegerType* n,
            const double*      X,
            const IntegerType* incX,
            double*            Y,
            const IntegerType* incY)
{
	dcopy_(n, X, incX, Y, incY);
}
