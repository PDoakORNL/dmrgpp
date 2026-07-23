// imethod3 was ported from a plain CPU nested loop to a Kokkos device
// kernel. imethod1 and imethod2 are untouched CPU implementations of the
// same mathematical operation (X += kron(op(A),op(B)) * Y), reachable
// through the same csr_kron_mult_method entry point, so they make a
// ready-made oracle: run the same random inputs through imethod1 and
// imethod3 and require identical (to floating-point tolerance) output.
#include "../csr_kron_mult.cpp"

#include <PsimagLite/Complex.h>
#include <PsimagLite/CrsMatrix.h>
#include <PsimagLite/Matrix.h>
#include <PsimagLite/MatrixNonOwned.h>

#include <Kokkos_Core.hpp>

#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <cmath>
#include <complex>
#include <random>

namespace {

template <typename ComplexOrRealType>
PsimagLite::Matrix<ComplexOrRealType>
randomDense(int nrow, int ncol, std::mt19937& rng, double density)
{
	std::uniform_real_distribution<double> value(-1.0, 1.0);
	std::uniform_real_distribution<double> keep(0.0, 1.0);

	PsimagLite::Matrix<ComplexOrRealType> m(nrow, ncol, ComplexOrRealType(0));
	for (int j = 0; j < ncol; ++j) {
		for (int i = 0; i < nrow; ++i) {
			if (keep(rng) > density)
				continue;

			if constexpr (PsimagLite::IsComplexNumber<ComplexOrRealType>::True) {
				m(i, j) = ComplexOrRealType(value(rng), value(rng));
			} else {
				m(i, j) = static_cast<ComplexOrRealType>(value(rng));
			}
		}
	}
	return m;
}

// Same nrow/ncol bookkeeping as csr_kron_mult_method itself: given
// op(A) and op(B), what shape must Y and X be.
struct KronShape {
	int nrowY;
	int ncolY;
	int nrowX;
	int ncolX;
};

KronShape computeShape(int nrowA, int ncolA, int nrowB, int ncolB, char transA, char transB)
{
	const bool isOpA = (transA == 'T') || (transA == 't') || (transA == 'C') || (transA == 'c');
	const bool isOpB = (transB == 'T') || (transB == 't') || (transB == 'C') || (transB == 'c');

	const int nrow1 = isOpA ? ncolA : nrowA;
	const int ncol1 = isOpA ? nrowA : ncolA;
	const int nrow2 = isOpB ? ncolB : nrowB;
	const int ncol2 = isOpB ? nrowB : ncolB;

	return KronShape { ncol2, ncol1, nrow2, nrow1 };
}

template <typename ComplexOrRealType>
void checkImethod3MatchesImethod1(char transA, char transB, int seed)
{
	using RealType = typename PsimagLite::Real<ComplexOrRealType>::Type;

	std::mt19937 rng(static_cast<unsigned>(seed));

	const int nrowA = 4;
	const int ncolA = 5;
	const int nrowB = 3;
	const int ncolB = 6;

	// density < 1 so imethod3's nnz(A)*nnz(B) pairwise loop actually
	// exercises sparsity, not just a dense matrix in CSR clothing.
	auto denseA = randomDense<ComplexOrRealType>(nrowA, ncolA, rng, 0.5);
	auto denseB = randomDense<ComplexOrRealType>(nrowB, ncolB, rng, 0.5);

	PsimagLite::CrsMatrix<ComplexOrRealType> a(denseA);
	PsimagLite::CrsMatrix<ComplexOrRealType> b(denseB);

	const KronShape shape = computeShape(nrowA, ncolA, nrowB, ncolB, transA, transB);

	auto denseY = randomDense<ComplexOrRealType>(shape.nrowY, shape.ncolY, rng, 1.0);

	PsimagLite::Matrix<ComplexOrRealType> xout1(shape.nrowX, shape.ncolX, ComplexOrRealType(0));
	PsimagLite::Matrix<ComplexOrRealType> xout3(shape.nrowX, shape.ncolX, ComplexOrRealType(0));

	PsimagLite::MatrixNonOwned<const ComplexOrRealType> yinRef(denseY);
	PsimagLite::MatrixNonOwned<ComplexOrRealType>       xout1Ref(xout1);
	PsimagLite::MatrixNonOwned<ComplexOrRealType>       xout3Ref(xout3);

	csr_kron_mult_method(1, transA, transB, a, b, yinRef, xout1Ref);
	csr_kron_mult_method(3, transA, transB, a, b, yinRef, xout3Ref);

	for (int i = 0; i < shape.nrowX; ++i) {
		for (int j = 0; j < shape.ncolX; ++j) {
			const RealType diff  = std::abs(xout1(i, j) - xout3(i, j));
			const RealType scale = RealType(1) + std::abs(xout1(i, j));
			INFO("i=" << i << " j=" << j << " imethod1=" << xout1(i, j)
			          << " imethod3=" << xout3(i, j));
			CHECK(diff < RealType(1e-10) * scale);
		}
	}
}

} // namespace

TEST_CASE("csr_kron_mult imethod3 matches imethod1, real", "[csr_kron_mult][imethod3]")
{
	const auto trans = GENERATE(std::pair<char, char> { 'N', 'N' },
	                            std::pair<char, char> { 'T', 'N' },
	                            std::pair<char, char> { 'N', 'T' },
	                            std::pair<char, char> { 'T', 'T' });
	checkImethod3MatchesImethod1<double>(trans.first, trans.second, 12345);
}

TEST_CASE("csr_kron_mult imethod3 matches imethod1, complex", "[csr_kron_mult][imethod3]")
{
	const auto trans = GENERATE(std::pair<char, char> { 'N', 'N' },
	                            std::pair<char, char> { 'C', 'N' },
	                            std::pair<char, char> { 'N', 'C' },
	                            std::pair<char, char> { 'C', 'C' },
	                            std::pair<char, char> { 'T', 'C' });
	checkImethod3MatchesImethod1<std::complex<double>>(trans.first, trans.second, 67890);
}

int main(int argc, char* argv[])
{
	Kokkos::ScopeGuard scopeGuard(argc, argv);
	return Catch::Session().run(argc, argv);
}
