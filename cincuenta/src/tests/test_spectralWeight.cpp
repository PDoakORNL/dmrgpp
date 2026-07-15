#include "ImpuritySolverBase.h"
#include <PsimagLite/Matsubaras.h>
#include <PsimagLite/PsimagLite.h>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <complex>

// spectralWeight computes lim_{omega->inf} -omega * Im(G(i*omega)).
// For the free propagator G(i*omega_n) = 1/(i*omega_n),
// Im(G) = -1/omega_n, so -omega_n * Im(G) = 1 exactly.
// The result must equal 1.0 regardless of beta or N.

TEST_CASE("spectralWeight of free propagator equals 1", "[spectralWeight]")
{
	using RealType          = double;
	using ComplexType       = std::complex<RealType>;
	using SolverBase        = Dmft::ImpuritySolverBase<ComplexType>;
	using MatsubarasType    = PsimagLite::Matsubaras<RealType>;
	using VectorComplexType = typename SolverBase::VectorComplexType;

	const RealType beta = 10.0;
	const SizeType N    = 100; // half the total; total = 2*N frequencies

	MatsubarasType m(beta, N, 0.0);

	VectorComplexType g(m.total());
	for (SizeType i = 0; i < m.total(); ++i) {
		const RealType omega = m.omega(i);
		// G(i*omega) = 1/(i*omega)
		g[i] = ComplexType(1) / (ComplexType(0, 1) * omega);
	}

	const RealType sw = SolverBase::spectralWeight(g, m);
	CHECK(sw == Catch::Approx(1.0).epsilon(1e-12));
}

TEST_CASE("spectralWeight scales linearly with prefactor", "[spectralWeight]")
{
	using RealType          = double;
	using ComplexType       = std::complex<RealType>;
	using SolverBase        = Dmft::ImpuritySolverBase<ComplexType>;
	using MatsubarasType    = PsimagLite::Matsubaras<RealType>;
	using VectorComplexType = typename SolverBase::VectorComplexType;

	const RealType beta = 5.0;
	const SizeType N    = 50;

	MatsubarasType m(beta, N, 0.0);

	VectorComplexType g(m.total());
	for (SizeType i = 0; i < m.total(); ++i) {
		const RealType omega = m.omega(i);
		// G = 0.7/(i*omega) → spectral weight = 0.7
		g[i] = ComplexType(0.7) / (ComplexType(0, 1) * omega);
	}

	const RealType sw = SolverBase::spectralWeight(g, m);
	CHECK(sw == Catch::Approx(0.7).epsilon(1e-12));
}
