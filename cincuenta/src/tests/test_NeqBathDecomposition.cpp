#include "KadanoffBaym.h"
#include "NeqBathDecomposition.h"
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <complex>
#include <vector>

using Dmft::KadanoffBaym;
using Dmft::NeqBathDecomposition;

// ── Helpers ───────────────────────────────────────────────────────────────

namespace {

// NeqBathDecomposition with a single dummy bath site (V=0, ε=0).
// This makes deltaMinusLesser = 0, so -iΔ⁺<(n,j) = -i * delta.lesser(n,j).
static NeqBathDecomposition<double> makeDecomp(SizeType rank, SizeType nT)
{
	PsimagLite::Vector<double>::Type bathParams = { 0.0, 0.0 };
	return NeqBathDecomposition<double>(rank,
	                                    /*beta=*/10.0,
	                                    /*mu=*/0.0,
	                                    bathParams,
	                                    nT,
	                                    /*nTau=*/4,
	                                    /*dt=*/0.2,
	                                    /*dtau=*/0.5);
}

// Fill delta.lesser so that -iΔ⁺<(n,j) = target[n][j].
// With zero first bath: iDeltaPlusLesser = i * delta.lesser,
// so setting delta.lesser(n,j) = i * target gives -iDeltaPlusLesser = target. ✓
static KadanoffBaym<double> makeDelta(SizeType nT, const std::vector<std::vector<double>>& target)
{
	const std::complex<double> I(0, 1);
	KadanoffBaym<double>       delta(nT, /*nTau=*/4, /*dt=*/0.2, /*dtau=*/0.5);
	for (SizeType n = 0; n <= nT; ++n)
		for (SizeType j = 0; j <= nT; ++j)
			delta.lesser(n, j) = I * target[n][j];
	return delta;
}

// Call update for n = 0..nT in sequence.
static void
runUpdates(NeqBathDecomposition<double>& decomp, const KadanoffBaym<double>& delta, SizeType nT)
{
	for (int n = 0; n <= static_cast<int>(nT); ++n)
		decomp.update(n, delta);
}

// Max |Σ_p V(n,p) conj(V(j,p)) − target[n][j]| over all n,j.
static double reconstructionError(NeqBathDecomposition<double>&           decomp,
                                  SizeType                                nT,
                                  const std::vector<std::vector<double>>& target)
{
	double maxErr = 0.0;
	for (SizeType n = 0; n <= nT; ++n) {
		for (SizeType j = 0; j <= nT; ++j) {
			std::complex<double> val(0);
			for (SizeType p = 0; p < decomp.rank(); ++p)
				val += decomp.Vplus(static_cast<int>(n), static_cast<int>(p))
				    * std::conj(
				           decomp.Vplus(static_cast<int>(j), static_cast<int>(p)));
			maxErr = std::max(maxErr, std::abs(val.real() - target[n][j]));
		}
	}
	return maxErr;
}

// Construct a physically-motivated PSD target matrix of rank r:
//   target[n][j] = Σ_{k=0}^{r-1} v[k][n] * v[k][j]
// with v[k][0] = 0 (degenerate at n=0, matching DMFT physics: Δ⁺<(0,j) = 0).
static std::vector<std::vector<double>> makePsdTarget(SizeType nT, int rank)
{
	const SizeType                   N = nT + 1;
	std::vector<std::vector<double>> v(rank, std::vector<double>(N, 0.0));
	for (int r = 0; r < rank; ++r) {
		const double decay = 0.3 + 0.2 * r;
		const double freq  = 0.5 * (r + 1);
		for (SizeType n = 1; n < N; ++n)
			v[r][n] = std::exp(-decay * n) * std::sin(freq * n);
	}

	std::vector<std::vector<double>> mat(N, std::vector<double>(N, 0.0));
	for (SizeType n = 0; n < N; ++n)
		for (SizeType j = 0; j < N; ++j)
			for (int r = 0; r < rank; ++r)
				mat[n][j] += v[r][n] * v[r][j];
	return mat;
}

} // namespace

// ── Seeding / structural tests ─────────────────────────────────────────────

TEST_CASE("n=0 is always skipped — V_(0,p) stays zero", "[NeqBathDecomposition][seeding]")
{
	const SizeType nT = 6;
	for (SizeType L : { 1u, 2u, 3u }) {
		auto target = makePsdTarget(nT, static_cast<int>(L));
		auto delta  = makeDelta(nT, target);
		auto decomp = makeDecomp(L, nT);
		runUpdates(decomp, delta, nT);

		for (SizeType p = 0; p < L; ++p)
			CHECK(std::abs(decomp.Vplus(0, static_cast<int>(p))) < 1e-15);
	}
}

TEST_CASE("L=2: column 0 is alive after fix", "[NeqBathDecomposition][seeding]")
{
	// Column 0 is seeded at n=1 (pivot row 1, not 0). With the old code
	// V_(0,0)≈0 caused a guard-to-zero and column 0 stayed dead (L_eff=1).
	const SizeType nT     = 6;
	auto           target = makePsdTarget(nT, 2);
	auto           delta  = makeDelta(nT, target);
	auto           decomp = makeDecomp(2, nT);
	runUpdates(decomp, delta, nT);

	CHECK(std::abs(decomp.Vplus(1, 0)) > 1e-6); // V_(1,0) non-zero
	CHECK(std::abs(decomp.Vplus(1, 1)) < 1e-15); // col 1 not yet seeded at n=1
	CHECK(std::abs(decomp.Vplus(2, 1)) > 1e-6); // V_(2,1) seeded at n=2
}

TEST_CASE("L=3: all three columns alive after seeding phase", "[NeqBathDecomposition][seeding]")
{
	const SizeType nT     = 8;
	auto           target = makePsdTarget(nT, 3);
	auto           delta  = makeDelta(nT, target);
	auto           decomp = makeDecomp(3, nT);
	runUpdates(decomp, delta, nT);

	CHECK(std::abs(decomp.Vplus(1, 0)) > 1e-6); // col 0 seeded at n=1
	CHECK(std::abs(decomp.Vplus(2, 1)) > 1e-6); // col 1 seeded at n=2
	CHECK(std::abs(decomp.Vplus(3, 2)) > 1e-6); // col 2 seeded at n=3

	// Columns not yet seeded at their seed step must be zero
	CHECK(std::abs(decomp.Vplus(1, 1)) < 1e-15);
	CHECK(std::abs(decomp.Vplus(1, 2)) < 1e-15);
	CHECK(std::abs(decomp.Vplus(2, 2)) < 1e-15);
}

// ── Reconstruction accuracy tests ─────────────────────────────────────────

TEST_CASE("L=1 reconstructs a rank-1 matrix exactly", "[NeqBathDecomposition][reconstruction]")
{
	const SizeType nT     = 10;
	auto           target = makePsdTarget(nT, 1);
	auto           delta  = makeDelta(nT, target);
	auto           decomp = makeDecomp(1, nT);
	runUpdates(decomp, delta, nT);

	// A rank-1 Cholesky of a rank-1 matrix should be exact up to floating point.
	CHECK(reconstructionError(decomp, nT, target) < 1e-10);
}

TEST_CASE("L=2 reconstructs a rank-2 matrix accurately", "[NeqBathDecomposition][reconstruction]")
{
	const SizeType nT     = 10;
	auto           target = makePsdTarget(nT, 2);
	auto           delta  = makeDelta(nT, target);
	auto           decomp = makeDecomp(2, nT);
	runUpdates(decomp, delta, nT);

	// The optimal update minimizes the residual in the column space; for an
	// exactly rank-2 input the reconstruction should be essentially exact.
	CHECK(reconstructionError(decomp, nT, target) < 1e-8);
}

TEST_CASE("L=1 has larger error than L=2 on a rank-2 matrix",
          "[NeqBathDecomposition][reconstruction]")
{
	// This is the regression test that would have caught the L_eff=1 bug:
	// if rank-2 runs secretly deliver rank-1 accuracy, this test fails.
	const SizeType nT     = 10;
	auto           target = makePsdTarget(nT, 2);
	auto           delta  = makeDelta(nT, target);

	auto decomp1 = makeDecomp(1, nT);
	runUpdates(decomp1, delta, nT);
	const double err1 = reconstructionError(decomp1, nT, target);

	auto decomp2 = makeDecomp(2, nT);
	runUpdates(decomp2, delta, nT);
	const double err2 = reconstructionError(decomp2, nT, target);

	// L=2 must be strictly better than L=1.
	CHECK(err2 < err1);
	// And L=1's residual must be measurable (not numerically zero).
	CHECK(err1 > 1e-4);
}

TEST_CASE("L=3 reconstructs a rank-3 matrix accurately", "[NeqBathDecomposition][reconstruction]")
{
	const SizeType nT     = 12;
	auto           target = makePsdTarget(nT, 3);
	auto           delta  = makeDelta(nT, target);
	auto           decomp = makeDecomp(3, nT);
	runUpdates(decomp, delta, nT);

	CHECK(reconstructionError(decomp, nT, target) < 1e-8);
}
