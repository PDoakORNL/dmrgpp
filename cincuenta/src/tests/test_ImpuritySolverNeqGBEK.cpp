#include "CincuentaInputCheck.h"
#include "ImpuritySolverNeqGBEK.h"
#include "TestMatrixUtils.h"
#include <catch2/catch_test_macros.hpp>

// ── Types ────────────────────────────────────────────────────────────────────

using ComplexType   = std::complex<double>;
using RealType      = double;
using SolverType    = Dmft::ImpuritySolverNeqGBEK<ComplexType>;
using ParamsType    = Dmft::ParamsNeqDmftSolver<ComplexType>;
using InputNgType   = PsimagLite::InputNg<Dmft::CincuentaInputCheck>;
using VectorComplex = typename PsimagLite::Vector<ComplexType>::Type;
using CrsMatrixType = PsimagLite::CrsMatrix<ComplexType>;
using WordType      = LanczosPlusPlus::LanczosGlobals::WordType;
using Acc           = Dmft::GBEKTestAccessor;

// ── Ainur config ─────────────────────────────────────────────────────────────
// Minimal GBEK config: 1 first-bath site, 2 time steps, rank-1 second bath.
// Equilibrium solver is ExactDiag; GBEK path invoked for neq propagation.
// DMRG fields are included because ParamsDmftSolver requires them even when
// the neq solver is GBEK.
static const std::string kConfig = "##Ainur1.0\n\n"
                                   "FicticiousBeta=10;\n"
                                   "ChemicalPotential=0.;\n"
                                   "Matsubaras=20;\n"
                                   "LatticeGf=\"energy,semicircular,4\";\n"
                                   "NumberOfBathPoints=1;\n"
                                   "DmftNumberOfIterations=1;\n"
                                   "DmftTolerance=1e-3;\n"
                                   "ImpuritySolver=\"exactdiag\";\n"
                                   "FitOptions=particleholesymmetric;\n"
                                   "MinParamsDelta=0.01;\n"
                                   "MinParamsDelta2=0.01;\n"
                                   "MinParamsTolerance=1e-4;\n"
                                   "MinParamsMaxIter=100;\n"
                                   "MinParamsVerbose=0;\n"
                                   "TargetElectronsUp=1;\n"
                                   "TargetElectronsDown=1;\n"
                                   "int ImpuritySite=0;\n"
                                   "real HubbardU=2.;\n"
                                   "HubbardUFinal=2.;\n"
                                   "RootOutputname=\"testGBEK\";\n"
                                   "InfiniteLoopKeptStates=20;\n"
                                   "matrix FiniteLoopsGs=[[@auto, 20, 0],[@auto, 20, 0]];\n"
                                   "real OmegaBegin=-4.;\n"
                                   "integer OmegaTotal=40;\n"
                                   "real OmegaStep=0.2;\n"
                                   "real OmegaDelta=0.2;\n"
                                   "integer TridiagSteps=100;\n"
                                   "real TridiagEps=1e-6;\n"
                                   "TruncationTolerance=\"1e-6,20\";\n"
                                   "CorrectionVectorEta=0.;\n"
                                   "GsWeight=0.1;\n"
                                   "matrix FiniteLoopsOmega=[[@auto, 20, 2],[@auto, 20, 2]];\n"
                                   "TmaxNeq=0.2;\n"
                                   "NtNeq=2;\n"
                                   "NeqDmftIter=1;\n"
                                   "NeqDmftTolerance=1e-4;\n"
                                   "NeqSolver=\"gbek\";\n"
                                   "NeqBathRank=1;\n"
                                   "BandwidthFinal=3.;\n";

// ── Shared fixture ────────────────────────────────────────────────────────────
// Constructed once (static) so the full solve() runs only on first access.
// Members are declared in construction order so lifetimes are nested correctly:
//   ioW -> io -> params -> solver.
struct GBEKFixture {
	InputNgType::Writeable ioW;
	InputNgType::Readable  io;
	ParamsType             params;
	SolverType             solver;

	GBEKFixture()
	    : ioW(Dmft::CincuentaInputCheck {}, kConfig)
	    , io(ioW)
	    , params(io)
	    , solver(params, io)
	{
		// V=0.5, epsilon=0 for the single first-bath site.
		solver.solve({ 0.5, 0.0 });
	}
};

static GBEKFixture& fixture()
{
	static GBEKFixture f;
	return f;
}

// ── Helpers ───────────────────────────────────────────────────────────────────

// Returns a test vMid with bathRank entries using known complex values.
static std::vector<ComplexType> makeTestVMid(SizeType bathRank)
{
	std::vector<ComplexType> v(bathRank);
	for (SizeType p = 0; p < bathRank; ++p)
		v[p] = ComplexType(0.3 + 0.1 * p, 0.2 - 0.05 * p);
	return v;
}

// Compare CSR SpMV against applyHext for a single (v, vMid) combination.
// Returns max|hv_csr - hv_ref|.
static RealType spMVDiff(const SolverType&               solver,
                         const VectorComplex&            v,
                         const std::vector<ComplexType>& vMid,
                         const CrsMatrixType&            csr,
                         const std::vector<WordType>&    upW,
                         const std::vector<WordType>&    dnW,
                         SizeType                        dim1)
{
	const SizeType dim = upW.size() * dnW.size();
	VectorComplex  hv_csr(dim, ComplexType(0));
	VectorComplex  hv_ref(dim, ComplexType(0));

	csr.matrixVectorProduct(hv_csr, v);
	Acc::applyHext(solver, v, hv_ref, vMid, upW, dnW, dim1);

	RealType maxDiff = 0;
	for (SizeType i = 0; i < dim; ++i)
		maxDiff = std::max(maxDiff, std::abs(hv_csr[i] - hv_ref[i]));
	return maxDiff;
}

// ── TC-1: applyHext is Hermitian ─────────────────────────────────────────────
// Validates the reference implementation before using it to check anything else.
// If this fails, the physics was broken before the CSR existed.
TEST_CASE("applyHext is Hermitian for N-1 and N+1 sectors", "[GBEK][applyHext]")
{
	auto&          f = fixture();
	auto&          s = f.solver;
	const SizeType r = Acc::bathRank(s);

	SECTION("N-1 sector, zero vMid (fixed part only)")
	{
		std::vector<ComplexType> zeroVMid(r, ComplexType(0));
		auto                     H = buildDenseFromApplyHext(
                    s, zeroVMid, Acc::upWordsNm1(s), Acc::dnWordsNm1(s), Acc::dim1Nm1(s));
		assertHermitian(H);
	}

	SECTION("N-1 sector, non-zero vMid")
	{
		auto vMid = makeTestVMid(r);
		auto H    = buildDenseFromApplyHext(
                    s, vMid, Acc::upWordsNm1(s), Acc::dnWordsNm1(s), Acc::dim1Nm1(s));
		assertHermitian(H);
	}

	SECTION("N+1 sector, non-zero vMid")
	{
		auto vMid = makeTestVMid(r);
		auto H    = buildDenseFromApplyHext(
                    s, vMid, Acc::upWordsNp1(s), Acc::dnWordsNp1(s), Acc::dim1Np1(s));
		assertHermitian(H);
	}
}

// ── TC-2: buildHextCSR produces a Hermitian matrix ───────────────────────────
// If TC-1 passes but TC-2 fails the bug is in buildHextCSR, not applyHext.
TEST_CASE("buildHextCSR produces Hermitian matrix for both sectors", "[GBEK][buildHextCSR]")
{
	auto& f = fixture();
	auto& s = f.solver;

	SECTION("N-1 sector after updateCSR with test vMid")
	{
		auto vMid = makeTestVMid(Acc::bathRank(s));
		Acc::updateCSR(s, Acc::csrNm1Mut(s), Acc::varNm1(s), vMid);
		auto H = buildDenseFromCSR(Acc::csrNm1(s));
		assertHermitian(H);
	}

	SECTION("N+1 sector after updateCSR with test vMid")
	{
		auto vMid = makeTestVMid(Acc::bathRank(s));
		Acc::updateCSR(s, Acc::csrNp1Mut(s), Acc::varNp1(s), vMid);
		auto H = buildDenseFromCSR(Acc::csrNp1(s));
		assertHermitian(H);
	}
}

// ── TC-3: updateCSR writes variable slots and preserves fixed slots ───────────
TEST_CASE("updateCSR modifies only variable slots", "[GBEK][updateCSR]")
{
	auto&          f    = fixture();
	auto&          s    = f.solver;
	const SizeType r    = Acc::bathRank(s);
	auto           vMid = makeTestVMid(r);

	// Work on N-1 sector.
	auto&          csr    = Acc::csrNm1Mut(s);
	const auto&    varVec = Acc::varNm1(s);
	const SizeType nnz    = csr.nonZeros();

	// Snapshot all values before update.
	std::vector<ComplexType> before(nnz);
	for (SizeType i = 0; i < nnz; ++i)
		before[i] = csr.getValue(i);

	Acc::updateCSR(s, csr, varVec, vMid);

	// Build set of variable indices for O(1) membership test.
	std::vector<bool> isVar(nnz, false);
	for (const auto& ve : varVec)
		isVar[ve.nnzIdx] = true;

	for (SizeType i = 0; i < nnz; ++i) {
		if (isVar[i]) {
			// Expected value computed from VarEntry.
			// (We check it via the aggregate SpMV test; here we just verify it
			// changed.) A variable slot set from a non-zero vMid should differ from its
			// zero build-time value unless the VarEntry happens to produce zero. Accept
			// either "changed" (normal) or "correctly zero" (edge case).
			(void)i; // checked in TC-4
		} else {
			// Fixed slots must be unchanged.
			INFO("Fixed slot " << i << " was modified by updateCSR");
			CHECK(csr.getValue(i) == before[i]);
		}
	}

	// Verify each variable slot has the value that the VarEntry formula predicts.
	for (const auto& ve : varVec) {
		const ComplexType Vp       = ve.isConj ? std::conj(vMid[ve.p]) : vMid[ve.p];
		const ComplexType expected = Vp * ComplexType(static_cast<RealType>(ve.sign));
		INFO("VarEntry nnzIdx=" << ve.nnzIdx << " p=" << ve.p << " isConj=" << ve.isConj
		                        << " sign=" << ve.sign);
		CHECK(std::abs(csr.getValue(ve.nnzIdx) - expected) < 1e-15);
	}
}

// ── TC-4: SpMV matches applyHext for multiple (v, vMid) combinations ─────────
TEST_CASE("CSR SpMV agrees with applyHext for multiple vectors and vMid", "[GBEK][SpMV][applyHext]")
{
	auto&          f = fixture();
	auto&          s = f.solver;
	const SizeType r = Acc::bathRank(s);

	// Work with the N-1 sector.
	const auto&    upW    = Acc::upWordsNm1(s);
	const auto&    dnW    = Acc::dnWordsNm1(s);
	const SizeType d      = upW.size() * dnW.size();
	const SizeType d1     = Acc::dim1Nm1(s);
	auto&          csr    = Acc::csrNm1Mut(s);
	const auto&    varVec = Acc::varNm1(s);

	REQUIRE(d > 0);

	// Three vMid cases.
	std::vector<std::vector<ComplexType>> vMidCases = {
	    std::vector<ComplexType>(r, ComplexType(0)),     // zero
	    makeTestVMid(r),                                  // complex test values
	    [r]() {                                           // purely imaginary
	        std::vector<ComplexType> v(r);
	        for (SizeType p = 0; p < r; ++p)
	            v[p] = ComplexType(0, 0.7 + 0.1 * p);
	        return v;
	    }()
	};

	// Three test vectors.
	VectorComplex e0(d, ComplexType(0));
	e0[0] = ComplexType(1);
	VectorComplex emid(d, ComplexType(0));
	emid[d / 2] = ComplexType(1);
	// Simple normalised linear combination.
	VectorComplex esum(d, ComplexType(0));
	if (d >= 2) {
		esum[0]     = ComplexType(1) / std::sqrt(2.0);
		esum[d - 1] = ComplexType(1) / std::sqrt(2.0);
	} else {
		esum[0] = ComplexType(1);
	}

	const std::vector<const VectorComplex*> vecs = { &e0, &emid, &esum };

	for (SizeType vi = 0; vi < vMidCases.size(); ++vi) {
		Acc::updateCSR(s, csr, varVec, vMidCases[vi]);
		for (SizeType ti = 0; ti < vecs.size(); ++ti) {
			const RealType diff
			    = spMVDiff(s, *vecs[ti], vMidCases[vi], csr, upW, dnW, d1);
			INFO("vMid case " << vi << ", test vec " << ti
			                  << ": max|CSR*v - applyHext(v)| = " << diff);
			CHECK(diff < 1e-12);
		}
	}
}

// ── TC-5: krylovExpmvCSR agrees with krylovExpmv ────────────────────────────
TEST_CASE("krylovExpmvCSR agrees with krylovExpmv to 1e-10", "[GBEK][Krylov]")
{
	auto&          f    = fixture();
	auto&          s    = f.solver;
	const SizeType r    = Acc::bathRank(s);
	auto           vMid = makeTestVMid(r);

	const auto&    upW    = Acc::upWordsNm1(s);
	const auto&    dnW    = Acc::dnWordsNm1(s);
	const SizeType d      = upW.size() * dnW.size();
	const SizeType d1     = Acc::dim1Nm1(s);
	auto&          csr    = Acc::csrNm1Mut(s);
	const auto&    varVec = Acc::varNm1(s);

	// dt from TmaxNeq/NtNeq = 0.2/2 = 0.1
	const RealType dt = 0.1;

	Acc::updateCSR(s, csr, varVec, vMid);

	auto checkPair = [&](const VectorComplex& psi, const char* label)
	{
		auto result_ref = Acc::krylovExpmv(s, psi, vMid, upW, dnW, d1, dt);
		auto result_csr = Acc::krylovExpmvCSR(s, psi, Acc::csrNm1(s), dt);

		REQUIRE(result_ref.size() == d);
		REQUIRE(result_csr.size() == d);

		RealType maxDiff = 0;
		for (SizeType i = 0; i < d; ++i)
			maxDiff = std::max(maxDiff, std::abs(result_csr[i] - result_ref[i]));

		INFO(label << ": max|krylovCSR - krylovRef| = " << maxDiff);
		CHECK(maxDiff < 1e-10);
	};

	// Test vector 1: standard basis e_0.
	VectorComplex e0(d, ComplexType(0));
	e0[0] = ComplexType(1);
	checkPair(e0, "e_0");

	// Test vector 2: PsiHist_[0] from the completed propagation.
	const auto& hist = Acc::psiHist(s);
	if (!hist.empty() && hist[0].size() == d) {
		checkPair(hist[0], "PsiHist[0]");
	}
}

// ── TC-6: krylovExpmvCSR preserves norm ──────────────────────────────────────
// A non-Hermitian H causes norm drift; this catches wrong-sign JW entries
// or numerical errors in the Krylov recurrence.
TEST_CASE("krylovExpmvCSR preserves vector norm to 1e-10", "[GBEK][Krylov][norm]")
{
	auto&          f    = fixture();
	auto&          s    = f.solver;
	const SizeType r    = Acc::bathRank(s);
	auto           vMid = makeTestVMid(r);

	const auto&    upW    = Acc::upWordsNm1(s);
	const auto&    dnW    = Acc::dnWordsNm1(s);
	const SizeType d      = upW.size() * dnW.size();
	auto&          csr    = Acc::csrNm1Mut(s);
	const auto&    varVec = Acc::varNm1(s);

	const RealType dt = 0.1;

	Acc::updateCSR(s, csr, varVec, vMid);

	// Unit vector e_0.
	VectorComplex e0(d, ComplexType(0));
	e0[0] = ComplexType(1);

	auto           result  = Acc::krylovExpmvCSR(s, e0, Acc::csrNm1(s), dt);
	const RealType normOut = vecNorm(result);
	INFO("||krylovExpmvCSR(e_0)|| = " << normOut);
	CHECK(std::abs(normOut - 1.0) < 1e-10);

	// Also check N+1 sector to exercise both paths.
	const auto&    upWp = Acc::upWordsNp1(s);
	const auto&    dnWp = Acc::dnWordsNp1(s);
	const SizeType dp   = upWp.size() * dnWp.size();
	auto&          csrP = Acc::csrNp1Mut(s);
	const auto&    varP = Acc::varNp1(s);

	Acc::updateCSR(s, csrP, varP, vMid);

	VectorComplex ep(dp, ComplexType(0));
	ep[0] = ComplexType(1);

	auto           resultP = Acc::krylovExpmvCSR(s, ep, Acc::csrNp1(s), dt);
	const RealType normP   = vecNorm(resultP);
	INFO("N+1 sector: ||krylovExpmvCSR(e_0)|| = " << normP);
	CHECK(std::abs(normP - 1.0) < 1e-10);
}

// ── TC-7: Gα/Gβ spin-seed averaging (GBEK Eq. 70) ──────────────────────────────
//
// Regression test for a real bug (found and fixed 2026-07-08): a
// spin-imbalanced base filling (nup != ndown, e.g. the true atomic-limit
// single-atom seed nup=1, ndown=0) produces a spin-polarized extended-Fock
// ground state on its own. Without averaging Galpha (impurity's extra
// electron up) and Gbeta (down), the up-channel occupation is frozen at
// n=1 for all time instead of oscillating around the physically-required
// n=1/2 -- exactly the failure this test guards against. The expected
// value 0.5 is not a tolerance-tuned regression anchor: it follows directly
// from GBEK Sec. VI-A (paramagnetic, particle-hole symmetric state at half
// filling, <n_sigma>=1/2 for all t) and from the exact self-consistency
// relation Delta(t,t)=v(t)^2 n(t) matching the paper's own Fig. 3 colorbar
// maximum of 0.5.
//
// nBath=0 (true atomic limit, Delta^-=0 exactly) isolates this from the
// separate second-bath Cholesky machinery: with an empty bathParams vector,
// solveLplus builds nsites_ext=1+2L with no first-bath sites at all, so any
// deviation from n=1/2 can only come from the (nup_ext,ndown_ext) seeding,
// not from Cholesky-approximation error.
TEST_CASE("Gimp(t,t) plateaus at 1/2 for a spin-imbalanced atomic-limit seed",
          "[GBEK][spin-averaging]")
{
	static const std::string kAtomicConfig
	    = "##Ainur1.0\n\n"
	      "FicticiousBeta=10;\n"
	      "ChemicalPotential=0.;\n"
	      "Matsubaras=20;\n"
	      "LatticeGf=\"energy,semicircular,4\";\n"
	      "NumberOfBathPoints=1;\n"
	      "DmftNumberOfIterations=1;\n"
	      "DmftTolerance=1e-3;\n"
	      "ImpuritySolver=\"exactdiag\";\n"
	      "FitOptions=particleholesymmetric;\n"
	      "MinParamsDelta=0.01;\n"
	      "MinParamsDelta2=0.01;\n"
	      "MinParamsTolerance=1e-4;\n"
	      "MinParamsMaxIter=100;\n"
	      "MinParamsVerbose=0;\n"
	      "TargetElectronsUp=1;\n"
	      "TargetElectronsDown=0;\n"
	      "int ImpuritySite=0;\n"
	      "real HubbardU=2.;\n"
	      "HubbardUFinal=2.;\n"
	      "RootOutputname=\"testGBEKAtomic\";\n"
	      "InfiniteLoopKeptStates=20;\n"
	      "matrix FiniteLoopsGs=[[@auto, 20, 0],[@auto, 20, 0]];\n"
	      "real OmegaBegin=-4.;\n"
	      "integer OmegaTotal=40;\n"
	      "real OmegaStep=0.2;\n"
	      "real OmegaDelta=0.2;\n"
	      "integer TridiagSteps=100;\n"
	      "real TridiagEps=1e-6;\n"
	      "TruncationTolerance=\"1e-6,20\";\n"
	      "CorrectionVectorEta=0.;\n"
	      "GsWeight=0.1;\n"
	      "matrix FiniteLoopsOmega=[[@auto, 20, 2],[@auto, 20, 2]];\n"
	      "TmaxNeq=0.2;\n"
	      "NtNeq=4;\n"
	      "NeqDmftIter=1;\n"
	      "NeqDmftTolerance=1e-4;\n"
	      "NeqSolver=\"gbek\";\n"
	      "NeqBathRank=2;\n"
	      "BandwidthFinal=4.;\n";

	InputNgType::Writeable ioW(Dmft::CincuentaInputCheck {}, kAtomicConfig);
	InputNgType::Readable  io(ioW);
	ParamsType             params(io);
	SolverType             solver(params, io);

	solver.solve({}); // empty bathParams: nBath=0, true atomic limit

	using KBType = SolverType::KBType;
	KBType slice(params.nT,
	             params.eqParams.nMatsubaras,
	             params.dt,
	             params.eqParams.ficticiousBeta
	                 / static_cast<RealType>(params.eqParams.nMatsubaras));

	for (int n = 0; n <= static_cast<int>(params.nT); ++n) {
		solver.computeGimp(slice, n);
		const ComplexType g
		    = slice.lesser(static_cast<SizeType>(n), static_cast<SizeType>(n));
		INFO("n=" << n << " G_imp(t,t)=" << g);
		CHECK(g.real() == Catch::Approx(0.0).margin(1e-10));
		CHECK(g.imag() == Catch::Approx(0.5).epsilon(1e-8));
	}
}
