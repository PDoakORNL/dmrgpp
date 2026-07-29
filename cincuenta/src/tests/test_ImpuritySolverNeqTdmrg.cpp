#include "CincuentaInputCheck.h"
#include "ImpuritySolverNeqExactDiag.h"
#include "ImpuritySolverNeqTdmrg.h"
#include "NeqDmftSolver.h"
#include "NeqLatticeGf.h"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <complex>
#include <iostream>
#include <string>

// Phase 1 (evolving-bath project) diagnostic test: verifies that chaining
// nT single-step restarts (ImpuritySolverNeqTdmrg::solveChainedColumn0)
// reproduces the SAME column-0 result as today's single monolithic
// Run3/Run5 (ImpuritySolverNeqTdmrg::solve()), before any new-column
// machinery is added. See /Users/epd/.claude/plans/fancy-painting-moon.md.
//
// Two variants, same geometry/bath, differing only in U:
//   U=0    -- matches the existing inputU0NeqTdmrg.ain regression input's
//             own rationale (tDMRG is exact here, no entanglement growth),
//             BUT this is also exactly the point flagged in project memory
//             (output3 test) as having DMRG degenerate-eigenvector
//             arch-dependence at U=0. A residual chained-vs-monolithic
//             mismatch at U=0 is therefore ambiguous: it could be a real
//             chaining/gauge bug, or it could be two independent DMRG runs
//             resolving a U=0 degeneracy differently.
//   U=0.5  -- small nonzero U lifts that degeneracy while keeping
//             truncation error small, disambiguating the two explanations.

using RealType          = double;
using ComplexType       = std::complex<RealType>;
using SolverType        = Dmft::ImpuritySolverNeqTdmrg<ComplexType>;
using ParamsType        = Dmft::ParamsNeqDmftSolver<ComplexType>;
using InputNgType       = PsimagLite::InputNg<Dmft::CincuentaInputCheck>;
using VectorRealType    = typename SolverType::VectorRealType;
using VectorComplexType = typename SolverType::VectorComplexType;
using ApplicationType   = PsimagLite::PsiApp;

// N=6 sites (5 bath + 1 impurity): TSPAdvanceEach=N-2=4, matching the
// convention already established in inputU0NeqTdmrg.ain. NtNeq kept tiny
// (2) to keep this fast -- this is a mechanics check, not a physics
// benchmark. %U% is substituted with the desired HubbardU/HubbardUFinal
// value (no quench in U either variant -- only U itself is varied between
// the two test cases, everything else identical).
// rootName (default "testTdmrgChain", matching every existing call site):
// pass a distinct bare-alphanumeric identifier for a TEST_CASE that would
// otherwise collide with another sharing the same RootOutputname-derived
// checkpoint file prefix -- confirmed as a REAL, reproducible cross-
// TEST_CASE contamination mechanism (not hypothetical). Must be bare
// alphanumeric: a period in the root name (e.g. embedding a raw "uValue"
// like "0.5") breaks log/checkpoint filename parsing outright (confirmed
// empirically -- do not parameterize this by uValue).
static std::string configWithU(const std::string& uValue,
                               const std::string& rootName = "testTdmrgChain")
{
	return "##Ainur1.0\n\n"
	       "FicticiousBeta=20;\n"
	       "ChemicalPotential=0.;\n"
	       "Matsubaras=200;\n"
	       "LatticeGf=\"energy,semicircular,4\";\n"
	       "NumberOfBathPoints=5;\n"
	       "DmftNumberOfIterations=1;\n"
	       "DmftTolerance=1e-6;\n"
	       "ImpuritySolver=\"exactdiag\";\n"
	       "FitOptions=particleholesymmetric;\n"
	       "MinParamsDelta=0.01;\n"
	       "MinParamsMaxIter=10000;\n"
	       "MinParamsDelta2=0.01;\n"
	       "MinParamsTolerance=1e-4;\n"
	       "MinParamsVerbose=0;\n"
	       "vector InitBathVector=[0.447, 0.447, 0.447, 1.0, 1.8];\n"
	       "int ImpuritySite=0;\n"
	       "real HubbardU="
	    + uValue
	    + ";\n"
	      "TargetElectronsUp=3;\n"
	      "TargetElectronsDown=3;\n"
	      "RootOutputname=\""
	    + rootName
	    + "\";\n"
	      "InfiniteLoopKeptStates=100;\n"
	      "matrix FiniteLoopsGs=[[@auto, 100, 0],[@auto, 100, 0]];\n"
	      "real OmegaBegin=-6.;\n"
	      "integer OmegaTotal=20;\n"
	      "real OmegaStep=0.3;\n"
	      "real OmegaDelta=0.1;\n"
	      "integer TridiagSteps=200;\n"
	      "real TridiagEps=1e-9;\n"
	      "TruncationTolerance=\"1e-10,100\";\n"
	      "CorrectionVectorEta=0.;\n"
	      "GsWeight=0.1;\n"
	      "matrix FiniteLoopsOmega=[[@auto, 100, 2],[@auto, 100, 2]];\n"
	      "HubbardUFinal="
	    + uValue
	    + ";\n"
	      "TmaxNeq=0.2;\n"
	      "NtNeq=2;\n"
	      "NeqDmftIter=1;\n"
	      "NeqDmftTolerance=0.001;\n"
	      "NeqSolver=\"tdmrg\";\n"
	      "matrix FiniteLoopsTdmrg=[\n"
	      "    [@auto, 100, 0],[@auto, 100, 0],\n"
	      "    [@auto, 100, 0],[@auto, 100, 0]];\n"
	      "TSPTimeSteps=5;\n"
	      "TSPAdvanceEach=4;\n";
}

static void runChainedVsMonolithicComparison(const std::string& uValue)
{
	int    argc    = 1;
	char   arg0[]  = "test_ImpuritySolverNeqTdmrg";
	char*  argv0[] = { arg0 };
	char** argv    = argv0;

	ApplicationType app("test_ImpuritySolverNeqTdmrg", &argc, &argv, 1);

	InputNgType::Writeable ioW(Dmft::CincuentaInputCheck {}, configWithU(uValue));
	InputNgType::Readable  io(ioW);
	ParamsType             params(io);
	SolverType             solver(params, app, io);

	// Symmetric 5-site bath, arbitrary non-degenerate values.
	const VectorRealType bathParams = { 0.3, 0.3, 0.3, 0.3, 0.3, -0.6, -0.3, 0.0, 0.3, 0.6 };

	solver.solve(bathParams);
	const auto& monolithic = solver.gimp();

	const auto chained = solver.solveChainedColumn0(bathParams);

	std::cout << "=== Phase 1 chained-vs-monolithic column-0 comparison (U=" << uValue
	          << ") ===\n";
	for (int n = 1; n <= static_cast<int>(params.nT); ++n) {
		const ComplexType monoRet = monolithic.retarded(n, 0);
		const ComplexType monoLes = monolithic.lesser(n, 0);

		auto itG = chained.ggt0.find(n);
		auto itL = chained.glt0.find(n);
		REQUIRE(itG != chained.ggt0.end());
		REQUIRE(itL != chained.glt0.end());

		const ComplexType chainedGgt = itG->second;
		const ComplexType chainedGlt = itL->second;
		const ComplexType chainedRet = chainedGgt - chainedGlt;

		std::cout << "n=" << n << " monolithic G^R(n,0)=" << monoRet
		          << " monolithic G^<(n,0)=" << monoLes
		          << " chained G^>(n,0)=" << chainedGgt
		          << " chained G^<(n,0)=" << chainedGlt
		          << " chained G^R(n,0)=" << chainedRet << "\n";

		// Tight tolerance: at U=0/0.5, TmaxNeq=0.2, NtNeq=2 this is a small
		// enough system that truncation error is negligible (m=100 is not
		// even close to saturated) -- any real disagreement here is a
		// chaining/gauge bug, not physics, so this should not be loosened.
		const RealType tol = 1e-6;
		CHECK(chainedRet.real() == Catch::Approx(monoRet.real()).margin(tol));
		CHECK(chainedRet.imag() == Catch::Approx(monoRet.imag()).margin(tol));
		CHECK(chainedGlt.real() == Catch::Approx(monoLes.real()).margin(tol));
		CHECK(chainedGlt.imag() == Catch::Approx(monoLes.imag()).margin(tol));
	}
}

TEST_CASE("ImpuritySolverNeqTdmrg chained column-0 matches monolithic solve() (U=0)",
          "[ImpuritySolverNeqTdmrg][Phase1]")
{
	runChainedVsMonolithicComparison("0.");
}

TEST_CASE("ImpuritySolverNeqTdmrg chained column-0 matches monolithic solve() (U=0.5, "
          "no U=0 degeneracy)",
          "[ImpuritySolverNeqTdmrg][Phase1]")
{
	runChainedVsMonolithicComparison("0.5");
}

// Phase 1 full-grid gate: computeFullGrid() must reproduce a genuine
// two-time G(t_n,t_j) grid, not just column 0. Compared against
// ImpuritySolverNeqExactDiag's own full-grid computeGimp (self-contained,
// no NeqDmftSolver driving needed -- see plan file) for all (n,j), j<n.
//
// The atomic-limit closed form (nBath=0) is NOT used here as the
// reference, unlike the plan's general suggestion: tDMRG's star geometry
// needs at least one real bath site (TSPAdvanceEach=nsites-2 is degenerate
// at nsites=1), so nBath=0 isn't a configuration tDMRG can actually run.
// ImpuritySolverNeqExactDiag's general (nonzero-bath) full grid is the
// right reference for this solver.
//
// Diagonal (n==j) IS included: computeFullGrid now fills it via a
// gauge-invariant equal-time measurement taken at each column's birth
// (see Column::ggtDiag/gltDiag).
static void runFullGridVsExactDiagComparison(const std::string& uValue)
{
	int    argc    = 1;
	char   arg0[]  = "test_ImpuritySolverNeqTdmrg";
	char*  argv0[] = { arg0 };
	char** argv    = argv0;

	ApplicationType app("test_ImpuritySolverNeqTdmrg", &argc, &argv, 1);

	InputNgType::Writeable ioW(Dmft::CincuentaInputCheck {}, configWithU(uValue));
	InputNgType::Readable  io(ioW);
	ParamsType             params(io);
	SolverType             tdmrg(params, app, io);

	using ExactDiagType = Dmft::ImpuritySolverNeqExactDiag<ComplexType>;
	using EdKBType      = typename ExactDiagType::KBType;
	InputNgType::Writeable ioW2(Dmft::CincuentaInputCheck {}, configWithU(uValue));
	InputNgType::Readable  io2(ioW2);
	ParamsType             paramsEd(io2);
	ExactDiagType          exactDiag(paramsEd, io2);

	// Symmetric 5-site bath, arbitrary non-degenerate values (same as the
	// column-0-only test above).
	const VectorRealType bathParams = { 0.3, 0.3, 0.3, 0.3, 0.3, -0.6, -0.3, 0.0, 0.3, 0.6 };

	exactDiag.solve(bathParams);

	// computeGimp fills a CALLER-owned grid (see NeqDmftSolver's own usage
	// pattern) -- it does not populate the solver's internal gimp_ (that's
	// a separate accessor, used only for Matsubara/left-mixing components
	// elsewhere). Build our own reference grid, same idiom as
	// test_ImpuritySolverNeqExactDiag.cpp's makeSlice.
	EdKBType reference(paramsEd.nT,
	                   paramsEd.eqParams.nMatsubaras,
	                   paramsEd.dt,
	                   paramsEd.eqParams.ficticiousBeta
	                       / static_cast<RealType>(paramsEd.eqParams.nMatsubaras));
	for (int n = 0; n <= static_cast<int>(paramsEd.nT); ++n)
		exactDiag.computeGimp(reference, n);

	const auto grid = tdmrg.computeFullGrid(bathParams);

	std::cout << "=== Phase 1 full-grid vs ExactDiag comparison (U=" << uValue << ") ===\n";
	const int      nT  = static_cast<int>(params.nT);
	const RealType tol = 1e-6;
	for (int n = 0; n <= nT; ++n) {
		for (int j = 0; j <= n; ++j) {
			const ComplexType refRet  = reference.retarded(n, j);
			const ComplexType refLes  = reference.lesser(n, j);
			const ComplexType gridRet = grid.retarded(n, j);
			const ComplexType gridLes = grid.lesser(n, j);

			std::cout << "n=" << n << " j=" << j << " ExactDiag G^R=" << refRet
			          << " G^<=" << refLes << " tDMRG G^R=" << gridRet
			          << " G^<=" << gridLes << "\n";

			CHECK(gridRet.real() == Catch::Approx(refRet.real()).margin(tol));
			CHECK(gridRet.imag() == Catch::Approx(refRet.imag()).margin(tol));
			CHECK(gridLes.real() == Catch::Approx(refLes.real()).margin(tol));
			CHECK(gridLes.imag() == Catch::Approx(refLes.imag()).margin(tol));
		}
	}
}

TEST_CASE("ImpuritySolverNeqTdmrg full grid matches ExactDiag full grid (U=0)",
          "[ImpuritySolverNeqTdmrg][Phase1][FullGrid]")
{
	runFullGridVsExactDiagComparison("0.");
}

TEST_CASE("ImpuritySolverNeqTdmrg full grid matches ExactDiag full grid (U=0.5)",
          "[ImpuritySolverNeqTdmrg][Phase1][FullGrid]")
{
	runFullGridVsExactDiagComparison("0.5");
}

// Phase 2 blocker-B gate: measureSecondBathOccupations must correctly seed
// L "occupied" (<n>~2) and L "empty" (<n>~0) second-bath sites via the
// eps-split GS potential, PROVIDED eps is scaled to the largest other
// coupling present (here, the single first-bath hopping V=0.5). A single
// first-bath site (nBath=1) is used, matching the standalone smoke test in
// build/tmp/eps_split_seed/ that first found this requirement empirically
// -- see project_tdmrg_evolving_bath memory / fancy-painting-moon.md plan.
static std::string configForSecondBathTest()
{
	return "##Ainur1.0\n\n"
	       "FicticiousBeta=20;\n"
	       "ChemicalPotential=0.;\n"
	       "Matsubaras=200;\n"
	       "LatticeGf=\"energy,semicircular,4\";\n"
	       "NumberOfBathPoints=1;\n"
	       "DmftNumberOfIterations=1;\n"
	       "DmftTolerance=1e-6;\n"
	       "ImpuritySolver=\"exactdiag\";\n"
	       "FitOptions=particleholesymmetric;\n"
	       "MinParamsDelta=0.01;\n"
	       "MinParamsMaxIter=10000;\n"
	       "MinParamsDelta2=0.01;\n"
	       "MinParamsTolerance=1e-4;\n"
	       "MinParamsVerbose=0;\n"
	       "vector InitBathVector=[0.5];\n"
	       "int ImpuritySite=0;\n"
	       "real HubbardU=0.;\n"
	       "TargetElectronsUp=1;\n"
	       "TargetElectronsDown=0;\n"
	       "RootOutputname=\"testTdmrgSecondBathSeed\";\n"
	       "InfiniteLoopKeptStates=60;\n"
	       "matrix FiniteLoopsGs=[[@auto, 60, 0],[@auto, 60, 0]];\n"
	       "real OmegaBegin=-6.;\n"
	       "integer OmegaTotal=20;\n"
	       "real OmegaStep=0.3;\n"
	       "real OmegaDelta=0.1;\n"
	       "integer TridiagSteps=200;\n"
	       "real TridiagEps=1e-9;\n"
	       "TruncationTolerance=\"1e-10,100\";\n"
	       "CorrectionVectorEta=0.;\n"
	       "GsWeight=0.1;\n"
	       "matrix FiniteLoopsOmega=[[@auto, 60, 2],[@auto, 60, 2]];\n"
	       "HubbardUFinal=0.;\n"
	       "TmaxNeq=0.2;\n"
	       "NtNeq=2;\n"
	       "NeqDmftIter=1;\n"
	       "NeqDmftTolerance=0.001;\n"
	       "NeqSolver=\"tdmrg\";\n"
	       "matrix FiniteLoopsTdmrg=[\n"
	       "    [@auto, 60, 0],[@auto, 60, 0],\n"
	       "    [@auto, 60, 0],[@auto, 60, 0]];\n"
	       "TSPTimeSteps=5;\n"
	       "TSPAdvanceEach=4;\n";
}

TEST_CASE("ImpuritySolverNeqTdmrg second-bath eps-split seeding requires eps "
          "scaled to the coupling strength",
          "[ImpuritySolverNeqTdmrg][Phase2][SecondBathSeeding]")
{
	int    argc    = 1;
	char   arg0[]  = "test_ImpuritySolverNeqTdmrg";
	char*  argv0[] = { arg0 };
	char** argv    = argv0;

	ApplicationType app("test_ImpuritySolverNeqTdmrg", &argc, &argv, 1);

	InputNgType::Writeable ioW(Dmft::CincuentaInputCheck {}, configForSecondBathTest());
	InputNgType::Readable  io(ioW);
	ParamsType             params(io);
	SolverType             solver(params, app, io);

	// nBath=1, hopping V=0.5, bathEps=0.
	const VectorRealType bathParams = { 0.5, 0.0 };
	const SizeType       L          = 1;

	SECTION("eps comparable to V: correct L-occupied/L-empty seeding")
	{
		const RealType eps = 5.0 * 0.5; // several times the first-bath hopping V
		const auto     occ = solver.measureSecondBathOccupations(bathParams, L, eps);
		REQUIRE(occ.size() == 2 * L);
		// occ[0] = intended-occupied site, occ[1] = intended-empty site.
		CHECK(occ[0] == Catch::Approx(2.0).margin(1e-4));
		CHECK(occ[1] == Catch::Approx(0.0).margin(1e-4));
	}

	SECTION("eps too small relative to V: seeding fails (documents the "
	        "calibration requirement, not a bug)")
	{
		const RealType eps = 0.01; // << V=0.5, matches the failing smoke-test case
		const auto     occ = solver.measureSecondBathOccupations(bathParams, L, eps);
		REQUIRE(occ.size() == 2 * L);
		// Does NOT match the intended 2/0 pattern -- this is the whole
		// point of the test: a fixed small eps is not safe in general.
		const bool seededCorrectly
		    = (std::abs(occ[0] - 2.0) < 1e-4) && (std::abs(occ[1] - 0.0) < 1e-4);
		CHECK_FALSE(seededCorrectly);
	}
}

// Phase 2 fan-out gate (advisor-recommended first check, per
// fancy-painting-moon.md): extend the geometry with 2L second-bath sites
// via computeFullGridWithInertSecondBath, holding Vplus=0 for EVERY step
// (not just n=0) -- an "inert spectator" test of the larger lattice, the
// +L electron counts, the eps-split GS, and the recalibrated
// TSPAdvanceEach, all exercised together, BEFORE any real second-bath
// physics or self-consistency is added. With Vplus=0 throughout, the 2L
// sites never couple to anything, so every (n,j) must still match today's
// UNEXTENDED computeFullGrid to the same 1e-6 tolerance as the Phase 1
// gate above -- any mismatch here means the fan-out mechanics broke on
// the larger geometry, not a physics difference.
TEST_CASE("ImpuritySolverNeqTdmrg extended-geometry (inert L=1 second bath) "
          "matches unextended full grid",
          "[ImpuritySolverNeqTdmrg][Phase2][FanOut]")
{
	int    argc    = 1;
	char   arg0[]  = "test_ImpuritySolverNeqTdmrg";
	char*  argv0[] = { arg0 };
	char** argv    = argv0;

	ApplicationType app("test_ImpuritySolverNeqTdmrg", &argc, &argv, 1);

	InputNgType::Writeable ioW(Dmft::CincuentaInputCheck {}, configWithU("0.5"));
	InputNgType::Readable  io(ioW);
	ParamsType             params(io);
	SolverType             solver(params, app, io);

	// Same symmetric 5-site bath used by the Phase 1 full-grid gate.
	const VectorRealType bathParams = { 0.3, 0.3, 0.3, 0.3, 0.3, -0.6, -0.3, 0.0, 0.3, 0.6 };
	const SizeType       L          = 1;
	// eps well above the largest coupling actually present (max|hopping|,
	// max|bathEps| here are 0.3 and 0.6) -- per the blocker-B calibration
	// rule, not a fixed small constant.
	const RealType eps = 3.0;

	const auto reference = solver.computeFullGrid(bathParams);
	const auto extended  = solver.computeFullGridWithInertSecondBath(bathParams, L, eps);

	const int      nT  = static_cast<int>(params.nT);
	const RealType tol = 1e-6;
	for (int n = 0; n <= nT; ++n) {
		for (int j = 0; j <= n; ++j) {
			const ComplexType refRet = reference.retarded(n, j);
			const ComplexType refLes = reference.lesser(n, j);
			const ComplexType extRet = extended.retarded(n, j);
			const ComplexType extLes = extended.lesser(n, j);

			CHECK(extRet.real() == Catch::Approx(refRet.real()).margin(tol));
			CHECK(extRet.imag() == Catch::Approx(refRet.imag()).margin(tol));
			CHECK(extLes.real() == Catch::Approx(refLes.real()).margin(tol));
			CHECK(extLes.imag() == Catch::Approx(refLes.imag()).margin(tol));
		}
	}
}

// ---- Phase 2: diagonal Connectors-independence check -----------------------
//
// Load-bearing claim for the self-consistent (evolving-Vplus) design: a
// column's equal-time diagonal G(born,born) is captured from the FIRST
// measurement occurrence in its first advanceColumn segment, which (per the
// advisor consult recorded in fancy-painting-moon.md) should be taken BEFORE
// that segment's own TimeEvolve acts -- i.e. before the segment's declared
// second-bath Connectors value has any chance to matter. This was previously
// confirmed only at Connectors=0 (computeFullGridWithInertSecondBath); this
// test exercises it with two DIFFERENT nonzero Connectors values and checks
// the diagonal is identical while the (necessarily Connectors-dependent)
// off-diagonal entry differs. If this fails, the byproduct-diagonal-capture
// design is unsafe for genuinely evolving Vplus and the fallback (measuring
// <n_imp> directly, per measureSecondBathOccupations) must be used instead.
TEST_CASE("ImpuritySolverNeqTdmrg column diagonal is independent of the "
          "second-bath Connectors value used for its first advance",
          "[ImpuritySolverNeqTdmrg][Phase2][ConnectorsIndependence]")
{
	int    argc    = 1;
	char   arg0[]  = "test_ImpuritySolverNeqTdmrg";
	char*  argv0[] = { arg0 };
	char** argv    = argv0;

	ApplicationType app("test_ImpuritySolverNeqTdmrg", &argc, &argv, 1);

	InputNgType::Writeable ioW(Dmft::CincuentaInputCheck {}, configWithU("0.5"));
	InputNgType::Readable  io(ioW);
	ParamsType             params(io);
	SolverType             solver(params, app, io);

	const VectorRealType bathParams = { 0.3, 0.3, 0.3, 0.3, 0.3, -0.6, -0.3, 0.0, 0.3, 0.6 };
	const SizeType       L          = 1;
	const RealType       eps        = 3.0;

	// Two distinct nonzero Connectors profiles for the 2L=2 second-bath
	// entries (occupied-site coupling, empty-site coupling), one real, one
	// complex -- deliberately different from each other and from the inert
	// (0) case already covered above.
	const VectorComplexType connectorsA = { ComplexType(0.3, 0.0), ComplexType(0.3, 0.0) };
	const VectorComplexType connectorsB = { ComplexType(0.7, 0.2), ComplexType(0.7, 0.2) };

	// Distinct rootSuffix per call: this test makes two back-to-back calls
	// against the same solver/root_, which otherwise write/read the same
	// checkpoint files -- observed to cause real cross-call contamination
	// when run as part of the full suite (see computeFullGridWithInertSecondBath's
	// doc comment).
	const auto gridA = solver.computeFullGridWithInertSecondBath(
	    bathParams, L, eps, connectorsA, "gridbathA_");
	const auto gridB = solver.computeFullGridWithInertSecondBath(
	    bathParams, L, eps, connectorsB, "gridbathB_");

	const RealType tol = 1e-6;

	// Diagonal entries: must match between the two runs (Connectors-
	// independent, per the hypothesis under test).
	const int nT = static_cast<int>(params.nT);
	for (int n = 0; n <= nT; ++n) {
		const ComplexType retA = gridA.retarded(n, n);
		const ComplexType retB = gridB.retarded(n, n);
		const ComplexType lesA = gridA.lesser(n, n);
		const ComplexType lesB = gridB.lesser(n, n);
		CHECK(retA.real() == Catch::Approx(retB.real()).margin(tol));
		CHECK(retA.imag() == Catch::Approx(retB.imag()).margin(tol));
		CHECK(lesA.real() == Catch::Approx(lesB.real()).margin(tol));
		CHECK(lesA.imag() == Catch::Approx(lesB.imag()).margin(tol));
	}

	// Off-diagonal (1,0): genuinely evolved under the segment's own
	// Connectors, so should DIFFER between the two profiles -- a sanity
	// check that this test isn't vacuously passing because everything is
	// insensitive to Connectors.
	if (nT >= 1) {
		const ComplexType offA = gridA.retarded(1, 0);
		const ComplexType offB = gridB.retarded(1, 0);
		const RealType    diff = std::abs(offA - offB);
		CHECK(diff > tol);
	}
}

// ---- Phase 2: NeqBathRank>0 wiring runs end-to-end -------------------------
//
// ImpuritySolverNeqTdmrg's constructor takes an extra ApplicationType&
// argument that NeqDmftSolver's generic impSolver_(params, io) construction
// call cannot supply -- confirmed this is exactly why cincuenta.cpp's own
// "tdmrg" branch never wraps ImpuritySolverNeqTdmrg in NeqDmftSolver at all
// (it drives ImpuritySolverNeqTdmrg directly). So this test cannot use
// NeqDmftSolver itself; it manually replicates NeqDmftSolver::solve()/
// timeStep()'s exact sequence (see NeqDmftSolver.h) using NeqLatticeGf
// directly, which needs no such extra argument. This is the correctness
// gate for Task 14 (solve()/computeGimp()/prepareTimeStep() dispatch for
// neqBathRank_>0) -- NOT the Phase 2 "full gate" against
// ImpuritySolverNeqGBEK (Task 15), which additionally requires resolving
// this constructor-signature mismatch (or writing an equally manual drive
// for GBEK too) and is deliberately left for that separate, later task.
TEST_CASE("ImpuritySolverNeqTdmrg NeqBathRank=1 self-consistent wiring runs "
          "end-to-end without crashing",
          "[ImpuritySolverNeqTdmrg][Phase2][SelfConsistentWiring]")
{
	int    argc    = 1;
	char   arg0[]  = "test_ImpuritySolverNeqTdmrg";
	char*  argv0[] = { arg0 };
	char** argv    = argv0;

	ApplicationType app("test_ImpuritySolverNeqTdmrg", &argc, &argv, 1);

	InputNgType::Writeable ioW(Dmft::CincuentaInputCheck {},
	                           configWithU("0.5", "testTdmrgChainSCW") + "NeqBathRank=1;\n");
	InputNgType::Readable  io(ioW);
	ParamsType             params(io);
	SolverType             solver(params, app, io);

	const VectorRealType bathParams = { 0.3, 0.3, 0.3, 0.3, 0.3, -0.6, -0.3, 0.0, 0.3, 0.6 };

	Dmft::NeqLatticeGf<ComplexType> latticeGf(params);
	Dmft::KadanoffBaym<ComplexType> gimp(
	    params.nT,
	    params.eqParams.nMatsubaras,
	    params.dt,
	    params.eqParams.ficticiousBeta / static_cast<RealType>(params.eqParams.nMatsubaras));

	solver.solve(bathParams);
	gimp.matsubara_t = solver.gimp().matsubara_t;
	gimp.matsubara_w = solver.gimp().matsubara_w;

	solver.computeGimp(gimp, 0);
	latticeGf.initialize(gimp);
	latticeGf.updateDelta(0, gimp);
	solver.prepareTimeStep(0, latticeGf.delta());

	const int nT = static_cast<int>(params.nT);
	for (int n = 1; n <= nT; ++n) {
		solver.computeGimp(gimp, n); // predictor
		for (SizeType iter = 0; iter < params.neqDmftIter; ++iter) {
			latticeGf.updateDelta(n, gimp);
			solver.prepareTimeStep(n, latticeGf.delta());
			solver.computeGimp(gimp, n); // corrector
		}
		latticeGf.advance(n);

		for (int j = 0; j <= n; ++j) {
			CHECK(std::isfinite(gimp.retarded(n, j).real()));
			CHECK(std::isfinite(gimp.retarded(n, j).imag()));
			CHECK(std::isfinite(gimp.lesser(n, j).real()));
			CHECK(std::isfinite(gimp.lesser(n, j).imag()));
		}
	}
}

// ---- Task 16: NeqDmftSolver<..., ImpuritySolverNeqTdmrg> constructibility --
//
// Confirms the new 3-arg NeqDmftSolver constructor overload (NeqDmftSolver.h)
// actually lets ImpuritySolverNeqTdmrg be driven through NeqDmftSolver's OWN
// solve()/timeStep() loop, exactly like ImpuritySolverNeqGBEK/ExactDiag/
// Lanczos already are -- superseding the manual replication in the
// "SelfConsistentWiring" test above (kept as-is; it still documents/
// exercises the underlying computeGimp/prepareTimeStep contract directly).
// This is what Task 15's full gate needs: both tDMRG and GBEK driven
// identically through NeqDmftSolver, not via two different hand-rolled
// harnesses.
TEST_CASE("ImpuritySolverNeqTdmrg can be driven through NeqDmftSolver's own "
          "solve()",
          "[ImpuritySolverNeqTdmrg][Phase2][NeqDmftSolverWiring]")
{
	int    argc    = 1;
	char   arg0[]  = "test_ImpuritySolverNeqTdmrg";
	char*  argv0[] = { arg0 };
	char** argv    = argv0;

	ApplicationType app("test_ImpuritySolverNeqTdmrg", &argc, &argv, 1);

	InputNgType::Writeable ioW(Dmft::CincuentaInputCheck {},
	                           configWithU("0.5", "testTdmrgChainNDS") + "NeqBathRank=1;\n");
	InputNgType::Readable  io(ioW);
	ParamsType             params(io);

	using TdmrgNeqSolverType = Dmft::NeqDmftSolver<ComplexType, Dmft::ImpuritySolverNeqTdmrg>;
	TdmrgNeqSolverType neqSolver(params, app, io);

	const VectorRealType bathParams = { 0.3, 0.3, 0.3, 0.3, 0.3, -0.6, -0.3, 0.0, 0.3, 0.6 };
	neqSolver.solve(bathParams);

	const auto& gimp = neqSolver.gimp();
	const int   nT   = static_cast<int>(params.nT);
	for (int n = 0; n <= nT; ++n) {
		for (int j = 0; j <= n; ++j) {
			CHECK(std::isfinite(gimp.retarded(n, j).real()));
			CHECK(std::isfinite(gimp.retarded(n, j).imag()));
			CHECK(std::isfinite(gimp.lesser(n, j).real()));
			CHECK(std::isfinite(gimp.lesser(n, j).imag()));
		}
	}
}
