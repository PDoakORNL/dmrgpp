#include "CincuentaInputCheck.h"
#include "ImpuritySolverNeqExactDiag.h"
#include "ImpuritySolverNeqGBEK.h"
#include "ImpuritySolverNeqTdmrg.h"
#include "NeqDmftSolver.h"
#include "NeqLatticeGf.h"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
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
static std::string configForSecondBathTest(const std::string& rootName = "testTdmrgSecondBathSeed")
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
	       "RootOutputname=\""
	    + rootName
	    + "\";\n"
	      "InfiniteLoopKeptStates=60;\n"
	      "matrix FiniteLoopsGs=[[@auto, 60, 0],[@auto, 60, 0],[@auto, 60, "
	      "0],[@auto, 60, 0]];\n"
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

// Same eps-split GS seeding mechanism as the L=1 test directly above, but at
// higher rank (L=3,4,5) -- item #2 of the fast-isolated-test set added
// 2026-08-05 to reduce reliance on slow multi-step self-consistent runs (see
// TDMRG_EVOLVING_BATH.md). Each case is a SINGLE throwaway GS run (no time
// evolution, no self-consistency loop, no NeqDmftSolver): the same
// nBath=1/V=0.5 geometry as configForSecondBathTest(), extended by
// measureSecondBathOccupations to nsitesExt = 2 + 2L sites, with eps = 5*V as
// established by the L=1 case's calibration rule. Confirms the intended
// product state -- L occupied auxiliary sites (<n>~2), L empty auxiliary
// sites (<n>~0), impurity singly-occupied -- is reached at each rank, not
// just L=1.
//
// FINDING (2026-08-05, use CHECK not REQUIRE below so all three L still run
// and report): this is NOT rank-independent. L=3 passes with the original
// FiniteLoopsGs (2 sweeps). L=4 needed FiniteLoopsGs bumped to 4 sweeps here
// to converge -- with only 2 it also failed. L=5, even with 4 sweeps AND eps
// pushed up to 60*V (30, vs the L=1 rule's 2.5), still fails: one electron
// migrates from the impurity/first-bath sector into the (nominally fully
// decoupled, since second-bath Connectors=0 at the GS stage) auxiliary
// manifold, and specific occ/empty sites end up swapped rather than merely
// partially mixed. This is a genuine rank-dependent breakdown of the
// eps-split seeding mechanism, not a test artifact -- more sweeps did NOT
// fix L=5, so it is not simply "needs more DMRG convergence effort" the way
// L=4 was. Plausibly related to the parallel NtNeq>=NeqBathRank/DM-eigs
// investigation's underlying concern (many near-degenerate decoupled
// orbitals at higher L), though the failure mode observed here (occupation
// migrating to the WRONG site, not a degenerate-DM truncation crash) is
// different enough that it should be treated as a separate, adjacent
// finding rather than assumed to be the same bug. Not investigated further
// here (see feedback_compute_discipline: don't chase an open physics
// question with more compute without a plan).
TEST_CASE("ImpuritySolverNeqTdmrg second-bath eps-split seeding holds at "
          "higher rank (L=3,4,5)",
          "[ImpuritySolverNeqTdmrg][Phase2][SecondBathSeeding]")
{
	int    argc    = 1;
	char   arg0[]  = "test_ImpuritySolverNeqTdmrg";
	char*  argv0[] = { arg0 };
	char** argv    = argv0;

	ApplicationType app("test_ImpuritySolverNeqTdmrg", &argc, &argv, 1);

	// GENERATE must come BEFORE building the config/root name: each L gets
	// its OWN RootOutputname (and measureSecondBathOccupations further
	// appends "secondbath_seed_gs" to it, but the geometry -- nsitesExt --
	// differs per L, so the prefix itself must differ too). Reusing one
	// prefix across differing lattice sizes is exactly the cross-TEST_CASE
	// checkpoint contamination mechanism documented above
	// configWithU()/runChainedVsMonolithicComparison -- confirmed real, not
	// hypothetical. Root names stay bare alphanumeric per that same
	// constraint (no periods).
	const SizeType L = GENERATE(SizeType(3), SizeType(4), SizeType(5));
	CAPTURE(L);
	const std::string rootName = "testTdmrgSecondBathSeedL" + ttos(L);

	InputNgType::Writeable ioW(Dmft::CincuentaInputCheck {}, configForSecondBathTest(rootName));
	InputNgType::Readable  io(ioW);
	ParamsType             params(io);
	SolverType             solver(params, app, io);

	// nBath=1, hopping V=0.5, bathEps=0 -- identical geometry to the L=1 case.
	const VectorRealType bathParams = { 0.5, 0.0 };
	const RealType       eps        = 5.0 * 0.5; // same calibration rule as the L=1 case

	const auto occ = solver.measureSecondBathOccupations(bathParams, L, eps);
	REQUIRE(occ.size() == 2 * L);
	RealType auxTotal = 0;
	for (SizeType p = 0; p < L; ++p) {
		CHECK(occ[p] == Catch::Approx(2.0).margin(1e-4));
		CHECK(occ[L + p] == Catch::Approx(0.0).margin(1e-4));
		auxTotal += occ[p] + occ[L + p];
	}
	// TargetElectronsUp=1, TargetElectronsDown=0 (base impurity sector) plus
	// L extra up- and L extra down-electrons for the eps-split (nupExt=L+1,
	// ndownExt=L, total 2L+1) -- if the 2L auxiliary sites correctly hold
	// 2L electrons between them (checked above), the remaining single
	// electron is confined to the impurity+first-bath sector's TOTAL
	// occupation. That is NOT the same as "the impurity is singly occupied"
	// at this nBath=1 config: the first-bath hopping V=0.5 delocalizes that
	// one electron across the impurity and bath site, so <n_imp> itself is
	// ~0.5, not 1 -- the "impurity singly occupied" claim is exact only at
	// nBath=0 (no bath site to share the electron with), which this
	// particular test config (inherited from the L=1 SecondBathSeeding
	// test, nBath=1) does not exercise.
	CHECK(auxTotal == Catch::Approx(2.0 * L).margin(1e-4));
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
	latticeGf.updateLambda(0, gimp);
	solver.prepareTimeStep(0, latticeGf.lambda());

	const int nT = static_cast<int>(params.nT);
	for (int n = 1; n <= nT; ++n) {
		solver.computeGimp(gimp, n); // predictor
		for (SizeType iter = 0; iter < params.neqDmftIter; ++iter) {
			latticeGf.updateLambda(n, gimp);
			solver.prepareTimeStep(n, latticeGf.lambda());
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

// ---- Task 29: nBath=0 atomic-limit self-consistent gate -------------------
//
// Same cross-check philosophy as Task 15's FullGate above, but with
// bathParams={} (literal nBath=0, Wolf et al.'s / GBEK's NeqAtomicLimit
// starting point -- Lambda^-===0) instead of a real 5-site bath. This is
// the acceptance gate for Task 29's dual-sector averaging + null-branch
// skip (see Column::particleNull, ScSector, solveSelfConsistent's doc
// comment). TargetElectronsUp=1/Down=0 matches the atomic-limit
// constraint (nup_+ndown_==1) ImpuritySolverNeqExactDiag::solve() enforces
// for its own nBath==0 bypass -- exactDiag_ is a member of both solvers
// here, so both sides hit that same bypass for their equilibrium/
// Matsubara piece.
//
// STATUS (2026-07-31): known RED -- SIGSEGV, not just a failing assertion.
// Caught cleanly by Catch2's signal handler (reports as one failed test
// case, does not abort the binary/suite), so safe to leave in the tracked
// suite, same as Task 15's own documented residual. Root cause is NOT the
// dual-sector/null-branch-skip machinery this gate was written to verify
// (that machinery is independently confirmed correct -- see
// fancy-painting-moon.md Task #29's "continued" section): column 0's GS
// birth run has dir0:Connectors=[0,0] at nBath=0 (no first bath, second
// bath always inert at birth time) -- a fully disconnected 3-site system,
// pathological for two-site DMRG's density-matrix truncation independent
// of which branch (particle/hole) is later applied. Confirmed via
// h5ls -r comparison that the GS checkpoint itself is structurally normal
// (identical /Def/FinalPsi layout to a working nBath>0 run); the crash is
// in the subsequent birth run reading it. Not yet resolved -- see the plan
// file for the current best next step (whether DMRG++ supports injecting
// an explicit product-state checkpoint, bypassing the variational GS run
// for this one case).
static std::string configAtomicLimit(const std::string& rootName)
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
	       "vector InitBathVector=[0.5, 0.0];\n"
	       "int ImpuritySite=0;\n"
	       "real HubbardU=2.;\n"
	       "TargetElectronsUp=1;\n"
	       "TargetElectronsDown=0;\n"
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
	      "HubbardUFinal=2.;\n"
	      "TmaxNeq=0.2;\n"
	      "NtNeq=2;\n"
	      "NeqDmftIter=1;\n"
	      "NeqDmftTolerance=0.001;\n"
	      "NeqSolver=\"tdmrg\";\n"
	      "matrix FiniteLoopsTdmrg=[\n"
	      "    [@auto, 100, 0],[@auto, 100, 0],\n"
	      "    [@auto, 100, 0],[@auto, 100, 0]];\n"
	      "TSPTimeSteps=5;\n"
	      "TSPAdvanceEach=1;\n"
	      "NeqBathRank=1;\n";
}

TEST_CASE("ImpuritySolverNeqTdmrg vs ImpuritySolverNeqGBEK: nBath=0 atomic "
          "limit, NeqBathRank=1 self-consistent bath evolution agree",
          "[ImpuritySolverNeqTdmrg][Phase2][Task29][AtomicLimitGate]")
{
	int    argc    = 1;
	char   arg0[]  = "test_ImpuritySolverNeqTdmrg";
	char*  argv0[] = { arg0 };
	char** argv    = argv0;

	ApplicationType app("test_ImpuritySolverNeqTdmrg", &argc, &argv, 1);

	const VectorRealType emptyBathParams; // nBath=0

	InputNgType::Writeable ioWT(Dmft::CincuentaInputCheck {},
	                            configAtomicLimit("testTdmrgAtomicGate"));
	InputNgType::Readable  io_t(ioWT);
	ParamsType             paramsT(io_t);
	using TdmrgNeqSolverType = Dmft::NeqDmftSolver<ComplexType, Dmft::ImpuritySolverNeqTdmrg>;
	TdmrgNeqSolverType tdmrgSolver(paramsT, app, io_t);
	tdmrgSolver.solve(emptyBathParams);

	InputNgType::Writeable ioWG(Dmft::CincuentaInputCheck {},
	                            configAtomicLimit("testTdmrgAtomicGateGBEK"));
	InputNgType::Readable  io_g(ioWG);
	ParamsType             paramsG(io_g);
	using GbekNeqSolverType = Dmft::NeqDmftSolver<ComplexType, Dmft::ImpuritySolverNeqGBEK>;
	GbekNeqSolverType gbekSolver(paramsG, io_g);
	gbekSolver.solve(emptyBathParams);

	const auto& gimpT = tdmrgSolver.gimp();
	const auto& gimpG = gbekSolver.gimp();

	const int      nT2 = static_cast<int>(paramsT.nT);
	const RealType tol = 1e-4;
	for (int n = 0; n <= nT2; ++n) {
		for (int j = 0; j <= n; ++j) {
			const ComplexType retT = gimpT.retarded(n, j);
			const ComplexType lesT = gimpT.lesser(n, j);
			const ComplexType retG = gimpG.retarded(n, j);
			const ComplexType lesG = gimpG.lesser(n, j);

			std::cout << "n=" << n << " j=" << j << " tDMRG G^R=" << retT
			          << " G^<=" << lesT << " GBEK G^R=" << retG << " G^<=" << lesG
			          << "\n";

			CHECK(retT.real() == Catch::Approx(retG.real()).margin(tol));
			CHECK(retT.imag() == Catch::Approx(retG.imag()).margin(tol));
			CHECK(lesT.real() == Catch::Approx(lesG.real()).margin(tol));
			CHECK(lesT.imag() == Catch::Approx(lesG.imag()).margin(tol));
		}
	}
}

// ---- Task 44: in-bounds NeqBathRank sanity check for d(t) --------------
//
// Per user direction (2026-08-04): NeqBathRank<3 is out of bounds (GBEK's
// own paper never validates rank<2/3; the committed NeqBathRank=1
// AtomicLimitGate test above is a known, out-of-bounds SIGABRT, not a
// real blocker). This test checks two things at NeqBathRank=3,4,5,
// nBath=0 (atomic limit): (1) does the self-consistent run complete at
// all (the NeqBathRank=1 crash's root cause -- disconnected-3-site GS --
// is not obviously rank-specific, so this is not assumed); (2) does d(t)
// look physically sane (finite, in [0,0.25]) and how does it move as
// rank increases (a first, cheap look at internal convergence -- NOT a
// tDMRG-vs-GBEK numeric comparison, which the standing comparison-
// methodology note rules out at matching rank).
// DIAGNOSTIC (temporary): same as configAtomicLimitRank but with an
// explicit NtNeq, so a caller can request NtNeq>=NeqBathRank -- the only
// configuration in which NeqBathDecomposition::update actually seeds every
// orbital-pair (row n seeds pair n-1, so pairs 0..NtNeq-1 get seeded; at
// NtNeq<NeqBathRank, pairs NtNeq..NeqBathRank-1 sit at exact zero coupling
// AND (per scPotTdmrg_, which never eps-splits auxiliary sites outside the
// one-off GS birth) exact zero potential-splitting for the whole run --
// an unbroken exact degeneracy, not merely "not yet distinguished". Not a
// permanent gate -- delete/inline once the degeneracy-vs-kept-states
// question is settled.
static std::string configAtomicLimitRankNtNeq(const std::string& rootName, int rank, int ntNeq)
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
	       "vector InitBathVector=[0.5, 0.0];\n"
	       "int ImpuritySite=0;\n"
	       "real HubbardU=2.;\n"
	       "TargetElectronsUp=1;\n"
	       "TargetElectronsDown=0;\n"
	       "RootOutputname=\""
	    + rootName
	    + "\";\n"
	      "NeqOutputPrefix=\""
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
	      "HubbardUFinal=2.;\n"
	      "TmaxNeq=0.2;\n"
	      "NtNeq="
	    + ttos(ntNeq)
	    + ";\n"
	      "NeqDmftIter=1;\n"
	      "NeqDmftTolerance=0.001;\n"
	      "NeqSolver=\"tdmrg\";\n"
	      "matrix FiniteLoopsTdmrg=[\n"
	      "    [@auto, 100, 0],[@auto, 100, 0],\n"
	      "    [@auto, 100, 0],[@auto, 100, 0]];\n"
	      "TSPTimeSteps=5;\n"
	      "TSPAdvanceEach=1;\n"
	      "NeqBathRank="
	    + ttos(rank) + ";\n";
}

// Task 44 follow-up: NtNeq=2 (< NeqBathRank) at rank>=3 leaves the
// unseeded orbital-pairs at exact zero coupling AND (per scPotTdmrg_,
// which never eps-splits auxiliary sites outside the one-off GS birth)
// exact zero potential-splitting for the whole run -- an unbroken exact
// degeneracy, not merely "not yet distinguished". This was confirmed to be
// the actual cause of a "printSumAndCheckEigs: DM eigs don't amount to
// one" failure seen at NeqBathRank=4/NtNeq=2 (preceded by a LAPACK gesdd
// info=1 non-convergence warning -- a classic symptom of clustered/
// degenerate singular values), NOT an inadequate kept-states cap: the SAME
// kept-states=100 cap that failed at NtNeq=2 passes cleanly once NtNeq is
// raised to fully seed every orbital-pair (see the NtNeq=rank tests
// below). Every NtNeq=2 test at rank>=3 -- the grouped {3,4,5} test, the
// isolated rank=5-alone test, and the isolated rank=4 investigation test
// used to find this -- has been removed: each one pins an ill-posed,
// permanently-degenerate configuration that was never physically
// meaningful to test in the first place, and would be permanently RED (or
// worse, flaky, depending on how badly gesdd struggles on a given run) for
// no useful reason. See TDMRG_EVOLVING_BATH.md / project memory for the
// full investigation.

// NeqBathRank=3 with NtNeq=3 (NtNeq>=NeqBathRank), the cheapest
// configuration in which every orbital-pair actually gets seeded at least
// once (pairs 0,1,2 via n=1,2,3) -- no permanently-inert, exactly-
// degenerate spectator sites remain for the whole run. Same kept-states
// cap (100) as the committed Task 44 {3,4,5} test above. This is the
// well-posed counterpart of that test and completes cleanly.
TEST_CASE("Task 44: tDMRG self-consistent d(t) at NeqBathRank=3, "
          "NtNeq=3 (fully seeded), atomic limit",
          "[ImpuritySolverNeqTdmrg][Task44][AtomicLimitDocc]")
{
	int    argc    = 1;
	char   arg0[]  = "test_ImpuritySolverNeqTdmrg";
	char*  argv0[] = { arg0 };
	char** argv    = argv0;

	ApplicationType      app("test_ImpuritySolverNeqTdmrg", &argc, &argv, 1);
	const VectorRealType emptyBathParams; // nBath=0

	const int         rank     = 3;
	const std::string rootName = "testTdmrgAtomicRank" + ttos(rank) + "FullSeed";

	InputNgType::Writeable ioW(Dmft::CincuentaInputCheck {},
	                           configAtomicLimitRankNtNeq(rootName, rank, rank));
	InputNgType::Readable  io(ioW);
	ParamsType             params(io);

	using TdmrgNeqSolverType = Dmft::NeqDmftSolver<ComplexType, Dmft::ImpuritySolverNeqTdmrg>;
	TdmrgNeqSolverType solver(params, app, io);
	solver.solve(emptyBathParams);

	const auto& gimp = solver.gimp();
	const int   nT2  = static_cast<int>(params.nT);
	for (int n = 0; n <= nT2; ++n) {
		for (int j = 0; j <= n; ++j) {
			CHECK(std::isfinite(gimp.retarded(n, j).real()));
			CHECK(std::isfinite(gimp.lesser(n, j).real()));
		}
	}

	solver.dumpGreenFunctions();
	std::ifstream doccFile(rootName + "-docc-energy");
	REQUIRE(doccFile.good());
	RealType t = 0, docc = 0, ekin = 0, eint = 0, etot = 0;
	while (doccFile >> t >> docc >> ekin >> eint >> etot) {
		std::cout << "rank=" << rank << " NtNeq=" << rank << " t=" << t << " docc=" << docc
		          << "\n";
		CHECK(std::isfinite(docc));
		CHECK(docc >= RealType(-1e-6));
		CHECK(docc <= RealType(0.25 + 1e-6));
	}
}

// NeqBathRank=4 with NtNeq=4 (fully seeded), same kept-states cap (100) as
// the NtNeq=2 rank=4 config that originally failed with
// "printSumAndCheckEigs: DM eigs don't amount to one". Confirms the
// degeneracy hypothesis at the rank where the original failure was
// observed: this completes cleanly under the SAME kept-states=100 cap
// that failed at NtNeq=2, so the original failure was an artifact of the
// permanently-degenerate (unseeded pairs 2,3) NtNeq=2 config, not a real
// kept-states-too-small or engine bug. Cost note: this run takes ~15
// minutes wall-clock (9-site lattice, 4 real time steps with genuine
// Lanczos refinement) -- worth tagging separately if CI runtime matters.
TEST_CASE("Task 44: tDMRG self-consistent d(t) at NeqBathRank=4, "
          "NtNeq=4 (fully seeded), atomic limit",
          "[ImpuritySolverNeqTdmrg][Task44][AtomicLimitDocc][Slow]")
{
	int    argc    = 1;
	char   arg0[]  = "test_ImpuritySolverNeqTdmrg";
	char*  argv0[] = { arg0 };
	char** argv    = argv0;

	ApplicationType      app("test_ImpuritySolverNeqTdmrg", &argc, &argv, 1);
	const VectorRealType emptyBathParams; // nBath=0

	const int         rank     = 4;
	const std::string rootName = "testTdmrgAtomicRank" + ttos(rank) + "FullSeed";

	InputNgType::Writeable ioW(Dmft::CincuentaInputCheck {},
	                           configAtomicLimitRankNtNeq(rootName, rank, rank));
	InputNgType::Readable  io(ioW);
	ParamsType             params(io);

	using TdmrgNeqSolverType = Dmft::NeqDmftSolver<ComplexType, Dmft::ImpuritySolverNeqTdmrg>;
	TdmrgNeqSolverType solver(params, app, io);
	solver.solve(emptyBathParams);

	const auto& gimp = solver.gimp();
	const int   nT2  = static_cast<int>(params.nT);
	for (int n = 0; n <= nT2; ++n) {
		for (int j = 0; j <= n; ++j) {
			CHECK(std::isfinite(gimp.retarded(n, j).real()));
			CHECK(std::isfinite(gimp.lesser(n, j).real()));
		}
	}

	solver.dumpGreenFunctions();
	std::ifstream doccFile(rootName + "-docc-energy");
	REQUIRE(doccFile.good());
	RealType t = 0, docc = 0, ekin = 0, eint = 0, etot = 0;
	while (doccFile >> t >> docc >> ekin >> eint >> etot) {
		std::cout << "rank=" << rank << " NtNeq=" << rank << " t=" << t << " docc=" << docc
		          << "\n";
		CHECK(std::isfinite(docc));
		CHECK(docc >= RealType(-1e-6));
		CHECK(docc <= RealType(0.25 + 1e-6));
	}
}

// ---- Advisor-recommended discriminator: does a REAL nonzero nBath (not
// nBath=0) at birth avoid the disconnected-lattice WFT crash entirely?
// nBath=2 gives the impurity real, nonzero hoppings from t=0 -- the
// lattice is never fully disconnected even while the NeqBathRank=3
// second-bath orbital-pairs are still being seeded. If this runs clean,
// the bug localizes to "impurity has zero connections at birth"
// specifically (a narrower, possibly cincuenta-side fix), not "the
// seeding window is structurally hostile to DMRG" in general.
TEST_CASE("Task 44 discriminator: nBath=2 (real bath), NeqBathRank=3 -- "
          "does a real bath at birth avoid the WFT crash",
          "[ImpuritySolverNeqTdmrg][Task44][RealBathDiscriminator]")
{
	int    argc    = 1;
	char   arg0[]  = "test_ImpuritySolverNeqTdmrg";
	char*  argv0[] = { arg0 };
	char** argv    = argv0;

	ApplicationType app("test_ImpuritySolverNeqTdmrg", &argc, &argv, 1);

	// nBath=2 -> nsites=3 (impurity+2 bath); ExactDiag's equilibrium
	// delegate requires half-filling (nup+ndown==nsites==3), so override
	// configWithU's default 3/3 (which was sized for its own 5-bath
	// config) down to 2/1.
	InputNgType::Writeable ioW(Dmft::CincuentaInputCheck {},
	                           configWithU("0.5", "testTdmrgRealBathDisc")
	                               + "NeqBathRank=3;\nTargetElectronsUp=2;\n"
	                                 "TargetElectronsDown=1;\n");
	InputNgType::Readable  io(ioW);
	ParamsType             params(io);

	using TdmrgNeqSolverType = Dmft::NeqDmftSolver<ComplexType, Dmft::ImpuritySolverNeqTdmrg>;
	TdmrgNeqSolverType solver(params, app, io);

	const VectorRealType bathParams = { 0.3, 0.3, -0.3, 0.3 }; // nBath=2
	solver.solve(bathParams);

	const auto& gimp = solver.gimp();
	const int   nT2  = static_cast<int>(params.nT);
	for (int n = 0; n <= nT2; ++n) {
		for (int j = 0; j <= n; ++j) {
			CHECK(std::isfinite(gimp.retarded(n, j).real()));
			CHECK(std::isfinite(gimp.lesser(n, j).real()));
		}
	}
}

// ---- DIAGNOSTIC (temporary): isolate whether fillSelfConsistentRow(gimp,2)
// is wrong standalone, or only when preceded by calls for n=0,1 within the
// same solver's lifetime (same reused chainRoot files). Not a permanent
// gate -- delete once the root cause is understood.
TEST_CASE("DIAGNOSTIC: computeGimp(gimp,2) called directly, skipping n=0,1",
          "[ImpuritySolverNeqTdmrg][Diagnostic]")
{
	int    argc    = 1;
	char   arg0[]  = "test_ImpuritySolverNeqTdmrg";
	char*  argv0[] = { arg0 };
	char** argv    = argv0;

	ApplicationType app("test_ImpuritySolverNeqTdmrg", &argc, &argv, 1);

	const VectorRealType bathParams = { 0.3, 0.3, 0.3, 0.3, 0.3, -0.6, -0.3, 0.0, 0.3, 0.6 };

	InputNgType::Writeable ioW(Dmft::CincuentaInputCheck {},
	                           configWithU("0.5", "testTdmrgChainDiag") + "NeqBathRank=1;\n");
	InputNgType::Readable  io(ioW);
	ParamsType             params(io);
	SolverType             solver(params, app, io);
	solver.solve(bathParams);

	Dmft::KadanoffBaym<ComplexType> gimp(
	    params.nT,
	    params.eqParams.nMatsubaras,
	    params.dt,
	    params.eqParams.ficticiousBeta / static_cast<RealType>(params.eqParams.nMatsubaras));

	solver.computeGimp(gimp, 2);

	std::cout << "DIAGNOSTIC standalone computeGimp(gimp,2): retarded(2,0)="
	          << gimp.retarded(2, 0) << " lesser(2,0)=" << gimp.lesser(2, 0)
	          << " retarded(2,1)=" << gimp.retarded(2, 1)
	          << " lesser(2,1)=" << gimp.lesser(2, 1) << "\n";
}

// ---- DIAGNOSTIC (temporary): does a column's SECOND advance actually
// respect a NEW Connectors value, or does it silently reuse the first
// advance's? Not a permanent gate -- delete once understood.
TEST_CASE("DIAGNOSTIC second advance connectors sensitivity",
          "[ImpuritySolverNeqTdmrg][Diagnostic]")
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

	const VectorComplexType connectorsStep1
	    = { ComplexType(0.26, 0.0), ComplexType(0.26, 0.0) };
	const VectorComplexType connectorsStep2A
	    = { ComplexType(0.26, 0.0), ComplexType(0.26, 0.0) };
	const VectorComplexType connectorsStep2B
	    = { ComplexType(0.52, -0.008), ComplexType(0.52, -0.008) };

	const auto resultA = solver.diagnosticSecondAdvanceConnectors(
	    bathParams, L, eps, connectorsStep1, connectorsStep2A, "diag2stepA_");
	const auto resultB = solver.diagnosticSecondAdvanceConnectors(
	    bathParams, L, eps, connectorsStep1, connectorsStep2B, "diag2stepB_");

	std::cout << "DIAGNOSTIC second-advance sensitivity: resultA=" << resultA
	          << " resultB=" << resultB << "\n";
}

// Positive control for the above: hold connectorsStep2 FIXED and instead vary
// connectorsStep1 (the first advance's Connectors). If the step-2 measurement
// tracks connectorsStep1 instead, that positively confirms "frozen at the
// FIRST advance's Hamiltonian" (not "frozen at zero" or "measurement is
// junk/independent of everything"). See advisor consult, 2026-07-29.
TEST_CASE("DIAGNOSTIC second advance connectors sensitivity -- positive control",
          "[ImpuritySolverNeqTdmrg][Diagnostic]")
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

	const VectorComplexType connectorsStep2Fixed
	    = { ComplexType(0.26, 0.0), ComplexType(0.26, 0.0) };
	const VectorComplexType connectorsStep1A
	    = { ComplexType(0.26, 0.0), ComplexType(0.26, 0.0) };
	const VectorComplexType connectorsStep1B
	    = { ComplexType(0.52, -0.008), ComplexType(0.52, -0.008) };

	const auto resultA = solver.diagnosticSecondAdvanceConnectors(
	    bathParams, L, eps, connectorsStep1A, connectorsStep2Fixed, "diag2stepC_");
	const auto resultB = solver.diagnosticSecondAdvanceConnectors(
	    bathParams, L, eps, connectorsStep1B, connectorsStep2Fixed, "diag2stepD_");

	std::cout << "DIAGNOSTIC positive control: resultA=" << resultA << " resultB=" << resultB
	          << "\n";
}

// Discriminating check (advisor consult, 2026-07-29): three chained advances
// with Connectors C1,C2,C3; harvest ggtRaw[3] (occurrence #1, a continuation
// segment) while varying ONLY C2. If the harvested value moves with C2, that
// confirms C2 reaches the harvest one step LATER (an off-by-one in which
// step's Connectors a continuation segment's occurrence #1 reflects, since
// occurrence #1 is a readout of the state the PREVIOUS segment produced) --
// not an engine freeze. See TDMRG_EVOLVING_BATH.md Link 9.
TEST_CASE("DIAGNOSTIC third advance -- does step-3 harvest track C2 (off-by-one check)",
          "[ImpuritySolverNeqTdmrg][Diagnostic]")
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

	const VectorComplexType c1  = { ComplexType(0.26, 0.0), ComplexType(0.26, 0.0) };
	const VectorComplexType c2A = { ComplexType(0.26, 0.0), ComplexType(0.26, 0.0) };
	const VectorComplexType c2B = { ComplexType(0.52, -0.008), ComplexType(0.52, -0.008) };
	const VectorComplexType c3  = { ComplexType(0.26, 0.0), ComplexType(0.26, 0.0) };

	const auto resultA = solver.diagnosticThirdAdvanceConnectors(
	    bathParams, L, eps, c1, c2A, c3, "diag3stepA_");
	const auto resultB = solver.diagnosticThirdAdvanceConnectors(
	    bathParams, L, eps, c1, c2B, c3, "diag3stepB_");

	std::cout << "DIAGNOSTIC third-advance off-by-one check: resultA=" << resultA
	          << " resultB=" << resultB << "\n";
}

// Two calibration attempts were tried here and abandoned (both refuted by
// data, not by inspection -- see TDMRG_EVOLVING_BATH.md Link 9 for the full
// account): (1) scaling advanceEach with border-gated firing unchanged --
// plateaued at the "2 advances/segment" value for every trial from 7-12,
// never reaching 1; (2) adding advanceUnrestricted (removing the border
// gate found in dmrg/Engine/ApplyOperatorExpression.h:769-773) and scaling
// advanceEach as a pure count -- plateaued at "2 steps total across 3
// calls" for trials 11-13 and "0 advances at all" for trials >=14, never
// converging on "1 advance/segment" for any value. advanceEach cannot
// express "exactly one advance per two-row segment" under either gating
// mode with this engine. A real fix needs either dmrg/Engine/ instrumentation
// or an outer-loop redesign around the 2-dt-per-segment time quantum --
// both are scope decisions, not something to keep scanning for here.

// Engine fix (option 1, user-approved 2026-07-29): dmrg/Engine's
// NonLocalForTargetingExpression::advanceInTimeOrNot hardcodes
// advanceOnlyAtBorder=true with no existing opt-out -- a genuinely
// different code path from ApplyOperatorExpression.h's border gate (the
// advanceUnrestricted flag above targeted the WRONG mechanism entirely,
// which is why it never converged). Added an opt-in maxAdvances cap
// (0=unlimited, the default -- every existing caller unaffected) on
// GroupOfOneTimeEvolutions::OneTimeEvolution, threaded through
// buildStepInput's TimeEvolve{...} string. This test verifies maxAdvances=1
// makes each segment fire exactly once: harvest imaginary parts across a
// 3-call static-Connectors chain should land at 1x, 2x, 3x a single
// dt-advance (not 2x, 4x, 6x as confirmed without the cap).
TEST_CASE("DIAGNOSTIC maxAdvances=1 gives exactly one advance per segment",
          "[ImpuritySolverNeqTdmrg][Diagnostic]")
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

	const VectorRealType    bathParams = { 0.3, 0.3, 0.3, 0.3, 0.3, -0.6, -0.3, 0.0, 0.3, 0.6 };
	const SizeType          L          = 1;
	const RealType          eps        = 3.0;
	const VectorComplexType c          = { ComplexType(0.26, 0.0), ComplexType(0.26, 0.0) };

	const auto result = solver.diagnosticThirdAdvanceConnectors(
	    bathParams, L, eps, c, c, c, "diagmaxadv_", 0, 1);

	std::cout << "DIAGNOSTIC maxAdvances=1 static-C step-3 harvest=" << result << "\n";
}

// With maxAdvances=1, harvest(3) should now depend on C3 (this call's OWN
// Connectors) -- the opposite of the pre-fix lag, where only C2 mattered.
TEST_CASE("DIAGNOSTIC maxAdvances=1 -- step-3 harvest now tracks its OWN Connectors (C3)",
          "[ImpuritySolverNeqTdmrg][Diagnostic]")
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

	const VectorRealType    bathParams = { 0.3, 0.3, 0.3, 0.3, 0.3, -0.6, -0.3, 0.0, 0.3, 0.6 };
	const SizeType          L          = 1;
	const RealType          eps        = 3.0;
	const VectorComplexType c1         = { ComplexType(0.26, 0.0), ComplexType(0.26, 0.0) };
	const VectorComplexType c2         = { ComplexType(0.26, 0.0), ComplexType(0.26, 0.0) };
	const VectorComplexType c3A        = { ComplexType(0.26, 0.0), ComplexType(0.26, 0.0) };
	const VectorComplexType c3B = { ComplexType(0.52, -0.008), ComplexType(0.52, -0.008) };

	const auto resultA = solver.diagnosticThirdAdvanceConnectors(
	    bathParams, L, eps, c1, c2, c3A, "diagmaxadvA_", 0, 1);
	const auto resultB = solver.diagnosticThirdAdvanceConnectors(
	    bathParams, L, eps, c1, c2, c3B, "diagmaxadvB_", 0, 1);

	std::cout << "DIAGNOSTIC maxAdvances=1 C3 sensitivity: resultA=" << resultA
	          << " resultB=" << resultB << "\n";
}

TEST_CASE("DIAGNOSTIC Link 12: does a born-mid-chain column's first advance "
          "since birth respond to its OWN Connectors?",
          "[ImpuritySolverNeqTdmrg][Diagnostic]")
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
	// Birth Connectors (used for column0's t0->t1 advance and column1's
	// birth) held fixed; column1's own advance-to-2 Connectors (its FIRST
	// advance since birth) differs between A and B. If column1's harvest
	// doesn't change, the birth checkpoint is killing Connectors
	// sensitivity on the very next segment -- a hard, narrow reproduction
	// of Link 12. If it does change, the bug is elsewhere (comparison
	// convention, applySignFlip, etc.), not this mechanism.
	const VectorComplexType c1  = { ComplexType(0.26, 0.0), ComplexType(0.26, 0.0) };
	const VectorComplexType c2A = { ComplexType(0.26, 0.0), ComplexType(0.26, 0.0) };
	const VectorComplexType c2B = { ComplexType(0.52, -0.008), ComplexType(0.52, -0.008) };

	const auto resultA = solver.diagnosticColumn1FirstAdvanceConnectors(
	    bathParams, L, eps, c1, c2A, "diagcol1A_", 0, 1);
	const auto resultB = solver.diagnosticColumn1FirstAdvanceConnectors(
	    bathParams, L, eps, c1, c2B, "diagcol1B_", 0, 1);

	std::cout << "DIAGNOSTIC Link12 column1-first-advance sensitivity: resultA=" << resultA
	          << " resultB=" << resultB << "\n";
}

TEST_CASE("DIAGNOSTIC Link 12: eps-split seeding check for the GATE's own "
          "bathParams/eps (5-site bath, not the 1-site config Task 8 validated)",
          "[ImpuritySolverNeqTdmrg][Diagnostic]")
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

	// Same bathParams as the FullGate test. maxCoupling = max(U=0.5,
	// hoppings=0.3, bathEps up to 0.6) = 0.6, so scEps_ = 5*0.6 = 3.0 --
	// reproduce that computation exactly rather than hardcoding 3.0, in
	// case solveSelfConsistent's formula and this diverge.
	const VectorRealType bathParams  = { 0.3, 0.3, 0.3, 0.3, 0.3, -0.6, -0.3, 0.0, 0.3, 0.6 };
	const SizeType       L           = 1;
	RealType             maxCoupling = 0.5;
	for (SizeType i = 0; i < 5; ++i) {
		maxCoupling = std::max(maxCoupling, std::abs(bathParams[i]));
		maxCoupling = std::max(maxCoupling, std::abs(bathParams[5 + i]));
	}
	const RealType eps = 5.0 * maxCoupling;

	const auto occ = solver.measureSecondBathOccupations(bathParams, L, eps);
	REQUIRE(occ.size() == 2 * L);
	std::cout << "DIAGNOSTIC Link12 gate-config eps-split seeding: eps=" << eps
	          << " occ[0](intended-occupied)=" << occ[0] << " occ[1](intended-empty)=" << occ[1]
	          << "\n";
}

// ---- Task 17: fillSelfConsistentRow's incremental design (persisting
// Column state across calls instead of rebuilding from scratch every
// call) rests on DmrgRunner being bit-deterministic given an identical
// input file and on-disk restart source -- the original always-recompute
// design never needed this property. Confirmed once before implementing
// (fancy-painting-moon.md, Task #17 scope), kept as permanent regression
// coverage since the shipped design depends on it going forward, not just
// at the point it was checked. Runs the SAME 3-advance chain twice
// (identical bathParams/L/eps/Connectors, only rootSuffix differs so the
// two runs don't collide on filenames) and requires the two results be
// EXACTLY equal.
TEST_CASE("DIAGNOSTIC Task 17: DmrgRunner is bit-deterministic given identical inputs",
          "[ImpuritySolverNeqTdmrg][Diagnostic]")
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

	const VectorRealType    bathParams = { 0.3, 0.3, 0.3, 0.3, 0.3, -0.6, -0.3, 0.0, 0.3, 0.6 };
	const SizeType          L          = 1;
	const RealType          eps        = 3.0;
	const VectorComplexType c1         = { ComplexType(0.26, 0.0), ComplexType(0.26, 0.0) };
	const VectorComplexType c2
	    = { ComplexType(0.524393, -0.0077008), ComplexType(0.524393, -0.0077008) };
	const VectorComplexType c3 = { ComplexType(0.4, 0.1), ComplexType(0.4, 0.1) };

	const auto resultA = solver.diagnosticThirdAdvanceConnectors(
	    bathParams, L, eps, c1, c2, c3, "diagdetA_", 0, 1);
	const auto resultB = solver.diagnosticThirdAdvanceConnectors(
	    bathParams, L, eps, c1, c2, c3, "diagdetB_", 0, 1);

	std::cout << "DIAGNOSTIC Task17 determinism: resultA=" << resultA << " resultB=" << resultB
	          << " diff=" << (resultA - resultB) << "\n";

	CHECK(resultA.real() == resultB.real());
	CHECK(resultA.imag() == resultB.imag());
}

// ---- Task 17 regression pin: the incremental fillSelfConsistentRow
// (O(nT^2*(1+neqDmftIter))) must keep reproducing the exact numbers the
// original "truncated batch recompute" (O(nT^3*neqDmftIter)) produced.
// The batch implementation itself was edited away during Task 17 (not
// kept as a standing reference -- see fancy-painting-moon.md, Task #17),
// so this pins the values that were manually confirmed bit-identical
// between the two implementations at every one of Task 17's five
// migration steps (see project_tdmrg_evolving_bath memory), independent
// of the separate GBEK cross-check (whose 1e-4 margin is far looser than
// what this test can catch). Same bathParams/config as the FullGate
// test, but tDMRG-only (no GBEK solve) -- keeps this fast and catches any
// future accidental drift in the incremental algorithm itself, not just
// gross disagreement with GBEK.
TEST_CASE("ImpuritySolverNeqTdmrg NeqBathRank=1 self-consistent: incremental "
          "algorithm matches pinned batch-recompute reference values",
          "[ImpuritySolverNeqTdmrg][Phase2][RegressionPin]")
{
	int    argc    = 1;
	char   arg0[]  = "test_ImpuritySolverNeqTdmrg";
	char*  argv0[] = { arg0 };
	char** argv    = argv0;

	ApplicationType app("test_ImpuritySolverNeqTdmrg", &argc, &argv, 1);

	const VectorRealType bathParams = { 0.3, 0.3, 0.3, 0.3, 0.3, -0.6, -0.3, 0.0, 0.3, 0.6 };

	InputNgType::Writeable ioW(Dmft::CincuentaInputCheck {},
	                           configWithU("0.5", "testTdmrgChainRegPin") + "NeqBathRank=1;\n");
	InputNgType::Readable  io(ioW);
	ParamsType             params(io);
	using TdmrgNeqSolverType = Dmft::NeqDmftSolver<ComplexType, Dmft::ImpuritySolverNeqTdmrg>;
	TdmrgNeqSolverType solver(params, app, io);
	solver.solve(bathParams);

	const auto& gimp = solver.gimp();

	// Values pinned from the last batch-recompute run before Task 17's
	// rewrite (and reconfirmed bit-identical after every migration step).
	// Tolerance (1e-6) is far tighter than the FullGate test's GBEK-
	// comparison margin (1e-4) -- this test isn't checking physics
	// against an independent reference, it's checking the incremental
	// rewrite against its own known-correct predecessor.
	const RealType tol = 1e-6;

	CHECK(gimp.retarded(0, 0).real() == Catch::Approx(0.0).margin(tol));
	CHECK(gimp.retarded(0, 0).imag() == Catch::Approx(-1.0).margin(tol));
	CHECK(gimp.lesser(0, 0).real() == Catch::Approx(0.0).margin(tol));
	CHECK(gimp.lesser(0, 0).imag() == Catch::Approx(0.5).margin(tol));

	CHECK(gimp.retarded(1, 0).real() == Catch::Approx(0.0).margin(tol));
	CHECK(gimp.retarded(1, 0).imag() == Catch::Approx(-0.996754).margin(tol));
	CHECK(gimp.lesser(1, 0).real() == Catch::Approx(-0.031695).margin(tol));
	CHECK(gimp.lesser(1, 0).imag() == Catch::Approx(0.498377).margin(tol));

	CHECK(gimp.retarded(1, 1).real() == Catch::Approx(0.0).margin(tol));
	CHECK(gimp.retarded(1, 1).imag() == Catch::Approx(-1.0).margin(tol));
	CHECK(gimp.lesser(1, 1).real() == Catch::Approx(0.0).margin(tol));
	CHECK(gimp.lesser(1, 1).imag() == Catch::Approx(0.5).margin(tol));

	CHECK(gimp.retarded(2, 0).real() == Catch::Approx(0.000118).margin(tol));
	CHECK(gimp.retarded(2, 0).imag() == Catch::Approx(-0.983645).margin(tol));
	CHECK(gimp.lesser(2, 0).real() == Catch::Approx(-0.063042).margin(tol));
	CHECK(gimp.lesser(2, 0).imag() == Catch::Approx(0.491821).margin(tol));

	CHECK(gimp.retarded(2, 1).real() == Catch::Approx(0.0).margin(tol));
	CHECK(gimp.retarded(2, 1).imag() == Catch::Approx(-0.994700).margin(tol));
	CHECK(gimp.lesser(2, 1).real() == Catch::Approx(-0.031684).margin(tol));
	CHECK(gimp.lesser(2, 1).imag() == Catch::Approx(0.497352).margin(tol));

	CHECK(gimp.retarded(2, 2).real() == Catch::Approx(0.0).margin(tol));
	CHECK(gimp.retarded(2, 2).imag() == Catch::Approx(-1.0).margin(tol));
	CHECK(gimp.lesser(2, 2).real() == Catch::Approx(0.0).margin(tol));
	CHECK(gimp.lesser(2, 2).imag() == Catch::Approx(0.500011).margin(tol));
}

// ---- Task 28: Wolf et al. (PRB 90, 235131) cosine hopping-ramp protocol ----
//
// Confirms the ramp mechanism itself (NeqLatticeGf::tStarAt, feeding
// updateLambda) needs NO impurity-solver-specific code: a pure NeqLatticeGf
// unit test with a synthetic constant gimp, no DMRG/ED involved at all.
// NeqAtomicLimit=1 forces t*=0 at n=0 (the paper's atomic-limit starting
// point -- "no impurity-bath correlations in the initial state"); QuenchShape
// "cosine"/QuenchDuration=0.25 then ramps t* toward BandwidthFinal/4 per the
// paper's v(t). Grep confirms ImpuritySolverNeqGBEK.h reads none of these
// three params either -- if GBEK needs zero solver-side ramp code, tDMRG
// needs zero too. See plan file Task #28 "What research found".
TEST_CASE("NeqLatticeGf cosine ramp: Lambda(n,n) starts at zero (atomic "
          "limit) and grows monotonically as t* ramps up",
          "[NeqLatticeGf][Task28][RampMechanics]")
{
	InputNgType::Writeable ioW(Dmft::CincuentaInputCheck {},
	                           configWithU("0.5", "testRampMechanics")
	                               + "QuenchShape=\"cosine\";\n"
	                                 "QuenchDuration=0.25;\n"
	                                 "BandwidthFinal=4;\n"
	                                 "NeqAtomicLimit=1;\n");
	InputNgType::Readable  io(ioW);
	ParamsType             params(io);

	Dmft::NeqLatticeGf<ComplexType> latticeGf(params);
	Dmft::KadanoffBaym<ComplexType> gimp(
	    params.nT,
	    params.eqParams.nMatsubaras,
	    params.dt,
	    params.eqParams.ficticiousBeta / static_cast<RealType>(params.eqParams.nMatsubaras));

	// Synthetic, constant-in-time atomic-limit-like G_imp: with a constant
	// gimp, any n-dependence seen in Lambda(n,n) can only come from the
	// ramp's own t*(n) scaling, not from the impurity solver.
	const int nT = static_cast<int>(params.nT);
	for (int n = 0; n <= nT; ++n) {
		for (int j = 0; j <= n; ++j) {
			gimp.retarded(n, j) = ComplexType(0, -1);
			gimp.lesser(n, j)   = ComplexType(0, 0.5);
		}
	}

	latticeGf.initialize(gimp);
	latticeGf.updateLambda(0, gimp);

	// NeqAtomicLimit forces t*(0)=0, so Lambda(0,0)=t*(0)^2*gimp(0,0) must
	// be exactly zero -- the paper's atomic-limit starting condition.
	CHECK(latticeGf.lambda().retarded(0, 0).real() == Catch::Approx(0.0).margin(1e-12));
	CHECK(latticeGf.lambda().retarded(0, 0).imag() == Catch::Approx(0.0).margin(1e-12));

	// TmaxNeq=0.2, NtNeq=2 (configWithU's fixed values) => dt=0.1, so
	// n=1,2 land at t=0.1,0.2, both still short of QuenchDuration=0.25 --
	// strictly inside the ramp, never reaching the tStarFinal_ plateau.
	RealType prevMag = 0;
	for (int n = 1; n <= nT; ++n) {
		latticeGf.updateLambda(n, gimp);
		const RealType mag = std::abs(latticeGf.lambda().retarded(n, n));
		CHECK(mag > prevMag);
		prevMag = mag;
	}
}

// ---- Task 28: config-only smoke test under a genuine, non-degenerate ------
//               cosine ramp (NOT the literal atomic limit)
//
// Reuses configWithU + the FullGate test's EXACT proven bathParams
// (nBath=5, symmetric nonzero hoppings, already known non-degenerate) --
// only the ramp params (QuenchShape/QuenchDuration/BandwidthFinal/
// NeqAtomicLimit) are added on top. This deliberately does NOT attempt
// GBEK/Wolf's literal "hoppings=0, decoupled first bath" atomic-limit
// starting point: reading GBEK's own production reference inputs
// (inputNeqAtomicLimitGBEKL3_fastramp.ain) revealed that NeqAtomicLimit=1
// makes GBEK BYPASS the equilibrium bath fit and construct a genuine
// nBath=0 (bare impurity) Fock space for the real neq computation --
// NumberOfBathPoints=1/TargetElectronsUp=1,Down=0 in that file are
// placeholder values whose result is discarded, not the real bath. tDMRG's
// star geometry cannot do literal nBath=0 (TSPAdvanceEach=nsites-2
// degenerates at nsites=1, already documented in Phase 1) -- so reaching
// GBEK's literal atomic limit is a SEPARATE, harder capability gap, not
// something "config-only" can close; it needs its own task (see plan file).
//
// Multiple attempts to approximate the atomic limit via "nBath>=1,
// hoppings=0" instead (matching this file's existing inert-second-bath
// precedent, just applied to the first bath) were tried and abandoned:
// symmetric per-spin target counts (Up=Down=3, nBath=5) hit a genuine
// ground-state degeneracy (the decoupled impurity's implicit level tied
// exactly with a bath eps value at the filled/empty boundary, confirmed via
// the advisor and by hand-computing the level ranking -- both solvers broke
// the tie differently, producing a same-magnitude opposite-sign Re[G^R]);
// a fully spin-polarized fix (Up=n,Down=0) hit ImpuritySolverNeqExactDiag's
// hard-enforced nup+ndown==nsites (documented gap, memory
// project_neq_ed_arbitrary_filling) and its Lehmann N+-1 sector needing an
// interior nup; a half-filled asymmetric fix (Up=2,Down=1,nBath=2) crashed
// with a near-zero-norm state in the eps-split GS run. Per advisor
// guidance, this was recognized as config-engineering against a
// structurally unreachable target rather than a config to keep searching
// for, and abandoned in favor of THIS test, which validates the actually
// achievable claim: the ramp mechanism itself, applied to a REAL (nonzero)
// first bath, needs no tDMRG-side code.
//
// BandwidthFinal=8 (not 4, matching the base LatticeGf's implied t*=1) is
// chosen so the ramp is genuinely non-trivial: NeqAtomicLimit=1 forces
// t*(0)=0, ramping via the cosine shape toward t*_f=0.25*8=2.0 -- a real
// start-to-finish change, not a same-value no-op.
TEST_CASE("ImpuritySolverNeqTdmrg runs under a cosine hopping ramp "
          "(NeqBathRank=1, nonzero first bath)",
          "[ImpuritySolverNeqTdmrg][Task28][RampSmoke]")
{
	int    argc    = 1;
	char   arg0[]  = "test_ImpuritySolverNeqTdmrg";
	char*  argv0[] = { arg0 };
	char** argv    = argv0;

	ApplicationType app("test_ImpuritySolverNeqTdmrg", &argc, &argv, 1);

	InputNgType::Writeable ioW(Dmft::CincuentaInputCheck {},
	                           configWithU("0.5", "testTdmrgRampSmoke")
	                               + "NeqBathRank=1;\n"
	                                 "QuenchShape=\"cosine\";\n"
	                                 "QuenchDuration=0.25;\n"
	                                 "BandwidthFinal=8;\n"
	                                 "NeqAtomicLimit=1;\n");
	InputNgType::Readable  io(ioW);
	ParamsType             params(io);

	using TdmrgNeqSolverType = Dmft::NeqDmftSolver<ComplexType, Dmft::ImpuritySolverNeqTdmrg>;
	TdmrgNeqSolverType neqSolver(params, app, io);

	// Same symmetric 5-site bath used throughout this file's Phase 1/2
	// gates -- already proven non-degenerate.
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

// ---- Task 28: tDMRG vs GBEK cross-check under a cosine hopping ramp -------
//
// Same cross-check philosophy as the FullGate test above, with the ramp
// params added on top of the identical config/bathParams -- the real
// acceptance gate for Task #28's central claim: if tDMRG's second-bath
// machinery already receives the ramped Lambda correctly (via the same
// solver-agnostic NeqLatticeGf path GBEK already uses), the two solvers
// must still agree to the tolerance the static-bath FullGate test
// established. See the RampSmoke test above for why this does NOT attempt
// GBEK/Wolf's literal nBath=0 atomic limit.
TEST_CASE("ImpuritySolverNeqTdmrg vs ImpuritySolverNeqGBEK agree under a "
          "cosine hopping ramp",
          "[ImpuritySolverNeqTdmrg][Task28][RampGate]")
{
	int    argc    = 1;
	char   arg0[]  = "test_ImpuritySolverNeqTdmrg";
	char*  argv0[] = { arg0 };
	char** argv    = argv0;

	ApplicationType app("test_ImpuritySolverNeqTdmrg", &argc, &argv, 1);

	const std::string rampExtra = "NeqBathRank=1;\n"
	                              "QuenchShape=\"cosine\";\n"
	                              "QuenchDuration=0.25;\n"
	                              "BandwidthFinal=8;\n"
	                              "NeqAtomicLimit=1;\n";

	// Same symmetric 5-site bath used throughout this file's Phase 1/2
	// gates -- already proven non-degenerate.
	const VectorRealType bathParams = { 0.3, 0.3, 0.3, 0.3, 0.3, -0.6, -0.3, 0.0, 0.3, 0.6 };

	InputNgType::Writeable ioWT(Dmft::CincuentaInputCheck {},
	                            configWithU("0.5", "testTdmrgRampGate") + rampExtra);
	InputNgType::Readable  io_t(ioWT);
	ParamsType             paramsT(io_t);
	using TdmrgNeqSolverType = Dmft::NeqDmftSolver<ComplexType, Dmft::ImpuritySolverNeqTdmrg>;
	TdmrgNeqSolverType tdmrgSolver(paramsT, app, io_t);
	tdmrgSolver.solve(bathParams);

	InputNgType::Writeable ioWG(Dmft::CincuentaInputCheck {},
	                            configWithU("0.5", "testTdmrgRampGateGBEK") + rampExtra);
	InputNgType::Readable  io_g(ioWG);
	ParamsType             paramsG(io_g);
	using GbekNeqSolverType = Dmft::NeqDmftSolver<ComplexType, Dmft::ImpuritySolverNeqGBEK>;
	GbekNeqSolverType gbekSolver(paramsG, io_g);
	gbekSolver.solve(bathParams);

	const auto& gimpT = tdmrgSolver.gimp();
	const auto& gimpG = gbekSolver.gimp();

	const int      nT2 = static_cast<int>(paramsT.nT);
	const RealType tol = 1e-4;
	for (int n = 0; n <= nT2; ++n) {
		for (int j = 0; j <= n; ++j) {
			const ComplexType retT = gimpT.retarded(n, j);
			const ComplexType lesT = gimpT.lesser(n, j);
			const ComplexType retG = gimpG.retarded(n, j);
			const ComplexType lesG = gimpG.lesser(n, j);

			std::cout << "[ramp] n=" << n << " j=" << j << " tDMRG G^R=" << retT
			          << " G^<=" << lesT << " GBEK G^R=" << retG << " G^<=" << lesG
			          << "\n";

			CHECK(retT.real() == Catch::Approx(retG.real()).margin(tol));
			CHECK(retT.imag() == Catch::Approx(retG.imag()).margin(tol));
			CHECK(lesT.real() == Catch::Approx(lesG.real()).margin(tol));
			CHECK(lesT.imag() == Catch::Approx(lesG.imag()).margin(tol));
		}
	}
}
