#include "CincuentaInputCheck.h"
#include "ImpuritySolverNeqTdmrg.h"
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

using RealType        = double;
using ComplexType     = std::complex<RealType>;
using SolverType      = Dmft::ImpuritySolverNeqTdmrg<ComplexType>;
using ParamsType      = Dmft::ParamsNeqDmftSolver<ComplexType>;
using InputNgType     = PsimagLite::InputNg<Dmft::CincuentaInputCheck>;
using VectorRealType  = typename SolverType::VectorRealType;
using ApplicationType = PsimagLite::PsiApp;

// N=6 sites (5 bath + 1 impurity): TSPAdvanceEach=N-2=4, matching the
// convention already established in inputU0NeqTdmrg.ain. NtNeq kept tiny
// (2) to keep this fast -- this is a mechanics check, not a physics
// benchmark. %U% is substituted with the desired HubbardU/HubbardUFinal
// value (no quench in U either variant -- only U itself is varied between
// the two test cases, everything else identical).
static std::string configWithU(const std::string& uValue)
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
	      "RootOutputname=\"testTdmrgChain\";\n"
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
