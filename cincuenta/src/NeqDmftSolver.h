#ifndef NEQ_DMFT_SOLVER_H
#define NEQ_DMFT_SOLVER_H

#include "CincuentaInputCheck.h"
#include "ImpuritySolverNeqExactDiag.h"
#include "KadanoffBaym.h"
#include "NeqLatticeGf.h"
#include "ParamsNeqDmftSolver.h"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace Dmft {

// Non-equilibrium DMFT self-consistency loop for an interaction quench U_i -> U_f.
//
// The loop advances the Kadanoff-Baym equations one real-time step at a time.
// At each step n:
//   1. Compute G_imp(n, j) from the impurity solver (fixed bath).
//   2. Update the hybridization Λ(n, j) = t*² G_imp(n, j) (Bethe self-consistency).
//   3. Advance the Weiss field G_0(n, j) via the Volterra integro-differential equation.
//
// ImpSolverTemplate selects the impurity solver:
//   ImpuritySolverNeqExactDiag  — full Lehmann (default, exact for small baths)
//   ImpuritySolverNeqLanczos    — truncated Lanczos Lehmann (larger baths)
//   ImpuritySolverNeqGBEK       — two-bath Cholesky scheme (GEBK PRB 88, 235106)
template <typename ComplexOrRealType,
          template <typename> class ImpSolverTemplate = ImpuritySolverNeqExactDiag>
class NeqDmftSolver {

public:

	using RealType          = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using ComplexType       = std::complex<RealType>;
	using VectorRealType    = typename PsimagLite::Vector<RealType>::Type;
	using VectorComplexType = typename PsimagLite::Vector<ComplexType>::Type;
	using KBType            = KadanoffBaym<ComplexOrRealType>;
	using KBDerivType       = KBDerivative<ComplexOrRealType>;
	using ParamsNeqType     = ParamsNeqDmftSolver<ComplexOrRealType>;
	using InputNgType       = PsimagLite::InputNg<CincuentaInputCheck>;
	using ImpSolverType     = ImpSolverTemplate<ComplexOrRealType>;
	using LatticeGfType     = NeqLatticeGf<ComplexOrRealType>;
	using ApplicationType   = PsimagLite::PsiApp;

	NeqDmftSolver(const ParamsNeqType& params, typename InputNgType::Readable& io)
	    : params_(params)
	    , impSolver_(params, io)
	    , latticeGf_(params)
	    , gimp_(params.nT,
	            params.eqParams.nMatsubaras,
	            params.dt,
	            params.eqParams.ficticiousBeta
	                / static_cast<RealType>(params.eqParams.nMatsubaras))
	{ }

	// Overload for impurity solvers whose constructor needs an extra
	// ApplicationType& (currently only ImpuritySolverNeqTdmrg, which drives
	// DmrgRunner subprocesses and so needs the PsiApp). Only instantiated
	// (and so only needs to compile) when ImpSolverType actually has a
	// matching 3-arg constructor -- member functions of a class template
	// are instantiated lazily on use, so this does not affect
	// ImpuritySolverNeqGBEK/ExactDiag/Lanczos, none of which have one.
	// Added to let ImpuritySolverNeqTdmrg be driven through NeqDmftSolver
	// like every other impurity solver (see fancy-painting-moon.md, Phase
	// 2 gate task) -- previously it could only be constructed directly,
	// bypassing NeqDmftSolver's self-consistency loop entirely.
	NeqDmftSolver(const ParamsNeqType&            params,
	              const ApplicationType&          app,
	              typename InputNgType::Readable& io)
	    : params_(params)
	    , impSolver_(params, app, io)
	    , latticeGf_(params)
	    , gimp_(params.nT,
	            params.eqParams.nMatsubaras,
	            params.dt,
	            params.eqParams.ficticiousBeta
	                / static_cast<RealType>(params.eqParams.nMatsubaras))
	{ }

	// Run the full neq-DMFT calculation.
	//   bathParams: equilibrium bath parameters {V_0..V_{nBath-1}, ε_0..ε_{nBath-1}}
	//               obtained from the preceding equilibrium DMFT run.
	void solve(const VectorRealType& bathParams)
	{
		std::cout << "NeqDmftSolver: running impurity solver setup" << std::endl;
		impSolver_.solve(bathParams);

		// Copy equilibrium Matsubara components from the solver's internal gimp —
		// computeGimp() only fills real-time (retarded/lesser/left-mixing) slices.
		gimp_.matsubara_t = impSolver_.gimp().matsubara_t;
		gimp_.matsubara_w = impSolver_.gimp().matsubara_w;

		// Populate t=0 (equilibrium) boundary conditions.
		impSolver_.computeGimp(gimp_, 0);
		latticeGf_.initialize(gimp_);
		// Seed the Cholesky decomposition for step 0.
		// prepareTimeStep(n) is called in the n>=1 loop; n=0 must be primed here.
		latticeGf_.updateLambda(0, gimp_);
		impSolver_.prepareTimeStep(0, latticeGf_.lambda());

		std::cout << "NeqDmftSolver: starting time propagation to t_max=" << params_.tMax
		          << " with nT=" << params_.nT << " steps" << std::endl;

		auto t_neq_start = std::chrono::steady_clock::now();
		for (int n = 1; n <= static_cast<int>(params_.nT); ++n) {
			auto t0 = std::chrono::steady_clock::now();
			timeStep(n);
			auto   t1      = std::chrono::steady_clock::now();
			double dt_wall = std::chrono::duration<double>(t1 - t0).count();
			// Format into a LOCAL stream, not std::cout directly: std::fixed/
			// setprecision are sticky on the stream object they're applied
			// to, and DmrgRunner's own in-situ measurement prints
			// (TargetingCommon::test()) write to this same global std::cout
			// (redirected per-run to different log files, but the format
			// state persists across redirects) -- setting precision(1) here
			// directly previously truncated every subsequent measurement log
			// in the process to 1 decimal digit, corrupting Task 15's gate
			// comparison. See TDMRG_EVOLVING_BATH.md Link 12.
			std::ostringstream dtStr;
			dtStr << std::fixed << std::setprecision(1) << dt_wall;
			std::cout << "  step " << n << " / " << params_.nT << "  (" << dtStr.str()
			          << " s)" << std::endl;
		}
		double neq_total
		    = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_neq_start)
		          .count();
		std::ostringstream neqTotalStr;
		neqTotalStr << std::fixed << std::setprecision(1) << neq_total;
		std::cout << "NeqDmftSolver: neq phase total " << neqTotalStr.str() << " s\n";
	}

	// Access the impurity GF (populated after solve()).
	const KBType& gimp() const { return gimp_; }

	// Access the Weiss field G_0 (populated after solve()).
	const KBType& g0() const { return latticeGf_.g0(); }

	// Write KB Green's functions to files.  When params_.neqOutputPrefix is set,
	// filenames are "{prefix}-green-retarded" etc.; otherwise "green-retarded" etc.
	void dumpGreenFunctions() const
	{
		const std::string& p = params_.neqOutputPrefix;
		gimp_.dump(p.empty() ? "green" : p + "-green");
		latticeGf_.g0().dump(p.empty() ? "weiss-green" : p + "-weiss-green");
		// Filename kept as "weiss-delta" (external tooling convention, ~20
		// gbek_reference/*.py scripts hardcode this path) even though the
		// accessor below is lambda() -- the quantity is GBEK/Wolf's Λ(t,t'),
		// not the equilibrium hybridization Δ; see NeqLatticeGf.h.
		latticeGf_.lambda().dump(p.empty() ? "weiss-delta" : p + "-weiss-delta");
		impSolver_.dumpPlusBath(p.empty() ? "plus-bath-lesser" : p + "-plus-bath-lesser");
		impSolver_.dumpV(p.empty() ? "cholesky-V" : p + "-cholesky-V");
		impSolver_.dumpDoccAndEnergy(p.empty() ? "docc-energy" : p + "-docc-energy");
	}

private:

	// Advance all KB components by one time step n.
	//
	// Self-consistency following GEBK Fig. 2(b) progressive scheme (PRB 88, 235106):
	//   Predictor: computeGimp uses V[n] from the previous step (extrapolation).
	//   Corrector iterations: updateLambda fills row n of Λ, prepareTimeStep updates
	//   the Cholesky bath V[n] from the complete row, then computeGimp re-evaluates.
	//   For ExactDiag/Lanczos prepareTimeStep is a no-op; NeqDmftIter correctors run.
	//   For GBEK L>0, NeqDmftIter=2-5 converges in 1-3 correctors.
	//   advance(n) runs after the final corrector so G^< is consistent with final V[n].
	void timeStep(int n)
	{
		// Predictor: G_imp(n,j) with V[n] inherited from the previous step.
		impSolver_.computeGimp(gimp_, n);

		for (SizeType iter = 0; iter < params_.neqDmftIter; ++iter) {
			// Λ(n, j) = t*² G_imp(n, j) — fills lambda row n
			latticeGf_.updateLambda(n, gimp_);

			// Update bath for step n using the now-complete lambda row n.
			// No-op for ExactDiag/Lanczos; updates Cholesky V[n] for GBEK.
			impSolver_.prepareTimeStep(n, latticeGf_.lambda());

			// Corrector: re-evaluate G_imp with the updated bath.
			// Always called so that advance(n) sees G^< computed with the
			// final V[n], not the penultimate one.
			impSolver_.computeGimp(gimp_, n);
		}

		// Advance G_0(n, j) via Volterra integro-differential equation.
		latticeGf_.advance(n);
	}

	const ParamsNeqType& params_;
	ImpSolverType        impSolver_;
	LatticeGfType        latticeGf_;
	KBType               gimp_; // local copy filled step by step
};

} // namespace Dmft
#endif // NEQ_DMFT_SOLVER_H
