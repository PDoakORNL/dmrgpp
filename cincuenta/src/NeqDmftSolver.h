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

namespace Dmft {

// Non-equilibrium DMFT self-consistency loop for an interaction quench U_i -> U_f.
//
// The loop advances the Kadanoff-Baym equations one real-time step at a time.
// At each step n:
//   1. Compute G_imp(n, j) from the impurity solver (fixed bath).
//   2. Update the hybridization Δ(n, j) = t*² G_imp(n, j) (Bethe self-consistency).
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
		latticeGf_.updateDelta(0, gimp_);
		impSolver_.prepareTimeStep(0, latticeGf_.delta());

		std::cout << "NeqDmftSolver: starting time propagation to t_max=" << params_.tMax
		          << " with nT=" << params_.nT << " steps" << std::endl;

		auto t_neq_start = std::chrono::steady_clock::now();
		for (int n = 1; n <= static_cast<int>(params_.nT); ++n) {
			auto t0 = std::chrono::steady_clock::now();
			timeStep(n);
			auto   t1      = std::chrono::steady_clock::now();
			double dt_wall = std::chrono::duration<double>(t1 - t0).count();
			std::cout << "  step " << n << " / " << params_.nT << "  (" << std::fixed
			          << std::setprecision(1) << dt_wall << " s)" << std::endl;
		}
		double neq_total
		    = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_neq_start)
		          .count();
		std::cout << "NeqDmftSolver: neq phase total " << std::fixed << std::setprecision(1)
		          << neq_total << " s\n";
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
		latticeGf_.delta().dump(p.empty() ? "weiss-delta" : p + "-weiss-delta");
		impSolver_.dumpPlusBath(p.empty() ? "plus-bath-lesser" : p + "-plus-bath-lesser");
	}

private:

	// Advance all KB components by one time step n.
	//
	// Self-consistency following GEBK Fig. 2(b) progressive scheme (PRB 88, 235106):
	//   Predictor: computeGimp uses V[n] from the previous step (extrapolation).
	//   Corrector iterations: updateDelta fills row n of Δ, prepareTimeStep updates
	//   the Cholesky bath V[n] from the complete row, then computeGimp re-evaluates.
	//   For ExactDiag/Lanczos prepareTimeStep is a no-op; NeqDmftIter correctors run.
	//   For GBEK L>0, NeqDmftIter=2-5 converges in 1-3 correctors.
	//   advance(n) runs after the final corrector so G^< is consistent with final V[n].
	void timeStep(int n)
	{
		// Predictor: G_imp(n,j) with V[n] inherited from the previous step.
		impSolver_.computeGimp(gimp_, n);

		for (SizeType iter = 0; iter < params_.neqDmftIter; ++iter) {
			// Δ(n, j) = t*² G_imp(n, j) — fills delta row n
			latticeGf_.updateDelta(n, gimp_);

			// Update bath for step n using the now-complete delta row n.
			// No-op for ExactDiag/Lanczos; updates Cholesky V[n] for GBEK.
			impSolver_.prepareTimeStep(n, latticeGf_.delta());

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
