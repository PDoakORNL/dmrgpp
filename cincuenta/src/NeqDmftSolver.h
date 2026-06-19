#ifndef NEQ_DMFT_SOLVER_H
#define NEQ_DMFT_SOLVER_H

#include "CincuentaInputCheck.h"
#include "ImpuritySolverNeqExactDiag.h"
#include "KadanoffBaym.h"
#include "NeqLatticeGf.h"
#include "ParamsNeqDmftSolver.h"
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
//   ImpuritySolverNeqLanczos       — truncated Lanczos Lehmann (larger baths)
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
		std::cout << "NeqDmftSolver: running impurity solver setup\n";
		impSolver_.solve(bathParams);

		// Copy equilibrium Matsubara components from the solver's internal gimp —
		// computeGimp() only fills real-time (retarded/lesser/left-mixing) slices.
		gimp_.matsubara_t = impSolver_.gimp().matsubara_t;
		gimp_.matsubara_w = impSolver_.gimp().matsubara_w;

		// Populate t=0 (equilibrium) boundary conditions.
		impSolver_.computeGimp(gimp_, 0);
		latticeGf_.initialize(gimp_);

		std::cout << "NeqDmftSolver: starting time propagation to t_max=" << params_.tMax
		          << " with nT=" << params_.nT << " steps\n";

		for (int n = 1; n <= static_cast<int>(params_.nT); ++n) {
			timeStep(n);
			if (n % 10 == 0 || n == static_cast<int>(params_.nT))
				std::cout << "  step " << n << " / " << params_.nT << "\n";
		}

		std::cout << "NeqDmftSolver: done\n";
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
	}

private:

	// Advance all KB components by one time step n.
	// Inner DMFT iterations are performed until convergence or maxIter is reached.
	void timeStep(int n)
	{
		// Step 1: compute G_imp(n, j) for j ≤ n from the Lehmann representation.
		impSolver_.computeGimp(gimp_, n);

		// Inner self-consistency loop (for future use with a self-consistent bath).
		// With the fixed-bath ExactDiag solver the loop always converges in 1 iter.
		for (SizeType iter = 0; iter < params_.neqDmftIter; ++iter) {

			// Step 2: Δ(n, j) = t*² G_imp(n, j)
			latticeGf_.updateDelta(n, gimp_);

			// Step 3: advance G_0(n, j) via Volterra integro-differential equation.
			latticeGf_.advance(n);

			// Convergence check placeholder.
			// A full self-consistent loop would compare the new G_0 against the
			// previous iteration and break when the difference falls below
			// params_.neqDmftError.  For the fixed-bath ExactDiag solver, one
			// iteration per time step is exact.
			break;
		}
	}

	const ParamsNeqType& params_;
	ImpSolverType        impSolver_;
	LatticeGfType        latticeGf_;
	KBType               gimp_; // local copy filled step by step
};

} // namespace Dmft
#endif // NEQ_DMFT_SOLVER_H
