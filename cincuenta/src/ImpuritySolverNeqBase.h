#ifndef IMPURITYSOLVER_NEQ_BASE_H
#define IMPURITYSOLVER_NEQ_BASE_H
#include "CincuentaInputCheck.h"
#include "InputNg.h"
#include "KadanoffBaym.h"
#include "Matsubaras.h"
#include "ModelParams.h"
#include "Vector.h"

namespace Dmft {

/*!
 * \brief ImpuritySolverNeqBase
 * Abstract interface for non-equilibrium impurity solvers.
 * Separates the two-time G_imp computation from the equilibrium ImpuritySolverBase
 * hierarchy (which returns a 1D Matsubara vector) to avoid polluting that interface.
 */
template <typename ComplexOrRealType> class ImpuritySolverNeqBase {

public:

	using RealType        = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using ComplexType     = std::complex<RealType>;
	using VectorRealType  = typename PsimagLite::Vector<RealType>::Type;
	using KBType          = KadanoffBaym<ComplexOrRealType>;
	using InputNgType     = PsimagLite::InputNg<CincuentaInputCheck>;
	using ModelParamsType = ModelParams<RealType>;

	virtual ~ImpuritySolverNeqBase() { }

	// Run the impurity solver for the given bath parameters.
	// For step-by-step solvers (ExactDiag, Lanczos): diagonalizes H(U_i)/H(U_f) so that
	// subsequent computeGimp(n) calls can fill time slices on demand.
	// For one-shot solvers (tDMRG): runs the full real-time propagation here.
	virtual void solve(const VectorRealType& bathParams) = 0;

	/*!
	 * \brief computeGimp
	 * Fill G_imp KB components for all (t_n, t_j) with j <= n, and G^{Left}(t_n, tau_j).
	 * Called once per time step inside the neq self-consistency loop.
	 *
	 * \param[in/out] gimp Two-time Green's function container to fill.
	 * \param[in] n Current time-step index.
	 */
	virtual void computeGimp(KBType& gimp, int n) const = 0;

	/*!
	 * \brief gimp
	 * Access the last-computed two-time impurity GF.
	 *
	 * \return Const reference to the KadanoffBaym Green's function.
	 */
	virtual const KBType& gimp() const = 0;
};

} // namespace Dmft
#endif // IMPURITYSOLVER_NEQ_BASE_H
