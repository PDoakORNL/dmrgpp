#ifndef IMPURITYSOLVER_NEQ_BASE_H
#define IMPURITYSOLVER_NEQ_BASE_H
#include "CincuentaInputCheck.h"
#include "InputNg.h"
#include "KadanoffBaym.h"
#include "Matsubaras.h"
#include "ModelParams.h"
#include "Vector.h"

namespace Dmft {

// Abstract interface for non-equilibrium impurity solvers.
// Separates the two-time G_imp computation from the equilibrium ImpuritySolverBase
// hierarchy (which returns a 1D Matsubara vector) to avoid polluting that interface.
template <typename ComplexOrRealType>
class ImpuritySolverNeqBase {

public:

	using RealType          = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using ComplexType       = std::complex<RealType>;
	using VectorRealType    = typename PsimagLite::Vector<RealType>::Type;
	using KBType            = KadanoffBaym<ComplexOrRealType>;
	using InputNgType       = PsimagLite::InputNg<CincuentaInputCheck>;
	using ModelParamsType   = ModelParams<RealType>;

	virtual ~ImpuritySolverNeqBase() {}

	// One-time setup: full diagonalization of H(U_i) and H(U_f) for the given bath.
	virtual void initialize(const VectorRealType& bathParams) = 0;

	// Fill G_imp KB components for all (t_n, t_j) with j <= n, and G^{Left}(t_n, tau_j).
	// Called once per time step inside the neq self-consistency loop.
	virtual void computeGimp(KBType& gimp, int n) const = 0;

	// Access the last-computed two-time impurity GF.
	virtual const KBType& gimp() const = 0;
};

} // namespace Dmft
#endif // IMPURITYSOLVER_NEQ_BASE_H
