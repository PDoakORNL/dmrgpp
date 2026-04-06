#ifndef IMPURITYSOLVER_BASE_H
#define IMPURITYSOLVER_BASE_H

#include "PsimagLite.h"
#include "Vector.h"

namespace Dmft {

template <typename ComplexOrRealType> class ImpuritySolverBase {

public:

	using RealType          = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using ComplexType       = std::complex<RealType>;
	using VectorRealType    = typename PsimagLite::Vector<RealType>::Type;
	using VectorComplexType = typename PsimagLite::Vector<ComplexType>::Type;
	using ApplicationType   = PsimagLite::PsiApp;

	virtual ~ImpuritySolverBase() { }

	// bathParams[0-nBath-1] ==> V ==> hoppings impurity --> bath
	// bathParams[nBath-...] ==> energies on each bath site
	virtual void solve(const VectorRealType& bathParams) = 0;

	virtual const VectorComplexType& gimp() const = 0;
};
}
#endif // IMPURITYSOLVER_BASE_H
