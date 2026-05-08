#ifndef ANDERSONFUNCTION_H
#define ANDERSONFUNCTION_H
#include "FunctionOfFrequency.h"
#include "Vector.h"

namespace Dmft {

template <typename ComplexOrRealType> class AndersonFunction {

public:

	using RealType                = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using VectorRealType          = typename PsimagLite::Vector<RealType>::Type;
	using FunctionOfFrequencyType = FunctionOfFrequency<ComplexOrRealType>;

	using FieldType = RealType;

	AndersonFunction(SizeType nBath, const RealType& mu)
	    : nBath_(nBath)
	    , mu_(mu)
	{ }

	ComplexOrRealType anderson(const VectorRealType& args, ComplexOrRealType iwn) const
	{
		return anderson(args, iwn, nBath_, mu_);
	}

	SizeType nBath() const { return nBath_; }

	// Returns \sum_{0<=j<nBath} V_j^2/(iwn - epsilon_j),
	// where the V_j are stored in the first half or args,
	// and the epsilon_j are stored in the last half or args
	static ComplexOrRealType anderson(const VectorRealType& args,
	                                  ComplexOrRealType     iwn,
	                                  SizeType              nBath,
	                                  const RealType&       mu)
	{
		assert(args.size() == 2 * nBath || args.size() == nBath);
		assert(PsimagLite::real(iwn) == 0);
		ComplexOrRealType sum = 0.0;
		for (SizeType i = 0; i < nBath; ++i) {
			const RealType valpha  = args[i];
			const RealType epsilon = (args.size() == nBath) ? 0 : args[i + nBath];
			sum += valpha * valpha / (iwn + mu - epsilon);
		}

		return sum;
	}

	static ComplexOrRealType squareOf(ComplexOrRealType x) { return x * x; }

	// For any 0 <= jnd < 2*nBath, this function returns the derivative of
	// the AndersonFunction above with respect to bath parameter jnd,
	// evaluated at the bath parameters args
	// The order in which bath parameters are stored is described
	// under AndersonFunction
	ComplexOrRealType
	andersonPrime(const VectorRealType& args, ComplexOrRealType iwn, SizeType jnd) const
	{
		assert(args.size() == 2 * nBath_ || args.size() == nBath_);
		assert(jnd < args.size());
		RealType valpha  = (jnd < nBath_) ? args[jnd] : args[jnd - nBath_];
		RealType epsilon = 0;
		if (args.size() == 2 * nBath_) {
			epsilon = (jnd < nBath_) ? args[jnd + nBath_] : args[jnd];
		}

		return (jnd < nBath_) ? 2.0 * valpha / (iwn + mu_ - epsilon)
		                      : squareOf(valpha / (iwn + mu_ - epsilon));
	}

private:

	SizeType nBath_; // Number of Bath sites
	RealType mu_; // the Chemical potential
};

}
#endif // ANDERSONFUNCTION_H
