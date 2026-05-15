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

	SizeType nBath() const { return nBath_; }

	const RealType& mu() const { return mu_; }

	// Returns \sum_{0<=j<nBath} V_j^2/(iwn - epsilon_j),
	// where the V_j are stored in the first half or args,
	// and the epsilon_j are stored in the last half or args
	ComplexOrRealType anderson(const VectorRealType& args, ComplexOrRealType iwn) const
	{
		assert(PsimagLite::real(iwn) == 0);
		ComplexOrRealType sum = 0.0;
		for (SizeType i = 0; i < nBath_; ++i) {
			const RealType valpha  = (args.size() == 2 * nBath_)
			     ? args[i]
			     : calcVsIfParticleHoleSymm(args, i);
			const RealType epsilon = (args.size() == 2 * nBath_)
			    ? args[i + nBath_]
			    : calcEpsilonIfParticleHoleSymm(args, i);
			sum += valpha * valpha / (iwn + mu_ - epsilon);
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
		assert(jnd < args.size());
		RealType valpha      = 0;
		RealType epsilon     = 0;
		bool     diff_valpha = false;
		if (args.size() == 2 * nBath_) {
			valpha      = (jnd < nBath_) ? args[jnd] : args[jnd - nBath_];
			epsilon     = (jnd < nBath_) ? args[jnd + nBath_] : args[jnd];
			diff_valpha = (jnd < nBath_); // diff with respect to valpha
		} else {
			assert(args.size() == (nBath_ & 1) ? nBath_ + (nBath_ - 1) / 2
			                                   : nBath_ + nBath_ / 2);
			SizeType knd = (jnd < nBath_) ? jnd : jnd - nBath_;
			valpha       = calcVsIfParticleHoleSymm(args, knd);
			epsilon      = calcEpsilonIfParticleHoleSymm(args, knd);
			diff_valpha  = (jnd < nBath_);
		}

		return (diff_valpha) ? 2.0 * valpha / (iwn + mu_ - epsilon)
		                     : squareOf(valpha / (iwn + mu_ - epsilon));
	}

private:

	RealType calcEpsilonIfParticleHoleSymm(const VectorRealType& args, SizeType ind) const
	{
		assert(ind < nBath_);
		SizeType unknowns = (nBath_ & 1) ? nBath_ + (nBath_ - 1) / 2 : nBath_ + nBath_ / 2;
		assert(args.size() == unknowns);
		SizeType offset = (nBath_ & 1) ? (nBath_ - 1) / 2 : nBath_ / 2;

		if (ind < offset) {
			return args[ind + nBath_];
		} else {
			if (nBath_ & 1) {
				return (ind == offset) ? 0 : -args[ind + offset];
			} else {
				return -args[ind + offset];
			}
		}
	}

	RealType calcVsIfParticleHoleSymm(const VectorRealType& args, SizeType ind) const
	{
		assert(ind < nBath_);
		assert(args.size() == (nBath_ & 1) ? nBath_ + (nBath_ - 1) / 2
		                                   : nBath_ + nBath_ / 2);

		return args[ind];
	}

	SizeType nBath_; // Number of Bath sites
	RealType mu_; // the Chemical potential
};

}
#endif // ANDERSONFUNCTION_H
