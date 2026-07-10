#ifndef ANDERSONFUNCTION_H
#define ANDERSONFUNCTION_H
#include "FunctionOfFrequency.h"
#include <PsimagLite/Vector.h>

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

	RealType mu() const { return mu_; }

	SizeType nBath() const { return nBath_; }

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
			     : calcVsIfParticleHoleSymm(args, i, nBath_);
			const RealType epsilon = (args.size() == 2 * nBath_)
			    ? args[i + nBath_]
			    : calcEpsilonIfParticleHoleSymm(args, i, nBath_);
			sum += valpha * valpha / (iwn + mu_ - epsilon);
		}

		return sum;
	}

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
			assert(args.size() == nBath_);
			SizeType offset  = (nBath_ & 1) ? (nBath_ + 1) / 2 : nBath_ / 2;
			SizeType Eoffset = (nBath_ & 1) ? (nBath_ - 1) / 2 : nBath_ / 2;
			SizeType knd     = (jnd < offset) ? jnd : jnd - offset;
			valpha           = calcVsIfParticleHoleSymm(args, knd, nBath_);
			epsilon          = calcEpsilonIfParticleHoleSymm(args, knd, nBath_);
			diff_valpha      = (jnd < offset);
			// Each V or epsilon parameter is shared by two mirror bath sites;
			// both contributions must be summed.
			if (diff_valpha) {
				if (knd < Eoffset) // paired baths: epsilon and -epsilon
					return ComplexOrRealType(2.0) * valpha
					    / (iwn + mu_ - epsilon)
					    + ComplexOrRealType(2.0) * valpha
					    / (iwn + mu_ + epsilon);
				return ComplexOrRealType(2.0) * valpha
				    / (iwn + mu_); // middle bath, epsilon=0
			}
			return squareOf(valpha / (iwn + mu_ - epsilon)) // epsilon derivative:
			    - squareOf(valpha
			               / (iwn + mu_ + epsilon)); // opposite sign for mirror bath
		}

		return (diff_valpha) ? 2.0 * valpha / (iwn + mu_ - epsilon)
		                     : squareOf(valpha / (iwn + mu_ - epsilon));
	}

	static ComplexOrRealType squareOf(ComplexOrRealType x) { return x * x; }

private:

	static RealType
	calcEpsilonIfParticleHoleSymm(const VectorRealType& args, SizeType ind, SizeType nBath)
	{
		assert(args.size() == nBath);
		SizeType offset      = (nBath & 1) ? (nBath - 1) / 2 : nBath / 2;
		SizeType one_or_zero = (nBath & 1);

		if (ind < offset) {
			return args[ind + offset + one_or_zero];
		} else {
			if (nBath & 1) {
				return (ind == offset) ? 0 : -args[ind];
			} else {
				return -args[ind];
			}
		}
	}

	static RealType
	calcVsIfParticleHoleSymm(const VectorRealType& args, SizeType ind, SizeType nBath)
	{
		assert(args.size() == nBath);
		SizeType offset = (nBath & 1) ? (nBath + 1) / 2 : nBath / 2;

		if (ind < offset) {
			return args[ind];
		} else {
			return args[ind - offset];
		}
	}

	SizeType nBath_; // Number of Bath sites
	RealType mu_; // the Chemical potential
};

}
#endif // ANDERSONFUNCTION_H
