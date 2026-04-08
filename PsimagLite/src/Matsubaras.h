#ifndef MATSUBARAS_H
#define MATSUBARAS_H
#include "FrequencyRange.hh"
#include "Vector.h"
#include <cassert>

namespace PsimagLite {

template <typename RealType_> class Matsubaras : public FrequencyRange<RealType_> {

public:

	using RealType       = RealType_;
	using VectorRealType = typename PsimagLite::Vector<RealType>::Type;

	Matsubaras(const RealType& fictBeta, SizeType nMatsubara, const RealType& delta)
	    : fictBeta_(fictBeta)
	    , delta_(delta)
	    , nMatsubara_(nMatsubara)
	    , matsubaras_(2 * nMatsubara)
	{
		for (SizeType i = 0; i < nMatsubara_; ++i) {
			matsubaras_[i + nMatsubara_] = M_PI * (2 * i + 1) / fictBeta_;
		}

		for (SizeType i = 0; i < nMatsubara_; ++i) {
			matsubaras_[i] = -matsubaras_[2 * nMatsubara_ - 1 - i];
		}
	}

	RealType omega(SizeType i) const final
	{
		assert(i < matsubaras_.size());
		return matsubaras_[i];
	}

	const RealType& fictitiousBeta() const { return fictBeta_; }

	SizeType total() const final { return matsubaras_.size(); }

	RealType delta() const final { return delta_; }

private:

	RealType       fictBeta_; // ficticious beta
	RealType       delta_;
	SizeType       nMatsubara_; // half the number of matsubaras
	VectorRealType matsubaras_; // wn starting at 0 with the most negative wn
};
}
#endif // MATSUBARAS_H
