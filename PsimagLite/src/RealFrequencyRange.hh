#ifndef REAL_FREQ_RANGE_H
#define REAL_FREQ_RANGE_H
#include "FrequencyRange.hh"
#include "Vector.h"
#include <cassert>

namespace PsimagLite {

template <typename RealType_> class RealFrequencyRange : public FrequencyRange<RealType_> {

public:

	using RealType       = RealType_;
	using VectorRealType = typename PsimagLite::Vector<RealType>::Type;

	RealFrequencyRange(const RealType& begin,
	                   const RealType& step,
	                   SizeType        total,
	                   const RealType& delta)
	    : begin_(begin)
	    , step_(step)
	    , delta_(delta)
	    , total_(total)
	{ }

	RealType omega(SizeType i) const { return begin_ + step_ * i; }

	SizeType total() const { return total_; }

	RealType delta() const { return delta_; }

private:

	RealType begin_;
	RealType step_;
	RealType delta_;
	SizeType total_;
};
}
#endif // REAL_FREQ_RANGE_H
