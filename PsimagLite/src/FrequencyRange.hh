#ifndef PSI_FREQ_RANGE_H
#define PSI_FREQ_RANGE_H
#include "Vector.h"
#include <cassert>

namespace PsimagLite {

template <typename RealType_> class FrequencyRange {
public:

	virtual ~FrequencyRange() { }

	virtual SizeType offset() const { return 0; }

	using RealType = RealType_;

	virtual RealType omega(SizeType i) const = 0;

	virtual SizeType total() const = 0;

	virtual RealType delta() const = 0;
};
}
#endif // PSI_FREQ_RANGE_H
