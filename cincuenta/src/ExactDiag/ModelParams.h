#ifndef CINC_MODELPARAMS_H
#define CINC_MODELPARAMS_H
#include "Geometry/Star.h"
#include "InputCheck.h"
#include "Matrix.h"
#include "Vector.h"

namespace Dmft {

template <typename ComplexOrRealType> struct ModelParams {

	using RealType       = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using VectorRealType = typename PsimagLite::Vector<RealType>::Type;
	using InputNgType    = PsimagLite::InputNg<Dmrg::InputCheck>;
	using StarType       = PsimagLite::Star<ComplexOrRealType, InputNgType::Readable>;

	ModelParams(const VectorRealType& bathParams, InputNgType::Readable& io)
	{
		io.readline(nsites_, "TotalNumberOfSites=");
		star_ = StarType(nsites_, io);

		SizeType bath = bathParams.size() / 2;
		assert((bathParams.size() & 1) == 0);
		assert(nsites_ == bath + 1);
		potentialV_.resize(nsites_, 0);

		// hoppings from the center to any other site
		hoppings_.resize(nsites_ - 1);
		SizeType center     = star_.CENTER;
		potentialV_[center] = 0; // potential at the impurity is set to zero
		for (SizeType i = 0; i < nsites_; ++i) {
			if (i == center) {
				continue;
			}

			SizeType handle = star_.handle(center, i);

			assert(handle < hoppings_.size());
			assert(handle + bath < bathParams.size());
			hoppings_[handle] = bathParams[handle];

			potentialV_[i] = bathParams[handle + bath];
		}
	}

	const PsimagLite::GeometryBase<ComplexOrRealType, InputNgType::Readable>& geometry() const
	{
		return star_;
	}

	SizeType numberOfSites() const { return nsites_; }

	const VectorRealType& potentialV() const { return potentialV_; }

	const VectorRealType& hoppings() const { return hoppings_; }

	SizeType impuritySite() const { return star_.CENTER; }

private:

	SizeType       nsites_;
	VectorRealType potentialV_;
	VectorRealType hoppings_;
	StarType       star_;
};

}
#endif // CINC_MODELPARAMS_H
