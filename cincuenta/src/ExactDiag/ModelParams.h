#ifndef CINC_MODELPARAMS_H
#define CINC_MODELPARAMS_H
#include "Matrix.h"
#include "Vector.h"

namespace Dmft {

template <typename RealType> struct ModelParams {

	using VectorRealType = typename PsimagLite::Vector<RealType>::Type;

	ModelParams(const VectorRealType& bathParams, SizeType center)
	    : center_site(center)
	{
		SizeType bath = bathParams.size() / 2;
		assert((bathParams.size() & 1) == 0);
		sites = bath + 1;
		potentialV.resize(sites);

		// hoppings from the center to any other site
		hoppings.resize(sites - 1);
		for (SizeType i = 0; i < sites; ++i) {
			if (i == center) {
				continue;
				potentialV[center] = 0;
			}

			SizeType offset = (i < center) ? 0 : 1;
			assert(i >= offset);

			hoppings[i]   = bathParams[i - offset];
			potentialV[i] = bathParams[i + bath - offset];
		}
	}

	SizeType       center_site;
	SizeType       sites;
	VectorRealType potentialV;
	VectorRealType hoppings;
};

}
#endif // CINC_MODELPARAMS_H
