#ifndef CLUSTERG__HPP
#define CLUSTERG__HPP
#include "AndersonFunction.h"
#include <vector>

namespace Dmft {

template <typename ComplexOrRealType> class ClusterG0 {
public:

	using AnderFunctionType = AndersonFunction<ComplexOrRealType>;
	using RealType          = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using VectorRealType    = std::vector<RealType>;

	ClusterG0(const AnderFunctionType& anderson_function)
	    : anderson_function_(anderson_function)
	{ }

	ComplexOrRealType g0cluster(const VectorRealType& args, const ComplexOrRealType& iwn) const
	{
		constexpr RealType epsilon0 = 0;
		RealType           mu       = anderson_function_.mu();
		ComplexOrRealType  tmp
		    = iwn + mu - epsilon0 - anderson_function_.anderson(args, iwn);
		return 1.0 / tmp;
	}

	ComplexOrRealType
	g0clusterPrime(const VectorRealType& args, const ComplexOrRealType& iwn, SizeType jnd) const
	{
		constexpr RealType epsilon0 = 0;
		RealType           mu       = anderson_function_.mu();
		ComplexOrRealType  tmp
		    = iwn + mu - epsilon0 - anderson_function_.anderson(args, iwn);
		ComplexOrRealType derivative_tmp
		    = -anderson_function_.andersonPrime(args, iwn, jnd);
		return -derivative_tmp / (tmp * tmp);
	}

private:

	const AnderFunctionType& anderson_function_;
};
}
#endif // CLUSTERG__HPP
