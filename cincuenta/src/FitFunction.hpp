#ifndef FITFUNCTION_HPP
#define FITFUNCTION_HPP
#include "AndersonFunction.h"
#include "FunctionOfFrequency.h"

namespace Dmft {

template <typename ComplexOrRealType> class FitFunction {

public:

	using AndersonFunctionType    = AndersonFunction<ComplexOrRealType>;
	using FunctionOfFrequencyType = FunctionOfFrequency<ComplexOrRealType>;
	using RealType                = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using VectorRealType          = typename AndersonFunctionType::VectorRealType;

	using FieldType = RealType;

	FitFunction(SizeType nBath, const FunctionOfFrequencyType& gammaG, const RealType& mu)
	    : anderson_function(nBath, mu)
	    , gammaG_(gammaG)
	{ }

	SizeType size() const { return 2 * anderson_function.nBath(); }

	// Returns \sum_n |Anderson(Valpha, eAlpha, iwn) - GammaG(iwn)|^2
	// See the AndersonFunction below
	RealType operator()(const VectorRealType& args) const
	{
		RealType       sum             = 0.0;
		const SizeType totalMatsubaras = gammaG_.totalMatsubaras();
		for (SizeType i = 0; i < totalMatsubaras; ++i) {
			const ComplexOrRealType iwn(0, gammaG_.omega(i));
			const ComplexOrRealType val
			    = anderson_function.anderson(args, iwn) - gammaG_(i);
			sum += PsimagLite::real(val * PsimagLite::conj(val));
		}

		return sum / totalMatsubaras;
	}

	// For each 0 <= j < 2*nBath, this function
	// returns the derivative of the function above with respect
	// to bath parameter j, evaluated at the bath parameters in src
	// and stores the result in dest[j]
	// for the order of bath parameters see AndersonFunction
	void df(VectorRealType& dest, const VectorRealType& src) const
	{
		const SizeType totalMatsubaras = gammaG_.totalMatsubaras();
		SizeType       nBath           = anderson_function.nBath();

		for (SizeType j = 0; j < 2 * nBath; ++j) {
			RealType sum = 0.0;
			for (SizeType i = 0; i < totalMatsubaras; ++i) {
				const ComplexOrRealType iwn(0, gammaG_.omega(i));
				const ComplexOrRealType val
				    = anderson_function.anderson(src, iwn) - gammaG_(i);

				const ComplexOrRealType valPrime
				    = anderson_function.andersonPrime(src, iwn, j);

				sum += PsimagLite::real(val * PsimagLite::conj(valPrime)
				                        + valPrime * PsimagLite::conj(val));
			}

			dest[j] = sum / totalMatsubaras;
		}
	}

private:

	AndersonFunctionType           anderson_function;
	const FunctionOfFrequencyType& gammaG_;
};
}
#endif // FITFUNCTION_HPP
