#ifndef FITFUNCTION_HPP
#define FITFUNCTION_HPP
#include "AndersonFunction.h"
#include "FunctionOfFrequency.h"
#include "ProgramGlobals.h"

namespace Dmft {

template <typename ComplexOrRealType> class FitFunction {

public:

	enum class Options
	{
		NONE,
		PARTICLE_HOLE_SYMM
	};

	using AndersonFunctionType    = AndersonFunction<ComplexOrRealType>;
	using FunctionOfFrequencyType = FunctionOfFrequency<ComplexOrRealType>;
	using RealType                = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using VectorRealType          = typename AndersonFunctionType::VectorRealType;

	using FieldType = RealType;

	// (FIXME document it here)
	FitFunction(SizeType                       nBath,
	            const FunctionOfFrequencyType& gammaG,
	            const RealType&                mu,
	            Options                        options)
	    : anderson_function(nBath, mu)
	    , gammaG_(gammaG)
	    , options_(options)
	    , unknowns_(2 * nBath)
	{
		if (options_ == Options::PARTICLE_HOLE_SYMM) {
			SizeType unknown_energies = (nBath & 1) ? (nBath - 1) / 2 : nBath / 2;
			unknowns_                 = nBath + unknown_energies;
		}
	}

	SizeType size() const { return unknowns_; }

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
		SizeType       n               = src.size();
		assert(dest.size() == n);

		for (SizeType j = 0; j < n; ++j) {
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

	static Options computeOptions(const std::string& str)
	{
		std::string optionslc = Dmrg::ProgramGlobals::toLower(str);
		if (optionslc.empty()) {
			return Options::NONE;
		} else if (optionslc == "particleholesymmetric") {
			return Options::PARTICLE_HOLE_SYMM;
		} else {
			throw std::runtime_error("FitFunction: unknown option" + str + "\n");
		}
	}

private:

	AndersonFunctionType           anderson_function;
	const FunctionOfFrequencyType& gammaG_;
	Options                        options_;
	SizeType                       unknowns_;
};
}
#endif // FITFUNCTION_HPP
