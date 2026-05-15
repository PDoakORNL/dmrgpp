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
	            const FunctionOfFrequencyType& g0,
	            const RealType&                mu,
	            Options                        options)
	    : anderson_function_(nBath, mu)
	    , g0_(g0)
	    , options_(options)
	    , unknowns_(2 * nBath)
	{
		if (options_ == Options::PARTICLE_HOLE_SYMM) {
			unknowns_ = nBath;
		}
	}

	SizeType size() const { return unknowns_; }

	// Returns \sum_n Weights_n |G0cluster(Valpha, eAlpha, iwn) - G0(iwn)|^2
	RealType operator()(const VectorRealType& args) const
	{
		RealType sum             = 0.0;
		SizeType totalMatsubaras = g0_.totalMatsubaras();
		for (SizeType i = 0; i < totalMatsubaras; ++i) {
			ComplexOrRealType iwn(0, g0_.omega(i));
			RealType          weight = 1.0 / g0_.omega(i);
			ComplexOrRealType val    = g0cluster(args, iwn) - g0_(i);
			sum += weight * PsimagLite::real(val * PsimagLite::conj(val));
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
		SizeType totalMatsubaras = g0_.totalMatsubaras();
		SizeType n               = src.size();
		assert(dest.size() == n);

		for (SizeType j = 0; j < n; ++j) {
			RealType sum = 0.0;
			for (SizeType i = 0; i < totalMatsubaras; ++i) {
				ComplexOrRealType iwn(0, g0_.omega(i));
				RealType          weight = 1.0 / g0_.omega(i);
				ComplexOrRealType val    = g0cluster(src, iwn) - g0_(i);

				ComplexOrRealType valPrime = g0clusterPrime(src, iwn, j);

				sum += weight
				    * (PsimagLite::real(val * PsimagLite::conj(valPrime)
				                        + valPrime * PsimagLite::conj(val)));
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

	AndersonFunctionType           anderson_function_;
	const FunctionOfFrequencyType& g0_;
	Options                        options_;
	SizeType                       unknowns_;
};
}
#endif // FITFUNCTION_HPP
