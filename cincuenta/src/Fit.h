#ifndef FIT_H
#define FIT_H
#include "AndersonFunction.h"
#include "FitFunction.hpp"
#include "MersenneTwister.h"
#include "MinParams.h"
#include "Minimizer.h"
#include "PsimagLite.h"

namespace Dmft {

template <typename ComplexOrRealType> class Fit {

public:

	using RealType                = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using VectorRealType          = std::vector<RealType>;
	using MinParamsType           = MinParams<RealType>;
	using AndersonFunctionType    = AndersonFunction<ComplexOrRealType>;
	using FitFunctionType         = FitFunction<ComplexOrRealType>;
	using FunctionOfFrequencyType = typename AndersonFunctionType::FunctionOfFrequencyType;
	using RngType                 = PsimagLite::MersenneTwister;
	using Options                 = typename FitFunctionType::Options;

	struct InitResults {

		InitResults(RealType ra_, RealType rb_, const VectorRealType& result_, bool reset_)
		    : ra(ra_)
		    , rb(rb_)
		    , result(result_)
		    , reset(reset_)
		{ }

		template <typename ReadableType>
		InitResults(ReadableType& io)
		    : ra(0)
		    , rb(0)
		    , reset(false)
		{
			try {
				io.readline(ra, "InitBathRa=");
			} catch (std::exception&) { }

			try {
				io.readline(rb, "InitBathRb=");
			} catch (std::exception&) { }

			try {
				io.read(result, "InitBathVector");
			} catch (std::exception&) { }

			try {
				int tmp = 0;
				io.readline(tmp, "InitBathReset=");
				reset = (tmp > 0);
			} catch (std::exception&) { }

			bool nonConstant = (result.size() > 0);
			bool isConstant  = (ra != 0 || rb != 0);
			if (isConstant && nonConstant)
				err("InitBath: Cannot have ra or rb and also a vector of init "
				    "results\n");

			if (!nonConstant && !isConstant)
				ra = 1;
		}

		RealType       ra;
		RealType       rb;
		VectorRealType result;
		bool           reset;
	};

	Fit(SizeType nBath, const MinParamsType& minParams, const InitResults& initResults)
	    : nBath_(nBath)
	    , minParams_(minParams)
	    , results_(2 * nBath_)
	    , rng_(1234)
	    , initResults_(initResults)
	{ }

	// Compute the optimized bath parameters and store them in vector gammaG
	// See AndersonFunction.h documentation for the fit function, and
	// for the order of storage of bath parameters
	void fit(const FunctionOfFrequencyType& gammaG, const RealType& mu, Options options)
	{
		FitFunctionType f(nBath_, gammaG, mu, options);

		VectorRealType results(f.size());
		setResults(results);

		PsimagLite::Minimizer<RealType, FitFunctionType> min(
		    f, minParams_.maxIter, minParams_.verbose);
		int iter = 0;
		if (minParams_.method == MinParamsType::Method::CONJUGATE_GRADIENT) {
			iter = min.conjugateGradient(
			    results, minParams_.delta, minParams_.delta2, minParams_.tolerance);
		} else {
			iter = min.simplex(results, minParams_.delta, minParams_.tolerance);
		}

		if (iter < 0)
			std::cerr << "No minimum found\n";

		assert(results.size() == nBath_ || results.size() == 2 * nBath_);
		if (results.size() == 2 * nBath_) {
			for (SizeType i = 0; i < 2 * nBath_; ++i) {
				results_[i] = results[i];
			}
		} else {
			// particle-hole symmetric case
			// copy Vs first, starting with the fitted ones...
			for (SizeType i = 0; i < number_of_fitted_Vs; ++i) {
				results_[i] = results[i];
			}

			// ...and then the rest of the Vs are a mirror
			// Works also for Nbath_ odd
			for (SizeType i = number_of_fitted_Vs; i < nBath_; ++i) {
				results_[i] = results[i - number_of_fitted_Vs];
			}

			// copy onsite energies, first the fitted ones...
			SizeType offset1 = nBath_ - number_of_fitted_Vs;
			for (SizeType i = number_of_fitted_Vs; i < nBath_; ++i) {
				results_[i + offset1] = results[i];
			}

			// ...and the rest are the opposites
			SizeType offset2 = 2 * (nbath - number_of_fitted_Vs);
			for (SizeType i = number_of_fitted_Vs; i < nBath_; ++i) {
				results_[i + offset2] = -results[i];
			}

			if (nBath_ & 1) {
				assert(nBath_ + offset2 + 1 == 2 * nBath_);
				results_[2 * nBath_ - 1] = 0; // last onsite if nbath is odd
			} else {
				assert(nBath_ + offset2 == 2 * nBath_);
			}
		}
	}

	const VectorRealType& result() const { return results_; }

	SizeType nBath() const { return nBath_; }

	static Options computeOptions(const std::string& options)
	{
		return FitFunctionType::computeOptions(options);
	}

private:

	void setResults(VectorRealType& results)
	{
		bool nonConstant = (initResults_.result.size() > 0);
		bool isConstant  = (initResults_.ra != 0 || initResults_.rb != 0);
		if (isConstant && nonConstant)
			err("InitResults: Cannot have ra or rb and also a vector of init "
			    "results\n");

		if (nonConstant && initResults_.result.size() != results.size())
			err(PsimagLite::String(
			        "InitResults: vector of init results has wrong size: ")
			    + "expected " + ttos(results.size()) + ", but found "
			    + ttos(initResults_.result.size()) + "\n");

		for (SizeType i = 0; i < results.size(); ++i)
			results[i] = (isConstant) ? initResults_.ra * rng_() + initResults_.rb
			                          : initResults_.result[i];
	}

	const SizeType       nBath_; // number of bath sites
	const MinParamsType& minParams_; // parameters for fitting algorithm
	VectorRealType       results_; // stores bath parameters
	RngType              rng_;
	const InitResults&   initResults_;
};
}
#endif // FIT_H
