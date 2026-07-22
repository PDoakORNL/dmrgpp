#ifndef PARAMSDMFTSOLVER_H
#define PARAMSDMFTSOLVER_H
#include "CincuentaInputCheck.h"
#include "MinParams.h"
#include <PsimagLite/Geometry/Star.h>
#include <PsimagLite/InputNg.h>
#include <PsimagLite/Vector.h>

namespace Dmft {

template <typename ComplexOrRealType_> struct ParamsDmftSolver {

	using ComplexOrRealType = ComplexOrRealType_;
	using RealType          = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using MinParamsType     = MinParams<RealType>;
	using InputNgType       = PsimagLite::InputNg<CincuentaInputCheck>;

	ParamsDmftSolver(typename InputNgType::Readable& io)
	    : echoInput(false)
	    , minParams(io)
	{
		io.readline(ficticiousBeta, "FicticiousBeta=");
		io.readline(mu, "ChemicalPotential=");
		io.readline(nMatsubaras, "Matsubaras=");
		io.readline(latticeGf, "LatticeGf=");
		io.readline(nBath, "NumberOfBathPoints=");
		io.readline(dmftIter, "DmftNumberOfIterations=");
		io.readline(dmftError, "DmftTolerance=");
		io.readline(impuritySolver, "ImpuritySolver=");

		try {
			io.readline(precision, "Precision=");
		} catch (std::exception&) { }

		try {
			io.readline(fit_options, "FitOptions=");
		} catch (std::exception&) { }
	}

	bool          echoInput;
	RealType      ficticiousBeta;
	RealType      mu;
	RealType      dmftError;
	SizeType      nMatsubaras;
	std::string   latticeGf;
	std::string   impuritySolver;
	std::string   fit_options;
	SizeType      nBath;
	SizeType      dmftIter;
	SizeType      precision;
	SizeType      center_site = PsimagLite::Star<ComplexOrRealType, int>::CENTER;
	MinParamsType minParams;
};
}
#endif // PARAMSDMFTSOLVER_H
