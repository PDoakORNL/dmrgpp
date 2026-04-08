#ifndef IMPURITYSOLVER_BASE_H
#define IMPURITYSOLVER_BASE_H

#include "InputNg.h"
#include "Matsubaras.h"
#include "ModelParams.h"
#include "PsimagLite.h"
#include "Vector.h"

namespace Dmft {

template <typename ComplexOrRealType> class ImpuritySolverBase {

public:

	using RealType          = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using ComplexType       = std::complex<RealType>;
	using VectorRealType    = typename PsimagLite::Vector<RealType>::Type;
	using VectorComplexType = typename PsimagLite::Vector<ComplexType>::Type;
	using ApplicationType   = PsimagLite::PsiApp;
	using ModelParamsType   = ModelParams<RealType>;
	using InputNgType       = PsimagLite::InputNg<Dmrg::InputCheck>;

	virtual ~ImpuritySolverBase() { }

	// bathParams[0-nBath-1] ==> V ==> hoppings impurity --> bath
	// bathParams[nBath-...] ==> energies on each bath site
	virtual void solve(const VectorRealType& bathParams) = 0;

	virtual const VectorComplexType& gimp() const = 0;

protected:

	static void writeGimpForDebugOnly(const std::string&                      file_out,
	                                  const std::vector<ComplexOrRealType>&   gimp,
	                                  const PsimagLite::Matsubaras<RealType>& matsubaras)
	{
		const SizeType n = gimp.size();
		std::ofstream  fout(file_out);
		if (!fout || !fout.good())
			err(std::string("Could not write to") + file_out + "\n");

		for (SizeType i = 0; i < n; ++i) {
			const ComplexType value = gimp[i];
			const RealType    omega = matsubaras.omega(i);
			fout << omega << " " << PsimagLite::real(value) << " "
			     << PsimagLite::imag(value) << "\n";
		}

		fout.close();
	}

	static std::string readAndModifyInput(const std::string&     gs_template,
	                                      const ModelParamsType& model_params)
	{
		PsimagLite::String data;
		InputNgType::Writeable::readFile(data, gs_template);
		PsimagLite::String data2 = addBathParams(data, model_params);
		return data2;
	}

	static std::string addBathParams(const std::string&     data,
	                                 const ModelParamsType& model_params)
	{
		const PsimagLite::String connectors = vectorToString(model_params.hoppings(), ",");
		const PsimagLite::String label      = "dir0:Connectors=[" + connectors + "];\n";
		const PsimagLite::String potentialV
		    = vectorToString(model_params.potentialV(), ",");
		const PsimagLite::String label2
		    = "potentialV=[" + potentialV + "," + potentialV + "];\n";

		return data + label + label2;
	}

private:

	static PsimagLite::String vectorToString(const VectorRealType& v,
	                                         const std::string&    separator)
	{
		PsimagLite::String buffer;
		SizeType           n = v.size();
		for (SizeType i = 0; i < n; ++i) {
			buffer += ttos(v[i]);
			if (i + 1 < n) {
				buffer += ",";
			}
		}

		return buffer;
	}
};
}
#endif // IMPURITYSOLVER_BASE_H
