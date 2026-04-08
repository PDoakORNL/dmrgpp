#ifndef IMPURITYSOLVER_EXACTD_H
#define IMPURITYSOLVER_EXACTD_H

#include "ContinuedFractionCollection.h"
#include "CrsMatrix.h"
#include "ImpuritySolverBase.h"
#include "InputCheck.h"
#include "LanczosPlusPlus/src/Engine/DefaultSymmetry.h"
#include "LanczosPlusPlus/src/Engine/Engine.h"
#include "LanczosPlusPlus/src/Engine/InternalProductOnTheFly.h"
#include "LanczosPlusPlus/src/Engine/LabeledOperator.h"
#include "LanczosPlusPlus/src/Engine/LanczosGlobals.h"
#include "LanczosPlusPlus/src/Engine/ModelSelector.h"
#include "Matsubaras.h"
#include "ParamsDmftSolver.h"
#include "PsimagLite.h"
#include "Vector.h"

namespace Dmft {

template <typename ComplexOrRealType>
class ImpuritySolverExactDiag : public ImpuritySolverBase<ComplexOrRealType> {

public:

	using BaseType             = ImpuritySolverBase<ComplexOrRealType>;
	using InputNgType          = typename BaseType::InputNgType;
	using ParamsDmftSolverType = ParamsDmftSolver<ComplexOrRealType>;
	using RealType             = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using ComplexType          = std::complex<RealType>;
	using VectorRealType       = typename PsimagLite::Vector<RealType>::Type;
	using VectorComplexType    = typename PsimagLite::Vector<ComplexType>::Type;
	using ApplicationType  = typename ImpuritySolverBase<ParamsDmftSolverType>::ApplicationType;
	using SparseMatrixType = PsimagLite::CrsMatrix<ComplexOrRealType>;
	using SolverParametersType = PsimagLite::ParametersForSolver<RealType>;
	using MatsubarasType       = PsimagLite::Matsubaras<RealType>;
	using GeometryType         = PsimagLite::Geometry<ComplexOrRealType,
	                                                  typename InputNgType::Readable,
	                                                  LanczosPlusPlus::LanczosGlobals>;
	using ModelSelectorType    = LanczosPlusPlus::
	    ModelSelector<ComplexOrRealType, GeometryType, typename InputNgType::Readable>;
	using ModelBaseType = LanczosPlusPlus::
	    ModelBase<ComplexOrRealType, GeometryType, typename InputNgType::Readable>;
	using LanzcosSymmetryType
	    = LanczosPlusPlus::DefaultSymmetry<GeometryType, typename ModelBaseType::BasisBaseType>;
	using InternalProductType
	    = LanczosPlusPlus::InternalProductOnTheFly<ModelBaseType, LanzcosSymmetryType>;
	using EngineType = LanczosPlusPlus::
	    Engine<ModelBaseType, LanczosPlusPlus::InternalProductOnTheFly, LanzcosSymmetryType>;
	using CollectionContFractionType = PsimagLite::ContinuedFractionCollection<RealType>;
	using ModelParamsType            = typename BaseType::ModelParamsType;

	ImpuritySolverExactDiag(const ParamsDmftSolverType& params,
	                        const ApplicationType&,
	                        typename InputNgType::Readable& io)
	    : params_(params)
	    , solverParams_(nullptr)
	    , matsubaras_(params.ficticiousBeta,
	                  params.nMatsubaras,
	                  0.) // last arg. is the real part
	    , gimp_(matsubaras_.total())
	    , io_(io)
	{ }

	~ImpuritySolverExactDiag()
	{
		delete solverParams_;
		solverParams_ = nullptr;
	}

	// bathParams[0-nBath-1] ==> V ==> hoppings impurity --> bath
	// bathParams[nBath-...] ==> energies on each bath site
	void solve(const VectorRealType& bathParams)
	{
		ModelParamsType model_params(bathParams, io_);

		PsimagLite::String data2
		    = BaseType::readAndModifyInput(params_.gsTemplate, model_params);

		// This will replaced by a LanzosRunner at some point
		// so that we don't repeat what's in lanczos.cpp main driver
		// and also honor and have all its features

		// std::cout << geometry;

		Dmrg::InputCheck                inputCheck;
		typename InputNgType::Writeable ioWriteable(inputCheck, data2);
		typename InputNgType::Readable  io(ioWriteable);

		GeometryType         geometry(io);
		ModelSelectorType    modelSelector(io, geometry);
		const ModelBaseType& modelPtr = modelSelector();
		EngineType           engine(modelPtr, io);

		RealType energy = engine.energies(0);
		std::cout << "Energy=" << energy << "\n";

		// We need a collection of cont. fraction because we have two in the formula
		// for the Green Function
		CollectionContFractionType       cfCollection(PsimagLite::FreqEnum::FREQ_MATSUBARA);
		LanczosPlusPlus::LabeledOperator OPERATOR_C
		    = LanczosPlusPlus::LabeledOperator::Label::OPERATOR_C;
		SizeType                      impurity_site = model_params.impuritySite();
		std::pair<SizeType, SizeType> spin_pair(0, 0); // spin symmetric, we choose up == 0
		std::pair<SizeType, SizeType> orb_pair(0, 0); // no orbitals
		std::vector<std::string>      vstr; // output
		engine.spectralFunction(cfCollection,
		                        vstr,
		                        OPERATOR_C,
		                        impurity_site,
		                        impurity_site,
		                        spin_pair,
		                        orb_pair);

		// compute gimp
		// computeGreenFunction(energy, gs, basis, mp);
		std::cout << "SumOfGimp=" << density() << "\n";
	}

	const VectorComplexType& gimp() const { return gimp_; }

private:

	ComplexType density() const
	{
		const SizeType n   = gimp_.size();
		ComplexType    sum = 0;
		for (SizeType i = 0; i < n; ++i)
			sum += gimp_[i];

		return sum;
	}

	void computeGreenFunction(const CollectionContFractionType& cf_collection)
	{
		typename CollectionContFractionType::PlotDataType   results;
		typename CollectionContFractionType::PlotParamsType params;
		params.beta               = params_.ficticiousBeta;
		params.numberOfMatsubaras = matsubaras_.size();

		cf_collection.plot(results, params);
	}

	const ParamsDmftSolverType&     params_;
	SolverParametersType*           solverParams_;
	MatsubarasType                  matsubaras_;
	VectorComplexType               gimp_;
	typename InputNgType::Readable& io_;
};
}
#endif // IMPURITYSOLVER_EXACTD_H
