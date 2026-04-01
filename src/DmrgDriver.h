#ifndef DMRGDRIVER_H
#define DMRGDRIVER_H
#define USE_PTHREADS_OR_NOT_NG
#include "BasisWithOperators.h"
#include "CanonicalExpression.h"
#include "ChebyshevSolver.h"
#include "CrsMatrix.h"
#include "DmrgSolver.h"
#include "InputCheck.h"
#include "InputNg.h"
#include "LanczosSolver.h"
#include "LeftRightSuper.h"
#include "MatrixVectorKron/MatrixVectorKron.h"
#include "MatrixVectorOnTheFly.h"
#include "MatrixVectorStored.h"
#include "ModelBase.h"
#include "ModelHelperLocal.h"
#include "ModelSelector.h"
#include "OperatorSpec.h"
#include "Operators.h"
#include "ParametersDmrgSolver.h"
#include "ProgramGlobals.h"
#include "Qn.h"
#include "SuperGeometry.h"
#include "TargetingBase.h"
#include "VectorWithOffset.h"
#include "VectorWithOffsets.h"

#ifndef USE_FLOAT
using RealType = double;
#else
using RealType = float;
#endif

struct OperatorOptions {

	enum class IntrospectEnum
	{
		EXPRESSION,
		MODEL_BASIS,
		MODEL_HAMILTONIAN
	};

	OperatorOptions()
	    : site(0)
	    , introspect(IntrospectEnum::EXPRESSION)
	    , enabled(false)
	{ }

	SizeType           site;
	PsimagLite::String opexpr;
	IntrospectEnum     introspect;
	bool               enabled;
};

template <typename ModelBaseType>
void operatorDriver(const ModelBaseType& model, const OperatorOptions& obsOptions)
{
	using ModelHelperType  = typename ModelBaseType::ModelHelperType;
	using OperatorsType    = typename ModelHelperType::OperatorsType;
	using OperatorType     = typename OperatorsType::OperatorType;
	using OperatorSpecType = Dmrg::OperatorSpec<ModelBaseType, OperatorType>;

	OperatorType       opC;
	const OperatorType opEmpty;
	OperatorSpecType   opSpec(model);
	int                site = -1;

	switch (obsOptions.introspect) {
	case OperatorOptions::IntrospectEnum::EXPRESSION:
		PsimagLite::CanonicalExpression<OperatorSpecType> canonicalExpression(opSpec);
		canonicalExpression(opC, obsOptions.opexpr, opEmpty, site);
		opC.write(std::cout);
		break;
	case OperatorOptions::IntrospectEnum::MODEL_BASIS:
		model.printBasis(obsOptions.site);
		break;
	case OperatorOptions::IntrospectEnum::MODEL_HAMILTONIAN:
		model.printTerms();
		break;
	}
}
#endif // DMRGDRIVER_H
