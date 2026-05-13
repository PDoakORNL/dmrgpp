#ifndef OPTIONSFORINTROSPECT_HH
#define OPTIONSFORINTROSPECT_HH
#include "PsimagLite.h"

namespace Dmrg {

struct OptionsForIntrospect {

	enum class IntrospectEnum
	{
		EXPRESSION,
		MODEL_BASIS,
		MODEL_HAMILTONIAN
	};

	OptionsForIntrospect()
	    : site(0)
	    , introspect(IntrospectEnum::EXPRESSION)
	{ }

	SizeType           site;
	PsimagLite::String opexpr;
	IntrospectEnum     introspect;
};
}
#endif // OPTIONSFORINTROSPECT_HH
