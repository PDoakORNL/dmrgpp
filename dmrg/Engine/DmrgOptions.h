#ifndef DMRG_OPTIONS_H
#define DMRG_OPTIONS_H
#include "ProgramGlobals.h"
#include <PsimagLite/PsimagLite.h>
#include <algorithm>
#include <cctype>
#include <numeric>

namespace Dmrg {

template <typename InputValidatorType> class DmrgOptions {

public:

	using CharType            = typename PsimagLite::String::value_type;
	using StringConstIterator = typename PsimagLite::String::const_iterator;
	using VectorStringType    = PsimagLite::Vector<PsimagLite::String>::Type;

	DmrgOptions(PsimagLite::String label, InputValidatorType& io)
	{
		PsimagLite::String tmp;
		io.readline(tmp, label);
		PsimagLite::split(vdata_, tmp, ",");
		lowerAll();
	}

	void operator+=(PsimagLite::String moreData)
	{
		VectorStringType vmore;
		PsimagLite::split(vmore, moreData, ",");
		vdata_.insert(vdata_.end(), vmore.begin(), vmore.end());
		lowerAll();
	}

	void write(PsimagLite::String label, PsimagLite::IoSerializer& ioSerializer) const
	{
		const PsimagLite::String tmp
		    = std::accumulate(vdata_.begin(), vdata_.end(), PsimagLite::String(","));
		ioSerializer.write(label, tmp);
	}

	bool isSet(PsimagLite::String what) const
	{
		what                                = ProgramGlobals::toLower(what);
		VectorStringType::const_iterator it = std::find(vdata_.begin(), vdata_.end(), what);
		return (it != vdata_.end());
	}

private:

	void lowerAll()
	{
		std::transform(vdata_.begin(),
		               vdata_.end(),
		               vdata_.begin(),
		               [](PsimagLite::String s) { return ProgramGlobals::toLower(s); });
	}

	VectorStringType vdata_;
};
}
#endif // DMRG_OPTIONS_H
