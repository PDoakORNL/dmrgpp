/*
Copyright (c) 2009-2015, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 3.0]
[by G.A., Oak Ridge National Laboratory]

UT Battelle Open Source Software License 11242008

OPEN SOURCE LICENSE

Subject to the conditions of this License, each
contributor to this software hereby grants, free of
charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), a
perpetual, worldwide, non-exclusive, no-charge,
royalty-free, irrevocable copyright license to use, copy,
modify, merge, publish, distribute, and/or sublicense
copies of the Software.

1. Redistributions of Software must retain the above
copyright and license notices, this list of conditions,
and the following disclaimer.  Changes or modifications
to, or derivative works of, the Software should be noted
with comments and the contributor and organization's
name.

2. Neither the names of UT-Battelle, LLC or the
Department of Energy nor the names of the Software
contributors may be used to endorse or promote products
derived from this software without specific prior written
permission of UT-Battelle.

3. The software and the end-user documentation included
with the redistribution, with or without modification,
must include the following acknowledgment:

"This product includes software produced by UT-Battelle,
LLC under Contract No. DE-AC05-00OR22725  with the
Department of Energy."

*********************************************************
DISCLAIMER

THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT OWNER, CONTRIBUTORS, UNITED STATES GOVERNMENT,
OR THE UNITED STATES DEPARTMENT OF ENERGY BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED
STATES DEPARTMENT OF ENERGY, NOR THE COPYRIGHT OWNER, NOR
ANY OF THEIR EMPLOYEES, REPRESENTS THAT THE USE OF ANY
INFORMATION, DATA, APPARATUS, PRODUCT, OR PROCESS
DISCLOSED WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*********************************************************

*/
/** \ingroup DMRG */
/*@{*/

/*! \file TargetQuantumElectrons.h
 *
 *
 */
#ifndef TargetQuantumElectrons_H
#define TargetQuantumElectrons_H
#include "Vector.h"
#include "ProgramGlobals.h"

namespace Dmrg {
//! Hubbard Model Parameters
template<typename RealType>
class TargetQuantumElectrons {

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

public:

	template<typename IoInputType>
	TargetQuantumElectrons(IoInputType& io, bool allowUpDown = true)
	    : twiceJ_(0)
	{
		PsimagLite::String  msg("TargetQuantumElectrons: ");
		bool hasTwiceJ = false;
		try {
			io.read(twiceJ_,"TargetSpinTimesTwo=");
			hasTwiceJ = true;
		} catch (std::exception&) {}

		SizeType ready = 0;
		if (allowUpDown) {
			SizeType electronsUp = 0;
			SizeType electronsDown = 0;
			try {
				io.read(electronsUp,"TargetElectronsUp=");
				io.read(electronsDown,"TargetElectronsDown=");
				totalElectrons_ = electronsUp + electronsDown;
				szPlusConst_ = electronsUp;
				ready++;
			} catch (std::exception&) {}
		}

		try {
			io.read(totalElectrons_,"TargetElectronsTotal=");
			io.read(szPlusConst_,"TargetSzPlusConst=");
			ready++;
		} catch (std::exception&) {}

		switch (ready) {
		case 2:
			msg += "Provide either up/down or total/sz but not both.\n";
			throw PsimagLite::RuntimeError(msg);

		case 0:
			msg += "Provide at least one of up/down or total/sz.\n";
			throw PsimagLite::RuntimeError(msg);
		}

		int tmp = 0;
		try {
			io.readline(tmp,"UseSu2Symmetry=");
		} catch (std::exception&) {}

		isSu2_ = (tmp > 0);

		if (isSu2_ && !hasTwiceJ) {
			msg += "Please provide TargetSpinTimesTwo when running with SU(2).\n";
			throw PsimagLite::RuntimeError(msg);
		}
	}

	template<typename SomeMemResolvType>
	SizeType memResolv(SomeMemResolvType&,
	                   SizeType,
	                   PsimagLite::String = "") const
	{
		return 0;
	}

	void setTargetNumbers(VectorSizeType& t,
	                      SizeType sites,
	                      SizeType totalSites,
	                      SizeType direction) const
	{
		t.resize((isSu2_) ? 3 : 2,0);

		if (direction == ProgramGlobals::INFINITE) {
			t[0] = static_cast<SizeType>(round(szPlusConst_*sites/totalSites));
			t[1] = static_cast<SizeType>(round(totalElectrons_*sites/totalSites));
		} else {
			t[0] = szPlusConst_;
			t[1] = totalElectrons_;
		}

		if (t.size() < 3) return;

		RealType jReal = twiceJ_*sites/static_cast<RealType>(totalSites);
		SizeType tmp = (direction == ProgramGlobals::INFINITE) ?
		            static_cast<SizeType>(round(jReal)) : twiceJ_;

		if (totalElectrons_%2==0) {
			if (tmp%2 != 0) tmp++;
		} else {
			if (tmp%2 == 0) tmp++;
		}

		t[2] = tmp;
	}

	void print(std::ostream& os) const
	{
		os<<"TargetElectronsTotal="<<totalElectrons_<<"\n";
		os<<"TargetSzPlusConst="<<szPlusConst_<<"\n";
		if (isSu2_)
			os<<"TargetSpinTimesTwo="<<twiceJ_<<"\n";
	}

private:

	bool isSu2_;
	SizeType totalElectrons_;
	SizeType szPlusConst_;
	SizeType twiceJ_;
};

//! Function that prints model parameters to stream os
template<typename RealTypeType>
std::ostream& operator<<(std::ostream &os,
                         const TargetQuantumElectrons<RealTypeType>& p)
{
	p.print(os);
	return os;
}
} // namespace Dmrg

/*@}*/
#endif

