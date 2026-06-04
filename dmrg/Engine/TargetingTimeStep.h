/*
Copyright (c) 2009-2014, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 5.]
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

#ifndef TARGETING_TIMESTEP_H
#define TARGETING_TIMESTEP_H

#include "BlockDiagonalMatrix.h"
#include "ParametersForSolver.h"
#include "PredicateAwesome.h"
#include "ProgramGlobals.h"
#include "ProgressIndicator.h"
#include "TargetParamsTimeStep.h"
#include "TargetingBase.h"
#include "TimeVectorsKrylov.h"
#include "TimeVectorsRungeKutta.h"
#include "TimeVectorsSuzukiTrotter.h"
#include <iostream>

namespace Dmrg {

template <typename LanczosSolverType_, typename VectorWithOffsetType_>
class TargetingTimeStep : public TargetingBase<LanczosSolverType_, VectorWithOffsetType_> {

	enum
	{
		BORDER_NEITHER,
		BORDER_LEFT,
		BORDER_RIGHT
	};

public:

	using LanczosSolverType      = LanczosSolverType_;
	using BaseType               = TargetingBase<LanczosSolverType, VectorWithOffsetType_>;
	using TargetingCommonType    = typename BaseType::TargetingCommonType;
	using CheckpointType         = typename BaseType::CheckpointType;
	using PairType               = std::pair<SizeType, SizeType>;
	using OptionsType            = typename BaseType::OptionsType;
	using MatrixVectorType       = typename BaseType::MatrixVectorType;
	using ModelType              = typename MatrixVectorType::ModelType;
	using RealType               = typename ModelType::RealType;
	using OperatorsType          = typename ModelType::OperatorsType;
	using ModelHelperType        = typename ModelType::ModelHelperType;
	using LeftRightSuperType     = typename ModelHelperType::LeftRightSuperType;
	using BasisWithOperatorsType = typename LeftRightSuperType::BasisWithOperatorsType;
	using VectorSizeType         = PsimagLite::Vector<SizeType>::Type;
	using WaveFunctionTransfType = typename BaseType::WaveFunctionTransfType;
	using VectorWithOffsetType   = typename WaveFunctionTransfType::VectorWithOffsetType;
	using ComplexOrRealType      = typename VectorWithOffsetType::value_type;
	using TargetVectorType       = typename VectorWithOffsetType::VectorType;
	using VectorRealType         = typename PsimagLite::Vector<RealType>::Type;
	using OperatorType           = typename BasisWithOperatorsType::OperatorType;
	using BasisType              = typename BasisWithOperatorsType::BasisType;
	using TargetParamsType       = TargetParamsTimeStep<ModelType>;
	using BlockType              = typename BasisType::BlockType;
	using TimeSerializerType     = typename TargetingCommonType::TimeSerializerType;
	using SparseMatrixType       = typename OperatorType::StorageType;
	using InputValidatorType     = typename ModelType::InputValidatorType;
	using QnType                 = typename BasisType::QnType;
	using StageEnumType          = typename TargetingCommonType::StageEnumType;

	TargetingTimeStep(const LeftRightSuperType&     lrs,
	                  const CheckpointType&         checkPoint,
	                  const WaveFunctionTransfType& wft,
	                  const QnType&,
	                  InputValidatorType& ioIn,
	                  PsimagLite::String  targeting)
	    : BaseType(lrs, checkPoint, wft, 0)
	    , tstStruct_(ioIn, targeting, checkPoint.model())
	    , wft_(wft)
	    , progress_(targeting)
	    , weight_(tstStruct_.times().size())
	    , tvEnergy_(tstStruct_.times().size(), 0.0)
	    , gsWeight_(tstStruct_.gsWeight())
	    , evolveGs_(tstStruct_.evolveGroundState())
	    , gsEvolvedIdx_(2 * tstStruct_.times().size() - 1)
	    , gsEvolvedCurrentBasis_(false)
	    , gsWftAppliedThisAdvance_(false)
	{
		if (!wft.isEnabled())
			err("TST needs an enabled wft\n");
		if (tstStruct_.sites() == 0)
			err("TST needs at least one TSPSite\n");

		RealType tau = tstStruct_.tau();
		RealType sum = 0;
		SizeType n   = tstStruct_.times().size();

		RealType factor = (n + 4.0) / (n + 2.0);
		factor *= (1.0 - gsWeight_);
		for (SizeType i = 0; i < n; i++) {
			tstStruct_.times()[i] = i * tau / (n - 1);
			weight_[i]            = factor / (n + 4);
			sum += weight_[i];
		}
		sum -= weight_[0];
		sum -= weight_[n - 1];
		weight_[0] = weight_[n - 1] = 2 * factor / (n + 4);
		sum += weight_[n - 1];
		sum += weight_[0];

		gsWeight_ = 1.0 - sum;
		sum += gsWeight_;
		assert(fabs(sum - 1.0) < 1e-5);

		this->common().aoeNonConst().initTimeVectors(tstStruct_, ioIn);

		if (evolveGs_) {
			// Register the extra slot index with TargetingCommon so that
			// "gsT" in dressed measurement labels resolves to this P-vector.
			// The slot itself is allocated by postCtor() → targetVectorsResize(targets()),
			// which picks up the +1 from the overridden targets() above.
			this->common().setEvolvedGsIndex(gsEvolvedIdx_);
		}
	}

	SizeType sites() const { return tstStruct_.sites(); }

	// With TSPEvolveGroundState=1 we append 2 extra Krylov slots for the GS:
	//   slot N   = "anchor" (GS at Krylov t=0, starting point for this advance)
	//   slot N+1 = evolved GS at the same Krylov time as P1 (= gsEvolvedIdx_)
	// N extra slots for GS Krylov: anchor (times[0]=0) + N-1 evolved states (times[1..N-1]).
	// The last slot (index 2N-1) accumulates a full tau step per advance, same as P-vectors.
	SizeType targets() const { return tstStruct_.times().size() * (evolveGs_ ? 2 : 1); }

	RealType weight(SizeType i) const
	{
		assert(!this->common().aoe().allStages(StageEnumType::DISABLED));
		// GS Krylov slots (N..2N-1) carry zero truncation weight.
		if (evolveGs_ && i >= tstStruct_.times().size())
			return 0.0;
		return weight_[i];
	}

	RealType gsWeight() const
	{
		if (this->common().aoe().allStages(StageEnumType::DISABLED))
			return 1.0;
		return gsWeight_;
	}

	bool includeGroundStage() const
	{
		if (!this->common().aoe().noStageIs(StageEnumType::DISABLED))
			return true;
		bool b = (fabs(gsWeight_) > 1e-6);
		return b;
	}

	void evolve(const VectorRealType&         energies,
	            ProgramGlobals::DirectionEnum direction,
	            const BlockType&              block1,
	            const BlockType&,
	            SizeType loopNumber)
	{
		assert(block1.size() > 0);
		SizeType site = block1[0];
		assert(energies.size() > 0);
		RealType Eg = energies[0];
		// Reset per-advance flags.
		gsEvolvedCurrentBasis_ = false;
		gsWftAppliedThisAdvance_ = false;
		evolveInternal(Eg, direction, block1, loopNumber);

		SizeType numberOfSites = this->lrs().super().block().size();

		if (site > 1 && site < numberOfSites - 2)
			return;

		if (direction == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) {
			if (site == 1)
				return;
		} else {
			if (site == numberOfSites - 2)
				return;
		}

		SizeType  x = (site == 1) ? 0 : numberOfSites - 1;
		BlockType block(1, x);
		evolveInternal(Eg, direction, block, loopNumber);
	}

	bool end() const
	{
		return (tstStruct_.maxTime() != 0
		        && this->common().aoe().timeVectors().time() >= tstStruct_.maxTime());
	}

	void read(typename TargetingCommonType::IoInputType& io, PsimagLite::String prefix)
	{
		this->common().readGSandNGSTs(io, prefix, "TimeStep");
		// Save |GS_i⟩ immediately after checkpoint load, before any DMRG sweep overwrites
		// psiConst.  This is the pre-quench GS with a sign fixed by the hdf5 file — the
		// same in both the particle and hole tDMRG runs, ensuring sign consistency in G^<.
		if (evolveGs_) {
			const auto& psi = this->common().aoe().psiConst();
			if (psi.size() > 0 && psi[0].size() > 0 && psi[0][0] != nullptr
			    && psi[0][0]->size() > 0) {
				gsInitial_ = *psi[0][0];
				// Note: we do NOT pre-seed targetVectors[gsEvolvedIdx_] here.
				// Between read() and the first advance, the DMRG changes the basis
				// at every site.  Seeding here would leave gsT in the checkpoint
				// basis, causing size mismatches in cocoon measurements.  The GS
				// Krylov block in evolveInternal() seeds it on the first advance
				// using the WFT to transform gsInitial_ to the current basis.
			}
		}
	}

	void write(const VectorSizeType&        block,
	           PsimagLite::IoSelector::Out& io,
	           PsimagLite::String           prefix) const
	{
		PsimagLite::OstringStream                     msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg << "Saving state...";
		progress_.printline(msgg, std::cout);

		this->common().write(io, block, prefix);
		this->common().writeNGSTs(io, prefix, block, "TimeStep");
	}

private:

	void evolveInternal(RealType                      Eg,
	                    ProgramGlobals::DirectionEnum direction,
	                    const BlockType&              block1,
	                    SizeType                      loopNumber)
	{
		if (direction == ProgramGlobals::DirectionEnum::INFINITE)
			return;
		VectorWithOffsetType phiNew;
		assert(block1.size() > 0);
		SizeType site = block1[0];

		// WFT the GS Krylov slot at the START of each DMRG advance, before any cocoon
		// or Krylov call.  This keeps targetVectors[gsEvolvedIdx_] in the current basis
		// by chaining through every WFT step from the checkpoint to the current advance.
		// gsWftAppliedThisAdvance_ prevents double-WFT on the corner sub-call.
		if (evolveGs_ && !gsWftAppliedThisAdvance_) {
			// Seed on very first call: copy gsInitial_ (checkpoint basis).
			// The WFT below then maps it to the current basis.
			if (this->common().aoe().targetVectors(gsEvolvedIdx_).size() == 0
			    && gsInitial_.size() > 0)
				this->tvNonConst(gsEvolvedIdx_) = gsInitial_;
			if (this->common().aoe().targetVectors(gsEvolvedIdx_).size() > 0)
				this->common().aoeNonConst().wftSome(site, gsEvolvedIdx_, gsEvolvedIdx_ + 1);
			gsWftAppliedThisAdvance_ = true;
		}

		SizeType numberOfSites    = this->lrs().super().block().size();
		bool     doBorderIfBorder = (site < 1 || site >= numberOfSites - 1);

		if (doBorderIfBorder) {
			if (loopNumber >= this->model().params().finiteLoop.size() - 1) {
				if (direction == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) {
					if (site >= numberOfSites - 1) {
						this->common().cocoon(block1, direction, false);
						return;
					}
				} else {
					if (site < 1) {
						this->common().cocoon(block1, direction, false);
						return;
					}
				}
			}
			this->common().cocoon(block1, direction, false);
		}

		this->common().aoeNonConst().getPhi(
		    &phiNew, Eg, direction, site, loopNumber, tstStruct_, nullptr);

		PairType startEnd(0, tstStruct_.times().size());
		bool     allOperatorsApplied
		    = (this->common().aoe().noStageIs(StageEnumType::DISABLED)
		       && this->common().aoe().noStageIs(StageEnumType::OPERATOR));

		VectorSizeType indices(startEnd.second - startEnd.first);
		for (SizeType i = 0; i < indices.size(); ++i)
			indices[i] = i + startEnd.first;

		static const bool isLastCall = true;
		this->common().aoeNonConst().calcTimeVectors(
		    indices,
		    Eg,
		    phiNew,
		    direction,
		    allOperatorsApplied,
		    false, // don't wft or advance indices[0]
		    block1,
		    isLastCall);

		// Krylov-evolve |GS_i⟩ → |Ψ_0(t)⟩ = e^{−iH_f t}|GS_i⟩ for "gsT" measurements.
		//
		// gsInitial_ is saved in read() from the checkpoint, before DMRG sweeps change
		// psiConst.  Both the particle and hole tDMRG runs load the same checkpoint, so
		// gsInitial_ carries the same sign in both runs — ensuring G^R = G^>−G^< is
		// purely imaginary as required by causality.
		//
		// Accumulation: WFT-transform targetVectors[gsEvolvedIdx_] from the previous
		// advance's basis to the current basis, then Krylov-evolve by one τ step.
		// A sign-flip check after WFT detects and corrects spurious −1 gauge factors.
		//
		// gsEvolvedCurrentBasis_ prevents the corner sub-call in evolve() from
		// re-WFT-ing a freshly computed vector (same-advance double-WFT → norm=0).
		if (evolveGs_ && allOperatorsApplied && !gsEvolvedCurrentBasis_) {
			const SizeType gsAnchorIdx = tstStruct_.times().size();
			const VectorWithOffsetType& gsEvolved
			    = this->common().aoe().targetVectors(gsEvolvedIdx_);
			// gsEvolved is already in the current DMRG basis: WFT was applied at the
			// start of evolveInternal() (gsWftAppliedThisAdvance_ flag).
			const VectorWithOffsetType psi0
			    = (gsEvolved.size() > 0) ? gsEvolved
			                             : *this->common().aoe().psiConst()[0][0];
			// Build indices {N, N+1, ..., 2N-1} so the Krylov loop uses
			// times[0..N-1], giving targetVectors[2N-1] a full tau evolution
			// per advance (matching the particle-sector accumulation rate).
			const SizeType N = tstStruct_.times().size();
			VectorSizeType gsIdx(N);
			for (SizeType j = 0; j < N; ++j)
				gsIdx[j] = gsAnchorIdx + j;
			this->common().aoeNonConst().calcTimeVectors(
			    gsIdx,
			    Eg,
			    psi0,
			    direction,
			    allOperatorsApplied,
			    false,
			    block1,
			    isLastCall);
			// Gauge sign correction: the WFT can flip the sign of the N-sector
			// block (a well-known DMRG gauge freedom).  Both the particle and hole
			// tDMRG runs experience the same flip (their N-sector DM is identical),
			// so a flip makes G^R → −G^R_correct for all subsequent time steps —
			// Max|dIm(G^R)| ≈ 2.  This correction is the same operation that
			// DMRG++ already applies to psiConst via <gs|penultimate> tracking.
			//
			// Reference: psiConst()[0][0] is always positively phased (DMRG++
			// Lanczos convention) and lives in the same N-sector.  The overlap
			// ⟨GS_f|Ψ_0(t)⟩ stays ≥ 0 for moderate t and quench (neq-DMFT regime);
			// for very strong quenches or very long t this could fail — document
			// with a comment if that regime is reached.
			{
				const VectorWithOffsetType& gsEvNew
				    = this->common().aoe().targetVectors(gsEvolvedIdx_);
				const auto& refPsi = this->common().aoe().psiConst();
				if (gsEvNew.size() > 0 && refPsi.size() > 0 && refPsi[0].size() > 0
				    && refPsi[0][0] != nullptr && refPsi[0][0]->size() > 0) {
					const ComplexOrRealType ov = *refPsi[0][0] * gsEvNew;
					if (PsimagLite::real(ov) < 0)
						this->tvNonConst(gsEvolvedIdx_) *= ComplexOrRealType(-1);
				}
			}
			// Mark the evolved GS as being in the current basis so the corner
			// sub-call from evolve() does not re-WFT or re-evolve it.
			gsEvolvedCurrentBasis_ = true;
		}


		this->common().cocoon(block1, direction, false);

		PsimagLite::String predicate = this->model().params().printHamiltonianAverage;
		const SizeType     center    = this->model().superGeometry().numberOfSites() / 2;
		PsimagLite::replaceAll(predicate, "c", ttos(center));
		PsimagLite::PredicateAwesome<> pAwesome(predicate);
		assert(block1.size() > 0);
		if (pAwesome.isTrue("s", block1[0]))
			printEnergies(); // in-situ

		const OptionsType& options = this->model().params().options;
		bool               normalizeTimeVectors
		    = (options.isSet("normalizeTimeVectors") || options.isSet("TargetingAncilla"));

		if (options.isSet("neverNormalizeVectors"))
			normalizeTimeVectors = false;

		if (normalizeTimeVectors)
			this->common().normalizeTimeVectors();

		this->common().printNormsAndWeights(gsWeight_, weight_);
	}

	void printEnergies() const
	{
		// Skip GS Krylov slots (indices >= N): they are for measurement only,
		// and tvEnergy_ has only N entries.
		const SizeType nMain = tstStruct_.times().size();
		for (SizeType i = 0; i < nMain; i++)
			printEnergies(this->tv(i), i);
	}

	void printEnergies(const VectorWithOffsetType& phi, SizeType whatTarget) const
	{
		for (SizeType ii = 0; ii < phi.sectors(); ii++) {
			SizeType i = phi.sector(ii);
			printEnergies(phi, whatTarget, i);
		}
	}

	void printEnergies(const VectorWithOffsetType& phi, SizeType whatTarget, SizeType i0) const
	{
		const SizeType p = this->lrs().super().findPartitionNumber(phi.offset(i0));
		typename ModelHelperType::Aux                 aux(p, BaseType::lrs());
		typename ModelType::HamiltonianConnectionType hc(
		    BaseType::lrs(),
		    ModelType::modelLinks(),
		    this->common().aoe().timeVectors().time(),
		    BaseType::model().superOpHelper(),
		    BaseType::model().ioIn());
		typename LanczosSolverType::MatrixType lanczosHelper(BaseType::model(), hc, aux);

		const SizeType   total = phi.effectiveSize(i0);
		TargetVectorType phi2(total);
		phi.extract(phi2, i0);
		TargetVectorType x(total);
		lanczosHelper.matrixVectorProduct(x, phi2);
		PsimagLite::OstringStream                     msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg << "Hamiltonian average at time=" << this->common().aoe().timeVectors().time();
		msg << " for target=" << whatTarget;
		ComplexOrRealType numerator = phi2 * x;
		ComplexOrRealType den       = phi2 * phi2;
		ComplexOrRealType division  = (PsimagLite::norm(den) < 1e-10) ? 0 : numerator / den;
		msg << " sector=" << i0 << " <phi(t)|H|phi(t)>=" << numerator;
		msg << " <phi(t)|phi(t)>=" << den << " " << division;
		progress_.printline(msgg, std::cout);
		tvEnergy_[whatTarget] = PsimagLite::real(division);
	}

	TargetParamsType              tstStruct_;
	const WaveFunctionTransfType& wft_;
	PsimagLite::ProgressIndicator progress_;
	VectorRealType                weight_;
	mutable VectorRealType        tvEnergy_;
	RealType                      gsWeight_;
	bool                          evolveGs_;              // TSPEvolveGroundState=1
	SizeType                      gsEvolvedIdx_;          // targetVectors slot for |Ψ_0(t)⟩
	bool                          gsEvolvedCurrentBasis_; // true once WFT+evolve done this advance
	VectorWithOffsetType          gsInitial_;             // |GS_i⟩ from checkpoint (pre-quench)
	bool                          gsWftAppliedThisAdvance_; // prevents double-WFT at corner
}; // class TargetingTimeStep
} // namespace Dmrg

#endif
