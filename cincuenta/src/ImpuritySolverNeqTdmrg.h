#ifndef IMPURITYSOLVER_NEQ_TDMRG_H
#define IMPURITYSOLVER_NEQ_TDMRG_H

#include "CmdLineOptions.hh"
#include "DmrgRunner.h"
#include "ImpuritySolverNeqBase.h"
#include "ImpuritySolverNeqExactDiag.h"
#include "KadanoffBaym.h"
#include "NeqBathDecomposition.h"
#include "ParamsNeqDmftSolver.h"
#include "PsimagLite.h"
#include "Vector.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <string>

namespace Dmft {

// Non-equilibrium impurity solver using DMRG++ time-dependent DMRG (tDMRG).
//
// Five-run approach using TargetingExpression with <P.last|P> gauge correction
// (from pr/dmrgpp-pvector-last-syntax).  This replaces the old TSPEvolveGroundState
// mechanism which did not fit the DMRG++ design.
//
// Run 1:  GS DMRG with H(U_i)  →  saves |GS_i⟩ to root_+gs.
// Run 2:  TargetingExpression restart, FiniteLoops flag=2 (keeps |GS_i⟩ as |gs⟩),
//         computes P0 = c'[0]*|gs⟩  →  saves to root_+particle.
// Run 3:  TargetingExpression tDMRG, restarts from root_+particle, flag=2:
//           P1 = TimeEvolve*|P0⟩ = e^{-iH_f t} c†|GS_i⟩
//           P2 = TimeEvolve*|gs⟩ = e^{-iH_f t} |GS_i⟩
//         In-situ: <P2|c|P1>   →  G^>(t,0)
//                  <P2.last|P2> →  consecutive GS overlap (gauge-phase correction)
// Run 4:  Same as Run 2 but with c[0]  →  saves to root_+hole.
// Run 5:  Same as Run 3 but restarts from root_+hole:
//         In-situ: <P1|c|P2>   →  G^<(t,0)
//                  <P2.last|P2> →  gauge-phase correction
//
// G^R(t,0) = G^>(t,0) - G^<(t,0)
//
// Validation: run with BandwidthFinal=BandwidthInitial and HubbardUFinal=HubbardU
// (no quench).  G^R(t,0) from tDMRG must agree with the ED solver to truncation error.
//
// Required input parameters:
//   TargetElectronsUp=, TargetElectronsDown=
//   RootOutputname=
//   InfiniteLoopKeptStates=
//   FiniteLoopsGs=        — loops for GS run (flag 0 OK)
//   FiniteLoopsTdmrg=     — loops for init/tDMRG runs (flag 2 enforced automatically;
//                           must contain enough sweeps for params_.nT time advances
//                           at TSPAdvanceEach sites each)
//   TSPTimeSteps=         — Krylov substep count (default 5, Gonzalo's preference)
//   TSPAdvanceEach=       — sites per time advance (use N-2 for N-site star; default 1)
template <typename ComplexOrRealType>
class ImpuritySolverNeqTdmrg : public ImpuritySolverNeqBase<ComplexOrRealType> {

public:

	using BaseType          = ImpuritySolverNeqBase<ComplexOrRealType>;
	using RealType          = typename BaseType::RealType;
	using ComplexType       = typename BaseType::ComplexType;
	using VectorRealType    = typename BaseType::VectorRealType;
	using VectorComplexType = typename PsimagLite::Vector<ComplexType>::Type;
	using KBType            = typename BaseType::KBType;
	using InputNgType       = typename BaseType::InputNgType;
	using ParamsNeqType     = ParamsNeqDmftSolver<ComplexOrRealType>;
	using DmrgRunnerType    = Dmrg::DmrgRunner<RealType>;
	using ApplicationType   = PsimagLite::PsiApp;
	using ExactDiagType     = ImpuritySolverNeqExactDiag<ComplexOrRealType>;
	using DecompType        = NeqBathDecomposition<ComplexOrRealType>;

	ImpuritySolverNeqTdmrg(const ParamsNeqType&            params,
	                       const ApplicationType&          app,
	                       typename InputNgType::Readable& io)
	    : params_(params)
	    , app_(app)
	    , io_(io)
	    , neqBathRank_(params.neqBathRank)
	    , exactDiag_(params, io)
	    , gimp_(params.nT,
	            params.eqParams.nMatsubaras,
	            params.dt,
	            params.eqParams.ficticiousBeta
	                / static_cast<RealType>(params.eqParams.nMatsubaras))
	{
		io.readline(nup_, "TargetElectronsUp=");
		io.readline(ndown_, "TargetElectronsDown=");
		io.readline(root_, "RootOutputname=");
		io.readline(infiniteLoops_, "InfiniteLoopKeptStates=");
		io.readline(finiteLoopsGs_, "FiniteLoopsGs=");

		try {
			io.readline(finiteLoopsTdmrg_, "FiniteLoopsTdmrg=");
		} catch (std::exception&) {
			finiteLoopsTdmrg_ = finiteLoopsGs_;
		}

		try {
			io.readline(tspTimeSteps_, "TSPTimeSteps=");
		} catch (std::exception&) {
			tspTimeSteps_ = 5;
		}

		try {
			io.readline(tspAdvanceEach_, "TSPAdvanceEach=");
		} catch (std::exception&) {
			tspAdvanceEach_ = 1;
		}
	}

	// Runs all five DMRG passes and fills the KB grid. NeqBathRank=0 only --
	// see solveSelfConsistent for the NeqBathRank>0 (evolving-bath) path,
	// which this dispatches to unchanged and untouched otherwise (byte-for-
	// byte the same code as before Phase 2).
	void solve(const VectorRealType& bathParams) override
	{
		if (neqBathRank_ > 0) {
			solveSelfConsistent(bathParams);
			return;
		}

		const SizeType nBath  = bathParams.size() / 2;
		const SizeType nsites = nBath + 1;

		VectorRealType hoppings(nBath), bathEps(nBath);
		for (SizeType i = 0; i < nBath; ++i) {
			hoppings[i] = bathParams[i];
			bathEps[i]  = bathParams[nBath + i];
		}

		VectorRealType potGS(nsites), potTdmrg(nsites);
		potGS[0]    = -RealType(0.5) * params_.uInitial;
		potTdmrg[0] = -RealType(0.5) * params_.uFinal;
		for (SizeType i = 0; i < nBath; ++i) {
			potGS[i + 1]    = bathEps[i];
			potTdmrg[i + 1] = bathEps[i];
		}

		// Run 1: GS DMRG for H(U_i)
		std::cout << "ImpuritySolverNeqTdmrg: GS DMRG U_i=" << params_.uInitial << "\n";
		{
			Dmrg::CmdLineOptions opts;
			opts.logfile = root_ + "gs.log";
			DmrgRunnerType runner(
			    app_,
			    buildGsInput(params_.uInitial, hoppings, potGS, nup_, ndown_, nsites),
			    opts);
			runner.doOneRun();
		}

		// Run 2: apply c'[0] to |GS_i⟩, save particle initial state
		std::cout << "ImpuritySolverNeqTdmrg: particle-init (c'|GS_i>)\n";
		{
			Dmrg::CmdLineOptions opts;
			opts.logfile = root_ + "particle_init.log";
			DmrgRunnerType runner(
			    app_,
			    buildParticleInitInput(
			        params_.uFinal, hoppings, potTdmrg, nup_, ndown_, nsites),
			    opts);
			runner.doOneRun();
		}

		// Run 3: particle tDMRG — G^>(t,0) and gauge
		const std::string tdmrgLog = root_ + "tdmrg.log";
		std::cout << "ImpuritySolverNeqTdmrg: particle tDMRG U_f=" << params_.uFinal
		          << " t_max=" << params_.tMax << "\n";
		{
			Dmrg::CmdLineOptions opts;
			opts.logfile              = tdmrgLog;
			opts.in_situ_measurements = "<P2|c|P1>,<P2.last|P2>";
			DmrgRunnerType runner(
			    app_,
			    buildTdmrgInput(
			        params_.uFinal, hoppings, potTdmrg, nup_, ndown_, nsites),
			    opts);
			runner.doOneRun();
		}

		// Run 4: apply c[0] to |GS_i⟩, save hole initial state
		std::cout << "ImpuritySolverNeqTdmrg: hole-init (c|GS_i>)\n";
		{
			Dmrg::CmdLineOptions opts;
			opts.logfile = root_ + "hole_init.log";
			DmrgRunnerType runner(
			    app_,
			    buildHoleInitInput(
			        params_.uFinal, hoppings, potTdmrg, nup_, ndown_, nsites),
			    opts);
			runner.doOneRun();
		}

		// Run 5: hole tDMRG — G^<(t,0) and gauge
		const std::string holeTdmrgLog = root_ + "tdmrg_hole.log";
		std::cout << "ImpuritySolverNeqTdmrg: hole tDMRG G^<(t,0)\n";
		{
			Dmrg::CmdLineOptions opts;
			opts.logfile              = holeTdmrgLog;
			opts.in_situ_measurements = "<P1|c|P2>,<P2.last|P2>";
			DmrgRunnerType runner(
			    app_,
			    buildHoleTdmrgInput(
			        params_.uFinal, hoppings, potTdmrg, nup_, ndown_, nsites),
			    opts);
			runner.doOneRun();
		}

		std::map<int, ComplexType> ggt0_at_step;
		parseTdmrgLog(tdmrgLog, ggt0_at_step);

		std::map<int, ComplexType> glt0_at_step;
		parseHoleTdmrgLog(holeTdmrgLog, glt0_at_step);

		fillKBGrid(ggt0_at_step, glt0_at_step);
	}

	void computeGimp(KBType& gimp, int n) const override
	{
		if (neqBathRank_ == 0) {
			gimp.retarded(n, 0) = gimp_.retarded(n, 0);
			gimp.lesser(n, n)   = gimp_.lesser(n, n);
			gimp.lesser(n, 0)   = gimp_.lesser(n, 0);
			return;
		}

		// Matsubara/left-mixing (equilibrium) components: delegated to
		// ExactDiag, exactly as ImpuritySolverNeqGBEK::computeGimp does --
		// tDMRG has no equilibrium solver of its own, and ExactDiag's
		// Lehmann-representation result is already what GBEK itself trusts
		// for this piece.
		exactDiag_.computeGimp(gimp, n);

		// Retarded/lesser (n,j), j=0..n: truncated-batch recompute -- rebuild
		// the fan-out from t=0 up through step n fresh, reading Vplus(k,p)
		// from decomp_'s CURRENT state. Correct and simple (no propagated-
		// through/rewind bookkeeping needed: a corrector's prepareTimeStep
		// call just refines decomp_ before the next computeGimp call, and
		// this recompute always uses whatever is current); cost is O(n^2)
		// per call, O(nT^3) total, acceptable at the tiny NtNeq (2-4) this
		// project's compute-discipline convention already mandates for
		// gate-style tests. See fancy-painting-moon.md, Phase 2, "truncated
		// batch recompute" (advisor consult).
		fillSelfConsistentRow(gimp, n);
	}

	const KBType& gimp() const override { return neqBathRank_ > 0 ? exactDiag_.gimp() : gimp_; }

	// Advance the Cholesky bath decomposition to step n. No propagated-
	// through invalidation needed (unlike GBEK): computeGimp's truncated-
	// batch recompute always rebuilds from scratch using decomp_'s current
	// state, so there is nothing to invalidate.
	void prepareTimeStep(int n, const KBType& delta) override
	{
		if (decomp_)
			decomp_->update(n, delta);
	}

	// ---- Phase 1 (evolving-bath project) diagnostic: chained column 0 -----
	//
	// Verification-only entry point, NOT wired into solve()/computeGimp()'s
	// existing dispatch. Recomputes G^>(t_n,0)/G^<(t_n,0) by chaining nT
	// single-step restarts (one DmrgRunner call per step, each restarting
	// from the previous step's checkpoint) instead of today's single
	// monolithic Run3/Run5 call spanning the whole trajectory. This is the
	// narrowest possible slice of the full multi-column two-time grid: it
	// only re-derives the SAME column (j=0) today's solve() already
	// computes, so its result can be diffed directly against
	// ImpuritySolverNeqTdmrg::solve()'s own gimp_ as a correctness check on
	// the chaining/restart machinery itself (RestartMappingTvs,
	// RestartSourceTvForPsi) before any new-column-birthing logic is added.
	// See /Users/epd/.claude/plans/fancy-painting-moon.md, Phase 1.
	struct ChainedResult {
		std::map<int, ComplexType> ggt0; // G^>(t_n, 0), gauge-corrected
		std::map<int, ComplexType> glt0; // G^<(t_n, 0), gauge-corrected
	};

	ChainedResult solveChainedColumn0(const VectorRealType& bathParams) const
	{
		const SizeType nBath  = bathParams.size() / 2;
		const SizeType nsites = nBath + 1;

		VectorRealType hoppings(nBath), bathEps(nBath);
		for (SizeType i = 0; i < nBath; ++i) {
			hoppings[i] = bathParams[i];
			bathEps[i]  = bathParams[nBath + i];
		}

		VectorRealType potGS(nsites), potTdmrg(nsites);
		potGS[0]    = -RealType(0.5) * params_.uInitial;
		potTdmrg[0] = -RealType(0.5) * params_.uFinal;
		for (SizeType i = 0; i < nBath; ++i) {
			potGS[i + 1]    = bathEps[i];
			potTdmrg[i + 1] = bathEps[i];
		}

		const std::string chainRoot = root_ + "chain_";

		{
			Dmrg::CmdLineOptions opts;
			opts.logfile = chainRoot + "gs.log";
			DmrgRunnerType runner(app_,
			                      buildGsInputAt(chainRoot + "gs",
			                                     params_.uInitial,
			                                     hoppings,
			                                     potGS,
			                                     nup_,
			                                     ndown_,
			                                     nsites),
			                      opts);
			runner.doOneRun();
		}

		std::string particleRoot = chainRoot + "particle0";
		{
			Dmrg::CmdLineOptions opts;
			opts.logfile = chainRoot + "particle_init.log";
			DmrgRunnerType runner(app_,
			                      buildInitInputAt(particleRoot,
			                                       chainRoot + "gs",
			                                       params_.uFinal,
			                                       hoppings,
			                                       potTdmrg,
			                                       nup_,
			                                       ndown_,
			                                       nsites,
			                                       "'"),
			                      opts);
			runner.doOneRun();
		}

		std::string holeRoot = chainRoot + "hole0";
		{
			Dmrg::CmdLineOptions opts;
			opts.logfile = chainRoot + "hole_init.log";
			DmrgRunnerType runner(app_,
			                      buildInitInputAt(holeRoot,
			                                       chainRoot + "gs",
			                                       params_.uFinal,
			                                       hoppings,
			                                       potTdmrg,
			                                       nup_,
			                                       ndown_,
			                                       nsites,
			                                       ""),
			                      opts);
			runner.doOneRun();
		}

		ChainedResult              result;
		std::map<int, ComplexType> gaugeP, gaugeH;
		int                        particleMapTv = 0, particleSrcTv = -1;
		int                        holeMapTv = 0, holeSrcTv = -1;

		for (int n = 1; n <= static_cast<int>(params_.nT); ++n) {
			{
				const std::string outRoot = chainRoot + "particle" + ttos(n);
				// A segment restarting straight from the *_init checkpoint
				// (accumulated time 0) reports its OWN measurement at
				// (local) time 0 first, then the wanted post-advance value
				// second -- take the LAST match. A segment restarting via
				// RestartSourceTvForPsi from a PREVIOUS chained step (which
				// already carries accumulated time forward) reports the
				// wanted value FIRST, then a further, unwanted one-more-
				// advance value second -- take the FIRST match. Confirmed
				// empirically by diffing against solve()'s monolithic
				// result; see project_tdmrg_evolving_bath memory.
				const bool           takeLast = (particleSrcTv < 0);
				Dmrg::CmdLineOptions opts;
				opts.logfile = chainRoot + "particle_step" + ttos(n) + ".log";
				opts.in_situ_measurements = "<P2|c|P1>,<P2.last|P2>";
				DmrgRunnerType runner(app_,
				                      buildStepInput(params_.uFinal,
				                                     hoppings,
				                                     potTdmrg,
				                                     nup_,
				                                     ndown_,
				                                     nsites,
				                                     particleRoot,
				                                     particleMapTv,
				                                     particleSrcTv,
				                                     outRoot),
				                      opts);
				runner.doOneRun();

				ComplexType ggt(0), gauge(0);
				parseSingleMeasurement(opts.logfile, "<P2|c|P1>", ggt, takeLast);
				if (parseSingleMeasurement(
				        opts.logfile, "<P2.last|P2>", gauge, takeLast))
					gaugeP[n] = gauge;
				result.ggt0[n] = ComplexType(0, -1) * ggt;

				particleRoot  = outRoot;
				particleMapTv = 1; // next step continues from this step's P1
				particleSrcTv = 2; // next step's |gs> seed is this step's P2
			}
			{
				const std::string outRoot = chainRoot + "hole" + ttos(n);
				const bool takeLast = (holeSrcTv < 0); // see particle branch above
				Dmrg::CmdLineOptions opts;
				opts.logfile = chainRoot + "hole_step" + ttos(n) + ".log";
				opts.in_situ_measurements = "<P1|c|P2>,<P2.last|P2>";
				DmrgRunnerType runner(app_,
				                      buildStepInput(params_.uFinal,
				                                     hoppings,
				                                     potTdmrg,
				                                     nup_,
				                                     ndown_,
				                                     nsites,
				                                     holeRoot,
				                                     holeMapTv,
				                                     holeSrcTv,
				                                     outRoot),
				                      opts);
				runner.doOneRun();

				ComplexType glt(0), gauge(0);
				parseSingleMeasurement(opts.logfile, "<P1|c|P2>", glt, takeLast);
				if (parseSingleMeasurement(
				        opts.logfile, "<P2.last|P2>", gauge, takeLast))
					gaugeH[n] = gauge;
				result.glt0[n] = ComplexType(0, 1) * glt;

				holeRoot  = outRoot;
				holeMapTv = 1;
				holeSrcTv = 2;
			}
		}

		applySignFlip(result.ggt0);
		for (auto& kv : result.ggt0) {
			auto it = gaugeP.find(kv.first);
			if (it != gaugeP.end() && std::abs(it->second) > RealType(1e-10))
				kv.second /= it->second;
		}

		applySignFlip(result.glt0);
		for (auto& kv : result.glt0) {
			auto it = gaugeH.find(kv.first);
			if (it != gaugeH.end() && std::abs(it->second) > RealType(1e-10))
				kv.second *= it->second;
		}

		// Note: applyGlobalPhase (existing helper) anchors at n=0, which this
		// chain never measures directly (Run "particle_init"/"hole_init" only
		// prepare states, they don't measure). Left unnormalized here
		// deliberately -- the Phase 1 test compares against the monolithic
		// path's OWN un-normalized-at-this-stage values, or compares
		// magnitudes, rather than assuming a shared phase convention.

		return result;
	}

	// ---- Phase 1 (evolving-bath project): genuine multi-column full grid --
	//
	// Generalizes solveChainedColumn0 above from "column 0 only" to a
	// genuine two-time G(t_n,t_j) grid: one MPS trajectory per insertion
	// time t_j (a "column"), chained forward in time, following the
	// architecture in /Users/epd/.claude/plans/fancy-painting-moon.md.
	// Still verification-only -- not wired into solve()/computeGimp()'s
	// existing dispatch.
	//
	// Static bath only (no NeqBathDecomposition wiring yet -- that is
	// Phase 2). Column 0 doubles as the perpetual seed source for birthing
	// every later column (no separate shared reference thread needed --
	// see plan file for why).
	//
	// KNOWN GAP, not silently dropped: does NOT fill the diagonal
	// gimp.lesser(n,n)/retarded(n,n). Getting the equal-time diagonal
	// requires either a zero-length TimeEvolve measurement or bundling an
	// extra advance into birth, both added mechanisms not yet validated
	// against the monolithic reference the way the off-diagonal chaining
	// below has been. Left as follow-up work; see project_tdmrg_evolving_bath
	// memory. NeqLatticeGf::updateDelta DOES read the diagonal (it feeds
	// NeqBathDecomposition's target `d` for the Cholesky column), so this
	// gap must be closed before Phase 2 can be correct, not just before it
	// is "nice to have".
	KBType computeFullGrid(const VectorRealType& bathParams) const
	{
		const SizeType nBath  = bathParams.size() / 2;
		const SizeType nsites = nBath + 1;

		VectorRealType hoppings(nBath), bathEps(nBath);
		for (SizeType i = 0; i < nBath; ++i) {
			hoppings[i] = bathParams[i];
			bathEps[i]  = bathParams[nBath + i];
		}

		VectorRealType potGS(nsites), potTdmrg(nsites);
		potGS[0]    = -RealType(0.5) * params_.uInitial;
		potTdmrg[0] = -RealType(0.5) * params_.uFinal;
		for (SizeType i = 0; i < nBath; ++i) {
			potGS[i + 1]    = bathEps[i];
			potTdmrg[i + 1] = bathEps[i];
		}

		const std::string chainRoot = root_ + "grid_";
		const int         nT        = static_cast<int>(params_.nT);

		{
			Dmrg::CmdLineOptions opts;
			opts.logfile = chainRoot + "gs.log";
			DmrgRunnerType runner(app_,
			                      buildGsInputAt(chainRoot + "gs",
			                                     params_.uInitial,
			                                     hoppings,
			                                     potGS,
			                                     nup_,
			                                     ndown_,
			                                     nsites),
			                      opts);
			runner.doOneRun();
		}

		std::vector<Column> columns;
		columns.reserve(static_cast<SizeType>(nT));
		columns.emplace_back();
		columns[0].born = 0;
		birthColumn(columns[0],
		            chainRoot + "gs",
		            -1,
		            chainRoot + "gs",
		            -1,
		            chainRoot,
		            "column0",
		            params_.uFinal,
		            hoppings,
		            potTdmrg,
		            nsites);

		for (int n = 1; n <= nT; ++n) {
			// 1. Advance column 0 -- both extends its own trajectory (giving
			//    G(n,0)) and produces the bare reference at t_n that step 2
			//    below needs to birth column n.
			advanceColumn(
			    columns[0], n, chainRoot, params_.uFinal, hoppings, potTdmrg, nsites);

			// 2. Birth column n from column 0's just-refreshed bare
			//    reference, unless n==nT (nothing would ever advance it).
			SizeType justBornIdx = columns.size(); // sentinel: none born this step
			if (n < nT) {
				columns.emplace_back();
				columns.back().born = n;
				birthColumn(columns.back(),
				            columns[0].particleRoot,
				            columns[0].particleSrcTv,
				            columns[0].holeRoot,
				            columns[0].holeSrcTv,
				            chainRoot,
				            "column" + ttos(n),
				            params_.uFinal,
				            hoppings,
				            potTdmrg,
				            nsites);
				justBornIdx = columns.size() - 1;
			}

			// 3. Advance every OTHER existing column (not column 0, already
			//    done in step 1; not the one just born in step 2, which has
			//    nothing to advance yet).
			for (SizeType idx = 1; idx < columns.size(); ++idx) {
				if (idx == justBornIdx)
					continue;
				advanceColumn(columns[idx],
				              n,
				              chainRoot,
				              params_.uFinal,
				              hoppings,
				              potTdmrg,
				              nsites);
			}
		}

		// Column nT is never birthed inside the loop above (birth only
		// happens when n < nT, since a column born at the very last step
		// would never be advanced further for OFF-diagonal purposes) --
		// but its diagonal G(nT,nT) is still wanted. Birth it from column
		// 0's FINAL state, then advance it ONCE (to a throwaway step
		// beyond nT, never written into fullGimp -- see the n=j+1..nT loop
		// below, which is empty for j=nT) purely to trigger the
		// diagonal-capture byproduct in advanceColumn's first-advance
		// path.
		if (nT > 0) {
			columns.emplace_back();
			columns.back().born = nT;
			birthColumn(columns.back(),
			            columns[0].particleRoot,
			            columns[0].particleSrcTv,
			            columns[0].holeRoot,
			            columns[0].holeSrcTv,
			            chainRoot,
			            "column" + ttos(nT),
			            params_.uFinal,
			            hoppings,
			            potTdmrg,
			            nsites);
			advanceColumn(columns.back(),
			              nT + 1,
			              chainRoot,
			              params_.uFinal,
			              hoppings,
			              potTdmrg,
			              nsites);
		}

		KBType fullGimp(params_.nT,
		                params_.eqParams.nMatsubaras,
		                params_.dt,
		                params_.eqParams.ficticiousBeta
		                    / static_cast<RealType>(params_.eqParams.nMatsubaras));

		for (auto& col : columns) {
			applySignFlip(col.ggtRaw);
			applySignFlip(col.gltRaw);

			const int j = col.born;

			// Diagonal G(j,j): captured once at birth, gauge-invariant by
			// construction (see Column::ggtDiag/gltDiag doc comment) --
			// no sign-flip/phase correction needed, unlike the
			// off-diagonal history below.
			{
				const ComplexType ggt   = ComplexType(0, -1) * col.ggtDiag;
				const ComplexType glt   = ComplexType(0, 1) * col.gltDiag;
				fullGimp.lesser(j, j)   = glt;
				fullGimp.retarded(j, j) = ggt - glt;
			}

			for (int n = j + 1; n <= nT; ++n) {
				auto itG = col.ggtRaw.find(n);
				auto itL = col.gltRaw.find(n);
				if (itG == col.ggtRaw.end() || itL == col.gltRaw.end())
					continue;
				const ComplexType ggt   = ComplexType(0, -1) * itG->second;
				const ComplexType glt   = ComplexType(0, 1) * itL->second;
				fullGimp.lesser(n, j)   = glt;
				fullGimp.retarded(n, j) = ggt - glt;
				if (j < n)
					fullGimp.lesser(j, n) = -std::conj(glt);
			}
		}

		return fullGimp;
	}

	// ---- Phase 2 (evolving-bath project): extended-geometry fan-out check --
	//
	// Same architecture as computeFullGrid above, but the star geometry is
	// extended with 2L second-bath sites (L "occupied", L "empty"), seeded
	// via the eps-split GS potential validated by
	// measureSecondBathOccupations. Vplus is fixed at 0 for EVERY step here
	// (not just n=0) -- an "inert spectator" test of the larger lattice,
	// the +L electron counts, the eps-split GS, and the recalibrated
	// TSPAdvanceEach in isolation, BEFORE any bath self-consistency (still
	// to come) is added. With Vplus=0 throughout, the 2L sites never
	// couple to anything and every (n,j) must match computeFullGrid's own
	// (unextended) output to the same tolerance -- this is the real test
	// that the fan-out mechanics survive the larger geometry, per the
	// advisor consult recorded in fancy-painting-moon.md (the TSPAdvanceEach
	// and diagonal-occurrence-count risks that motivated this specific
	// check were verified empirically in build/tmp/advanceeach_check/
	// before this method was written).
	// fixedConnectors (default empty): when non-empty (size 2L), used as the
	// second-bath Connectors for EVERY step instead of the all-zero inert
	// default -- added to empirically test whether Column::ggtDiag/gltDiag
	// (captured from a column's first advanceColumn call, see that method's
	// doc comment) actually depends on the segment's own declared Connectors
	// value, which matters for the self-consistent (evolving Vplus) case
	// where the diagonal for row n must be available before Vplus(n+1) is
	// known. See fancy-painting-moon.md, Phase 2, "diagonal Connectors-
	// independence check".
	// rootSuffix (default "gridbath_"): distinguishes the output/checkpoint
	// file prefix across multiple calls made from the SAME test process
	// (e.g. two calls with different fixedConnectors back-to-back) --
	// callers making more than one call against the same root_ MUST pass
	// distinct suffixes, or the second call's DmrgRunner-written checkpoint
	// files collide with the first's (observed as a real, reproducible
	// cross-call contamination: a "column diagonal is independent of
	// Connectors" test calling this twice with the default suffix saw its
	// second run's off-diagonal spuriously equal the first's, only when
	// run as part of the full suite where a THIRD call -- from a different
	// TEST_CASE sharing the same default root_ -- had already written to
	// the same files earlier in the same process; passing a distinct
	// suffix per call eliminates this at the source rather than chasing
	// the exact HDF5/checkpoint reuse mechanism).
	KBType computeFullGridWithInertSecondBath(const VectorRealType&    bathParams,
	                                          SizeType                 L,
	                                          RealType                 eps,
	                                          const VectorComplexType& fixedConnectors
	                                          = VectorComplexType(),
	                                          const std::string& rootSuffix = "gridbath_") const
	{
		const SizeType nBath     = bathParams.size() / 2;
		const SizeType nsites    = nBath + 1;
		const SizeType nsitesExt = nsites + 2 * L;

		VectorRealType hoppings(nBath), bathEps(nBath);
		for (SizeType i = 0; i < nBath; ++i) {
			hoppings[i] = bathParams[i];
			bathEps[i]  = bathParams[nBath + i];
		}

		VectorRealType potGS(nsitesExt, RealType(0)), potTdmrg(nsitesExt, RealType(0));
		potGS[0]    = -RealType(0.5) * params_.uInitial;
		potTdmrg[0] = -RealType(0.5) * params_.uFinal;
		for (SizeType i = 0; i < nBath; ++i) {
			potGS[i + 1]    = bathEps[i];
			potTdmrg[i + 1] = bathEps[i];
		}
		// eps-split GS-only seeding (occupied sites first, empty sites
		// second, matching NeqBathDecomposition's/eqHybDecomp's own
		// convention); potTdmrg's second-bath entries stay 0 (the true,
		// always-eps=0 evolving-bath potential).
		for (SizeType p = 0; p < L; ++p)
			potGS[nsites + p] = -eps;
		for (SizeType p = 0; p < L; ++p)
			potGS[nsites + L + p] = eps;

		SecondBathExt secondBath;
		secondBath.active      = true;
		secondBath.nup         = nup_ + L;
		secondBath.ndown       = ndown_ + L;
		secondBath.connectors  = fixedConnectors.empty()
		     ? VectorComplexType(2 * L, ComplexType(0))
		     : fixedConnectors;
		secondBath.advanceEach = nsitesExt - 2;

		const std::string chainRoot = root_ + rootSuffix;
		const int         nT        = static_cast<int>(params_.nT);

		{
			Dmrg::CmdLineOptions opts;
			opts.logfile = chainRoot + "gs.log";
			DmrgRunnerType runner(app_,
			                      buildGsInputAt(chainRoot + "gs",
			                                     params_.uInitial,
			                                     hoppings,
			                                     potGS,
			                                     secondBath.nup,
			                                     secondBath.ndown,
			                                     nsitesExt,
			                                     secondBath.connectors),
			                      opts);
			runner.doOneRun();
		}

		std::vector<Column> columns;
		columns.reserve(static_cast<SizeType>(nT));
		columns.emplace_back();
		columns[0].born = 0;
		birthColumn(columns[0],
		            chainRoot + "gs",
		            -1,
		            chainRoot + "gs",
		            -1,
		            chainRoot,
		            "column0",
		            params_.uFinal,
		            hoppings,
		            potTdmrg,
		            nsitesExt,
		            secondBath);

		for (int n = 1; n <= nT; ++n) {
			advanceColumn(columns[0],
			              n,
			              chainRoot,
			              params_.uFinal,
			              hoppings,
			              potTdmrg,
			              nsitesExt,
			              secondBath);

			SizeType justBornIdx = columns.size();
			if (n < nT) {
				columns.emplace_back();
				columns.back().born = n;
				birthColumn(columns.back(),
				            columns[0].particleRoot,
				            columns[0].particleSrcTv,
				            columns[0].holeRoot,
				            columns[0].holeSrcTv,
				            chainRoot,
				            "column" + ttos(n),
				            params_.uFinal,
				            hoppings,
				            potTdmrg,
				            nsitesExt,
				            secondBath);
				justBornIdx = columns.size() - 1;
			}

			for (SizeType idx = 1; idx < columns.size(); ++idx) {
				if (idx == justBornIdx)
					continue;
				advanceColumn(columns[idx],
				              n,
				              chainRoot,
				              params_.uFinal,
				              hoppings,
				              potTdmrg,
				              nsitesExt,
				              secondBath);
			}
		}

		if (nT > 0) {
			columns.emplace_back();
			columns.back().born = nT;
			birthColumn(columns.back(),
			            columns[0].particleRoot,
			            columns[0].particleSrcTv,
			            columns[0].holeRoot,
			            columns[0].holeSrcTv,
			            chainRoot,
			            "column" + ttos(nT),
			            params_.uFinal,
			            hoppings,
			            potTdmrg,
			            nsitesExt,
			            secondBath);
			advanceColumn(columns.back(),
			              nT + 1,
			              chainRoot,
			              params_.uFinal,
			              hoppings,
			              potTdmrg,
			              nsitesExt,
			              secondBath);
		}

		KBType fullGimp(params_.nT,
		                params_.eqParams.nMatsubaras,
		                params_.dt,
		                params_.eqParams.ficticiousBeta
		                    / static_cast<RealType>(params_.eqParams.nMatsubaras));

		for (auto& col : columns) {
			applySignFlip(col.ggtRaw);
			applySignFlip(col.gltRaw);

			const int j = col.born;

			{
				const ComplexType ggt   = ComplexType(0, -1) * col.ggtDiag;
				const ComplexType glt   = ComplexType(0, 1) * col.gltDiag;
				fullGimp.lesser(j, j)   = glt;
				fullGimp.retarded(j, j) = ggt - glt;
			}

			for (int n = j + 1; n <= nT; ++n) {
				auto itG = col.ggtRaw.find(n);
				auto itL = col.gltRaw.find(n);
				if (itG == col.ggtRaw.end() || itL == col.gltRaw.end())
					continue;
				const ComplexType ggt   = ComplexType(0, -1) * itG->second;
				const ComplexType glt   = ComplexType(0, 1) * itL->second;
				fullGimp.lesser(n, j)   = glt;
				fullGimp.retarded(n, j) = ggt - glt;
				if (j < n)
					fullGimp.lesser(j, n) = -std::conj(glt);
			}
		}

		return fullGimp;
	}

	// ---- Phase 2 (evolving-bath project): second-bath seeding check --------
	//
	// GBEK's second bath is 2L extra star-geometry sites (L "occupied", L
	// "empty") coupled to the impurity with amplitude Vplus(n,p), which is
	// exactly 0 at n=0 (NeqBathDecomposition skips n=0) -- so at the point
	// where tDMRG's GS run would diagonalize, these 2L sites are EXACTLY
	// decoupled from the rest of the system. DMRG's Lanczos GS solver then
	// faces an exactly degenerate manifold across every occupation
	// combination among them, with no natural mechanism (unlike GBEK's own
	// explicit ED basis-vector construction) to land on the specific
	// L-occupied/L-empty configuration NeqBathDecomposition expects.
	//
	// Fix: break the degeneracy with a small potential split at the GS run
	// only (-eps on the intended-occupied sites, +eps on the intended-empty
	// sites); the true evolving-bath potential is 0 for all these sites at
	// every subsequent time step, so this split is purely a GS-seeding
	// device, not part of the physical model.
	//
	// IMPORTANT, confirmed empirically (see project_tdmrg_evolving_bath
	// memory / fancy-painting-moon.md plan): eps must be scaled to the
	// largest OTHER coupling actually present (first-bath hoppings, U),
	// not just placed above the Lanczos convergence floor. A too-small eps
	// lets the coupled impurity/first-bath subsystem's OWN delocalization
	// energy gain outcompete the split, pulling electrons away from the
	// intended second-bath occupation pattern instead of leaving them
	// there. Callers must choose eps accordingly (e.g. several times
	// max(|hoppings|, |bathEps|, U)) and verify with this very method --
	// do not assume a fixed eps works across parameter regimes.
	//
	// Returns the measured <n_p> for each of the 2L second-bath sites,
	// ordered [occupied_0..occupied_{L-1}, empty_0..empty_{L-1}] to match
	// NeqBathDecomposition's/eqHybDecomp's own convention (first half
	// occupied, second half unoccupied). Verification-only: builds one
	// extra, throwaway GS run; not wired into solve()/computeFullGrid.
	VectorRealType measureSecondBathOccupations(const VectorRealType& bathParams,
	                                            SizeType              L,
	                                            RealType              eps) const
	{
		const SizeType nBath     = bathParams.size() / 2;
		const SizeType nsites    = nBath + 1;
		const SizeType nsitesExt = nsites + 2 * L;

		VectorRealType hoppings(nBath), bathEps(nBath);
		for (SizeType i = 0; i < nBath; ++i) {
			hoppings[i] = bathParams[i];
			bathEps[i]  = bathParams[nBath + i];
		}

		// Extended Connectors: first-bath hoppings unchanged, second-bath
		// hoppings all zero (Vplus(0,p)=0 by construction at t=0).
		VectorRealType hoppingsExt(nsitesExt - 1, RealType(0));
		for (SizeType i = 0; i < nBath; ++i)
			hoppingsExt[i] = hoppings[i];

		// Extended potentialV: impurity/first-bath unchanged, then L
		// occupied sites at -eps, L empty sites at +eps.
		VectorRealType potExt(nsitesExt, RealType(0));
		potExt[0] = -RealType(0.5) * params_.uInitial;
		for (SizeType i = 0; i < nBath; ++i)
			potExt[i + 1] = bathEps[i];
		for (SizeType p = 0; p < L; ++p)
			potExt[nsites + p] = -eps;
		for (SizeType p = 0; p < L; ++p)
			potExt[nsites + L + p] = eps;

		const SizeType nupExt   = nup_ + L;
		const SizeType ndownExt = ndown_ + L;

		const std::string outRoot = root_ + "secondbath_seed_gs";

		std::string insitu;
		for (SizeType p = 0; p < 2 * L; ++p) {
			if (p > 0)
				insitu += ",";
			insitu += "<gs|n[" + ttos(nsites + p) + "]|gs>";
		}

		Dmrg::CmdLineOptions opts;
		opts.logfile              = outRoot + ".log";
		opts.in_situ_measurements = insitu;
		DmrgRunnerType runner(app_,
		                      buildGsInputAt(outRoot,
		                                     params_.uInitial,
		                                     hoppingsExt,
		                                     potExt,
		                                     nupExt,
		                                     ndownExt,
		                                     nsitesExt),
		                      opts);
		runner.doOneRun();

		VectorRealType occ(2 * L);
		for (SizeType p = 0; p < 2 * L; ++p) {
			const std::string label = "<gs|n[" + ttos(nsites + p) + "]|gs>";
			ComplexType       val(0);
			parseSingleMeasurementAtSite(opts.logfile, label, nsites + p, val, true);
			occ[p] = std::real(val);
		}
		return occ;
	}

private:

	// ---- Phase 1 multi-column state -----------------------------------------
	//
	// One "column" = one insertion time t_j = born*dt, tracked per channel
	// (particle/hole) as: the current checkpoint root, which TV index to
	// restart-map for the NEXT segment's P0 (mappedTv) and which TV index to
	// use as RestartSourceTvForPsi for the NEXT segment's |gs> seed (srcTv;
	// -1 means "use restartRoot's own natural single state", true right
	// after birth), and the raw measurement history keyed by outer step n.
	// applySignFlip (a discrete, magnitude-based heuristic, unrelated to the
	// phase-division gauge correction removed below) still needs a column's
	// own full, chronologically-ordered history, hence deferred to
	// computeFullGrid's final pass rather than applied per-step.
	//
	// No <P2.last|P2>-based phase-gauge correction is stored here (removed
	// -- see advanceColumn's comment for the empirical finding: it was a
	// no-op in every Phase 1 test and actively wrong once tested against an
	// eps-split extended geometry).
	struct Column {
		int                        born = 0;
		std::string                particleRoot, holeRoot;
		int                        particleMapTv = 0, particleSrcTv = -1;
		int                        holeMapTv = 0, holeSrcTv = -1;
		std::map<int, ComplexType> ggtRaw; // keyed by n, n > born
		std::map<int, ComplexType> gltRaw;
		// Equal-time diagonal G(born,born), captured once at birth (see
		// birthColumn) -- gauge-invariant by construction, no sign-flip/
		// phase correction needed (unlike ggtRaw/gltRaw above).
		ComplexType ggtDiag = 0, gltDiag = 0;
	};

	// Birth a new column at time col.born. The particle and hole branches
	// restart from SEPARATE sources -- (particleSourceRoot, particleSourceTv)
	// and (holeSourceRoot, holeSourceTv) -- because a source column's
	// particle-channel and hole-channel checkpoint chains are independently
	// propagated (see ImpuritySolverNeqTdmrg's class doc: Run2->Run3 vs
	// Run4->Run5), not interchangeable. No measurement, no time advance --
	// exactly today's buildParticleInitInput/buildHoleInitInput pattern (see
	// buildInitInputAt), just parameterized to restart from an arbitrary
	// column's current checkpoint instead of always from the t=0 gs run.
	// After birth, col.particleRoot/holeRoot hold single-TV checkpoints, so
	// the column's own first subsequent advanceColumn call is structurally
	// identical to column 0's very first step (mapTv=0, srcTv=-1) --
	// exactly the case solveChainedColumn0 already validated needs
	// takeLast=true.
	// Bundles the per-call second-bath (evolving-bath) extension for
	// birthColumn/advanceColumn, used only when neqBathRank_ > 0.
	// Default-constructed (active==false) reproduces today's behavior
	// exactly: nup_/ndown_ members and tspAdvanceEach_ as read from the
	// input file, no second-bath Connectors appended. A dedicated struct
	// (rather than more trailing optional parameters) avoids sentinel-value
	// ambiguity for nup/ndown/advanceEach, none of which have a safe
	// "unused" value of their own.
	struct SecondBathExt {
		bool              active = false;
		SizeType          nup    = 0;
		SizeType          ndown  = 0;
		VectorComplexType connectors; // size 2L, duplicated occ/empty per p
		SizeType          advanceEach = 0; // nsitesExt - 2, recalibrated
	};

	void birthColumn(Column&               col,
	                 const std::string&    particleSourceRoot,
	                 int                   particleSourceTv,
	                 const std::string&    holeSourceRoot,
	                 int                   holeSourceTv,
	                 const std::string&    chainRoot,
	                 const std::string&    tag,
	                 RealType              uFinal,
	                 const VectorRealType& hoppings,
	                 const VectorRealType& potTdmrg,
	                 SizeType              nsites,
	                 const SecondBathExt&  secondBath = SecondBathExt()) const
	{
		const SizeType nup   = secondBath.active ? secondBath.nup : nup_;
		const SizeType ndown = secondBath.active ? secondBath.ndown : ndown_;

		col.particleRoot = chainRoot + tag + "_particle_birth";
		{
			Dmrg::CmdLineOptions opts;
			opts.logfile = col.particleRoot + ".log";
			DmrgRunnerType runner(app_,
			                      buildInitInputAt(col.particleRoot,
			                                       particleSourceRoot,
			                                       uFinal,
			                                       hoppings,
			                                       potTdmrg,
			                                       nup,
			                                       ndown,
			                                       nsites,
			                                       "'",
			                                       particleSourceTv,
			                                       particleSourceTv >= 0,
			                                       secondBath.connectors),
			                      opts);
			runner.doOneRun();
		}
		col.particleMapTv = 0;
		col.particleSrcTv = -1;

		col.holeRoot = chainRoot + tag + "_hole_birth";
		{
			Dmrg::CmdLineOptions opts;
			opts.logfile = col.holeRoot + ".log";
			DmrgRunnerType runner(app_,
			                      buildInitInputAt(col.holeRoot,
			                                       holeSourceRoot,
			                                       uFinal,
			                                       hoppings,
			                                       potTdmrg,
			                                       nup,
			                                       ndown,
			                                       nsites,
			                                       "",
			                                       holeSourceTv,
			                                       holeSourceTv >= 0,
			                                       secondBath.connectors),
			                      opts);
			runner.doOneRun();
		}
		col.holeMapTv = 0;
		col.holeSrcTv = -1;
	}

	// Advance column col by one step, to outer step n (n = col's
	// last-advanced step + 1, tracked implicitly by the caller's loop
	// order). Fills col.ggtRaw[n]/gltRaw[n] (raw, not yet gauge-corrected --
	// see computeFullGrid). takeLast mirrors solveChainedColumn0's
	// empirically-validated rule: true only for a column's own first
	// advance since birth (particleSrcTv/holeSrcTv < 0), false thereafter.
	void advanceColumn(Column&               col,
	                   int                   n,
	                   const std::string&    chainRoot,
	                   RealType              uFinal,
	                   const VectorRealType& hoppings,
	                   const VectorRealType& potTdmrg,
	                   SizeType              nsites,
	                   const SecondBathExt&  secondBath = SecondBathExt()) const
	{
		const SizeType nup   = secondBath.active ? secondBath.nup : nup_;
		const SizeType ndown = secondBath.active ? secondBath.ndown : ndown_;

		const std::string tag = "j" + ttos(col.born) + "_n" + ttos(n);
		{
			const bool           takeLast = (col.particleSrcTv < 0);
			const std::string    outRoot  = chainRoot + tag + "_particle";
			Dmrg::CmdLineOptions opts;
			opts.logfile              = outRoot + ".log";
			opts.in_situ_measurements = "<P2|c|P1>";
			DmrgRunnerType runner(app_,
			                      buildStepInput(uFinal,
			                                     hoppings,
			                                     potTdmrg,
			                                     nup,
			                                     ndown,
			                                     nsites,
			                                     col.particleRoot,
			                                     col.particleMapTv,
			                                     col.particleSrcTv,
			                                     outRoot,
			                                     secondBath.connectors,
			                                     secondBath.advanceEach),
			                      opts);
			runner.doOneRun();

			ComplexType ggt(0);
			parseSingleMeasurement(opts.logfile, "<P2|c|P1>", ggt, takeLast);
			// No <P2.last|P2>-based gauge correction (removed -- see
			// fancy-painting-moon.md, "gauge-correction finding"): <P2|c|P1>'s
			// bra and ket both descend from the SAME loaded |gs> reference in
			// every segment (whether this is a column's first advance since
			// birth or a later one restarted via RestartSourceTvForPsi from
			// its own prior segment), so any arbitrary phase that reference
			// carries cancels in the bra-ket already. Dividing by
			// <P2.last|P2> anyway was a no-op in every Phase 1 test (it
			// measured exactly 1 there) but actively wrong once tested
			// against an eps-split extended geometry (Phase 2): confirmed
			// empirically that <P2.last|P2> is then a large,
			// sweep-position-independent and step-independent phase (DMRG's
			// own arbitrary basis choice for a now more-degenerate reduced
			// density matrix, not a physical dynamical phase), and dividing
			// by it introduced a spurious rotation into an otherwise-correct
			// raw measurement. Both extended-GS and unextended-GS sector
			// energies were checked to decompose cleanly (differ by exactly
			// -2*eps*L), confirming the impurity block itself is undisturbed
			// and the phase is pure basis-choice arbitrariness, not physics.
			col.ggtRaw[n] = ggt;

			// A column's own FIRST advance since birth (takeLast==true) is
			// restarting from a plain, single-TV birth checkpoint with NO
			// RestartSourceTvForPsi involved -- structurally identical to
			// column 0's very own first advance. The FIRST occurrence in
			// THAT specific case is the equal-time diagonal G(born,born)
			// (confirmed empirically: matches ImpuritySolverNeqExactDiag's
			// G(0,0) exactly for column 0 -- see project_tdmrg_evolving_bath
			// memory). Capture it here instead of via a separate, fragile
			// birth-time measurement.
			if (takeLast)
				parseSingleMeasurement(
				    opts.logfile, "<P2|c|P1>", col.ggtDiag, false);

			col.particleRoot  = outRoot;
			col.particleMapTv = 1;
			col.particleSrcTv = 2;
		}
		{
			const bool           takeLast = (col.holeSrcTv < 0);
			const std::string    outRoot  = chainRoot + tag + "_hole";
			Dmrg::CmdLineOptions opts;
			opts.logfile              = outRoot + ".log";
			opts.in_situ_measurements = "<P1|c|P2>";
			DmrgRunnerType runner(app_,
			                      buildStepInput(uFinal,
			                                     hoppings,
			                                     potTdmrg,
			                                     nup,
			                                     ndown,
			                                     nsites,
			                                     col.holeRoot,
			                                     col.holeMapTv,
			                                     col.holeSrcTv,
			                                     outRoot,
			                                     secondBath.connectors,
			                                     secondBath.advanceEach),
			                      opts);
			runner.doOneRun();

			ComplexType glt(0);
			parseSingleMeasurement(opts.logfile, "<P1|c|P2>", glt, takeLast);
			// No gauge correction -- see the matching comment in the
			// particle branch above.
			col.gltRaw[n] = glt;

			if (takeLast)
				parseSingleMeasurement(
				    opts.logfile, "<P1|c|P2>", col.gltDiag, false);

			col.holeRoot  = outRoot;
			col.holeMapTv = 1;
			col.holeSrcTv = 2;
		}
	}

	// ---- Phase 2 (evolving-bath project): self-consistent solve() setup ---
	//
	// One-time (per solve() call) setup only: delegates the equilibrium/
	// Matsubara solve to ExactDiag (mirroring ImpuritySolverNeqGBEK::solve,
	// which does the same for the same reason -- tDMRG has no equilibrium
	// solver of its own), constructs decomp_, and stores the extended-
	// geometry parameters fillSelfConsistentRow needs. Deliberately does NOT
	// run any DMRG passes itself -- fillSelfConsistentRow rebuilds the GS
	// run and the whole column fan-out from scratch on every computeGimp
	// call (see that method's doc comment for why this is correct, not just
	// simple).
	void solveSelfConsistent(const VectorRealType& bathParams)
	{
		exactDiag_.solve(bathParams);

		const RealType beta = params_.eqParams.ficticiousBeta;
		const RealType mu   = 0;
		decomp_             = std::make_unique<DecompType>(
                    neqBathRank_,
                    beta,
                    mu,
                    bathParams,
                    params_.nT,
                    params_.eqParams.nMatsubaras,
                    params_.dt,
                    params_.eqParams.ficticiousBeta
                        / static_cast<RealType>(params_.eqParams.nMatsubaras));

		const SizeType nBath  = bathParams.size() / 2;
		const SizeType nsites = nBath + 1;
		scL_                  = neqBathRank_;
		scNsitesExt_          = nsites + 2 * scL_;
		scNup_                = nup_ + scL_;
		scNdown_              = ndown_ + scL_;

		scHoppings_.resize(nBath);
		scBathEps_.resize(nBath);
		for (SizeType i = 0; i < nBath; ++i) {
			scHoppings_[i] = bathParams[i];
			scBathEps_[i]  = bathParams[nBath + i];
		}

		scPotTdmrg_.assign(scNsitesExt_, RealType(0));
		scPotTdmrg_[0] = -RealType(0.5) * params_.uFinal;
		for (SizeType i = 0; i < nBath; ++i)
			scPotTdmrg_[i + 1] = scBathEps_[i];
		// Second-bath entries of scPotTdmrg_ stay 0: the true evolving-bath
		// potential is eps=0 for every time-evolution segment (see
		// measureSecondBathOccupations' doc comment) -- eps only applies to
		// the one-off GS seeding run rebuilt fresh in fillSelfConsistentRow.

		// eps calibration rule confirmed empirically (blocker B, see
		// fancy-painting-moon.md): must scale to the largest OTHER coupling
		// present, not a fixed small constant.
		RealType maxCoupling = std::max(params_.uInitial, params_.uFinal);
		for (SizeType i = 0; i < nBath; ++i) {
			maxCoupling = std::max(maxCoupling, std::abs(scHoppings_[i]));
			maxCoupling = std::max(maxCoupling, std::abs(scBathEps_[i]));
		}
		if (maxCoupling == RealType(0))
			maxCoupling = RealType(1);
		scEps_ = RealType(5) * maxCoupling;

		scChainRoot_ = root_ + "sc_";
	}

	// ---- Phase 2 (evolving-bath project): self-consistent row fill --------
	//
	// Fills gimp's retarded/lesser (n,j) for j=0..n (and the anti-Hermitian
	// transpose lesser(j,n) for j<n) by rebuilding the FULL fan-out from
	// t=0 through step n fresh -- same Column/birthColumn/advanceColumn
	// mechanics as computeFullGridWithInertSecondBath, but bounded by n
	// (not params_.nT) and reading the second-bath Connectors at each outer
	// step k from decomp_->Vplus(k-1,p)/Vplus(k,p) (midpoint-averaged,
	// mirroring ImpuritySolverNeqGBEK::computeVMid) instead of a fixed 0.
	//
	// Deliberately does NOT persist columns_ or any progress marker across
	// calls: NeqDmftSolver's corrector loop (see NeqDmftSolver::timeStep)
	// calls prepareTimeStep(n,...)/computeGimp(gimp,n) repeatedly for the
	// SAME n as decomp_ is refined, and a full recompute from decomp_'s
	// CURRENT state is trivially correct for that without any propagated-
	// through/rewind bookkeeping. Cost is O(n^2) segments per call, O(nT^3)
	// total over a run -- accepted for now (see the advisor consult
	// recorded in fancy-painting-moon.md): incremental/persistent-state
	// optimization is deferred, to be built later against this method's own
	// passing gate as its regression net, not attempted alongside it.
	//
	// The newest column (born at step n, never naturally advanced within
	// this bounded recompute since there is no step n+1 here) still needs
	// its diagonal G(n,n): mirrors computeFullGrid's original handling of
	// column nT exactly -- birth it, then one extra "throwaway" advance
	// purely to trigger the diagonal-capture byproduct (see advanceColumn),
	// discarding the resulting off-diagonal value. This throwaway advance's
	// own Connectors value provably does not matter (confirmed empirically:
	// see the "column diagonal is independent of the second-bath Connectors
	// value" test) -- it reuses whatever this step's own vMid already is,
	// rather than needing Vplus(n+1,*), which is not yet determined at this
	// point (prepareTimeStep(n+1,...) has not run).
	void fillSelfConsistentRow(KBType& gimp, int n) const
	{
		const SizeType nBath  = scHoppings_.size();
		const SizeType nsites = nBath + 1;

		VectorRealType potGS(scNsitesExt_, RealType(0));
		potGS[0] = -RealType(0.5) * params_.uInitial;
		for (SizeType i = 0; i < nBath; ++i)
			potGS[i + 1] = scBathEps_[i];
		for (SizeType p = 0; p < scL_; ++p)
			potGS[nsites + p] = -scEps_;
		for (SizeType p = 0; p < scL_; ++p)
			potGS[nsites + scL_ + p] = scEps_;

		auto vMidConnectors = [this](int k)
		{
			VectorComplexType c(2 * scL_);
			for (SizeType p = 0; p < scL_; ++p) {
				const ComplexType vPrev
				    = decomp_->Vplus(k - 1, static_cast<int>(p));
				const ComplexType vCurr = decomp_->Vplus(k, static_cast<int>(p));
				const ComplexType vMid  = RealType(0.5) * (vPrev + vCurr);
				c[p]                    = vMid;
				c[scL_ + p]             = vMid;
			}
			return c;
		};

		{
			Dmrg::CmdLineOptions opts;
			opts.logfile = scChainRoot_ + "gs.log";
			DmrgRunnerType runner(
			    app_,
			    buildGsInputAt(scChainRoot_ + "gs",
			                   params_.uInitial,
			                   scHoppings_,
			                   potGS,
			                   scNup_,
			                   scNdown_,
			                   scNsitesExt_,
			                   VectorComplexType(2 * scL_, ComplexType(0))),
			    opts);
			runner.doOneRun();
		}

		std::vector<Column> columns;
		columns.reserve(static_cast<SizeType>(n) + 1);
		columns.emplace_back();
		columns[0].born = 0;
		birthColumn(columns[0],
		            scChainRoot_ + "gs",
		            -1,
		            scChainRoot_ + "gs",
		            -1,
		            scChainRoot_,
		            "column0",
		            params_.uFinal,
		            scHoppings_,
		            scPotTdmrg_,
		            scNsitesExt_,
		            SecondBathExt { true,
		                            scNup_,
		                            scNdown_,
		                            VectorComplexType(2 * scL_, ComplexType(0)),
		                            scNsitesExt_ - 2 });

		VectorComplexType lastConnectors(2 * scL_, ComplexType(0));
		for (int k = 1; k <= n; ++k) {
			lastConnectors = vMidConnectors(k);
			SecondBathExt secondBath {
				true, scNup_, scNdown_, lastConnectors, scNsitesExt_ - 2
			};

			advanceColumn(columns[0],
			              k,
			              scChainRoot_,
			              params_.uFinal,
			              scHoppings_,
			              scPotTdmrg_,
			              scNsitesExt_,
			              secondBath);

			SizeType justBornIdx = columns.size();
			if (k < n) {
				columns.emplace_back();
				columns.back().born = k;
				birthColumn(columns.back(),
				            columns[0].particleRoot,
				            columns[0].particleSrcTv,
				            columns[0].holeRoot,
				            columns[0].holeSrcTv,
				            scChainRoot_,
				            "column" + ttos(k),
				            params_.uFinal,
				            scHoppings_,
				            scPotTdmrg_,
				            scNsitesExt_,
				            secondBath);
				justBornIdx = columns.size() - 1;
			}

			for (SizeType idx = 1; idx < columns.size(); ++idx) {
				if (idx == justBornIdx)
					continue;
				advanceColumn(columns[idx],
				              k,
				              scChainRoot_,
				              params_.uFinal,
				              scHoppings_,
				              scPotTdmrg_,
				              scNsitesExt_,
				              secondBath);
			}
		}

		if (n > 0) {
			SecondBathExt secondBath {
				true, scNup_, scNdown_, lastConnectors, scNsitesExt_ - 2
			};
			columns.emplace_back();
			columns.back().born = n;
			birthColumn(columns.back(),
			            columns[0].particleRoot,
			            columns[0].particleSrcTv,
			            columns[0].holeRoot,
			            columns[0].holeSrcTv,
			            scChainRoot_,
			            "column" + ttos(n),
			            params_.uFinal,
			            scHoppings_,
			            scPotTdmrg_,
			            scNsitesExt_,
			            secondBath);
			advanceColumn(columns.back(),
			              n + 1,
			              scChainRoot_,
			              params_.uFinal,
			              scHoppings_,
			              scPotTdmrg_,
			              scNsitesExt_,
			              secondBath);
		}

		for (auto& col : columns) {
			applySignFlip(col.ggtRaw);
			applySignFlip(col.gltRaw);

			const int j = col.born;
			if (j != n)
				continue; // only need row n's own diagonal from this column

			const ComplexType ggtD = ComplexType(0, -1) * col.ggtDiag;
			const ComplexType gltD = ComplexType(0, 1) * col.gltDiag;
			gimp.lesser(j, j)      = gltD;
			gimp.retarded(j, j)    = ggtD - gltD;
		}

		for (auto& col : columns) {
			const int j = col.born;
			if (j == n)
				continue;
			auto itG = col.ggtRaw.find(n);
			auto itL = col.gltRaw.find(n);
			if (itG == col.ggtRaw.end() || itL == col.gltRaw.end())
				continue;
			const ComplexType ggt = ComplexType(0, -1) * itG->second;
			const ComplexType glt = ComplexType(0, 1) * itL->second;
			gimp.lesser(n, j)     = glt;
			gimp.retarded(n, j)   = ggt - glt;
			gimp.lesser(j, n)     = -std::conj(glt);
		}
	}

	// ---- Input construction ------------------------------------------------

	std::string buildGsInput(RealType              U,
	                         const VectorRealType& hoppings,
	                         const VectorRealType& potV,
	                         SizeType              nup,
	                         SizeType              ndown,
	                         SizeType              nsites) const
	{
		std::string s = "##Ainur1.0\n\n";
		s += geomHeader(nsites, U);
		s += "SolverOptions=twositedmrg,geometryallinsystem;\n";
		s += "Version=neqTdmrg;\n";
		s += "OutputFile=" + root_ + "gs;\n";
		s += "InfiniteLoopKeptStates=" + ttos(infiniteLoops_) + ";\n";
		s += "FiniteLoops=" + finiteLoopsGs_ + ";\n";
		s += "TargetElectronsUp=" + ttos(nup) + ";\n";
		s += "TargetElectronsDown=" + ttos(ndown) + ";\n";
		s += "dir0:Connectors=" + buildConnectorsStr(hoppings) + ";\n";
		s += "potentialV=" + buildPotentialVStr(potV) + ";\n";
		return s;
	}

	// Run 2: apply c'[0] to |GS_i⟩.
	// FiniteLoops flag=2 prevents DMRG from re-optimising |gs⟩ to |GS_f⟩,
	// keeping it as the WFT-transformed |GS_i⟩ throughout.
	std::string buildParticleInitInput(RealType              U_f,
	                                   const VectorRealType& hoppings,
	                                   const VectorRealType& potV,
	                                   SizeType              nup,
	                                   SizeType              ndown,
	                                   SizeType              nsites) const
	{
		std::string s = "##Ainur1.0\n\n";
		s += geomHeader(nsites, U_f);
		s += "SolverOptions=twositedmrg,geometryallinsystem,TargetingExpression,restart;\n";
		s += "Version=neqTdmrg;\n";
		s += "OutputFile=" + root_ + "particle;\n";
		s += "InfiniteLoopKeptStates=" + ttos(infiniteLoops_) + ";\n";
		s += "FiniteLoops=" + enforceFlag2(finiteLoopsGs_) + ";\n";
		s += "TargetElectronsUp=" + ttos(nup) + ";\n";
		s += "TargetElectronsDown=" + ttos(ndown) + ";\n";
		s += "dir0:Connectors=" + buildConnectorsStr(hoppings) + ";\n";
		s += "potentialV=" + buildPotentialVStr(potV) + ";\n";
		s += "RestartFilename=" + root_ + "gs;\n";
		s += "GsWeight=0.1;\n";
		s += "string P0=\"c'[0]*|gs>\";\n";
		return s;
	}

	// Run 3: particle tDMRG. P1=e^{-iH_f t}c†|GS_i⟩, P2=e^{-iH_f t}|GS_i⟩.
	// <P2|c|P1> = <GS_i(t)|c|Φ_+(t)⟩ → G^>(t,0) after multiplying by -i.
	// <P2.last|P2> = <GS_i(t-Δt)|GS_i(t)⟩ → N-sector gauge phase correction.
	std::string buildTdmrgInput(RealType              U_f,
	                            const VectorRealType& hoppings,
	                            const VectorRealType& potV,
	                            SizeType              nup,
	                            SizeType              ndown,
	                            SizeType              nsites) const
	{
		std::string s = "##Ainur1.0\n\n";
		s += geomHeader(nsites, U_f);
		s += "SolverOptions=twositedmrg,geometryallinsystem,TargetingExpression,restart,"
		     "usecomplex;\n";
		s += "Version=neqTdmrg;\n";
		s += "OutputFile=" + root_ + "tdmrg;\n";
		s += "InfiniteLoopKeptStates=" + ttos(infiniteLoops_) + ";\n";
		s += "FiniteLoops=" + enforceFlag2(finiteLoopsTdmrg_) + ";\n";
		s += "TargetElectronsUp=" + ttos(nup) + ";\n";
		s += "TargetElectronsDown=" + ttos(ndown) + ";\n";
		s += "dir0:Connectors=" + buildConnectorsStr(hoppings) + ";\n";
		s += "potentialV=" + buildPotentialVStr(potV) + ";\n";
		s += "RestartFilename=" + root_ + "particle;\n";
		s += "RestartMappingTvs=[0, -1, -1];\n";
		s += "GsWeight=0.1;\n";
		s += "string P0=|P0>;\n";
		s += "string P1=\"TimeEvolve{tau=" + ttos(params_.dt) + ",steps="
		    + ttos(tspTimeSteps_) + ",advanceEach=" + ttos(tspAdvanceEach_) + "}*|P0>\";\n";
		s += "string P2=\"TimeEvolve{tau=" + ttos(params_.dt) + ",steps="
		    + ttos(tspTimeSteps_) + ",advanceEach=" + ttos(tspAdvanceEach_) + "}*|gs>\";\n";
		return s;
	}

	// Run 4: apply c[0] to |GS_i⟩.  Mirror of Run 2 with annihilation operator.
	std::string buildHoleInitInput(RealType              U_f,
	                               const VectorRealType& hoppings,
	                               const VectorRealType& potV,
	                               SizeType              nup,
	                               SizeType              ndown,
	                               SizeType              nsites) const
	{
		std::string s = "##Ainur1.0\n\n";
		s += geomHeader(nsites, U_f);
		s += "SolverOptions=twositedmrg,geometryallinsystem,TargetingExpression,restart;\n";
		s += "Version=neqTdmrg;\n";
		s += "OutputFile=" + root_ + "hole;\n";
		s += "InfiniteLoopKeptStates=" + ttos(infiniteLoops_) + ";\n";
		s += "FiniteLoops=" + enforceFlag2(finiteLoopsGs_) + ";\n";
		s += "TargetElectronsUp=" + ttos(nup) + ";\n";
		s += "TargetElectronsDown=" + ttos(ndown) + ";\n";
		s += "dir0:Connectors=" + buildConnectorsStr(hoppings) + ";\n";
		s += "potentialV=" + buildPotentialVStr(potV) + ";\n";
		s += "RestartFilename=" + root_ + "gs;\n";
		s += "GsWeight=0.1;\n";
		s += "string P0=\"c[0]*|gs>\";\n";
		return s;
	}

	// Run 5: hole tDMRG.  P1=e^{-iH_f t}c|GS_i⟩, P2=e^{-iH_f t}|GS_i⟩.
	// <P1|c|P2> = <φ_h(t)|c|GS_i(t)⟩ → G^<(t,0) after multiplying by +i.
	std::string buildHoleTdmrgInput(RealType              U_f,
	                                const VectorRealType& hoppings,
	                                const VectorRealType& potV,
	                                SizeType              nup,
	                                SizeType              ndown,
	                                SizeType              nsites) const
	{
		std::string s = "##Ainur1.0\n\n";
		s += geomHeader(nsites, U_f);
		s += "SolverOptions=twositedmrg,geometryallinsystem,TargetingExpression,restart,"
		     "usecomplex;\n";
		s += "Version=neqTdmrg;\n";
		s += "OutputFile=" + root_ + "tdmrg_hole;\n";
		s += "InfiniteLoopKeptStates=" + ttos(infiniteLoops_) + ";\n";
		s += "FiniteLoops=" + enforceFlag2(finiteLoopsTdmrg_) + ";\n";
		s += "TargetElectronsUp=" + ttos(nup) + ";\n";
		s += "TargetElectronsDown=" + ttos(ndown) + ";\n";
		s += "dir0:Connectors=" + buildConnectorsStr(hoppings) + ";\n";
		s += "potentialV=" + buildPotentialVStr(potV) + ";\n";
		s += "RestartFilename=" + root_ + "hole;\n";
		s += "RestartMappingTvs=[0, -1, -1];\n";
		s += "GsWeight=0.1;\n";
		s += "string P0=|P0>;\n";
		s += "string P1=\"TimeEvolve{tau=" + ttos(params_.dt) + ",steps="
		    + ttos(tspTimeSteps_) + ",advanceEach=" + ttos(tspAdvanceEach_) + "}*|P0>\";\n";
		s += "string P2=\"TimeEvolve{tau=" + ttos(params_.dt) + ",steps="
		    + ttos(tspTimeSteps_) + ",advanceEach=" + ttos(tspAdvanceEach_) + "}*|gs>\";\n";
		return s;
	}

	// ---- Phase 1 diagnostic: parameterized single-step builders -----------
	// Mirror buildGsInput/buildParticleInitInput/buildHoleInitInput/
	// buildTdmrgInput/buildHoleTdmrgInput above, but with an explicit output
	// root (so a chained run doesn't clobber solve()'s own checkpoint files)
	// and, for buildStepInput, an explicit restart source and a SINGLE
	// FiniteLoops row (one time advance) instead of the whole-trajectory
	// finiteLoopsTdmrg_ matrix.

	// secondBathConnectors (default empty): when non-empty, its entries are
	// appended (as complex literals) after hoppings' real entries in
	// dir0:Connectors -- the evolving-bath second-bath couplings, already
	// duplicated once per occupied-p/empty-p site pair by the caller (see
	// fancy-painting-moon.md's Phase 2 architecture section for why no
	// extra conjugation is needed). Empty (the default) reproduces today's
	// NeqBathRank=0 output exactly (buildConnectorsStrWithSecondBath with
	// an empty second argument formats identically to buildConnectorsStr).
	std::string buildGsInputAt(const std::string&       outRoot,
	                           RealType                 U,
	                           const VectorRealType&    hoppings,
	                           const VectorRealType&    potV,
	                           SizeType                 nup,
	                           SizeType                 ndown,
	                           SizeType                 nsites,
	                           const VectorComplexType& secondBathConnectors
	                           = VectorComplexType()) const
	{
		std::string s = "##Ainur1.0\n\n";
		s += geomHeader(nsites, U);
		s += "SolverOptions=twositedmrg,geometryallinsystem";
		s += hasNonzeroImag(secondBathConnectors) ? ",usecomplex;\n" : ";\n";
		s += "Version=neqTdmrg;\n";
		s += "OutputFile=" + outRoot + ";\n";
		s += "InfiniteLoopKeptStates=" + ttos(infiniteLoops_) + ";\n";
		s += "FiniteLoops=" + finiteLoopsGs_ + ";\n";
		s += "TargetElectronsUp=" + ttos(nup) + ";\n";
		s += "TargetElectronsDown=" + ttos(ndown) + ";\n";
		s += "dir0:Connectors="
		    + buildConnectorsStrWithSecondBath(hoppings, secondBathConnectors) + ";\n";
		s += "potentialV=" + buildPotentialVStr(potV) + ";\n";
		return s;
	}

	// opChar: "'" for the particle branch (c'[0]), "" for the hole branch (c[0]).
	// sourceTv: RestartSourceTvForPsi to select which TV in restartRoot's
	//   checkpoint becomes this run's |gs> reference; -1 (default, omit the
	//   key) is correct when restartRoot has a single natural state to
	//   restore (the plain GS run, or any *_init checkpoint). Non-negative
	//   is needed when birthing a NEW column from an ALREADY-ADVANCED
	//   column's multi-TV checkpoint (see birthColumn/computeFullGrid).
	// complexSource: must match whether restartRoot's own checkpoint was
	//   written with "usecomplex" (buildStepInput always sets it) or not
	//   (the plain GS run / this function's own output never does) --
	//   Checkpoint validation rejects a complex/real mismatch ("Previous
	//   run was complex and this one is not"). false (default) is correct
	//   for column 0's original birth from the plain GS run; true is
	//   needed when birthing a later column from an already-advanced
	//   (hence complex) column's checkpoint.
	std::string buildInitInputAt(const std::string&       outRoot,
	                             const std::string&       restartRoot,
	                             RealType                 U_f,
	                             const VectorRealType&    hoppings,
	                             const VectorRealType&    potV,
	                             SizeType                 nup,
	                             SizeType                 ndown,
	                             SizeType                 nsites,
	                             const std::string&       opChar,
	                             int                      sourceTv      = -1,
	                             bool                     complexSource = false,
	                             const VectorComplexType& secondBathConnectors
	                             = VectorComplexType()) const
	{
		// This run's own OUTPUT needs usecomplex whenever ITS Hamiltonian is
		// complex, not only when the restart SOURCE was (complexSource) --
		// in practice these always coincide for how birthColumn is actually
		// invoked (Vplus(0,p)=0 is real, and by the time Vplus is complex
		// the birth source is already an advanced, complex column), but
		// checking directly rather than relying on that coincidence.
		const bool hasComplexSecondBath = hasNonzeroImag(secondBathConnectors);

		std::string s = "##Ainur1.0\n\n";
		s += geomHeader(nsites, U_f);
		s += "SolverOptions=twositedmrg,geometryallinsystem,TargetingExpression,restart";
		s += (complexSource || hasComplexSecondBath) ? ",usecomplex;\n" : ";\n";
		s += "Version=neqTdmrg;\n";
		s += "OutputFile=" + outRoot + ";\n";
		s += "InfiniteLoopKeptStates=" + ttos(infiniteLoops_) + ";\n";
		s += "FiniteLoops=" + enforceFlag2(finiteLoopsGs_) + ";\n";
		s += "TargetElectronsUp=" + ttos(nup) + ";\n";
		s += "TargetElectronsDown=" + ttos(ndown) + ";\n";
		s += "dir0:Connectors="
		    + buildConnectorsStrWithSecondBath(hoppings, secondBathConnectors) + ";\n";
		s += "potentialV=" + buildPotentialVStr(potV) + ";\n";
		s += "RestartFilename=" + restartRoot + ";\n";
		if (sourceTv >= 0)
			s += "RestartSourceTvForPsi=" + ttos(sourceTv) + ";\n";
		s += "GsWeight=0.1;\n";
		s += "string P0=\"c" + opChar + "[0]*|gs>\";\n";
		return s;
	}

	// One single time-advance segment, restarting from restartRoot.
	// mappedP0Tv: old-run TV index to seed this run's P0 from (RestartMappingTvs).
	// sourceTvForPsi: old-run TV index to seed this run's |gs> reference from
	//   (RestartSourceTvForPsi); -1 (omit the key) means "use whatever |gs>
	//   already is in restartRoot" -- correct only for the very first segment
	//   (n=1), which restarts from the *_init checkpoint's untouched |gs>.
	// advanceEachOverride: 0 (default) means "use tspAdvanceEach_ as read
	// from the input file", correct for the NeqBathRank=0 geometry it was
	// calibrated against. A nonzero value is required once the geometry is
	// extended with 2L second-bath sites (nsitesExt = nsites+2L) -- the
	// N-2 convention must be recomputed from nsitesExt, NOT left at the
	// original (smaller-lattice) value. Confirmed empirically
	// (build/tmp/advanceeach_check/, L=1): the stale value still gave
	// exactly one time-advance per segment for L=1 by coincidence/slack,
	// but this is not guaranteed for larger L -- always pass the
	// recalibrated nsitesExt-2 explicitly when a second bath is present.
	std::string buildStepInput(RealType                 U_f,
	                           const VectorRealType&    hoppings,
	                           const VectorRealType&    potV,
	                           SizeType                 nup,
	                           SizeType                 ndown,
	                           SizeType                 nsites,
	                           const std::string&       restartRoot,
	                           int                      mappedP0Tv,
	                           int                      sourceTvForPsi,
	                           const std::string&       outRoot,
	                           const VectorComplexType& secondBathConnectors
	                           = VectorComplexType(),
	                           SizeType advanceEachOverride = 0) const
	{
		const SizeType advanceEach
		    = (advanceEachOverride == 0) ? tspAdvanceEach_ : advanceEachOverride;

		std::string s = "##Ainur1.0\n\n";
		s += geomHeader(nsites, U_f);
		s += "SolverOptions=twositedmrg,geometryallinsystem,TargetingExpression,restart,"
		     "usecomplex;\n";
		s += "Version=neqTdmrg;\n";
		s += "OutputFile=" + outRoot + ";\n";
		s += "InfiniteLoopKeptStates=" + ttos(infiniteLoops_) + ";\n";
		// Two rows (there-and-back sweep), not one: a single row's sweep is
		// not guaranteed to reach site 0, where the in-situ measurement this
		// segment exists to produce actually lives (confirmed empirically --
		// see project_tdmrg_evolving_bath memory, "Phase 1 narrow-slice
		// progress" -- a single row left site-0 measurements missing
		// entirely from the log starting at the second chained segment).
		s += "FiniteLoops=[[@auto, " + ttos(infiniteLoops_) + ", 2],[@auto, "
		    + ttos(infiniteLoops_) + ", 2]];\n";
		s += "TargetElectronsUp=" + ttos(nup) + ";\n";
		s += "TargetElectronsDown=" + ttos(ndown) + ";\n";
		s += "dir0:Connectors="
		    + buildConnectorsStrWithSecondBath(hoppings, secondBathConnectors) + ";\n";
		s += "potentialV=" + buildPotentialVStr(potV) + ";\n";
		s += "RestartFilename=" + restartRoot + ";\n";
		s += "RestartMappingTvs=[" + ttos(mappedP0Tv) + ", -1, -1];\n";
		if (sourceTvForPsi >= 0)
			s += "RestartSourceTvForPsi=" + ttos(sourceTvForPsi) + ";\n";
		s += "GsWeight=0.1;\n";
		s += "string P0=|P0>;\n";
		s += "string P1=\"TimeEvolve{tau=" + ttos(params_.dt) + ",steps="
		    + ttos(tspTimeSteps_) + ",advanceEach=" + ttos(advanceEach) + "}*|P0>\";\n";
		s += "string P2=\"TimeEvolve{tau=" + ttos(params_.dt) + ",steps="
		    + ttos(tspTimeSteps_) + ",advanceEach=" + ttos(advanceEach) + "}*|gs>\";\n";
		return s;
	}

	// ---- Log parsing -------------------------------------------------------

	// Parse particle tDMRG log for <P2|c|P1> and <P2.last|P2>.
	// Applies gauge correction: ggt0 /= gauge_overlap.
	// Measurement format (site 0, from in_situ_measurements):
	//   "<site> (<re>,<im>) <time> <label>"
	void parseTdmrgLog(const std::string& logfile, std::map<int, ComplexType>& ggt0_at_step)
	{
		std::ifstream fin(logfile);
		if (!fin || !fin.good()) {
			err("ImpuritySolverNeqTdmrg: cannot open tDMRG log '" + logfile + "'\n");
			return;
		}

		std::map<int, ComplexType> gauge_at_step;

		std::string line;
		while (std::getline(fin, line)) {
			SizeType    site = 0;
			std::string valStr, label;
			RealType    t = 0;
			if (!parseMeasurementLine(line, site, valStr, t, label))
				continue;
			if (site != 0)
				continue;

			RealType re = 0, im = 0;
			if (!parseComplex(valStr, re, im))
				continue;

			const int n = static_cast<int>(std::round(t / params_.dt));
			if (n < 0 || n > static_cast<int>(params_.nT))
				continue;

			if (label == "<P2|c|P1>")
				ggt0_at_step[n] = ComplexType(re, im);
			else if (label == "<P2.last|P2>")
				gauge_at_step[n] = ComplexType(re, im);
		}

		if (ggt0_at_step.empty())
			std::cerr << "ImpuritySolverNeqTdmrg: WARNING: no G^> measurements in '"
			          << logfile << "'\n";

		applySignFlip(ggt0_at_step);

		// Divide by gauge: corrects N-sector phase drift between consecutive advances.
		for (auto& kv : ggt0_at_step) {
			auto it = gauge_at_step.find(kv.first);
			if (it != gauge_at_step.end() && std::abs(it->second) > RealType(1e-10))
				kv.second /= it->second;
		}

		applyGlobalPhase(ggt0_at_step);
	}

	// Parse hole tDMRG log for <P1|c|P2> and <P2.last|P2>.
	// Applies gauge correction: glt0 *= gauge_overlap.
	void parseHoleTdmrgLog(const std::string& logfile, std::map<int, ComplexType>& glt0_at_step)
	{
		std::ifstream fin(logfile);
		if (!fin || !fin.good()) {
			std::cerr << "ImpuritySolverNeqTdmrg: WARNING: cannot open hole log '"
			          << logfile << "'\n";
			return;
		}

		std::map<int, ComplexType> gauge_at_step;

		std::string line;
		while (std::getline(fin, line)) {
			SizeType    site = 0;
			std::string valStr, label;
			RealType    t = 0;
			if (!parseMeasurementLine(line, site, valStr, t, label))
				continue;
			if (site != 0)
				continue;

			RealType re = 0, im = 0;
			if (!parseComplex(valStr, re, im))
				continue;

			const int n = static_cast<int>(std::round(t / params_.dt));
			if (n < 0 || n > static_cast<int>(params_.nT))
				continue;

			if (label == "<P1|c|P2>")
				glt0_at_step[n] = ComplexType(re, im);
			else if (label == "<P2.last|P2>")
				gauge_at_step[n] = ComplexType(re, im);
		}

		if (glt0_at_step.empty()) {
			std::cerr << "ImpuritySolverNeqTdmrg: WARNING: no G^< measurements in '"
			          << logfile << "'\n";
			return;
		}

		applySignFlip(glt0_at_step);

		// Multiply by gauge: corrects N-sector phase for hole-run bra.
		for (auto& kv : glt0_at_step) {
			auto it = gauge_at_step.find(kv.first);
			if (it != gauge_at_step.end() && std::abs(it->second) > RealType(1e-10))
				kv.second *= it->second;
		}

		applyGlobalPhase(glt0_at_step);
	}

	void fillKBGrid(const std::map<int, ComplexType>& ggt0_at_step,
	                const std::map<int, ComplexType>& glt0_at_step)
	{
		const int         nT = static_cast<int>(params_.nT);
		const ComplexType pI = ComplexType(0, 1);
		const ComplexType mI = ComplexType(0, -1);

		for (int n = 0; n <= nT; ++n) {
			auto itG = ggt0_at_step.find(n);
			if (itG == ggt0_at_step.end())
				continue;
			const ComplexType ggt = mI * itG->second; // G^> = -i<P2|c|P1>

			auto itL = glt0_at_step.find(n);
			if (itL != glt0_at_step.end()) {
				const ComplexType glt = pI * itL->second; // G^< = +i<P1|c|P2>
				gimp_.lesser(n, 0)    = glt;
				gimp_.retarded(n, 0)  = ggt - glt;
			} else {
				gimp_.retarded(n, 0) = ggt;
			}
		}
	}

	// ---- Utilities ---------------------------------------------------------

	// Detect and correct sign flips from DMRG gauge jumps at sweep reversals.
	// A true sign flip satisfies |M(t)+M(t-dt)|² ≪ |M(t)-M(t-dt)|².
	static void applySignFlip(std::map<int, ComplexType>& m)
	{
		if (m.size() < 2)
			return;
		auto it   = m.begin();
		auto prev = it;
		++it;
		for (; it != m.end(); ++it, ++prev) {
			const ComplexType& a         = prev->second;
			const RealType     sum_norm  = std::norm(a + it->second);
			const RealType     diff_norm = std::norm(a - it->second);
			if (sum_norm < RealType(0.1) * diff_norm)
				it->second = -it->second;
		}
	}

	// Remove residual global phase: rotate so that the t=0 value is real.
	static void applyGlobalPhase(std::map<int, ComplexType>& m)
	{
		auto it0 = m.find(0);
		if (it0 == m.end())
			return;
		const ComplexType z0 = it0->second;
		if (std::abs(z0) < RealType(1e-10))
			return;
		const ComplexType phase_inv = ComplexType(std::abs(z0)) / z0;
		for (auto& kv : m)
			kv.second *= phase_inv;
	}

	// In-situ measurement log line: "<site> (<re>,<im>) <time> <label>"
	static bool parseMeasurementLine(const std::string& line,
	                                 SizeType&          site,
	                                 std::string&       valStr,
	                                 RealType&          t,
	                                 std::string&       label)
	{
		std::istringstream iss(line);
		if (!(iss >> site))
			return false;
		if (!(iss >> valStr))
			return false;
		if (!(iss >> t))
			return false;
		if (!(iss >> label))
			return false;
		return true;
	}

	static bool parseComplex(const std::string& s, RealType& re, RealType& im)
	{
		if (s.empty())
			return false;
		if (s[0] == '(') {
			const std::string body  = s.substr(1, s.size() > 2 ? s.size() - 2 : 0);
			auto              comma = body.find(',');
			if (comma == std::string::npos)
				return false;
			try {
				re = std::stod(body.substr(0, comma));
				im = std::stod(body.substr(comma + 1));
				return true;
			} catch (...) {
				return false;
			}
		}
		try {
			re = std::stod(s);
			im = RealType(0);
			return true;
		} catch (...) {
			return false;
		}
	}

	// Phase 1 diagnostic: a single-step segment's own FiniteLoops has TWO
	// rows (see buildStepInput), so its log has TWO site-0 measurement
	// lines for a given label. Which one is the wanted post-advance value
	// depends on how this segment was restarted (confirmed empirically by
	// diffing against ImpuritySolverNeqTdmrg::solve()'s monolithic result;
	// see project_tdmrg_evolving_bath memory):
	//   - restarting straight from the *_init checkpoint (accumulated time
	//     0): the FIRST match is the trivial t=0 value, the SECOND is the
	//     wanted one -- takeLast=true.
	//   - restarting via RestartSourceTvForPsi from a PREVIOUS chained step
	//     (which already carries accumulated time forward): the FIRST
	//     match is the wanted value, the SECOND is one advance too far --
	//     takeLast=false.
	static bool parseSingleMeasurement(const std::string& logfile,
	                                   const std::string& label,
	                                   ComplexType&       outVal,
	                                   bool               takeLast)
	{
		return parseSingleMeasurementAtSite(logfile, label, 0, outVal, takeLast);
	}

	// Generalizes parseSingleMeasurement's hardcoded site==0 filter (valid
	// only for measurements anchored at the impurity). Needed for
	// measureSecondBathOccupations, whose <gs|n[p]|gs> labels are anchored
	// at the second-bath sites (site index p >= nsites), not site 0.
	static bool parseSingleMeasurementAtSite(const std::string& logfile,
	                                         const std::string& label,
	                                         SizeType           targetSite,
	                                         ComplexType&       outVal,
	                                         bool               takeLast)
	{
		std::ifstream fin(logfile);
		if (!fin || !fin.good())
			return false;

		bool        found = false;
		std::string line;
		while (std::getline(fin, line)) {
			SizeType    site = 0;
			std::string valStr, lbl;
			RealType    t = 0;
			if (!parseMeasurementLine(line, site, valStr, t, lbl))
				continue;
			if (site != targetSite || lbl != label)
				continue;

			RealType re = 0, im = 0;
			if (!parseComplex(valStr, re, im))
				continue;

			outVal = ComplexType(re, im);
			found  = true;
			if (!takeLast)
				return true;
		}
		return found;
	}

	// Replace FiniteLoops flag 0 with flag 2 to prevent |gs⟩ re-optimisation.
	// The "2" (fast-WFT) flag keeps |gs⟩ = WFT-transformed |GS_i⟩, which is
	// required so that P2 = TimeEvolve*|gs⟩ evolves the pre-quench ground state.
	static std::string enforceFlag2(const std::string& loops)
	{
		std::string            result = loops;
		std::string::size_type pos    = 0;
		while ((pos = result.find(", 0]", pos)) != std::string::npos) {
			result.replace(pos, 4, ", 2]");
			pos += 4;
		}
		return result;
	}

	// Common geometry/model header block for all five runs.
	std::string geomHeader(SizeType nsites, RealType U) const
	{
		std::string s;
		s += "TotalNumberOfSites=" + ttos(nsites) + ";\n";
		s += "NumberOfTerms=1;\n";
		s += "DegreesOfFreedom=1;\n";
		s += "GeometryKind=star;\n";
		s += "GeometryOptions=none;\n";
		s += "hubbardU=" + buildHubbardUStr(U, nsites) + ";\n";
		s += "Model=HubbardOneBand;\n";
		return s;
	}

	static std::string buildHubbardUStr(RealType U, SizeType nsites)
	{
		std::string s = "[" + ttos(U);
		for (SizeType i = 1; i < nsites; ++i)
			s += ", 0.";
		return s + "]";
	}

	static std::string buildConnectorsStr(const VectorRealType& hoppings)
	{
		std::string s = "[";
		for (SizeType i = 0; i < hoppings.size(); ++i) {
			if (i > 0)
				s += ",";
			s += ttos(hoppings[i]);
		}
		return s + "]";
	}

	// Ainur's complex literal syntax is "<real>i<imag>" (e.g. "0.5i-0.3" ->
	// 0.5-0.3i), NOT "(re,im)" or "a+bi" -- confirmed against the live
	// AinurComplex::toComplex parser (see fancy-painting-moon.md Phase 2
	// blocker-A). A zero imaginary part is written as a plain real literal
	// (no "i" token at all), matching what buildConnectorsStr already emits
	// for the (always-real) first bath.
	// True if any entry has nonzero imaginary part -- used to decide whether
	// a run's SolverOptions needs "usecomplex". Shared by buildGsInputAt and
	// buildInitInputAt; buildGsInputAt originally omitted this check
	// entirely (a real bug: it fed a complex-literal-bearing Connectors
	// string to a non-"usecomplex" run, which fails Ainur parsing with
	// e.g. "atof received a non-digit 0.7i0.2" -- found via the diagonal
	// Connectors-independence test, fancy-painting-moon.md Phase 2).
	static bool hasNonzeroImag(const VectorComplexType& v)
	{
		for (SizeType p = 0; p < v.size(); ++p)
			if (std::imag(v[p]) != RealType(0))
				return true;
		return false;
	}

	static std::string formatComplexLiteral(ComplexType v)
	{
		const RealType re = std::real(v);
		const RealType im = std::imag(v);
		if (im == RealType(0))
			return ttos(re);
		return ttos(re) + "i" + ttos(im);
	}

	// Second-bath-aware Connectors: first-bath hoppings (always real, as
	// today) followed by the second-bath couplings Vplus(n,p), which are
	// genuinely complex in general (Cholesky factor of Lambda^gtrless --
	// see NeqBathDecomposition) and so need the complex literal syntax
	// above. Used only when neqBathRank_ > 0; the NeqBathRank=0 path
	// keeps using the plain buildConnectorsStr unchanged.
	static std::string
	buildConnectorsStrWithSecondBath(const VectorRealType&    firstBathHoppings,
	                                 const VectorComplexType& secondBathHoppings)
	{
		std::string s     = "[";
		bool        first = true;
		for (SizeType i = 0; i < firstBathHoppings.size(); ++i) {
			if (!first)
				s += ",";
			s += ttos(firstBathHoppings[i]);
			first = false;
		}
		for (SizeType i = 0; i < secondBathHoppings.size(); ++i) {
			if (!first)
				s += ",";
			s += formatComplexLiteral(secondBathHoppings[i]);
			first = false;
		}
		return s + "]";
	}

	static std::string buildPotentialVStr(const VectorRealType& potV)
	{
		std::string inner;
		for (SizeType i = 0; i < potV.size(); ++i) {
			if (i > 0)
				inner += ",";
			inner += ttos(potV[i]);
		}
		return "[" + inner + "," + inner + "]";
	}

	// ---- Member data -------------------------------------------------------
	const ParamsNeqType&            params_;
	const ApplicationType&          app_;
	typename InputNgType::Readable& io_;
	const SizeType                  neqBathRank_;
	ExactDiagType                   exactDiag_;
	KBType                          gimp_;
	SizeType                        nup_   = 0;
	SizeType                        ndown_ = 0;
	std::string                     root_;
	SizeType                        infiniteLoops_ = 0;
	std::string                     finiteLoopsGs_;
	std::string                     finiteLoopsTdmrg_;
	SizeType                        tspTimeSteps_   = 5;
	SizeType                        tspAdvanceEach_ = 1;

	// ---- Phase 2 (evolving-bath project): self-consistency state ----------
	// Populated by solve() only when neqBathRank_ > 0; unused (default-
	// constructed/empty) on the NeqBathRank=0 path, which stays on the
	// original code above unchanged.
	std::unique_ptr<DecompType> decomp_;
	VectorRealType              scHoppings_, scBathEps_, scPotTdmrg_;
	SizeType                    scNsitesExt_ = 0;
	SizeType                    scL_         = 0;
	SizeType                    scNup_ = 0, scNdown_ = 0;
	RealType                    scEps_ = 0;
	std::string                 scChainRoot_;
};

} // namespace Dmft
#endif // IMPURITYSOLVER_NEQ_TDMRG_H
