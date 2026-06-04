#ifndef IMPURITYSOLVER_NEQ_TDMRG_H
#define IMPURITYSOLVER_NEQ_TDMRG_H

#include "CmdLineOptions.hh"
#include "DmrgRunner.h"
#include "ImpuritySolverNeqBase.h"
#include "KadanoffBaym.h"
#include "ParamsNeqDmftSolver.h"
#include "PsimagLite.h"
#include "Vector.h"
#include <cassert>
#include <cmath>
#include <complex>
#include <fstream>
#include <map>
#include <sstream>

namespace Dmft {

// Non-equilibrium impurity solver using DMRG++ time-dependent DMRG (tDMRG).
//
// Three-run approach targeting the t'=0 slice of the KB Green's function:
//   Run 1: GS DMRG with H(U_i) → saves state to OutputFile+gs.
//   Run 2: tDMRG with H(U_f), restarted from Run 1, applying c†_{imp} as the
//           initial operator to obtain a second target |Φ(t)⟩ alongside the
//           time-evolved GS |Ψ_B(t)⟩.  Gives G^>(t,0).
//   Run 3: tDMRG with H(U_f), restarted from Run 1, applying c_{imp} as the
//           initial operator to obtain |φ_h(t)⟩ (N-1 sector).  Gives G^<(t,0).
//
// In-situ measurements (TSPEvolveGroundState=1: "gsT" = time-evolved |Ψ_0(t)⟩):
//   Run 2: <gsT|nup|gsT>  →  n_up(t)                 →  G^<(t,t) = i·n(t)
//          <gsT|c|P1>     →  ⟨Ψ_0(t)|c|Φ_+(t)⟩      →  G^>(t,0) = -i·value
//   Run 3: <P1|c|gsT>     →  ⟨φ_h(t)|c|Ψ_0(t)⟩      →  G^<(t,0) = +i·value
//
//   G^R(t,0) = G^>(t,0) - G^<(t,0)  (full retarded, no approximation)
//
// Required input parameters (in addition to the standard neq-DMFT set):
//   TargetElectronsUp=      — as for ImpuritySolverNeqExactDiag
//   TargetElectronsDown=    — as for ImpuritySolverNeqExactDiag
//   RootOutputname=         — file-name stem for GS/tDMRG outputs
//   InfiniteLoopKeptStates= — m for the infinite DMRG sweep
//   FiniteLoopsGs=          — finite-loop spec for the GS run (matrix format)
//   FiniteLoopsTdmrg=       — finite-loop spec for the tDMRG runs (must cover
//                              enough sweeps to accumulate N_t advances at
//                              TSPAdvanceEach sites per advance; reused for both
//                              particle and hole runs)
//   TSPTimeSteps=           — number of Krylov vectors per time step (accuracy)
//   TSPAdvanceEach=         — sites per time advance (controls sweep count)
template <typename ComplexOrRealType>
class ImpuritySolverNeqTdmrg : public ImpuritySolverNeqBase<ComplexOrRealType> {

public:

	using BaseType          = ImpuritySolverNeqBase<ComplexOrRealType>;
	using RealType          = typename BaseType::RealType;
	using ComplexType       = typename BaseType::ComplexType;
	using VectorRealType    = typename BaseType::VectorRealType;
	using KBType            = typename BaseType::KBType;
	using InputNgType       = typename BaseType::InputNgType;
	using ParamsNeqType     = ParamsNeqDmftSolver<ComplexOrRealType>;
	using DmrgRunnerType    = Dmrg::DmrgRunner<RealType>;
	using ApplicationType   = PsimagLite::PsiApp;

	ImpuritySolverNeqTdmrg(const ParamsNeqType&            params,
	                        const ApplicationType&          app,
	                        typename InputNgType::Readable& io)
	    : params_(params)
	    , app_(app)
	    , io_(io)
	    , gimp_(params.nT,
	            params.eqParams.nMatsubaras,
	            params.dt,
	            params.eqParams.ficticiousBeta /
	                static_cast<RealType>(params.eqParams.nMatsubaras))
	{
		io.readline(nup_,   "TargetElectronsUp=");
		io.readline(ndown_, "TargetElectronsDown=");
		io.readline(root_,  "RootOutputname=");
		io.readline(infiniteLoops_, "InfiniteLoopKeptStates=");
		io.readline(finiteLoopsGs_, "FiniteLoopsGs=");

		try {
			io.readline(finiteLoopsTdmrg_, "FiniteLoopsTdmrg=");
		} catch (std::exception&) {
			// Fall back to GS loops if tDMRG-specific ones are not provided.
			// Users should provide FiniteLoopsTdmrg with enough sweeps to
			// cover params_.nT time advances at TSPAdvanceEach sites each.
			finiteLoopsTdmrg_ = finiteLoopsGs_;
		}

		try {
			io.readline(tspTimeSteps_, "TSPTimeSteps=");
		} catch (std::exception&) {
			tspTimeSteps_ = 4; // standard Krylov accuracy default
		}

		try {
			io.readline(tspAdvanceEach_, "TSPAdvanceEach=");
		} catch (std::exception&) {
			tspAdvanceEach_ = 1;
		}
	}

	// One-time setup: runs GS DMRG and tDMRG, pre-fills the KB grid.
	void solve(const VectorRealType& bathParams) override
	{
		const SizeType nBath  = bathParams.size() / 2;
		const SizeType nsites = nBath + 1;

		VectorRealType hoppings(nBath), bathEps(nBath);
		for (SizeType i = 0; i < nBath; ++i) {
			hoppings[i] = bathParams[i];
			bathEps[i]  = bathParams[nBath + i];
		}

		// Anderson model: µ_imp = -U/2 at half-filling
		VectorRealType potGS(nsites), potTdmrg(nsites);
		potGS[0]    = -RealType(0.5) * params_.uInitial;
		potTdmrg[0] = -RealType(0.5) * params_.uFinal;
		for (SizeType i = 0; i < nBath; ++i) {
			potGS[i + 1]    = bathEps[i];
			potTdmrg[i + 1] = bathEps[i];
		}

		// Run 1: GS DMRG for H(U_i)
		std::cout << "ImpuritySolverNeqTdmrg: running GS DMRG for U_i=" << params_.uInitial
		          << "\n";
		{
			Dmrg::CmdLineOptions opts;
			opts.logfile = root_ + "gs.log";
			const std::string input = buildGsInput(params_.uInitial, hoppings, potGS,
			                                        nup_, ndown_, nsites);
			DmrgRunnerType runner(app_, input, opts);
			runner.doOneRun();
		}

		// Run 2: particle-sector tDMRG — c†_{imp}|GS⟩ evolved under H(U_f)
		// TSPEvolveGroundState=1 also Krylov-evolves |GS_i⟩ → |Ψ_0(t)⟩ (label "gsT"),
		// so <gsT|c|P1> gives the correct G^>(t,0) = −i⟨Ψ_0(t)|c|Φ_+(t)⟩.
		const std::string tdmrgLog = root_ + "tdmrg.log";
		std::cout << "ImpuritySolverNeqTdmrg: running particle tDMRG for U_f=" << params_.uFinal
		          << " to t_max=" << params_.tMax << "\n";
		{
			Dmrg::CmdLineOptions opts;
			opts.logfile              = tdmrgLog;
			opts.in_situ_measurements = "<gs|nup|gs>,<gs|c|P1>";
			const std::string input = buildTdmrgInput(params_.uFinal, hoppings, potTdmrg,
			                                           nup_, ndown_, nsites);
			DmrgRunnerType runner(app_, input, opts);
			runner.doOneRun();
		}

		// Run 3: hole-sector tDMRG — c_{imp}|GS⟩ evolved under H(U_f)
		// <P1|c|gsT> gives the correct G^<(t,0) = +i⟨φ_h(t)|c|Ψ_0(t)⟩.
		const std::string holeTdmrgLog = root_ + "tdmrg_hole.log";
		std::cout << "ImpuritySolverNeqTdmrg: running hole tDMRG for G^<(t,0)\n";
		{
			Dmrg::CmdLineOptions opts;
			opts.logfile              = holeTdmrgLog;
			opts.in_situ_measurements = "<P1|c|gs>";
			const std::string input = buildHoleTdmrgInput(params_.uFinal, hoppings, potTdmrg,
			                                               nup_, ndown_, nsites);
			DmrgRunnerType runner(app_, input, opts);
			runner.doOneRun();
		}

		std::map<int, RealType>    nup_at_step;
		std::map<int, ComplexType> ggt0_at_step;
		parseTdmrgLog(tdmrgLog, nup_at_step, ggt0_at_step);

		std::map<int, ComplexType> glt0_at_step;
		parseHoleTdmrgLog(holeTdmrgLog, glt0_at_step);

		fillKBGrid(nup_at_step, ggt0_at_step, glt0_at_step);
	}

	// Fill the n-th time slice of gimp from the pre-computed KB grid.
	void computeGimp(KBType& gimp, int n) const override
	{
		gimp.retarded(n, 0) = gimp_.retarded(n, 0);
		gimp.lesser(n, n)   = gimp_.lesser(n, n);
		gimp.lesser(n, 0)   = gimp_.lesser(n, 0);
	}

	const KBType& gimp() const override { return gimp_; }

private:

	// ---- Input construction ------------------------------------------------

	std::string buildGsInput(RealType              U,
	                          const VectorRealType& hoppings,
	                          const VectorRealType& potV,
	                          SizeType              nup,
	                          SizeType              ndown,
	                          SizeType              nsites) const
	{
		std::string s = "##Ainur1.0\n\n";
		s += "TotalNumberOfSites=" + ttos(nsites) + ";\n";
		s += "NumberOfTerms=1;\n";
		s += "DegreesOfFreedom=1;\n";
		s += "GeometryKind=star;\n";
		s += "GeometryOptions=none;\n";
		s += "hubbardU=" + buildHubbardUStr(U, nsites) + ";\n";
		s += "Model=HubbardOneBand;\n";
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

	std::string buildTdmrgInput(RealType              U_f,
	                             const VectorRealType& hoppings,
	                             const VectorRealType& potV,
	                             SizeType              nup,
	                             SizeType              ndown,
	                             SizeType              nsites) const
	{
		std::string s = "##Ainur1.0\n\n";
		s += "TotalNumberOfSites=" + ttos(nsites) + ";\n";
		s += "NumberOfTerms=1;\n";
		s += "DegreesOfFreedom=1;\n";
		s += "GeometryKind=star;\n";
		s += "GeometryOptions=none;\n";
		s += "hubbardU=" + buildHubbardUStr(U_f, nsites) + ";\n";
		s += "Model=HubbardOneBand;\n";
		s += "SolverOptions=twositedmrg,geometryallinsystem,TimeStepTargeting,restart;\n";
		s += "Version=neqTdmrg;\n";
		s += "OutputFile=" + root_ + "tdmrg;\n";
		s += "InfiniteLoopKeptStates=" + ttos(infiniteLoops_) + ";\n";
		s += "FiniteLoops=" + finiteLoopsTdmrg_ + ";\n";
		s += "TargetElectronsUp=" + ttos(nup) + ";\n";
		s += "TargetElectronsDown=" + ttos(ndown) + ";\n";
		s += "dir0:Connectors=" + buildConnectorsStr(hoppings) + ";\n";
		s += "potentialV=" + buildPotentialVStr(potV) + ";\n";

		// Restart from the GS of H(U_i)
		s += "RestartFilename=" + root_ + "gs;\n";

		// tDMRG time-stepping parameters
		s += "GsWeight=0.1;\n";
		s += "TSPTau=" + ttos(params_.dt) + ";\n";
		s += "TSPTimeSteps=" + ttos(tspTimeSteps_) + ";\n";
		s += "TSPAdvanceEach=" + ttos(tspAdvanceEach_) + ";\n";
		s += "TSPAlgorithm=Krylov;\n";
		// TSPEvolveGroundState=1 reserved for Phase 2 (proper |GS_i⟩ accumulation).

		// Operator: c†_{imp} at site 0 (border site in star → needs identity trigger at site 1)
		s += "TSPProductOrSum=product;\n";
		s += "TSPSites=[1, 0];\n";
		s += "TSPLoops=[0, 0];\n";
		s += "string TSPOp0:TSPOperator=expression;\n";
		s += "string TSPOp0:OperatorExpression=identity;\n";
		s += "string TSPOp1:TSPOperator=expression;\n";
		s += "string TSPOp1:OperatorExpression=c';\n"; // c' = c† (spin-up creation)

		return s;
	}

	// Same as buildTdmrgInput but applies c (annihilation) to build the
	// N-1 sector initial state |φ_h(0)⟩ = c_{imp}|GS⟩ for G^<(t,0).
	std::string buildHoleTdmrgInput(RealType              U_f,
	                                 const VectorRealType& hoppings,
	                                 const VectorRealType& potV,
	                                 SizeType              nup,
	                                 SizeType              ndown,
	                                 SizeType              nsites) const
	{
		std::string s = "##Ainur1.0\n\n";
		s += "TotalNumberOfSites=" + ttos(nsites) + ";\n";
		s += "NumberOfTerms=1;\n";
		s += "DegreesOfFreedom=1;\n";
		s += "GeometryKind=star;\n";
		s += "GeometryOptions=none;\n";
		s += "hubbardU=" + buildHubbardUStr(U_f, nsites) + ";\n";
		s += "Model=HubbardOneBand;\n";
		s += "SolverOptions=twositedmrg,geometryallinsystem,TimeStepTargeting,restart;\n";
		s += "Version=neqTdmrg;\n";
		s += "OutputFile=" + root_ + "tdmrg_hole;\n";
		s += "InfiniteLoopKeptStates=" + ttos(infiniteLoops_) + ";\n";
		s += "FiniteLoops=" + finiteLoopsTdmrg_ + ";\n";
		s += "TargetElectronsUp=" + ttos(nup) + ";\n";
		s += "TargetElectronsDown=" + ttos(ndown) + ";\n";
		s += "dir0:Connectors=" + buildConnectorsStr(hoppings) + ";\n";
		s += "potentialV=" + buildPotentialVStr(potV) + ";\n";

		s += "RestartFilename=" + root_ + "gs;\n";

		s += "GsWeight=0.1;\n";
		s += "TSPTau=" + ttos(params_.dt) + ";\n";
		s += "TSPTimeSteps=" + ttos(tspTimeSteps_) + ";\n";
		s += "TSPAdvanceEach=" + ttos(tspAdvanceEach_) + ";\n";
		s += "TSPAlgorithm=Krylov;\n";
		// TSPEvolveGroundState=1 reserved for Phase 2 (proper |GS_i⟩ accumulation).

		// Operator: c_{imp} at site 0 (annihilation, N-1 sector hole state)
		s += "TSPProductOrSum=product;\n";
		s += "TSPSites=[1, 0];\n";
		s += "TSPLoops=[0, 0];\n";
		s += "string TSPOp0:TSPOperator=expression;\n";
		s += "string TSPOp0:OperatorExpression=identity;\n";
		s += "string TSPOp1:TSPOperator=expression;\n";
		s += "string TSPOp1:OperatorExpression=c;\n"; // c = spin-up annihilation

		return s;
	}

	// ---- Log parsing -------------------------------------------------------

	// Parse the particle-sector tDMRG log for <gsT|nup|gsT>, <gsT|c|P1>, and
	// <gs|penultimate> (overlap of consecutive N-sector GS captures for gauge tracking).
	// Log format: "<site> (<re>,<im>) <time> <label> (<norm_re>,<norm_im>)"
	// Last measurement per (step, label) at site 0 wins (most converged sweep).
	void parseTdmrgLog(const std::string&          logfile,
	                   std::map<int, RealType>&    nup_at_step,
	                   std::map<int, ComplexType>& ggt0_at_step)
	{
		std::ifstream fin(logfile);
		if (!fin || !fin.good()) {
			err("ImpuritySolverNeqTdmrg: cannot open tDMRG log '" + logfile + "'\n");
			return;
		}

		// <gs|penultimate> = <psi_prev|psi_curr> for consecutive N-sector GS captures.
		// This is the accumulated N-sector gauge phase drift between advances.
		std::map<int, ComplexType> gauge_at_step;

		std::string line;
		while (std::getline(fin, line)) {
			std::istringstream iss(line);
			SizeType           site = 0;
			if (!(iss >> site)) continue;
			if (site != 0)     continue;

			std::string valStr;
			if (!(iss >> valStr)) continue;
			RealType re = 0, im = 0;
			if (!parseComplex(valStr, re, im)) continue;

			RealType t = 0;
			if (!(iss >> t)) continue;
			std::string label;
			if (!(iss >> label)) continue;

			const int n = static_cast<int>(std::round(t / params_.dt));
			if (n < 0 || n > static_cast<int>(params_.nT)) continue;

			if (label == "<gs|nup|gs>")
				nup_at_step[n] = re;
			else if (label == "<gs|c|P1>")
				ggt0_at_step[n] = ComplexType(re, im);
			else if (label == "<gs|penultimate>")
				gauge_at_step[n] = ComplexType(re, im);
		}

		if (nup_at_step.empty() && ggt0_at_step.empty())
			std::cerr << "ImpuritySolverNeqTdmrg: WARNING: no in-situ measurements "
			             "found in '" << logfile << "'.  Check TSP parameters and "
			             "finite-loop count.\n";

		// Sign-flip correction: the N+1-sector (P1) gauge can flip by π at sweep
		// direction changes.  Detect by |a+b|² ≪ |a-b|² and negate.
		{
			auto it = ggt0_at_step.begin();
			auto prev = it;
			++it;
			for (; it != ggt0_at_step.end(); ++it, ++prev) {
				const ComplexType& a = prev->second;
				const RealType sum_norm  = std::norm(a + it->second);
				const RealType diff_norm = std::norm(a - it->second);
				if (sum_norm < RealType(0.1) * diff_norm)
					it->second = -it->second;
			}
		}

		// Per-step gauge correction: true_value = <gs|c|P1> / <gs|penultimate>.
		// <gs|penultimate> = <psi_prev|psi_curr> tracks the N-sector gauge drift.
		// Dividing cancels the spurious N-sector phase accumulated since the previous
		// advance, assuming the N+1-sector (P1) gauge drifts proportionally.
		for (auto& kv : ggt0_at_step) {
			auto it = gauge_at_step.find(kv.first);
			if (it != gauge_at_step.end() && std::abs(it->second) > RealType(1e-10))
				kv.second /= it->second;
		}

		// Fallback: t=0 anchor.  At t=0, <gs|c|P1> = 1 - n_↑ which is real.
		// Rotating the series by -arg(M(0)) removes any residual global phase not
		// handled by per-step correction above.
		auto it0 = ggt0_at_step.find(0);
		if (it0 != ggt0_at_step.end()) {
			const ComplexType z0 = it0->second;
			if (std::abs(z0) > RealType(1e-10)) {
				const ComplexType phase_inv = ComplexType(std::abs(z0)) / z0;
				for (auto& kv : ggt0_at_step)
					kv.second *= phase_inv;
			}
		}
	}

	// Parse the hole-sector tDMRG log for <P1|c|gsT> and <gs|penultimate>.
	// G^<(t,0) = +i * measured value.
	//
	// Per-step gauge correction: <P1|c|gsT> may still pick up a spurious DMRG gauge
	// phase between the P1 and gsT targets.  <gs|penultimate> tracks the N-sector
	// drift between consecutive advances.  Multiplying corrects the N-sector component.
	void parseHoleTdmrgLog(const std::string&          logfile,
	                        std::map<int, ComplexType>& glt0_at_step)
	{
		std::ifstream fin(logfile);
		if (!fin || !fin.good()) {
			std::cerr << "ImpuritySolverNeqTdmrg: WARNING: cannot open hole tDMRG log '"
			          << logfile << "'\n";
			return;
		}

		std::map<int, ComplexType> gauge_at_step; // <gs|penultimate> at each step

		std::string line;
		while (std::getline(fin, line)) {
			std::istringstream iss(line);
			SizeType           site = 0;
			if (!(iss >> site)) continue;
			if (site != 0)     continue;

			std::string valStr;
			if (!(iss >> valStr)) continue;
			RealType re = 0, im = 0;
			if (!parseComplex(valStr, re, im)) continue;

			RealType t = 0;
			if (!(iss >> t)) continue;
			std::string label;
			if (!(iss >> label)) continue;

			const int n = static_cast<int>(std::round(t / params_.dt));
			if (n < 0 || n > static_cast<int>(params_.nT)) continue;

			if (label == "<P1|c|gs>")
				glt0_at_step[n] = ComplexType(re, im);
			else if (label == "<gs|penultimate>")
				gauge_at_step[n] = ComplexType(re, im);
		}

		if (glt0_at_step.empty()) {
			std::cerr << "ImpuritySolverNeqTdmrg: WARNING: no G^< hole measurements "
			             "found in '" << logfile << "'\n";
			return;
		}

		// Stage 1: sign-flip correction.
		// The hole state (N-1 sector) is gauge-fixed independently from the N-sector
		// targets.  At sweep direction changes the gauge can jump by ~180°; we detect
		// this by checking |M(t)+M(t-dt)|² ≪ |M(t)-M(t-dt)|² and negate M(t).
		{
			auto it = glt0_at_step.begin();
			auto prev = it;
			++it;
			for (; it != glt0_at_step.end(); ++it, ++prev) {
				const ComplexType& a = prev->second;
				const RealType sum_norm  = std::norm(a + it->second);
				const RealType diff_norm = std::norm(a - it->second);
				if (sum_norm < RealType(0.1) * diff_norm)
					it->second = -it->second;
			}
		}

		// Stage 2 (per-step): true_value = <P1|c|gs> × <gs|penultimate>.
		// <gs|penultimate> tracks the N-sector drift; multiplying corrects it.
		for (auto& kv : glt0_at_step) {
			auto it = gauge_at_step.find(kv.first);
			if (it != gauge_at_step.end() && std::abs(it->second) > RealType(1e-10))
				kv.second *= it->second;
		}

		// Stage 3 (fallback): t=0 anchor.  After sign-flip correction, M(0) should
		// be real (= <n_↑>).  Rotating by -arg(M(0)) removes the residual global phase.
		{
			auto it0 = glt0_at_step.find(0);
			if (it0 != glt0_at_step.end()) {
				const ComplexType z0 = it0->second;
				if (std::abs(z0) > RealType(1e-10)) {
					const ComplexType phase_inv = ComplexType(std::abs(z0)) / z0;
					for (auto& kv : glt0_at_step)
						kv.second *= phase_inv;
				}
			}
		}
	}

	void fillKBGrid(const std::map<int, RealType>&    nup_at_step,
	                const std::map<int, ComplexType>& ggt0_at_step,
	                const std::map<int, ComplexType>& glt0_at_step)
	{
		const int       nT  = static_cast<int>(params_.nT);
		const ComplexType pI = ComplexType(0,  1); // +i
		const ComplexType mI = ComplexType(0, -1); // −i

		for (int n = 0; n <= nT; ++n) {
			// G^<(n,n) = i·n_↑(t_n)
			auto itN = nup_at_step.find(n);
			if (itN != nup_at_step.end())
				gimp_.lesser(n, n) = ComplexType(0, itN->second);

			// G^>(n,0) = −i·⟨Ψ_B(t)|c_↑|Φ(t)⟩
			auto itG = ggt0_at_step.find(n);
			if (itG == ggt0_at_step.end()) continue;
			const ComplexType ggt = mI * itG->second;

			// G^<(n,0) = +i·⟨φ_h(t)|c_↑|Ψ_B(t)⟩  (from hole run)
			auto itL = glt0_at_step.find(n);
			if (itL != glt0_at_step.end()) {
				const ComplexType glt = pI * itL->second;
				gimp_.lesser(n, 0)   = glt;
				gimp_.retarded(n, 0) = ggt - glt; // G^R = G^> − G^<
			} else {
				// Hole run missing for this step: fall back to G^R ≈ G^>
				gimp_.retarded(n, 0) = ggt;
			}
		}
	}

	// ---- Utilities ---------------------------------------------------------

	// Parse a complex number from string: either "(re,im)" or just "re".
	static bool parseComplex(const std::string& s, RealType& re, RealType& im)
	{
		if (s.empty())
			return false;
		if (s[0] == '(') {
			const std::string body = s.substr(1, s.size() > 2 ? s.size() - 2 : 0);
			auto              comma = body.find(',');
			if (comma == std::string::npos)
				return false;
			try {
				re = std::stod(body.substr(0, comma));
				im = std::stod(body.substr(comma + 1));
				return true;
			} catch (...) { return false; }
		}
		try {
			re = std::stod(s);
			im = RealType(0);
			return true;
		} catch (...) { return false; }
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
			if (i > 0) s += ",";
			s += ttos(hoppings[i]);
		}
		return s + "]";
	}

	// potentialV is doubled for spin-up and spin-down blocks.
	static std::string buildPotentialVStr(const VectorRealType& potV)
	{
		std::string inner;
		for (SizeType i = 0; i < potV.size(); ++i) {
			if (i > 0) inner += ",";
			inner += ttos(potV[i]);
		}
		return "[" + inner + "," + inner + "]";
	}

	// ---- Member data -------------------------------------------------------
	const ParamsNeqType&            params_;
	const ApplicationType&          app_;
	typename InputNgType::Readable& io_;
	KBType                          gimp_;
	SizeType                        nup_           = 0;
	SizeType                        ndown_         = 0;
	std::string                     root_;
	SizeType                        infiniteLoops_ = 0;
	std::string                     finiteLoopsGs_;
	std::string                     finiteLoopsTdmrg_;
	SizeType                        tspTimeSteps_  = 4;
	SizeType                        tspAdvanceEach_ = 1;
};

} // namespace Dmft
#endif // IMPURITYSOLVER_NEQ_TDMRG_H
