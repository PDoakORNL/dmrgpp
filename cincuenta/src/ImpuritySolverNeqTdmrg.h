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
// Two-run approach targeting the t'=0 slice of the KB Green's function:
//   Run 1: GS DMRG with H(U_i) → saves state to OutputFile+gs.
//   Run 2: tDMRG with H(U_f), restarted from Run 1, applying c†_{imp} as the
//           initial operator to obtain a second target |Φ(t)⟩ alongside the
//           time-evolved GS |Ψ_B(t)⟩.
//
// In-situ measurements during Run 2:
//   <gs|nup|gs>  →  n_up(t) at impurity site  →  G^<(t,t) = i·n(t)
//   <gs|c|P1>    →  ⟨Ψ_B(t)|c|Φ(t)⟩          →  G^>(t,0) = -i·⟨gs|c|P1⟩
//   G^R(t,0)     ≈  G^>(t,0)  (G^< off-diagonal not yet computed)
//
// Required input parameters (in addition to the standard neq-DMFT set):
//   TargetElectronsUp=      — as for ImpuritySolverNeqExactDiag
//   TargetElectronsDown=    — as for ImpuritySolverNeqExactDiag
//   RootOutputname=         — file-name stem for GS/tDMRG outputs
//   InfiniteLoopKeptStates= — m for the infinite DMRG sweep
//   FiniteLoopsGs=          — finite-loop spec for the GS run (matrix format)
//   FiniteLoopsTdmrg=       — finite-loop spec for the tDMRG run (must cover
//                              enough sweeps to accumulate N_t advances at
//                              TSPAdvanceEach sites per advance)
//   TSPTimeSteps=           — number of Krylov vectors per time step (accuracy)
//   TSPAdvanceEach=         — sites per time advance (controls sweep count)
//
// Limitation (MVP): G^<(t,0) for t>0 and G^{Left}(t,τ) are not computed here;
// those components are zero in gimp_.  The t'=0 slice is sufficient to compare
// G^R(t,0) and n(t) against ImpuritySolverNeqExactDiag and Tsuji's IPT results.
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
	void initialize(const VectorRealType& bathParams) override
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

		// Run 2: tDMRG for H(U_f), restarted from GS, with c†_{imp} as initial operator
		const std::string tdmrgLog = root_ + "tdmrg.log";
		std::cout << "ImpuritySolverNeqTdmrg: running tDMRG for U_f=" << params_.uFinal
		          << " to t_max=" << params_.tMax << "\n";
		{
			Dmrg::CmdLineOptions opts;
			opts.logfile             = tdmrgLog;
			opts.in_situ_measurements = "<gs|nup|gs>,<gs|c|P1>";
			const std::string input = buildTdmrgInput(params_.uFinal, hoppings, potTdmrg,
			                                           nup_, ndown_, nsites);
			DmrgRunnerType runner(app_, input, opts);
			runner.doOneRun();
		}

		parseTdmrgLog(tdmrgLog);
	}

	// Fill the n-th time slice of gimp from the pre-computed KB grid.
	// Only G^R(n,0) and the diagonal G^<(n,n) are populated; other components
	// that require the hole-sector tDMRG run remain zero for this MVP.
	void computeGimp(KBType& gimp, int n) const override
	{
		gimp.retarded(n, 0) = gimp_.retarded(n, 0);
		gimp.lesser(n, n)   = gimp_.lesser(n, n);
		// Matsubara: not computed by tDMRG; leave as zeros.
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

	// ---- Log parsing -------------------------------------------------------

	// Parse the tDMRG log file written by DmrgRunner to extract time-series
	// in-situ measurements and fill the KB grid.
	//
	// Log format (printed by TargetingCommon::test()):
	//   "<site> (<re>,<im>) <time> <label> (<norm_re>,<norm_im>)"
	//
	// Lines that cannot be parsed as this pattern are silently skipped.
	// We keep the last measurement per (time_step, label) pair at site 0
	// (impurity), which corresponds to the most-swept (most converged) sweep.
	void parseTdmrgLog(const std::string& logfile)
	{
		std::ifstream fin(logfile);
		if (!fin || !fin.good()) {
			err("ImpuritySolverNeqTdmrg: cannot open tDMRG log '" + logfile + "'\n");
			return;
		}

		// Maps from integer time step index n → measurement value at site 0.
		// Later entries overwrite earlier ones (last-wins = most converged).
		std::map<int, RealType>    nup_at_step;
		std::map<int, ComplexType> ggt0_at_step; // <gs|c|P1> at impurity

		std::string line;
		while (std::getline(fin, line)) {
			// Parse: "<site> <value> <time> <label> <norm>"
			std::istringstream iss(line);
			SizeType           site = 0;
			if (!(iss >> site))
				continue;
			if (site != 0)
				continue; // impurity site only

			std::string valStr;
			if (!(iss >> valStr))
				continue;

			RealType re = 0, im = 0;
			if (!parseComplex(valStr, re, im))
				continue;

			RealType t = 0;
			if (!(iss >> t))
				continue;

			std::string label;
			if (!(iss >> label))
				continue;

			// Convert physical time to time-step index, rounding to nearest integer.
			const int n = static_cast<int>(std::round(t / params_.dt));
			if (n < 0 || n > static_cast<int>(params_.nT))
				continue;

			if (label == "<gs|nup|gs>") {
				nup_at_step[n] = re;
			} else if (label == "<gs|c|P1>") {
				ggt0_at_step[n] = ComplexType(re, im);
			}
		}

		if (nup_at_step.empty() && ggt0_at_step.empty()) {
			std::cerr << "ImpuritySolverNeqTdmrg: WARNING: no in-situ measurements "
			             "found in '" << logfile << "'.  Check TSP parameters and "
			             "finite-loop count.\n";
			return;
		}

		fillKBGrid(nup_at_step, ggt0_at_step);
	}

	void fillKBGrid(const std::map<int, RealType>&    nup_at_step,
	                const std::map<int, ComplexType>& ggt0_at_step)
	{
		const int       nT    = static_cast<int>(params_.nT);
		const ComplexType mI   = ComplexType(0, -1); // −i

		for (int n = 0; n <= nT; ++n) {
			// Diagonal G^<(n,n) = i · n(t_n)  (n total = 2 × n_up by SU(2) symmetry)
			auto itN = nup_at_step.find(n);
			if (itN != nup_at_step.end()) {
				const RealType ntot = RealType(2) * itN->second;
				gimp_.lesser(n, n)  = ComplexType(0, ntot);
			}

			// G^>(n,0) = −i · <gs|c|P1>  (spin-up channel; SU(2) gives full GF via ×2)
			// G^R(n,0) ≈ G^>(n,0) — placeholder until G^<(n,0) is computed.
			auto itG = ggt0_at_step.find(n);
			if (itG != ggt0_at_step.end()) {
				const ComplexType ggt0 = mI * itG->second * RealType(2);
				gimp_.retarded(n, 0)  = ggt0;
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
