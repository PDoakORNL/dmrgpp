#ifndef PARAMS_NEQ_DMFT_SOLVER_H
#define PARAMS_NEQ_DMFT_SOLVER_H
#include "CincuentaInputCheck.h"
#include "ParamsDmftSolver.h"

namespace Dmft {

// Parameters for the non-equilibrium DMFT extension.
// Contains all equilibrium parameters (via composition) plus the
// interaction-quench and real-time grid parameters.
template <typename ComplexOrRealType>
struct ParamsNeqDmftSolver {

	using RealType        = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using InputNgType     = PsimagLite::InputNg<CincuentaInputCheck>;
	using EqParamsType    = ParamsDmftSolver<ComplexOrRealType>;

	ParamsNeqDmftSolver(typename InputNgType::Readable& io)
	    : eqParams(io)
	{
		io.readline(uInitial, "HubbardU=");
		io.readline(uFinal,   "HubbardUFinal=");
		io.readline(tMax,     "TmaxNeq=");
		io.readline(nT,       "NtNeq=");
		dt = tMax / static_cast<RealType>(nT);

		try {
			io.readline(neqDmftIter, "NeqDmftIter=");
		} catch (std::exception&) {
			neqDmftIter = eqParams.dmftIter;
		}

		try {
			io.readline(neqDmftError, "NeqDmftTolerance=");
		} catch (std::exception&) {
			neqDmftError = eqParams.dmftError;
		}
	}

	// Equilibrium DMFT parameters (beta, mu, nBath, etc.)
	EqParamsType eqParams;

	// Interaction quench: U_i -> U_f at t = 0
	RealType uInitial = 0;
	RealType uFinal   = 0;

	// Real-time grid
	RealType tMax = 0;
	SizeType nT   = 0;
	RealType dt   = 0;

	// Inner DMFT convergence at each time step
	SizeType neqDmftIter  = 10;
	RealType neqDmftError = 1e-4;
};

} // namespace Dmft
#endif // PARAMS_NEQ_DMFT_SOLVER_H
