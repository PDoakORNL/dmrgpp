#ifndef PARAMS_NEQ_DMFT_SOLVER_H
#define PARAMS_NEQ_DMFT_SOLVER_H
#include "CincuentaInputCheck.h"
#include "ParamsDmftSolver.h"

namespace Dmft {

/*!
 * \brief ParamsNeqDmftSolver
 * Parameters for the non-equilibrium DMFT extension.
 * Contains all equilibrium parameters (via composition) plus the
 * interaction-quench and real-time grid parameters.
 */
template <typename ComplexOrRealType> struct ParamsNeqDmftSolver {

	using RealType     = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using InputNgType  = PsimagLite::InputNg<CincuentaInputCheck>;
	using EqParamsType = ParamsDmftSolver<ComplexOrRealType>;

	ParamsNeqDmftSolver(typename InputNgType::Readable& io)
	    : eqParams(io)
	{
		io.readline(uInitial, "HubbardU=");
		io.readline(uFinal, "HubbardUFinal=");
		io.readline(tMax, "TmaxNeq=");
		io.readline(nT, "NtNeq=");
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

		try {
			io.readline(bandwidthFinal, "BandwidthFinal=");
		} catch (std::exception&) { }

		try {
			io.readline(neqOutputPrefix, "NeqOutputPrefix=");
		} catch (std::exception&) { }
	}

	EqParamsType eqParams; ///< Equilibrium DMFT parameters (beta, mu, nBath, etc.)

	RealType uInitial = 0; ///< Interaction quench: U_i -> U_f at t = 0
	RealType uFinal   = 0; ///<

	RealType tMax = 0; ///< Real-time grid
	SizeType nT   = 0; ///<
	RealType dt   = 0; ///<

	SizeType neqDmftIter  = 10; ///< Inner DMFT convergence at each time step
	RealType neqDmftError = 1e-4; ///<

	/// Hopping quench: Bethe lattice bandwidth for t > 0.
	/// 0 (default) means no quench — use the equilibrium bandwidth from LatticeGf.
	RealType bandwidthFinal = 0;

	/// Optional prefix for output Green's function files.
	/// Empty (default) → "green-retarded" etc.  Non-empty → "{prefix}-green-retarded" etc.
	std::string neqOutputPrefix = "";
};

} // namespace Dmft
#endif // PARAMS_NEQ_DMFT_SOLVER_H
