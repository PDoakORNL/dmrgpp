#ifndef IMPURITYSOLVER_NEQ_LANCZOS_H
#define IMPURITYSOLVER_NEQ_LANCZOS_H

#include "CincuentaInputCheck.h"
#include "ImpuritySolverNeqBase.h"
#include "InputCheck.h"
#include "KadanoffBaym.h"
#include "LanczosPlusPlus/src/Engine/DefaultSymmetry.h"
#include "LanczosPlusPlus/src/Engine/InternalProductOnTheFly.h"
#include "LanczosPlusPlus/src/Engine/InternalProductStored.h"
#include "LanczosPlusPlus/src/Engine/LabeledOperator.h"
#include "LanczosPlusPlus/src/Engine/LanczosGlobals.h"
#include "LanczosPlusPlus/src/Engine/ModelSelector.h"
#include "LanczosSolver.h"
#include "Matrix.h"
#include "ParametersForSolver.h"
#include "ParamsNeqDmftSolver.h"
#include "PsimagLite.h"
#include "Vector.h"
#include <cassert>
#include <cmath>
#include <complex>
#include <memory>

namespace Dmft {

// Non-equilibrium impurity solver using truncated Lanczos diagonalization.
//
// Builds an approximate Lehmann representation keeping only the nStates_
// lowest eigenstates per particle-number sector.  This extends
// ImpuritySolverNeqExactDiag to problems where full diagonalization is
// infeasible: Lanczos scales as O(k * D) where k << D is the number of
// kept states, while fullDiag scales as O(D^3).
//
// Physics: interaction quench U_i -> U_f at t=0, fixed bath.
// Green's functions are computed from the truncated Lehmann representation.
//
// Required extra input parameter:
//   NstatesNeq=<int>  — number of Lanczos eigenstates to keep per sector.
//
// Accuracy improves monotonically with NstatesNeq.  For NstatesNeq >= sector
// dimension the result is identical to ImpuritySolverNeqExactDiag.
template <typename ComplexOrRealType>
class ImpuritySolverNeqLanczos : public ImpuritySolverNeqBase<ComplexOrRealType> {

public:

	using BaseType    = ImpuritySolverNeqBase<ComplexOrRealType>;
	using RealType    = typename BaseType::RealType;
	using ComplexType = typename BaseType::ComplexType;
	using VectorRealType    = typename BaseType::VectorRealType;
	using KBType            = typename BaseType::KBType;
	using InputNgType       = typename BaseType::InputNgType;
	using VectorComplexType = typename PsimagLite::Vector<ComplexType>::Type;
	using MatrixType        = PsimagLite::Matrix<ComplexOrRealType>;
	using MatrixComplexType = PsimagLite::Matrix<ComplexType>;
	using ParamsNeqType     = ParamsNeqDmftSolver<ComplexOrRealType>;

	using DmrgInputReadable = typename PsimagLite::InputNg<Dmrg::InputCheck>::Readable;
	using GeometryType = PsimagLite::
	    Geometry<ComplexOrRealType, DmrgInputReadable, LanczosPlusPlus::LanczosGlobals>;
	using ModelSelectorType =
	    LanczosPlusPlus::ModelSelector<ComplexOrRealType, GeometryType, DmrgInputReadable>;
	using ModelBaseType =
	    LanczosPlusPlus::ModelBase<ComplexOrRealType, GeometryType, DmrgInputReadable>;
	using BasisBaseType         = typename ModelBaseType::BasisBaseType;
	using DefaultSymmetryType   = LanczosPlusPlus::DefaultSymmetry<GeometryType, BasisBaseType>;
	using InternalProductOnTheFlyType =
	    LanczosPlusPlus::InternalProductOnTheFly<ModelBaseType, DefaultSymmetryType>;
	using InternalProductStoredType =
	    LanczosPlusPlus::InternalProductStored<ModelBaseType, DefaultSymmetryType>;
	using LabeledOperatorType = LanczosPlusPlus::LabeledOperator;
	using WordType            = LanczosPlusPlus::LanczosGlobals::WordType;
	using PairIntType         = LanczosPlusPlus::LanczosGlobals::PairIntType;

	using VectorType         = typename PsimagLite::Vector<ComplexOrRealType>::Type;
	using VectorVectorType   = typename PsimagLite::Vector<VectorType>::Type;
	using ParametersSolverType = PsimagLite::ParametersForSolver<RealType>;
	using LanczosSolverType    = PsimagLite::LanczosSolver<InternalProductOnTheFlyType>;

	ImpuritySolverNeqLanczos(const ParamsNeqType&            params,
	                       typename InputNgType::Readable& io)
	    : params_(params)
	    , io_(io)
	    , nup_(0)
	    , ndown_(0)
	    , nT_(params.nT)
	    , nTau_(params.eqParams.nMatsubaras)
	    , dtau_(params.eqParams.ficticiousBeta /
	            static_cast<RealType>(params.eqParams.nMatsubaras))
	    , nStates_(0)
	    , gimp_(params.nT,
	            params.eqParams.nMatsubaras,
	            params.dt,
	            params.eqParams.ficticiousBeta /
	                static_cast<RealType>(params.eqParams.nMatsubaras))
	    , E0_pre_(0)
	{
		io.readline(nup_,     "TargetElectronsUp=");
		io.readline(ndown_,   "TargetElectronsDown=");
		io.readline(nStates_, "NstatesNeq=");
		if (nStates_ == 0)
			err("ImpuritySolverNeqLanczos: NstatesNeq must be > 0\n");
	}

	// One-time setup: truncated Lanczos diagonalization of H_pre and H_post,
	// then build the approximate Lehmann representation.
	void initialize(const VectorRealType& bathParams) override
	{
		const SizeType nBath   = bathParams.size() / 2;
		const SizeType nsites  = nBath + 1;
		const SizeType impSite = 0;

		assert(nup_ > 0);

		VectorRealType hoppings(nBath), bathEps(nBath);
		for (SizeType i = 0; i < nBath; ++i) {
			hoppings[i] = bathParams[i];
			bathEps[i]  = bathParams[nBath + i];
		}

		VectorRealType potPre(nsites, 0), potPost(nsites, 0);
		potPre[impSite]  = -RealType(0.5) * params_.uInitial;
		potPost[impSite] = -RealType(0.5) * params_.uFinal;
		for (SizeType i = 0; i < nBath; ++i) {
			potPre[i + 1]  = bathEps[i];
			potPost[i + 1] = bathEps[i];
		}

		// ---- Pre-quench: N and N+1 sectors --------------------------------
		VectorRealType energiesPreN;
		MatrixType     eigvecsPreN, eigvecsPreN1;
		{
			const std::string inputStr = buildLanczosInput(
			    params_.uInitial, nup_, ndown_, hoppings, potPre, nsites);
			Dmrg::InputCheck ic;
			typename PsimagLite::InputNg<Dmrg::InputCheck>::Writeable ioW(ic, inputStr);
			DmrgInputReadable ioR(ioW);
			GeometryType      geom(ioR);
			ModelSelectorType ms(ioR, geom);
			const ModelBaseType& model = ms();

			diagWithBasisTruncated(model, model.basis(), geom, energiesPreN, eigvecsPreN);
			E0_pre_ = energiesPreN[0];

			std::unique_ptr<BasisBaseType> bN1pre(model.createBasis(nup_ + 1, ndown_));
			diagWithBasisTruncated(model, *bN1pre, geom, energiesN1_pre_, eigvecsPreN1);

			buildF(model.basis(), *bN1pre, eigvecsPreN, eigvecsPreN1, impSite);
		}

		// ---- Post-quench: N, N+1, and N-1 sectors -------------------------
		MatrixType eigvecsPostN, eigvecsPostN1, eigvecsPostNm1;
		{
			const std::string inputStr = buildLanczosInput(
			    params_.uFinal, nup_, ndown_, hoppings, potPost, nsites);
			Dmrg::InputCheck ic;
			typename PsimagLite::InputNg<Dmrg::InputCheck>::Writeable ioW(ic, inputStr);
			DmrgInputReadable ioR(ioW);
			GeometryType      geom(ioR);
			ModelSelectorType ms(ioR, geom);
			const ModelBaseType& model = ms();

			diagWithBasisTruncated(model, model.basis(), geom, energiesN_post_, eigvecsPostN);

			std::unique_ptr<BasisBaseType> bN1post(model.createBasis(nup_ + 1, ndown_));
			diagWithBasisTruncated(model, *bN1post, geom, energiesN1_post_, eigvecsPostN1);

			std::unique_ptr<BasisBaseType> bNm1post(model.createBasis(nup_ - 1, ndown_));
			diagWithBasisTruncated(model, *bNm1post, geom, energiesNm1_post_, eigvecsPostNm1);

			// b_n = <n^N_post | GS_pre>
			// Outer: kept post-quench N states; inner: full Fock basis
			const SizeType dimN_full   = eigvecsPostN.rows();
			const SizeType nKeptN_post = eigvecsPostN.cols();
			b_.resize(nKeptN_post, ComplexType(0));
			for (SizeType n = 0; n < nKeptN_post; ++n) {
				ComplexType s = 0;
				for (SizeType i = 0; i < dimN_full; ++i)
					s += std::conj(ComplexType(eigvecsPostN(i, n)))
					   * ComplexType(eigvecsPreN(i, 0));
				b_[n] = s;
			}

			buildPhiPsi(model.basis(), *bN1post, *bNm1post,
			            eigvecsPostN, eigvecsPostN1, eigvecsPostNm1, impSite);
		}

		buildON1(eigvecsPreN1, eigvecsPostN1);
		buildChi();
		buildMatsubara();
	}

	// Fill the n-th time-slice of gimp (retarded, lesser, left-mixing).
	void computeGimp(KBType& gimp, int n) const override
	{
		for (int j = 0; j <= n; ++j)
			gimp.retarded(n, j) = gGreater(n, j) - gLesser(n, j);

		for (int j = 0; j <= n; ++j)
			gimp.lesser(n, j) = gLesser(n, j);
		for (int j = 0; j < n; ++j)
			gimp.lesser(j, n) = -std::conj(gimp.lesser(n, j));

		for (SizeType j = 0; j <= nTau_; ++j)
			gimp.left_mixing(n, j) = gLeft(n, j);
	}

	const KBType& gimp() const override { return gimp_; }

private:

	// Truncated Lanczos diagonalization of the Hamiltonian in the given basis.
	// Returns up to min(nStates_, dim) lowest eigenpairs.
	// eigvecs has shape (dim, nKept): eigvecs(j, i) = component j of eigenvector i.
	// Falls back to fullDiag (via InternalProductStored) if Lanczos fails or
	// the sector is smaller than nStates_.
	void diagWithBasisTruncated(const ModelBaseType& model,
	                            const BasisBaseType& basis,
	                            const GeometryType&  geom,
	                            VectorRealType&      eigs,
	                            MatrixType&          eigvecs)
	{
		const SizeType dim   = basis.size();
		const SizeType nKeep = std::min(nStates_, dim);

		DefaultSymmetryType         rs(basis, geom, "");
		InternalProductOnTheFlyType ham(model, basis, rs);

		VectorType       initial(dim);
		PsimagLite::fillRandom(initial);
		VectorVectorType zs(nKeep, VectorType(dim));
		eigs.resize(nKeep);

		ParametersSolverType lparams;
		lparams.lotaMemory = true;
		lparams.options    = "reortho"; // full re-orthogonalization keeps eigenvectors orthonormal
		LanczosSolverType lanczosSolver(ham, lparams);

		try {
			lanczosSolver.computeAllStatesBelow(eigs, zs, initial, nKeep);
			// diag() inside computeAllStatesBelow resizes eigs to the full
			// tridiagonal matrix dimension (Lanczos steps), not nKeep.
			// Truncate to the nKeep lowest eigenvalues, which are in eigs[0..nKeep-1].
			eigs.resize(nKeep);
		} catch (std::exception&) {
			// Sector is small or Lanczos failed: use fullDiag instead.
			InternalProductStoredType hamS(model, basis, rs);
			VectorRealType eigs2(dim);
			MatrixType     fm;
			hamS.fullDiag(eigs2, fm);
			eigs.resize(nKeep);
			zs.resize(nKeep, VectorType(dim));
			for (SizeType i = 0; i < nKeep; ++i) {
				eigs[i] = eigs2[i];
				for (SizeType j = 0; j < dim; ++j)
					zs[i][j] = fm(j, i);
			}
		}

		// Pack into column-major matrix: eigvecs(j, i) = i-th eigenvector at index j
		// eigs.size() == nKeep; zs[0..nKeep-1] are valid.
		eigvecs.resize(dim, eigs.size());
		for (SizeType i = 0; i < eigs.size(); ++i)
			for (SizeType j = 0; j < dim; ++j)
				eigvecs(j, i) = zs[i][j];
	}

	// f_[l] = <l^{N+1}_pre | c†_{imp,up} | GS_pre>
	// l runs over nKeptN1_pre kept eigenstates; inner sums over full Fock basis.
	void buildF(const BasisBaseType& basisN,
	            const BasisBaseType& basisN1,
	            const MatrixType&    eigvecsN,   // (dimN_full, nKeptN_pre)
	            const MatrixType&    eigvecsN1,  // (dimN1_full, nKeptN1_pre)
	            SizeType             impSite)
	{
		const SizeType dimN_full  = basisN.size();
		const SizeType dimN1_full = basisN1.size();
		const SizeType nKeptN1    = eigvecsN1.cols();

		const LabeledOperatorType opCdag(LabeledOperatorType::Label::OPERATOR_CDAGGER);
		const SizeType            spinUp = LanczosPlusPlus::LanczosGlobals::SPIN_UP;
		const SizeType            spinDn = LanczosPlusPlus::LanczosGlobals::SPIN_DOWN;

		MatrixType cdagFock(dimN1_full, dimN_full);
		for (SizeType m = 0; m < dimN_full; ++m) {
			WordType    ket1 = basisN(m, spinUp);
			WordType    ket2 = basisN(m, spinDn);
			PairIntType bra  = basisN1.getBraIndex(ket1, ket2, opCdag, impSite, spinUp, 0);
			if (bra.first < 0)
				continue;
			int sign = basisN.doSignGf(ket1, ket2, impSite, spinUp, 0);
			cdagFock(static_cast<SizeType>(bra.first), m) = ComplexOrRealType(sign);
		}

		// tmp[k] = sum_m cdagFock(k,m) * GS_pre[m]
		VectorComplexType tmp(dimN1_full, ComplexType(0));
		for (SizeType k = 0; k < dimN1_full; ++k)
			for (SizeType m = 0; m < dimN_full; ++m)
				tmp[k] += ComplexType(cdagFock(k, m)) * ComplexType(eigvecsN(m, 0));

		// f_[l] = sum_k conj(eigvecsN1(k,l)) * tmp[k]
		f_.resize(nKeptN1, ComplexType(0));
		for (SizeType l = 0; l < nKeptN1; ++l)
			for (SizeType k = 0; k < dimN1_full; ++k)
				f_[l] += std::conj(ComplexType(eigvecsN1(k, l))) * tmp[k];
	}

	// Build Phi_(k,n) and Psi_(k,n) for all time steps n = 0..nT.
	// Inner sums over Fock basis; outer indices over kept eigenstates.
	void buildPhiPsi(const BasisBaseType& basisN,
	                 const BasisBaseType& basisN1,
	                 const BasisBaseType& basisNm1,
	                 const MatrixType&    eigvecsN,    // (dimN_full,   nKeptN_post)
	                 const MatrixType&    eigvecsN1,   // (dimN1_full,  nKeptN1_post)
	                 const MatrixType&    eigvecsNm1,  // (dimNm1_full, nKeptNm1_post)
	                 SizeType             impSite)
	{
		const SizeType dimN_full   = basisN.size();
		const SizeType dimN1_full  = basisN1.size();
		const SizeType dimNm1_full = basisNm1.size();
		const SizeType nKeptN      = eigvecsN.cols();
		const SizeType nKeptN1     = eigvecsN1.cols();
		const SizeType nKeptNm1    = eigvecsNm1.cols();
		const SizeType nSteps      = nT_ + 1;
		const SizeType spinUp      = LanczosPlusPlus::LanczosGlobals::SPIN_UP;
		const SizeType spinDn      = LanczosPlusPlus::LanczosGlobals::SPIN_DOWN;
		const LabeledOperatorType opCdag(LabeledOperatorType::Label::OPERATOR_CDAGGER);
		const LabeledOperatorType opC(LabeledOperatorType::Label::OPERATOR_C);

		// -- Particle channel (c†) ------------------------------------------
		// Mcdag(k,m) = <k^{N+1}_post | c†_imp | m^N_post>
		// k in 0..nKeptN1-1, m in 0..nKeptN-1
		MatrixComplexType Mcdag(nKeptN1, nKeptN);
		{
			MatrixType cdagFock(dimN1_full, dimN_full);
			for (SizeType m = 0; m < dimN_full; ++m) {
				WordType    ket1 = basisN(m, spinUp);
				WordType    ket2 = basisN(m, spinDn);
				PairIntType bra  = basisN1.getBraIndex(ket1, ket2, opCdag, impSite, spinUp, 0);
				if (bra.first < 0)
					continue;
				int sign = basisN.doSignGf(ket1, ket2, impSite, spinUp, 0);
				cdagFock(static_cast<SizeType>(bra.first), m) = ComplexOrRealType(sign);
			}
			// tmp(k, m) = sum_a conj(eigvecsN1(a,k)) * cdagFock(a,m)
			// k in 0..nKeptN1-1, a in 0..dimN1_full-1, m in 0..dimN_full-1
			MatrixComplexType tmp(nKeptN1, dimN_full);
			for (SizeType k = 0; k < nKeptN1; ++k)
				for (SizeType m = 0; m < dimN_full; ++m)
					for (SizeType a = 0; a < dimN1_full; ++a)
						tmp(k, m) += std::conj(ComplexType(eigvecsN1(a, k)))
						           * ComplexType(cdagFock(a, m));
			// Mcdag(k,m) = sum_b tmp(k,b) * eigvecsN(b,m)
			// b in 0..dimN_full-1
			for (SizeType k = 0; k < nKeptN1; ++k)
				for (SizeType m = 0; m < nKeptN; ++m)
					for (SizeType b = 0; b < dimN_full; ++b)
						Mcdag(k, m) += tmp(k, b) * ComplexType(eigvecsN(b, m));
		}

		// Phi_(k, n) = sum_m Mcdag(k,m) * b_m * exp(i * omega_km * t_n)
		Phi_.resize(nKeptN1, nSteps);
		for (SizeType k = 0; k < nKeptN1; ++k) {
			for (SizeType n = 0; n < nSteps; ++n) {
				const RealType tn = n * params_.dt;
				ComplexType    s  = 0;
				for (SizeType m = 0; m < nKeptN; ++m) {
					const RealType phase = (energiesN1_post_[k] - energiesN_post_[m]) * tn;
					s += Mcdag(k, m) * b_[m]
					   * ComplexType(std::cos(phase), std::sin(phase));
				}
				Phi_(k, n) = s;
			}
		}

		// -- Hole channel (c) -----------------------------------------------
		// Mc(k,m) = <k^{N-1}_post | c_imp | m^N_post>
		MatrixComplexType Mc(nKeptNm1, nKeptN);
		{
			MatrixType cFock(dimNm1_full, dimN_full);
			for (SizeType m = 0; m < dimN_full; ++m) {
				WordType    ket1 = basisN(m, spinUp);
				WordType    ket2 = basisN(m, spinDn);
				PairIntType bra  = basisNm1.getBraIndex(ket1, ket2, opC, impSite, spinUp, 0);
				if (bra.first < 0)
					continue;
				int sign = basisN.doSignGf(ket1, ket2, impSite, spinUp, 0);
				cFock(static_cast<SizeType>(bra.first), m) = ComplexOrRealType(sign);
			}
			MatrixComplexType tmp(nKeptNm1, dimN_full);
			for (SizeType k = 0; k < nKeptNm1; ++k)
				for (SizeType m = 0; m < dimN_full; ++m)
					for (SizeType a = 0; a < dimNm1_full; ++a)
						tmp(k, m) += std::conj(ComplexType(eigvecsNm1(a, k)))
						           * ComplexType(cFock(a, m));
			for (SizeType k = 0; k < nKeptNm1; ++k)
				for (SizeType m = 0; m < nKeptN; ++m)
					for (SizeType b = 0; b < dimN_full; ++b)
						Mc(k, m) += tmp(k, b) * ComplexType(eigvecsN(b, m));
		}

		// Psi_(k, n) = sum_m Mc(k,m) * b_m * exp(i * omega_km * t_n)
		Psi_.resize(nKeptNm1, nSteps);
		for (SizeType k = 0; k < nKeptNm1; ++k) {
			for (SizeType n = 0; n < nSteps; ++n) {
				const RealType tn = n * params_.dt;
				ComplexType    s  = 0;
				for (SizeType m = 0; m < nKeptN; ++m) {
					const RealType phase = (energiesNm1_post_[k] - energiesN_post_[m]) * tn;
					s += Mc(k, m) * b_[m]
					   * ComplexType(std::cos(phase), std::sin(phase));
				}
				Psi_(k, n) = s;
			}
		}
	}

	// O_N1_(k, l) = <k^{N+1}_post | l^{N+1}_pre>
	// k in 0..nKeptN1_post-1, l in 0..nKeptN1_pre-1
	// inner sum over full N+1 Hilbert space
	void buildON1(const MatrixType& eigvecsPreN1, const MatrixType& eigvecsPostN1)
	{
		const SizeType dimN1_full = eigvecsPostN1.rows();
		const SizeType nKeptPre   = eigvecsPreN1.cols();
		const SizeType nKeptPost  = eigvecsPostN1.cols();

		O_N1_.resize(nKeptPost, nKeptPre);
		for (SizeType k = 0; k < nKeptPost; ++k)
			for (SizeType l = 0; l < nKeptPre; ++l) {
				ComplexType s = 0;
				for (SizeType i = 0; i < dimN1_full; ++i)
					s += std::conj(ComplexType(eigvecsPostN1(i, k)))
					   * ComplexType(eigvecsPreN1(i, l));
				O_N1_(k, l) = s;
			}
	}

	// chi_(k, j) = sum_l O_N1(k,l) * f_l * exp(Omega_l * tau_j)
	void buildChi()
	{
		const SizeType dimN1post = O_N1_.rows();
		const SizeType dimN1pre  = energiesN1_pre_.size();
		const SizeType nTauSteps = nTau_ + 1;

		chi_.resize(dimN1post, nTauSteps);
		for (SizeType k = 0; k < dimN1post; ++k) {
			for (SizeType j = 0; j < nTauSteps; ++j) {
				const RealType tau = j * dtau_;
				ComplexType    s   = 0;
				for (SizeType l = 0; l < dimN1pre; ++l) {
					const RealType OmegaL = energiesN1_pre_[l] - E0_pre_;
					s += O_N1_(k, l) * f_[l] * std::exp(OmegaL * tau);
				}
				chi_(k, j) = s;
			}
		}
	}

	// G^M(tau) and G^M(i*omega_k) from the pre-quench N+1 Lehmann representation.
	void buildMatsubara()
	{
		const SizeType dimN1pre  = energiesN1_pre_.size();
		const SizeType nTauSteps = nTau_ + 1;
		const SizeType nFreqs    = nTau_;
		const RealType beta      = params_.eqParams.ficticiousBeta;

		for (SizeType j = 0; j < nTauSteps; ++j) {
			const RealType tau = j * dtau_;
			ComplexType    gm  = 0;
			for (SizeType l = 0; l < dimN1pre; ++l) {
				const RealType OmegaL = energiesN1_pre_[l] - E0_pre_;
				gm -= RealType(std::norm(f_[l])) * std::exp(-OmegaL * tau);
			}
			gimp_.matsubara_t[j] = gm;
		}

		for (SizeType k = 0; k < nFreqs; ++k) {
			const int      kInt   = static_cast<int>(k);
			const int      nInt   = static_cast<int>(nFreqs);
			const RealType omegaK = RealType(2 * kInt - nInt + 1) * M_PI / beta;
			ComplexType    gw     = 0;
			for (SizeType l = 0; l < dimN1pre; ++l) {
				const RealType OmegaL = energiesN1_pre_[l] - E0_pre_;
				gw += RealType(std::norm(f_[l]))
				    / (ComplexType(0, omegaK) - OmegaL);
			}
			gimp_.matsubara_w[k] = gw;
		}
	}

	// G^>(t_n, t_j) = -i sum_k conj(Phi_k[n]) * Phi_k[j]
	ComplexType gGreater(int n, int j) const
	{
		const SizeType nKeptN1 = Phi_.rows();
		ComplexType    s       = 0;
		for (SizeType k = 0; k < nKeptN1; ++k)
			s += std::conj(Phi_(k, static_cast<SizeType>(n)))
			   * Phi_(k, static_cast<SizeType>(j));
		return ComplexType(0, -1) * s;
	}

	// G^<(t_n, t_j) = i sum_k conj(Psi_k[j]) * Psi_k[n]
	ComplexType gLesser(int n, int j) const
	{
		const SizeType nKeptNm1 = Psi_.rows();
		ComplexType    s        = 0;
		for (SizeType k = 0; k < nKeptNm1; ++k)
			s += std::conj(Psi_(k, static_cast<SizeType>(j)))
			   * Psi_(k, static_cast<SizeType>(n));
		return ComplexType(0, 1) * s;
	}

	// G^{Left}(t_n, tau_j) = -i sum_k chi_k[j] * conj(Phi_k[n])
	ComplexType gLeft(int n, SizeType j) const
	{
		const SizeType nKeptN1 = chi_.rows();
		ComplexType    s       = 0;
		for (SizeType k = 0; k < nKeptN1; ++k)
			s += chi_(k, j) * std::conj(Phi_(k, static_cast<SizeType>(n)));
		return ComplexType(0, -1) * s;
	}

	// Build LanczosPlusPlus Anderson-model input string (same as ExactDiag version).
	static std::string buildLanczosInput(RealType              U,
	                                     SizeType              nup,
	                                     SizeType              ndown,
	                                     const VectorRealType& hoppings,
	                                     const VectorRealType& potV,
	                                     SizeType              nsites)
	{
		std::string uStr = "[" + ttos(U);
		for (SizeType i = 1; i < nsites; ++i)
			uStr += ", 0.";
		uStr += "]";

		std::string connStr = "[";
		for (SizeType i = 0; i < hoppings.size(); ++i) {
			if (i > 0) connStr += ",";
			connStr += ttos(hoppings[i]);
		}
		connStr += "]";

		std::string potStr = "[";
		for (SizeType i = 0; i < nsites; ++i) {
			if (i > 0) potStr += ",";
			potStr += ttos(potV[i]);
		}
		potStr += ",";
		for (SizeType i = 0; i < nsites; ++i) {
			if (i > 0) potStr += ",";
			potStr += ttos(potV[i]);
		}
		potStr += "]";

		std::string s = "##Ainur1.0\n\n";
		s += "TotalNumberOfSites=" + ttos(nsites) + ";\n";
		s += "NumberOfTerms=1;\n";
		s += "DegreesOfFreedom=1;\n";
		s += "GeometryKind=star;\n";
		s += "GeometryOptions=none;\n";
		s += "hubbardU=" + uStr + ";\n";
		s += "Model=HubbardOneBand;\n";
		s += "SolverOptions=twositedmrg,geometryallinsystem,hd5dontprint;\n";
		s += "Version=templateForDMFT;\n";
		s += "OutputFile=neqDmrgDummy;\n";
		s += "InfiniteLoopKeptStates=1;\n";
		s += "FiniteLoops=0 0 0;\n";
		s += "TargetElectronsUp=" + ttos(nup) + ";\n";
		s += "TargetElectronsDown=" + ttos(ndown) + ";\n";
		s += "dir0:Connectors=" + connStr + ";\n";
		s += "potentialV=" + potStr + ";\n";
		return s;
	}

	// ---- Member variables --------------------------------------------------
	const ParamsNeqType&            params_;
	typename InputNgType::Readable& io_;
	SizeType                        nup_;
	SizeType                        ndown_;
	SizeType                        nT_;
	SizeType                        nTau_;
	RealType                        dtau_;
	SizeType                        nStates_; // Lanczos states to keep per sector

	// Pre-quench N+1 spectrum and operator amplitudes
	VectorRealType    energiesN1_pre_;
	RealType          E0_pre_;
	VectorComplexType f_;

	// Post-quench spectra (sizes = nKept per sector)
	VectorRealType energiesN_post_, energiesN1_post_, energiesNm1_post_;

	// Quench overlaps: b_n = <n^N_post | GS_pre>
	VectorComplexType b_;

	// Overlap matrix: O_N1_(k,l) = <k^{N+1}_post | l^{N+1}_pre>
	MatrixComplexType O_N1_;

	// Time-dependent amplitudes (rows = kept states, cols = time steps)
	MatrixComplexType Phi_; // particle channel
	MatrixComplexType Psi_; // hole channel
	MatrixComplexType chi_; // imaginary-time factor for G^{Left}

	KBType gimp_;
};

} // namespace Dmft
#endif // IMPURITYSOLVER_NEQ_LANCZOS_H
