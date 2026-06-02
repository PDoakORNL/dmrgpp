#ifndef IMPURITYSOLVER_NEQ_EXACTDIAG_H
#define IMPURITYSOLVER_NEQ_EXACTDIAG_H

#include "CincuentaInputCheck.h"
#include "ImpuritySolverNeqBase.h"
#include "InputCheck.h"
#include "KadanoffBaym.h"
#include "LanczosPlusPlus/src/Engine/DefaultSymmetry.h"
#include "LanczosPlusPlus/src/Engine/InternalProductStored.h"
#include "LanczosPlusPlus/src/Engine/LabeledOperator.h"
#include "LanczosPlusPlus/src/Engine/LanczosGlobals.h"
#include "LanczosPlusPlus/src/Engine/ModelSelector.h"
#include "Matrix.h"
#include "ParamsNeqDmftSolver.h"
#include "PsimagLite.h"
#include "Vector.h"
#include <cassert>
#include <cmath>
#include <complex>
#include <memory>

namespace Dmft {

// Non-equilibrium impurity solver using LanczosPlusPlus full diagonalization.
//
// Physics: interaction quench U_i -> U_f at t=0, bath parameters fixed at
// their equilibrium values.  Green's functions are computed from the exact
// Lehmann representation on the Keldysh-Matsubara contour.
//
// initialize(bathParams) must be called once with the equilibrium bath
// parameters before the time-stepping loop.  computeGimp(gimp, n) then fills
// the n-th time slice of gimp for each real-time step n = 0..nT.
template <typename ComplexOrRealType>
class ImpuritySolverNeqExactDiag : public ImpuritySolverNeqBase<ComplexOrRealType> {

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
	    LanczosPlusPlus::LanczosModelBase<ComplexOrRealType, GeometryType, DmrgInputReadable>;
	using BasisBaseType         = typename ModelBaseType::BasisBaseType;
	using DefaultSymmetryType   = LanczosPlusPlus::DefaultSymmetry<GeometryType, BasisBaseType>;
	using InternalProductStoredType =
	    LanczosPlusPlus::InternalProductStored<ModelBaseType, DefaultSymmetryType>;
	using LabeledOperatorType = LanczosPlusPlus::LabeledOperator;
	using WordType            = LanczosPlusPlus::LanczosGlobals::WordType;
	using PairIntType         = LanczosPlusPlus::LanczosGlobals::PairIntType;

	ImpuritySolverNeqExactDiag(const ParamsNeqType&            params,
	                            typename InputNgType::Readable& io)
	    : params_(params)
	    , io_(io)
	    , nup_(0)
	    , ndown_(0)
	    , nT_(params.nT)
	    , nTau_(params.eqParams.nMatsubaras)
	    , dtau_(params.eqParams.ficticiousBeta /
	            static_cast<RealType>(params.eqParams.nMatsubaras))
	    , gimp_(params.nT,
	            params.eqParams.nMatsubaras,
	            params.dt,
	            params.eqParams.ficticiousBeta /
	                static_cast<RealType>(params.eqParams.nMatsubaras))
	    , E0_pre_(0)
	{
		io.readline(nup_,   "TargetElectronsUp=");
		io.readline(ndown_, "TargetElectronsDown=");
	}

	// Diagonalize H_pre and H_post for the given bath,
	// precompute all time-dependent spectral amplitudes.
	void solve(const VectorRealType& bathParams) override
	{
		const SizeType nBath   = bathParams.size() / 2;
		const SizeType nsites  = nBath + 1;
		const SizeType impSite = 0; // star center

		assert(nup_ > 0); // need N-1 sector

		VectorRealType hoppings(nBath), bathEps(nBath);
		for (SizeType i = 0; i < nBath; ++i) {
			hoppings[i] = bathParams[i];
			bathEps[i]  = bathParams[nBath + i];
		}

		// Anderson impurity at half filling: mu_imp = -U/2
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

			diagWithBasis(model, model.basis(), geom, energiesPreN, eigvecsPreN);
			E0_pre_ = energiesPreN[0];

			std::unique_ptr<BasisBaseType> bN1pre(model.createBasis(nup_ + 1, ndown_));
			diagWithBasis(model, *bN1pre, geom, energiesN1_pre_, eigvecsPreN1);

			buildF(model.basis(), *bN1pre, eigvecsPreN, eigvecsPreN1, impSite);

			// N-1 sector: hole contributions to G^M
			std::unique_ptr<BasisBaseType> bNm1pre(model.createBasis(nup_ - 1, ndown_));
			MatrixType eigvecsNm1_pre;
			diagWithBasis(model, *bNm1pre, geom, energiesNm1_pre_, eigvecsNm1_pre);
			buildHPre(model.basis(), *bNm1pre, eigvecsPreN, eigvecsNm1_pre, impSite);
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

			diagWithBasis(model, model.basis(), geom, energiesN_post_, eigvecsPostN);

			std::unique_ptr<BasisBaseType> bN1post(model.createBasis(nup_ + 1, ndown_));
			diagWithBasis(model, *bN1post, geom, energiesN1_post_, eigvecsPostN1);

			std::unique_ptr<BasisBaseType> bNm1post(model.createBasis(nup_ - 1, ndown_));
			diagWithBasis(model, *bNm1post, geom, energiesNm1_post_, eigvecsPostNm1);

			// b_n = <n^N_post | GS_pre>  (GS_pre = col 0 of eigvecsPreN)
			const SizeType dimN = eigvecsPostN.rows();
			b_.resize(dimN, ComplexType(0));
			for (SizeType n = 0; n < dimN; ++n) {
				ComplexType s = 0;
				for (SizeType i = 0; i < dimN; ++i)
					s += std::conj(ComplexType(eigvecsPostN(i, n)))
					   * ComplexType(eigvecsPreN(i, 0));
				b_[n] = s;
			}

			buildPhiPsi(model.basis(), *bN1post, *bNm1post,
			            eigvecsPostN, eigvecsPostN1, eigvecsPostNm1, impSite);
		}

		// O^{N+1}_{kl} = <k^{N+1}_post | l^{N+1}_pre>
		buildON1(eigvecsPreN1, eigvecsPostN1);

		// chi_k(tau_j) = sum_l O_{kl} f_l exp(Omega_l tau_j)
		buildChi();

		// Populate the Matsubara components of gimp_ from pre-quench spectrum
		buildMatsubara();
	}

	// Fill the n-th time-slice of gimp (retarded, lesser, left-mixing).
	// Precondition: initialize() has been called.
	void computeGimp(KBType& gimp, int n) const override
	{
		// Retarded: G^R(n,j) = G^>(n,j) - G^<(n,j),  j <= n
		for (int j = 0; j <= n; ++j)
			gimp.retarded(n, j) = gGreater(n, j) - gLesser(n, j);

		// Lesser: G^<(n,j) for j <= n, then fill column n by anti-Hermiticity
		for (int j = 0; j <= n; ++j)
			gimp.lesser(n, j) = gLesser(n, j);
		for (int j = 0; j < n; ++j)
			gimp.lesser(j, n) = -std::conj(gimp.lesser(n, j));

		// Left-mixing: G^{Left}(n, tau_j) for j = 0..nTau
		for (SizeType j = 0; j <= nTau_; ++j)
			gimp.left_mixing(n, j) = gLeft(n, j);
	}

	const KBType& gimp() const override { return gimp_; }

private:

	// Build a minimal LanczosPlusPlus input string for the Anderson model
	// with explicit U, electron numbers, bath parameters.
	static std::string buildLanczosInput(RealType              U,
	                                     SizeType              nup,
	                                     SizeType              ndown,
	                                     const VectorRealType& hoppings,
	                                     const VectorRealType& potV,
	                                     SizeType              nsites)
	{
		// hubbardU vector: U at site 0 (impurity), 0 elsewhere
		std::string uStr = "[" + ttos(U);
		for (SizeType i = 1; i < nsites; ++i)
			uStr += ", 0.";
		uStr += "]";

		// star geometry connectors
		std::string connStr = "[";
		for (SizeType i = 0; i < hoppings.size(); ++i) {
			if (i > 0) connStr += ",";
			connStr += ttos(hoppings[i]);
		}
		connStr += "]";

		// potentialV: doubled for spin-up and spin-down blocks
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
		s += "OutputFile=neqExactDiagDummy;\n";
		s += "InfiniteLoopKeptStates=1;\n";
		s += "FiniteLoops=0 0 0;\n";
		s += "TargetElectronsUp=" + ttos(nup) + ";\n";
		s += "TargetElectronsDown=" + ttos(ndown) + ";\n";
		s += "dir0:Connectors=" + connStr + ";\n";
		s += "potentialV=" + potStr + ";\n";
		return s;
	}

	// Full diagonalization of the Hamiltonian restricted to basis.
	static void diagWithBasis(const ModelBaseType& model,
	                          const BasisBaseType& basis,
	                          const GeometryType&  geom,
	                          VectorRealType&      eigs,
	                          MatrixType&          eigvecs)
	{
		DefaultSymmetryType       rs(basis, geom, "");
		InternalProductStoredType ham(model, basis, rs);
		const SizeType            dim = basis.size();
		eigs.resize(dim);
		eigvecs.resize(dim, dim);
		ham.fullDiag(eigs, eigvecs);
	}

	// f_[l] = <l^{N+1}_pre | c†_{imp,up} | GS_pre>
	// Uses the pre-quench N and N+1 bases and eigenvectors.
	void buildF(const BasisBaseType& basisN,
	            const BasisBaseType& basisN1,
	            const MatrixType&    eigvecsN,
	            const MatrixType&    eigvecsN1,
	            SizeType             impSite)
	{
		const SizeType dimN  = basisN.size();
		const SizeType dimN1 = basisN1.size();
		const LabeledOperatorType opCdag(LabeledOperatorType::Label::OPERATOR_CDAGGER);
		const SizeType            spinUp = LanczosPlusPlus::LanczosGlobals::SPIN_UP;
		const SizeType            spinDn = LanczosPlusPlus::LanczosGlobals::SPIN_DOWN;

		// Fock-basis matrix: cdagFock(k,m) = <k^{N+1}|c†|m^N> in Fock repr.
		MatrixType cdagFock(dimN1, dimN);
		for (SizeType m = 0; m < dimN; ++m) {
			WordType    ket1 = basisN(m, spinUp);
			WordType    ket2 = basisN(m, spinDn);
			PairIntType bra  = basisN1.getBraIndex(ket1, ket2, opCdag, impSite, spinUp, 0);
			if (bra.first < 0)
				continue;
			int sign = basisN.doSignGf(ket1, ket2, impSite, spinUp, 0);
			cdagFock(static_cast<SizeType>(bra.first), m) = ComplexOrRealType(sign);
		}

		// Two-step transform: f = eigvecsN1† * cdagFock * eigvecsN[:, 0]
		// f_l = sum_{k,m} conj(eigvecsN1(k,l)) * cdagFock(k,m) * eigvecsN(m,0)
		// Step 1: tmp(k) = sum_m cdagFock(k,m) * eigvecsN(m,0)
		VectorComplexType tmp(dimN1, ComplexType(0));
		for (SizeType k = 0; k < dimN1; ++k)
			for (SizeType m = 0; m < dimN; ++m)
				tmp[k] += ComplexType(cdagFock(k, m)) * ComplexType(eigvecsN(m, 0));

		// Step 2: f_l = sum_k conj(eigvecsN1(k,l)) * tmp(k)
		f_.resize(dimN1, ComplexType(0));
		for (SizeType l = 0; l < dimN1; ++l)
			for (SizeType k = 0; k < dimN1; ++k)
				f_[l] += std::conj(ComplexType(eigvecsN1(k, l))) * tmp[k];
	}

	// h_pre_[k] = <k^{N-1}_pre | c_imp_up | GS^N_pre>  (pre-quench hole amplitudes)
	void buildHPre(const BasisBaseType& basisN,
	               const BasisBaseType& basisNm1,
	               const MatrixType&    eigvecsN,
	               const MatrixType&    eigvecsNm1,
	               SizeType             impSite)
	{
		const SizeType dimN   = basisN.size();
		const SizeType dimNm1 = basisNm1.size();
		const SizeType spinUp = LanczosPlusPlus::LanczosGlobals::SPIN_UP;
		const SizeType spinDn = LanczosPlusPlus::LanczosGlobals::SPIN_DOWN;
		const LabeledOperatorType opC(LabeledOperatorType::Label::OPERATOR_C);

		MatrixType cFock(dimNm1, dimN);
		for (SizeType m = 0; m < dimN; ++m) {
			WordType    ket1 = basisN(m, spinUp);
			WordType    ket2 = basisN(m, spinDn);
			PairIntType bra  = basisNm1.getBraIndex(ket1, ket2, opC, impSite, spinUp, 0);
			if (bra.first < 0)
				continue;
			int sign = basisN.doSignGf(ket1, ket2, impSite, spinUp, 0);
			cFock(static_cast<SizeType>(bra.first), m) = ComplexOrRealType(sign);
		}

		VectorComplexType tmp(dimNm1, ComplexType(0));
		for (SizeType k = 0; k < dimNm1; ++k)
			for (SizeType m = 0; m < dimN; ++m)
				tmp[k] += ComplexType(cFock(k, m)) * ComplexType(eigvecsN(m, 0));

		h_pre_.resize(dimNm1, ComplexType(0));
		for (SizeType l = 0; l < dimNm1; ++l)
			for (SizeType k = 0; k < dimNm1; ++k)
				h_pre_[l] += std::conj(ComplexType(eigvecsNm1(k, l))) * tmp[k];
	}

	// Build Phi_(k, n) and Psi_(k, n):
	//   Phi_(k,n) = sum_m Mcdag(k,m) b_m exp(i(E^{N+1}_k - E^N_m) t_n)
	//   Psi_(k,n) = sum_m Mc(k,m)   b_m exp(i(E^{N-1}_k - E^N_m) t_n)
	// where Mcdag and Mc are c†/c matrix elements in the post-quench eigenbasis.
	void buildPhiPsi(const BasisBaseType& basisN,
	                 const BasisBaseType& basisN1,
	                 const BasisBaseType& basisNm1,
	                 const MatrixType&    eigvecsN,
	                 const MatrixType&    eigvecsN1,
	                 const MatrixType&    eigvecsNm1,
	                 SizeType             impSite)
	{
		const SizeType dimN   = basisN.size();
		const SizeType dimN1  = basisN1.size();
		const SizeType dimNm1 = basisNm1.size();
		const SizeType nSteps = nT_ + 1;
		const SizeType spinUp = LanczosPlusPlus::LanczosGlobals::SPIN_UP;
		const SizeType spinDn = LanczosPlusPlus::LanczosGlobals::SPIN_DOWN;
		const LabeledOperatorType opCdag(LabeledOperatorType::Label::OPERATOR_CDAGGER);
		const LabeledOperatorType opC(LabeledOperatorType::Label::OPERATOR_C);

		// -- Particle channel (c†) -------------------------------------------
		// Mcdag(k,m) = eigvecsN1† * cdagFock * eigvecsN
		MatrixComplexType Mcdag(dimN1, dimN);
		{
			MatrixType cdagFock(dimN1, dimN);
			for (SizeType m = 0; m < dimN; ++m) {
				WordType    ket1 = basisN(m, spinUp);
				WordType    ket2 = basisN(m, spinDn);
				PairIntType bra  = basisN1.getBraIndex(ket1, ket2, opCdag, impSite, spinUp, 0);
				if (bra.first < 0)
					continue;
				int sign = basisN.doSignGf(ket1, ket2, impSite, spinUp, 0);
				cdagFock(static_cast<SizeType>(bra.first), m) = ComplexOrRealType(sign);
			}
			// tmp(k, m) = sum_a conj(eigvecsN1(a,k)) * cdagFock(a,m)
			MatrixComplexType tmp(dimN1, dimN);
			for (SizeType k = 0; k < dimN1; ++k)
				for (SizeType m = 0; m < dimN; ++m)
					for (SizeType a = 0; a < dimN1; ++a)
						tmp(k, m) += std::conj(ComplexType(eigvecsN1(a, k)))
						           * ComplexType(cdagFock(a, m));
			// Mcdag(k,m) = sum_b tmp(k,b) * eigvecsN(b,m)
			for (SizeType k = 0; k < dimN1; ++k)
				for (SizeType m = 0; m < dimN; ++m)
					for (SizeType b = 0; b < dimN; ++b)
						Mcdag(k, m) += tmp(k, b) * ComplexType(eigvecsN(b, m));
		}

		// Phi_(k, n) = sum_m Mcdag(k,m) * b_m * exp(i * omega_km * t_n)
		Phi_.resize(dimN1, nSteps);
		for (SizeType k = 0; k < dimN1; ++k) {
			for (SizeType n = 0; n < nSteps; ++n) {
				const RealType tn = n * params_.dt;
				ComplexType    s  = 0;
				for (SizeType m = 0; m < dimN; ++m) {
					const RealType phase = (energiesN1_post_[k] - energiesN_post_[m]) * tn;
					s += Mcdag(k, m) * b_[m]
					   * ComplexType(std::cos(phase), std::sin(phase));
				}
				Phi_(k, n) = s;
			}
		}

		// -- Hole channel (c) ------------------------------------------------
		// Mc(k,m) = eigvecsNm1† * cFock * eigvecsN
		MatrixComplexType Mc(dimNm1, dimN);
		{
			MatrixType cFock(dimNm1, dimN);
			for (SizeType m = 0; m < dimN; ++m) {
				WordType    ket1 = basisN(m, spinUp);
				WordType    ket2 = basisN(m, spinDn);
				PairIntType bra  = basisNm1.getBraIndex(ket1, ket2, opC, impSite, spinUp, 0);
				if (bra.first < 0)
					continue;
				int sign = basisN.doSignGf(ket1, ket2, impSite, spinUp, 0);
				cFock(static_cast<SizeType>(bra.first), m) = ComplexOrRealType(sign);
			}
			MatrixComplexType tmp(dimNm1, dimN);
			for (SizeType k = 0; k < dimNm1; ++k)
				for (SizeType m = 0; m < dimN; ++m)
					for (SizeType a = 0; a < dimNm1; ++a)
						tmp(k, m) += std::conj(ComplexType(eigvecsNm1(a, k)))
						           * ComplexType(cFock(a, m));
			for (SizeType k = 0; k < dimNm1; ++k)
				for (SizeType m = 0; m < dimN; ++m)
					for (SizeType b = 0; b < dimN; ++b)
						Mc(k, m) += tmp(k, b) * ComplexType(eigvecsN(b, m));
		}

		// Psi_(k, n) = sum_m Mc(k,m) * b_m * exp(i * omega_km * t_n)
		Psi_.resize(dimNm1, nSteps);
		for (SizeType k = 0; k < dimNm1; ++k) {
			for (SizeType n = 0; n < nSteps; ++n) {
				const RealType tn = n * params_.dt;
				ComplexType    s  = 0;
				for (SizeType m = 0; m < dimN; ++m) {
					const RealType phase = (energiesNm1_post_[k] - energiesN_post_[m]) * tn;
					s += Mc(k, m) * b_[m]
					   * ComplexType(std::cos(phase), std::sin(phase));
				}
				Psi_(k, n) = s;
			}
		}
	}

	// O_N1_(k, l) = <k^{N+1}_post | l^{N+1}_pre>
	void buildON1(const MatrixType& eigvecsPreN1, const MatrixType& eigvecsPostN1)
	{
		const SizeType dim = eigvecsPostN1.rows();
		assert(eigvecsPreN1.rows() == dim);
		O_N1_.resize(dim, dim);
		for (SizeType k = 0; k < dim; ++k)
			for (SizeType l = 0; l < dim; ++l) {
				ComplexType s = 0;
				for (SizeType i = 0; i < dim; ++i)
					s += std::conj(ComplexType(eigvecsPostN1(i, k)))
					   * ComplexType(eigvecsPreN1(i, l));
				O_N1_(k, l) = s;
			}
	}

	// chi_(k, j) = sum_l O_N1(k,l) * f_l * exp(-Omega_l * (beta - tau_j))
	// G^{Left} involves only the particle sector; no hole term needed.
	// Omega_l = E^{N+1}_{l,pre} - E^N_{0,pre}
	void buildChi()
	{
		const SizeType dimN1post = O_N1_.rows();
		const SizeType dimN1pre  = energiesN1_pre_.size();
		const SizeType nTauSteps = nTau_ + 1;
		const RealType beta      = nTau_ * dtau_;

		chi_.resize(dimN1post, nTauSteps);
		for (SizeType k = 0; k < dimN1post; ++k) {
			for (SizeType j = 0; j < nTauSteps; ++j) {
				const RealType tau = j * dtau_;
				ComplexType    s   = 0;
				for (SizeType l = 0; l < dimN1pre; ++l) {
					const RealType OmegaL = energiesN1_pre_[l] - E0_pre_;
					// Use exp(-OmegaL*(beta-tau)) so chi decays from tau=beta
					// toward tau=0, matching G^{Left}(0,tau) = -i G^M(beta-tau).
					s += O_N1_(k, l) * f_[l] * std::exp(-OmegaL * (beta - tau));
				}
				chi_(k, j) = s;
			}
		}
	}

	// Populate gimp_.matsubara_t and gimp_.matsubara_w from the pre-quench
	// Lehmann representation including both particle (N+1) and hole (N-1) sectors.
	//
	// G^M(tau) = -sum_l |f_l|^2 exp(-Omega_l^p tau)         [particle, T=0]
	//            -sum_k |h_k|^2 exp(-Omega_k^h (beta - tau)) [hole, finite beta]
	// G^M(iw)  =  sum_l |f_l|^2 / (iw - Omega_l^p)
	//           + sum_k |h_k|^2 / (iw + Omega_k^h)
	// where Omega^p = E^{N+1} - E^N_0 > 0 and Omega^h = E^N_0 - E^{N-1} > 0.
	// Both sectors are needed for a PH-symmetric purely imaginary G^M at mu=0.
	void buildMatsubara()
	{
		const SizeType dimN1pre  = energiesN1_pre_.size();
		const SizeType dimNm1pre = energiesNm1_pre_.size();
		const SizeType nTauSteps = nTau_ + 1;
		const SizeType nFreqs    = nTau_;
		const RealType beta      = params_.eqParams.ficticiousBeta;

		for (SizeType j = 0; j < nTauSteps; ++j) {
			const RealType tau = j * dtau_;
			ComplexType    gm  = 0;
			for (SizeType l = 0; l < dimN1pre; ++l) {
				const RealType OmegaL = energiesN1_pre_[l] - E0_pre_;
				gm -= std::norm(f_[l]) * std::exp(-OmegaL * tau);
			}
			for (SizeType l = 0; l < dimNm1pre; ++l) {
				const RealType OmegaH = E0_pre_ - energiesNm1_pre_[l];
				if (OmegaH <= RealType(0))
					continue; // E_k^{N-1} >= E_0^N: unreachable from GS, skip
				gm -= std::norm(h_pre_[l]) * std::exp(-OmegaH * (beta - tau));
			}
			gimp_.matsubara_t[j] = gm;
		}

		for (SizeType k = 0; k < nFreqs; ++k) {
			const int    kInt = static_cast<int>(k);
			const int    nInt = static_cast<int>(nFreqs);
			const RealType omegaK =
			    RealType(2 * kInt - nInt + 1) * M_PI / beta;
			ComplexType gw = 0;
			for (SizeType l = 0; l < dimN1pre; ++l) {
				const RealType OmegaL = energiesN1_pre_[l] - E0_pre_;
				gw += RealType(std::norm(f_[l]))
				    / (ComplexType(0, omegaK) - OmegaL);
			}
			for (SizeType l = 0; l < dimNm1pre; ++l) {
				const RealType OmegaH = E0_pre_ - energiesNm1_pre_[l];
				if (OmegaH <= RealType(0))
					continue; // E_k^{N-1} >= E_0^N: unreachable from GS, skip
				gw += RealType(std::norm(h_pre_[l]))
				    / (ComplexType(0, omegaK) + OmegaH);
			}
			gimp_.matsubara_w[k] = gw;
		}
	}

	// G^>(t_n, t_j) = -i sum_k conj(Phi_k[n]) * Phi_k[j]
	ComplexType gGreater(int n, int j) const
	{
		const SizeType dimN1 = Phi_.rows();
		ComplexType    s     = 0;
		for (SizeType k = 0; k < dimN1; ++k)
			s += std::conj(Phi_(k, static_cast<SizeType>(n)))
			   * Phi_(k, static_cast<SizeType>(j));
		return ComplexType(0, -1) * s;
	}

	// G^<(t_n, t_j) = i sum_k conj(Psi_k[j]) * Psi_k[n]
	ComplexType gLesser(int n, int j) const
	{
		const SizeType dimNm1 = Psi_.rows();
		ComplexType    s      = 0;
		for (SizeType k = 0; k < dimNm1; ++k)
			s += std::conj(Psi_(k, static_cast<SizeType>(j)))
			   * Psi_(k, static_cast<SizeType>(n));
		return ComplexType(0, 1) * s;
	}

	// G^{Left}(t_n, tau_j) = +i sum_k chi_k[j] * conj(Phi_k[n])
	// Sign is +i (not -i) because chi uses exp(-OmegaL*(beta-tau)); together
	// they satisfy G^{Left}(0,tau) = -i G^M(beta-tau).
	ComplexType gLeft(int n, SizeType j) const
	{
		const SizeType dimN1 = chi_.rows();
		ComplexType    s     = 0;
		for (SizeType k = 0; k < dimN1; ++k)
			s += chi_(k, j) * std::conj(Phi_(k, static_cast<SizeType>(n)));
		return ComplexType(0, 1) * s;
	}

	// ---- Member variables --------------------------------------------------
	const ParamsNeqType&            params_;
	typename InputNgType::Readable& io_;
	SizeType                        nup_;
	SizeType                        ndown_;
	SizeType                        nT_;
	SizeType                        nTau_;
	RealType                        dtau_;
	KBType                          gimp_;

	// Pre-quench N+1 spectrum and operator amplitudes
	VectorRealType    energiesN1_pre_;
	RealType          E0_pre_;
	VectorComplexType f_;      // f_l = <l^{N+1}_pre | c†_imp | GS_pre>

	// Pre-quench N-1 spectrum and hole amplitudes (for complete G^M)
	VectorRealType    energiesNm1_pre_;
	VectorComplexType h_pre_;  // h_k = <k^{N-1}_pre | c_imp | GS_pre>

	// Post-quench spectra
	VectorRealType energiesN_post_, energiesN1_post_, energiesNm1_post_;

	// Quench overlaps: b_n = <n^N_post | GS_pre>
	VectorComplexType b_;

	// Overlap matrix between pre and post N+1 sectors
	MatrixComplexType O_N1_;

	// Precomputed time-dependent amplitudes (rows=sector states, cols=time steps)
	MatrixComplexType Phi_; // particle channel
	MatrixComplexType Psi_; // hole channel
	MatrixComplexType chi_; // imaginary-time factor for G^{Left}
};

} // namespace Dmft
#endif // IMPURITYSOLVER_NEQ_EXACTDIAG_H
