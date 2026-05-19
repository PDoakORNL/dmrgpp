#ifndef IMPURITYSOLVER_EQ_LANCZOS_H
#define IMPURITYSOLVER_EQ_LANCZOS_H

#include "CincuentaInputCheck.h"
#include "ImpuritySolverBase.h"
#include "InputCheck.h"
#include "LanczosPlusPlus/src/Engine/LabeledOperator.h"
#include "LanczosPlusPlus/src/Engine/LanczosGlobals.h"
#include "LanczosPlusPlus/src/Engine/ModelSelector.h"
#include "LanczosImpurityUtils.h"
#include "Matrix.h"
#include "ParamsDmftSolver.h"
#include "PsimagLite.h"
#include "Vector.h"
#include <complex>
#include <memory>

namespace Dmft {

// Equilibrium impurity solver using truncated Lanczos diagonalization.
//
// Replaces ImpuritySolverDmrg for the equilibrium DMFT self-consistency loop.
// Rather than spawning one DMRG correction-vector sub-process per frequency,
// this solver diagonalizes H once per DMFT iteration with LanczosPlusPlus
// and evaluates G(z) = sum_l |f_l|^2 / (z - Omega_l) analytically at all
// Matsubara or real frequencies from the resulting Lehmann poles.
//
// Accuracy is controlled by NstatesEq (Lanczos states kept per sector).
// Setting NstatesEq=0 falls back to full diagonalization (exact for small systems).
//
// Required input parameters (beyond the standard equilibrium set):
//   TargetElectronsUp=<int>
//   TargetElectronsDown=<int>
//   NstatesEq=<int>   (0 = full diag; >0 = truncated Lanczos)
template <typename ComplexOrRealType>
class ImpuritySolverEqLanczos : public ImpuritySolverBase<ComplexOrRealType> {

public:

	using BaseType             = ImpuritySolverBase<ComplexOrRealType>;
	using RealType             = typename BaseType::RealType;
	using ComplexType          = std::complex<RealType>;
	using VectorRealType       = typename BaseType::VectorRealType;
	using VectorComplexType    = typename BaseType::VectorComplexType;
	using ApplicationType      = typename BaseType::ApplicationType;
	using InputNgType          = typename BaseType::InputNgType;
	using ParamsDmftSolverType = ParamsDmftSolver<ComplexOrRealType>;
	using MatrixType           = PsimagLite::Matrix<ComplexOrRealType>;

	using Utils               = LanczosImpurityUtils<ComplexOrRealType>;
	using DmrgInputReadable   = typename Utils::DmrgInputReadable;
	using GeometryType        = typename Utils::GeometryType;
	using ModelSelectorType =
	    LanczosPlusPlus::ModelSelector<ComplexOrRealType, GeometryType, DmrgInputReadable>;
	using ModelBaseType   = typename Utils::ModelBaseType;
	using BasisBaseType   = typename Utils::BasisBaseType;
	using LabeledOperatorType = LanczosPlusPlus::LabeledOperator;
	using WordType            = LanczosPlusPlus::LanczosGlobals::WordType;
	using PairIntType         = LanczosPlusPlus::LanczosGlobals::PairIntType;

	ImpuritySolverEqLanczos(const ParamsDmftSolverType&     params,
	                        const ApplicationType&,
	                        typename InputNgType::Readable& io)
	    : BaseType(params.ficticiousBeta, params.nMatsubaras, io)
	    , io_(io)
	    , nup_(0)
	    , ndown_(0)
	    , nStates_(0)
	    , eta_(0.1)
	    , freq_enum_(PsimagLite::FreqEnum::MATSUBARA)
	{
		io.readline(nup_,     "TargetElectronsUp=");
		io.readline(ndown_,   "TargetElectronsDown=");
		io.readline(nStates_, "NstatesEq=");
		try {
			io.readline(eta_, "OmegaDelta=");
		} catch (std::exception&) {}
	}

	// Diagonalize H(U) with the current bath parameters, build the Lehmann
	// representation, and fill gimp_ at all Matsubara or real frequencies.
	void solve(const VectorRealType& bathParams,
	           PsimagLite::FreqEnum  freq_enum,
	           SizeType) override
	{
		const SizeType nBath  = bathParams.size() / 2;
		const SizeType nsites = nBath + 1;
		const SizeType impSite = 0;

		VectorRealType hoppings(nBath), bathEps(nBath);
		for (SizeType i = 0; i < nBath; ++i) {
			hoppings[i] = bathParams[i];
			bathEps[i]  = bathParams[nBath + i];
		}

		RealType U = 0;
		io_.readline(U, "HubbardU=");

		VectorRealType potV(nsites, 0);
		potV[impSite] = -RealType(0.5) * U;
		for (SizeType i = 0; i < nBath; ++i)
			potV[i + 1] = bathEps[i];

		const std::string inputStr =
		    Utils::buildLanczosInput(U, nup_, ndown_, hoppings, potV, nsites);

		Dmrg::InputCheck ic;
		typename PsimagLite::InputNg<Dmrg::InputCheck>::Writeable ioW(ic, inputStr);
		DmrgInputReadable ioR(ioW);
		GeometryType      geom(ioR);
		ModelSelectorType ms(ioR, geom);
		const ModelBaseType& model = ms();

		VectorRealType energiesN;
		MatrixType     eigvecsN;
		Utils::diagWithBasisTruncated(model, model.basis(), geom, nStates_,
		                              energiesN, eigvecsN);
		const RealType E0 = energiesN[0];

		VectorRealType energiesN1;
		MatrixType     eigvecsN1;
		std::unique_ptr<BasisBaseType> bN1(model.createBasis(nup_ + 1, ndown_));
		Utils::diagWithBasisTruncated(model, *bN1, geom, nStates_,
		                              energiesN1, eigvecsN1);

		VectorComplexType f;
		buildF(model.basis(), *bN1, eigvecsN, eigvecsN1, impSite, f);

		VectorRealType energiesNm1;
		MatrixType     eigvecsNm1;
		std::unique_ptr<BasisBaseType> bNm1(model.createBasis(nup_ - 1, ndown_));
		Utils::diagWithBasisTruncated(model, *bNm1, geom, nStates_,
		                              energiesNm1, eigvecsNm1);

		VectorComplexType h;
		buildH(model.basis(), *bNm1, eigvecsN, eigvecsNm1, impSite, h);

		fillGimp(f, energiesN1, h, energiesNm1, E0, freq_enum);
		freq_enum_ = freq_enum;
	}

	const VectorComplexType& gimp()     const override { return gimp_; }
	PsimagLite::FreqEnum     freqEnum() const override { return freq_enum_; }

private:

	// f[l] = <l^{N+1} | c†_{imp,up} | GS_N>
	void buildF(const BasisBaseType& basisN,
	            const BasisBaseType& basisN1,
	            const MatrixType&    eigvecsN,   // (dimN,  1) — only GS needed
	            const MatrixType&    eigvecsN1,  // (dimN1, nKept)
	            SizeType             impSite,
	            VectorComplexType&   f)
	{
		const SizeType dimN_full  = basisN.size();
		const SizeType dimN1_full = basisN1.size();
		const SizeType nKeptN1    = eigvecsN1.cols();

		const LabeledOperatorType opCdag(LabeledOperatorType::Label::OPERATOR_CDAGGER);
		const SizeType spinUp = LanczosPlusPlus::LanczosGlobals::SPIN_UP;

		// cdagFock(k, m) = <k^{N+1}| c†_{imp,up} |m^N>  in Fock basis
		MatrixType cdagFock(dimN1_full, dimN_full);
		for (SizeType m = 0; m < dimN_full; ++m) {
			WordType    ket1 = basisN(m, spinUp);
			WordType    ket2 = basisN(m, LanczosPlusPlus::LanczosGlobals::SPIN_DOWN);
			PairIntType bra  = basisN1.getBraIndex(ket1, ket2, opCdag, impSite, spinUp, 0);
			if (bra.first < 0) continue;
			int sign = basisN.doSignGf(ket1, ket2, impSite, spinUp, 0);
			cdagFock(static_cast<SizeType>(bra.first), m) = ComplexOrRealType(sign);
		}

		// tmp[k] = sum_m cdagFock(k,m) * GS[m]
		VectorComplexType tmp(dimN1_full, ComplexType(0));
		for (SizeType k = 0; k < dimN1_full; ++k)
			for (SizeType m = 0; m < dimN_full; ++m)
				tmp[k] += ComplexType(cdagFock(k, m)) * ComplexType(eigvecsN(m, 0));

		// f[l] = sum_k conj(eigvecsN1(k, l)) * tmp[k]
		f.assign(nKeptN1, ComplexType(0));
		for (SizeType l = 0; l < nKeptN1; ++l)
			for (SizeType k = 0; k < dimN1_full; ++k)
				f[l] += std::conj(ComplexType(eigvecsN1(k, l))) * tmp[k];
	}

	// h[k] = <k^{N-1} | c_{imp,up} | GS_N>
	void buildH(const BasisBaseType& basisN,
	            const BasisBaseType& basisNm1,
	            const MatrixType&    eigvecsN,
	            const MatrixType&    eigvecsNm1,
	            SizeType             impSite,
	            VectorComplexType&   h)
	{
		const SizeType dimN_full   = basisN.size();
		const SizeType dimNm1_full = basisNm1.size();
		const SizeType nKeptNm1    = eigvecsNm1.cols();

		const LabeledOperatorType opC(LabeledOperatorType::Label::OPERATOR_C);
		const SizeType spinUp = LanczosPlusPlus::LanczosGlobals::SPIN_UP;

		// cFock(k, m) = <k^{N-1}| c_{imp,up} |m^N>  in Fock basis
		MatrixType cFock(dimNm1_full, dimN_full);
		for (SizeType m = 0; m < dimN_full; ++m) {
			WordType    ket1 = basisN(m, spinUp);
			WordType    ket2 = basisN(m, LanczosPlusPlus::LanczosGlobals::SPIN_DOWN);
			PairIntType bra  = basisNm1.getBraIndex(ket1, ket2, opC, impSite, spinUp, 0);
			if (bra.first < 0) continue;
			int sign = basisN.doSignGf(ket1, ket2, impSite, spinUp, 0);
			cFock(static_cast<SizeType>(bra.first), m) = ComplexOrRealType(sign);
		}

		// tmp[k] = sum_m cFock(k,m) * GS[m]
		VectorComplexType tmp(dimNm1_full, ComplexType(0));
		for (SizeType k = 0; k < dimNm1_full; ++k)
			for (SizeType m = 0; m < dimN_full; ++m)
				tmp[k] += ComplexType(cFock(k, m)) * ComplexType(eigvecsN(m, 0));

		// h[l] = sum_k conj(eigvecsNm1(k, l)) * tmp[k]
		h.assign(nKeptNm1, ComplexType(0));
		for (SizeType l = 0; l < nKeptNm1; ++l)
			for (SizeType k = 0; k < dimNm1_full; ++k)
				h[l] += std::conj(ComplexType(eigvecsNm1(k, l))) * tmp[k];
	}

	// Fill gimp_ from the Lehmann representation, scaled to match LanczosPlusPlus convention.
	//
	// Physical formula:
	//   G(z) = sum_l |f_l|^2 / (z - Omega_l^p)   [particle, N+1 sector]
	//        + sum_k |h_k|^2 / (z - Omega_k^h)   [hole,     N-1 sector]
	// where Omega_l^p = E_{N+1,l} - E0 > 0 (UHB) and Omega_k^h = E0 - E_{N-1,k} < 0 (LHB).
	//
	// The 4x prefactor matches the LanczosPlusPlus diagonal-GF convention used by
	// ImpuritySolverExactDiag: getModifiedState accumulates c†_isite + c†_jsite for
	// isite==jsite, doubling the vector and giving weight = ||2c†|GS>||^2 = 4*||c†|GS>||^2
	// per sector instead of the physical ||c†|GS>||^2.  Both the DMFT Weiss-field fit and
	// the self-energy update 1/G_imp are calibrated for this 4x-scaled convention.
	void fillGimp(const VectorComplexType& f,
	              const VectorRealType&    energiesN1,
	              const VectorComplexType& h,
	              const VectorRealType&    energiesNm1,
	              RealType                 E0,
	              PsimagLite::FreqEnum     freq_enum)
	{
		const SizeType nParticle = f.size();
		const SizeType nHole     = h.size();
		const RealType scale     = 4.0;

		if (freq_enum == PsimagLite::FreqEnum::MATSUBARA) {
			const SizeType n = this->matsubaras().total();
			gimp_.resize(n);
			for (SizeType i = 0; i < n; ++i) {
				const RealType wn = this->matsubaras().omega(i);
				ComplexType    gw(0);
				for (SizeType l = 0; l < nParticle; ++l) {
					const RealType OmegaP = energiesN1[l] - E0;
					gw += RealType(std::norm(f[l])) / (ComplexType(0, wn) - OmegaP);
				}
				for (SizeType k = 0; k < nHole; ++k) {
					const RealType OmegaH = E0 - energiesNm1[k];
					gw += RealType(std::norm(h[k])) / (ComplexType(0, wn) - OmegaH);
				}
				gimp_[i] = scale * gw;
			}
		} else {
			const SizeType n = this->realFreqRange().total();
			gimp_.resize(n);
			for (SizeType i = 0; i < n; ++i) {
				const RealType omega = this->realFreqRange().omega(i);
				ComplexType    gw(0);
				for (SizeType l = 0; l < nParticle; ++l) {
					const RealType OmegaP = energiesN1[l] - E0;
					gw += RealType(std::norm(f[l])) / (ComplexType(omega, eta_) - OmegaP);
				}
				for (SizeType k = 0; k < nHole; ++k) {
					const RealType OmegaH = E0 - energiesNm1[k];
					gw += RealType(std::norm(h[k])) / (ComplexType(omega, eta_) - OmegaH);
				}
				gimp_[i] = scale * gw;
			}
		}
	}

	typename InputNgType::Readable& io_;
	SizeType                        nup_;
	SizeType                        ndown_;
	SizeType                        nStates_;
	RealType                        eta_;
	VectorComplexType               gimp_;
	PsimagLite::FreqEnum            freq_enum_;
};

} // namespace Dmft
#endif // IMPURITYSOLVER_EQ_LANCZOS_H
