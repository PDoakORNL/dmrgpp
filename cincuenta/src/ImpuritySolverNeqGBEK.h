#ifndef IMPURITYSOLVER_NEQ_GBEK_H
#define IMPURITYSOLVER_NEQ_GBEK_H

#include "CincuentaInputCheck.h"
#include "ImpuritySolverNeqBase.h"
#include "ImpuritySolverNeqExactDiag.h"
#include "KadanoffBaym.h"
#include "LanczosPlusPlus/src/Engine/DefaultSymmetry.h"
#include "LanczosPlusPlus/src/Engine/InputCheck.h"
#include "LanczosPlusPlus/src/Engine/InternalProductStored.h"
#include "LanczosPlusPlus/src/Engine/LabeledOperator.h"
#include "LanczosPlusPlus/src/Engine/LanczosGlobals.h"
#include "LanczosPlusPlus/src/Engine/ModelSelector.h"
#include "LanczosSolver.h"
#include "Matrix.h"
#include "NeqBathDecomposition.h"
#include "ParametersForSolver.h"
#include "ParamsNeqDmftSolver.h"
#include "PsimagLite.h"
#include "Vector.h"

#include "CrsMatrix.h"

#include <algorithm>
#include <cassert>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

// Non-equilibrium impurity solver implementing the GBEK two-bath scheme:
//   Gramsch, Balzer, Eckstein, Kollar — PRB 88, 235106 (2013).
//
// The hybridization Δ is split as Δ = Δ⁻ + Δ⁺:
//   Δ⁻ : first bath — equilibrium memory from fixed {V_α, ε_α}
//   Δ⁺ : second bath — neq dynamics via rank-L Cholesky of i·Δ⁺_<
//
// For L=0 (first bath only), delegates entirely to ImpuritySolverNeqExactDiag
// (exact Lehmann representation), which is the correct L=0 limit.
//
// For L>0, the extended Fock space includes 2L second-bath sites:
//   L "empty" sites (initially unoccupied, ε=0) coupled to impurity via V^+_{n,p}
//   L "occupied" sites (initially doubly occupied, ε=0) coupled to impurity via V^+_{n,p}
// The many-body state is propagated step by step under the time-dependent H(t_n)
// using a Krylov Lanczos matrix exponential.
//
// Green's functions are computed as inner products (GBEK Eqs. 72-74):
//   G^<(t_j, t_n) = +i <Ψ(t_j) | Ψ(t_n)>   [N-1 sector]
//   G^>(t_n, t_j) = -i <Φ(t_n) | Φ(t_j)>   [N+1 sector]
// G^{Left}(t_n, τ) is reused from ImpuritySolverNeqExactDiag (Matsubara-based).
//
// Gα/Gβ spin-seed averaging (GBEK Eq. 70): whenever the base filling is
// spin-imbalanced (nup_ != ndown_ -- e.g. the NeqAtomicLimit single-atom
// seed nup=1, ndown=0), a single extended-Fock-space ground state is
// spin-polarized and does not represent the paramagnetic, particle-hole
// symmetric physics on its own. The paper's remedy: run the WHOLE
// extended-system calculation twice -- once seeded as if the impurity's
// extra electron were up (system alpha, nup_ext=nup+L, ndown_ext=ndown+L),
// once as if it were down (system beta, nup_ext=ndown+L, ndown_ext=nup+L)
// -- and average the resulting (always up-channel-probed) Green's
// functions. When nup_==ndown_ (e.g. the near-atomic path's half-filled
// multi-site bath), alpha and beta are the identical configuration, so
// this reduces to a no-op and only one calculation is actually run.
namespace Dmft {

template <typename ComplexOrRealType>
class ImpuritySolverNeqGBEK : public ImpuritySolverNeqBase<ComplexOrRealType> {

public:

	using BaseType          = ImpuritySolverNeqBase<ComplexOrRealType>;
	using RealType          = typename BaseType::RealType;
	using ComplexType       = typename BaseType::ComplexType;
	using VectorRealType    = typename BaseType::VectorRealType;
	using KBType            = typename BaseType::KBType;
	using InputNgType       = typename BaseType::InputNgType;
	using ParamsNeqType     = ParamsNeqDmftSolver<ComplexOrRealType>;
	using ExactDiagType     = ImpuritySolverNeqExactDiag<ComplexOrRealType>;
	using DecompType        = NeqBathDecomposition<ComplexOrRealType>;
	using VectorComplexType = typename PsimagLite::Vector<ComplexType>::Type;
	using MatrixComplexType = PsimagLite::Matrix<ComplexType>;
	using MatrixType        = PsimagLite::Matrix<ComplexOrRealType>;

	// LanczosPlusPlus / PsimagLite types (parallel to ImpuritySolverNeqExactDiag)
	using LppInputReadable =
	    typename PsimagLite::InputNg<LanczosPlusPlus::InputCheck>::Readable;
	using GeometryType = PsimagLite::
	    Geometry<ComplexOrRealType, LppInputReadable, LanczosPlusPlus::LanczosGlobals>;
	using ModelSelectorType
	    = LanczosPlusPlus::ModelSelector<ComplexOrRealType, GeometryType, LppInputReadable>;
	using ModelBaseType
	    = LanczosPlusPlus::LanczosModelBase<ComplexOrRealType, GeometryType, LppInputReadable>;
	using BasisBaseType       = typename ModelBaseType::BasisBaseType;
	using DefaultSymmetryType = LanczosPlusPlus::DefaultSymmetry<GeometryType, BasisBaseType>;
	using InternalProductStoredType
	    = LanczosPlusPlus::InternalProductStored<ModelBaseType, DefaultSymmetryType>;
	using LabeledOperatorType    = LanczosPlusPlus::LabeledOperator;
	using WordType               = LanczosPlusPlus::LanczosGlobals::WordType;
	using PairIntType            = LanczosPlusPlus::LanczosGlobals::PairIntType;
	using LanczosSolverForGSType = PsimagLite::LanczosSolver<InternalProductStoredType>;
	using CrsMatrixComplexType   = PsimagLite::CrsMatrix<ComplexType>;

	ImpuritySolverNeqGBEK(const ParamsNeqType& params, typename InputNgType::Readable& io)
	    : bathRank_(params.neqBathRank)
	    , params_(params)
	    , exactDiag_(params, io)
	    , nup_(0)
	    , ndown_(0)
	    , nBath_(0)
	    , nsites_ext_(0)
	{
		if (bathRank_ > 0) {
			io.readline(nup_, "TargetElectronsUp=");
			io.readline(ndown_, "TargetElectronsDown=");
		}
	}

	// Initialise the impurity solver and bath decomposition.
	// For L=0: delegates to ExactDiag (Lehmann representation).
	// For L>0: also builds the extended Fock space and seeds time-propagated states.
	void solve(const VectorRealType& bathParams) override
	{
		const RealType beta = params_.eqParams.ficticiousBeta;
		const RealType mu   = 0;

		decomp_ = std::make_unique<DecompType>(
		    bathRank_,
		    beta,
		    mu,
		    bathParams,
		    params_.nT,
		    params_.eqParams.nMatsubaras,
		    params_.dt,
		    params_.eqParams.ficticiousBeta
		        / static_cast<RealType>(params_.eqParams.nMatsubaras));

		exactDiag_.solve(bathParams);

		if (bathRank_ > 0)
			solveLplus(bathParams);
	}

	// Advance the Cholesky decomposition to step n and (for L>0) invalidate
	// the propagated state so the next computeGimp call re-propagates.
	void prepareTimeStep(int n, const KBType& delta) override
	{
		if (decomp_)
			decomp_->update(n, delta);
		if (bathRank_ > 0 && n > 0) {
			// V[n] just changed; invalidate cached PsiHist[n] and beyond, for
			// both spin configurations.
			sectorAlpha_.propagatedThrough
			    = std::min(sectorAlpha_.propagatedThrough, n - 1);
			if (!sameConfig_)
				sectorBeta_.propagatedThrough
				    = std::min(sectorBeta_.propagatedThrough, n - 1);
		}
	}

	// Fill G_imp KB components for all (t_n, t_j) with j <= n.
	// For L=0: delegates to ExactDiag (Lehmann).
	// For L>0: ExactDiag fills G^{Left}; GBEK inner products overwrite G^< and G^R.
	void computeGimp(KBType& gimp, int n) const override
	{
		// Always call ExactDiag for G^{Left} and Matsubara components
		exactDiag_.computeGimp(gimp, n);
		if (bathRank_ == 0)
			return;

		// Ensure states are propagated up to step n (lazy, idempotent)
		ensurePropagated(n);

		// Overwrite G^< and G^R with GBEK inner-product formulas
		for (int j = 0; j <= n; ++j) {
			const SizeType    sn    = static_cast<SizeType>(n);
			const SizeType    sj    = static_cast<SizeType>(j);
			const ComplexType gLess = gLesserGBEK(n, j);
			const ComplexType gGrtr = gGreaterGBEK(n, j);
			gimp.lesser(sn, sj)     = gLess;
			gimp.retarded(sn, sj)   = gGrtr - gLess;
		}
		// Anti-Hermitian: G^<(t_j, t_n) = -conj(G^<(t_n, t_j)) for j < n
		for (int j = 0; j < n; ++j)
			gimp.lesser(static_cast<SizeType>(j), static_cast<SizeType>(n))
			    = -std::conj(
			        gimp.lesser(static_cast<SizeType>(n), static_cast<SizeType>(j)));
	}

	// Matsubara components come from ExactDiag in both L=0 and L>0 cases.
	const KBType& gimp() const override { return exactDiag_.gimp(); }

	const DecompType* decomposition() const { return decomp_.get(); }

	void dumpPlusBath(const std::string& filename) const override
	{
		if (decomp_)
			decomp_->dumpPlusBath(filename);
	}

	// Raw Cholesky factor V_(n,p), for row-by-row comparison against an
	// independent offline trace (see NeqBathDecomposition::dumpV).
	void dumpV(const std::string& filename) const override
	{
		if (decomp_)
			decomp_->dumpV(filename);
	}

private:

	// Position and update info for one second-bath variable entry in the CSR matrix.
	struct VarEntry {
		int      nnzIdx; // index into CrsMatrix values array (matches setValues arg type)
		SizeType p; // second-bath rank index
		bool     isConj; // true → store conj(vMid[p])*sign, false → vMid[p]*sign
		int      sign; // Jordan-Wigner sign (±1)
	};

	// ========== Gα/Gβ per-configuration extended-Fock state ==========
	//
	// Everything that depends on WHICH spin the impurity's extra electron
	// (nup_ext,ndown_ext) is seeded with. System alpha uses
	// (nup+L, ndown+L) (impurity's extra electron is up); system beta
	// swaps to (ndown+L, nup+L) (as if it were down). Both are always
	// probed with the up-spin operator (see seedState/gLesserGBEKSector),
	// matching GBEK Eq. 70's Gα_σ/Gβ_σ for σ=up.
	struct ExtendedSector {
		SizeType                               dim1Nm1 = 0, dim1Np1 = 0;
		std::vector<WordType>                  upWordsNm1, dnWordsNm1;
		std::vector<WordType>                  upWordsNp1, dnWordsNp1;
		CrsMatrixComplexType                   csrNm1, csrNp1;
		std::vector<VarEntry>                  varNm1, varNp1;
		mutable std::vector<VectorComplexType> PsiHist, PhiHist;
		mutable int                            propagatedThrough = -1;
	};

	// ========== L>0 setup ==========

	// Build extended Fock space, diagonalise pre-quench H, seed PsiHist/PhiHist
	// for both spin configurations (alpha always; beta only if nup_ != ndown_,
	// since otherwise it is identical to alpha -- see class docstring).
	void solveLplus(const VectorRealType& bathParams)
	{
		const SizeType L = bathRank_;
		nBath_           = bathParams.size() / 2;
		nsites_ext_      = nBath_ + 1 + 2 * L;

		firstBathHop_.resize(nBath_);
		firstBathEps_.resize(nBath_);
		for (SizeType i = 0; i < nBath_; ++i) {
			firstBathHop_[i] = bathParams[i];
			firstBathEps_[i] = bathParams[nBath_ + i];
		}

		// Post-quench on-site potentials (ε=0 for second bath)
		potPost_.assign(nsites_ext_, RealType(0));
		potPost_[0] = -RealType(0.5) * params_.uFinal;
		for (SizeType i = 0; i < nBath_; ++i)
			potPost_[i + 1] = firstBathEps_[i];

		sameConfig_ = (nup_ == ndown_);

		buildSector(sectorAlpha_, nup_ + L, ndown_ + L, bathParams);
		if (!sameConfig_)
			buildSector(sectorBeta_, ndown_ + L, nup_ + L, bathParams);
	}

	// Build one spin configuration's extended Fock space, diagonalise its
	// pre-quench Hamiltonian, and seed PsiHist[0]/PhiHist[0] -- shared logic
	// for both system alpha and system beta (GBEK Eq. 70).
	void buildSector(ExtendedSector&       sector,
	                 SizeType              nupExt,
	                 SizeType              ndownExt,
	                 const VectorRealType& bathParams)
	{
		const SizeType L = bathRank_;

		// Extended hoppings (real; second bath starts at 0, updated later via Cholesky)
		VectorRealType hopExt(nBath_ + 2 * L, RealType(0));
		for (SizeType i = 0; i < nBath_; ++i)
			hopExt[i] = firstBathHop_[i];

		// Pre-quench extended potential: large ±ε for second bath to fix initial occupation
		const RealType bigEps = RealType(500);
		VectorRealType potPre(nsites_ext_, RealType(0));
		potPre[0] = -RealType(0.5) * params_.uInitial;
		for (SizeType i = 0; i < nBath_; ++i)
			potPre[i + 1] = firstBathEps_[i];
		for (SizeType p = 0; p < L; ++p)
			potPre[nBath_ + 1 + p] = +bigEps; // empty
		for (SizeType p = 0; p < L; ++p)
			potPre[nBath_ + 1 + L + p] = -bigEps; // occupied

		// Build extended pre-quench model
		const std::string inputPre = buildLanczosInput(
		    params_.uInitial, nupExt, ndownExt, hopExt, potPre, nsites_ext_);

		LanczosPlusPlus::InputCheck                                          ic;
		typename PsimagLite::InputNg<LanczosPlusPlus::InputCheck>::Writeable ioW(ic,
		                                                                         inputPre);
		LppInputReadable                                                     ioR(ioW);
		GeometryType                                                         geom(ioR);
		ModelSelectorType                                                    ms(ioR, geom);
		const ModelBaseType&                                                 model = ms();

		// Diagonalise N sector → pre-quench GS.
		// Full diag is feasible for L=1 (dim~4900); for L>=2 (dim~63k+) use Lanczos GS.
		VectorRealType energiesN;
		MatrixType     eigvecsN;
		if (nsites_ext_ <= 8) {
			diagWithBasis(model, model.basis(), geom, energiesN, eigvecsN);
		} else {
			eigvecsN = lanczosGS(model, geom);
		}

		// Build N-1 and N+1 sector bases for the extended system
		std::unique_ptr<BasisBaseType> bNm1(model.createBasis(nupExt - 1, ndownExt));
		std::unique_ptr<BasisBaseType> bNp1(model.createBasis(nupExt + 1, ndownExt));

		// Build fast-lookup tables from sorted Fock words
		buildFockLookup(*bNm1, sector.upWordsNm1, sector.dnWordsNm1, sector.dim1Nm1);
		buildFockLookup(*bNp1, sector.upWordsNp1, sector.dnWordsNp1, sector.dim1Np1);

		// Build CSR sparse matrices for H_ext and validate against applyHext
		buildHextCSR(sector.upWordsNm1,
		             sector.dnWordsNm1,
		             sector.dim1Nm1,
		             sector.csrNm1,
		             sector.varNm1);
		buildHextCSR(sector.upWordsNp1,
		             sector.dnWordsNp1,
		             sector.dim1Np1,
		             sector.csrNp1,
		             sector.varNp1);
		checkApplyHextEquivalence(sector);

		// Allocate history arrays
		const SizeType nSteps = params_.nT + 1;
		sector.PsiHist.assign(nSteps, VectorComplexType(bNm1->size(), ComplexType(0)));
		sector.PhiHist.assign(nSteps, VectorComplexType(bNp1->size(), ComplexType(0)));

		// Seed step 0: c_{imp,↑}|GS_pre_ext⟩ (N-1) and c†_{imp,↑}|GS_pre_ext⟩ (N+1)
		seedState(sector, eigvecsN, model.basis(), *bNm1, *bNp1);
		sector.propagatedThrough = 0;
	}

	// Apply c_{imp,↑} and c†_{imp,↑} to the pre-quench GS to seed PsiHist[0]/PhiHist[0].
	// Always the up-spin operator, regardless of whether this sector is
	// system alpha or system beta -- see class docstring.
	void seedState(ExtendedSector&      sector,
	               const MatrixType&    eigvecsN,
	               const BasisBaseType& basisN,
	               const BasisBaseType& basisNm1,
	               const BasisBaseType& basisNp1)
	{
		const SizeType            dimN   = basisN.size();
		const SizeType            spinUp = LanczosPlusPlus::LanczosGlobals::SPIN_UP;
		const SizeType            spinDn = LanczosPlusPlus::LanczosGlobals::SPIN_DOWN;
		const LabeledOperatorType opC(LabeledOperatorType::Label::OPERATOR_C);
		const LabeledOperatorType opCdag(LabeledOperatorType::Label::OPERATOR_CDAGGER);
		const SizeType            impSite = 0;

		for (SizeType m = 0; m < dimN; ++m) {
			const ComplexType amp(eigvecsN(m, 0));
			if (std::abs(amp) < RealType(1e-15))
				continue;

			const WordType ket1 = basisN(m, spinUp);
			const WordType ket2 = basisN(m, spinDn);

			// c_{imp,↑} |m⟩ → N-1 sector
			const PairIntType braC
			    = basisNm1.getBraIndex(ket1, ket2, opC, impSite, spinUp, 0);
			if (braC.first >= 0) {
				const int sign = basisN.doSignGf(ket1, ket2, impSite, spinUp, 0);
				sector.PsiHist[0][static_cast<SizeType>(braC.first)]
				    += ComplexType(sign) * amp;
			}

			// c†_{imp,↑} |m⟩ → N+1 sector
			const PairIntType braCdag
			    = basisNp1.getBraIndex(ket1, ket2, opCdag, impSite, spinUp, 0);
			if (braCdag.first >= 0) {
				const int sign = basisN.doSignGf(ket1, ket2, impSite, spinUp, 0);
				sector.PhiHist[0][static_cast<SizeType>(braCdag.first)]
				    += ComplexType(sign) * amp;
			}
		}
	}

	// Extract sorted up/dn Fock words from a basis for O(log n) index lookup.
	static void buildFockLookup(const BasisBaseType&   basis,
	                            std::vector<WordType>& upWords,
	                            std::vector<WordType>& dnWords,
	                            SizeType&              dim1)
	{
		const SizeType dim    = basis.size();
		const SizeType spinUp = LanczosPlusPlus::LanczosGlobals::SPIN_UP;
		const SizeType spinDn = LanczosPlusPlus::LanczosGlobals::SPIN_DOWN;

		if (dim == 0) {
			dim1 = 0;
			return;
		}

		// Determine dim1: period of the up-spin word (it cycles with period = #distinct up
		// words)
		const WordType firstUp = basis(0, spinUp);
		dim1                   = 1;
		while (dim1 < dim && basis(dim1, spinUp) != firstUp)
			++dim1;

		upWords.resize(dim1);
		for (SizeType x = 0; x < dim1; ++x)
			upWords[x] = basis(x, spinUp);

		const SizeType dim2 = dim / dim1;
		dnWords.resize(dim2);
		for (SizeType y = 0; y < dim2; ++y)
			dnWords[y] = basis(y * dim1, spinDn);
	}

	// ========== Precomputed CSR sparse Hamiltonian ==========

	// Build the CSR sparse matrix for H_ext in one sector.
	// Fixed entries (diagonal + first-bath hops) go directly into matrix values.
	// Variable entries (second-bath hops, time-dependent) are stored as zero initially;
	// their positions are recorded in varEntries for O(1) update at each time step.
	// Asserts that no (row,col) slot is simultaneously fixed and variable.
	void buildHextCSR(const std::vector<WordType>& upWords,
	                  const std::vector<WordType>& dnWords,
	                  SizeType                     dim1,
	                  CrsMatrixComplexType&        csr,
	                  std::vector<VarEntry>&       varEntries) const
	{
		const SizeType dim2   = dnWords.size();
		const SizeType dim    = dim1 * dim2;
		const SizeType spinUp = LanczosPlusPlus::LanczosGlobals::SPIN_UP;

		struct TripletEntry {
			SizeType    row, col;
			ComplexType fixedVal;
			bool        isVar;
			SizeType    p;
			bool        isConj;
			int         sign;
		};

		std::vector<TripletEntry> triplets;
		triplets.reserve(dim * (1 + 2 * (2 * nBath_ + 4 * bathRank_)));

		for (SizeType m = 0; m < dim; ++m) {
			const SizeType y    = m / dim1;
			const SizeType x    = m % dim1;
			const WordType up_m = upWords[x];
			const WordType dn_m = dnWords[y];

			// Diagonal term
			RealType diag = 0;
			if ((up_m & WordType(1)) && (dn_m & WordType(1)))
				diag += params_.uFinal;
			for (SizeType i = 0; i < nsites_ext_; ++i) {
				if ((up_m >> i) & WordType(1))
					diag += potPost_[i];
				if ((dn_m >> i) & WordType(1))
					diag += potPost_[i];
			}
			if (diag != RealType(0))
				triplets.push_back({ m, m, ComplexType(diag), false, 0, false, 1 });

			for (SizeType spin = 0; spin < 2; ++spin) {
				const WordType w       = (spin == spinUp) ? up_m : dn_m;
				const WordType w_other = (spin == spinUp) ? dn_m : up_m;

				auto tryHop = [&](SizeType    from,
				                  SizeType    to,
				                  bool        isVar,
				                  SizeType    p,
				                  bool        isConj,
				                  ComplexType fixedT)
				{
					if (!((w >> from) & WordType(1)))
						return;
					if ((w >> to) & WordType(1))
						return;
					const WordType new_w
					    = (w & ~(WordType(1) << from)) | (WordType(1) << to);
					const WordType new_up = (spin == spinUp) ? new_w : w_other;
					const WordType new_dn = (spin == spinUp) ? w_other : new_w;
					const SizeType m_new  = lookupIndex(
                                            new_up, new_dn, upWords, dnWords, dim1, dim2);
					if (m_new == std::numeric_limits<SizeType>::max())
						return;
					const int         sign = computeJWSign(w, from, to);
					const ComplexType val
					    = isVar ? ComplexType(0) : fixedT * ComplexType(sign);
					triplets.push_back(
					    { m_new, m, val, isVar, p, isConj, sign });
				};

				// First bath (fixed, real hoppings V_α)
				for (SizeType a = 0; a < nBath_; ++a) {
					const ComplexType Va(firstBathHop_[a]);
					tryHop(0, a + 1, false, 0, false, Va);
					tryHop(a + 1, 0, false, 0, false, Va);
				}

				// Second bath (variable: hop ↔ empty-site, ↔ occupied-site)
				for (SizeType p = 0; p < bathRank_; ++p) {
					const SizeType emptyS = 1 + nBath_ + p;
					const SizeType occS   = 1 + nBath_ + bathRank_ + p;
					tryHop(0, emptyS, true, p, false, ComplexType(0));
					tryHop(emptyS, 0, true, p, true, ComplexType(0));
					tryHop(0, occS, true, p, false, ComplexType(0));
					tryHop(occS, 0, true, p, true, ComplexType(0));
				}
			}
		}

		// Sort by (row, col) so we can build row-major CSR.
		std::sort(triplets.begin(),
		          triplets.end(),
		          [](const TripletEntry& a, const TripletEntry& b)
		          { return a.row < b.row || (a.row == b.row && a.col < b.col); });

		// Assert no fixed/var collision at the same (row,col) slot.
		for (SizeType i = 1; i < triplets.size(); ++i) {
			if (triplets[i].row == triplets[i - 1].row
			    && triplets[i].col == triplets[i - 1].col)
				assert(false
				       && "buildHextCSR: fixed and variable hop collide at same "
				          "(row,col)");
		}

		varEntries.clear();
		csr.resize(dim, dim);
		csr.reserve(triplets.size());

		int      nnzCounter = 0;
		SizeType tripletIdx = 0;
		for (SizeType row = 0; row < dim; ++row) {
			csr.setRow(row, nnzCounter);
			while (tripletIdx < triplets.size() && triplets[tripletIdx].row == row) {
				const TripletEntry& t = triplets[tripletIdx];
				csr.pushCol(t.col);
				csr.pushValue(t.fixedVal);
				if (t.isVar)
					varEntries.push_back({ nnzCounter, t.p, t.isConj, t.sign });
				++nnzCounter;
				++tripletIdx;
			}
		}
		csr.setRow(dim, nnzCounter);
	}

	// Overwrite second-bath (variable) CSR entries with current midpoint hoppings.
	// All other entries retain their fixed values from buildHextCSR.
	void updateCSR(CrsMatrixComplexType&           csr,
	               const std::vector<VarEntry>&    varEntries,
	               const std::vector<ComplexType>& vMid) const
	{
		for (const VarEntry& ve : varEntries) {
			const ComplexType Vp = ve.isConj ? std::conj(vMid[ve.p]) : vMid[ve.p];
			csr.setValues(ve.nnzIdx, Vp * ComplexType(ve.sign));
		}
	}

	// ========== Time propagation ==========

	// Ensure PsiHist_[0..n] and PhiHist_[0..n] are populated with the current V values.
	void ensurePropagated(int n) const
	{
		ensureSectorPropagated(sectorAlpha_, n);
		if (!sameConfig_)
			ensureSectorPropagated(sectorBeta_, n);
	}

	void ensureSectorPropagated(ExtendedSector& sector, int n) const
	{
		for (int k = sector.propagatedThrough + 1; k <= n; ++k)
			propagateOneStep(sector, k);
		if (n > sector.propagatedThrough)
			sector.propagatedThrough = n;
	}

	// Propagate step n from step n-1 using the midpoint Hamiltonian.
	// The midpoint second-bath hoppings V_mid come from the shared Cholesky
	// decomposition (decomp_) -- identical for both spin configurations,
	// since there is only one physical second bath.
	void propagateOneStep(ExtendedSector& sector, int n) const
	{
		assert(n > 0 && n < static_cast<int>(sector.PsiHist.size()));

		// Midpoint second-bath hoppings: V_mid = (V[n-1] + V[n]) / 2
		std::vector<ComplexType> vMid(bathRank_);
		for (SizeType p = 0; p < bathRank_; ++p) {
			const ComplexType vPrev = decomp_->Vplus(n - 1, static_cast<int>(p));
			const ComplexType vCurr = decomp_->Vplus(n, static_cast<int>(p));
			vMid[p]                 = RealType(0.5) * (vPrev + vCurr);
		}

		updateCSR(sector.csrNm1, sector.varNm1, vMid);
		updateCSR(sector.csrNp1, sector.varNp1, vMid);

		sector.PsiHist[n]
		    = krylovExpmvCSR(sector.PsiHist[n - 1], sector.csrNm1, params_.dt);
		sector.PhiHist[n]
		    = krylovExpmvCSR(sector.PhiHist[n - 1], sector.csrNp1, params_.dt);
	}

	// Krylov (Lanczos) approximation to exp(-i H dt) |psi⟩.
	// Uses the Lanczos recursion to project H into a small Krylov subspace,
	// then computes the matrix exponential there.
	VectorComplexType krylovExpmv(const VectorComplexType&        psi,
	                              const std::vector<ComplexType>& vMid,
	                              const std::vector<WordType>&    upWords,
	                              const std::vector<WordType>&    dnWords,
	                              SizeType                        dim1,
	                              RealType                        dt,
	                              int                             m = 40) const
	{
		const SizeType dim   = psi.size();
		const RealType norm0 = std::sqrt(std::real(innerProduct(psi, psi)));
		if (norm0 < RealType(1e-14))
			return psi;

		// Lanczos iteration to build K_m(H, psi/norm0)
		std::vector<VectorComplexType> Q;
		Q.reserve(m + 1);
		Q.push_back(psi);
		scaleVec(Q[0], ComplexType(RealType(1) / norm0));

		std::vector<RealType> alpha, beta;
		VectorComplexType     Hq(dim);
		int                   mUsed = 0;

		for (int j = 0; j < m; ++j) {
			applyHext(Q[j], Hq, vMid, upWords, dnWords, dim1);

			alpha.push_back(std::real(innerProduct(Q[j], Hq)));
			VectorComplexType w = Hq;
			axpy(w, -ComplexType(alpha[j]), Q[j]);
			if (j > 0)
				axpy(w, -ComplexType(beta[j - 1]), Q[j - 1]);

			const RealType b = std::sqrt(std::real(innerProduct(w, w)));
			mUsed            = j + 1;
			if (b < RealType(1e-10)) {
				std::cout << "  krylov early-exit at j=" << j << " / " << m << "\n";
				break; // invariant subspace
			}

			beta.push_back(b);
			Q.push_back(w);
			scaleVec(Q.back(), ComplexType(RealType(1) / b));
		}

		// Build complex tridiagonal T (real symmetric, stored as complex for
		// PsimagLite::diag)
		MatrixComplexType Tc(mUsed, mUsed, ComplexType(0));
		for (int j = 0; j < mUsed; ++j)
			Tc(j, j) = ComplexType(alpha[j]);
		for (int j = 0; j < mUsed - 1; ++j) {
			Tc(j + 1, j) = ComplexType(beta[j]);
			Tc(j, j + 1) = ComplexType(beta[j]);
		}

		// Diagonalise T_m: Tc gets overwritten with eigenvectors (columns)
		VectorRealType eigsTm;
		PsimagLite::diag(Tc, eigsTm, 'V');

		// Compute exp(-i T_m dt) e_1 in the eigenbasis and transform back
		// expCoeff[j] = Σ_k conj(Tc[0,k]) * exp(-i λ_k dt) * Tc[j,k]
		VectorComplexType expCoeff(mUsed, ComplexType(0));
		for (int k = 0; k < mUsed; ++k) {
			const ComplexType phase  = std::exp(ComplexType(0, -eigsTm[k] * dt));
			const ComplexType factor = std::conj(Tc(0, k)) * phase;
			for (int j = 0; j < mUsed; ++j)
				expCoeff[j] += factor * Tc(j, k);
		}

		// Result: norm0 * Σ_j expCoeff[j] * Q[j]
		VectorComplexType result(dim, ComplexType(0));
		for (int j = 0; j < mUsed; ++j) {
			const ComplexType c = norm0 * expCoeff[j];
			for (SizeType i = 0; i < dim; ++i)
				result[i] += c * Q[j][i];
		}
		return result;
	}

	// Krylov approximation using precomputed CSR sparse matrix.
	// Caller must call updateCSR before this to load current vMid values.
	// Identical algorithm to krylovExpmv; SpMV replaces applyHext.
	VectorComplexType krylovExpmvCSR(const VectorComplexType&    psi,
	                                 const CrsMatrixComplexType& csr,
	                                 RealType                    dt,
	                                 int                         m = 40) const
	{
		const SizeType dim   = psi.size();
		const RealType norm0 = std::sqrt(std::real(innerProduct(psi, psi)));
		if (norm0 < RealType(1e-14))
			return psi;

		std::vector<VectorComplexType> Q;
		Q.reserve(m + 1);
		Q.push_back(psi);
		scaleVec(Q[0], ComplexType(RealType(1) / norm0));

		std::vector<RealType> alpha, beta;
		VectorComplexType     Hq(dim, ComplexType(0));
		int                   mUsed = 0;

		for (int j = 0; j < m; ++j) {
			std::fill(Hq.begin(), Hq.end(), ComplexType(0));
			csr.matrixVectorProduct(Hq, Q[j]);

			alpha.push_back(std::real(innerProduct(Q[j], Hq)));
			VectorComplexType w = Hq;
			axpy(w, -ComplexType(alpha[j]), Q[j]);
			if (j > 0)
				axpy(w, -ComplexType(beta[j - 1]), Q[j - 1]);

			const RealType b = std::sqrt(std::real(innerProduct(w, w)));
			mUsed            = j + 1;
			if (b < RealType(1e-10)) {
				std::cout << "  krylov CSR early-exit at j=" << j << " / " << m
				          << "\n";
				break;
			}

			beta.push_back(b);
			Q.push_back(w);
			scaleVec(Q.back(), ComplexType(RealType(1) / b));
		}

		MatrixComplexType Tc(mUsed, mUsed, ComplexType(0));
		for (int j = 0; j < mUsed; ++j)
			Tc(j, j) = ComplexType(alpha[j]);
		for (int j = 0; j < mUsed - 1; ++j) {
			Tc(j + 1, j) = ComplexType(beta[j]);
			Tc(j, j + 1) = ComplexType(beta[j]);
		}

		VectorRealType eigsTm;
		PsimagLite::diag(Tc, eigsTm, 'V');

		VectorComplexType expCoeff(mUsed, ComplexType(0));
		for (int k = 0; k < mUsed; ++k) {
			const ComplexType phase  = std::exp(ComplexType(0, -eigsTm[k] * dt));
			const ComplexType factor = std::conj(Tc(0, k)) * phase;
			for (int j = 0; j < mUsed; ++j)
				expCoeff[j] += factor * Tc(j, k);
		}

		VectorComplexType result(dim, ComplexType(0));
		for (int j = 0; j < mUsed; ++j) {
			const ComplexType c = norm0 * expCoeff[j];
			for (SizeType i = 0; i < dim; ++i)
				result[i] += c * Q[j][i];
		}
		return result;
	}

	// One-shot equivalence check: assert max|applyHext(v) - CSR*v| < tol for both sectors.
	// Called once in solveLplus after buildHextCSR to validate the CSR construction.
	// Uses a non-zero artificial vMid to exercise the variable second-bath entries.
	void checkApplyHextEquivalence(ExtendedSector& sector) const
	{
		if (sector.upWordsNm1.empty() || sector.upWordsNp1.empty())
			return;

		// Artificial non-zero vMid so variable entries are non-trivially tested.
		std::vector<ComplexType> vMid(bathRank_);
		for (SizeType p = 0; p < bathRank_; ++p)
			vMid[p] = ComplexType(RealType(0.3) + RealType(0.1) * p,
			                      RealType(0.2) - RealType(0.05) * p);

		auto checkSector = [&](const std::string&           name,
		                       const std::vector<WordType>& upWords,
		                       const std::vector<WordType>& dnWords,
		                       SizeType                     dim1,
		                       CrsMatrixComplexType&        csr,
		                       const std::vector<VarEntry>& varEntries)
		{
			const SizeType dim = upWords.size() * dnWords.size();
			// Normalised unit vector as test input.
			VectorComplexType testV(
			    dim, ComplexType(RealType(1) / std::sqrt(RealType(dim))));

			updateCSR(csr, varEntries, vMid);

			VectorComplexType hvRef(dim, ComplexType(0));
			applyHext(testV, hvRef, vMid, upWords, dnWords, dim1);

			VectorComplexType hvCsr(dim, ComplexType(0));
			csr.matrixVectorProduct(hvCsr, testV);

			RealType maxDiff = 0;
			for (SizeType i = 0; i < dim; ++i)
				maxDiff = std::max(maxDiff, std::abs(hvRef[i] - hvCsr[i]));

			std::cout << "  checkApplyHextEquivalence [" << name << "]: dim=" << dim
			          << "  max|applyHext - CSR*v| = " << maxDiff << "\n";
			assert(maxDiff < RealType(1e-12)
			       && "applyHextCSR disagrees with applyHext — CSR build error");
		};

		checkSector("N-1",
		            sector.upWordsNm1,
		            sector.dnWordsNm1,
		            sector.dim1Nm1,
		            sector.csrNm1,
		            sector.varNm1);
		checkSector("N+1",
		            sector.upWordsNp1,
		            sector.dnWordsNp1,
		            sector.dim1Np1,
		            sector.csrNp1,
		            sector.varNp1);
	}

	// ========== Fock-space Hamiltonian matrix-vector product ==========

	// Apply the post-quench extended Hamiltonian H_ext(V^+) to state v → Hv.
	// H = U_f n↑n↓|imp + Σ_i potPost_[i]*(n↑+n↓)|_i
	//   + Σ_{α,σ} V_α (c†_0 c_{α+1} + h.c.)             [first bath, real]
	//   + Σ_{p,σ} V^+_p (c†_0 c_{ep} + h.c.)            [second bath empty sites]
	//   + Σ_{p,σ} V^+_p (c†_0 c_{op} + h.c.)            [second bath occupied sites]
	void applyHext(const VectorComplexType&        v,
	               VectorComplexType&              hv,
	               const std::vector<ComplexType>& gbekHop,
	               const std::vector<WordType>&    upWords,
	               const std::vector<WordType>&    dnWords,
	               SizeType                        dim1) const
	{
		const SizeType dim    = v.size();
		const SizeType spinUp = LanczosPlusPlus::LanczosGlobals::SPIN_UP;
		const SizeType dim2   = (dim1 > 0) ? (dim / dim1) : 0;

		std::fill(hv.begin(), hv.end(), ComplexType(0));

		for (SizeType m = 0; m < dim; ++m) {
			if (std::abs(v[m]) < RealType(1e-15))
				continue;

			const SizeType    y    = m / dim1;
			const SizeType    x    = m % dim1;
			const WordType    up_m = upWords[x];
			const WordType    dn_m = dnWords[y];
			const ComplexType vm   = v[m];

			// --- Diagonal terms ---
			RealType diag = 0;
			if ((up_m & WordType(1)) && (dn_m & WordType(1)))
				diag += params_.uFinal; // Hubbard U at impurity (site 0)
			for (SizeType i = 0; i < nsites_ext_; ++i) {
				if ((up_m >> i) & WordType(1))
					diag += potPost_[i];
				if ((dn_m >> i) & WordType(1))
					diag += potPost_[i];
			}
			hv[m] += ComplexType(diag) * vm;

			// --- Hopping terms for each spin ---
			for (SizeType spin = 0; spin < 2; ++spin) {
				const WordType w       = (spin == spinUp) ? up_m : dn_m;
				const WordType w_other = (spin == spinUp) ? dn_m : up_m;

				// Apply c†_to c_from (hop from → to) on spin word w
				auto applyHop = [&](SizeType from, SizeType to, ComplexType t)
				{
					if (!((w >> from) & WordType(1)))
						return; // source unoccupied
					if ((w >> to) & WordType(1))
						return; // target occupied

					const WordType new_w
					    = (w & ~(WordType(1) << from)) | (WordType(1) << to);
					const WordType new_up = (spin == spinUp) ? new_w : w_other;
					const WordType new_dn = (spin == spinUp) ? w_other : new_w;

					const SizeType m_new = lookupIndex(
					    new_up, new_dn, upWords, dnWords, dim1, dim2);
					if (m_new == std::numeric_limits<SizeType>::max())
						return;

					const int sign = computeJWSign(w, from, to);
					hv[m_new] += t * ComplexType(sign) * vm;
				};

				// First bath: site 0 ↔ site α+1 (real hoppings V_α)
				for (SizeType a = 0; a < nBath_; ++a) {
					const ComplexType Va(firstBathHop_[a]);
					applyHop(0, a + 1, Va);
					applyHop(a + 1, 0, Va); // real: conj = same
				}

				// Second bath: site 0 ↔ empty sites and occupied sites (complex
				// V^+_p)
				for (SizeType p = 0; p < bathRank_; ++p) {
					const SizeType    emptyS = 1 + nBath_ + p;
					const SizeType    occS   = 1 + nBath_ + bathRank_ + p;
					const ComplexType Vp     = gbekHop[p];
					applyHop(0, emptyS, Vp);
					applyHop(emptyS, 0, std::conj(Vp));
					applyHop(0, occS, Vp);
					applyHop(occS, 0, std::conj(Vp));
				}
			}
		}
	}

	// ========== Fock-space index lookup ==========

	// Find the basis index for (up_word, dn_word) using pre-sorted arrays.
	// Returns std::numeric_limits<SizeType>::max() if not found (shouldn't happen for valid
	// hops).
	static SizeType lookupIndex(WordType                     up,
	                            WordType                     dn,
	                            const std::vector<WordType>& upWords,
	                            const std::vector<WordType>& dnWords,
	                            SizeType                     dim1,
	                            SizeType /*dim2*/)
	{
		const auto it = std::lower_bound(upWords.begin(), upWords.end(), up);
		if (it == upWords.end() || *it != up)
			return std::numeric_limits<SizeType>::max();
		const SizeType x = static_cast<SizeType>(it - upWords.begin());

		const auto jt = std::lower_bound(dnWords.begin(), dnWords.end(), dn);
		if (jt == dnWords.end() || *jt != dn)
			return std::numeric_limits<SizeType>::max();
		const SizeType y = static_cast<SizeType>(jt - dnWords.begin());

		return x + y * dim1;
	}

	// Jordan-Wigner sign for c†_{to} c_{from} applied to spin word w.
	// Assumes bit 'from' is set and bit 'to' is unset in w.
	static int computeJWSign(WordType w, SizeType from, SizeType to)
	{
		// P_from: # bits at positions 0..from-1 in w
		int parity = 0;
		for (SizeType k = 0; k < from; ++k)
			parity += static_cast<int>((w >> k) & WordType(1));
		// After destroying at 'from':
		const WordType w1 = w ^ (WordType(1) << from);
		// P_to: # bits at positions 0..to-1 in w1
		for (SizeType k = 0; k < to; ++k)
			parity += static_cast<int>((w1 >> k) & WordType(1));
		return (parity & 1) ? -1 : 1;
	}

	// ========== Green's function formulas ==========

	// G^<(t_n, t_j) = +i <Ψ(t_j) | Ψ(t_n)>, one spin configuration.
	static ComplexType gLesserGBEKSector(const ExtendedSector& sector, int n, int j)
	{
		const SizeType dim = sector.PsiHist[0].size();
		ComplexType    s(0);
		for (SizeType k = 0; k < dim; ++k)
			s += std::conj(sector.PsiHist[static_cast<SizeType>(j)][k])
			    * sector.PsiHist[static_cast<SizeType>(n)][k];
		return ComplexType(0, 1) * s;
	}

	// G^>(t_n, t_j) = -i <Φ(t_n) | Φ(t_j)>, one spin configuration.
	static ComplexType gGreaterGBEKSector(const ExtendedSector& sector, int n, int j)
	{
		const SizeType dim = sector.PhiHist[0].size();
		ComplexType    s(0);
		for (SizeType k = 0; k < dim; ++k)
			s += std::conj(sector.PhiHist[static_cast<SizeType>(n)][k])
			    * sector.PhiHist[static_cast<SizeType>(j)][k];
		return ComplexType(0, -1) * s;
	}

	// Gα/Gβ average (GBEK Eq. 70): G_up = (1/2)(Gα_up + Gβ_up). When
	// nup_==ndown_, system beta is identical to system alpha and was never
	// built, so this reduces to just system alpha's value.
	ComplexType gLesserGBEK(int n, int j) const
	{
		if (sameConfig_)
			return gLesserGBEKSector(sectorAlpha_, n, j);
		return RealType(0.5)
		    * (gLesserGBEKSector(sectorAlpha_, n, j)
		       + gLesserGBEKSector(sectorBeta_, n, j));
	}

	ComplexType gGreaterGBEK(int n, int j) const
	{
		if (sameConfig_)
			return gGreaterGBEKSector(sectorAlpha_, n, j);
		return RealType(0.5)
		    * (gGreaterGBEKSector(sectorAlpha_, n, j)
		       + gGreaterGBEKSector(sectorBeta_, n, j));
	}

	// ========== Vector helpers ==========

	static void axpy(VectorComplexType& w, ComplexType c, const VectorComplexType& v)
	{
		for (SizeType i = 0; i < w.size(); ++i)
			w[i] += c * v[i];
	}

	static void scaleVec(VectorComplexType& v, ComplexType c)
	{
		for (SizeType i = 0; i < v.size(); ++i)
			v[i] *= c;
	}

	static ComplexType innerProduct(const VectorComplexType& a, const VectorComplexType& b)
	{
		ComplexType s(0);
		for (SizeType i = 0; i < a.size(); ++i)
			s += std::conj(a[i]) * b[i];
		return s;
	}

	// ========== LanczosPlusPlus infrastructure (mirrors ExactDiag) ==========

	// Minimal LanczosPlusPlus input string for a star-geometry Anderson model.
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
			if (i > 0)
				connStr += ",";
			connStr += ttos(hoppings[i]);
		}
		connStr += "]";

		std::string potStr = "[";
		for (SizeType i = 0; i < nsites; ++i) {
			if (i > 0)
				potStr += ",";
			potStr += ttos(potV[i]);
		}
		potStr += ",";
		for (SizeType i = 0; i < nsites; ++i) {
			if (i > 0)
				potStr += ",";
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
		s += "OutputFile=neqGbekDummy;\n";
		s += "InfiniteLoopKeptStates=1;\n";
		s += "FiniteLoops=0 0 0;\n";
		s += "TargetElectronsUp=" + ttos(nup) + ";\n";
		s += "TargetElectronsDown=" + ttos(ndown) + ";\n";
		s += "dir0:Connectors=" + connStr + ";\n";
		s += "potentialV=" + potStr + ";\n";
		return s;
	}

	// Full diagonalisation of a sector of the Hamiltonian.
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

	// Lanczos ground state for sectors too large for full diagonalisation (L>=2).
	// Returns a dim×1 matrix whose column 0 is the normalised GS eigenvector.
	static MatrixType lanczosGS(const ModelBaseType& model, const GeometryType& geom)
	{
		const BasisBaseType&      basis = model.basis();
		const SizeType            dim   = basis.size();
		DefaultSymmetryType       rs(basis, geom, "");
		InternalProductStoredType ham(model, basis, rs);

		PsimagLite::ParametersForSolver<RealType> lparams;
		lparams.steps      = static_cast<SizeType>(std::min(dim, SizeType(200)));
		lparams.tolerance  = 1e-12;
		lparams.lotaMemory = true;

		LanczosSolverForGSType lanczos(ham, lparams);

		typename LanczosSolverForGSType::VectorType initVec(
		    dim, ComplexOrRealType(1.0 / std::sqrt(RealType(dim))));
		RealType                                    gsEnergy = 0;
		typename LanczosSolverForGSType::VectorType gsVec(dim, ComplexOrRealType(0));

		lanczos.computeOneState(gsEnergy, gsVec, initVec, 0);

		MatrixType result(dim, 1);
		for (SizeType i = 0; i < dim; ++i)
			result(i, 0) = gsVec[i];
		return result;
	}

	// ========== Member variables ==========

	SizeType                    bathRank_;
	const ParamsNeqType&        params_;
	ExactDiagType               exactDiag_;
	std::unique_ptr<DecompType> decomp_;

	// Electron counts (read from io for L>0)
	SizeType nup_, ndown_;

	// Extended system geometry
	SizeType nBath_; // # first-bath sites
	SizeType nsites_ext_; // total sites = nBath+1+2L

	// First bath parameters (extracted from bathParams in solveLplus)
	VectorRealType firstBathHop_;
	VectorRealType firstBathEps_;

	// Post-quench on-site potentials for the extended system (shared by both
	// spin configurations -- U_final and first-bath epsilons don't depend on
	// which spin the impurity's extra electron carries).
	VectorRealType potPost_;

	// Gα (nup_ext=nup+L, ndown_ext=ndown+L) and Gβ (swapped) extended-Fock
	// state (GBEK Eq. 70). sameConfig_ is true when nup_==ndown_, in which
	// case sectorBeta_ is never built and is identical to sectorAlpha_.
	mutable ExtendedSector sectorAlpha_, sectorBeta_;
	bool                   sameConfig_ = false;

	friend struct GBEKTestAccessor;
};

} // namespace Dmft
#endif // IMPURITYSOLVER_NEQ_GBEK_H
