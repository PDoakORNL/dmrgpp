#ifndef NEQ_BATH_DECOMPOSITION_H
#define NEQ_BATH_DECOMPOSITION_H
#include "KadanoffBaym.h"
#include "Svd.h"
#include "Vector.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>

// Bath decomposition for neq-DMFT following Gramsch, Balzer, Eckstein, Kollar
// PRB 88, 235106 (2013).
//
// Splits the hybridization Δ = Δ⁻ + Δ⁺:
//   Δ⁻ : first bath — memory of the initial equilibrium state,
//         computed analytically from the free propagators of the equilibrium {V_α, ε_α}.
//   Δ⁺ : second bath — neq dynamics accumulated for t>0,
//         represented via a rank-L low-rank Cholesky of i·Δ⁺_<.
//
// At each time step n (in order 0,1,2,...) call update(n, delta) after
// delta has been filled up to row n.  The resulting second-bath hoppings
// are available via Vplus(n, p).
//
// bathParams layout: [V_0 .. V_{N-1},  ε_0 .. ε_{N-1}]
namespace Dmft {

template <typename ComplexOrRealType> class NeqBathDecomposition {

public:

	using RealType          = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using ComplexType       = std::complex<RealType>;
	using VectorRealType    = typename PsimagLite::Vector<RealType>::Type;
	using VectorComplexType = typename PsimagLite::Vector<ComplexType>::Type;
	using MatrixComplexType = PsimagLite::Matrix<ComplexType>;
	using KBType            = KadanoffBaym<ComplexOrRealType>;

	NeqBathDecomposition(SizeType              rank,
	                     RealType              beta,
	                     RealType              mu,
	                     const VectorRealType& bathParams,
	                     SizeType              nT,
	                     SizeType              nTau,
	                     RealType              dt,
	                     RealType              dtau)
	    : rank_(rank)
	    , beta_(beta)
	    , mu_(mu)
	    , nT_(nT)
	    , nTau_(nTau)
	    , dt_(dt)
	    , dtau_(dtau)
	    , V_(nT + 1, std::max(rank, SizeType(1)), ComplexType(0))
	{
		const SizeType nBath = bathParams.size() / 2;
		hoppings_.resize(nBath);
		bathEps_.resize(nBath);
		for (SizeType i = 0; i < nBath; ++i) {
			hoppings_[i] = bathParams[i];
			bathEps_[i]  = bathParams[nBath + i];
		}
	}

	// Incremental Cholesky update for time step n.
	// delta.lesser(n, j) for j <= n must be filled before calling.
	//
	// n=0 is skipped: at t=0 the system is in equilibrium so Δ⁺<(0,0)=0.
	// Skipping ensures V_(0,p)=0 exactly. The standard Cholesky phase then
	// seeds column p from pivot row p+1 (n=1..L), giving true rank-L behavior.
	// Without the skip col 0 of V stays zero (L_eff=1) because V_(0,0)≈0 from
	// the bath-fit residual causes V_(1,0)=Δ⁺<(1,0)/V_(0,0) to blow up.
	void update(int n, const KBType& delta)
	{
		if (rank_ == 0)
			return;
		if (n == 0)
			return;

		const int L = static_cast<int>(rank_);

		if (n <= L) {
			// Standard Cholesky: n=1 seeds col 0, n=2 seeds col 1, ..., n=L seeds col
			// L-1. Pivot for column p lives at row p+1 (row 0 is degenerate).
			const int col = n - 1;
			for (int p = 0; p < col; ++p)
				V_(n, p) = offDiagElement(n, p, delta);

			RealType d = -std::real(iDeltaPlusLesser(n, n, delta));
			for (int k = 0; k < col; ++k)
				d -= std::norm(V_(n, k));
			V_(n, col) = std::sqrt(std::max(d, RealType(0)));
		} else {
			// Optimal update: minimize ||Q q - a||² where Q_{kp} = conj(V_(k,p)),
			// k = 1..L (reference rows 1..L; row 0 is degenerate and excluded).
			choleskyOptimalUpdate(n, delta);
		}
	}

	// V^+(t_n, p): second-bath hopping for orbital p = 0..rank_-1
	ComplexType Vplus(int n, int p) const
	{
		assert(p >= 0 && static_cast<SizeType>(p) < rank_);
		return V_(n, p);
	}

	// V^-(t_n, α): first-bath hopping V_α exp(−i ε_α t_n)
	ComplexType Vminus(int n, int alpha) const
	{
		const ComplexType I(0, 1);
		return static_cast<ComplexType>(hoppings_[alpha])
		    * std::exp(-I * static_cast<ComplexType>(bathEps_[alpha] * n * dt_));
	}

	// Δ⁻_<(t_n, t_j): first-bath lesser hybridization
	ComplexType deltaMinusLesser(int n, int j) const
	{
		const ComplexType I(0, 1);
		ComplexType       result(0);
		for (SizeType a = 0; a < hoppings_.size(); ++a) {
			const RealType V2    = hoppings_[a] * hoppings_[a];
			const RealType eps   = bathEps_[a];
			const RealType fermi = fermiFunc(eps);
			result += V2 * I * fermi
			    * std::exp(-I * static_cast<ComplexType>(eps * (n - j) * dt_));
		}
		return result;
	}

	// Δ⁻_R(t_n, t_j): first-bath retarded hybridization  (n >= j)
	ComplexType deltaMinusRetarded(int n, int j) const
	{
		assert(n >= j);
		const ComplexType I(0, 1);
		ComplexType       result(0);
		for (SizeType a = 0; a < hoppings_.size(); ++a) {
			const RealType V2  = hoppings_[a] * hoppings_[a];
			const RealType eps = bathEps_[a];
			result += V2 * (-I)
			    * std::exp(-I * static_cast<ComplexType>(eps * (n - j) * dt_));
		}
		return result;
	}

	// Δ⁻_⌐(t_n, τ_j): first-bath left-mixing hybridization
	// Uses g_α^M(τ) = −f(ε_α) exp(−ε_α τ) for the free bath-site Matsubara GF
	ComplexType deltaMinusLeftMixing(int n, int j) const
	{
		const ComplexType I(0, 1);
		ComplexType       result(0);
		for (SizeType a = 0; a < hoppings_.size(); ++a) {
			const RealType V2    = hoppings_[a] * hoppings_[a];
			const RealType eps   = bathEps_[a];
			const RealType fermi = fermiFunc(eps);
			const RealType tau   = j * dtau_;
			result += V2 * std::exp(-I * static_cast<ComplexType>(eps * n * dt_))
			    * (-fermi) * static_cast<ComplexType>(std::exp(-eps * tau));
		}
		return result;
	}

	SizeType rank() const { return rank_; }
	SizeType nFirst() const { return hoppings_.size(); }

	// Reconstruct -iΔ⁺_<(t_n, t_j) = Σ_p V(n,p) conj(V(j,p)) and write to file.
	// File columns: n  j  Re(-iΔ⁺_<)  Im(-iΔ⁺_<)
	void dumpPlusBath(const std::string& filename) const
	{
		std::ofstream f(filename);
		for (SizeType n = 0; n <= nT_; ++n) {
			for (SizeType j = 0; j <= nT_; ++j) {
				ComplexType val(0);
				for (SizeType p = 0; p < rank_; ++p)
					val += V_(n, p) * std::conj(V_(j, p));
				f << n * dt_ << " " << j * dt_ << " " << val.real() << " "
				  << val.imag() << "\n";
			}
		}
	}

private:

	// Fermi function f(ε) = 1 / (1 + exp(β(ε − μ)))
	RealType fermiFunc(RealType eps) const
	{
		const RealType x = beta_ * (eps - mu_);
		if (x > 500)
			return 0;
		if (x < -500)
			return 1;
		return 1 / (1 + std::exp(x));
	}

	// i Δ⁺_<(t_n, t_j): helper with antisymmetry for j > n
	ComplexType iDeltaPlusLesser(int n, int j, const KBType& delta) const
	{
		const ComplexType I(0, 1);
		ComplexType       dless;
		if (j <= n)
			dless = static_cast<ComplexType>(delta.lesser(n, j));
		else
			dless = -std::conj(static_cast<ComplexType>(delta.lesser(j, n)));
		return I * (dless - deltaMinusLesser(n, j));
	}

	// Standard Cholesky off-diagonal element V_(n, p) in the standard phase.
	// Pivot for column p is at row p+1 (row 0 is degenerate and skipped).
	ComplexType offDiagElement(int n, int p, const KBType& delta) const
	{
		ComplexType val = -iDeltaPlusLesser(n, p + 1, delta);
		for (int k = 0; k < p; ++k)
			val -= V_(n, k) * std::conj(V_(p + 1, k));
		const ComplexType pivot = V_(p + 1, p);
		if (std::abs(pivot) < 1e-14)
			return ComplexType(0);
		return val / std::conj(pivot);
	}

	// Optimal update for n > rank_: solve the JOINT minimization from GBEK
	// Eq. 63,
	//   F(q) = 2||Q_s q - a||^2 + |q^H q - d|^2,
	// where Q_s is the s x L matrix of ALL previously-determined rows 1..s
	// (s = n-1), growing by one row every step, and d = -iΔ⁺_<(t_n,t_n) is
	// the target diagonal value.
	//
	// Two bugs, fixed in sequence (see
	// cincuenta/TestSuite/gbek_reference/gbek_cholesky.py for the
	// independent Python implementation both were caught against):
	//
	//  1. Fixed reference window: an earlier version of this function
	//     hardcoded the reference set to the first L seeding rows forever,
	//     instead of growing it to all s = n-1 rows. This discards almost
	//     all previously-built structure once n grows large relative to L,
	//     causing far-too-fast (near-total) collapse of the reconstructed
	//     hybridization.
	//  2. Missing diagonal constraint: after fixing (1), this function
	//     still solved only the linear part of Eq. 63 (Q^H Q q = Q^H a),
	//     ignoring the q^H q ~ d term -- NOT separable from the linear fit,
	//     since F(q) couples them. This systematically undershoots the
	//     diagonal by a smaller but still real amount (verified: ~5-15% on
	//     realistic full-rank targets, growing worse with n).
	//
	// a[k] = -iΔ⁺_<(t_n, t_{k+1}) for k = 0..s-1 (row k+1), one entry per
	// row of Q_s.
	void choleskyOptimalUpdate(int n, const KBType& delta)
	{
		const int L = static_cast<int>(rank_);
		const int s = n - 1; // number of previously-determined rows (1..n-1)

		VectorComplexType a(s);
		for (int k = 0; k < s; ++k)
			a[k] = -iDeltaPlusLesser(n, k + 1, delta);
		// Gram matrix QtQ[p][p'] = Σ_{k=1}^{s} V_(k,p) conj(V_(k,p'))
		MatrixComplexType QtQ(L, L, ComplexType(0));
		for (int p = 0; p < L; ++p) {
			for (int pp = 0; pp < L; ++pp) {
				ComplexType sum(0);
				for (int k = 1; k <= s; ++k)
					sum += V_(k, p) * std::conj(V_(k, pp));
				QtQ(p, pp) = sum;
			}
		}

		// Rhs Qta[p] = Σ_{k=1}^{s} V_(k,p) a[k-1]
		VectorComplexType Qta(L, ComplexType(0));
		for (int p = 0; p < L; ++p)
			for (int k = 1; k <= s; ++k)
				Qta[p] += V_(k, p) * a[k - 1];

		const RealType    d = -std::real(iDeltaPlusLesser(n, n, delta));
		VectorComplexType q = solveOptimalUpdateJoint(QtQ, Qta, d, L);
		for (int p = 0; p < L; ++p)
			V_(n, p) = q[p];
	}

	// Solve GBEK Eq. 63's joint objective F(q) = 2||Qq-a||^2 + |q^H q-d|^2
	// exactly. Setting the Wirtinger gradient dF/dq* to zero gives (the
	// factors of 2 cancel completely -- easy to get wrong; a draft of this
	// fix kept a stray factor of 2 and got answers close to, but measurably
	// different from, the true minimum found by numerically minimizing F(q)
	// directly with BFGS in Python):
	//
	//   [QtQ + (mu - d) I] q(mu) = Qta,   mu = q(mu)^H q(mu)
	//
	// This is ONE nonlinear SCALAR equation in mu (not a joint multivariate
	// optimization): for fixed mu, q(mu) solves a small L x L linear
	// system; g(mu) = ||q(mu)||^2 - mu = 0 is then a 1D root-find, solved
	// here by bisection after bracketing -- no nonlinear-optimization
	// library needed.
	//
	// At mu = d, this reduces exactly to the old (bug 2) linear-only
	// solution: [QtQ] q = Qta. Verified ||q(d)||^2 <= d in every case
	// checked; g is decreasing as mu decreases below d, so the root lies
	// in [0, d] and bisection there is robust.
	static VectorComplexType solveOptimalUpdateJoint(const MatrixComplexType& QtQ,
	                                                 const VectorComplexType& Qta,
	                                                 RealType                 d,
	                                                 int                      L)
	{
		auto qOf = [&](RealType mu)
		{
			MatrixComplexType A = QtQ;
			for (int p = 0; p < L; ++p)
				A(p, p) += ComplexType(mu - d, 0);
			return solveLinear(A, Qta, L);
		};
		auto gOf = [&](RealType mu)
		{
			VectorComplexType q     = qOf(mu);
			RealType          norm2 = 0;
			for (int p = 0; p < L; ++p)
				norm2 += std::norm(q[p]);
			return norm2 - mu;
		};

		RealType muHi = d;
		if (gOf(muHi) >= 0)
			return qOf(muHi);

		RealType muLo = muHi;
		RealType step = std::max(d, RealType(1.0)) * 0.5;
		for (int iter = 0; iter < 60; ++iter) {
			muLo = std::max(muLo - step, RealType(0.0));
			if (gOf(muLo) >= 0 || muLo <= 0.0)
				break;
			step *= 1.5;
		}

		for (int iter = 0; iter < 100; ++iter) {
			const RealType muMid = 0.5 * (muLo + muHi);
			if (gOf(muMid) >= 0)
				muLo = muMid;
			else
				muHi = muMid;
			if (muHi - muLo < 1e-14 * std::max(d, RealType(1.0)))
				break;
		}
		return qOf(muLo);
	}

	// Least-squares solve of A q = b (A is L×L, Hermitian in practice since
	// it is QtQ + scalar*I, but solved generally) via SVD-based pseudo-
	// inverse, truncating singular values below rcond*max(s) -- matching
	// numpy.linalg.lstsq's default behavior.
	//
	// This replaced a plain Gaussian-elimination-with-partial-pivoting
	// solve (skip-if-|pivot|<1e-14) that was found, via independent
	// cross-check against cincuenta/TestSuite/gbek_reference/gbek_cholesky.py
	// on the true atomic-limit target (Delta^- === 0 exactly, an extreme,
	// genuinely rank-deficient seed), to silently return large-but-finite
	// garbage on NEAR-singular (not exactly singular) systems -- exactly
	// the "near-singular Gram matrix blowup" failure mode documented in
	// gbek_cholesky.py's module docstring, which was fixed on the Python
	// side (switched to np.linalg.lstsq) but never ported here. The two
	// implementations agreed on the near-atomic (W_i=0.1) target because
	// that target's early rows are not exactly singular, only degenerate;
	// the true atomic limit's Delta_<(0,0)=0 exactly is a harsher seed and
	// exposed the gap.
	static VectorComplexType
	solveLinear(const MatrixComplexType& A, const VectorComplexType& b, int L)
	{
		MatrixComplexType            u(A);
		VectorRealType               s;
		MatrixComplexType            vt;
		PsimagLite::Svd<ComplexType> svd;
		svd('S', u, s, vt);

		// rcond=1e-10 matches gbek_cholesky.py::_solve_optimal_update's
		// np.linalg.lstsq(A, Qta, rcond=1e-10) exactly -- deliberately far
		// looser than machine-epsilon truncation, since the near-singular
		// directions here are not noise but genuinely poorly-determined
		// combinations that must be dropped, not merely rounding error.
		const RealType sMax
		    = (s.empty()) ? RealType(0) : *std::max_element(s.begin(), s.end());
		const RealType rcond  = RealType(1e-10);
		const RealType thresh = rcond * sMax;

		// x = V * diag(1/s_i, truncated) * U^H * b
		VectorComplexType uhB(L, ComplexType(0));
		for (int i = 0; i < L; ++i)
			for (int k = 0; k < L; ++k)
				uhB[i] += std::conj(u(k, i)) * b[k];

		VectorComplexType x(L, ComplexType(0));
		for (int i = 0; i < L; ++i) {
			if (s[i] <= thresh)
				continue;
			const ComplexType coeff = uhB[i] / s[i];
			for (int c = 0; c < L; ++c)
				x[c] += std::conj(vt(i, c)) * coeff;
		}
		return x;
	}

	SizeType          rank_;
	RealType          beta_, mu_, dt_, dtau_;
	SizeType          nT_, nTau_;
	VectorRealType    hoppings_; // V_α from equilibrium bath
	VectorRealType    bathEps_; // ε_α from equilibrium bath
	MatrixComplexType V_; // second-bath hoppings V_[n][p], (nT+1) × rank_
};

} // namespace Dmft
#endif // NEQ_BATH_DECOMPOSITION_H
