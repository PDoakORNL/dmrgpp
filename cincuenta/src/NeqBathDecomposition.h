#ifndef NEQ_BATH_DECOMPOSITION_H
#define NEQ_BATH_DECOMPOSITION_H
#include "KadanoffBaym.h"
#include "Vector.h"
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

	// Optimal update for n >= rank_: solve normal equations Q†Q q = Q†a
	// Q_{kp} = conj(V_{kp}) for k=0..L-1, a_k = iΔ⁺_<(t_n, t_k)
	void choleskyOptimalUpdate(int n, const KBType& delta)
	{
		const int L = static_cast<int>(rank_);

		// Reference rows 1..L (row 0 is degenerate and excluded from both Gram and rhs).
		// Build rhs vector a[k] = -iΔ⁺_<(t_n, t_{k+1}) for k = 0..L-1 (row k+1)
		VectorComplexType a(L);
		for (int k = 0; k < L; ++k)
			a[k] = -iDeltaPlusLesser(n, k + 1, delta);
		// Gram matrix G[p][p'] = Σ_{k=1}^{L} V_(k,p) conj(V_(k,p'))
		MatrixComplexType G(L, L, ComplexType(0));
		for (int p = 0; p < L; ++p) {
			for (int pp = 0; pp < L; ++pp) {
				ComplexType sum(0);
				for (int k = 1; k <= L; ++k)
					sum += V_(k, p) * std::conj(V_(k, pp));
				G(p, pp) = sum;
			}
		}

		// Rhs b[p] = Σ_{k=1}^{L} V_(k,p) a[k-1]
		VectorComplexType b(L, ComplexType(0));
		for (int p = 0; p < L; ++p)
			for (int k = 1; k <= L; ++k)
				b[p] += V_(k, p) * a[k - 1];

		VectorComplexType q = solveLinear(G, b, L);
		for (int p = 0; p < L; ++p)
			V_(n, p) = q[p];
	}

	// Gaussian elimination with partial pivoting for an L×L system Ax = b.
	static VectorComplexType solveLinear(MatrixComplexType A, VectorComplexType b, int L)
	{
		// Forward elimination
		for (int col = 0; col < L; ++col) {
			// Find pivot
			int      pivot  = col;
			RealType maxVal = std::abs(A(col, col));
			for (int row = col + 1; row < L; ++row) {
				RealType v = std::abs(A(row, col));
				if (v > maxVal) {
					maxVal = v;
					pivot  = row;
				}
			}
			// Swap rows
			if (pivot != col) {
				for (int c = 0; c < L; ++c)
					std::swap(A(col, c), A(pivot, c));
				std::swap(b[col], b[pivot]);
			}
			if (std::abs(A(col, col)) < 1e-14)
				continue; // singular, skip
			// Eliminate
			for (int row = col + 1; row < L; ++row) {
				const ComplexType factor = A(row, col) / A(col, col);
				for (int c = col; c < L; ++c)
					A(row, c) -= factor * A(col, c);
				b[row] -= factor * b[col];
			}
		}
		// Back substitution
		VectorComplexType x(L, ComplexType(0));
		for (int row = L - 1; row >= 0; --row) {
			if (std::abs(A(row, row)) < 1e-14)
				continue;
			ComplexType sum = b[row];
			for (int c = row + 1; c < L; ++c)
				sum -= A(row, c) * x[c];
			x[row] = sum / A(row, row);
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
