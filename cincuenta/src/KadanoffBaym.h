#ifndef KADANOFFBAYM_H
#define KADANOFFBAYM_H
#include "Matrix.h"
#include "Vector.h"
#include <algorithm>
#include <cassert>
#include <fstream>
#include <string>

// Translated and adapted from sample code by N. Tsuji accompanying:
//   H. Aoki, N. Tsuji, M. Eckstein, M. Kollar, T. Oka, P. Werner,
//   "Nonequilibrium dynamical mean-field theory and its applications",
//   Rev. Mod. Phys. 86, 779 (2014).  https://doi.org/10.1103/RevModPhys.86.779
//
// Per the authors: the sample programs may be used and modified for
// non-commercial purposes, but their use must be acknowledged in
// publications with a citation to the review article above.
//
// Changes from the original noneq-dmft/cxx (2013): merged green.h/green.cxx
// into a single header, templated on ComplexOrRealType, and replaced bare
// STL containers with PsimagLite Matrix/Vector types.

namespace Dmft {

template <typename ComplexOrRealType> class KBDerivative;

/*!
 * \brief KadanoffBaym
 * Two-time Kadanoff-Baym Green's function G(t,t') on a uniform real-time
 * grid 0..N_t and imaginary-time grid 0..N_tau.
 *
 * Convention (following Aoki et al. 2014):
 *   retarded(n,j)    = G^R(t_n, t_j),       j <= n
 *   left_mixing(n,j) = G^{Left}(t_n, tau_j)
 *   lesser(n,j)      = G^{<}(t_n, t_j)
 *   matsubara_w[k]   = G^M(i*omega_k),  k = 0..N_tau-1
 *   matsubara_t[j]   = G^M(tau_j),      j = 0..N_tau
 *     (matsubara_t has N_tau+1 entries so that the periodicity access
 *      matsubara_t[N_tau + k - j] = G^M(beta - (tau_j - tau_k)) is valid)
 */
template <typename ComplexOrRealType> class KadanoffBaym {

public:

	using RealType          = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using ComplexType       = std::complex<RealType>;
	using VectorComplexType = typename PsimagLite::Vector<ComplexType>::Type;
	using MatrixComplexType = PsimagLite::Matrix<ComplexType>;

	KadanoffBaym() = default;

	KadanoffBaym(SizeType nT, SizeType nTau, RealType dt, RealType dtau)
	    : matsubara_w(nTau)
	    , matsubara_t(nTau + 1)
	    , retarded(nT + 1, nT + 1)
	    , left_mixing(nT + 1, nTau + 1)
	    , lesser(nT + 1, nT + 1)
	    , nT_(nT)
	    , nTau_(nTau)
	    , dt_(dt)
	    , dtau_(dtau)
	{ }

	SizeType nT() const { return nT_; }
	SizeType nTau() const { return nTau_; }
	RealType dt() const { return dt_; }
	RealType dtau() const { return dtau_; }

	// Matsubara parts — populated during equilibrium initialisation
	VectorComplexType matsubara_w; ///< G^M(i*omega_k)
	VectorComplexType matsubara_t; ///< G^M(tau_j)

	// Real-time Kadanoff-Baym components
	MatrixComplexType retarded; ///< G^R(t_n, t_j)
	MatrixComplexType left_mixing; ///< G^{Left}(t_n, tau_j)
	MatrixComplexType lesser; ///< G^{<}(t_n, t_j)

	/*!
	 * \brief convolute
	 * C(t,t') = (A*B)(t,t') for t=n*dt or t'=n*dt.
	 * Integrals via trapezoidal rule.
	 */
	void convolute(int n, const KadanoffBaym& A, const KadanoffBaym& B)
	{
		const int         Ntau = static_cast<int>(nTau_);
		VectorComplexType AxB(std::max(Ntau + 1, n + 1));

		// C^R(t_n, t_j) = int_{t_j}^{t_n} ds A^R(t_n,s) B^R(s,t_j)
		for (int j = 0; j <= n; ++j) {
			retarded(n, j) = 0;
			for (int k = j; k <= n; ++k)
				AxB[k] = A.retarded(n, k) * B.retarded(k, j);
			retarded(n, j) += dt_ * trapezoid(AxB, j, n);
		}

		// C^{Left}(t_n, tau_j)
		//   = int_0^beta dtau A^{Left}(t_n,tau) B^M(tau,tau_j)
		//   + int_0^{t_n} ds  A^R(t_n,s) B^{Left}(s,tau_j)
		for (int j = 0; j <= Ntau; ++j) {
			left_mixing(n, j) = 0;
			for (int k = 0; k <= j; ++k)
				AxB[k] = A.left_mixing(n, k) * B.matsubara_t[Ntau + k - j];
			left_mixing(n, j) += -dtau_ * trapezoid(AxB, 0, j);
			for (int k = j; k <= Ntau; ++k)
				AxB[k] = A.left_mixing(n, k) * B.matsubara_t[k - j];
			left_mixing(n, j) += dtau_ * trapezoid(AxB, j, Ntau);
			for (int k = 0; k <= n; ++k)
				AxB[k] = A.retarded(n, k) * B.left_mixing(k, j);
			left_mixing(n, j) += dt_ * trapezoid(AxB, 0, n);
		}

		// C^{<}(t_n, t_j), j < n:
		//   = -i int_0^beta dtau A^{Left}(t_n,tau) [B^{Left}(t_j,beta-tau)]^*
		//   + int_0^{t_j} ds A^{<}(t_n,s) [B^R(t_j,s)]^*
		//   + int_0^{t_n} ds A^R(t_n,s)  B^{<}(s,t_j)
		const ComplexType I(0, 1);
		for (int j = 0; j <= n - 1; ++j) {
			lesser(n, j) = 0;
			for (int k = 0; k <= Ntau; ++k)
				AxB[k] = A.left_mixing(n, k)
				    * PsimagLite::conj(B.left_mixing(j, Ntau - k));
			lesser(n, j) += -I * dtau_ * trapezoid(AxB, 0, Ntau);
			for (int k = 0; k <= j; ++k)
				AxB[k] = A.lesser(n, k) * PsimagLite::conj(B.retarded(j, k));
			lesser(n, j) += dt_ * trapezoid(AxB, 0, j);
			for (int k = 0; k <= n; ++k)
				AxB[k] = A.retarded(n, k) * B.lesser(k, j);
			lesser(n, j) += dt_ * trapezoid(AxB, 0, n);
		}

		// C^{<}(t_i, t_n), i <= n  (fill the n-th column)
		for (int i = 0; i <= n; ++i) {
			lesser(i, n) = 0;
			for (int k = 0; k <= Ntau; ++k)
				AxB[k] = A.left_mixing(i, k)
				    * PsimagLite::conj(B.left_mixing(n, Ntau - k));
			lesser(i, n) += -I * dtau_ * trapezoid(AxB, 0, Ntau);
			for (int k = 0; k <= n; ++k)
				AxB[k] = A.lesser(i, k) * PsimagLite::conj(B.retarded(n, k));
			lesser(i, n) += dt_ * trapezoid(AxB, 0, n);
			for (int k = 0; k <= i; ++k)
				AxB[k] = A.retarded(i, k) * B.lesser(k, n);
			lesser(i, n) += dt_ * trapezoid(AxB, 0, i);
		}
	}

	/*!
	 * \brief volterra_int
	 * Solves G = G0 + K*G (Volterra integral, 2nd-order) at index n.
	 * G^R diagonal G(n,n) = G0(n,n) is the free propagator.
	 */
	void volterra_int(int n, const KadanoffBaym& G0, const KadanoffBaym& K)
	{
		const int         Ntau = static_cast<int>(nTau_);
		VectorComplexType KxG(std::max(Ntau + 1, n + 1));
		const ComplexType I(0, 1);

		// Retarded
		retarded(n, n) = G0.retarded(n, n);
		for (int j = 0; j <= n - 1; ++j) {
			retarded(n, j) = G0.retarded(n, j);
			for (int l = j; l <= n - 1; ++l)
				KxG[l] = K.retarded(n, l) * retarded(l, j);
			retarded(n, j) += dt_ * trapezoid_half_edge(KxG, j, n - 1);
			retarded(n, j) /= ComplexType(1) - RealType(0.5) * dt_ * K.retarded(n, n);
		}

		// Left-mixing
		for (int j = 0; j <= Ntau; ++j) {
			left_mixing(n, j) = G0.left_mixing(n, j);
			for (int l = 0; l <= j; ++l)
				KxG[l] = K.left_mixing(n, l) * matsubara_t[Ntau + l - j];
			left_mixing(n, j) += -dtau_ * trapezoid(KxG, 0, j);
			for (int l = j; l <= Ntau; ++l)
				KxG[l] = K.left_mixing(n, l) * matsubara_t[l - j];
			left_mixing(n, j) += dtau_ * trapezoid(KxG, j, Ntau);
			for (int l = 0; l <= n - 1; ++l)
				KxG[l] = K.retarded(n, l) * left_mixing(l, j);
			left_mixing(n, j) += dt_ * trapezoid_half_edge(KxG, 0, n - 1);
			left_mixing(n, j)
			    /= ComplexType(1) - RealType(0.5) * dt_ * K.retarded(n, n);
		}

		// Lesser — G^<(t_n, t_j) for j < n
		for (int j = 0; j <= n - 1; ++j) {
			lesser(n, j) = G0.lesser(n, j);
			for (int l = 0; l <= Ntau; ++l)
				KxG[l] = K.left_mixing(n, l)
				    * PsimagLite::conj(left_mixing(j, Ntau - l));
			lesser(n, j) += -I * dtau_ * trapezoid(KxG, 0, Ntau);
			for (int l = 0; l <= j; ++l)
				KxG[l] = K.lesser(n, l) * PsimagLite::conj(retarded(j, l));
			lesser(n, j) += dt_ * trapezoid(KxG, 0, j);
			for (int l = 0; l <= n - 1; ++l)
				KxG[l] = K.retarded(n, l) * lesser(l, j);
			lesser(n, j) += dt_ * trapezoid_half_edge(KxG, 0, n - 1);
			lesser(n, j) /= ComplexType(1) - RealType(0.5) * dt_ * K.retarded(n, n);
		}
		// G^{<}(t_i, t_n) = -[G^{<}(t_n, t_i)]^*  (Hermite conjugate)
		for (int i = 0; i <= n - 1; ++i)
			lesser(i, n) = -PsimagLite::conj(lesser(n, i));
		// Diagonal G^<(t_n, t_n)
		lesser(n, n) = G0.lesser(n, n);
		for (int l = 0; l <= Ntau; ++l)
			KxG[l] = K.left_mixing(n, l) * PsimagLite::conj(left_mixing(n, Ntau - l));
		lesser(n, n) += -I * dtau_ * trapezoid(KxG, 0, Ntau);
		for (int l = 0; l <= n; ++l)
			KxG[l] = K.lesser(n, l) * PsimagLite::conj(retarded(n, l));
		lesser(n, n) += dt_ * trapezoid(KxG, 0, n);
		for (int l = 0; l <= n - 1; ++l)
			KxG[l] = K.retarded(n, l) * lesser(l, n);
		lesser(n, n) += dt_ * trapezoid_half_edge(KxG, 0, n - 1);
		lesser(n, n) /= ComplexType(1) - RealType(0.5) * dt_ * K.retarded(n, n);
	}

	/*!
	 * \brief volterra_intdiff
	 * Solves [i d/dt - h(t)] G(t,t') = delta(t,t') + (K*G)(t,t')
	 * at index n using 2nd-order implicit Runge-Kutta (n >= 1).
	 * Updates both *this (G values) and G_der_new (derivative at t_n).
	 */
	void volterra_intdiff(int                              n,
	                      const VectorComplexType&         h,
	                      const KadanoffBaym&              K,
	                      KBDerivative<ComplexOrRealType>& G_der,
	                      KBDerivative<ComplexOrRealType>& G_der_new)
	{
		assert(n >= 1);
		const int         Ntau = static_cast<int>(nTau_);
		VectorComplexType KxG(std::max(Ntau + 1, n + 1));
		const ComplexType I(0, 1);

		// Retarded
		retarded(n, n)        = -I;
		G_der_new.retarded[n] = -I * h[n] * retarded(n, n);
		for (int j = 0; j <= n - 1; ++j) {
			retarded(n, j)
			    = retarded(n - 1, j) + RealType(0.5) * dt_ * G_der.retarded[j];
			for (int l = j; l <= n - 1; ++l)
				KxG[l] = K.retarded(n, l) * retarded(l, j);
			G_der_new.retarded[j] = -I * dt_ * trapezoid_half_edge(KxG, j, n - 1);
			retarded(n, j) += RealType(0.5) * dt_ * G_der_new.retarded[j];
			retarded(n, j) /= ComplexType(1) + RealType(0.5) * I * dt_ * h[n]
			    + RealType(0.25) * I * dt_ * dt_ * K.retarded(n, n);
			G_der_new.retarded[j] += -I * h[n] * retarded(n, j)
			    - RealType(0.5) * I * dt_ * K.retarded(n, n) * retarded(n, j);
		}

		// Left-mixing
		for (int j = 0; j <= Ntau; ++j) {
			left_mixing(n, j)
			    = left_mixing(n - 1, j) + RealType(0.5) * dt_ * G_der.left_mixing[j];
			for (int l = 0; l <= j; ++l)
				KxG[l] = K.left_mixing(n, l) * matsubara_t[Ntau + l - j];
			G_der_new.left_mixing[j] = I * dtau_ * trapezoid(KxG, 0, j);
			for (int l = j; l <= Ntau; ++l)
				KxG[l] = K.left_mixing(n, l) * matsubara_t[l - j];
			G_der_new.left_mixing[j] += -I * dtau_ * trapezoid(KxG, j, Ntau);
			for (int l = 0; l <= n - 1; ++l)
				KxG[l] = K.retarded(n, l) * left_mixing(l, j);
			G_der_new.left_mixing[j] += -I * dt_ * trapezoid_half_edge(KxG, 0, n - 1);
			left_mixing(n, j) += RealType(0.5) * dt_ * G_der_new.left_mixing[j];
			left_mixing(n, j) /= ComplexType(1) + RealType(0.5) * I * dt_ * h[n]
			    + RealType(0.25) * I * dt_ * dt_ * K.retarded(n, n);
			G_der_new.left_mixing[j] += -I * h[n] * left_mixing(n, j)
			    - RealType(0.5) * I * dt_ * K.retarded(n, n) * left_mixing(n, j);
		}

		// Lesser — G^<(t_n, t_j) for j < n
		for (int j = 0; j <= n - 1; ++j) {
			lesser(n, j) = lesser(n - 1, j) + RealType(0.5) * dt_ * G_der.lesser[j];
			for (int l = 0; l <= Ntau; ++l)
				KxG[l] = K.left_mixing(n, l)
				    * PsimagLite::conj(left_mixing(j, Ntau - l));
			G_der_new.lesser[j] = -I * (-I) * dtau_ * trapezoid(KxG, 0, Ntau);
			for (int l = 0; l <= j; ++l)
				KxG[l] = K.lesser(n, l) * PsimagLite::conj(retarded(j, l));
			G_der_new.lesser[j] += -I * dt_ * trapezoid(KxG, 0, j);
			for (int l = 0; l <= n - 1; ++l)
				KxG[l] = K.retarded(n, l) * lesser(l, j);
			G_der_new.lesser[j] += -I * dt_ * trapezoid_half_edge(KxG, 0, n - 1);
			lesser(n, j) += RealType(0.5) * dt_ * G_der_new.lesser[j];
			lesser(n, j) /= ComplexType(1) + RealType(0.5) * I * dt_ * h[n]
			    + RealType(0.25) * I * dt_ * dt_ * K.retarded(n, n);
			G_der_new.lesser[j] += -I * h[n] * lesser(n, j)
			    - RealType(0.5) * I * dt_ * K.retarded(n, n) * lesser(n, j);
		}
		// Hermite conjugate: G^<(t_i, t_n) = -[G^<(t_n, t_i)]^*
		for (int i = 0; i <= n - 1; ++i)
			lesser(i, n) = -PsimagLite::conj(lesser(n, i));

		// d/dt G^<(t_{n-1}, t_n) — needed for the diagonal step below
		ComplexType G_der_lesser = -I * h[n - 1] * lesser(n - 1, n);
		for (int l = 0; l <= Ntau; ++l)
			KxG[l]
			    = K.left_mixing(n - 1, l) * PsimagLite::conj(left_mixing(n, Ntau - l));
		G_der_lesser += -I * (-I) * dtau_ * trapezoid(KxG, 0, Ntau);
		for (int l = 0; l <= n; ++l)
			KxG[l] = K.lesser(n - 1, l) * PsimagLite::conj(retarded(n, l));
		G_der_lesser += -I * dt_ * trapezoid(KxG, 0, n);
		for (int l = 0; l <= n - 1; ++l)
			KxG[l] = K.retarded(n - 1, l) * lesser(l, n);
		G_der_lesser += -I * dt_ * trapezoid(KxG, 0, n - 1);

		// Diagonal G^<(t_n, t_n)
		lesser(n, n) = lesser(n - 1, n) + RealType(0.5) * dt_ * G_der_lesser;
		for (int l = 0; l <= Ntau; ++l)
			KxG[l] = K.left_mixing(n, l) * PsimagLite::conj(left_mixing(n, Ntau - l));
		G_der_new.lesser[n] = -I * (-I) * dtau_ * trapezoid(KxG, 0, Ntau);
		for (int l = 0; l <= n; ++l)
			KxG[l] = K.lesser(n, l) * PsimagLite::conj(retarded(n, l));
		G_der_new.lesser[n] += -I * dt_ * trapezoid(KxG, 0, n);
		for (int l = 0; l <= n - 1; ++l)
			KxG[l] = K.retarded(n, l) * lesser(l, n);
		G_der_new.lesser[n] += -I * dt_ * trapezoid_half_edge(KxG, 0, n - 1);
		lesser(n, n) += RealType(0.5) * dt_ * G_der_new.lesser[n];
		lesser(n, n) /= ComplexType(1) + RealType(0.5) * I * dt_ * h[n]
		    + RealType(0.25) * I * dt_ * dt_ * K.retarded(n, n);
		G_der_new.lesser[n] += -I * h[n] * lesser(n, n)
		    - RealType(0.5) * I * dt_ * K.retarded(n, n) * lesser(n, n);
	}

	/*!
	 * \brief dump
	 * Write all four KB components to plain-text files matching the noneq-dmft
	 * reference format (Aoki et al. 2014).  Files are named <prefix>-retarded, etc.
	 * matsubara_t writes "tau Re Im" (reference writes only Re; imaginary part
	 * is near zero for single-band systems with time-reversal symmetry).
	 */
	void dump(const std::string& prefix) const
	{
		const int Nt  = static_cast<int>(nT_);
		const int Nta = static_cast<int>(nTau_);

		// retarded: causal upper triangle, t >= t'
		{
			std::ofstream ofs(prefix + "-retarded");
			ofs.precision(10);
			for (int n = 0; n <= Nt; ++n) {
				const RealType t = n * dt_;
				for (int j = 0; j <= n; ++j) {
					const RealType tp = j * dt_;
					ofs << std::fixed << t << " " << tp << " "
					    << retarded(n, j).real() << " " << retarded(n, j).imag()
					    << "\n";
				}
			}
		}

		// lesser: full (nT+1) x (nT+1) matrix
		{
			std::ofstream ofs(prefix + "-lesser");
			ofs.precision(10);
			for (int n = 0; n <= Nt; ++n) {
				const RealType t = n * dt_;
				for (int j = 0; j <= Nt; ++j) {
					const RealType tp = j * dt_;
					ofs << std::fixed << t << " " << tp << " "
					    << lesser(n, j).real() << " " << lesser(n, j).imag()
					    << "\n";
				}
			}
		}

		// left_mixing: (nT+1) x (nTau+1)
		{
			std::ofstream ofs(prefix + "-left-mixing");
			ofs.precision(10);
			for (int n = 0; n <= Nt; ++n) {
				const RealType t = n * dt_;
				for (int j = 0; j <= Nta; ++j) {
					const RealType tau = j * dtau_;
					ofs << std::fixed << t << " " << tau << " "
					    << left_mixing(n, j).real() << " "
					    << left_mixing(n, j).imag() << "\n";
				}
			}
		}

		// matsubara_t: equilibrium imaginary-time GF, nTau+1 points
		{
			std::ofstream ofs(prefix + "-matsubara-t");
			ofs.precision(10);
			for (int j = 0; j <= Nta; ++j) {
				const RealType tau = j * dtau_;
				ofs << std::fixed << tau << " " << matsubara_t[j].real() << " "
				    << matsubara_t[j].imag() << "\n";
			}
		}
	}

private:

	SizeType nT_   = 0;
	SizeType nTau_ = 0;
	RealType dt_   = 0;
	RealType dtau_ = 0;

	/*!
	 * \brief trapezoid
	 * Trapezoidal weights: w_i = w_j = 1/2, w_k = 1 for i < k < j.
	 * Returns zero when i == j (empty interval).
	 */
	static ComplexType trapezoid(const VectorComplexType& f, int i, int j)
	{
		if (j == i)
			return ComplexType(0);
		ComplexType s = RealType(0.5) * (f[i] + f[j]);
		for (int k = i + 1; k < j; ++k)
			s += f[k];
		return s;
	}

	/*!
	 * \brief trapezoid_half_edge
	 * Half-edge trapezoid: w_i = 1/2, w_k = 1 for k > i.
	 */
	static ComplexType trapezoid_half_edge(const VectorComplexType& f, int i, int j)
	{
		ComplexType s = RealType(0.5) * f[i];
		for (int k = i + 1; k <= j; ++k)
			s += f[k];
		return s;
	}
};

/*!
 * \brief KBDerivative
 * Stores d/dt G(t_n, t_j) for fixed n at the current time step.
 * Required by volterra_intdiff for the Runge-Kutta predictor-corrector.
 */
template <typename ComplexOrRealType> class KBDerivative {

public:

	using RealType          = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using ComplexType       = std::complex<RealType>;
	using VectorComplexType = typename PsimagLite::Vector<ComplexType>::Type;

	KBDerivative() = default;

	KBDerivative(SizeType nT, SizeType nTau)
	    : retarded(nT + 2)
	    , left_mixing(nTau + 1)
	    , lesser(nT + 2)
	{ }

	/*!
	 * \brief update
	 * Copy derivatives from src for all j up to n, and all imaginary-time indices.
	 */
	void update(SizeType n, SizeType nTau, const KBDerivative& src)
	{
		for (SizeType i = 0; i <= n; ++i) {
			retarded[i] = src.retarded[i];
			lesser[i]   = src.lesser[i];
		}
		for (SizeType j = 0; j <= nTau; ++j)
			left_mixing[j] = src.left_mixing[j];
	}

	VectorComplexType retarded;
	VectorComplexType left_mixing;
	VectorComplexType lesser;
};

} // namespace Dmft
#endif // KADANOFFBAYM_H
