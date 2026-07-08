#include "CincuentaInputCheck.h"
#include "ImpuritySolverNeqExactDiag.h"
#include "KadanoffBaym.h"
#include "LanczosPlusPlus/src/Engine/InputCheck.h"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <complex>
#include <iostream>
#include <string>

using RealType       = double;
using ComplexType    = std::complex<RealType>;
using SolverType     = Dmft::ImpuritySolverNeqExactDiag<ComplexType>;
using ParamsType     = Dmft::ParamsNeqDmftSolver<ComplexType>;
using InputNgType    = PsimagLite::InputNg<Dmft::CincuentaInputCheck>;
using KBType         = typename SolverType::KBType;
using VectorRealType = typename SolverType::VectorRealType;

// Variables pre-declared by CincuentaInputCheck::import() (which chains
// LanczosPlusPlus::InputCheck::import() and Dmrg::InputCheck::import()) must be
// assigned without a type keyword; re-declaring with "integer"/"real" triggers
// Ainur's "Already declared".  HubbardU uses "real" because import() only knows
// the lowercase-vector form "hubbardU".
static const std::string kConfig = "##Ainur1.0\n\n"
                                   "FicticiousBeta=10;\n"
                                   "ChemicalPotential=0.;\n"
                                   "Matsubaras=20;\n"
                                   "LatticeGf=\"energy,semicircular,4\";\n"
                                   "NumberOfBathPoints=1;\n"
                                   "DmftNumberOfIterations=1;\n"
                                   "DmftTolerance=1e-6;\n"
                                   "ImpuritySolver=\"exactdiag\";\n"
                                   "MinParamsDelta=0.01;\n"
                                   "MinParamsDelta2=0.01;\n"
                                   "MinParamsTolerance=1e-4;\n"
                                   "MinParamsMaxIter=100;\n"
                                   "MinParamsVerbose=0;\n"
                                   "real HubbardU=0.;\n"
                                   "HubbardUFinal=0.;\n"
                                   "TmaxNeq=1.0;\n"
                                   "NtNeq=4;\n"
                                   "TargetElectronsUp=1;\n"
                                   "TargetElectronsDown=1;\n";

/*!
 * \brief Construct a zero-initialised KadanoffBaym grid sized to match the solved params.
 * \param[in] solverParams Solver parameters supplying grid dimensions, time step, and beta.
 * \return A fresh zero-initialised Keldysh grid ready for computeGimp.
 */
static KBType makeSlice(const ParamsType& solverParams)
{
	return KBType(solverParams.nT,
	              solverParams.eqParams.nMatsubaras,
	              solverParams.dt,
	              solverParams.eqParams.ficticiousBeta
	                  / static_cast<RealType>(solverParams.eqParams.nMatsubaras));
}

// ---- buildLanczosInput -------------------------------------------------

TEST_CASE("buildLanczosInput produces valid LanczosPlusPlus input string", "[buildLanczosInput]")
{
	// Use distinct, non-trivial values so each key-value mapping can be verified
	// independently rather than matching an accidental coincidence.
	const RealType       U        = 2.0;
	const SizeType       nup      = 1;
	const SizeType       ndown    = 1;
	const VectorRealType hoppings = { 0.5 };
	const VectorRealType potV     = { -0.5, 0.3 }; // impurity site, bath site
	const SizeType       nsites   = 2;

	const std::string inputStr
	    = SolverType::buildLanczosInput(U, nup, ndown, hoppings, potV, nsites);

	std::cout << "=== buildLanczosInput output ===\n"
	          << inputStr << "================================\n";

	// Parse the string via LanczosPlusPlus::InputCheck to verify it is valid
	// LanczosPlusPlus input and that every argument reaches the correct key.
	using LppInputNg = PsimagLite::InputNg<LanczosPlusPlus::InputCheck>;
	LanczosPlusPlus::InputCheck lppInputCheck;
	LppInputNg::Writeable       lppIoW(lppInputCheck, inputStr);
	LppInputNg::Readable        lppIo(lppIoW);

	SizeType totalSites = 0;
	lppIo.readline(totalSites, "TotalNumberOfSites=");
	CHECK(totalSites == nsites);

	SizeType targetUp = 0;
	lppIo.readline(targetUp, "TargetElectronsUp=");
	CHECK(targetUp == nup);

	SizeType targetDown = 0;
	lppIo.readline(targetDown, "TargetElectronsDown=");
	CHECK(targetDown == ndown);

	// hubbardU: U at impurity (index 0), zero at every bath site.
	VectorRealType hubbardUVec;
	lppIo.read(hubbardUVec, "hubbardU");
	REQUIRE(hubbardUVec.size() == nsites);
	CHECK(hubbardUVec[0] == Catch::Approx(U));
	CHECK(hubbardUVec[1] == Catch::Approx(0.0));

	// dir0:Connectors: bath hopping amplitudes in order.
	VectorRealType connectors;
	lppIo.read(connectors, "dir0:Connectors");
	REQUIRE(connectors.size() == hoppings.size());
	CHECK(connectors[0] == Catch::Approx(hoppings[0]));

	// potentialV: spin-up block (sites 0..nsites-1) then spin-down block (sites
	// nsites..2*nsites-1).
	VectorRealType potentialV;
	lppIo.read(potentialV, "potentialV");
	REQUIRE(potentialV.size() == 2 * nsites);
	for (SizeType siteIdx = 0; siteIdx < nsites; ++siteIdx) {
		CHECK(potentialV[siteIdx] == Catch::Approx(potV[siteIdx]));
		CHECK(potentialV[siteIdx + nsites] == Catch::Approx(potV[siteIdx]));
	}
}

// ---- Keldysh GF boundary conditions -----------------------------------

TEST_CASE("ImpuritySolverNeqExactDiag Keldysh GF at U=0 half-filling", "[GreensFunctions]")
{
	// U=0, 1-bath site, half-filling: nup=ndown=1, hopping=0.5, epsilon_bath=0.
	// All GF properties below follow analytically from the free propagator.
	InputNgType::Writeable ioW(Dmft::CincuentaInputCheck {}, kConfig);
	InputNgType::Readable  io(ioW);
	ParamsType             params(io);
	SolverType             solver(params, io);

	const VectorRealType bathParams = { 0.5, 0.0 };
	solver.solve(bathParams);

	SECTION("G^R(t,t) = -i (parameter-independent boundary condition)")
	{
		KBType slice = makeSlice(params);
		solver.computeGimp(slice, 0);
		CHECK(slice.retarded(0, 0).real() == Catch::Approx(0.0).epsilon(1e-12));
		CHECK(slice.retarded(0, 0).imag() == Catch::Approx(-1.0).epsilon(1e-12));
	}

	SECTION("Im[G^<(0,0)] = 0.5 (half-filling occupancy)")
	{
		KBType slice = makeSlice(params);
		solver.computeGimp(slice, 0);
		CHECK(slice.lesser(0, 0).imag() == Catch::Approx(0.5).epsilon(1e-8));
	}

	SECTION("Re[G^M(tau=0)] = -0.5 (Matsubara normalisation at half-filling)")
	{
		CHECK(solver.gimp().matsubara_t[0].real() == Catch::Approx(-0.5).epsilon(1e-6));
	}

	SECTION("G^R is causal: retarded(0,j) = 0 for j > 0")
	{
		KBType slice = makeSlice(params);
		for (int n = 0; n <= static_cast<int>(params.nT); ++n)
			solver.computeGimp(slice, n);
		CHECK(slice.retarded(0, 1).real() == Catch::Approx(0.0).epsilon(1e-15));
		CHECK(slice.retarded(0, 1).imag() == Catch::Approx(0.0).epsilon(1e-15));
		CHECK(slice.retarded(0, 2).real() == Catch::Approx(0.0).epsilon(1e-15));
		CHECK(slice.retarded(0, 2).imag() == Catch::Approx(0.0).epsilon(1e-15));
	}
}

// ---- NeqAtomicLimit (nBath=0) closed-form bypass -----------------------
//
// Expected values below were generated by
// cincuenta/TestSuite/gbek_reference/atomic_limit_reference.py (U=2, dt=0.04,
// FicticiousBeta=0.8, Matsubaras=20 so the Matsubara grid coincides with the
// real-time grid), which independently derives these closed forms from the
// T=0 Lehmann sum and cross-checks them to machine precision against a
// direct 4-state exact diagonalization -- see that file's docstring,
// especially the note that G^Left decays at rate U/4, not the naively
// expected gap U/2.
static const std::string kConfigAtomic = "##Ainur1.0\n\n"
                                         "FicticiousBeta=0.8;\n"
                                         "ChemicalPotential=0.;\n"
                                         "Matsubaras=20;\n"
                                         "LatticeGf=\"energy,semicircular,4\";\n"
                                         "NumberOfBathPoints=1;\n"
                                         "DmftNumberOfIterations=1;\n"
                                         "DmftTolerance=1e-6;\n"
                                         "ImpuritySolver=\"exactdiag\";\n"
                                         "MinParamsDelta=0.01;\n"
                                         "MinParamsDelta2=0.01;\n"
                                         "MinParamsTolerance=1e-4;\n"
                                         "MinParamsMaxIter=100;\n"
                                         "MinParamsVerbose=0;\n"
                                         "real HubbardU=2.;\n"
                                         "HubbardUFinal=2.;\n"
                                         "TmaxNeq=0.8;\n"
                                         "NtNeq=20;\n"
                                         "TargetElectronsUp=1;\n"
                                         "TargetElectronsDown=0;\n";

TEST_CASE("ImpuritySolverNeqExactDiag NeqAtomicLimit (nBath=0) closed form",
          "[GreensFunctions][AtomicLimit]")
{
	InputNgType::Writeable ioW(Dmft::CincuentaInputCheck {}, kConfigAtomic);
	InputNgType::Readable  io(ioW);
	ParamsType             params(io);
	SolverType             solver(params, io);

	const VectorRealType emptyBathParams; // nBath=0 triggers the atomic-limit bypass
	solver.solve(emptyBathParams);

	KBType slice = makeSlice(params);
	for (int n = 0; n <= static_cast<int>(params.nT); ++n)
		solver.computeGimp(slice, n);

	SECTION("G^R(t_n,t_j) = -i*cos(U*(t_n-t_j)/2)")
	{
		CHECK(slice.retarded(0, 0).real() == Catch::Approx(0.0).epsilon(1e-12));
		CHECK(slice.retarded(0, 0).imag() == Catch::Approx(-1.0).epsilon(1e-12));
		CHECK(slice.retarded(1, 0).imag()
		      == Catch::Approx(-0.9992001066609779).epsilon(1e-10));
		CHECK(slice.retarded(5, 0).imag()
		      == Catch::Approx(-0.9800665778412416).epsilon(1e-10));
		CHECK(slice.retarded(20, 0).imag()
		      == Catch::Approx(-0.6967067093471654).epsilon(1e-10));
		CHECK(slice.retarded(20, 20).imag() == Catch::Approx(-1.0).epsilon(1e-12));
	}

	SECTION("G^<(t_n,t_j) = (i/2)*exp(i*U*(t_n-t_j)/2)")
	{
		CHECK(slice.lesser(0, 0).real() == Catch::Approx(0.0).epsilon(1e-12));
		CHECK(slice.lesser(0, 0).imag() == Catch::Approx(0.5).epsilon(1e-12));
		CHECK(slice.lesser(1, 0).real()
		      == Catch::Approx(-0.01999466709331708).epsilon(1e-10));
		CHECK(slice.lesser(1, 0).imag()
		      == Catch::Approx(0.49960005333048896).epsilon(1e-10));
		CHECK(slice.lesser(20, 0).real()
		      == Catch::Approx(-0.3586780454497614).epsilon(1e-10));
		CHECK(slice.lesser(20, 0).imag()
		      == Catch::Approx(0.3483533546735827).epsilon(1e-10));
	}

	SECTION("G^M(tau) = -(1/2)*exp(-(U/2)*tau)")
	{
		CHECK(solver.gimp().matsubara_t[0].real() == Catch::Approx(-0.5).epsilon(1e-10));
		CHECK(solver.gimp().matsubara_t[1].real()
		      == Catch::Approx(-0.4803947195761616).epsilon(1e-10));
		CHECK(solver.gimp().matsubara_t[5].real()
		      == Catch::Approx(-0.4093653765389909).epsilon(1e-10));
		CHECK(solver.gimp().matsubara_t[20].real()
		      == Catch::Approx(-0.22466448205861078).epsilon(1e-10));
	}

	SECTION("G^Left(t_n,tau_j) = -(i/2)*exp(i*U*t_n/4)*exp(-(U/4)*tau_j)")
	{
		CHECK(slice.left_mixing(0, 0).real() == Catch::Approx(0.0).epsilon(1e-12));
		CHECK(slice.left_mixing(0, 0).imag() == Catch::Approx(-0.5).epsilon(1e-12));
		CHECK(slice.left_mixing(1, 0).real()
		      == Catch::Approx(0.00999933334666654).epsilon(1e-10));
		CHECK(slice.left_mixing(1, 0).imag()
		      == Catch::Approx(-0.4999000033332889).epsilon(1e-10));
		CHECK(slice.left_mixing(10, 10).real()
		      == Catch::Approx(0.08132834540766957).epsilon(1e-10));
		CHECK(slice.left_mixing(10, 10).imag()
		      == Catch::Approx(-0.4012053236712601).epsilon(1e-10));
		CHECK(slice.left_mixing(20, 20).real()
		      == Catch::Approx(0.1305174605717285).epsilon(1e-10));
		CHECK(slice.left_mixing(20, 20).imag()
		      == Catch::Approx(-0.3087028239508228).epsilon(1e-10));
	}
}
