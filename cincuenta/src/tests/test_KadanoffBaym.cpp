#include "KadanoffBaym.h"
#include <PsimagLite/PsimagLite.h>
#include <catch2/catch_test_macros.hpp>

using Dmft::KadanoffBaym;
using Dmft::KBDerivative;

TEST_CASE("KadanoffBaym dimensions after construction", "[KadanoffBaym]")
{
	const SizeType       nT = 5, nTau = 3;
	KadanoffBaym<double> kb(nT, nTau, 0.1, 0.5);

	CHECK(kb.matsubara_w.size() == nTau);
	CHECK(kb.matsubara_t.size() == nTau + 1);
	CHECK(SizeType(kb.retarded.rows()) == nT + 1);
	CHECK(SizeType(kb.retarded.cols()) == nT + 1);
	CHECK(SizeType(kb.left_mixing.rows()) == nT + 1);
	CHECK(SizeType(kb.left_mixing.cols()) == nTau + 1);
	CHECK(SizeType(kb.lesser.rows()) == nT + 1);
	CHECK(SizeType(kb.lesser.cols()) == nT + 1);
}

TEST_CASE("KadanoffBaym accessors return constructor arguments", "[KadanoffBaym]")
{
	KadanoffBaym<double> kb(4, 6, 0.05, 0.25);
	CHECK(kb.nT() == 4);
	CHECK(kb.nTau() == 6);
	CHECK(kb.dt() == 0.05);
	CHECK(kb.dtau() == 0.25);
}

TEST_CASE("KadanoffBaym default-constructed has empty containers", "[KadanoffBaym]")
{
	KadanoffBaym<double> kb;
	CHECK(kb.nT() == 0);
	CHECK(kb.nTau() == 0);
	CHECK(kb.matsubara_w.size() == 0);
	CHECK(kb.matsubara_t.size() == 0);
}

TEST_CASE("KBDerivative dimensions after construction", "[KadanoffBaym]")
{
	const SizeType       nT = 7, nTau = 4;
	KBDerivative<double> der(nT, nTau);

	CHECK(der.retarded.size() == nT + 2);
	CHECK(der.lesser.size() == nT + 2);
	CHECK(der.left_mixing.size() == nTau + 1);
}

TEST_CASE("KBDerivative default-constructed has empty containers", "[KadanoffBaym]")
{
	KBDerivative<double> der;
	CHECK(der.retarded.size() == 0);
	CHECK(der.lesser.size() == 0);
	CHECK(der.left_mixing.size() == 0);
}
