#include "GetBraOrKet.h"
#include "PsimagLite.h"
#include <catch2/catch_test_macros.hpp>

using G = PsimagLite::GetBraOrKet;

TEST_CASE("GetBraOrKet ket vs bra detection", "[GetBraOrKet]")
{
	REQUIRE(G("|gs>").isKet() == true);
	REQUIRE(G("<gs|").isKet() == false);
	REQUIRE(G("|P3>").isKet() == true);
	REQUIRE(G("<P3|").isKet() == false);
}

TEST_CASE("GetBraOrKet ground state", "[GetBraOrKet]")
{
	G gs("|gs>");
	REQUIRE(gs.isPvector() == false);
	REQUIRE(gs.isRvector() == false);
	REQUIRE(gs.isLastKrylov() == false);
}

TEST_CASE("GetBraOrKet P-vector index parsing", "[GetBraOrKet]")
{
	SECTION("P0")
	{
		G p0("|P0>");
		REQUIRE(p0.isPvector() == true);
		REQUIRE(p0.pIndex() == 0);
		REQUIRE(p0.isLastKrylov() == false);
	}

	SECTION("P2")
	{
		G p2("|P2>");
		REQUIRE(p2.isPvector() == true);
		REQUIRE(p2.pIndex() == 2);
		REQUIRE(p2.isLastKrylov() == false);
	}

	SECTION("P12 multi-digit")
	{
		G p12("|P12>");
		REQUIRE(p12.isPvector() == true);
		REQUIRE(p12.pIndex() == 12);
		REQUIRE(p12.isLastKrylov() == false);
	}

	SECTION("bra form <P2|")
	{
		G p2bra("<P2|");
		REQUIRE(p2bra.isPvector() == true);
		REQUIRE(p2bra.pIndex() == 2);
		REQUIRE(p2bra.isKet() == false);
	}
}

TEST_CASE("GetBraOrKet .last suffix", "[GetBraOrKet]")
{
	SECTION("ket |P2.last>")
	{
		G p2last("|P2.last>");
		REQUIRE(p2last.isPvector() == true);
		REQUIRE(p2last.pIndex() == 2);
		REQUIRE(p2last.isLastKrylov() == true);
		REQUIRE(p2last.isKet() == true);
	}

	SECTION("bra <P2.last|")
	{
		G p2lastBra("<P2.last|");
		REQUIRE(p2lastBra.isPvector() == true);
		REQUIRE(p2lastBra.pIndex() == 2);
		REQUIRE(p2lastBra.isLastKrylov() == true);
		REQUIRE(p2lastBra.isKet() == false);
	}

	SECTION("P0.last")
	{
		G p0last("|P0.last>");
		REQUIRE(p0last.pIndex() == 0);
		REQUIRE(p0last.isLastKrylov() == true);
	}

	SECTION("regular P2 does not have .last set")
	{
		G p2("|P2>");
		REQUIRE(p2.isLastKrylov() == false);
	}
}
