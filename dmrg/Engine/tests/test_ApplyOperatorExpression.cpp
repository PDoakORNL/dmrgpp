#include "LastKrylovSlots.h"
#include "PsimagLite.h"
#include <catch2/catch_test_macros.hpp>

// LastKrylovSlots is the registry extracted from ApplyOperatorExpression that
// maps a P-vector's primary slot index to the last Krylov slot allocated for
// its time evolution.  ApplyOperatorExpression delegates its
// registerLastKrylovSlot / getLastKrylovSlot methods to this class.

using Dmrg::LastKrylovSlots;

TEST_CASE("LastKrylovSlots returns -1 for unregistered index", "[LastKrylovSlots]")
{
	LastKrylovSlots reg;
	REQUIRE(reg.getSlot(0) == -1);
	REQUIRE(reg.getSlot(5) == -1);
}

TEST_CASE("LastKrylovSlots registers and retrieves a single slot", "[LastKrylovSlots]")
{
	LastKrylovSlots reg;
	reg.registerSlot(2, 7);
	REQUIRE(reg.getSlot(2) == 7);
}

TEST_CASE("LastKrylovSlots handles multiple P-vectors independently", "[LastKrylovSlots]")
{
	LastKrylovSlots reg;
	reg.registerSlot(0, 4);
	reg.registerSlot(1, 8);
	reg.registerSlot(3, 12);

	REQUIRE(reg.getSlot(0) == 4);
	REQUIRE(reg.getSlot(1) == 8);
	REQUIRE(reg.getSlot(2) == -1); // never registered
	REQUIRE(reg.getSlot(3) == 12);
}

TEST_CASE("LastKrylovSlots overwrites on re-registration", "[LastKrylovSlots]")
{
	LastKrylovSlots reg;
	reg.registerSlot(0, 4);
	reg.registerSlot(0, 9);
	REQUIRE(reg.getSlot(0) == 9);
}
