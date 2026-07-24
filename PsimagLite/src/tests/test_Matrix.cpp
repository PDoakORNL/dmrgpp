#include <PsimagLite/Matrix.h>
#include <catch2/catch_test_macros.hpp>

using MatrixType = PsimagLite::Matrix<int>;

// Locks down resize(nrow, ncol, val)'s documented contract: NOT
// element-preserving across a change in nrow. Found via a real bug: a
// process-lifetime cache (LanczosPlusPlus::Combinatorial) relied on this
// overload leaving untouched cells at their old value across repeated
// resizes to different sizes, and instead got a different, unrelated old
// cell's value reinterpreted at the new (row, col) position.
TEST_CASE("Matrix resize(nrow,ncol,val) is not element-preserving across nrow changes",
          "[Matrix][resize]")
{
	MatrixType m(3, 3);
	for (SizeType i = 0; i < 3; ++i)
		for (SizeType j = 0; j < 3; ++j)
			m(i, j) = static_cast<int>(10 * i + j);

	// (0,1) is 1 before the resize.
	CHECK(m(0, 1) == 1);

	// Shrink to 2x2 with nrow changing (3 -> 2). Storage is column-major
	// (flat index = i + j*nrow_), so this is a raw flat-array truncation,
	// not a reshape: new (0,1)'s flat index is 0 + 1*2 = 2, which under
	// the OLD nrow=3 layout was flat index 2, i.e. old (2,0) = 20 -- a
	// cell unrelated to (0,1) by row/col identity.
	m.resize(2, 2, -1);
	CHECK(m(0, 1) == 20); // old (2,0), NOT the fill value -1 and NOT old (0,1)
}

// Locks down resize(nrow, ncol)'s documented contract: element-preserving
// and shape-aware. Every (i,j) present in both shapes keeps its value;
// newly-added (i,j) are value-initialized (0 for int). This is the
// resize to use whenever a reshape must not silently reinterpret data
// from a different logical cell -- e.g. Combinatorial's cache, fixed by
// switching to this overload.
TEST_CASE("Matrix resize(nrow,ncol) preserves data at true (row,col) identity", "[Matrix][resize]")
{
	MatrixType m(3, 3);
	for (SizeType i = 0; i < 3; ++i)
		for (SizeType j = 0; j < 3; ++j)
			m(i, j) = static_cast<int>(10 * i + j);

	SECTION("Shrinking nrow keeps each surviving (i,j) at its own value")
	{
		m.resize(2, 2);
		CHECK(m(0, 0) == 0);
		CHECK(m(0, 1) == 1);
		CHECK(m(1, 0) == 10);
		CHECK(m(1, 1) == 11);
	}

	SECTION("Growing nrow keeps each old (i,j) at its own value; new cells are "
	        "value-initialized")
	{
		m.resize(5, 4);
		for (SizeType i = 0; i < 3; ++i)
			for (SizeType j = 0; j < 3; ++j)
				CHECK(m(i, j) == static_cast<int>(10 * i + j));
		CHECK(m(3, 0) == 0);
		CHECK(m(4, 3) == 0);
		CHECK(m(0, 3) == 0);
	}

	SECTION("Repeated resize to different nrow never resurrects an unrelated old cell")
	{
		// This is the exact shape of the real bug: grow, then shrink back
		// down to a size that was NOT the original size, and check a cell
		// that was never explicitly written keeps its value-initialized
		// default instead of inheriting some other cell's old content.
		MatrixType fresh(1, 1);
		fresh.resize(12, 12); // grow: (0,1) is new -> 0
		CHECK(fresh(0, 1) == 0);
		fresh(8, 0) = 99; // an unrelated cell, far from (0,1)
		fresh.resize(8, 8); // shrink: (0,1) must stay 0, not become 99
		CHECK(fresh(0, 1) == 0);
	}
}
