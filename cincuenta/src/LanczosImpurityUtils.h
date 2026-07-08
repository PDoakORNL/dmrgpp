#ifndef LANCZOS_IMPURITY_UTILS_H
#define LANCZOS_IMPURITY_UTILS_H

#include "InputCheck.h"
#include "LanczosPlusPlus/src/Engine/DefaultSymmetry.h"
#include "LanczosPlusPlus/src/Engine/InternalProductOnTheFly.h"
#include "LanczosPlusPlus/src/Engine/InternalProductStored.h"
#include "LanczosPlusPlus/src/Engine/LanczosGlobals.h"
#include "LanczosPlusPlus/src/Engine/LanczosModelBase.hpp"
#include "LanczosSolver.h"
#include "Matrix.h"
#include "ParametersForSolver.h"
#include "PsimagLite.h"
#include "Vector.h"

namespace Dmft {

// Shared Lanczos diagonalization and input-building utilities used by both
// ImpuritySolverEqLanczos (equilibrium) and ImpuritySolverNeqLanczos (NEQ).
template <typename ComplexOrRealType> struct LanczosImpurityUtils {

	using RealType       = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using VectorRealType = typename PsimagLite::Vector<RealType>::Type;
	using MatrixType     = PsimagLite::Matrix<ComplexOrRealType>;

	using DmrgInputReadable = typename PsimagLite::InputNg<Dmrg::InputCheck>::Readable;
	using GeometryType      = PsimagLite::
	    Geometry<ComplexOrRealType, DmrgInputReadable, LanczosPlusPlus::LanczosGlobals>;
	using ModelBaseType
	    = LanczosPlusPlus::LanczosModelBase<ComplexOrRealType, GeometryType, DmrgInputReadable>;
	using BasisBaseType       = typename ModelBaseType::BasisBaseType;
	using DefaultSymmetryType = LanczosPlusPlus::DefaultSymmetry<GeometryType, BasisBaseType>;
	using InternalProductOnTheFlyType
	    = LanczosPlusPlus::InternalProductOnTheFly<ModelBaseType, DefaultSymmetryType>;
	using InternalProductStoredType
	    = LanczosPlusPlus::InternalProductStored<ModelBaseType, DefaultSymmetryType>;
	using VectorType           = typename PsimagLite::Vector<ComplexOrRealType>::Type;
	using VectorVectorType     = typename PsimagLite::Vector<VectorType>::Type;
	using ParametersSolverType = PsimagLite::ParametersForSolver<RealType>;
	using LanczosSolverType    = PsimagLite::LanczosSolver<InternalProductOnTheFlyType>;

	// Truncated Lanczos diagonalization of the Hamiltonian in the given basis.
	// Returns up to min(nStates, dim) lowest eigenpairs.
	// eigvecs has shape (dim, nKept): eigvecs(j, i) = component j of eigenvector i.
	// If nStates == 0 or nStates >= dim, falls back to full diagonalization.
	// Also falls back to fullDiag if Lanczos throws.
	static void diagWithBasisTruncated(const ModelBaseType& model,
	                                   const BasisBaseType& basis,
	                                   const GeometryType&  geom,
	                                   SizeType             nStates,
	                                   VectorRealType&      eigs,
	                                   MatrixType&          eigvecs)
	{
		const SizeType dim   = basis.size();
		const SizeType nKeep = (nStates == 0 || nStates >= dim) ? dim : nStates;

		DefaultSymmetryType rs(basis, geom, "");

		if (nKeep == dim) {
			InternalProductStoredType hamS(model, basis, rs);
			VectorRealType            eigs2(dim);
			MatrixType                fm;
			hamS.fullDiag(eigs2, fm);
			eigs.resize(nKeep);
			eigvecs.resize(dim, nKeep);
			for (SizeType i = 0; i < nKeep; ++i) {
				eigs[i] = eigs2[i];
				for (SizeType j = 0; j < dim; ++j)
					eigvecs(j, i) = fm(j, i);
			}
			return;
		}

		InternalProductOnTheFlyType ham(model, basis, rs);
		VectorType                  initial(dim);
		PsimagLite::fillRandom(initial);
		VectorVectorType zs(nKeep, VectorType(dim));
		eigs.resize(nKeep);

		ParametersSolverType lparams;
		lparams.lotaMemory = true;
		lparams.options    = "reortho";
		LanczosSolverType lanczosSolver(ham, lparams);

		try {
			lanczosSolver.computeAllStatesBelow(eigs, zs, initial, nKeep);
			eigs.resize(nKeep);
		} catch (std::exception&) {
			InternalProductStoredType hamS(model, basis, rs);
			VectorRealType            eigs2(dim);
			MatrixType                fm;
			hamS.fullDiag(eigs2, fm);
			eigs.resize(nKeep);
			for (SizeType i = 0; i < nKeep; ++i) {
				eigs[i] = eigs2[i];
				for (SizeType j = 0; j < dim; ++j)
					zs[i][j] = fm(j, i);
			}
		}

		eigvecs.resize(dim, eigs.size());
		for (SizeType i = 0; i < eigs.size(); ++i)
			for (SizeType j = 0; j < dim; ++j)
				eigvecs(j, i) = zs[i][j];
	}

	// Build LanczosPlusPlus input string for star-geometry Anderson model.
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

		// A single site (no bath) has no bonds at all; the Ainur vector
		// grammar cannot parse an empty "[]" literal, so NumberOfTerms=0 is
		// used instead to skip the Connectors read entirely rather than
		// trying to serialize a zero-length vector.
		const bool  hasHoppings = (hoppings.size() > 0);
		std::string connStr     = "[";
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
		s += "NumberOfTerms=" + std::string(hasHoppings ? "1" : "0") + ";\n";
		s += "DegreesOfFreedom=1;\n";
		s += "GeometryKind=star;\n";
		s += "GeometryOptions=none;\n";
		s += "hubbardU=" + uStr + ";\n";
		s += "Model=HubbardOneBand;\n";
		s += "SolverOptions=twositedmrg,geometryallinsystem,hd5dontprint;\n";
		s += "Version=templateForDMFT;\n";
		s += "OutputFile=neqDmrgDummy;\n";
		s += "InfiniteLoopKeptStates=1;\n";
		s += "FiniteLoops=0 0 0;\n";
		s += "TargetElectronsUp=" + ttos(nup) + ";\n";
		s += "TargetElectronsDown=" + ttos(ndown) + ";\n";
		if (hasHoppings)
			s += "dir0:Connectors=" + connStr + ";\n";
		s += "potentialV=" + potStr + ";\n";
		return s;
	}
};

} // namespace Dmft
#endif // LANCZOS_IMPURITY_UTILS_H
