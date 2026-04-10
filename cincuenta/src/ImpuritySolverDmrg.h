#ifndef IMPURITYSOLVER_DMRG_H
#define IMPURITYSOLVER_DMRG_H

#include "CmdLineOptions.hh"
#include "DmrgRunner.h"
#include "Geometry/Star.h"
#include "ImpuritySolverBase.h"
#include "InputCheck.h"
#include "LanczosSolver.h"
#include "ManyOmegas.h"
#include "Matsubaras.h"
#include "ModelParams.h"
#include "ParamsDmftSolver.h"
#include "ProcOmegas.h"
#include "PsiBase64.h"
#include "PsimagLite.h"
#include "Vector.h"

namespace Dmft {

template <typename ComplexOrRealType>
class ImpuritySolverDmrg : public ImpuritySolverBase<ComplexOrRealType> {

	enum class DmrgType
	{
		TYPE_0,
		TYPE_1
	};

public:

	using BaseType             = ImpuritySolverBase<ComplexOrRealType>;
	using ParamsDmftSolverType = ParamsDmftSolver<ComplexOrRealType>;
	using RealType             = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using ComplexType          = std::complex<RealType>;
	using VectorRealType       = typename PsimagLite::Vector<RealType>::Type;
	using VectorComplexType    = typename PsimagLite::Vector<ComplexType>::Type;
	using DmrgRunnerType       = Dmrg::DmrgRunner<RealType>;
	using ApplicationType      = PsimagLite::PsiApp;
	using MatsubarasType       = PsimagLite::Matsubaras<RealType>;
	using ManyOmegasType       = Dmrg::ManyOmegas<RealType, MatsubarasType>;
	using ProcOmegasType       = Dmrg::ProcOmegas<RealType, MatsubarasType>;
	using ModelParamsType      = typename BaseType::ModelParamsType;
	using InputNgType          = typename BaseType::InputNgType;

	ImpuritySolverDmrg(const ParamsDmftSolverType&     params,
	                   const ApplicationType&          app,
	                   typename InputNgType::Readable& io)
	    : params_(params)
	    , app_(app)
	    , io_(io)
	{ }

	// bathParams[0-nBath-1] ==> V ==> hoppings impurity --> bath
	// bathParams[nBath-...] ==> energies on each bath site
	void solve(const VectorRealType& bathParams)
	{
		ModelParamsType model_params(bathParams, io_);
		SizeType        mpiRank = PsimagLite::MPI::commRank(PsimagLite::MPI::COMM_WORLD);

		if (mpiRank == 0) {
			PsimagLite::String data2
			    = BaseType::readAndModifyInput(params_.gsTemplate, model_params);
			PsimagLite::String insitu = "<gs|nup|gs>";

			Dmrg::CmdLineOptions cmdline_options;
			cmdline_options.in_situ_measurements = "<gs|nup|gs>";
			cmdline_options.logfile              = "-";

			DmrgRunnerType runner(app_, data2, cmdline_options);
			runner.doOneRun();
		}

		PsimagLite::MPI::barrier(PsimagLite::MPI::COMM_WORLD);

		PsimagLite::String data3;
		InputNgType::Writeable::readFile(data3, params_.omegaTemplate);
		PsimagLite::String data4 = BaseType::addBathParams(data3, model_params);

		SizeType impurity_site = model_params.impuritySite();

		doType(DmrgType::TYPE_0, data4, impurity_site, mpiRank);

		doType(DmrgType::TYPE_1, data4, impurity_site, mpiRank);

		if (mpiRank == 0) {
			scaleGimp();

			std::cerr << "Sum of Gimp=" << density() << "\n";
			MatsubarasType matsubaras(params_.ficticiousBeta, params_.nMatsubaras, 0.);
			BaseType::writeGimpForDebugOnly("gimp_dmrg.txt", gimp_, matsubaras);
		}

		PsimagLite::MPI::barrier(PsimagLite::MPI::COMM_WORLD);
	}

	const VectorComplexType& gimp() const { return gimp_; }

private:

	/*
	 * This business of communicating by strings is far from ideal,
	 * but DMRG++ doesn't have an internal API right now
	 */
	static std::string
	addTypeAndObs(DmrgType t, SizeType impurity_site, PsimagLite::String data)
	{
		const std::string obsTc = (t == DmrgType::TYPE_0) ? "c'" : "c";
		return data + addType(t) + addObs(impurity_site, obsTc);
	}

	/* The type added is because there are two terms that need to be computed
	 * by separate runs and the terms differ by a sign
	 */
	static std::string addType(DmrgType t)
	{

		const SizeType tt = (t == DmrgType::TYPE_0) ? 0 : 1;
		return "DynamicDmrgType=" + ttos(tt) + ";\n";
	}

	/*
	 * The observable is just one for the Green function: c (the destruction operator)
	 * and we apply at the impurity. Now, if the impurity is site 0 (a border site)
	 * we need to add a trigger at site 1, that is, and operator at site 1: the identity,
	 * that is we multiply by the identity in this case.
	 */
	static std::string addObs(SizeType impurity_site, const std::string& obs)
	{
		std::string str = "TSPProductOrSum=product;\n";
		SizeType    ind = 0;
		if (impurity_site > 0) {
			str += "TSPSites=[" + ttos(impurity_site)
			    + "];\n"
			      " TSPLoops=[0]\n";
			ind = 0;
		} else {
			str += "TSPSites=[1, 0];\n"
			       "TSPLoops=[0, 0];\n"
			       "string TSPOp0:TSPOperator=expression;\n"
			       "string TSPOp0:OperatorExpression=identity;\n"
			       "string TSPOp1:TSPOperator=expression;\n";
			ind = 1;
		}

		str += "string TSPOp" + ttos(ind) + ":OperatorExpression=" + obs + ";\n";
		return str;
	}

	void doType(DmrgType t, PsimagLite::String data, SizeType impurity_site, SizeType mpiRank)
	{
		PsimagLite::String obs   = (t == DmrgType::TYPE_0) ? "c" : "c'";
		PsimagLite::String data2 = addTypeAndObs(t, impurity_site, data);

		PsimagLite::Matsubaras<RealType> matsubaras(params_.ficticiousBeta,
		                                            params_.nMatsubaras,
		                                            0.); // last argument is real part

		ManyOmegasType manyOmegas(data2, matsubaras, app_);

		const bool               dryrun   = false;
		const PsimagLite::String rootname = "dmftDynamics";
		Dmrg::CmdLineOptions     cmdline_options;
		cmdline_options.in_situ_measurements = "<gs|" + obs + "|P2>,<gs|" + obs + "|P3>";
		manyOmegas.run(dryrun, rootname, cmdline_options);

		if (mpiRank != 0)
			return;

		const PsimagLite::String rootIname   = "input";
		const PsimagLite::String rootOname   = "OUTPUT";
		const bool               skipFourier = true;

		Dmrg::InputCheck                inputCheck;
		typename InputNgType::Writeable ioW(inputCheck, data2);
		typename InputNgType::Readable  io(ioW);
		ProcOmegasType                  procOmegas(
                    io, params_.precision, skipFourier, rootIname, rootOname, matsubaras);

		procOmegas.run();

		readGimp(rootOname, matsubaras, t);
	}

	void readGimp(PsimagLite::String filename, const MatsubarasType& matsubaras, DmrgType t)
	{
		std::ifstream fin(filename);
		if (!fin || !fin.good() || fin.bad())
			err("readGimp: Could not open " + filename + "\n");

		if (gimp_.size() == 0)
			gimp_.resize(matsubaras.total());

		SizeType ind = 0;
		while (!fin.eof()) {
			RealType val = 0;
			fin >> val;
			SizeType n = 0;
			fin >> n;

			bool     centerSeen = false;
			RealType val1       = 0;
			RealType val2       = 0;

			for (SizeType i = 0; i <= params_.center_site; ++i) {
				SizeType site = 0;
				fin >> site;

				fin >> val2;

				fin >> val1;

				if (site == params_.center_site) {
					centerSeen = true;
					break;
				}
			}

			if (!centerSeen)
				err("Internal error: center " + ttos(params_.center_site)
				    + " not seen, freq id = " + ttos(ind) + "\n");

			if (t == DmrgType::TYPE_0) {
				gimp_[ind] = ComplexType(-val1, -val2);
			} else {
				gimp_[ind] += ComplexType(-val1, -val2);
			}

			++ind;

			if (ind >= gimp_.size())
				break;

			if (n == 1)
				continue;

			const SizeType tmp = params_.center_site + 1;
			n -= tmp;

			for (SizeType i = 0; i < 3 * n; ++i)
				fin >> val1;
		}

		if (ind < gimp_.size())
			err("readGimp: Not all values computed\n");
	}

	void scaleGimp()
	{
		const SizeType n           = gimp_.size();
		ComplexType    pre_density = density();
		const RealType factor
		    = std::sqrt(std::real(1.0 / pre_density / std::conj(pre_density)));
		for (SizeType i = 0; i < n; ++i) {
			gimp_[i] *= factor;
		}
	}

	ComplexType density() const
	{
		const SizeType n   = gimp_.size();
		ComplexType    sum = 0;
		for (SizeType i = 0; i < n; ++i)
			sum += gimp_[i];

		return sum;
	}

	const ParamsDmftSolverType&     params_;
	const ApplicationType&          app_;
	VectorComplexType               gimp_;
	typename InputNgType::Readable& io_;
};
}
#endif // IMPURITYSOLVER_DMRG_H
