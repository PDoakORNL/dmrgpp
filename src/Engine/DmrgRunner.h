#ifndef DMRGRUNNER_H
#define DMRGRUNNER_H
#include "BasisWithOperators.h"
#include "CmdLineOptions.hh"
#include "CrsMatrix.h"
#include "DmrgSolver.h"
#include "InputCheck.h"
#include "LeftRightSuper.h"
#include "MatrixVectorKron/MatrixVectorKron.h"
#include "MatrixVectorOnTheFly.h"
#include "MatrixVectorStored.h"
#include "ModelHelperLocal.h"
#include "ModelSelector.h"
#include "ParametersDmrgSolver.h"
#include "ProgramGlobals.h"
#include "PsimagLite.h"
#include "Qn.h"
#include "RedirectOutput.hh"
#include "RunFinished.h"
#include "SuperGeometry.h"
#include "VectorWithOffset.h"

namespace Dmrg {

template <typename RealType> class DmrgRunner {

public:

	using InputNgType              = PsimagLite::InputNg<InputCheck>;
	using ParametersDmrgSolverType = ParametersDmrgSolver<RealType, InputNgType::Readable, Qn>;
	using ApplicationType          = PsimagLite::PsiApp;

	DmrgRunner(const ApplicationType& application)
	    : application_(application)
	{ }

	void doOneRun(const PsimagLite::String& data, const CmdLineOptions& cmd_line) const
	{
		InputCheck             inputCheck;
		InputNgType::Writeable ioWriteable(inputCheck, data);
		InputNgType::Readable  io(ioWriteable);

		ParametersDmrgSolverType dmrgSolverParams(io, cmd_line.solver_options);

		adjustDmrgSolverParams(dmrgSolverParams,
		                       cmd_line.in_situ_measurements,
		                       cmd_line.precision,
		                       cmd_line.number_of_threads);

		if (endDueToClobber(dmrgSolverParams.options.isSet("noClobber"),
		                    cmd_line.logfile)) {
			return;
		}

		dealWithConsoleOutput(
		    dmrgSolverParams, cmd_line.logfile, cmd_line.unbuffered_output);

		constexpr bool PRINT_HEADER = true;
		application_.base64encode(std::cout, data, PRINT_HEADER);

		application_.printCmdLine(std::cout);

		setCorrectThreading(dmrgSolverParams, io);

		bool isComplex = (dmrgSolverParams.options.isSet("useComplex")
		                  || dmrgSolverParams.options.isSet("TimeStepTargeting"));

		if (isComplex)
			doOneRun2<std::complex<RealType>>(dmrgSolverParams, io);
		else
			doOneRun2<RealType>(dmrgSolverParams, io);
	}

	void printOutputChange(PsimagLite::String logfile, PsimagLite::String data) const
	{
		std::cerr << Provenance::logo(application_.name());
		std::cerr << "Standard output sent to ";
		std::cerr << logfile << "\n";
		std::cerr.flush();

		application_.printCmdLine(std::cout);
		application_.base64encode(std::cout, data, true);
	}

	const ApplicationType& application() const { return application_; }

private:

	template <typename ComplexOrRealType>
	void doOneRun2(const ParametersDmrgSolverType& dmrgSolverParams,
	               InputNgType::Readable&          io) const
	{
		using SuperGeometryType
		    = SuperGeometry<ComplexOrRealType, InputNgType::Readable, ProgramGlobals>;
		using VectorWithOffsetType   = VectorWithOffset<ComplexOrRealType, Qn>;
		using MySparseMatrixComplex  = PsimagLite::CrsMatrix<ComplexOrRealType>;
		using BasisType              = Basis<MySparseMatrixComplex>;
		using BasisWithOperatorsType = BasisWithOperators<BasisType>;
		using LeftRightSuperType     = LeftRightSuper<BasisWithOperatorsType, BasisType>;
		using ModelHelperType        = ModelHelperLocal<LeftRightSuperType>;
		using ModelBaseType          = ModelBase<ModelHelperType,
		                                         ParametersDmrgSolverType,
		                                         InputNgType::Readable,
		                                         SuperGeometryType>;

		if (dmrgSolverParams.options.isSet("MatrixVectorStored")) {
			doOneRun3<MatrixVectorStored<ModelBaseType>>(dmrgSolverParams, io);
		} else if (dmrgSolverParams.options.isSet("MatrixVectorOnTheFly")) {
			doOneRun3<MatrixVectorOnTheFly<ModelBaseType>>(dmrgSolverParams, io);
		} else {
			doOneRun3<MatrixVectorKron<ModelBaseType>>(dmrgSolverParams, io);
		}
	}

	template <typename MatrixVectorType>
	void doOneRun3(const ParametersDmrgSolverType& dmrgSolverParams,
	               InputNgType::Readable&          io) const
	{
		using ComplexOrRealType = typename MatrixVectorType::ComplexOrRealType;
		if (dmrgSolverParams.options.isSet("vectorwithoffsets")) {
			typedef VectorWithOffsets<ComplexOrRealType, Qn> VectorWithOffsetType;
			doOneRun4<MatrixVectorType, VectorWithOffsetType>(dmrgSolverParams, io);
		} else {
			typedef VectorWithOffset<ComplexOrRealType, Qn> VectorWithOffsetType;
			doOneRun4<MatrixVectorType, VectorWithOffsetType>(dmrgSolverParams, io);
		}
	}
	template <typename MatrixVectorType, typename VectorWithOffsetType>
	void doOneRun4(const ParametersDmrgSolverType& dmrgSolverParams,
	               InputNgType::Readable&          io) const
	{
		using ComplexOrRealType = typename MatrixVectorType::ComplexOrRealType;
		using SuperGeometryType
		    = SuperGeometry<ComplexOrRealType, InputNgType::Readable, ProgramGlobals>;

		SuperGeometryType geometry(io);
		if (dmrgSolverParams.options.isSet("printgeometry"))
			std::cout << geometry;

		using SolverType    = PsimagLite::LanczosSolver<MatrixVectorType>;
		using ModelBaseType = typename SolverType::MatrixType::ModelType;

		//! Setup the Model
		ModelSelector<ModelBaseType> modelSelector(dmrgSolverParams.model);
		ModelBaseType&               model = modelSelector(dmrgSolverParams, io, geometry);

		//! Setup the dmrg solver
		using DmrgSolverType = DmrgSolver<SolverType, VectorWithOffsetType>;
		DmrgSolverType dmrgSolver(model, io);

		//! Calculate observables:
		dmrgSolver.main(geometry);

		std::cout.flush();
	}

	void dealWithConsoleOutput(const ParametersDmrgSolverType& dmrgSolverParams,
	                           const std::string&              logfile,
	                           bool                            unbuffered) const
	{
		if (logfile == "-") {
			return;
		}

		PsimagLite::RedirectOutput::setAppName(application_.name(),
		                                       Provenance::logo(application_.name()));

		std::ios_base::openmode open_mode
		    = (dmrgSolverParams.autoRestart) ? std::ofstream::app : std::ofstream::out;

		PsimagLite::RedirectOutput::doIt(logfile, open_mode, unbuffered);
	}

	void adjustDmrgSolverParams(ParametersDmrgSolverType& dmrgSolverParams,
	                            const std::string&        insitu,
	                            SizeType                  precision,
	                            SizeType                  threadsInCmd) const
	{
		if (threadsInCmd > 0)
			dmrgSolverParams.nthreads = threadsInCmd;

		if (precision > 0)
			dmrgSolverParams.precision = precision;

		dmrgSolverParams.insitu = insitu;

		if (dmrgSolverParams.options.isSet("minimizeDisk")) {
			dmrgSolverParams.options += ",noSaveWft,noSaveStacks,noSaveData";
		}

		// The below doesn't change dmrgSolverParams
		if (dmrgSolverParams.options.isSet("hd5DontPrint"))
			PsimagLite::IoNg::dontPrintDebug();

		if (dmrgSolverParams.autoRestart) {
			std::cout << "\nAutoRestart possible\n";
		}
	}

	static void setCorrectThreading(const ParametersDmrgSolverType& dmrgSolverParams,
	                                InputNgType::Readable&          io)
	{
		typedef PsimagLite::Concurrency ConcurrencyType;

		SizeType threadsStackSize = 0;
		try {
			io.readline(threadsStackSize, "ThreadsStackSize=");
		} catch (std::exception&) { }

		constexpr bool setAffinities = false; // no longer supported, I mean, deleted

		PsimagLite::CodeSectionParams codeSection(dmrgSolverParams.nthreads,
		                                          dmrgSolverParams.nthreads2,
		                                          setAffinities,
		                                          threadsStackSize);
		ConcurrencyType::setOptions(codeSection);
	}

	static bool endDueToClobber(bool is_no_clobber_set, const std::string& logfile)
	{
		if (logfile == "-")
			return false;

		RunFinished runFinished(is_no_clobber_set);
		if (runFinished.OK(logfile.c_str())) {
			runFinished.printTermination(logfile.c_str());
			return true;
		}

		return false;
	}

	const ApplicationType& application_;
};
}
#endif // DMRGRUNNER_H
