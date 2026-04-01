#include "BlockDiagonalMatrix.h"
#include "CmdLineOptions.hh"
#include "Concurrency.h"
#include "DmrgRunner.h"
#include "Io/IoNg.h"
#include "OptionsForIntrospect.hh"
#include "Provenance.h"

typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

using namespace Dmrg;

typedef PsimagLite::Concurrency ConcurrencyType;

void printLicense(const PsimagLite::PsiApp& app, const OperatorOptions& options)
{
	if (!ConcurrencyType::root() || options.enabled)
		return;

	std::cout << ProgramGlobals::license;
	Provenance provenance;
	std::cout << provenance;
	std::cout << Provenance::logo(app.name()) << "\n";
	app.checkMicroArch(std::cout, Provenance::compiledMicroArch());
}

void usageOperator()
{
	std::cerr << "USAGE is operator -f filename -e canonical_operator_expression\n";
	std::cerr << "Deprecated options are: -l label [-d dof] [-s site] [-t]\n";
}

int main(int argc, char** argv)
{
	PsimagLite::PsiApp application("DMRG++::instrospect", &argc, &argv, 1);
	PsimagLite::String filename = "";
	int                opt      = 0;
	OperatorOptions    operator_options;
	CmdLineOptions     cmdline_options;
	cmdline_options.solver_options = ",introspect";
	PsimagLite::String strUsage(application.name());
	strUsage += " -f filename [-p precision] [-s site] [-V] -e expression [-H | -B]";
	PsimagLite::String sOptions;
	bool               versionOnly = false;
	/* PSIDOC OperatorDriver
	 The arguments to the \verb!operator! executable are as follows.
	\begin{itemize}
	 \item[-f] [Mandatory, String] Input to use. The Model= line is
	very important in input.inp.

	\item[-e] [Mandatory unless -H or -B, String] OperatorExpression; see manual

	\item[-s] [Optional, Integer] \emph{Deprecated. Use -e.}
	Site for operator.
	Meaningful only for Models where
	the Hilbert space depends on the site (different kinds of atoms).
	Defaults to 0.

	\item[-B] [Optional] Prints the basis and all operators for the model

	\item[-H] [Optional] Prints the Hamiltonian terms for the model
	\end{itemize}
	 */
	while ((opt = getopt(argc, argv, "f:s:p:e:s:HB")) != -1) {
		switch (opt) {
		case 'f':
			filename = optarg;
			break;
		case 's':
			operator_options.site = atoi(optarg);
			break;
		case 'p':
			cmdline_options.precision = atoi(optarg);
			std::cout.precision(cmdline_options.precision);
			std::cerr.precision(cmdline_options.precision);
			break;
		case 'e':
			operator_options.introspect = OperatorOptions::IntrospectEnum::EXPRESSION;
			operator_options.opexpr     = optarg;
			break;
		case 'o':
			sOptions += optarg;
			break;
		case 'S':
			cmdline_options.number_of_threads = atoi(optarg);
			break;
		case 'B':
			operator_options.introspect = OperatorOptions::IntrospectEnum::MODEL_BASIS;
			break;
		case 'H':
			operator_options.introspect
			    = OperatorOptions::IntrospectEnum::MODEL_HAMILTONIAN;
			break;
		case 'U':
			cmdline_options.unbuffered_output = true;
			break;
		case 'V':
			versionOnly             = true;
			cmdline_options.logfile = "-";
			break;
		default:
			InputCheck::usageMain(strUsage);
			return 1;
		}
	}

	// sanity checks here
	if (filename == "" && !versionOnly) {
		InputCheck::usageMain(strUsage);
		return 1;
	}

	if (!cmdline_options.logfile.empty()) {
		throw std::runtime_error("operator does not support -l logfile. "
		                         "Did you mean -e?\n");
	}

	if (optind < argc) {
		std::cerr << "WARNING: Garbage at end of command line will be ignored\n";
	}

	// print license
	if (versionOnly) {
		printLicense(application, operator_options);
		return 0;
	}

	printLicense(application, operator_options);

	DmrgRunner<RealType> dmrg_runner(application);

	PsimagLite::String data;
	PsimagLite::InputNg<InputCheck>::Writeable::readFile(data, filename);

	dmrg_runner.doOneRun(data, cmdline_options, operator_options);
}
