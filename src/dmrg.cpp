#include "BlockDiagonalMatrix.h"
#include "CmdLineOptions.hh"
#include "Concurrency.h"
#include "DmrgDriver.h"
#include "DmrgRunner.h"
#include "Io/IoNg.h"
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
	PsimagLite::PsiApp application("DMRG++", &argc, &argv, 1);
	PsimagLite::String filename = "";
	int                opt      = 0;
	OperatorOptions    options;
	CmdLineOptions     cmdline_options;
	PsimagLite::String strUsage(application.name());
	if (PsimagLite::basename(argv[0]) == "operator")
		options.enabled = true;
	strUsage += " -f filename [-k] [-p precision] [-o solverOptions] [-V] [whatToMeasure]";
	PsimagLite::String sOptions;
	bool               versionOnly = false;
	/* PSIDOC DmrgDriver
There is a single input file that is passed as the
argument to \verb!-f!, like so
\begin{lstlisting}
	./dmrg -f input.inp whatToMeasure
\end{lstlisting}
where \verb!whatToMeasure! is optional. The command line arguments
to the main dmrg driver are the following.
	  \begin{itemize}
	  \item[-f] {[}Mandatory, String{]} Input to use.
	  \item[-o] {[}Optional, String{]} Extra options for SolverOptions
	  \item[-p] [Optional, Integer] Digits of precision for printing.
	  \item[whatToMeasure] {[}Optional, String{]} What to measure in-situ.
	  This is a comma-separated list of braket specifications.
	  Braket specifications can be bare or dressed, and are explained elsewhere.
	  \item[-l] {[}Optional, String{]} Without this option std::cout is redirected
	  to a file.
	  This option with the string ``?'' prints name of such log file.
	  This option with the string ``-'' writes std::cout to terminal.
	  In other cases, string is the name of the file to redirect std::cout to.
	 \item[-k] [Optional] Keep untar files
	 \item[-U] [Optional] Make cout output unbuffered
	 \item[-S] [Optional, number] Ignore the Threads= line if present in the input,
	  and run with Threads=number
	 \item[-V] [Optional] Print version and exit
	  \end{itemize}
	 */
	/* PSIDOC OperatorDriver
	 The arguments to the \verb!operator! executable are as follows.
	\begin{itemize}
	 \item[-f] [Mandatory, String] Input to use. The Model= line is
	very important in input.inp.

	\item[-e] [Mandatory unless -l, String] OperatorExpression; see manual

	\item[-s] [Optional, Integer] \emph{Deprecated. Use -e.}
	Site for operator.
	Meaningful only for Models where
	the Hilbert space depends on the site (different kinds of atoms).
	Defaults to 0.

	\item[-B] [Optional] Prints the basis and all operators for the model

	\item[-H] [Optional] Prints the Hamiltonian terms for the model
	\end{itemize}
	 */
	while ((opt = getopt(argc, argv, "f:s:l:p:e:o:S:kBHUV")) != -1) {
		switch (opt) {
		case 'f':
			filename = optarg;
			break;
		case 's':
			options.site = atoi(optarg);
			break;
		case 'l':
			cmdline_options.logfile = optarg;
			break;
		case 'p':
			cmdline_options.precision = atoi(optarg);
			std::cout.precision(cmdline_options.precision);
			std::cerr.precision(cmdline_options.precision);
			break;
		case 'e':
			options.introspect = OperatorOptions::IntrospectEnum::EXPRESSION;
			options.opexpr     = optarg;
			break;
		case 'o':
			sOptions += optarg;
			break;
		case 'S':
			cmdline_options.number_of_threads = atoi(optarg);
			break;
		case 'B':
			options.introspect = OperatorOptions::IntrospectEnum::MODEL_BASIS;
			break;
		case 'H':
			options.introspect = OperatorOptions::IntrospectEnum::MODEL_HAMILTONIAN;
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

	if (options.enabled && !cmdline_options.logfile.empty()) {
		throw std::runtime_error("operator does not support -l logfile. "
		                         "Did you mean -e?\n");
	}

	cmdline_options.in_situ_measurements = (optind < argc) ? argv[optind] : "";

	if (cmdline_options.logfile == "?") { // query only
		std::cout << ProgramGlobals::coutName(filename, application.name()) << "\n";
		return 0;
	} else if (cmdline_options.logfile.empty()) {
		cmdline_options.logfile = ProgramGlobals::coutName(filename, application.name());
	}

	// print license
	if (versionOnly) {
		printLicense(application, options);
		return 0;
	}

	printLicense(application, options);

	DmrgRunner<RealType> dmrg_runner(application);

	PsimagLite::String data;
	PsimagLite::InputNg<InputCheck>::Writeable::readFile(data, filename);

	dmrg_runner.doOneRun(data, cmdline_options);
}
