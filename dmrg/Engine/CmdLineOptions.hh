#ifndef CMDLINEOPTIONS_HH
#define CMDLINEOPTIONS_HH
#include <string>

namespace Dmrg {

struct CmdLineOptions {
	std::string  solver_options;
	std::string  in_situ_measurements;
	std::string  logfile;
	unsigned int precision         = 12;
	unsigned int number_of_threads = 0;
	bool         unbuffered_output = false;
};
}
#endif // CMDLINEOPTIONS_HH
