#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <regex>
#include <string>

int main(int argc, char* argv[])
{
	// Optional flag: --first  => stop at the first regex match (default: use last match)
	bool use_first = false;
	int  arg_start = 1;
	if (argc >= 2 && std::string(argv[1]) == "--first") {
		use_first = true;
		arg_start = 2;
	}

	if (argc - arg_start != 4) {
		std::cerr << "Usage: " << argv[0]
		          << " [--first] pattern expected tolerance file\n";
		return 1;
	}

	std::string pattern { argv[arg_start] };
	double      expected  = std::stod(argv[arg_start + 1]);
	double      tolerance = std::stod(argv[arg_start + 2]);
	std::string filename { argv[arg_start + 3] };

	std::regex    re { pattern };
	std::smatch   match;
	std::ifstream file { filename };
	if (!file) {
		std::cerr << "Error: cannot open file '" << filename << "'\n";
		return 1;
	}

	std::string line;
	double      found_float = 0.0;
	bool        matched     = false;

	while (std::getline(file, line)) {
		std::string temp_line = line;
		while (std::regex_search(temp_line, match, re)) {
			found_float = std::stod(match[1].str());
			matched     = true;
			if (use_first)
				break;
			temp_line = temp_line.substr(match.position() + match.length());
		}
		if (matched && use_first)
			break;
	}

	if (!matched) {
		std::cerr << "Regex not found in file\n";
		return 1;
	}

	double diff = std::abs(found_float - expected);
	if (diff <= tolerance) {
		std::cout << "Value " << found_float << " within tolerance of " << tolerance
		          << " of expected value " << expected << " diff: " << diff << '\n';
		return 0;
	} else {
		std::cout << "FAIL! Value " << found_float << " differs from expected value "
		          << expected << " by " << diff << '\n';
		return 1;
	}
}
