#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <regex>
#include <string>

int main(int argc, char* argv[])
{
	if (argc != 5) {
		std::cerr << "Usage: " << argv[0] << " pattern expected tolerance file\n";
		return 1;
	}

	std::string pattern { argv[1] };
	double      expected  = std::stod(argv[2]);
	double      tolerance = std::stod(argv[3]);
	std::string filename { argv[4] };

	std::regex    re { pattern };
	std::smatch   match;
	std::ifstream file { filename };
	if (!file) {
		std::cerr << "Error: cannot open file '" << filename << "'\n";
		return 1;
	}

	std::string line;
	double      last_found_float = 0.0;
	bool        matched          = false;

	while (std::getline(file, line)) {
		std::string temp_line = line;
		while (std::regex_search(temp_line, match, re)) {
			matched          = true;
			last_found_float = std::stod(match[1].str());
			temp_line        = temp_line.substr(match.position() + match.length());
		}
	}

	if (!matched) {
		std::cerr << "Regex not found in file\n";
		return 1;
	}

	double diff = std::abs(last_found_float - expected);
	if (diff <= tolerance) {
		std::cout << "Value " << last_found_float << " within tolerance of " << tolerance
		          << " of expected value " << expected << " diff: " << diff << '\n';
		return 0;
	} else {
		std::cout << "FAIL! Value " << last_found_float << " differs from expected value "
		          << expected << " by " << diff << '\n';
		return 1;
	}
}
