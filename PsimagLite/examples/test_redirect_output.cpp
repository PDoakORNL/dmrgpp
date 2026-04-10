#include "RedirectOutput.hh"

int main(int argc, char* argv[])
{
	PsimagLite::RedirectOutput redirect_output;

	if (argc == 2) {
		redirect_output.doIt("test", std::ofstream::out, true);
	}

	std::cout << "Hello!\n";
}
