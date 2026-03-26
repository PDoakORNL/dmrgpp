#ifndef COOKINPUTEXPRESSION_HH
#define COOKINPUTEXPRESSION_HH
#include "InputCheck.h"
#include "InputNg.h"
#include "Matrix.h"
#include "PsimagLite.h"
#include <string>
#include <vector>

namespace Dmrg {

template <typename ComplexOrRealType> class CookInputExpression {
public:

	using RealType             = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using InputNgType          = PsimagLite::InputNg<Dmrg::InputCheck>;
	using InputNgValidatorType = InputNgType::Readable;

	CookInputExpression(InputNgValidatorType& io)
	    : io_(io)
	{ }

	std::string operator()(const std::string& expr)
	{
		std::string label = "!readTable";
		if (expr.substr(0, label.size()) == label) {
			std::string str = expr.substr(label.size(), std::string::npos);

			// delete spaces if any
			str.erase(std::remove_if(str.begin(), str.end(), isspace), str.end());

			// delete parens if any
			SizeType length = str.length();
			if (length > 0 && str[0] == '(' && str[length - 1] == ')') {
				str.erase(0, 1);
				if (length > 1) {
					str.erase(length - 2, 1);
				}
			}

			// split on comma
			std::vector<std::string> args;
			PsimagLite::split(args, str, ",");
			if (args.size() != 2) {
				err("readTable expects two arguments\n");
			}

			PsimagLite::Matrix<RealType> matrix;
			io_.read(matrix, args[0]);
			RealType value = findValueFor(matrix, PsimagLite::atof(args[1]));
			return ttos(value);
		} else {
			return expr;
		}
	}

private:

	static RealType findValueFor(const PsimagLite::Matrix<RealType>& matrix, const RealType& t)
	{
		SizeType rows = matrix.rows();
		if (matrix.cols() != 2) {
			err("findValueFor(): not a table\n");
		}

		for (SizeType i = 0; i < rows; ++i) {
			if (matrix(i, 0) == t) {
				return matrix(i, 1);
			}
		}

		throw std::runtime_error("Value not found in table\n");
	}

	InputNgValidatorType& io_;
};
}
#endif // COOKINPUTEXPRESSION_HH
