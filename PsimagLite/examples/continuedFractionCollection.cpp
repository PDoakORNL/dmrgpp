// BEGIN LICENSE BLOCK
/*
Copyright (c) 2009 , UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 1.0.0]

*********************************************************
THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED.

Please see full open source license included in file LICENSE.
*********************************************************

*/
// END LICENSE BLOCK

#include "ContinuedFractionCollection.h"
#include "ContinuedFraction.h"
#include "Io/IoSimple.h"
#include "Matsubaras.h"
#include "RealFrequencyRange.hh"
#include "TridiagonalMatrix.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

using namespace PsimagLite;
typedef double                                RealType;
typedef ContinuedFraction<RealType>           ContinuedFractionType;
typedef ContinuedFractionCollection<RealType> ContinuedFractionCollectionType;

void usage(const char* progName)
{
	std::cerr << "Usage: " << progName << " -f file  -b omega1";
	std::cerr << " -t total -s omegaStep -d delta -B beta -m matsubaras\n";
}

template <typename PlotParamsType>
void plotAll(const ContinuedFractionCollectionType& cfCollection, const PlotParamsType& params)
{
	ContinuedFractionCollectionType::PlotDataType v;
	cfCollection.plot(v, params);
	std::cout.precision(12);
	for (SizeType x = 0; x < v.size(); x++) {
		std::cout << v[x].first << " " << PsimagLite::imag(v[x].second);
		std::cout << " " << PsimagLite::real(v[x].second) << "\n";
	}
}

template <typename PlotParamsType>
void plotOneByOne(const ContinuedFractionCollectionType& cfCollection, const PlotParamsType& params)
{

	std::cout.precision(12);
	for (SizeType i = 0; i < cfCollection.size(); i++) {
		ContinuedFractionCollectionType::PlotDataType v;
		cfCollection.plotOne(i, v, params);
		for (SizeType x = 0; x < v.size(); x++) {
			std::cout << v[x].first << " " << PsimagLite::imag(v[x].second);
			std::cout << " " << PsimagLite::real(v[x].second) << "\n";
		}
	}
}

int main(int argc, char* argv[])
{
	int      opt        = 0;
	String   file       = "";
	RealType wbegin     = 0;
	SizeType total      = 0;
	RealType wstep      = 0;
	RealType delta      = 0;
	RealType beta       = 0.0;
	SizeType matsubaras = 0;
	bool     oneByOne   = false;
	while ((opt = getopt(argc, argv, "f:b:t:s:d:m:B:1")) != -1) {
		switch (opt) {
		case 'f':
			file = optarg;
			break;
		case 'b':
			wbegin = atof(optarg);
			break;
		case 't':
			total = atoi(optarg);
			break;
		case 's':
			wstep = atof(optarg);
			break;
		case 'd':
			delta = atof(optarg);
			break;
		case '1':
			oneByOne = true;
			break;
		case 'm':
			matsubaras = atoi(optarg);
			break;
		case 'B':
			beta = atof(optarg);
			break;
		default:
			usage(argv[0]);
			return 1;
		}
	}
	// sanity checks:
	bool real1                = (total < 0 || wstep <= 0 || delta <= 0);
	bool imag1                = (beta <= 0 || matsubaras == 0);
	bool not_valid_matsubaras = (beta > 0 && matsubaras == 0);
	not_valid_matsubaras      = not_valid_matsubaras || (beta <= 0 && matsubaras > 0);

	if (file == "" || (real1 & imag1) || not_valid_matsubaras || beta < 0) {
		usage(argv[0]);
		return 1;
	}

	IoSimple::In                    io(file);
	ContinuedFractionCollectionType cfCollection(io);

	FrequencyRange<RealType>* params        = nullptr;
	bool                      is_matsubaras = (beta > 0 && matsubaras > 0);
	if (is_matsubaras) {
		params = new PsimagLite::Matsubaras<RealType>(beta, matsubaras, delta);
	} else {
		params = new PsimagLite::RealFrequencyRange<RealType>(wbegin, wstep, total, delta);
	}

	if (!oneByOne)
		plotAll(cfCollection, *params);
	else
		plotOneByOne(cfCollection, *params);

	delete params;
	params = nullptr;
}
