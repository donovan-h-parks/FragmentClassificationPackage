//=======================================================================
// Author: Donovan Parks
//
// Copyright 2013 Donovan Parks
//
// This file is part of NaiveBayes.
//
// NaiveBayes is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// NaiveBayes is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with NaiveBayes.  If not, see <http://www.gnu.org/licenses/>.
//=======================================================================

#include "stdafx.h"

#include "FastaIO.hpp"
#include "KmerModel.hpp"
#include "Utils.hpp"

struct Parameters
{
	bool bShowHelp, bShowVersion, bShowContactInfo;
	std::string inputFile, outputFile;
	int kmerSize;
};

void help()
{
	std::cout << "Usage: [options] -i <input-file> -o <output-file>" << std::endl;
	std::cout << std::endl;
	std::cout << "Required parameters:" << std::endl;
	std::cout << "  <input-file>  FASTA file containing sequences for null model." << std::endl;
	std::cout << "  <output-file> Output file for null model." << std::endl;
	std::cout << std::endl;
	std::cout << "Optional parameters:" << std::endl;
	std::cout << "  --help        Print help message." << std::endl;
	std::cout << "  --version     Print version information." << std::endl;
	std::cout << "  --contact     Print contact information." << std::endl;
	std::cout << "  -n <integer>  Desired oligonucleotide length (default = 8)." << std::endl;
	std::cout << std::endl;
	std::cout << "Typical usage:" << std::endl;
	std::cout << "  nb-null-model -s sequences.fna -m null_model.txt"  << std::endl << std::endl;
}

bool parseCommandLine(int argc, char* argv[], Parameters& parameters)
{
	// set default values
	parameters.bShowHelp = false;
	parameters.bShowContactInfo = false;
	parameters.bShowVersion = false;
	parameters.kmerSize = 8;

	// parse parameters
	int p = 1;
	while(p < argc)
	{
		if(strcmp(argv[p], "-n") == 0)
		{
			parameters.kmerSize = atoi(argv[p+1]);
			p += 2;
		}
		else if(strcmp(argv[p], "-i") == 0)
		{
			parameters.inputFile = argv[p+1];
			p += 2;
		}
		else if(strcmp(argv[p], "-o") == 0)
		{
			parameters.outputFile = argv[p+1];
			p += 2;
		}
		else if(strcmp(argv[p], "--help") == 0)
		{
			parameters.bShowHelp = true;
			p += 1;
		}
		else if(strcmp(argv[p], "--version") == 0)
		{
			parameters.bShowVersion = true;
			p += 1;
		}
		else if(strcmp(argv[p], "--contact") == 0)
		{
			parameters.bShowContactInfo = true;
			p += 1;
		}
		else
		{
			std::cout << "Unrecognized parameter: " << argv[p] << std::endl << std::endl;
			return false;
		}
	}

	return true;
}

int main(int argc, char* argv[])
{
	// Parse command-line arguments
	Parameters parameters;
	bool bParsed = parseCommandLine(argc, argv, parameters);

	if(!bParsed || parameters.bShowHelp || argc == 1) 
	{			
		help();	
    return 0;
  }
	else if(parameters.bShowVersion)
	{
		std::cout << "nb-null-model v1.1 by Donovan Parks, Norm MacDonald, and Rob Beiko." << std::endl;
		return 0;
	}
	else if(parameters.bShowContactInfo)
	{
		std::cout << "Comments, suggestions, and bug reports can be sent to Donovan Parks (donovan.parks@gmail.com)." << std::endl;
		return 0;
	}
	else if(parameters.inputFile.empty() || parameters.outputFile.empty())
	{
		std::cout << "Must specify both an input (-i) and output file (-o)." << std::endl << std::endl;
		help();
		return 0;
	}

	// train null model
	std::cout << "Building null model..." << std::endl;

	KmerModel kmerModel(parameters.kmerSize);
	kmerModel.name("Null Model");	

	FastaIO fastaIO;
	bool bOK = fastaIO.open(parameters.inputFile);
	if(!bOK)
	{
		std::cerr << "Error opening file: " << parameters.inputFile << std::endl;
		return -1;
	}

	while(true)
	{
		SeqInfo seqInfo;
		bool bNextSeq = fastaIO.nextSeq(seqInfo);
		if(!bNextSeq)
			break;

		bool bOK = kmerModel.constructModel(seqInfo);
		if(!bOK)
		{
			std::cerr << "[Error] Failed to build model." << std::endl;
			return -1;
		}
	}	

	kmerModel.calculateConditionalProbabilities();
	kmerModel.write(parameters.outputFile);

	return 0;
}

