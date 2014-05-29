//=======================================================================
// Author: Donovan Parks
//
// Copyright 2010 Donovan Parks
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
	std::string sequenceFile, outputDir;
	int kmerSize;
};

void help()
{
	std::cout << "Naive Bayes Train v1.0.5" << std::endl;
	std::cout << std::endl;
	std::cout << "Usage: [options] -s <sequence-file> -m <model-dir>" << std::endl;
	std::cout << std::endl;
	std::cout << "Required parameters:" << std::endl;
	std::cout << "  <sequence-file>  File listing path to each FASTA file for which a model shoud be built." << std::endl;
	std::cout << "  <model-dir>      Directory to store models." << std::endl;
	std::cout << std::endl;
	std::cout << "Optional parameters:" << std::endl;
	std::cout << "  --help        Print help message." << std::endl;
	std::cout << "  --version     Print version information." << std::endl;
	std::cout << "  --contact     Print contact information." << std::endl;
	std::cout << "  -n <integer>  Desired oligonucleotide length (default = 8)." << std::endl;
	std::cout << std::endl;
	std::cout << "Typical usage:" << std::endl;
	std::cout << "  nb-train -s sequences.txt -m ./models/"  << std::endl << std::endl;
}

bool parseCommandLine(int argc, char* argv[], Parameters& parameters)
{
	// set default values
	parameters.bShowHelp = false;
	parameters.bShowContactInfo = false;
	parameters.bShowVersion = false;
	parameters.kmerSize = 10;

	// parse parameters
	int p = 1;
	while(p < argc)
	{
		if(strcmp(argv[p], "-n") == 0)
		{
			parameters.kmerSize = atoi(argv[p+1]);
			p += 2;
		}
		else if(strcmp(argv[p], "-s") == 0)
		{
			parameters.sequenceFile = argv[p+1];
			p += 2;
		}
		else if(strcmp(argv[p], "-m") == 0)
		{
			parameters.outputDir = argv[p+1];
			
			char lastChar = parameters.outputDir[strlen(parameters.outputDir.c_str())-1];
			if(lastChar != '/' && lastChar != '\\')
				parameters.outputDir += '/';

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
		std::cout << "Naive Bayes Train v1.0.4 by Donovan Parks, Norm MacDonald, and Rob Beiko." << std::endl;
		return 0;
	}
	else if(parameters.bShowContactInfo)
	{
		std::cout << "Comments, suggestions, and bug reports can be sent to Donovan Parks (donovan.parks@gmail.com)." << std::endl;
		return 0;
	}
	else if(parameters.outputDir.empty() || parameters.sequenceFile.empty())
	{
		std::cout << "Must specify sequence file (-s) and output directory (-m)." << std::endl << std::endl;
		help();
		return 0;
	}

	// train model for each sequence in the sequence file
	std::cout << "Training models..." << std::endl;
	uint numModels = 0;
	std::ifstream modelStream(parameters.sequenceFile.c_str(), std::ios::in);
	while(!modelStream.eof())
	{
		std::string line;
		getline(modelStream, line);

		line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());

		if(line.empty())
			continue;

		std::string modelName = line.substr(line.find_last_of('/')+1, std::string::npos);
		modelName = modelName.substr(0, modelName.find_last_of('.'));

		std::cout << "  Processing model " << modelName << std::endl;

		KmerModel kmerModel(parameters.kmerSize);
		kmerModel.name(modelName);

		numModels++;		
	
		FastaIO fastaIO;
		
		bool bOK = fastaIO.open(line);
		if(!bOK)
		{
			std::cerr << "Error opening file: " << line << std::endl;
			return -1;
		}

		while(true)
		{
			SeqInfo seqInfo;
			bool bNextSeq = fastaIO.nextSeq(seqInfo);
			if(!bNextSeq)
				break;

			bOK = kmerModel.constructModel(seqInfo);
			if(!bOK)
			{
				std::cerr << "Error building model." << std::endl;
				return -1;
			}
		}	

		kmerModel.calculateConditionalProbabilities();
		kmerModel.write(parameters.outputDir + modelName + ".txt");
	}

	std::cout << std::endl;
	std::cout << "Number of models: " << numModels << std::endl;

	return 0;
}

