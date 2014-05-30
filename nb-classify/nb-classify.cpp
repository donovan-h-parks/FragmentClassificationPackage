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
#include "KmerCalculator.hpp"
#include "KmerModel.hpp"
#include "Utils.hpp"

struct Parameters
{
	bool bShowHelp, bShowVersion, bShowContactInfo;
	std::string queryFile, modelFile, resultsFile, tempExtension;
	int batchSize, topModels, verbose;
};

void help()
{
	std::cout << "Naive Bayes Classification v1.0.5" << std::endl;
	std::cout << std::endl;
	std::cout << "  Usage: [options] -q <query-file> -m <model-file> -r <results-file>" << std::endl;
	std::cout << std::endl;
	std::cout << "Required parameters:" << std::endl;
	std::cout << "  <query-file>    Multi-FASTA file containing query fragments to classify." << std::endl;
	std::cout << "  <model-file>    File indicating models to use for classification." << std::endl;
	std::cout << "  <results-file>  File to write classification results to." << std::endl;
	std::cout << std::endl;
	std::cout << "Optional parameters:" << std::endl;
	std::cout << "  --help        Print help message." << std::endl;
	std::cout << "  --version     Print version information." << std::endl;
	std::cout << "  --contact     Print contact information." << std::endl;
	std::cout << "  -b <integer>  Number of fragments to classify at a time (default = 50000)." << std::endl;	
	std::cout << "  -t <integer>  Log likelihood of the top T models will be returned. If you" << std::endl;  
	std::cout << "                  wish to have the log likelihood of all models in the" << std::endl;
	std::cout << "                  results file set T = 0 (default = 0)." << std::endl;
	std::cout << "  -v <integer>  Level of output information (default = 1)." << std::endl;	
	std::cout << "  -e <string>   Extension to add to temporary files (default = txt)." << std::endl;	
	std::cout << std::endl;
	std::cout << "Typical usage:" << std::endl;
	std::cout << "  nb-classify -q test.fasta -m models.txt -r nb_results.txt" << std::endl << std::endl;
}

bool parseCommandLine(int argc, char* argv[], Parameters& parameters)
{
	// set default values
	parameters.bShowHelp = false;
	parameters.bShowContactInfo = false;
	parameters.bShowVersion = false;
	parameters.batchSize = 50000;
	parameters.topModels = 0;
	parameters.verbose = 1;
	parameters.tempExtension = "txt";

	// parse parameters
	int p = 1;
	while(p < argc)
	{
		if(strcmp(argv[p], "-q") == 0)
		{
			parameters.queryFile = argv[p+1];
			p += 2;
		}
		else if(strcmp(argv[p], "-m") == 0)
		{
			parameters.modelFile = argv[p+1];
			p += 2;
		}
		else if(strcmp(argv[p], "-r") == 0)
		{
			parameters.resultsFile = argv[p+1];
			p += 2;
		}
		else if(strcmp(argv[p], "-b") == 0)
		{
			parameters.batchSize = atoi(argv[p+1]);
			p += 2;
		}
		else if(strcmp(argv[p], "-t") == 0)
		{
			parameters.topModels = atoi(argv[p+1]);
			p += 2;
		}
		else if(strcmp(argv[p], "-v") == 0)
		{
			parameters.verbose = atoi(argv[p+1]);
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
		else if(strcmp(argv[p], "-e") == 0)
		{
			parameters.tempExtension = argv[p+1];
			p += 2;
		}
		else
		{
			std::cout << "Unrecognized parameter: " << argv[p] << std::endl << std::endl;
			return false;
		}
	}

	return true;
}

struct TopModel
{
	TopModel(): modelNum(0), logLikelihood(-1e10) {} 
	TopModel(uint _modelNum, float _logLikelihood): modelNum(_modelNum), logLikelihood(_logLikelihood) {}

	uint modelNum;
	float logLikelihood;
};

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
		std::cout << "Naive Bayes Classify v1.0.5 by Donovan Parks, Norm MacDonald, and Rob Beiko." << std::endl;
		return 0;
	}
	else if(parameters.bShowContactInfo)
	{
		std::cout << "Comments, suggestions, and bug reports can be sent to Donovan Parks (donovan.parks@gmail.com)." << std::endl;
		return 0;
	}
	else if(parameters.queryFile.empty() || parameters.modelFile.empty() || parameters.resultsFile.empty())
	{
		std::cout << "Must specify query (-q), model (-m), and result (-r) file." << std::endl << std::endl;
		help();
		return 0;
	}

	bool bRecordAllModels = false;
	if(parameters.topModels <= 0)
	{
		bRecordAllModels = true;
		parameters.topModels = 0;
	}

	// Get model k-mer length
	if(parameters.verbose >= 1)
		std::cout << "Determining n-mer length..." << std::endl;

	std::ifstream tempStream(parameters.modelFile.c_str(), std::ios::in);	
	if(tempStream.fail())
	{
		std::cout << "Failed to open model file: " << parameters.modelFile << std::endl << std::endl;
		return -1;
	}
	std::string line;
	std::getline(tempStream, line);	
	KmerModel tempModel(line);
	uint kmerLength = tempModel.kmerLength();
	if(parameters.verbose >= 1)
		std::cout << "  n-mer length: " << kmerLength << std::endl << std::endl;
	
	// Read query fragments

	if(parameters.verbose >= 1)
		std::cout << "Reading query fragments..." << std::endl;

	char* buffer = NULL;
	std::vector<SeqInfo> querySeqs;
	FastaIO fastaIO;
	bool bSuccess = fastaIO.readSeqs(parameters.queryFile, querySeqs, buffer, parameters.verbose);
	if(!bSuccess)
	{
		std::cout << "Failed to open query fragment file: " << parameters.queryFile << std::endl;
		return -1;
	}
	if(parameters.verbose >= 1)
		std::cout << "  Number of query fragments: " << querySeqs.size() << std::endl << std::endl;

	// Classify query fragments in batches in order to keep memory requirements within reason (~ 1GB)
	if(parameters.verbose >= 1)
		std::cout << "Processing query fragments in batches of " << parameters.batchSize << "." << std::endl << std::endl;

	KmerCalculator kmerCalculator(kmerLength);
	for(uint batchNum = 0; batchNum < ceil(double(querySeqs.size()) / parameters.batchSize); ++batchNum)
	{
		if(parameters.verbose >= 1)
			std::cout << "Batch #" << (batchNum+1) << std::endl;

		// get k-mers for each query fragment
		if(parameters.verbose >= 1)
			std::cout << "  Calculating n-mers in query fragment: " << std::endl;	

		std::vector< std::vector<uint> > queryKmerProfiles;
		queryKmerProfiles.reserve(parameters.batchSize);
		for(uint seqIndex = batchNum*parameters.batchSize; 
					seqIndex < std::min(ulong(querySeqs.size()), ulong(batchNum+1)*parameters.batchSize); 
					++seqIndex)
		{
			if(parameters.verbose >= 3)
				std::cout << querySeqs.at(seqIndex).seqId << std::endl;
			else if (seqIndex % 5000 == 0 && parameters.verbose >= 1)
				std::cout << "." << std::flush;

			std::vector<uint> profile;
			kmerCalculator.extractForwardKmers(querySeqs.at(seqIndex), profile);
			queryKmerProfiles.push_back(profile);
		}
		if(parameters.verbose >= 1)
			std::cout << std::endl;

		// apply each model to each query sequence
		if(parameters.verbose >= 1)
			std::cout << "  Applying models to query sequences: " << std::endl;

		std::ifstream modelStream(parameters.modelFile.c_str(), std::ios::in);

		uint modelNum = 0;

		std::vector<std::string> modelNames;
		std::vector< std::list<TopModel> > topModelsPerFragment(queryKmerProfiles.size());		
		std::vector< std::vector<float> > modelLogLikelihoods;
		while(!modelStream.eof())
		{
			std::string line;
			std::getline(modelStream, line);

			if(line.empty())
				break;

			if(modelNum % 200 == 0 && parameters.verbose >= 1)
				std::cout << " " << modelNum << std::flush;
			
			KmerModel kmerModel(line);
			modelNames.push_back(kmerModel.name());
			if(parameters.verbose >= 2)
			{
				kmerModel.printModelInfo(std::cout);
				std::cout << std::endl;
			}

			ulong size = 0;
			if(bRecordAllModels)
				size = queryKmerProfiles.size();		
			std::vector<float> logLikelihoods(size);
			for(uint seqIndex = 0; seqIndex < queryKmerProfiles.size(); ++seqIndex)
			{
				SeqInfo querySeqInfo = querySeqs[seqIndex + batchNum*parameters.batchSize];	
				float logLikelihood = kmerModel.classify(querySeqInfo, queryKmerProfiles[seqIndex]);

				// record models with highest log likelihood
				if(bRecordAllModels)
				{
					logLikelihoods[seqIndex] = logLikelihood;
				}
				else
				{
					std::list<TopModel> topModels = topModelsPerFragment.at(seqIndex);

					if(topModels.size() == 0)
						topModels.push_front(TopModel(modelNum, logLikelihood));

					std::list<TopModel>::iterator it;
					bool bInserted = false;
					for(it = topModels.begin(); it != topModels.end(); it++)
					{
						if(logLikelihood > it->logLikelihood)
						{
							topModels.insert(it, TopModel(modelNum, logLikelihood));
							bInserted = true;
							break;
						}
					}

					if((int)topModels.size() < parameters.topModels && !bInserted)
						topModels.push_back(TopModel(modelNum, logLikelihood));
					else if((int)topModels.size() > parameters.topModels)
							topModels.pop_back();

					topModelsPerFragment.at(seqIndex) = topModels;
				}					
			}

			if(bRecordAllModels)
				modelLogLikelihoods.push_back(logLikelihoods);	
		
			modelNum++;
		}
		if(parameters.verbose >= 1)
			std::cout << std::endl;

		// write out classification
		if(parameters.verbose >= 1)
			std::cout << "  Writing out classification results." << std::endl << std::endl;

		std::stringstream outputTempResults;
		outputTempResults << "./batch_" << batchNum << "." << parameters.tempExtension;		
		std::ofstream fout(outputTempResults.str().c_str(), std::ios::out);	
		if(fout.fail())
		{
			std::cout << "Failed to write temporary results file: " << outputTempResults.str() << std::endl;
			return -1;
		}

		// check if all model results are to be written out
		if(bRecordAllModels)
		{		
			if(batchNum == 0)
			{
				fout << "Fragment Id" << "\t" << "Length" << "\t" << "Valid n-mers";
				for(uint modelIndex = 0; modelIndex < modelNames.size(); ++modelIndex)
					fout << "\t" << modelNames[modelIndex];
				fout << std::endl;
			}

			for(uint seqIndex = 0; seqIndex < queryKmerProfiles.size(); ++seqIndex)
			{
				SeqInfo querySeqInfo = querySeqs.at(seqIndex + batchNum*parameters.batchSize);

				fout << querySeqInfo.seqId << "\t" << querySeqInfo.length << "\t" << querySeqInfo.validKmers;

				for(uint modelIndex = 0; modelIndex < modelNames.size(); ++modelIndex)
					fout << "\t" << modelLogLikelihoods[modelIndex][seqIndex];
				fout << std::endl;
			}
		}
		else
		{
			for(uint seqIndex = 0; seqIndex < queryKmerProfiles.size(); ++seqIndex)
			{
				SeqInfo querySeqInfo = querySeqs.at(seqIndex + batchNum*parameters.batchSize);

				fout << querySeqInfo.seqId << "\t" << querySeqInfo.length << "\t" << querySeqInfo.validKmers;

				std::list<TopModel>::iterator it;
				for(it = topModelsPerFragment.at(seqIndex).begin(); it != topModelsPerFragment.at(seqIndex).end(); it++)
					fout << "\t" << modelNames[it->modelNum] << "\t" << it->logLikelihood;
			
				fout << std::endl;
			}
		}

		fout.close();
	}
	
	// free memory allocated to hold query fragment data
	delete[] buffer;

	// Concatenate result files
	if(parameters.verbose >= 1)
		std::cout << "Building results file: ";

	std::ofstream resultsStream(parameters.resultsFile.c_str(), std::ios::out | std::ios::binary);
	for(uint batchNum = 0; batchNum < ceil(double(querySeqs.size()) / parameters.batchSize); ++batchNum)
	{
		if(parameters.verbose >= 1)
			std::cout << "." << std::flush;

		std::stringstream tempResultFile;
		tempResultFile << "./batch_" << batchNum  << "." << parameters.tempExtension;		
		std::ifstream tempStream(tempResultFile.str().c_str(), std::ios::binary);
		if(tempStream.fail() || tempStream.bad())
		{
			std::cout << "Failed to open file: " << tempResultFile.str() << std::endl;
			return -1;
		}

		// calculate size of file
		tempStream.seekg(0, std::ios::end);
		ulong fileSize = tempStream.tellg();
		tempStream.seekg(0, std::ios::beg);
		
		// write out data in reasonable sized chunks
		ulong chunkSize = 64*1024*1024;
		
		// allocate memory for reading file
		char* tempBuffer = new char[chunkSize];
		if(tempBuffer == NULL)
		{
			std::cout << std::endl << "Failed to allocate memory required by file: " << tempResultFile.str() << std::endl;
			return -1;
		}
		
		for(uint chunk = 0; chunk < ceil(float(fileSize) / chunkSize); ++chunk)
		{
			ulong currentChunkSize = std::min(chunkSize, fileSize - chunk*chunkSize);

			// read file into buffer
			tempStream.read(tempBuffer, currentChunkSize);
			if(tempStream.fail() || tempStream.bad())
			{
				std::cout << std::endl << "Failed to read data from " << tempResultFile.str() << std::endl;
				return -1;
			}

			resultsStream.write(tempBuffer, currentChunkSize);
			resultsStream.flush();
		}

		tempStream.close();
		delete[] tempBuffer;
	}
	resultsStream.close();

	if(parameters.verbose >= 1)
	{
		std::cout << std::endl;
		std::cout << "Done." << std::endl;
	}

	for(uint batchNum = 0; batchNum < ceil(double(querySeqs.size()) / parameters.batchSize); ++batchNum)
	{
		std::stringstream filename;
		filename << "./batch_" << batchNum  << "." << parameters.tempExtension;
		std::remove(filename.str().c_str());
	}
	
	return 0;
}

