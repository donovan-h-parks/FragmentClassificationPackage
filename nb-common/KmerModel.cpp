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

#include "KmerModel.hpp"

using namespace std;

KmerModel::KmerModel(uint wordLength): m_kmerCalculator(new KmerCalculator(wordLength)), m_wordLength(wordLength)
{
	// store log probability of each kmer number
	m_logConditionalProb = new float[m_kmerCalculator->numPossibleWords()];
	memset(m_logConditionalProb, 0, m_kmerCalculator->numPossibleWords()*sizeof(float));
}

KmerModel::KmerModel(const std::string& modelFile)
{
	read(modelFile);
}

KmerModel::~KmerModel()
{
	delete m_kmerCalculator;
	delete[] m_logConditionalProb;
}

bool KmerModel::constructModel(SeqInfo& seqInfo)
{
	if(seqInfo.length < m_wordLength)
	{
		std::cerr << "Sequence is shorter then k-mer length: " << seqInfo.seqId << "." << std::endl;
		return false;
	}

	// ensure all sequences are from the same strain
	if(m_modelInfo.taxonomy.strain == "")
	{
		// initialize model
		m_modelInfo.taxonomy = seqInfo.taxonomy;
		m_modelInfo.numWords = 0;
		m_modelInfo.numSeqs = 0;
	}
	else
	{
		if(seqInfo.taxonomy.strain != m_modelInfo.taxonomy.strain)
		{
			std::cerr << "Sequence " << seqInfo.seqId << " is from strain " << seqInfo.taxonomy.strain; 
			std::cerr << ", expecting sequences from strain " << m_modelInfo.taxonomy.strain << "." << std::endl;
			return false;
		}
	}

	// calculate kmer profile for sequence
	std::vector<uint> kmerVector;
	m_kmerCalculator->extractKmers(seqInfo, kmerVector);
	m_modelInfo.numWords += seqInfo.validKmers;
	m_modelInfo.numSeqs++;

	// record kmers
	for(ulong i = 0; i < seqInfo.validKmers; ++i)
		m_logConditionalProb[kmerVector.at(i)]++;

	return true;
}


void KmerModel::calculateConditionalProbabilities()
{
	// calculate log conditional probabilities
	for(uint i = 0; i < m_kmerCalculator->numPossibleWords(); ++i)
	{
		float count = m_logConditionalProb[i];
		m_logConditionalProb[i] = log((count + 1.0f) / (m_modelInfo.numWords + m_kmerCalculator->numPossibleWords())); 
	}
}

float KmerModel::classify(SeqInfo& seqInfo) 
{
	std::vector<uint> kmers;
	m_kmerCalculator->extractForwardKmers(seqInfo, kmers);

	return classify(seqInfo, kmers); 
}


float KmerModel::classify(const SeqInfo& seqInfo, const std::vector<uint>& profile) 
{
	float logProb = 0.0f;
	for(ulong i = 0; i < seqInfo.validKmers; ++i)
		logProb += m_logConditionalProb[profile.at(i)];

	return logProb;
}

void KmerModel::write(const std::string& filename) const
{
	std::ofstream fout(filename.c_str(), std::ios::out | std::ios::binary);
	if(!fout.is_open())
	{
		std::cout << "Failed to write model to file: " << filename << std::endl;
		return;
	}

	fout.write((char*)&m_wordLength, sizeof(uint));

	size_t size = m_modelInfo.name.size();
	fout.write((char*)&size, sizeof(size_t));
	fout.write(m_modelInfo.name.c_str(), m_modelInfo.name.size());

	size = m_modelInfo.taxonomy.superKingdom.size();
	fout.write((char*)&size, sizeof(size_t));
	fout.write(m_modelInfo.taxonomy.superKingdom.c_str(), m_modelInfo.taxonomy.superKingdom.size());

	size = m_modelInfo.taxonomy.phylum.size();
	fout.write((char*)&size, sizeof(size_t));
	fout.write(m_modelInfo.taxonomy.phylum.c_str(), m_modelInfo.taxonomy.phylum.size());

	size = m_modelInfo.taxonomy.taxonomicClass.size();
	fout.write((char*)&size, sizeof(size_t));
	fout.write(m_modelInfo.taxonomy.taxonomicClass.c_str(), m_modelInfo.taxonomy.taxonomicClass.size());

	size = m_modelInfo.taxonomy.order.size();
	fout.write((char*)&size, sizeof(size_t));
	fout.write(m_modelInfo.taxonomy.order.c_str(), m_modelInfo.taxonomy.order.size());

	size = m_modelInfo.taxonomy.family.size();
	fout.write((char*)&size, sizeof(size_t));
	fout.write(m_modelInfo.taxonomy.family.c_str(), m_modelInfo.taxonomy.family.size());

	size = m_modelInfo.taxonomy.genus.size();
	fout.write((char*)&size, sizeof(size_t));
	fout.write(m_modelInfo.taxonomy.genus.c_str(), m_modelInfo.taxonomy.genus.size());

	size = m_modelInfo.taxonomy.species.size();
	fout.write((char*)&size, sizeof(size_t));
	fout.write(m_modelInfo.taxonomy.species.c_str(), m_modelInfo.taxonomy.species.size());

	size = m_modelInfo.taxonomy.strain.size();
	fout.write((char*)&size, sizeof(size_t));
	fout.write(m_modelInfo.taxonomy.strain.c_str(), m_modelInfo.taxonomy.strain.size());

	fout.write((char*)m_logConditionalProb, sizeof(float)*m_kmerCalculator->numPossibleWords());

	fout.close();
}

void KmerModel::read(const std::string& filename)
{
	std::ifstream fin(filename.c_str(), std::ios::in | std::ios::binary);
	if(!fin.is_open())
	{
		std::cout << "Failed to read model from file: " << filename << "." << std::endl;
		return;
	}

	fin.read((char*)&m_wordLength, sizeof(uint));

	size_t size;
	char buffer[1024];

	fin.read((char*)&size, sizeof(size_t));
	fin.read(buffer, size);
	buffer[size] = 0;
	m_modelInfo.name = buffer;

	fin.read((char*)&size, sizeof(size_t));
	fin.read(buffer, size);
	buffer[size] = 0;
	m_modelInfo.taxonomy.superKingdom = buffer;

	fin.read((char*)&size, sizeof(size_t));
	fin.read(buffer, size);
	buffer[size] = 0;
	m_modelInfo.taxonomy.phylum = buffer;

	fin.read((char*)&size, sizeof(size_t));
	fin.read(buffer, size);
	buffer[size] = 0;
	m_modelInfo.taxonomy.taxonomicClass = buffer;

	fin.read((char*)&size, sizeof(size_t));
	fin.read(buffer, size);
	buffer[size] = 0;
	m_modelInfo.taxonomy.order = buffer;

	fin.read((char*)&size, sizeof(size_t));
	fin.read(buffer, size);
	buffer[size] = 0;
	m_modelInfo.taxonomy.family = buffer;

	fin.read((char*)&size, sizeof(size_t));
	fin.read(buffer, size);
	buffer[size] = 0;
	m_modelInfo.taxonomy.genus = buffer;

	fin.read((char*)&size, sizeof(size_t));
	fin.read(buffer, size);
	buffer[size] = 0;
	m_modelInfo.taxonomy.species = buffer;

	fin.read((char*)&size, sizeof(size_t));
	fin.read(buffer, size);
	buffer[size] = 0;
	m_modelInfo.taxonomy.strain = buffer;

	m_kmerCalculator = new KmerCalculator(m_wordLength);
	m_logConditionalProb = new float[m_kmerCalculator->numPossibleWords()];

	fin.read((char*)m_logConditionalProb, sizeof(float) * m_kmerCalculator->numPossibleWords());

	fin.close();
}

void KmerModel::printModelInfo(std::ostream& out) const
{
	out << "Model name: " << m_modelInfo.name << std::endl;
	out << "Model taxonomy: " << m_modelInfo.taxonomy.taxonomyStr() << std::endl;
	out << "N-mer length: " << m_wordLength << std::endl;
}