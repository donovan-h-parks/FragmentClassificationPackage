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

#ifndef KMER_MODEL
#define KMER_MODEL

#include "stdafx.h"

#include "KmerCalculator.hpp"

class KmerModel
{
public:
	static const byte INVALID_NT_CHARACTER = 255;

public:
	KmerModel(uint wordLength);
	KmerModel(const std::string& modelFile);

	~KmerModel();

	bool constructModel(SeqInfo& seqInfo);

	void calculateConditionalProbabilities();

	float classify(SeqInfo& seqInfo);
	float classify(const SeqInfo& seqInfo, const std::vector<uint>& profile);

	void write(const std::string& filename) const;
	
	TaxonomyModel taxonomy() const { return m_modelInfo.taxonomy; }

	void name(const std::string& name) { m_modelInfo.name = name; }
	std::string name() const { return m_modelInfo.name; }
	
	uint kmerLength() const { return m_wordLength; }

	void printModelInfo(std::ostream& out) const;

private:
	void read(const std::string& filename);

private:
	struct ModelInfo
	{
		ModelInfo(): numWords(0), numSeqs(0) {}

		ModelInfo(const std::string _name, TaxonomyModel _taxonomy, uint _numWords, uint _numSeqs)
			: name(_name), taxonomy(_taxonomy), numWords(_numWords), numSeqs(_numSeqs) {}

		std::string name;

		TaxonomyModel taxonomy;

		ulong numWords;
		ulong numSeqs;
	};

private:
	ModelInfo m_modelInfo;

	KmerCalculator* m_kmerCalculator;

	float* m_logConditionalProb;

	uint m_wordLength;
};

#endif