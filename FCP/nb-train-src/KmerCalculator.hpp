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

#ifndef KMER_CALCULATOR
#define KMER_CALCULATOR

#include "stdafx.h"

class KmerCalculator
{
public:
	static const byte INVALID_NT_CHARACTER = 255;

public:
	KmerCalculator(uint wordLength);
	~KmerCalculator();

	ulong numPossibleWords() const { return m_numPossibleWords; }
	ulong topMultiplier() const { return m_topMultiplier; }

	byte ntValue(char c) const { return m_ntValues[c]; }

	void extractKmers(SeqInfo& seqInfo, std::vector<uint>& kmerValues);
	void extractForwardKmers(SeqInfo& seqInfo, std::vector<uint>& kmerValues);

	void baseFrequencies(SeqInfo& seqInfo, std::vector<float>& baseFrequencies, ulong& numValidBases);

private:
	byte* m_ntValues;
	byte* m_ntReverseValues;
	uint* m_kmerWordCount;

	uint m_wordLength;

	ulong m_numPossibleWords;
	ulong m_topMultiplier;

	uint m_bitShift;
};

#endif