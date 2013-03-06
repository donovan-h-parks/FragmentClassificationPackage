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

#include "KmerCalculator.hpp"

using namespace std;

enum ENCODING_TYPES { NUCLEOTIDE, YR_ENCODING };
const int gENCODING = NUCLEOTIDE;	

KmerCalculator::KmerCalculator(uint wordLength): m_wordLength(wordLength)
{
	if(gENCODING == NUCLEOTIDE)
	{
		m_bitShift = 2;

		m_numPossibleWords = 1 << (2*m_wordLength);

		m_topMultiplier = 1 << (2 * (m_wordLength-1));		// equivalent to pow(4, m_wordLength-1)

		m_kmerWordCount = new uint[m_numPossibleWords];

		// map nucleotides to values
		m_ntValues = new byte[256];
		memset(m_ntValues, INVALID_NT_CHARACTER, 256*sizeof(byte));

		m_ntValues['a'] = m_ntValues['A'] = 0;
		m_ntValues['c'] = m_ntValues['C'] = 1;
		m_ntValues['g'] = m_ntValues['G'] = 2;
		m_ntValues['t'] = m_ntValues['T'] = 3;
		m_ntValues['u'] = m_ntValues['U'] = 3;

		m_ntReverseValues = new byte[256];
		memset(m_ntReverseValues, INVALID_NT_CHARACTER, 256*sizeof(byte));

		m_ntReverseValues['a'] = m_ntReverseValues['A'] = 3;
		m_ntReverseValues['c'] = m_ntReverseValues['C'] = 2;
		m_ntReverseValues['g'] = m_ntReverseValues['G'] = 1;
		m_ntReverseValues['t'] = m_ntReverseValues['T'] = 0;
		m_ntReverseValues['u'] = m_ntReverseValues['U'] = 0;
	}
	else if(gENCODING == YR_ENCODING)
	{
		m_bitShift = 1;

		m_numPossibleWords = 1 << (m_wordLength);

		m_topMultiplier = 1 << (m_wordLength-1);		// equivalent to pow(2, m_wordLength-1)

		m_kmerWordCount = new uint[m_numPossibleWords];

		// map nucleotides to values
		m_ntValues = new byte[256];
		memset(m_ntValues, INVALID_NT_CHARACTER, 256*sizeof(byte));

		m_ntValues['a'] = m_ntValues['A'] = 1;
		m_ntValues['c'] = m_ntValues['C'] = 0;
		m_ntValues['g'] = m_ntValues['G'] = 1;
		m_ntValues['t'] = m_ntValues['T'] = 0;
		m_ntValues['u'] = m_ntValues['U'] = 0;

		m_ntReverseValues = new byte[256];
		memset(m_ntReverseValues, INVALID_NT_CHARACTER, 256*sizeof(byte));

		m_ntReverseValues['a'] = m_ntReverseValues['A'] = 0;
		m_ntReverseValues['c'] = m_ntReverseValues['C'] = 1;
		m_ntReverseValues['g'] = m_ntReverseValues['G'] = 0;
		m_ntReverseValues['t'] = m_ntReverseValues['T'] = 1;
		m_ntReverseValues['u'] = m_ntReverseValues['U'] = 1;
	}
}

KmerCalculator::~KmerCalculator()
{
	delete[] m_kmerWordCount;
	delete[] m_ntValues;
}

void KmerCalculator::extractForwardKmers(SeqInfo& seqInfo, std::vector<uint>& kmerValues)
{	
	seqInfo.validKmers = 0;
	if(seqInfo.length < m_wordLength)
		return;

	kmerValues.reserve(seqInfo.length - m_wordLength + 1);

	// calculate kmer number for initial window
	const char* seq = seqInfo.seq;
	ulong word = 0;
	int indexOfInvalidCharacter = -1;
	for(ulong i = 0; i < m_wordLength; ++i)
	{
		byte value = m_ntValues[seq[i]];
		if(value == INVALID_NT_CHARACTER)
			indexOfInvalidCharacter = m_wordLength;

		word = (word << m_bitShift) + value;

		indexOfInvalidCharacter = max(-1, indexOfInvalidCharacter - 1);
	}

	// record word
	if(indexOfInvalidCharacter == -1)
	{
		kmerValues.push_back(word);
		seqInfo.validKmers += 1;
	}

	// calculate new kmer number for each 1 nt shift of the window
	for(ulong i = m_wordLength; i < seqInfo.length; ++i)
	{
		// subtract value for nt leaving window
		word -= m_topMultiplier * m_ntValues[seq[i - m_wordLength]];

		// add value of nt entering window
		byte value = m_ntValues[seq[i]];

		if(value == INVALID_NT_CHARACTER)
			indexOfInvalidCharacter = m_wordLength;

		word = (word << m_bitShift) + value;

		indexOfInvalidCharacter = max(-1, indexOfInvalidCharacter - 1);

		// record word
		if(indexOfInvalidCharacter == -1)
		{
			kmerValues.push_back(word);
			seqInfo.validKmers += 1;
		}		
	}
}

void KmerCalculator::extractKmers(SeqInfo& seqInfo, std::vector<uint>& kmerValues)
{	
	kmerValues.reserve(2 * (seqInfo.length - m_wordLength + 1));

	// calculate kmer number for initial window
	const char* seq = seqInfo.seq;
	ulong word = 0;
	ulong reverseWord = 0;
	int indexOfInvalidCharacter = -1;
	seqInfo.validKmers = 0;
	for(ulong i = 0; i < m_wordLength; ++i)
	{
		byte value = m_ntValues[seq[i]];
		byte reverseValue = m_ntReverseValues[seq[m_wordLength-i-1]];

		if(value == INVALID_NT_CHARACTER)
			indexOfInvalidCharacter = m_wordLength;

		word = (word << m_bitShift) + value;
		reverseWord = (reverseWord << m_bitShift) + reverseValue;

		indexOfInvalidCharacter = max(-1, indexOfInvalidCharacter - 1);
	}

	// record word
	if(indexOfInvalidCharacter == -1)
	{
		kmerValues.push_back(word);
		kmerValues.push_back(reverseWord);

		seqInfo.validKmers += 2;
	}

	// calculate new kmer number for each 1 nt shift of the window
	for(ulong i = m_wordLength; i < seqInfo.length; ++i)
	{
		// subtract value for nt leaving window
		word -= m_topMultiplier * m_ntValues[seq[i - m_wordLength]];
		reverseWord -= m_ntReverseValues[seq[i - m_wordLength]];
		
		// add value of nt entering window
		byte value = m_ntValues[seq[i]];
		byte reverseValue = m_ntReverseValues[seq[i]];

		if(value == INVALID_NT_CHARACTER)
			indexOfInvalidCharacter = m_wordLength;

		word = (word << m_bitShift) + value;
		reverseWord = (reverseWord >> m_bitShift) + m_topMultiplier * reverseValue;

		indexOfInvalidCharacter = max(-1, indexOfInvalidCharacter - 1);

		// record word
		if(indexOfInvalidCharacter == -1)
		{
			kmerValues.push_back(word);
			kmerValues.push_back(reverseWord);

			seqInfo.validKmers += 2;
		}		
	}
}

void KmerCalculator::baseFrequencies(SeqInfo& seqInfo, std::vector<float>& baseFrequencies, ulong& numValidBases)
{
	// get base frequencies
	ulong atCount = 0;
	ulong gcCount = 0;
	for(uint i = 0; i < seqInfo.length; ++i)
	{
		if(seqInfo.seq[i] == 'A' || seqInfo.seq[i] == 'T')
		{
			atCount++;
			numValidBases += 2;
		}
		else if (seqInfo.seq[i] == 'C' || seqInfo.seq[i] == 'G')
		{
			gcCount++;	
			numValidBases += 2;
		}
	}

	// A => 0
	// C => 1
	// G => 2
	// T => 3
	baseFrequencies[0] += atCount;
	baseFrequencies[1] += gcCount;
	baseFrequencies[2] += gcCount;
	baseFrequencies[3] += atCount;
}

