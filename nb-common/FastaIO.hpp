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

#ifndef FASTA_IO
#define FASTA_IO

#include "stdafx.h"

class FastaIO 
{
public:
	static const long BUFFER_SIZE = 500 * 1024 * 1024;	// 100 Mbyte buffer

public:
	FastaIO();
	FastaIO(const std::string& filename);

	~FastaIO();

	bool open(const std::string& filename);

	bool nextSeq(SeqInfo& seqInfo, bool bKeepAlignment = true);

	float percentageProcessed() const;

	void setToFirstSeq();

	bool readSeqs(const std::string& filename, std::vector<SeqInfo>& seqInfo, char* buffer, uint verbose = 0, bool bKeepAlignment = true);

private:
	void nextSeq(SeqInfo& seqInfo, const bool bKeepAlignment, char* buffer, ulong& bufferPos, const ulong bytesInBuffer);

	bool fillBuffer();

	void initialize();
	void releaseMemory();

private:

	std::ifstream m_fileStream;

	ulong m_fileSize;
	ulong m_bytesRead;

	char* m_buffer;
	ulong m_bufferPos;
	ulong m_bytesInBuffer;
};

#endif