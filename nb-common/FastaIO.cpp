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

FastaIO::FastaIO()
{
	m_buffer = NULL;

	initialize();
}

FastaIO::FastaIO(const std::string& filename) 
{ 
	m_buffer = NULL;

	initialize();

	open(filename); 
}

void FastaIO::initialize()
{
	m_fileSize = 0;
	m_bytesRead = 0;

	m_bufferPos = 0;
	m_bytesInBuffer = 0;
}

void FastaIO::setToFirstSeq() 
{ 
	m_bytesRead = 0;
	m_bufferPos = 0;
	m_bytesInBuffer = 0;

	m_fileStream.clear();
}

FastaIO::~FastaIO() 
{
	releaseMemory();
	m_fileStream.close();
}

void FastaIO::releaseMemory()
{
	if(m_buffer != NULL)
		delete[] m_buffer;
	
	m_buffer = NULL;
}

bool FastaIO::open(const std::string& filename)
{
	initialize();

	m_fileStream.open(filename.c_str(), std::ios::in | std::ios::binary);
	if(!m_fileStream.is_open())
		return false;

	// calculate size of file
	m_fileStream.seekg(0, std::ios::end);
	m_fileSize = m_fileStream.tellg();
	m_fileStream.seekg(0, std::ios::beg);

	// allocate memory for reading file
	releaseMemory();
	m_buffer = new char[BUFFER_SIZE];

	return true;
}

bool FastaIO::nextSeq(SeqInfo& seqInfo, bool bKeepAlignment)
{
	if (m_bytesRead - m_bytesInBuffer + m_bufferPos == m_fileSize)
			return false;	// finished processing all sequences

	// determine if more sequences must be read from file
	if(m_bufferPos == m_bytesInBuffer)
	{
		if(!fillBuffer())
			return false;
	}

	nextSeq(seqInfo, bKeepAlignment, m_buffer, m_bufferPos, m_bytesInBuffer);
	
	return true;	
}

void FastaIO::nextSeq(SeqInfo& seqInfo, const bool bKeepAlignment, char* buffer, ulong& bufferPos, const ulong bytesInBuffer)
{
	// extract sequences id and null terminate
	char* seqIdStart = (char *)(buffer + bufferPos + 1);
	char* headerEnd = (char *)strchr(seqIdStart, '\n');
	char* seqIdEnd = (char *)memchr(seqIdStart, ' ', headerEnd-seqIdStart);
	
	if(seqIdEnd != NULL)
		seqIdEnd[0] = 0;	// null terminate at first space
	
	if(headerEnd == NULL)
	{
		// perhaps lines are terminated with just a carriage return (i.e., old mac style)
		headerEnd = (char *)strchr(seqIdStart, '\r');
	}

	// check if line is terminated with a \r\n which is typical under Microsoft Windows
	if(headerEnd[-1] == '\r')
		headerEnd[-1] = 0;
	else
		headerEnd[0] = 0;
		
	seqInfo.seqId = seqIdStart;

	// extract sequence (removing any whitespace characters and optionally any alignment characters)
	char *seqStart = headerEnd + 1;

	char *seqEnd = seqStart;
	bufferPos = (ulong) (seqStart - buffer);
	while(bufferPos < bytesInBuffer)
	{
		char c = buffer[bufferPos];
		if(c == '>')
			break;

		++bufferPos;
		if (isspace(c))
			continue;

		if(!bKeepAlignment && (c == '-' || c == '.'))
			continue;

		*seqEnd++ = c;
	}

	seqInfo.seq = seqStart;
	seqInfo.length = seqEnd - seqStart;
	
	// null terminate sequence data
	seqEnd[0] = 0;	
}

bool FastaIO::fillBuffer()
{
	m_bufferPos = 0;

	long bytesToRead = m_fileSize - m_bytesRead;
	if (bytesToRead > BUFFER_SIZE)
		bytesToRead = BUFFER_SIZE;

	m_fileStream.seekg(m_bytesRead, std::ios::beg);
	if(m_fileStream.fail() || m_fileStream.bad())
	{
		std::cerr << "Error reading FASTA file (seeking error)." << std::endl;
		return false;
	}

	m_fileStream.read(m_buffer, bytesToRead);
	if(m_fileStream.fail() || m_fileStream.bad())
	{
		std::cerr << "Error reading FASTA file." << std::endl;
		return false;
	}
	else if(m_buffer[0] != '>')
	{
		std::cerr << "Invalid FASTA file format." << std::endl;
		return false;
	}
	
	m_bytesInBuffer = bytesToRead;
	if(m_bytesRead + bytesToRead != m_fileSize)
	{
		// determine number of bytes until last complete sequence
		int i = 0;
		for(i = bytesToRead-1; i >= 0; i--)
		{
			if(m_buffer[i] == '>')
				break;
		}

		m_bytesInBuffer = i;
	}

	m_bytesRead += m_bytesInBuffer;		

	return true;
}

float FastaIO::percentageProcessed() const
{
	return (m_bytesRead + m_bufferPos - m_bytesInBuffer) * 100.0f / m_fileSize;
}

bool FastaIO::readSeqs(const std::string& filename, std::vector<SeqInfo>& seqInfo, char* buffer, uint verbose, bool bKeepAlignment)
{
	m_fileStream.open(filename.c_str(), std::ios::in | std::ios::binary);
	if(!m_fileStream.is_open())
		return false;

	// calculate size of file
	m_fileStream.seekg(0, std::ios::end);
	m_fileSize = m_fileStream.tellg();
	m_fileStream.seekg(0, std::ios::beg);

	// allocate memory for reading file
	buffer = new char[m_fileSize];

	// read file into buffer
	m_fileStream.read(buffer, m_fileSize);
	if(m_fileStream.fail() || m_fileStream.bad())
		return false;

	// determine position of each sequence in buffer
	ulong bufferPos = 0;
	while(bufferPos != m_fileSize)
	{
		SeqInfo info;
		nextSeq(info, bKeepAlignment, buffer, bufferPos, m_fileSize);
		seqInfo.push_back(info);

		if(verbose >= 3)
			std::cout << "Read fragment: " << info.seqId << std::endl;
	}

	return true;
}