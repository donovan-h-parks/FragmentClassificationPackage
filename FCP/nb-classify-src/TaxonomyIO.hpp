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

#ifndef TAXONOMY_IO
#define TAXONOMY_IO

#include "stdafx.h"

class TaxonomyIO 
{
public:
	TaxonomyIO();
	TaxonomyIO(const std::string& filename);

	~TaxonomyIO();

	bool open(const std::string& filename);

	uint numCategories(TAXONOMIC_RANK rank);

	std::string category(const std::string& seqId, TAXONOMIC_RANK rank);

	uint numTaxa(const std::string& taxa);

	TaxonomyModel taxonomy(const std::string& seqId) const;

	std::vector<uint> numChildren(TaxonomyModel queryModel);

private:
	void parse();
	void releaseMemory();

private:
	std::ifstream m_fileStream;

	std::map<std::string, TaxonomyModel> m_seqToTaxonomyMap;
	std::map<std::string, uint> m_taxaCountMap;

	char* m_buffer;
};

#endif