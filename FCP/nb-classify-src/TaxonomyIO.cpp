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

#include "TaxonomyIO.hpp"

TaxonomyIO::TaxonomyIO()
{
	m_buffer = NULL;
}

TaxonomyIO::TaxonomyIO(const std::string& filename) 
{ 
	m_buffer = NULL;

	open(filename); 
}

TaxonomyIO::~TaxonomyIO() 
{
	releaseMemory();
}

void TaxonomyIO::releaseMemory()
{
	if(m_buffer != NULL)
		delete[] m_buffer;
	
	m_buffer = NULL;
}

bool TaxonomyIO::open(const std::string& filename)
{
	m_fileStream.open(filename.c_str(), std::ios::in | std::ios::binary);
	if(!m_fileStream.is_open())
		return false;

	// calculate size of file
	m_fileStream.seekg(0, std::ios::end);
	long fileSize = m_fileStream.tellg();
	m_fileStream.seekg(0, std::ios::beg);

	// read file into buffer
	releaseMemory();
	m_buffer = new char[fileSize];
	m_fileStream.read(m_buffer, fileSize);
	if(m_fileStream.fail() || m_fileStream.bad())
		return false;
	
	// parse file
	parse();

	return true;
}

void TaxonomyIO::parse()
{
	char* curPos = m_buffer;

	while(true)
	{
		// read sequence id and taxonomy string
		char* seqIdStart = curPos;
		char* seqIdEnd = (char *)strchr(curPos, '\t');
		if(seqIdEnd == NULL)
			break;

		char* taxonomyStrStart = seqIdEnd + 1;
		char* taxonomyStrEnd = (char *)strchr(seqIdEnd, '\n');
		if(taxonomyStrEnd == NULL)
			break;

		curPos = taxonomyStrEnd+1;

		// null terminate these strings
		seqIdEnd[0] = 0;
		taxonomyStrEnd[0] = 0;
		
		// parse taxonomy information
		TaxonomyModel taxonomyModel;

		char* categoryStart = taxonomyStrStart;
		int modelLevel = 0;
		for(char* ch = taxonomyStrStart; ch < taxonomyStrEnd; ++ch)
		{
			if(*ch == ';')
			{
				*ch = 0;
				
				std::map<std::string, uint>::iterator it = m_taxaCountMap.find(categoryStart);
				if(it == m_taxaCountMap.end())
					m_taxaCountMap[categoryStart] = 1;
				else
					m_taxaCountMap[categoryStart] = it->second + 1;

				if(modelLevel == 0)
					taxonomyModel.superKingdom = categoryStart;
				else if(modelLevel == 1)
					taxonomyModel.phylum = categoryStart;
				else if(modelLevel == 2)
					taxonomyModel.taxonomicClass = categoryStart;
				else if(modelLevel == 3)
					taxonomyModel.order = categoryStart;
				else if(modelLevel == 4)
					taxonomyModel.family = categoryStart;
				else if(modelLevel == 5)
					taxonomyModel.genus = categoryStart;
				else if(modelLevel == 6)
					taxonomyModel.species = categoryStart;
				else if(modelLevel == 7)
					taxonomyModel.strain = categoryStart;

				categoryStart = ch + 1;
				modelLevel++;
			}
		}

		m_seqToTaxonomyMap[seqIdStart] = taxonomyModel;
	}
}

uint TaxonomyIO::numCategories(TAXONOMIC_RANK rank)
{
	std::set<std::string> categories;

	std::map<std::string, TaxonomyModel>::iterator it;
	for(it = m_seqToTaxonomyMap.begin(); it != m_seqToTaxonomyMap.end(); ++it)
	{
		TaxonomyModel model = it->second;
		if(rank == SUPERKINGDOM_RANK)
			categories.insert(model.superKingdom);
		else if(rank == PHYLUM_RANK)
			categories.insert(model.phylum);
		else if(rank == CLASS_RANK)
			categories.insert(model.taxonomicClass);
		else if(rank == ORDER_RANK)
			categories.insert(model.order);
		else if(rank == FAMILY_RANK)
			categories.insert(model.family);
		else if(rank == GENUS_RANK)
			categories.insert(model.genus);	
		else if(rank == SPECIES_RANK)
			categories.insert(model.species);	
		else if(rank == STRAIN_RANK)
			categories.insert(model.strain);	
	}

	return categories.size() - categories.count("");
}

std::string TaxonomyIO::category(const std::string& seqId, TAXONOMIC_RANK rank)
{
	std::map<std::string, TaxonomyModel>::iterator it;

	it = m_seqToTaxonomyMap.find(seqId);
	if(it != m_seqToTaxonomyMap.end())
	{
		if(rank == SUPERKINGDOM_RANK)
			return it->second.superKingdom;
		else if(rank == PHYLUM_RANK)
			return it->second.phylum;
		else if(rank == CLASS_RANK)
			return it->second.taxonomicClass;
		else if(rank == ORDER_RANK)
			return it->second.order;
		else if(rank == FAMILY_RANK)
			return it->second.family;
		else if(rank == GENUS_RANK)
			return it->second.genus;
		else if(rank == SPECIES_RANK)
			return it->second.species;
		else if(rank == STRAIN_RANK)
			return it->second.strain;
		else
			return "";
	}

	std::cout << "Taxonomic information for sequence is missing: " << seqId << std::endl;
	return "";
}

TaxonomyModel TaxonomyIO::taxonomy(const std::string& seqId) const 
{ 
	std::string id = seqId.substr(0, seqId.find_last_of('-'));
	std::map<std::string, TaxonomyModel>::const_iterator it = m_seqToTaxonomyMap.find(id);
	if(it != m_seqToTaxonomyMap.end())
		return it->second; 

	std::cerr << "Sequence has no taxonomic information!" << std::endl;

	return TaxonomyModel();
}

uint TaxonomyIO::numTaxa(const std::string& taxa)
{
	std::map<std::string, uint>::iterator it = m_taxaCountMap.find(taxa);
	if(it != m_taxaCountMap.end())
		return it->second;

	return 0;
}	

std::vector<uint> TaxonomyIO::numChildren(TaxonomyModel queryModel)
{
	std::set<std::string> phylum;
	std::set<std::string> classes;
	std::set<std::string> order;
	std::set<std::string> family;
	std::set<std::string> genus;
	std::set<std::string> species;
	std::set<std::string> strain;

	std::map<std::string, TaxonomyModel>::iterator it;
	for(it = m_seqToTaxonomyMap.begin(); it != m_seqToTaxonomyMap.end(); ++it)
	{
		TaxonomyModel model = it->second;

		if(queryModel.category(SUPERKINGDOM_RANK) == model.category(SUPERKINGDOM_RANK))
			phylum.insert(model.category(PHYLUM_RANK));

		if(queryModel.category(PHYLUM_RANK) == model.category(PHYLUM_RANK))
			classes.insert(model.category(CLASS_RANK));

		if(queryModel.category(CLASS_RANK) == model.category(CLASS_RANK))
			order.insert(model.category(ORDER_RANK));

		if(queryModel.category(ORDER_RANK) == model.category(ORDER_RANK))
			family.insert(model.category(FAMILY_RANK));

		if(queryModel.category(FAMILY_RANK) == model.category(FAMILY_RANK))
			genus.insert(model.category(GENUS_RANK));

		if(queryModel.category(GENUS_RANK) == model.category(GENUS_RANK))
			species.insert(model.category(SPECIES_RANK));

		if(queryModel.category(SPECIES_RANK) == model.category(SPECIES_RANK))
			strain.insert(model.category(STRAIN_RANK));
	}

	std::vector<uint> counts;
	counts.push_back(phylum.size());
	counts.push_back(classes.size());
	counts.push_back(order.size());
	counts.push_back(family.size());
	counts.push_back(genus.size());
	counts.push_back(species.size());
	counts.push_back(strain.size());

	return counts;
}