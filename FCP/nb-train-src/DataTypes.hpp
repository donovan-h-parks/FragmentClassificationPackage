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

#ifndef DATA_TYPES
#define DATA_TYPES

typedef unsigned char byte;
typedef unsigned int uint;
typedef unsigned long ulong;

const uint NUM_TAXONOMIC_RANKS = 7;
enum TAXONOMIC_RANK { SUPERKINGDOM_RANK, PHYLUM_RANK, CLASS_RANK, ORDER_RANK, FAMILY_RANK, GENUS_RANK, SPECIES_RANK, STRAIN_RANK };

struct TaxonomyModel
{
	std::string taxonomyStr() const
	{
		return superKingdom + ";" + phylum + ";" + taxonomicClass + ";" + order 
						+ ";" + family + ";" + genus + ";" + species + ";" + strain + ";";
	}

	std::string category(TAXONOMIC_RANK rank) const
	{
		if(rank == SUPERKINGDOM_RANK)
			return superKingdom;
		else if(rank == PHYLUM_RANK)
			return phylum;
		else if(rank == CLASS_RANK)
			return taxonomicClass;
		else if(rank == ORDER_RANK)
			return order;
		else if(rank == FAMILY_RANK)
			return family;
		else if(rank == GENUS_RANK)
			return genus;
		else if(rank == SPECIES_RANK)
			return species;
		else if(rank == STRAIN_RANK)
				return strain;
		else
			return "";
	}

	std::string superKingdom;
	std::string phylum;
	std::string taxonomicClass;
	std::string order;
	std::string family;
	std::string genus;
	std::string species;
	std::string strain;
};

struct SeqInfo
{
	TaxonomyModel taxonomy;

	const char* seqId;

	const char* seq;		
	ulong length;
	ulong validKmers;
};

#endif