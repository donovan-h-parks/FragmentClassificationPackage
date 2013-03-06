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

#include "Utils.hpp"

std::string numberToStr(int number)
{
	std::stringstream out;
	out << number;
	
	return out.str();
}

std::string numberToStr(uint number)
{
	std::stringstream out;
	out << number;
	
	return out.str();
}

std::string numberToStr(float number)
{
	std::stringstream out;
	out << number;
	
	return out.str();
}