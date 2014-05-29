#!/usr/bin/env python
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

# This script classifies query fragments using BLASTN

import sys, os
import fileinput

if len(sys.argv) != 4:
	print 'BLASTN v1.0 by Donovan Parks, Norm MacDonald, and Rob Beiko'
	print 'Usage: python BLASTN.py <blastn-path> <query-file> <results-file>'
	print ''
	print 'Required parameters:'
	print '  <blastn-path>   Full path to blastn.'
	print '  <query-file>    Multi-FASTA file containing query fragments to classify.'
	print '  <results-file>  File to write classification results to.'
	print ''
	print 'Typical usage:'
	print '  python BLASTN.py /path/to/blastn test.fasta blastn_results.txt'
	print ''
	exit()

blastnEXE = sys.argv[1]
queryFile = sys.argv[2]
resultsFile = sys.argv[3]

print 'Running BLASTN...'

os.system(blastnEXE + ' -query ' + queryFile + ' -db ./blast_data/BacteriaAndArchaeaGenomesDB -evalue 10 -outfmt 7 -task blastn -out ' + resultsFile)

print 'Done.'