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

import sys
import os
import fileinput
import platform

if len(sys.argv) != 11:
	print 'AddModel v1.0.1 by Donovan Parks, Norm MacDonald, and Rob Beiko'
	print ''
	print 'Usage: python AddModel.py <N> <sequence-file> <domain> <phylum> <class> <order> <family> <genus> <species> <strain>'
	print ''
	print 'Required parameters:'
	print '  <N>              Desired n-mer length (must be the same for all models).'
	print '  <sequence-file>  Multi-FASTA file containing sequence data to build model from.'
	print '  <domain>         Domain of model.'
	print '  ...'
	print '  <stain>          Strain of model.'
	print ''
	print 'Typical usage:'
	print '  python AddModel.py 8 ./training/custom/MySeqData.fasta Bacteria Proteobacteria'
	print '            Betaproteobacteria Burkholderiales Burkholderiaceae Ralstonia '
	print '            "Ralstonia pickettii" "Ralstonia pickettii 12D"'
	print ''
	exit()

N = sys.argv[1]
fastaFile = sys.argv[2]
domain = sys.argv[3]
phylum = sys.argv[4]
classRank = sys.argv[5]
order = sys.argv[6]
family = sys.argv[7]
genus = sys.argv[8]
species = sys.argv[9]
strain = sys.argv[10]

taxonomy = [domain, phylum, classRank, order, family, genus, species, strain]

# get id for each sequence in new sequence file
seqIds = []
for line in fileinput.input([fastaFile]): 
	if line[0] == '>':
		id = line[1:].strip()
		seqIds.append(id)

# add new sequence to taxonomy file
print 'Adding sequence to taxonomy file (taxonomy.txt)...'
fout = open('taxonomy.txt', 'a')
for id in seqIds:
	fout.write(id + '\t')
	
	for r in xrange(0, len(taxonomy)):
		fout.write(taxonomy[r] + ';')
	fout.write('\n')
fout.close()

# add new sequence to sequence file
print 'Adding sequence to sequence file (sequence.txt)...'
fout = open('./training/sequences.txt', 'a')
fout.write('.' + fastaFile + '\n')
fout.close()

# write temporary sequence file
fout = open('./training/__temp__.txt', 'w')
fout.write(fastaFile + '\n')
fout.close()

# build new model
print 'Building model for new sequence...'
if platform.system() == 'Windows':
	os.system('nb-train-windows.exe -n %s -s ./training/__temp__.txt -m ./models/genomes/' % N)
else: # assume the system can build the executable from source
	os.system('./nb-train -n %s -s ./training/__temp__.txt -m ./models/genomes/' % N)
	
# removing temporary sequence file
os.remove('./training/__temp__.txt')

# add new model to model file
print 'Adding new model to model file (models.txt)...'
fout = open('./models/models.txt', 'a')

filename = fastaFile[fastaFile.rfind('/')+1:fastaFile.rfind('.')]
fout.write('./models/genomes/' + filename + '.txt' + '\n')
fout.close()

print 'Done.'

