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

# This scrip downloads all complete sequenced bacterial and
# archaeal genomes in the NCBI RefSeq database. Due to the size
# of this file (> 3.5GB) it can take several hours to download.
# Taxonomic information is also obtained by downloading a portion
# of the NCBI Taxonomy database.

import os
import sys
import platform
import tarfile
import urllib
import fileinput
from ftplib import FTP

#***
PROTEIN_FILE = ''
#***

# Write a file to disk obtained via FTP
class FtpWriter:
	def __init__(self, file):
		self.f = open(file, 'wb')
		self.count = 0

	def __call__(self, block):
		self.f.write(block)
		
		if self.count % 1000 == 0:
			print '.',
			sys.stdout.flush()
			
		self.count += 1
		
	def close(self):
		self.f.close()

# Download bacterial and archaeal genomes in NCBI RefSeq
def DownloadGenomes(genomeFile):
	bDownload = True
	if os.path.exists('./' + genomeFile):
		bValidResponse = False
		while not bValidResponse:
			response = raw_input('NCBI genome file ' + genomeFile + ' already exists. Would you like to download the latest version [Y/N]? ')
			if response[0] == 'Y' or response[0] == 'y':
				bDownload = True
				bValidResponse = True
			elif response[0] == 'N' or response[0] == 'n':
				bDownload = False
				bValidResponse = True
			
	if bDownload:
		ncbiFTP = 'ftp.ncbi.nih.gov'
		genomeDir = '/genomes/Bacteria'
		
		# Connect to NBCI's FTP site using an anonymous account
		print 'Connecting to NCBI FTP site (' + ncbiFTP + ')...'
		ftp = FTP(ncbiFTP)
		print ftp.login()
		print '\n'

		# Change to directory containing bacterial and archaeal genomes
		print 'Changing to directory ' + genomeDir
		print ftp.cwd(genomeDir)
		print '\n'

		# Download bacterial and archaeal genomes
		print 'Downloading bacterial and archaeal genomes (' + genomeFile + ')...'
		print '  This file is ~4GB and may take awhile to download.'
		ftpWriter = FtpWriter(genomeFile)
		msg = ftp.retrbinary('RETR ' + genomeFile, ftpWriter, 32*1024*1024)
		print '\n'
		print msg
		ftpWriter.close()

		ftp.quit()
	
# Download NCBI taxonomy database
def DownloadTaxonomy(taxonomyDump):
	bDownload = True
	if os.path.exists('./' + taxonomyDump):
		bValidResponse = False
		while not bValidResponse:
			response = raw_input('NCBI taxonomy file ' + taxonomyDump + ' already exists. Would you like to download the latest version [Y/N]? ')
			if response[0] == 'Y' or response[0] == 'y':
				bDownload = True
				bValidResponse = True
			elif response[0] == 'N' or response[0] == 'n':
				bDownload = False
				bValidResponse = True
				
	if bDownload:
		ncbiFTP = 'ftp.ncbi.nih.gov'
		taxonomyDir = '/pub/taxonomy'
		
		# Connect to NBCI's FTP site using an anonymous account
		print 'Connecting to NCBI FTP site (' + ncbiFTP + ')...'
		ftp = FTP(ncbiFTP)
		print ftp.login()
		print '\n'

		# Change to directory containing taxonomy files
		print 'Changing to directory ' + taxonomyDir
		print ftp.cwd(taxonomyDir)
		print '\n'

		# Download taxonomy files
		print 'Downloading taxonomy database files...'
		print '  It may take a few minutes to download these files.'
		
		ftpWriter = FtpWriter(taxonomyDump)
		msg = ftp.retrbinary('RETR ' + taxonomyDump, ftpWriter, 32*1024*1024)
		print '\n'
		print msg
		ftpWriter.close()
		
		ftp.quit()

# Decompress genome file
def DecompressGenomes(genomeFile):
	tar = tarfile.open(genomeFile, 'r:gz')
	tar.extractall('./ncbi_genomes/')
	tar.close()
	
# Decompress taxonomy files
def DecompressTaxonomy(taxonomyDump):
	tar = tarfile.open(taxonomyDump, 'r:gz')
	tar.extractall('./taxonomy/')
	tar.close()
	
# Get full taxonomy of all prokaryotes
def BuildTaxonomyFile():
	# read taxon Id number of all contigs
	print 'Extracting taxon Id from each contig...'
	assessionToTaxonId = {}
	accessionToSource = {}
	genomeDirs = os.listdir('./ncbi_genomes/')
	for dir in genomeDirs:
		if not os.path.isdir('./ncbi_genomes/' + dir):
			continue
			
		for filename in os.listdir('./ncbi_genomes/' + dir):
			if os.path.isdir('./ncbi_genomes/' + dir + '/' + filename):
				continue
				
			accession = filename.split('.')[0]
			for line in fileinput.input(['./ncbi_genomes/' + dir + '/' + filename]): 
				if 'SOURCE' in line:
					source = line[len('SOURCE'):].strip()
					accessionToSource[accession] = source.replace('/', '_')
				if '/db_xref="taxon:' in line:
					taxonId = line.split(':')[1]
					taxonId = int(taxonId[0:taxonId.rfind('\"')])
					assessionToTaxonId[accession] = taxonId
					fileinput.close()
					
					break
	
	print 'Number of contigs: ' + str(len(assessionToTaxonId))

	# extract taxonomy of each contig
	print 'Extracting taxonomy of each contig...'

	nodeIdToName = {}
	for line in fileinput.input(['./taxonomy/names.dmp']): 
		lineSplit = line.split('|')
		id = int(lineSplit[0])
		name = lineSplit[1].strip()
		type = lineSplit[3].strip()
		
		if type == 'scientific name':
			nodeIdToName[id] = name
		
	taxonIdToNode = {}
	for line in fileinput.input(['./taxonomy/nodes.dmp']): 
		lineSplit = line.split('|')
		taxonId = int(lineSplit[0])
		parentId = int(lineSplit[1])
		rank = lineSplit[2].strip()
		
		taxonIdToNode[taxonId] = [rank, parentId]
				
	ranks = ['strain', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']
	fout = open('taxonomy.txt', 'w')
	missingTaxonIds = {}
	for assession in assessionToTaxonId:
		taxonId = assessionToTaxonId[assession]
		source = accessionToSource[assession]

		taxonomy = ['','','','','','','','']
		rankIndex = 0
		
		while True:
			if taxonId not in nodeIdToName:
				bKnownTaxon = False
				break
			elif nodeIdToName[taxonId] == 'root':
				bKnownTaxon = True
				break
				
			node = taxonIdToNode[taxonId]
			if node[0] in ranks:
				while rankIndex < ranks.index(node[0]):
					if rankIndex != 0:
						taxonomy[rankIndex] = nodeIdToName[taxonId] + ' (' + ranks[rankIndex] + ')'
					else:
						taxonomy[rankIndex] = source
					rankIndex += 1
					
				taxonomy[ranks.index(node[0])] = nodeIdToName[taxonId]
				rankIndex += 1
			
			taxonId = node[1]
		
		if bKnownTaxon:
			fout.write(assession + '\t')
			for r in xrange(7, -1, -1):
				fout.write(taxonomy[r] + ';')
			fout.write('\n')
		else:
			missingTaxonIds[taxonId] = source

	fout.close()
	
	return missingTaxonIds

# create genome-level input files
def CreateStrainSeqFiles():
	# determine genome of each sequence
	assessionToGenome = {}
	for line in fileinput.input(['taxonomy.txt']): 
		lineSplit = line.split('\t')
		
		seqId = lineSplit[0]
		category = lineSplit[1].split(';')[7]

		assessionToGenome[seqId] = category
		
	# creat required directories
	if not os.path.exists('./training'):
		os.makedirs('./training')
		os.makedirs('./training/sequences')
		os.makedirs('./training/custom')
	
	# remove any previously created models
	for assession in assessionToGenome:
		genome = assessionToGenome[assession]
		genomeFile = genome.replace(' ', '_')
		genomeFile = genomeFile.replace(':', '_')
		genomeFile += '.fasta'
		
		if os.path.exists('./training/sequences/' + genomeFile):
			os.remove('./training/sequences/' + genomeFile)
		
	#***
	if PROTEIN_FILE:
		fprotein = open(PROTEIN_FILE,'w')
	#***
		
	# convert genbank files to fasta files
	genomeDirs = os.listdir('./ncbi_genomes/')
	for dir in genomeDirs:
		if not os.path.isdir('./ncbi_genomes/' + dir):
			continue
			
		for filename in os.listdir('./ncbi_genomes/' + dir):
			if os.path.isdir('./ncbi_genomes/' + dir + '/' + filename):
				continue
				
			fullFilename = './ncbi_genomes/' + dir + '/' + filename
			
			# read sequence data from genbank file
			data = open(fullFilename).read()
			origin = data.rfind('ORIGIN')
			start = data.find('1', origin)
			end = data.find('//', origin)
			
			seqLines = data[start:end].split('\n')
			
			seq = ''
			for line in seqLines:
				subseq = line.split()
				seq += ''.join(subseq[1:])
			
			# write fasta file
			assession = filename.split('.')[0]
			
			if assession in assessionToGenome:
				genome = assessionToGenome[assession]
			
				print assession

				genomeFile = genome.replace(' ', '_')
				genomeFile = genomeFile.replace(':', '_')
				
				fout = open('./training/sequences/' + genomeFile + '.fasta', 'a')
				fout.write('>' + assession + '\n')
				
				index = 0
				while index+60 < len(seq):
					fout.write(seq[index:index+60] + '\n')
					index += 60
				fout.write(seq[index:] + '\n')
				fout.close()
			
			#***
			if PROTEIN_FILE:
				plen = len('/protein_id="')
				tlen = len('/translation="')
				cdsPosition = data.find('CDS')
				while cdsPosition >=0:
					proteinIdPosition = data.find('/protein_id',cdsPosition)
					proteinId = data[proteinIdPosition+plen:data.find('"',proteinIdPosition+plen+1)]
					translationPosition = data.find('/translation',cdsPosition)
					translation = data[translationPosition+tlen:data.find('"',translationPosition+tlen+1)].replace(' ','')
					fprotein.write('>%s\t%s\t%s\n%s\n'%(assession,proteinId,genome,translation))
					cdsPosition = data.find('CDS',cdsPosition+1)
			#***
			
	
	# create training file for genome models
	trainingSet = open('./training/sequences.txt', 'w')
	for filename in os.listdir('./training/sequences/'):
		trainingSet.write('./training/sequences/' + filename + '\n')
	trainingSet.close()
	
	#***
	if PROTEIN_FILE:
		fprotein.close()
	
	
# Build Naive Bayes models
def BuildNaiveBayesModels():
	# creat required directories
	if not os.path.exists('./models'):
		os.makedirs('./models')
		os.makedirs('./models/genomes')
		
	# build stain-level models
	if platform.system() == 'Windows':
		print 'Building genome-level models...'
		os.system('nb-train-windows.exe -n 8 -s ./training/sequences.txt -m ./models/genomes/')
	else: # assume the system can build the executable from source
		print 'Compiling nb-train...'
		os.chdir('./nb-train-src')
		os.system('make')
		os.chdir('..')
		os.system('cp ./nb-train-src/nb-train .')
		
		print 'Compiling nb-classify...'
		os.chdir('./nb-classify-src')
		os.system('make')
		os.chdir('..')
		os.system('cp ./nb-classify-src/nb-classify .')
		
		print 'Building genome-level models...'
		os.system('./nb-train -n 8 -s ./training/sequences.txt -m ./models/genomes/')

	# create model file for classifying query fragments
	modelFile = open('./models/models.txt', 'w')
	for line in fileinput.input(['./training/sequences.txt']): 
		genome = line[line.rfind('/')+1:line.rfind('.')]
		modelFile.write('./models/genomes/' + genome + '.txt' + '\n')
	modelFile.close()

genomeFile = 'all.gbk.tar.gz'
taxonomyDump = 'taxdump.tar.gz'

print 'FCP_install v1.0.5'
print ''
print 'This script is maintained by Donovan Parks (donovan.parks@gmail.com), Norm MacDonald, and Rob Beiko.'
print ''

#***
for flagi, flag in enumerate(sys.argv):
	if flag.lower() == '--protein':
		PROTEIN_FILE = sys.argv[flagi+1]
#***

print '\n'
print 'Downloading bacterial and archaeal genomes from NCBI:'
DownloadGenomes(genomeFile)
print '\n----------------------------------------------------------\n'

print 'Decompressing genomes:'
DecompressGenomes(genomeFile)
print '\n----------------------------------------------------------\n'

print 'Downloading NCBI taxonomy database:'
DownloadTaxonomy(taxonomyDump)
print '\n----------------------------------------------------------\n'

print 'Decompressing taxonomy files:'
DecompressTaxonomy(taxonomyDump)
print '\n----------------------------------------------------------\n'

print 'Building taxonomy file for genomes:'
missingTaxonIds = BuildTaxonomyFile()
print '\n----------------------------------------------------------\n'

print 'Creating input sequence file for each genome:'
CreateStrainSeqFiles()
print '\n----------------------------------------------------------\n'

print 'Building Naive Bayes models for each genomes:'
BuildNaiveBayesModels()
print '\n----------------------------------------------------------\n'

if len(missingTaxonIds) != 0:
	print '\n-- Warnings --------------------------------------------\n'
	for taxonId in missingTaxonIds:
		source = missingTaxonIds[taxonId]
		print "No NCBI taxon information for taxon id '" + str(taxonId) + "'."
		print " This taxon id is associated with the following source: " + source
		print " No model was build for this taxon."
		print ''
	print 'Installation complete (' + str(len(missingTaxonIds)) + ' warning, 0 errors).'
else:
	print 'Installation complete (0 warnings, 0 errors). '
