# Classify fragments using NB-BL

MIN_SCORE = -1e100
TAXONOMY_TXT = 'taxonomy.txt'

import sys
import fileinput
import math

if len(sys.argv) != 4:
	print 'NB-BL v1.1 by Donovan Parks, Norm MacDonald, and Rob Beiko'
	print '           with optimizations by Jose R Valverde (CNB/CSIC)'
	print ''
	print 'Usage: python NB-BL.py <nb-results> <blastn-results> <results-file>'
	print ''
	print 'Required parameters:'
	print '  <blastn-path>   Results of NB classifier with T=0.'
	print '  <query-file>    Results of BLASTN classifier.'
	print '  <results-file>  File to write classification results to.'
	print ''
	print 'Typical usage:'
	print '  python NB-BL.py nb_results.txt blastn_results.txt nb-bl_results.txt'
	print ''
	exit()
	
nbResults = sys.argv[1]
blastnResults = sys.argv[2]
resultsFile = sys.argv[3]

# read in complete taxonomy of each strain/genome
def loadTaxonomy(taxonomy_file):
	# read taxonomy.txt and prepare accesstionToTaxonomy and 
	# strainToTaxonomy dictionaries
	try:
		tax = open(taxonomy_file, "r")
	except:
		print "Cannot open", taxonomy_file

	strainToTaxonomy = {}
	accessionToTaxonomy = {}

	for line in tax:
		lineSplit = line.split('\t')
		
		accession = lineSplit[0]
		taxonomy = lineSplit[1].split(';')
		
		accessionToTaxonomy[accession] = taxonomy
		
		strain = taxonomy[7].replace(':', '_').replace(' ', '_')
		strainToTaxonomy[strain] = taxonomy
		
	tax.close()
	return (accessionToTaxonomy, strainToTaxonomy)
	

accessionToTaxonomy, strainToTaxonomy = loadTaxonomy(TAXONOMY_TXT)

#### JR ###
# This could be even more efficient:
#
#	open blast results
#	open NB results
#	while there are sequences to be processed:
#   		while there are blast results for the same current sequence
#   			read all blast results for the current sequence
#   			read all NB results for this same current sequence
#   			consolidate them as below
#   			write our results
#
# Or, alternately, it may be easier to move blast processing
# down, i. e. as one reads NB results and finds a new sequence,
# read in all its corresponding blast results and do the combination.
# In a napkin sketch, sleight-of-hand, quick-and-dirty way it could
# look like
#
#	read one blastline
#	while there are NB lines:
#   		each line corresponds to a different fragment/read/sequence
#   		read in blast lines corresponding to this fragment/read/sequence
#   		assuming NB and blast are in the same order:
#   		if blastline.split('\t')[0] == fragmentId:
#   			# there are blast results for this NBid
#   			while blastline.split('\t')[0] == fragmentId:
#   				process blast line as below
#   				read one blastline
#   			# at this point blastline points to the next blast hit
#   		else:
#   			pass # there are no blast hits for this NB line/fragment/read/seq
#   		do NB+Blast combination as below
#
#### JR ### TO BE DONE (eventually)


# determine smallest blastn E-value for each strain (genome) for each fragment
#
print 'Processing blastn results...'

def loadBlastnEvalues(blastnResults, accessionToTaxonomy):
	# blastnEvalues is a dictionary containing the minimum e-value for each
	# strain matched by each sequence. I.e. it is a dictionary where each
	# entry corresponds to a query sequence, and each value is in turn a
	# dictionary of strains matched and their minimum e-value
	# This shouldn't be too big (several GB, approx 1K * N_sequences).
	blastnEvalues = {}

	try:
		bl = open(blastnResults, "r")
	except:
		print "Cannot open", blastnResults

	#	A blastn output file contains sequence hits and scores
	#	A given sequence may match more than one target from the same
	#		original species, as well as more than one species
	#	For each query sequence
	#		For each species it matches
	#			store only the best matching score
	#
	#	We don't want just the best overall score because when later
	#	it is combined with a LogLikelihood score, the combination
	#	may result in a differnt ranking order.
	for line in  bl:
		# ignore comments
		if line[0] == '#':
			continue
			
		# parse line and extract ID, accno and e-value
		lineSplit = line.split('\t')
		queryId = lineSplit[0]
		accession = lineSplit[1]
		evalue = float(lineSplit[10])
		
		# lookup taxonomy from accno and get strain
		taxonomy = accessionToTaxonomy[accession]
		strain = accessionToTaxonomy[accession][7]
		
		# get key associated to ID or an empty dictionary {} instead
		#   The first time blastnEvalues is empty
		evalues = blastnEvalues.get(queryId, {})
		
		# if this strain has already been seen
		if strain in evalues:
			# update its evalue to the minimum of (old, new)
			if evalue < evalues[strain]:
				evalues[strain] = evalue
		else:
			# otherwise, just enter its evalue
			evalues[strain] = evalue
		
		# update blastnEvalues	
		blastnEvalues[queryId] = evalues
		
	bl.close()
	return blastnEvalues
	

# use blastn results to associate each sequence with all the
#   species it matches and the best score obtained for each (of 
#   possibly many) species match
blastnEvalues = loadBlastnEvalues(blastnResults, accessionToTaxonomy)


#
# get NB likelihood for each fragment
#
print 'Processing NB results, combining with BL and writing out results...'

# prepare output file
try:
	fout = open(resultsFile, 'w')
except:
	print "Cannot open", resultsFile

# print header
fout.write('Fragment Id\tTaxonomic classification\n')

# nbLogLikelihoods is a dictionary containing the LogLikelihood
# of each query sequence being each given strain for all strains
nbLogLikelihoods = {}
bHeaderLine = True
strains = []

# Format of NB file is
#	FragmentID	Length	Valid_n-mers	Species1	Species2 ... SpeciesN
#	Each fragment uses one line with the scores for all tested species
try:
	nb = open(nbResults, "r")
except:
	print "Cannot open", nbResults

# work fragment by fragment
#	This reduces memory footprint to proportional to the number of
#	species considered plus the blastnEvalues dictionary instead
#	of the number of fragments BY the number of species (which would
#	give a huge number) plus the blastnEvalues dictionary
#
#	Cost could be further reduced if we sorted BlastN results and
#	NB results by fragmentID prior to matching, so that we would
#	be able to open BlastN as well and work on only one fragment
#	with both at a time (instead of keeping all BlastN results in
#	memory).
#
#	That should be easy to do with a grep(1)/sort(1) combination. We
#	could get rid of comments in blast since we need not worry about
#	them.
for line in nb:

	lineSplit = line.split('\t')
	
	# The first line is special: it contains the header
	#	with the names of all tested species (from column 4 [3] on)
	# Retrieve the name of all the species tested by NB
	if bHeaderLine:
		for i in xrange(3, len(lineSplit)):
			strains.append(lineSplit[i].strip())
		bHeaderLine = False
		continue
	
	# we expect each line to correspong to a different fragmentId
	fragmentId = lineSplit[0]
	
	# obtain likelihoods (they are all in the same line)
	logLikelihoods = {}
	for i in xrange(3, len(lineSplit)):
		strain = strains[i-3]
		logLikelihoods[strain] = float(lineSplit[i])
		# the following could be added here as well saving loops
		
	# Combine with NB-BL to obtain actual scores
	scores = {}
	for strain in logLikelihoods:
		logLikelihood = logLikelihoods[strain]
		
		# Here we could read in the corresponding blast results on
		# the fly instead of keeping all of them in memory.
		# Assuming we have already read a blast line:
		#
		# 	if blastline.split('\t')[0] == fragmentId:
		#   		while blast.split('\t')[0] == fragmentID:
		#   			blastnEvalues = {}
		#			... process blast line as above			
		#   			read blastline
		
		if (fragmentId not in blastnEvalues) or (strain not in blastnEvalues[fragmentId]):
			# just take NB score
			scores[strain] = logLikelihood
		else:
			evalue = blastnEvalues[fragmentId][strain]
			
			if evalue == 0:
				scores[strain] = logLikelihood + 10000	# big bonus to score
			else:
				scores[strain] = logLikelihood + 12*(4 - math.log(evalue)) # seems like this should be log10, but PhymmBL uses the natural log
		
	# find top score for this fragment
	topScore = MIN_SCORE
	for strain in scores:
		if scores[strain] > topScore:
			topTaxonomy = strainToTaxonomy[strain]
			topScore = scores[strain]
			
	# print out top score for this fragment
	fout.write(fragmentId + '\t')
	for r in xrange(0, 8):
		fout.write(topTaxonomy[r] + ';')
		fout.write("\t")
		for r in xrange(0, 8):
			fout.write("%f;"%(topScore))
	fout.write('\n')
	
nb.close()
fout.close()

print 'Done.'

