# GC Skew Analysis Script
# C. Breuer 9/2021
#
# Script calculates the GCS using multiple calculations: 
#     By gene:
#         - all nucleotides ("whole gene") or only nucleotides in codon positions 1,2 or 3
#     By sliding window:
#         - across the full genome using a sliding window of adjustable length
# 
# Multiple Outputs:
#     One Gene Coordinates flat file with the original data plus columns for the various whole-gene GC skew values
#     Five Wiggle Files:
#         - Sliding window across the whole genome
#         - Whole gene:  1) all nucleotides, 2-4) Only codon positions 1,2, or 3
# 
# Notes:
#     - For whole gene GC skew values, the value is relative to the leading strand of the replication fork.
#     - This is where ter information is used - nucleotide 1 is assumed to be part of the origin, and ter is given
#       by the user. Note: ter coordinates can be looked up on the DoriC website: http://tubic.org/doric/public/index.php
#     - As a result, whole-gene GC skew values post-ter will have the opposite sign (-/+) as the sliding window value.
#     - Discussion of script function is at the bottom. 
#     - File/folder names and sliding window hyperparameters can be changed in the analyzeGCskew() function

import os
#################################################################################################
def getGenomesAndTers(GTfile, inputDataFolder):
	'''
	- Input is the genome/ter coordinates file name
	- The genome names in this file are used as a guide - script will then look for a <genome>.fasta and
	  a <genome>_coords.txt file specifying the start/stop/strand information for each gene feature.
	- Output is a list of genomes to process (from first column in the coords file)
	'''
	genomes = []
	try:
		filepath = inputDataFolder + "/" + GTfile
		with open(filepath, "r") as crd_file:
			for line in crd_file:
				splt = line.split("\t")

				# Start and Stop coordinates for each feature
				genomes.append((splt[0], splt[1]))
	except FileNotFoundError:
		print()
		print("# Error: Coult not find file path " + inputDataFolder + "/" + GTfile)
		print("# Exiting.")
		print()
		exit()
	return genomes

#################################################################################################
def get_coordinates(coords_file, ter, inFolder2):
	"""
	- Input is coordinates file
	- Output is coordinates list
	"""

	coordinates_list = []
	linenum = 0
	with open(inFolder2 + "/" + coords_file, "r") as cf:
		coordFile =  cf.readlines()
		for line in coordFile:
			error = "OK"
			linenum += 1

			if linenum == 1:
				continue 
			else:
				try:
					splt1 = line.split("\t")
					try:
						start = int(splt1[0])
						end = int(splt1[1])
						if start > end:						# The first gene sometimes spans the origin. This throws off the calc. switch to base pair 1.
							start = 1
					except ValueError:						# Some coordinate lines will be "Error". Throw these out.
						error = "Errant data point"

					#ori = splt1[17].strip("\n")
					sys = splt1[5].strip("\n")
					strand = splt1[2]
					name = splt1[4]
					func = splt1[3]
					if start < ter:
						reverse_comp = "no"
					else:
						reverse_comp = "yes"
				except IndexError:
					print("Index Error: Check that all columns are present in coordinates file.")
					print("File name: ", coords_file)
					print()
					exit()
					
				if error == "OK":
					coordinates_list.append([start, end, reverse_comp, name, sys, strand, func])
				else:
					continue
	return coordinates_list, coordFile

################################################################################################################################
def helper_RC(sequence):
	'''Returns the reverse-complement of the sequence. In string out string.'''

	comp = ""
	rc = ""
	sequence_len = len(sequence) +1

	# Complement the sequence
	for j in sequence:
		if j == "A" or j == "a":
			comp = comp + "T"
		elif j == "C" or j == "c":
			comp = comp + "G"
		elif j == "T" or j == "t":
			comp = comp + "A"
		elif j == "G" or j == "g":
			comp = comp + "C"
		else:
			comp = comp + "A"	# Needed for rare "N", "R", "K", etc. laying around. Can't just delete them for
								# without messing up subsequent functions so script just adds in an "A". For most
								# situations this shouldn't (meaningfully) affect anything.
			print("Error during RC process: Bad Base -- ", j)
			print('Replacing bad base with "A"')
			# print("Found in sequence: ", sequence)
			continue

	for k in range(1,sequence_len):
		l = k * -1
		rc = rc + comp[l]

	return rc
#################################################################################################
def compress(ntstring, strand):
	'''
	Input: whole gene nt sequence
	Output: string comprised only of the 1st, 2nd, or 3rd codon position nucleotides 
	'''
	outseq3 = ""
	outseq2 = ""
	outseq1 = ""

	
	if strand == "+":
		codon_position1 = 1
		codon_position2 = 2
		codon_position3 = 0
	elif strand == "-":
		codon_position1 = 0
		codon_position2 = 2
		codon_position3 = 1
	
	ctr = 0
	for n in ntstring:
		ctr += 1
		current_spot = ctr %3
		try:
			if current_spot == codon_position3:
				outseq3 = outseq3 + n
			elif current_spot == codon_position2:
				outseq2 = outseq2 + n
			elif current_spot == codon_position1:
				outseq1 = outseq1 + n
		except UnboundLocalError:
			print("# Error. Exiting.")
			exit()

	return outseq1, outseq2, outseq3
#################################################################################################
def analyze_sequence(crd, gnm_data):
	'''
	- Input is a set of coordinates
	- Output is nt sequence for that set of coordinates
	'''
	START = crd[0]-1
	STOP = crd[1]
	sequence_string = gnm_data[START:STOP]					# gives TOP STRAND Seq at specified coordinates
	pastTer = crd[2]										# Gene is located past ter: yes or no
	
	c1, c2, c3 = compress(sequence_string, crd[5])			# Reduces sequence to only codon positions 1, 2, or 3
	all_sequences = [sequence_string, c1, c2, c3]
	out_data = []

	for entry in all_sequences:
		nts = [0,0,0,0, pastTer]
		for nt in entry:
			if nt == "A" or nt == "a":
				nts[0] += 1
			elif nt == "C" or nt == "c":
				nts[1] += 1
			elif nt == "G" or nt == "g":
				nts[2] += 1
			else: # number of Ts
				nts[3] += 1
		out_data.append(nts)
	return out_data

#################################################################################################
def GCSanalysis(cds, genomedata, original_lines, genome_nm, originalDataLines, dataFolder):
	'''
	- Main analysis function, called by main overall function analyzeGCskew
	- Input: 	0) [[startGene1, endGene1, reverse_compGene1, nameGene1, sysGene1, strandGene1, funcGene1], [...]]
				1) full genome fasta sequence
				2) all coordinate file lines as list
				3) genome name as string, no file suffix
	- Outputs:
				Each is a list of lists: [[startGene1, stopGene1, GCskewGene1], [startGene2...
				1) Whole gene GC skew
				2) Codon position 1 GC skew
				3) Codon position 2 GC skew
				4) Codon position 3 GC skew
	'''
	with open(dataFolder + "/" + genome_nm + "_GCS123.txt", "w") as outfile:	# New coordinates format file with GC skew column
		outfile.write(original_lines[0].strip("\n") + "\tWholeGene_GCS\tCP1_GCS\tCP2_GCS\tCP3_GCS\n" )

		forwig = []
		forwigCP1 = []
		forwigCP2 = []
		forwigCP3 = []

		HO = [0,0,0,0]	# A, C, G, T number overall
		CD = [0,0,0,0]

		place_counter = 0
		for entry in cds:
			place_counter += 1									# Coordinates line number index
			new_coord_line = originalDataLines[place_counter].strip("\n")

			# Calculate GC skew
			# Data is #A,#C,#G,#T in the top strand for each gene and if it's before/after Ter
			analysis_out = analyze_sequence(entry, genomedata)	# Number A,C,G,T in each gene [[A,C,G,T, PastTer], [...
			switch = 0

			# For every gene...
			# Genes after ter need to have their GC skew reversed (value * -1)
			# Analyze out returns a list of lists. 1 list per gene, with the number of A,C,G,T, and pastTer for each
			# of the various GC skew calculations (whole gene, CP1, CP2, CP3). 
			# analysis_out = [[#A_gene1_whole, #C_gene1_whole, #G_gene1_whole, #T_gene1_whole, pastTer], [#A_gene1_CP1, #C_gene1_CP1...]]
			
			for v in analysis_out:
				switch += 1
				if switch == 1:
					skew = get_GCskew(v)
					# Record data
					forwig.append([entry[0], entry[1], skew])
				elif switch == 2:
					skew = get_GCskew(v)
					# Record data
					forwigCP1.append([entry[0], entry[1], skew])
				elif switch == 3:
					skew = get_GCskew(v)
					# Record data
					forwigCP2.append([entry[0], entry[1], skew])
				elif switch == 4:
					skew = get_GCskew(v)
					# Record data
					forwigCP3.append([entry[0], entry[1], skew])

				# Make the new outfile line by appending each GC skew to the line
				new_coord_line = new_coord_line + "\t"+ str(skew)
			
			# Write the new data line
			outfile.write(new_coord_line + "\n")

	return forwig, forwigCP1, forwigCP2, forwigCP3

#################################################################################################
def get_GCskew(inlist):
	'''
	In: List [#A,#C,#G,#T, past ter (Y/N)]
	Out: GC skew value
	'''
	ntA = inlist[0] 
	ntC = inlist[1]
	ntG = inlist[2]
	ntT = inlist[3]
	past_ter = inlist[4]
	try:
		s = (ntG-ntC)/float(ntG+ntC)
	except ZeroDivisionError:
		s = 0.1							# Rare, but occasional error yields a non-inverter data point
										# This is a qualitative determination, based on typical quatitative cutoff.
	if past_ter == "yes":				# Reverses all GC skew values after Ter
		s = s*-1

	return s

#################################################################################################
def calcGCS(g, c):
	'''
	- Calculates the GC skew
	- Input is the number of Gs and Cs in the string
	- Output is the GC skew, or a specified value in rare cases of zero division error
	'''
	try:
		return (g-c)/float(g+c)
	except ZeroDivisionError:
		if g > 1:
			return 0.1					# Arbitrary value due to no Cs but >1 G
		else:
			return 0

#################################################################################################
def load_fasta(target, inFolder):
	"""
	- Input is fasta
	- Output is string of full DNA sequence
	"""
	
	seq = ""
	line_num = 0

	with open(inFolder + "/" + target, "r") as fastaFile:
		for line in fastaFile:
			line_num += 1

			#ignore first line, get name though
			if line_num == 1:
				sys_name1 = line.strip(">")
				sys_name = sys_name1.strip("\n")
			else:
				seq = seq + line.strip("\n")

	return seq, sys_name

#################################################################################################################
def GC_skew_slide(seq, window): 
	"""
	- Calculates GC skew (G-C)/(G+C) over the specified sliding window length
	- Returns a list of ratios (floats), controlled by the length of the sequence and the size of the window. 
	- Does NOT look at any ambiguous nucleotides.
	- Returns list of tuples: (nt position which is the middle nt in the window, GC skew)
	""" 
	values = [] 
	s = seq[:window]
	number_g = s.count('G') + s.count('g') 
	number_c = s.count('C') + s.count('c') 
	window_middle = int(window/2)
	gcs = calcGCS(number_g, number_c)
	values.append((window_middle, gcs))	
	window_back = 1
	window_front = window
	window_middle = int(window/2)

	for i in range(len(seq)-window):
		window_middle += 1
		front_nt = seq[window_front + i]
		if front_nt == "g" or front_nt == "G":
			number_g += 1
		elif front_nt == "c" or front_nt == "C":
			number_c +=1
		
		back_nt = seq[window_back +i]
		if back_nt == "g" or back_nt == "G":
			number_g -= 1
		elif back_nt == "c" or back_nt == "C":
			number_c -=1
		values.append((window_middle, calcGCS(number_g, number_c)))
	return values

#################################################################################################
def create_wig(genome_len, gcdata, genomename, outfilename, gsn, dataFolder):
	"""
	- Inputs are genome length, GCskew data, genome name, outfile name, genome systematic name
		- GC skew data are [start, stop, GCSvalue]
		- For each nt in the range of the gene, the GCS value is recorded to the wiggle file
	- Outputs a new file with all nucleotide positions annotated.
	"""

	with open(dataFolder + "/" + outfilename + ".wig", "w") as out:
		out.write("track\nvariableStep chrom=" + gsn + " span=1\n")
		wigdata = [] # list of tuples
		
		for i in range(1,genome_len):
			wigdata.append([i,0])

		for j in gcdata:	
			s = j[0]
			t = j[1]
			v = j[2]
			for k in range(s-1,t):
				try:
					wigdata[k][1] = v
				except IndexError:
					continue

		for l in wigdata:
			out.write(str(l[0]) + "\t" + str(l[1]) + "\n")
	return None

#################################################################################################################
def print_wig_from_slide(vals, fname, ws, dataFolder, systematicName):
	"""
	- Writes sliding window GC skew data to an outfile
	- Input is GC skew values [nucleotide position, GC skew] 
	- Output is a wiggle file 
	"""
	with open(dataFolder + "/" + fname + "_GCSslide_" + str(ws) + ".wig", "w") as gcsk_out:
		gcsk_out.write("track\nvariableStep chrom=" + systematicName + " span=1\n")

		for entry in vals:
			gcsk_out.write(str(entry[0]) + "\t" + str(entry[1]) + "\n")
	return None

#################################################################################################################
def fileChecker(gnmNames, dataFldr):
	'''
	- Input is a list of genome names, folder where data are kept
	- Prints "Files are ok" or "Files not found" and exits if a coordinate and fasta files isn't found for
	  each of the specified genomes in the data directory.
	- Output is None
	'''
	missing = []
	missingFiles = "no"
	if not os.path.isdir("./" + dataFldr):
		print()
		print('Error: Data file folder "' + dataFldr + '" does not exist. Exiting')
		print()
		exit()
	else:
		files = os.listdir(dataFldr)
		for crdFileTuple in gnmNames:
			genomeName = crdFileTuple[0]
			coord_file = genomeName + "_coords.txt"
			fastagenome = genomeName + ".fasta"
			if coord_file in files:
				if fastagenome in files:
					continue
				else:
					missingFiles = "yes"
					missing.append(fastagenome)
			else:
				missingFiles = "yes"
				missing.append(coord_file)

		if missingFiles == "yes":
			print("#")
			print("# Error: The following are missing:")
			for x in missing:
				print("#     " + x)
			print("# Add the missing files and re-start the program.")
			print("# Exiting.")
			print("")
			exit()

#################################################################################################
def initialize(header):
	'''
	- Starts the program
	- Returns either additional information or starts the analysis
	'''
	
	# Print the program header
	if header == "y" or header == "Y":
		print()
		print()
		print(70*"#")
		print(70*"#")
		print("#")
		print("# GC Skew Analysis Tool v.5")
		print("# ")
		print("# Written by C. Breuer, Updated: 9/2021")
		print("#")
		print("# Reference: Gene inversion potentiates bacterial evolvability and virulence")
		print("#            Nature Communications 2018, PMID: 30405125")
		print("#")
		print(70*"#")
		print(70*"#")
		print()

	print()
	print("     Minimum Input Requirements:")
	print()
	print("     1. A folder named 'Data' containing:")
	print('     2. The file "ter_coords.txt" (See help menu.)')
	print()
	print("        For each genome:")
	print("     3. <genome_abbreviation>.fasta")
	print("     4. <genome_abbreviation>_coords.txt")
	print()
	print()
	answer = input("     Continue? (Y/N/Help): ")

	if answer == "Y" or answer == "y":
		print()
		print(70*"#")

		# Run the program
		return analyzeGCskew()
	
	# Print help menu, then re-initialize the program without the header
	elif answer == "h" or answer == "H" or answer == "Help" or answer == "help":
		print()
		print()
		print(70*"#")
		print("# GC Skew Analysis Help:")
		print("#")
		print("# Script Function:")
		print("#      - The script will perform multiple GC Skew analyses on all genome files in the specified")
		print("#        directory. (Default is the working directory.")
		print("#")
		print("#      - Analyses include calculation of the GC Skew over a sliding window, and over whole gene regions.")
		print("#")
		print("#      - For whole gene regions, the GC skew will be calculated using:")
		print("#              - The full gene sequence")
		print("#              - Only nucleotides in codon position 1")
		print("#              - Only nucleotides in codon position 2")
		print("#              - Only nucleotides in codon position 3")
		print("#")
		print("#       - Note that the gene-based GC skew values will be calculated relative to the leading arm of")
		print("#                the chromosome, meaning that all values after ter will be reversed (x*-1).")
		print("#")
		print("# Target files:")
		print("#     - Target genome files are specified by the names in the 'ter_coords.txt' file which")
		print("#       must be a tab-delimited text file with each line formatted as:")
		print("#       <genome_name><tab><ter_coordinate>")
		print("#")
		print("#     - For every genome in the first column, the script will look for an associated .fasta and")
		print("#       coordinates file in the target directory.")
		print("#")
		print("#     - File naming conventions are key.")
		print("#")
		print("#     - Also note that the systematic name of each genome, recorded in the first line of the .fasta")
		print("#       file will be written into the .wig files. This name must match the genome name in the")
		print("#       visualization software (e.g. MochiView) or you may run into issues looking at the data.")
		print(70*"#")
		print()
		return initialize("no")
	else:
		print()
		print("Ok. Exiting.")
		print()
		exit()
#################################################################################################################
def analyzeGCskew():
	'''
	- Main function
	- User can customize outfolder names, ter_coords file name, etc.
	- Implements GC skew analysis, largely by calling the analysis function
	'''
	
	# Hyperparameters (customize as needed):
	inputFolder = "Data"
	outfolder = "GCS_AnalysisOut"
	genomesFile = "ter_coords.txt"
	sliding_window_size = 200

	# Check on the existencen of a folder called "data" and a make a new folder if needed
	if not os.path.isdir("./" + outfolder):
		os.mkdir(outfolder)

	# Input is ter_coords.txt file. Change this is you change the file naming convention
	# Yields [(genome1, ter1), (genome2...
	all_genomes = getGenomesAndTers(genomesFile, inputFolder)

	# Makes sure all files are there
	fileChecker(all_genomes, inputFolder)
	print("#")

	print("# The following genomes will be analyzed:")
	for g in all_genomes:
		print("# " + g[0])
	print("#")
	print(70*"#")

	# Iterates through all genomes in ter coordinates file	
	for genome in all_genomes:													# Process each genome
		print("#")
		print("# Currently processing: ", genome[0])
		#################################################################################################
		# Collect all genome information in memory
		Gname = genome[0]														# Genome name, and file prefix for files
		ter_position = int(genome[1])											# Ter coordinate
		coordinates_file = Gname + "_coords.txt"								# Coords view file name
		fastagenome = Gname + ".fasta"											# Fasta file name										
		print("# Loading fasta file...")
		genomesequence, gnm_sys_name = load_fasta(fastagenome, inputFolder)					# This genome's sequence and genome systematic name
		genome_length = len(genomesequence) -1									# This genome's length
		#################################################################################################
		# Print GCS whole genome, Sliding window
		print("# Calculating GC skew...")
		sliding_data = GC_skew_slide(genomesequence, sliding_window_size)
		print("# Writing wiggle file #1 (sliding window for full genome)...")
		print_wig_from_slide(sliding_data, Gname, sliding_window_size, outfolder, gnm_sys_name)
		#################################################################################################
		# Calculate the GC skew for individual genes, and write the data to a new coordinates file
		
		# Collect gene coordinates from file
		print("# Collecting gene coordinates for ", genome[0])
		cds1, coord_data = get_coordinates(coordinates_file, ter_position, inputFolder)
		
		# Analysis function, executes other functions
		All_data, CP1_data, CP2_data, CP3_data = GCSanalysis(cds1, genomesequence, 						
			coord_data, Gname, coord_data, outfolder) 															# Creates new coordinates file copy with additional column for GC skew
		
		print("# Writing wiggle files showing the by-gene GC skew value using various calculations:")
		print("#      Wiggle file #2 (whole-gene)")
		create_wig(genome_length, All_data, Gname, (Gname +"_GCS_genes"), gnm_sys_name, outfolder)							# Creates By-Gene GC skew wiggle file
		print("#      Wiggle file #3 (CP1)")
		create_wig(genome_length, CP1_data, Gname, (Gname +"_GCS_CP1_genes"), gnm_sys_name, outfolder)						# Creates By-Gene GC skew wiggle file
		print("#      Wiggle file #4 (CP2)")
		create_wig(genome_length, CP2_data, Gname, (Gname +"_GCS_CP2_genes"), gnm_sys_name, outfolder)						# Creates By-Gene GC skew wiggle file
		print("#      Wiggle file #5 (CP3)")
		create_wig(genome_length, CP3_data, Gname, (Gname +"_GCS_CP3_genes"), gnm_sys_name, outfolder)						# Creates By-Gene GC skew wiggle file
	print("#")
	print(70*"#")
	print("#")
	print("# All processing has completed.")
	print("#")
	print(70*"#")
	print()
	exit()
#################################################################################################
# Function calls 
#################################################################################################
initialize("Y")
# Script design:
# User customization of hyperparameters should happen within the analyzeGCskew() function.
# Function initialize() calls analyzeGCskew() which calls:
#    1. File checker
#    2. load_fasta
#    3. GC_skew_slide - calculates the sliding window GC skew profile for the whole genome
#    4. print_wig_from_slide - writes the sliding window skew data to a .wig file 
#	 5. get_coordinates() - loads the coordinate information for all features
#    6. GCSanalysis():
#             - calculates the GCS for all features (genes) using all methods (CP1/2/3/all)
#             - records a copy of the original feature file data with additional columns for GCS data
#	 7. create_wig - called for each method, cretes a new by-gene wiggle file
