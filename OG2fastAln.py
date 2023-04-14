#!/usr/bin/python3
import shutil
import sys
import os
import fnmatch
import fileinput
import re
import sys
import os
from os import listdir
from os.path import isfile, join
from Bio import SeqIO
from Bio.SeqIO import FastaIO
import string
from Bio import AlignIO
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Align.Applications import MafftCommandline
import argparse
import itertools
import progressbar
import time
from alive_progress import alive_bar
from time import sleep


###################################################################################################
#### Script from Guillaume Cossard:
### v 1.0
### 05/03/2023

###################################################################################################
# Format of the ctl file output if --ctl option -- Modify to adapt to the CodeML files needed
###################################################################################################
CTL_text = '''seqfile = {0}.phy    ### sequence data file name
treefile = Tree_4sp_Gracilaria_Unrooted.txt ### tree structure file name

outfile = {1}.mlc          ### main result file name
noisy = 0   ### 0,1,2,3,9: how much rubbish on the screen
verbose = 0   ### 1: detailed output, 0: concise output
runmode = -2   ### 0: user tree;  1: semi-automatic;  2: automatic
		### 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

seqtype = 1   ### 1:codons; 2:AAs; 3:codons-->AAs
CodonFreq = 2   ### 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
clock = 0   ### 0: no clock, unrooted tree, 1: clock, rooted tree
aaDist = 0   ### 0:equal, +:geometric; -:linear, (1-5:G1974,Miyata,c,p,v)
model = 0

NSsites = 0 1 2
		### 0:one w; 1:NearlyNeutral; 2:PositiveSelection; 3:discrete;
		### 4:freqs; 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;10:3normal
icode = 0   ### 0:standard genetic code; 1:mammalian mt; 2-10:see below
Mgene = 0   ### 0:rates, 1:separate; 2:pi, 3:kappa, 4:all

fix_kappa = 0   ### 1: kappa fixed, 0: kappa to be estimated
kappa = 2   ### initial or fixed kappa
fix_omega = 0   ### 1: omega or omega_1 fixed, 0: estimate
omega = 0.5  ### initial or fixed omega, for codons or codon-based AAs
ncatG = 10   ### # of categories in the dG or AdG models of rates

getSE = 0   ### 0: dont want them, 1: want S.E.s of estimates
RateAncestor = 0   ### (0,1,2): rates (alpha>0) or ancestral states (1 or 2)

Small_Diff = .45e-6
cleandata = 1  ### remove sites with ambiguity data (1:yes, 0:no)?
fix_blength = 0  ### 0: ignore, -1: random, 1: initial, 2: fixed
'''


def parseArgs():
	parser = argparse.ArgumentParser(
		description='Generate phylip alignment (.phy) trimmed with Gblocks (as well as .ctl option file for codeml if asked) \
		from the list of Single copy Orthogroups (output of OrthoFinder) \
		Order of genome in SampleFile has to be the same as in the Orthogroups header \
		Running time ~ 65min for 12,000 genes \
		Requires : 	translatorx_vLocal.pl (runs with perlbrew) ; MAFFT ; Fasta2Phylip.py ; GBlocks (And check python modules) ; RmEmptyseqFromFasta.sh\n\n\n Gblocks directory has to be in the PATH to permit Gblocks command to be found ! --> No cleaned alignements (_clean_ali.fasta) are gonna be output from TranslatorX otherwise')
	parser.add_argument('-SC', '--SingleCopyOrthologs', required=True,
                     type=str, help='list of all single copy orthogroups, one column file. Can include OG not shared among all species')
	parser.add_argument('-OG', '--orthogroups', required=True,
                     type=str, help='Orthogroups.csv (csv-format from OrthoFinder) - Table listing genes of each species for each orthogroup. Could be csv or tsv file')
	parser.add_argument('-g', '--genomes', required=True, type=str,
                     help="txt file that indicate species, genome path and phylip abreviation to apply (HAS TO BE In THE SAME ORDER AS IN ORTHOGROUPS file). With header")
	parser.add_argument('-aa', '--AApath', required=False, type=str,
                     help="Full path to the folder containing Orthogroups sequences in aa (from Orthofinder)")
	parser.add_argument('-c', '--ctl', required=False,
                     help="use if ctl file have to be provided (CodeML input)", action="store_true")  # The keyword “action” is being given the value “store_true” which means that if the option is specifed, then assign the value “True”
	return parser


def main(args):
	args, leftovers = parseArgs().parse_known_args(args)
	list_Ortho = args.SingleCopyOrthologs
	Groups = args.orthogroups
	# filter empty lines and removes header = nb species
	nb_sp = sum(1 for line in open(args.genomes, 'r') if line.rstrip())-1
	print(str(nb_sp))
	Samplefile = args.genomes
	if args.AApath is None:
		aaOGPAth = "./"
	else:
		aaOGPAth = args.AApath
	SCO_list = SCO2list(list_Ortho)
	SPlist = SP2list(nb_sp, Groups)
	# print(SPlist)
	SC_groupDict = Ortho2dict(nb_sp, Groups, SPlist, SCO_list)
	# print(SC_groupDict["OG0021638"])
	GenomeList = Genomes2dict(nb_sp, Samplefile)
	AbrevList = Abrev2List(Samplefile)
	dirname = 'Aln_Ortho'
# write SCorthologs fasta & per species fasta with OG names and Species gene names in "Aln_Ortho" directory"
	writeSC(nb_sp, SC_groupDict, GenomeList, AbrevList, SPlist, dirname)
	TransXAlign(dirname, aaOGPAth, AbrevList)
	if args.ctl is None:
		print("All alignments trimmed - end of script")
		pass
	else:
		PAMLctl()
		print("All alignments trimmed - CTL files ready for CodeML -- end of script")


def SP2list(nb_sp, Groups):
	with open(Groups, 'r') as grouping:
		Firstline = grouping.readline()
		List_Genelist = [list()] * nb_sp  # list of list of nb_sp elements
		# print(List_Genelist)
		SPlist = [str()] * nb_sp
		for i in range(nb_sp):  #  from 0 to 7
			SPlist[i] = Firstline.strip().split('\t')[i+1]
			print("\nSpecies {0} = {1}\n".format(str(i+1), str(SPlist[i])))
		return (SPlist)


def SCO2list(list_Ortho):
	SClist = []
	with open(list_Ortho, 'r') as SC:
		for line in SC:  #  from 0 to 7
			line = line.strip()
			if line not in SClist:
				SClist.append(line)
	return (SClist)


def Ortho2dict(nb_sp, Groups, SPlist, SCO_list):
	groupDict = {}
	SC_groupDict = {}
	with open(Groups, 'r') as grouping:
		Firstline = grouping.readline()
		List_Genelist = [list()] * nb_sp  # list of list of nb_sp elements
		for line in grouping:
			OGROUP = line.strip().split('\t')[0]
			#  from 1 to 8 with increment of 1 (because no need the first column, so i=0) ##replace 9 by nb_sp+1 ?
			for i in range(1, nb_sp+1, 1):
				# print(i)
				List_Genelist[i-1] = line.split('\t')[i].strip().split(", ")
				if OGROUP not in groupDict:
					groupDict[OGROUP] = {}
					groupDict[OGROUP][SPlist[i-1]] = List_Genelist[(i-1)]
				else:
					groupDict[OGROUP][SPlist[i-1]] = List_Genelist[(i-1)]
				if OGROUP not in SC_groupDict and OGROUP in SCO_list:
					# NbGenes = len(line.strip().split('\t'))
					SC_groupDict[OGROUP] = {}
					SC_groupDict[OGROUP][SPlist[i-1]] = List_Genelist[i-1]
				elif OGROUP in SC_groupDict and OGROUP in SCO_list:
					SC_groupDict[OGROUP][SPlist[i-1]] = List_Genelist[i-1]
				# print(List_Genelist[(i-1)])
	# return(groupDict)
	return (SC_groupDict)


def Genomes2dict(nb_sp, Samplefile):
	print(nb_sp)
	GenomeList = [dict()] * nb_sp
	FileList = []
	with open(Samplefile, 'r') as Samples:
		Samples.readline()
		for line in Samples:
			if line.rstrip():  #  if line is not empty
				print(line.strip().split('\t'))
				FileList.append(str(line.strip().split('\t')[1]))
		# print(len(FileList))
		for i in range(nb_sp):  #  from 0 to 7
			# print(i)
			SPdict = SeqIO.to_dict(SeqIO.parse(open(FileList[i], 'r'), 'fasta'))
			GenomeList[i] = SPdict
		return (GenomeList)


def Abrev2List(Samplefile):
	AbrevList = []
	with open(Samplefile, 'r') as Samples:
		Samples.readline()
		for line in Samples:
			if line.rstrip():  #  if line is not empty
				print(line.strip().split('\t'))
				AbrevList.append(str(line.strip().split('\t')[2]))
	return (AbrevList)


def writeSC(nb_sp, SC_groupDict, GenomeList, AbrevList, SPlist, dirname):
	dir = dirname
	try:
		os.mkdir(dir)
	except OSError:
		print('\nTry to create Folder "{0}" but it already exists.'.format(dir))
		REP = input('remove Folder :: YES[y] or NO[n] ?')
		if REP == "n":
			print("Change the name of the Folder to continue")
			dir = input('Enter a new folder name ? ')
			os.mkdir(dir)
		else:
			shutil.rmtree(dir)
			os.mkdir(dir)

	#  creates list of x empty element, x = nb of species
	List_outSEQ = [str()] * nb_sp
	List_outOGSEQ = [str()] * nb_sp
	ListGeneSeq = [str()] * nb_sp

	for i in range(nb_sp):  # 0 to 7
		List_outSEQ[i] = open('SCOrtho_GeneName_{0}.fasta'.format(AbrevList[i]), 'w')
		List_outOGSEQ[i] = open('SCOrtho_OGname_{0}.fasta'.format(AbrevList[i]), 'w')
	cnt = 0
	for k in SC_groupDict:
		Err = open('Errors_OrthoAln.txt', 'w')
		OUTORTHO = open("SC-Ortho_{0}.fasta".format(k), "w")
		# print(SC_groupDict[k])
		for i in range(nb_sp):  # 0 to 7
			ListGeneSeq[i] = str(SC_groupDict[k][SPlist[i]][0])
			if len(ListGeneSeq[i]) == 0:
				print("WARNING :: NO sequence of Orthogroup {0} for species : {1}".format(
					k, AbrevList[i]))
				Err.write(
					'NO sequence of Orthogroup {0} for species : {1}"\n'.format(k, AbrevList[i]))
				List_outSEQ[i].write('>{0}\n{1}\n'.format(str(k), ""))
				List_outOGSEQ[i].write('>{0}\n{1}\n'.format(str(k), ""))
			elif ListGeneSeq[i] not in GenomeList[i]:
				print("WARNING :: sequence  {0} of Orthogroup {1} not in genome fasta file of {2}".format(
					str(ListGeneSeq[i]), k, AbrevList[i]))
				Err.write('sequence  {0} of Orthogroup {1} not in genome fasta file of {2}\n'.format(
					str(ListGeneSeq[i]), k, AbrevList[i]))
				List_outSEQ[i].write('>{0}\n{1}\n'.format(str(k), ""))
				List_outOGSEQ[i].write('>{0}\n{1}\n'.format(str(k), ""))
			elif ListGeneSeq[i] in GenomeList[i]:
				List_outSEQ[i].write('>{0}\n{1}\n'.format(
					str(ListGeneSeq[i]), GenomeList[i][ListGeneSeq[i]].seq))
				List_outOGSEQ[i].write('>{0}\n{1}\n'.format(
					str(k), GenomeList[i][ListGeneSeq[i]].seq))
				# remove last 3 character of codon stop as translatorX does not read the (*) of Stop codons in AA fasta files
				OUTORTHO.write('>{0}\n{1}\n'.format(
					AbrevList[i], GenomeList[i][ListGeneSeq[i]].seq[:-3]))
		cnt += 1
		shutil.move("SC-Ortho_{0}.fasta".format(k), dir)

	for j in range(nb_sp):  # Create perSP fasta without empty sequences (necessary for RNAseq mapping with kallisto) with OGnames
		List_outOGSEQ[j].close()
		os.system('RmEmptyseqFromFasta.sh {0}' .format(
			"SCOrtho_OGname_{0}.fasta".format(AbrevList[j])))
	print(
		'\nNumber of Single-copy Orthologs written in alignment fasta = {0}\n'.format(str(cnt)))


def TransXAlign(dirname, aaOGPAth, AbrevList):
	dir = dirname
	os.chdir(dir)
	Path = os.getcwd()
	aaOGPATH = aaOGPAth
	print("\n\nCurrent directory is : {0}\n".format(Path))
	os.mkdir("Clean_nucl_Align")
	onlyfiles = [f for f in listdir(Path) if isfile(join(Path, f))]
	count = 0
	stop_codon_list = ["TAA", "TAG", "TGA"]
	StopcodonFile = open('Stop_codons_in_alignments.txt', 'w')
	NbseqwithStopCodons = 0
	with alive_bar(len(onlyfiles), bar='smooth', spinner='elements', title="Processed files") as bar:
		for file in onlyfiles:
			if file.startswith('SC-Ortho_') and file.endswith('.fasta'):
				aaOG = file.split('SC-Ortho_')[1].split('.fasta')[0]
				aaOGname = str(aaOG+'.fa')
				aaOGfile = join(aaOGPATH, aaOGname)
				# FirstAAdict=SeqIO.to_dict(SeqIO.parse(open(aaOGfile, 'r'), 'fasta'))
				# FirstAA = open(aaOGfile, 'r')
				cntsp = 0

				# Block part to update if wanting to input the ALIGNEMENT of aa directly into TranslatorX (used in Gblocks in .fasta format). Could be done with pasting in the block after deinterleave de fasta:
				# import subprocess
				# in_file = 'AA_{0}.fa'.format(aaOG)
				# subprocess.call(["mafft", "--out", "Aln_AA_{0}.fasta".format(aaOG), in_file])
				###

				# with open('temp_AA_{0}.fa'.format(aaOG), 'w') as AAout: ### rename sequences by abreviated species name from Sample file (correspondance nt - aa)
				# 	for line in FirstAA:
				# 		if line.startswith('>'):
				# 			AAout.write('>'+AbrevList[cntsp]+'\n')
				# 			cntsp += 1
				# 		else:
				# 			AAout.write(line)
				# 	FirstAA.close()
				# 	AAout.close()
				# aaAbvOGfile=str('temp_AA_{0}.fa'.format(aaOG))
				# ## code from deinterleave_fasta.py in Scripts
				# linenum = 0
				# infile = open(aaAbvOGfile, 'r')
				# with open('AA_{0}.fa'.format(aaOG), 'w') as outfile:
				# 	for line in infile:
				# 		line = line.strip()
				# 		if len(line) > 0:
				# 			if line[0] == '>':
				# 				if linenum == 0:
				# 					outfile.write(line + '\n')
				# 					linenum += 1
				# 				else:
				# 					outfile.write('\n' + line + '\n')
				# 			else:
				# 				outfile.write(line)
				# outfile.close()
				if os.stat("%s" % (file)).st_size != 0:
					# -p F  means MAFFT is used to align AA sequences.  Input aligned aa sequence (from mafft alignement exectued before, by adding -a Aln_AA_{0}.fasta) with at the end (.format(aaOG))
					os.system(
						'translatorx_vLocal.pl -i {0} -o TransX_{0} -t T -p M -g "-t=p" 1>>stdOUT_TransX.txt 2>>stderr_TransX.txt'.format(file))
					# this is the name of alignement files out of translatorX
					Newfile = "TransX_{0}.nt_cleanali.fasta".format(file)
					count += 1
					if isfile(Newfile):
						OpFile = open(Newfile, 'r')
						XX = SeqIO.parse(OpFile, 'fasta')
						for seq in XX:  # Check stop codons to record them in output text file as theywill need to get removed for phylip conversion
							if len(seq.seq) % 3 != 0:
								print("WARNING :: Sequence {0}, in file {1}, has a length that's not a multiple of 3 !".format(str(Newfile, seq.id)))
								continue
							else:
								n = len(seq.seq)
								k = 0
								# store the starting positions of the codons
								found_codon_positions = []
								while k < n-2:  # While loop to check for stop codons in the alignment, replace by Ns
									# extract a three-nucleotide subsequence
									possible_codon = seq.seq[k:k+3]
									# print(possible_codon)
									if possible_codon in stop_codon_list:
										found_codon_positions.append(k)
									k += 1
									NbseqwithStopCodons += 1
								# print('found stop codons at indices {0}'.format(found_codon_positions))
								StopcodonFile.write(str(Newfile) + ':\n' + str(seq.id) + ':\n' +
								                    'found stop codons at indices {0}'.format(found_codon_positions) + '\n\n\n')
						OpFile.close()
					shutil.move(Newfile, "Clean_nucl_Align")
					# print("Done {0}/{1} alignments ... \n".format(str(count), str(len(onlyfiles))))
					# bar2.update(count)
					bar()
				else:
					pass
			for f in listdir(Path):
				# and not f.endswith('aa_based_codon_coloured.html')
				if f.startswith("TransX") and isfile(f) and not f.endswith('nt_cleanali.fasta') and not f.endswith('_TransX.txt'):
					os.remove(f)


def PAMLctl():
	print('\nProceeding to "PAML-ready" phylip conversion and codeml configuration files (.ctl)')
	os.chdir("Clean_nucl_Align")
	Path = os.getcwd()
	onlyfiles = [f for f in listdir(Path) if isfile(join(Path, f))]
	totfiles = len(onlyfiles)
	numfiles = 0
	with alive_bar(len(onlyfiles), bar='filling', spinner="dots_waves", title="Processed files") as bar:
		for file in onlyfiles:
			if file.endswith(".fasta"):
				numfiles += 1
				os.system("Fasta2NOSTOPPhylip.py {0} {1}.phy".format(file, file))
				bar()
				# print("{}/{} files converted to Phylip with replaced stop codons".format(numfiles, totfiles))
	count = 0
	onlyfiles = [f for f in listdir(Path) if isfile(join(Path, f))]
	# print(onlyfiles)
	# print(listdir(Path))
	with alive_bar(len(onlyfiles), bar='filling', spinner="elements", title="Processed files") as bar3:
		for file in onlyfiles:
			if file.endswith(".phy"):
				# os.system("phylip2paml.pl {0}".format(file))
				ctlname = file.split(".phy")[0]
				Outclt = open("{0}.ctl".format(ctlname), "w")
				Outclt.write(CTL_text.format(ctlname, ctlname))
				count += 1
				# print("Processed {0}/{1} files ...".format(str(count), totfiles))
			bar3()

if __name__ == "__main__":
	args = None
	if len(sys.argv) == 1:
		args = ["--help"]
	main(args)