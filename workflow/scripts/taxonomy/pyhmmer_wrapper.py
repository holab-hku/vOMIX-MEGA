import os
import sys
import psutil
import argparse
import time

import pyhmmer
from pyhmmer.easel import SequenceFile
from pyhmmer.plan7 import HMMFile
from pyhmmer.hmmer import hmmsearch
from pyhmmer.hmmer import hmmscan




def run_pyhmmer(proteins, hmmdb, scan, tblout, domtblout, bitcutoff, cores_n):
	if cores_n != 0:
		cpu_count = cores_n
	else:
		cpu_count = psutil.cpu_count()
	cpu_total = psutil.cpu_count()
	print("Total number of CPUs:", cpu_total)
	print("Number of CPUs in use:", cpu_count)

	available_memory = psutil.virtual_memory().available
	proteins_size = os.stat(proteins).st_size
	database_size = os.stat(hmmdb).st_size
	input_size = proteins_size + database_size
	
	# There are 2 different strategies for memory allocation 
	# 1) If input and database only take less than 10% of the total memory of your local machine
	# 2) If input only take less that a 10% of the total memory of the local machine
	# 3) If input takes more that a 10% of the total memory of the local machine

	with HMMFile(hmmdb) as hmm_file:
		hmms = list(hmm_file)
		n_hmms = len(hmms)	

		with SequenceFile(proteins, digital=True) as seq_file:
			if input_size < available_memory * 0.1:
				print("\nEnough available memory!")
				print("Pre-fetching targets into memory...\n")
				print(f"Total Input size: {input_size/(1024**3):.2f}G")
				print(f"Total available memory: {available_memory/(1024**3):.2f}G")
				print(f"Percentage memory used: {input_size/available_memory:.2f} %")
				seqs = seq_file.read_block()
			else:
				seqs = seq_file

			print("\nPerforming pyhmmer hmmsearch...")
			t1 = time.time()
			all_hits_list = list(hmmsearch(hmms, seqs, cpus=cores_n, Z=n_hmms, bit_cutoffs=bitcutoff))
			
			time_in_seconds = time.time() - t1
			hours = time_in_seconds // 3600
			minutes = (time_in_seconds % 3600) // 60
			seconds = (time_in_seconds % 3600) % 60

			print(f"Hmmsearch of {n_hmms} HMMs on {cpu_count} CPUs took {hours} hours, {minutes} minutes, {seconds:.0} seconds")
				
			print("Writing results into output file...")
			t1 = time.time()
			with open(domtblout, "wb") as f:
				all_hits_list[0].write(f, format="domains", header=True)
			with open(domtblout, "ab") as f:
				for hits in all_hits_list[1:]:
					hits.write(f, format="domains", header=False)
			print(f"Writing the file took {time.time() - t1:.3} seconds")

	
			if tblout is not None:
				with open(tblout, "wb") as f:
					all_hits.write(f, format="targets", header=True)

		#	if domtblout is not None:
		#		with open(domtblout, "wb") as f:
		#			all_hits.write(f, format="domains", header=True)

	
	


if __name__ == "__main__":
	parser = argparse.ArgumentParser()

	parser.description = """A wrapper script of pyhmmer which takes into account memory of the system and either parses faa input into memory and pre-fetches, greatly improving speed."""
					
	parser.add_argument(
		"-i",
		"--proteins",
		dest="proteins_faa",
		help="path to input protein fasta file", 
		required=True)

	parser.add_argument(
                "-db",
                "--hmmdb",
                dest="hmm_db",
                help="path to HMM database with .hmm ending", 
		required=True)

	parser.add_argument(
                "-t",
                "--tblout",
                dest="tbl_out",
                help="saves the output file of domains equivelant to the --tblout flag in hmmer.",
		default=None)

	parser.add_argument(
		"-d", 
                "--domtblout",
                dest="domtbl_out",
                help="saves the output file of domains equivelant to the --domtblout flag in hmmer.", 
		default=None)
	
	parser.add_argument(
		"-s",
		"--hmmscan",
		action="store_true", 
		default=False, 
                dest="hmm_scan",
                help="if the flag is on, the program will run hmmscan instead of hmmsearch[default] and swap Z scores")

	parser.add_argument(
                "-b",
                "--bit_cutoff",
                dest="bit_cutoff",
                help="choose bit cutoff as ['gathering', 'noise', 'trusted']. None by default",
                default=None)

	parser.add_argument(
                "-c",
                "--cores",
                dest="n_cores",
		type=int, 
                help="Number of cores to use. By defualt is 0 which means pyhmmer chooses the optimal number.",
                default=0)

	
	args = parser.parse_args()
	
	if args.tbl_out is None and args.domtbl_out is None:
		print("At least one type of output must be selected.")
		print("Set the output file path with --tblout or --domtblout.")
		print("Run --help or -h to find out more.")
		sys.exit(1)

	run_pyhmmer(args.proteins_faa, args.hmm_db, args.hmm_scan, 
		args.tbl_out, args.domtbl_out, args.bit_cutoff, args.n_cores)
