#!/usr/bin/env python

import os 
import sys 
import subprocess
import argparse
import json
import time
import re
import warnings
from io import StringIO
import datetime

import pandas as pd
from Bio import Entrez

from rich.console import Console
from rich.progress import Progress
from rich.layout import Layout
from rich.panel import Panel
console = Console()


def validate_samples(samples):
	"""
    	Performs checks on the samples provided in the sample_list.tsv
    	to either verify their existence locally or on the SRA.

   	 If you are using local files, you can either:
   	 1) Provide full file paths in the R1 and R2 columns of sample_list.tsv [only R1 if single-end]
	 2) Place the fastq files in the config['datadir'] path with <sample>_{1,2}.fastq.gz naming format.
      	 If a file has no paired {2} label, it is assumed to be single-end.

    	"""
	console.print(Panel.fit("""[dim italic]Validating availability of local and remote SRA raw sequences.\n\nIf you are using local files, you can either:\n1) Provide full file paths in the R1 and R2 columns of sample_list.tsv [only R1 if single-end]\n2) Place the fastq files in the config['datadir'] path with <sample>_{1,2}.fastq.gz naming format.\nCo-assemblies and mix-assemblies can be setup by writing the same <assembly> column for different samples.\n""", title="Sample Validation", subtitle="In Progress..."))
	with Progress(transient=True) as progress:
		task = progress.add_task("[cyan]Validating...", total=len(samples))

		for sample, items in samples.items():

			# check if files exist locally
			R1path = samples[sample]["R1"]
			R2path = samples[sample]["R2"]
			if os.path.exists(R1path) and os.path.exists(R2path):
				R1size = os.path.getsize(R1path) / (1024 ** 3)
				R2size = os.path.getsize(R2path) / (1024 ** 3)

				sizegb = round((R1size + R2size), 2)
				console.print("[dim] {accs} pre-downloaded [{size_gb}GB][/dim]".format(accs=sample, size_gb=sizegb))
				progress.update(task, advance=1)
				time.sleep(0.5)
				continue

			# check if it exists in SRA if not present locally
			acc = items['accession']

			try:
				handle = Entrez.efetch(db="sra", id=acc, retmax=1000, rettype="full", retmode="xml")
				record = handle.read()
				if not isinstance(record, str):
					record = record.decode("utf-8")
				match = re.search(r'size="(\d+)"', record)
				sizebyte = match.group(1)
				sizegb = round(int(sizebyte) /pow(1024, 3), 2)
				console.print("[dim] {accs} downloading [{size_gb}GB][/dim]".format(accs=acc, size_gb=sizegb))
				
			except Exception as e:
				console.print(Panel.fit(f"Sample '{sample}' [dim]with accession '{acc}' could not be processed for the following reasons:\n\n1) The Accession Could noound by efetch in NCBI's Entrez Direct\n2) The accession is not a valid run\n3) The NCBI server is limiting the number of requests\n4) The file is not in correct fasta gunzipped format with a 'fasta.gz' ending.\n\nPlease ensure that it is a valid SRA accession or reduce the number of requests for sample validation. Alternetaively you could change config['email'] parameter.If you intend to use locally stored fastq files, make sure your sample list contains the column 'R1' for single-end and 'R2' for paired-end files.", title="Run Error", subtitle="SRA Accession Not Parsed"))
				sys.exit(1)

			progress.update(task, advance=1)
			time.sleep(0.5)

		progress.stop_task(task)


def parse_sample_list(f, datadir, outdir, email, time):
	"""
	Parse the sample list. Each sample is stored as a dictionary in the samples{} dictionary.
	samples{sample_name} will have the following information:
	
	samples[sample_name] = {'R1': 'path to R1',
				'R2': 'path to R2',
				'accession': 'accession id'}
	"""

	if not os.path.exists(datadir):
		os.makedirs(datadir)
	if not datadir.endswith(os.sep):
		datadir = os.path.join(datadir, '')
	
	samples = {}
	
	###################
	# PARSE SAMPLE DF #
	###################

	df = pd.read_csv(f, comment='#', header=0, sep='\t', index_col=False, dtype=str)
	df['sample_id'] = df['sample_id'].fillna(df['accession'])
	df = df.map(lambda x: x.strip() if isinstance(x, str) else x) # strip white space
	df['sample_id'] = df['sample_id'].astype(str)
	df = df.replace(r'^\s*$', float('nan'), regex=True)
	if 'assembly' not in df.columns:
		    df['assembly'] = float('nan')
	
	# iterate through df and if sample_id is missing, replace is with SRA accession
	# also add assembly as sample_id if no assembly is given (single-sample assembly)
	for i, row in df.iterrows():
		if pd.isna(row['sample_id']) and not pd.isna(row['accession']):
			df.loc[i, 'sample_id'] = row['accession']

		elif pd.isna(row['sample_id']) and pd.isna(row['accession']):
			raise ValueError("Column with both empty sample_id and SRA accession, please provide at least one")
	
	# if there are no assemblies set, make it the same as sample_id (single sample assembly)
	df['assembly'] = [row['assembly'] if not pd.isna(row['assembly']) else row['sample_id'] for _, row in df.iterrows()]
	# if there are no R1 or R2s, generate them
	df['R1'] = [row['R1'] if not pd.isna(row['R1']) else '{dir}{s}_1.fastq.gz'.format(dir= datadir, s=row['sample_id']) for _, row in df.iterrows()]
	df['R2'] = [row['R2'] if not pd.isna(row['R2']) else '{dir}{s}_2.fastq.gz'.format(dir = datadir, s=row['sample_id']) for _, row in df.iterrows()]

	# FUTURE WARNINGN: Setting an item of incompatible dtype is deprecated and will raise an error in a future version of pandas. Value '' has dtype incompatible with float64, please explicitly cast to a compatible dtype first.
	warnings.simplefilter(action='ignore', category=FutureWarning)
	df.fillna('', inplace=True)
	# set unique names for the file index
	df.set_index('sample_id', inplace=True)
	
	df.index[df.index.duplicated()]	
	if df.index.duplicated().any():
		console.print(Panel.fit("ValueError on df.set_index('sample_id'). Values in the sample_id column may be non-unique. Please check {}".format(f), title="Value Error", subtitle="Duplicate 'sample_id' Names"))
		raise ValueError()

	# If it's just a one-column file, expand it by assuming that the
	# SRA accessions are in the first column
	if (df.shape[1] == 0):
		if ('accession' not in df.columns):
			df = df.assign(accession = pd.Series(df.index, index = df.index))
			df.index.name = "sample_id"
	
	# Remove duplicates rows and throw a warning
	duplicaterow = df[df.duplicated()]
	duplicateid = df[df.index.duplicated()]
	

	if not (duplicaterow.empty and duplicateid.empty):

		duprow = duplicaterow.index.tolist()
		dupid = duplicateid.index.tolist()
		duplist = duprow + dupid

		console.print(Panel.fit("Duplicate rows or SRA accessions found.\nPlease check your sample_list.tsv file.\nWarning list:\n{}\n\nAt the moment having the same file in different assemblies is not supported, but we are working on it for future versions".format(duplist), title="Sample ID Error", subtitle="Duplicate Sample IDs"))
		sys.exit(1)


	######################
	# SETUP DICTIONARIES #
	######################

	# setup co-assemblies if needed
	assemblies = {}
	for assembly in df.assembly.unique():
		df_filt = df[df['assembly'] == assembly]
		R1s = df_filt['R1'].tolist()
		R2s = df_filt['R2'].tolist()

		try:
			accessions = df_filt['accession'].tolist()
		except (NameError, AttributeError, KeyError):
			accessions = []

		sample_ids = df_filt.index.tolist()

		assemblies[assembly] = {'R1': R1s, 'R2': R2s, 
				'sample_id' : sample_ids, 
				'accession' : accessions}
	
	for sample_id in df.index:
		try:
			R1 = df.loc[sample_id, 'R1']
			R2 = df.loc[sample_id, 'R2']
		except KeyError:
			R1f = '{}_1.fastq.gz'.format(sample_id)
			R1 = '{dir}/{f}'.format(f=R1f, dir=datadir)
			R2f = '{}_2.fastq.gz'.format(sample_id)
			R2 = '{dir}/{f}'.format(f=R2f, dir=datadir)
 		
		if 'accession' in df.columns:
			accession = df.loc[sample_id, 'accession']
		else:
			accession = ''

		try:
			assemblyid = df.loc[sample_id, 'assembly']
		except KeyError:
			print(f"Error: sample_id '{sample_id}' not found in the DataFrame.")
			sys.exit()

		samples[sample_id] = {'R1': R1,
				'R2': R2,
				'accession': accession, 
				'assembly': assemblyid}
	

	# save log of input files
	logdir = os.path.join(outdir, (".vomix/log/vomix" + time))
	
	with open(os.path.join(logdir,  "sample.json"), "w") as samplelog:
		json.dump(samples, samplelog)
	with open(os.path.join(logdir,  "assemblies.json"), "w") as assemblylog:
		json.dump(assemblies, assemblylog)
	
	# check if samples.json already exists or has changed
	# if samples.json is identical, assembly.json is also
	# presumed to be identical
	
	samplejson = os.path.join(outdir, ".vomix/samples.json")
	assemblyjson = os.path.join(outdir, ".vomix/assemblies.json")

	if os.path.exists(samplejson):
		with open(samplejson, "r") as sampleold:
			samples_old = json.load(sampleold)
			if samples_old == samples:
				console.print(Panel.fit("""[bold]Warning[/bold]: [dim]{json} already exists and is identical to the sample list provided {fi}. Skipping validation. If you would like to redo sample validation, run 'rm {json}' and try again.""".format(json = samplejson, fi=f), title="Warning", subtitle="Sample List JSON"))
				return samples, assemblies
				sys.exit()

				
	# If not exited already, validate and write
	Entrez.email = email
	validate_samples(samples)
	
	with open(samplejson, "w") as sampleout:
		json.dump(samples, sampleout)
	with open(assemblyjson, "w") as assemblyout:
		json.dump(assemblies, assemblyout)

	
	return samples, assemblies

