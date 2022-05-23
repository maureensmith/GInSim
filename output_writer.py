
import pandas as pd
from pathlib import Path
import random
from yaml import dump
import datetime
strptime = datetime.datetime.strptime
import uuid
import numpy as np
import itertools as it



class writer:
	'''
	Class containing input and output information, reading and writing functions.
	'''

	def __init__(self, outputpath, file_prefix, init_date="2020-01-01"):
		self.fasta_path = outputpath + "/fasta"
		self.table_path = outputpath + "/table"
		self.config_path = outputpath + "/config"
		self.ref_path = outputpath + "/reference"
		self.init_date = init_date
		self.file_prefix = file_prefix

		try:
			Path(self.fasta_path).mkdir(parents=True, exist_ok=True)
		except FileExistsError:
			print("Creation of the directory %s failed" % self.fasta_path)
			exit()
		try:
			Path(self.table_path).mkdir(parents=True, exist_ok=True)
		except FileExistsError:
			print("Creation of the directory %s failed" % self.fasta_path)
			exit()
		try:
			Path(self.config_path).mkdir(parents=True, exist_ok=True)
		except FileExistsError:
			print("Creation of the directory %s failed" % self.config_path)
			exit()
		try:
			Path(self.ref_path).mkdir(parents=True, exist_ok=True)
		except FileExistsError:
			print("Creation of the directory %s failed" % self.ref_path)
			exit()

	def write_reference(self, sequence):
		'''
		Write initial sequence as reference file.
		'''
		outputfile = self.ref_path + "/" + self.file_prefix +".fasta"
		try:
			file = open(outputfile, "w+")
			print("--- Write reference into file " + outputfile + "---")
			file.write(">initialSequence" + "\n" + sequence + "\n")
		except OSError:
			print("Writing of reference file %s failed" % outputfile)
			exit()

	def write_fasta(self, header_prefix, species_dict, file_suffix = None, initialSequence = None, p_shred=None, p_dismiss=None):
		'''
		Write the simulated sequences into one fasta file.
		Returns a table with dates and counts.
		'''
		outputfile = self.fasta_path + "/" + self.file_prefix + "_" + file_suffix + ".fasta"
		rows_list = []

		start_time = strptime(self.init_date, "%Y-%m-%d").date()
		t = 0
		try:
			file = open(outputfile, "w+")

			print("--- Write sequences into file " + outputfile + "---")

			for spec_dict in species_dict:
				date = start_time + datetime.timedelta(days=t)

				sampled_N = 0
				for seq, num in spec_dict.items():
					sampled_N += num
					for _ in range(num):
						#shred sequence, if paramters are given
						shredded_seq = shred_seq(seq, p_shred, p_dismiss)
						#write each part as separate read into fasta
						id=uuid.uuid4()
						for i,s in enumerate(shredded_seq):
							# add unique id to sequence
							header = header_prefix +"|"+ str(id) + "_" + str(i) + "|" + str(date.strftime("%Y-%m-%d"))
							file.write(header + "\n" + s + "\n")

				nMut=sum(spec_dict.values())-spec_dict.get(initialSequence, 0)
				# faster way to add rows to a df
				dict1 = {"t": t, "date": date, "sampled_N": sampled_N, "haplotypes": len(spec_dict), "numMut": nMut}
				rows_list.append(dict1)

				t += 1

		except OSError:
			print("Writing of fasta file %s failed" % outputfile)
			exit()
		return pd.DataFrame(rows_list)


	def write_table(self, table, file_suffix = None):
	    '''
		Write the table with true and sampled sequence counts for each time step
		'''
	    outputfile = self.table_path + "/" + self.file_prefix + "_" + file_suffix + ".tsv"
	    print("--- Write table into file " + outputfile + "---")
	    table.to_csv(outputfile, sep="\t", header=True, index=False)

	def write_config_yaml(self, file_suffix=None):
		'''
		Create the config.yaml for calling the snakemake pipeline
		'''
		file_name = self.file_prefix + "_" + file_suffix
		config_file = self.config_path +  "/config_" + file_name + ".yaml"
		fasta_file = self.fasta_path +  "/" + file_name + ".fasta"
		table_file = self.table_path +  "/" + file_name + ".tsv"
		ref_file =  self.ref_path + "/"+ self.file_prefix+ ".fasta"

		config_dict={"samples":  fasta_file,
		            # "reported_cases" : [table_file, "\t", "date", "sampled_N", "%Y-%m-%d"],
		             "reported_cases": [table_file, "\t", "date", "true_N", "%Y-%m-%d"],
		             "consensus": ref_file,
		             "number_per_bin": [],
		             "days_per_bin": [2, 4, 6, 8, 10],
		             "min_bin_size": 15,
		             "min_days_span": 2,
		             "max_days_span": 20,
		             "freq_cutoff": 0,
		             "group": file_name,
					 "R0": "n"
		             }


		with open(config_file, 'w') as file:
			dump(config_dict, file)


def get_slices(seq, l, p_dismiss=None):
	it1, it2 = it.tee(l)

	next(it2)
	for start, end in zip(it1, it2):
		include_subread = True
		# dismiss each snippet with a certain probability
		if p_dismiss is not None:
			include_subread = random.random() > p_dismiss

		if include_subread:
			yield seq[start:end]

def shred_seq(seq, p_shred=None, p_dismiss=None):
	print("whole sequence ")
	print(seq)
	if p_shred is not None:
		# draw number of cuts
		n_cut = np.random.poisson((len(seq) - 1) * p_shred)
		print("Number of cuts " + str(n_cut))
		print(len(seq))
		cutpos = random.sample(range(1, len(seq) - 1), n_cut)
		cutpos.sort()
		print(cutpos)
		cut_positions = [0] + cutpos + [len(seq)]
		print("Cut positions")
		print(cut_positions)
		sub_seqs = list(get_slices(seq, cut_positions, p_dismiss))
	else:
		sub_seqs = [seq]

	print("sub sequence ")
	for x in sub_seqs:
		print(x)

	return sub_seqs
