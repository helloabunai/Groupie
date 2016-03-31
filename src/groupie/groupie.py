#!/usr/bin/python
## Alastair Maxwell -- University of Glasgow
## alastair.maxwell@glasgow.ac.uk

import argparse
import sys
import os
import glob
import subprocess
import csv
import datetime
import logging
import multiprocessing

class clr:

	def __init__(self):

		pass

	purple = '\033[95m'
	cyan = '\033[96m'
	darkcyan = '\033[36m'
	blue = '\033[94m'
	green = '\033[92m'
	yellow = '\033[93m'
	red = '\033[91m'

	bold = '\033[1m'
	underline = '\033[4m'
	end = '\033[0m'

class GroupieException(Exception):
	pass

class integrity:

	def __init__(self):

		self.script_dir = sys.path[0]

	@staticmethod
	def general_test(parser):

		## Argument checker
		if not len(sys.argv) > 1:
			parser.print_help()
			sys.exit(0)

		in_dir = ''
		out_dir = ''

		arguments = parser.parse_args()
		if arguments.sam is not None:
			if type(arguments.sam) is list:
				in_dir = arguments.sam[0]
			if type(arguments.sam) is str:
				in_dir = arguments.sam

		if arguments.out is not None:
			if type(arguments.out) is list:
				out_dir = arguments.out[0]
			if type(arguments.out) is str:
				out_dir = arguments.out

		directories = {'in_dir': in_dir, 'out_dir': out_dir}
		return directories

	@staticmethod
	def input_test(input_dir):

		## Input data test method
		## List of potential files to be batched
		group_targets = []

		## Checks for input validity
		if not os.path.exists(input_dir):
			raise GroupieException("ERROR: Specified input directory does not exist!")

		## If the input dir specified is for a single file..
		## Print for UXP, check file ends with correct format
		## Append to list, return to main
		if os.path.isfile(input_dir):
			print '{}{}{}{}'.format(clr.bold, 'groupie__', clr.end, ' Grouping an individual file...')
			if input_dir.endswith(".sam"):
				genotype_target = os.path.abspath(input_dir)
				print '{}{}{}{}{}'.format(clr.bold, clr.green, 'groupie__', clr.end, ' Input directory/files OK!')
				return genotype_target

		## If the input dir specified is for a folder of files..
		## Print for UXP, check all files in folder for format
		## Append to list, return to main
		if os.path.isdir(input_dir):
			print '{}{}{}{}'.format(clr.bold, 'groupie__', clr.end, ' Grouping multiple files...')

			for candidate in glob.glob(os.path.join(input_dir, '*')):
				if candidate.endswith(".sam"):
					group_targets.append(os.path.abspath(candidate))
			print '{}{}{}{}{}'.format(clr.bold, clr.green, 'groupie__', clr.end, ' Input directory/files OK!')
			return group_targets

	@staticmethod
	def output_test(out_dir_root):

		## Ensures root output is a real directory
		## Generates folder name based on date (for run ident)
		date = datetime.date.today().strftime('%d-%m-%Y')
		time = datetime.datetime.now().strftime('%H%M%S')
		today = date + '-' + time

		## If the user specified root doesn't exist, make it
		## Then make the run directory for datetime
		if not os.path.exists(out_dir_root):
			print '{}{}{}{}'.format(clr.bold, 'groupie__', clr.end, ' Creating output root... ')
			os.mkdir(out_dir_root)
		run_dir = out_dir_root + '/GroupieRun' + today
		print '{}{}{}{}'.format(clr.bold, 'groupie__', clr.end, ' Creating instance run directory.. ')
		os.mkdir(run_dir)

		## Inform user it's all gonna be okaaaayyyy
		print '{}{}{}{}{}'.format(clr.bold, clr.green, 'groupie__', clr.end, ' Output directories OK!')
		return run_dir

	@staticmethod
	def thread_test(parser):

		cpu_threads = 0
		arguments = parser.parse_args()
		if arguments.cpu is not None:
			if type(arguments.cpu) is list:
				cpu_threads = arguments.cpu[0]
			if type(arguments.cpu) is int:
				cpu_threads = arguments.cpu

		if int(cpu_threads) > multiprocessing.cpu_count():
			raise GroupieException("ERROR: Specified CPU threads greater than available on system!")
		else:
			print '{}{}{}{}'.format(clr.bold, 'groupie__', clr.end, ' Utilising ' + str(cpu_threads) + ' threads...')
			return str(cpu_threads)

	@staticmethod
	def library_test():

		which_samtools = subprocess.Popen(['which', 'samtools'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
		samtools_dir = which_samtools.communicate()[0]
		which_samtools.wait()

		if not 'samtools' in samtools_dir:
			raise GroupieException('Samtools Binary missing from $PATH, "which samtools == null". Please install samtools.')
		else:
			print '{}{}{}{}{}'.format(clr.bold, clr.green, 'groupie__', clr.end, ' Libraries OK!')
			return True

	@staticmethod
	def csv_cleanup(input_outdir, sample_name):

		raw_distro = os.path.join(input_outdir, 'raw_repeatdistro.csv')
		cln_distro = os.path.join(input_outdir, 'cln_repeatdistro.csv')

		with open(raw_distro) as csv_file:

			file_root = str(sample_name.split('.')[0]).split('/')[-1]
			repeat_values = []
			data_string = ''

			for line in csv_file.readlines()[:-1]:
				values = line.split('\t')
				data_string += values[0] + ',' + values[1] +',' + values[2] + ',0\n'

		filestring = str(len(repeat_values))+',3,' + file_root + '\n'
		filestring += data_string
		cleaned_file = open(cln_distro, 'w')
		cleaned_file.write(filestring)
		cleaned_file.close()

class groupie:

	def __init__(self):

		print '{}{}{}{}{}'.format('\n', clr.bold, 'groupie__', clr.end, ' Welcome to Groupie!')
		print '{}{}{}{}{}'.format(clr.bold, 'groupie__', clr.end, ' alastair.maxwell@glasgow.ac.uk', '\n')

		self.input_dir = ''
		self.output_dir = ''

		defpath = os.path.join(os.path.expanduser('~'),'Groupie')
		self.parser = argparse.ArgumentParser(description='Groupie: Simple repeat count distribution batching script.')
		self.parser.add_argument('-sam', help='Path to a *.sam file, or a folder of *.sam files for input.', nargs=1)
		self.parser.add_argument('-rng', help='Numerical range to batch targets reads into.', default=10, type=int, nargs=1)
		self.parser.add_argument('-cpu', help='Number of processor threads to use in supported mutli-threaded libraries.', nargs=1, default=multiprocessing.cpu_count())
		self.parser.add_argument('-out', help='Path to direct output to.', nargs=1, default=defpath)
		self.arguments = self.parser.parse_args()
		self.group_range = self.arguments.rng

		##
		## Argument checks etc
		self.data_dirs = integrity.general_test(self.parser)
		self.input_dir = self.data_dirs['in_dir']
		self.output_dir = self.data_dirs['out_dir']

		self.group_targets = integrity.input_test(self.input_dir)
		self.run_dir = integrity.output_test(self.output_dir)
		self.cpu_threads = integrity.thread_test(self.parser)
		integrity.library_test()

		##
		## str = single file, no loop
		## lst = dir of files, loop(x)
		if type(self.group_targets) is str:
			self.single_flow()
		if type(self.group_targets) is list:
			self.multi_flow()

	def single_flow(self):

		print '\n'.lstrip()

		##
		## Prepare Data
		print '{}{}{}{}'.format(clr.bold, 'groupie__', clr.end, ' Processing repeat distributions... ')
		self.extract_distro(self.group_targets)
		self.group_distro(self.group_targets)
		print '{}{}{}{}{}'.format('\n',clr.bold, 'groupie__', clr.end, ' Output saved to ' + self.run_dir)
		print '{}{}{}{}'.format(clr.bold, 'groupie__', clr.end, ' Exiting...')

	def multi_flow(self):

		##
		## Prepare Data
		print '{}{}{}{}'.format(clr.bold, 'groupie__', clr.end, ' Processing repeat distributions... ')
		for i in range(0, len(self.group_targets)):
			try:
				self.extract_distro(self.group_targets[i])
				self.group_distro(self.group_targets[i])
			except:
				continue
		print '{}{}{}{}{}'.format('\n',clr.bold, 'groupie__', clr.end, ' Output saved to ' + self.run_dir)
		print '{}{}{}{}'.format(clr.bold, 'groupie__', clr.end, ' Exiting...')

	def extract_distro(self, samfile):

		"""
		Function to extract repeat distribution from aligned sequence in sam format
		Program ends up here if the user specified sam input files
		Distribution produced from this section is then utilised for machine learning
		"""

		## Generate a folder for this file's outputs
		## and all associated data
		filename = samfile.split('/')[-1].split('.')[0]
		input_outdir = os.path.join(self.run_dir, filename)
		os.makedirs(input_outdir)

		## Samtools view and sort, index and idxstats
		## End product = extracted repeat distribution for this file
		samtools_view = subprocess.Popen(['samtools', 'view','-bS','-@', self.cpu_threads,samfile],stdout=subprocess.PIPE)
		samtools_sort = subprocess.Popen(['samtools', 'sort','-@', self.cpu_threads,'-', os.path.join(input_outdir,'sorted_assembly')],stdin=samtools_view.stdout)
		samtools_view.wait(); samtools_sort.wait()

		samtools_indx = subprocess.Popen(['samtools', 'index',os.path.join(input_outdir,'sorted_assembly.bam')],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
		samtools_indx.wait()

		read_expr_dir = os.path.join(input_outdir, 'raw_repeatdistro.csv')
		read_expr_file = open(read_expr_dir, 'w')

		samtools_idxs = subprocess.Popen(['samtools', 'idxstats',os.path.join(input_outdir,'sorted_assembly.bam')],stdout=read_expr_file)
		samtools_idxs.wait(); read_expr_file.close()

		##
		## Clean the raw distribution up
		## So it interfaces with rest of program OK
		integrity.csv_cleanup(input_outdir, filename)

	def group_distro(self, samfile):

		filename = samfile.split('/')[-1].split('.')[0]
		input_outdir = os.path.join(self.run_dir, filename)
		distro_path = os.path.join(input_outdir, 'cln_repeatdistro.csv')

		grouped_expressions = []
		grouped_results = []
		grouped_path = os.path.join(input_outdir, 'grouped_distribution.csv')


		with open(distro_path, 'rb') as csvfile:
			reader = csv.reader(csvfile, delimiter=',')

			data = list(reader)
			row_count = len(data[2:])

			if not row_count % self.group_range == 0:
				raise GroupieException('ERROR: Specified group range does not evenly divide into row count of input CSV.')

			i=2
			while True:
				lst_block = data[i:i+self.group_range]
				grouped_expressions.append(lst_block)
				i+= self.group_range
				if i >= row_count:
					break

		for block in grouped_expressions:

			subrange = []
			subreads = []

			for subline in block:
				range_figure = subline[0].split('_')[3]
				subrange.append(range_figure)
				read_indice = subline[2]
				subreads.append(int(read_indice))

			range_start = subrange[0]
			range_end = subrange[-1]
			range_string = '{}{}{}'.format(str(range_start),'-',str(range_end))
			read_total = sum(subreads)

			to_csv = [range_string, read_total]
			grouped_results.append(to_csv)

		with open(grouped_path, 'wb') as fp:
			a = csv.writer(fp, delimiter=',')
			a.writerows(grouped_results)

def main():

	try:
		groupie()
	except:
		logdir = os.path.join(os.path.expanduser('~'),'Groupie','system_logs')
		rundate = datetime.date.today().strftime('%d-%m-%Y')
		runtime = datetime.datetime.now().time().strftime('%H-%M')

		if not os.path.isdir(logdir): os.makedirs(logdir)
		filen = rundate + '_' + runtime + '.log'; logfi = os.path.join(logdir, filen)
		logging.basicConfig(filename=logfi,level=logging.DEBUG)
		logging.exception('Main thread exception handler. Check error log.')
		print '{}{}{}{}{}'.format(clr.bold, clr.red, 'groupie__', clr.end, ' Fatal Exception.\n\tLogged to ' + logdir + '\n\tExiting Groupie.. ')
		sys.exit(2)

if __name__ == '__main__':
	main()