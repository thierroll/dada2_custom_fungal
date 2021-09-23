#! /usr/bin/env python2
import sys,os,argparse
import yingseqtools as y
import ref_file_settings as ref

def main():
	#def strip_barcode_addons(seqfile1,seqfile2,pdiffs=0,separate_barfile=False,test=False):
	parser = argparse.ArgumentParser(description='YT: Runs blastn on a fasta file. Creates output file containing taxonomy info')
	parser.add_argument('fastafile',help='Sequence file to be queried with BLASTN.')
	parser.add_argument('-db',choices=ref.blastdb_list,help='NCBI database to be used. Default is refseq_rna',default='refseq_rna')
	parser.add_argument('-taxoutfile',help='Filename for taxonomy output. Default is [fastafile].blastn.txt',default=None)
	parser.add_argument('-evalue',type=float,help='E-value threshold. Default 0.001',default=10)
	parser.add_argument('-max_target_seqs',type=int,help='Number of hits to return per sequence. Default is 50',default=50)   #changed to 50 as otherwise some Saccharomyces sequences will not be identified
	args = parser.parse_args()
	y.blastn(**vars(args))

if __name__ == '__main__':
	main()


