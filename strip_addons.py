#! /usr/bin/env python2.7
import sys,os,argparse
import yingseqtools as y

def main():
	#def strip_barcode_addons(seqfile1,seqfile2,pdiffs=0,separate_barfile=False,test=False):
	parser = argparse.ArgumentParser(description='YT: Trims off barcode addons. Output will be named reads1.fastq and reads2.fastq. Sequences that do not match both forward and reverse are stored in scrap.fastq.')
	parser.add_argument('seqfile1',help='Forward fastq file.')
	parser.add_argument('seqfile2',help='Reverse fastq file.')
	parser.add_argument('-pdiffs',type=int,help='Number of primer bp mismatches allowed (default is 0).',metavar='<int>',default=0)
	parser.add_argument('-remove_bar_primer',action='store_true',help='In addition to stripping addons, remove barcode and primer, and store barcode separately in barcodes.fastq. Default=false, where barcode and primer stay on the sequences.')
	parser.add_argument('-fw_primer',type=str,help='Forward primer. Default is 16S V4: AYTGGGYDTAAAGNG')
	parser.add_argument('-rev_primer',type=str,help='Reverse primer. Default is 16S V4: CCGTCAATTYHTTTRAGT.')
	parser.add_argument('-test',action='store_true',help='Test mode. Only does first 300 seqs.')
	args = parser.parse_args()
	y.strip_barcode_addons(**vars(args))

if __name__ == '__main__':
	main()


