#! /usr/bin/env python2.7
import sys,os,glob,shutil,fnmatch,random,itertools,gzip
import regex as re
from operator import add
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline
import ref_file_settings as ref

strip_barcode = ref.strip_barcode
#strip_barcode = ref.file['strip_barcode']


#strip barcode addons at beginning of fastq files
def strip_barcode_addons(seqfile1,seqfile2,pdiffs=0,remove_bar_primer=False,fw_primer=None,rev_primer=None,test=False):
	if fw_primer is None:
		fw_primer = 'AYTGGGYDTAAAGNG'
	if rev_primer is None:
		rev_primer = 'CCGTCAATTYHTTTRAGT'

	print 'primer',fw_primer,rev_primer
	#can be fasta, fastq, or fastq.gz
	seq_pattern = re.compile('.(fasta|fastq)(.gz)?$')
	assert seq_pattern.search(seqfile1) is not None, 'YTError:'+seqfile1+'needs to end in .fasta, .fastq, or .fastq.gz!'
	assert seq_pattern.search(seqfile2) is not None, 'YTError:'+seqfile2+'needs to end in .fasta, .fastq, or .fastq.gz!'
	seq_format = seq_pattern.search(seqfile1).group(1)
	seq_format2 = seq_pattern.search(seqfile2).group(1)
	gzipped = seq_pattern.search(seqfile1).group(2) is not None
	gzipped2 = seq_pattern.search(seqfile2).group(2) is not None
	assert seq_format==seq_format2 and gzipped==gzipped2, 'YTError: These 2 sequence files do not match!'
	print 'Trimming off barcode addons...'
	#fw_primer = 'AYTGGGYDTAAAGNG'
	#rev_primer = 'CCGTCAATTYHTTTRAGT'
	#fw_primer_regex = re.compile('(A[CT]TGGG[CT][AGT]TAAAG[ACGT]G){e<='+str(pdiffs)+'}')
	#rev_primer_regex = re.compile('(CCGTCAATT[CT][ACT]TTT[AG]AGT){e<='+str(pdiffs)+'}')
	fw_primer_regex = re.compile(primer_pattern(fw_primer,pdiffs=pdiffs))
	rev_primer_regex = re.compile(primer_pattern(rev_primer,pdiffs=pdiffs))
	len_barcode = 12 #golay barcodes
	max_len_addon = 8 #addon is 1-8 bp
	max_len_fw = len_barcode + max_len_addon + len(fw_primer)
	max_len_rev = len_barcode + max_len_addon + len(rev_primer)
	trim_seqfile1 = 'reads1.fastq'
	trim_seqfile2 = 'reads2.fastq'
	scrap_seqfile = 'scrap.fastq'
	bar_seqfile = 'barcodes.fastq'
	if gzipped: #gzipped
		seq1 = gzip.open(seqfile1,'rb')
		seq2 = gzip.open(seqfile2,'rb')
	else: #not gzipped
		seq1 = open(seqfile1,'rb')
		seq2 = open(seqfile2,'rb')
	newseq1 = open(trim_seqfile1,'wb')
	newseq2 = open(trim_seqfile2,'wb')
	badseq = open(scrap_seqfile,'wb')
	if remove_bar_primer:
		barcode = open(bar_seqfile,'wb')
	n_forward_reverse = 0 #counts forward>reverse seqs (p1-p2)
	n_reverse_forward = 0 #counts reverse>forward seqs (p2-p1)
	n_noP1 = 0 #counts forward non-hits
	n_noP2 = 0 #counts reverse non-hits
	i=0
	for (forward,reverse) in itertools.izip(SeqIO.parse(seq1,seq_format),SeqIO.parse(seq2,seq_format)):
		assert forward.id==reverse.id, 'YTError: sequence IDs do not match!'
		#these will erase the second half of the header (id=first part,description=whole line)
		#for compatibility with usearch
		forward.description=forward.id
		reverse.description=reverse.id
		#look for forward primer in seq1 and reverse primer in seq2
		r_front = fw_primer_regex.search(str(forward.seq[:max_len_fw]))
		if r_front is not None: #forward primer hit seq1, look for rev primer in seq2.
			r_end = rev_primer_regex.search(str(reverse.seq[:max_len_rev]))
			if r_end is not None: #seq1 is forward, seq2 is reverse
				primerstart1 = r_front.start()
				barstart1 = max(primerstart1-len_barcode,0)
				seqstart1 = primerstart1+len(fw_primer)
				primerstart2 = r_end.start()				
				barstart2 = max(primerstart2-len_barcode,0)
				seqstart2 = primerstart2+len(rev_primer)
				if remove_bar_primer:
					SeqIO.write(forward[seqstart1:],newseq1,seq_format)
					SeqIO.write(reverse[seqstart2:],newseq2,seq_format)
					bar1 = forward[barstart1:primerstart1]
					bar2 = reverse[barstart2:primerstart2]
					bar = bar1 + bar2
					bar.description = bar1.description
					SeqIO.write(bar,barcode,seq_format)
				else:
					SeqIO.write(forward[barstart1:],newseq1,seq_format)
					SeqIO.write(reverse[barstart2:],newseq2,seq_format)
				n_forward_reverse += 1
			else: #seq1 is forward but seq2 nohit, ERROR
				forward.id = forward.id+'|noP2'
				reverse.id = reverse.id+'|noP2'
				SeqIO.write(forward,badseq,seq_format)
				SeqIO.write(reverse,badseq,seq_format)
				n_noP2 += 1
		else: #first primer nohit, look for rev primer in seq1
			r_front = rev_primer_regex.search(str(forward.seq[:max_len_rev]))
			if r_front is not None: #rev primer hit seq1, look for forward primer in seq2.
				r_end = fw_primer_regex.search(str(reverse.seq[:max_len_fw]))
				if r_end is not None: #seq1 is reverse primer, seq2 is forward primer.
					primerstart1 = r_front.start()
					barstart1 = max(primerstart1-len_barcode,0)
					seqstart1 = primerstart1+len(rev_primer)
					primerstart2 = r_end.start()				
					barstart2 = max(primerstart2-len_barcode,0)
					seqstart2 = primerstart2+len(fw_primer)
					if remove_bar_primer:
						#SeqIO.write(forward[seqstart1:],newseq1,seq_format)
						#SeqIO.write(reverse[seqstart2:],newseq2,seq_format)
						SeqIO.write(forward[seqstart1:],newseq2,seq_format)
						SeqIO.write(reverse[seqstart2:],newseq1,seq_format)
						bar1 = forward[barstart1:primerstart1]
						bar2 = reverse[barstart2:primerstart2]
						#bar = bar1 + bar2
						bar = bar2 + bar1
						#bar.description = bar1.description
						bar.description = bar2.description
						SeqIO.write(bar,barcode,seq_format)
					else:
						#SeqIO.write(forward[barstart1:],newseq1,seq_format)
						#SeqIO.write(reverse[barstart2:],newseq2,seq_format)
						SeqIO.write(forward[barstart1:],newseq2,seq_format)
						SeqIO.write(reverse[barstart2:],newseq1,seq_format)
					n_reverse_forward += 1
				else: #rev primer hit seq1, but seq2 is nohit for forward primer.
					forward.id = forward.id+'|noP2'
					reverse.id = reverse.id+'|noP2'
					SeqIO.write(forward,badseq,seq_format)
					SeqIO.write(reverse,badseq,seq_format)
					n_noP2 += 1
			else: #seq1 is nohit for both primers
				forward.id = forward.id+'|noP1'
				reverse.id = reverse.id+'|noP1'
				SeqIO.write(forward,badseq,seq_format)
				SeqIO.write(reverse,badseq,seq_format)
				n_noP1 += 1
		if test and n_reverse_forward+n_forward_reverse>=300:
			break
	seq1.close()
	seq2.close()
	newseq1.close()
	newseq2.close()
	badseq.close()
	if remove_bar_primer:
		barcode.close()
	totalseqs = n_forward_reverse+n_reverse_forward+n_noP1+n_noP2
	print 'forward-reverse primer:'+str(n_forward_reverse)+'('+str(round(float(n_forward_reverse)/totalseqs*100,1))+'%)'
	print 'reverse-forward primer: '+str(n_reverse_forward)+'('+str(round(float(n_reverse_forward)/totalseqs*100,1))+'%)'
	print 'no P1 match: '+str(n_noP1)+'('+str(round(float(n_noP1)/totalseqs*100,1))+'%)'
	print 'no P2 match: '+str(n_noP2)+'('+str(round(float(n_noP2)/totalseqs*100,1))+'%)'
	print 'Total Seqs: '+str(totalseqs)+' (100%)'

def blastn(fastafile,db='refseq_rna',taxoutfile=None,evalue=10,max_target_seqs=50,outfmt='qseqid qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send nident mismatch gapopen gaps ppos frames qframe sframe sblastnames sskingdoms stitle salltitles sstrand qcovs qcovhsp evalue bitscore score length pident sscinames scomnames staxids'):
	#command line equivalent: blastn -query fastafile -db refseq_rna -evalue 10 -max_target_seqs 50 -out taxoutfile -outfmt "6 qseqid sseqid evalue bitscore score length pident sscinames"
	#performs blastn on fasta file
	names_file = ref.names_file
	nodes_file = ref.nodes_file
	merged_file = ref.merged_file

	assert re.search('.(fasta|fna)$',fastafile), 'YTError: sequence file does not end in fasta!'
	#assert fastafile.endswith('fasta'), 'YTError: sequence file does not end in fasta!'
	assert db in ref.blastdb_list, 'YTError, blast db not found locally!'
	if taxoutfile is None:
		outfile = fastafile + '.blastn.' + db + '.txt'
	else:
		outfile = taxoutfile
	outfile0 = outfile + '.temp'
	logfile = outfile + '.log'
	print 'Running blast....'
	blastn_cline = NcbiblastnCommandline(query=fastafile, db=db, evalue=evalue, max_target_seqs=max_target_seqs, outfmt='"6 '+outfmt+'"', out=outfile0)
	#print blastn_cline #look at it, will print: blastn -out "C:\Users\taury\Desktop\Patric Files\cbotoutput2.xml" -outfmt 5 -query "C:\Users\taury\Desktop\Patric Files\cbot.txt" -db refseq_genomic -evalue 0.001
	print 'Running blastn:', fastafile
	stdout, stderr = blastn_cline() #run it (might take a while)
	#read in taxonomy
	print 'Calculating full taxonomy....'
	tax = {}
	with open(nodes_file,'rb') as nodes:
		for line in nodes:
			values = line.split('\t|\t')
			taxid = long(values[0])
			parentid = long(values[1])
			rank = values[2]
			division = values[4]
			tax[str(taxid)] = {'parentid':parentid,'rank':rank,'division':division}
	with open(names_file,'rb') as names:
		for line in names:
			values = line.split('\t|\t')
			taxid = long(values[0])
			taxname = values[1]
			taxnameclass = re.sub('\t\|\n','',values[3])
			if taxnameclass=='scientific name':
				tax[str(taxid)]['taxname'] = taxname

	with open(merged_file,'rb') as merged:
		for line in merged:
			values = line.split('\t|\t')
			oldtaxid = values[0]
			newtaxid = re.sub('\t\|\n','',values[1])
			try:
				tax[str(oldtaxid)] = tax[str(newtaxid)]
			except KeyError as e:
				print 'YTError: Encountered merged taxid ',oldtaxid,', but unable to find new taxid, ',newtaxid,'.'
				raise

	#put header on file.
	def getlineage(taxid,abbrev=False):
		if abbrev:
			full_lineage = getlineage(taxid).split('|')
			#levels = ['[superkingdom]','[phylum]','[class]','[order]','[family]','[genus]','[species]']
			levels = ['[superkingdom]','[kingdom]','[subkingdom]','[superphylum]','[phylum]','[subphylum]','[superclass]','[class]','[subclass]','[infraclass]','[superorder]','[order]','[suborder]','[infraorder]','[parvorder]','[superfamily]','[family]','[subfamily]','[tribe]','[subtribe]','[genus]','[subgenus]','[species group]','[species subgroup]','[species]','[subspecies]','[varietas]','[forma]']
			abbrev_lineage = []
			for lvl in levels:
				taxhits = [t for t in full_lineage if t.endswith(lvl)]
				if len(taxhits)==0:
					abbrev_lineage.append('xxxx'+lvl)
				elif len(taxhits)==1:
					abbrev_lineage.append(taxhits[0])
				else:
					raise Exception('YTError: in getlineage, more than one hit for a level:',taxid,taxhits,full_lineage) 
					#abbrev_lineage.append(taxhits[0])
			return('|'.join(abbrev_lineage))
		
		try:
			taxon = tax[str(taxid)]
		except KeyError as e:
			print 'YTError: taxid',taxid,'not found! Consider updating the taxdmp files.'
			raise
		name = taxon['taxname'] + '[' + taxon['rank'] + ']'
		parentid = taxon['parentid']
		rank = taxon['rank']
		if long(taxid) == 1:
			return name 
		else:
			return getlineage(parentid) + '|' + name
	fields = outfmt.split(' ')
	assert 'staxids' in fields, 'YTError: could not find [staxids] field in output!'
	taxid_pos = long([i for i,f in enumerate(fields) if f=='staxids'][0])
	header = '\t'.join(fields)
	out0 = open(outfile0,'rb')
	out = open(outfile,'wb')
	out.write(header+'\ttaxonomy\tfulltaxonomy\n')
	for line in out0:
		values = re.sub('\n$','',line).split('\t')
		taxid = values[taxid_pos]
		taxid = taxid.split(';')[0]
		shorttax = getlineage(taxid,abbrev=True)
		fulltax = getlineage(taxid)
		values.append(shorttax)
		values.append(fulltax)
		newline = '\t'.join(values) + '\n'
		out.write(newline)
	out0.close()
	out.close()
	os.remove(outfile0)
	with open(logfile,'wb') as log:
		log.write('Input Sequence File: '+fastafile+'\n')
		log.write('BLASTN Results: '+outfile+'\n')
		log.write('Database (db): '+db+'\n')
		log.write('E-value threshold: '+str(evalue)+'\n')
		log.write('Number of hits/seq (max_target_seqs): '+str(max_target_seqs)+'\n')
	print 'Done!'