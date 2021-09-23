# **Code accompanying manuscript on DADA2 customization**

## **File descriptions**

*yingseqtools.py*: Python file which is called upon by specific commands used below
*preparation_from_raw_reads.R*: This script removes primers, sorts forward and reverse reads into separate files by using *strip_addons.py* and then removes primer run-in via cutadapt
*process_dada2.R* runs dada2 on the pre-processed files (filtering variables maxEE and truncQ can be changed in this script)
*dada2_pipeline.R* combines denoised files into one phyloseq object, runs taxonomic annotation (you can chose UNITE or NT for a BLAST-based approach, or use RDP), and includes a tracing file
*refseq.txt* and *refseqs.fasta* include expected sequences from our mock community
*test_dada2_variables.R* script used to test different values for maxEE, truncQ, minOverlap, and maxMismatch in the DADA2 pipeline
*Qualityprofile.R* script used to compute and visualize quality per cycle per read and to compute expected error per read

## **Preparation of UNITE database for BLAST**

To convert a UNITE database to a BLAST database:

```
cat unite_2020_general.fasta | perl -ne 's/[^\x00-\x7F]+/ /g; print;' >unite_2020_general_no_ascii.fasta
makeblastdb -in unite_2020_general_no_ascii.fasta  -blastdb_version 5 -taxid_map tax_map.txt -title "UNITE with taxid" -dbtype nucl
```