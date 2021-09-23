**Code accompanying manuscript on DADA2 customization**

**File description**

**Preparation of UNITE database for BLAST**
To convert a UNITE database (e.g. the full UNITE + INSD database) to a BLAST database:

*cat unite_2020_general.fasta | perl -ne 's/[^\x00-\x7F]+/ /g; print;' >unite_2020_general_no_ascii.fasta*
*makeblastdb -in unite_2020_general_no_ascii.fasta  -blastdb_version 5 -taxid_map tax_map.txt -title "UNITE with taxid" -dbtype nucl*