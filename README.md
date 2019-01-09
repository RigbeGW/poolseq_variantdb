# Variant db from poolSeq

This python script creates a proteogenomic database from Pooled sequencing experiments by interrogating sorted bam alignment files to include more than one version(due to mutations)
of a tryptic peptide. This way it reduces the database size since only those variant sequences containing mutations in a certain order as seen in the alignment file and not all combinations of variant sequences are used.
In addition to the tryptic peptides containing mutations, their 'wild type' protein sequences are written to the database. Also the proteins that do not contain mutations are recorded.

Requirements:

It uses samtools,the pyteomics and biopython packages.

The following flowchart summarizes the steps for the creation of the proteogenomic database from Pool-seq.


<img src="https://github.com/RigbeGW/poolseq_variantdb/blob/master/variantdbFlowchart.png" width="473" height="811">

The script can be invoked from the command line as:

```
usage: Poolseq_variantdb.py [-h] --genomeID GENOMEID --genbankFile GENBANKFILE
                            --bamFile BAMFILE --SnpEffTextOutputFile
                            SNPEFFTEXTOUTPUTFILE
                            [--geneticCodeID GENETICCODEID]
                            [--fastaFile FASTAFILE] [--poolID POOLID]

Create Pool-seq based variant proteome database

optional arguments:
  -h, --help            show this help message and exit
  --genomeID GENOMEID   the ID for the reference genome that is found in the
                        bam
  --genbankFile GENBANKFILE
                        the path to genbank file for the reference genome
  --bamFile BAMFILE     the path to the sorted alignment bam file for
                        retrieving reads from
  --SnpEffTextOutputFile SNPEFFTEXTOUTPUTFILE
                        the path to the text output file from SnpEff
  --geneticCodeID GENETICCODEID
                        the genetic code id of your specieshttps: //
                        www.ncbi.nlm.nih.gov / Taxonomy / Utils / wprintgc.cgi
  --fastaFile FASTAFILE
                        the path to the fasta output file
  --poolID POOLID       Qualifiers to add to the fasta header

```
