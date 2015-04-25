
## Github repo version of Genome/src 

This is a collection of unix command line program and script utitlities for processing DNA, RNA and protein sequences.  Many of these read naturally from files specified as command line arguments or from stdin when no arguments or given, and typically write output on standard out.   They are intended to be usable in unix pipes.   

Most of these read and write sequence data in fasta format.

### Perl scripts

* **fasta2md5** - reads fasta sequences and writes the md5 hex hashcode along with the sequence header (DNA, RNA, protein)
* **sel**   - (select) reads fasta sequences and filters according to various criteria (DNA, RNA, protein)
* **seq_len** - (sequence length) reads fasta sequences and prints out their lengths (DNA, RNA, protein)
* **trans.pl** - (translate) translates DNA or RNA sequences into protein sequences (fasta format)

### compiled C programs

* **nt** - nucleotide frequencies of DNA or RNA sequences
* **kmers** - nucleotide k-mer counts
* **prosearch** - search DNA for binding motifs specified by patterns
