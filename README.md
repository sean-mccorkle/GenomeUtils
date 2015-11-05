
## Github repo version of Genome/src 

This is a collection of unix command line program and script utilities mostly for processing DNA, RNA and protein sequences.  Many of these read naturally from files specified as command line arguments or from stdin when no arguments or given, and typically write output on standard out.   They are intended to be usable in unix pipes.   

Most of these read and write sequence data in fasta format.

### Perl scripts

* **codon_freqs** - reports codon frequences of DNA & RNA sequences
* **connected_subgraphs** Reads pairwise edges (node node) and outputs connected subgraphs
* **extract_sequences_list** Read list of ids, extract those fasta sequences with those ids in header
* **fasta2md5** - read fasta sequences and write the md5 hex hashcode along with the sequence header (DNA, RNA, protein)
* **kmers_by_frame** - print out kmer counts in each reading frame (DNA, RNA - can work on proteins too).
* **new_overlapper** - reads two (sorted) tables of chromosome coordinates and reports overlaps between them
* **ngrams** - reads fasta-format sequences, prints counts of ngrams (DNA, RNA, protein)
* **numseqs** - report count of numbers of fasta sequences (DNA, RNA, protein)
* **orfs.pl** - find open reading frames in DNA sequences
* **poly_a** - report long polynucleotide streches in DNA and RNA sequences
* **pp** - fasta format sequence pretty-printer (DNA, RNA, protein)
* **rc** - reverse complement DNA fasta sequences
* **sel**   - (select) reads fasta sequences and filters according to various criteria (DNA, RNA, protein)
* **seq_len** - (sequence length) reads fasta sequences and prints out their lengths (DNA, RNA, protein)
* **seq_range** - (range of sequences filter) reads fasta sequences, prints out only those in specified numeric range  
* **splt** - breaks files of multiple fasta sequences into separate files (DNA, RNA, protein)
* **trans.pl** - (translate) translates DNA or RNA sequences into protein sequences (fasta format)

### compiled C programs

* **intervals** - extract multiple sequence intervals from a large DNA sequence (chromosome sized)
* **kmers** - nucleotide k-mer counts
* **nt** - nucleotide frequencies of DNA or RNA sequences
* **prosearch** - search DNA for binding motifs specified by patterns
