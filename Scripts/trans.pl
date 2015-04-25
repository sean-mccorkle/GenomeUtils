#!/usr/bin/perl
#
# Program:      trans.pl
# Programmer:   Sean R. McCorkle, Biology Dept. Brookhaven National Laboratory
# Language:     perl
# Description:  Read DNA sequences, translate into protein (amino acid) 
#               sequences
#
# Usage:        trans.pl [-f<frame>]  [<seqfile> ...]
#
#               where <seq-file>'s are in FASTA format.  Stdin is scanned
#               if no files are specified.
#
# Options:
#               -f <frame>    translate frame <frame> which is one of 
#                             F1, F2, F3, R1, R2, R3 or a commma-separated
#                             list of these for more than one, or "all"
#                             which stands for "F1,F2,F3,R1,R2,R3"
#                             (defaults to F1) 
#
# Output:
#               fasta format protein sequences.   * indicates a stop codon,
#               X indicates a non-translatable codon (i.e. a codon which
#               is short (at the end) or a codon which contains a letter 
#               other than A, C, G, or T
#
# TODO:
#               1) add "R1-3", etc as allowable frames
#               2) use strict
#
# $Id: trans.pl,v 1.1 2013/02/27 16:49:54 seanmccorkle Exp seanmccorkle $
#

use  Getopt::Std;

getopts( "f:");

                               ###########
                               # Globals #
                               ###########

$frame = $opt_f ? check_frame( $opt_f ) : "F1";

$line_len = 50;                         # for fasta output

                          ######################
                          # Codon translations #
                          ######################

%trans = (
          'AAA' => 'K', 'AAC' => 'N', 'AAG' => 'K', 'AAT' => 'N', 'ACA' => 'T',
          'ACC' => 'T', 'ACG' => 'T', 'ACT' => 'T', 'AGA' => 'R', 'AGC' => 'S',
          'AGG' => 'R', 'AGT' => 'S', 'ATA' => 'I', 'ATC' => 'I', 'ATG' => 'M',
          'ATT' => 'I', 'CAA' => 'Q', 'CAC' => 'H', 'CAG' => 'Q', 'CAT' => 'H',
          'CCA' => 'P', 'CCC' => 'P', 'CCG' => 'P', 'CCT' => 'P', 'CGA' => 'R',
          'CGC' => 'R', 'CGG' => 'R', 'CGT' => 'R', 'CTA' => 'L', 'CTC' => 'L',
          'CTG' => 'L', 'CTT' => 'L', 'GAA' => 'E', 'GAC' => 'D', 'GAG' => 'E',
          'GAT' => 'D', 'GCA' => 'A', 'GCC' => 'A', 'GCG' => 'A', 'GCT' => 'A',
          'GGA' => 'G', 'GGC' => 'G', 'GGG' => 'G', 'GGT' => 'G', 'GTA' => 'V',
          'GTC' => 'V', 'GTG' => 'V', 'GTT' => 'V', 'TAA' => '*', 'TAC' => 'Y',
          'TAG' => '*', 'TAT' => 'Y', 'TCA' => 'S', 'TCC' => 'S', 'TCG' => 'S',
          'TCT' => 'S', 'TGA' => '*', 'TGC' => 'C', 'TGG' => 'W', 'TGT' => 'C',
          'TTA' => 'L', 'TTC' => 'F', 'TTG' => 'L', 'TTT' => 'F',
         );


                            ################
                            # Main Program #
                            ################

while ( <> )
   {
    if ( /^>(.*)/ )
       {
        $new_hdr = $1;
        translate_and_print( $seq, $hdr ) if ( defined( $hdr ) );
        $hdr = $new_hdr;
        $seq = "";
       }
    else
       {
        chomp;
        s/\s//g;
        $seq .= $_;
       }
   }

translate_and_print( $seq, $hdr ) if ( defined( $hdr ) );

                             ###############
                             # Subroutines #
                             ###############

sub  translate_and_print
   {
    my ( $seq, $hdr ) = @_;
    my $dna;

    foreach my $fr ( split( /,/, $frame ) )
       {
        $fr =~ /^([FR])([123])$/ || die "bad frame \"$frame\"\n";
        my ($dir,$off) = ($1,$2-1);

        # reframe DNA sequence as needed

        if ( $dir eq 'F' )
           { $dna = uc( substr( $seq, $off ) ); }
        else
           { $dna = uc( substr( rc( $seq ), $off ) ); }
    
        # truncate length to a multiple of 3 if necessary

        if ( (length( $dna ) % 3) > 0 )
           {
            $dna = substr( $dna, 0, 3 * int( length( $dna ) / 3) );
           }

        $prot = translate( $dna );
        my $newhdr = $hdr;
        $newhdr .= "-$fr";
        print_fasta( $newhdr, $prot );
       }
   }

sub  translate
   {
    my $dna = shift;
    my $l = length( $dna );
    die "bad length $l\n" if ( $l % 3 );
    my $prot = "";
    for ( my $j = 0;  $j < $l;  $j += 3 )
       {
        my $codon = substr( $dna, $j, 3 );
        $prot .= $trans{$codon} ? $trans{$codon} : 'X';
       }
    return( $prot );
   }

sub  print_fasta
   {
    my ( $hdr, $seq ) = @_;
    my $i;

    print ">$hdr\n";
    my $n = length( $seq );                       # calculate length once
    for ( $i = 0; $i < $n; $i++ )                 # for each character
       {
        print substr( $seq, $i, 1 );              # print it out
        if ( ( ($i + 1) % $line_len ) == 0 )      # every 50 characters
            { print "\n"; }                       #   put a newline
       }
    print "\n" if ( ( $i % $line_len ) != 0 );    # last newline ensures tidy 
                                                  #output
   }

#
# rc( $seq ) returns the reverse complement of $seq
# 
sub  rc
   {
    my $r = reverse( split( '', shift ) );
    $r =~ tr/acgtmrwsykvhdbnACGTMRWSYKVHDBN/tgcakywsrmbdhvnTGCAKYWSRMBDHVN/;
    return( $r );
   }

#
# check_frame( $f )  verifies that $f is upper or lower case F1, F2, F3, R1,
#                    R2 or R3 or a comma-separated lists of these, or "all".   
#                    If so, returns the uppercase form, otherwise dies with error message.
#                    ("ALL" is convered to "F1,F2,F3,R1,R2,R3")
#
sub  check_frame
   {
    my $f = uc( shift );
    $f = "F1,F2,F3,R1,R2,R3" if ( $f eq "ALL" );
    ($f =~ /^[FR][123](,[FR][123])*$/ ) 
      || die "bad frame specifier \"$f\"; " .
         "must be one of F1, F2, F3, R1, R2 or R3, or a comma-separated list of these\n";
    return( $f );
   }
