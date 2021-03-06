#!/usr/bin/perl
#
# Program:      orfs.pl
# Programmer:   Sean R. McCorkle
# Language:     perl
# Description:  Read DNA sequences and report ORFs in all reading frames, 
#               both directions
#
# Usage:        orfs.pl [-l<n>] [-W] [-s cod[,...]]  [<seqfile> ...]
#
#               where <seq-file>'s are in FASTA format.  Stdin is scanned
#               if no files are specified.
#
# Options:
#               -l <n>        report ORF only if length meets or exceeds length
#                             of <n> (default 300)
#               -s cod[,...]  use specified codons for starts, rather than 
#                             ATG, GTG and TTG (comma-separated)
#               -W            assume sequence is a subsequence - report 
#                             possible ORFs on 5' ends without a start (assume 
#                             a start occured past the end of the sequence), 
#                             and reportm posssible orfs on 3' ends without 
#                             stops.
#
# Output:
#               fasta DNA sequences for each orf.  Input sequence Header 
#               is prepended by a string describing the orf, 
#
#                  <frome>-<start pos>-<length nt>-<start codon>-<stop codon>
#
#               where frame is one of F1, F2, F3, R1, R2, or R3
#               for example
#
#                       F2-33497-678-TTG-TGA
#               or
#                       R3-4192-1029-ATG-TGA
#                           
# TODO:   
#               . add options to control output format.  maybe table output?
#
# $Id: orfs.pl,v 1.4 2009/09/10 12:51:22 mccorkle Exp mccorkle $
#

use  Getopt::Std;

getopts( "l:s:W");

                               ###########
                               # Globals #
                               ###########

$nt_thresh = $opt_l ? $opt_l : 300;     # report ORF only if this length 
                                        # (in nucleotides) is exceeded

@starts = $opt_s ? check_codons( $opt_s ) : ( "ATG", "GTG", "TTG" );
@stops  = ( "TAA", "TAG", "TGA" );

map {$signal{$_} = 'start';} @starts;   # make a single codon lookup table
map {$signal{$_} = 'stop';}  @stops;    # that reports 'start' or 'stop'

$line_len = 50;                         # for fasta output


                            ################
                            # Main Program #
                            ################

while ( <> )
   {
    if ( /^>(.*)/ )
       {
        $new_hdr = $1;
        report_orfs( $seq, $hdr ) if ( defined( $hdr ) );
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

report_orfs( $seq, $hdr ) if ( defined( $hdr ) );

                             ###############
                             # Subroutines #
                             ###############


#
# find and report all ORFs in $seq, all reading frames, both directions
#
sub  report_orfs
   {
    my ( $seq, $hdr ) = @_;
    
    for ( my $i = 0; $i < 3; $i++ )
        { scan_frame( $seq, $hdr, $i, "F" ); }
    my $rseq = rc( $seq );
    for ( my $i = 0; $i < 3; $i++ )
        { scan_frame( $rseq, $hdr, $i, "R" ); }
   }

#
#  looks for occurances of a start followed by a stop in one direction
#  ($dir == "F" or "R") and one reading frame ($off = 0, 1, 2) only
#
sub  scan_frame
   {
    my ( $seq, $hdr, $off, $dir ) = @_;

    my $len = length( $seq );
    my ($state,$startpos) = $opt_W ? ("inside",$off) : ("outside", -1000000);

    my $j;
    for ( $j = $off; $j + 3 <= $len; $j += 3 )
       {
        $codon = uc( substr( $seq, $j, 3 ) );
        $sig = $signal{$codon};
        if ( $sig eq 'start' )
            {
             if ( $state eq "outside" )
                {
                 $startpos = $j;
                 $state = "inside";
                }
            }
        elsif ( $sig eq 'stop' )
            {
             if ( $state eq "inside" )
                {
                 $state = "outside";
                 output_orf( $seq, $hdr, $off, $dir, $startpos, $j );
                }
            }
       }
    if ( $opt_W && $state eq "inside" )
       {
        $state = "outside";
        output_orf( $seq, $hdr, $off, $dir, $startpos, $j );
       }
   }


sub  output_orf
   {
    my ( $seq, $hdr, $off, $dir, $a, $b ) = @_;
    my $len_nt = $b - $a;
    my $orfhdr;

    if ( $len_nt >= $nt_thresh )
       {
        if ( $opt_W && $a < 3 )
           { $start_codon = "???"; }
        else
           { $start_codon = substr( $seq, $a, 3 ); }

        if ( $opt_W && $b >= length( $seq ) - 3 )
           { $stop_codon = "???";  }
        else
           { $stop_codon = substr( $seq, $b, 3 ); }

        if ( $dir eq "F" )
           {
            $orfhdr = "$dir" . ($off+1) . "-" . ($a+1) . 
                        "-$len_nt-$start_codon-$stop_codon";
           }
        else
           {
            my $s = length( $seq ) - ($a + $len_nt);
            $orfhdr = "$dir" . ($off+1) . "-" . ($s+1) . 
                        "-$len_nt-$start_codon-$stop_codon";
           }
        $orfhdr .= " $hdr";
        print_fasta( $orfhdr, substr( $seq, $a, $len_nt ) );
       }
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
# check_codons( $s ) - verifies that $s is a comma-separted list of DNA 
#                      triplets, and then returns them as a list.  Dies with 
#                      error message if check fails.
#
sub  check_codons
   {
    my $s = shift;
    ($s =~ /^[ACGT]{3}(,[ACGT]{3})*$/) ||
       die "bad codon list \"$s\"; " . 
            "must be comma-separated list of DNA triplets, ala ATG,TTG,GTG\n";
    return( split( /,/, $s ) );
   }
