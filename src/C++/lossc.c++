// Program:      lossc
// Programmers:  Lawrence Brenner, initial development Aug, 2001.
//               updated by Sean McCorkle, Mar 2002
// Language:     C++, with STL
//
// Description:  Lots of short sequence comparisons.
//
// Usage:        lossc [-hmV] [-v<n>] <thresh> <data seqs> [<db seqs>]
//  
//               <thresh> is an integer indicated the maximum edit distance
//                        considered for a match 
//               <data seqs> is a file of short sequences and descriptions
//               <db seqs>  is a file of short sequences and descriptions
//
//                 -m  show misses as well as hits
//                 -v<n>  verbosity level n (1, 2, ...) higher means more mess.
//                 -h     print usage message, then exit
//                 -V     print version, then exit
//
// $Id: lossc.c++,v 0.5 2002/05/02 21:28:45 mccorkle Exp mccorkle $
//
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <stl.h>
#include <stdlib.h>
#ifndef SOLARIS        /* Sun's are really a pain in the neck sometimes */
#include <unistd.h>
#endif


static char const rcsid[] = "$Id: lossc.c++,v 0.5 2002/05/02 21:28:45 mccorkle Exp mccorkle $";

#define  LARGE 1000000

const int  max_num_seqs = 500000;    // maximum number of sequences
const int  max_seq_size  = 100;      // max length of the sequences is maxtrix
                                     // size.
int        min_seq_length;           // set by main()
int        misses;                   // set by -m option in GetArgs
int        verbosity;                // set by -v option in GetArgs

static  int   EdDist[max_seq_size][max_seq_size]; //The Smith-Waterman matrix

//
// usage() - print usage message and bomb out
//
void usage( void )
   {
    cout << "Usage:  lossc [-hmV] [-v<n>] <thresh> <seqfile> [<seqfile2>]\n";
    exit( 1 );
   }

//
// spacer( c, n ) - return a string of n c's
//
string spacer( char c, int n )
   {
    string b( n, c );
    return( b );
   }

//
// GetArgs() - parse args for options and set global values accordingly,
// and then extract edit_threshold and filenames and return them.
//

void GetArgs( int     argc,         // input- argc from main( argc, argv )
              char   *argv[],       // input- argv from main( argc, argv )
              int    &edit_thresh,  // output- edit distance threshold
              string &file1,        // output- mandatory filename
              string &file2         // output- optional 2nd file or "" if none
            )
   {
    extern char *optarg;
    extern int   optind;
    int          c;

    misses = 0;
    verbosity = 0;

    while ( (c = getopt( argc, argv, "hmv:V")) != -1 )
        switch ( c )
           {
            case 'm':  misses = 1;
                       break;
            case 'v':  verbosity = atoi( optarg );
                       break;
            case 'V':  cout << "lossc, $Revision: 0.5 $\n";
                       break;
            case 'h':  
            default:
                       usage();
           }
    argc -= optind;
    argv += optind;

    if ( argc < 2 )
        usage();

    edit_thresh = atoi(argv[0]);
    file1 = argv[1];
    if ( argc > 2 )
        file2 = argv[2];
    else
        file2 = "";
   }


// note to self: put in lots of input checks, including max seq length
// check.

void  Load_Sequences(  string         filename,   // input- short seq filename
                       vector<string> &sequences,  // output- array of seqs
                       vector<string> &descrip,    // output- corres. descrips
                       int            &count       // 
                    )
 {
  ifstream           in; 
  string             line;
  int                i = 0;
  int                k;

  in.open( filename.c_str() );
  if ( ! in )
    {
      cerr << "Can not open "<< filename.c_str() << "\n";
      exit;
    }

  while ( ! in.eof() )
    {
      getline(in, line);
      k = line.find( ' ' );
      if ( verbosity > 1 )
          cout << "k " << k << " line: [" << line << "] " 
               << line.length() << "\n";
      if ( k > 0 )
        {
         sequences[i] = line.substr( 0, k );
         descrip[i] = line.substr( k+1, line.length() );
         if ( verbosity > 1 )
             cout << "input: [" << sequences[i] << "] [" << descrip[i] 
                  << "]\n";
         i++;
	}
      else if ( line.length() > 0 )
        {
         sequences[i] = line;
         descrip[i] = ' ';
         if ( verbosity > 1 )
             cout << "input: [" << sequences[i] << "] [" << descrip[i] 
                  << "]\n";
         i++;
        }
    }
    in.close();
    count = i; //i - 1; - why was this here?
    if ( verbosity == 1 )
        cout << count << " sequences in file " << filename << "\n";

    sequences.resize( count );
    descrip.resize( count );
} 

//
// min_length() returns the minium length of the seqs in the vector
//
int   min_length( vector<string> seqs )
   {
    vector<string>::iterator s;
    int                      min = max_seq_size + 1;

    for ( s = seqs.begin();  s != seqs.end();  s++ )
        if ( (*s).length() < min )
            min = (*s).length();
    if ( verbosity > 0 )
        cout << "Min length " << min << "\n";
    return( min );
   }

//
// trunc_lengths() truncates all seqs down to size len
//
void  trunc_lengths( vector<string> &seqs,  int len )
   {
    vector<string>::iterator s;
    
    for ( s = seqs.begin();  s != seqs.end();  s++ )
        if ( (*s).length() > len )
           {
            if ( verbosity > 1 )
                cout << "    truncating " << *s;
	    (*s).erase( ((*s).begin() + len), (*s).end() );
            if ( verbosity > 1 )
                cout << " to " << *s << "\n";
           }
   }

               ///////////////////////////////
               // Edit distance calculation //
               ///////////////////////////////
//
// Fills the very top row and left column with increasing scores
//
void  Initialize_Matrix( void )
   {
    int row;
    int col;

    for ( row = 0; row < max_seq_size; row++ ) //Sets the values of r,0 
        EdDist[row][0] = row;

    for ( col = 1; col < max_seq_size; col++ ) //Sets the values of 0,c 
        EdDist[0][col] = col;
   }


//
// UpdateMin() - selects the minimum score at EdDist[row,col], which may 
// necessitate checking seq1[row-1] and seq2[col-1].  Returns the minimum
// value in final
//
void  UpdateMin( int     row,   // input- row number for evaluation
                 int     col,   // input- col number for evaluation
                 string  seq1,  // input- 1st short sequence
                 string  seq2,  // input- 2nd short sequence
                 int    &final    // output- minimum value
               )
   {
    int           b = 0;          //Used in the substitution, match function.
  

    if ( seq1[row - 1] == seq2[col - 1] ) //Checks for matches.
        b = 0;
    else                                     // Checks for substituitons.
        b = 1;
  
    final = EdDist[row - 1][col - 1] + b;    // - 1;
  
    if ( EdDist[row - 1][col] + 1 < final )  // Checks for insertions.
        final = EdDist[row - 1][col] + 1;
    if ( EdDist[row][col - 1]  + 1 < final ) // Checks for deletions.
        final  = EdDist[row][col - 1] + 1; 

    EdDist[row][col] = final; //++final;
   }

// This version of Edit_Distance counts on the previous values of
// EdDist[r][c] from the previous sequence comparsion to be present,
// up to the subsquare determined by start_diag.  

// returns:  dist - edit distance (or threshold + 1 if aborted)
//           last_diag - matrix diag indicates where calculation stopped 
//                       (< seq.length() if aborted early because threshold
//                        exceeded)

void Edit_Distance( string seq1,        // input- first short sequence 
		    string seq2,        // input- 2nd short sequence 
		    int    threshold,   // input- abort if dist > threshold
		    int    start_diag,  // input- recalculate mat @ this column
		    int    &last_diag,  // output- last diag calc. before abort
		    int    &dist )      // output- resultant distance
   {
    int  eta = 0;                // Minimum edit distance changes.
    int  min;                    // The lowest calculated value at one position
    int  r, c;                   // row and column counters
 
    int  rows = seq1.length();   // number of rows in the matrix
    int  cols = seq2.length();   // number of columns
    int  diag;                   // Variable that analyzes only the revelant diagonals.

    // Begin filling the matrix EdDist[r][c] in subsquares, starting
    // in the upper right corner, or actually at the specified start position

    for ( diag = start_diag; ((diag <= rows) && (eta <= threshold)); diag++ )
       {
        eta = EdDist[0][diag];
        for (  r = 1;  r < diag;  r++  )         // do the column
	   {
	    UpdateMin( r, diag, seq1, seq2, min );
	    if ( min < eta )
	        eta = min;
           }

        for (  c = 1;  c <= diag;  c++  )        // now do the row
	   {
	    UpdateMin( diag, c, seq1, seq2, min );
            if ( min < eta )
               eta = min;
           }
       }
 
    if ( eta > threshold )
        dist = threshold + 1;    
    else
        dist = EdDist[rows][cols];

    last_diag = diag;
   }


//
// first_diff( a, b ) returns the index of the first position where
// strings a & b differ.  (Note: we assume both strings are the same length 
// - is this always going to be true?
//
int  first_diff( string a, string b )
   {
    int i = 0; 

    while ( ( i < a.length() ) && ( a[i] == b[i] ) )  
        i++;
    return i;
   }


void  Cross_Compare( int             edit_thresh,
                     vector <string> tags, 
                     vector <string> tag_descrip, 
                     int             tcount, 
                     vector <string> database, 
                     vector <string> db_descrip, 
                     int             dcount
                   )
   {
    int             t;              // database sequence index
    int             d;              // counter for the database.
    int             ed;             // edit distance between the sequences.
    int             same;           // difference between the two sequences.
    int             last = LARGE;  // last col position calculated in matrix.
    int             best_ed;

    for ( t = 0; t < tcount; t++ ) //Checks the tags to the database.
       {
        best_ed = edit_thresh + 1;
        for ( d = 0; d < dcount; d++ )
           { 
            if ( d == 0 )
                same = 0;
            else
                same = first_diff( database[d - 1], database[d] ); 

            // if we aborted early the previous time, and this sequence
            // is the same BEYOND that point, then there's no need to 
            // to bother - this one doesn't match in the threshold either.

            if ( same < last )
               {
                Edit_Distance( tags[t], database[d], edit_thresh, same, 
                               last, ed );
                if ( ed <= edit_thresh )
                   {
                    cout << tags[t] << " " << database[d] << " " << ed <<
                         " " << tag_descrip[t] << "|" << db_descrip[d] << "\n";
                    best_ed = ed;
                   }
               }
           }
        if ( misses && ( best_ed > edit_thresh ) )  // no hits, print miss output
            cout << tags[t] << " " << spacer( '-', min_seq_length )
                 << " X " <<tag_descrip[t]<<"\n";

       }   
   }


void  Self_Compare( int             edit_thresh,
                    vector <string> tags, 
                    vector <string> tag_descrip, 
                    int             tcount )
   {
    int             i, j;
    int             same;           // difference between the two sequences.
    int             last = LARGE;  // last col position calculated in matrix.
    int             ed;
    int             best_ed;

    if ( verbosity > 1 )
        cout << "Self check\n";

    for ( i = 0; i < tcount; i++ )
       {
        best_ed = edit_thresh + 1;
        for ( j = i+1; j < tcount; j++ )
           {
            if ( j == i + 1 )
                same = 0;
            else
                same = first_diff( tags[j - 1], tags[j] ); 

            // if we aborted early the previous time, and this sequence
            // is the same BEYOND that point, then there's no need to 
            // to bother - this one doesn't match in the threshold either.

            if ( same < last )
               {
                Edit_Distance( tags[i], tags[j], edit_thresh, same, 
                               last, ed );
                if ( ed <= edit_thresh )
                   {
                    cout << tags[i] << " " << tags[j] << " " << ed << " " <<
                         tag_descrip[i] << "|" << tag_descrip[j] << "\n";
                    best_ed = ed;
                   }
               }
           }
        if ( misses && ( best_ed > edit_thresh ) )  // no hits, print miss output
            cout << tags[i] << " " << spacer( '-', min_seq_length ) 
                 << " X " <<tag_descrip[i]<<"\n";
       }
   }

                           //////////////////
                           // Main Program //
                           //////////////////


int  main( int argc, char *argv[] )
  
   {
    int             dcount; //The counter for the database.
    int             tcount; //The counter for the tags.
    int             edit_thresh; //How much of an edit distance you want.
    int             minl2;
    vector <string> tags(max_num_seqs); //The vector of tags that are to be checked against the databse.
    vector <string> database(max_num_seqs); //The vector of sequences that are in the databse.
    vector <string> tag_descrip(max_num_seqs);
    vector <string> db_descrip(max_num_seqs);
    string          file1;
    string          file2;

    GetArgs( argc, argv, edit_thresh, file1, file2 );

    Initialize_Matrix();

    Load_Sequences( file1, tags, tag_descrip, tcount );    // read data tags
    min_seq_length = min_length( tags );
    if ( file2 != "" ) // read database tags
       {
        Load_Sequences( file2, database, db_descrip, dcount ); 
        minl2 = min_length( database );
        if ( minl2 < min_seq_length ) 
            min_seq_length = minl2;
        trunc_lengths( tags, min_seq_length );        
        trunc_lengths( database, min_seq_length );        
        sort( database.begin(), database.end() );  // sort alphabetically
        Cross_Compare( edit_thresh, tags, tag_descrip, tcount, 
                     database, db_descrip, dcount );
       }
    else
       {
        trunc_lengths( tags, min_seq_length );        
        sort( tags.begin(), tags.end() );
        Self_Compare( edit_thresh, tags, tag_descrip, tcount );
       }

   }



