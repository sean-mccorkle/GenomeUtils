/* File:         trie.h                                                      */
/* Programmer:   Sean R. McCorkle                                            */
/* Language:     C                                                           */
/* Description:  Header file for trie library: data structures,              */
/*               prototypes, etc.                                            */
/*                                                                           */
/* $Id: trie.h,v 0.1 2000/08/07 14:22:12 mccorkle Exp mccorkle $  */
/*****************************************************************************/


         /**********************************************************/
         /* A trie is composed of nodes of type TR_NODE and stores */
         /* information in nodes of type TR_LEAF, which is just a  */
         /* linked list of pointers to any arbitrary record        */
         /**********************************************************/

struct tr_leaf {
                void           *data;
                struct tr_leaf *next;
               };
typedef struct tr_leaf  TR_LEAF;

struct tr_node {
                char            c;
                TR_LEAF        *leaves;
                int             count;   
                struct tr_node *children;
                struct tr_node *next;
               };
typedef struct tr_node  TR_NODE;



            /**************************************************/
            /* public use function and procedure declarations */
            /**************************************************/

void     add_trie_string( TR_NODE *t, char *s, void *data );
TR_NODE *create_trie( void );
void     DFScount( TR_NODE *t );
int      DFSdepth( void );
char    *DFSpath( void );
void     DFSprocess( TR_NODE *t, void (*proc)(TR_LEAF *) );
void     print_trie( TR_NODE *t, int indent, void (*print_data)(void *) );
void     DFSprintleaves( TR_NODE *t, void (*print_data)(void *) );
TR_NODE *search_trie( TR_NODE *t, char *s, int (*compare)( char, char ) );






