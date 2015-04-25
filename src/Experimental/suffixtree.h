/* File:         suffixtree.h                                                */
/* Programmer:   Sean R. McCorkle                                            */
/* Language:     C                                                           */
/* Description:  Header file for suffixtree library: data structures,        */
/*               prototypes, etc.                                            */
/*                                                                           */
/* $Id$  */
/*****************************************************************************/

          /******************************************************/
          /* A suffix tree is composed of nodes of type ST_NODE */
          /******************************************************/

struct st_node {
                char   c;
                int    loc;
                struct st_node *children;
                struct st_node *next;
               };
typedef struct st_node  ST_NODE;


               /********************************************/
               /* function and procedure type declarations */
               /********************************************/


ST_NODE *new_st_node( char c, int loc, ST_NODE *children, ST_NODE *next );
ST_NODE *linear_suffix_tree( char *s, int loc );
ST_NODE *search_st_list( ST_NODE *l, char c );
void     add_to_suffix_tree( ST_NODE *t, char *s, int loc );
ST_NODE *create_suffix_tree( char *s );
void     print_suffix_tree( ST_NODE *t, int indent );
void     enumerate_st_hits( ST_NODE *t );
void     search_suffix_tree( ST_NODE *t, char *s );





