/* File:         trie.h                                                      */
/* Programmer:   Sean R. McCorkle                                            */
/* Language:     C                                                           */
/* Description:  Header file for trie library: data structures,              */
/*               prototypes, etc.                                            */
/*                                                                           */
/* $Id$  */
/*****************************************************************************/

typedef enum { FORWARD = 'F', REVERSE = 'R' } DIRECTION;

         /**********************************************************/
         /* A trie is composed of nodes of type TR_NODE and stores */
         /* information in nodes of type TR_LEAF                   */
         /**********************************************************/

struct tr_leaf {
                int             string_id;
                int             pos;
                DIRECTION       dir;
                struct tr_leaf *next;
               };
typedef struct tr_leaf  TR_LEAF;

struct tr_node {
                char            c;
                TR_LEAF        *leaf;
                struct tr_node *children;
                struct tr_node *next;
               };
typedef struct tr_node  TR_NODE;



            /**************************************************/
            /* public use function and procedure declarations */
            /**************************************************/

void      print_trie( TR_NODE *t, int indent );
void      add_trie_string( TR_NODE *t, char *s, int string_id, 
                           int pos, DIRECTION dir );






