#ifndef _TYPES_H
#define _TYPES_H

#include "defs.h"

typedef struct {
  char   atom_name[5];
  char   atom_number[6];
  char   residue_name[4];
  char   residue_number[5];
  double pos[3];
  int    type;
  int    residue_type;
  int    hphil;
  int    act;
  int    class;   /*PCS-06Jun00*/
  char   icode;   /*RSK-11Nov01*/
  char   chainID; /*RSK-11Nov01*/
  double    rad;
} pdb_atom_t;

typedef pdb_atom_t  *pdb_atom_pt;

typedef struct {
  double pos[3];
  int    type;
} point_t;

typedef point_t  *point_pt;

/* data structure describing a bit string */
typedef struct {
  unsigned char      *bits;                 /* bits, stored in a byte vector */
  unsigned int       len;               /* number of bits in this bit string */
} bitstring_t;

typedef bitstring_t  *bitstring_pt;

/* atom entry in mol2 file */
typedef struct {
  int                number;              /* number of this atom in the file */
  char               name[5];                       /* atom name in the file */
  double             pos[3];                      /* coordinates of the atom */
  char               type_str[7];           /* string defining the atom type */
  int                type;                               /* type of the atom */
  int                orbit;       /* orbit or whatever follows the . in type */
  int                hyd;              /* DONOR, ACCEPTOR, DONEPTOR, NOTHING */
  bitstring_pt       fragments;             /* fragments the atom is part of */
} atom_t;

typedef atom_t  *atom_pt;

/* bond description */
typedef struct {
  int                number;        /* number of bond given in the mol2 file */
  int                atom1;              /* first atom (index in atom field) */
  int                atom2;                                   /* second atom */
  int                type;                            /* SINGLE, DOUBLE, ... */
  char               type_str[5];         /* string describing the bond type */
} bond_t;

typedef bond_t  *bond_pt;

/* structure describing a complete molecule */
typedef struct {
  char               name[MAX_LEN_MOL2_COMPOUND_NAME];      /* compound name */
  atom_pt            atoms;                           /* vector of all atoms */
  bond_pt            bonds;                           /* vector of all bonds */
  int                *flexible_bonds;  /* list of indices of rotatable bonds */
  int                **neighbors;      /* adjacency list for structure graph */
  int                **bond_to_neighbor;             /* list of bond indices */
  int                *number_of_neighbors;   /* number of neighbors of atoms */
  float              **carbon_ring_centers;     /* positions of ring centers */
  int                *carbon_ring_atom;      /* index of an atom on the ring */
  int                number_of_atoms;     /* number of atoms in the compound */
  int                number_of_bonds;               /* total number of bonds */
  int                number_of_flexible_bonds;            /* rotatable bonds */
  int                number_of_carbon_rings;       /* number of carbon rings */
  float              score;                 /* score assigned to this ligand */
} molecule_t;

typedef molecule_t  *molecule_pt;

/* entry in a rule list in the DOCK rule syntax */
typedef struct rlist {
  int                number;                      /* minimal number of atoms */
  int                type;                                  /* type of atoms */
  int                orbit;                                /* orbit of atoms */
  struct rlist       *required;                /* list of required neighbors */
  struct rlist       *prohibited;            /* list of prohibited neighbors */
} rule_list_t;

typedef rule_list_t  *rule_list_pt;

/* header of a rule list, descibes base atom, i.e. donor/acceptor or atom
   agjacent to flexible bond */
typedef struct {
  int                type;                                   /* type of atom */
  int                orbit;                                 /* orbit of atom */
  rule_list_pt       required[MAX_RULES_PER_BONDED_ATOM];  /* req. neighbors */
  rule_list_pt       prohibited[MAX_RULES_PER_BONDED_ATOM]; /* prohb. neigh. */
} rule_t;

typedef rule_t  *rule_pt;

/* definition of a flexible bond, based on DOCK syntax */
typedef struct {
  char               name[32];           /* identifier for this type of bond */
  rule_t             atom[2];                /* rules for the adjacent atoms */
} flex_bond_defn_t;

typedef flex_bond_defn_t  *flex_bond_defn_pt;

/* definition of the rules for hydrogen-bond donors and acceptors */
typedef struct {
  rule_t             donor_rules[MAX_NUMBER_OF_HYD_RULES];    /* donor rules */
  rule_t             acceptor_rules[MAX_NUMBER_OF_HYD_RULES];  /* acc. rules */
  int                number_of_donor_rules;               /* number of rules */
  int                number_of_acceptor_rules;            /* number of rules */
} hyd_defn_t;

typedef hyd_defn_t  *hyd_defn_pt;  

typedef double              transform_matrix_t[3][4];
typedef transform_matrix_t  *transform_matrix_pt;

struct _spoint
{
    point_t tpnt;
    int oindx[100];             /** list of clustered point **/
    int nindx;                  /** number of clustered points **/
    struct _spoint *next;
};

/** list of closest pair of _spoints **/
struct _cpair
{
    struct _spoint *pnts[2];
    double dist;
};


#endif
