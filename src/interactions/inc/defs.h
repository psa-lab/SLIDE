#ifndef _DEFS_H
#define _DEFS_H

/** uncomment the next line to create SLIDE demo version **/
/**#define DEMO_VERSION**/

#define  MAX_NUMBER_OF_MOL2_ATOMS      10000
#define  MAX_NUMBER_OF_MOL2_BONDS      10000
#define  MAX_NUMBER_OF_FLEXIBLE_BONDS  500
#define  MAX_NUMBER_OF_CARBON_RINGS     50
#define  MAX_MOL2_LINELENGTH           128
#define  MAX_LEN_MOL2_COMPOUND_NAME    128
#define  MAX_RULES_PER_BONDED_ATOM       5
#define  MAX_NUMBER_OF_FLEX_BOND_RULES  40
#define  MAX_NUMBER_OF_HYD_RULES        10
#define  MAX_NEIGHBORS                   4
#define  MAX_CARBON_RING_SIZE            8
#define  MAX_NUMBER_OF_CLUSTERING_POINTS 10000

#define  ACCEPTOR    10
#define  DONOR       11
#define  DONEPTOR    12
#define  HYDROPHOB   13
#define  NOTHING     42

#define  FALSE    0
#define  TRUE     1

#define  X        0
#define  Y        1
#define  Z        2

#define  ANY     99

#define  NO_FLEX_BOND 0
#define  FLEX         1
#define  FLEX_REVERSE 2

#define  FIXED        0
#define  ROTATABLE    1

#define  FAILURE      0
#define  SUCCESS      1

#define  UNKNOWN     -1
#define  SINGLE       1
#define  DOUBLE       2
#define  TRIPLE       3
#define  QUADRUPLE    4
#define  AROMATIC     5
#define  AMIDE        6
#define  DELOCALIZED  7
#define  CYCLE_BOND  10

#define  NO_CYCLE     0
#define  CYCLE        1
#define  UNVISITED    6
#define  VISITED     20
#define  CYCLE_START 30

#define  NO_SUBSTITUENTS 0
#define  HPHOB_SUBSTITUENT 1
#define  HPHIL_SUBSTITUENT 2

#define  NO_HPHIL_NEIGHBORS 0
#define  HPHIL_NEIGHBORS 1

#define  IN       10
#define  OUT      20

#define RING_THRESHOLD 1.6
#define OVERALL_THRESHOLD 1.8
#define FINAL_THRESHOLD 0.8

#define  CA        1
#define  CB        2
#define  C         3
#define  O         4
#define  N         5
#define  CD        6
#define  CD1       7
#define  CD2       8
#define  CE        9
#define  CE1      10
#define  CE2      11
#define  CE3      12
#define  CG       13
#define  CG1      14
#define  CG2      15
#define  CH2      16
#define  CZ       17
#define  CZ2      18
#define  CZ3      19
#define  ND1      20
#define  ND2      21
#define  NE       22
#define  NE1      23
#define  NE2      24
#define  NH1      25
#define  NH2      26
#define  NZ       27
#define  OD1      28
#define  OD2      29
#define  OE1      30
#define  OE2      31
#define  OG       32
#define  OG1      33
#define  OXT      35
#define  SD       36
#define  SG       37
#define  S        40
#define  P        41
#define  H        42
#define  SI       43
#define  CO       44
#define  AG       45
#define  NI       46
#define  BR       47
#define  TL       48
#define  OS       49
#define  CL       50
#define  LA       51
#define  K        52
#define  FE       53
#define  DUMMY    54    /* dummy atom */
#define  DU       54    /* dummy atom */
#define  MO       55
#define  I        56
#define  AS       57
#define  U        58
#define  F        59
#define  RE       60
#define  RU       61
#define  ER       62
#define  LI       63
#define  MN       64
#define  PD       65
#define  CU       66
#define  AL       67
#define  MG       68
#define  NA       69
#define  CR       70
#define  CS       71
#define  TI       72
#define  ZN       73
#define  RB       74
#define  V        75
#define  RH       76
#define  GE       77
#define  ZR       78
#define  GA       79
#define  Y_       80
/***** Added by PCS -- 22-Mar-00 *****/
#define  OC       81  /* Too Close Oxygen - Very small radius */
/*************************************/

#define  SP1       1    /*  .1   */
#define  SP2       2    /*  .2   */
#define  SP3       3    /*  .3   */
/*       O         4        .o   */
#define  SP4       5    /*  .4   */
#define  AR        6    /*  .ar  */
#define  CO2       7    /*  .co2 */
#define  AM        8    /*  .am  */
#define  PL3       9    /*  .pl3 */
#define  CAT      10    /*  .cat */
#define  OH       11    /*  .oh  */
#define  O2       12    /*  .o2  */
#define  TH       13    /*  .th  */
#endif

#define MIN_DOUBLE 0.00000001
