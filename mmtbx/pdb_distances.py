import sys
import numpy


##########################################################################
##########################################################################
##########################################################################
#SECTION: IMPORTANT NOTES                                 BEGINNING
##########################################################################
##########################################################################
##########################################################################

#########################################
#Syntax: "phenix.python pdb_distances_######.py pdb-file.pdb > something.log" (where ###### is a date when the program was last edited)
#########################################

###########################
#######040610 PROGRAM FLOW:
        #"MASTER_Basepairs_bonds" is a list that carries pre-compiled basepair information regarding a denomination for the type of basepair interaction, residues invoved, number of bonds, and atoms involved in the bonds
                #Loaded in "SECTION: LIST DEFINITION"
                #"MASTER_Basepairs_bonds" contains as many lines as basepairs
        #FUNCTION "run(args)" is the program created by Ralf W. Grosse-Kunstle and provides a list of ATOM-to-ATOM distances calculated from a .pdb file. The rest of the program was created around UNCTION "run(args)"
                #Called from "Sub SECTION Calling run(args)" within "SECTION MAIN"
        #"First_List" Contains the filtered output from FUNCTION "run(args)". Complete list of lines carrying ATOM to ATOM distances below a rough cutoff of 5A.
                #The line "pair_asu_table = xray_structure.pair_asu_table(    distance_cutoff=5.0)" in FUNCTION "run(args)" carries the cutoff value used in FUNCTION "run(args)"

        #"MASTER_Basepairs" receives from "First_List" all the lines whose ATOM to ATOM distances are consistent with hydrogen bonds in basepair interactions, as defined by "MASTER_Basepairs_bonds".
                #The transfer of information from "First_List" to "MASTER_Basepairs" is performed in "Section: First sorting of basepair candidates" of FUNCTION "program".
                        #FUNCTION "program" is called from "Sub SECTION MAIN LOOP"
                #"MASTER_Basepairs" contains as many first-level sublists as lines in "MASTER_Basepairs_bonds" (see below). Therefore, all the information in every second-level sublist is consistent with the geometry of the basepair specified by the equivalent first-level sublist in "MASTER_Basepairs_bonds"
                #The information in each second-level sublist of "MASTER_Basepairs" merely involves pairs of ATOMs whose distance to one another is below a new CUTOFF value, dynamically assigned by the program at every run of "Sub SECTION MAIN LOOP". Therefore, in order to assign a 2-bond basepair at a particular CUTOFF value, there would have to be 2 second-level sublists somewhere in the "MASTER_Basepairs" first-level sublist corresponding to that basepair. Three lines would then be required for a three-bond basepair, although the assignment of these basepairs with only two bonds is allowed
        #"MASTER_Basepairs_summary" has the same number of first-level sublists as lines in "MASTER_Basepairs_bonds" and first-level sublists in "MASTER_Basepairs". "MASTER_Basepairs_summary" receives the content of "MASTER_Basepairs" and places it into a new set of second level sublists, each of which corresponds to  a single basepair candidate. Information for all possible bonds found relevant to the pair of residues in the basepair candidate are added onto these second-level sublists
                #All operations for the assigmnent of basepairs are carried out with "MASTER_Basepairs_summary"
                #Many second-level sublists in "MASTER_Basepairs_summary" will not be assigned to any basepairs
        #"new_list_end" will be used to filter out the information of "MASTER_Basepairs_summary" into a single-level list carrying only assigned basepairs. These basepairs will be ordered by RESIDUE (from smallest res# to largest res#)
                #Loaded in "SECTION ORDERING OUTPUT"
        #"run_cutoff_LISTS" will carry the same information as "new_list_end" but ordered by CUTOFF at which the basepair was identified. First level lists are the CUTOFF levels, according to the values stored in "run_cutoff". "run_cutoff_LISTS" will be used for STATISTICS
                #Loaded in SECTION STATISTICS
                        #"run_cutoff" is loaded in "Sub SECTION MAIN LOOP"
#######040610 PROGRAM FLOW:
###########################

#Note 012110:
  #By testing this script with different cutoffs for detection of basepairs, (see 'CUTOFF'), I have realized that different subsets of basepairs are identified when the resolution is changed from 3.2 to 3.6. At low cuttoff, only high-confidence basepairs are identified and very few of them get tagged as involved in MULTIPLE CONTACTS. As the cutoff is increased, the number of basepairs involved in MULTIPLE CONTACTS increases and so does the number of misidentified basepairs. This suggests that the best way to maximize the number of identified basepairs should be to start the run at low cuttoff, let's say 3.2 and perform successive runs at increasing cuttoffs, until 3.6 is reached. If the basepairs identified with high confidence are fixed, new basepairs can be identified with confidence at higher cutoffs by extension of preformed helices
##########################################################################
##########################################################################
##########################################################################
#SECTION: IMPORTANT NOTES                                 END
##########################################################################
##########################################################################
##########################################################################


##########################################################################
##########################################################################
##########################################################################
#SECTION: LIST DEFINITION                                 BEGINNING
##########################################################################
##########################################################################
##########################################################################
#"First_List" Contains the filtered output from FUNCTION "run(args)". Listof ATOM to ATOM distances below a rough cutoff of 5A

#"MASTER_Basepairs" master list: This list will contain lines from "First_List" that carry ATOM to ATOM distances that are consistent with basepair geometry, as established in "Basepair Lists SECTION" below.
        #"MASTER_Basepairs" will have as many sublists as established in "Basepair Lists SECTION" below. The basepair denomination for these sublists is (NOTE that this denomination will be carried separately in "MASTER_Basepairs_schemes" and also in "MASTER_Basepairs_summary[i][0]").
        #"MASTER_Basepairs" will have as many sub-sublists as basepair candidates are identified for a particular basepair scheme
###i = 0 in MASTER_Basepairs (I_AA)
###i = 1 in MASTER_Basepairs (II_AA)
###i = 2 in MASTER_Basepairs (III_GG)
###i = 3 in MASTER_Basepairs (IV_GG)
###i = 4 MASTER_Basepairs (V_AA)
###i = 5 in MASTER_Basepairs (VI_GG)
###i = 6 in MASTER_Basepairs (VII_GG)
###i = 7 in MASTER_Basepairs (VIII_AG)
###i = 8 in MASTER_Basepairs (VIII_GA)
###i = 9 in MASTER_Basepairs (IX_AG)
###i = 10 in MASTER_Basepairs (IX_GA)
###i = 11 in MASTER_Basepairs (X_AG)
###i = 12 in MASTER_Basepairs (X_GA)
###i = 13 in MASTER_Basepairs (XI_AG)
###i = 14 in MASTER_Basepairs (XI_GA)
###i = 15 in MASTER_Basepairs (XII_UU)
###i = 16 in MASTER_Basepairs (XIII_UU)
###i = 17 in MASTER_Basepairs (XIV_CC)
###i = 18 in MASTER_Basepairs (XV_CC)
###i = 19 in MASTER_Basepairs (XVII_CU)
###i = 20 in MASTER_Basepairs (XVII_UC)
###i = 21 in MASTER_Basepairs (XVIII_CU)
###i = 22 in MASTER_Basepairs (XVIII_UC)
###i = 23 in MASTER_Basepairs (XIX_CG_WC)
###i = 24 in MASTER_Basepairs (XIX_GC_WC)
###i = 25 in MASTER_Basepairs (XX_AU_WC)
###i = 26 in MASTER_Basepairs (XX_UA_WC)
###i = 27 in MASTER_Basepairs (XXI_AU)
###i = 28 in MASTER_Basepairs (XXI_UA)
###i = 29 in MASTER_Basepairs (XXII_CG)
###i = 30 in MASTER_Basepairs (XXII_GC)
###i = 31 in MASTER_Basepairs (XXIII_AU)
###i = 32 in MASTER_Basepairs (XXIII_UA)
###i = 33 in MASTER_Basepairs (XXIV_AU)
###i = 34 in MASTER_Basepairs (XXIV_UA)
###i = 35 in MASTER_Basepairs (XXV_AC)
###i = 36 in MASTER_Basepairs (XXV_CA)
###i = 37 in MASTER_Basepairs (XXVI_AC)
###i = 38 in MASTER_Basepairs (XXVI_CA)
###i = 39 in MASTER_Basepairs (XXVII_GU)
###i = 40 in MASTER_Basepairs (XXVII_UG)
###i = 41 in MASTER_Basepairs (XXVIII_GU)
###i = 42 in MASTER_Basepairs (XXVIII_UG)
###i = 44 in MASTER_Basepairs (XXX_CA)

MASTER_Basepairs = []
MASTER_Basepairs_excluded = [] #Will carry all the lines excluded from MASTER_Basepairs because they exceed CUTOFF value. Will be important while attempting to recover missing bonds. Will be important while attempting to recover missing bonds in "SECTION STATISTICS"
#THE BASEPAIR LISTS SECTION: The individual basepair lists in THE BASEPAIR LISTS SECTION start with a roman numeral that defines the basepair according to Saenger.
        #For basepairs formed with bases of different identitiy, two lists are created, depending on what the first base of the pair is (the one with the lower resid # ----> atom_i.resid()).
        #For every basepair in , a short list will also be created that contains:
                #1) number of bonds ---> n,
                #2) resid i (atom_i.resid())
                #3) resid j (atom_j.resid())
                #5) bond atom 1 for bond 1,
                #6) bond atom 2 for bond 1, ...
                #x-1) bond atom 1 for bond n,
                #x) bond atom 2 for bond n.
                #These lists will be appended to MASTER_Basepairs_bonds, which should also have a length of 43 like MASTER_Basepairs. Therefore, the relative positions for all basepair searching parameters in MASTER_Basepairs_bonds can be easily accessed by using their position in MASTER_Basepairs

MASTER_Basepairs_bonds = []


########################### Basepair Lists SECTION #################################
##### bonds[0] = # of bonds for current geometry
##### bonds[1] = Resid name for base 1
##### bonds[2] = Resid name for base 2
##### bonds[3] = atom name in bond #1 for base #1
##### bonds[4] = atom name in bond #1 for base #2
##### bonds[5] = atom name in bond #2 for base #1
##### bonds[6] = atom name in bond #2 for base #2
##### bonds[7] = atom name in bond #3 for base #1, NA if only 2 bonds for current geometry
##### bonds[8] = atom name in bond #3 for base #2, NA if only 2 bonds for current geometry
##### bonds[9] = average P-P distance for current geometry, NA if not determined
##### bonds[10] = standard deviation P-P distance for current geometry, NA if not determined
##### bonds[11] = average C1-C1 distance for current geometry, NA if not determined
##### bonds[12] = standard deviation C1-C1 distance for current geometry, NA if not determined


##### Homo purine
###i = 0 in MASTER_Basepairs (I_AA)
I_AA = []
#I Homo purine, Base-pairing pattern AA: AA_2
#Bond   A       A       Length Ave      Length Std      Attribute
#1      N1      N6      2.92    0.14    T
#2      N6      N1      3.06    0.14    T
bonds = [2, "A", "A", "N1", "N6", "N6", "N1", "NA", "NA", "NA", "NA", "NA", "NA"]
MASTER_Basepairs_bonds.append(bonds)
###i = 1 in MASTER_Basepairs (II_AA)
II_AA = []
#II Base-pairing pattern AA: AA_9
#Bond   A       A       Length Ave      Length Std      Attribute
#1      N6      N7      3.01    0.20    T
#2      N7      N6      2.89    0.21    T
bonds = [2, "A", "A", "N6", "N7", "N7", "N6", "NA", "NA", "NA", "NA", "NA", "NA"]
MASTER_Basepairs_bonds.append(bonds)
###i = 2 in MASTER_Basepairs (III_GG)
III_GG = []
#III Base-pairing pattern GG: GG_90
#Bond   G       G       Length Ave      Length Std      Attribute
#1      N1      O6      2.85    0.12    T
#2      O6      N1      2.88    0.12    T
bonds = [2, "G", "G", "N1", "O6", "O6", "N1", "NA", "NA", "NA", "NA", "NA", "NA"]
MASTER_Basepairs_bonds.append(bonds)
###i = 3 in MASTER_Basepairs (IV_GG)
IV_GG = []
#IV Base-pairing pattern GG: GG_109
#Bond   G       G       Length Ave      Length Std      Attribute
#1      N2      N3      2.96    0.29    T
#2      N3      N2      3.23    0.15    T
bonds = [2, "G", "G", "N2", "N3", "N3", "N2", "NA", "NA", "NA", "NA", "NA", "NA"]
MASTER_Basepairs_bonds.append(bonds)
###i = 4 MASTER_Basepairs (V_AA)
V_AA = [] #V Base-pairing pattern AA: AA_4
#Bond   A       A       Length Ave      Length Std      Attribute
#1      N1      N6      2.93    0.18    T
#2      N6      N7      3.07    0.17    T
#STATISTICAL DATA obtained from 2J02.pdb, 5 basepairs
bonds = [2, "A", "A", "N1", "N6", "N6", "N7", "NA", "NA", 12.63, 0.964, 12.20, 0.220]
MASTER_Basepairs_bonds.append(bonds)
###i = 5 in MASTER_Basepairs (VI_GG)
VI_GG = [] #VI Base-pairing pattern GG: GG_16
#Bond   G       G       Length Ave      Length Std      Attribute
#1      N1      O6      2.88    0.16    T
#2      N2      N7      2.91    0.15    T
bonds = [2, "G", "G", "N1", "O6", "N2", "N7", "NA", "NA", "NA", "NA", "NA", "NA"]
MASTER_Basepairs_bonds.append(bonds)
###i = 6 in MASTER_Basepairs (VII_GG)
VII_GG = [] #VII Base-pairing pattern GG: GG_21
#Bond   G       G       Length Ave      Length Std      Attribute
#1      N1      N7      2.93    0.17    T
#2      N2      O6      2.76    0.24    T
bonds = [2, "G", "G", "N1", "N7", "N2", "O6", "NA", "NA", "NA", "NA", "NA", "NA"]
MASTER_Basepairs_bonds.append(bonds)

#####Hetero Purine
###i = 7 in MASTER_Basepairs (VIII_AG)
VIII_AG = [] #VIII (AG Imino). Base-pairing pattern AG: AG_3
#Bond   A       G       Length Ave      Length Std      Attribute
#1      N1      N1      2.88    0.15    T
#2      N6      O6      2.95    0.21    TI
#STATISTICAL DATA obtained from 2J02.pdb, 11 basepairs. Combined VIII_AG and VIII_GA. Needs to be revised
bonds = [2, "A", "G", "N1", "N1", "N6", "O6", "NA", "NA", 19.598, 1.085, 12.795, 0.222]
MASTER_Basepairs_bonds.append(bonds)
###i = 8 in MASTER_Basepairs (VIII_GA)
VIII_GA = [] #VIII (GA Imino). Base-pairing pattern GA: GA_3
#Bond   G       A       Length Ave      Length Std      Attribute
#1      N1      N1      2.88    0.15    T
#2      O6      N6      2.95    0.21    T
#STATISTICAL DATA obtained from 2J02.pdb, 11 basepairs. Combined VIII_AG and VIII_GA. Needs to be revised
bonds = [2, "G", "A", "N1", "N1", "O6", "N6", "NA", "NA", 19.598, 1.085, 12.795, 0.222]
MASTER_Basepairs_bonds.append(bonds)
###i = 9 in MASTER_Basepairs (IX_AG)
IX_AG = [] #IX Base-pairing pattern AG: AG_81
#Bond   A       G       Length Ave      Length Std      Attribute
#1      N7      N1      3.08    0.08    T
#2      N6      O6      2.70    0.09    T
bonds = [2, "A", "G", "N7", "N1", "N6", "O6", "NA", "NA", "NA", "NA", "NA", "NA"]
MASTER_Basepairs_bonds.append(bonds)
###i = 10 in MASTER_Basepairs (IX_GA)
IX_GA = [] #IX Base-pairing pattern GA: GA_81
#Bond   G       A       Length Ave      Length Std      Attribute
#1      N1      N7      3.08    0.08    T
#2      O6      N6      2.70    0.09    T
bonds = [2, "G", "A", "N1", "N7", "O6", "N6", "NA", "NA", "NA", "NA", "NA", "NA"]
MASTER_Basepairs_bonds.append(bonds)
###i = 11 in MASTER_Basepairs (X_AG)
X_AG = [] #X Base-pairing pattern AG: AG_38
#Bond   A       G       Length Ave      Length Std      Attribute
#1      N1      N2      3.02    0.19    T
#2      N6      N3      3.15    0.11    T
bonds = [2, "A", "G", "N1", "N2", "N6", "N3", "NA", "NA", "NA", "NA", "NA", "NA"]
MASTER_Basepairs_bonds.append(bonds)
###i = 12 in MASTER_Basepairs (X_GA)
X_GA = [] #X Base-pairing pattern GA: GA_38
#Bond   G       A       Length Ave      Length Std      Attribute
#1      N2      N1      3.02    0.19    T
#2      N3      N6      3.15    0.11    T
bonds = [2, "G", "A", "N2", "N1", "N3", "N6", "NA", "NA", "NA", "NA", "NA", "NA"]
MASTER_Basepairs_bonds.append(bonds)
###i = 13 in MASTER_Basepairs (XI_AG)
XI_AG = [] #XI (AG Sheared). Base-pairing pattern AG: AG_24
#Bond   A       G       Length Ave      Length Std      Attribute
#1      N7      N2      3.075           0.236           T
#2      N6      N3      3.132           0.227           T
#STATISTICAL DATA obtained from 2J02.pdb, 10 basepairs
bonds = [2, "A", "G", "N7", "N2", "N6", "N3", "NA", "NA", 16.04, 1.262, 9.512, 0.399]
MASTER_Basepairs_bonds.append(bonds)
###i = 14 in MASTER_Basepairs (XI_GA)
XI_GA = [] #XI (GA Sheared). Base-pairing pattern GA: GA_24
#Bond   G       A       Length Ave      Length Std      Attribute
#1      N2      N7      2.966           0.248           T
#2      N3      N6      3.201           0.225           T
#STATISTICAL DATA obtained from 2J02.pdb, 9 basepairs
bonds = [2, "G", "A", "N2", "N7", "N3", "N6", "NA", "NA", 15.45, 1.321, 9.523, 0.367]
MASTER_Basepairs_bonds.append(bonds)

#####Homo pyrimidine
###i = 15 in MASTER_Basepairs (XII_UU)
XII_UU = [] #XII Base-pairing pattern UU: UU_20
#Bond   U       U       Length Ave      Length Std      Attribute
#1      N3      O4      2.98    0.06    T
#2      O4      N3      2.74    0.05    T
bonds = [2, "U", "U", "N3", "O4", "O4", "N3", "NA", "NA", "NA", "NA", "NA", "NA"]
MASTER_Basepairs_bonds.append(bonds)
###i = 16 in MASTER_Basepairs (XIII_UU)
XIII_UU = [] #XIII Base-pairing pattern UU: UU_11
#Bond   U       U       Length Ave      Length Std      Attribute
#1      N3      O2      2.70    0.08    T
#2      O2      N3      2.77    0.15    T
bonds = [2, "U", "U", "N3", "O2", "O2", "N3", "NA", "NA", "NA", "NA", "NA", "NA"]
MASTER_Basepairs_bonds.append(bonds)
###i = 17 in MASTER_Basepairs (XIV_CC)
XIV_CC = [] #XIV Base-pairing pattern CC: CC_29
#Bond   C       C       Length Ave      Length Std      Attribute
#1      N3      N4      3.34    0.00    T
#2      N4      N3      2.09    0.00    T
bonds = [2, "C", "C", "N3", "N4", "N4", "N3", "NA", "NA", "NA", "NA", "NA", "NA"]
MASTER_Basepairs_bonds.append(bonds)
###i = 18 in MASTER_Basepairs (XV_CC)
XV_CC = [] #Base-pairing pattern CC: CC_6
#Bond   C       C       Length Ave      Length Std      Attribute
#1      N3      N3      2.92    0.04    P
#2      N4      O2      2.83    0.01    T
#3      O2      N4      3.12    0.07    T
bonds = [3, "C", "C", "N3", "N3", "N4", "O2", "O2", "N4", "NA", "NA", "NA", "NA"]
MASTER_Basepairs_bonds.append(bonds)

#####Hetero pyrimidine
###i = 19 in MASTER_Basepairs (XVII_CU)
XVII_CU = [] #Base-pairing pattern CU: CU_35
#Bond   C       U       Length Ave      Length Std      Attribute
#1      N3      N3      2.98    0.05    T
#2      N4      O2      2.91    0.07    T
#3      O2      O4      3.18    0.06    P
bonds = [3, "C", "U", "N3", "N3", "N4", "O2", "O2", "O4", "NA", "NA", "NA", "NA"]
MASTER_Basepairs_bonds.append(bonds)
###i = 20 in MASTER_Basepairs (XVII_UC)
XVII_UC = [] #Base-pairing pattern UC: UC_35
#Bond   U       C       Length Ave      Length Std      Attribute
#1      N3      N3      2.98    0.05    T
#2      O2      N4      2.91    0.07    T
#3      O4      O2      3.18    0.06    P
bonds = [3, "U", "C", "N3", "N3", "O2", "N4", "O4", "O2", "NA", "NA", "NA", "NA"]
MASTER_Basepairs_bonds.append(bonds)
###i = 21 in MASTER_Basepairs (XVIII_CU)
XVIII_CU = [] #Base-pairing pattern CU: CU_36
#Bond   C       U       Length Ave      Length Std      Attribute
#1      N3      N3      3.16    0.06    T
#2      N4      O4      3.19    0.07    T
bonds = [2, "C", "U", "N3", "N3", "N4", "O4", "NA", "NA", "NA", "NA", "NA", "NA"]
MASTER_Basepairs_bonds.append(bonds)
###i = 22 in MASTER_Basepairs (XVIII_UC)
XVIII_UC = [] #Base-pairing pattern UC: UC_36
#Bond   U       C       Length Ave      Length Std      Attribute
#1      N3      N3      3.16    0.06    T
#2      O4      N4      3.19    0.07    T
bonds = [2, "U", "C", "N3", "N3", "O4", "N4", "NA", "NA", "NA", "NA", "NA", "NA"]
MASTER_Basepairs_bonds.append(bonds)

#####Purine pyrimidine
###i = 23 in MASTER_Basepairs (XIX_CG_WC)
XIX_CG_WC = [] #XIX (Watson-Crick CG)
#Bond    C       G       Length Ave      Length Std      Attribute
#1       N3      N1      3.039           0.195           T
#2       O2      N2      2.931           0.300           T
#3       N4      O6      3.075           0.310           T
#STATISTICAL DATA obtained from 2J02.pdb, 143 basepairs
#bonds = [3, "C", "G", "N3", "N1", "O2", "N2", "N4", "O6", 18.58, 0.900, 10.69, 0.241]
#STATISTICAL DATA obtained from 4TNA.pdb, 6 basepairs
bonds = [3, "C", "G", "N3", "N1", "O2", "N2", "N4", "O6", 18.62, 1.145, 10.87, 0.436]
MASTER_Basepairs_bonds.append(bonds)
###i = 24 in MASTER_Basepairs (XIX_GC_WC)
XIX_GC_WC = [] #XIX (Watson-Crick GC)
#Bond    G       C       Length Ave      Length Std      Attribute
#1       N1      N3      3.005           0.206           T
#2       N2      O2      2.922           0.267           T
#3       O6      N4      3.029           0.326           T
#STATISTICAL DATA obtained from 2J02.pdb, 157 basepairs
#bonds = [3, "G", "C", "N1", "N3", "N2", "O2", "O6", "N4", 18.49, 0.742, 10.64, 0.257]
#STATISTICAL DATA obtained from 4TNA.pdb, 7 basepairs
bonds = [3, "G", "C", "N1", "N3", "N2", "O2", "O6", "N4", 17.89, 1.160, 10.69, 0.390]
MASTER_Basepairs_bonds.append(bonds)
###i = 25 in MASTER_Basepairs (XX_AU_WC)
XX_AU_WC = [] #XX (Watson-Crick AU). Base-pairing pattern AU: AU_2
#Bond   A       U       Length Ave      Length Std      Attribute
#1      N1      N3      2.927           0.216           T
#2      N6      O4      3.032           0.241           T
#STATISTICAL DATA obtained from 2J02.pdb, 21 basepairs
#bonds = [2, "A", "U", "N1", "N3", "N6", "O4", "NA", "NA", 18.68, 0.269, 10.67, 0.27]
#STATISTICAL DATA obtained from 4TNA.pdb, 7 basepairs (AU and UA combined)
bonds = [2, "A", "U", "N1", "N3", "N6", "O4", "NA", "NA", 18.46, 0.258, 10.726, 0.270]
MASTER_Basepairs_bonds.append(bonds)
###i = 26 in MASTER_Basepairs (XX_UA_WC)
XX_UA_WC = [] #XX (Watson-Crick UA). Base-pairing pattern UA: UA_2
#Bond   U       A       Length Ave      Length Std      Attribute
#1      N3      N1      2.974           0.188           T
#2      O4      N6      3.140           0.269           T
#STATISTICAL DATA obtained from 2J02.pdb, 30 basepairs
#bonds = [2, "U", "A", "N3", "N1", "O4", "N6", "NA", "NA", 18.60, 0.816, 10.55, 0.298]
#STATISTICAL DATA obtained from 4TNA.pdb, 7 basepairs (AU and UA combined)
bonds = [2, "A", "U", "N3", "N1", "O4", "N6", "NA", "NA", 18.46, 0.258, 10.726, 0.270]
MASTER_Basepairs_bonds.append(bonds)
###i = 27 in MASTER_Basepairs (XXI_AU)
XXI_AU = [] #XXI (AU Reversed Watson-Crick). Base-pairing pattern AU: AU_30
#Bond   A       U       Length Ave      Length Std      Attribute
#1      N1      N3      2.84    0.13    T
#2      N6      O2      2.94    0.17    T
#STATISTICAL DATA obtained from 2J02.pdb, 5 combined basepairs of the types XXI_AU and XXI_UA. Should be rechecked
bonds = [2, "A", "U", "N1", "N3", "N6", "O2", "NA", "NA", 16.73, 2.117, 10.97, 0.109]
MASTER_Basepairs_bonds.append(bonds)
###i = 28 in MASTER_Basepairs (XXI_UA)
XXI_UA = [] #XXI (UA Reversed Watson-Crick). Base-pairing pattern UA: UA_30
#Bond   U       A       Length Ave      Length Std      Attribute
#1      N3      N1      2.84    0.13    T
#2      O2      N6      2.94    0.17    T
#STATISTICAL DATA obtained from 2J02.pdb, 5 combined basepairs of the types XXI_AU and XXI_UA. Should be rechecked
bonds = [2, "U", "A", "N3", "N1", "O2", "N6", "NA", "NA", 16.73, 2.117, 10.97, 0.109]
MASTER_Basepairs_bonds.append(bonds)
###i = 29 in MASTER_Basepairs (XXII_CG)
XXII_CG = [] #XXII (CG Reversed Watson-Crick). Base-pairing pattern CG: CG_31
#Bond   C       G       Length Ave      Length Std      Attribute
#1      O2      N1      2.80    0.17    T
#2      N3      N2      2.86    0.18    T
bonds = [2, "C", "G", "O2", "N1", "N3", "N2", "NA", "NA", "NA", "NA", "NA", "NA"]
MASTER_Basepairs_bonds.append(bonds)
###i = 30 in MASTER_Basepairs (XXII_GC)
XXII_GC = [] #XXII (GC Reversed Watson-Crick). Base-pairing pattern GC: GC_31
#Bond   G       C       Length Ave      Length Std      Attribute
#1      N1      O2      2.80    0.17    T
#2      N2      N3      2.86    0.18    T
bonds = [2, "G", "C", "N1", "O2", "N2", "N3", "NA", "NA", "NA", "NA", "NA", "NA"]
MASTER_Basepairs_bonds.append(bonds)
###i = 31 in MASTER_Basepairs (XXIII_AU)
XXIII_AU = [] #XXIII (AU Hoogsteen). Base-pairing pattern AU: AU_29
#Bond   A       U       Length Ave      Length Std      Attribute
#1      N6      O4      3.05    0.15    T
#2      N7      N3      2.96    0.15    T
#STATISTICAL DATA obtained from 2J02.pdb, 4 basepairs of the type XXIII_UA. Should be rechecked
bonds = [2, "A", "U", "N6", "O4", "N7", "N3", "NA", "NA", 10.89, 0.436, 8.208, 0.260]
MASTER_Basepairs_bonds.append(bonds)
###i = 32 in MASTER_Basepairs (XXIII_UA)
XXIII_UA = [] #XXIII (UA Hoogsteen). Base-pairing pattern UA: UA_29
#Bond   U       A       Length Ave      Length Std      Attribute
#1      O4      N6      3.05    0.15    T
#2      N3      N7      2.96    0.15    T
#STATISTICAL DATA obtained from 2J02.pdb, 4 basepairs.
bonds = [2, "U", "A", "O4", "N6", "N3", "N7", "NA", "NA", 10.89, 0.436, 8.208, 0.260]
MASTER_Basepairs_bonds.append(bonds)
###i = 33 in MASTER_Basepairs (XXIV_AU)
XXIV_AU = [] #XXIV (AU Reversed Hoogsteen). Base-pairing pattern AU: AU_14
#Bond   A       U       Length Ave      Length Std      Attribute
#1      N6      O2      2.838           0.191           T
#2      N7      N3      2.943           0.048           T
#STATISTICAL DATA obtained from 2J02.pdb,7 basepairs
bonds = [2, "A", "U", "N6", "O2", "N7", "N3", "NA", "NA", 12.24, 1.077, 9.766, 0.205]
MASTER_Basepairs_bonds.append(bonds)
###i = 34 in MASTER_Basepairs (XXIV_UA)
XXIV_UA = [] #XXIV (UA Reversed Hoogsteen). Base-pairing patterMASTER_Basepairs_bonds.append(bonds)
#Bond   U       A       Length Ave      Length Std      Attribute
#1      O2      N6      3.002 SD        0.260           T
#2      N3      N7      2.935           0.249           T
#STATISTICAL DATA obtained from 2J02.pdb,11 basepairs
bonds = [2, "U", "A", "O2", "N6", "N3", "N7", "NA", "NA", 12.20, 1.556, 9.689, 0.334]
MASTER_Basepairs_bonds.append(bonds)
###i = 35 in MASTER_Basepairs (XXV_AC)
XXV_AC = [] #XXV (AC Reversed Hoogsteen). Base-pairing pattern AC: AC_7
#Bond   A       C       Length Ave      Length Std      Attribute
#1      N6      N3      3.355           0.198           T
#2      N7      N4      3.263           0.213           T
#STATISTICAL DATA obtained from 2J02.pdb, 3 basepairs(but compare to XXV_CA)
bonds = [2, "A", "C", "N6", "N3", "N7", "N4", "NA", "NA", 13.02, 0.814, 11.26, 0.228]
MASTER_Basepairs_bonds.append(bonds)
###i = 36 in MASTER_Basepairs (XXV_CA)
XXV_CA = [] #XXV (CA Reversed Hoogsteen). Base-pairing pattern CA: CA_7
#Bond   C       A       Length Ave      Length Std      Attribute
#1      N3      N6      3.067           0.202           T
#2      N4      N7      2.980           0.213           T
#STATISTICAL DATA obtained from 2J02.pdb, 5 basepairs
bonds = [2, "C", "A", "N3", "N6", "N4", "N7", "NA", "NA", 13.10, 1.713, 11.24, 0.149]
MASTER_Basepairs_bonds.append(bonds)
###i = 37 in MASTER_Basepairs (XXVI_AC)
XXVI_AC = [] #XXVI (AC Reversed Wobble). Base-pairing pattern AC: AC_42
#Bond   A       C       Length Ave      Length Std      Attribute
#1      N1      N4      3.01    0.08    T
#2      N6      N3      3.09    0.14    T
bonds = [2, "A", "C", "N1", "N4", "N6", "N3", "NA", "NA", "NA", "NA", "NA", "NA"]
MASTER_Basepairs_bonds.append(bonds)
###i = 38 in MASTER_Basepairs (XXVI_CA)
XXVI_CA = [] #XXVI (CA Reversed Wobble). Base-pairing pattern CA: CA_42
#Bond   C       A       Length Ave      Length Std      Attribute
#1      N4      N1      3.01    0.08    T
#2      N3      N6      3.09    0.14    T
bonds = [2, "C", "A", "N4", "N1", "N3", "N6", "NA", "NA", "NA", "NA", "NA", "NA"]
MASTER_Basepairs_bonds.append(bonds)
###i = 39 in MASTER_Basepairs (XXVII_GU)
XXVII_GU = [] #XXVII (Reversed GU Wobble). Base-pairing pattern GU: GU_29
#Bond   G       U       Length Ave      Length Std      Attribute
#1      N1      O4      2.99    0.07    T
#2      O6      N3      2.99    0.15    T
bonds = [2, "G", "U", "N1", "O4", "O6", "N3", "NA", "NA", "NA", "NA", "NA", "NA"]
MASTER_Basepairs_bonds.append(bonds)
###i = 40 in MASTER_Basepairs (XXVII_UG)
XXVII_UG = [] #XXVII (Reversed UG Wobble). Base-pairing pattern UG: UG_29
#Bond   U       G       Length Ave      Length Std      Attribute
#1      O4      N1      2.99    0.07    T
#2      N3      O6      2.99    0.15    T
bonds = [2, "U", "G", "O4", "N1", "N3", "O6", "NA", "NA", "NA", "NA", "NA", "NA"]
MASTER_Basepairs_bonds.append(bonds)
###i = 41 in MASTER_Basepairs (XXVIII_GU)
XXVIII_GU = [] #XXVIII (GU Wobble). Base-pairing pattern GU: GU_1
#Bond   G       U       Length Ave      Length Std      Attribute
#1      N1      O2      2.885           0.238           T
#2      O6      N3      2.846           0.213           T
#STATISTICAL DATA obtained from 2J02.pdb, 21 basepairs
bonds = [2, "G", "U", "N1", "O2", "O6", "N3", "NA", "NA", 18.45, 0.896, 10.44, 0.245]
MASTER_Basepairs_bonds.append(bonds)
###i = 42 in MASTER_Basepairs (XXVIII_UG)
XXVIII_UG = [] #XXVIII (UG Wobble). Base-pairing pattern UG: UG_1
#Bond   U       G       Length Ave      Length Std      Attribute
#1      O2      N1      2.932           0.214           T
#2      N3      O6      2.846           0.201           T
#STATISTICAL DATA obtained from 2J02.pdb, 22 basepairs
bonds = [2, "U", "G", "O2", "N1", "N3", "O6", "NA", "NA", 18.38, 0.807, 10.49, 0.213]
MASTER_Basepairs_bonds.append(bonds)

####ADDITIONAL BASEPAIRS, NOT IN SAENGER'S COMPILATION
###i = 43 in MASTER_Basepairs (XXIX_AC)
XXIX_AC = [] #XXIX (AC Wobble). Base-pairing pattern AC
#Bond   A       C       Length Ave      Length Std      Attribute
#1      N6      N3      unk             unk
#2      N1      O2      unk             unk
bonds = [2, "A", "C", "N6", "N3", "N1", "O2", "NA", "NA", "NA", "NA", "NA", "NA"]
MASTER_Basepairs_bonds.append(bonds)
###i = 44 in MASTER_Basepairs (XXX_CA)
XXIX_CA = [] #XXIX (CA Wobble). Base-pairing pattern CA
#Bond   C       A       Length Ave      Length Std      Attribute
#1      N3      N6      unk             unk
#2      O2      N1      unk             unk
bonds = [2, "A", "C", "N6", "N3", "N1", "O2", "NA", "NA", "NA", "NA", "NA", "NA"]
MASTER_Basepairs_bonds.append(bonds)
###i = 45 in MASTER_Basepairs (XXXI_GC)
XXX_GC = [] #XXX Base-pair between positions G1034 and C1028 of T. thermophilus 16S rRNA (E. coli numbering)
#Bond   G       C       Length Ave      Length Std      Attribute
#1      N2      N4      3.1             unk
#2      N3      N3      3.2             unk
#Not totally sure it is a real basepair
bonds = [2, "G", "C", "N2", "N4", "N3", "N3", "NA", "NA", "NA", "NA", "NA", "NA"]
MASTER_Basepairs_bonds.append(bonds)
###i = 46 in MASTER_Basepairs (XXXII_CG)
XXX_CG = [] #XXX Base-pair between positions G1034 and C1028 of T. thermophilus 16S rRNA (E. coli numbering)
#Bond   C       G       Length Ave      Length Std      Attribute
#1      N4      N2      3.1             unk
#2      N3      N3      3.2             unk
#Not totally sure it is a real basepair
bonds = [2, "C", "G", "N4", "N2", "N3", "N3", "NA", "NA", "NA", "NA", "NA", "NA"]
MASTER_Basepairs_bonds.append(bonds)
###i = 47 in MASTER_Basepairs (XXXIII_GG)
XXXI_GG = [] #XXXI Base-pair between positions G1030A and C1031 of T. thermophilus 16S rRNA (E. coli numbering)
#Bond   G       G       Length Ave      Length Std      Attribute
#1      N2      N1      2.9             unk
#2      N3      N2      3.4             unk
bonds = [2, "G", "G", "N2", "N1", "N3", "N2", "NA", "NA", "NA", "NA", "NA", "NA"]
MASTER_Basepairs_bonds.append(bonds)


########################### Basepair Lists SECTION #################################

#MASTER_Basepairs will be loaded after the sequential reading of First_List. This list is the output of FUNCTION "run(args)" with the lines carrying possible basepair bond information given the CUTOFF (search CUTOFF)
        #every MASTER_Basepair[i] below will be assigned as many [j] sublists as possible bonds are identified. MASTER_Basepair[i][j] will be composed of
                #k=0, atom_i.resid()
                #k=1, resname_i
                #k=2, atmname_i
                #k=3, atom_j.resid()
                #k=4, resname_j
                #k=5, atmname_j
                #k=6, distance


MASTER_Basepairs_excluded = [[], [], [], [], [], [], [], [], [], [],  [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []]

#MASTER_Basepairs_schemes only carries strings with the names of the basepairing Schemes
MASTER_Basepairs_schemes = ["I_AA", "II_AA", "III_GG", "IV_GG", "V_AA", "VI_GG", "VII_GG", "VIII_AG", "VIII_GA", "IX_AG",  "IX_GA", "X_AG", "X_GA", "XI_AG", "XI_GA", "XII_UU", "XIII_UU", "XIV_CC", "XV_CC", "XVII_CU", "XVII_UC", "XVIII_CU", "XVIII_UC", "XIX_CG_WC", "XIX_GC_WC", "XX_AU_WC", "XX_UA_WC", "XXI_AU", "XXI_UA", "XXII_CG", "XXII_GC", "XXIII_AU", "XXIII_UA", "XXIV_AU", "XXIV_UA", "XXV_AC",  "XXV_CA", "XXVI_AC", "XXVI_CA", "XXVII_GU", "XXVII_UG", "XXVIII_GU", "XXVIII_UG", "XXIX_AC", "XXIX_CA", "XXX_GC", "XXX_CG", "XXXI_GG"]

#MASTER_Basepairs_summary: The lines in MASTER_Basepairs, carrying possible basepairing bonds, will be sorted into this list, once they have been confirmed as possible basepairs. Same format as MASTER_Basepairs
MASTER_Basepairs_summary = [[], [], [], [], [], [], [], [], [], [],  [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []]


##########################################################################
##########################################################################
##########################################################################
#SECTION: LIST DEFINITION                                       END
##########################################################################
##########################################################################
##########################################################################


##########################################################################
##########################################################################
##########################################################################
#SECTION FUNCTIONS                                         BEGINNING
##########################################################################
##########################################################################
##########################################################################

################################################
################################################
#FUNCTION "run(args)"                BEGINNING
################################################
################################################
#FUNCTION run(args) was created by Ralf W. Grosse-Kunstle and provides a list of ATOM-to-ATOM distances calculated from a .pdb file.


def pair_sym_table_as_antons_master(
      master,
      unit_cell,
      pdb_atoms,
      sites_frac,
      pair_sym_table,
      reindexing_array,
      omit_symmetry_interactions=True):
  for table_i_seq,pair_sym_dict in enumerate(pair_sym_table):
    i_seq = reindexing_array[table_i_seq]
    site_i = sites_frac[i_seq]
    atom_i = pdb_atoms[i_seq]
    resname_i = atom_i.resname
    atmname_i = atom_i.name
    for table_j_seq,sym_ops in pair_sym_dict.items():
      j_seq = reindexing_array[table_j_seq]
      site_j = sites_frac[j_seq]
      atom_j = pdb_atoms[j_seq]
      resname_j = atom_j.resname
      atmname_j = atom_j.name
      for sym_op in sym_ops:
        if (omit_symmetry_interactions and not sym_op.is_unit_mx()):
          continue
        site_ji = sym_op * site_j
        distance = unit_cell.distance(site_i, site_ji)
        if atom_i.resid() != atom_j.resid():
           master.append([
             atom_i.resid(),
             resname_i,
             atmname_i,
             atom_j.resid(),
             resname_j,
             atmname_j,
             distance])



def run(args):
  assert len(args) == 1
  import iotbx.pdb
  pdb_inp = iotbx.pdb.input(file_name=args[0])
  crystal_symmetry = pdb_inp.crystal_symmetry()
  sites_cart = pdb_inp.atoms().extract_xyz()
  if (crystal_symmetry is not None):
    unit_cell = crystal_symmetry.unit_cell()
  else:
    from cctbx import crystal, uctbx
    unit_cell = uctbx.non_crystallographic_unit_cell(sites_cart=sites_cart)
    crystal_symmetry = crystal.symmetry(
      unit_cell=unit_cell, space_group_symbol="P1")
  sites_frac = unit_cell.fractionalization_matrix() * sites_cart
  #
  from cctbx.array_family import flex
  pair_sym_table = crystal_symmetry.special_position_settings() \
    .pair_asu_table(
      distance_cutoff=5.0,
      sites_frac=sites_frac).extract_pair_sym_table()
  pdb_atoms = pdb_inp.atoms_with_labels()
  master = [args[0]]
  pair_sym_table_as_antons_master(
    master=master,
    unit_cell=unit_cell,
    pdb_atoms=pdb_atoms,
    sites_frac=sites_frac,
    pair_sym_table=pair_sym_table,
    reindexing_array=flex.size_t_range(sites_frac.size()))
  #
  p_selection = pdb_inp.atoms().extract_name() == " P  "
  p_pair_sym_table = crystal_symmetry.special_position_settings() \
    .pair_asu_table(
      distance_cutoff=25.0,
      sites_frac=sites_frac.select(p_selection)).extract_pair_sym_table()
  pair_sym_table_as_antons_master(
    master=master,
    unit_cell=unit_cell,
    pdb_atoms=pdb_atoms,
    sites_frac=sites_frac,
    pair_sym_table=p_pair_sym_table,
    reindexing_array=p_selection.iselection())
  #
  c1_selection_1 = pdb_inp.atoms().extract_name() ==  " C1*"
  c1_pair_sym_table = crystal_symmetry.special_position_settings() \
    .pair_asu_table(
      distance_cutoff=22.0,
      sites_frac=sites_frac.select(c1_selection_1)).extract_pair_sym_table()
  pair_sym_table_as_antons_master(
    master=master,
    unit_cell=unit_cell,
    pdb_atoms=pdb_atoms,
    sites_frac=sites_frac,
    pair_sym_table=c1_pair_sym_table,
    reindexing_array=c1_selection_1.iselection())

  c1_selection_2 = pdb_inp.atoms().extract_name() ==  " C1'"
  c1_pair_sym_table = crystal_symmetry.special_position_settings() \
    .pair_asu_table(
      distance_cutoff=22.0,
      sites_frac=sites_frac.select(c1_selection_2)).extract_pair_sym_table()
  pair_sym_table_as_antons_master(
    master=master,
    unit_cell=unit_cell,
    pdb_atoms=pdb_atoms,
    sites_frac=sites_frac,
    pair_sym_table=c1_pair_sym_table,
    reindexing_array=c1_selection_2.iselection())

  return master
################################################
################################################
#FUNCTION "run(args)"                   END
################################################
################################################



################################################
################################################
#FUNCTION "CONTROL"               BEGINNING
################################################
################################################
#The bases identified with this function cannot be assigned to a single basepair. The lists of bases involved in MULTIPLE CONTACTS (sublists with more than 2 bases in list 'control', as determined in this FUNCTION, can be later deconvouluted by using the CONTINUITY OF HELICITY CRITERION or the ELIMINATION principle. Should the program not be able to deconvolute sublists with more than 2 bases in list 'control', the user will be prompted to check the assignement of these bases in the structure after displaying it in pymol
#Called from "FUNCTION Section: LOOKING FOR BASES INVOLVED IN MULTIPLE CONTACTS"
#FUNCTION "CONTROL" is executed before any basepair is assigned

def CONTROL(control, resid1, resid2, base1, base2):
    count = 0 #see avoiding repetitions
    collect = []
    position_control = ['none', 'none']
    if len(control) == 0: #First couple of bases in the 'control' list
        collect =[resid1, base1, resid2, base2]
    else:   #The 'control' list has already been appended. The next lines will determine which of the 2 bases from the current basepair candidate have been already identified in other possible basepairs
        for j in range (len(control)):
           if (resid1 in control[j]):
              count = count + 1
              position_control[0] = j  #see avoiding repetitions
           if (resid2 in control[j]):
              count = count + 1
              position_control[1] = j  #see avoiding repetitions
    if (count == 0) and (len(control) > 0):
        collect = [resid1, base1, resid2, base2]
#avoiding repetitions
#Both resid1 and resid2 from MASTER_Basepairs_summary[i][j] could be involved in a base pair interaction, other than the one described by MASTER_Basepairs_summary[i][j]. The following lines will avoid repetitions in the 'control' lines
    if count > 0:
        collect = [resid1, base1, resid2, base2]
        if ('none' not in position_control) and (position_control[1] < position_control[0]):
           position_control = [position_control[1], position_control[0]]
        for i in range (len(position_control)):
           if (position_control[i] != 'none'):
               for j in range (len(control[position_control[i]])/2):
                   if (resid1 not in control[position_control[i]][j*2]) and (resid2 not in control[position_control[i]][j*2]):
                      collect.append(control[position_control[i]][j*2])
                      collect.append(control[position_control[i]][(j*2)+1])
               control[position_control[i]] = []
    if len(collect) > 0:
        control.append(collect)
    return control
################################################
################################################
#FUNCTION "CONTROL"                   END
################################################
################################################


################################################
################################################
#FUNCTION "CONVERT"
################################################
################################################
#Converts a list composed of space-containing strings of the form ' resid # ', ' resid ', ' resid # ', ' resid ', etc into a string of the type 'resid' and 'resid #'
def CONVERT(list):
                   count = 2
                   resid = []
                   for i in range (len(list)):
                       resid.append([])
                   resid_count = 0
#                   print "resid", resid, "list", list
                   while count <= len(resid):
                       if list[resid_count] != []:
                           m=re.search('\d+', list[resid_count])
                           resid[resid_count] = m.group()
                           n=re.search('\w', list[resid_count + 1])
                           resid[resid_count + 1] = n.group()
                           count = count + 2
                           resid_count = resid_count + 2
                       else:
                           count = count + 2
                   return resid
################################################
################################################
#FUNCTION "CONVERT"
################################################
################################################


################################################
################################################
#FUNCTION "PYMOL_OUTPUT"                BEGINNING
################################################
################################################
def PYMOL_OUTPUT(list, h, pdb_file_main, From, file):
#list = new_list_end[b] or
#h = position within the original list from which LIST1 derives
          line = []
          if (h%5)/4 == 1:
             color = "blue"
          elif (h%5)/3 == 1:
             color = "cyan"
          elif (h%5)/2 == 1:
             color = "forest"
          elif (h%5)/1 == 1:
             color = "purple"
          else:
             color = "yellow"
          for k in range (len(file)):
              if file[k] != []:
                 for a in range (1, 3):
                     if From == 'basepairs':
                         m=re.search('\d+', list[2 * a])
                         resid = m.group()
                         RESID = str(resid)
                         line = [["color " + color + ", /" + pdb_file_main + "/*/*/" + RESID + '\n',"hide sticks, /" + pdb_file_main + "/*/*/" + RESID + '\n'],["color " + color + ", /" + pdb_file_main + "/*/*/" + RESID + '\n', "show sticks, /" + pdb_file_main + "/*/*/" + RESID + '\n']]
                         for l in range (len(line[k])):
                             file[k].write(line[k][l])
                     else:
                         m=re.search('\d+', list[2 * a])
                         resid = m.group()
                         RESID = str(resid)
                         line = "color " + color + ", /" + pdb_file_main + "/*/*/" + RESID + '\n'
                         file[k].write(line)
                         line1 = "show sticks, /" + pdb_file_main + "/*/*/" + RESID + '\n'
                         file[k].write(line1)

###############################################
################################################
#FUNCTION "PYMOL_OUTPUT"                   END
################################################
################################################

###############################################
################################################
#FUNCTION "CLOSE"                   END
################################################
################################################
def CLOSE(file):
   for h in range (len(file)):
      if file[h] != []:
          file[h].close()

###############################################
################################################
#FUNCTION "CLOSE"                   END
################################################
################################################


################################################
################################################
#FUNCTION "SEARCH"
################################################
################################################
def SEARCH(string, search):
   m=re.search(search, string)
   SEARCH = m.group()
   return SEARCH
################################################
################################################
#FUNCTION "SEARCH"
################################################
################################################

################################################
################################################
#FUNCTION "to_int"
################################################
################################################
#Searched in Internet. Converts strings with #s to integers
def to_int(in_str):
    """Converts a string to an integer"""
    out_num = 0

    if in_str[0] == "-":
        multiplier = -1
        in_str = in_str[1:]
    else:
        multiplier = 1
    for x in range(0, len(in_str)):
        out_num = out_num * 10 + ord(in_str[x]) - ord('0')
    return out_num * multiplier
################################################
################################################
#FUNCTION "to_int"
################################################
################################################

################################################
################################################
#FUNCTION "to_string"
################################################
################################################
#Searched in Internet. Converts #s to strings
def to_string(in_int):
    """Converts an integer to a string"""
    out_str = ""
    prefix = ""
    if in_int < 0:
        prefix = "-"
        in_int = -in_int
    while in_int / 10 != 0:
        out_str = chr(ord('0') + in_int % 10) + out_str
        in_int = in_int / 10
    out_str = chr(ord('0') + in_int % 10) + out_str
    return prefix + out_str
################################################
################################################
#FUNCTION "to_string"
################################################
################################################


################################################
################################################
#FUNCTION "ADD"
################################################
################################################
def ADD(SEARCH1, SEARCH2):
   A = SEARCH1
   A = to_string(A)
   B = SEARCH2
   B = to_string(B)
   a = SEARCH1 + 1
   a = to_string(a)
   b = SEARCH2 - 1
   b = to_string(b)
   c = SEARCH1 - 1
   c = to_string(c)
   d = SEARCH2 + 1
   d = to_string(d)
   return[a, b, c,d]
################################################
################################################
#FUNCTION "ADD"
################################################
################################################


################################################
################################################
#FUNCTION "SEQUENTIAL_READOUT"   BEGINNING
################################################
################################################
######Sequential readout of lines from lists 'control' and 'MASTER_Basepairs_summary'. Not necessary for execution
#LIST1 = MASTER_Basepairs_summary
#LIST2 = MASTER_Basepairs_schemes
#LIST3 = control
def SEQUENTIAL_READOUT(LIST1, LIST2, LIST3, print_control, line2):
    control_control = 1 #Initializing. This variable signals the presence of bases with possible multiple contacts
    print "\n##################################################\n##### SEQUENTIAL OUTPUT READOUT        BEGINNING"
    if print_control == 0:
        print "##### LINES CARRYING POSSIBLE BASEPAIR INFORMATION",
        print line2
        for i in range (len(LIST1)):
            line1 = "Basepairing scheme: " + LIST2[i]
            print line1
            for k in range (len(LIST1[i])):
                line1 = ''
                for l in range (len(LIST1[i][k])):
                    a = str(LIST1[i][k][l])
                    line1 = line1 + a + ' '
                print line1
        print "#####LINES CARRYING POSSIBLE BASEPAIR INFORMATION",
        print line2
        print "##################################################"
#OUTPUT regarding MULTIPLE CONTACTS
    print "#################################################"
    print "####LIST OF BASES WITH MULTIPLE POSSIBLE CONTACTS",
    print line2
    count1 = 0
    if LIST3 != []:
        line3 = ''
        list1 = []
        list = []
        count = 0
        a = ""
        for i in range (len(LIST3)):
            for j in range (len(LIST3[i])/2):
               if LIST3[i][j*2] != []:
                  count = count + 1
                  list.append(LIST3[i][j*2])
                  list.append(LIST3[i][(j*2)+1])
            if len(list) > 4: #Only sublists with more than two ATOMs should be printed
               count1 = count1 + 1
               list = CONVERT(list)
               for k in range (len(list)):
                  a = a + list[k]
                  if k%2 != 0:
                      a = a + ' '
               line3 = a
               list1.append(line3)
               a = ""
            list = []
        line1 = ""
        line2 = ""
        line3 = ""
    if count1 == 0:
        line1 = "NO BASES DETECTED WITH MULTIPLE POSSIBLE CONTACTS\n#################################################"
        print line1
        control_control = 0 #Bases with multiple contacts = 0. This will avoid unnecessary runs of ELIMINATION
    else:
        print "##### GROUPS OF BASES DETECTED WITH MULTIPLE POSSIBLE CONTACTS = ", count1
        for i in range (len(list1)):
            print "GROUP", i + 1, "  ",
            print list1[i]
        print "####LIST OF BASES WITH MULTIPLE POSSIBLE CONTACTS\n#################################################"
    print "##################################################\n###### SEQUENTIAL OUTPUT READOUT          END\n##################################################\n"

    return control_control
#    import sys
#    sys.exit()
######Sequential readout of lines from lists 'control' and 'MASTER_Basepairs_summary'. Not necessary for execution
################################################
################################################
#FUNCTION "SEQUENTIAL_READOUT"   END
################################################
################################################



################################################
################################################
#FUNCTION "FIRST_APPEND"   BEGINNING
################################################
################################################
#This FUNCTION will be used to compare MASTER_Basepairs_summary to 'control' (already loaded), to determine what to do with base pairs involved in MULTIPLE CONTACTS accroding to the criteria chosen.
        #If the base pair is found to be involved in MULTIPLE CONTACTS, a '1' will be appended to its line in MASTER_Basepairs_summary
                #Having a '1' will not affect the assignment of 3-base basepairs for which all bonds have been identified
        #If the base pair is found not to be involved in MULTIPLE CONTACTS, a '0' will be appended to its line in MASTER_Basepairs_summary

def FIRST_APPEND(LIST1, LIST2, max_cutoff, cutoff):
#LIST1 = MASTER_Basepairs_summary[l][m]
#LIST2 = control

               for n in range (len(LIST2)):
                   append = '0'
                   if (len(LIST2[n]) > 4) and ((LIST1[2] in LIST2[n]) or (LIST1[4] in LIST2[n])):
                      append = '1'
                   if '1' in append:
#Base pair candidate is involved in MULTIPLE CONTACTS. Should be appended as 'no' unless it is a 3-base basepair with all bonds identified
                      if (LIST1[1] == 3) and (LIST1[6] == 3):
                          LIST1.append('yes')
                          LIST1.append('1')
                      if (LIST1[1] == 2) and (LIST1[6] == 2):
                          LIST1.append('no')
                          LIST1.append('1')
                      if (LIST1[1] == 3) and (LIST1[6] == 2):
                          LIST1.append('no')
                          LIST1.append('1')
                   elif (LIST1[2] in LIST2[n]) or (LIST1[4] in LIST2[n]):
#Base pair candidate is not involved in MULTIPLE CONTACTS. Should be appended as 'yes' unless it is a 3-base basepair with only 2 bonds identified. These basepair candidates will only be appended as 'yes' when the MAX_CUTOFF is reached
                      if (LIST1[1] == LIST1[6]):
                          LIST1.append('yes')
                          LIST1.append('0')
                      if (LIST1[1] == 3 and LIST1[6] == 2) and (cutoff != max_cutoff):
                          LIST1.append('no')
                          LIST1.append('0')
                      elif (LIST1[1] == 3 and LIST1[6] == 2) and (cutoff == max_cutoff):
                          LIST1.append('yes')
                          LIST1.append('0')

#               print "LIST1\n", LIST1
               return LIST1
################################################
################################################
#FUNCTION "FIRST_APPEND"   END
################################################
################################################


################################################
################################################
#FUNCTION "FUNCTION: BOND_PLACING"   BEGINING
################################################
################################################
#Called from FUNCTION "program"
def BOND_PLACING(LIST1, LIST2, LIST3):
#LIST1 = MASTER_Basepairs_bonds[i]
#LIST2 = MASTER_Basepairs_summary[i][k]
#LIST3 = atom
    #place the atoms involved in the new bond just detected (currently in MASTER_Basepairs[i][j][2] and MASTER_Basepairs[i][j][5]. Distance is in MASTER_Basepairs[i][j][6]) into the proper location within MASTER_Basepairs_summary[i][k] (See MASTER_Basepairs_summary SCHEME)  by first checking the location of these atoms in MASTER_Basepairs_bonds. There are several possibilities
#FOLLOW A LINE
#                                   if ('  31 ' in LIST2) and ('  39 ' in LIST2):
#                                      print "entering BOND_PLACING", "LIST1", LIST1, "\nLIST2", LIST2, "\nLIST3", LIST3
#                                      print "entering BOND_PLACING"
#FOLLOW A LINE

                                   if LIST1[3] == LIST3[0] and LIST1[4] == LIST3[1]:
                                       LIST2[7] = LIST3[0]
                                       LIST2[8] = LIST3[1]
                                       LIST2[9] = LIST3[2]
                                       LIST2[6] = LIST2[6] + 1
                                   elif LIST1[5] == LIST3[0] and LIST1[6] == LIST3[1]:
                                       LIST2[10] = LIST3[0]
                                       LIST2[11] = LIST3[1]
                                       LIST2[12] = LIST3[2]
                                       LIST2[6] = LIST2[6] + 1
                                   elif LIST1[7] == LIST3[0] and LIST1[8] == LIST3[1]:
                                       LIST2[13] = LIST3[0]
                                       LIST2[14] = LIST3[1]
                                       LIST2[15] = LIST3[2]
                                       LIST2[6] = LIST2[6] + 1
#FOLLOW A LINE
#                                   if ('  31 ' in LIST2) and ('  39 ' in LIST2):
#                                      print "exiting BOND_PLACING", "LIST1", LIST1, "\nLIST2", LIST2, "\nLIST3", LIST3
#                                      print "exiting BOND_PLACING" 
#FOLLOW A LINE
                                   return LIST2
################################################
################################################
#FUNCTION "FUNCTION: BOND_PLACING"   BEGINING
################################################
################################################

################################################
################################################
#FUNCTION "C1-C1 DISTANCE CRITERION"    END
################################################
################################################
#IMPORTANT NOTE 050310: the C1* to C1* may be a much better criterion, as it seems to be more stable than P-P distance
def C1_C1_DISTANCE(LIST1, LIST2, A0, A1, A2, A3):
#LIST1 = MASTER_Basepairs_summary[i][j]
#LIST2 = MASTER_Basepairs_bonds[i]
#A0 = convert[0]
#A1 = convert[1]
#A2 = convert[2]
#A3 = convert[3]
    print "\n...using the C1\'-C1\' distance criterion"
    if LIST1[16] != []:
       diffP_P = LIST1[16] - LIST2[9]
    else:
       LIST1[16] = LIST2[9] #To avoid crashes when a phosphate is missing, like at the 5' end of the molecule, an average value is given. This will solve problems as long as there is an average value in MASTER_Basepairs_bonds, so it could create problems for rare basepairs, but one would not expect such basepairs at the beginning of helices.  
       diffP_P = "NA"
    diffC1_C1 = LIST1[17] - LIST2[11]
    line = "The basepair formed by residues " + A0 + A1 + ":" + A2 + A3 + " displays a " + LIST1[0] + " geometry with a P-P distance of " + str(LIST1[16]) + " and a C1\'-C1\' distance of " + str(LIST1[17]) + ".\n     Empirically determined average C1\'-C1\' distance for this geometry = " + str(LIST2[11]) + " + SD = " + str(LIST2[12]) + "\n     Empirically determined average P-P distance for this geometry = " + str(LIST2[9]) + " + SD = " + str(LIST2[10])
    print line
    if abs(diffC1_C1) > 3 * LIST2[12]:
        LIST1.append('REMOVED')
        print "The C1\'-C1\' distance for this candidate basepair is more than 3 times the recorded standard deviation for this geometry. The basepair will be prevented from further processing by apending it as 'REMOVED'"
    LIST1[18:18] = [abs(diffC1_C1)]
    return LIST1
################################################
################################################
#FUNCTION "C1-C1 DISTANCE CRITERION"    END
################################################
################################################


################################################
################################################
#FUNCTION "Loop"             BEGINNING
################################################
################################################

def loop(MASTER_Basepairs_summary, control, loop_control, CUTOFF, positions):

#The first part of the loop will attempt to expand the initial assignment of legitimate basepairs by using the CONTINUITY OF HELICITY CRITERION. Basepairs for which only two bonds have been identified and tagged as 'no' in "SECTION: LOOKING FOR BASES INVOLVED IN MULTIPLE CONTACTS" will be included in this search. More info in FUNCTION Section: Final assignment of Basepairs.

#The second part of the loop will check whether any basepairs identified by the first part of the loop are involved in MULTIPLE CONTACTS. If so, they will be removed from this list and the existence of 'lone pairs' that could indicate additional basepairs will be checked. Assignment of 'lone pairs' as legitimate basepairs will be done if a 2-bond candidate basepair carrying the 2 basesof the lone pair exists. The criterion is called 'ELIMINATION'

   new_basepairs = 0
   CUTOFF_str = str(CUTOFF)
   run = 0 #This variable will count the number of runs of the loop
   count_found = 1 #This variable will ensure that the loop runs until no more additional 'no'-tagged basepairs can be identified as legitimate. Will be reinitialized to 0 within the loop. In addition this variable will keep track of how many additional basepairs are found in each run of the loop
   total_count = 0 #will keep track of the total amount of basepairs identified in the loop



#First part of LOOP:
   while count_found > 0:
        run = run + 1
#        if run == 1:
        count_found = 0 #Needs to be reinitialized at every run of the loop. Will keep track of all basepairs found by first part of the loop.
        run_str = to_string(run)
        line1 = "##################################\n\"Attempting to find additional basepairs by the CONTINUITY OF HELICITY CRITERION\"       Attempt # " + run_str + "\n##################################"
        print line1
        transient = []
        #First the positions of adjacent basepairs to those involved in MULTIPLE CONTACTS, are calculated and outputed
        for i in range (len(MASTER_Basepairs_summary)):
            for j in range (len(MASTER_Basepairs_summary[i])):
                MH_yes = 'y'
                if ((MASTER_Basepairs_summary[i][j][1] == 3) and (MASTER_Basepairs_summary[i][j][6] == 2) and (MAX_CUTOFF_str[0:3] not in CUTOFF_str[0:3])) or ('REMOVED' in MASTER_Basepairs_summary[i][j]):
                    MH_yes = 'n' #Three-bond-basepairs for which only two bonds have been identified are prevented from entering the rest of the loop lines until MAX_CUTOFF_str[0:3]  == CUTOFF_str[0:3]
                if len(MASTER_Basepairs_summary[i]) > 0:
                    if ('no' in MASTER_Basepairs_summary[i][j]) and ('yes' not in MASTER_Basepairs_summary[i][j]) and MH_yes == 'y':
                        m=re.search('\d+', MASTER_Basepairs_summary[i][j][2])
                        search1 = m.group()
                        SEARCH1 = to_int(search1) #convert string to integers
                        m=re.search('\d+', MASTER_Basepairs_summary[i][j][4])
                        search2 = m.group()
                        SEARCH2 = to_int(search2) #convert string to integers
                        list = [MASTER_Basepairs_summary[i][j][2], MASTER_Basepairs_summary[i][j][3], MASTER_Basepairs_summary[i][j][4], MASTER_Basepairs_summary[i][j][5]]
                        convert = CONVERT(list)
                        add = ADD(SEARCH1, SEARCH2)
                        line1 = "\nBasepair tagged for its involvement in MULTIPLE CONTACTS has been detected!!!!!!.\n  The basepair is formed by residues " + convert[0] + convert[1] + " and " + convert[2] + convert[3] + "\n    Attempting to determine whether it could be used to expand the length of an identified helical segment.\n      Searching for the following possible basepair candidates: " + add[0] + " and " + add[1] + " or " + add[2] + " and " + add[3]
                        print line1

#Second, the whole MASTER_Basepairs_summary is searched for the existence of either one of the adjacent basepairs calculated above. Done in FUNCTION HELICITY
                        where_from = 'Loop' #Used in SECTION HELICITY to determine what to do whether the program is coming from 'Loop' or 'CONTINUITY_HELICITY_CRITERION', or from '    print "AFTER LIST1", LIST1REMOVAL'
                        helicity = HELICITY(MASTER_Basepairs_summary, add[0], add[1], add[2], add[3], where_from)
                        count_helicity = 0 #Control output in this region
                        for t in range (0,2):
                            if helicity[t] != []:
                                count_helicity = count_helicity + 1  #Signals entering this 'if'
                                str(helicity[t][0])
                                str(helicity[t][1])

#One or the 2 adjacent basepairs is found
                                if count_helicity == 1: #One of the 2 possible adjacent basepairs is found. This is enought to assign the basepair whose assignemnt is in question by 'MH' criterion
                                    count_found = count_found + 1 #Will keep tract of how many basepairs are found. Serves also to keep the loop going
                                    line1  = "     #####\n     !!!!Found an adjacent basepair to residues " + convert[0] + convert[1] + " and "  + convert[2] + convert[3] +  "\n     This basepair is formed by residues " + helicity[t][0] + " and " +  helicity[t][1] + "\n     #####\n    Assigning " + convert[0] + convert[1] + " and "  + convert[2] + convert[3] + " as a legitimate basepair by the CONTINUITY OF HELICITY CRITERION\n"
                                    print line1
                                    MASTER_Basepairs_summary[i][j].append('yes')
                                    MASTER_Basepairs_summary[i][j].append('MH')  #'MH'Maximization of Helicity TAG
                                    transient = [MASTER_Basepairs_summary[i][j][2], MASTER_Basepairs_summary[i][j][4]]
                                    positions.append(transient)
                                    CUTOFF_str = repr(CUTOFF) #Converts 'float' to 'str'
                                    MASTER_Basepairs_summary[i][j].append(CUTOFF_str[0:3])
                                if count_helicity == 2: #The second adjacent basepair has been found. The basepair whose assignemnt is in question has already bee assigned when 'count_helicity == 1'. The program will only output the new finding
                                    line1  = "     #####\n     !!!!Found a second adjacent basepair to residues " + convert[0] + convert[1] + " and "  + convert[2] + convert[3] +  "\n     This basepair is formed by residues " + helicity[t][0] + " and " +  helicity[t][1] + "\n     #####\n"
                                    print line1

#Neither one of the possible adjacent basepairs was found by the CONTINUITY OF HELICITY CRITERION
                        if count_helicity == 0:
                            line1  = "                     Neither basepair was found!!!!!"
                            print line1

#Time to recap
        a = str(count_found)
        b = str(run_str )
        total_count = total_count + count_found
        line1 = "###### Number of basepairs identified by the CONTINUITY OF HELICITY CRITERION = " + a + "      Attempt # " + run_str
        print line1

#Second part of LOOP
        #Attempting to find additional basepairs by the 'ELIMINATION' criterion. running REMOVAL in 'semi-bulk' mode by using a terminal 'MH' as a tag.
        removal = []
        if (total_count > 0):
            mode_removal = 1 #"Run FUNCTION REMOVAL" in 'semi-bulk' mode
            CUTOFF_str = str(CUTOFF)
            print "\n##############################\n REMOVING the newly identified basepairs from the list of bases with possible mutliple contacts.", run_number_str, " Distance Cutoff = ", CUTOFF_str[0:3], "\n##############################\n "
            removal = REMOVAL(MASTER_Basepairs_summary, control, convert[0], convert[2], CUTOFF, mode_removal, positions)
            MASTER_Basepairs_summary = removal[0]
            control = removal[1]
            new_basepairs = removal[2]
            count_found = new_basepairs #If additional basepairs have been assigned by 'ELIMINATION' the Loop should run again
            line1 = "\n\n##### Semi-Bulk Removal"
            print line1
            a = str(new_basepairs)
            line1 = a + " additional basepairs have been identified by ELIMINATION"
            print line1
            line1 = "##### Semi-Bulk Removal"
            print line1

        a = ''

#END OF LOOP
   return [MASTER_Basepairs_summary, control] #The original MASTER_Basepairs_summary and control lists are returned after modification

################################################
################################################
#FUNCTION "Loop"                        END
################################################
################################################


################################################
################################################
#FUNCTION "REMOVAL"                    BEGINNING
################################################
################################################
#Will remove newly identified basepairs from 'control' list, now in LIST2 within FUNCTION REMOVAL, carrying basepairs involved in MULTIPLE CONTACTS.
#After the newly identified basepairs are removed and the removed LIST2 sublists loaded in 'control_removed', the latterlist  is checked for the presence of 'lone pairs' that could represent additional basepairs
#Can be called from:
#1) "FUNCTION Section: Final assignment of Basepairs".
#2) "FUNCTION "Loop""

     #Uses FUNCTION HELICITY and

     #SUBFUNCTION "DEFINING_A_B"
          #Determines what 'a' and 'b' values are used in REMOVAL,
     #SUBFUNCTION "ACTUAL_REMOVAL"
          #Where lists LIST1 and LIST2 will be compared.


#'base_a' and 'base_b' = residues 1 and 2 of the current basepair (only if run in 'semi-bulk' mode)
#Two ways of executing this FUNCTION.
    #In 'semi-bulk mode' ('mode_removal' = 1). When called from "FUNCTION "Loop". It should only remove basepairs tagged as 'MH'
    #In 'bulk' mode ('mode_removal' = 0). For the whole LIST1 = MASTER_Basepairs_summary. When called from "FUNCTION Section: Final assignment of Basepairs"
       #Controlled by the argument 'mode_removal', which may be '1' ('semi-bulk' mode) or '0' ('bulk' mode). Peformed in SUBFUNCTION "DEFINING_A_B" above
       #Bases from legitimate basepairs that had been included in 'control' (LIST2), will be removed from this list. Removed 'control' sublists will be loaded into 'control_removed'.
       #In a second step, the remaining bases, loaded into 'control_removed', will be searched for the presence of 'lone pairs' that may be indicative of additional basepairs in "Searching for 'lone pairs' left over in 'control'"


#Identifying basepairs involved in MULTIPLE CONTACTS.
   #Most of the searching is performed in  ACTUAL_REMOVAL

#LIST1 = MASTER_Basepairs_summary
#LIST2 = control
def REMOVAL(LIST1, LIST2, base_a, base_b, CUTOFF, mode_removal, positions):
            CUTOFF_str = repr(CUTOFF) #Converts 'float' to 'str'
            new_basepairs = 0
            ACTUAL = []
            output_control = 0
            list = [] #Will carry the reminder of LIST2 sublists, once after "removal". Thenthe contents of "list" will be appendedto "control_removed". Has to be defined here and at the beginning of SUBFUNCTION "ACTUAL_REMOVAL"
            control_removed = [] #Will carry the control sublists, now in LIST2 within FUNCTION REMOVAL, after basepairs, identified as legitimate, have been removed

            for t in range (len(LIST2)):
               if len(LIST2[t]) > 4:

                   for u in range (len(LIST2[t])/2):
                      if LIST2[t][u*2] != []:
#                          print "LIST2[t][u*2]", LIST2[t][u*2]
                          for a in range (len(positions)):
#                              print "positions[a]", positions[a]
                              if LIST2[t][u*2] in positions[a]:
                                  line1 = "Residue" + LIST2[t][u*2] + "has been removed from the list of residues involved in MULTIPLE CONTACTS"
#                                  print line1
                                  LIST2[t][u*2] = []
                                  LIST2[t][(u*2)+1] = []


#Once all basepairs have been removed, the modified LIST2 ('control'), will be searched for lines with removed residues. The remaining residues will be transferred to 'control_removed' to search for 'lone pairs'.

            for l in range (len(LIST2)):
                 list = []
                 if LIST2 != []:
                    if [] in LIST2[l]:
#                       print "LIST2[l]", LIST2[l]
                       for m in range (len(LIST2[l])/2):
                          if LIST2[l][m*2] != []:
                             list.append(LIST2[l][m*2])
                             list.append(LIST2[l][(m*2)+1])
                       if list != []:
                          list.append(l) #This will indicate which line of LIST2 ('control') contains possible lone pairs
#                          print "Appending control_removed"
                          control_removed.append(list)
#                          print "control_removed", control_removed

#Searching for 'lone pairs' left over in 'control'
#Now the program searches through 'control_removed' to determine whether any "lone pairs" have been left after the removal process. This event should be indicated by the length of  the list 'control_removed', loaded above
            for i in range (len(control_removed)):
#               print "checking whether removed 'control' lines may carry \"lone pairs\" of bases which could correspond to actual basepairs"
               if len(control_removed[i]) == 5: #resid #1, resid1, resid #2, resid2, 't'or location in LIST2 (control), as calculated above
                   LIST2[control_removed[i][4]] = []
                   convert = CONVERT(control_removed[i])
                   line1 = "\n#####\nFound a lone pair of bases after removing newly identified basepairs. Bases are:\n" + convert[0] + convert[1] + " and " + convert[2] + convert[3] + "\nMatching against basepair list\n##### "
                   print line1
                   where_from = 'REMOVAL' #Used in SECTION HELICITY to determine what to do whether the program is coming from 'Loop' or 'CONTINUITY_HELICITY_CRITERION', or from 'REMOVAL'
                   a = "iiiiiiii" #Just to fill up the list of 'HELICITY''s arguments
                   b = "iiiiiiii" #Just to fill up the list of 'HELICITY''s arguments
                   helicity = HELICITY(LIST1, control_removed[i][0], control_removed[i][2], a, b, where_from)
                   line = ""
                   if helicity[2] != []:
                       if (LIST1[helicity[2][0]][helicity[2][1]][1] == 3 and LIST1[helicity[2][0]][helicity[2][1]][6] == 2) and (MAX_CUTOFF_str[0:3] not in CUTOFF_str[0:3]):
#candidate 3-base basepairs for which only 2 bases have been identified will be excluded from final assignment by 'ELIMINATION' until MAX_CUTOFF is reached.
                          line = "However, this is a 3-base basepairs for which only 2 bases have been so far identified, final assignment will await until final distance CUTTOFF is used"
                       else:
                          LIST1[helicity[2][0]][helicity[2][1]].append('yes')
                          LIST1[helicity[2][0]][helicity[2][1]].append('ELI')  #Assignment by 'Elimination' TAG.
                          CUTOFF_str = repr(CUTOFF) #Converts 'float' to 'str'
                          LIST1[helicity[2][0]][helicity[2][1]].append(CUTOFF_str[0:3])
                          new_basepairs = new_basepairs + 1
                       line1 = "     #####\n     Found!!! A basepair between residues " + convert[0] + convert[1] + " and " + convert[2] + convert[3] + " has been identified by the 'ELIMINATION' criterion\n     #####\n"
                   else:
                       line1 = "No basepair was identified"
                   print line1
                   print line
            return [LIST1, LIST2, new_basepairs]

################################################
################################################
#FUNCTION "REMOVAL"                    END
################################################
################################################


################################################
################################################
#FUNCTION "program"               BEGINNING
################################################
################################################

#This is the main function of the program. It started as the program itself but had to be converted to a function to allow its repetitive iteration
        #Called from 'MAIN'

#First_List = Contains the filtered output from FUNCTION "run(args)". List of ATOM to ATOM distances below a rough cutoff of 5A
def program(First_List, MASTER_Basepairs_summary, CUTOFF, pdb_file_main, control, run_number_str, CUTOFF_str, run_number, positions):

    First_List_smaller = []
    First_List_P = []
    First_List_C1 = []
    for a in range (len(First_List)):
       if First_List[a] != []:
          if (First_List[a][len(First_List[a]) - 1] < 5):
              First_List_smaller.append(First_List[a])
          if ' P  ' in First_List[a][2] and ' P  ' in First_List[a][5]:
              First_List_P.append(First_List[a])
          C1 = [' C1\'', ' C1*']
          if First_List[a][2] in C1 and First_List[a][5] in C1:
#FOLLOW A LINE
#              if ('  31 ' in First_List[a][0] and '  39 ' in First_List[a][3]):
#                 print "First_List[a] ", First_List[a]
#FOLLOW A LINE

              First_List_C1.append(First_List[a])
    First_List = First_List_smaller
    MASTER_Basepairs = [[], [], [], [], [], [], [], [], [], [],  [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []]

#####################################
#FUNCTION "Program". Section: First sorting of basepair candidates          BEGINNING
#####################################
    line1 = "####################################\nFirst sorting of basepair candidates: Run # " + run_number_str + " Distance Cutoff = " + CUTOFF_str[0:3] + "\n#####################################"
    print line1
    collect = []
    for i in range (1, len(First_List)):
#FOLLOW A LINE
#        if ('  31 ' in First_List[i][0] and '  39 ' in First_List[i][3]):
#            print "First_List[i] ", First_List[i]
#FOLLOW A LINE
        if First_List[i][0] != First_List[i][3]:
            for j in range(len(MASTER_Basepairs)):
               #Matching the bases to those of MASTER_Basepairs_bonds
               if (MASTER_Basepairs_bonds[j][1] in First_List[i][1]) and (MASTER_Basepairs_bonds[j][2] in First_List[i][4]):
                   count = 0
                   transient = MASTER_Basepairs_bonds[j][3:9]
#FOLLOW A LINE
#                   if ('  31 ' in First_List[i][0] and '  39 ' in First_List[i][3]):
#                         print "transient", transient, "\nMASTER_Basepairs_bonds[j]", MASTER_Basepairs_bonds[j]
#FOLLOW A LINE
                   #Matching the ATOMS to those of MASTER_Basepairs_bonds
                   for k in range ((len(transient)/2)):
#                       print "transient", transient
                       if (transient[2*k] in First_List[i][2]) and (transient[2*k + 1] in First_List[i][5]) and (First_List[i][6] > PREV_CUTOFF) and (First_List[i][6] <= CUTOFF):
                           collect = [First_List[i][0], First_List[i][1], First_List[i][2], First_List[i][3], First_List[i][4], First_List[i][5], First_List[i][6]]
                           MASTER_Basepairs[j].append(collect)
#FOLLOW A LINE
#                           if ('  31 ' in First_List[i][0] and '  39 ' in First_List[i][3]):
#                               print "\nMASTER_Basepairs_schemes[j]", MASTER_Basepairs_schemes[j]
#                               print "MASTER_Basepairs[j][len(MASTER_Basepairs[j])-1]", MASTER_Basepairs[j][len(MASTER_Basepairs[j])-1]
#FOLLOW A LINE

                       elif (transient[2*k] in First_List[i][2]) and (transient[2*k + 1] in First_List[i][5]) and (First_List[i][6] > CUTOFF) and (MAX_CUTOFF_str[0:3] in CUTOFF_str[0:3]):
                           collect = [First_List[i][0], First_List[i][1], First_List[i][2], First_List[i][3], First_List[i][4], First_List[i][5], First_List[i][6], MASTER_Basepairs_schemes[j]]
                           MASTER_Basepairs_excluded[j].append(collect)
#FOLLOW A LINE
#                           if ('  31 ' in First_List[i][0] and '  39 ' in First_List[i][3]):
#                               for q in range (len(MASTER_Basepairs_excluded[j])):
#                                  print "MASTER_Basepairs_excluded[j][q]", MASTER_Basepairs_excluded[j][q]
#                               print "MASTER_Basepairs_excluded[j][len(MASTER_Basepairs_excluded[j])-1]", MASTER_Basepairs_excluded[j][len(MASTER_Basepairs_excluded[j])-1]
#FOLLOW A LINE
                       count = count + 1
                       collect = []

#At this point the list MASTER_Basepairs carries all the lines that are consistent with basepair interactions based only on the identity of the bases involved and the CUTOFF value. To learn which of these correspond to actual basepairs, the program has to determine how many lines per every two identified base candidates are contained in the corresponding basepairing Scheme sub-list. Done in "Second sorting of basepair candidates"
#####################################
#       FUNCTION "Program". Section: First sorting of basepair candidates                        END
#####################################


#####################################
#       FUNCTION "Program". Section: Second sorting of basepair candidates               BEGINNING
#####################################
#To find out how many lines of MASTER_Basepairs correspond to the same candidate basepair, the program will go over every MASTER_Basepairs sublist and assign it to a candidate basepair. This information will be appended to MASTER_Basepairs_summary, which can be accessed like the other MASTER_ lists.
        #MASTER_Basepairs_summary SCHEME: The sub-lists in MASTER_Basepairs_summary will be organized by basepairing scheme [i] and then by basepair candidate [j]. It will carry the following values:
                #position 0: "Basepairing Scheme" = MASTER_Basepairs_schemes[i]
                #position 1: "number of bonds expected" = MASTER_Basepairs_bonds[i][0],
                #position 2: "Base # for base 1" = MASTER_Basepairs[i][0][0],
                #position 3: "Base identity for base 1" = MASTER_Basepairs[i][0][1],
                #position 4: "Base # for base 2" = MASTER_Basepairs[i][0][3],
                #position 5: "Base identity for base 2" = MASTER_Basepairs[i][0][4],
                #position 6: "number of bonds found" = calculated below,
                #position 7: "atom 1 for bond 1 according to MASTER_Basepairs_bonds[i][3]" MASTER_Basepairs[i][0][2],

                #position 8: "atom 2 for bond 1 according to MASTER_Basepairs_bonds[i][4]" MASTER_Basepairs[i][0][5],
                #position 9: "Distance for H-bond 1 between atoms 1 and 2" MASTER_Basepairs[i+x][0][6],
                #position 10: "atom 1 for bond 2 according to MASTER_Basepairs_bonds[i][5]" MASTER_Basepairs[i+x][0][2],
                #position 11: "atom 2 for bond 2 according to MASTER_Basepairs_bonds[i][6]" MASTER_Basepairs[i+x][0][5],
                #position 12: "Distance for H-bond 2 between atoms 1 and 2" MASTER_Basepairs[i+x][0][6],
                #position 13: "atom 1 for bond 3 according to MASTER_Basepairs_bonds[i][7]" MASTER_Basepairs[i+x+2][0][2],
                #position 14: "atom 2 for bond 3 according to MASTER_Basepairs_bonds[i][8]" MASTER_Basepairs[i+x+s][0][5]]
                #position 15: "Distance for H-bond 3 between atoms 1 and 2" MASTER_Basepairs[i+x][0][6],
                        #NOTE: The latter three strings will only be present in three-bond basepairing schemes
                #position 16: "Phosphate-to-phosphate 'P-to-P' distance". Carried in 'First_List_P' and collected in 'FUNCTION Program. Section: LOOKING FOR BASES WITH POSSIBLE MULTIPLE CONTACTS'
                #position 17: "C1'-to-C1'" distance". Carried in 'First_List_C1' and collected in 'FUNCTION Program. Section: LOOKING FOR BASES WITH POSSIBLE MULTIPLE CONTACTS'
    line1 = "####################################\nSecond sorting of basepair candidates: Run # " + run_number_str + " Distance Cutoff = " + CUTOFF_str[0:3] + "\n#####################################"
    print line1
    for i in range(len(MASTER_Basepairs)): #'i' will correspond to the same level in both MASTER_Basepairs and MASTER_Basepairs_summary, the level of basepair scheme
        if len(MASTER_Basepairs[i]) > 0:
            for j in range (len(MASTER_Basepairs[i])):
#Detection of lines carrying consecutive bases. These should be not allowed to enter the assignment section, as they can give rise to FALSE basepairs
                resid = [[],[]]
                atom = [[],[],[]]
                m=re.search('\d+', MASTER_Basepairs[i][j][0])
                resid[0] = int(m.group())
                m=re.search('\d+', MASTER_Basepairs[i][j][3])
                resid[1] = int(m.group())
                m=re.search('\S+', MASTER_Basepairs[i][j][2])
                atom[0] = str(m.group())
                m=re.search('\S+', MASTER_Basepairs[i][j][5])
                atom[1] = str(m.group())
                atom[2] = MASTER_Basepairs[i][j][6]
                if (resid[1] - resid[0] > 1) or (resid[1] - resid[0] < -1):

#FOLLOW A LINE
#                    if ('  31 ' in MASTER_Basepairs[i][j][0]) and ('  39 ' in MASTER_Basepairs[i][j][3]):
#                        print "\nBefore Second sorting: MASTER_Basepairs[i][j]", MASTER_Basepairs_schemes[i], MASTER_Basepairs[i][j], "atom", atom
#                        print "i", i, "j", j, "MASTER_Basepairs_schemes[i]", MASTER_Basepairs_schemes[i], "MASTER_Basepairs_bonds[i]", MASTER_Basepairs_bonds[i]

#FOLLOW A LINE
                    if (len(MASTER_Basepairs_summary[i]) > 0):
                        found = 0 #Will control that the program does not execute the next 'elif' if it has executed this 'if'
                        for k in range (len(MASTER_Basepairs_summary[i])):

#This basepair has already been detected, new bond for the basepair. Here we should distinguish between RUN 1 and later RUNs, as during the latter the program will run into pre-loaded MASTER_Basepairs_summary[i][k] lines
                           if (MASTER_Basepairs[i][j][0] in MASTER_Basepairs_summary[i][k][2]) and (MASTER_Basepairs[i][j][3] in MASTER_Basepairs_summary[i][k][4]):
                              found = 1 #Will control that the program does not execute the next 'elif' if it has executed this 'if'
                             # if (MASTER_Basepairs[i][j][2] not in MASTER_Basepairs_summary[i][k][7]) and (MASTER_Basepairs[i][j][2] not in MASTER_Basepairs_summary[i][k][10]) and (MASTER_Basepairs[i][j][2] not in MASTER_Basepairs_summary[i][k][13]):
                              if (atom[0] != MASTER_Basepairs_summary[i][k][7]) and (atom[0] != MASTER_Basepairs_summary[i][k][10]) and (atom[0] !=  MASTER_Basepairs_summary[i][k][13]):
########Checking that MASTER_Basepairs[i][j] and MASTER_Basepairs_summary[i][k] contain the same residues but that the former brings new ATOMS
#                                  MASTER_Basepairs_summary[i][k][6] = MASTER_Basepairs_summary[i][k][6] + 1  #One new bond added
                              #Bond ATOMS and distances will be placed at their proper position in "FUNCTION BOND_PLACING"
                                  MASTER_Basepairs_summary[i][k] = BOND_PLACING(MASTER_Basepairs_bonds[i], MASTER_Basepairs_summary[i][k], atom)
#FOLLOW A LINE
#                                  if ('  31 ' in MASTER_Basepairs[i][j][0]) and ('  39 ' in MASTER_Basepairs[i][j][3]):
#                                      print "new bond for the basepair: MASTER_Basepairs_summary[i][k]", MASTER_Basepairs_summary[i][k]
#FOLLOW A LINE


#New basepair detected for a basepairing Scheme for which other basepairs have been identified
                           elif (found == 0) and (k == len(MASTER_Basepairs_summary[i]) - 1):
                               collect = [MASTER_Basepairs_schemes[i], MASTER_Basepairs_bonds[i][0], MASTER_Basepairs[i][j][0], MASTER_Basepairs[i][j][1], MASTER_Basepairs[i][j][3], MASTER_Basepairs[i][j][4], 0, [], [], [], [], [], [], [], [], [], [], []]
                               MASTER_Basepairs_summary[i].append(collect)
                       #Bond ATOMS and distances will be placed at their proper position in "FUNCTION BOND_PLACING"
                               MASTER_Basepairs_summary[i][len(MASTER_Basepairs_summary[i])-1] = BOND_PLACING(MASTER_Basepairs_bonds[i], MASTER_Basepairs_summary[i][len(MASTER_Basepairs_summary[i])-1], atom)
#FOLLOW A LINE
#                               if ('  31 ' in MASTER_Basepairs[i][j][0]) and ('  39 ' in MASTER_Basepairs[i][j][3]):
#                                  print "New basepair detected for a basepairing Scheme: MASTER_Basepairs_summary[i][k]", MASTER_Basepairs_summary[i][k]
#FOLLOW A LINE

#First basepair for this basepairing Scheme
                    else: ######
                        collect = [MASTER_Basepairs_schemes[i], MASTER_Basepairs_bonds[i][0], MASTER_Basepairs[i][j][0], MASTER_Basepairs[i][j][1], MASTER_Basepairs[i][j][3], MASTER_Basepairs[i][j][4], 0, [], [], [], [], [], [], [], [], [], [], []]
                        MASTER_Basepairs_summary[i].append(collect)
                        #Bond ATOMS and distances will be placed at their proper position in "FUNCTION BOND_PLACING"
                        MASTER_Basepairs_summary[i][0] = BOND_PLACING(MASTER_Basepairs_bonds[i], MASTER_Basepairs_summary[i][0], atom)
#FOLLOW A LINE
#                        if ('  31 ' in MASTER_Basepairs[i][j][0]) and ('  39 ' in MASTER_Basepairs[i][j][3]):
#                            print "First basepair for this basepairing Scheme: MASTER_Basepairs_summary[i][0]", MASTER_Basepairs_summary[i][0]
#FOLLOW A LINE
#####################################
#       FUNCTION "Program". Section: Second sorting of basepair candidates                END
#####################################

########################################################
#       FUNCTION "Program". Section: LOOKING FOR BASES WITH POSSIBLE MULTIPLE CONTACTS
########################################################
#The following code is designed to identify bases from basepair candidates carried in MASTER_Basepairs_summary, which are involved in MULTIPLE CONTACTS. These bases will not be assigned to base pairs, except if they belong to 3-bond basepairs for which all 3 bonds have been identified.
    #DEFINITION of MULTIPLE CONTACTS: These are ATOM-TO-ATOM contacts identified by 'iotbx.pdb' which are below the CUTOFF and that do not necessarily reflect basepairing interactions but obscure the assignment of such interactions. Indeed, a large part of all the contacts identified in "SECTION: First sorting of basepair candidates" are not from actual base pairs but from stacked bases. In addition, some of these multiple interactions could be used to detect base triples (NOT CLEAR HOW TO DO THIS YET).

#Due to the importance of these MULTIPLE CONTACTS, this section is devoted to find them. After detecting the existence of MULTIPLE CONTACTS, the program will attempt to remove as  many bases as possible from them by using the "CONTINUITY OF HELICITY CRITERION" in the next section.

    #Searching for these bases has to be done before the assignment of legitimate basepairs, so that the complete list of multiply involved bases can be used during this assignment in SECTION: Final assignment of Basepairs
   #In addition to the code contained in this section, function 'CONTROL' is in charge of performing this search. The list of bases involved in MULTIPLE CONTACTS will be loaded into 'control' by function 'CONTROL'.


    print   "#####################################\nSearching for bases with possible MULTIPLE CONTACTS Run #", run_number_str, " Distance Cutoff = ", CUTOFF_str[0:3], "\n#####################################"
    for i in range (len(MASTER_Basepairs_summary)):
        if len(MASTER_Basepairs_summary[i]) > 0:
            for j in range (len(MASTER_Basepairs_summary[i])):
                append = '0' #Initializing. Will control what to append
                if MASTER_Basepairs_summary[i][j] != [] and MASTER_Basepairs_summary[i][j][16] == []:
                    for r in range (len(First_List_P)):
                        if First_List_P[r][0] in MASTER_Basepairs_summary[i][j] and First_List_P[r][3] in MASTER_Basepairs_summary[i][j][4]:
                           MASTER_Basepairs_summary[i][j][16] = First_List_P[r][6]
                    C1_compare = 22
                    for s in range (len(First_List_C1)):
                        if First_List_C1[s][0] in MASTER_Basepairs_summary[i][j] and First_List_C1[s][3] in MASTER_Basepairs_summary[i][j][4]:
                            if First_List_C1[s][6] < C1_compare:
                                C1_compare = First_List_C1[s][6]
                    MASTER_Basepairs_summary[i][j][17] = C1_compare

#Now that the C1-C1 distances have been added to the basepair candidates, this distance can be used to perform a first discrimination by weeding out basepairs with C1-C1 distances that are more than 3X the precalculated average C1-Ca. This will be done in FUNCTION C1_C1_DISTANCE. Candidate basepairs with 'bad' C1-C1 distances will be appended as 'REMOVED'
#                    for i in range (len(MASTER_Basepairs_summary)):
#                        for j in range (len(MASTER_Basepairs_summary[i])):
                    if (MASTER_Basepairs_bonds[i][9] != 'NA'):
                        list = [MASTER_Basepairs_summary[i][j][2], MASTER_Basepairs_summary[i][j][3], MASTER_Basepairs_summary[i][j][4], MASTER_Basepairs_summary[i][j][5]]
                        convert = CONVERT(list)
                        MASTER_Basepairs_summary[i][j] = C1_C1_DISTANCE(MASTER_Basepairs_summary[i][j], MASTER_Basepairs_bonds[i], convert[0], convert[1], convert[2], convert[3])


                if (MASTER_Basepairs_summary[i][j][6] >= 2) and ('REMOVED' not in MASTER_Basepairs_summary[i][j]):

#Calling FUNCTION "CONTROL" to perform an initial assessment of whether the identified bases might be involved in MULTIPLE CONTACTS
                    control = CONTROL(control, MASTER_Basepairs_summary[i][j][2], MASTER_Basepairs_summary[i][j][4], MASTER_Basepairs_summary[i][j][3], MASTER_Basepairs_summary[i][j][5])

######ASSIGNING '1' or '0' to candidate basepairs(see also explanation for DISTANCE criterion)######
#Once the 'control' list has been loaded, MASTER_Basepairs_summary has to be compared to 'control' in order to properly tag all basepairs candidates prior to their assignment as legitimate basepairs.
        #Performed in #This FUNCTION will be used to compare MASTER_Basepairs_summary to 'control' (already loaded), to determine what to do with bases with possible multiple contacts according to the criteria chosen.
                #If the base pair is found to be possibly involved in MULTIPLE CONTACTS, a '1' will be appended to its line in MASTER_Basepairs_summary
                        #Having a '1' will not affect the assignment of 3-base basepairs for which all bonds have been identified
                #If the base pair is found not to be involved in MULTIPLE CONTACTS, a '0' will be appended to its line in MASTER_Basepairs_summary

    for l in range (len(MASTER_Basepairs_summary)):
        for m in range (len(MASTER_Basepairs_summary[l])):
            if ('D' not in MASTER_Basepairs_summary[l][m]) and ('MH' not in MASTER_Basepairs_summary[l][m]) and ('ELI' not in MASTER_Basepairs_summary[l][m]):
               first_append = FIRST_APPEND(MASTER_Basepairs_summary[l][m], control, MAX_CUTOFF_str[0:3], CUTOFF_str[0:3])
               MASTER_Basepairs_summary[l][m] = first_append
######ASSIGNING '1' or '0' to candidate basepairs(see also explanation for DISTANCE criterion)######

#Printing output.
    line1 = "AFTER SEARCHING FOR BASES WITH POSSIBLE MULTIPLE CONTACTS: Run # " + run_number_str + " Distance Cutoff = " + CUTOFF_str[0:3]
    print_control = 1
    control_control = SEQUENTIAL_READOUT(MASTER_Basepairs_summary, MASTER_Basepairs_schemes, control, print_control, line1) #control_control = 0 signals that no bases with possible multiple contacts have been found. Will avoid unnecesary runs of some parts of the program
########################################################
#       FUNCTION "Program". Section: LOOKING FOR BASES WITH POSSIBLE MULTIPLE CONTACTS                 END
########################################################



################################################
#       FUNCTION "Program". Section: Final assignment of Basepairs BEGINNING
################################################

#Performed by following three different criteria:
        #First: DISTANCE criterion
        #Second: CONTINUITY OF HELICITY CRITERION
        #Third: ELIMINATION criterion

#1) DISTANCE CRITERION (Performed in #######SUBSECTION Initial Assignment):
                #First: The numeral '0' (not involved in MULTIPLE CONTACTS), or '1' (involved in MULTIPLE CONTACTS), will be appended to all candidate basepairs after comparing MASTER_Basepairs_summary to 'control'.
                #Second: All possible legitimate basepairs will be identified based on the DISTANCE criterion. These include:
                        #3-bond basepairs for which all 3 bonds have been detected, regardless of whether they carry a '0' or '1' tag. A 'yes' TAG will be appended to these basepairs. In addition, a 'D' (Distance) TAG
                        #2-bond basepairs which carry a '0', (i.e not involved in MULTIPLE CONTACTS). A 'yes' tag will be appended to these basepairs. In addition, a 'D' (Distance) TAG will be appended to these basepairs.
                        #3-bond basepairs for which only 2 bonds have been identified, which carry a '0', (i.e not involved in MULTIPLE CONTACTS). A 'yes' tag will be appended to these basepairs but only when MAX_CUTOFF_str[0:3] in CUTOFF_str[0:3]. In addition, a 'D' (Distance) TAG will be appended to these basepairs.
                                #The reason to keep these basepairs unassigned until the last run is that the missing third bond may be identified at MAX_CUTTOF, making their assignment stronger
                        #All other basepairs will be tagged as 'no'. This tag can be reversed by the next two CRITERIA
                        #After one round of "Initial Assignment" 'no'-tagged basepairs (only 2-bond basepairs) will be run through FUNCTION "REMOVAL". This will allow the identification of furhter basepairs through the "ELIMINATION" criterion.

#2)"ELIMINATION" criterion. As important as detecting basepairs possibly involved in MULTIPLE CONTACTS, is to remove the ones assigned as legitimate from this list ('control'). This is done by FUNCTION REMOVAL as follows.
                #The removal of newly assigned, legitimate basepairs would have left the length of these sub-lists intact but the strings corresponding to the newly assigned basepairs will carry now '[]'. Therefore if only 2 bases are left in some of these sublists after the removal of newly assigned basepairs from them, these leftover bases may identify a lone pair that has been removed of multiple interactions.
                #Any new basepair identified by this method will be appended a 'ELI'(Elimination) TAG at position (16 + ('run_number' - 1)*4). For STATISTICS.


#2) CONTINUITY OF HELICITY CRITERION
        #A 'no' string  can be reversed in Loop "Attempting to find additional basepairs by the CONTINUITY OF HELICITY CRITERION" if an adjacent basepair is found. This would  assume that the continuity of helical segments is a strong parameter while determining the likelihood of a given basepair
                #The involvement of many basepairs in MULTIPLE CONTACTS (see definition of MULTIPLE CONTACTS), which makes it impossible to accurately assign basepairs with less than 3 bonds, severely limits the ability of using distance data to accurately assign base pairs. At medium cutoffs about only half the basepairs can be identified.
#This proportion increases with lower cutoff distances (more restricitive search) but then, other basepairs with bonds longer than the cutoff are missed. Therefore, if no other criterion was used only 3-bond basepairs and a few other isolated basepairs, not involved in mutiple interactions could be identified by this program.
                #To improve this, I have written, "SUBSECTION Loop", containing the loop named "Attempting to find additional basepairs by the CONTINUITY OF HELICITY CRITERION". This loop will attempt to expand the initial assignment of legitimate basepairs by including basepairs for which only two bonds have been identified and tagged as 'no' in "SECTION: LOOKING FOR BASES INVOLVED IN MULTIPLE CONTACTS".
                #It works as follows:
                        #If an immediately adjacent legitimate basepair can be found for the current 'no'-tagged basepair candidate, the current basepair would be retagged as a 'yes', or legitimate, basepair. This is based on the idea that retagging the current basepair would expand an already established helical element
                        #The loop "Attempting to find additional basepairs by the CONTINUITY OF HELICITY CRITERION" will run until no more additional basepairs are found. This will ensure that every newly identified legitimate basepair can be used to increase the chances that other basepairs can be assigned as legitimate.
                        #Any additional basepairs identified in this section will be removed from the 'control' list carrying basepairs involved in MULTIPLE CONTACTS.
                        #NOTE:The actual loop is contained within FUNCTION "Loop", which is called from this section
                        #Any new basepair identified by the loop will be appended a 'MH'(Maximization of Helicity) TAG at position (16 + ('run_number' - 1)*4). For STATISTICS
                        #After one round of "Initial Assignment" 'no' tagged basepairs will be run through FUNCTION "REMOVAL". This will allow the identification of furhter basepairs through the "ELIMINATION" criterion.


#TAGS APPENDED: Every run, the following TAGS will be appended at the end of every line in MASTER_Basepairs_summary
               # 1) Position 13: TAG = 'none', globally added before MASTER_Basepairs_summary is searched for bases involved in MULTIPLE CONTACTS by comparing to list 'control'
               # 2) Position 13: TAG = 'yes' or 'no', which mean that the basepair candidate can be accurately assigned or not, depending on whether it is involved in MULTIPLE CONTACTS with other possible basepairs. For 3-bond basepairs for which all 3 bonds have been identified, a 'no' tag will be bypassed and the basepair will be considered unambiguously assigned in SECTION "Final assignment of Basepairs". The 'no' tag will then be reversed
               # 3) Position 14: TAG = can take the values '1' or '0', depending on whether the basepair was detected in MULTIPLE CONTACTS or not. For STATISTICS purposes only.
               # 5) Position 15: TAG = CUTOFF
               # 6) Position 16: TAG = can take the values 'D' (Distance), 'MH'(Maximization of Helicity), or 'ELI' (Elimination), to keep a record of how the basepair was identified (see SECTION: Final assignment of Basepairs). For STATISTICS only.
              # 7 and further: Position 17 through 20 (4 positions/run) will be appended in Run 2, and so on
                     #The actual position of the TAGs in every run can be calculated as: (13 + ('run_number' - 1)*4), in which the initial value for 'run_number' is = 1
#######SUBSECTION Initial Assignment

    convert = []
    new_basepairs = 0
    print   "##############################\n##### Initial Basepair Assignment. Run #", run_number_str, " Distance Cutoff = ", CUTOFF_str[0:3], "PREV_CUTOFF", PREV_CUTOFF, "\n##############################"
    transient = []
    if (MAX_CUTOFF_str[0:3] in CUTOFF_str[0:3]):
        print "MAX_CUTOFF_str[0:3] ", MAX_CUTOFF_str[0:3], "CUTOFF_str[0:3] ", CUTOFF_str[0:3]
    if CUTOFF_str[0:3] in MAX_CUTOFF_str[0:3]:
        print "CUTOFF_str[0:3] ", CUTOFF_str[0:3], "MAX_CUTOFF_str[0:3] ", MAX_CUTOFF_str[0:3]
    for i in range (len(MASTER_Basepairs_summary)):
        print "Basepairing scheme:", MASTER_Basepairs_schemes[i]
        for j in range (len(MASTER_Basepairs_summary[i])):
                if (len(MASTER_Basepairs_summary[i]) > 0) and ('D' not in MASTER_Basepairs_summary[i][j]) and ('MH' not in MASTER_Basepairs_summary[i][j]) and ('ELI' not in MASTER_Basepairs_summary[i][j]) and ('REMOVED' not in MASTER_Basepairs_summary[i][j]):

#basepairs with 3 bonds, all identified
                    if (MASTER_Basepairs_summary[i][j][1] == 3) and (MASTER_Basepairs_summary[i][j][6] == 3):
                        #if (('D' not in MASTER_Basepairs_summary[i][j]) and ('MH' not in MASTER_Basepairs_summary[i][j]) and ('ELI' not in MASTER_Basepairs_summary[i][j]) and ('REMOVED' not in MASTER_Basepairs_summary[i][j])):
                            list = [MASTER_Basepairs_summary[i][j][2], MASTER_Basepairs_summary[i][j][3], MASTER_Basepairs_summary[i][j][4], MASTER_Basepairs_summary[i][j][5]]
                            convert = CONVERT(list)
                            print  "----->High confidence 3-bond basepair formed by residues " + convert[1] + convert[0] + " and " + convert[3] + convert[2]
                            length_1 = len(MASTER_Basepairs_summary[i][j]) - 1
                            MASTER_Basepairs_summary[i][j].append('D') #Distance 'D' TAG
                            MASTER_Basepairs_summary[i][j].append(CUTOFF_str[0:3])
                            transient = [MASTER_Basepairs_summary[i][j][2], MASTER_Basepairs_summary[i][j][4]]
                            positions.append(transient)
                            new_basepairs = new_basepairs + 1

                            if ('no' in MASTER_Basepairs_summary[i][j]):
                                print "    The bases in this basepair may be participating in an interaction with other bases. However, the detection of 3 bonds, makes its assignment as a 3-bond basepair highly probable."

#basepairs with 3 bonds for which only 2 bonds were identified
                    elif ((MASTER_Basepairs_summary[i][j][1] == 3) and (MASTER_Basepairs_summary[i][j][6] == 2) and (MAX_CUTOFF_str[0:3] in CUTOFF_str[0:3])):
#These basepairs are excluded from 'DISTANCE' and 'ELIMINATION' criteria (MAX_CUTOFF_str[0:3] in CUTOFF_str[0:3]), so they should not arrive here with either 'D' or 'ELI'
                        if ('yes' in  MASTER_Basepairs_summary[i][j]):
                            list = [MASTER_Basepairs_summary[i][j][2], MASTER_Basepairs_summary[i][j][3], MASTER_Basepairs_summary[i][j][4], MASTER_Basepairs_summary[i][j][5]]
                            convert = CONVERT(list)
                            print "----->At least 2 bonds identified out of 3 expected for the base pair formed by residues " + convert[1] + convert[0] + " and " + convert[3] + convert[2]
                            MASTER_Basepairs_summary[i][j].append('D') #Distance 'D' TAG
                            MASTER_Basepairs_summary[i][j].append(CUTOFF_str[0:3])
                            transient = [MASTER_Basepairs_summary[i][j][2], MASTER_Basepairs_summary[i][j][4]]
                            positions.append(transient)
                            new_basepairs = new_basepairs + 1

#basepairs with 2 bonds
                    elif (MASTER_Basepairs_summary[i][j][1] == 2) and (MASTER_Basepairs_summary[i][j][6] == 2):
                        #if (('D' not in MASTER_Basepairs_summary[i][j]) and ('MH' not in MASTER_Basepairs_summary[i][j]) and ('ELI' not in MASTER_Basepairs_summary[i][j]) and ('REMOVED' not in MASTER_Basepairs_summary[i][j])):
                            if ('yes' in  MASTER_Basepairs_summary[i][j]):
                                list = [MASTER_Basepairs_summary[i][j][2], MASTER_Basepairs_summary[i][j][3], MASTER_Basepairs_summary[i][j][4], MASTER_Basepairs_summary[i][j][5]]
                                convert = CONVERT(list)
                                print "----->2-bond basepair identified formed by residues " + convert[1] + convert[0] + " and " + convert[3] + convert[2]
                                MASTER_Basepairs_summary[i][j].append('D') #Distance 'D' TAG
                                MASTER_Basepairs_summary[i][j].append(CUTOFF_str[0:3])
                                transient = [MASTER_Basepairs_summary[i][j][2], MASTER_Basepairs_summary[i][j][4]]
                                positions.append(transient)
                                new_basepairs = new_basepairs + 1
    print   "##### END of Initial Basepair Assignment. Run #", run_number_str, " Distance Cutoff = ", CUTOFF_str[0:3], "\n##############################"

#The newly identified basepairs should be removed from 'control' (list of bases involved in MULTIPLE CONTACTS)
    if new_basepairs > 0 and control_control == 1:
        print "\n##############################\n##### Searching for additional basepairs by the ELIMINATION criterion. Run # ", run_number_str, " Distance Cutoff = ", CUTOFF_str[0:3]
        mode_removal = 0
        removal = REMOVAL(MASTER_Basepairs_summary, control, convert[0], convert[2], CUTOFF, mode_removal, positions)
        MASTER_Basepairs_summary = removal[0]
        control = removal[1]
        new_basepairs = removal[2]

#Printing output.
        a = str(new_basepairs)
        line1 = "##### Number of additional basepairs identified by ELIMINATION = " + a + "  Run # " + run_number_str +  " Distance Cutoff = " + CUTOFF_str[0:3]
        print line1
        print "############################## "
        line2 = "AFTER Bulk Removal: Run # " + run_number_str + " Distance Cutoff = " + CUTOFF_str[0:3] #This line will be outputed in FUNCTION "SEQUENTIAL_READOUT"
        print_control = 0
        SEQUENTIAL_READOUT(MASTER_Basepairs_summary, MASTER_Basepairs_schemes, control, print_control, line2)
    elif new_basepairs == 0 and control_control == 1:
        print "##############################\n NO newly identified basepairs to be removed from the list BASES WITH POSSIBLE MULTIPLE CONTACTS.\n ", "Run # ", run_number_str, " Distance Cutoff = ", CUTOFF_str[0:3]
        print "############################## "
    else:
        print "##############################\n NO BASES WITH POSSIBLE MULTIPLE CONTACTS have been identified.\nAssignement of furhter additional basepairs by the ELIMINATION criterion not possible\n ", "Run # ", run_number_str, " Distance Cutoff = ", CUTOFF_str[0:3]
        print "############################## "

######SUBSECTION Initial Assignment

######SUBSECTION Loop
#######Loop "Attempting to find additional basepairs by the CONTINUITY OF HELICITY CRITERION"
#see more info on this loop above

#FIRST LOOP RUN, searching for basepairs that would increase helicity

    if control_control == 1:
        print "\n##############################\n##### SEARCHING FOR ADDITIONAL BASEPAIRS BY USING THE CONTINUITY OF HELICITY CRITERIUM /nRun #", run_number_str, " Distance Cutoff = ", CUTOFF_str[0:3]

        loop_output = []
        loop_control = 1 #Allow multiple runs of the loop

#calling FUNCTION "Loop"
        loop_output = loop(MASTER_Basepairs_summary, control, loop_control, CUTOFF, positions)
        MASTER_Basepairs_summary = loop_output[0]
        control = loop_output[1]
#Calling FUNCTION "SEQUENTIAL_READOUT"
        line2 = "AFTER SEARCHING FOR ADDITIONAL BASEPAIRS BY USING THE CONTINUITY OF HELICITY CRITERIUM: Run # " + run_number_str + " Distance Cutoff = " + CUTOFF_str[0:3]
        print_control = 0
        SEQUENTIAL_READOUT(MASTER_Basepairs_summary, MASTER_Basepairs_schemes, control, print_control, line2)
        print "\n##############################"
    else:
        print "\n##############################\nBYPASSING THE SEARCH FOR ADDITIONAL BASEPAIRS BY USING THE CONTINUITY OF HELICITY CRITERIUM /nRun #", run_number_str, " Distance Cutoff = ", CUTOFF_str[0:3], "\n##############################"

################################################
#       FUNCTION "Program". Section: Final assignment of Basepairs       END
################################################

################################################
################################################
#FUNCTION "program"                     END
################################################
################################################

##########################################################################
##########################################################################
##########################################################################
#SECTION FUNCTIONS                                                 END
##########################################################################
##########################################################################
##########################################################################


##########################################################################
##########################################################################
##########################################################################
#SECTION HELICITY                        BEGINNING
##########################################################################
##########################################################################
##########################################################################
#Searches for matches to a potential basepair candidate within MASTER_Basepairs_summary (named LIST within this FUNCTION)
def HELICITY(LIST, a0, a1, a2, a3, where_from):
                        a0 = str(a0)
                        a1 = str(a1)
                        a2 = str(a2)
                        a3 = str(a3)
                        line1 = "            Searching for possible basepairs between residues " + a0 + " and " +  a1

                        print line1
                        if ('iiiiiiiiiiiiiii' in a2) and ('iiiiiiiiiiiiiii' in a3):
                              line1 = line1 + a2 + " and " + a3 + "\n"

#Called from FUNCTION "Loop"
#LIST = MASTER_Basepairs_summary
                        match_found1 = []
                        match_found2 = []
                        match_position = []
                        for r in range (len(LIST)):
                            for s in range (len(LIST[r])):
                                m=re.search('\d+', LIST[r][s][2])
                                match1 = m.group()
                                m=re.search('\d+', LIST[r][s][4])
                                match2 = m.group()
#                                print "LIST[r][s]", LIST[r][s], "\n    match1", match1, "match2", match2, "a0", a0, "a1", a1, "a2", a2, "a3", a3
                                if (('Loop' in where_from) and ('REMOVED' not in LIST[r][s]) and ('yes' in LIST[r][s]) and ((match1 in a0) and (match2 in a1) or (match2 in a0) and (match1 in a1))):
                                   match_found1 = [a0, a1]
                                   match_position = [r, s]
#                                   print "match found: a0", a0, "a1", a1
                                elif (('Loop' in where_from) and ('REMOVED' not in LIST[r][s]) and ('yes' in LIST[r][s]) and ((match1 in a2) and (match2 in a3) or (match1 in a3) and (match2 in a2))):
                                   match_found2 = [a2, a3]
                                   match_position = [r, s]
#                                   print "match found: a2", a2, "a3", a3
                                elif ('REMOVAL' in where_from) and ('REMOVED' not in LIST[r][s]) and ('no' in LIST[r][s]) and ((match1 in a0) and (match2 in a1) or (match2 in a0) and (match1 in a1)):
                                   match_found1 = [a0, a1]
                                   match_position = [r, s]
#                                   print "where_from ", where_from, "LIST[r][s] ", LIST[r][s],  "\na0", a0, "a1", a1, "match1 ", match1, "match2 ", match2

                        return[match_found1, match_found2, match_position]
##########################################################################
##########################################################################
##########################################################################
#SECTION HELICITY                                 END
##########################################################################
##########################################################################
##########################################################################




##########################################################################
##########################################################################
##########################################################################
#SECTION MAIN                                    BEGINNING
##########################################################################
##########################################################################
##########################################################################

##########################################################################
#Sub SECTION Calling run(args)
##########################################################################
#FUNCTION run(args) was created by Ralf W. Grosse-Kunstle and provides a list of ATOM-to-ATOM distances calculated from a .pdb file.

First_List = run(args=sys.argv[1:]) #Will carry the output of FUNCTION run(args)
pdb_file = First_List[0] #Name of the file read in
First_List[0] = []
##########################################################################
#Sub SECTION Calling run(args)
##########################################################################

##########################################################################
#Sub SECTION file name
##########################################################################
import os
dir =  os.getcwd()

import re
ma=re.compile(r'\/+')
iterator = ma.finditer(pdb_file)
aaa=[]
for match in iterator:
   aaa.append(match.span())
directory = pdb_file[0:aaa[len(aaa)-1][1]]
FILE = pdb_file[aaa[len(aaa)-1][1]:len(pdb_file)]
mm=re.search(r'\.', FILE)
START = mm.start()
pdb_file_main = FILE[0:START]

##########################################################################
#Sub SECTION file name
##########################################################################

##########################################################################
#Sub SECTION MAIN LOOP                        BEGINNING
##########################################################################
run_number = 1 #In addition to signaling the Run #, this variable will allow differential execution of certain events depending on whether the program is in 'first run' at the lowest CUTOFF, or later runs at higher CUTOFFs
run_cutoff = [] #Will register all CUTOFF values used in each run

CUTOFF = 3.0
CUTOFF_str = repr(CUTOFF) #Converts 'float' to 'str'
run_cutoff.append(run_number)
run_cutoff.append(CUTOFF_str[0:3])
run_cutoff_LISTS = [] #Will carry all the assigned basepairs, ordered by CUTOFF value at which they have been assigned. Used in STATISTICS
MAX_CUTOFF = 3.6
PREV_CUTOFF = 0.0 #Initializing
increment = 0.2
        #ARGUMENTS
                #First_List = Contains the filtered output from FUNCTION "run(args)". Listof ATOM to ATOM distances below a rough cutoff of 5A
                #CUTOFF = Defined above. CUTOFF to assign ATOM-to-ATOM distances that are compatible with basepair formation. CUTOFF has a Dynamic value, increasing by the value of 'increment' at every new execution of FUNCTION "Program"
                #pdb_file_main = main part of the input .pdb file withou the extension
                #file = name of the .pml file for PYMOL
                #file1 = name of the output file for the run
positions = [] #Will store the base #s of the basepairs identified for their easy removal in FUNCTION REMOVAL
while CUTOFF < MAX_CUTOFF:
    if run_number > 1:
        CUTOFF = CUTOFF + increment
        CUTOFF_str = repr(CUTOFF) #Converts 'float' to 'str'
        run_cutoff.append(run_number)
        run_cutoff.append(CUTOFF_str[0:3])
#    else:
    control = [] #The list of bases involved in MULTIPLE CONTACTS will be loaded into 'control' by function 'CONTROL'. See FUNCTION Section: LOOKING FOR BASES INVOLVED IN MULTIPLE CONTACTS:

    print "\n############################################################################################\n############################################################################################\n###### Run # ", run_number, "DISTANCE CUTOFF = ", CUTOFF, "\n############################################################################################\n############################################################################################\n"
    run_number_str = to_string(run_number)
    MAX_CUTOFF_str = str(MAX_CUTOFF) #Converts 'float' to 'str'

#Executing FUNCTION "program", main part of the program:
    #It started as the program itself but had to be converted to a function to allow its repetitive iteration
    program(First_List, MASTER_Basepairs_summary, CUTOFF, pdb_file_main, control, run_number_str, CUTOFF_str, run_number, positions)

    run_number = run_number + 1
    PREV_CUTOFF = CUTOFF
##########################################################################
#Sub SECTION MAIN LOOP                          END
##########################################################################

##########################################################################
##########################################################################
##########################################################################
#SECTION MAIN                                   END
##########################################################################
##########################################################################
##########################################################################


##########################################################################
##########################################################################
##########################################################################
#SECTION DETECTION OF FALSE BASEPAIRS                  BEGININNG
##########################################################################
##########################################################################
##########################################################################

#I have observed that in addition to predicting basepairs, the CONTINUITY OF HELICITY CRITERION could be used to detect some false onse. For example, I have observed that false 3-bond basepairs could be erroneously assigned by bases that are closely stacked to real basepairs. This is not the situation for 2-bond basepairs candidates and 3-bond basepair candidates, for which only 2 bonds have been detected, because these are excluded from being falsely assigned by the MULTIPLE INTERACTION CRITERION. Interestinlgy, I have observed that after assignment of additional basepairs by the CONTINUITY OF HELICITY CRITERION, new 2-bond basepairs may be assigned which share a base with the false basepair. Therefore, a 'MH'-tagged basepair could be used to detect a false 3-bond basepair, with which it shares a base, and to remove it from the final basepair list. However, before this is done, the non-'MH'-tagged basepair should be run through checked for its possible involvement in a helical element (see below). In addition, this observation also indicates that in other cases, in which the exclusion of a erroneously assigned basepair is not as stright forward (for example, 2 D-tagged, 3-bond basepairs that may share one base), subjecting the two basepairs to the "CONTINUITY OF HELICITY CRITERION" may be used to tell the legitimate basepair from the false one.


#I have also observed that some false Distance-identified 2-bond basepairs can escape the MULTIPLE INTERACTION CRITERION if they lie adjacent to a 3-bond basepair for which all bonds have been identified. This is not surprising, as the latter basepairs are excluded from re-entering the pool of bases involved in  MULTIPLE CONTACTS. In this case, running the CONTINUITY OF HELICITY CRITERION should identify the false basepair. However, the possibility of allowing these basepairs to enter the pool of bases involved in  MULTIPLE CONTACTS should be considered.
        #Note 021010: Indeed, allowing 3-bond basepairs for which all bonds have been identified to enter the  pool of bases involved in  MULTIPLE CONTACTS, prevents some of these bases from being erroneoulsy assigned as legitimate basepairs.


####FUNCTIONING
#In the MAIN part of the section MASTER_Basepairs_summary_bak will be compared to MASTER_Basepairs_summary to detect false basepairs as described above. Several criteria will be cyclically used at different levels of discrimination stringency to decide which one of the debated basepairs is LEGITIMATE and which one is FALSE. Initially a list of base-sharing basepairs is defined. There are 2 categories:
        #Basepairs which share a single base. These are the most abundant. They are removed from the list first.
        #Basepairs which share both bases. In this case the 2 bases clearly form a basepair, only the right geometry has to be inferred. These are removed from the list once the basepairs which share a single base have been processed.


        #"C1-C1 and P-P DISTANCE CRITERION": Performed only at HIGH STRINGENCY and before any of the other criteria is used. Will allow initial assignement of FALSE basepairs independently of the rest of the loop.
                #Performed in SECTION FUNCTION "C1_C1 DISTANCE CRITERION". As of 050310 version, I have added statistical information regarding the C1-C1 and P-P distance. The information was empirically generated by the program and then fed back into into it by means of the list "MASTER_Basepairs_bonds". Historically, while this information was not necessary for the prediction of basepairs in tRNA (4TNA.pdb), 16S rRNA was more of a challenge, in particular regions in which the structure is a little bit disordered, the program has problems telling real basepairs from false ones due to the noise brought about by stacking interactions. The quality of this information should be improved as more structures are processed by the program, specially large ones. IMPORTANT: only really well refined structures should be used to get statistical data, as this information is highly sensitive to the quality of the structure!!!!!
        #"CONTINUITY_HELICITY_CRITERION CRITERION": Initially this criterion was the only one used to identify FALSE basepairs and worked well enough with tRNA. However, tests with 16S rRNA made it clear that this criterion by itself was not enough.
                #The debated basepairs are sent to SECTION FUNCTION "CONTINUITY_HELICITY_CRITERION" to check whether they could form part of a helical element. The debated basepairs will be given a number from "0" to "2", based on the number of helical elements that are continuous with the proposed basepair.
                        #Score
        #"CUTOFF_CRITERION": Performed in SECTION FUNCTION "CUTOFF_CRITERION". Not very valuable and possibly targeted for removal after implementing the "C1-C1 DISTANCE CRITERION".
        #"BOND_CRITERION": Perofrmed in SECTION FUNCTION "BOND_CRITERION". Not very valuable and possibly targeted for removal after implementing the "C1-C1 DISTANCE CRITERION".

#At HIGH STRINGENCY, basepairs with this distance falling beyond the C1-C1 distance average +/- standard deviation are not considered further. Other criteria not considered at this stringency
#At LOW STRINGENCY, the C1-C1 distance average +/- standard deviation values could be used to weigh the rest of the scores, coming from the other discrimination criteria


#The recorded score after submitting the debated basepairs to all the above criteria determines which one of the them is defined as LEGITIMATE and which one as FALSE. A 'LEGITIMATE' tag will be appended to this basepair and a 'FALSE' tag to the other
             #If both basepairs have equal scores associated to them, the 'yes' tag of both will be replaced with a 'no' and an 'UNDETERMINED' tag will be appended to both

#################################
#####SECTION FUNCTION "CONTINUITY_HELICITY_CRITERION"          BEGINNING
###### The following lines are adapted from FUNCTION LOOP. It was too complicated to run again the program through FUNCTION LOOP, so i have adapted the relevant lines of FUNCTION LOOP to check the amount of helix continuity afforded by each of the possible false basepairs.
#Called below

def CONTINUITY_HELICITY_CRITERION(MASTER_Basepairs_summary_line, MASTER_Basepairs_summary, CUTOFF):
                        m=re.search('\d+', MASTER_Basepairs_summary_line[2])
                        search1 = m.group()
                        SEARCH1 = to_int(search1) #convert string to integers
                        m=re.search('\d+', MASTER_Basepairs_summary_line[4])
                        search2 = m.group()
                        SEARCH2 = to_int(search2) #convert string to integers
                        list = [MASTER_Basepairs_summary_line[2], MASTER_Basepairs_summary_line[3], MASTER_Basepairs_summary_line[4], MASTER_Basepairs_summary_line[5]]
                        convert = CONVERT(list)
                        add = ADD(SEARCH1, SEARCH2)
                        line = "\n...using the CONTINUITY OF HELICITY CRITERION criterion for " + convert[0] + convert[1] + " and " + convert[2] + convert[3]
                        print line
                        line1 = "\nChecking whether the basepair formed by residues " + convert[0] + convert[1] + " and " + convert[2] + convert[3] + " could be involved in a helical element\n    Searching for the following possible basepair candidates: " + add[0] + " and " + add[1] + " or " + add[2] + " and " + add[3]
                        print line1

#the whole MASTER_Basepairs_summary is searched for the existence of either one of the adjacent basepairs calculated above. Done in FUNCTION HELICITY
                        where_from = 'Loop' #Used in SECTION HELICITY to determine what to do whether the program is coming from 'Loop' or 'CONTINUITY_HELICITY_CRITERION', or from 'REMOVAL'
                        helicity = HELICITY(MASTER_Basepairs_summary, add[0], add[1], add[2], add[3], where_from)
                        count_helicity = 0 #Control output in this region
                        for t in range (0,2):
                            if helicity[t] != []:
                                count_helicity = count_helicity + 1  #Signals entering this 'if'
                                str(helicity[t][0])
                                str(helicity[t][1])
#                                count_found = count_found + 1 #Will keep tract of how many basepairs are found. Serves also to keep the loop going

#One or the 2 adjacent basepairs is found
                                if count_helicity == 1: #One of the 2 possible adjacent basepairs is found. This is enought to assign the basepair whose assignemnt is in question by 'MH' criterion
#                                    MASTER_Basepairs_summary_line.append('MH')  #'MH'Maximization of Helicity TAG
                                    line1  = "     #####\n     !!!!Found an adjacent basepair to residues " + convert[0] + convert[1] + " and "  + convert[2] + convert[3] +  "\n     This basepair is formed by residues " + helicity[t][0] + " and " +  helicity[t][1] + "\n     #####\n"
                                    print line1
                                if count_helicity == 2: #The second adjacent basepair has been found. The basepair whose assignemnt is in question has already bee assigned when 'count_helicity == 1'. The program will only output the new finding
                                    line1  = "     #####\n     !!!!Found a second adjacent basepair to residues " + convert[0] + convert[1] + " and "  + convert[2] + convert[3] +  "\n     This basepair is formed by residues " + helicity[t][0] + " and " +  helicity[t][1] + "     \n     #####\n"
                                    print line1
#Neither one of the possible adjacent basepairs was found by the CONTINUITY OF HELICITY CRITERION
                        if count_helicity == 0:
                            line1  = "                     Neither basepair was found!!!!!"
                            print line1
                        return [MASTER_Basepairs_summary_line, count_helicity]
#####SECTION FUNCTION "CONTINUITY_HELICITY_CRITERION"                END
#################################


#################################
#####SECTION FUNCTION "CUTOFF_CRITERION"          BEGINNING
def CUTOFF_CRITERION(LIST1, LIST2, A0, A1, A2, A3, A4, A5, A6, A7):
#LIST1 = MASTER_Basepairs_summary[i][j]
#LIST2 = MASTER_Basepairs_summary[l][m]
#A0 = convert[0]
#A1 = convert[1]
#A2 = convert[2]
#A3 = convert[3]
#A4 = convert[4]
#A5 = convert[5]
#A6 = convert[6]
#A7 = convert[7]
    score = [0, 0]
    print "\n...using the CUTOFF criterion"
    #The CUTOFF string has to be converted back into an 'float' number
    float_1 = float(LIST1[len(LIST1) - 1][0]) + float(LIST1[len(LIST1) - 1][2])/10
    float_2 = float(LIST2[len(LIST2) - 1][0]) + float(LIST2[len(LIST2) - 1][2])/10
    line = "\n     The basepair formed by residues " + A0 + A1 + ":" + A2 + A3 + " and with a " + LIST1[0] + " geometry was assigned at a CUTOFF value of " + LIST1[len(LIST1) - 1]  + " whereas the basepair formed by residues " + A4 + A5 + ":" + A6 + A7 + "and with a " + LIST2[0] + " geometry was assigned at a CUTOFF value of " + LIST2[len(LIST2) - 1]
    if float_1 < float_2:
        print line
        line1 = "          therefore the " + A0 + A1 + ":" + A2 + A3 + " basepair, with a " + LIST1[0] + " geometry, fares better than the " + A4 + A5 + ":" + A6 + A7 + " basepair, with a " +  LIST2[0] + " geometry by the CUTOFF criterion"
        print line1
        score = [2, 0]

    elif float_1 > float_2:
        print line
        line1 = "          therefore the " + A4 + A5 + ":" + A6 + A7 + " basepair, with a " + LIST2[0] + " geometry, fares better than the " + A0 + A1 + ":" + A2 + A3 + " basepair, with a " +  LIST1[0] + " geometry by the CUTOFF criterion"
        print line1
        score = [0 , 2]
    else:
        print "\n     Both basepairs have been assigned at the same CUTOFF value of " + LIST2[len(LIST2) - 1]

    return score

#####SECTION FUNCTION "CUTOFF_CRITERION"          END
#################################

#################################
#####SECTION FUNCTION "BOND_CRITERION"          BEGINNING
def BOND_CRITERION(LIST1, LIST2, A0, A1, A2, A3, A4, A5, A6, A7):
#LIST1 = MASTER_Basepairs_summary[i][j]
#LIST2 = MASTER_Basepairs_summary[l][m]
#A0 = convert[0]
#A1 = convert[1]
#A2 = convert[2]
#A3 = convert[3]
#A4 = convert[4]
#A5 = convert[5]
#A6 = convert[6]
#A7 = convert[7]
     score = [0, 0]
                #number of bonds criterion. The larger this number, the more chances that the basepair is correctly assigned
     print "\n...using the BOND criterion"
     if LIST1[6] > LIST2[6]:
         line = "     therefore the " + A0 + A1 + ":" + A2 + A3 + " basepair fares better under the 'NUMBER-OF-BONDS CRITERION' than the " + A4 + A5 + ":" + A6 + A7 + "basepair."
         score[0] = 1
         score[1] = 0
     elif LIST1[6] < LIST2[6]:
         line = "     therefore the " + A4 + A5 + ":" + A6 + A7 + " basepair fares better under the 'NUMBER-OF-BONDS CRITERION' than the " + A0 + A1 + ":" + A2 + A3 + "basepair."
         score[1] = 1
         score[0] = 0
     else:
         line = "     therefore the 'NUMBER-OF-BONDS CRITERION' cannot be used to sort out basepair legitimacy"
     A = str(LIST1[6])
     B = str(LIST2[6])
     line1 = A + " bonds have been identified for the " + A0 + A1 + ":" + A2 + A3 + " basepair and " + B + " bonds have been identified for the " + A4 + A5 + ":" + A6 + A7 + " basepair."
     print line1
     print line

     return score

#####SECTION FUNCTION "BOND_CRITERION"          END
#################################

#################################
#####SECTION FUNCTION "DECISSIONS"              BEGINNING
def DECISSION(LIST1, LIST2, A0, A1, A2, A3, A4, A5, A6, A7, removed, undetermined, two_base_shared):
#LIST1 and LIST2 can either be MASTER_Basepairs_summary[i][j] or MASTER_Basepairs_summary[l][m] depending on their scores
#A0 = convert[0]
#A1 = convert[1]
#A2 = convert[2]
#A3 = convert[3]
#A4 = convert[4]
#A5 = convert[5]
#A6 = convert[6]
#A7 = convert[7]
    if undetermined == 'n':
        if two_base_shared == "n":
           line = "The basepair formed by residues " + A4 + A5 + ":" + A6 + A7 + " and with geometry " + LIST2[0] + " will be assigned as 'FALSE'. Please inspect this basepair and its surroundings for the presence of additional misassigned basepairs"
           LIST2.append('FALSE')
        else:
           line = "The basepair formed by residues " + A4 + A5 + ":" + A6 + A7 + " and with geometry " + LIST2[0] + " will be assigned as 'WRONG'. Please inspect this basepair and its surroundings for the presence of additional misassigned basepairs"
           LIST2.append('WRONG')
        print line
        line = "The basepair formed by residues " + A0 + A1 + ":" + A2 + A3 + " and with geometry " + LIST1[0] + " will be assigned as 'LEGITIMATE'"
        print line
        print "################"
        LIST1[19:19] = ['LEGITIMATE']
    else:
        line = "\nBasepairs, " + convert[0] + convert[1] + ":" + convert[2] + convert[3] + " with geometry " + MASTER_Basepairs_summary[i][j][0] + " and " + convert[4] + convert[5] + ":" + convert[6] + convert[7] +  " with geometry " + MASTER_Basepairs_summary[l][m][0] + " cannot be deconvoluted. Their assignment as legitimate basepairs remains 'UNDETERMINED'"
        print line
        print "################"
        LIST1.append('UNDETERMINED')
        LIST2.append('UNDETERMINED')
    removed = removed + 1

    return [LIST1, LIST2, removed]
#####SECTION FUNCTION "DECISSIONS"              END
#################################

######################################################
#####MAIN part of SECTION DETECTION OF FALSE BASEPAIRS
#MASTER_Basepairs_summary_bak = MASTER_Basepairs_summary
print "SEARCHING FOR ERRONEOUSLY ASSIGNED BASEPAIRS"

CONTINUITY_HELICITY_CRITERION_output = [] #Will collect the output of 'FUNCTION CONTINUITY_HELICITY_CRITERION'
while_control = 'y' #Will end the next 'while' loop when == 'n'
run_while = 1
stringency = 'HIGH' #Will signal whether the next 'while' loop will be executed with 'HIGH' or 'LOW' stringency
single_base_shared = 'n' #When 'y', it will signal that deconvolution of basepairs sharing one base can start
two_base_shared = 'n' #When 'y', it will signal that deconvolution of basepairs sharing two bases can start. Requires prior deconvolution of basepairs sharing a single base

###############
#'FALSE-BASEPAIR-ASSIGNMENT LOOP'             BEGINNING
MODE = "Single base sharing"
while while_control == 'y':

    removed = 0 #This variable will keep track of how many basepairs sharing bases are deconvoluted in this loop
    if two_base_shared == 'y':
       line = "deconvoluting candidate basepairs sharing both bases\n"
       run_while = 1
    else:
       line = "deconvoluting candidate basepairs sharing one base\n"
    print "#######################################################\n#####SEARCHING FOR ERRONEOUSLY ASSIGNED BASEPAIRS: STRINGENCY ", stringency, " Run =", run_while, "\nMODE =", MODE, "\n#######################################################\n"
#, line, "\nremoved", removed, " while_control", while_control, " two_base_shared", two_base_shared, " single_base_shared", single_base_shared, " stringency", stringency
    for i in range (len(MASTER_Basepairs_summary)):
        if (len(MASTER_Basepairs_summary[i]) > 0):
            for j in range (len(MASTER_Basepairs_summary[i])):


                if (('REMOVED' not in MASTER_Basepairs_summary[i][j]) and ('FALSE' not in MASTER_Basepairs_summary[i][j]) and ('WRONG' not in MASTER_Basepairs_summary[i][j]) and ('UNDETERMINED' not in MASTER_Basepairs_summary[i][j])) and (('D' in MASTER_Basepairs_summary[i][j]) or ('MH' in MASTER_Basepairs_summary[i][j]) or ('ELI' in MASTER_Basepairs_summary[i][j])):
                    for l in range (len(MASTER_Basepairs_summary)):
                        if (len(MASTER_Basepairs_summary[l]) > 0):
                            for m in range (len(MASTER_Basepairs_summary[l])):
                                 if (('REMOVED' not in MASTER_Basepairs_summary[i][j]) and ('FALSE' not in MASTER_Basepairs_summary[l][m]) and ('WRONG' not in MASTER_Basepairs_summary[l][m]) and ('UNDETERMINED' not in MASTER_Basepairs_summary[l][m])) and (('D' in MASTER_Basepairs_summary[l][m]) or ('MH' in MASTER_Basepairs_summary[l][m]) or ('ELI' in MASTER_Basepairs_summary[l][m])):
                                     list = [MASTER_Basepairs_summary[i][j][2], MASTER_Basepairs_summary[i][j][3], MASTER_Basepairs_summary[i][j][4], MASTER_Basepairs_summary[i][j][5], MASTER_Basepairs_summary[l][m][2], MASTER_Basepairs_summary[l][m][3], MASTER_Basepairs_summary[l][m][4], MASTER_Basepairs_summary[l][m][5]]
                                     convert = CONVERT(list)
                                     SCORE = [[],[],[]] #To keep track of the performance of the debated basepairs in each of the three CRITERIA that will be used to sort out their legitimacy. FORMAT: 'SCORE[0] will be used for CONTINUITY_HELICITY_CRITERION', 'SCORE[1] will be used for CUTOFF_CRITERION', and 'SCORE[2] will be used for BOND_CRITERION'
                                     score = [[],[]] #'score[0]' and 'score[1]' will carry a value representative of how the debated basepairs fare under all criteria used in the loop+
                                     single_base_shared = 'n' ##When 'y', it will signal that deconvolution of basepairs sharing one base can start
                                     if (MASTER_Basepairs_summary[i][j][2] in MASTER_Basepairs_summary[l][m]) and (MASTER_Basepairs_summary[i][j][4] not in MASTER_Basepairs_summary[l][m]):
                                        single_base_shared = 'y'
#                                        print "single_base_shared = 'y' MASTER_Basepairs_summary[i][j]", MASTER_Basepairs_summary[i][j], "\nMASTER_Basepairs_summary[l][m]", MASTER_Basepairs_summary[l][m]
                                     if (MASTER_Basepairs_summary[i][j][2] not in MASTER_Basepairs_summary[l][m]) and (MASTER_Basepairs_summary[i][j][4] in MASTER_Basepairs_summary[l][m]):
                                        single_base_shared = 'y'
#                                        print "single_base_shared = 'y' MASTER_Basepairs_summary[i][j]", MASTER_Basepairs_summary[i][j], "\nMASTER_Basepairs_summary[l][m]", MASTER_Basepairs_summary[l][m]

######CASE 1: Only one base is shared. Most of the times this is the case.
#In this case there are three criteria that could be used to sort out which of the two is 'LEGITIMATE' and which one is 'FALSE'.
                                     if single_base_shared == 'y':
                                        line = "\n################\nBasepair, " + convert[0] + convert[1] + ":" + convert[2] + convert[3] + " with geometry " + MASTER_Basepairs_summary[i][j][0] + " and basepair " + convert[4] + convert[5] + ":" + convert[6] + convert[7] + " with geometry " + MASTER_Basepairs_summary[l][m][0] + " share one base!!!!!!\nAttempting to deconvolute these two basepairs with " + stringency + " stringency by"
                                        print line

        #CRITERION #1: CONTINUITY_HELICITY_CRITERION: Checking helicity of the basepair. This is done in 'FUNCTION CONTINUITY_HELICITY_CRITERION', where the basepair will be assigned a number between '0' and '2' depending on the number of helical elements that are continuous with it. HIGHEST PRIORITY!!
                                        CONTINUITY_HELICITY_CRITERION_output = CONTINUITY_HELICITY_CRITERION(MASTER_Basepairs_summary[i][j], MASTER_Basepairs_summary, CUTOFF)
                                        MASTER_Basepairs_summary[i][j] = CONTINUITY_HELICITY_CRITERION_output[0]
                                        SCORE[0].append(CONTINUITY_HELICITY_CRITERION_output[1] * 4)
                                        CONTINUITY_HELICITY_CRITERION_output = CONTINUITY_HELICITY_CRITERION(MASTER_Basepairs_summary[l][m], MASTER_Basepairs_summary, CUTOFF)
                                        MASTER_Basepairs_summary[l][m] = CONTINUITY_HELICITY_CRITERION_output[0]
                                        SCORE[0].append(CONTINUITY_HELICITY_CRITERION_output[1] * 4)
                                        print "Number of helical elements continuous with basepair, ", convert[0], convert[1], ":", convert[2], convert[3], " = ", int(SCORE[0][0]/4)
                                        print "Number of helical elements continuous with basepair, ", convert[4], convert[5], ":", convert[6], convert[7], " = ", int(SCORE[0][1]/4)

#                                        print "AFTER CONTINUITY_HELICITY_CRITERION SCORE", SCORE, "\nscore", score
        #CRITERION #2: CUTOFF_CRITERION: the smaller the CUTOFF value, the more chances that the basepair is correctly assigned. MEDIUM PRIORITY.
                                   #The cutoff string is at MASTER_Basepairs_summary[i][j][len(MASTER_Basepairs_summary[i][j]) - 1] and MASTER_Basepairs_summary[l][m][len(MASTER_Basepairs_summary[l][m]) - 1]

                                        SCORE[1] = CUTOFF_CRITERION(MASTER_Basepairs_summary[i][j], MASTER_Basepairs_summary[l][m], convert[0], convert[1], convert[2], convert[3], convert[4], convert[5], convert[6], convert[7])
#                                        print "AFTER CUTOFF_CRITERION SCORE", SCORE, "\nscore", score
        #CRITERION #3: BOND_CRITERION: The larger the number of bonds between two bases, the more chances that the basepair is correctly assigned. Should have the lowest priority. LOWEST PRIORITY.
                                        SCORE[2] = BOND_CRITERION(MASTER_Basepairs_summary[i][j], MASTER_Basepairs_summary[l][m], convert[0], convert[1], convert[2], convert[3], convert[4], convert[5], convert[6], convert[7])

####CASE 2: both bases are shared, the distances between ATOMS for these bases is compatible with more than 1 basepair geometry. SHOULD ONLY BE PERFORMED AFTER ALL BASEPAIRS WHICH SHARE A SINGLE BASE HAVE BEEN SORTED OUT
        #In this case there are two criteria that could be used to sort out which of the two geometries is 'LEGITIMATE' and which one is 'FALSE':
                                     if (two_base_shared == 'y') and (MASTER_Basepairs_summary[i][j][2] == MASTER_Basepairs_summary[l][m][2]) and (MASTER_Basepairs_summary[i][j][4] ==  MASTER_Basepairs_summary[l][m][4]) and not ((i == l) and (j == m)):
                                         SCORE[0] = [0, 0]
        #CRITERION #1: CUTOFF_CRITERION: the smaller the CUTOFF value, the more chances that the basepair is correctly assigned.
                                   #The cutoff string is at MASTER_Basepairs_summary[i][j][len(MASTER_Basepairs_summary[i][j]) - 1] and MASTER_Basepairs_summary[l][m][len(MASTER_Basepairs_summary[l][m]) - 1]
#                                         print "MASTER_Basepairs_summary[i][j]", MASTER_Basepairs_summary[i][j], "/nMASTER_Basepairs_summary[l][m]", MASTER_Basepairs_summary[l][m]
                                         line = "\n################\nTwo possible geometries for the basepair formed by " + convert[0] + convert[1] + ":" + convert[2] + convert[3] + ", i.e.: " + MASTER_Basepairs_summary[i][j][0] + " and " + MASTER_Basepairs_summary[l][m][0]
                                         print line
                                         SCORE[1] = CUTOFF_CRITERION(MASTER_Basepairs_summary[i][j], MASTER_Basepairs_summary[l][m], convert[0], convert[1], convert[2], convert[3], convert[4], convert[5], convert[6], convert[7])
        #CRITERION #2: BOND_CRITERION: The larger the number of bonds between two bases, the more chances that the basepair is correctly assigned. Should have the lowest priority
                                         SCORE[2] = BOND_CRITERION(MASTER_Basepairs_summary[i][j], MASTER_Basepairs_summary[l][m], convert[0], convert[1], convert[2], convert[3], convert[4], convert[5], convert[6], convert[7])

#####DECISIONS
#The 'FALSE-BASEPAIR-ASSIGNMENT LOOP' will run at highest stringency until no more FALSE basepairs can be removed, i.e, when 'removed' = 0 and 'while_control' = 'n'.
        #Highest stringency means that in order to assign a pair of debated basepairs as LEGITIMATE and FALSE, the former has to fare better than the latter under all three CRITERIA
        #Therefore, to be assigned as LEGITIMATE and FALSE, one of the debated basepairs must have a score of '3' and the other a score of '0'
                                     if SCORE != [[],[],[]]:
#                                         print "\nENTERING DECISSIONS removed", removed, " while_control", while_control, " two_base_shared", two_base_shared, " single_base_shared", single_base_shared, " stringency", stringency, "score", score
                                         print "CALCULATING SCORE", SCORE,
                                         score[0] = SCORE[0][0] + SCORE[1][0] + SCORE[2][0]
                                         score[1] = SCORE[0][1] + SCORE[1][1] + SCORE[2][1]
                                         print "CALCULATING score", score
                                         decission = []
                                         if (stringency == 'HIGH'):
                                             undetermined = 'n'
                                             if score[0] >= 7:
                                                 decission = DECISSION(MASTER_Basepairs_summary[i][j], MASTER_Basepairs_summary[l][m], convert[0], convert[1], convert[2], convert[3], convert[4], convert[5], convert[6], convert[7], removed, undetermined, two_base_shared)
                                             elif score[1] >= 7:
                                                 decission = DECISSION(MASTER_Basepairs_summary[l][m], MASTER_Basepairs_summary[i][j], convert[4], convert[5], convert[6], convert[7], convert[0], convert[1], convert[2], convert[3], removed, undetermined, two_base_shared)
#                                             print "END DECISSIONS with stringency HIGH: removed", removed, "MASTER_Basepairs_summary[i][j] ", MASTER_Basepairs_summary[i][j], "\nMASTER_Basepairs_summary[l][m] ", MASTER_Basepairs_summary[l][m]
                                         elif (stringency == 'LOW'):
                                             if score[0] > score[1]:
#                                                 print "score[0] > score[1]"
                                                 undetermined = 'n'
                                                 decission = DECISSION(MASTER_Basepairs_summary[i][j], MASTER_Basepairs_summary[l][m], convert[0], convert[1], convert[2], convert[3], convert[4], convert[5], convert[6], convert[7], removed, undetermined, two_base_shared)
                                             elif score[0] < score[1]:
#                                                 print "score[0] < score[1]"
                                                 undetermined = 'n'
                                                 decission = DECISSION(MASTER_Basepairs_summary[l][m], MASTER_Basepairs_summary[i][j], convert[4], convert[5], convert[6], convert[7], convert[0], convert[1], convert[2], convert[3], removed, undetermined, two_base_shared)
                                             else:
                                                 undetermined = 'y'
                                                 decission = DECISSION(MASTER_Basepairs_summary[l][m], MASTER_Basepairs_summary[i][j], convert[0], convert[1], convert[2], convert[3], convert[4], convert[5], convert[6], convert[7], removed, undetermined, two_base_shared)

                                         if decission != []:
                                             MASTER_Basepairs_summary[i][j] = decission[0]
                                             MASTER_Basepairs_summary[l][m] = decission[1]
                                             removed = decission[2]
#                                         print "END DECISSIONS with stringency LOW: removed", removed, "MASTER_Basepairs_summary[i][j] ", MASTER_Basepairs_summary[i][j], "\nMASTER_Basepairs_summary[l][m] ", MASTER_Basepairs_summary[l][m], "\nscore[1]", score[1], "score[0]", score[0]


#Loop control:
    print "Number of succesful deconvolutions =", removed
    if two_base_shared == 'y':
       while_control = 'n' #Ends the 'while' loop
    if stringency == 'LOW' and removed == 0:
       two_base_shared = 'y' #Now the program will search for basepairs with both bases shared and then exit the loop
       MODE = "Two base sharing"
    if stringency == 'HIGH' and removed == 0:
       stringency = 'LOW'
       run_while = 0
#Loop control:

    run_while = run_while + 1


#'FALSE-BASEPAIR-ASSIGNMENT LOOP'             END

#####MAIN part of SECTION DETECTION OF FALSE BASEPAIRS
######################################################
print_control = 0
line1 = " AFTER DETECTION OF FALSE BASEPAIRS"
SEQUENTIAL_READOUT(MASTER_Basepairs_summary, MASTER_Basepairs_schemes, control, print_control, line1)
##########################################################################
##########################################################################
##########################################################################
#SECTION DETECTION OF FALSE BASEPAIRS                  END
##########################################################################
##########################################################################
##########################################################################




##########################################################################
##########################################################################
##########################################################################
#SECTION ORDERING OUTPUT                             BEGINNING
##########################################################################
##########################################################################
##########################################################################

print "\n##########################################################################\nORDERING OUTPUT\n##########################################################################\n"
new_list_end = [] #Will collect all legitimate basepairs from MASTER_Basepairs_summary, from largest first residue to smallest
max = []
yes_count = 1 #This will allow to enter the next 'while'where it will be made '0'
while (yes_count > 0):
    yes_count = 0 #Initializing, will count the # of 'yes' lines in left in MASTER_Basepairs_summary after they are sequentially removed by this loop
    for i in range (len(MASTER_Basepairs_summary)):
       if (len(MASTER_Basepairs_summary[i]) > 0):
          for j in range (len(MASTER_Basepairs_summary[i])):
              if len(MASTER_Basepairs_summary[i][j]) > 0:
                   if 'yes' in MASTER_Basepairs_summary[i][j]:
                      if yes_count == 0: #Initializing 'min' and 'max'
                          count = 1
                          m=re.search('\d+', MASTER_Basepairs_summary[i][j][2])
                          search1 = m.group()
                          compare_to = to_int(search1)
                          max = [compare_to, i, j]
                          yes_count = 1
                      else:
                          m=re.search('\d+', MASTER_Basepairs_summary[i][j][2])
                          search1 = m.group()
                          transient = to_int(search1)
                          if transient >= max[0]:
                              max = [transient, i, j]
                          yes_count = yes_count + 1
    if yes_count > 0:
            new_list_end.append(MASTER_Basepairs_summary[max[1]][max[2]])
            MASTER_Basepairs_summary[max[1]][max[2]] = []

for i in range (len(new_list_end)):
   print new_list_end[i]
##########################################################################
##########################################################################
##########################################################################
#SECTION ORDERING OUTPUT                             END
##########################################################################
##########################################################################
##########################################################################


##########################################################################
##########################################################################
##########################################################################
#SECTION STATISTICS                                  BEGINNING
##########################################################################
##########################################################################
##########################################################################

#Now the list of baseipairs has been ordered and contained in 'new_list_end'. The list will be outputed in several different ways, under two main categories: text output and PYMOL files.
##################################
#####SECTION FUNCTION "ADD_SPACES"
def ADD_SPACES(STRING, NUM):
#STRING = new_list_end[count_for][l]
#NUM = key[l]

    result = NUM - len(STRING)
    space = " "
    STRING = space*result + STRING
    return STRING
#####SECTION FUNCTION "ADD_SPACES"
##################################

###################################
#####SECTION FUNCTION "DIFF_CALC"
def DIFF_CALC(ARR1, ARR2, LIST1, LIST2):
#ARR1 = new_list_end_floats[b,:]
#ARR2 = STATS[c,:]
#LIST1 = can be 'new_list_end[b]',falsely_assigned[b], wrongly_assigned[b]
#LIST2 = MASTER_Basepairs_bonds[c]
#    print "ARR1", ARR1, "ARR2", ARR2, "\LIST1", LIST1
    diff_1 = str(ARR1[0] - ARR2[0])
    diff_2 = str(ARR1[1] - ARR2[2])
    diff_3 = str(ARR1[2] - ARR2[4])
#    diff_P_dist = str(ARR1[3] - ARR2[6])
#    diff_C1_dist = str(ARR1[4] - ARR2[8])
    if LIST2[9] != "NA":
        diff_P_dist = str(ARR1[3] - LIST2[9])
        diff_C1_dist = str(ARR1[4] - LIST2[11])
    else:
        diff_P_dist = str(ARR1[3] - ARR2[6])
        diff_C1_dist = str(ARR1[4] - ARR2[8])
    list = [diff_1[:6], diff_2[:6], diff_3[:6], diff_P_dist[:6], diff_C1_dist[:6]]
    key = 5
    for l in range (len(list)):
        if list[l] == '0.0':
           list[l] = '0.000'
        if list[l][0] != "-":
           list[l] = ' ' + list[l][0:5]
        if len(list[l]) < key:
           list[l] = ADD_SPACES(list[l], key)
#To format the distances so that the ones carrying a '*' do not jutt out
    if len(LIST1[9]) == 5:
        LIST1[9] = LIST1[9] + ' '
    if len(LIST1[12]) == 5:
        LIST1[12] = LIST1[12] + ' '
    if len(LIST1[15]) == 5:
        LIST1[15] = LIST1[15] + ' '
    print LIST1[0], "  ", LIST1[2], " ", LIST1[3], " ", LIST1[4], " ", LIST1[5], "    ", LIST1[1], "    ", LIST1[6], "  ", LIST1[7], " ", LIST1[8], LIST1[9], list[0], LIST1[16], list[3], LIST1[17], list[4]
    print "                                                     ", LIST1[10], " ", LIST1[11], LIST1[12], list[1]
    print "                                                     ", LIST1[13], " ", LIST1[14], LIST1[15], list[2]

    return list
#####SECTION FUNCTION "DIFF_CALC"
##################################


key = [] #Will carry information to format the output
        #OUTPUT FIELDS
                #position 0: "Basepairing Scheme"                                     STRING of 9 characters
key.append(9)
                #position 1: "number of bonds expected"                               STRING of 1 character
key.append(1)
                #position 2: "Base # for base 1"                                      STRING of 4 characters
key.append(4)
                #position 3: "Base identity for base 1"                               STRING of 3 characters
key.append(3)
                #position 4: "Base # for base 2"                                      STRING of 4 characters
key.append(4)
                #position 5: "Base identity for base 2"                               STRING of 3 characters
key.append(3)
                #position 6: "number of bonds found"                                  STRING of 1 character
key.append(1)
                #position 7: "atom 1 for bond 1"                                      STRING of 2 characters
key.append(2)
                #position 8: "atom 2 for bond 1"                                      STRING of 2 characters
key.append(2)
                #position 9: distance for H-bond 1                                    STRING of 4 characters
key.append(4)
                #position 10: "atom 1 for bond 2"                                     STRING of 2 characters
key.append(2)
                #position 11: "atom 2 for bond 2"                                     STRING of 2 characters
key.append(2)
                #position 12: distance for H-bond 2                                   STRING of 4 characters
key.append(4)
                #position 13: "atom 1 for bond 3"                                     STRING of 2 characters
key.append(2)
                #position 14: "atom 2 for bond 3"                                     STRING of 2 characters
key.append(2)
                #position 15: distance for H-bond 2                                   STRING of 4 characters
key.append(4)
                #position 16 "P-P bond distance"                                      STRING of 4 characters
key.append(5)
                #position 17 "Information regarding involvement in MULTIPLE CONTACTS", carried in new_list_end[count_for][len(new_list_end[count_for])-3], or at new_list_end[count_for][len(new_list_end[count_for])-4] if the basepair was assigned by 'MH' or 'ELI'. The information in this string is either '0' (not involved in MULTIPLE CONTACTS) or '1' (involved in MULTIPLE CONTACTS) but will be switched to 'NO' or 'YES', which will make more sense during statistics                                                                       STRING of 3 characters
key.append(3)
                #position 18 "Method used to assign the basepair", carried in new_list_end[count_for][len(new_list_end[count_for])-2]                                              STRING of 3 characters
key.append(3)
                #position 19 "CUTOFF", carried in new_list_end[count_for][len(new_list_end[count_for])-1]                                              STRING of 3 characters
key.append(3)

###########################FORMATTING 'new_list_end'############################################
print "######################################\n##########ORDERED RAW OUTPUT##########        BEGINNING\n######################################"
print "############### '*' notes bonds with length larger tan maximum distance cutoff ###############"
list = []
count_for = len(new_list_end) - 1 #Controld the next 'while'
new_list_end_floats = numpy.zeros(5 * len(new_list_end)).reshape(len(new_list_end), 5) #Will have as many lines as 'new_list_end', so it can be accessed just like 'new_list_end, but will carry only the bond- and P-distances as float, so they can be used in statistical operations
while count_for > -1:
           line1 = ''

####Finding Missing Bonds
#To allow STATISTICS regarding bond distances in the assigned basepairs, the missing bonds for 3-bond basepairs with only 2 bonds identified need to be found.
           missing = []
           if (new_list_end[count_for][1] == 3) and (new_list_end[count_for][6] == 2):
# and ('FALSE' not in new_list_end[count_for]) and ('UNDETERMINED' not in new_list_end[count_for]):
               print "Basepair formed by residues", new_list_end[count_for][2], " and", new_list_end[count_for][4], "is missing one bond. Searching for missing bond"
               missing_control = "0"
#               print "new_list_end[count_for]", new_list_end[count_for]
          #Let's first find the basepair scheme
               for i in range (len(MASTER_Basepairs_schemes)):
                   if MASTER_Basepairs_schemes[i] == new_list_end[count_for][0]:
#                       print "MASTER_Basepairs_schemes[i] ", MASTER_Basepairs_schemes[i], "len(MASTER_Basepairs_schemes)", len(MASTER_Basepairs_schemes), "MASTER_Basepairs_bonds[i] ", MASTER_Basepairs_bonds[i], "len(MASTER_Basepairs_bonds)", len(MASTER_Basepairs_bonds), "len(MASTER_Basepairs_excluded)", len(MASTER_Basepairs_excluded)
#                       for q in range (len(MASTER_Basepairs_excluded[i])):
#                           print MASTER_Basepairs_excluded[i][q]
          #Let's now find the location of the missing bonds in new_list_end[count_for] and in MASTER_Basepairs_bonds[i]. Append the location and the identity of the participating ATOMS to the list 'missing'
                       for j in range (len(new_list_end[count_for])):
                           if new_list_end[count_for][j] == []:
                               if j == 7:
                                  missing.append(j)
                                  missing.append(MASTER_Basepairs_bonds[i][3])
                                  missing.append(j+1)
                                  missing.append(MASTER_Basepairs_bonds[i][4])
                                  break
                               elif j == 10:
                                  missing.append(j)
                                  missing.append(MASTER_Basepairs_bonds[i][5])
                                  missing.append(j+1)
                                  missing.append(MASTER_Basepairs_bonds[i][6])
                                  break
                               elif j == 13:
                                  missing.append(j)
                                  missing.append(MASTER_Basepairs_bonds[i][7])
                                  missing.append(j+1)
                                  missing.append(MASTER_Basepairs_bonds[i][8])
                                  break
          #Let's now find the line in MASTER_Basepairs that carries the missing bond
#                       print "missing", missing
                       for k in range (len(MASTER_Basepairs_excluded[i])):
#                           print "MASTER_Basepairs_excluded[i][k]", MASTER_Basepairs_excluded[i][k]
                           if (new_list_end[count_for][2] == MASTER_Basepairs_excluded[i][k][0]) and (new_list_end[count_for][4] == MASTER_Basepairs_excluded[i][k][3]):
                               if (missing[1] in MASTER_Basepairs_excluded[i][k][2]) and (missing[3] in MASTER_Basepairs_excluded[i][k][5]):
                                  new_list_end[count_for][missing[0]] = MASTER_Basepairs_excluded[i][k][2]
                                  new_list_end[count_for][missing[2]] = MASTER_Basepairs_excluded[i][k][5]
                                  new_list_end[count_for][missing[2]+1] = MASTER_Basepairs_excluded[i][k][6]
                                  missing_control = "1"
                       if missing_control == "0": #This will prevent the program from crashing if a bond remains missing. Following these basepairs can lead to the identification of new basepair geometries
                             new_list_end[count_for][missing[0]] = "1"
                             new_list_end[count_for][missing[2]] = "1"
                             new_list_end[count_for][missing[2]+1] = "11"
                             #The newly found bond will be marked for output a few lines below
####Finding Missing Bonds



#DELETING UNWANTED SPACES and SAVING THE DISTANCE AS FLOATS IN 'new_list_end_floats' FOR FUTURE MATHEMATICAL OPERATIONS
        #For basepairs with only 2 bonds identified, positions 13, 14, and 15 will have []. Let's replace them with a string carrying '0', or '00' for new_list_end[count_for][15], so 'FUNCTION CONVERT' will work properly (see also next explanation)
#           print "\nnew_list_end[count_for]", new_list_end[count_for]
#           print "count_for", count_for, "new_list_end_floats[count_for]", new_list_end_floats[count_for]
           if "11" != new_list_end[count_for][15] and "11" != new_list_end[count_for][12] and "11" != new_list_end[count_for][9]:
               if new_list_end[count_for][13] == []: #Basepairs with only 2 bonds
                   print "new_list_end[count_for]", new_list_end[count_for]
                   new_list_end[count_for][13] = '0'
                   new_list_end[count_for][14] = '0'
                   new_list_end[count_for][15] = '00'
                   new_list_end_floats[count_for,:] = new_list_end_floats[count_for,:] + [new_list_end[count_for][9], new_list_end[count_for][12], 0, new_list_end[count_for][16], new_list_end[count_for][17]]
               else:
                   print "new_list_end[count_for]", new_list_end[count_for]
                   new_list_end_floats[count_for,:] = new_list_end_floats[count_for,:] + [new_list_end[count_for][9], new_list_end[count_for][12], new_list_end[count_for][15], new_list_end[count_for][16], new_list_end[count_for][17]]
                   a3 = str(new_list_end[count_for][15])
                   new_list_end[count_for][15] = a3[:5]
               a1 = str(new_list_end[count_for][9])
               new_list_end[count_for][9] = a1[:5]
               a2 = str(new_list_end[count_for][12])
               new_list_end[count_for][12] = a2[:5]
               a3 = str(new_list_end[count_for][16])
               if new_list_end[count_for][16] < 10:
                   new_list_end[count_for][16] = ' ' + a3[:5]
               if new_list_end[count_for][16] >= 10:
                   new_list_end[count_for][16] = a3[:6]
               a4 = str(new_list_end[count_for][17])
               new_list_end[count_for][17] = a4[:6]
        #positions 2, 4, 7 and 8, 10 and 11, 13 and 14,  have lefthand spaces. They need to be removed. Strings "new_list_end[count_for][7]" and above are introduced twice to make 'FUNCTION CONVERT' work properly. First time through "FUNCTION CONVERT" the 'number' part of the string is extracted and the second time the 'word' part of the string is extracted

               list = [new_list_end[count_for][2], new_list_end[count_for][3], new_list_end[count_for][4], new_list_end[count_for][5], new_list_end[count_for][7], new_list_end[count_for][7], new_list_end[count_for][8], new_list_end[count_for][8], new_list_end[count_for][10], new_list_end[count_for][10], new_list_end[count_for][11], new_list_end[count_for][11], new_list_end[count_for][13], new_list_end[count_for][13], new_list_end[count_for][14], new_list_end[count_for][14]]

               convert = CONVERT(list)
               new_list_end[count_for][2] = convert[0]
               new_list_end[count_for][4] = convert[2]
               new_list_end[count_for][7] = convert[5] + convert[4]
               new_list_end[count_for][8] = convert[7] + convert[6]
               new_list_end[count_for][10] = convert[9] + convert[8]
               new_list_end[count_for][11] = convert[11] + convert[10]
               new_list_end[count_for][13] = convert[13] + convert[12]
               new_list_end[count_for][14] = convert[15] + convert[14]
               if new_list_end[count_for][13] == '00':
                   new_list_end[count_for][13] = '--'
                   new_list_end[count_for][14] = '--'
                   new_list_end[count_for][15] = '  -- '
               if new_list_end[count_for][7] == '11' or new_list_end[count_for][10] == '11' or new_list_end[count_for][13] == '11':
                   new_list_end[count_for][13] = 'NF'
                   new_list_end[count_for][14] = 'NF'
                   new_list_end[count_for][15] = '  NF '

#DELETING UNWANTED SPACES

#ADDING DESIRED SPACES
#Adding spaces for format. In conjuction with 'SECTION FUNCTION ADD_SPACES'
           #Formatting most of the info carried in 'new_list_end'
           for l in range (0,17): #For most of the info carried in 'new_list_end'
              a = str(new_list_end[count_for][l])
              if len(a) < key[l]:
                 new_list_end[count_for][l] = ADD_SPACES(a, key[l])
              else:
                 new_list_end[count_for][l] = a
              if missing != [] and l == missing[2]+1 and missing_control == "1":
#                 print "[missing[2]+1] ", [missing[2]+1], "and l", l
                 new_list_end[count_for][missing[2]+1] = new_list_end[count_for][missing[2]+1] + '*'
#                 print "AFTER append new_list_end[count_for][missing[2]+1]", new_list_end[count_for][missing[2]+1]
              line1 = line1 + new_list_end[count_for][l] + ' '
           #Formatting the method ('MH', 'D', 'ELI'), carried in new_list_end[count_for][len(new_list_end[count_for])-3], after appending the 'MULTIPLE CONTACTS' information, right above
           if ('FALSE' in new_list_end[count_for]) or ('WRONG' in new_list_end[count_for]):
               new_list_end[count_for][len(new_list_end[count_for])-3] = ADD_SPACES(new_list_end[count_for][len(new_list_end[count_for])-3], key[17])
               line1 = line1 + new_list_end[count_for][len(new_list_end[count_for])-3] + ' '
           else:
               new_list_end[count_for][len(new_list_end[count_for])-2] = ADD_SPACES(new_list_end[count_for][len(new_list_end[count_for])-2], key[17])
               line1 = line1 + new_list_end[count_for][len(new_list_end[count_for])-2] + ' '
           #Formatting the cutoff, carried in new_list_end[count_for][len(new_list_end[count_for])-2], after appending the 'MULTIPLE CONTACTS' information, right above
           if ('FALSE' in new_list_end[count_for]) or ('WRONG' in new_list_end[count_for]):
                line1 = line1 + new_list_end[count_for][len(new_list_end[count_for])-2] + ' '
           else:
                line1 = line1 + new_list_end[count_for][len(new_list_end[count_for])-1] + ' '
           #Formatting the MULTIPLE CONTACTS information ('0' or '1'), carried in new_list_end[count_for][len(new_list_end[count_for])-3], or at new_list_end[count_for][len(new_list_end[count_for])-4] if the basepair was assigned by 'MH' or 'ELI'. To make it simpler later, the '0' or '1' will be switched to 'NO' or 'YES' and the formatted string will be appended at the end of 'new_list_end[count_for]'
           if ('FALSE' in new_list_end[count_for]) or ('WRONG' in new_list_end[count_for]) or ('UNDETERMINED' in new_list_end[count_for]):
               place = 4
               if ('FALSE' in new_list_end[count_for]):
                  ending = ' FALSE'
               elif ('WRONG' in new_list_end[count_for]):
                  ending = 'WRONG'
               elif ('UNDETERMINED' in new_list_end[count_for]):
                  ending = ' UNDETERMINED'
           else:
               place = 3
               ending = ''
           if '  D' in new_list_end[count_for]:
               if '0' in str(new_list_end[count_for][len(new_list_end[count_for])- place]):
                  a = 'NO'
               if '1' in str(new_list_end[count_for][len(new_list_end[count_for])- place]):
                  a = 'YES'
               new_list_end[count_for].append(ADD_SPACES(a, key[16]))
           else:
               if '0' in str(new_list_end[count_for][len(new_list_end[count_for])- place -1]):
                  a = 'NO'
               if '1' in str(new_list_end[count_for][len(new_list_end[count_for])- place -1]):
                  a = 'YES'
               new_list_end[count_for].append(ADD_SPACES(a, key[16]))
           line1 = line1 + new_list_end[count_for][len(new_list_end[count_for])-1] + ending
           print line1
#           print new_list_end[count_for]
           count_for = count_for - 1
print "############### '*' notes bonds with length larger tan maximum distance cutoff ###############"
print "######################################\n##########ORDERED RAW OUTPUT##########       END\n######################################"
#ADDING DESIRED SPACES

#print "AFTER loading new_list_end_floats \nlen(new_list_end)", len(new_list_end), "new_list_end_floats.ndim", new_list_end_floats.ndim, "new_list_end_floats.shape", new_list_end_floats.shape, "count_for", count_for
###########################FORMATTING 'new_list_end'############################################

#################Transferring the list of ORDERED and FORMATTED basepairs to run_cutoff_LISTS, in which they will be further order by CUTOFF. Also 'FALSE', 'LEGITIMATE', and 'UNDETERMINED' basepairs get transferred to corresponding lists
LEGITIMATE_tagged = [] #Will collect basepairs with a 'LEGITIMATE' tag from MASTER_Basepairs_summary, from largest first residue to smallest
falsely_assigned = [] #Will collect FALSE basepairs from MASTER_Basepairs_summary, from largest first residue to smallest
wrongly_assigned = [] #Will collect basepairs with 'WRONG' geometry from MASTER_Basepairs_summary, from largest first residue to smallest
undetermined = [] #Will collect 'UNDETERMINED' basepairs from MASTER_Basepairs_summary, from largest first residue to smallest
missing_bond = [] #Will collect basepairs with missing bonds. Tagged as " miss"
file = [[], []]
pymol0 = directory + "/" + pdb_file_main + "_all-basepairs_REVERSED" + "_PYMOL" + "_script.pml"
file[0] = open(pymol0, 'w')
line0 = "load " + pdb_file + '\n' + "color white, /" + pdb_file_main + "/*" + '\n' + "set stick_ball, on" + '\n' + "set stick_ball_ratio, 2.00000" + '\n' + "show sticks" + '\n'
file[0].write(line0)
pymol1 = directory + "/" + pdb_file_main + "_all-basepairs" + "_PYMOL" + "_script.pml"
file[1] = open(pymol1, 'w')
line1 = "load " + pdb_file + '\n' + "color white, /" + pdb_file_main + "/*" + '\n' + "set stick_ball, on" + '\n' + "set stick_ball_ratio, 2.00000" + '\n'
file[1].write(line1)
From = 'basepairs'

basepair_count = 0
for h in range (len(run_cutoff)/2):
   run_cutoff_LISTS.append([])
for b in range (len(new_list_end)):
   if 'REMOVED' not in new_list_end[b] and 'FALSE' not in new_list_end[b] and 'WRONG' not in new_list_end[b] and 'UNDETERMINED' not in new_list_end[b] and " miss" not in new_list_end[b]:
       PYMOL_OUTPUT(new_list_end[b], b, pdb_file_main, From, file)
   if ('LEGITIMATE' in new_list_end[b]):
       LEGITIMATE_tagged.append(new_list_end[b])
#      LEGITIMATE_tagged[len(LEGITIMATE_tagged)-1].append(b) #For easy future access to the corresponding line in 'new_list_end'
   elif ('FALSE' in new_list_end[b]):
      falsely_assigned.append(new_list_end[b])
      falsely_assigned[len(falsely_assigned)-1].append(b) #For easy future access to the corresponding line in 'new_list_end'
   elif ('WRONG' in new_list_end[b]):
      wrongly_assigned.append(new_list_end[b])
      wrongly_assigned[len(wrongly_assigned)-1].append(b) #For easy future access to the corresponding line in 'new_list_end'
   elif ('UNDETERMINED' in new_list_end[b]):
      undetermined.append(new_list_end[b])
      undetermined[len(undetermined)-1].append(b) #For easy future access to the corresponding line in 'new_list_end'
   elif(" miss" in new_list_end[b]):
      missing_bond.append(new_list_end[b])
      missing_bond[len(missing_bond)-1].append(b) #For easy future access to the corresponding line in 'new_list_end'
   for h in range (len(run_cutoff_LISTS)):
      if run_cutoff[h*2 + 1] in new_list_end[b] and ('REMOVED' not in new_list_end[b]) and ('FALSE' not in new_list_end[b]) and ('WRONG' not in new_list_end[b]) and ('UNDETERMINED' not in new_list_end[b]) and " miss" not in new_list_end[b]:
           run_cutoff_LISTS[h].append(new_list_end[b])
           basepair_count = basepair_count + 1

CLOSE(file)

#################Transferring the list of ORDERED and FORMATTED basepairs to run_cutoff_LISTS, in which they will be further order by CUTOFF. Also 'FALSE', 'LEGITIMATE', and 'UNDETERMINED' basepairs get transferred to corresponding lists


#########################
########STATISTICS OUTPUT              BEGINNING
print "\n\n#######################################################\n##################### STATISTICS ######################\n#######################################################\n"
print "###### A total of", basepair_count, "basepairs has been identified"
##### 1) PRINTING OUTPUT BY CUTOFF    BEGINNING
file = [[], []]
for h in range (len(run_cutoff)/2):
   print "\n\nBASEPAIRS FULLY IDENTIFIED AT DISTANCE CUTOFF = ", run_cutoff[h*2 + 1]
   print           "==============================================================================="
   print "   Scheme   Resid  Base  Resid  Base Expect. Found ATOM ATOM  DIST. Method M.I."
   print           "             #1     #1    #2     #2   Bonds  Bonds  #1   #2"
   print           "==============================================================================="

   #For PYMOL output
   pymol1 = directory + "/" + pdb_file_main + "_basepairs-at-cutoff_" + run_cutoff[h*2 + 1] + "_PYMOL" + "_script.pml"
   file[1] = open(pymol1, 'w')
   line1 = "load " + pdb_file + '\n' + "color white, /" + pdb_file_main + "/*" + '\n' + "set stick_ball, on" + '\n' + "set stick_ball_ratio, 2.00000" + '\n'
   file[1].write(line1)
   #For PYMOL output

   for a in range (len(run_cutoff_LISTS[h])):

      PYMOL_OUTPUT(run_cutoff_LISTS[h][a], a, pdb_file_main, From, file)

      print run_cutoff_LISTS[h][a][0], "  ", run_cutoff_LISTS[h][a][2], " ", run_cutoff_LISTS[h][a][3], " ", run_cutoff_LISTS[h][a][4], " ", run_cutoff_LISTS[h][a][5], "    ", run_cutoff_LISTS[h][a][1], "    ", run_cutoff_LISTS[h][a][6], "  ", run_cutoff_LISTS[h][a][7], " ", run_cutoff_LISTS[h][a][8], run_cutoff_LISTS[h][a][9], "  ", run_cutoff_LISTS[h][a][len(run_cutoff_LISTS[h][a])-3], run_cutoff_LISTS[h][a][len(run_cutoff_LISTS[h][a])-1]
      print "                                                     ", run_cutoff_LISTS[h][a][10], " ", run_cutoff_LISTS[h][a][11], run_cutoff_LISTS[h][a][12]
      print "                                                     ", run_cutoff_LISTS[h][a][13], " ", run_cutoff_LISTS[h][a][14], run_cutoff_LISTS[h][a][15]
   print           "==============================================================================="
   print "Number of basepairs identified at CUTOFF =", run_cutoff[h*2 + 1], "is", len(run_cutoff_LISTS[h])
print "############### notes bonds with length larger tan maximum distance cutoff ###############"
CLOSE(file)
##### 1) PRINTING OUTPUT BY CUTOFF    END

##### 2) LIST OF LEGITIMATE BASEPAIRS
if LEGITIMATE_tagged != []:
    file = [[], []]
    pymol1 = directory + "/" + pdb_file_main + "_LEGITIMATED-basepairs" + "_PYMOL" + "_script.pml"
    file[1] = open(pymol1, 'w')
    line1 = "load " + pdb_file + '\n' + "color white, /" + pdb_file_main + "/*" + '\n' + "set stick_ball, on" + '\n' + "set stick_ball_ratio, 2.00000" + '\n'
    file[1].write(line1)
    color = 'RED'
    From = ''

    print "\n\nLIST OF LEGITIMATE BASEPAIRS FOUND TO SHARE BASES WITH FALSELY ASSIGNED BASEPAIRS"
    print           "==============================================================================="
    print "   Scheme   Resid  Base  Resid  Base Expect. Found ATOM ATOM  DIST. Method M.I."
    print           "             #1     #1    #2     #2   Bonds  Bonds  #1   #2"
    print           "==============================================================================="
    for h in range (len(LEGITIMATE_tagged)):
       PYMOL_OUTPUT(LEGITIMATE_tagged[h], h, pdb_file_main, From, file)
       print LEGITIMATE_tagged[h][0], "  ", LEGITIMATE_tagged[h][2], " ", LEGITIMATE_tagged[h][3], " ", LEGITIMATE_tagged[h][4], " ", LEGITIMATE_tagged[h][5], "    ", LEGITIMATE_tagged[h][1], "    ", LEGITIMATE_tagged[h][6], "  ", LEGITIMATE_tagged[h][7], " ", LEGITIMATE_tagged[h][8], LEGITIMATE_tagged[h][9], "  ", LEGITIMATE_tagged[h][len(LEGITIMATE_tagged[h])-4], LEGITIMATE_tagged[h][len(LEGITIMATE_tagged[h])-2]
       print "                                                     ", LEGITIMATE_tagged[h][10], " ", LEGITIMATE_tagged[h][11], LEGITIMATE_tagged[h][12]
       print "                                                     ", LEGITIMATE_tagged[h][13], " ", LEGITIMATE_tagged[h][14], LEGITIMATE_tagged[h][15]
    print           "==============================================================================="
    print "############### '*' notes bonds with length larger tan maximum distance cutoff ###############"
    print "Number of LEGITIMATE BASEPAIRS FOUND TO SHARE BASES WITH FALSELY ASSIGNED BASEPAIRS =", len(LEGITIMATE_tagged)
    CLOSE(file)

##### 3) LIST OF FALSE BASEPAIRS
if falsely_assigned != []:
    file = [[], []]
    pymol1 = directory + "/" + pdb_file_main + "_FALSE-basepairs" + "_PYMOL" + "_script.pml"
    file[1] = open(pymol1, 'w')
    line1 = "load " + pdb_file + '\n' + "color white, /" + pdb_file_main + "/*" + '\n' + "set stick_ball, on" + '\n' + "set stick_ball_ratio, 2.00000" + '\n'
    file[1].write(line1)
    color = 'RED'
    From = ''

    print "\n\nLIST OF FALSE BASEPAIRS"
    print           "==============================================================================="
    print "   Scheme   Resid  Base  Resid  Base Expect. Found ATOM ATOM  DIST. Method M.I."
    print           "             #1     #1    #2     #2   Bonds  Bonds  #1   #2"
    print           "==============================================================================="
    for h in range (len(falsely_assigned)):
       PYMOL_OUTPUT(falsely_assigned[h], h, pdb_file_main, From, file)
       print falsely_assigned[h][0], "  ", falsely_assigned[h][2], " ", falsely_assigned[h][3], " ", falsely_assigned[h][4], " ", falsely_assigned[h][5], "    ", falsely_assigned[h][1], "    ", falsely_assigned[h][6], "  ", falsely_assigned[h][7], " ", falsely_assigned[h][8], falsely_assigned[h][9], "  ", falsely_assigned[h][len(falsely_assigned[h])-5], falsely_assigned[h][len(falsely_assigned[h])-2]
       print "                                                     ", falsely_assigned[h][10], " ", falsely_assigned[h][11], falsely_assigned[h][12]
       print "                                                     ", falsely_assigned[h][13], " ", falsely_assigned[h][14], falsely_assigned[h][15]
    print           "==============================================================================="
    print "############### '*' notes bonds with length larger tan maximum distance cutoff ###############"
    print "Number of FALSELY ASSIGNED BASEPAIRS =", len(falsely_assigned)
    CLOSE(file)

##### 4) LIST OF BASEPAIRS WITH WRONG GEOMETRY
if wrongly_assigned != []:
    file = [[], []]
    pymol1 = directory + "/" + pdb_file_main + "_WRONG-basepairs" + "_PYMOL" + "_script.pml"
    file[1] = open(pymol1, 'w')
    line1 = "load " + pdb_file + '\n' + "color white, /" + pdb_file_main + "/*" + '\n' + "set stick_ball, on" + '\n' + "set stick_ball_ratio, 2.00000" + '\n'
    file[1].write(line1)
    color = 'RED'
    From = ''

    print "\n\nLIST OF BASEPAIRS WITH WRONG GEOMETRY"
    print           "==============================================================================="
    print "   Scheme   Resid  Base  Resid  Base Expect. Found ATOM ATOM  DIST. Method M.I."
    print           "             #1     #1    #2     #2   Bonds  Bonds  #1   #2"
    print           "==============================================================================="
    for h in range (len(wrongly_assigned)):
       PYMOL_OUTPUT(wrongly_assigned[h], h, pdb_file_main, From, file)
       print wrongly_assigned[h][0], "  ", wrongly_assigned[h][2], " ", wrongly_assigned[h][3], " ", wrongly_assigned[h][4], " ", wrongly_assigned[h][5], "    ", wrongly_assigned[h][1], "    ", wrongly_assigned[h][6], "  ", wrongly_assigned[h][7], " ", wrongly_assigned[h][8], wrongly_assigned[h][9], "  ", wrongly_assigned[h][len(wrongly_assigned[h])-5], wrongly_assigned[h][len(wrongly_assigned[h])-2]
       print "                                                     ", wrongly_assigned[h][10], " ", wrongly_assigned[h][11], wrongly_assigned[h][12]
       print "                                                     ", wrongly_assigned[h][13], " ", wrongly_assigned[h][14], wrongly_assigned[h][15]
    print           "==============================================================================="
    print "############### '*' notes bonds with length larger tan maximum distance cutoff ###############"
    print "Number of BASEPAIRS WITH WRONG GEOMETRY =", len(wrongly_assigned)
    CLOSE(file)


##### 5) LIST OF UNDETERMINED BASEPAIRS BEGINNING
if undetermined != []:
    file = [[], []]
    pymol1 = directory + "/" + pdb_file_main + "_UNDETERMINED-basepairs" + "_PYMOL" + "_script.pml"
    file[1] = open(pymol1, 'w')
    line1 = "load " + pdb_file + '\n' + "color white, /" + pdb_file_main + "/*" + '\n' + "set stick_ball, on" + '\n' + "set stick_ball_ratio, 2.00000" + '\n'
    file[1].write(line1)
    color = 'RED'
    From = ''
    print "\n\nLIST OF BASE-SHARING BASEPAIRS THAT CANNOT BE DECONVOLUTED"
    print           "==============================================================================="
    print "   Scheme   Resid  Base  Resid  Base Expect. Found ATOM ATOM  DIST. Method M.I."
    print           "             #1     #1    #2     #2   Bonds  Bonds  #1   #2"
    print           "==============================================================================="
    for h in range (len(undetermined)):
      PYMOL_OUTPUT(undetermined[h], h, pdb_file_main, From, file)
      print undetermined[h][0], "  ", undetermined[h][2], " ", undetermined[h][3], " ", undetermined[h][4], " ", undetermined[h][5], "    ", undetermined[h][1], "    ", undetermined[h][6], "  ", undetermined[h][7], " ", undetermined[h][8], undetermined[h][9], "  ", undetermined[h][len(undetermined[h])-4], undetermined[h][len(undetermined[h])-1]
      print "                                                     ", undetermined[h][10], " ", undetermined[h][11], undetermined[h][12]
      print "                                                     ", undetermined[h][13], " ", undetermined[h][14], undetermined[h][15]
    print           "==============================================================================="
    print "Number of BASE-SHARING BASEPAIRS THAT CANNOT BE DECONVOLUTED =", len(undetermined)
    print "############### '*' notes bonds with length larger tan maximum distance cutoff ###############"
    CLOSE(file)
##### 5) LIST OF UNDETERMINED BASEPAIRS END

##### 6) LIST OF MISSING-BOND BASEPAIRS BEGINNING
if missing_bond != []:
    file = [[], []]
    pymol1 = directory + "/" + pdb_file_main + "_MISSING-BOND-basepairs" + "_PYMOL" + "_script.pml"
    file[1] = open(pymol1, 'w')
    line1 = "load " + pdb_file + '\n' + "color white, /" + pdb_file_main + "/*" + '\n' + "set stick_ball, on" + '\n' + "set stick_ball_ratio, 2.00000" + '\n'
    file[1].write(line1)
    color = 'RED'
    From = ''
    print "\n\nLIST OF MISSING-BOND BASEPAIRS"
    print           "==============================================================================="
    print "   Scheme   Resid  Base  Resid  Base Expect. Found ATOM ATOM  DIST. Method M.I."
    print           "             #1     #1    #2     #2   Bonds  Bonds  #1   #2"
    print           "==============================================================================="
    for h in range (len(missing_bond)):
      PYMOL_OUTPUT(missing_bond[h], h, pdb_file_main, From, file)
      print missing_bond[h][0], "  ", missing_bond[h][2], " ", missing_bond[h][3], " ", missing_bond[h][4], " ", missing_bond[h][5], "    ", missing_bond[h][1], "    ", missing_bond[h][6], "  ", missing_bond[h][7], " ", missing_bond[h][8], missing_bond[h][9], "  ", missing_bond[h][len(missing_bond[h])-4], missing_bond[h][len(missing_bond[h])-1]
      print "                                                     ", missing_bond[h][10], " ", missing_bond[h][11], missing_bond[h][12]
      print "                                                     ", missing_bond[h][13], " ", missing_bond[h][14], missing_bond[h][15]
    print           "==============================================================================="
    print "Number of MISSING-BOND BASEPAIRS =", len(missing_bond)
    print "############### 'NF' notes Not-found bond ###############"
    CLOSE(file)
##### 6) LIST OF MISSING-BOND BASEPAIRS END



##### 7) LIST OF BASEPAIRS BY GEOMETRY BEGINNIG
###BOND STATISTICS. Important arrays:
        # new_list_end_floats = numpy.zeros(4 * len(new_list_end)).reshape(len(new_list_end), 4) #Will have as many lines as 'new_list_end', so it can be accessed just like 'new_list_end, but will carry only the bond- and P-distances as float, so they can be used in statistical operations
                #Defined in Subsection "DELETING UNWANTED SPACES and SAVING THE DISTANCE AS FLOATS" above
        # STATS, see definition below


list = []
STATS = numpy.zeros(11 * len(MASTER_Basepairs_schemes)).reshape(len(MASTER_Basepairs_schemes), 11) #can be accessed with MASTER_Basepairs_schemes indexes. Will store:
                #STATS[a,0] = AVERAGE for bond-1's length for the GEOMETRY specified by the 'MASTER_Basepairs_schemes[a]'
                #STATS[a,1] = STANDARD DEVIATION for bond-1's length for the GEOMETRY specified by the 'MASTER_Basepairs_schemes[a]'
                #STATS[a,2] = AVERAGE for bond-2's length for the GEOMETRY specified by the 'MASTER_Basepairs_schemes[a]'
                #STATS[a,3] = STANDARD DEVIATION for bond-2's length for the GEOMETRY specified by the 'MASTER_Basepairs_schemes[a]'
                #STATS[a,4] = AVERAGE for bond-3's length for the GEOMETRY specified by the 'MASTER_Basepairs_schemes[a]'
                #STATS[a,5] = STANDARD DEVIATION for bond-3's length for the GEOMETRY specified by the 'MASTER_Basepairs_schemes[a]'
                #STATS[a,6] = AVERAGE for P-P distance for the GEOMETRY specified by the 'MASTER_Basepairs_schemes[a]'
                #STATS[a,7] = STANDARD DEVIATION for P-P distance for the GEOMETRY specified by the 'MASTER_Basepairs_schemes[a]'
                #STATS[a,8] = AVERAGE for C1'-C1' distance for the GEOMETRY specified by the 'MASTER_Basepairs_schemes[a]'
                #STATS[a,9] = STANDARD DEVIATION for C1'-C1' distance for the GEOMETRY specified by the 'MASTER_Basepairs_schemes[a]'
                #STATS[a,] = N, nuber of basepairs found of the GEOMETRY specified by the 'MASTER_Basepairs_schemes[a]'
for a in range (len(MASTER_Basepairs_schemes)):
    SUM = []
    counter = 0
#Calculating mean and standard deviation
#    if "XXVIII_GU" in MASTER_Basepairs_schemes[a]:
#        print "MASTER_Basepairs_schemes[a]", MASTER_Basepairs_schemes[a]
    for e in range (len(new_list_end)):
#        print "OUTSIDE \"for e\" new_list_end[e]", new_list_end[e]
        if  MASTER_Basepairs_schemes[a] in new_list_end[e][0] and 'REMOVED' not in new_list_end[e] and 'FALSE' not in new_list_end[e] and 'WRONG' not in new_list_end[e] and 'UNDETERMINED' not in new_list_end[e] and 'miss' not in new_list_end[e]:
#            print "IN \"for e\" new_list_end[e]", new_list_end[e]
            SUM.append(new_list_end_floats[e,:])
            counter = counter + 1
    if SUM != []:
        SUM = numpy.array(SUM)
#        if "XXVIII_GU" in MASTER_Basepairs_schemes[a]:
#            print "SUM", SUM, "\nSUM.ndim", SUM.ndim, "SUM.shape", SUM.shape, "counter", counter, "SUM.shape[0]",SUM.shape[0]
        SUM.reshape(counter, 5)
#        print "STATS[a,]", STATS[a,], "\nSUM[0,]", SUM[0,]
#        print "numpy.mean(SUM[:,0])", numpy.mean(SUM[:,0])
        STATS[a,] = [numpy.mean(SUM[:,0]), numpy.std(SUM[:,0]), numpy.mean(SUM[:,1]), numpy.std(SUM[:,1]), numpy.mean(SUM[:,2]), numpy.std(SUM[:,2]), numpy.mean(SUM[:,3]), numpy.std(SUM[:,3]), numpy.mean(SUM[:,4]), numpy.std(SUM[:,4]), SUM.shape[0]]
#        if "XXVIII_GU" in MASTER_Basepairs_schemes[a]:
#            print "STATS[a,]", STATS[a,], "STATS[a,0]", STATS[a,0]
#        list = [str(STATS[a,0]), str(STATS[a,1]), str(STATS[a,2]), str(STATS[a,3]), str(STATS[a,4]), str(STATS[a,5]), str(STATS[a,6]), str(STATS[a,7]), str(STATS[a,8])]
#        print "Bond1: MEAN =", list[0][:5], "SD =", list[1][:5], "Bond2: MEAN = ", list[2][:5], "SD =", list[3][:5], "Bond3: MEAN =", list[4][:5], "SD =", list[5][:5], "P-to-P distance =", list[6][:5], "SD =", list[7][:5], "Number of basepairs with this geometry =", list[8][:1]

#SUMMARY STATISTICS BY BASEPAIR
for a in range (len(MASTER_Basepairs_schemes)):
    print MASTER_Basepairs_schemes[a],
    list = [str(STATS[a,0]), str(STATS[a,1]), str(STATS[a,2]), str(STATS[a,3]), str(STATS[a,4]), str(STATS[a,5]), str(STATS[a,6]), str(STATS[a,7]), str(STATS[a,8]), str(STATS[a,9]), str(STATS[a,10])]
    if list[0] != '0.0':
        print "STATISTICAL VALUES CALCULATED WITHIN STRUCTURE"
        print "Bond1: MEAN =", list[0][:5], "SD =", list[1][:5], "Bond2: MEAN = ", list[2][:5], "SD =", list[3][:5], "Bond3: MEAN =", list[4][:5], "SD =", list[5][:5], "P-to-P distance =", list[6][:5], "SD =", list[7][:5], "C1\'-to-C1\' distance =", list[8][:5], "SD =", list[9][:5], "No Basepairs =", list[7][:1]
    else:
        print "No basepairs of this geometry were found"

#WHOLE LIST
list = []
print "\n\nLIST OF BASEPAIRS ORDERED BY GEOMETRY"
for c in range (len(MASTER_Basepairs_schemes)):
    geometry_counter = 0 #Will count all basepairs within 'new_list_end[b]' with 'MASTER_Basepairs_schemes[c]' geometry
    geometry_counter_real = 0 #Will count only non-FALSE, non-UNDETERMINED basepairs within 'new_list_end[b]' with 'MASTER_Basepairs_schemes[c]' geometry
#TRUE BASEPAIRS
    for b in range (len(new_list_end)):
        m=re.search('\S+', new_list_end[b][0])
        scheme = str(m.group())
        if MASTER_Basepairs_schemes[c] == scheme:
           geometry_counter = geometry_counter + 1
           if geometry_counter == 1:
              print           "======================================================================================================="
              print "   Scheme   Resid  Base  Resid  Base Expect. Found  ATOM ATOM  DIST.  DIFF.  DIST   DIFF.  DIST   DIFF."
              print           "             #1     #1    #2     #2   Bonds  Bonds   #1   #2    BOND         P-P          C1'-C1'"
              print           "======================================================================================================="
           if 'REMOVED' not in new_list_end[b] and 'FALSE' not in new_list_end[b] and 'WRONG' not in new_list_end[b] and 'UNDETERMINED' not in new_list_end[b]:
               geometry_counter_real = geometry_counter_real + 1
               list = DIFF_CALC(new_list_end_floats[b,:], STATS[c,:], new_list_end[b], MASTER_Basepairs_bonds[c])

    if geometry_counter_real > 0 and 'REMOVED' not in new_list_end[b] and 'FALSE' not in new_list_end[b] and 'WRONG' not in new_list_end[b] and 'UNDETERMINED' not in new_list_end[b]:
         list = [str(STATS[c,0]), str(STATS[c,1]), str(STATS[c,2]), str(STATS[c,3]), str(STATS[c,4]), str(STATS[c,5]), str(STATS[c,6]), str(STATS[c,7]), str(STATS[c,8]), str(STATS[c,9])]
         print "PRECALCULATED STATISTICAL DATA:", "\nP-to-P distance =", MASTER_Basepairs_bonds[c][9], " SD =", MASTER_Basepairs_bonds[c][10], "\nC1\'-to-C1\' distance =", MASTER_Basepairs_bonds[c][11], " SD =", MASTER_Basepairs_bonds[c][12]
         print "STATISTICAL DATA FROM CURRENT STRUCTURE:", "\nBond1: MEAN =", list[0][:5], "SD =", list[1][:5], "\nBond2: MEAN =", list[2][:5], " SD =", list[3][:5], "\nBond3: MEAN =", list[4][:5], " SD =", list[5][:5], "\nP-to-P distance =", list[6][:5], " SD =", list[7][:5], "\nC1\'-to-C1\' distance =", list[8][:5], " SD =", list[9][:5]
         print "NUMBER OF BASEPAIRS WITH", MASTER_Basepairs_schemes[c], "geometry =", geometry_counter_real, "\n"

#FALSE BASEPAIRS
    false_counter = 0
    for b in range (len(falsely_assigned)):
        if MASTER_Basepairs_schemes[c] in falsely_assigned[b][0]:
            if b == 0:
               print "LIST FALSE BASEPAIRS WITH", MASTER_Basepairs_schemes[c], "geometry:"
            false_counter = false_counter + 1
            list = DIFF_CALC(new_list_end_floats[falsely_assigned[b][len(falsely_assigned[b])-1],:], STATS[c,:], falsely_assigned[b], MASTER_Basepairs_bonds[c])
    if false_counter > 0:
        print "NUMBER OF FALSE BASEPAIRS WITH", MASTER_Basepairs_schemes[c], "geometry =", false_counter, "\n"

#UNDETERMINED BASEPAIRS
    undetermined_counter = 0
    for b in range (len(undetermined)):
        if MASTER_Basepairs_schemes[c] in undetermined[b][0]:
            if b == 0:
                print "LIST OF UNDETERMINED BASEPAIRS WITH", MASTER_Basepairs_schemes[c], "geometry:"
            undetermined_counter = undetermined_counter + 1
            list = DIFF_CALC(new_list_end_floats[undetermined[b][len(undetermined[b])-1],:], STATS[c,:], undetermined[b], MASTER_Basepairs_bonds[c])
    if undetermined_counter > 0:
        print "NUMBER OF UNDETERMINED BASEPAIRS WITH", MASTER_Basepairs_schemes[c], "geometry =", undetermined_counter, "\n"



print           "==============================================================================="
#### 7) LIST OF BASEPAIRS BY GEOMETRY END




print "\n\n#######################################################\n##################### STATISTICS ######################\n#######################################################\n"
#########################################################################
##########################################################################
##########################################################################
#SECTION STATISTICS                                  END
##########################################################################
##########################################################################
##########################################################################


##########################################################################
##########################################################################
##########################################################################
#SECTION FINAL OUTPUT
##########################################################################
##########################################################################
##########################################################################


#OUTPUT REGARDING BASES WITH POSSIBLE MULTIPLE CONTACTS
#print "\n\nThe following group of bases have been identified as possible basepairs but cannot be assigned to a single basepair interaction\nPlease visually check the assignment in Pymol.\n##### These bases will be displayed in RED #####."
file = [[], []]
pymol1 = directory + "/" + pdb_file_main + "_UNASSIGNABLE-basepairs" + "_PYMOL" + "_script.pml"
file[1] = open(pymol1, 'w')
line1 = "load " + pdb_file + '\n' + "color white, /" + pdb_file_main + "/*" + '\n' + "set stick_ball, on" + '\n' + "set stick_ball_ratio, 2.00000" + '\n'
file[1].write(line1)

for i in range (len(control)):
    line1 = [] #will carry messages for printing out
    if len(control[i]) > 4:
       processed_control = [] #Will shuttle control[i] sub-lists with more than 4 elements. These are arrays of candidate basepairs possibly involved in MULTIPLE CONTACTS that could not be assigned as legitimate basepairs
       list = [] #Will carry formatted basepairs, so that they can be properly outputed
       for j in range (len(control[i])):
           if control[i][j] != []:
              processed_control.append(control[i][j])
       for u in range (len(processed_control)): #To format the 2 strings for the current basepair, resid and resid #, into residresid #
           if u%2 == 0:
              search = '\d+'
           else:
              search = '\w+'
           formatted = SEARCH(processed_control[u], search)
           list.append(formatted)
       for v in range (len(list)/2):
           line = "color red, /" + pdb_file_main + "/*/*/" + list[v*2] + '\n'
           file[1].write(line)
           line = "show sticks, /" + pdb_file_main + "/*/*/" + list[v*2] + '\n'
           file[1].write(line)
           line1.append(list[v*2] + list[(v*2)+1])


print "########################################################\n#################### END OF PROGRAM ####################\n########################################################"

import sys
sys.exit()

#STILL TO DO:
#v.041910 #The following basepairs in Thermus 16S rRNA cannot be assigned.
        #Between A978 and G1316 of the X_AG geometry.
        #No helical elements continuous with the basepair
        #Involved in multiple interactions with A1318 via A978
        #Distances between A`1318/N1 and A`978/C8 = 4.3, A`1318/C2 and A`978/C4 = 4.2
                #These distances could be used to deconvolute these basepairs. May not work in other cases, though
        #Between G1139 and G1142 of the VI_GG geometry
        #No helical elements continuous with the basepair
        #Involved in multiple interactions with G1134 via G1142
        #Distance between G`1142/N9 and G`1134/N2 = 3.7
                #This distance
        #Problems with residue G1030A in Thermus 16S rRNA (2J02)

        #SOLVED: Does not detect basepair between C490 and C444. This one should be a three bond basepair with 2 bonds identified
        #Could it recognize base to backbone interactions?:
                #G481 06 to A451 O2* ---> 3.3A
                #G481 N2 to A451 O1P ---> 3.3A
                #G481 N1 to A451 O1P ---> 3.1A
                #Note that A451 is also pretty close to G481 backbone
                #See also the interaction between G485's N1 and N2 with A448's O2P

        #note that G447 and A487 form a sheared basepair in which the G's N3 and the A's N6 are more than 3.7A apart. The distance from A's N6 to 447's O2* is 3.0A, though. Is this a common theme in sheared basepairs?
                #Sheared between A441 and G493, from A's N6 to G's O2* 2.8A, whereas from A's N6 to G's N3 3.7A
        #SOLVED 051210: Wrong helical element running from C456 to C458 and from G474 to G476. At maximum CUTOFF of 3.6 only two basepairs should have been identified
                #3-bond base pair with all bonds identified between C456 and G475
                #3-bond base pair with  only two bonds identified between C455 and G476

        #051510:Possible new basepair between A441 and G493 (involving A`441/N6 and G`493/O2*, and A`441/N7 and G`493/N2)
        #051610 (BUILD pdb_distances_051210.py): Basepair between 978 and 1316 of the X_AG geometry cannot be deconvoluted due to a stacking noise from A 1318, which is detected as II_AA. C1-C1 distance should be able to take care of this.
                # You clicked /2J02_16SrRNA//A/G`1316/N2 -> (pk1)
                # You clicked /2J02_16SrRNA//A/A`978/N1 -> (pk1)
                # You clicked /2J02_16SrRNA//A/G`1316/N3 -> (pk1)
                # You clicked /2J02_16SrRNA//A/A`978/N6 -> (pk1)
        #051610 (SOLVED) (BUILD pdb_distances_051210.py): Basepair between G1022 and C1007 (C1-C1 = 11.2) cannot be deconvoluted from stacking noise between C1008 and G1022 (C1-C1 = 8.3). C1-C1 distances to sort out this noise before attempting deconvolution by MH??????
        #Basepairs C444:G490 and G445:C489 are not identified. Instead a basepair between C444 and C489 appears to be identified

        # Base-pair between positions G1030A and C1031 of T. thermophilus 16S rRNA not recognized, probably because of consecutive bases involved
###i = 47 in MASTER_Basepairs (XXXIII_GG)
#XXXI_GG = [] #XXXIII Base-pair between positions G1030A and C1031 of T. thermophilus 16S rRNA (E. coli numbering)
#Bond   G       G       Length Ave      Length Std      Attribute
#1      N2      N1      2.9             unk
#2      N3      N2      3.4             unk
#bonds = [2, "G", "G", "N2", "N1", "N3", "N2", "NA", "NA", "NA", "NA", "NA", "NA"]
#1GMASTER_Basepairs_bonds.append(bonds)

        #TRIPLE involving bases 64, 68, and 101. the program detects the a basepair between 69 and 101 and appends that between 68 and 64 as FALSE
#v.042110
#From Ralf 050310
#    Hi Anton,

#    > but how is 'neighboring' defined? in terms of distance? in terms of residue
#    > number?

#    Simply the relative positions in the list of residues.
#    Maybe I should have said "consecutive residues".
#    I was thinking that could tell you if you have a base-pair interaction
#    or a stacking interaction.
#    But I'm at the limit of my biochemistry knowledge here.

#    > how would this distinguish between atoms forming a basepair and
#    > stacking atoms? unless you have some other way to tell them apart
#    > neighboring bases within 5A can be part of a stacking interaction, or even
#    > be part of a neighboring helix.

#    From the hierarchical iotbx.pdb objects you can easily tell if
#    atoms are in different chains.

#    > thanks for the explanation on the objects. however, I am far from
#    > understanding how you build them. What is the name of the type of
#    > hierarchical objects that you use? so i can read something more about them.

#    We have thousands of objects of this kind. It is impossible to write
#    documentation for all of them. The way to work with an object (after
#    somebody told you they have the information you need) is to insert
#    "print dir(atom)" or "help(atom)" into the script and run it. It will
#    show you all the things the object contains. We use long variable
#    names to make them as self-explanatory as possible.

#    > maybe, I could try to finish the program the way I am writing it and later
#    > on translate it into objects, like you suggested. i need to learn how they
#    > work first.

#    Yes, sounds like a good idea.
#    But what you do here...

#              master.append([
#                atom_i.resid(),
#                resname_i,
#                atmname_i,
#                atom_j.resid(),
#                resname_j,
#                atmname_j,
#                distance])

#    "throw away everything I don't understand" and "create a list
#    of anonymous numbers", is causing a lot of problems. I think it will
#    pay to play with the atom objects a little, and to start right here
#    making meaningful objects where each number has a name someone else
#    can understand. A simple start would be

#     class atom_pair:

#       def __init__(self, atom_i, atom_j, distance):
#         self.atom_i = atom_i
#         self.atom_j = atom_j
#         self.distance = distance

#       def atoms_are_in_same_chain(self):
#         # returns True or False

#       def atoms_are_in_consecutive_residues(self)
#         # returns True or False

#    and

#     master.append(atom_pair(atom_i, atom_j, distance))

#    Method names like above would tell me what you need and I could help
#    filling in the details.

#    Ralf
