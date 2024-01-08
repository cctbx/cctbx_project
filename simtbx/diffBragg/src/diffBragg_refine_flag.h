#ifndef SIMTBX_DIFFBRAGG_REFINE_FLAG
#define SIMTBX_DIFFBRAGG_REFINE_FLAG
/*
enum refine_id {
    ROTX_ID = (1u << 0),
    ROTY_ID = (1u << 1),
    ROTZ_ID = (1u << 2),
    UCELL_A_ID = (1u << 3),
    UCELL_B_ID = (1u << 4),
    UCELL_C_ID = (1u << 5),
    UCELL_ALPHA_ID = (1u << 6),
    UCELL_BETA_ID = (1u << 7),
    UCELL_GAMMA_ID = (1u << 8),
    NCELLS_ID = (1u << 9),
    PANELZ_ID = (1u << 10),
    FCELL_ID = (1u << 11),
    LAMBDA_OFFSET_ID = (1u << 12),
    LAMBDA_SCALE_ID = (1u << 13),
    PANEL_ROTO_ID = (1u << 14),
    PANELX_ID = (1u << 15),
    PANELY_ID = (1u << 16),
    PANEL_ROTF_ID = (1u << 17),
    PANEL_ROTS_ID = (1u << 18),
    ETA_ID = (1u << 19),
    NCELLS_OFFDIAG_ID = (1u << 21),
    F_PRIME_F_DPRIME_ID = (1u << 22),
    DIFFUSE_ID = (1u << 23),
};
*/

enum refine_flag : uint32_t {
    REFINE_DIFFUSE = (1u << 0),
    REFINE_FP_FDP = (1u << 1),
    REFINE_UMAT1 = (1u << 2),
    REFINE_UMAT2 = (1u << 3),
    REFINE_UMAT3 = (1u << 4),
    REFINE_BMAT1 = (1u << 5),
    REFINE_BMAT2 = (1u << 6),
    REFINE_BMAT3 = (1u << 7),
    REFINE_BMAT4 = (1u << 8),
    REFINE_BMAT5 = (1u << 9),
    REFINE_BMAT6 = (1u << 10),
    REFINE_NCELLS1 = (1u << 11),
    REFINE_NCELLS2 = (1u << 12),
    REFINE_NCELLS3 = (1u << 13),
    REFINE_NCELLS_DEF = (1u << 14),
    REFINE_PANEL_ORIGIN1 = (1u << 15),
    REFINE_PANEL_ORIGIN2 = (1u << 16),
    REFINE_PANEL_ORIGIN3 = (1u << 17),
    REFINE_PANEL_ROT1 = (1u << 18),
    REFINE_PANEL_ROT2 = (1u << 19),
    REFINE_PANEL_ROT3 = (1u << 20),
    REFINE_FCELL = (1u << 21),
    REFINE_ETA = (1u << 22),
    REFINE_LAMBDA1 = (1u << 23),
    REFINE_LAMBDA2 = (1u << 24),
    REFINE_ICELL = (1u << 25)
};

constexpr unsigned int REFINE_BMAT = REFINE_BMAT1 | REFINE_BMAT2 | REFINE_BMAT3 | REFINE_BMAT4 | REFINE_BMAT5 | REFINE_BMAT6;
constexpr unsigned int REFINE_UMAT = REFINE_UMAT1 | REFINE_UMAT2 | REFINE_UMAT3;
constexpr unsigned int REFINE_NCELLS = REFINE_NCELLS1 | REFINE_NCELLS2 | REFINE_NCELLS3;
constexpr unsigned int REFINE_PANEL_ORIGIN = REFINE_PANEL_ORIGIN1 | REFINE_PANEL_ORIGIN2 | REFINE_PANEL_ORIGIN3;
constexpr unsigned int REFINE_PANEL_ROT = REFINE_PANEL_ROT1 | REFINE_PANEL_ROT2 | REFINE_PANEL_ROT3;
constexpr unsigned int REFINE_LAMBDA = REFINE_LAMBDA1 | REFINE_LAMBDA2;

#endif
