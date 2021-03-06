MODULE MVELOCIDADES
    REAL(8), DIMENSION(:), ALLOCATABLE:: VEL_X, VEL_Y, W_X, W_Y
END MODULE MVELOCIDADES

MODULE MVARIABGEN
    REAL(8), DIMENSION(:, :), ALLOCATABLE:: U, U1, U2, RHS, RHS1, RHS2, RHS3, UN
END MODULE MVARIABGEN

MODULE MVARIABLES
    REAL(8), DIMENSION(:), ALLOCATABLE:: P, T, RHO, E, RMACH
END MODULE MVARIABLES

MODULE MESTABILIZACION
    !REAL(8) CTE
    REAL(8) , DIMENSION(:), ALLOCATABLE:: SHOC, T_SUGN1, T_SUGN2, T_SUGN3
END MODULE MESTABILIZACION

!MODULE MNEWMARK
    !REAL(8), DIMENSION(:), ALLOCATABLE:: DXPOS, DYPOS
!END MODULE MNEWMARK

!MODULE MATRICES
    !INTEGER NDIM
    !REAL(8), DIMENSION(:, :), ALLOCATABLE:: MASA, C, K
!END MODULE MATRICES

!MODULE MAT2
    !REAL(8), DIMENSION(:), ALLOCATABLE:: DISN, DISS, VELN, VELS, ACCN, ACCS, R, G, G1, F, CTEN
!END MODULE MAT2

MODULE TIMERS
    integer rate, start_t, end_t, sub_start, sub_end, sub_rate, m_start, m_end, m_rate
    real calcrhs_t, cuarto_t, output_t, masas_t, deriv_t, laplace_t, normales_t, forces_t, newmark_t, grad_t, transf_t, estab_t
	real spmv_t, residuo_t, fuente_t, total_t
END MODULE TIMERS

!MODULE DATOS_REFINAMIENTO
!    REAL(8)ETA_REFIN, HHMAX_REFIN, HHMIN_REFIN
!END MODULE DATOS_REFINAMIENTO
!
!MODULE DATOS_ENTRADA
!    INTEGER NGAS
!    REAL(8) GAMA, UINF, VINF, TINF, RHOINF, MACHINF, PINF, CINF
!    REAL(8) FR, FMU, FGX, FGY, QH, FK, FCv
!END MODULE DATOS_ENTRADA
!
!MODULE MNORMALES
!    INTEGER   NNORMV
!    INTEGER , DIMENSION(:, :), ALLOCATABLE:: IVN
!    INTEGER , DIMENSION(:), ALLOCATABLE:: INORMV_NODE
!    REAL(8) , DIMENSION(:), ALLOCATABLE:: RNORMV_VALUEX, RNORMV_VALUEY, RNX, RNY, BAUX
!END MODULE MNORMALES
!
!MODULE MALLOCAR
!    INTEGER NNOD, NELEM, NFIXRHO, NFIXVI, NFIXV, NELNORM, NFIXT, NSETS, NMASTER
!    INTEGER NSLAVE, NFIX_MOVE, NMOVE
!END MODULE MALLOCAR

!MODULE MVARIABFIX
!    INTEGER , DIMENSION(:), ALLOCATABLE:: IFIXRHO_NODE, IFIXV_NODE, IFIXT_NODE
!    REAL(8) , DIMENSION(:), ALLOCATABLE:: RFIXRHO_VALUE, RFIXV_VALUEX
!    REAL(8) , DIMENSION(:), ALLOCATABLE:: RFIXV_VALUEY, RFIXT_VALUE
!END MODULE MVARIABFIX
!
!MODULE MGEOMETRIA
!    INTEGER , DIMENSION(:, :), ALLOCATABLE:: N
!    REAL(8) , DIMENSION(:, :), ALLOCATABLE:: DNX, DNY
!    REAL(8) , DIMENSION(:), ALLOCATABLE:: X, Y, HHX, HHY, AREA, HH, M
!END MODULE MGEOMETRIA

!MODULE MFUERZAS
!    INTEGER , DIMENSION(:, :, :), ALLOCATABLE:: ISET
!    INTEGER , DIMENSION(:), ALLOCATABLE:: IPER_MASTER, IPER_SLAVE, IPER_AUX
!	integer ielem_sets(10)
!END MODULE MFUERZAS

!MODULE MMOVIMIENTO
!    INTEGER , DIMENSION(:), ALLOCATABLE:: IFM, I_M, ILAUX
!END MODULE MMOVIMIENTO

!MODULE MLAPLACE
!    INTEGER , DIMENSION(:, :), ALLOCATABLE:: INDEL
!    INTEGER , DIMENSION(:), ALLOCATABLE:: IND, NN1, NN2
!    REAL(8) , DIMENSION(:), ALLOCATABLE:: ADIAG, RAUX, S
!END MODULE MLAPLACE

!MODULE MPRINTRES
!    CHARACTER*4 RHOCHAR, VEL2CHAR, MACHCHAR, PRESCHAR, TEMPCHAR, ENERCHAR, POSCHAR
!END MODULE MPRINTRES
