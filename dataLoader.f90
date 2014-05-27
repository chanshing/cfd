module InputData !Contiene datos de euler-1.dat
    real(8) U_inf, V_inf, T_inf, RHO_inf, MACH_inf, P_inf, C_inf
    real(8) FR, FMU, FGX, FGY, QH, FK, FCv, GAMA
	real(8) XREF, YREF
	real(8) FSAFE, CTE
    real(8) ETA_REFIN, HHMAX_REFIN, HHMIN_REFIN
	integer IRESTART, MAXITER, IPRINT, MOVIE, ITLOCAL, MOVING
    integer NGAS
    character(4) RHOCHAR, VEL2CHAR, MACHCHAR, PRESCHAR, TEMPCHAR, ENERCHAR, POSCHAR
	character(80) FILENAME 
contains
subroutine readInputData
	implicit none
    open (1, FILE = 'EULER.DAT', STATUS = 'OLD')
    read (1, '(A)') FILENAME
    CLOSE(1)

    !CCCC-------------------------------------------------------!CCCC
    !CCCC----> LECTURA DE LOS DATOS GENERALES DEL PROBLEMA <----!CCCC
    !CCCC-------------------------------------------------------!CCCC
    open(1, FILE = trim(FILENAME)//'-1.dat', STATUS = 'OLD')
    !CCCC  ----->
    read(1, *)
    read(1, *) IRESTART, MAXITER, IPRINT, MOVIE, ITLOCAL
    !CCCC  ----->
    read(1, *)
    read(1, *) FSAFE, U_inf, V_inf, MACH_inf, T_inf, RHO_inf, P_inf
    !CCCC  ----->
    read(1, *)
    read(1, *) FMU, FGX, FGY, QH
    !CCCC  ----->
    read(1, *)
    read(1, *) FK, FR, FCv, GAMA, NGAS
    !CCCC  -----> CTE DE SHOCK-CAPTURING
    read(1, *)
    read(1, *) CTE
    !CCCC  -----> ALE
    read(1, *)
    read(1, *)
    read(1, *) MOVING, XREF, YREF
    !CCCC  -----> IMPRESION DE RESULTADOS SEGUN AGUSTIN
    read(1, *)
    read(1, *)
    read(1, *) RHOCHAR, VEL2CHAR, MACHCHAR, PRESCHAR, TEMPCHAR, ENERCHAR, POSCHAR
    !CCCC  -----> REFINANMIENTO
    read(1, *)
    read(1, *)
    read(1, *)
    read(1, *)
    read(1, *) ETA_REFIN, HHMAX_REFIN, HHMIN_REFIN
    CLOSE(1)
	
    CTE = 1.D0/CTE
    if(T_inf.EQ.0.D0) T_inf = P_inf/(FR*RHO_inf)
    if(P_inf.EQ.0.) P_inf = RHO_inf*FR*T_inf
    if(RHO_inf.EQ.0.D0) RHO_inf = P_inf/(FR*T_inf)
    C_inf = dsqrt(GAMA*FR*T_inf)
    if(dsqrt(U_inf**2.D0 + V_inf**2.D0) .EQ. 0.D0) U_inf = C_inf*MACH_inf
end subroutine readInputData
end module InputData

module MeshData !Contiene mvariabfix, mgeometria, mallocar, mfuerzas, mmovimiento & mnormales
	! = = = = = = = = = = Geometria  = = = = = = = = = ==
    integer, dimension(:, :), allocatable :: inpoel
    real(8), dimension(:, :), allocatable :: dNx, dNy
    real(8), dimension(:), allocatable :: X, Y, HHX, HHY, area, HH, M
    integer npoin, nelem, nFixRho, nFixVi, nFixV, nelNorm, nFixT, nsets, nFix_move, nmove
    integer NSLAVE, NMASTER 
	! = = = = = = = = = = Fixed  = = = = = = = = = ==
    integer, dimension(:), allocatable :: IFIXRHO_NODE, IFIXV_NODE, IFIXT_NODE
    real(8), dimension(:), allocatable :: RFIXRHO_VALUE, RFIXV_VALUEX
    real(8), dimension(:), allocatable :: RFIXV_VALUEY, RFIXT_VALUE
	! = = = = = = = = = = Normales  = = = = = = = = = =
    integer, dimension(:, :), allocatable :: IVN
	! = = = = = = = = = = Set Fuerzas  = = = = = = = = = ==
    integer, dimension(:, :, :), allocatable :: ISET
    integer, dimension(:), allocatable :: IPER_MASTER, IPER_SLAVE, IPER_AUX
	integer IELEM_SETS(10)
	integer NELE_SET, NSET_NUMB
	! = = = = = = = = = = MOVIMIENTO  = = = = = = = = = =
    integer, dimension(:), allocatable :: IFM, I_M, ILAUX
contains
subroutine loadMeshData
	use InputData
	!integer inod, ierror, ifixrho, ifixv, ielem, iel, im, imaster, ifixt
	!integer inel, twall, ifixt_a, innod, iset1, iset2, isetnumb, isets
	!real(8) dummy 

    open(1, FILE = trim(FILENAME)//'.dat', STATUS = 'OLD')

    read(1, *)      !NUMERO DE NODOS   !NUMERO DE ELEMENTOS
    read(1, *)       npoin              , nelem
    read(1, *)
    read(1, *) NFIXRHO, NFIXVI, NFIXV, NELNORM, NFIXT, NSETS, NMASTER, NSLAVE, NFIX_MOVE, NMOVE

    CALL ALLOC
  
    read(1, *)
    read(1, *)
    read(1, *)
    read(1, *)
    !CCCC-----> COORDENADAS DE LOS NODOS
    write(*, *) ' reading coordinates....'
    IERROR = 1
    do INOD = 1, npoin
        read(1, *, ERR = 27) INNOD, X(INNOD), Y(INNOD)
    end do
    write(*, *) ' ready '
    write(*, '(A, I5//)') ' TOTAL NODOS LEIDOS:', INNOD
  
    !CCCC----> CONECTIVIDADES DE LOS ELEMENTOS
    write(*, *) ' reading elements....'
    IERROR = 2
    read(1, *)
    do IELEM = 1, nelem
        read(1, *, ERR = 27) IEL, inpoel(1, IEL), inpoel(2, IEL), inpoel(3, IEL)
    end do
    write(*, *) ' ready '
    write(*, '(A, I6//)') ' TOTAL ELEMENTOS LEIDOS:', IEL
  
    !CCCC----> PUNTOS CON DENSIDAD PRESCRITA
    write(*, *) 'reading fix density points....'
    IERROR = 3
    read(1, *)
    do IFIXRHO = 1, NFIXRHO
        read(1, *, ERR = 27) IFIXRHO_NODE(IFIXRHO), RFIXRHO_VALUE(IFIXRHO)
        if (RFIXRHO_VALUE(IFIXRHO).LT.0) THEN
            RFIXRHO_VALUE(IFIXRHO) = 1.225d0
        ELSE
            RFIXRHO_VALUE(IFIXRHO) = RFIXRHO_VALUE(IFIXRHO)*RHO_inf
        end if
    end do
    write(*, *) ' ready '
    write(*, '(A, I5//)') ' TOTAL NODOS CON DENSIDAD IMPUESTA:', NFIXRHO
  
    !CCCC----> NODOS CON VELOCIDAD FIJA
    write(*, *) ' reading fix velocities (INFLOW)....'
    IERROR = 4
    read(1, *)
    do IFIXV = 1, NFIXVI
        read(1, *, ERR = 27) IFIXV_NODE(IFIXV), RFIXV_VALUEX(IFIXV) &
            , RFIXV_VALUEY(IFIXV)
        RFIXV_VALUEX(IFIXV) = RFIXV_VALUEX(IFIXV)*U_inf
        RFIXV_VALUEY(IFIXV) = RFIXV_VALUEY(IFIXV)*V_inf
    end do
    write(*, *) ' ready '
    write(*, '(A, I5//)') ' TOTAL ELEMENTOS CON VELOCIDAD IMPUESTA:', NFIXVI
  
    !CCCC----> NODOS CON VELOCIDAD FIJA (NO SLIP)
    !CCCC----> CALCULA LA TEMPERATURA DE ESTANCAMIENTO PARA IMPONERLA EN LAS
    !CCCC----> PAREDES DONDE SE PRESCRIBE LA CONDICION DE "NO SLIP"
    TWALL = T_inf*(1.D0 + (GAMA-1)/2.D0*MACH_inf*MACH_inf)
    IFIXT_A = 0
    write(*, *) ' reading fix velocities (NO SLIP)....'
    IERROR = 4
    read(1, *)
    do IFIXV = NFIXVI + 1, NFIXVI + NFIXV
        read(1, *, ERR = 27) IFIXV_NODE(IFIXV), RFIXV_VALUEX(IFIXV) &
            , RFIXV_VALUEY(IFIXV)
        RFIXV_VALUEX(IFIXV) = 0.D0
        RFIXV_VALUEY(IFIXV) = 0.D0
        IFIXT_A = IFIXT_A + 1
        IFIXT_NODE(IFIXT_A) = IFIXV_NODE(IFIXV)
        RFIXT_VALUE(IFIXT_A) = TWALL
    end do
    write(*, *) ' ready '
    NFIXV = NFIXV + NFIXVI
    write(*, '(A, I5//)') ' TOTAL ELEMENTOS CON VELOCIDAD IMPUESTA:', NFIXV
  
    !CCCC----> VELOCIDAD NORMAL NULA
    write(*, *) ' reading elements with normal velocity prescribe....'
    IERROR = 5
    read(1, *)
    do INEL = 1, NELNORM
        read(1, *, ERR = 27) IVN(1, INEL), IVN(2, INEL)
    end do
    write(*, *) ' ready '
    write(*, '(A, I6/)')' TOTAL ELEMENTOS CON VELOCIDAD NORMAL IMPUESTA:', NELNORM
  
    !CCCC----> TEMPERATURA PRESCRITA
    write(*, *) ' reading fix temperature nodes....'
    IERROR = 8
    read(1, *)
    do IFIXT = 1, NFIXT
        read(1, *, ERR = 27) IFIXT_NODE(IFIXT + IFIXT_A), RFIXT_VALUE(IFIXT + IFIXT_A)
        RFIXT_VALUE(IFIXT + IFIXT_A) = RFIXT_VALUE(IFIXT + IFIXT_A)*T_inf
    end do
    write(*, *) ' ready '
    NFIXT = NFIXT + IFIXT_A
    write(*, '(A, I6/)') ' TOTAL NODOS COM TEMPERATURA IMPUESTA:', NFIXT
  
    !CCCC----> LECTURA DE LOS SETS PARA EL
    !CCCC----> CALCULO DE LAS FUERZAS
    write(*, *) ' reading sets nodes....'
    IERROR = 9
	IELEM_SETS = 0
    NSET_NUMB = 0
    read(1, *)
    do ISETS = 1, NSETS
        read(1, *, ERR = 27) NELE_SET, ISET1, ISET2, ISETNUMB
        if (ISETNUMB.GT.NSET_NUMB) NSET_NUMB = ISETNUMB
        IELEM_SETS(ISETNUMB) = IELEM_SETS(ISETNUMB) + 1
        ISET(1, ISETNUMB, IELEM_SETS(ISETNUMB)) = ISET1
        ISET(2, ISETNUMB, IELEM_SETS(ISETNUMB)) = ISET2
        ISET(3, ISETNUMB, IELEM_SETS(ISETNUMB)) = NELE_SET
    end do
    write(*, *) ' ready '
    write(*, '(A, I6)') ' NUMERO DE SETS:', NSET_NUMB
    write(*, '(A, I6//)') ' NUMERO DE ELEMENTOS DE LOS SETS:', NSETS
  
    !CCCC----> PERIODICAL MASTER
    write(*, *) ' reading nodes with periodical master and slave condition....'
    IERROR = 10
    if(NMASTER.NE.NSLAVE)THEN
        write(*, *)'ERROR NODOS MASTER DISTINTO NODOS SLAVE'
        stop
    end if
    read(1, *)
    do IMASTER = 1, NMASTER
        read(1, *, ERR = 27)IPER_MASTER(IMASTER)
    end do
    read(1, *)
    do IMASTER = 1, NMASTER
        read(1, *, ERR = 27)IPER_SLAVE(IMASTER)
    end do
    write(*, *) ' ready '
    write(*, '(A, I6/)')' TOTAL NODOS CON CONDICION PERIODICA:', NMASTER + NSLAVE
  
    !CCCC----> NODOS CON MOVIMIENTO FIJO
    write(*, *) ' reading fix movement ....'
    read(1, *)
    do IM = 1, NFIX_MOVE
        read(1, *, ERR = 27) IFM(IM), dummy
    end do
    write(*, *) ' ready '
    write(*, '(A, I5//)') ' TOTAL NODOS CON MOVIMIENTO FIJO:', NFIX_MOVE

    !CCCC----> NODOS CON MOVIMIENTO
    !CCCC--------------------------
    write(*, *) ' reading movement ....'
    read(1, *)
    do IM = 1, NMOVE
        read(1, *, ERR = 27) I_M(IM), dummy
    end do
    write(*, *) ' ready '
    write(*, '(A, I5//)') ' TOTAL NODOS CON MOVIMIENTO FIJO:', NMOVE

	close(1)
	IERROR = 0
27  if (IERROR.NE.0) THEN
        CALL READ_ERROR(IERROR)
    end if
end subroutine loadMeshData

subroutine read_error(IERROR)
    implicit real(8) (a-h, o-z)
    character(10) errType(5)

    data errType /'NODOS', 'ELEMENTOS', 'NO_TRAC', 'FIX_VEL', 'NORM_VEL'/
    if (IERROR.ne.0) then
        write(*, '(A)') 'ERROR EN LA LECTURA DE ', errType(IERROR)
        stop
    end if
end subroutine read_error 

subroutine ALLOC
    !use DATOS_REFINAMIENTO
    !use DATOS_ENTRADA
    !use MNORMALES
    !use MALLOCAR
    !use MVARIABFIX
    !use MGEOMETRIA
	!use MeshData
    !use MFUERZAS
    !use MMOVIMIENTO
	!use InputData
    use MVELOCIDADES
    use MVARIABGEN
    use MVARIABLES
    use MESTABILIZACION
    use MNEWMARK
    ALLOCATE(VEL_X(npoin), VEL_Y(npoin), W_X(npoin), W_Y(npoin))
    ALLOCATE(U(4, npoin), U1(4, npoin), U2(4, npoin), RHS(4, npoin), RHS1(4, npoin), UN(4, npoin))
    ALLOCATE(RHS2(4, npoin), RHS3(4, npoin))
    !ALLOCATE(RNORMV_VALUEX(npoin), RNORMV_VALUEY(npoin), IVN(2, NELNORM))
	allocate(ivn(2, nelnorm))
    !ALLOCATE(INORMV_NODE(npoin), RNX(npoin), RNY(npoin), BAUX(npoin))
    ALLOCATE(IFIXRHO_NODE(NFIXRHO), RFIXRHO_VALUE(NFIXRHO), RFIXV_VALUEX(NFIXVI + NFIXV))
    ALLOCATE(RFIXV_VALUEY(NFIXVI + NFIXV), IFIXV_NODE(NFIXVI + NFIXV))
    ALLOCATE(RFIXT_VALUE(NFIXV + NFIXT), IFIXT_NODE(NFIXT + NFIXV))
    ALLOCATE(X(npoin), Y(npoin), inpoel(3, NELEM), HHX(NELEM), HHY(NELEM), DNX(3, NELEM))
    ALLOCATE(DNY(3, NELEM), AREA(NELEM), HH(NELEM), M(npoin))
    ALLOCATE(P(npoin), T(npoin), RHO(npoin), E(npoin), RMACH(npoin))
    !ALLOCATE(S(npoin*15), NN1(npoin*15), NN2(npoin*15), IND(npoin), INDEL(10, npoin))
    !ALLOCATE(ADIAG(npoin), RAUX(npoin))
    ALLOCATE(ISET(3, 10, NSETS), IPER_MASTER(NMASTER), IPER_SLAVE(NMASTER))
    ALLOCATE(IPER_AUX(NMASTER), IFM(NFIX_MOVE), I_M(NMOVE), ILAUX(npoin))
    ALLOCATE(SHOC(NELEM), T_SUGN1(NELEM), T_SUGN2(NELEM), T_SUGN3(NELEM))
    ALLOCATE(DXPOS(NMOVE), DYPOS(NMOVE))
end subroutine ALLOC
end module MeshData
