!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCC                                                       CCC  
!CCC                   NAVIER-STOKES 2D                    CCC     
!CCC                                                       CCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
PROGRAM NSComp2D
    !use DATOS_REFINAMIENTO
    !use DATOS_ENTRADA
    use Mnormales
    !use MALLOCAR
    use MVELOCIDADES
    use MVARIABGEN
    !use MVARIABFIX
    !use MGEOMETRIA
    use MVARIABLES
    !use MFUERZAS
    !use MMOVIMIENTO
    use MESTABILIZACION
    !use MPRINTRES
    use MNEWMARK
    use MATRICES
    use MAT2
    use TIMERS
	use MelementSurrPoint
	use MpointSurrPoint
	use Mlaplace
	use MbiconjGrad
	use InputData
	use MeshData
    implicit real(8) (A-H, O-Z)

    integer BANDERA, NESTAB
    ! CCCC SUAVIZADO DE LOS ELEMENTOS
    integer, dimension(:, :), allocatable:: SMOOTH_SIM
    real(8) FX(10), FY(10), RM(10), F_VX(10), F_VY(10)
    real(8), dimension(:), allocatable:: GAMM
    real(8) ER(4), ERR(4)
    ! CCCC LOCAL TIME STEP
    real(8), dimension(:), allocatable:: DTL, DT
    ! CCCC NODOS VECINOS Y LAPLACIANO
    real(8), dimension(:), allocatable:: B
    real(8), dimension(:), allocatable:: XPOS, YPOS, POS_AUX
    ! CCCC REFINAMIENTO
    real(8), dimension(:), allocatable:: HH_NEW
    ! CCCC SUAVIZADO SMOOTHING
    logical, dimension(:), allocatable:: SMOOTH_FIX

	cuarto_t=0; calcrhs_t = 0; output_t=0; masas_t=0; deriv_t=0; laplace_t=0; normales_t=0; fuente_t=0
	forces_t=0; newmark_t=0; transf_t=0; grad_t=0; residuo_t=0; spmv_t=0

    !CCCC----> TIEMPO DE CPU
	call system_clock(m_start, m_rate)

    open(22, FILE='FORCES_HISTORY', STATUS='UNKNOWN')

	call readInputData
	call loadMeshData

    ALLOCATE(SMOOTH_SIM(2, npoin))
    ALLOCATE(GAMM(npoin), DTL(NELEM), DT(NELEM), B(npoin))
    ALLOCATE(XPOS(npoin), YPOS(npoin), POS_AUX(npoin))
    ALLOCATE(HH_NEW(npoin), SMOOTH_FIX(npoin))

    !CCCC------------------------------------!CCCC
    !CCCC----->  COMIENZO DEL CALCULO  <-----!CCCC
    !CCCC------------------------------------!CCCC
    !CCCC----> ASIGNA VALOR A LAS VARIABLES PRIMITIVAS Y REALIZA UN RESTART
    !CCCC----> SI IRESTART.EQ.1
    !CCCC------------------------------------------------------------
  
    GAMM = GAMA
    call RESTART(GAMM)

    if(NGAS.NE.1) GAMM=GAMA
 
    !CCCC-----> CALCULO DE LOS NODOS CON PERIODICIDAD
    !CCCC---------------------------------------------------------
    if(NMASTER.NE.0.AND.NSLAVE.NE.0)THEN
        call PERIODIC
    end if

    !CCCC-----> CALCULO DE LAS NORMALES
    !CCCC---------------------------------------------------------
    call NORMALES

    !CCCC----> CALCULO DE LAS DERIVADAS, EL AREA Y LA LONGITUD CARAC.
    !CCCC------------------------------------------------------------
    call DERIV(HMIN)

    !CCCC----> CALCULO DE LAS MATRICES DE MASAS LUMPED
    !CCCC  ----------------------------------------------------------
    call MASAS

    !CCCC----> CALCULO DE LOS NODOS VECINOS Y EL LAPLACIANO
    !CCCC  ------------------------------------------------------
	call elementSurrPoint(inpoel, nelem, npoin)
	call pointSurrPoint(inpoel, nelem, npoin)

	call laplace(inpoel, area, dNx, dNy, nelem, npoin)
	npos = size(lap_sparse)

    !!$ ! FIJO QUE NODOS SE MUEVEN Y CUALES NO..
    NNMOVE=NFIX_MOVE+NMOVE
    KK=0
    do I=1, NMOVE
        KK=KK+1
        N1=I_M(I)
        ILAUX(KK)=N1
    end do
  
    do I=1, NFIX_MOVE
        KK=KK+1
        N1=IFM(I)
        ILAUX(KK)=N1
    end do

    !!$ ! NODOS QUE NO SON SUAVIZADOS
    SMOOTH_FIX(1:npoin)=.TRUE.
    SMOOTH_SIM(1:2, 1:npoin)=1
    do I=1, NMOVE
        N1=I_M(I)
        SMOOTH_FIX(N1)=.FALSE.
    end do
     
    do I=1, NFIX_MOVE
        N1=IFM(I)
        SMOOTH_FIX(N1)=.FALSE.
    end do

    !CCCC----> ABRE EL ARCHIVO DE CONVERGENCIA
    !CCCC----------------------------------------------
    open(7, FILE=trim(filename)//'.cnv', STATUS='UNKNOWN')

    !CCCC -----------------------------------------CCCC
    !CCCC ----> COMIENZO DE LAS ITERACIONES <----  CCCC
    !CCCC -----------------------------------------CCCC
    TIME = 0.D0
    ISAL = 0
    ITER = 0
    DTMIN = 0.D0
    ITERPRINT = 0
    FLUX1 = DABS(FGX+FGY+QH)          !VARIABLE PARA CALCULAR LOS TERMINOS FUENTES
    RMACH = (U_inf**2+V_inf**2)/DSQRT(GAMA*FR*T_inf)
    !CCCC---->  ORDEN DE INTEGRACION DE RUNGE-KUTTA
    NRK = 4
  
    ! ****DIMENSION DEL SISTEMA****
    NDIM=2
    call ALLOC1

    !CCCC----> DEFINO VELOCIDAD DE LA MALLA
    W_X=-0.D0 ; W_Y=0.D0
    YPOS=0.D0 ; XPOS=0.D0 ; IMOOTH=0

    write(*, '(A, I3, A)') '****-------> RUNGE-KUTTA DE', NRK, '  ORDEN <-------****'
    write(*, *)''
  
    !CCCC----> ABRE EL ARCHIVO DONDE SE IMPRIMEN LOS RESULTADOS
    if (MOVIE.EQ.1)THEN
        open(2, FILE=trim(filename)//'.flavia.res', STATUS='UNKNOWN')
    end if
  
    BANDERA=1
    NESTAB=1
  
    !!CONDICIONES INICIALES NEWMARK!!!
    DISN(1)=0.D0; DISN(2)=0.D0
    VELN(1)=0.D0; VELN(2)=0.D0
    ACCN=0.D0
    F(1)=0.d0
    F(2)=0.d0
    pi=dacos(-1.d0)
    YPOSR=DISN(1)
    ALPHA=DISN(2)
   
    do WHILE (ITER.LT.MAXITER.AND.ISAL.EQ.0)
 
        ITER=ITER+1

        !CCCC  ----> CALCULO DEL dT
        !CCCC  --------------------
     
        call DELTAT(DTMIN, DT)
      
        if (BANDERA.EQ.1) THEN
            DTMIN1=DTMIN
            BANDERA=2
        end if
        PORC=DABS((DTMIN-DTMIN1)/DTMIN)
        if (100.D0*PORC.LE.1.D0) THEN
            DTMIN=DTMIN1
        else
            DTMIN1=DTMIN
            BANDERA=2
        end if

        !CCCC----> FUNCION PARA ESTABLECER LOCAL-TIME-STEP
        if (ITLOCAL.NE.0) THEN
            DTFACT = 1.D0 - DEXP(-ITER*4.6D0/ITLOCAL)
            DTL= DTMIN*DTFACT+DT*(1.D0-DTFACT)
        else
            DTL=DTMIN
        end if
     
        TIME=TIME+DTMIN
     
        U1 = U
    
        if (BANDERA.LE.4) THEN
            call RK(DTMIN, NRK, BANDERA, GAMM, dtl)
        else
            call ADAMSB(DTMIN, NESTAB, GAMM, dtl)
        end if
   
        !CCCC-------------------------------------------------CCCC
        !CCCC  ----> CALCULO DE LOS ERRORES (RESIDUOS) <----  CCCC
        !CCCC  ----> PARA CONTROLAR LA CONVERGENCIA    <----  CCCC
        !CCCC-------------------------------------------------CCCC
        ER = 0.D0; ERR = 0.D0
        do INOD=1, npoin
            ER(:) = ER(:) + (U(:, INOD) - U1(:, INOD))**2
            ERR(:) = ERR(:) + U1(:, INOD)**2.D0
        end do

        U = U1
     
		call system_clock(start_t, rate)
        write(7, '(I7, 4E14.6)') ITER, TIME, DSQRT(ER(1)/ERR(1)), DSQRT(ER(4)/ERR(4))
		call system_clock(end_t)
		output_t = output_t + real(end_t-start_t)/real(rate)

		call system_clock(start_t, rate)
        call FORCES(NSETS, NSET_NUMB, IELEM_SETS, npoin &
            , X, Y, P, FX, FY, RM, X_REF, Y_REF)
		call system_clock(end_t)
		forces_t = forces_t+ real(end_t-start_t)/real(rate)
     
        F(1) = DCOS(DISN(2))*FY(1) - DSIN(DISN(2))*FX(1)
        F(2) = RM(1) + 0.1*DSIN(DISN(2))*FX(1) - DCOS(DISN(2))*0.1*FY(1)

        ALPHAV = DISN(2); YPOSRV = DISN(1)
		call system_clock(start_t, rate)
        call NEWMARK_METHOD(DTMIN)
		call system_clock(end_t)
		newmark_t = newmark_t + real(end_t-start_t)/real(rate)
     
        ALPHA = DISN(2) - ALPHAV; YPOSR = DISN(1) - YPOSRV
        !open(17, FILE='DESPLAZAMIENTO', STATUS='UNKNOWN')
        !write(17, *)ITER, DISN(1), DISN(2)

		call system_clock(start_t, rate)
        call TRANSF(ALPHA, YPOSR, XREF, YREF)
		call system_clock(end_t)
		transf_t = transf_t + real(end_t-start_t)/real(rate)
      
        POS_AUX(1:NMOVE) = DXPOS(1:NMOVE)
         !NODOS SIN MOVIMIENTO
        POS_AUX(1 + NMOVE : NFIX_MOVE + NMOVE) = 0.D0
      
        B = 0.D0
		call system_clock(start_t, rate)
		call biconjGrad(lap_sparse, lap_idx, lap_rowptr, lap_diag, xpos, b, pos_aux(1:nnmove), ilaux(1:nnmove), &
				npoin, nnmove)
		call system_clock(end_t)
		grad_t = grad_t + real(end_t-start_t)/real(rate)
             
        X = X + XPOS
        W_X = XPOS/DTMIN
     
        POS_AUX(1:NMOVE) = DYPOS(1:NMOVE)
     
        B = 0.D0
		call system_clock(start_t, rate)
		call biconjGrad(lap_sparse, lap_idx, lap_rowptr, lap_diag, ypos, b, pos_aux(1:nnmove), ilaux(1:nnmove), &
				npoin, nnmove)
		call system_clock(end_t)
		grad_t = grad_t + real(end_t-start_t)/real(rate)
	
        Y = Y + YPOS
        W_Y = YPOS/DTMIN
     
        !CCCC--------------------------CCCC
        !CCCC  ----> IMPRESION <----   CCCC
        !CCCC--------------------------CCCC
        ITERPRINT = ITERPRINT + 1
        if (ITERPRINT.EQ.IPRINT.OR.ITER.EQ.MAXITER) THEN
            !CCCC----> Fusible para que corte si error>1.d2
            FUSIBLE=MAXVAL(DSQRT(ER/ERR))
        
            if(FUSIBLE.GT.1.D2)THEN
                print*, '      ERROR CONVERGENCIA'
                print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*, '    *****  OVERFLOW  *****'
                print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                STOP
            end if
        
            write(*, '(A, 3X)')'CCCC-----------------------------------------CCCC '
            write(*, '(A, 3X)')'CCCC  ----> INFORMACION DE LA CORRIDA <----  CCCC '
            write(*, '(A, 3X/)')'CCCC-----------------------------------------CCCC '
            write(*, '(A, I6/)')'PASOS EJECUTADOS:', ITER
            write(*, '(A, E12.4)')'TIEMPO ACUMULADO:', TIME
            write(*, '(A, E12.4)')'PASO DE TIEMPO:', DTMIN
            write(*, '(A, F7.2/)')'NUMERO DE MACH MAXIMO:', MAXVAL(RMACH)
        
            write(*, '(5X, A)')'******   ERRORES   *******'
            write(*, '(A, X, E12.4)')'Continuidad', DSQRT(ER(1)/ERR(1))
            write(*, '(A, X, E12.4)')'Momento u  ', DSQRT(ER(2)/ERR(2))
            write(*, '(A, X, E12.4)')'Momento v  ', DSQRT(ER(3)/ERR(3))
            write(*, '(A, X, E12.4, /)')'Energia    ', DSQRT(ER(4)/ERR(4))
        
            call PRINTFLAVIA(GAMM, RHO, VEL_X-W_X, VEL_Y-W_Y, P, T, E, RMACH, X, Y &
                , npoin, ITER)
        
            ITERPRINT=0
        
            if(FMU.NE.0.D0)THEN
                call FORCE_VISC(NELEM, npoin, IFF &
                    , IELEM_SETS, ISET, inpoel, NSET_NUMB, U_inf, V_inf, RHO_inf, T_inf &
                    , X, Y, P, T, VEL_X, VEL_Y, DNX, DNY, FMU, RHO &
                    , F_VX, F_VY)
            end if
        
            open(33, FILE='FORCES', STATUS='UNKNOWN')
        
            do ISET_NUMB=1, NSET_NUMB
                write(33, '(A, I2)') 'SET NUMERO', ISET_NUMB
                write(33, '(A, E14.5)') 'FUERZA EN X:', FX(ISET_NUMB)
                write(33, '(A, E14.5/)') 'FUERZA EN Y:', FY(ISET_NUMB)
                !write(33, '(A, E14.5/)') 'MOMENTO:', RM(ISET_NUMB)
                write(33, '(A, E14.5)') 'FUERZA VISCOSA EN X:', F_VX(ISET_NUMB)
                write(33, '(A, E14.5/)') 'FUERZA VISCOSA EN Y:', F_VY(ISET_NUMB)
                write(33, '(A, E14.5)') 'FUERZA TOTAL EN X:', FX(ISET_NUMB)+F_VX(ISET_NUMB)
                write(33, '(A, E14.5/)') 'FUERZA TOTAL EN Y:', FY(ISET_NUMB)+F_VY(ISET_NUMB)
            end do
            close(33)
        
            call PRINTREST(ITER, npoin, U, T, GAMM, TIME, filename)
        
        end if

        BANDERA=BANDERA+1

        !CCCC  ----> CALCULOS DE LOS TERMINOS PARA CADA ITERACION
        if(MOVING.EQ.1)THEN
            !CCCC-----> CALCULO DE LAS NORMALES
            !CCCC---------------------------------------------------------
			call system_clock(start_t, rate)
            call NORMALES
			call system_clock(end_t)
			normales_t= normales_t+ real(end_t-start_t)/real(rate)

            !CCCC----> CALCULO DE LAS DERIVADAS, EL AREA Y LA LONGITUD CARAC.
            !CCCC------------------------------------------------------------
			call system_clock(start_t, rate)
            call DERIV(HMIN)
			call system_clock(end_t)
			deriv_t= deriv_t+ real(end_t-start_t)/real(rate)

            !CCCC----> CALCULO DE LAS MATRICES DE MASAS LUMPED
            !CCCC  ----------------------------------------------------------
			call system_clock(start_t, rate)
            call MASAS
			call system_clock(end_t)
			masas_t= masas_t+ real(end_t-start_t)/real(rate)

            !CCCC----> CALCULO DE LOS NODOS VECINOS Y EL LAPLACIANO
            !CCCC  ------------------------------------------------------
			call system_clock(start_t, rate)
			call laplace(inpoel, area, dNx, dNy, nelem, npoin)
			call system_clock(end_t)
			laplace_t = laplace_t + real(end_t-start_t)/real(rate)
        end if
    end do

    !CCCC  -----> CIERRA EL ARCHIVO DE RESULTADOS..
    if (MOVIE.EQ.1)THEN
        close(2)
    end if
    !CCCC-----------------------------------CCCC
    !CCCC  ----->    FIN DEL LOOP   <-----  CCCC
    !CCCC-----------------------------------CCCC
  
    !CCCC----> CIERRA EL ARCHIVO DE LA CONVERGENCIA
    !CCCC------------------------------------------
    write(*, *)'****---------------------****'
    write(*, *)'****----> REFINANDO <----****'
    write(*, *)'****---------------------****'

	print *, "cuarto_t: ", cuarto_t
	print *, "estab_t: ", estab_t
	print *, "calcrhs_t : ", calcrhs_t
	print *, "fuente_t: ", fuente_t
	print *, "output_t: ", output_t
	print *, "masas_t: ", masas_t
	print *, "deriv_t: ", deriv_t
	print *, "normales_t: ", normales_t
	print *, "laplace_t: ", laplace_t
	print *, "forces_t: ", forces_t
	print *, "newmark_t: ", newmark_t
	print *, "transf_t: ", transf_t
	print *, "grad_t: ", grad_t
	!print *, "spmv_t: ", spmv_t

    !!$  !call NEW_SIZE(npoin, NELEM, inpoel, DNX, DNY, AREA, M, HH, U &
    !!$  !     , HH_NEW)
    !do I=1, npoin
    !   RMACH(I)=DSQRT(VEL_X(I)**2+VEL_Y(I)**2)/DSQRT(GAMA*FR*T(I))
    !end do
    !hh_new=hhmax_refin

   
    !call ESTIMADOR_ERR(npoin, NELEM, inpoel, DNX, DNY, HH, M, AREA, 3.5D0, RMACH, P)
    !do I=1, npoin
    !   if(HH_NEW(I).GT.P(I))HH_NEW(I)=P(I)
    !end do

    !call ESTIMADOR_ERR(npoin, NELEM, inpoel, DNX, DNY, HH, M, AREA, 4.D0, RHO, P)
    !do I=1, npoin
    !   if(HH_NEW(I).GT.P(I))HH_NEW(I)=P(I)
    !end do

    !HHMIN_REFIN=0.01D0
    !call ESTIMADOR_ERR(npoin, NELEM, inpoel, DNX, DNY, HH, M, AREA, 6.D0, RMACH, P)
    !do I=1, npoin
    !   if(HH_NEW(I).GT.P(I))HH_NEW(I)=P(I)
    !end do

    ! ACOMODO LOS TAMANOS SOBRE EL/LOS CUERPOS
    !do ISET_NUMB=1, NSET_NUMB
    !   do II=1, IELEM_SETS(ISET_NUMB)
            
    !      N1=ISET(1, ISET_NUMB, II)
    !      N2=ISET(2, ISET_NUMB, II)
    !      IELEM=ISET(3, ISET_NUMB, II)
    !      HH_NEW(N1)=.01D0 !hhmin_refin!DSQRT(2.D0*AREA(IELEM))
    !      HH_NEW(N2)=.01D0 !hhmin_refin!DSQRT(2.D0*AREA(IELEM))
    !   end do
    !end do
  
    !do I=1, NELEM
    !   AR=AREA(I)
    !   CC=DSQRT(1.5*2.D0*AR)
    !   do J=1, 3
    !      N1=inpoel(J, I)
    !      if(HH_NEW(N1).GT.CC) HH_NEW(N1)=HHMIN_REFIN
    !   end do
    !end do


    !open(12, FILE='remeshing.msh', STATUS='UNKNOWN')
  
    !write(12, '(A)') 'BackgroundMesh V 1.0'
    !write(12, '(A)') 'MESH dimension 2 ElemType Triangle Nnode 3'
    !write(12, '(A)') 'Coordinates'
    !do INOD=1, npoin
    !   write(12, '(I7, 3E14.6)') INOD, X(INOD), Y(INOD)
    !end do
    !write(12, '(A)') 'end Coordinates'
  
    !write(12, '(A)') 'Elements'
    !do IELEM=1, NELEM
    !   write(12, '(5I7)') IELEM, inpoel(1, IELEM), inpoel(2, IELEM), inpoel(3, IELEM)
    !end do
    !write(12, '(A)') 'end Elements'
  
    !write(12, '(A)') 'DesiredSize Nodes'
  
    !do I=1, npoin
    !   write(12, '(I7, E14.6)') I, HH_NEW(I)
    !end do
    !write(12, '(A)') 'end DesiredSize Nodes'
  
    !close(12)

    !CCCC  ----> TIEMPO DE CPU
	call system_clock(m_end)
	print *, "TIEMPO TOTAL: ", real(m_end - m_start)/real(m_rate)
    close(7)
    write(*, *)
    write(*, *)'****-------> FIN DEL CALCULO <-------****'
    write(*, *)
end PROGRAM NSComp2D

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCC---->              RESTART           <----CCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
subroutine RESTART(GAMM) 
    !use DATOS_ENTRADA
	use InputData
	use MeshData
    use MVARIABGEN
    use MVELOCIDADES
    use MVARIABLES
    implicit real(8) (A-H, O-Z)
    real(8) GAMM(npoin)
  
    if (IRESTART.NE.1) THEN
        !CCCC----> SI IOLDSOL .NE. 1 INICIALIZA LAS VARIABLES
        !CCCC----> CON LOS VALORES DEL INFINITO
        !CCCC------------------------------------------------
        RHOAMB=1.225D0 ; PAMB=101325.D0 ; TAMB=288.3D0 ;UAMB=U_inf ; VAMB=0.D0
        do INOD=1, npoin
        
            !CCCC---> VARIABLES CONSERVATIVAS
            U(1, INOD)=RHOAMB
            U(2, INOD)=RHOAMB*UAMB
            U(3, INOD)=RHOAMB*VAMB
            ENERGIA=PAMB/((GAMM(INOD)-1.D0)*RHOAMB)+.5d0*(UAMB**2.D0+VAMB**2.D0)
            U(4, INOD)=ENERGIA*RHOAMB
            VEL_X(INOD)=UAMB
            VEL_Y(INOD)=VAMB
            T(INOD)=TAMB
        end do
    else

        !CCCC----> SI IOLDSOL .EQ. 1 HACE UN RESTART
        !CCCC---------------------------------------
        open(1, FILE=trim(filename)//'.RST', FORM='UNFORMATTED', STATUS='UNKNOWN')
        read(1) ITER_OLD, TIME
        do INOD=1, npoin
            read(1) (U(J, INOD), J=1, 4), T(INOD), GAMM(INOD)
        end do
     
        close(1)

    end if
    return
end subroutine RESTART

subroutine PERIODIC
    !use MGEOMETRIA
	use MeshData
    implicit real(8) (A-H, O-Z)

    do I=1, NMASTER
        N1=IPER_MASTER(I)
        do J=1, NMASTER
            N2=IPER_SLAVE(J)
            XX=DABS(X(N1)-X(N2))
            if(XX.LT.1.D-6)THEN
                IPER_AUX(I)=N2
            end if
        end do
    end do
    return
end subroutine PERIODIC

!CCCC  ----> RUTINA DE INTERPOLACION PARA GAS EN EQUILIBRIO TERMOQUIMICO
subroutine TGAS(E, RHO, P, A, T, GAMM)
    implicit real(8) (A-H, O-Z)

    Y2=DLOG10(RHO/1.292D0)
    Z2=DLOG10(E/78408.4D0)
    if(Y2.GT.-.50D0) GO TO 11
    if(Y2.GT.-4.50D0)GO TO 6
    if(Z2.GT..65D0)  GO TO 1
    GAMM=1.4D0
    SNDSQ=E*.560D0
    GO TO 18
1   if(Z2.GT.1.50D0) GO TO 2
    GAMM=1.46543D0+(.007625D0+.000292D0*Y2)*Y2-(.254500D0+.017244D0*Y2)*Z2 &
        +(.355907D0+.015422D0*Y2-.163235D0*Z2)*Z2*Z2
    GAME=2.304D0*(-.25450D0-.017244D0*Y2+(.711814D0+.030844D0*Y2-.489705D0*Z2)*Z2)
    GAMR=2.304D0*(.007625D0+(-.017244D0+.015422D0*Z2)*Z2+.000584D0*Y2)
    A1=-.000954D0
    A2=.171187D0
    A3=.004567D0
    GO TO 17
2   if(Z2.GT.2.20D0) GO TO 3
    GAS1=2.02636D0+.0584931D0*Y2
    GAS2=.454886D0+.027433D0*Y2
    GAS3=.165265D0+.014275D0*Y2
    GAS4=.136685D0+.010071D0*Y2
    GAS5=.058493D0-.027433D0*Z2
    GAS6=-.014275D0+.010071D0*Z2
    GAS7=DEXP((-10.D0+3.D0*Y2)*(Z2-.023D0*Y2-2.025D0))
    DERE=-30.0D0
    DERR=0.285D0
    A1=.008737D0
    A2=.184842D0
    A3=-.302441D0
    GO TO 15
3   if(Z2.GT.3.05D0) GO TO 4
    GAS1=1.60804D0+.034791D0*Y2
    GAS2=.188906D0+.010927D0*Y2
    GAS3=.124117D0+.007277D0*Y2
    GAS4=.069839D0+.003985D0*Y2
    GAS5=.034791D0-.010927D0*Z2
    GAS6=-.07277D0+.03985D0*Z2
    GAS7=DEXP(-30.0D0*(Z2+.007D0*Y2-2.691D0))
    DERE=-30.0D0
    DERR=0.21D0
    A1=.017884D0
    A2=.153672D0
    A3=-.930224D0
    GO TO 15
4   if(Z2.GT.3.38D0) GO TO 5
    GAS1=1.25672D0+.007073D0*Y2
    GAS2=.039228D0-.000491D0*Y2
    GAS3=-.721798D0-.073753D0*Y2
    GAS4=-.198942D0-.021539D0*Y2
    GAS5=.007073D0+.000491D0*Z2
    GAS6=.073753D0-.021539D0*Z2
    GAS7=DEXP(0.425D0*Y2-50.0D0*Z2+166.7D0)
    DERE=-50.D0
    DERR=0.325D0
    A1=.002379D0
    A2=.217959D0
    A3=.005943D0
    GO TO 15
5   GAMM=-84.0327D0+(-.331761D0+.001153D0*Y2)*Y2+(72.2066D0+.491914D0*Y2)*Z2 &
        +(-20.3559D0-.070617D0*Y2+1.90979D0*Z2)*Z2*Z2
    GAME=2.304D0*( 72.2066D0+.491914D0*Y2+(-40.7118D0-.141234D0*Y2+5.72937D0*Z2)*Z2)
    GAMR=2.304D0*(-.831761D0+.002306D0*Y2+(.491914D0-.070617D0*Z2)*Z2)
    A1=.006572D0
    A2=.183396D0
    A3=-.135960D0
    GO TO 17
6   if(Z2.GT..65D0) GO TO 7
    GAMM=1.4D0
    SNDSQ=E*.560D0
    GO TO 18
7   if(Z2.GT.1.54D0) GO TO 8
    GAS1=1.44813D0+.001292D0*Y2
    GAS2=.073510D0+.001948D0*Y2
    GAS3=-.054745D0+.013705D0*Y2
    GAS4=-.055473D0+.021874D0*Y2
    GAS5=.001292D0-.001948D0*Z2
    GAS6=-.013705D0+.021874D0*Z2
    GAS7=DEXP(-10.0D0*(Z2-1.42D0))
    DERE=-1.D0
    DERR=0.D0
    A1=-.001973D0
    A2=.185233D0
    A3=-.059952D0
    GO TO 15
8   if(Z2.GT.2.22D0) GO TO 9
    GAS1=1.73158D0+.003902D0*Y2
    GAS2=.272846D0-.006237D0*Y2
    GAS3=-.041419D0-.037475D0*Y2
    GAS4=.016984D0-.018038D0*Y2
    GAS5=.003902D0+.006237D0*Z2
    GAS6=.037475D0-.018038D0*Z2
    GAS7=DEXP((-10.0D0+3.0D0*Y2)*(Z2-.023D0*Y2-2.025D0))
    DERE=3.D0*Y2-10.D0
    DERR=3.D0*Z2+12.15D0*Y2-20.325D0
    A1=-.013027D0
    A2=.07427D0
    A3=.012889D0
    GO TO 15
9   if(Z2.GT.2.90D0) GO TO 10
    GAS1=1.59350D0+.075324D0*Y2
    GAS2=.176186D0+.026072D0*Y2
    GAS3=.200838D0+.058536D0*Y2
    GAS6=.099687D0+.025287D0*Y2 !VIENE POR ACA
    GAS5=.075324D0-.026072D0*Z2
    GAS6=-.058536D0+.025287D0*Z2
    GAS7=DEXP((-10.D0+5.D0*Y2)*(Z2-2.7D0))
    DERE=5.D0*Y2-10.D0
    DERR=5.D0*Z2-13.5D0
    A1=.004342D0
    A2=.212192D0
    A3=-.001293D0
    GO TO 15
10  GAS1=1.12688D0-.025957D0*Y2
    GAS2=-.013602D0-.013772D0*Y2
    GAS3=.127737D0+.087942D0*Y2
    GAS4=.043104D0+.023547D0*Y2
    GAS5=-.025957D0+.013772D0*Z2
    GAS6=-.087942D0+.023547D0*Z2
    !C      GAS7=EXP(-20.0D0*Z2+(4.0D0*Z2-13.2D0)*Y2+66.D0)
    GAS7=DEXP((-20.D0+4.D0*Y2)*(Z2-3.3D0))
    DERE=-20.D0+4.D0*Y2
    DERR=4.D0*Z2-13.2D0
    A1=.006348D0
    A2=.209716D0
    A3=-.006001D0
    GO TO 15
11  if(Z2.GT..65D0) GO TO 12
    GAMM=1.4D0
    SNDSQ=E*.560D0
    GO TO 18
12  if(Z2.GT.1.68D0) GO TO 13
    GAS1=1.45510D0-.000102D0*Y2
    GAS2=.081537D0-.000166D0*Y2
    GAS3=-.128647D0+.049454D0*Y2
    GAS4=-.101036D0+.033518D0*Y2
    GAS5=-.000102D0+.000166D0*Z2
    GAS6=-.049454D0+.033518D0*Z2
    GAS7=DEXP(-15.D0*(Z2-1.420D0))
    DERE=-15.D0
    DERR=0.D0
    A1=.00045D0
    A2=.203892D0
    A3=.101797D0
    GO TO 15
13  if(Z2.GT.2.46D0) GO TO 14
    GAS1=1.59608D0-.042426D0*Y2
    GAS2=.192840D0-.029353D0*Y2
    GAS3=.019430D0-.005954D0*Y2
    GAS4=.026097D0-.006164D0*Y2
    GAS5=-.042426D0+.029353D0*Z2
    GAS6=.005954D0-.006164D0*Z2
    GAS7=DEXP(-15.D0*(Z2-2.050D0))
    DERE=-15.D0
    DERR=0.D0
    A1=-.006609D0
    A2=.127637D0
    A3=.297037D0
    GO TO 15
14  GAS1=1.54363D0-.049071D0*Y2
    GAS2=.153562D0-.029209D0*Y2
    GAS3=.324907D0+.077599D0*Y2
    GAS4=.142408D0+.022071D0*Y2
    GAS5=-.049071D0+.029209D0*Z2
    GAS6=-.077599D0+.022071D0*Z2
    GAS7=DEXP(-10.0D0*(Z2-2.708D0))
    DERE=-10.D0
    DERR=0.D0
    A1=-.000081D0
    A2=.226601D0
    A3=.170922D0
15  GAS10=1.D0/(1.+GAS7)
16  GAS8=GAS3-GAS4*Z2
    GAS8=GAS3-GAS4*Z2
    GAS9=GAS8*GAS7*GAS10*GAS10
    GAMM=GAS1-GAS2*Z2-GAS8*GAS10
    GAME=2.304D0*(-GAS2+GAS4*GAS10+GAS9*DERE)
    GAMR=2.304D0*(GAS5+GAS6*GAS10+GAS9*DERR)
17  SNDSQ=E*(A1+(GAMM-1.D0)*(GAMM+A2*GAME)+A3*GAMR)
18  A=DSQRT(SNDSQ)
    P=RHO*E*(GAMM-1.D0)
    X2=DLOG10(P/1.0134D5)
    Y2=Y2+.0231264D0
    Z3=X2-Y2
    if(Y2.GT.-.50D0) GO TO 29
    if(Y2.GT.-4.50D0)GO TO 24
    if(Z3.GT..30D0)  GO TO 19
    T=P/(287.D0*RHO)
    return
19  if(Z3.GT.1.0D0) GO TO 20
    T=10.D0**(.2718D0+.00074D0*Y2+(.990136D0-.004947D0*Y2)*Z3+(.990717D0 &
        +.175194D0*Y2-(.982407D0+.159233D0*Y2)*Z3)/(1.D0+DEXP(-20.D0*(Z3-.88D0))))
    GO TO 32
20  if(Z3.GT.1.35D0) GO TO 21
    T=10.D0**(1.39925D0+.167780D0*Y2+(-.143168D0-.159234D0*Y2)*Z3+(-.027614D0 &
        -.090761D0*Y2+(.307036D0+.121621D0*Y2)*Z3)/(1.D0+DEXP(-20.D0*(Z3-1.17D0))))
    GO TO 32
21  if(Z3.GT.1.79D0) GO TO 22
    T=10.D0**(1.11401D0+.002221D0*Y2+(.351875D0+.017246D0*Y2)*Z3+(-1.15099D0 &
        -.173555D0*Y2+(.673342D0-.088399D0*Y2)*Z3)/(1.D0+DEXP(-20.D0*(Z3-1.56D0))))
    GO TO 32

22  if(Z3.GT.2.47D0) GO TO 23
    T=10.D0**(1.01722D0-.017918D0*Y2+(.473523D0+.025456D0*Y2)*Z3+(-2.17978D0 &
        -.334716D0*Y2+(.898619D0+.127386D0*Y2)*Z3)/(1.D0+DEXP(-20.D0*(Z3-2.22D0))))
    GO TO 32
23  T=10.D0**(-45.0871D0-9.00504D0*Y2+(35.8685D0+6.79222D0*Y2)*Z3-(6.77699D0 &
        +1.2737D0*Y2)*Z3*Z3+(-.064705D0+.025325D0*Z3)*Y2*Y2)
    GO TO 32
24  if(Z3.GT..48D0) GO TO 25
    T=P/(287.D0*RHO)
    return
25  if(Z3.GT..9165D0) GO TO 26
    T=10.D0**(.284312D0+.987912D0*Z3+.001644D0*Y2)
    GO TO 32
26  if(Z3.GT.1.478D0) GO TO 27
    T=10.D0**(.502071D0-.01299D0*Y2+(.774818D0+.025397D0*Y2)*Z3+(.009912D0 &
        -.150527D0*Y2+(-.000385D0+.105734D0*Y2)*Z3)/(1.D0+DEXP(-15.D0*(Z3-1.28D0))))
    GO TO 32
27  if(Z3.GT.2.176D0) GO TO 28
    T=10.D0**(1.02294D0+.021535D0*Y2+(.427212D0+.0069D0*Y2)*Z3+(-.427823D0 &
        -.211991D0*Y2+(.257096D0+.101192D0*Y2)*Z3)/(1.D0+DEXP(-12.D0*(Z3-1.778D0))))
    GO TO 32
28  T=10.D0**(1.47540D0+.12962D0*Y2+(.254154D0-.046411D0*Y2)*Z3+(-.221229D0 &
        -.057077D0*Y2+(.158116D0+.03043D0*Y2)*Z3)/(1.D0+DEXP(5.D0*Y2*(Z3-2.40D0))))
    GO TO 32
29  if(Z3.GT..48D0) GO TO 30
    T=P/(287.D0*RHO)
    return
30  if(Z3.GT.1.07D0) GO TO 31
    T=10.D0**(.279268D0+.992172D0*Z3)
    GO TO 32
31  T=10.D0**(.233261D0-.056383D0*Y2+(1.19783D0+.063121D0*Y2-.165985D0*Z3)*Z3+ &
        (-.814535D0+.099233D0*Y2+(.602385D0-.067428D0*Y2-.093991D0*Z3)*Z3)/(1.D0+DEXP(( &
        5.D0*Y2-20.D0)*(Z3-1.78D0))))
32  T=T*151.777778D0
    return
end

!CCCC----------------------------------------------------!CCCC
!CCCC  ----> Imprime resultados en formato FLAVIA <----  !CCCC
!CCCC----------------------------------------------------!CCCC 
subroutine PRINTFLAVIA(GAMM, RHO, VEL_X, VEL_Y, P, T, E, RMACH, XPOS, YPOS &
    , npoin, ITER)
     
    !use MPRINTRES
	use InputData
    implicit real(8) (A-H, O-Z)
  
    INTEGER npoin
  
    real(8) RHO(npoin), VEL_X(npoin), VEL_Y(npoin)
    real(8) P(npoin), T(npoin), E(npoin), GAMM(npoin), RMACH(npoin)
    real(8) XPOS(npoin), YPOS(npoin), U(4, npoin)
  
    !CCCC----> IMPRESION DE RESULTADOS
    !CCCC-----------------------------
    if (MOVIE.EQ.0)THEN
        open(2, FILE=trim(filename)//'.flavia.res', STATUS='UNKNOWN')
    end if
  
    !CCCC----> ESCRITURA DE VELOCIDADES
    !CCCC------------------------------
    if(VEL2CHAR.EQ.'.si.')THEN
        write(2, '(A15, 5(I6, 2X))') 'VELOCITY', 2, ITER, 2, 1, 1
        write(2, '(A)') 'VEL_X'
        write(2, '(A)') 'VEL_Y'
  
        do INOD=1, npoin
            write(2, '(I6, 3E13.4)') INOD, VEL_X(INOD), VEL_Y(INOD)
        end do
    end if
  
    !CCCC----> ESCRITURA DE POSICIONES
    !CCCC------------------------------
    if(POSCHAR.EQ.'.si.')THEN
        write(2, '(A15, 5(I6, 2X))') 'POSITION', 2, ITER, 2, 1, 1
        write(2, '(A)') 'X'
        write(2, '(A)') 'Y'
  
        do INOD=1, npoin
            write(2, '(I6, 3E16.6)') INOD, XPOS(INOD), YPOS(INOD)
        end do
    end if
  
    !CCCC----> ESCRITURA DE LAS DENSIDADES
    !CCCC---------------------------------
    if(RHOCHAR.EQ.'.si.')THEN
        write(2, '(A15, 5(I6, 2X))') 'DENSITY', 2, ITER, 1, 1, 1
        write(2, '(A)') 'DENSITY'
  
        do INOD=1, npoin
            write(2, '(I6, 1E13.4)') INOD, RHO(INOD)
        end do
    end if
  
    !CCCC----> ESCRITURA DE LAS PRESIONES
    !CCCC--------------------------------
    if(PRESCHAR.EQ.'.si.')THEN
        write(2, '(A15, 5(I6, 2X))') 'PRESSURE', 2, ITER, 1, 1, 1
        write(2, '(A)') 'PRESSURE'
  
        do INOD=1, npoin
            write(2, '(I6, 1E16.3)') INOD, P(INOD)
        end do
    end if
  
    !CCCC  -----> ESCRITURA DE LAS TEMPERATURAS
    !CCCC  ------------------------------------
    if(TEMPCHAR.EQ.'.si.')THEN
        write(2, '(A15, 5(I6, 2X))') 'TEMPERATURE', 2, ITER, 1, 1, 1
        write(2, '(A)') 'TEMPERATURE'
  
        do INOD=1, npoin
            write(2, '(I6, 1E13.3)') INOD, T(INOD)
        end do
    end if
  
    !CCCC----> ESCRITURA DEL NUMERO DE MACH
    !CCCC-----------------------------------
    if(MACHCHAR.EQ.'.si.')THEN
        write(2, '(A15, 5(I6, 2X))') 'Mach_Number', 2, ITER, 1, 1, 1
        write(2, '(A)') 'Mach_Number'
        !write(2, '(A)') 'Mach_Number_Tgas'
  
        do INOD=1, npoin
            VEL=DSQRT(VEL_X(INOD)*VEL_X(INOD)+VEL_Y(INOD)*VEL_Y(INOD))
            VC=DSQRT(GAMM(INOD)*FR*T(INOD))
            write(2, '(I6, 2F11.2)') INOD, VEL/VC!, RMACH(INOD)
     
        end do
    end if
  
    !CCCC----> ESCRITURA DE LA ENERGIA INTERNA
    !CCCC-----------------------------------
    if(ENERCHAR.EQ.'.si.')THEN
        write(2, '(A15, 5(I6, 2X))') 'Internal_Energy', 2, ITER, 1, 1, 1
        write(2, '(A)') 'Internal_Energy'
  
        do INOD=1, npoin
            write(2, '(I6, 1E16.8)') INOD, E(INOD)
        end do
    end if
  
    !CCCC----> ESCRITURA DE GAMA
    !CCCC-----------------------------------
    !write(2, '(A15, 5I6)') 'GAMA', 2, ITER, 1, 1, 1
    !write(2, '(A)') 'GAMA'
  
    !do INOD=1, npoin
    !   write(2, '(I6, 1E16.8)') INOD, GAMM(INOD)
    !end do

    if (MOVIE.EQ.0)THEN
        close(2)
    end if
  
    return
end subroutine PRINTFLAVIA

!CCCC-----------------------------!CCCC
!CCCC  CALCULO DE LAS FUERZAS EN  !CCCC
!CCCC  CADA UNO DE LOS SETS       !CCCC
!CCCC-----------------------------!CCCC
subroutine FORCES(NSETS, NSET_NUMB, IELEM_SETS, npoin &
    , X, Y, P &
    , FX, FY, RM, X_ROT, Y_ROT)
    !use MFUERZAS
	use MeshData, only: ISET, IPER_MASTER, IPER_SLAVE, IPER_AUX
    implicit real(8) (A-H, O-Z)
  
    INTEGER IELEM_SETS(10)
  
    real(8) X(npoin), Y(npoin)
    real(8) P(npoin)
    real(8) FX(10), FY(10), RM(10)
  
    FX(1:NSET_NUMB)=0.D0
    FY(1:NSET_NUMB)=0.D0
    RM(1:NSET_NUMB)=0.D0
  
    do ISET_NUMB=1, NSET_NUMB
        do II=1, IELEM_SETS(ISET_NUMB)
            N1=ISET(1, ISET_NUMB, II)
            N2=ISET(2, ISET_NUMB, II)
        
            D_PRESS=(P(N1)+P(N2))/2.D0
            RLX=X(N1)-X(N2)
            RLY=Y(N2)-Y(N1)
        
            DFX=D_PRESS*RLY
            DFY=D_PRESS*RLX
        
            FX(ISET_NUMB)=FX(ISET_NUMB)+DFX
            FY(ISET_NUMB)=FY(ISET_NUMB)+DFY
        
            XC=(X(N1)+X(N2))/2.D0
            YC=(Y(N2)+Y(N1))/2.D0
        
            RM(ISET_NUMB)=RM(ISET_NUMB)+DFY*(XC-X_ROT)-DFX*(YC-Y_ROT)
        
        end do
    end do
    return
end subroutine FORCES
      
!CCCC-----------------------------------!CCCC
!CCCC  CALCULO DE LAS FUERZAS VISCOSAS  !CCCC
!CCCC-----------------------------------!CCCC
subroutine FORCE_VISC(NELEM, npoin, IFF &
    , IELEM_SETS, ISET, inpoel, NSET_NUMB, U_inf, V_inf, RHO_inf, T_inf &
    , X, Y, P, T, VEL_X, VEL_Y, DNX, DNY, RMU, RHO &
    , F_VX, F_VY)
      

    implicit real(8) (A-H, O-Z)

    INTEGER IELEM_SETS(10), ISET(3, 10, IFF), inpoel(3, NELEM)

    real(8) X(npoin), Y(npoin), P(npoin), VEL_X(npoin), VEL_Y(npoin), RHO(npoin)
    real(8) DNX(3, NELEM), DNY(3, NELEM), T(npoin)
    real(8) F_VX(10), F_VY(10)
      
    F_VX=0.D0 ; F_VY=0.D0
    open(1, FILE='SKIN.DAT', STATUS='UNKNOWN')
    do ISET_NUMB=1, NSET_NUMB

        do II=1, IELEM_SETS(ISET_NUMB)
            
            NN1=ISET(1, ISET_NUMB, II)
            NN2=ISET(2, ISET_NUMB, II)
            IELEM=ISET(3, ISET_NUMB, II)

            TEMP=(T(NN1)+T(NN2))/2.D0
            FMU= 1.716d-5*162.6/(TEMP-110.55)*(TEMP/273.15)**.75D0 !SUTHERLAND
            RLY=-(X(NN2)-X(NN1))
            RLX=Y(NN2)-Y(NN1)
            RMOD=DSQRT(RLX*RLX+RLY*RLY)

            DX=-RLY
            DY=RLX

            RLX=RLX/RMOD
            RLY=RLY/RMOD
            
            DUX=0.D0 ; DUY=0.D0 ; DVX=0.D0 ; DVY=0.D0
            PRESS=0.D0
            do JJ=1, 3
                NN=inpoel(JJ, IELEM)
                !CCCC  ----> DERIVADA DE VEL_X
                DUX=DUX+DNX(JJ, IELEM)*VEL_X(NN)
                DUY=DUY+DNY(JJ, IELEM)*VEL_X(NN)
                !CCCC  ----> DERIVADA DE VEL_Y
                DVX=DVX+DNX(JJ, IELEM)*VEL_Y(NN)
                DVY=DVY+DNY(JJ, IELEM)*VEL_Y(NN)
                PRESS=PRESS+P(NN)
            end do

            PRESS=PRESS/3.D0
            !CCCC  ----> TENSOR 2D
            TXX=-FMU*2.D0/3.D0*(DUX+DVY)+2.D0*FMU*DUX
            TXY=FMU*(DUY+DVX)
            
            TYX=TXY
            TYY=-FMU*2.D0/3.D0*(DUX+DVY)+2.D0*FMU*DVY

            !CCCC  ----> CALCULO PARA SACAR EL SKIN FRICTION
            TTX=(-PRESS+TXX)*RLX+TXY*RLY
            TTY=TYX*RLX+(-PRESS+TYY)*RLY
            TMOD=-TTX*RLY+TTY*RLX
            UU=U_inf*U_inf+V_inf*V_inf
            SKIN=TMOD/(.5d0*RHO_inf*UU)

            F_VX(ISET_NUMB)=F_VX(ISET_NUMB)+((-PRESS+TXX)*RLX+TXY*RLY)*DX
            F_VY(ISET_NUMB)=F_VY(ISET_NUMB)+(TYX*RLX+(-PRESS+TYY)*RLY)*DY
            write(1, *)SKIN, (X(NN2)+X(NN1))/2.D0
        end do
    end do
    close(1)
    return
end

!CCCC-----------------------------------------------!CCCC
!CCCC  ----> Imprime resultados para RESTART <----  !CCCC
!CCCC-----------------------------------------------!CCCC 
subroutine PRINTREST(ITER, npoin, U, T, GAMM, TIME, filename)

    implicit real(8) (A-H, O-Z)
      
    real(8) U(4, npoin), T(npoin), GAMM(npoin)

    character(80) filename 
      
    open(1, FILE=trim(filename)//'.RST', FORM='UNFORMATTED', STATUS='UNKNOWN')

    write(1) ITER, TIME
      
    do INOD=1, npoin
        write(1) (U(J, INOD), J=1, 4), T(INOD), GAMM(INOD)
    end do
      
    close(1)
      
    return
end

!CCCC---------------------------CCCC
!CCCC----> PRESSURE SWITCH <----CCCC
!CCCC---------------------------CCCC
subroutine NEWPRES(npoin, NPOS, NN1, NN2, P, PS)

    implicit real(8) (A-H, O-Z)
    INTEGER(4) NN1(NPOS), NN2(NPOS)
    real(8)S1(npoin), S2(npoin), PS(npoin), P(npoin)

    S1=0.D0
    S2=0.D0
  
    do IPOS=1, NPOS
        N1=NN1(IPOS)
        N2=NN2(IPOS)
        AUX=P(N1)-P(N2)
        S1(N1)=S1(N1)+AUX
        S2(N1)=S2(N1)+DABS(AUX)
    end do
  
    do I=1, npoin
        PS(I)=DABS(S1(I))/(S2(I))
        if(S2(I).LT.5.D-2)PS(I)=0.D0
        if(PS(I).LT..2D0)PS(I)=0.D0 !MODIFICAR SI ES NECESARIO
    end do
    return
end subroutine NEWPRES

!CCCC-----> CALCULA EL TAMA\D1O OPTIMO PARA REFINAMIENTO ADAPTATIVO
!CCCC------------------------------------------------------------
subroutine NEW_SIZE(npoin, NELEM, inpoel, DNX, DNY, AREA, M, HH, U &
    , AUX)
    !use DATOS_REFINAMIENTO
	use InputData
    implicit real*8 (A-H, O-Z)
  
    INTEGER inpoel(3, NELEM)
    real(8) DNX(3, NELEM), DNY(3, NELEM), AREA(NELEM), HH(NELEM), U(4, npoin)
    real(8) HH_NEW(NELEM)
    real(8) SXX(npoin), SXY(npoin), SYY(npoin), AUX(npoin), M(npoin)
  
    SXX=0.D0;  SYY=0.D0;  SXY=0.D0;   AUX=0.D0

    do IELEM=1, NELEM
        N1=inpoel(1, IELEM) ; N2=inpoel(2, IELEM) ; N3=inpoel(3, IELEM)
     
        RNX1=DNX(1, IELEM) ; RNX2=DNX(2, IELEM) ; RNX3=DNX(3, IELEM)
        RNY1=DNY(1, IELEM) ; RNY2=DNY(2, IELEM) ; RNY3=DNY(3, IELEM)
     
        !CCCC  ----> DERIVADA DE RHO.VEL_X
        DUX= RNX1*U(2, N1)+RNX2*U(2, N2)+RNX3*U(2, N3)
        DUY= RNY1*U(2, N1)+RNY2*U(2, N2)+RNY3*U(2, N3)
        !CCCC  ----> DERIVADA DE RHO.VEL_Y
        DVX= RNX1*U(3, N1)+RNX2*U(3, N2)+RNX3*U(3, N3)
        DVY= RNY1*U(3, N1)+RNY2*U(3, N2)+RNY3*U(3, N3)
     
        SIGMXX=2*DUX;   SIGMYY=2*DVY;    SIGMXY=DUY+DVX
     
        do I=1, 3
            N1=inpoel(I, IELEM)
            SXX(N1)=SXX(N1)+SIGMXX*AREA(IELEM)
            SYY(N1)=SYY(N1)+SIGMYY*AREA(IELEM)
            SXY(N1)=SXY(N1)+SIGMXY*AREA(IELEM)
            AUX(N1)=AUX(N1)+AREA(IELEM)
        end do
    end do
  
    do INOD=1, npoin
        SXX(INOD)=SXX(INOD)/AUX(INOD)
        SYY(INOD)=SYY(INOD)/AUX(INOD)
        SXY(INOD)=SXY(INOD)/AUX(INOD)
    end do
  
  
    ESIG=0.D0;   SIG=0.D0
    do IELEM=1, NELEM
        N1=inpoel(1, IELEM) ; N2=inpoel(2, IELEM) ; N3=inpoel(3, IELEM)
     
        RNX1=DNX(1, IELEM) ; RNX2=DNX(2, IELEM) ; RNX3=DNX(3, IELEM)
        RNY1=DNY(1, IELEM) ; RNY2=DNY(2, IELEM) ; RNY3=DNY(3, IELEM)
     
        !CCCC  ----> DERIVADA DE RHO.VEL_X
        DUX= RNX1*U(2, N1)+RNX2*U(2, N2)+RNX3*U(2, N3)
        DUY= RNY1*U(2, N1)+RNY2*U(2, N2)+RNY3*U(2, N3)
        !CCCC  ----> DERIVADA DE RHO.VEL_Y
        DVX= RNX1*U(3, N1)+RNX2*U(3, N2)+RNX3*U(3, N3)
        DVY= RNY1*U(3, N1)+RNY2*U(3, N2)+RNY3*U(3, N3)
     
        SIGMXX=2*DUX;   SIGMYY=2*DVY;    SIGMXY=DUY+DVX
     
        ESIGXX=SIGMXX-(SXX(N1)+SXX(N2)+SXX(N3))/3.D0
        SIGXX=(SXX(N1)+SXX(N2)+SXX(N3))/3.D0
     
        ESIGYY=SIGMYY-(SYY(N1)+SYY(N2)+SYY(N3))/3.D0
        SIGYY=(SYY(N1)+SYY(N2)+SYY(N3))/3.D0
     
        ESIGXY=SIGMXY-(SXY(N1)+SXY(N2)+SXY(N3))/3.D0
        SIGXY=(SXY(N1)+SXY(N2)+SXY(N3))/3.D0
         
        ESIG=ESIG+(ESIGXX*ESIGXX+ESIGYY*ESIGYY+ESIGXY*ESIGXY)*AREA(IELEM)
        SIG=SIG+(SIGXX*SIGXX+SIGYY*SIGYY+SIGXY*SIGXY)*AREA(IELEM)
    end do
  
    ERR_AVERAGE=ETA_REFIN*DSQRT((SIG+ESIG)/DFLOAT(NELEM))

    do IELEM=1, NELEM
     
        N1=inpoel(1, IELEM) ; N2=inpoel(2, IELEM) ; N3=inpoel(3, IELEM)
     
        RNX1=DNX(1, IELEM) ; RNX2=DNX(2, IELEM) ; RNX3=DNX(3, IELEM)
        RNY1=DNY(1, IELEM) ; RNY2=DNY(2, IELEM) ; RNY3=DNY(3, IELEM)
     
        !CCCC  ----> DERIVADA DE RHO.VEL_X
        DUX= RNX1*U(2, N1)+RNX2*U(2, N2)+RNX3*U(2, N3)
        DUY= RNY1*U(2, N1)+RNY2*U(2, N2)+RNY3*U(2, N3)
        !CCCC  ----> DERIVADA DE RHO.VEL_Y
        DVX= RNX1*U(3, N1)+RNX2*U(3, N2)+RNX3*U(3, N3)
        DVY= RNY1*U(3, N1)+RNY2*U(3, N2)+RNY3*U(3, N3)
     
        SIGMXX=2*DUX;   SIGMYY=2*DVY;    SIGMXY=DUY+DVX
     
        ESIGXX=SIGMXX-(SXX(N1)+SXX(N2)+SXX(N3))/3.D0
     
        ESIGYY=SIGMYY-(SYY(N1)+SYY(N2)+SYY(N3))/3.D0
         
        ESIGXY=SIGMXY-(SXY(N1)+SXY(N2)+SXY(N3))/3.D0
         
        ESIG=(ESIGXX*ESIGXX+ESIGYY*ESIGYY+ESIGXY*ESIGXY)*AREA(IELEM)
     
        EPS_I=ESIG/ERR_AVERAGE
        HH_NEW(IELEM)=HH(IELEM)/EPS_I
     
    end do

    AUX=0.D0
    do I=1, NELEM
        AR=AREA(I)
        do J=1, 3
            N1=inpoel(J, I)
            AUX(N1)=AUX(N1)+HH_NEW(I)*AR
        end do
    end do
    do I=1, npoin
        AUX(I)=AUX(I)/M(I)
        if(AUX(I).GT.HHMAX_REFIN) AUX(I)=HHMAX_REFIN
        if(AUX(I).LT.HHMIN_REFIN) AUX(I)=HHMIN_REFIN
    end do
  
    return
end subroutine NEW_SIZE

!CCCC-----> ESTIMADOR DE ERROR
!CCCC-------------------------------
! INPUT: VAR= VARIABLE DE ENTRADA A LA CUAL CALCULALOS EL ERROR
! OUTPUT: AUX= ERROR ESTIMADO SUAVIZADO A LOS NODOS
subroutine ESTIMADOR_ERR(npoin, NELEM, inpoel, DNX, DNY, HH, M, AREA, BETA, VAR, AUX)
    !use DATOS_REFINAMIENTO
	use InputData
    implicit real(8) (A-H, O-Z)
    INTEGER inpoel(3, NELEM)
    real(8)DNX(3, NELEM), DNY(3, NELEM), HH(NELEM), VAR(npoin)
    real(8) TITA(NELEM), AUX(npoin), AREA(NELEM), M(npoin)
    TITA=0.D0 ; TITA_PROM=0.D0
    do IELEM=1, NELEM
        N1=inpoel(1, IELEM) ; N2=inpoel(2, IELEM) ; N3=inpoel(3, IELEM)
     
        RNX1=DNX(1, IELEM) ; RNX2=DNX(2, IELEM) ; RNX3=DNX(3, IELEM)
        RNY1=DNY(1, IELEM) ; RNY2=DNY(2, IELEM) ; RNY3=DNY(3, IELEM)
     
        !CCCC  ----> DERIVADA DE LA VARIABLE
        DUX= RNX1*VAR(N1)+RNX2*VAR(N2)+RNX3*VAR(N3)
        DUY= RNY1*VAR(N1)+RNY2*VAR(N2)+RNY3*VAR(N3)

        TITA(IELEM)=DSQRT((DABS(DUX)+DABS(DUY))*.5D0)*HH(IELEM)
        TITA_PROM=TITA_PROM+TITA(IELEM)
    end do
    TITA_PROM=TITA_PROM/NELEM

    TITA_SD=0.D0
    do I=1, NELEM
        TITA_SD=TITA_SD+(TITA(I)-TITA_PROM)**2
    end do
    TITA_SD=DSQRT(TITA_SD/NELEM)

    !CCCC----> SUAVIZO EN A LOS NODOS
    AUX=0.D0
    do I=1, NELEM
        AR=AREA(I)
        do J=1, 3
            N1=inpoel(J, I)
            AUX(N1)=AUX(N1)+TITA(I)*AR
        end do
    end do
    do I=1, npoin
        AUX(I)=AUX(I)/M(I)
    end do

    !CCCC----> APLICO FILTRO PASA ALTO
    !BETA=3.5D0 ! MODIFICO A OJO
    COMPARADOR=TITA_PROM+BETA*TITA_SD
    do I=1, npoin
        if(AUX(I).GT.COMPARADOR)THEN
            AUX(I)=HHMIN_REFIN
        else
            AUX(I)=HHMAX_REFIN
        end if
    end do

    return
end subroutine ESTIMADOR_ERR

!CCCC---- SMOOTHING DE LOS ELEMENTOS
!CCCC-------------------------------
subroutine SMOOTH_MESH(npoin, NELEM, inpoel, X, Y &
    , SMOOTH_FIX, SMOOTH_SIM)
  
    implicit real(8) (A-H, O-Z)
  
    INTEGER inpoel(3, NELEM), SMOOTH_SIM(2, npoin)
    INTEGER NELINOD(npoin), NE(20, npoin)
  
    real(8) X(npoin), Y(npoin)
    real(8) MU_0, MU_X, MU_Y, MU_Z, RMU1, RMUMIN, GG, GGI, GAM1, GAMMIN
    real(8) MU(100), GX(100), GY(100), XX(4), YY(4)
  
    LOGICAL SMOOTH_FIX(npoin)
  
    write(*, *) 'suavizado inicial  1'
  
    do INOD=1, npoin
        NELINOD(INOD)=0
    end do
  
    do IELEM=1, NELEM
        do IN=1, 3
            INOD=inpoel(IN, IELEM)
            NELINOD(INOD)=NELINOD(INOD)+1
            NE(NELINOD(INOD), INOD)=IELEM
        end do
    end do
  
    DELTA=1.D-5
    !CCCC   COMIENZO DE ITERACIONES DE OPTIMIZACION
  
    do IT=1, 1000
     
     
        RMM=1.D10
        ERR=0.D0; xmax=0
        NUMEL=0
        do INOD=1, npoin
        
            if (SMOOTH_FIX(INOD)) THEN
           
                RMUMIN=1.D10
                do IEL=1, NELINOD(INOD)
                    IELEM=NE(IEL, INOD)
              
                    do IN=1, 3
                        N1=inpoel(IN, IELEM)
                        XX(IN)=X(N1)
                        YY(IN)=Y(N1)
                        if (N1.EQ.INOD) IN_P=IN
                    end do
              
                    MU_0=RMU1(XX, YY)
              
                    if (RMM.GT.MU_0) THEN
                        RMM=MU_0
                        IELMIN=IELEM
                    end if
              
                    MU(IEL)=MU_0
                    if (MU_0.LT.RMUMIN) THEN
                        RMUMIN=MU_0
                        IELM=IEL
                    end if
              
                    XX(IN_P)=X(INOD)+DELTA
                    MU_X=RMU1(XX, YY)
                    XX(IN_P)=X(INOD)
                    GX(IEL)=(MU_X-MU_0)/DELTA
              
                    YY(IN_P)=Y(INOD)+DELTA
                    MU_Y=RMU1(XX, YY)
                    YY(IN_P)=Y(INOD)
                    GY(IEL)=(MU_Y-MU_0)/DELTA
                end do
           
                if (RMUMIN.GT.0.8D0) CYCLE
                NUMEL=NUMEL+1
           
                GMX=GX(IELM)
                GMY=GY(IELM)
           
                GAM1=0.D0
                GAMMIN=1.D10
                do IEL=1, NELINOD(INOD)
                    GG=GMX*GMX+GMY*GMY
                    if (IEL.NE.IELM) THEN
                        GGI=GX(IEL)*GMX+GY(IEL)*GMY
                        if (GGI.LT.0) THEN
                            GAM1=(MU(IEL)-RMUMIN)/(GG-GGI)
                            if (GAMMIN.GT.GAM1.AND.GAM1.NE.0.) GAMMIN=GAM1
                        end if
                    end if
                end do
           
                GAMMIN=GAMMIN/8
                if (GAMMIN.GT.100) GAMMIN=1.D-6
           
                X1=X(INOD)+GMX*GAMMIN*SMOOTH_SIM(1, INOD)
                Y1=Y(INOD)+GMY*GAMMIN*SMOOTH_SIM(2, INOD)
           
                ERR=ERR+(X(INOD)-X1)**2+(Y(INOD)-Y1)**2
           
                X(INOD)=X1
                Y(INOD)=Y1
           
            end if
        end do
     
        if (NUMEL.EQ.0) EXIT
     
    end do
  
    write(*, *) 'ERROR FINAL:', ERR, '  PASOS:', IT-1, NUMEL
  
    return
end subroutine SMOOTH_MESH

FUNCTION RMU1(XX, YY)
    implicit real(8) (A-H, O-Z)
    real(8) XX(4), YY(4)
  
    AREA=XX(2)*YY(3)+XX(3)*YY(1)+XX(1)*YY(2) &
        -(XX(2)*YY(1)+XX(3)*YY(2)+XX(1)*YY(3))
  
    RL21=(XX(2)-XX(1))**2.D0+(YY(2)-YY(1))**2.D0
    RL32=(XX(3)-XX(2))**2.D0+(YY(3)-YY(2))**2.D0
    RL13=(XX(1)-XX(3))**2.D0+(YY(1)-YY(3))**2.D0
  
    RL=RL21+RL32+RL13
  
    RMU1=2.D0*DSQRT(3.D0)*AREA/RL
  
    return
end FUNCTION RMU1

!--------------------------------------------------------------------------       
!                 Newmark's Direct Integration Method
!--------------------------------------------------------------------------
! Code written by : - Ing. German Weht                                    
!                     Resercher & PhD Student                              
!                     Departamento de Mecanica Aeronautica                     
!                     Instituto Universitario Aeronautico           
!                     Cordoba, Argentina  
! E-mail : gerweht@gmail.com
!-------------------------------------------------------------------------
! PURPOSE
!        Solve linear structural differential equations 
!               *** Newmark time integration ***
!        The system of ordinary differential equations
!                [M]{ACC}+[C]{VEL}+[K]{DIS}={F}
!
!        
! INPUT
!        [M] :       System Mass              [inpoel, inpoel]
!        [C] :       System Damping           [inpoel, inpoel]
!        [K] :       System Stiffness         [inpoel, inpoel]
!        [F] :       Externally Applied Load  [inpoel]
!        [DISS] :    Initial Position         [inpoel]
!        [VELS] :    Initial Velocity         [inpoel]
!        [ACCS] :    Initial Aceleration      [inpoel]
!        [T] :       Initial Time              
!        [TMAX] :    Finish Time
!        [H] :       Time step 
!        [G] :       Auxiliaries              [inpoel]
!        [G1] :      Auxiliaries              [inpoel]
! OUTPUT
!       [DIS]:       Displacemente            [inpoel]
!       [VEL]:       Velocity                 [inpoel]
!       [ACC]:       Acceleration             [inpoel]
!
! The options include changing the value of the "gamma" and "beta"
! coefficient which appear in the formulation of the method. By default
! these values are set to gamma = 1/2 and beta = 1/4.
!
!
!-------------------------------------------------------------------------
subroutine NEWMARK_METHOD(H)
    use MATRICES
    use MAT2
    implicit real(8)(A-H, O-Z)
    !real(8), ALLOCATABLE:: M(:, :), C(:, :), K(:, :)
    !real(8), ALLOCATABLE:: DIS(:), DISS(:), VEL(:), VELS(:), ACC(:), ACCS(:), R(:), G(:), G1(:), F(:), CTE(:)
    PI=DACOS(-1.D0)
    GR=1.d0!9.81D0 ! gravedad


  

  
    ! **** DEFINICION DE CONSTANTES DEL METODO POR DEFAULT GAMA=1/2 BETA=1/4 ****
    GAMA=.5D0 ; BETA=.25D0

    A1 = 1.D0/(BETA*H*H) ; A2 = 1.D0/(BETA*H) ; A3 = 1.D0/(2.D0*BETA) - 1.D0 ; A4 = (1.D0 - GAMA)*H ; A5 = GAMA*H
    A1D = GAMA/(BETA*H) ; A2D = GAMA/BETA - 1.D0 ; A3D = 0.5D0*H*(GAMA/BETA - 2.D0)
  

  
    !ALLOCATE(M(NDIM, NDIM), C(NDIM, NDIM), K(NDIM, NDIM))
    !ALLOCATE(DIS(NDIM), DISS(NDIM), VEL(NDIM), VELS(NDIM), ACC(NDIM), ACCS(NDIM))
    !ALLOCATE(R(NDIM), G(NDIM), G1(NDIM), F(NDIM), CTE(NDIM))
  

    !***** MATRICES *****
    ! MATRIZ DE MASA
    MASA(1, 1)=51.5d0   ; MASA(1, 2)=-51.5d0*.0429d0
    MASA(2, 1)=-51.5d0*.0429d0   ; MASA(2, 2)=2.275d0

    ! MATRIZ DE AMORTIGUAMIENTO
    C(1, 1)=32.358D0/GR ;  C(1, 2)=0.D0
    C(2, 1)=0.D0   ;  C(2, 2)=5.718D0/GR

    ! MATRIZ DE RIGUIDEZ
    K(1, 1)=50828.463D0/GR   ; K(1, 2)=0.D0
    K(2, 1)=0.D0   ; K(2, 2)=35923.241D0/GR

    ! VECTOR DE FUERZAS EXTERNAS
  
 
    ! CONDICIONES INICIALES
    !DIS(1)=0.05D0; DIS(2)=PI/8.D0
    !VEL=0.D0
    !ACC=0.D0

    ! Matrix with [Keff]=[K]+a1[M]+a1D[C]
    K=K+A1*MASA+A1D*C
  
    !CCCC  -----> TRIANGULAR SUPERIOR DE LA MATRIZ K, SI ES CONSTANTE
    do I=1, NDIM-1
        do J=I+1, NDIM
            CTEN(J)=K(J, I)/K(I, I)
            do l=I+1, NDIM
                K(J, l)=K(J, l)-CTEN(J)*K(I, l)
            end do
        end do
    end do
  
    !do WHILE(T.LT.TMAX)
    !Guardo el resultado del paso anterior
    DISS = DISN
    VELS = VELN
    ACCS = ACCN
    R=0.D0
    ! Vector of effective loading forces at time t+dt
    ! R = F + [M]*(a1*dis + a2*vel + a3*acc) + [C]*(a1d*dis + a2d*vel + a3d*acc)
    ! Sumo ambos vectores, todas operaciones vectoriales
     
    ! Primera parte [M]*(a1*dis + a2*vel + a3*acc)
    G= A1*DISN+A2*VELN+A3*ACCN
    call MATVEC(MASA, G, G1, NDIM)
    R=F+G1

    ! Segunda parte [C]*(a1d*dis + a2d*vel + a3d*acc)
    G=A1D*DISN+A2D*VELN+A3D*ACCN
    call MATVEC(C, G, G1, NDIM)
    R=R+G1

    ! Resuelvo el sistema [Keff]q=Feff y obtengo los desplazamientos en t+dt
    call GAUSS(K, R, DISN, CTEN, NDIM)
     
    ! Aceleraciones el t+dt
    ACCN = A1*(DISN - DISS) - A2*VELS - A3*ACCS
    ! Velocidad a t+dt
    VELN = VELS + A4*ACCS + A5*ACCN
    TIEMPO=TIEMPO+H
     !write(1, *) T, DIS, VEL
 ! end do

end subroutine NEWMARK_METHOD

! SUBRUTINA PRODUCTO MATRIZ VECTOR
subroutine MATVEC(A, X, AUX, NDIM)
    implicit real(8)(A-H, O-Z)
    real(8)A(NDIM, NDIM), X(NDIM), AUX(NDIM)
  
    AUX=0.D0
    do I=1, NDIM
        do J=1, NDIM
            AUX(I)=AUX(I)+A(I, J)*X(J)
        end do
    end do
  
    return
end subroutine MATVEC

! METODO DE GAUSS JORDAN PARA SISTEMAS LINEALES
subroutine GAUSS(A, B, X, CTE, NDIM)
    implicit real(8) (A-H, O-Z)
    real(8) B(NDIM), X(NDIM), A(NDIM, NDIM), CTE(NDIM)
  
    !CCCC  -----> TRIANGULAR SUPERIOR AL VECTOR TERMINO INDEPENDIENTE
    do I=1, NDIM-1
        do J=I+1, NDIM
            B(J)=B(J)-CTE(J)*B(I)
        end do
    end do
     
    !CCCC  -----> RETORSUSTITUCION
    do I=NDIM, 1, -1
        S=0.D0
        do J=I+1, NDIM
            S=S+A(I, J)*X(J)
        end do
        X(I)=(B(I)-S)/A(I, I)
    end do
  
end subroutine GAUSS

subroutine ALLOC1
    use MATRICES
    use MAT2
    ALLOCATE(MASA(NDIM, NDIM), C(NDIM, NDIM), K(NDIM, NDIM))
    ALLOCATE(DISN(NDIM), DISS(NDIM), VELN(NDIM), VELS(NDIM), ACCN(NDIM), ACCS(NDIM))
    ALLOCATE(R(NDIM), G(NDIM), G1(NDIM), F(NDIM), CTEN(NDIM))
end subroutine ALLOC1

!subroutine DISTANCIA(X, Y, npoin, DISTX, DISTY, ILAUX, NMOVE, YPOSR, XREF, YREF)

!implicit real(8) (A-H, O-Z)
!real(8) X(npoin), Y(npoin), DISTX(NMOVE), DISTY(NMOVE)
!integer ILAUX(npoin)

!do I=1, NMOVE

!end do

!end subroutine

subroutine TRANSF(ALPHA, YPOSR, XREF, YREF)
    !use MGEOMETRIA
	use MeshData
    !use MALLOCAR
    !use MMOVIMIENTO
    use MNEWMARK

    implicit real(8) (A-H, O-Z)

    do I=1, NMOVE

        DISTX = X(ILAUX(I)) - XREF
        DISTY = Y(ILAUX(I)) - YREF+YPOSR

        DXPOS(I) = DCOS(ALPHA)*DISTX + DSIN(ALPHA)*DISTY - DISTX
        DYPOS(I) = -DSIN(ALPHA)*DISTX + DCOS(ALPHA)*DISTY - DISTY+YPOSR
    end do

end subroutine
