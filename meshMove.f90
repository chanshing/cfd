#define timer(func, store) call system_clock(start_t, rate); call func; call system_clock(end_t); store  = store + real(end_t - start_t)/real(rate);
module MeshMove
    real(8), dimension(:), allocatable :: b, pos_aux
    real(8), dimension(:), allocatable  :: dxpos, dypos
	real(8) yposr, alpha
    real(8), dimension(:, :), allocatable :: masa, c, k
    real(8), dimension(:), allocatable :: disn, diss, veln, vels, accn, accs, r, g, g1, f, cten
    real(8) fx(10), fy(10), rm(10), f_vx(10), f_vy(10)
    real(8), dimension(:), allocatable :: xpos, ypos

	private

	public :: setNewmarkCondition, fluidStructure
	public :: fx, fy, rm, f_vx, f_vy, xpos, ypos
contains
subroutine setNewmarkCondition
	implicit none

	call allocateNewmark

    !condiciones iniciales newmark
    DISN(1) = 0.d0; DISN(2) = 0.d0
    VELN(1) = 0.d0; VELN(2) = 0.d0
    ACCN = 0.d0
    F(1) = 0.d0
    F(2) = 0.d0
!    pi = dacos(-1.d0)
    YPOSR = DISN(1)
    ALPHA = DISN(2)
end subroutine setNewmarkCondition

subroutine fluidStructure(dtmin)
	use MeshData
	use Mlaplace
	use mvelocidades
	use BiconjGrad
	use mvariables
	use mvariabgen
	use InputData, only: XREF, YREF, NDIM
	use timers
	implicit none
	integer ipoin
	real(8) dtmin, alphav, yposrv

	if(.not. allocated(xpos)) allocate(xpos(npoin))
	if(.not. allocated(ypos)) allocate(ypos(npoin))
	if(.not. allocated(pos_aux)) allocate(pos_aux(npoin))
	if(.not. allocated(dxpos)) allocate(dxpos(nmove))
	if(.not. allocated(dypos)) allocate(dypos(nmove))
	if(.not. allocated(b)) allocate(b(npoin))

        timer(FORCES(NSETS, NSET_NUMB, IELEM_SETS, npoin, X, Y, P, FX, FY, RM, XREF, YREF), forces_t)
!        call FORCES(NSETS, NSET_NUMB, IELEM_SETS, npoin, X, Y, P, FX, FY, RM, XREF, YREF)
     
        F(1) = dcos(DISN(2))*FY(1) - dsin(DISN(2))*FX(1)
        F(2) = RM(1) + 0.1d0*dsin(DISN(2))*FX(1) - dcos(DISN(2))*0.1d0*FY(1)

        ALPHAV = DISN(2)
		YPOSRV = DISN(1)

!        call NEWMARK_METHOD(DTMIN)
        timer(NEWMARK_METHOD(DTMIN), newmark_t)
     
        ALPHA = DISN(2) - ALPHAV
		YPOSR = DISN(1) - YPOSRV

!        call TRANSF(ALPHA, YPOSR, XREF, YREF)
        timer(TRANSF(ALPHA, YPOSR, XREF, YREF), transf_t)
      
        !NODOS SIN MOVIMIENTO
		!$omp parallel do private(ipoin)
		do ipoin = 1 + NMOVE, NFIX_MOVE + NMOVE
			POS_AUX(ipoin) = 0.d0
		end do
		!$omp end parallel do
      
		!$omp parallel do private(ipoin)
		do ipoin = 1, NMOVE
			POS_AUX(ipoin) = DXPOS(ipoin)
		end do
		!$omp end parallel do

		!$omp parallel do private(ipoin)
		do ipoin = 1, npoin
			B(ipoin) = 0.d0
		end do
		!$omp end parallel do

!		call biCG(lap_sparse, lap_idx, lap_rowptr, lap_diag, xpos, b, pos_aux(1:nnmove), ilaux(1:nnmove), npoin, nnmove)
		timer(biCG(lap_sparse, lap_idx, lap_rowptr, lap_diag, xpos, b, pos_aux(1:nnmove), ilaux(1:nnmove), npoin, nnmove), grad_t)
             
		!$omp parallel do private(ipoin)
		do ipoin = 1, npoin
			X(ipoin) = X(ipoin) + XPOS(ipoin)
			W_X(ipoin) = XPOS(ipoin)/DTMIN
		end do
		!$omp end parallel do
     
		!$omp parallel do private(ipoin)
		do ipoin = 1, NMOVE
			POS_AUX(ipoin) = DYPOS(ipoin)
		end do
		!$omp end parallel do
     
		!$omp parallel do private(ipoin)
		do ipoin = 1, npoin
			B(ipoin) = 0.d0
		end do
		!$omp end parallel do

!		call biCG(lap_sparse, lap_idx, lap_rowptr, lap_diag, ypos, b, pos_aux(1:nnmove), ilaux(1:nnmove), npoin, nnmove)
		timer(biCG(lap_sparse, lap_idx, lap_rowptr, lap_diag, ypos, b, pos_aux(1:nnmove), ilaux(1:nnmove), npoin, nnmove), grad_t)
	
		!$omp parallel do private(ipoin)
		do ipoin = 1, npoin
			Y(ipoin) = Y(ipoin) + YPOS(ipoin)
			W_Y(ipoin) = YPOS(ipoin)/DTMIN
		end do
		!$omp end parallel do
end subroutine

subroutine allocateNewmark
    !use MATRICES
    !use MAT2
	use InputData, only: NDIM
    ALLOCATE(MASA(NDIM, NDIM), C(NDIM, NDIM), K(NDIM, NDIM))
    ALLOCATE(DISN(NDIM), DISS(NDIM), VELN(NDIM), VELS(NDIM), ACCN(NDIM), ACCS(NDIM))
    ALLOCATE(R(NDIM), G(NDIM), G1(NDIM), F(NDIM), CTEN(NDIM))
end subroutine allocateNewmark

subroutine FORCES(NSETS, NSET_NUMB, IELEM_SETS, npoin &
    , X, Y, P , FX, FY, RM, X_ROT, Y_ROT)
!cccc  calculo de las fuerzas en  cccc
!cccc  cada uno de los sets       cccc
    !use MFUERZAS
	use MeshData, only: ISET, IPER_MASTER, IPER_SLAVE, IPER_AUX
    implicit real(8) (A-H, O-Z)
  
    INTEGER IELEM_SETS(10)
  
    real(8) X(npoin), Y(npoin)
    real(8) P(npoin)
    real(8) FX(10), FY(10), RM(10)
  
    FX(1:NSET_NUMB)=0.d0
    FY(1:NSET_NUMB)=0.d0
    RM(1:NSET_NUMB)=0.d0
  
    do ISET_NUMB=1, NSET_NUMB
        do II=1, IELEM_SETS(ISET_NUMB)
            N1=ISET(1, ISET_NUMB, II)
            N2=ISET(2, ISET_NUMB, II)
        
            D_PRESS=(P(N1)+P(N2))/2.d0
            RLX=X(N1)-X(N2)
            RLY=Y(N2)-Y(N1)
        
            DFX=D_PRESS*RLY
            DFY=D_PRESS*RLX
        
            FX(ISET_NUMB)=FX(ISET_NUMB)+DFX
            FY(ISET_NUMB)=FY(ISET_NUMB)+DFY
        
            XC=(X(N1)+X(N2))/2.d0
            YC=(Y(N2)+Y(N1))/2.d0
        
            RM(ISET_NUMB)=RM(ISET_NUMB)+DFY*(XC-X_ROT)-DFX*(YC-Y_ROT)
        
        end do
    end do
    return
end subroutine FORCES

subroutine NEWMARK_METHOD(H)
!--------------------------------------------------------------------------       
!                 Newmark's Direct Integration Method
!--------------------------------------------------------------------------
! Code written by : - Ing. German Weht                                    
!                     Researcher & PhD Student                              
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
    !use MATRICES
    !use MAT2
	use InputData, only: NDIM
    implicit real(8)(A-H, O-Z)
    !real(8), ALLOCATABLE:: M(:, :), C(:, :), K(:, :)
    !real(8), ALLOCATABLE:: DIS(:), DISS(:), VEL(:), VELS(:), ACC(:), ACCS(:), R(:), G(:), G1(:), F(:), CTE(:)
    PI=DACOS(-1.d0)
    GR=1.d0!9.81D0 ! gravedad


  

  
    ! **** DEFINICION DE CONSTANTES DEL METODO POR DEFAULT GAMA=1/2 BETA=1/4 ****
    GAMA=.5D0 ; BETA=.25D0

    A1 = 1.d0/(BETA*H*H) ; A2 = 1.d0/(BETA*H) ; A3 = 1.d0/(2.d0*BETA) - 1.d0 ; A4 = (1.d0 - GAMA)*H ; A5 = GAMA*H
    A1D = GAMA/(BETA*H) ; A2D = GAMA/BETA - 1.d0 ; A3D = 0.5D0*H*(GAMA/BETA - 2.d0)
  

  
    !ALLOCATE(M(NDIM, NDIM), C(NDIM, NDIM), K(NDIM, NDIM))
    !ALLOCATE(DIS(NDIM), DISS(NDIM), VEL(NDIM), VELS(NDIM), ACC(NDIM), ACCS(NDIM))
    !ALLOCATE(R(NDIM), G(NDIM), G1(NDIM), F(NDIM), CTE(NDIM))
  

    !***** MATRICES *****
    ! MATRIZ DE MASA
    MASA(1, 1)=51.5d0   ; MASA(1, 2)=-51.5d0*.0429d0
    MASA(2, 1)=-51.5d0*.0429d0   ; MASA(2, 2)=2.275d0

    ! MATRIZ DE AMORTIGUAMIENTO
    C(1, 1)=32.358D0/GR ;  C(1, 2)=0.d0
    C(2, 1)=0.d0   ;  C(2, 2)=5.718D0/GR

    ! MATRIZ DE RIGUIDEZ
    K(1, 1)=50828.463D0/GR   ; K(1, 2)=0.d0
    K(2, 1)=0.d0   ; K(2, 2)=35923.241D0/GR

    ! VECTOR DE FUERZAS EXTERNAS
  
 
    ! CONDICIONES INICIALES
    !DIS(1)=0.05D0; DIS(2)=PI/8.d0
    !VEL=0.d0
    !ACC=0.d0

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
    R=0.d0
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

subroutine MATVEC(A, X, AUX, NDIM)
! SUBRUTINA PRODUCTO MATRIZ VECTOR
    implicit real(8)(A-H, O-Z)
    real(8)A(NDIM, NDIM), X(NDIM), AUX(NDIM)
  
    AUX=0.d0
    do I=1, NDIM
        do J=1, NDIM
            AUX(I)=AUX(I)+A(I, J)*X(J)
        end do
    end do
  
    return
end subroutine MATVEC

subroutine GAUSS(A, B, X, CTE, NDIM)
! METODO DE GAUSS JORDAN PARA SISTEMAS LINEALES
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
        S=0.d0
        do J=I+1, NDIM
            S=S+A(I, J)*X(J)
        end do
        X(I)=(B(I)-S)/A(I, I)
    end do
  
end subroutine GAUSS

subroutine TRANSF(ALPHA, YPOSR, XREF, YREF)
    !use MGEOMETRIA
    !use MALLOCAR
    !use MMOVIMIENTO
	use MeshData
    !use MNEWMARK
    implicit real(8) (A-H, O-Z)

    do I=1, NMOVE

        DISTX = X(ILAUX(I)) - XREF
        DISTY = Y(ILAUX(I)) - YREF+YPOSR

        DXPOS(I) = dcos(ALPHA)*DISTX + dsin(ALPHA)*DISTY - DISTX
        DYPOS(I) = -dsin(ALPHA)*DISTX + dcos(ALPHA)*DISTY - DISTY+YPOSR
    end do

end subroutine
end module MeshMove
