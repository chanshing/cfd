module Mnormales
    integer, private :: nnormvp
    integer, dimension(:), allocatable, private :: normvpIdx
    real(8), dimension(:), allocatable, private :: n_x, n_y, num_l_x, num_l_y, den_l
contains
subroutine normales
    !use mgeometria
	use meshdata, only: ivn, nelnorm, x, y, npoin
	implicit none
	real(8) l_edge_x, l_edge_y, l_edge, l_x, l_y, l_nrm
	integer ielem, ipoin

	if(.not.allocated(normvpIdx)) allocate(normvpIdx(npoin))
	if(.not.allocated(n_x)) allocate(n_x(npoin))
	if(.not.allocated(n_y)) allocate(n_y(npoin))
	if(.not.allocated(num_l_x)) allocate(num_l_x(npoin))
	if(.not.allocated(num_l_y)) allocate(num_l_y(npoin))
	if(.not.allocated(den_l)) allocate(den_l(npoin))

    num_l_x = 0.d0 
    num_l_y = 0.d0 
    den_l = 0.d0 
    nnormvp = 0 

    !$omp parallel &
    !$omp private(l_edge_x, l_edge_y, l_edge, ielem)
    !$omp do
    do ielem = 1, nelnorm
        l_edge_x = y(ivn(2, ielem)) - y(ivn(1, ielem))
        l_edge_y = -(x(ivn(2, ielem)) - x(ivn(1, ielem)))
        l_edge = dsqrt(l_edge_x*l_edge_x + l_edge_y*l_edge_y)
     
        !$omp atomic
        num_l_x(ivn(1, ielem)) = num_l_x(ivn(1, ielem)) + l_edge_x
        !$omp atomic
        num_l_y(ivn(1, ielem)) = num_l_y(ivn(1, ielem)) + l_edge_y
        !$omp atomic
        den_l(ivn(1, ielem)) = den_l(ivn(1, ielem)) + l_edge
        !$omp atomic
        num_l_x(ivn(2, ielem)) = num_l_x(ivn(2, ielem)) + l_edge_x
        !$omp atomic
        num_l_y(ivn(2, ielem)) = num_l_y(ivn(2, ielem)) + l_edge_y
        !$omp atomic
        den_l(ivn(2, ielem)) = den_l(ivn(2, ielem)) + l_edge
    end do
    !$omp end do
    !$omp end parallel

    do ipoin = 1, npoin
        if (den_l(ipoin) > 1.d-6) then
            l_x = num_l_x(ipoin)/den_l(ipoin)
            l_y = num_l_y(ipoin)/den_l(ipoin)
            l_nrm = dsqrt(l_x*l_x + l_y*l_y)

            if (l_nrm > 0.2d0) then
                nnormvp = nnormvp + 1
                normvpIdx(nnormvp) = ipoin
                n_x(nnormvp) = l_x/l_nrm
                n_y(nnormvp) = l_y/l_nrm
            end if
        end if
    end do
end subroutine normales

subroutine normalvel
    !use mallocar
    !use mnormales
	!use meshdata
    use mvelocidades
	implicit none
	integer inormvp, ipoin 
	real(8) vx, vy, wx, wy, rnx1, rny1, utp

    do inormvp = 1, nnormvp
        ipoin = normvpIdx(inormvp)

        vx = vel_x(ipoin); vy = vel_y(ipoin)
        wx = w_x(ipoin); wy = w_y(ipoin)

        rnx1 = n_x(inormvp); rny1 = n_y(inormvp)
     
        utp = -rny1*(vx - wx) + rnx1*(vy - wy)
     
        vel_x(ipoin) = -rny1*utp + wx
        vel_y(ipoin) = rnx1*utp + wy
    end do
end subroutine normalvel
end module Mnormales

subroutine deriv(hmin) !<----- PARA QUE SE USA HMIN?
!------------------------------------------
! Calcula dNx, dNy y area de cada elemento
!------------------------------------------
    !use MGEOMETRIA
	use MeshData, only: X, Y, inpoel, area, HH, HHX, HHY, dNx, dNy, nelem
	implicit none
	integer ielem
	real(8) hmin
	real(8) X_loc(3), Y_loc(3), a

	!$OMP PARALLEL &
	!$OMP PRIVATE(X_loc, Y_loc, ielem, a)
	!$OMP do
    do ielem = 1,nelem
        X_loc(:) = X(inpoel(:, ielem))
        Y_loc(:) = Y(inpoel(:, ielem))
     
        area(ielem) = (X_loc(2)*Y_loc(3) + X_loc(3)*Y_loc(1) + X_loc(1)*Y_loc(2)&
				- (X_loc(2)*Y_loc(1) + X_loc(3)*Y_loc(2) + X_loc(1)*Y_loc(3)))/2.d0
		a = area(ielem)

        dNx(1, ielem) = (Y_loc(2) - Y_loc(3))/(2.D0*a)
        dNx(2, ielem) = (Y_loc(3) - Y_loc(1))/(2.D0*a)
        dNx(3, ielem) = (Y_loc(1) - Y_loc(2))/(2.D0*a)
        dNy(1, ielem) = (X_loc(3) - X_loc(2))/(2.D0*a)
        dNy(2, ielem) = (X_loc(1) - X_loc(3))/(2.D0*a)
        dNy(3, ielem) = (X_loc(2) - X_loc(1))/(2.D0*a)
     
        HH(ielem) = dsqrt(area(ielem))
        HHX(ielem) = abs(min( X_loc(3) - X_loc(2), X_loc(1) - X_loc(3), X_loc(2) - X_loc(1) ))
        HHY(ielem) = abs(min( Y_loc(3) - Y_loc(2), Y_loc(1) - Y_loc(3), Y_loc(2) - Y_loc(1) ))
    end do
	!$OMP end do
	!$OMP end PARALLEL

	hmin = minval(hh) 
	if(hmin > 1.d10) hmin = 1.d10
end subroutine deriv

subroutine MASAS
!----------------------------------------
! Ensambla la matriz de masas condensada
!----------------------------------------
    !use MGEOMETRIA
	use MeshData, only: M, inpoel, area, nelem
	implicit none
	integer ielem
  
    M = 0.d0
    do ielem = 1,nelem
        M(inpoel(:, ielem)) = M(inpoel(:, ielem)) + area(ielem)/3.d0
    end do
end subroutine MASAS

subroutine deltat(dtmin, dt)
!-------------------------------
!		Calcula delta t
!-------------------------------
    !use MALLOCAR
    !use DATOS_ENTRADA
	use MeshData, only: inpoel, nelem, hh
	use InputData, only: FSAFE, fr, gama
    use MVELOCIDADES
    use MVARIABLES
	implicit none
    real(8) DT(nelem), dtmin
	integer ielem, i, ipoin
	real(8) T_iel, vumax, vvmax, vu, vv, vel, dtelem, cota, vc

    DTMIN = 1.D20
  
    do ielem = 1, nelem
        VUMAX = 0.D0
        VVMAX = 0.D0

        T_iel = (T(inpoel(1, ielem)) + T(inpoel(2, ielem)) + T(inpoel(3, ielem)))/3.D0
     
        !GM = GAMA-1.D0 !(GAMM(inpoel(1, ielem))+GAMM(inpoel(2, ielem))+GAMM(inpoel(3, ielem)))/3.D0
        VC = DSQRT(GAMA*FR*T_iel)
     
        !CCCC ----> CALCULO DE LA MAXIMA VELOCIDAD ELEMENTAL
        do i = 1, 3
            ipoin = inpoel(i, ielem)
            VU = dabs(VEL_X(ipoin) - W_X(ipoin))
            VV = dabs(VEL_Y(ipoin) - W_Y(ipoin))
        
            if (VU.GT.VUMAX) VUMAX = VU
            if (VV.GT.VVMAX) VVMAX = VV
        end do
     
        VEL = dsqrt(VUMAX**2 + VVMAX**2)
     
        DTELEM = .5D0*HH(ielem)/(VC+VEL)
        DT(ielem) = DTELEM*FSAFE
        if (DTELEM < DTMIN) DTMIN = DTELEM*FSAFE
    end do
  
    COTA = 10.D0*DTMIN !EL VALOR 10 ESTA PUESTO A OJO #MODIFICAR SI ES NECESARIO#
    do ielem = 1, nelem
        if(DT(ielem) > COTA) DT(ielem) = COTA
    end do
    return
end subroutine deltat

subroutine CUARTO_ORDEN(U, U_n, FR, gamm)
!---------------------------------------------
!		Calculo de la proyeccion U_n
!---------------------------------------------
    !use MALLOCAR
    !use MGEOMETRIA
	use MeshData
    implicit real(8) (A-H,O-Z)
    real(8) U(4,npoin), U_n(4,npoin), GAMM(npoin)
  
    real(8) A1(3,4,3), A2(3,4,3), Adv(4,3), U_loc(4,3), Ux(4), Uy(4), Nx(3), Ny(3), U_n_tmp(4,3)
    real(8) sp(3,3) !Funciones de forma
	real(8) c(3), vx(3), vy(3), V_sq(3), TEMP(3), e(3)
    INTEGER ipoin(3)

    U_n = 0.D0
	sp(:,1) = (/ .5D0, .5D0, 0.D0 /)
	sp(:,2) = (/ 0.D0, .5D0, .5D0 /)
    sp(:,3) = (/ .5D0, 0.D0, .5D0 /)

    !$OMP PARALLEL &
    !$OMP PRIVATE(ielem,ipoin,gama,TEMP,FMU,Nx,Ny,Ux,Uy,AR,&
    !$OMP U_loc,vx,vy,V_sq,e,c,A1,A2,Adv,U_n_tmp,i)

    !$OMP do
    do ielem = 1,nelem
        U_n_tmp = 0.D0
        ipoin = inpoel(:, ielem)

        gama = (GAMM(ipoin(1))+GAMM(ipoin(2))+GAMM(ipoin(3)))/3.D0

        Nx = dNx(:, ielem)
        Ny = dNy(:, ielem)

        Ux(:) = U(:,ipoin(1))*Nx(1) + U(:,ipoin(2))*Nx(2) + U(:,ipoin(3))*Nx(3)
        Uy(:) = U(:,ipoin(1))*Ny(1) + U(:,ipoin(2))*Ny(2) + U(:,ipoin(3))*Ny(3)

        AR = area(ielem)/3.D0
        
        !CCCC  ----> INTEGRO LAS VARIABLES EN LOS PUNTOS DE GAUSS
        U_loc(:,1) = sp(1,1)*U(:,ipoin(1)) + sp(2,1)*U(:,ipoin(2)) + sp(3,1)*U(:,ipoin(3))
        U_loc(:,2) = sp(1,2)*U(:,ipoin(1)) + sp(2,2)*U(:,ipoin(2)) + sp(3,2)*U(:,ipoin(3))
        U_loc(:,3) = sp(1,3)*U(:,ipoin(1)) + sp(2,3)*U(:,ipoin(2)) + sp(3,3)*U(:,ipoin(3))

        !CCCC  ----> DEFINO VARIABLES PRIMITIVAS
        vx(:) = U_loc(2,:)/U_loc(1,:)
        vy(:) = U_loc(3,:)/U_loc(1,:)
        e(:) = U_loc(4,:)/U_loc(1,:)
        V_sq(:) = vx*vx+vy*vy
        temp(:) = (gama-1.D0)/FR*(e-.5D0*V_sq) !FR = CTE. UNIVERSAL DE LOS GASES
        c(:) = dsqrt(gama*fr*temp(:))

		! = = = GAUSS QUAD = ==
		do i = 1,3
        !A1(0,:,i) = (/ 0.D0, 1.D0, 0.D0, 0.D0 /)
        A1(1,:,i) = (/ (gama-1.D0)/2.D0*V_sq(i)-vx(i)*vx(i),&
					(3.D0-gama)*vx(i),&
					-(gama-1.D0)*vy(i),&
					(gama-1.D0) /)
        A1(2,:,i) = (/ -vx(i)*vy(i),&
					vy(i),&
					vx(i),&
					0.D0 /)
        A1(3,:,i) = (/ ((gama-1.D0)*V_sq(i)-gama*e(i))*vx(i),&
					gama*e(i)-(gama-1.D0)/2.D0*V_sq(i)-(gama-1.D0)*vx(i)*vx(i),&
					-(gama-1.D0)*vx(i)*vy(i),&
					gama*vx(i) /)
        !A2(0,:,i) = (/ 0.D0, 0.D0, 1.D0, 0.D0 /)
        A2(1,:,i) = (/ -vx(i)*vy(i),&
					vy(i),&
					vx(i),&
					0.D0 /)
        A2(2,:,i) = (/ (gama-1.D0)/2.D0*V_sq(i)-vy(i)*vy(i),&
					-(gama-1.D0)*vx(i),&
					(3.D0-gama)*vy(i),&
					(gama-1.D0) /)
        A2(3,:,i) = (/ ((gama-1.D0)*V_sq(i)-gama*e(i))*vy(i),&
					-(gama-1.D0)*vx(i)*vy(i),&
					gama*e(i)-(gama-1.D0)/2.D0*V_sq(i)-(gama-1.D0)*vy(i)*vy(i),&
					gama*vy(i) /)

        Adv(1,i) = Ux(2) + Uy(3)
        Adv(2:4,i) = A1(:,1,i)*Ux(1) + A1(:,2,i)*Ux(2) + A1(:,3,i)*Ux(3) + A1(:,4,i)*Ux(4)&
					+A2(:,1,i)*Uy(1) + A2(:,2,i)*Uy(2) + A2(:,3,i)*Uy(3) + A2(:,4,i)*Uy(4)
		end do
		! = = = end GAUSS QUAD = ==

        U_n_tmp(:,1) = Adv(:,1)*sp(1,1)*AR + Adv(:,2)*sp(1,2)*AR + Adv(:,3)*sp(1,3)*AR
        U_n_tmp(:,2) = Adv(:,1)*sp(2,1)*AR + Adv(:,2)*sp(2,2)*AR + Adv(:,3)*sp(2,3)*AR
        U_n_tmp(:,3) = Adv(:,1)*sp(3,1)*AR + Adv(:,2)*sp(3,2)*AR + Adv(:,3)*sp(3,3)*AR

		do i = 1,3
		!$OMP ATOMIC
        U_n(1,ipoin(i)) = U_n(1,ipoin(i)) + U_n_tmp(1,i)
		!$OMP ATOMIC
        U_n(2,ipoin(i)) = U_n(2,ipoin(i)) + U_n_tmp(2,i)
		!$OMP ATOMIC
        U_n(3,ipoin(i)) = U_n(3,ipoin(i)) + U_n_tmp(3,i)
		!$OMP ATOMIC
        U_n(4,ipoin(i)) = U_n(4,ipoin(i)) + U_n_tmp(4,i)
		end do

    end do
    !$OMP end do

    !$OMP do
    do inod = 1,npoin
        U_n(:,inod) = -U_n(:,inod)/M(inod)
    end do
    !$OMP end do
    !$OMP end PARALLEL

    return
end subroutine CUARTO_ORDEN

subroutine ESTAB(U,T,GAMA,FR,RMU &
    ,DTMIN,RHO_inf,T_inf,U_inf,V_inf, GAMM)
!-------------------------------------------
! Calculo de los terminos de estabilizacion
!-------------------------------------------
    !use MALLOCAR
    !use MGEOMETRIA
	use MeshData
    use MVELOCIDADES
    use MESTABILIZACION
	use InputData, only: cte
    implicit real(8) (A-H,O-Z)

    real(8) U(4,npoin), T(npoin), GAMM(npoin)
	integer IPOIN(3)
  
    !CCCC  ----> CTES DE CALCULO PARA SHOCK-CAPTURING
    CC = DSQRT(GAMA*FR*T_inf)
    VEL2 = DSQRT(U_inf*U_inf + V_inf*V_inf)
  
    !$OMP PARALLEL do &
    !$OMP PRIVATE(i, ielem,IPOIN,GM,TAU,H_RGNE,H_RGN,H_JGN,RHO_ELEM,VX,VY,wx,wy,VEL2,DRX,DRY,DR2,&
    !$OMP DTX,DTY,DT2,DUX,DUY,DU2,RTX,RTY,RJX,RJY,RUX,RUY,TEMP,C,FMU,TERM_1,TERM_2,H_RGN1,H_RGN2,tmp)
    do ielem = 1,nelem
        IPOIN = inpoel(:, ielem)
        GM = (GAMM(IPOIN(1)) + GAMM(IPOIN(2)) + GAMM(IPOIN(3)))/3.D0
        TAU = 0.D0
        H_RGNE = 0.D0
        H_RGN = 0.D0
        H_JGN = 0.D0
     
        !CCCC  ----> VARIABLES ELEMENTALES
        RHO_ELEM = (U(1,IPOIN(1)) + U(1,IPOIN(2)) + U(1,IPOIN(3)))/3.D0
        VX = (VEL_X(IPOIN(1)) + VEL_X(IPOIN(2)) + VEL_X(IPOIN(3)))/3.D0
        VY = (VEL_Y(IPOIN(1)) + VEL_Y(IPOIN(2)) + VEL_Y(IPOIN(3)))/3.D0
     
        !PARTES NUEVAS
        wx = (W_X(IPOIN(1)) + W_X(IPOIN(2)) + W_X(IPOIN(3)))/3.D0
        wy = (W_Y(IPOIN(1)) + W_Y(IPOIN(2)) + W_Y(IPOIN(3)))/3.D0
        VX = VX-wx ; VY = VY-wy
        VEL2 = DSQRT(VX*VX+VY*VY)
     
        !CCCC  ----> DERIVADA DE RHO
        DRX = U(1,IPOIN(1))*dNx(1, ielem) + U(1,IPOIN(2))*dNx(2, ielem) + U(1,IPOIN(3))*dNx(3, ielem)
        DRY = U(1,IPOIN(1))*dNy(1, ielem) + U(1,IPOIN(2))*dNy(2, ielem) + U(1,IPOIN(3))*dNy(3, ielem)
        DR2 = DSQRT(DRX*DRX+DRY*DRY)+1.D-20
        !CCCC  ----> DERIVADA DE TEMPERATURA
        DTX = T(IPOIN(1))*dNx(1, ielem) + T(IPOIN(2))*dNx(2, ielem) + T(IPOIN(3))*dNx(3, ielem)
        DTY = T(IPOIN(1))*dNy(1, ielem) + T(IPOIN(2))*dNy(2, ielem) + T(IPOIN(3))*dNy(3, ielem)
        DT2 = DSQRT(DTX*DTX+DTY*DTY)+1.D-20
        !CCCC  ----> DERIVADA DE LA VELOCIDAD
        DUX = VEL2*dNx(1, ielem)+VEL2*dNx(2, ielem)+VEL2*dNx(3, ielem)
        DUY = VEL2*dNy(1, ielem)+VEL2*dNy(2, ielem)+VEL2*dNy(3, ielem)
        DU2 = DSQRT(DUX*DUX+DUY*DUY)+1.D-20
     
        !CCCC  ----> VECTOR UNIDAD THETA
        RTX = DTX/DT2
        RTY = DTY/DT2
        !CCCC  ----> VECTOR UNIDAD J
        RJX = DRX/DR2
        RJY = DRY/DR2
        !CCCC  ----> VECTOR UNIDAD VELOCIDAD
        RUX = DUX/DU2
        RUY = DUY/DU2
     
        TEMP = (T(IPOIN(1))+T(IPOIN(2))+T(IPOIN(3)))/3.D0
        C = DSQRT(GM*FR*TEMP)
     
        FMU= RMU*162.6/(TEMP-110.55)*(temp/273.15)**.75D0 !SUTHERLAND

        do i = 1,3
            TERM_1 = dabs(VX*dNx(i, ielem)+VY*dNy(i, ielem))
            TERM_2 = dabs(RJX*dNx(i, ielem)+RJY*dNy(i, ielem))
        
            H_RGN1 = dabs(RTX*dNx(i, ielem)+RTY*dNy(i, ielem)) !CALCULO PARA ECU. ENERGIA
            H_RGN2 = dabs(RUX*dNx(i, ielem)+RUY*dNy(i, ielem)) !CALCULO PARA ECU. MOMENTO
        
            TAU = TAU+TERM_1+TERM_2*C
        
            H_RGNE = H_RGNE+H_RGN1
            H_RGN = H_RGN+H_RGN2
            H_JGN = H_JGN+TERM_2
        end do
     
        TAU = 1.D0/TAU
        H_RGNE = 2.D0/H_RGNE
        H_RGN = 2.D0/H_RGN

        if(H_RGN.GT.1.D1) H_RGN = 0.D0

        H_JGN = 2.D0/H_JGN

        if(H_JGN.GT.1.D1) H_JGN = 0.D0

        SHOC(ielem) = (DR2*H_JGN/RHO_ELEM + (DR2*H_JGN/RHO_ELEM)**2)*.5D0*C**2*H_JGN/(2.D0*C)
     
        tmp = (1.D0/TAU**2.D0 +(2.D0/DTMIN)**2.D0)**(-.5D0)
        T_SUGN1(ielem)= tmp
        T_SUGN2(ielem)= tmp
        T_SUGN3(ielem)= tmp

     	!if(RMU.NE.0.D0)THEN
     	!   TAU_SUNG3=   H_RGN**2.D0/(4.D0*FMU/RHO_inf)
     	!   TAU_SUNG3_E= H_RGNE**2.D0/(4.D0*FMU/RHO_inf)
        
     	!   T_SUGN2(ielem) = (RESUMEN+1.D0/TAU_SUNG3**2.D0)**(-.5D0)
     	!   T_SUGN3(ielem) = (RESUMEN+1.D0/TAU_SUNG3_E**2.D0)**(-.5D0)
     	!end if
     
    end do
    !$OMP end PARALLEL do
  
    return
end subroutine ESTAB
      
!CCCC------------------------------------------------------------!CCCC
!CCCC  CALCULO DE LA INTEGRAL DE LOS TERMINOS DE ESTABILIZACION  !CCCC
!CCCC------------------------------------------------------------!CCCC
subroutine calcRHS(U, U_n, rhs, P, FR, RMU, FK, FCV, T_inf, dtl, gamm)
    !use MALLOCAR
    !use MGEOMETRIA
	use MeshData
    use MESTABILIZACION
	use InputData, only: cte
    implicit real(8) (A-H,O-Z)
    real(8) U(4,npoin),P(npoin),U_n(4,npoin),T(npoin), dtl(nelem), gamm(npoin)
    real(8) rhs(4,npoin)
    real(8) A1(3,4,3), A2(3,4,3), Adv(4,3), AA1(4), AA2(4), Adv_phi(4,3), U_loc(4,3), Ux(4), Uy(4), Nx(3), Ny(3)
    real(8) tau(4), choq(3), phi_loc(4,3), rhs_tmp(4,3)
    real(8) sp(3,3) !Funciones de forma
	real(8) c(3), vx(3), vy(3), V_sq(3), temp(3), e(3)
    INTEGER ipoin(3)

	sp(:,1) = (/ .5D0, .5D0, 0.D0 /)
	sp(:,2) = (/ 0.D0, .5D0, .5D0 /)
    sp(:,3) = (/ .5D0, 0.D0, .5D0 /)

    !$OMP PARALLEL &
    !$OMP PRIVATE(ielem,ipoin,gama,temp,FMU,Nx,Ny,Ux,Uy,AR,HLONG,HLONGX,HLONGY,tau,ALFA_MU,i,&
    !$OMP U_loc,phi_loc,vx,vy,V_sq,e,C,A1,A2,Adv,AA1,AA2,Adv_phi,choq,rhs_tmp)

    !$OMP do
    do ielem = 1,nelem
        rhs_tmp = 0.D0
        ipoin = inpoel(:, ielem)

        gama = (GAMM(ipoin(1))+GAMM(ipoin(2))+GAMM(ipoin(3)))/3.D0
        !temp = (T(ipoin(1))+T(ipoin(2))+T(ipoin(3)))/3.D0
        !FMU= 1.716D-5*162.6/(temp-110.55)*(temp/273.15)**.75D0     !SUTHERLAND

        Nx = dNx(:, ielem)
        Ny = dNy(:, ielem)

        Ux(:) = U(:,ipoin(1))*Nx(1) + U(:,ipoin(2))*Nx(2) + U(:,ipoin(3))*Nx(3)
        Uy(:) = U(:,ipoin(1))*Ny(1) + U(:,ipoin(2))*Ny(2) + U(:,ipoin(3))*Ny(3)

        AR = area(ielem)*dtl(ielem)/3.D0
        !CCCC  ----> LONG. CARACTERISTICA
        HLONGX = HHX(ielem)
        HLONGY = HHY(ielem)
        HLONG = DSQRT(area(ielem))
        !CCCC  ----> ESTAB. TEZDUyAR
        tau(1) = T_SUGN1(ielem)
        tau(2) = T_SUGN2(ielem)
        tau(3) = T_SUGN2(ielem)
        tau(4) = T_SUGN3(ielem)
        ALFA_MU = SHOC(ielem)
		choq(1) = ALFA_MU*AR*CTE
        choq(2) = ALFA_MU*AR*CTE
        choq(3) = ALFA_MU*AR*CTE

        !CCCC  ----> INTEGRO LAS VARIABLES EN LOS PUNTOS DE GAUSS
        U_loc(:,1) = sp(1,1)*U(:,ipoin(1)) + sp(2,1)*U(:,ipoin(2)) + sp(3,1)*U(:,ipoin(3))
        U_loc(:,2) = sp(1,2)*U(:,ipoin(1)) + sp(2,2)*U(:,ipoin(2)) + sp(3,2)*U(:,ipoin(3))
        U_loc(:,3) = sp(1,3)*U(:,ipoin(1)) + sp(2,3)*U(:,ipoin(2)) + sp(3,3)*U(:,ipoin(3))

        phi_loc(:,1) = sp(1,1)*U_n(:,ipoin(1)) + sp(2,1)*U_n(:,ipoin(2)) + sp(3,1)*U_n(:,ipoin(3))
        phi_loc(:,2) = sp(1,2)*U_n(:,ipoin(1)) + sp(2,2)*U_n(:,ipoin(2)) + sp(3,2)*U_n(:,ipoin(3))
        phi_loc(:,3) = sp(1,3)*U_n(:,ipoin(1)) + sp(2,3)*U_n(:,ipoin(2)) + sp(3,3)*U_n(:,ipoin(3))

        !CCCC  ----> DEFINO VARIABLES PRIMITIVAS
        vx(:) = U_loc(2,:)/U_loc(1,:)
        vy(:) = U_loc(3,:)/U_loc(1,:)
        e(:) = U_loc(4,:)/U_loc(1,:)
        V_sq(:) = vx*vx+vy*vy
        temp(:) = (gama-1.D0)/FR*(e-.5D0*V_sq) !FR = CTE. UNIVERSAL DE LOS GASES
        c(:) = DSQRT(gama*FR*temp(:))

		! = = = GAUSS QUAD = ==
		do i = 1,3 
        !A1(0,:,i) = (/ 0.D0, 1.D0, 0.D0, 0.D0 /)
        A1(1,:,i) = (/ (gama-1.D0)/2.D0*V_sq(i)-vx(i)*vx(i),&
					(3.D0-gama)*vx(i),&
					-(gama-1.D0)*vy(i),&
					(gama-1.D0) /)
        A1(2,:,i) = (/ -vx(i)*vy(i),&
					vy(i),&
					vx(i),&
					0.D0 /)
        A1(3,:,i) = (/ ((gama-1.D0)*V_sq(i)-gama*e(i))*vx(i),&
					gama*e(i)-(gama-1.D0)/2.D0*V_sq(i)-(gama-1.D0)*vx(i)*vx(i),&
					-(gama-1.D0)*vx(i)*vy(i),&
					gama*vx(i) /)
        !A2(0,:,i) = (/ 0.D0, 0.D0, 1.D0, 0.D0 /)
        A2(1,:,i) = (/ -vx(i)*vy(i),&
					vy(i),&
					vx(i),&
					0.D0 /)
        A2(2,:,i) = (/ (gama-1.D0)/2.D0*V_sq(i)-vy(i)*vy(i),&
					-(gama-1.D0)*vx(i),&
					(3.D0-gama)*vy(i),&
					(gama-1.D0) /)
        A2(3,:,i) = (/ ((gama-1.D0)*V_sq(i)-gama*e(i))*vy(i),&
					-(gama-1.D0)*vx(i)*vy(i),&
					gama*e(i)-(gama-1.D0)/2.D0*V_sq(i)-(gama-1.D0)*vy(i)*vy(i),&
					gama*vy(i) /)

        Adv(1,i) = Ux(2) + Uy(3)
        Adv(2:4,i) = A1(:,1,i)*Ux(1) + A1(:,2,i)*Ux(2) + A1(:,3,i)*Ux(3) + A1(:,4,i)*Ux(4)&
					+A2(:,1,i)*Uy(1) + A2(:,2,i)*Uy(2) + A2(:,3,i)*Uy(3) + A2(:,4,i)*Uy(4)
		end do 
		! = = = end GAUSS QUAD = ==

		Adv_phi = Adv + phi_loc

        AA1(1) = (Adv_phi(2,1) + Adv_phi(2,2) + Adv_phi(2,3))*AR
        AA1(2:4) = (A1(:,1,1)*Adv_phi(1,1) + A1(:,2,1)*Adv_phi(2,1) + A1(:,3,1)*Adv_phi(3,1) + A1(:,4,1)*Adv_phi(4,1) + &
				   A1(:,1,2)*Adv_phi(1,2) + A1(:,2,2)*Adv_phi(2,2) + A1(:,3,2)*Adv_phi(3,2) + A1(:,4,2)*Adv_phi(4,2) + &
				   A1(:,1,3)*Adv_phi(1,3) + A1(:,2,3)*Adv_phi(2,3) + A1(:,3,3)*Adv_phi(3,3) + A1(:,4,3)*Adv_phi(4,3))*AR
        AA2(1) = (Adv_phi(3,1) + Adv_phi(3,2) + Adv_phi(3,3))*AR
        AA2(2:4) = (A2(:,1,1)*Adv_phi(1,1) + A2(:,2,1)*Adv_phi(2,1) + A2(:,3,1)*Adv_phi(3,1) + A2(:,4,1)*Adv_phi(4,1) + &
				   A2(:,1,2)*Adv_phi(1,2) + A2(:,2,2)*Adv_phi(2,2) + A2(:,3,2)*Adv_phi(3,2) + A2(:,4,2)*Adv_phi(4,2) + &
				   A2(:,1,3)*Adv_phi(1,3) + A2(:,2,3)*Adv_phi(2,3) + A2(:,3,3)*Adv_phi(3,3) + A2(:,4,3)*Adv_phi(4,3))*AR

        rhs_tmp(:,1) = (Nx(1)*AA1 + Ny(1)*AA2)*tau + &
       		(Adv(:,1)*sp(1,1) + Adv(:,2)*sp(1,2) + Adv(:,3)*sp(1,3))*AR + 3*(Nx(1)*Ux + Ny(1)*Uy)*choq(1)
        rhs_tmp(:,2) = (Nx(2)*AA1 + Ny(2)*AA2)*tau + &
            (Adv(:,1)*sp(2,1) + Adv(:,2)*sp(2,2) + Adv(:,3)*sp(2,3))*AR + 3*(Nx(2)*Ux + Ny(2)*Uy)*choq(2)
        rhs_tmp(:,3) = (Nx(3)*AA1 + Ny(3)*AA2)*tau + &
            (Adv(:,1)*sp(3,1) + Adv(:,2)*sp(3,2) + Adv(:,3)*sp(3,3))*AR + 3*(Nx(3)*Ux + Ny(3)*Uy)*choq(3)

		do i = 1,3
		!$OMP ATOMIC
        rhs(1,ipoin(i)) = rhs(1,ipoin(i)) + rhs_tmp(1,i)
		!$OMP ATOMIC
        rhs(2,ipoin(i)) = rhs(2,ipoin(i)) + rhs_tmp(2,i)
		!$OMP ATOMIC
        rhs(3,ipoin(i)) = rhs(3,ipoin(i)) + rhs_tmp(3,i)
		!$OMP ATOMIC
        rhs(4,ipoin(i)) = rhs(4,ipoin(i)) + rhs_tmp(4,i)
		end do

    ENDDO
    !$OMP end do
    !$OMP end PARALLEL
    return
end subroutine calcRHS

!CCCC-----------------------------------!CCCC
!CCCC  ----> Impone la velocidad <----  !CCCC
!CCCC-----------------------------------!CCCC  
subroutine FIXVEL
    !use MALLOCAR
    !use MVARIABFIX
    use MVELOCIDADES
	use MeshData
    implicit real(8) (A-H,O-Z)
   
    do IFIXV = 1,NFIXV
        N1 = IFIXV_NODE(IFIXV)
        VEL_X(N1) = RFIXV_VALUEX(IFIXV)
        VEL_Y(N1) = RFIXV_VALUEY(IFIXV)
    end do
end subroutine FIXVEL

!CCCC--------------------------------------!CCCC
!CCCC     FIJA LOS VALORES DE T Y RHO      !CCCC
!CCCC     EN LOS CONTORNOS                 !CCCC
!CCCC--------------------------------------!CCCC
subroutine FIX(FR,GAMM)
    !use MALLOCAR
    !use MVARIABFIX
    use MVELOCIDADES
    use MVARIABLES
	use MeshData
    implicit real(8) (A-H,O-Z)
    real(8) GAMM(npoin)
   
    do IFIXRHO = 1,NFIXRHO
        INOD = IFIXRHO_NODE(IFIXRHO)
        RHO(INOD) = RFIXRHO_VALUE(IFIXRHO)
    end do
  
    do IFIXT = 1,NFIXT
        INOD = IFIXT_NODE(IFIXT)
        GM = GAMM(INOD)-1.D0
        T(INOD) = RFIXT_VALUE(IFIXT)
        E(INOD) = T(INOD)*FR/GM+.5D0*(VEL_X(INOD)**2.D0+VEL_Y(INOD)**2.D0)
    end do
end subroutine FIX

subroutine RK(DTMIN, NRK, BANDERA, GAMM, dtl)
    !use DATOS_REFINAMIENTO
    !use DATOS_ENTRADA
    !use MALLOCAR
    !use MVARIABFIX
    !use MGEOMETRIA
    !use MFUERZAS
    !use MMOVIMIENTO
	use InputData
    use MNORMALES
    use MVELOCIDADES
    use MVARIABGEN
	use MeshData
    use MVARIABLES
    use MESTABILIZACION
    use TIMERS
    implicit real(8)(A-H,O-Z)
	integer BANDERA
	real(8) GAMM(npoin)
	real(8) dtl(nelem)

    do IRK = 1,NRK
        RK_FACT = 1.D0/(NRK+1-IRK)
     
        !CCCC  ----> SOLO CALCULO UNA VEZ EL TERMINO DE ESTABILIZACION
        if(IRK.EQ.1)THEN
          
			call system_clock(start_t,rate)
        	call CUARTO_ORDEN(U1,UN,FR,gamm)
			call system_clock(end_t)
			cuarto_t = cuarto_t + real(end_t-start_t)/real(rate)
         
			call system_clock(start_t,rate)
            call ESTAB(U,T,GAMA,FR,RMU &
                ,DTMIN,RHO_inf,T_inf,U_inf,V_inf,GAMM)
			call system_clock(end_t)
			estab_t= estab_t+ real(end_t-start_t)/real(rate)
        
        end if
     
        !CCCC ------> CALCULO EL PRESSURE-SWICHT
        !CCCC --------------------------------------------------------
        !call NEWPRES(npoin,NPOS,NN1,NN2,dabs(RHO),PR)
        !call NEWPRES(npoin,NPOS,NN1,NN2,dabs(RHO*RMACH),PM)
        !call NEWPRES(npoin,NPOS,NN1,NN2,dabs(rmach),PS)
        PS = 1.0d0
     
        RHS = 0.D0
    
		call system_clock(start_t,rate)
        call calcRHS(U,UN,RHS,P,FR,RMU,FK,FCV,T_inf,dtl, gamm)
		call system_clock(end_t)
		calcrhs_t = calcrhs_t + real(end_t - start_t)/real(rate)
     
        !CCCC  ----> CALCULO DE LOS TERMINOS FUENTES
		call system_clock(start_t,rate)
        call FUENTE(dtl)
		call system_clock(end_t)
		fuente_t = fuente_t + real(end_t - start_t)/real(rate)
     
        !CCCC ----> INTEGRADOR TEMPORAL
        !$OMP PARALLEL do
        do INOD = 1,npoin
            U1(:,INOD) = U(:,INOD)-RK_FACT/M(INOD)*RHS(:,INOD)
        end do
        !$OMP end PARALLEL do

        !CCCC----------------------------------------------
        !CCCC----> PASA A LA VARIABLE PRIMARIA PARA APLICAR
        !CCCC----> LAS CONDICIONES DE CONTORNO
        !CCCC----------------------------------------------
        if(NGAS.EQ.1) GO TO 112
        do INOD = 1,npoin
            RHO(INOD) = U1(1,INOD)
            VEL_X(INOD) = U1(2,INOD)/RHO(INOD)
            VEL_Y(INOD) = U1(3,INOD)/RHO(INOD)
            E(INOD) = U1(4,INOD)/RHO(INOD)
            VEL2 = (VEL_X(INOD)**2.D0+VEL_Y(INOD)**2.D0)
            P(INOD) = RHO(INOD)*(GAMM(INOD)-1.D0)*(E(INOD)-.5D0*VEL2)
            T(INOD) = P(INOD)/(RHO(INOD)*FR)
            RMACH(INOD) = DSQRT(VEL2)/T(INOD)
        end do
112 CONTINUE
     
    !CCCC----> CASO PARA AIRE EN EQUILIBRIO
    !CCCC----------------------------------------------
    if(NGAS.NE.0)THEN
        do INOD = 1,npoin
            RHO(INOD) = U1(1,INOD)
            E(INOD) = U1(4,INOD)/RHO(INOD)
            VEL_X(INOD) = U1(2,INOD)/RHO(INOD)
            VEL_Y(INOD) = U1(3,INOD)/RHO(INOD)
            VEL2 = VEL_X(INOD)**2.D0+VEL_Y(INOD)**2.D0
            P(INOD) = RHO(INOD)*(GAMM(INOD)-1.D0)*(E(INOD)-.5D0*VEL2)
            T(INOD) = P(INOD)/(RHO(INOD)*FR)
            if(IRK.EQ.NRK)THEN
                call TGAS(E(INOD)-.5D0*VEL2,RHO(INOD),PGAS,AGAS,TGASi,GAMI)
                GAMM(INOD) = GAMI
                P(INOD) = PGAS
                T(INOD) = TGASi
                RMACH(INOD) = DSQRT(VEL2)/AGAS
            end if
        end do
    end if
     
    !!$        !CCCC----------------------------------------------
    !!$        !CCCC----> PASA INFORMACION DE PERIODICIDAD
    !!$
    !!$        if(NMASTER.NE.0)THEN
    !!$           do IMASTER = 1,NMASTER
    !!$              N1 = IPER_MASTER(IMASTER)
    !!$              N2 = IPER_AUX(IMASTER)
    !!$
    !!$              RRHO = (RHO(N1)+RHO(N2))/2.D0
    !!$              RHO(N1) = RRHO ; RHO(N2) = RRHO
    !!$
    !!$              RVELX = (VEL_X(N1)+VEL_X(N2))/2.D0
    !!$              VEL_X(N1) = RVELX ; VEL_X(N2) = RVELX
    !!$
    !!$              RVELY = (VEL_Y(N1)+VEL_Y(N2))/2.D0
    !!$              VEL_Y(N1) = RVELY ; VEL_Y(N2) = RVELY
    !!$
    !!$              RE = (E(N1)+E(N2))/2.D0
    !!$              E(N1) = RE ; E(N2) = RE
    !!$
    !!$              RP = (P(N1)+P(N2))/2.D0
    !!$              P(N1) = RP ; P(N2) = RP
    !!$
    !!$              RT = (T(N1)+T(N2))/2.D0
    !!$              T(N1) = RT ; T(N2) = RT
    !!$           end do
    !!$        end if

    !CCCC---------------------------------------CCCC
    !CCCC  ----> MOVIMIENTO DE MALLA     <----  CCCC
    !CCCC---------------------------------------CCCC

    !DELTAX = -dtmin*500.d0
    !do i = 1,NMOVE
    !   POS_AUX(i) = DELTAX
    !end do
    ! !NODOS SIN MOVIMIENTO
    !do i = NMOVE+1,NFIX_MOVE+NMOVE
    !   POS_AUX(i) = 0.D0
    !end do
      
    !B = 0.D0
    !call GRADCONJ2(S,XPOS,B,NN1,NN2,npoin,NPOS &
    !     ,ILAUX,POS_AUX,NNMOVE,RAUX,ADIAG &
    !     ,PK,APK,Z)
             
    !W_X = XPOS/DTMIN
    !X = X+XPOS
    !YPOS = Y
    !!$        IMOOTH = IMOOTH+1
    !!$        if(IMOOTH.EQ.IPRINT)THEN
    !!$           IMOOTH = 0
    !!$           call SMOOTH_MESH(npoin,nelem,inpoel,XPOS,YPOS &
    !!$                ,SMOOTH_FIX,SMOOTH_SIM)
    !!$           W_X = (X-XPOS)/DTMIN
    !!$           W_Y = (Y-YPOS)/DTMIN
    !!$        end if
    !       X = XPOS
    !       Y = YPOS

    !CCCC---------------------------------------CCCC
    !CCCC  ----> CONDICIONES DE CONTORNO <----  CCCC
    !CCCC---------------------------------------CCCC
     
    	!CCCC----> VELOCIDADES IMPUESTAS
    	!CCCC---------------------------
		call FIXVEL
       
    	!CCCC----> CORRECCION DE LAS VELOCIDADES NORMALES
    	!CCCC--------------------------------------------
    	call NORMALVEL
         
    	!CCCC----> VALORES IMPUESTOS
    	!CCCC-----------------------
    	call FIX(FR,GAMM)
    
		!$OMP PARALLEL do 
		do INOD = 1,npoin
        	U1(1,INOD) = RHO(INOD)
        	U1(2,INOD) = VEL_X(INOD)*RHO(INOD)
        	U1(3,INOD) = VEL_Y(INOD)*RHO(INOD)
        	U1(4,INOD) = E(INOD)*RHO(INOD)
		end do
		!$OMP end PARALLEL do

	end do

	if (BANDERA.EQ.2) THEN
    	RHS3 = RHS
	else if(BANDERA.EQ.3) THEN
    	RHS2 = RHS
	else if(BANDERA.EQ.4) THEN
    	RHS1 = RHS
	end if
end subroutine RK

!CCCC---->       ADAMS BASHFORTH
subroutine ADAMSB(DTMIN, NESTAB, GAMM, dtl)
    !use DATOS_REFINAMIENTO
    !use DATOS_ENTRADA
    !use MGEOMETRIA
	use InputData
    use MVELOCIDADES
    use MVARIABGEN
	use MeshData
    use MVARIABLES
    use MESTABILIZACION
    use TIMERS
	use Mnormales
    implicit real(8) (A-H,O-Z)
	integer NESTAB
	real(8) GAMM(npoin)
	real(8) dtl(nelem)

       !CCCC  ----> SOLO CALCULO UNA VEZ EL TERMINO DE ESTABILIZACION
    if (NESTAB.EQ.4) NESTAB = 1
    if (NESTAB.EQ.2) THEN
		call system_clock(start_t,rate)
        call CUARTO_ORDEN(U1,UN,FR,gamm)
		call system_clock(end_t)
		cuarto_t = cuarto_t + real(end_t-start_t)/real(rate)

        call ESTAB(U,T,GAMA,FR,RMU &
            ,DTMIN,RHO_inf,T_inf,U_inf,V_inf, GAMM)
        
    end if
    NESTAB = NESTAB+1
    !CCCC ------> CALCULO EL PRESSURE-SWICHT
    !CCCC --------------------------------------------------------
       !call NEWPRES(npoin,NPOS,NN1,NN2,dabs(RHO),PR)
       !call NEWPRES(npoin,NPOS,NN1,NN2,dabs(RHO*RMACH),PM)
       !call NEWPRES(npoin,NPOS,NN1,NN2,dabs(rmach),PS)
    PS = 1.d0
        
    RHS = 0.D0
        
	call system_clock(start_t,rate)
    call calcRHS(U,UN,RHS,P,FR,RMU,FK,FCV,T_inf,dtl, gamm)
	call system_clock(end_t)
	calcrhs_t = calcrhs_t + real(end_t - start_t)/real(rate)

    !CCCC  ----> CALCULO DE LOS TERMINOS FUENTES

	call system_clock(start_t,rate)
    call FUENTE(dtl)
	call system_clock(end_t)
	fuente_t = fuente_t + real(end_t - start_t)/real(rate)
       
        !CCCC ----> INTEGRADOR TEMPORAL
    !$OMP PARALLEL do PRIVATE(RL)
    do INOD = 1,npoin
        RL = 24.D0*M(INOD)
        U1(:,INOD) = U(:,INOD)-(55.D0*RHS(:,INOD)-59.D0*RHS1(:,INOD)+37.D0*RHS2(:,INOD)-9.D0*RHS3(:,INOD))/RL
    end do
    !$OMP end PARALLEL do

    RHS1 = RHS
    RHS2 = RHS1
    RHS3 = RHS2

    !CCCC----------------------------------------------
    !CCCC----> PASA A LA VARIABLE PRIMARIA PARA APLICAR
    !CCCC----> LAS CONDICIONES DE CONTORNO
    !CCCC----------------------------------------------
    if(NGAS.EQ.1) GO TO 111
    do INOD = 1,npoin
        RHO(INOD) = U1(1,INOD)
        VEL_X(INOD) = U1(2,INOD)/RHO(INOD)
        VEL_Y(INOD) = U1(3,INOD)/RHO(INOD)
        E(INOD) = U1(4,INOD)/RHO(INOD)
        VEL2 = (VEL_X(INOD)**2.D0+VEL_Y(INOD)**2.D0)
        P(INOD) = RHO(INOD)*(GAMM(INOD)-1.D0)*(E(INOD)-.5D0*VEL2)
        T(INOD) = P(INOD)/(RHO(INOD)*FR)
        RMACH(INOD) = DSQRT(VEL2)/T(INOD)
    end do
111 CONTINUE

    !CCCC----> CASO PARA AIRE EN EQUILIBRIO
    !CCCC----------------------------------------------
    if(NGAS.NE.0)THEN
        do INOD = 1,npoin
            RHO(INOD) = U1(1,INOD)
            E(INOD) = U1(4,INOD)/RHO(INOD)
            VEL_X(INOD) = U1(2,INOD)/RHO(INOD)
            VEL_Y(INOD) = U1(3,INOD)/RHO(INOD)
            VEL2 = VEL_X(INOD)**2.D0+VEL_Y(INOD)**2.D0
            P(INOD) = RHO(INOD)*(GAMM(INOD)-1.D0)*(E(INOD)-.5D0*VEL2)
            T(INOD) = P(INOD)/(RHO(INOD)*FR)
              
            call TGAS(E(INOD)-.5D0*VEL2,RHO(INOD),PGAS,AGAS,TGASi,GAMI)
            GAMM(INOD) = GAMI
            P(INOD) = PGAS
            T(INOD) = TGASi
            RMACH(INOD) = DSQRT(VEL2)/AGAS
              
        end do
    end if

    !!$        !CCCC----------------------------------------------
    !!$        !CCCC----> PASA INFORMACION DE PERIODICIDAD
    !!$
    !!$        if(NMASTER.NE.0)THEN
    !!$           do IMASTER = 1,NMASTER
    !!$              N1 = IPER_MASTER(IMASTER)
    !!$              N2 = IPER_AUX(IMASTER)
    !!$
    !!$              RRHO = (RHO(N1)+RHO(N2))/2.D0
    !!$              RHO(N1) = RRHO ; RHO(N2) = RRHO
    !!$
    !!$              RVELX = (VEL_X(N1)+VEL_X(N2))/2.D0
    !!$              VEL_X(N1) = RVELX ; VEL_X(N2) = RVELX
    !!$
    !!$              RVELY = (VEL_Y(N1)+VEL_Y(N2))/2.D0
    !!$              VEL_Y(N1) = RVELY ; VEL_Y(N2) = RVELY
    !!$
    !!$              RE = (E(N1)+E(N2))/2.D0
    !!$              E(N1) = RE ; E(N2) = RE
    !!$
    !!$              RP = (P(N1)+P(N2))/2.D0
    !!$              P(N1) = RP ; P(N2) = RP
    !!$
    !!$              RT = (T(N1)+T(N2))/2.D0
    !!$              T(N1) = RT ; T(N2) = RT
    !!$           end do
    !!$        end if

    !CCCC---------------------------------------CCCC
    !CCCC  ----> MOVIMIENTO DE MALLA     <----  CCCC
    !CCCC---------------------------------------CCCC
 
     !DELTAX = -dtmin*500.d0
     !do i = 1,NMOVE
     !  POS_AUX(i) = DELTAX
     !end do
     !
    !NODOS SIN MOVIMIENTO
    ! do i = NMOVE+1,NFIX_MOVE+NMOVE
    !    POS_AUX(i) = 0.D0
    ! end do
    ! B = 0.D0

    !call GRADCONJ2(S,XPOS,B,NN1,NN2,npoin,NPOS &
    !     ,ILAUX,POS_AUX,NNMOVE,RAUX,ADIAG &
    !      ,PK,APK,Z)

    !W_X = XPOS/DTMIN
    !X = X+XPOS
     !YPOS = Y
    !        IMOOTH = IMOOTH+1
    !        if(IMOOTH.EQ.IPRINT)THEN
    !           IMOOTH = 0
    !           call SMOOTH_MESH(npoin,nelem,inpoel,XPOS,YPOS &
    !                ,SMOOTH_FIX,SMOOTH_SIM)
    !           W_X = (X-XPOS)/DTMIN
    !           W_Y = (Y-YPOS)/DTMIN
    !        end if
    !X = XPOS
    !Y = YPOS

    !CCCC---------------------------------------CCCC
    !CCCC  ----> CONDICIONES DE CONTORNO <----  CCCC
    !CCCC---------------------------------------CCCC
        
    !CCCC----> VELOCIDADES IMPUESTAS
    !CCCC---------------------------
    call FIXVEL
       
    !CCCC----> CORRECCION DE LAS VELOCIDADES NORMALES
    !CCCC--------------------------------------------
    call NORMALVEL
         
    !CCCC----> VALORES IMPUESTOS
    !CCCC-----------------------
    call FIX(FR,GAMM)
    
    !$OMP PARALLEL do
    do INOD = 1,npoin
        U1(1,INOD) = RHO(INOD)
        U1(2,INOD) = VEL_X(INOD)*RHO(INOD)
        U1(3,INOD) = VEL_Y(INOD)*RHO(INOD)
        U1(4,INOD) = E(INOD)*RHO(INOD)
    end do
	!OMP end PARALLEL
end subroutine ADAMSB

!CCCC-----------------------------------CCCC
!CCCC  CALCULO DE LOS TERMINOS FUENTES  CCCC
!CCCC-----------------------------------CCCC
subroutine FUENTE(dtl)
    !use MALLOCAR
    !use MGEOMETRIA
	use MeshData
    use MVELOCIDADES
    use MVARIABGEN
    !use DATOS_ENTRADA
	use InputData

    implicit real(8) (A-H, O-Z)
  
    real(8) rhs_tmp(4, 3), Ux(4), Uy(4), Nx(3), Ny(3), dtl(nelem)
	real(8) sp(3, 3)
	real(8) wx(3), wy(3)
    integer ipoin(3)
  
	sp(:, 1) = (/ .5D0, .5D0, 0.D0 /)
	sp(:, 2) = (/ 0.D0, .5D0, .5D0 /)
    sp(:, 3) = (/ .5D0, 0.D0, .5D0 /)
  
    !$OMP PARALLEL &
    !$OMP PRIVATE(ielem, Ux, Uy, AR, wx, wy, rhs_tmp, ipoin, i)

    !$OMP do
    do ielem = 1, nelem
        rhs_tmp = 0.d0
        ipoin = inpoel(:, ielem)
 
        Ux = U(:, ipoin(1))*dNx(1, ielem) + U(:, ipoin(2))*dNx(2, ielem) + U(:, ipoin(3))*dNx(3, ielem)
        Uy = U(:, ipoin(1))*dNy(1, ielem) + U(:, ipoin(2))*dNy(2, ielem) + U(:, ipoin(3))*dNy(3, ielem)

        AR = AREA(ielem)*dtl(ielem)/3.D0
     
		! = = = GAUSS POINT 1 = ==
        wx(1) = sp(1, 1)*W_X(ipoin(1)) + sp(2, 1)*W_X(ipoin(2)) + sp(3, 1)*W_X(ipoin(3))
        wy(1) = sp(1, 1)*W_Y(ipoin(1)) + sp(2, 1)*W_Y(ipoin(2)) + sp(3, 1)*W_Y(ipoin(3))
		! = = = GAUSS POINT 2 = ==
        wx(2) = sp(1, 2)*W_X(ipoin(1)) + sp(2, 2)*W_X(ipoin(2)) + sp(3, 2)*W_X(ipoin(3))
        wy(2) = sp(1, 2)*W_Y(ipoin(1)) + sp(2, 2)*W_Y(ipoin(2)) + sp(3, 2)*W_Y(ipoin(3))
		! = = = GAUSS POINT 3 = ==
        wx(3) = sp(1, 3)*W_X(ipoin(1)) + sp(2, 3)*W_X(ipoin(2)) + sp(3, 3)*W_X(ipoin(3))
        wy(3) = sp(1, 3)*W_Y(ipoin(1)) + sp(2, 3)*W_Y(ipoin(2)) + sp(3, 3)*W_Y(ipoin(3))

        rhs_tmp(:, 1) = -AR*(sp(1, 1)*(Ux*wx(1) + Uy*wy(1)) + sp(1, 2)*(Ux*wx(2) + Uy*wy(2)) + sp(1, 3)*(Ux*wx(3) + Uy*wy(3)))
        rhs_tmp(:, 2) = -AR*(sp(2, 1)*(Ux*wx(1) + Uy*wy(1)) + sp(2, 2)*(Ux*wx(2) + Uy*wy(2)) + sp(2, 3)*(Ux*wx(3) + Uy*wy(3)))
        rhs_tmp(:, 3) = -AR*(sp(3, 1)*(Ux*wx(1) + Uy*wy(1)) + sp(3, 2)*(Ux*wx(2) + Uy*wy(2)) + sp(3, 3)*(Ux*wx(3) + Uy*wy(3)))

		do i = 1, 3
		!$OMP ATOMIC
        rhs(1, ipoin(i)) = rhs(1, ipoin(i)) + rhs_tmp(1, i)
		!$OMP ATOMIC
        rhs(2, ipoin(i)) = rhs(2, ipoin(i)) + rhs_tmp(2, i)
		!$OMP ATOMIC
        rhs(3, ipoin(i)) = rhs(3, ipoin(i)) + rhs_tmp(3, i)
		!$OMP ATOMIC
        rhs(4, ipoin(i)) = rhs(4, ipoin(i)) + rhs_tmp(4, i)
		end do

    end do
    !$OMP end do
    !$OMP end PARALLEL
end subroutine FUENTE
