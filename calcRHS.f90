module calcRHS_mod
	implicit none
contains
	subroutine calcRHS(rhs, U, theta, dNx, dNy, area, shoc, dtl, t_sugn1, t_sugn2, t_sugn3, inpoel, nelem, npoin)
		use InputData, only: Cv => FCV, lambda_ref => FK, mu_ref => FMU, gamma0 => gama, T_inf, cte
		use Mvariables, only: T
		integer, intent(in) :: npoin, nelem, inpoel(3,nelem)
		real*8, intent(inout) :: rhs(4,npoin)
		real*8, intent(in) :: U(4,npoin), theta(4,npoin), dNx(3,nelem), dNy(3,nelem)
		real*8, intent(in), dimension(nelem) :: area, dtl, shoc, t_sugn1, t_sugn2, t_sugn3
		!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		integer :: ielem, ipoi1, ipoi2, ipoi3, i, k
		real*8 :: rho, v1, v2, e, V_sq
		real*8 :: tau1, tau2, tau3, nu
		real*8, dimension(4) :: U_k, theta_k, Ux, Uy
		real*8, dimension(3) :: Nx, Ny
		real*8, dimension(4,3) :: rhs_tmp
		real*8, parameter, dimension(3,3) :: N = &
			reshape((/&
			0.d0, .5d0, .5d0, &
			.5d0, 0.d0, .5d0, &
			.5d0, .5d0, 0.d0 &
			/), (/3,3/))
		REAL*8, dimension(1:4) :: AiUi
		REAL*8, dimension(1:4) :: AiUi_theta
		REAL*8, dimension(1:4) :: A1AiUi_theta
		REAL*8, dimension(1:4) :: A2AiUi_theta
		REAL*8 :: lambda, mu, T_avg 
		REAL*8, dimension(2:4) :: K1jUj
		REAL*8, dimension(2:4) :: K2jUj
		!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		!$OMP PARALLEL DO DEFAULT(PRIVATE) &
		!$OMP SHARED(rhs, U, theta, dNx, dNy, area, shoc, dtl, &
		!$OMP t_sugn1, t_sugn2, t_sugn3, inpoel, nelem, npoin, &
		!$OMP T, Cv, lambda_ref, mu_ref, gamma0, T_inf, cte, N)
		do ielem = 1, nelem
			rhs_tmp(:,:) = 0.d0
			ipoi1 = inpoel(1,ielem)
			ipoi2 = inpoel(2,ielem)
			ipoi3 = inpoel(3,ielem)
			Nx(1:3) = dNx(1:3,ielem)
			Ny(1:3) = dNy(1:3,ielem)
			Ux(1:4) = U(1:4,ipoi1)*Nx(1) + U(1:4,ipoi2)*Nx(2) + U(1:4,ipoi3)*Nx(3)
			Uy(1:4) = U(1:4,ipoi1)*Ny(1) + U(1:4,ipoi2)*Ny(2) + U(1:4,ipoi3)*Ny(3)
			tau1 = t_sugn1(ielem)
			tau2 = t_sugn2(ielem)
			tau3 = t_sugn3(ielem)
			nu = shoc(ielem)*cte
			T_avg = (T(ipoi1) + T(ipoi2) + T(ipoi3))/3.d0
			mu = mu_ref*(T_avg/T_inf)**1.5d0*(T_inf + 110)/(T_avg + 110)
			lambda = lambda_ref*(T_avg/T_inf)**1.5d0*(T_inf + 194)/(T_avg + 194)
				do k = 1, 3
					U_k(1:4) = &
						N(1,k)*U(1:4,ipoi1) +&
						N(2,k)*U(1:4,ipoi2) +&
						N(3,k)*U(1:4,ipoi3)
					theta_k(1:4) = &
						N(1,k)*theta(1:4,ipoi1) +&
						N(2,k)*theta(1:4,ipoi2) +&
						N(3,k)*theta(1:4,ipoi3)
					rho = U_k(1)
					v1 = U_k(2)/rho
					v2 = U_k(3)/rho
					e = U_k(4)/rho
					V_sq = v1**2 + v2**2
					!******************************************************************************
					!*                    Code generated with sympy 0.7.5-git                     *
					!*                                                                            *
					!*              See http://www.sympy.org/ for more information.               *
					!*                                                                            *
					!*                       This file is part of 'ns2DComp'                       *
					!******************************************************************************
					AiUi(1) = Ux(2) + Uy(3)
					AiUi(2) = (1.0d0/2.0d0)*Ux(1)*(V_sq*(gamma0 - 1) - 2*v1**2) - Ux(2)*v1*( &
						gamma0 - 3) - Ux(3)*v2*(gamma0 - 1) + Ux(4)*(gamma0 - 1) - Uy(1)*v1*v2 + &
						Uy(2)*v2 + Uy(3)*v1
					AiUi(3) = -Ux(1)*v1*v2 + Ux(2)*v2 + Ux(3)*v1 + (1.0d0/2.0d0)*Uy(1)*(V_sq*( &
						gamma0 - 1) - 2*v2**2) - Uy(2)*v1*(gamma0 - 1) - Uy(3)*v2*(gamma0 - 3) + &
						Uy(4)*(gamma0 - 1)
					AiUi(4) = Ux(1)*v1*(V_sq*(gamma0 - 1) - e*gamma0) - 1.0d0/2.0d0*Ux(2)*(V_sq &
						*(gamma0 - 1) - 2*e*gamma0 + 2*v1**2*(gamma0 - 1)) - Ux(3)*v1*v2*( &
						gamma0 - 1) + Ux(4)*gamma0*v1 + Uy(1)*v2*(V_sq*(gamma0 - 1) - e*gamma0) - &
						Uy(2)*v1*v2*(gamma0 - 1) - 1.0d0/2.0d0*Uy(3)*(V_sq*(gamma0 - 1) - 2*e* &
						gamma0 + 2*v2**2*(gamma0 - 1)) + Uy(4)*gamma0*v2

					AiUi_theta(1:4) = +theta_k(1:4) + AiUi(1:4)

					A1AiUi_theta(1) = AiUi_theta(2)
					A1AiUi_theta(2) = v1*(-gamma0 + 3)*AiUi_theta(2) - v2*(gamma0 - 1)* &
						AiUi_theta(3) + (gamma0 - 1)*AiUi_theta(4) + ((1.0d0/2.0d0)* &
						V_sq*(gamma0 - 1) - v1**2)*AiUi_theta(1)
					A1AiUi_theta(3) = -v1*v2*AiUi_theta(1) + v1*AiUi_theta(3) + v2* &
						AiUi_theta(2)
					A1AiUi_theta(4) = gamma0*v1*AiUi_theta(4) - v1*v2*(gamma0 - 1)* &
						AiUi_theta(3) + v1*(V_sq*(gamma0 - 1) - e*gamma0)*AiUi_theta(1) &
						+ (-1.0d0/2.0d0*V_sq*(gamma0 - 1) + e*gamma0 - v1**2*(gamma0 - 1 &
						))*AiUi_theta(2)
					A2AiUi_theta(1) = AiUi_theta(3)
					A2AiUi_theta(2) = -v1*v2*AiUi_theta(1) + v1*AiUi_theta(3) + v2* &
						AiUi_theta(2)
					A2AiUi_theta(3) = -v1*(gamma0 - 1)*AiUi_theta(2) + v2*(-gamma0 + 3)* &
						AiUi_theta(3) + (gamma0 - 1)*AiUi_theta(4) + ((1.0d0/2.0d0)* &
						V_sq*(gamma0 - 1) - v2**2)*AiUi_theta(1)
					A2AiUi_theta(4) = gamma0*v2*AiUi_theta(4) - v1*v2*(gamma0 - 1)* &
						AiUi_theta(2) + v2*(V_sq*(gamma0 - 1) - e*gamma0)*AiUi_theta(1) &
						+ (-1.0d0/2.0d0*V_sq*(gamma0 - 1) + e*gamma0 - v2**2*(gamma0 - 1 &
						))*AiUi_theta(3)

					rhs_tmp(1:4,1) = rhs_tmp(1:4,1) + N(1,k)*AiUi(1:4) + &
						tau1*(Nx(1)*A1AiUi_theta(1:4) + Ny(1)*A2AiUi_theta(1:4)) + &
						nu*(Nx(1)*Ux(1:4) + Ny(1)*Uy(1:4))
					rhs_tmp(1:4,2) = rhs_tmp(1:4,2) + N(2,k)*AiUi(1:4) + &
						tau2*(Nx(2)*A1AiUi_theta(1:4) + Ny(2)*A2AiUi_theta(1:4)) + &
						nu*(Nx(2)*Ux(1:4) + Ny(2)*Uy(1:4))
					rhs_tmp(1:4,3) = rhs_tmp(1:4,3) + N(3,k)*AiUi(1:4) + &
						tau3*(Nx(3)*A1AiUi_theta(1:4) + Ny(3)*A2AiUi_theta(1:4)) + &
						nu*(Nx(3)*Ux(1:4) + Ny(3)*Uy(1:4))

					if(mu_ref > tiny(0d0)) then
						! K1jUj(1) = 0
						K1jUj(2) = (2.0d0/3.0d0)*mu*(-2*Ux(1)*v1 + 2*Ux(2) + Uy(1)*v2 - Uy(3))/rho
						K1jUj(3) = mu*(-Ux(1)*v2 + Ux(3) - Uy(1)*v1 + Uy(2))/rho
						K1jUj(4) = (1.0d0/3.0d0)*(Cv*mu*(-Uy(1)*v1*v2 + 3*Uy(2)*v2 - 2*Uy(3)*v1) &
							- Ux(1)*(Cv*mu*(3*V_sq + v1**2) - 3*lambda*(V_sq - e)) + Ux(2)*v1*(4*Cv*mu &
							- 3*lambda) + 3*Ux(3)*v2*(Cv*mu - lambda) + 3*Ux(4)*lambda)/(Cv*rho)

						! K2jUj(1) = 0
						K2jUj(2) = mu*(-Ux(1)*v2 + Ux(3) - Uy(1)*v1 + Uy(2))/rho
						K2jUj(3) = (2.0d0/3.0d0)*mu*(Ux(1)*v1 - Ux(2) - 2*Uy(1)*v2 + 2*Uy(3))/rho
						K2jUj(4) = (1.0d0/3.0d0)*(Cv*mu*(-Ux(1)*v1*v2 - 2*Ux(2)*v2 + 3*Ux(3)*v1) &
							- Uy(1)*(Cv*mu*(3*V_sq + v2**2) - 3*lambda*(V_sq - e)) + 3*Uy(2)*v1*(Cv*mu &
							- lambda) + Uy(3)*v2*(4*Cv*mu - 3*lambda) + 3*Uy(4)*lambda)/(Cv*rho)

						rhs_tmp(2:4,1) = rhs_tmp(2:4,1) + (Nx(1)*K1jUj(2:4) + Ny(1)*K2jUj(2:4))
						rhs_tmp(2:4,2) = rhs_tmp(2:4,2) + (Nx(2)*K1jUj(2:4) + Ny(2)*K2jUj(2:4))
						rhs_tmp(2:4,3) = rhs_tmp(2:4,3) + (Nx(3)*K1jUj(2:4) + Ny(3)*K2jUj(2:4))
					end if
					!******************************************************************************
				end do

				rhs_tmp(:,:) = rhs_tmp(:,:)*area(ielem)*dtl(ielem)/3.d0

			do i = 1, 4
				!$OMP ATOMIC
				RHS(i,ipoi1) = RHS(i,ipoi1) + rhs_tmp(i,1)
				!$OMP ATOMIC
				RHS(i,ipoi2) = RHS(i,ipoi2) + rhs_tmp(i,2)
				!$OMP ATOMIC
				RHS(i,ipoi3) = RHS(i,ipoi3) + rhs_tmp(i,3)
			end do
		end do
		!$OMP END PARALLEL DO

	end subroutine calcRHS

end module calcRHS_mod
