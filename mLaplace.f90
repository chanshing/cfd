module Mlaplace
	integer, dimension(:), allocatable :: lap_idx, lap_rowptr
	real(8), dimension(:), allocatable :: lap_sparse, lap_diag
	real*8, parameter :: TWOSQRT3 = 3.46410161513775d0
	private initialize, mu 
contains 
subroutine laplace(inpoel, area, dNx, dNy, nelem, npoin)
!------------------------------------------------------
!	Ensambla el vector sparse del laplaciano, lap_sparse
!------------------------------------------------------
	use PointNeighbor, only: esup1, esup2
	use MeshData, only: X, Y
	implicit none
	integer, intent(in) :: nelem, npoin
	integer, intent(in) :: inpoel(3,nelem)
    real*8, intent(in) :: dNx(3,nelem), dNy(3,nelem), area(nelem)
	!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	real*8 :: a, q, X3(3), Y3(3)
	integer :: ipoin, jpoin, kpoin, i, j, k, iesup, ielem
	logical, save :: isFirstCall = .true.
	!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if(isFirstCall) then
		call initialize(inpoel, nelem, npoin)
		isFirstCall = .false.
	end if
	!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	!$OMP PARALLEL DO PRIVATE(ipoin)
	do ipoin = 1, lap_rowptr(npoin + 1)
		lap_sparse(ipoin) = 0.d0
	end do
	!$OMP END PARALLEL DO

	!$OMP PARALLEL DO PRIVATE(ipoin, iesup, ielem, a, i, jpoin, j, kpoin, k)
	do ipoin = 1, npoin
		do iesup = esup2(ipoin) + 1, esup2(ipoin + 1)
			ielem = esup1(iesup); a = area(ielem)
			X3 = X(inpoel(:,ielem))
			Y3 = Y(inpoel(:,ielem))
			q = 1/(mu(X3,Y3)**2)
			do i = 1, 3
				jpoin = inpoel(i, ielem)
				if (jpoin == ipoin) then !Encuentra el nodo en cuestion
					do j = 1, 3
						kpoin = inpoel(j, ielem)
						do k = lap_rowptr(ipoin) + 1, lap_rowptr(ipoin + 1)
							if (lap_idx(k) == kpoin) then !Encuentra la ubicacion donde guardar
								lap_sparse(k) = lap_sparse(k) + (dNx(i, ielem)*dNx(j, ielem) + dNy(i, ielem)*dNy(j, ielem))*q
							end if
						end do
					end do
				end if
			end do
		end do
		lap_diag(ipoin) = lap_sparse(lap_rowptr(ipoin) + 1)
	end do
	!$OMP END PARALLEL DO
end subroutine laplace

subroutine initialize(inpoel, nelem, npoin)
!----------------------------------------------------------------
!	Inicializa variables necesarias. Arma lap_rowptr & lap_idx
!----------------------------------------------------------------
	use PointNeighbor, only: psup1, psup2, getPsup
	implicit none
	integer, intent(in) :: nelem, npoin
	integer, intent(in) :: inpoel(3, nelem)
	!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	integer :: i, ipoin, ipsup

	if(.not.allocated(psup1).or..not.allocated(psup2)) &
		call getPsup(inpoel, nelem, npoin)

	allocate(lap_sparse(size(psup1)+npoin))
	allocate(lap_idx(size(psup1)+npoin))
	allocate(lap_rowptr(npoin+1))
	allocate(lap_diag(npoin))

	lap_sparse = 0.d0
	lap_rowptr(1) = 0

	do i = 2, npoin+1
		lap_rowptr(i) = psup2(i) + i - 1
	end do

	do ipoin = 1, npoin
		lap_idx(lap_rowptr(ipoin)+1) = ipoin
		ipsup = psup2(ipoin)+1
		do i = lap_rowptr(ipoin)+2,lap_rowptr(ipoin+1)
			lap_idx(i) = psup1(ipsup)
			ipsup = ipsup + 1
		end do
	end do
end subroutine initialize 

pure real*8 function mu(X3, Y3)
	real*8, intent(in) :: X3(:), Y3(:)
	!%%%%%%%%%%%%%%%%%%%%
	real*8 :: l1, l2, l3, l, area
	!%%%%%%%%%%%%%%%%%%%%
	area = X3(2)*Y3(3)+X3(3)*Y3(1)+X3(1)*Y3(2) &
		-(X3(2)*Y3(1)+X3(3)*Y3(2)+X3(1)*Y3(3))
	l1 = (X3(3)-X3(2))**2 + (Y3(3)-Y3(2))**2
	l2 = (X3(1)-X3(3))**2 + (Y3(1)-Y3(3))**2
	l3 = (X3(2)-X3(1))**2 + (Y3(2)-Y3(1))**2
	l = l1 + l2 + l3
	mu = TWOSQRT3*area/l
end function mu
end module Mlaplace
