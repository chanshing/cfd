module Mlaplace
	integer, dimension(:), allocatable :: lap_idx, lap_rowptr
	real(8), dimension(:), allocatable :: lap_sparse, lap_diag
	private initialize 
contains 
subroutine laplace(inpoel, area, dNx, dNy, nelem, npoin)
!------------------------------------------------------
!	Ensambla el vector sparse del laplaciano, lap_sparse
!------------------------------------------------------
	use PointNeighbor, only: esup1, esup2
	implicit none
	integer nelem, npoin
	integer inpoel(3,nelem)
    real(8) dNx(3,nelem), dNy(3,nelem), area(nelem), a
	integer ipoin, iesup, ielem, i, jpoin, j, kpoin, k

	if(.not.allocated(lap_sparse)) call initialize(inpoel, nelem, npoin)

	!$omp parallel do private(ipoin)
	do ipoin = 1, lap_rowptr(npoin + 1)
		lap_sparse(ipoin) = 0.d0
	end do
	!$omp end parallel do

	!$OMP PARALLEL DO PRIVATE(ipoin, iesup, ielem, a, i, jpoin, j, kpoin, k)
	do ipoin = 1, npoin
		do iesup = esup2(ipoin) + 1, esup2(ipoin + 1)
			ielem = esup1(iesup); a = area(ielem)
			do i = 1, 3
				jpoin = inpoel(i, ielem)
				if (jpoin == ipoin) then !Encuentra el nodo en cuestion
					do j = 1, 3
						kpoin = inpoel(j, ielem)
						do k = lap_rowptr(ipoin) + 1, lap_rowptr(ipoin + 1)
							if (lap_idx(k) == kpoin) then !Encuentra la ubicacion donde guardar
								lap_sparse(k) = lap_sparse(k) + (dNx(i, ielem)*dNx(j, ielem) + dNy(i, ielem)*dNy(j, ielem))/a**2
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
	integer nelem, npoin
	integer inpoel(3, nelem)
	integer i, ipoin, i_psup1

	if(.not.allocated(psup1)) call getPsup(inpoel, nelem, npoin)

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
		i_psup1 = psup2(ipoin)+1
		do i = lap_rowptr(ipoin)+2,lap_rowptr(ipoin+1)
			lap_idx(i) = psup1(i_psup1)
			i_psup1 = i_psup1 + 1
		end do
	end do
end subroutine initialize 
end module Mlaplace
