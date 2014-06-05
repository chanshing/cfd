module MsparseLaplace
	integer, dimension(:), allocatable :: lap_idx, lap_rowptr
	real(8), dimension(:), allocatable :: lap_sparse, lap_diag
end module

subroutine laplace2(inpoel, area, dNx, dNy, nelem, npoin)
	use MsparseLaplace
	use MelementSurrPoint
	integer nelem, npoin
	integer inpoel(3,nelem)
    real(8) dNx(3,nelem), dNy(3,nelem), area(nelem), a

	lap_sparse = 0.d0
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
end subroutine

subroutine initLapSparse(npoin)
	use MpointSurrPoint
	use MsparseLaplace
	integer npoin
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
end subroutine