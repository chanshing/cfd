module PointNeighbor
	integer, dimension(:), allocatable :: esup1, esup2 
	integer, dimension(:), allocatable :: psup1, psup2, lpoin
contains
subroutine getEsup(inpoel,nelem,npoin)
	implicit none
	integer nelem, npoin
	integer inpoel(3,nelem) 
	integer ielem, ipoi1, ipoin, istor, i
	allocate(esup2(npoin+1))
	esup2 = 0

	!contar la cantidad de elementos vecinos de cada nodo (*histogram pattern*)
	do ielem=1,nelem
		do i=1,3
			ipoi1 = inpoel(i,ielem) + 1
			esup2(ipoi1) = esup2(ipoi1) + 1
		end do
	end do
	
	!reshuffle para obtener esup2 (*scan pattern*)
	do ipoin = 2,npoin+1
		esup2(ipoin) = esup2(ipoin) + esup2(ipoin-1)
	end do

	allocate(esup1(esup2(npoin+1))) 
		
	!obtener el array esup1 de los elementos vecinos
	do ielem=1,nelem
		do i=1,3
			ipoin = inpoel(i,ielem)
			istor = esup2(ipoin) + 1 
			esup2(ipoin) = istor     
			esup1(istor) = ielem
		end do
	end do

	!restaurar esup2
	do ipoin=npoin+1,2,-1
		esup2(ipoin) = esup2(ipoin-1)
	end do
	esup2(1) = 0
end subroutine getEsup

subroutine getPsup(inpoel, nelem, npoin)
	implicit none
	integer nelem, npoin
	integer inpoel(3,nelem)
	integer ielem, jpoin, ipoin, istor, iesup, i

	if(.not.allocated(esup1)) call getEsup(inpoel, nelem, npoin)

	allocate(psup2(npoin + 1))
	allocate(lpoin(npoin))
	lpoin = 0; psup2(1) = 0; istor = 0

	!calcular total de nodos vecinos
	do ipoin=1,npoin
		do iesup=esup2(ipoin)+1, esup2(ipoin+1)
			ielem=esup1(iesup)
			do i=1,3
				jpoin=inpoel(i,ielem)
				if(jpoin /= ipoin .and. lpoin(jpoin) /= ipoin) then
					istor = istor + 1
					lpoin(jpoin) = ipoin
				end if
			enddo
		end do
		psup2(ipoin+1) = istor
	end do

	allocate(psup1(istor))

	!obtener array psup1
	lpoin=0; istor=0
	do ipoin=1,npoin
		do iesup=esup2(ipoin)+1, esup2(ipoin+1)
			ielem=esup1(iesup)
			do i=1,3
				jpoin=inpoel(i,ielem)
				if(jpoin /= ipoin .and. lpoin(jpoin) /= ipoin) then
					istor = istor+1
					psup1(istor) = jpoin
					lpoin(jpoin) = ipoin
				end if
			enddo
		end do
	end do

	deallocate(lpoin)
end subroutine getPsup
end module PointNeighbor
