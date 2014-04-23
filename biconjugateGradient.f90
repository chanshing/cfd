module MbiconjGrad
	real(8), allocatable, dimension(:), private :: y, p, r, z
contains !allocBiCG & biconjGrad
subroutine allocBiCG(npoin)
	integer npoin
	allocate(y(npoin), p(npoin), r(npoin), z(npoin))
end subroutine allocBiCG

subroutine biconjGrad(spMtx, spIdx, spRowptr, diagMtx, x, b, x_fix, x_fixIdx, nMtxSize, npoin, nfix)
!-------------------------------------------------------------------------------------------------------------
!						 Calcula x dado A.x = b, mediante gradiente biconjugado
!-------------------------------------------------------------------------------------------------------------
	use TIMERS
	implicit none
	integer nMtxSize, npoin, nfix, k
	integer spIdx(nMtxSize), spRowptr(npoin+1)
	integer x_fixIdx(nfix)
	real(8) spMtx(nMtxSize), x(npoin), b(npoin), diagMtx(npoin)
	real(8)	x_fix(nfix)
	real(8) tol, err_old, err_new, alfa, beta, dot

	k = 0
	tol = 1.d-10
	
	x(x_fixIdx(:)) = x_fix(:)
	call SpMV(spMtx, spIdx, spRowptr, x, y, npoin, nMtxSize)
	y(x_fixIdx(:)) = x(x_fixIdx(:))*1.d30
	call sumArray(b, y, r, -1.d0, npoin) 
	r(x_fixIdx(:)) = 0.d0

	if(sum(r*r) < tol) return

	call divArray(r, diagMtx, p, npoin)
	call dotProduct(r, p, err_new, npoin)
	call SpMV(spMtx, spIdx, spRowptr, p, y, npoin, nMtxSize)
	y(x_fixIdx(:)) = p(x_fixIdx(:))*1.d30
	call dotProduct(p, y, dot, npoin)
	alfa = err_new/dot
	call sumArray(x, p, x, alfa, npoin) 
	err_old = err_new

	do while(dabs(err_old) > tol .and. k < 1000)
		k = k + 1
		call sumArray(r, y, r, -alfa, npoin) 
		call divArray(r, diagMtx, z, npoin)
		call dotProduct(r, z, err_new, npoin)
		beta = err_new/err_old
		call sumArray(z, p, p, beta, npoin) 
		call SpMV(spMtx, spIdx, spRowptr, p, y, npoin, nMtxSize)
		y(x_fixIdx(:)) = p(x_fixIdx(:))*1.d30
		call dotProduct(p, y, dot, npoin)
		alfa = err_new/dot
		call sumArray(x, p, x, alfa, npoin) 
		err_old = err_new
	end do
end subroutine biconjGrad
end module MbiconjGrad

subroutine divArray(x, y, z, n)
!---------------------------------
!		Calcula z = x/y
!---------------------------------
	implicit none
	integer n, i
	real(8) x(n), y(n), z(n)

	!$OMP PARALLEL DO
	do i = 1, n
		z(i) = x(i)/y(i)	
	end do
	!$OMP END PARALLEL DO
end subroutine divArray

subroutine sumArray(x, y, z, alfa, n)
!-----------------------------------
!		Calcula z = x + alfa.y
!-----------------------------------
	implicit none
	integer n, i
	real(8) x(n), y(n), z(n), alfa

	!$OMP PARALLEL DO
	do i = 1, n
		z(i) = x(i) + alfa*y(i)
	end do
	!$OMP END PARALLEL DO
end subroutine sumArray

subroutine dotProduct(x, y, res, n)
!-----------------------------------
!		Calcula res = <x,y>
!-----------------------------------
	implicit none
	integer n, i
	real(8) x(n), y(n), res
	res = 0.d0

	!$OMP PARALLEL DO REDUCTION(+:res)
	do i = 1, n
		res = res + x(i)*y(i)
	end do
	!$OMP END PARALLEL DO
end subroutine dotProduct
