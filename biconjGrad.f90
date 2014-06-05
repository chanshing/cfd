module BiconjGrad
	real(8), allocatable, dimension(:) :: y, p, r, z

	private

	public biCG
contains 
subroutine biCG(spMtx, spIdx, spRowptr, diagMtx, x, b, x_fix, x_fixIdx, npoin, nfix)
!--------------------------------------------------------------------------------------------------
!					 Calcula x dado A.x = b, mediante gradiente biconjugado
!--------------------------------------------------------------------------------------------------
	implicit none
	integer npoin, nfix, k
	integer spRowptr(npoin + 1), spIdx(spRowptr(npoin + 1))
	integer x_fixIdx(nfix)
	real(8) spMtx(spRowptr(npoin + 1)), x(npoin), b(npoin), diagMtx(npoin)
	real(8)	x_fix(nfix)
	real(8) tol, err_old, err_new, alfa, beta, py

	if(.not.allocated(y)) allocate(y(npoin))
	if(.not.allocated(p)) allocate(p(npoin))
	if(.not.allocated(r)) allocate(r(npoin))
	if(.not.allocated(z)) allocate(z(npoin))

	k = 0
	tol = 1.d-10
	
	call copy1(nfix, 1.d0, x_fix, x_fixIdx, npoin, x)
	call SpMV(spMtx, spIdx, spRowptr, x, y, npoin, spRowptr(npoin+1))
!	call mkl_dcsrgemv('n', npoin, spMtx, spRowptr, spIdx, x, y)
	call copy2(npoin, 1.d30, x, nfix, x_fixIdx, y)
	call vecsum(npoin, -1.d0, y, b, r)
	call assign2(npoin, r, nfix, x_fixIdx, 0.d0)

	if(vecdot(npoin, r, r) < tol) return

	call vecdiv(npoin, r, diagMtx, p)
	err_new = vecdot(npoin, r, p)
	call SpMV(spMtx, spIdx, spRowptr, p, y, npoin, spRowptr(npoin+1))
!	call mkl_dcsrgemv('n', npoin, spMtx, spRowptr, spIdx, p, y)
	call copy2(npoin, 1.d30, p, nfix, x_fixIdx, y)
	py = vecdot(npoin, p, y)
	alfa = err_new/py
	call vecsum(npoin, alfa, p, x, x)
	err_old = err_new

	do while(dabs(err_old) > tol .and. k < 1000)
		k = k + 1
		call vecsum(npoin, -alfa, y, r, r)
		call vecdiv(npoin, r, diagMtx, z)
		err_new = vecdot(npoin, r, z)
		beta = err_new/err_old
		call vecsum(npoin, beta, p, z, p)
		call SpMV(spMtx, spIdx, spRowptr, p, y, npoin, spRowptr(npoin+1))
!		call mkl_dcsrgemv('n', npoin, spMtx, spRowptr, spIdx, p, y)
		call copy2(npoin, 1.d30, p, nfix, x_fixIdx, y)
		py = vecdot(npoin, p, y)
		alfa = err_new/py
		call vecsum(npoin, alfa, p, x, x)
		err_old = err_new
	end do
end subroutine biCG

subroutine vecdiv(n, x, y, z)
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
end subroutine vecdiv

subroutine vecsum(n, alfa, x, y, z)
!-----------------------------------
!		Calcula z = alfa*x + y
!-----------------------------------
	implicit none
	integer n, i
	real(8) x(n), y(n), z(n), alfa

	!$OMP PARALLEL DO
	do i = 1, n
		z(i) = alfa*x(i) + y(i)
	end do
	!$OMP END PARALLEL DO
end subroutine vecsum

subroutine assign(n, y, scal)
!-----------------------------
!		y(:) = scal
!-----------------------------
	implicit none
	integer n, i
	real(8) scal, y(n)

	!$OMP PARALLEL DO PRIVATE(i)
	do i = 1, n
		y(i) = scal
	end do
	!$OMP END PARALLEL DO
end subroutine

subroutine assign2(n, y, m, idx, scal)
	implicit none
	integer n, m, i
	real(8) y(n), scal
	integer idx(m)

	!$OMP PARALLEL DO PRIVATE(i)
	do i = 1, m
		y(idx(i)) = scal
	end do
	!$OMP END PARALLEL DO
end subroutine

subroutine copy1(m, alfa, x, idx, n, y)
!-----------------------------------------
!			y(idx(:)) = alfa*x
!-----------------------------------------
	implicit none
	integer n, m, i
	real(8) x(m), y(n), alfa
	integer idx(m)
	!$OMP PARALLEL DO PRIVATE(i)
	do i=1, m
		y(idx(i)) = alfa*x(i)
	end do
	!$OMP END PARALLEL DO
end subroutine

subroutine copy2(n, alfa, x, m, idx, y)
!-----------------------------------------
!			y(idx(:)) = alfa*x(idx(:))
!-----------------------------------------
	implicit none
	integer n, m, i
	real(8) x(n), y(n), alfa
	integer idx(m)

	!$OMP PARALLEL DO PRIVATE(i)
	do i=1, m
		y(idx(i)) = alfa*x(idx(i))
	end do
	!$OMP END PARALLEL DO
end subroutine copy2

real(8) function vecdot(n, x, y)
!-----------------------------------
!	Devuelve producto punto <x,y>
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

	vecdot = res
end function vecdot

subroutine SpMV(spMtx, spIdx, spRowptr, v, y, npoin, npos)
!---------------------------------------------------------
!			Multiplicacion sparse y = M.v
!---------------------------------------------------------
	implicit none
	integer npoin, npos
	integer i, j
	real(8) spMtx(npos), v(npoin), y(npoin), dot
	integer spIdx(npos), spRowptr(npoin+1)

	!$OMP PARALLEL DO PRIVATE(i, j, dot)
	do i = 1, npoin
		dot = 0.d0
		do j = spRowptr(i)+1, spRowptr(i+1)
			dot = dot + spMtx(j)*v(spIdx(j))
		end do
		y(i) = dot
	end do
	!$OMP END PARALLEL DO 
end subroutine
end module BiconjGrad
