module MbiconjGrad
	real(8), allocatable, dimension(:) :: y_BiCG, p_BiCG, r_BiCG, z_BiCG
end module

subroutine allocBiCG(npoin)
	use MbiconjGrad
	allocate(y_BiCG(npoin), p_BiCG(npoin), r_BiCG(npoin), z_BiCG(npoin))
end subroutine

subroutine biconjGrad(spMtx, spIdx, spRowptr, diagMtx, x, b, x_fix, x_fixIdx, nMtxSize, npoin, nfix)
!-------------------------------------------------------------------------------------------------------------
!						 Calcula x dado A.x = b, mediante gradiente biconjugado
!-------------------------------------------------------------------------------------------------------------
	use TIMERS
	use MbiconjGrad
	integer spIdx(nMtxSize), spRowptr(npoin+1)
	integer x_fixIdx(nfix)
	real(8) spMtx(nMtxSize), x(npoin), b(npoin), diagMtx(npoin)
	real(8)	x_fix(nfix)
	real(8) tol, err_old, err_new, alfa, beta

	k = 0
	tol = 1.d-10
	
	x(x_fixIdx(:)) = x_fix(:)
	call system_clock(sub_start, sub_rate)
	call SpMV(spMtx, spIdx, spRowptr, x, y_BiCG, npoin, nMtxSize)
	y_BiCG(x_fixIdx(:)) = x(x_fixIdx(:))*1.d30
	call system_clock(sub_end)
	spmv_t = spmv_t + real(sub_end - sub_start)/real(sub_rate)
	r_BiCG = b - y_BiCG
	r_BiCG(x_fixIdx(:)) = 0.d0

	if(sum(r_BiCG*r_BiCG) < tol) return

	p_BiCG = r_BiCG/diagMtx
	err_new = sum(r_BiCG*p_BiCG)
	call system_clock(sub_start, sub_rate)
	call SpMV(spMtx, spIdx, spRowptr, p_BiCG, y_BiCG, npoin, nMtxSize)
	y_BiCG(x_fixIdx(:)) = p_BiCG(x_fixIdx(:))*1.d30
	call system_clock(sub_end)
	spmv_t = spmv_t + real(sub_end - sub_start)/real(sub_rate)
	alfa = err_new/sum(p_BiCG*y_BiCG)
	x = x + alfa*p_BiCG
	err_old = err_new

	do while(dabs(err_old) > tol .and. k < 1000)
		k = k + 1
		r_BiCG = r_BiCG - alfa*y_BiCG
		z_BiCG = r_BiCG/diagMtx
		err_new = sum(r_BiCG*z_BiCG)
		beta = err_new/err_old
		p_BiCG = z_BiCG + beta*p_BiCG
		call system_clock(sub_start, sub_rate)
		call SpMV(spMtx, spIdx, spRowptr, p_BiCG, y_BiCG, npoin, nMtxSize)
		y_BiCG(x_fixIdx(:)) = p_BiCG(x_fixIdx(:))*1.d30
		call system_clock(sub_end)
		spmv_t = spmv_t + real(sub_end - sub_start)/real(sub_rate)
		alfa = err_new/sum(p_BiCG*y_BiCG)
		x = x + alfa*p_BiCG
		err_old = err_new
	end do

end subroutine
