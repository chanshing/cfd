module gcl_mod
	implicit none
	private
	real*8, dimension(:), allocatable :: tot1, tot2
	real*8, dimension(:), allocatable :: area_old, W_x_old, W_y_old
	public :: main, putW, putArea
	contains 
		subroutine main(M, W_x, W_y, dNx, dNy, area, inpoel, dt)
			real*8, intent(inout), dimension(:) :: M
			real*8, intent(in), dimension(:,:) :: dNx, dNy
			real*8, intent(in), dimension(:) :: W_x, W_y, area
			integer, intent(in), dimension(:,:) :: inpoel
			real*8, intent(in) :: dt
			!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			real*8 :: divW, divW_old
			real*8, dimension(3) :: WW_x, WW_y, WW_x_old, WW_y_old
			integer :: ielem, nelem, npoin
			!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			nelem = size(inpoel,2)
			npoin = size(M)
			if(.not.allocated(tot1)) allocate(tot1(npoin))
			if(.not.allocated(tot2)) allocate(tot2(npoin))
			if(.not.allocated(area_old) .or.&
				.not.allocated(W_x_old) .or.&
				.not.allocated(W_y_old)) stop 'Faltan valores (GCL)'
			!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			tot1 = 0.d0
			tot2 = 0.d0
			do ielem = 1, nelem
			! v_x(:) = U(2,inpoel(:,ielem))/U(1,inpoel(:,ielem))
			! v_y(:) = U(3,inpoel(:,ielem))/U(1,inpoel(:,ielem))
			WW_x_old(:) = W_x_old(inpoel(:,ielem))
			WW_y_old(:) = W_y_old(inpoel(:,ielem))
			WW_x(:) = W_x(inpoel(:,ielem))
			WW_y(:) = W_y(inpoel(:,ielem))
			! divV = dNx(1,ielem)*v_x(1) + dNx(2,ielem)*v_x(2) + dNx(3,ielem)*v_x(3)&
			! 	+ dNy(1,ielem)*v_y(1) + dNy(2,ielem)*v_y(2) + dNy(3,ielem)*v_y(3)
			divW = dNx(1,ielem)*WW_x(1) + dNx(2,ielem)*WW_x(2) + dNx(3,ielem)*WW_x(3)&
				+ dNy(1,ielem)*WW_x(1) + dNy(2,ielem)*WW_x(2) + dNy(3,ielem)*WW_x(3)
			divW_old = dNx(1,ielem)*WW_x_old(1) + dNx(2,ielem)*WW_x_old(2) + dNx(3,ielem)*WW_x_old(3)&
				+ dNy(1,ielem)*WW_x_old(1) + dNy(2,ielem)*WW_x_old(2) + dNy(3,ielem)*WW_x_old(3)
			tot1(inpoel(:,ielem)) = tot1(inpoel(:,ielem)) + divW*area(ielem)/3.d0
			tot2(inpoel(:,ielem)) = tot2(inpoel(:,ielem)) + divW_old*area_old(ielem)/3.d0
			end do
			M = M + dt*(tot1 + tot2)/2.d0
		end subroutine
		subroutine putW(W_x, W_y)
			real*8, dimension(:), intent(in) :: W_x, W_y
			!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			if(.not.allocated(W_x_old)) allocate(W_x_old(size(W_x)))
			if(.not.allocated(W_y_old)) allocate(W_y_old(size(W_y)))
			!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			W_x_old = W_x
			W_y_old = W_y
		end subroutine
		subroutine putArea(area)
			real*8, dimension(:), intent(in) :: area
			!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			if(.not.allocated(area_old)) allocate(area_old(size(area)))
			!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			area_old = area
		end subroutine
end module
