module smoothing_mod
	implicit none
	private
	real*8, parameter :: TWOSQRT3 = 3.46410161513775d0
	integer, parameter :: MAX_ELEM_PER_NODE = 20
	integer, parameter :: NITER = 100, MITER = 2, NTRY = 4 
	real*8, parameter :: TOL_METRIC = .9d0
	real*8, parameter :: FACTOR_TOL_DIST = 1d-3
	real*8, parameter :: FACTOR_DELTA = 1d-2
	real*8, parameter :: FACTOR_PLUS = 1.d0
	real*8, parameter :: FACTOR_STEP = 1d-2
	real*8 :: H_MIN_GLOBAL
	public :: smoothing
contains
	subroutine smoothing(X, Y, npoin, inpoel, nelem, smoothable, fix)
		integer, intent(in) :: npoin, nelem, inpoel(3,nelem)
		real*8, intent(inout) :: X(npoin), Y(npoin)
		logical, intent(inout) :: smoothable(npoin)
		logical, intent(in) :: fix(npoin)
		!%%%%%%%%%%%%%%%
		integer :: iter, ipoin, min_idx
		real*8 :: d_max, TOL_DIST
		real*8, dimension(:), allocatable, save :: mu_vec, gx, gy
		logical, save :: isInitialized = .false.
		!%%%%%%%%%%%%%%%
		if(.not.isInitialized) then
			isInitialized = .true.
			allocate(mu_vec(MAX_ELEM_PER_NODE), gx(MAX_ELEM_PER_NODE), gy(MAX_ELEM_PER_NODE))
		end if
		!%%%%%%%%%%%%%%%
		call check_smoothable(npoin, nelem, inpoel, X, Y, smoothable)
		TOL_DIST = FACTOR_TOL_DIST*H_MIN_GLOBAL
		do iter = 1, NITER
		d_max = 0.d0
		do ipoin = 1, npoin
		if(smoothable(ipoin) .and. .not. fix(ipoin)) then
			smoothable(ipoin) = .false. 
			mu_vec = 0.d0
			gx = 0.d0
			gy = 0.d0
			call get_mu_vec(mu_vec, ipoin, npoin, inpoel, nelem, X, Y, min_idx)
			if(mu_vec(min_idx) < TOL_METRIC) then
				call move_ipoin(X, Y, ipoin, npoin, inpoel, nelem, gx, gy, mu_vec, min_idx, d_max, smoothable)
			end if
		end if
		end do
		if(d_max < TOL_DIST) exit
		end do
		print*, '=============SMOOTHING============='
		print*, '# iteraciones:', iter
		print*, 'Min Distortion Metric:', get_mesh_quality(npoin, nelem, X, Y, inpoel)
		print*, 'H_MIN_GLOBAL:', H_MIN_GLOBAL
		print*, '==================================='
	end subroutine smoothing

	subroutine move_ipoin(X, Y, ipoin, npoin, inpoel, nelem, gx, gy, mu_vec, min_idx, d_max, smoothable)
		integer, intent(in) :: ipoin, npoin, nelem, inpoel(3,nelem)
		integer, intent(inout) :: min_idx
		real*8, intent(inout) :: gx(:), gy(:), mu_vec(:)
		real*8, intent(inout) :: X(npoin), Y(npoin)
		real*8, intent(inout) :: d_max
		logical, intent(inout) :: smoothable(npoin)
		!%%%%%%%%%%%%%%%%%%%%
		real*8 :: x_old, y_old, x0, y0, mu_min
		real*8 :: dx, dy, step, d_move
		integer :: i, j, min_idx_new
		logical :: accepted
		!%%%%%%%%%%%%%%%%%%%%
		x0 = X(ipoin)
		y0 = Y(ipoin)
		do i = 1, MITER
		mu_min = mu_vec(min_idx)
		x_old = X(ipoin)
		y_old = Y(ipoin)
		call get_g(gx, gy, ipoin, npoin, inpoel, nelem, X, Y, mu_vec)
		step =  get_step(gx, gy, mu_vec, min_idx, ipoin, npoin, X, Y, nelem, inpoel)
		!TRY NEW POSITION
		accepted = .false.
		do j = 1, NTRY
		dx = step*gx(min_idx)
		dy = step*gy(min_idx)
		X(ipoin) = X(ipoin) + dx
		Y(ipoin) = Y(ipoin) + dy
		call get_mu_vec(mu_vec, ipoin, npoin, inpoel, nelem, X, Y, min_idx_new)
		if(mu_vec(min_idx_new) > mu_min*FACTOR_PLUS) then
			min_idx = min_idx_new
			accepted = .true.
			exit
		else
			step = .5d0*step
		end if
		end do
		!IF FAILED, RESTORE & EXIT
		if(.not.accepted) then
			X(ipoin) = x_old
			Y(ipoin) = y_old
			exit
		end if
		end do
		d_move = (X(ipoin) - x0)**2 + (Y(ipoin) - y0)**2
		!TRACK LARGEST DISTANCE MOVED
		if(d_move > d_max) d_max = d_move
		!IF MOVED, UPDATE LIST
		if(d_move > tiny(1d0)) call update_list(ipoin, smoothable)
	end subroutine move_ipoin

	subroutine get_mu_vec(mu_vec, ipoin, npoin, inpoel, nelem, X, Y, min_idx)
		use PointNeighbor, only: esup1, esup2
		integer, intent(in) :: ipoin, npoin, nelem, inpoel(3,nelem)
		real*8, intent(in) :: X(npoin), Y(npoin)
		real*8, intent(out) :: mu_vec(:)
		integer, optional, intent(out) :: min_idx
		!%%%%%%%%%%%%%%%%%%%%
		integer :: iesup, inode, jpoin
		real*8 :: X_loc(3), Y_loc(3)
		real*8 :: mu_min, mu1
		!%%%%%%%%%%%%%%%%%%%%
		mu_min = 1
		do iesup = esup2(ipoin) + 1, esup2(ipoin + 1)
		do inode = 1, 3
		jpoin = inpoel(inode, esup1(iesup))
		X_loc(inode) = X(jpoin)
		Y_loc(inode) = Y(jpoin)
		end do
		mu1 = mu(X_loc, Y_loc)
		mu_vec(iesup - esup2(ipoin)) = mu1
		if(present(min_idx)) then
			if(mu1 < mu_min) then
				mu_min = mu1
				min_idx = iesup - esup2(ipoin)
			end if
		end if
		end do
	end subroutine get_mu_vec

	subroutine get_g(gx, gy, ipoin, npoin, inpoel, nelem, X, Y, mu_vec)
		integer, intent(in) :: ipoin, npoin, nelem, inpoel(3,nelem)
		real*8, intent(in) :: mu_vec(:)
		real*8, intent(inout) :: X(npoin), Y(npoin) 
		real*8, intent(out) :: gx(:), gy(:)
		!%%%%%%%%%%%%%%%%%%%%
		real*8 :: x_old, y_old
		real*8 :: DELTA
		!%%%%%%%%%%%%%%%%%%%%
		DELTA = H_MIN_GLOBAL*FACTOR_DELTA
		!--- GX ---
		x_old = X(ipoin)
		X(ipoin) = X(ipoin) + DELTA
		call get_mu_vec(gx, ipoin, npoin, inpoel, nelem, X, Y)
		gx = (gx - mu_vec)/DELTA
		X(ipoin) = x_old
		!--- GY ---
		y_old = Y(ipoin)
		Y(ipoin) = Y(ipoin) + DELTA
		call get_mu_vec(gy, ipoin, npoin, inpoel, nelem, X, Y)
		gy = (gy - mu_vec)/DELTA
		Y(ipoin) = y_old
	end subroutine get_g

	function get_step(gx, gy, mu_vec, min_idx, ipoin, npoin, X, Y, nelem, inpoel)
		real*8 :: get_step
		real*8, intent(in) :: gx(:), gy(:), mu_vec(:)
		integer, intent(in) :: min_idx, ipoin, npoin, nelem, inpoel(3,nelem)
		real*8, intent(in) :: X(npoin), Y(npoin)
		!%%%%%%%%%%%%%%%%%%%%
		real*8 :: g2, gg, step1
		real*8 :: gx_min, gy_min, mu_min
		integer :: i
		!%%%%%%%%%%%%%%%%%%%%
		get_step = get_h_min(ipoin, npoin, X, Y, nelem, inpoel)*FACTOR_STEP
		gx_min = gx(min_idx)
		gy_min = gy(min_idx)
		mu_min = mu_vec(min_idx)
		g2 = gx_min**2 + gy_min**2
		do i = 1, size(mu_vec)
		gg = gx_min*gx(i) + gy_min*gy(i)
		if(gg < 0) then
			step1 = (mu_vec(i) - mu_min)/(g2 - gg)
			if(step1 < get_step) get_step = step1
		end if
		end do
	end function get_step

	function get_h_min(ipoin, npoin, X, Y, nelem, inpoel)
		use PointNeighbor, only: esup1, esup2
		real*8 :: get_h_min
		integer, intent(in) :: ipoin, npoin, nelem, inpoel(3,nelem)
		real*8, intent(in) :: X(npoin), Y(npoin)
		!%%%%%%%%%%%%%%%%%%%% 
		integer :: iesup
		real*8 :: X_loc(3), Y_loc(3)
		real*8 :: h1
		!%%%%%%%%%%%%%%%%%%%% 
		get_h_min = 1
		do iesup = esup2(ipoin) + 1, esup2(ipoin + 1)
		X_loc(:) = X(inpoel(:,esup1(iesup)))
		Y_loc(:) = Y(inpoel(:,esup1(iesup)))
		h1 = h(X_loc, Y_loc)
		if(h1 < get_h_min) get_h_min = h1
		end do
	end function

	function mu(X_loc, Y_loc)
		real*8, intent(in) :: X_loc(:), Y_loc(:)
		!%%%%%%%%%%%%%%%%%%%%
		real*8 :: l1, l2, l3, l, area, mu
		!%%%%%%%%%%%%%%%%%%%%
		area = X_loc(2)*Y_loc(3)+X_loc(3)*Y_loc(1)+X_loc(1)*Y_loc(2) &
			-(X_loc(2)*Y_loc(1)+X_loc(3)*Y_loc(2)+X_loc(1)*Y_loc(3))
		l1 = (X_loc(3)-X_loc(2))**2 + (Y_loc(3)-Y_loc(2))**2
		l2 = (X_loc(1)-X_loc(3))**2 + (Y_loc(1)-Y_loc(3))**2
		l3 = (X_loc(2)-X_loc(1))**2 + (Y_loc(2)-Y_loc(1))**2
		l = l1 + l2 + l3
		mu = TWOSQRT3*area/l
	end function mu

	function h(X_loc, Y_loc)
		real*8, intent(in) :: X_loc(:), Y_loc(:)
		!%%%%%%%%%%%%%%%%%%%%
		real*8 :: d1, d2, d3, d, area, h
		!%%%%%%%%%%%%%%%%%%%%
		area = X_loc(2)*Y_loc(3) + X_loc(3)*Y_loc(1) + X_loc(1)*Y_loc(2) &
			-(X_loc(2)*Y_loc(1) + X_loc(3)*Y_loc(2) + X_loc(1)*Y_loc(3))
		d1 = dabs(X_loc(3)-X_loc(2)) + dabs(Y_loc(3)-Y_loc(2))
		d2 = dabs(X_loc(1)-X_loc(3)) + dabs(Y_loc(1)-Y_loc(3))
		d3 = dabs(X_loc(2)-X_loc(1)) + dabs(Y_loc(2)-Y_loc(1))
		d = d1 + d2 + d3
		!H ES MENOR O IGUAL A 1/4 LA ALTURA MINIMA
		h = dabs(area)/d
	end function h

	subroutine update_list(ipoin, smoothable)
		use PointNeighbor, only: psup1, psup2
		integer, intent(in) :: ipoin
		logical, intent(inout) :: smoothable(:)
		!%%%%%%%%%%%%%%%%%%%%
		integer :: ipsup
		!%%%%%%%%%%%%%%%%%%%%
		do ipsup = psup2(ipoin) + 1, psup2(ipoin + 1)
		smoothable(psup1(ipsup)) = .true.
		end do
	end subroutine update_list

	function get_mesh_quality(npoin, nelem, X, Y, inpoel)
		real*8 :: get_mesh_quality
		integer, intent(in) :: npoin, nelem, inpoel(3,nelem)
		real*8, intent(in) :: X(npoin), Y(npoin)
		!%%%%%%%%%%%%%%%%%%%% 
		integer :: ielem
		real*8 :: X_loc(3), Y_loc(3)
		real*8 :: mu1
		!%%%%%%%%%%%%%%%%%%%% 
		get_mesh_quality = 1
		do ielem = 1, nelem
		X_loc(:) = X(inpoel(:,ielem))
		Y_loc(:) = Y(inpoel(:,ielem))
		mu1 = mu(X_loc, Y_loc)
		if(mu1 < get_mesh_quality) get_mesh_quality = mu1
		end do
	end function get_mesh_quality

	subroutine check_smoothable(npoin, nelem, inpoel, X, Y, smoothable)
		integer, intent(in) :: npoin, nelem, inpoel(3,nelem)
		real*8, intent(in) :: X(npoin), Y(npoin) 
		logical, intent(out) :: smoothable(npoin)
		!%%%%%%%%%%%%%%%
		integer :: ielem, inode, jpoin
		real*8 :: X_loc(3), Y_loc(3)
		real*8 :: h1
		!%%%%%%%%%%%%%%%
		H_MIN_GLOBAL = 1
		do ielem = 1, nelem
		do inode = 1, 3
		jpoin = inpoel(inode,ielem)
		X_loc(inode) = X(jpoin)
		Y_loc(inode) = Y(jpoin)
		end do
		if(mu(X_loc, Y_loc) < TOL_METRIC) then
			smoothable(inpoel(:,ielem)) = .true.
		end if
		h1 = h(X_loc, Y_loc)
		if(h1 < H_MIN_GLOBAL) H_MIN_GLOBAL = h1
		end do
	end subroutine check_smoothable
end module smoothing_mod
