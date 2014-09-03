module smoothing_mod
	implicit none
	save
	private
	real*8, parameter :: TWOSQRT3 = 3.46410161513775d0
	integer, parameter :: MAX_ELEM_PER_NODE = 20
	integer, parameter :: NITER = 100, MITER = 2, NTRY = 8
	real*8, parameter :: TOL_METRIC = .85D0
	real*8, parameter :: FACTOR_TOL_DIST = 1d-2
	real*8, parameter :: FACTOR_DELTA = 1d-2
	real*8, parameter :: FACTOR_PLUS = 1.d0
	real*8, parameter :: FACTOR_STEP = 3.d0
	!%%%%%%%%%%%%%%%
	integer :: NELEM, NPOIN
	real*8 :: H_MIN_GLOBAL
	logical, dimension(:), allocatable :: smoothable
	real*8, dimension(:), allocatable :: tol_move
	logical, save :: SWITCH_TOL_MOVE = .false.
	public :: smoothing, setTol_move
contains
	subroutine smoothing(X, Y, inpoel, fixed, npoin0, nelem0, dX, dY)
		integer, intent(in) :: npoin0, nelem0, inpoel(3,nelem0)
		real*8, intent(inout) :: X(npoin0), Y(npoin0)
		logical, intent(in) :: fixed(npoin0)
		real*8, intent(in), optional :: dX(npoin), dY(npoin)
		!%%%%%%%%%%%%%%%
		integer :: iter, ipoin, min_idx
		real*8 :: d_max, TOL_DIST
		real*8, allocatable, dimension(:), save :: mu_vec
		logical, save :: isFirstCall = .true.
		integer, save :: counter = 0
		integer, parameter :: PRINT_INTERVAL = 500
		!%%%%%%%%%%%%%%%
		counter = counter + 1
		if(isFirstCall) then
			NELEM = nelem0
			NPOIN = npoin0
			allocate(mu_vec(MAX_ELEM_PER_NODE))
			allocate(smoothable(NPOIN))
			counter = PRINT_INTERVAL
			isFirstCall = .false.
		end if
		!%%%%%%%%%%%%%%%
		call checkMesh(inpoel, X, Y)
		if(.not. any(smoothable == .true.)) return
		if(present(dX).and.present(dY)) call setTol_move(dX, dY)
		do iter = 1, 2
		call laplacianSmoothing(X, Y, fixed, inpoel)
		end do
		TOL_DIST = FACTOR_TOL_DIST*H_MIN_GLOBAL
		do iter = 1, NITER
		d_max = 0.d0
		do ipoin = 1, NPOIN
		if(smoothable(ipoin) .and. .not. fixed(ipoin)) then
			smoothable(ipoin) = .false. 
			call getMu_vec(mu_vec, ipoin, inpoel, X, Y, min_idx)
			if(mu_vec(min_idx) < TOL_METRIC) then
				call moveIpoin(X, Y, ipoin, inpoel, mu_vec, min_idx, d_max)
			end if
		end if
		end do
		if(d_max < TOL_DIST) exit
		end do
		if(counter == PRINT_INTERVAL) then
			counter = 0
			print*, '=============SMOOTHING============='
			print*, '# iteraciones:', iter
			print*, 'Min Distortion Metric:', getMeshQuality(X, Y, inpoel)
			print*, 'H_MIN_GLOBAL:', H_MIN_GLOBAL
			print*, '==================================='
		end if
	end subroutine smoothing

	subroutine laplacianSmoothing(X, Y, fixed, inpoel)
		use PointNeighbor!, only: psup1, psup2
		integer, intent(in) :: inpoel(3,NELEM)
		real*8, intent(inout) :: X(NPOIN), Y(NPOIN)
		logical, intent(in) :: fixed(NPOIN)
		!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		integer :: ipoin, ipsup, n
		real*8 :: X_new, Y_new, X_old, Y_old, dX, dY, mu_old
		!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if(.not.allocated(psup1).or..not.allocated(psup2)) call getPsup(inpoel, NELEM, NPOIN)
		do ipoin = 1, NPOIN
		if(smoothable(ipoin).and..not.fixed(ipoin)) then
			X_new = 0
			Y_new = 0
			mu_old = getMu_min(X, Y, ipoin, inpoel)
			X_old = X(ipoin)
			Y_old = Y(ipoin)
			n = psup2(ipoin + 1) - psup2(ipoin)
			do ipsup = psup2(ipoin) + 1, psup2(ipoin + 1)
			X_new = X_new + X(psup1(ipsup))
			Y_new = Y_new + Y(psup1(ipsup))
			end do
			X_new = X_new/n
			Y_new = Y_new/n
			if(SWITCH_TOL_MOVE) then
				dX = X_new - X(ipoin)
				dY = Y_new - Y(ipoin)
				X(ipoin) = X(ipoin) + dX/(dabs(dX) + dabs(dY))*tol_move(ipoin)/2*FACTOR_STEP
				Y(ipoin) = Y(ipoin) + dY/(dabs(dX) + dabs(dY))*tol_move(ipoin)/2*FACTOR_STEP
			else
				X(ipoin) = X_new
				Y(ipoin) = Y_new
			end if
			!IF WORSE, RESTORE
			if(getMu_min(X, Y, ipoin, inpoel) < mu_old) then
				X(ipoin) = X_old
				Y(ipoin) = Y_old
			end if
		end if
		end do
	end subroutine laplacianSmoothing

	subroutine moveIpoin(X, Y, ipoin, inpoel, mu_vec, min_idx, d_max)
		integer, intent(in) :: ipoin, inpoel(3,NELEM)
		integer, intent(inout) :: min_idx
		real*8, intent(inout) :: mu_vec(MAX_ELEM_PER_NODE)
		real*8, intent(inout) :: X(NPOIN), Y(NPOIN)
		real*8, intent(inout) :: d_max
		!%%%%%%%%%%%%%%%%%%%%
		real*8, allocatable, dimension(:), save :: gx, gy
		real*8 :: x_old, y_old, x0, y0, mu_min
		real*8 :: dx, dy, step, d_move
		integer :: i, j, min_idx_new
		logical :: accepted
		logical, save :: isFirstCall = .true.
		!%%%%%%%%%%%%%%%%%%%%
		if(isFirstCall) then
			allocate(gx(MAX_ELEM_PER_NODE), gy(MAX_ELEM_PER_NODE))
			isFirstCall = .false.
		end if
		!%%%%%%%%%%%%%%%%%%%%
		x0 = X(ipoin)
		y0 = Y(ipoin)
		do i = 1, MITER
		mu_min = mu_vec(min_idx)
		x_old = X(ipoin)
		y_old = Y(ipoin)
		call getG(gx, gy, ipoin, inpoel, X, Y, mu_vec)
		step =  getStep(gx, gy, mu_vec, min_idx, ipoin, X, Y, inpoel)
		!TRY NEW POSITION
		accepted = .false.
		do j = 1, NTRY
		dx = step*gx(min_idx)
		dy = step*gy(min_idx)
		X(ipoin) = X(ipoin) + dx
		Y(ipoin) = Y(ipoin) + dy
		call getMu_vec(mu_vec, ipoin, inpoel, X, Y, min_idx_new)
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
		if(d_move > tiny(1.d0)) call update_list(ipoin, inpoel)
	end subroutine moveIpoin

	subroutine getMu_vec(mu_vec, ipoin, inpoel, X, Y, min_idx)
		use PointNeighbor!, only: esup1, esup2
		integer, intent(in) :: ipoin, inpoel(3,NELEM)
		real*8, intent(in) :: X(NPOIN), Y(NPOIN)
		real*8, intent(out) :: mu_vec(MAX_ELEM_PER_NODE)
		integer, optional, intent(out) :: min_idx
		!%%%%%%%%%%%%%%%%%%%%
		integer :: iesup
		real*8 :: X_loc(3), Y_loc(3)
		real*8 :: mu_min, mu1
		!%%%%%%%%%%%%%%%%%%%%
		if(.not.allocated(esup1) .or. .not.allocated(esup2)) call getEsup(inpoel, NELEM, NPOIN)
		if(present(min_idx)) min_idx = 1
		mu_min = 1
		do iesup = esup2(ipoin) + 1, esup2(ipoin + 1)
		X_loc(:) = X(inpoel(:, esup1(iesup)))
		Y_loc(:) = Y(inpoel(:, esup1(iesup)))
		mu1 = mu(X_loc, Y_loc)
		mu_vec(iesup - esup2(ipoin)) = mu1
		if(present(min_idx)) then
			if(mu1 < mu_min) then
				mu_min = mu1
				min_idx = iesup - esup2(ipoin)
			end if
		end if
		end do
	end subroutine getMu_vec

	real*8 function getMu_min(X, Y, ipoin, inpoel)
		use PointNeighbor
		integer, intent(in) :: ipoin
		real*8, intent(in) :: X(NPOIN), Y(NPOIN)
		integer, intent(in) :: inpoel(3,NELEM)
		!%%%%%%%%%%%%%%%%%%%%
		integer :: iesup
		real*8 :: X_loc(3), Y_loc(3)
		real*8 :: mu1
		!%%%%%%%%%%%%%%%%%%%%
		if(.not.allocated(esup1) .or. .not.allocated(esup2)) call getEsup(inpoel, NELEM, NPOIN)
		getMu_min = 1
		do iesup = esup2(ipoin) + 1, esup2(ipoin + 1)
		X_loc(:) = X(inpoel(:, esup1(iesup)))
		Y_loc(:) = Y(inpoel(:, esup1(iesup)))
		mu1 = mu(X_loc, Y_loc)
		if(mu1 <= getMu_min) getMu_min = mu1
		end do
	end function

	subroutine getG(gx, gy, ipoin, inpoel, X, Y, mu_vec)
		integer, intent(in) :: ipoin, inpoel(3,NELEM)
		real*8, intent(in) :: mu_vec(MAX_ELEM_PER_NODE)
		real*8, intent(inout) :: X(NPOIN), Y(NPOIN) 
		real*8, intent(out) :: gx(:), gy(:)
		!%%%%%%%%%%%%%%%%%%%%
		real*8 :: x_old, y_old
		real*8 :: DELTA
		!%%%%%%%%%%%%%%%%%%%%
		DELTA = H_MIN_GLOBAL*FACTOR_DELTA
		!--- GX ---
		x_old = X(ipoin)
		X(ipoin) = X(ipoin) + DELTA
		call getMu_vec(gx, ipoin, inpoel, X, Y)
		gx = (gx - mu_vec)/DELTA
		X(ipoin) = x_old
		!--- GY ---
		y_old = Y(ipoin)
		Y(ipoin) = Y(ipoin) + DELTA
		call getMu_vec(gy, ipoin, inpoel, X, Y)
		gy = (gy - mu_vec)/DELTA
		Y(ipoin) = y_old
	end subroutine getG

	real*8 function getStep(gx, gy, mu_vec, min_idx, ipoin, X, Y, inpoel)
		use PointNeighbor!, only: esup1, esup2
		real*8, intent(in) :: gx(:), gy(:), mu_vec(MAX_ELEM_PER_NODE)
		integer, intent(in) :: min_idx, ipoin, inpoel(3,NELEM)
		real*8, intent(in) :: X(NPOIN), Y(NPOIN)
		!%%%%%%%%%%%%%%%%%%%%
		real*8 :: g2, gg, step1
		real*8 :: gx_min, gy_min, mu_min
		integer :: i
		!%%%%%%%%%%%%%%%%%%%%
		if(.not.allocated(esup1) .or. .not.allocated(esup2)) call getEsup(inpoel, NELEM, NPOIN)
		gx_min = gx(min_idx)
		gy_min = gy(min_idx)
		mu_min = mu_vec(min_idx)
		g2 = gx_min**2 + gy_min**2
		if(SWITCH_TOL_MOVE) then
			getStep = tol_move(ipoin)/((dabs(gx_min) + dabs(gy_min))*NITER)
		else
			getStep = getH_min(ipoin, X, Y, inpoel)*FACTOR_STEP/(dabs(gx_min) + dabs(gy_min))
		end if
		do i = 1, esup2(ipoin + 1) - esup2(ipoin)
		gg = gx_min*gx(i) + gy_min*gy(i)
		if(gg < 0) then
			step1 = (mu_vec(i) - mu_min)/(g2 - gg)
			if(step1 < getStep) getStep = step1
		end if
		end do
	end function getStep

	real*8 function getH_min(ipoin, X, Y, inpoel)
		use PointNeighbor!, only: esup1, esup2
		integer, intent(in) :: ipoin, inpoel(3,NELEM)
		real*8, intent(in) :: X(NPOIN), Y(NPOIN)
		!%%%%%%%%%%%%%%%%%%%% 
		integer :: iesup
		real*8 :: X_loc(3), Y_loc(3)
		real*8 :: h1
		!%%%%%%%%%%%%%%%%%%%% 
		if(.not.allocated(esup1) .or. .not.allocated(esup2)) call getEsup(inpoel, NELEM, NPOIN)
		getH_min = 1
		do iesup = esup2(ipoin) + 1, esup2(ipoin + 1)
		X_loc(:) = X(inpoel(:,esup1(iesup)))
		Y_loc(:) = Y(inpoel(:,esup1(iesup)))
		h1 = h(X_loc, Y_loc)
		if(h1 < getH_min) getH_min = h1
		end do
	end function getH_min

	subroutine update_list(ipoin, inpoel)
		use PointNeighbor!, only: psup1, psup2
		integer, intent(in) :: ipoin, inpoel(3,NELEM)
		!%%%%%%%%%%%%%%%%%%%%
		integer :: ipsup
		!%%%%%%%%%%%%%%%%%%%%
		if(.not.allocated(psup1) .or. .not.allocated(psup2)) call getPsup(inpoel, NELEM, NPOIN)
		forall (ipsup = psup2(ipoin) + 1 : psup2(ipoin + 1))
			smoothable(psup1(ipsup)) = .true.
		end forall
	end subroutine update_list

	pure real*8 function mu(X_loc, Y_loc)
		real*8, intent(in) :: X_loc(:), Y_loc(:)
		!%%%%%%%%%%%%%%%%%%%%
		real*8 :: l1, l2, l3, l, area
		!%%%%%%%%%%%%%%%%%%%%
		area = X_loc(2)*Y_loc(3)+X_loc(3)*Y_loc(1)+X_loc(1)*Y_loc(2) &
			-(X_loc(2)*Y_loc(1)+X_loc(3)*Y_loc(2)+X_loc(1)*Y_loc(3))
		l1 = (X_loc(3)-X_loc(2))**2 + (Y_loc(3)-Y_loc(2))**2
		l2 = (X_loc(1)-X_loc(3))**2 + (Y_loc(1)-Y_loc(3))**2
		l3 = (X_loc(2)-X_loc(1))**2 + (Y_loc(2)-Y_loc(1))**2
		l = l1 + l2 + l3
		mu = TWOSQRT3*area/l
	end function mu

	pure real*8 function h(X_loc, Y_loc)
		real*8, intent(in) :: X_loc(:), Y_loc(:)
		!%%%%%%%%%%%%%%%%%%%%
		real*8 :: d1, d2, d3, d, area
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

	pure real*8 function getMeshQuality(X, Y, inpoel)
		integer, intent(in) :: inpoel(3,NELEM)
		real*8, intent(in) :: X(NPOIN), Y(NPOIN)
		!%%%%%%%%%%%%%%%%%%%% 
		integer :: ielem
		real*8 :: X_loc(3), Y_loc(3)
		real*8 :: mu1
		!%%%%%%%%%%%%%%%%%%%% 
		getMeshQuality = 1
		do ielem = 1, nelem
		X_loc(:) = X(inpoel(:,ielem))
		Y_loc(:) = Y(inpoel(:,ielem))
		mu1 = mu(X_loc, Y_loc)
		if(mu1 < getMeshQuality) getMeshQuality = mu1
		end do
	end function getMeshQuality

	subroutine checkMesh(inpoel, X, Y)
		integer, intent(in) :: inpoel(3,NELEM)
		real*8, intent(in) :: X(NPOIN), Y(NPOIN) 
		!%%%%%%%%%%%%%%%%%%%%%%%%%
		integer :: ielem
		real*8 :: X_loc(3), Y_loc(3)
		real*8 :: h1
		!%%%%%%%%%%%%%%%%%%%%%%%%%
		H_MIN_GLOBAL = 1
		smoothable(:) = .false.
		do ielem = 1, NELEM
		X_loc(:) = X(inpoel(:,ielem))
		Y_loc(:) = Y(inpoel(:,ielem))
		if(mu(X_loc, Y_loc) < TOL_METRIC) then
			smoothable(inpoel(:,ielem)) = .true.
		end if
		h1 = h(X_loc, Y_loc)
		if(h1 < H_MIN_GLOBAL) H_MIN_GLOBAL = h1
		end do
	end subroutine checkMesh

	subroutine setTol_move(dX, dY)
		real*8, dimension(:), intent(in) :: dX, dY
		!%%%%%%%%%%%%%%%%%%%%%%%%%
		integer :: ipoin, npoin
		!%%%%%%%%%%%%%%%%%%%%%%%%%
		npoin = size(dX)
		if(.not.allocated(tol_move)) allocate(tol_move(npoin))
		!%%%%%%%%%%%%%%%%%%%%%%%%%
		SWITCH_TOL_MOVE = .true.
		!$OMP PARALLEL DO PRIVATE(ipoin)
		do ipoin = 1, npoin
		tol_move(ipoin) = dabs(dX(ipoin)) + dabs(dY(ipoin))
		end do
		!$OMP END PARALLEL DO
	end subroutine setTol_move
end module smoothing_mod
