program calcRHS
	use omp_lib
	parameter (nelem=2000000, npoin=1000000, loopsize=100)
	integer start_t, end_t, rate, start_mp 
	integer, allocatable, dimension(:,:) :: n
	real(8), allocatable, dimension(:,:) :: dnx, dny, u1, un, rhs
	real(8), allocatable, dimension(:) :: t_sugn1, t_sugn2, t_sugn3, dtl, ps, area, hhx, hhy, shoc
	real(8), allocatable, dimension(:) :: p, t, gamm
	real(8) fr, rmu, fk, fcv, tinf, cte
	allocate(n(3,nelem))
	allocate(p(npoin), t(npoin), gamm(npoin))
	allocate(dnx(3,nelem), dny(3,nelem)); allocate(u1(4,npoin), un(4,npoin)); allocate(rhs(4,npoin))
	allocate(t_sugn1(nelem), t_sugn2(nelem), t_sugn3(nelem)); allocate(dtl(nelem)); allocate(ps(nelem))
	allocate(area(nelem)); allocate(hhx(nelem), hhy(nelem)); allocate(shoc(nelem))
	fr=0; fcv=0; rmu=0; cte=1; tinf=0

	call fill_inpoel(n, nelem, npoin); call fill_3xnelem(dnx, nelem); call fill_3xnelem(dny, nelem)
	call fill_nelem(t_sugn1, nelem); call fill_nelem(t_sugn2, nelem); call fill_nelem(t_sugn3, nelem)
	call fill_nelem(shoc, nelem); call fill_nelem(dtl, nelem); call fill_nelem(area, nelem)
	call fill_nelem(hhx, nelem); call fill_nelem(hhy, nelem); call fill_u(u1, npoin)
	call fill_u(un, npoin); call fill_npoin(p, npoin); call fill_npoin(gamm, npoin)
	call fill_npoin(t, npoin); call fill_npoin(ps, npoin)

	if(.FALSE.) then 
	rhs=0.d0
	!$ start_mp = omp_get_wtime()
	do i=1,loopsize
		call todo(npoin,nelem,n,dnx,dny,area,hhx,hhy &
             ,t_sugn1,t_sugn2,t_sugn3,shoc,ps,dtl &
             ,u1,un,rhs,p,gamm,fr,rmu,fk,fcv,tinf,cte)
	end do
	!$ print *, "Serial version: ", omp_get_wtime()-start_mp, "sec"
	call checkit(rhs,npoin,c)
	print *, "check... ", c
	print *,
	rhs=0.d0

	rhs=0.d0
	!$ start_mp = omp_get_wtime()
	do i=1,loopsize
		call todo_mp_atomic_only(npoin,nelem,n,dnx,dny,area,hhx,hhy &
             ,t_sugn1,t_sugn2,t_sugn3,shoc,ps,dtl &
             ,u1,un,rhs,p,gamm,fr,rmu,fk,fcv,tinf,cte)
	end do
	!$ print *, "OMP atomic version: ", omp_get_wtime()-start_mp, "sec"
	call checkit(rhs,npoin,c)
	print *, "check... ", c
	print *,

	rhs=0.d0
	!$ start_mp = omp_get_wtime()
	do i=1,loopsize
		call todo_mp_critical_fp(npoin,nelem,n,dnx,dny,area,hhx,hhy &
             ,t_sugn1,t_sugn2,t_sugn3,shoc,ps,dtl &
             ,u1,un,rhs,p,gamm,fr,rmu,fk,fcv,tinf,cte)
	end do
	!$ print *, "OMP private+critical_FP version: ", omp_get_wtime()-start_mp, "sec"
	call checkit(rhs,npoin,c)
	print *, "check... ", c
	print *,

	rhs=0.d0
	!$ start_mp = omp_get_wtime()
	do i=1,loopsize
		call todo_mp_critical_only(npoin,nelem,n,dnx,dny,area,hhx,hhy &
             ,t_sugn1,t_sugn2,t_sugn3,shoc,ps,dtl &
             ,u1,un,rhs,p,gamm,fr,rmu,fk,fcv,tinf,cte)
	end do
	!$ print *, "OMP critical version: ", omp_get_wtime()-start_mp, "sec"
	call checkit(rhs,npoin,c)
	print *, "check... ", c
	print *,

	rhs=0.d0
	!$ start_mp = omp_get_wtime()
	do i=1,loopsize
		call todo_mp_atomic(npoin,nelem,n,dnx,dny,area,hhx,hhy &
             ,t_sugn1,t_sugn2,t_sugn3,shoc,ps,dtl &
             ,u1,un,rhs,p,gamm,fr,rmu,fk,fcv,tinf,cte)
	end do
	!$ print *, "OMP private+atomic version: ", omp_get_wtime()-start_mp, "sec"
	call checkit(rhs,npoin,c)
	print *, "check... ", c

	rhs=0.d0
	!$ start_mp = omp_get_wtime()
	do i=1,loopsize
		call todo_mp_noAtomic(npoin,nelem,n,dnx,dny,area,hhx,hhy &
             ,t_sugn1,t_sugn2,t_sugn3,shoc,ps,dtl &
             ,u1,un,rhs,p,gamm,fr,rmu,fk,fcv,tinf,cte)
	end do
	!$ print *, "OMP noAtomic version: ", omp_get_wtime()-start_mp, "sec"
	call checkit(rhs,npoin,c)
	print *, "check... ", c
	print *,
	end if

	rhs=0.d0
	!$ start_mp = omp_get_wtime()
	do i=1,loopsize
		call todo_mp_critical(npoin,nelem,n,dnx,dny,area,hhx,hhy &
             ,t_sugn1,t_sugn2,t_sugn3,shoc,ps,dtl &
             ,u1,un,rhs,p,gamm,fr,rmu,fk,fcv,tinf,cte)
	end do
	!$ print *, "OMP private+critical version: ", omp_get_wtime()-start_mp, "sec"
	call checkit(rhs,npoin,c)
	print *, "check... ", c
	print *,

	rhs=0.d0
	!$ start_mp = omp_get_wtime()
	do i=1,loopsize
		call todo_mp_critical_only_outside(npoin,nelem,n,dnx,dny,area,hhx,hhy &
             ,t_sugn1,t_sugn2,t_sugn3,shoc,ps,dtl &
             ,u1,un,rhs,p,gamm,fr,rmu,fk,fcv,tinf,cte)
	end do
	!$ print *, "OMP critical outside version: ", omp_get_wtime()-start_mp, "sec"
	call checkit(rhs,npoin,c)
	print *, "check... ", c
	print *,

	rhs=0.d0
	!$ start_mp = omp_get_wtime()
	do i=1,loopsize
		call todo_mp_atomic_only_outside(npoin,nelem,n,dnx,dny,area,hhx,hhy &
             ,t_sugn1,t_sugn2,t_sugn3,shoc,ps,dtl &
             ,u1,un,rhs,p,gamm,fr,rmu,fk,fcv,tinf,cte)
	end do
	!$ print *, "OMP atomic outside version: ", omp_get_wtime()-start_mp, "sec"
	call checkit(rhs,npoin,c)
	print *, "check... ", c
	print *,

end program calcRHS

subroutine todo(nnod,nelem,n,dnx,dny,area,hhx,hhy &
     ,t_sugn1,t_sugn2,t_sugn3,shoc,ps,dtl &
     ,u,un,rhs,p,gamm,fr,rmu,fk,fcv,tinf,cte)

  implicit real(8) (a-h,o-z)

  integer n(3,nelem)

  real(8) dnx(3,nelem),dny(3,nelem),u(4,nnod),p(nnod),un(4,nnod),t(nnod),gamm(nnod)
  real(8) area(nelem),rhs(4,nnod),hhx(nelem),hhy(nelem)

  real(8) alf(3),bet(3)
  real(8) t_sugn1(nelem),t_sugn2(nelem),t_sugn3(nelem),shoc(nelem)
  real(8) ps(nnod),dtl(nelem)

  !c---->     rho=u(1)
  !c---->     rho*vel_x=u(2)
  !c---->     rho*vel_y=u(3)
  !c---->     rho*et=u(4)

  data alf/.5d0,.5d0,0.d0/
  data bet/0.d0,.5d0,.5d0/

  ngauss=3    !ptos de gauss donde voy a intergrar

  do ielem=1,nelem

     n1=n(1,ielem) ; n2=n(2,ielem) ; n3=n(3,ielem)

     gama=(gamm(n1)+gamm(n2)+gamm(n3))/3.d0
     gm=gama-1.d0
     temp=(t(n1)+t(n2)+t(n3))/3.d0
     fmu= 1.716d-5*162.6/(temp-110.55)*(temp/273.15)**.75d0     !sutherland

     rnx1=dnx(1,ielem) ; rnx2=dnx(2,ielem) ; rnx3=dnx(3,ielem)
     rny1=dny(1,ielem) ; rny2=dny(2,ielem) ; rny3=dny(3,ielem)

     !cccc  ----> derivada de rho
     drx= rnx1*u(1,n1)+rnx2*u(1,n2)+rnx3*u(1,n3)
     dry= rny1*u(1,n1)+rny2*u(1,n2)+rny3*u(1,n3)
     !cccc  ----> derivada de rho.vel_x
     drux= rnx1*u(2,n1)+rnx2*u(2,n2)+rnx3*u(2,n3)
     druy= rny1*u(2,n1)+rny2*u(2,n2)+rny3*u(2,n3)
     !cccc  ----> derivada de rho.vel_y
     drvx= rnx1*u(3,n1)+rnx2*u(3,n2)+rnx3*u(3,n3)
     drvy= rny1*u(3,n1)+rny2*u(3,n2)+rny3*u(3,n3)
     !cccc  ----> derivada de rho.et
     drex= rnx1*u(4,n1)+rnx2*u(4,n2)+rnx3*u(4,n3)
     drey= rny1*u(4,n1)+rny2*u(4,n2)+rny3*u(4,n3)

     ar=area(ielem)*dtl(ielem)/3.d0
     !cccc  ----> long. caracteristica
     hlongx=hhx(ielem)
     hlongy=hhy(ielem)
     hlong=dsqrt(area(ielem))
     !cccc  ----> estab. tezduyar
     tau1=t_sugn1(ielem)
     tau2=t_sugn2(ielem)
     tau3=t_sugn3(ielem)
     alfa_mu=shoc(ielem)

     do j=1,ngauss

        rn1=1.d0-alf(j)-bet(j)
        rn2=alf(j)
        rn3=bet(j)

        !cccc  ----> integro las variables en los puntos de gauss
        u1= rn1*u(1,n1)+rn2*u(1,n2)+rn3*u(1,n3)
        u2= rn1*u(2,n1)+rn2*u(2,n2)+rn3*u(2,n3)
        u3= rn1*u(3,n1)+rn2*u(3,n2)+rn3*u(3,n3)
        u4= rn1*u(4,n1)+rn2*u(4,n2)+rn3*u(4,n3)

        fi_1= rn1*un(1,n1)+rn2*un(1,n2)+rn3*un(1,n3)
        fi_2= rn1*un(2,n1)+rn2*un(2,n2)+rn3*un(2,n3)
        fi_3= rn1*un(3,n1)+rn2*un(3,n2)+rn3*un(3,n3)
        fi_4= rn1*un(4,n1)+rn2*un(4,n2)+rn3*un(4,n3)
        !cccc  ----> defino variables primitivas
        vx=u2/u1
        vy=u3/u1
        et=u4/u1
        rmod2=vx*vx+vy*vy
        temp=gm/fr*(et-.5d0*rmod2) !fr=cte. universal de los gases

        c=dsqrt(gama*fr*temp)
        
        !cccc  ----> definicion de las matrices a1 y a2
        !cccc  ----> a1
        a1_21= gm/2.d0*rmod2-vx*vx ; a1_22=(3.d0-gama)*vx ; a1_23=-gm*vy ; a1_24=gm
        a1_31= -vx*vy ; a1_32=vy ; a1_33=vx
        a1_41=(gm*rmod2-gama*et)*vx ; a1_42=gama*et-gm/2.d0*rmod2-gm*vx*vx ; a1_43=-gm*vx*vy ; a1_44=gama*vx
        !cccc  ----> a2
        a2_21=-vx*vy ; a2_22=vy ; a2_23=vx
        a2_31=gm/2.d0*rmod2-vy*vy ; a2_32=-gm*vx ; a2_33=(3.d0-gama)*vy ; a2_34=gm
        a2_41=(gm*rmod2-gama*et)*vy ; a2_42=-gm*vx*vy ; a2_43=gama*et-gm/2.d0*rmod2-gm*vy*vy ; a2_44=gama*vy
       
        !cccc----> terminos de estabilizacion y captura de choque
        arr1=ar*tau1
        arr2=ar*tau2
        arr3=ar*tau3

        !cccc---->  calculo de mua y sus componentes  <----!cccc
        v11=c+dabs(vx)+dabs(vy)

        !cccc  ----> armo el pressure switch
        choq1=alfa_mu*ar*ps(n1)*cte
        choq2=alfa_mu*ar*ps(n2)*cte
        choq3=alfa_mu*ar*ps(n3)*cte

        !cccc  ----> multiplico por partes para simplificar el asunto
        !cccc  ----> 'a' por las derivadas
        auxa1=                  drux                                       +                    drvy
        auxa2=a1_21*drx + a1_22*drux + a1_23*drvx + a1_24*drex + a2_21*dry + a2_22*druy + a2_23*drvy
        auxa3=a1_31*drx + a1_32*drux + a1_33*drvx +              a2_31*dry + a2_32*druy + a2_33*drvy + a2_34*drey
        auxa4=a1_41*drx + a1_42*drux + a1_43*drvx + a1_44*drex + a2_41*dry + a2_42*druy + a2_43*drvy + a2_44*drey

        auxa11=auxa1+fi_1
        auxa22=auxa2+fi_2
        auxa33=auxa3+fi_3
        auxa44=auxa4+fi_4

        !cccc  ----> lo anterior por 'a' transpuesta
        aa1=                     auxa22
        aa2=a1_21*auxa11 + a1_22*auxa22 + a1_23*auxa33 + a1_24*auxa44
        aa3=a1_31*auxa11 + a1_32*auxa22 + a1_33*auxa33
        aa4=a1_41*auxa11 + a1_42*auxa22 + a1_43*auxa33 + a1_44*auxa44
        aa5=                                    auxa33
        aa6=a2_21*auxa11 + a2_22*auxa22 + a2_23*auxa33
        aa7=a2_31*auxa11 + a2_32*auxa22 + a2_33*auxa33 + a2_34*auxa44
        aa8=a2_41*auxa11 + a2_42*auxa22 + a2_43*auxa33 + a2_44*auxa44

		!cccc  ----> ensamble del right hand side del nodo n1
        rhs(1,n1)=rhs(1,n1)+(rnx1*aa1+rny1*aa5)*arr1 +auxa1*rn1*ar +(rnx1*drx+rny1*dry)*choq1
        rhs(2,n1)=rhs(2,n1)+(rnx1*aa2+rny1*aa6)*arr2 +auxa2*rn1*ar +(rnx1*drux+rny1*druy)*choq1
        rhs(3,n1)=rhs(3,n1)+(rnx1*aa3+rny1*aa7)*arr2 +auxa3*rn1*ar +(rnx1*drvx+rny1*drvy)*choq1
        rhs(4,n1)=rhs(4,n1)+(rnx1*aa4+rny1*aa8)*arr3 +auxa4*rn1*ar +(rnx1*drex+rny1*drey)*choq1
        !cccc  ----> ensamble del right hand side del nodo n2
        rhs(1,n2)=rhs(1,n2)+(rnx2*aa1+rny2*aa5)*arr1 +auxa1*rn2*ar +(rnx2*drx+rny2*dry)*choq2
        rhs(2,n2)=rhs(2,n2)+(rnx2*aa2+rny2*aa6)*arr2 +auxa2*rn2*ar +(rnx2*drux+rny2*druy)*choq2
        rhs(3,n2)=rhs(3,n2)+(rnx2*aa3+rny2*aa7)*arr2 +auxa3*rn2*ar +(rnx2*drvx+rny2*drvy)*choq2
        rhs(4,n2)=rhs(4,n2)+(rnx2*aa4+rny2*aa8)*arr3 +auxa4*rn2*ar +(rnx2*drex+rny2*drey)*choq2
        !cccc  ----> ensamble del right hand side del nodo n3
        rhs(1,n3)=rhs(1,n3)+(rnx3*aa1+rny3*aa5)*arr1 +auxa1*rn3*ar +(rnx3*drx+rny3*dry)*choq3
        rhs(2,n3)=rhs(2,n3)+(rnx3*aa2+rny3*aa6)*arr2 +auxa2*rn3*ar +(rnx3*drux+rny3*druy)*choq3
        rhs(3,n3)=rhs(3,n3)+(rnx3*aa3+rny3*aa7)*arr2 +auxa3*rn3*ar +(rnx3*drvx+rny3*drvy)*choq3
        rhs(4,n3)=rhs(4,n3)+(rnx3*aa4+rny3*aa8)*arr3 +auxa4*rn3*ar +(rnx3*drex+rny3*drey)*choq3
     end do
  end do
  return
end subroutine todo

subroutine todo_mp_atomic_only(nnod,nelem,n,dnx,dny,area,hhx,hhy &
    ,t_sugn1,t_sugn2,t_sugn3,shoc,ps,dtl &
    ,u,un,rhs,p,gamm,fr,rmu,fk,fcv,tinf,cte)
	use omp_lib
    implicit real(8) (a-h,o-z)
    integer n(3,nelem), ipoin(3)
    real(8) dnx(3,nelem),dny(3,nelem),u(4,nnod),p(nnod),un(4,nnod),t(nnod),gamm(nnod)
    real(8) area(nelem),rhs(4,nnod),hhx(nelem),hhy(nelem)
	real (8) a1(4,4), a2(4,4), ux(4), uy(4), u_loc(4), phi_loc(4), aux(4), aa(8), aux_phi(4), nx(3), ny(3)
    real(8) alf(3),bet(3)
    real(8) t_sugn1(nelem),t_sugn2(nelem),t_sugn3(nelem),shoc(nelem)
    real(8) ps(nnod),dtl(nelem)

    data alf/.5d0,.5d0,0.d0/
    data bet/0.d0,.5d0,.5d0/

    ngauss=3    !ptos de gauss donde voy a intergrar

	!$omp parallel &
	!$omp private(ielem,ipoin,gama,gm,temp,fmu,nx,ny,ux,uy,ar,hlong,hlongx,hlongy,tau1,tau2,tau3,alfa_mu,i,j,&
	!$omp rn1,rn2,rn3,u_loc,phi_loc,vx,vy,rmod2,et,c,a1,a2,aux,aa,aux_phi,choq1,choq2,choq3,arr1,arr2,arr3)

	!$omp do
    do ielem=1,nelem

		ipoin = n(:,ielem)

        gama=(gamm(ipoin(1))+gamm(ipoin(2))+gamm(ipoin(3)))/3.d0
        gm=gama-1.d0
        temp=(t(ipoin(1))+t(ipoin(2))+t(ipoin(3)))/3.d0
        fmu= 1.716d-5*162.6/(temp-110.55)*(temp/273.15)**.75d0     !sutherland

		nx = dnx(:,ielem)
		ny = dny(:,ielem)

		ux(:) = u(:,ipoin(1))*nx(1) + u(:,ipoin(2))*nx(2) + u(:,ipoin(3))*nx(3)
		uy(:) = u(:,ipoin(1))*ny(1) + u(:,ipoin(2))*ny(2) + u(:,ipoin(3))*ny(3)

        ar=area(ielem)*dtl(ielem)/3.d0
        !cccc  ----> long. caracteristica
        hlongx=hhx(ielem)
        hlongy=hhy(ielem)
        hlong=dsqrt(area(ielem))
        !cccc  ----> estab. tezduyar
        tau1=t_sugn1(ielem)
        tau2=t_sugn2(ielem)
        tau3=t_sugn3(ielem)
        alfa_mu=shoc(ielem)

        do j=1,ngauss

            rn1=1.d0-alf(j)-bet(j)
            rn2=alf(j)
            rn3=bet(j)

            !cccc  ----> integro las variables en los puntos de gauss
			u_loc = rn1*u(:,ipoin(1)) + rn2*u(:,ipoin(2)) + rn3*u(:,ipoin(3))
			phi_loc = rn1*un(:,ipoin(1)) + rn2*un(:,ipoin(2)) + rn3*un(:,ipoin(3))		

            !cccc  ----> defino variables primitivas
            vx=u_loc(2)/u_loc(1)
            vy=u_loc(3)/u_loc(1)
            et=u_loc(4)/u_loc(1)
            rmod2=vx*vx+vy*vy
            temp=gm/fr*(et-.5d0*rmod2) !fr=cte. universal de los gases

            c=dsqrt(gama*fr*temp)

            !cccc  ----> definicion de las matrices a1 y a2

            a1(:,1) = (/ 0.d0				  , 1.d0				          , 0.d0		  , 0.d0                       /)
            a1(:,2) = (/ gm/2.d0*rmod2-vx*vx  , (3.d0-gama)*vx                , -gm*vy        , gm                         /)
            a1(:,3) = (/ -vx*vy 			  , vy 			                  , vx	          , 0.d0                       /)
            a1(:,4) = (/ (gm*rmod2-gama*et)*vx, gama*et-gm/2.d0*rmod2-gm*vx*vx, -gm*vx*vy     , gama*vx                    /)

            a2(:,1) = (/ 0.d0				  , 0.d0						  , 1.d0	      , 0.d0                       /)
            a2(:,2) = (/ -vx*vy               , vy                            , vx            , 0.d0 			           /)
            a2(:,3) = (/ gm/2.d0*rmod2-vy*vy  , -gm*vx                        , (3.d0-gama)*vy, gm  		               /)
            a2(:,4) = (/ (gm*rmod2-gama*et)*vy, -gm*vx*vy                     , gama*et-gm/2.d0*rmod2-gm*vy*vy,    gama*vy /)

            !cccc----> terminos de estabilizacion y captura de choque
            arr1=ar*tau1
            arr2=ar*tau2
            arr3=ar*tau3

            !cccc  ----> armo el pressure switch
            choq1=alfa_mu*ar*ps(ipoin(1))*cte
            choq2=alfa_mu*ar*ps(ipoin(2))*cte
            choq3=alfa_mu*ar*ps(ipoin(3))*cte

            !cccc  ----> multiplico por partes para simplificar el asunto
            !cccc  ----> 'a' por las derivadas
        	aux(1)=                  ux(2)                                       +                    uy(3)
        	aux(2)=a1(1,2)*ux(1) + a1(2,2)*ux(2) + a1(3,2)*ux(3) + a1(4,2)*ux(4) + a2(1,2)*uy(1) + a2(2,2)*uy(2) + a2(3,2)*uy(3)
            aux(3)=a1(1,3)*ux(1) + a1(2,3)*ux(2) + a1(3,3)*ux(3) +                 a2(1,3)*uy(1) + a2(2,3)*uy(2) + a2(3,3)*uy(3)&
				+ a2(4,3)*uy(4)
        	aux(4)=a1(1,4)*ux(1) + a1(2,4)*ux(2) + a1(3,4)*ux(3) + a1(4,4)*ux(4) + a2(1,4)*uy(1) + a2(2,4)*uy(2) + a2(3,4)*uy(3)&
				+ a2(4,4)*uy(4)

			aux_phi = aux + phi_loc

            !cccc  ----> lo anterior por 'a' transpuesta
            aa(1)=                     aux_phi(2)
            aa(2)=a1(1,2)*aux_phi(1) + a1(2,2)*aux_phi(2) + a1(3,2)*aux_phi(3) + a1(4,2)*aux_phi(4)
            aa(3)=a1(1,3)*aux_phi(1) + a1(2,3)*aux_phi(2) + a1(3,3)*aux_phi(3)
            aa(4)=a1(1,4)*aux_phi(1) + a1(2,4)*aux_phi(2) + a1(3,4)*aux_phi(3) + a1(4,4)*aux_phi(4)
            aa(5)=                                   	    aux_phi(3)
            aa(6)=a2(1,2)*aux_phi(1) + a2(2,2)*aux_phi(2) + a2(3,2)*aux_phi(3)
            aa(7)=a2(1,3)*aux_phi(1) + a2(2,3)*aux_phi(2) + a2(3,3)*aux_phi(3) + a2(4,3)*aux_phi(4)
            aa(8)=a2(1,4)*aux_phi(1) + a2(2,4)*aux_phi(2) + a2(3,4)*aux_phi(3) + a2(4,4)*aux_phi(4)

            !cccc  ----> ensamble del right hand side del nodo ipoin(1)
			!$omp atomic
            rhs(1,ipoin(1))=rhs(1,ipoin(1))+(nx(1)*aa(1)+ny(1)*aa(5))*arr1 +aux(1)*rn1*ar&
					+(nx(1)*ux(1)+ny(1)*uy(1))*choq1
			!$omp atomic
            rhs(2,ipoin(1))=rhs(2,ipoin(1))+(nx(1)*aa(2)+ny(1)*aa(6))*arr2 +aux(2)*rn1*ar&
					+(nx(1)*ux(2)+ny(1)*uy(2))*choq1
			!$omp atomic
            rhs(3,ipoin(1))=rhs(3,ipoin(1))+(nx(1)*aa(3)+ny(1)*aa(7))*arr2 +aux(3)*rn1*ar&
					+(nx(1)*ux(3)+ny(1)*uy(3))*choq1
			!$omp atomic
            rhs(4,ipoin(1))=rhs(4,ipoin(1))+(nx(1)*aa(4)+ny(1)*aa(8))*arr3 +aux(4)*rn1*ar&
					+(nx(1)*ux(4)+ny(1)*uy(4))*choq1
            !cccc  ----> ensamble del right hand side del nodo ipoin(2)
			!$omp atomic
            rhs(1,ipoin(2))=rhs(1,ipoin(2))+(nx(2)*aa(1)+ny(2)*aa(5))*arr1 +aux(1)*rn2*ar&
					+(nx(2)*ux(1)+ny(2)*uy(1))*choq2
			!$omp atomic
            rhs(2,ipoin(2))=rhs(2,ipoin(2))+(nx(2)*aa(2)+ny(2)*aa(6))*arr2 +aux(2)*rn2*ar&
					+(nx(2)*ux(2)+ny(2)*uy(2))*choq2
			!$omp atomic
            rhs(3,ipoin(2))=rhs(3,ipoin(2))+(nx(2)*aa(3)+ny(2)*aa(7))*arr2 +aux(3)*rn2*ar&
					+(nx(2)*ux(3)+ny(2)*uy(3))*choq2
			!$omp atomic
            rhs(4,ipoin(2))=rhs(4,ipoin(2))+(nx(2)*aa(4)+ny(2)*aa(8))*arr3 +aux(4)*rn2*ar&
					+(nx(2)*ux(4)+ny(2)*uy(4))*choq2
            !cccc  ----> ensamble del right hand side del nodo ipoin(3)
			!$omp atomic
            rhs(1,ipoin(3))=rhs(1,ipoin(3))+(nx(3)*aa(1)+ny(3)*aa(5))*arr1 +aux(1)*rn3*ar&
					+(nx(3)*ux(1)+ny(3)*uy(1))*choq3
			!$omp atomic
            rhs(2,ipoin(3))=rhs(2,ipoin(3))+(nx(3)*aa(2)+ny(3)*aa(6))*arr2 +aux(2)*rn3*ar&
					+(nx(3)*ux(2)+ny(3)*uy(2))*choq3
			!$omp atomic
            rhs(3,ipoin(3))=rhs(3,ipoin(3))+(nx(3)*aa(3)+ny(3)*aa(7))*arr2 +aux(3)*rn3*ar& 
					+(nx(3)*ux(3)+ny(3)*uy(3))*choq3
			!$omp atomic
            rhs(4,ipoin(3))=rhs(4,ipoin(3))+(nx(3)*aa(4)+ny(3)*aa(8))*arr3 +aux(4)*rn3*ar& 
					+(nx(3)*ux(4)+ny(3)*uy(4))*choq3
        enddo
    enddo
	!$omp end do
	!$omp end parallel

    return
end subroutine todo_mp_atomic_only

subroutine todo_mp_atomic_only_outside(nnod,nelem,n,dnx,dny,area,hhx,hhy &
    ,t_sugn1,t_sugn2,t_sugn3,shoc,ps,dtl &
    ,u,un,rhs,p,gamm,fr,rmu,fk,fcv,tinf,cte)
	use omp_lib
    implicit real(8) (a-h,o-z)
    integer n(3,nelem), ipoin(3)
    real(8) dnx(3,nelem),dny(3,nelem),u(4,nnod),p(nnod),un(4,nnod),t(nnod),gamm(nnod)
    real(8) area(nelem),rhs(4,nnod),hhx(nelem),hhy(nelem)
	real(8) a1(3,4), a2(3,4), aux(4), aa(8), aux_phi(4), nx(3), ny(3), ux(4), uy(4), u_loc(4)
	real(8) arr(4), tau(4), choq(3), phi_loc(4), rhs_temp(4,3)
    real(8) alf(3),bet(3)
    real(8) t_sugn1(nelem),t_sugn2(nelem),t_sugn3(nelem),shoc(nelem)
    real(8) ps(nnod),dtl(nelem)

    data alf/.5d0,.5d0,0.d0/
    data bet/0.d0,.5d0,.5d0/

    ngauss=3    !ptos de gauss donde voy a intergrar

	!$omp parallel &
	!$omp private(ielem,ipoin,gama,gm,temp,fmu,nx,ny,ux,uy,ar,hlong,hlongx,hlongy,tau,alfa_mu,i,j,&
	!$omp rn1,rn2,rn3,u_loc,phi_loc,vx,vy,rmod2,et,c,a1,a2,aux,aa,aux_phi,choq,arr,rhs_temp)

	!$omp do
    do ielem=1,nelem
		rhs_temp = 0.d0
		ipoin = n(:,ielem)

        gama=(gamm(ipoin(1))+gamm(ipoin(2))+gamm(ipoin(3)))/3.d0
        gm=gama-1.d0
        temp=(t(ipoin(1))+t(ipoin(2))+t(ipoin(3)))/3.d0
        fmu= 1.716d-5*162.6/(temp-110.55)*(temp/273.15)**.75d0     !sutherland

		nx = dnx(:,ielem)
		ny = dny(:,ielem)

		ux(:) = u(:,ipoin(1))*nx(1) + u(:,ipoin(2))*nx(2) + u(:,ipoin(3))*nx(3)
		uy(:) = u(:,ipoin(1))*ny(1) + u(:,ipoin(2))*ny(2) + u(:,ipoin(3))*ny(3)

        ar=area(ielem)*dtl(ielem)/3.d0
        !cccc  ----> long. caracteristica
        hlongx=hhx(ielem)
        hlongy=hhy(ielem)
        hlong=dsqrt(area(ielem))
        !cccc  ----> estab. tezduyar
        tau(1)=t_sugn1(ielem)
        tau(2)=t_sugn2(ielem)
        tau(3)=t_sugn2(ielem)
        tau(4)=t_sugn3(ielem)
        alfa_mu=shoc(ielem)

        do j=1,ngauss

            rn1=1.d0-alf(j)-bet(j)
            rn2=alf(j)
            rn3=bet(j)

            !cccc  ----> integro las variables en los puntos de gauss
			u_loc = rn1*u(:,ipoin(1)) + rn2*u(:,ipoin(2)) + rn3*u(:,ipoin(3))
			phi_loc = rn1*un(:,ipoin(1)) + rn2*un(:,ipoin(2)) + rn3*un(:,ipoin(3))		

            !cccc  ----> defino variables primitivas
            vx=u_loc(2)/u_loc(1)
            vy=u_loc(3)/u_loc(1)
            et=u_loc(4)/u_loc(1)
            rmod2=vx*vx+vy*vy
            temp=gm/fr*(et-.5d0*rmod2) !fr=cte. universal de los gases

            c=dsqrt(gama*fr*temp)

            !cccc  ----> definicion de las matrices a1 y a2
            !a1(:,1) = (/ 0.d0				  , 1.d0				          , 0.d0		  , 0.d0                       /)
            a1(1,:) = (/ gm/2.d0*rmod2-vx*vx  , (3.d0-gama)*vx                , -gm*vy        , gm                         /)
            a1(2,:) = (/ -vx*vy 			  , vy 			                  , vx	          , 0.d0                       /)
            a1(3,:) = (/ (gm*rmod2-gama*et)*vx, gama*et-gm/2.d0*rmod2-gm*vx*vx, -gm*vx*vy     , gama*vx                    /)

            !a2(:,1) = (/ 0.d0				  , 0.d0						  , 1.d0	      , 0.d0                       /)
            a2(1,:) = (/ -vx*vy               , vy                            , vx            , 0.d0 			           /)
            a2(2,:) = (/ gm/2.d0*rmod2-vy*vy  , -gm*vx                        , (3.d0-gama)*vy, gm  		               /)
            a2(3,:) = (/ (gm*rmod2-gama*et)*vy, -gm*vx*vy                     , gama*et-gm/2.d0*rmod2-gm*vy*vy,    gama*vy /)

            !cccc----> terminos de estabilizacion y captura de choque
			arr = ar*tau

            !cccc  ----> armo el pressure switch
            choq(1)=alfa_mu*ar*ps(ipoin(1))*cte
            choq(2)=alfa_mu*ar*ps(ipoin(2))*cte
            choq(3)=alfa_mu*ar*ps(ipoin(3))*cte

            !cccc  ----> multiplico por partes para simplificar el asunto
            !cccc  ----> 'a' por las derivadas
        	aux(1)=                  ux(2)                                       +                    uy(3)
			aux(2:4) = a1(:,1)*ux(1) + a1(:,2)*ux(2) + a1(:,3)*ux(3) + a1(:,4)*ux(4) + a2(:,1)*uy(1) + a2(:,2)*uy(2)&
				+ a2(:,3)*uy(3) + a2(:,4)*uy(4)

			aux_phi = aux + phi_loc

            !cccc  ----> lo anterior por 'a' transpuesta
			aa(1) = aux_phi(2)
			aa(2:4) = a1(:,1)*aux_phi(1) + a1(:,2)*aux_phi(2) + a1(:,3)*aux_phi(3) + a1(:,4)*aux_phi(4)
			aa(5) = aux_phi(3)
			aa(6:8) = a2(:,1)*aux_phi(1) + a2(:,2)*aux_phi(2) + a2(:,3)*aux_phi(3) + a2(:,4)*aux_phi(4)

            !cccc  ----> ensamble del right hand side del nodo ipoin(1)
			rhs_temp(:,1) = rhs_temp(:,1) + (nx(1)*aa(1:4) + ny(1)*aa(5:8))*arr(:) + &
					aux(:)*rn1*ar + (nx(1)*ux(:) + ny(1)*uy(:))*choq(1)
			rhs_temp(:,2) = rhs_temp(:,2) + (nx(2)*aa(1:4) + ny(2)*aa(5:8))*arr(:) + &
					aux(:)*rn2*ar + (nx(2)*ux(:) + ny(2)*uy(:))*choq(2)
			rhs_temp(:,3) = rhs_temp(:,3) + (nx(3)*aa(1:4) + ny(3)*aa(5:8))*arr(:) + &
					aux(:)*rn3*ar + (nx(3)*ux(:) + ny(3)*uy(:))*choq(3)
        enddo
		do i=1,4
			!$omp atomic
			rhs(i,ipoin(1)) = rhs(i,ipoin(1)) + rhs_temp(i,1)
			!$omp atomic
			rhs(i,ipoin(2)) = rhs(i,ipoin(2)) + rhs_temp(i,2)
			!$omp atomic
			rhs(i,ipoin(3)) = rhs(i,ipoin(3)) + rhs_temp(i,3)
		end do
    enddo
	!$omp end do
	!$omp end parallel
    return
end subroutine todo_mp_atomic_only_outside

subroutine todo_mp_critical(nnod,nelem,n,dnx,dny,area,hhx,hhy &
    ,t_sugn1,t_sugn2,t_sugn3,shoc,ps,dtl &
    ,u,un,rhs,p,gamm,fr,rmu,fk,fcv,tinf,cte)
	use omp_lib
    implicit real(8) (a-h,o-z)
    integer n(3,nelem), ipoin(3)
    real(8) dnx(3,nelem),dny(3,nelem),u(4,nnod),p(nnod),un(4,nnod),t(nnod),gamm(nnod)
    real(8) area(nelem),rhs(4,nnod),hhx(nelem),hhy(nelem)
	real(8) a1(3,4), a2(3,4), aux(4), aa(8), aux_phi(4), nx(3), ny(3), ux(4), uy(4), u_loc(4)
	real(8) arr(4), tau(4), choq(3), phi_loc(4)
	real(8), dimension(:,:), allocatable :: prhs
    real(8) alf(3),bet(3)
    real(8) t_sugn1(nelem),t_sugn2(nelem),t_sugn3(nelem),shoc(nelem)
    real(8) ps(nnod),dtl(nelem)

    data alf/.5d0,.5d0,0.d0/
    data bet/0.d0,.5d0,.5d0/
	allocate(prhs(4,nnod))

    ngauss=3    !ptos de gauss donde voy a intergrar

	!$omp parallel &
	!$omp private(ielem,ipoin,gama,gm,temp,fmu,nx,ny,ux,uy,ar,hlong,hlongx,hlongy,tau,alfa_mu,i,j,&
	!$omp rn1,rn2,rn3,u_loc,phi_loc,vx,vy,rmod2,et,c,a1,a2,aux,aa,aux_phi,choq,arr,prhs)
	prhs=0.d0

	!$omp do
    do ielem=1,nelem
		ipoin = n(:,ielem)

        gama=(gamm(ipoin(1))+gamm(ipoin(2))+gamm(ipoin(3)))/3.d0
        gm=gama-1.d0
        temp=(t(ipoin(1))+t(ipoin(2))+t(ipoin(3)))/3.d0
        fmu= 1.716d-5*162.6/(temp-110.55)*(temp/273.15)**.75d0     !sutherland

		nx = dnx(:,ielem)
		ny = dny(:,ielem)

		ux(:) = u(:,ipoin(1))*nx(1) + u(:,ipoin(2))*nx(2) + u(:,ipoin(3))*nx(3)
		uy(:) = u(:,ipoin(1))*ny(1) + u(:,ipoin(2))*ny(2) + u(:,ipoin(3))*ny(3)

        ar=area(ielem)*dtl(ielem)/3.d0
        !cccc  ----> long. caracteristica
        hlongx=hhx(ielem)
        hlongy=hhy(ielem)
        hlong=dsqrt(area(ielem))
        !cccc  ----> estab. tezduyar
        tau(1)=t_sugn1(ielem)
        tau(2)=t_sugn2(ielem)
        tau(3)=t_sugn2(ielem)
        tau(4)=t_sugn3(ielem)
        alfa_mu=shoc(ielem)

        do j=1,ngauss

            rn1=1.d0-alf(j)-bet(j)
            rn2=alf(j)
            rn3=bet(j)

            !cccc  ----> integro las variables en los puntos de gauss
			u_loc = rn1*u(:,ipoin(1)) + rn2*u(:,ipoin(2)) + rn3*u(:,ipoin(3))
			phi_loc = rn1*un(:,ipoin(1)) + rn2*un(:,ipoin(2)) + rn3*un(:,ipoin(3))		

            !cccc  ----> defino variables primitivas
            vx=u_loc(2)/u_loc(1)
            vy=u_loc(3)/u_loc(1)
            et=u_loc(4)/u_loc(1)
            rmod2=vx*vx+vy*vy
            temp=gm/fr*(et-.5d0*rmod2) !fr=cte. universal de los gases

            c=dsqrt(gama*fr*temp)

            !cccc  ----> definicion de las matrices a1 y a2
            !a1(:,1) = (/ 0.d0				  , 1.d0				          , 0.d0		  , 0.d0                       /)
            a1(1,:) = (/ gm/2.d0*rmod2-vx*vx  , (3.d0-gama)*vx                , -gm*vy        , gm                         /)
            a1(2,:) = (/ -vx*vy 			  , vy 			                  , vx	          , 0.d0                       /)
            a1(3,:) = (/ (gm*rmod2-gama*et)*vx, gama*et-gm/2.d0*rmod2-gm*vx*vx, -gm*vx*vy     , gama*vx                    /)

            !a2(:,1) = (/ 0.d0				  , 0.d0						  , 1.d0	      , 0.d0                       /)
            a2(1,:) = (/ -vx*vy               , vy                            , vx            , 0.d0 			           /)
            a2(2,:) = (/ gm/2.d0*rmod2-vy*vy  , -gm*vx                        , (3.d0-gama)*vy, gm  		               /)
            a2(3,:) = (/ (gm*rmod2-gama*et)*vy, -gm*vx*vy                     , gama*et-gm/2.d0*rmod2-gm*vy*vy,    gama*vy /)

            !cccc----> terminos de estabilizacion y captura de choque
			arr = ar*tau

            !cccc  ----> armo el pressure switch
            choq(1)=alfa_mu*ar*ps(ipoin(1))*cte
            choq(2)=alfa_mu*ar*ps(ipoin(2))*cte
            choq(3)=alfa_mu*ar*ps(ipoin(3))*cte

            !cccc  ----> multiplico por partes para simplificar el asunto
            !cccc  ----> 'a' por las derivadas
        	aux(1)=                  ux(2)                                       +                    uy(3)
			aux(2:4) = a1(:,1)*ux(1) + a1(:,2)*ux(2) + a1(:,3)*ux(3) + a1(:,4)*ux(4) + a2(:,1)*uy(1) + a2(:,2)*uy(2)&
				+ a2(:,3)*uy(3) + a2(:,4)*uy(4)

			aux_phi = aux + phi_loc

            !cccc  ----> lo anterior por 'a' transpuesta
			aa(1) = aux_phi(2)
			aa(2:4) = a1(:,1)*aux_phi(1) + a1(:,2)*aux_phi(2) + a1(:,3)*aux_phi(3) + a1(:,4)*aux_phi(4)
			aa(5) = aux_phi(3)
			aa(6:8) = a2(:,1)*aux_phi(1) + a2(:,2)*aux_phi(2) + a2(:,3)*aux_phi(3) + a2(:,4)*aux_phi(4)

            !cccc  ----> ensamble del right hand side del nodo ipoin(1)
			prhs(:,ipoin(1)) = prhs(:,ipoin(1)) + (nx(1)*aa(1:4) + ny(1)*aa(5:8))*arr(:) + &
					aux(:)*rn1*ar + (nx(1)*ux(:) + ny(1)*uy(:))*choq(1)
			prhs(:,ipoin(2)) = prhs(:,ipoin(2)) + (nx(2)*aa(1:4) + ny(2)*aa(5:8))*arr(:) + &
					aux(:)*rn2*ar + (nx(2)*ux(:) + ny(2)*uy(:))*choq(2)
			prhs(:,ipoin(3)) = prhs(:,ipoin(3)) + (nx(3)*aa(1:4) + ny(3)*aa(5:8))*arr(:) + &
					aux(:)*rn3*ar + (nx(3)*ux(:) + ny(3)*uy(:))*choq(3)
        enddo
    enddo
	!$omp end do

	!$omp critical
			rhs = rhs + prhs
	!$omp end critical

	!$omp end parallel
    return
end subroutine todo_mp_critical

subroutine todo_mp_critical_only(nnod,nelem,n,dnx,dny,area,hhx,hhy &
    ,t_sugn1,t_sugn2,t_sugn3,shoc,ps,dtl &
    ,u,un,rhs,p,gamm,fr,rmu,fk,fcv,tinf,cte)
	use omp_lib
    implicit real(8) (a-h,o-z)
    integer n(3,nelem), ipoin(3)
    real(8) dnx(3,nelem),dny(3,nelem),u(4,nnod),p(nnod),un(4,nnod),t(nnod),gamm(nnod)
    real(8) area(nelem),rhs(4,nnod),hhx(nelem),hhy(nelem)
	real(8) a1(3,4), a2(3,4), aux(4), aa(8), aux_phi(4), nx(3), ny(3), ux(4), uy(4), u_loc(4)
	real(8) arr(4), tau(4), choq(3), phi_loc(4)
    real(8) alf(3),bet(3)
    real(8) t_sugn1(nelem),t_sugn2(nelem),t_sugn3(nelem),shoc(nelem)
    real(8) ps(nnod),dtl(nelem)

    data alf/.5d0,.5d0,0.d0/
    data bet/0.d0,.5d0,.5d0/

    ngauss=3    !ptos de gauss donde voy a intergrar

	!$omp parallel &
	!$omp private(ielem,ipoin,gama,gm,temp,fmu,nx,ny,ux,uy,ar,hlong,hlongx,hlongy,tau,alfa_mu,i,j,&
	!$omp rn1,rn2,rn3,u_loc,phi_loc,vx,vy,rmod2,et,c,a1,a2,aux,aa,aux_phi,choq,arr)

	!$omp do
    do ielem=1,nelem
		ipoin = n(:,ielem)

        gama=(gamm(ipoin(1))+gamm(ipoin(2))+gamm(ipoin(3)))/3.d0
        gm=gama-1.d0
        temp=(t(ipoin(1))+t(ipoin(2))+t(ipoin(3)))/3.d0
        fmu= 1.716d-5*162.6/(temp-110.55)*(temp/273.15)**.75d0     !sutherland

		nx = dnx(:,ielem)
		ny = dny(:,ielem)

		ux(:) = u(:,ipoin(1))*nx(1) + u(:,ipoin(2))*nx(2) + u(:,ipoin(3))*nx(3)
		uy(:) = u(:,ipoin(1))*ny(1) + u(:,ipoin(2))*ny(2) + u(:,ipoin(3))*ny(3)

        ar=area(ielem)*dtl(ielem)/3.d0
        !cccc  ----> long. caracteristica
        hlongx=hhx(ielem)
        hlongy=hhy(ielem)
        hlong=dsqrt(area(ielem))
        !cccc  ----> estab. tezduyar
        tau(1)=t_sugn1(ielem)
        tau(2)=t_sugn2(ielem)
        tau(3)=t_sugn2(ielem)
        tau(4)=t_sugn3(ielem)
        alfa_mu=shoc(ielem)

        do j=1,ngauss

            rn1=1.d0-alf(j)-bet(j)
            rn2=alf(j)
            rn3=bet(j)

            !cccc  ----> integro las variables en los puntos de gauss
			u_loc = rn1*u(:,ipoin(1)) + rn2*u(:,ipoin(2)) + rn3*u(:,ipoin(3))
			phi_loc = rn1*un(:,ipoin(1)) + rn2*un(:,ipoin(2)) + rn3*un(:,ipoin(3))		

            !cccc  ----> defino variables primitivas
            vx=u_loc(2)/u_loc(1)
            vy=u_loc(3)/u_loc(1)
            et=u_loc(4)/u_loc(1)
            rmod2=vx*vx+vy*vy
            temp=gm/fr*(et-.5d0*rmod2) !fr=cte. universal de los gases

            c=dsqrt(gama*fr*temp)

            !cccc  ----> definicion de las matrices a1 y a2
            !a1(:,1) = (/ 0.d0				  , 1.d0				          , 0.d0		  , 0.d0                       /)
            a1(1,:) = (/ gm/2.d0*rmod2-vx*vx  , (3.d0-gama)*vx                , -gm*vy        , gm                         /)
            a1(2,:) = (/ -vx*vy 			  , vy 			                  , vx	          , 0.d0                       /)
            a1(3,:) = (/ (gm*rmod2-gama*et)*vx, gama*et-gm/2.d0*rmod2-gm*vx*vx, -gm*vx*vy     , gama*vx                    /)

            !a2(:,1) = (/ 0.d0				  , 0.d0						  , 1.d0	      , 0.d0                       /)
            a2(1,:) = (/ -vx*vy               , vy                            , vx            , 0.d0 			           /)
            a2(2,:) = (/ gm/2.d0*rmod2-vy*vy  , -gm*vx                        , (3.d0-gama)*vy, gm  		               /)
            a2(3,:) = (/ (gm*rmod2-gama*et)*vy, -gm*vx*vy                     , gama*et-gm/2.d0*rmod2-gm*vy*vy,    gama*vy /)

            !cccc----> terminos de estabilizacion y captura de choque
			arr = ar*tau

            !cccc  ----> armo el pressure switch
            choq(1)=alfa_mu*ar*ps(ipoin(1))*cte
            choq(2)=alfa_mu*ar*ps(ipoin(2))*cte
            choq(3)=alfa_mu*ar*ps(ipoin(3))*cte

            !cccc  ----> multiplico por partes para simplificar el asunto
            !cccc  ----> 'a' por las derivadas
        	aux(1)=                  ux(2)                                       +                    uy(3)
			aux(2:4) = a1(:,1)*ux(1) + a1(:,2)*ux(2) + a1(:,3)*ux(3) + a1(:,4)*ux(4) + a2(:,1)*uy(1) + a2(:,2)*uy(2)&
				+ a2(:,3)*uy(3) + a2(:,4)*uy(4)

			aux_phi = aux + phi_loc

            !cccc  ----> lo anterior por 'a' transpuesta
			aa(1) = aux_phi(2)
			aa(2:4) = a1(:,1)*aux_phi(1) + a1(:,2)*aux_phi(2) + a1(:,3)*aux_phi(3) + a1(:,4)*aux_phi(4)
			aa(5) = aux_phi(3)
			aa(6:8) = a2(:,1)*aux_phi(1) + a2(:,2)*aux_phi(2) + a2(:,3)*aux_phi(3) + a2(:,4)*aux_phi(4)

            !cccc  ----> ensamble del right hand side del nodo ipoin(1)
	!$omp critical
			rhs(:,ipoin(1)) = rhs(:,ipoin(1)) + (nx(1)*aa(1:4) + ny(1)*aa(5:8))*arr(:) + &
					aux(:)*rn1*ar + (nx(1)*ux(:) + ny(1)*uy(:))*choq(1)
			rhs(:,ipoin(2)) = rhs(:,ipoin(2)) + (nx(2)*aa(1:4) + ny(2)*aa(5:8))*arr(:) + &
					aux(:)*rn2*ar + (nx(2)*ux(:) + ny(2)*uy(:))*choq(2)
			rhs(:,ipoin(3)) = rhs(:,ipoin(3)) + (nx(3)*aa(1:4) + ny(3)*aa(5:8))*arr(:) + &
					aux(:)*rn3*ar + (nx(3)*ux(:) + ny(3)*uy(:))*choq(3)
	!$omp end critical
        enddo
    enddo
	!$omp end do
	!$omp end parallel
    return
end subroutine todo_mp_critical_only

subroutine todo_mp_critical_only_outside(nnod,nelem,n,dnx,dny,area,hhx,hhy &
    ,t_sugn1,t_sugn2,t_sugn3,shoc,ps,dtl &
    ,u,un,rhs,p,gamm,fr,rmu,fk,fcv,tinf,cte)
	use omp_lib
    implicit real(8) (a-h,o-z)
    integer n(3,nelem), ipoin(3)
    real(8) dnx(3,nelem),dny(3,nelem),u(4,nnod),p(nnod),un(4,nnod),t(nnod),gamm(nnod)
    real(8) area(nelem),rhs(4,nnod),hhx(nelem),hhy(nelem)
	real(8) a1(3,4), a2(3,4), aux(4), aa(8), aux_phi(4), nx(3), ny(3), ux(4), uy(4), u_loc(4)
	real(8) arr(4), tau(4), choq(3), phi_loc(4), rhs_temp(4,3)
    real(8) alf(3),bet(3)
    real(8) t_sugn1(nelem),t_sugn2(nelem),t_sugn3(nelem),shoc(nelem)
    real(8) ps(nnod),dtl(nelem)

    data alf/.5d0,.5d0,0.d0/
    data bet/0.d0,.5d0,.5d0/

    ngauss=3    !ptos de gauss donde voy a intergrar

	!$omp parallel &
	!$omp private(ielem,ipoin,gama,gm,temp,fmu,nx,ny,ux,uy,ar,hlong,hlongx,hlongy,tau,alfa_mu,i,j,&
	!$omp rn1,rn2,rn3,u_loc,phi_loc,vx,vy,rmod2,et,c,a1,a2,aux,aa,aux_phi,choq,arr,rhs_temp)

	!$omp do
    do ielem=1,nelem
		rhs_temp = 0.d0
		ipoin = n(:,ielem)

        gama=(gamm(ipoin(1))+gamm(ipoin(2))+gamm(ipoin(3)))/3.d0
        gm=gama-1.d0
        temp=(t(ipoin(1))+t(ipoin(2))+t(ipoin(3)))/3.d0
        fmu= 1.716d-5*162.6/(temp-110.55)*(temp/273.15)**.75d0     !sutherland

		nx = dnx(:,ielem)
		ny = dny(:,ielem)

		ux(:) = u(:,ipoin(1))*nx(1) + u(:,ipoin(2))*nx(2) + u(:,ipoin(3))*nx(3)
		uy(:) = u(:,ipoin(1))*ny(1) + u(:,ipoin(2))*ny(2) + u(:,ipoin(3))*ny(3)

        ar=area(ielem)*dtl(ielem)/3.d0
        !cccc  ----> long. caracteristica
        hlongx=hhx(ielem)
        hlongy=hhy(ielem)
        hlong=dsqrt(area(ielem))
        !cccc  ----> estab. tezduyar
        tau(1)=t_sugn1(ielem)
        tau(2)=t_sugn2(ielem)
        tau(3)=t_sugn2(ielem)
        tau(4)=t_sugn3(ielem)
        alfa_mu=shoc(ielem)

        do j=1,ngauss

            rn1=1.d0-alf(j)-bet(j)
            rn2=alf(j)
            rn3=bet(j)

            !cccc  ----> integro las variables en los puntos de gauss
			u_loc = rn1*u(:,ipoin(1)) + rn2*u(:,ipoin(2)) + rn3*u(:,ipoin(3))
			phi_loc = rn1*un(:,ipoin(1)) + rn2*un(:,ipoin(2)) + rn3*un(:,ipoin(3))		

            !cccc  ----> defino variables primitivas
            vx=u_loc(2)/u_loc(1)
            vy=u_loc(3)/u_loc(1)
            et=u_loc(4)/u_loc(1)
            rmod2=vx*vx+vy*vy
            temp=gm/fr*(et-.5d0*rmod2) !fr=cte. universal de los gases

            c=dsqrt(gama*fr*temp)

            !cccc  ----> definicion de las matrices a1 y a2
            !a1(:,1) = (/ 0.d0				  , 1.d0				          , 0.d0		  , 0.d0                       /)
            a1(1,:) = (/ gm/2.d0*rmod2-vx*vx  , (3.d0-gama)*vx                , -gm*vy        , gm                         /)
            a1(2,:) = (/ -vx*vy 			  , vy 			                  , vx	          , 0.d0                       /)
            a1(3,:) = (/ (gm*rmod2-gama*et)*vx, gama*et-gm/2.d0*rmod2-gm*vx*vx, -gm*vx*vy     , gama*vx                    /)

            !a2(:,1) = (/ 0.d0				  , 0.d0						  , 1.d0	      , 0.d0                       /)
            a2(1,:) = (/ -vx*vy               , vy                            , vx            , 0.d0 			           /)
            a2(2,:) = (/ gm/2.d0*rmod2-vy*vy  , -gm*vx                        , (3.d0-gama)*vy, gm  		               /)
            a2(3,:) = (/ (gm*rmod2-gama*et)*vy, -gm*vx*vy                     , gama*et-gm/2.d0*rmod2-gm*vy*vy,    gama*vy /)

            !cccc----> terminos de estabilizacion y captura de choque
			arr = ar*tau

            !cccc  ----> armo el pressure switch
            choq(1)=alfa_mu*ar*ps(ipoin(1))*cte
            choq(2)=alfa_mu*ar*ps(ipoin(2))*cte
            choq(3)=alfa_mu*ar*ps(ipoin(3))*cte

            !cccc  ----> multiplico por partes para simplificar el asunto
            !cccc  ----> 'a' por las derivadas
        	aux(1)=                  ux(2)                                       +                    uy(3)
			aux(2:4) = a1(:,1)*ux(1) + a1(:,2)*ux(2) + a1(:,3)*ux(3) + a1(:,4)*ux(4) + a2(:,1)*uy(1) + a2(:,2)*uy(2)&
				+ a2(:,3)*uy(3) + a2(:,4)*uy(4)

			aux_phi = aux + phi_loc

            !cccc  ----> lo anterior por 'a' transpuesta
			aa(1) = aux_phi(2)
			aa(2:4) = a1(:,1)*aux_phi(1) + a1(:,2)*aux_phi(2) + a1(:,3)*aux_phi(3) + a1(:,4)*aux_phi(4)
			aa(5) = aux_phi(3)
			aa(6:8) = a2(:,1)*aux_phi(1) + a2(:,2)*aux_phi(2) + a2(:,3)*aux_phi(3) + a2(:,4)*aux_phi(4)

            !cccc  ----> ensamble del right hand side del nodo ipoin(1)
			rhs_temp(:,1) = rhs_temp(:,1) + (nx(1)*aa(1:4) + ny(1)*aa(5:8))*arr(:) + &
					aux(:)*rn1*ar + (nx(1)*ux(:) + ny(1)*uy(:))*choq(1)
			rhs_temp(:,2) = rhs_temp(:,2) + (nx(2)*aa(1:4) + ny(2)*aa(5:8))*arr(:) + &
					aux(:)*rn2*ar + (nx(2)*ux(:) + ny(2)*uy(:))*choq(2)
			rhs_temp(:,3) = rhs_temp(:,3) + (nx(3)*aa(1:4) + ny(3)*aa(5:8))*arr(:) + &
					aux(:)*rn3*ar + (nx(3)*ux(:) + ny(3)*uy(:))*choq(3)
        enddo
	!$omp critical
	rhs(:,ipoin(1)) = rhs(:,ipoin(1)) + rhs_temp(:,1)
	rhs(:,ipoin(2)) = rhs(:,ipoin(2)) + rhs_temp(:,2)
	rhs(:,ipoin(3)) = rhs(:,ipoin(3)) + rhs_temp(:,3)
	!$omp end critical
    enddo
	!$omp end do
	!$omp end parallel
    return
end subroutine todo_mp_critical_only_outside

subroutine todo_mp_critical_fp(nnod,nelem,n,dnx,dny,area,hhx,hhy &
    ,t_sugn1,t_sugn2,t_sugn3,shoc,ps,dtl &
    ,u,un,rhs,p,gamm,fr,rmu,fk,fcv,tinf,cte)
	use omp_lib
    implicit real(8) (a-h,o-z)
    integer n(3,nelem), ipoin(3)
    real(8) dnx(3,nelem),dny(3,nelem),u(4,nnod),p(nnod),un(4,nnod),t(nnod),gamm(nnod)
    real(8) area(nelem),rhs(4,nnod),hhx(nelem),hhy(nelem)
	real(8) a1(3,4), a2(3,4), aux(4), aa(8), aux_phi(4), nx(3), ny(3), ux(4), uy(4), u_loc(4)
	real(8) arr(4), tau(4), choq(3), phi_loc(4)
	real(8), dimension(:,:), allocatable :: prhs
    real(8) alf(3),bet(3)
    real(8) t_sugn1(nelem),t_sugn2(nelem),t_sugn3(nelem),shoc(nelem)
    real(8) ps(nnod),dtl(nelem)

    data alf/.5d0,.5d0,0.d0/
    data bet/0.d0,.5d0,.5d0/

	allocate(prhs(4,nnod))
	prhs=0.d0

    ngauss=3    !ptos de gauss donde voy a intergrar

	!$omp parallel &
	!$omp firstprivate(ielem,ipoin,gama,gm,temp,fmu,nx,ny,ux,uy,ar,hlong,hlongx,hlongy,tau,alfa_mu,i,j,&
	!$omp rn1,rn2,rn3,u_loc,phi_loc,vx,vy,rmod2,et,c,a1,a2,aux,aa,aux_phi,choq,arr,prhs)

	!$omp do
    do ielem=1,nelem
		ipoin = n(:,ielem)

        gama=(gamm(ipoin(1))+gamm(ipoin(2))+gamm(ipoin(3)))/3.d0
        gm=gama-1.d0
        temp=(t(ipoin(1))+t(ipoin(2))+t(ipoin(3)))/3.d0
        fmu= 1.716d-5*162.6/(temp-110.55)*(temp/273.15)**.75d0     !sutherland

		nx = dnx(:,ielem)
		ny = dny(:,ielem)

		ux(:) = u(:,ipoin(1))*nx(1) + u(:,ipoin(2))*nx(2) + u(:,ipoin(3))*nx(3)
		uy(:) = u(:,ipoin(1))*ny(1) + u(:,ipoin(2))*ny(2) + u(:,ipoin(3))*ny(3)

        ar=area(ielem)*dtl(ielem)/3.d0
        !cccc  ----> long. caracteristica
        hlongx=hhx(ielem)
        hlongy=hhy(ielem)
        hlong=dsqrt(area(ielem))
        !cccc  ----> estab. tezduyar
        tau(1)=t_sugn1(ielem)
        tau(2)=t_sugn2(ielem)
        tau(3)=t_sugn2(ielem)
        tau(4)=t_sugn3(ielem)
        alfa_mu=shoc(ielem)

        do j=1,ngauss

            rn1=1.d0-alf(j)-bet(j)
            rn2=alf(j)
            rn3=bet(j)

            !cccc  ----> integro las variables en los puntos de gauss
			u_loc = rn1*u(:,ipoin(1)) + rn2*u(:,ipoin(2)) + rn3*u(:,ipoin(3))
			phi_loc = rn1*un(:,ipoin(1)) + rn2*un(:,ipoin(2)) + rn3*un(:,ipoin(3))		

            !cccc  ----> defino variables primitivas
            vx=u_loc(2)/u_loc(1)
            vy=u_loc(3)/u_loc(1)
            et=u_loc(4)/u_loc(1)
            rmod2=vx*vx+vy*vy
            temp=gm/fr*(et-.5d0*rmod2) !fr=cte. universal de los gases

            c=dsqrt(gama*fr*temp)

            !cccc  ----> definicion de las matrices a1 y a2
            !a1(:,1) = (/ 0.d0				  , 1.d0				          , 0.d0		  , 0.d0                       /)
            a1(1,:) = (/ gm/2.d0*rmod2-vx*vx  , (3.d0-gama)*vx                , -gm*vy        , gm                         /)
            a1(2,:) = (/ -vx*vy 			  , vy 			                  , vx	          , 0.d0                       /)
            a1(3,:) = (/ (gm*rmod2-gama*et)*vx, gama*et-gm/2.d0*rmod2-gm*vx*vx, -gm*vx*vy     , gama*vx                    /)

            !a2(:,1) = (/ 0.d0				  , 0.d0						  , 1.d0	      , 0.d0                       /)
            a2(1,:) = (/ -vx*vy               , vy                            , vx            , 0.d0 			           /)
            a2(2,:) = (/ gm/2.d0*rmod2-vy*vy  , -gm*vx                        , (3.d0-gama)*vy, gm  		               /)
            a2(3,:) = (/ (gm*rmod2-gama*et)*vy, -gm*vx*vy                     , gama*et-gm/2.d0*rmod2-gm*vy*vy,    gama*vy /)

            !cccc----> terminos de estabilizacion y captura de choque
			arr = ar*tau

            !cccc  ----> armo el pressure switch
            choq(1)=alfa_mu*ar*ps(ipoin(1))*cte
            choq(2)=alfa_mu*ar*ps(ipoin(2))*cte
            choq(3)=alfa_mu*ar*ps(ipoin(3))*cte

            !cccc  ----> multiplico por partes para simplificar el asunto
            !cccc  ----> 'a' por las derivadas
        	aux(1)=                  ux(2)                                       +                    uy(3)
			aux(2:4) = a1(:,1)*ux(1) + a1(:,2)*ux(2) + a1(:,3)*ux(3) + a1(:,4)*ux(4) + a2(:,1)*uy(1) + a2(:,2)*uy(2)&
				+ a2(:,3)*uy(3) + a2(:,4)*uy(4)

			aux_phi = aux + phi_loc

            !cccc  ----> lo anterior por 'a' transpuesta
			aa(1) = aux_phi(2)
			aa(2:4) = a1(:,1)*aux_phi(1) + a1(:,2)*aux_phi(2) + a1(:,3)*aux_phi(3) + a1(:,4)*aux_phi(4)
			aa(5) = aux_phi(3)
			aa(6:8) = a2(:,1)*aux_phi(1) + a2(:,2)*aux_phi(2) + a2(:,3)*aux_phi(3) + a2(:,4)*aux_phi(4)

            !cccc  ----> ensamble del right hand side del nodo ipoin(1)
			prhs(:,ipoin(1)) = prhs(:,ipoin(1)) + (nx(1)*aa(1:4) + ny(1)*aa(5:8))*arr(:) + &
					aux(:)*rn1*ar + (nx(1)*ux(:) + ny(1)*uy(:))*choq(1)
			prhs(:,ipoin(2)) = prhs(:,ipoin(2)) + (nx(2)*aa(1:4) + ny(2)*aa(5:8))*arr(:) + &
					aux(:)*rn2*ar + (nx(2)*ux(:) + ny(2)*uy(:))*choq(2)
			prhs(:,ipoin(3)) = prhs(:,ipoin(3)) + (nx(3)*aa(1:4) + ny(3)*aa(5:8))*arr(:) + &
					aux(:)*rn3*ar + (nx(3)*ux(:) + ny(3)*uy(:))*choq(3)
        enddo
    enddo
	!$omp end do

	!$omp critical
			rhs = rhs + prhs
	!$omp end critical
	!$omp end parallel

    return
end subroutine todo_mp_critical_fp

subroutine todo_mp_atomic(nnod,nelem,n,dnx,dny,area,hhx,hhy &
    ,t_sugn1,t_sugn2,t_sugn3,shoc,ps,dtl &
    ,u,un,rhs,p,gamm,fr,rmu,fk,fcv,tinf,cte)
	use omp_lib
    implicit real(8) (a-h,o-z)
    integer n(3,nelem), ipoin(3)
    real(8) dnx(3,nelem),dny(3,nelem),u(4,nnod),p(nnod),un(4,nnod),t(nnod),gamm(nnod)
    real(8) area(nelem),rhs(4,nnod),hhx(nelem),hhy(nelem)
	real(8) a1(3,4), a2(3,4), aux(4), aa(8), aux_phi(4), nx(3), ny(3), ux(4), uy(4), u_loc(4)
	real(8) arr(4), tau(4), choq(3), phi_loc(4)
	real(8), dimension(:,:), allocatable :: prhs
    real(8) alf(3),bet(3)
    real(8) t_sugn1(nelem),t_sugn2(nelem),t_sugn3(nelem),shoc(nelem)
    real(8) ps(nnod),dtl(nelem)

    data alf/.5d0,.5d0,0.d0/
    data bet/0.d0,.5d0,.5d0/

	allocate(prhs(4,nnod))

    ngauss=3    !ptos de gauss donde voy a intergrar

	!$omp parallel &
	!$omp private(ielem,ipoin,gama,gm,temp,fmu,nx,ny,ux,uy,ar,hlong,hlongx,hlongy,tau,alfa_mu,i,j,&
	!$omp rn1,rn2,rn3,u_loc,phi_loc,vx,vy,rmod2,et,c,a1,a2,aux,aa,aux_phi,choq,arr,prhs)

	prhs = 0.d0

	!$omp do
    do ielem=1,nelem
		ipoin = n(:,ielem)

        gama=(gamm(ipoin(1))+gamm(ipoin(2))+gamm(ipoin(3)))/3.d0
        gm=gama-1.d0
        temp=(t(ipoin(1))+t(ipoin(2))+t(ipoin(3)))/3.d0
        fmu= 1.716d-5*162.6/(temp-110.55)*(temp/273.15)**.75d0     !sutherland

		nx = dnx(:,ielem)
		ny = dny(:,ielem)

		ux(:) = u(:,ipoin(1))*nx(1) + u(:,ipoin(2))*nx(2) + u(:,ipoin(3))*nx(3)
		uy(:) = u(:,ipoin(1))*ny(1) + u(:,ipoin(2))*ny(2) + u(:,ipoin(3))*ny(3)

        ar=area(ielem)*dtl(ielem)/3.d0
        !cccc  ----> long. caracteristica
        hlongx=hhx(ielem)
        hlongy=hhy(ielem)
        hlong=dsqrt(area(ielem))
        !cccc  ----> estab. tezduyar
        tau(1)=t_sugn1(ielem)
        tau(2)=t_sugn2(ielem)
        tau(3)=t_sugn2(ielem)
        tau(4)=t_sugn3(ielem)
        alfa_mu=shoc(ielem)

        do j=1,ngauss

            rn1=1.d0-alf(j)-bet(j)
            rn2=alf(j)
            rn3=bet(j)

            !cccc  ----> integro las variables en los puntos de gauss
			u_loc = rn1*u(:,ipoin(1)) + rn2*u(:,ipoin(2)) + rn3*u(:,ipoin(3))
			phi_loc = rn1*un(:,ipoin(1)) + rn2*un(:,ipoin(2)) + rn3*un(:,ipoin(3))		

            !cccc  ----> defino variables primitivas
            vx=u_loc(2)/u_loc(1)
            vy=u_loc(3)/u_loc(1)
            et=u_loc(4)/u_loc(1)
            rmod2=vx*vx+vy*vy
            temp=gm/fr*(et-.5d0*rmod2) !fr=cte. universal de los gases

            c=dsqrt(gama*fr*temp)

            !cccc  ----> definicion de las matrices a1 y a2
            !a1(:,1) = (/ 0.d0				  , 1.d0				          , 0.d0		  , 0.d0                       /)
            a1(1,:) = (/ gm/2.d0*rmod2-vx*vx  , (3.d0-gama)*vx                , -gm*vy        , gm                         /)
            a1(2,:) = (/ -vx*vy 			  , vy 			                  , vx	          , 0.d0                       /)
            a1(3,:) = (/ (gm*rmod2-gama*et)*vx, gama*et-gm/2.d0*rmod2-gm*vx*vx, -gm*vx*vy     , gama*vx                    /)

            !a2(:,1) = (/ 0.d0				  , 0.d0						  , 1.d0	      , 0.d0                       /)
            a2(1,:) = (/ -vx*vy               , vy                            , vx            , 0.d0 			           /)
            a2(2,:) = (/ gm/2.d0*rmod2-vy*vy  , -gm*vx                        , (3.d0-gama)*vy, gm  		               /)
            a2(3,:) = (/ (gm*rmod2-gama*et)*vy, -gm*vx*vy                     , gama*et-gm/2.d0*rmod2-gm*vy*vy,    gama*vy /)

            !cccc----> terminos de estabilizacion y captura de choque
			arr = ar*tau

            !cccc  ----> armo el pressure switch
            choq(1)=alfa_mu*ar*ps(ipoin(1))*cte
            choq(2)=alfa_mu*ar*ps(ipoin(2))*cte
            choq(3)=alfa_mu*ar*ps(ipoin(3))*cte

            !cccc  ----> multiplico por partes para simplificar el asunto
            !cccc  ----> 'a' por las derivadas
        	aux(1)=                  ux(2)                                       +                    uy(3)
			aux(2:4) = a1(:,1)*ux(1) + a1(:,2)*ux(2) + a1(:,3)*ux(3) + a1(:,4)*ux(4) + a2(:,1)*uy(1) + a2(:,2)*uy(2)&
				+ a2(:,3)*uy(3) + a2(:,4)*uy(4)

			aux_phi = aux + phi_loc

            !cccc  ----> lo anterior por 'a' transpuesta
			aa(1) = aux_phi(2)
			aa(2:4) = a1(:,1)*aux_phi(1) + a1(:,2)*aux_phi(2) + a1(:,3)*aux_phi(3) + a1(:,4)*aux_phi(4)
			aa(5) = aux_phi(3)
			aa(6:8) = a2(:,1)*aux_phi(1) + a2(:,2)*aux_phi(2) + a2(:,3)*aux_phi(3) + a2(:,4)*aux_phi(4)

            !cccc  ----> ensamble del right hand side del nodo ipoin(1)
			prhs(:,ipoin(1)) = prhs(:,ipoin(1)) + (nx(1)*aa(1:4) + ny(1)*aa(5:8))*arr(:) + &
					aux(:)*rn1*ar + (nx(1)*ux(:) + ny(1)*uy(:))*choq(1)
			prhs(:,ipoin(2)) = prhs(:,ipoin(2)) + (nx(2)*aa(1:4) + ny(2)*aa(5:8))*arr(:) + &
					aux(:)*rn2*ar + (nx(2)*ux(:) + ny(2)*uy(:))*choq(2)
			prhs(:,ipoin(3)) = prhs(:,ipoin(3)) + (nx(3)*aa(1:4) + ny(3)*aa(5:8))*arr(:) + &
					aux(:)*rn3*ar + (nx(3)*ux(:) + ny(3)*uy(:))*choq(3)
        enddo
    enddo
	!$omp end do

	do i=1,nnod
		do j=1,4
	!$omp atomic
			rhs(j,i) = rhs(j,i) + prhs(j,i)
		end do
	end do
	!$omp end parallel

    return
end subroutine todo_mp_atomic

subroutine todo_mp_noAtomic(nnod,nelem,n,dnx,dny,area,hhx,hhy &
    ,t_sugn1,t_sugn2,t_sugn3,shoc,ps,dtl &
    ,u,un,rhs,p,gamm,fr,rmu,fk,fcv,tinf,cte)
	use omp_lib
    implicit real(8) (a-h,o-z)
    integer n(3,nelem), ipoin(3)
    real(8) dnx(3,nelem),dny(3,nelem),u(4,nnod),p(nnod),un(4,nnod),t(nnod),gamm(nnod)
    real(8) area(nelem),rhs(4,nnod),hhx(nelem),hhy(nelem)
	real(8) a1(3,4), a2(3,4), aux(4), aa(8), aux_phi(4), nx(3), ny(3), ux(4), uy(4), u_loc(4)
	real(8) arr(4), tau(4), choq(3), phi_loc(4)
	real(8), dimension(:,:,:), allocatable :: brhs
    real(8) alf(3),bet(3)
    real(8) t_sugn1(nelem),t_sugn2(nelem),t_sugn3(nelem),shoc(nelem)
    real(8) ps(nnod),dtl(nelem)

    data alf/.5d0,.5d0,0.d0/
    data bet/0.d0,.5d0,.5d0/

    ngauss=3    !ptos de gauss donde voy a intergrar
	allocate(brhs(4,3,nelem))
	brhs=0.d0

	!$omp parallel &
	!$omp private(ielem,ipoin,gama,gm,temp,fmu,nx,ny,ux,uy,ar,hlong,hlongx,hlongy,tau,alfa_mu,i,j,&
	!$omp rn1,rn2,rn3,u_loc,phi_loc,vx,vy,rmod2,et,c,a1,a2,aux,aa,aux_phi,choq,arr)

	!$omp do
    do ielem=1,nelem
		ipoin = n(:,ielem)

        gama=(gamm(ipoin(1))+gamm(ipoin(2))+gamm(ipoin(3)))/3.d0
        gm=gama-1.d0
        temp=(t(ipoin(1))+t(ipoin(2))+t(ipoin(3)))/3.d0
        fmu= 1.716d-5*162.6/(temp-110.55)*(temp/273.15)**.75d0     !sutherland

		nx = dnx(:,ielem)
		ny = dny(:,ielem)

		ux(:) = u(:,ipoin(1))*nx(1) + u(:,ipoin(2))*nx(2) + u(:,ipoin(3))*nx(3)
		uy(:) = u(:,ipoin(1))*ny(1) + u(:,ipoin(2))*ny(2) + u(:,ipoin(3))*ny(3)

        ar=area(ielem)*dtl(ielem)/3.d0
        !cccc  ----> long. caracteristica
        hlongx=hhx(ielem)
        hlongy=hhy(ielem)
        hlong=dsqrt(area(ielem))
        !cccc  ----> estab. tezduyar
        tau(1)=t_sugn1(ielem)
        tau(2)=t_sugn2(ielem)
        tau(3)=t_sugn2(ielem)
        tau(4)=t_sugn3(ielem)
        alfa_mu=shoc(ielem)

        do j=1,ngauss

            rn1=1.d0-alf(j)-bet(j)
            rn2=alf(j)
            rn3=bet(j)

            !cccc  ----> integro las variables en los puntos de gauss
			u_loc = rn1*u(:,ipoin(1)) + rn2*u(:,ipoin(2)) + rn3*u(:,ipoin(3))
			phi_loc = rn1*un(:,ipoin(1)) + rn2*un(:,ipoin(2)) + rn3*un(:,ipoin(3))		

            !cccc  ----> defino variables primitivas
            vx=u_loc(2)/u_loc(1)
            vy=u_loc(3)/u_loc(1)
            et=u_loc(4)/u_loc(1)
            rmod2=vx*vx+vy*vy
            temp=gm/fr*(et-.5d0*rmod2) !fr=cte. universal de los gases

            c=dsqrt(gama*fr*temp)

            !cccc  ----> definicion de las matrices a1 y a2
            !a1(:,1) = (/ 0.d0				  , 1.d0				          , 0.d0		  , 0.d0                       /)
            a1(1,:) = (/ gm/2.d0*rmod2-vx*vx  , (3.d0-gama)*vx                , -gm*vy        , gm                         /)
            a1(2,:) = (/ -vx*vy 			  , vy 			                  , vx	          , 0.d0                       /)
            a1(3,:) = (/ (gm*rmod2-gama*et)*vx, gama*et-gm/2.d0*rmod2-gm*vx*vx, -gm*vx*vy     , gama*vx                    /)

            !a2(:,1) = (/ 0.d0				  , 0.d0						  , 1.d0	      , 0.d0                       /)
            a2(1,:) = (/ -vx*vy               , vy                            , vx            , 0.d0 			           /)
            a2(2,:) = (/ gm/2.d0*rmod2-vy*vy  , -gm*vx                        , (3.d0-gama)*vy, gm  		               /)
            a2(3,:) = (/ (gm*rmod2-gama*et)*vy, -gm*vx*vy                     , gama*et-gm/2.d0*rmod2-gm*vy*vy,    gama*vy /)

            !cccc----> terminos de estabilizacion y captura de choque
			arr = ar*tau

            !cccc  ----> armo el pressure switch
            choq(1)=alfa_mu*ar*ps(ipoin(1))*cte
            choq(2)=alfa_mu*ar*ps(ipoin(2))*cte
            choq(3)=alfa_mu*ar*ps(ipoin(3))*cte

            !cccc  ----> multiplico por partes para simplificar el asunto
            !cccc  ----> 'a' por las derivadas
        	aux(1)=                  ux(2)                                       +                    uy(3)
			aux(2:4) = a1(:,1)*ux(1) + a1(:,2)*ux(2) + a1(:,3)*ux(3) + a1(:,4)*ux(4) + a2(:,1)*uy(1) + a2(:,2)*uy(2)&
				+ a2(:,3)*uy(3) + a2(:,4)*uy(4)

			aux_phi = aux + phi_loc

            !cccc  ----> lo anterior por 'a' transpuesta
			aa(1) = aux_phi(2)
			aa(2:4) = a1(:,1)*aux_phi(1) + a1(:,2)*aux_phi(2) + a1(:,3)*aux_phi(3) + a1(:,4)*aux_phi(4)
			aa(5) = aux_phi(3)
			aa(6:8) = a2(:,1)*aux_phi(1) + a2(:,2)*aux_phi(2) + a2(:,3)*aux_phi(3) + a2(:,4)*aux_phi(4)

            !cccc  ----> ensamble del right hand side del nodo ipoin(1)
			brhs(:,1,ielem) = brhs(:,1,ielem) + (nx(1)*aa(1:4) + ny(1)*aa(5:8))*arr(:) + &
					aux(:)*rn1*ar + (nx(1)*ux(:) + ny(1)*uy(:))*choq(1)
			brhs(:,2,ielem) = brhs(:,2,ielem) + (nx(2)*aa(1:4) + ny(2)*aa(5:8))*arr(:) + &
					aux(:)*rn2*ar + (nx(2)*ux(:) + ny(2)*uy(:))*choq(2)
			brhs(:,3,ielem) = brhs(:,3,ielem) + (nx(3)*aa(1:4) + ny(3)*aa(5:8))*arr(:) + &
					aux(:)*rn3*ar + (nx(3)*ux(:) + ny(3)*uy(:))*choq(3)
        enddo
    enddo
	!$omp end do
	!$omp master
		do i=1,nelem
			rhs(:,n(1,i)) = brhs(:,1,i)
			rhs(:,n(2,i)) = brhs(:,2,i)
			rhs(:,n(3,i)) = brhs(:,3,i)
		end do
	!$omp end master
	!$omp end parallel
    return
end subroutine todo_mp_noAtomic

subroutine fill_3xnelem(x, nelem)
	real(8) x(3,nelem)

	do ielem=1,nelem
		x(1,ielem)=1
		x(2,ielem)=2
		x(3,ielem)=3
	end do
end subroutine fill_3xnelem

subroutine fill_inpoel(x, nelem, npoin)
	integer c
	integer x(3,nelem)
	c=0

	do ielem=1, nelem
		x(1,ielem)= mod(c,npoin)+1.d0
		x(2,ielem)= mod(c+1,npoin)+1.d0
		x(3,ielem)= mod(c+2,npoin)+1.d0
		c = c+3
	end do
end subroutine fill_inpoel

subroutine fill_u(x, npoin)
	real(8) x(4,npoin)

	do ipoin=1,npoin
		x(1,ipoin)= 1.d0
		x(2,ipoin)= 2.d0
		x(3,ipoin)= 3.d0
		x(4,ipoin)= 4.d0
	end do
end subroutine fill_u

subroutine fill_npoin(x, npoin)
	real(8) x(npoin)
	do ipoin=1,npoin
		x(ipoin)=mod(ipoin,100)+1.d0
	end do
end subroutine fill_npoin

subroutine fill_nelem(x, nelem)
	real(8) x(nelem)
	do ielem=1,nelem
		x(ielem)=(mod(ielem,10)+1.d0)/10.d0
	end do
end subroutine fill_nelem

subroutine checkit(rhs,npoin,c)
	real(8) rhs(4,npoin)
	c = sum(sum(rhs/10e6, DIM=1))
end subroutine checkit
