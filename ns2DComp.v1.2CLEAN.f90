MODULE DATOS_REFINAMIENTO
  REAL(8)ETA_REFIN,HHMAX_REFIN,HHMIN_REFIN
END MODULE DATOS_REFINAMIENTO

MODULE DATOS_ENTRADA
  REAL(8) FR,GAMA,UINF,VINF,TINF,RHOINF,MACHINF,PINF,CINF
END MODULE DATOS_ENTRADA

PROGRAM NSComp2D
  USE DATOS_REFINAMIENTO
  USE DATOS_ENTRADA
  IMPLICIT REAL(8) (A-H,O-Z)
  
  PARAMETER (NELE=250000,NODOS=130000,IFF=70000)
  
  INTEGER N(3,NELE),IVN(2,IFF)
  INTEGER IFIXRHO_NODE(IFF),IFIXV_NODE(IFF),IFIXT_NODE(IFF)
  INTEGER IELEM_SETS(10),ISET(3,10,IFF),INORMV_NODE(IFF)    
  INTEGER IPER_MASTER(IFF),IPER_SLAVE(IFF),IPER_AUX(IFF)
  !Nodos vecinos y laplaciano
  INTEGER IND(NODOS),INDEL(10,NODOS),IAUX(50),NN1(NODOS*15),NN2(NODOS*15) 
  INTEGER ILAUX(NODOS)
  !Bloque de moviemiento fijo y movil
  INTEGER IFM(IFF),I_M(IFF)
  !Suavizado de los elementos
  INTEGER SMOOTH_SIM(2,NODOS)

  REAL(8) RFIXRHO_VALUE(IFF),RFIXV_VALUEX(IFF),RFIXV_VALUEY(IFF)
  REAL(8) RNORMV_VALUEX(IFF),RNORMV_VALUEY(IFF)
  REAL(8) RFIXT_VALUE(IFF)
  
  REAL(8) X(NODOS),Y(NODOS),HHX(NELE),HHY(NELE)
  REAL(8) DNX(3,NELE),DNY(3,NELE),AREA(NELE),HH(NELE)
  
  REAL(8) VEL_X(NODOS),VEL_Y(NODOS),T(NODOS),P(NODOS),E(NODOS),RHO(NODOS)
  REAL(8) U(4,NODOS),U1(4,NODOS),UN(4,NODOS),RHS(4,NODOS),SHOC(NELE)
  REAL(8) M(NODOS),T_SUGN1(NELE),T_SUGN2(NELE),T_SUGN3(NELE)
  
  REAL(8) FX(10),FY(10),RM(10),F_VX(10),F_VY(10),GAMM(NODOS),RMACH(NODOS)
  REAL(8) ER(4),ERR(4),PS(NODOS)
  !Local time step
  REAL(8) DTL(NELE),DT(NELE)
  !Nodos vecinos y laplaciano
  REAL(8)RAUX(NODOS),ADIAG(NODOS),S(NODOS*15),B(NODOS),RES(NODOS)
  REAL(8) XPOS(NODOS),YPOS(NODOS),PK(NODOS),APK(NODOS),Z(NODOS),POS_AUX(NODOS)
  !Refinamiento
  REAL(8) HH_NEW(NODOS)
  !Velocidad de la malla
  REAL(8) W_X(NODOS),W_Y(NODOS)
  !Suavizado smoothing
  LOGICAL SMOOTH_FIX(NODOS)

  REAL(4) ZZ(2),ETIME 
  CHARACTER FILE*80
  
  !Tiempo de cpu
  HITE=ETIME(ZZ)
  OPEN(22,FILE='FORCES_HISTORY',STATUS='UNKNOWN') 
  OPEN (1,FILE='EULER.DAT',STATUS='OLD')
  READ (1,'(A)') FILE
  CLOSE(1)
  ILONG=LONG_FILE(FILE)
  !CCCC-------------------------------------------------------!CCCC
  !CCCC----> LECTURA DE LOS DATOS GENERALES DEL PROBLEMA <----!CCCC
  !CCCC-------------------------------------------------------!CCCC
  NMASTER=0 ; NSLAVE=0 ; ITLOCAL=0
  
  !Lee datos generales del problema
  OPEN(1,FILE=FILE(1:ILONG)//'-1.dat',STATUS='OLD')
  READ(1,*)
  READ(1,*) IRESTART,MAXITER,IPRINT,MOVIE,ITLOCAL
  READ(1,*)
  READ(1,*) FSAFE,UINF,VINF,MACHINF,TINF,RHOINF,PINF
  READ(1,*)
  READ(1,*) FMU,FGX,FGY,QH
  READ(1,*)
  READ(1,*) FK,FR,FCv,GAMA,NGAS
  !Cte de shock-capturing
  READ(1,*)
  READ(1,*) CTE
  !Refinamiento
  READ(1,*)
  READ(1,*)
  READ(1,*)
  READ(1,*)
  READ(1,*)ETA_REFIN,HHMAX_REFIN,HHMIN_REFIN
  CLOSE(1)
  CTE=1.D0/CTE

  !Multiples formas de entrada de datos
  CALL VARIABLES 

  OPEN (1,FILE=FILE(1:ILONG)//'.dat',STATUS='OLD') 
  !Coordenadas de los nodos 
  READ(1,*)      !Numero de nodos   !Numero de elementos  
  READ(1,*)       NNOD              ,NELEM 
  READ(1,*) 
  READ(1,*) NFIXRHO,NFIXVI,NFIXV,NELNORM,NFIXT,NSETS,NMASTER,NSLAVE,NFIX_MOVE,NMOVE
  
  IF (NNOD.GT.NODOS) STOP 'HACER MAS GRANDE "NODOS"' 
  IF (NELEM.GT.NELE) STOP 'HACER MAS GRANDE "NELE "' 
  IF (NFIXRHO.GT.IFF.OR.NFIXT.GT.IFF.OR.NELNORM.GT.IFF)  &
       STOP 'HACER MAS GRANDE "IFF"'   
  
  DO I=1,4 
     READ(1,*) 
  END DO
  
  !Coordenadas de los nodos 
  WRITE(*,*) ' reading coordinates....' 
  IERROR=1 
  DO INOD=1,NNOD 
     READ(1,*,ERR=27) INNOD,X(INNOD),Y(INNOD) 
  END DO
  WRITE(*,*) ' ready ' 
  WRITE(*,'(A,I5//)') ' TOTAL NODOS LEIDOS:',INNOD 
  
  !Conectividades de los elementos 
  WRITE(*,*) ' reading elements....' 
  IERROR=2 
  READ(1,*) 
  DO IELEM=1,NELEM 
     READ(1,*,ERR=27) IEL,N(1,IEL),N(2,IEL),N(3,IEL) 
  END DO
  WRITE(*,*) ' ready ' 
  WRITE(*,'(A,I6//)') ' TOTAL ELEMENTOS LEIDOS:',IEL 
  
  !Puntos con densidad prescrita
  WRITE(*,*) 'reading fix density points....'
  IERROR=3
  READ(1,*)
  DO IFIXRHO=1,NFIXRHO
     READ(1,*,ERR=27) IFIXRHO_NODE(IFIXRHO),RFIXRHO_VALUE(IFIXRHO)
     
     IF (RFIXRHO_VALUE(IFIXRHO).LT.0) THEN
        RFIXRHO_VALUE(IFIXRHO)=1.225d0
     ELSE
        RFIXRHO_VALUE(IFIXRHO)=RFIXRHO_VALUE(IFIXRHO)*RHOINF
     END IF
     
  END DO
  WRITE(*,*) ' ready '
  WRITE(*,'(A,I5//)') ' TOTAL NODOS CON DENSIDAD IMPUESTA:',NFIXRHO
  
  !Nodos con velocidad fija 
  WRITE(*,*) ' reading fix velocities (INFLOW)....' 
  IERROR=4 
  READ(1,*) 
  DO IFIXV=1,NFIXVI
     READ(1,*,ERR=27) IFIXV_NODE(IFIXV),RFIXV_VALUEX(IFIXV) &
          ,RFIXV_VALUEY(IFIXV)
     RFIXV_VALUEX(IFIXV)=RFIXV_VALUEX(IFIXV)*UINF
     RFIXV_VALUEY(IFIXV)=RFIXV_VALUEY(IFIXV)*VINF
  END DO
  WRITE(*,*) ' ready ' 
  WRITE(*,'(A,I5//)') ' TOTAL ELEMENTOS CON VELOCIDAD IMPUESTA:',NFIXVI 
  
  !Nodos con velocidad fija (no slip)
  !Calcula la temperatura de estancamiento para imponerla en las 
  !paredes donde se prescribe la condicion de "no slip"
  TWALL=TINF*(1.D0+(GAMA-1)/2.D0*MACHINF*MACHINF)
  IFIXT_A=0
  
  WRITE(*,*) ' reading fix velocities (NO SLIP)....' 
  IERROR=4 
  READ(1,*) 
  DO IFIXV=NFIXVI+1,NFIXVI+NFIXV
     READ(1,*,ERR=27) IFIXV_NODE(IFIXV),RFIXV_VALUEX(IFIXV) &
          ,RFIXV_VALUEY(IFIXV)
     RFIXV_VALUEX(IFIXV)=0.D0
     RFIXV_VALUEY(IFIXV)=0.D0
     IFIXT_A=IFIXT_A+1
     IFIXT_NODE(IFIXT_A)=IFIXV_NODE(IFIXV)
     RFIXT_VALUE(IFIXT_A)=TWALL
  END DO
  WRITE(*,*) ' ready ' 
  NFIXV=NFIXV+NFIXVI
  WRITE(*,'(A,I5//)') ' TOTAL ELEMENTOS CON VELOCIDAD IMPUESTA:',NFIXV
  
  !Velocidad normal nula 
  WRITE(*,*) ' reading elements with normal velocity prescribe....' 
  IERROR=5 
  READ(1,*) 
  DO INEL=1,NELNORM 
     READ(1,*,ERR=27) IVN(1,INEL),IVN(2,INEL)
  END DO
  WRITE(*,*) ' ready ' 
  WRITE(*,'(A,I6/)')' TOTAL ELEMENTOS CON VELOCIDAD NORMAL IMPUESTA:',NELNORM 
  
  !Temperatura prescrita
  WRITE(*,*) ' reading fix temperature nodes....'
  IERROR=8
  READ(1,*)
  DO IFIXT=1,NFIXT 
     READ(1,*,ERR=27) IFIXT_NODE(IFIXT+IFIXT_A),RFIXT_VALUE(IFIXT+IFIXT_A)
     RFIXT_VALUE(IFIXT+IFIXT_A)=RFIXT_VALUE(IFIXT+IFIXT_A)*TINF
  END DO
  WRITE(*,*) ' ready ' 
  NFIXT=NFIXT+IFIXT_A
  WRITE(*,'(A,I6/)') ' TOTAL NODOS COM TEMPERATURA IMPUESTA:',NFIXT  
  
  !Lectura de los sets para el calculo de las fuerzas 
  WRITE(*,*) ' reading sets nodes....'
  IERROR=9
  DO ISETS=1,10
     IELEM_SETS(ISETS)=0
  END DO
  NSET_NUMB=0
  READ(1,*)
  DO ISETS=1,NSETS
     READ(1,*,ERR=27) NELE_SET,ISET1,ISET2,ISETNUMB
     IF (ISETNUMB.GT.NSET_NUMB) NSET_NUMB=ISETNUMB
     IELEM_SETS(ISETNUMB)=IELEM_SETS(ISETNUMB)+1
     ISET(1,ISETNUMB,IELEM_SETS(ISETNUMB))=ISET1
     ISET(2,ISETNUMB,IELEM_SETS(ISETNUMB))=ISET2
     ISET(3,ISETNUMB,IELEM_SETS(ISETNUMB))=NELE_SET
  END DO
  WRITE(*,*) ' ready ' 
  WRITE(*,'(A,I6)') ' NUMERO DE SETS:',NSET_NUMB
  WRITE(*,'(A,I6//)') ' NUMERO DE ELEMENTOS DE LOS SETS:',NSETS
  
  !Periodical master 
  WRITE(*,*) ' reading nodes with periodical master and slave condition....' 
  IERROR=10
  IF(NMASTER.NE.NSLAVE)THEN
     WRITE(*,*)'ERROR NODOS MASTER DISTINTO NODOS SLAVE'
     STOP
  END IF
  READ(1,*) 
  DO IMASTER=1,NMASTER
     READ(1,*,ERR=27)IPER_MASTER(IMASTER)
  END DO
  READ(1,*) 
  DO IMASTER=1,NMASTER
     READ(1,*,ERR=27)IPER_SLAVE(IMASTER)
  END DO
  WRITE(*,*) ' ready ' 
  WRITE(*,'(A,I6/)')' TOTAL NODOS CON CONDICION PERIODICA:',NMASTER+NSLAVE 
  
  !Nodos con movimiento fijo 
  WRITE(*,*) ' reading fix movement ....' 
  READ(1,*) 
  DO IM=1,NFIX_MOVE
     READ(1,*,ERR=27) IFM(IM),walberto
  END DO
  WRITE(*,*) ' ready ' 
  WRITE(*,'(A,I5//)') ' TOTAL NODOS CON MOVIMIENTO FIJO:',NFIX_MOVE  

  !Nodos con movimiento 
  WRITE(*,*) ' reading movement ....' 
  READ(1,*) 
  DO IM=1,NMOVE
     READ(1,*,ERR=27) I_M(IM),walberto
  END DO
  WRITE(*,*) ' ready ' 
  WRITE(*,'(A,I5//)') ' TOTAL NODOS CON MOVIMIENTO FIJO:',NMOVE  

  IERROR=0 
  !Control de errores en la lectura 
27 IF (IERROR.NE.0) THEN 
     CALL READ_ERROR(IERROR,IEL) 
  END IF
 
  !CCCCCCCCCCCCCCCCC                                  CCCCCCCCCCCCCCCCC 
  !CCCCCCCCCCCCCCCCC       COMIENZO DEL CALCULO       CCCCCCCCCCCCCCCCC
  !CCCCCCCCCCCCCCCCC                                  CCCCCCCCCCCCCCCCC

  !Asigna valor a las variables primitivas y realiza un restart si IRESTART == 1
  GAMM=GAMA
  CALL RESTART(NNOD,FILE,IRESTART &
       ,GAMM &
       ,U,VEL_X,VEL_Y,T)
  IF(NGAS.NE.1) GAMM=GAMA
 
  !Calculo de los nodos con periodicidad
  IF(NMASTER.NE.0.AND.NSLAVE.NE.0)THEN
     CALL PERIODIC(NNOD,IFF,NMASTER,NSLAVE &
          ,IPER_MASTER,IPER_SLAVE,X &
          ,IPER_AUX)
  END IF

  !Calculo de las normales 
  CALL NORMALES(VEL_X,VEL_Y,X,Y,IVN,NELNORM,NNOD,IFF &
       ,NNORMV,INORMV_NODE,RNORMV_VALUEX,RNORMV_VALUEY)
   
  !Calculo de las derivadas, el area y la longitud carac.   
  CALL DERIV(HMIN,HHX,HHY,X,Y,DNX,DNY,AREA,HH,N,NELEM,NNOD)

  !Calculo de las matrices de masas lumped
  CALL MASAS(M,AREA,N,NELEM,NNOD)

  !Calculo de los nodos vecinos y el laplaciano
  CALL LAPLACE(NNOD,NELEM,N,AREA,DNX,DNY,S,NN1,NN2,IAUX &
       ,NPOS,IND,INDEL,ADIAG,RAUX)
 
!!$ ! Fijo que nodos se mueven y cuales no..
  NNMOVE=NFIX_MOVE+NMOVE
  KK=0
  DO I=1,NMOVE
     KK=KK+1
     N1=I_M(I)
     ILAUX(KK)=N1
  END DO
  
  DO I=1,NFIX_MOVE
     KK=KK+1
     N1=IFM(I)
     ILAUX(KK)=N1
  END DO

!!$ ! Nodos que no son suavizados
     SMOOTH_FIX(1:NNOD)=.TRUE.
     SMOOTH_SIM(1:2,1:NNOD)=1
     DO I=1,NMOVE
        N1=I_M(I)
        SMOOTH_FIX(N1)=.FALSE.
     END DO
     
     DO I=1,NFIX_MOVE
        N1=IFM(I)
        SMOOTH_FIX(N1)=.FALSE.
     END DO
  
!!$!Movimiento en y
!!$  KK=0
!!$  YPOS=Y ; B=0.D0 ; POS_AUX=0.D0
!!$  DO I=1,NMOVE
!!$     KK=KK+1
!!$     POS_AUX(KK)=DELTAX
!!$  END DO
!!$  DO I=1,NFIX_MOVE
!!$     KK=KK+1
!!$     POS_AUX(KK)=0.D0
!!$  END DO
  
!!$  CALL GRADCONJ2(S,YPOS,B,NN1,NN2,NNOD,NPOS &
!!$       ,ILAUX,POS_AUX,NNMOVE,RAUX,ADIAG &
!!$       ,PK,APK,Z) 
   
  !CCCC----> ABRE EL ARCHIVO DE CONVERGENCIA
  !CCCC----------------------------------------------
  OPEN(7,FILE=FILE(1:ILONG)//'.cnv',STATUS='UNKNOWN')

  !CCCC -----------------------------------------CCCC
  !CCCC ----> COMIENZO DE LAS ITERACIONES <----  CCCC
  !CCCC -----------------------------------------CCCC
  TIME=0.D0
  ISAL=0
  ITER=0
  DTMIN=0.D0
  ITERPRINT=0
  FLUX1=DABS(FGX+FGY+QH) !Variable para calcular los terminos fuentes
  RMACH=(UINF**2+VINF**2)/DSQRT(GAMA*FR*TINF)
  !Orden de integracion de RUNGE-KUTTA    
  NRK=4

  !Defino velocidad de la malla
  W_X=0.D0 ; W_Y=0.D0
  MOVING=0 ; YPOS=0.D0 ; IMOOTH=0

  WRITE(*,'(A,I3,A)') '****-------> RUNGE-KUTTA DE',NRK,'  ORDEN <-------****'
  WRITE(*,*)''
  
  !Abre el archivo donde se imprimen los resultados
  IF (MOVIE.EQ.1)THEN
     OPEN(2,FILE=FILE(1:ILONG)//'.flavia.res',STATUS='UNKNOWN')
  END IF
   
  DO WHILE (ITER.LT.MAXITER.AND.ISAL.EQ.0)
     
     ITER=ITER+1
     !CCCC  ----> CALCULOS DE LOS TERMINOS PARA CADA ITERACION
     IF(MOVING.EQ.1)THEN
        !CCCC-----> CALCULO DE LAS NORMALES 
        !CCCC---------------------------------------------------------
        CALL NORMALES(VEL_X,VEL_Y,X,Y,IVN,NELNORM,NNOD,IFF &
             ,NNORMV,INORMV_NODE,RNORMV_VALUEX,RNORMV_VALUEY)
        !CCCC----> CALCULO DE LAS DERIVADAS, EL AREA Y LA LONGITUD CARAC.   
        !CCCC------------------------------------------------------------
        CALL DERIV(HMIN,HHX,HHY,X,Y,DNX,DNY,AREA,HH,N,NELEM,NNOD)
        
        !CCCC----> CALCULO DE LAS MATRICES DE MASAS LUMPED
        !CCCC  ----------------------------------------------------------
        CALL MASAS(M,AREA,N,NELEM,NNOD)
        
        !CCCC----> CALCULO DE LOS NODOS VECINOS Y EL LAPLACIANO
        !CCCC  ------------------------------------------------------
        CALL LAPLACE(NNOD,NELEM,N,AREA,DNX,DNY,S,NN1,NN2,IAUX &
             ,NPOS,IND,INDEL,ADIAG,RAUX)
     END IF

     !CCCC  ----> CALCULO DEL dT
     !CCCC  --------------------
     
     CALL DELTAT(NNOD,NELEM,N,FSAFE,FR,GAMA &
          ,T,U,W_X,W_Y,DTMIN,HH,DT)

     !CCCC----> FUNCION PARA ESTABLECER LOCAL-TIME-STEP
     IF (ITLOCAL.NE.0) THEN
        DTFACT = 1.D0 - DEXP(-ITER*4.6D0/ITLOCAL)
        DO IELEM=1,NELEM
           DTL(IELEM)= DTMIN*DTFACT+DT(IELEM)*(1.D0-DTFACT)
        END DO
     ELSE
        DO IELEM=1,NELEM
           DTL(IELEM)=DTMIN
        END DO
     END IF
     
     TIME=TIME+DTMIN
     
     U1=U
        
     DO IRK=1,NRK
        
        RK_FACT=1.D0/(NRK+1-IRK)
        
        !CCCC  ----> SOLO CALCULO UNA VEZ EL TERMINO DE ESTABILIZACION
        
        IF(IRK.EQ.1)THEN
           CALL CUARTO_ORDEN(NNOD,NELEM,N,DNX,DNY,AREA,M &
                ,U1,UN,GAMM,FR)
          
           CALL ESTAB(NNOD,NELEM,N,DNX,DNY,U,T &
                ,VEL_X,VEL_Y,W_X,W_Y,GAMA,FR,FMU &
                ,DTMIN,RHOINF,TINF,UINF,VINF &
                ,T_SUGN1,T_SUGN2,T_SUGN3,SHOC,GAMM)
        END IF
                
        !CCCC ------> CALCULO EL PRESSURE-SWICHT
        !CCCC --------------------------------------------------------
        CALL NEWPRES(NNOD,NPOS,NN1,NN2,DABS(RMACH),PS)
        RHS=0.D0

        CALL TODO(NNOD,NELEM,N,DNX,DNY,AREA,HHX,HHY &
             ,T_SUGN1,T_SUGN2,T_SUGN3,SHOC,PS,DTL &
             ,U1,UN,RHS,P,GAMM,FR,FMU,FK,FCV,TINF,CTE)
        
        !CCCC  ----> CALCULO DE LOS TERMINOS FUENTES

        CALL FUENTE(NNOD,NELEM,N,AREA,DNX,DNY,W_X,W_Y,DTL &
             ,U,RHS,FGX,FGY,QH)
        
        !CCCC ----> INTEGRADOR TEMPORAL
        DO INOD=1,NNOD
           DO J=1,4
              U1(J,INOD)=U(J,INOD)-RK_FACT/M(INOD)*RHS(J,INOD)
           END DO
        END DO

        !CCCC----------------------------------------------
        !CCCC----> PASA A LA VARIABLE PRIMARIA PARA APLICAR 
        !CCCC----> LAS CONDICIONES DE CONTORNO 
        !CCCC---------------------------------------------- 
        IF(NGAS.EQ.1) GO TO 111
        DO INOD=1,NNOD
           RHO(INOD)=U1(1,INOD)
           VEL_X(INOD)=U1(2,INOD)/RHO(INOD)
           VEL_Y(INOD)=U1(3,INOD)/RHO(INOD)
           E(INOD)=U1(4,INOD)/RHO(INOD)
           VEL2=(VEL_X(INOD)**2.D0+VEL_Y(INOD)**2.D0)
           P(INOD)=RHO(INOD)*(GAMM(INOD)-1.D0)*(E(INOD)-.5D0*VEL2)
           T(INOD)=P(INOD)/(RHO(INOD)*FR)
           RMACH(INOD)=DSQRT(VEL2)/T(INOD)
        END DO
111     CONTINUE
        
        !CCCC----> CASO PARA AIRE EN EQUILIBRIO
        !CCCC----------------------------------------------
        IF(NGAS.NE.0)THEN
           DO INOD=1,NNOD
              RHO(INOD)=U1(1,INOD)
              E(INOD)=U1(4,INOD)/RHO(INOD)
              VEL_X(INOD)=U1(2,INOD)/RHO(INOD)
              VEL_Y(INOD)=U1(3,INOD)/RHO(INOD)
              VEL2=VEL_X(INOD)**2.D0+VEL_Y(INOD)**2.D0
              P(INOD)=RHO(INOD)*(GAMM(INOD)-1.D0)*(E(INOD)-.5D0*VEL2)
              T(INOD)=P(INOD)/(RHO(INOD)*FR)
              IF(IRK.EQ.NRK)THEN
                 CALL TGAS(E(INOD)-.5D0*VEL2,RHO(INOD),PGAS,AGAS,TGASi,GAMI)
                 GAMM(INOD)=GAMI
                 P(INOD)=PGAS
                 T(INOD)=TGASi
                 RMACH(INOD)=DSQRT(VEL2)/AGAS
              END IF
           END DO
        END IF
        

        !CCCC---------------------------------------CCCC 
        !CCCC  ----> CONDICIONES DE CONTORNO <----  CCCC 
        !CCCC---------------------------------------CCCC 
        
        !CCCC----> VELOCIDADES IMPUESTAS 
        !CCCC--------------------------- 
        CALL FIXVEL(NNOD,VEL_X,VEL_Y &
             ,NFIXV,IFIXV_NODE,RFIXV_VALUEX,RFIXV_VALUEY)

        !CCCC----> CORRECCION DE LAS VELOCIDADES NORMALES 
        !CCCC-------------------------------------------- 
        CALL NORMALVEL(NNOD,VEL_X,VEL_Y,W_X,W_Y,RHO,INORMV_NODE &
             ,RNORMV_VALUEX,RNORMV_VALUEY,NNORMV) 

        !CCCC----> VALORES IMPUESTOS 
        !CCCC----------------------- 
        CALL FIX(NNOD,RHO,VEL_X,VEL_Y,T,E &
             ,NFIXRHO,IFIXRHO_NODE,RFIXRHO_VALUE &
             ,NFIXT,IFIXT_NODE,RFIXT_VALUE &
             ,FR,GAMM,GAMA)
        
        DO INOD=1,NNOD
           U1(1,INOD)=RHO(INOD)
           U1(2,INOD)=VEL_X(INOD)*RHO(INOD)
           U1(3,INOD)=VEL_Y(INOD)*RHO(INOD)
           U1(4,INOD)=E(INOD)*RHO(INOD)
        END DO 

     END DO
     
     !CCCC-------------------------------------------------CCCC 
     !CCCC  ----> CALCULO DE LOS ERRORES (RESIDUOS) <----  CCCC 
     !CCCC  ----> PARA CONTROLAR LA CONVERGENCIA    <----  CCCC   
     !CCCC-------------------------------------------------CCCC          
     ER=0.D0; ERR=0.D0
     DO INOD=1,NNOD
        DO J=1,4
           ER(J)=ER(J)+(U(J,INOD)-U1(J,INOD))**2.D0
           ERR(J)=ERR(J)+U1(J,INOD)**2.D0
           U(J,INOD)=U1(J,INOD)
        END DO
     END DO
     
     
     WRITE(7,'(I7,4E14.6)') ITER,TIME,DSQRT(ER(1)/ERR(1)),DSQRT(ER(4)/ERR(4))
     
     
     !CCCC--------------------------CCCC 
     !CCCC  ----> IMPRESION <----   CCCC 
     !CCCC--------------------------CCCC 
     ITERPRINT=ITERPRINT+1
     
     IF (ITERPRINT.EQ.IPRINT.OR.ITER.EQ.MAXITER) THEN
        WRITE(*,'(A,I5,3X,A,2X,E15.4)') 'Impresion en el paso:' ,ITER,'ERROR:',MAXVAL(DSQRT(ER/ERR))
        
        CALL PRINTFLAVIA(FR,GAMM,RHO,VEL_X-W_X,VEL_Y-W_Y,P,T,E,RMACH,XPOS,YPOS &
             ,NNOD,ITER,MOVIE,FILE,ILONG)
        ITERPRINT=0
        !CALL PRINTFLAVIA(FR,GAMM,U,XPOS,YPOS,RMACH,NNOD,ITER,MOVIE,FILE,ILONG)
        !CALL FORCES(NSETS,NSET_NUMB,IELEM_SETS,ISET,NNOD,IFF &
        !     ,X,Y,P &
        !     ,FX,FY,RM)
        
        IF(FMU.NE.0.D0)THEN
           CALL FORCE_VISC(NELEM,NNOD,IFF &
                ,IELEM_SETS,ISET,N,NSET_NUMB,UINF,VINF,RHOINF,TINF &
                ,X,Y,P,T,VEL_X,VEL_Y,DNX,DNY,FMU,RHO &
                ,F_VX,F_VY)      
        END IF
        
        OPEN(33,FILE='FORCES',STATUS='UNKNOWN')           
        
        DO ISET_NUMB=1,NSET_NUMB
           WRITE(33,'(A,I2)') 'SET NUMERO',ISET_NUMB
           WRITE(33,'(A,E14.5)') 'FUERZA EN X:',FX(ISET_NUMB)
           WRITE(33,'(A,E14.5/)') 'FUERZA EN Y:',FY(ISET_NUMB)
           ! WRITE(33,'(A,E14.5/)') 'MOMENTO:',RM(ISET_NUMB)
           WRITE(33,'(A,E14.5)') 'FUERZA VISCOSA EN X:',F_VX(ISET_NUMB)
           WRITE(33,'(A,E14.5/)') 'FUERZA VISCOSA EN Y:',F_VY(ISET_NUMB)
           WRITE(33,'(A,E14.5)') 'FUERZA TOTAL EN X:',FX(ISET_NUMB)+F_VX(ISET_NUMB)
           WRITE(33,'(A,E14.5/)') 'FUERZA TOTAL EN Y:',FY(ISET_NUMB)+F_VY(ISET_NUMB)
        END DO
        CLOSE(33)
        
        CALL PRINTREST(ITER,NNOD,U,T,GAMM,TIME,FILE,ILONG)
        
     END IF    
     
  END DO

  !Cierra el archivo de resultados..
  IF (MOVIE.EQ.1)THEN
     CLOSE(2)
  END IF

  !CCCC-----------------------------------CCCC 
  !CCCC  ----->    FIN DEL LOOP   <-----  CCCC 
  !CCCC-----------------------------------CCCC 
  
  !Cierra el archivo de la convergencia  
  WRITE(*,*)'****---------------------****'
  WRITE(*,*)'****----> REFINANDO <----****'
  WRITE(*,*)'****---------------------****'

!!$  !CALL NEW_SIZE(NNOD,NELEM,N,DNX,DNY,AREA,M,HH,U &
!!$  !     ,HH_NEW)
  DO I=1,NNOD
     RMACH(I)=DSQRT(VEL_X(I)**2+VEL_Y(I)**2)/DSQRT(GAMA*FR*T(I))
  END DO
  hh_new=hhmax_refin
  
  CALL ESTIMADOR_ERR(NNOD,NELEM,N,DNX,DNY,HH,M,AREA,3.5D0,RMACH,P)
  DO I=1,NNOD
     IF(HH_NEW(I).GT.P(I))HH_NEW(I)=P(I)
  END DO

  CALL ESTIMADOR_ERR(NNOD,NELEM,N,DNX,DNY,HH,M,AREA,4.D0,RHO,P)
  DO I=1,NNOD
     IF(HH_NEW(I).GT.P(I))HH_NEW(I)=P(I)
  END DO

  HHMIN_REFIN=0.01D0
  CALL ESTIMADOR_ERR(NNOD,NELEM,N,DNX,DNY,HH,M,AREA,6.D0,RMACH,P)
  DO I=1,NNOD
     IF(HH_NEW(I).GT.P(I))HH_NEW(I)=P(I)
  END DO

  !Acomodo los tamanos sobre el/los cuerpos
  DO ISET_NUMB=1,NSET_NUMB
     DO II=1,IELEM_SETS(ISET_NUMB)
            
        N1=ISET(1,ISET_NUMB,II)
        N2=ISET(2,ISET_NUMB,II)
        IELEM=ISET(3,ISET_NUMB,II)
        HH_NEW(N1)=.01D0 !hhmin_refin!DSQRT(2.D0*AREA(IELEM))
        HH_NEW(N2)=.01D0 !hhmin_refin!DSQRT(2.D0*AREA(IELEM))
     END DO
  END DO
  
  !DO I=1,NELEM
  !   AR=AREA(I)
  !   CC=DSQRT(1.5*2.D0*AR)
  !   DO J=1,3
  !      N1=N(J,I)
  !      IF(HH_NEW(N1).GT.CC) HH_NEW(N1)=HHMIN_REFIN
  !   END DO
  !END DO

  OPEN(12,FILE='remeshing.msh',STATUS='UNKNOWN')
  
  WRITE(12,'(A)') 'BackgroundMesh V 1.0'
  WRITE(12,'(A)') 'MESH dimension 2 ElemType Triangle Nnode 3'
  WRITE(12,'(A)') 'Coordinates'
  DO INOD=1,NNOD
     WRITE(12,'(I7,3E14.6)') INOD,X(INOD),Y(INOD)
  END DO
  WRITE(12,'(A)') 'End Coordinates'
  
  WRITE(12,'(A)') 'Elements'
  DO IELEM=1,NELEM
     WRITE(12,'(5I7)') IELEM,N(1,IELEM),N(2,IELEM),N(3,IELEM)
  END DO
  WRITE(12,'(A)') 'End Elements'
  
  WRITE(12,'(A)') 'DesiredSize Nodes'
  DO I=1,NNOD
     WRITE(12,'(I7,E14.6)') I,HH_NEW(I)
  END DO
  WRITE(12,'(A)') 'End DesiredSize Nodes'          
  
  CLOSE(12)

  !Tiempo de cpu
  HITE=ETIME(ZZ)
  CLOSE(7)
  WRITE(*,*)
  WRITE(*,'(3X,A)')'TIEMPO DE CPU'
  WRITE(*,'(X,F10.4,2X,A)')HITE/60.D0,'Minutos'
  WRITE(*,*)''
  WRITE(*,*)'****-------> FIN DEL CALCULO <-------****'
  WRITE(*,*)

END PROGRAM NSComp2D

!Restart
SUBROUTINE RESTART(NNOD,FILE,IRESTART &
     ,GAMM &
     ,U,VEL_X,VEL_Y,T)

  USE DATOS_ENTRADA
  IMPLICIT REAL(8) (A-H,O-Z)
  
  REAL(8) U(4,NNOD),VEL_X(NNOD),VEL_Y(NNOD),T(NNOD),GAMM(NNOD)
  
  CHARACTER FILE*80
  
  ILONG=LONG_FILE(FILE)
  
  IF (IRESTART.NE.1) THEN   
     !CCCC----> SI IOLDSOL .NE. 1 INICIALIZA LAS VARIABLES
     !CCCC----> CON LOS VALORES DEL INFINITO
     !CCCC------------------------------------------------
     DO INOD=1,NNOD
        
        !CCCC---> VARIABLES CONSERVATIVAS
        U(1,INOD)=RHOINF
        U(2,INOD)=RHOINF*UINF
        U(3,INOD)=RHOINF*VINF
        ENERGIA=PINF/((GAMM(INOD)-1.D0)*RHOINF)+.5*(UINF**2.D0+VINF**2.D0)
        U(4,INOD)=ENERGIA*RHOINF
        VEL_X(INOD)=UINF
        VEL_Y(INOD)=VINF
        T(INOD)=TINF
     END DO
  ELSE                     
     !CCCC----> SI IOLDSOL .EQ. 1 HACE UN RESTART
     !CCCC---------------------------------------
     OPEN(1,FILE=FILE(1:ILONG)//'.RST',FORM='UNFORMATTED',STATUS='UNKNOWN')
     
     READ(1) ITER_OLD,TIME
     
     DO INOD=1,NNOD
        READ(1) (U(J,INOD),J=1,4),T(INOD),GAMM(INOD)
     END DO
     
     CLOSE(1)
     
     
  END IF
  
  RETURN
END SUBROUTINE RESTART

SUBROUTINE PERIODIC(NNOD,IFF,NMASTER,NSLAVE &
     ,IPER_MASTER,IPER_SLAVE,X &
     ,IPER_AUX)
  
  IMPLICIT REAL(8) (A-H,O-Z)
  
  INTEGER IPER_MASTER(IFF),IPER_SLAVE(IFF),IPER_AUX(IFF)
  REAL(8) X(NNOD)

  DO I=1,NMASTER
     N1=IPER_MASTER(I)
     DO J=1,NMASTER
        N2=IPER_SLAVE(J)
        XX=DABS(X(N1)-X(N2))
        IF(XX.LT.1.D-6)THEN
           IPER_AUX(I)=N2
        END IF
     END DO
  END DO
  RETURN
END SUBROUTINE PERIODIC

!Calcular normales en los nodos
SUBROUTINE NORMALES(RNX,RNY,X,Y,IVN,NELNORM,NNOD,IFF &
     ,NNORMV,INORMV_NODE,RNORMV_VALUEX,RNORMV_VALUEY)
  
  IMPLICIT REAL(8) (A-H,O-Z)
  
  INTEGER IVN(2,IFF),INORMV_NODE(IFF)
  
  REAL(8) X(NNOD),Y(NNOD)
  REAL(8) RNX(NNOD),RNY(NNOD),BAUX(NNOD)
  REAL(8) RNORMV_VALUEX(IFF),RNORMV_VALUEY(IFF)
  REAL(8) NX,NY
  
  DO INOD=1,NNOD
     RNX(INOD)=0.D0
     RNY(INOD)=0.D0
     BAUX(INOD)=0.D0
  END DO
  
  !CCCCC NORMALES DE LAS PAREDES   !CCCCC
  DO I=1,NELNORM
     
     N1=IVN(1,I)
     N2=IVN(2,I)
     
     NX=(Y(N2)-Y(N1))
     NY=-(X(N2)-X(N1))
     
     RMOD=DSQRT(NX*NX+NY*NY)
     NX=NX/RMOD
     NY=NY/RMOD
     
     DO IN=1,2
        RNX(IVN(IN,I))=RNX(IVN(IN,I))+NX*RMOD
        RNY(IVN(IN,I))=RNY(IVN(IN,I))+NY*RMOD
        BAUX(IVN(IN,I))=BAUX(IVN(IN,I))+RMOD
     END DO
  END DO
  
  NNORMV=0
  DO INOD=1,NNOD
     IF (BAUX(INOD).GT.1.D-6) THEN
        
        RX=RNX(INOD)/BAUX(INOD)
        RY=RNY(INOD)/BAUX(INOD)
        RMOD=DSQRT(RX*RX+RY*RY)
        IF (RMOD.GT.0.2D0) THEN
           NNORMV=NNORMV+1
           INORMV_NODE(NNORMV)=INOD
           RNORMV_VALUEX(NNORMV)=RX/RMOD
           RNORMV_VALUEY(NNORMV)=RY/RMOD
        END IF
     END IF
  END DO
  
  RETURN
END SUBROUTINE NORMALES

!Calcular shape function derivatives & elements area
SUBROUTINE DERIV(HMIN,HHX,HHY,X,Y,DNX,DNY,AREA,HH,N,NELEM,NNOD)
  IMPLICIT REAL(8) (A-H,O-Z)
  
  INTEGER N(3,NELEM)
  
  REAL(8) X(NNOD),Y(NNOD),AREA(NELEM)
  REAL(8) DNX(3,NELEM),DNY(3,NELEM),HH(NELEM),HHX(NELEM),HHY(NELEM)
  
  HMIN=1.D10
  DO IELEM=1,NELEM
     N1=N(1,IELEM)
     N2=N(2,IELEM)
     N3=N(3,IELEM)
     
     AR=X(N2)*Y(N3)+X(N3)*Y(N1)+X(N1)*Y(N2)-(X(N2)*Y(N1)+X(N3)*Y(N2)+X(N1)*Y(N3))
     AR=AR/2.D0
     AREA(IELEM)=AR
     
     DNX(1,IELEM)=(Y(N2)-Y(N3))/(2.D0*AR)
     DNX(2,IELEM)=(Y(N3)-Y(N1))/(2.D0*AR)
     DNX(3,IELEM)=(Y(N1)-Y(N2))/(2.D0*AR)        
     DNY(1,IELEM)=(X(N3)-X(N2))/(2.D0*AR)
     DNY(2,IELEM)=(X(N1)-X(N3))/(2.D0*AR)
     DNY(3,IELEM)=(X(N2)-X(N1))/(2.D0*AR)
     
     HH(IELEM)=DSQRT(AREA(IELEM))
     IF (HH(IELEM).LT.HMIN) HMIN=HH(IELEM)
     ATA1=MIN( X(N3)-X(N2) , X(N1)-X(N3) , X(N2)-X(N1) )
     ATA2=MIN( Y(N3)-Y(N2) , Y(N1)-Y(N3) , Y(N2)-Y(N1) ) 
     HHX(IELEM)=ABS(ATA1)
     HHY(IELEM)=ABS(ATA2)
  END DO
  
  RETURN
END SUBROUTINE DERIV

!Calcular matriz condensada
SUBROUTINE MASAS(M,AREA,N,NELEM,NNOD)

  IMPLICIT REAL(8) (A-H,O-Z)
  INTEGER N(3,NELEM)
  REAL(8) M(NNOD),AREA(NELEM)
  
  M=0.D0 
  
  DO IELEM=1,NELEM
     
     N1=N(1,IELEM)
     N2=N(2,IELEM)
     N3=N(3,IELEM)
     
     AR=AREA(IELEM)/3.D0
     
     M(N1)=M(N1)+AR
     M(N2)=M(N2)+AR
     M(N3)=M(N3)+AR 
     
  END DO
  
  RETURN
END SUBROUTINE MASAS

!Calcular delta t
SUBROUTINE DELTAT(NNOD,NELEM,N,FSAFE,FR,GAMA &
     ,T,U,W_X,W_Y,DTMIN,HH,DT)
  
  IMPLICIT REAL(8) (A-H,O-Z)
  
  INTEGER N(3,NELEM)
  REAL(8) HH(NELEM),U(4,NNOD)
  REAL(8) T(NNOD),DT(NELEM)
  REAL(8) W_X(NNOD),W_Y(NNOD)
  DTMIN=1.D20
  
  DO IELEM=1,NELEM
     VUMAX=0.D0
     VVMAX=0.D0

     T_IEL=(T(N(1,IELEM))+T(N(2,IELEM))+T(N(3,IELEM)))/3.D0
     
     !GM=GAMA-1.D0 !(GAMM(N(1,IELEM))+GAMM(N(2,IELEM))+GAMM(N(3,IELEM)))/3.D0     
     VC=DSQRT(GAMA*FR*T_IEL)
     
     !CCCC ----> CALCULO DE LA MAXIMA VELOCIDAD ELEMENTAL
     DO I=1,3
        N1=N(I,IELEM)
        VU=DABS(U(2,N1)/U(1,N1)-W_X(N1))
        VV=DABS(U(3,N1)/U(1,N1)-W_Y(N1))
        
        IF (VU.GT.VUMAX) VUMAX=VU
        IF (VV.GT.VVMAX) VVMAX=VV            
     END DO
     
     VEL=(VUMAX**2.D0+VVMAX**2.D0)**.5D0         
     
     DTELEM=.5D0*HH(IELEM)/(VC+VEL)
     DT(IELEM)=DTELEM*FSAFE
     IF (DTELEM.LT.DTMIN) DTMIN=DTELEM*FSAFE
     
  END DO
  
  COTA=10.D0*DTMIN !EL VALOR 10 ESTA PUESTO A OJO #MODIFICAR SI ES NECESARIO#
  DO IELEM=1,NELEM
     IF(DT(IELEM).GT.COTA) DT(IELEM)=COTA
  END DO
  RETURN
END SUBROUTINE DELTAT

!Calcular la integral de los terminos de estabilizacion.
SUBROUTINE CUARTO_ORDEN(NNOD,NELEM,N,DNX,DNY,AREA,M &
     ,U,UN,GAMM,FR)
  
  IMPLICIT REAL(8) (A-H,O-Z)
  
  INTEGER N(3,NELEM)
  
  REAL(8) DNX(3,NELEM),DNY(3,NELEM),U(4,NNOD),UN(4,NNOD)
  REAL(8) AREA(NELEM),M(NNOD),GAMM(NNOD)
  REAL(8) ALF(3),BET(3)
  
  !C---->     RHO=U(1)
  !C---->     RHO*VEL_X=U(2)
  !C---->     RHO*VEL_Y=U(3)
  !C---->     RHO*ET=U(4)
  
  DATA ALF/.5D0,.5D0,0.D0/
  DATA BET/0.D0,.5D0,.5D0/
  
  NGAUSS=3    !PTOS DE GAUSS DONDE VOY A INTERGRAR
  
  UN=0.D0
  
  DO IELEM=1,NELEM
     
     N1=N(1,IELEM)
     N2=N(2,IELEM)
     N3=N(3,IELEM)
     
     GAMA=(GAMM(N1)+GAMM(N2)+GAMM(N3))/3.D0
     GM=GAMA-1.D0
     RNX1=DNX(1,IELEM) ; RNX2=DNX(2,IELEM) ; RNX3=DNX(3,IELEM)
     RNY1=DNY(1,IELEM) ; RNY2=DNY(2,IELEM) ; RNY3=DNY(3,IELEM)
     
     !CCCC  ----> DERIVADA DE RHO           
     DRX= RNX1*U(1,N1)+RNX2*U(1,N2)+RNX3*U(1,N3)
     DRY= RNY1*U(1,N1)+RNY2*U(1,N2)+RNY3*U(1,N3)
     !CCCC  ----> DERIVADA DE RHO.VEL_X           
     DRUX= RNX1*U(2,N1)+RNX2*U(2,N2)+RNX3*U(2,N3)
     DRUY= RNY1*U(2,N1)+RNY2*U(2,N2)+RNY3*U(2,N3)
     !CCCC  ----> DERIVADA DE RHO.VEL_Y           
     DRVX= RNX1*U(3,N1)+RNX2*U(3,N2)+RNX3*U(3,N3)
     DRVY= RNY1*U(3,N1)+RNY2*U(3,N2)+RNY3*U(3,N3)
     !CCCC  ----> DERIVADA DE RHO.ET           
     DREX= RNX1*U(4,N1)+RNX2*U(4,N2)+RNX3*U(4,N3)
     DREY= RNY1*U(4,N1)+RNY2*U(4,N2)+RNY3*U(4,N3)
    
     
     AR=AREA(IELEM)/3.D0
     DO J=1,NGAUSS
        
        RN1=1.D0-ALF(J)-BET(J)
        RN2=ALF(J)
        RN3=BET(J)
        
        !CCCC  ----> INTEGRO LAS VARIABLES PRIMITIVAS EN LOS PUNTOS DE GAUSS
        U1= RN1*U(1,N1)+RN2*U(1,N2)+RN3*U(1,N3)
        U2= RN1*U(2,N1)+RN2*U(2,N2)+RN3*U(2,N3)
        U3= RN1*U(3,N1)+RN2*U(3,N2)+RN3*U(3,N3)
        U4= RN1*U(4,N1)+RN2*U(4,N2)+RN3*U(4,N3)
        
        !CCCC  ----> DEFINO VARIABLES PURAS
        VX=U2/U1
        VY=U3/U1
        ET=U4/U1 
        RMOD2=VX*VX+VY*VY
        TEMP=GM/FR*(ET-.5D0*RMOD) !FR=CTE. UNIVERSAL DE LOS GASES
        C=DSQRT(GAMA*FR*TEMP)
        
        !CCCC  ----> DEFINICION DE LAS MATRICES A1 Y A2
        !CCCC  ----> A1 
        A1_21= GM/2.D0*RMOD2-VX*VX ; A1_22=(3.D0-GAMA)*VX ; A1_23=-GM*VY ; A1_24=GM
        A1_31= -VX*VY ; A1_32=VY ; A1_33=VX
        A1_41=(GM*RMOD2-GAMA*ET)*VX ; A1_42=GAMA*ET-GM/2.D0*RMOD2-GM*VX*VX ; A1_43=-GM*VX*VY ; A1_44=GAMA*VX
        !CCCC  ----> A2
        A2_21=-VX*VY ; A2_22=VY ; A2_23=VX
        A2_31=GM/2.D0*RMOD2-VY*VY ; A2_32=-GM*VX ; A2_33=(3.D0-GAMA)*VY ; A2_34=GM
        A2_41=(GM*RMOD2-GAMA*ET)*VY ; A2_42=-GM*VX*VY ; A2_43=GAMA*ET-GM/2.D0*RMOD2-GM*VY*VY ; A2_44=GAMA*VY 
        
        !CCCC------------------------------------------------------------------------------------------------------
        !CCCC  ----> MULTIPLICO POR PARTES PARA SIMPLIFICAR EL ASUNTO
        !CCCC  ----> 'A' POR LAS DERIVADAS
        AUXA1=          +       DRUX                                       +                    DRVY
        AUXA2=A1_21*DRX + A1_22*DRUX + A1_23*DRVX + A1_24*DREX + A2_21*DRY + A2_22*DRUY + A2_23*DRVY
        AUXA3=A1_31*DRX + A1_32*DRUX + A1_33*DRVX +              A2_31*DRY + A2_32*DRUY + A2_33*DRVY + A2_34*DREY
        AUXA4=A1_41*DRX + A1_42*DRUX + A1_43*DRVX + A1_44*DREX + A2_41*DRY + A2_42*DRUY + A2_43*DRVY + A2_44*DREY
        
        !CCCC  ----> ENSAMBLE DEL RIGHT HAND SIDE DEL NODO N1
        UN(1,N1)=UN(1,N1)+AUXA1*RN1*AR
        UN(2,N1)=UN(2,N1)+AUXA2*RN1*AR 
        UN(3,N1)=UN(3,N1)+AUXA3*RN1*AR 
        UN(4,N1)=UN(4,N1)+AUXA4*RN1*AR 
        !CCCC  ----> ENSAMBLE DEL RIGHT HAND SIDE DEL NODO N2
        UN(1,N2)=UN(1,N2)+AUXA1*RN2*AR 
        UN(2,N2)=UN(2,N2)+AUXA2*RN2*AR 
        UN(3,N2)=UN(3,N2)+AUXA3*RN2*AR 
        UN(4,N2)=UN(4,N2)+AUXA4*RN2*AR 
        !CCCC  ----> ENSAMBLE DEL RIGHT HAND SIDE DEL NODO N3
        UN(1,N3)=UN(1,N3)+AUXA1*RN3*AR 
        UN(2,N3)=UN(2,N3)+AUXA2*RN3*AR 
        UN(3,N3)=UN(3,N3)+AUXA3*RN3*AR 
        UN(4,N3)=UN(4,N3)+AUXA4*RN3*AR 
        
     END DO
  END DO
  
  DO INOD=1,NNOD
     DO J=1,4
        UN(J,INOD)=-UN(J,INOD)/M(INOD)
     END DO
  END DO
  
  RETURN
END SUBROUTINE CUARTO_ORDEN

!Calculo de los terminos de estabilizacion (Tezduyar)
SUBROUTINE ESTAB(NNOD,NELEM,N,DNX,DNY,U,T &
     ,VEL_X,VEL_Y,W_X,W_Y,GAMA,FR,RMU &
     ,DTMIN,RHOINF,TINF,UINF,VINF &
     ,T_SUGN1,T_SUGN2,T_SUGN3,SHOC,GAMM)
  
  IMPLICIT REAL(8) (A-H,O-Z)
  
  INTEGER N(3,NELEM)
  
  REAL(8) DNX(3,NELEM),DNY(3,NELEM),U(4,NNOD),T(NNOD),GAMM(NNOD)
  REAL(8) T_SUGN1(NELEM),T_SUGN2(NELEM),T_SUGN3(NELEM)
  REAL(8) VEL_X(NNOD),VEL_Y(NNOD),SHOC(NELEM)
  REAL(8) W_X(NNOD),W_Y(NNOD)
  
  !CCCC  ----> CTES DE CALCULO PARA SHOCK-CAPTURING      
  CC=DSQRT(GAMA*FR*TINF)
  VEL2=DSQRT(UINF*UINF+VINF*VINF)  
  
  DO IELEM=1,NELEM
     
     N1=N(1,IELEM)
     N2=N(2,IELEM)
     N3=N(3,IELEM)
     GM=(GAMM(N1)+GAMM(N2)+GAMM(N3))/3.D0
     TAU=0.D0
     H_RGNE=0.D0
     H_RGN=0.D0
     H_JGN=0.D0
     
     !CCCC  ----> VARIABLES ELEMENTALES
     RHO_ELEM=(U(1,N1)+U(1,N2)+U(1,N3))/3.D0
     VX=(VEL_X(N1)+VEL_X(N2)+VEL_X(N3))/3.D0
     VY=(VEL_Y(N1)+VEL_Y(N2)+VEL_Y(N3))/3.D0
     
     !PARTES NUEVAS
     WX=(W_X(N1)+W_X(N2)+W_X(N3))/3.D0
     WY=(W_Y(N1)+W_Y(N2)+W_Y(N3))/3.D0
     VX=VX-WX ; VY=VY-WY
     VEL2=DSQRT(VX*VX+VY*VY)
     
     !CCCC  ----> DERIVADA DE RHO           
     DRX= U(1,N1)*DNX(1,IELEM)+U(1,N2)*DNX(2,IELEM)+U(1,N3)*DNX(3,IELEM)
     DRY= U(1,N1)*DNY(1,IELEM)+U(1,N2)*DNY(2,IELEM)+U(1,N3)*DNY(3,IELEM)
     DR2=DSQRT(DRX*DRX+DRY*DRY)+1.D-20
     !CCCC  ----> DERIVADA DE TEMPERATURA
     DTX= T(N1)*DNX(1,IELEM)+T(N2)*DNX(2,IELEM)+T(N3)*DNX(3,IELEM)
     DTY= T(N1)*DNY(1,IELEM)+T(N2)*DNY(2,IELEM)+T(N3)*DNY(3,IELEM)
     DT2=DSQRT(DTX*DTX+DTY*DTY)+1.D-20
     !CCCC  ----> DERIVADA DE LA VELOCIDAD 
     DUX= VEL2*DNX(1,IELEM)+VEL2*DNX(2,IELEM)+VEL2*DNX(3,IELEM)
     DUY= VEL2*DNY(1,IELEM)+VEL2*DNY(2,IELEM)+VEL2*DNY(3,IELEM)
     
     DU2=DSQRT(DUX*DUX+DUY*DUY)+1.D-20
     
     !CCCC  ----> VECTOR UNIDAD THETA
     RTX=DTX/DT2
     RTY=DTY/DT2
     !CCCC  ----> VECTOR UNIDAD J
     RJX=DRX/DR2
     RJY=DRY/DR2
     !CCCC  ----> VECTOR UNIDAD VELOCIDAD 
     RUX=DUX/DU2
     RUY=DUY/DU2  
     
     TEMP=(T(N1)+T(N2)+T(N3))/3.D0
     C=DSQRT(GM*FR*TEMP)
     
     FMU= 1.716d-5*162.6/(TEMP-110.55)*(TEMP/273.15)**.75D0 !SUTHERLAND
     DO I=1,3
        
        TERM_1=DABS(VX*DNX(I,IELEM)+VY*DNY(I,IELEM))
        TERM_2=DABS(RJX*DNX(I,IELEM)+RJY*DNY(I,IELEM))
        
        H_RGN1=DABS(RTX*DNX(I,IELEM)+RTY*DNY(I,IELEM)) !CALCULO PARA ECU. ENERGIA
        H_RGN2=DABS(RUX*DNX(I,IELEM)+RUY*DNY(I,IELEM)) !CALCULO PARA ECU. MOMENTO
        
        TAU=TAU+TERM_1+TERM_2*VEL2
        
        H_RGNE=H_RGNE+H_RGN1
        H_RGN=H_RGN+H_RGN2
        H_JGN=H_JGN+TERM_2
     END DO
     
     TAU=1.D0/TAU
     H_RGNE=2.D0/H_RGNE
     
     H_RGN=2.D0/H_RGN
     IF(H_RGN.GT.1.D3)H_RGN=0.D0
     H_JGN=2.D0/H_JGN       
     IF(H_JGN.GT.1.D3)H_JGN=0.D0
     TR=H_JGN/4.D0
     TR1=DR2*H_JGN/RHO_ELEM
     ZZZ=H_JGN/(2.D0*C)
     SHOC(IELEM)=(DSQRT(TR1)+TR1**2)*.5D0*C**2*ZZZ
     
     RESUMEN=1.D0/TAU**2.D0 +(2.D0/DTMIN)**2.D0
     RRR=RESUMEN**(-.5D0)  
     T_SUGN1(IELEM)=RRR
     T_SUGN2(IELEM)=RRR
     T_SUGN3(IELEM)=RRR
     
     IF(FMU.NE.0.D0)THEN
        TAU_SUNG3=   H_RGN**2.D0/(4.D0*FMU/RHOINF)
        TAU_SUNG3_E= H_RGNE**2.D0/(4.D0*FMU/RHOINF)
        
        T_SUGN2(IELEM)=(RESUMEN+1.D0/TAU_SUNG3**2.D0)**(-.5D0)
        T_SUGN3(IELEM)=(RESUMEN+1.D0/TAU_SUNG3_E**2.D0)**(-.5D0)
     END IF
     
  END DO
  
  RETURN
END SUBROUTINE ESTAB
      
!Calcular la integral de los terminos de estabilizacion.
SUBROUTINE TODO(NNOD,NELEM,N,DNX,DNY,AREA,HHX,HHY &
     ,T_SUGN1,T_SUGN2,T_SUGN3,SHOC,PS,DTL &
     ,U,UN,RHS,P,GAMM,FR,RMU,FK,FCV,TINF,CTE)
  
  IMPLICIT REAL(8) (A-H,O-Z)
  
  INTEGER N(3,NELEM)
  
  REAL(8) DNX(3,NELEM),DNY(3,NELEM),U(4,NNOD),P(NNOD),UN(4,NNOD),T(NNOD),GAMM(NNOD)
  REAL(8) AREA(NELEM),RHS(4,NNOD),HHX(NELEM),HHY(NELEM)
  
  REAL(8) ALF(3),BET(3)
  REAL(8) T_SUGN1(NELEM),T_SUGN2(NELEM),T_SUGN3(NELEM),SHOC(NELEM)
  REAL(8) PS(NNOD),DTL(NELEM)
  
  !C---->     RHO=U(1)
  !C---->     RHO*VEL_X=U(2)
  !C---->     RHO*VEL_Y=U(3)
  !C---->     RHO*ET=U(4)
  
  DATA ALF/.5D0,.5D0,0.D0/
  DATA BET/0.D0,.5D0,.5D0/
  
  NGAUSS=3    !PTOS DE GAUSS DONDE VOY A INTERGRAR
  
  DO IELEM=1,NELEM
     
     N1=N(1,IELEM) ; N2=N(2,IELEM) ; N3=N(3,IELEM)
     
     GAMA=(GAMM(N1)+GAMM(N2)+GAMM(N3))/3.D0
     GM=GAMA-1.D0
     TEMP=(T(N1)+T(N2)+T(N3))/3.D0
     FMU= 1.716d-5*162.6/(TEMP-110.55)*(TEMP/273.15)**.75D0     !SUTHERLAND
     
     RNX1=DNX(1,IELEM) ; RNX2=DNX(2,IELEM) ; RNX3=DNX(3,IELEM)
     RNY1=DNY(1,IELEM) ; RNY2=DNY(2,IELEM) ; RNY3=DNY(3,IELEM)
     
     !CCCC  ----> DERIVADA DE RHO           
     DRX= RNX1*U(1,N1)+RNX2*U(1,N2)+RNX3*U(1,N3)
     DRY= RNY1*U(1,N1)+RNY2*U(1,N2)+RNY3*U(1,N3)
     !CCCC  ----> DERIVADA DE RHO.VEL_X           
     DRUX= RNX1*U(2,N1)+RNX2*U(2,N2)+RNX3*U(2,N3)
     DRUY= RNY1*U(2,N1)+RNY2*U(2,N2)+RNY3*U(2,N3)
     !CCCC  ----> DERIVADA DE RHO.VEL_Y           
     DRVX= RNX1*U(3,N1)+RNX2*U(3,N2)+RNX3*U(3,N3)
     DRVY= RNY1*U(3,N1)+RNY2*U(3,N2)+RNY3*U(3,N3)
     !CCCC  ----> DERIVADA DE RHO.ET           
     DREX= RNX1*U(4,N1)+RNX2*U(4,N2)+RNX3*U(4,N3)
     DREY= RNY1*U(4,N1)+RNY2*U(4,N2)+RNY3*U(4,N3)
     
     AR=AREA(IELEM)*DTL(IELEM)/3.D0
     !AR=AREA(IELEM)/3.D0
     !CCCC  ----> LONG. CARACTERISTICA
     HLONGX=HHX(IELEM)
     HLONGY=HHY(IELEM)
     HLONG=DSQRT(AREA(IELEM))
     !CCCC  ----> ESTAB. TEZDUYAR
     TAU1=T_SUGN1(IELEM)
     TAU2=T_SUGN2(IELEM)
     TAU3=T_SUGN3(IELEM)
     ALFA_MU=SHOC(IELEM)
         
     DO J=1,NGAUSS
        
        RN1=1.D0-ALF(J)-BET(J)
        RN2=ALF(J)
        RN3=BET(J)
        
        !CCCC  ----> INTEGRO LAS VARIABLES EN LOS PUNTOS DE GAUSS
        U1= RN1*U(1,N1)+RN2*U(1,N2)+RN3*U(1,N3)
        U2= RN1*U(2,N1)+RN2*U(2,N2)+RN3*U(2,N3)
        U3= RN1*U(3,N1)+RN2*U(3,N2)+RN3*U(3,N3)
        U4= RN1*U(4,N1)+RN2*U(4,N2)+RN3*U(4,N3)

        FI_1= RN1*UN(1,N1)+RN2*UN(1,N2)+RN3*UN(1,N3)
        FI_2= RN1*UN(2,N1)+RN2*UN(2,N2)+RN3*UN(2,N3)
        FI_3= RN1*UN(3,N1)+RN2*UN(3,N2)+RN3*UN(3,N3)
        FI_4= RN1*UN(4,N1)+RN2*UN(4,N2)+RN3*UN(4,N3)
        !CCCC  ----> DEFINO VARIABLES PRIMITIVAS
        VX=U2/U1
        VY=U3/U1
        ET=U4/U1
        RMOD2=VX*VX+VY*VY
        TEMP=GM/FR*(ET-.5D0*RMOD2) !FR=CTE. UNIVERSAL DE LOS GASES
        
        C=DSQRT(GAMA*FR*TEMP)
        
        !CCCC  ----> VARIABLES COMPACTAS
        FMU1=FMU-(FK/FCV)   ! FK=CONDUCTIVIDAD TERMICA; FCV=CALOR ESPECIFICO A VOL. CTE
        FMU43=4.D0/3.D0*FMU ! FMU=VISCOSIDAD
        FMU23=2.D0/3.D0*FMU
        RU1=AR/U1           ! RU1= FACTOR COMUN DEL JACOBIANO DE LAS MATRICES K
        
        !CCCC  ----> DEFINICION DE LAS MATRICES A1 Y A2

        !CCCC  ----> A1
        A1_21= GM/2.D0*RMOD2-VX*VX ; A1_22=(3.D0-GAMA)*VX ; A1_23=-GM*VY ; A1_24=GM
        A1_31= -VX*VY ; A1_32=VY ; A1_33=VX
        A1_41=(GM*RMOD2-GAMA*ET)*VX ; A1_42=GAMA*ET-GM/2.D0*RMOD2-GM*VX*VX ; A1_43=-GM*VX*VY ; A1_44=GAMA*VX
        !CCCC  ----> A2
        A2_21=-VX*VY ; A2_22=VY ; A2_23=VX
        A2_31=GM/2.D0*RMOD2-VY*VY ; A2_32=-GM*VX ; A2_33=(3.D0-GAMA)*VY ; A2_34=GM
        A2_41=(GM*RMOD2-GAMA*ET)*VY ; A2_42=-GM*VX*VY ; A2_43=GAMA*ET-GM/2.D0*RMOD2-GM*VY*VY ; A2_44=GAMA*VY 
       
        !CCCC----------------------------------------
        !CCCC  ----> NO CALCULO LOS TERMINOS VISCOSOS CUANDO MU=0.0
        IF(FMU.GT.0.D0)THEN
           !CCCC  ----> DEFINICION DE LAS MATRICES K11, K12, K21 Y K22
           
           !CCCC  ----> K11
           RK11_21=-FMU43*VX ; RK11_22=FMU43
           RK11_31=-FMU*VY ; RK11_33=FMU 
           RK11_41=FK/FCV*(.5D0*RMOD2-FCV*TEMP)-FMU*RMOD2-FMU/3.D0*VX*VX ; RK11_42=(FMU/3.D0+FMU1)*VX
           RK11_43=FMU1*VY ; RK11_44=FK/FCV
           !CCCC  ----> K12
           RK12_21=FMU23*VY ; RK12_23=-FMU23
           RK12_31=-FMU*VX ; RK12_32=FMU
           RK12_41=-FMU/3.D0*VX*VY ; RK12_42=FMU*VY ; RK12_43=-FMU23*VX
           !CCCC  ----> K21
           RK21_21=-FMU*VY ; RK21_23=FMU
           RK21_31=FMU23*VX ; RK21_32=-FMU23
           RK21_41=-FMU/3.D0*VX*VY ; RK21_42=-FMU23*VY ; RK21_43=FMU*VX
           !CCCC  ----> K22
           RK22_21=-FMU*VX ; RK22_22=FMU
           RK22_31=-FMU43*VY ; RK22_33=FMU43
           RK22_41=FK/FCV*(.5D0*RMOD2-FCV*TEMP)-FMU*RMOD2-FMU/3.D0*VY*VY ; RK22_42=FMU1*VX
           RK22_43=(FMU/3.D0+FMU1)*VY ; RK22_44=FK/FCV

           !CCCC  ----> TERMINOS VISCOSOS
           VV2=RK11_21*DRX + RK11_22*DRUX + RK12_21*DRY + RK12_23*DRVY
           VV3=RK11_31*DRX + RK11_33*DRVX + RK12_31*DRY + RK12_32*DRUY
           VV4=RK11_41*DRX + RK11_42*DRUX + RK11_43*DRVX + RK11_44*DREX + RK12_41*DRY + RK12_42*DRUY + RK12_43*DRVY
           
           VV6=RK21_21*DRX + RK21_23*DRVX + RK22_21*DRY + RK22_22*DRUY
           VV7=RK21_31*DRX + RK21_32*DRUX + RK22_31*DRY + RK22_33*DRVY
           VV8=RK21_41*DRX + RK21_42*DRUX + RK21_43*DRVX + RK22_41*DRY + RK22_42*DRUY + RK22_43*DRVY + RK22_44*DREY
           
           !CCCC  ----> ENSAMBLE DEL RIGHT HAND SIDE DEL NODO N1
           
           RHS(2,N1)=RHS(2,N1)+RU1*(RNX1*VV2+RNY1*VV6)
           RHS(3,N1)=RHS(3,N1)+RU1*(RNX1*VV3+RNY1*VV7) 
           RHS(4,N1)=RHS(4,N1)+RU1*(RNX1*VV4+RNY1*VV8)
           !CCCC  ----> ENSAMBLE DEL RIGHT HAND SIDE DEL NODO N2
           
           RHS(2,N2)=RHS(2,N2)+RU1*(RNX2*VV2+RNY2*VV6)
           RHS(3,N2)=RHS(3,N2)+RU1*(RNX2*VV3+RNY2*VV7)
           RHS(4,N2)=RHS(4,N2)+RU1*(RNX2*VV4+RNY2*VV8)
           !CCCC  ----> ENSAMBLE DEL RIGHT HAND SIDE DEL NODO N3
           
           RHS(2,N3)=RHS(2,N3)+RU1*(RNX3*VV2+RNY3*VV6)
           RHS(3,N3)=RHS(3,N3)+RU1*(RNX3*VV3+RNY3*VV7)
           RHS(4,N3)=RHS(4,N3)+RU1*(RNX3*VV4+RNY3*VV8)    
        END IF
        
        !CCCC-----------------------------------------------------------------------------------------------------------------
        !CCCC----> TERMINOS DE ESTABILIZACION Y CAPTURA DE CHOQUE                        
        ARR1=AR*TAU1
        ARR2=AR*TAU2
        ARR3=AR*TAU3
        
        !CCCC---->************************************<----!CCCC 
        !CCCC---->  CALCULO DE MUa Y SUS COMPONENTES  <----!CCCC
        !CCCC---->************************************<----!CCCC
        V11=C+DABS(VX)+DABS(VY) 
        
        !CCCC  ----> ARMO EL PRESSURE SWITCH
        
        CHOQ1=alfa_mu*AR*PS(N1)*CTE
        CHOQ2=alfa_mu*AR*PS(N2)*CTE
        CHOQ3=alfa_mu*AR*PS(N3)*CTE
        !CCCC----------------------------------------------------------------------------------------------------
        
        
        !CCCC----------------------------------------------------------------------------------------------------
        !CCCC  ----> MULTIPLICO POR PARTES PARA SIMPLIFICAR EL ASUNTO
        !CCCC  ----> 'A' POR LAS DERIVADAS
        AUXA1=                  DRUX                                       +                    DRVY
        AUXA2=A1_21*DRX + A1_22*DRUX + A1_23*DRVX + A1_24*DREX + A2_21*DRY + A2_22*DRUY + A2_23*DRVY
        AUXA3=A1_31*DRX + A1_32*DRUX + A1_33*DRVX +              A2_31*DRY + A2_32*DRUY + A2_33*DRVY + A2_34*DREY
        AUXA4=A1_41*DRX + A1_42*DRUX + A1_43*DRVX + A1_44*DREX + A2_41*DRY + A2_42*DRUY + A2_43*DRVY + A2_44*DREY
        
        AUXA11=AUXA1+FI_1
        AUXA22=AUXA2+FI_2
        AUXA33=AUXA3+FI_3
        AUXA44=AUXA4+FI_4
        !CCCC  ----> LO ANTERIOR POR 'A' TRANSPUESTA

        AA1=                     AUXA22  
        AA2=A1_21*AUXA11 + A1_22*AUXA22 + A1_23*AUXA33 + A1_24*AUXA44                  
        AA3=A1_31*AUXA11 + A1_32*AUXA22 + A1_33*AUXA33 
        AA4=A1_41*AUXA11 + A1_42*AUXA22 + A1_43*AUXA33 + A1_44*AUXA44
        AA5=                                    AUXA33  
        AA6=A2_21*AUXA11 + A2_22*AUXA22 + A2_23*AUXA33 
        AA7=A2_31*AUXA11 + A2_32*AUXA22 + A2_33*AUXA33 + A2_34*AUXA44
        AA8=A2_41*AUXA11 + A2_42*AUXA22 + A2_43*AUXA33 + A2_44*AUXA44
        
       
        !CCCC  ----> ENSAMBLE DEL RIGHT HAND SIDE DEL NODO N1
        RHS(1,N1)=RHS(1,N1)+(RNX1*AA1+RNY1*AA5)*ARR1 +AUXA1*RN1*AR +(RNX1*DRX+RNY1*DRY)*CHOQ1
        RHS(2,N1)=RHS(2,N1)+(RNX1*AA2+RNY1*AA6)*ARR2 +AUXA2*RN1*AR +(RNX1*DRUX+RNY1*DRUY)*CHOQ1
        RHS(3,N1)=RHS(3,N1)+(RNX1*AA3+RNY1*AA7)*ARR2 +AUXA3*RN1*AR +(RNX1*DRVX+RNY1*DRVY)*CHOQ1
        RHS(4,N1)=RHS(4,N1)+(RNX1*AA4+RNY1*AA8)*ARR3 +AUXA4*RN1*AR +(RNX1*DREX+RNY1*DREY)*CHOQ1
        !CCCC  ----> ENSAMBLE DEL RIGHT HAND SIDE DEL NODO N2
        RHS(1,N2)=RHS(1,N2)+(RNX2*AA1+RNY2*AA5)*ARR1 +AUXA1*RN2*AR +(RNX2*DRX+RNY2*DRY)*CHOQ2
        RHS(2,N2)=RHS(2,N2)+(RNX2*AA2+RNY2*AA6)*ARR2 +AUXA2*RN2*AR +(RNX2*DRUX+RNY2*DRUY)*CHOQ2
        RHS(3,N2)=RHS(3,N2)+(RNX2*AA3+RNY2*AA7)*ARR2 +AUXA3*RN2*AR +(RNX2*DRVX+RNY2*DRVY)*CHOQ2
        RHS(4,N2)=RHS(4,N2)+(RNX2*AA4+RNY2*AA8)*ARR3 +AUXA4*RN2*AR +(RNX2*DREX+RNY2*DREY)*CHOQ2
        !CCCC  ----> ENSAMBLE DEL RIGHT HAND SIDE DEL NODO N3
        RHS(1,N3)=RHS(1,N3)+(RNX3*AA1+RNY3*AA5)*ARR1 +AUXA1*RN3*AR +(RNX3*DRX+RNY3*DRY)*CHOQ3
        RHS(2,N3)=RHS(2,N3)+(RNX3*AA2+RNY3*AA6)*ARR2 +AUXA2*RN3*AR +(RNX3*DRUX+RNY3*DRUY)*CHOQ3
        RHS(3,N3)=RHS(3,N3)+(RNX3*AA3+RNY3*AA7)*ARR2 +AUXA3*RN3*AR +(RNX3*DRVX+RNY3*DRVY)*CHOQ3
        RHS(4,N3)=RHS(4,N3)+(RNX3*AA4+RNY3*AA8)*ARR3 +AUXA4*RN3*AR +(RNX3*DREX+RNY3*DREY)*CHOQ3
     END DO
  END DO
  
  RETURN
END SUBROUTINE TODO

!Interpolacion para gas en equilibrio
      SUBROUTINE TGAS(E,RHO,P,A,T,GAMM)
      IMPLICIT REAL(8) (A-H,O-Z)

      Y2=DLOG10(RHO/1.292D0)
      Z2=DLOG10(E/78408.4D0)
      IF(Y2.GT.-.50D0) GO TO 11
      IF(Y2.GT.-4.50D0)GO TO 6
      IF(Z2.GT..65D0)  GO TO 1
      GAMM=1.4D0
      SNDSQ=E*.560D0
      GO TO 18
 1    IF(Z2.GT.1.50D0) GO TO 2
      GAMM=1.46543D0+(.007625D0+.000292D0*Y2)*Y2-(.254500D0+.017244D0*Y2)*Z2 &
          +(.355907D0+.015422D0*Y2-.163235D0*Z2)*Z2*Z2
      GAME=2.304D0*(-.25450D0-.017244D0*Y2+(.711814D0+.030844D0*Y2-.489705D0*Z2)*Z2)
      GAMR=2.304D0*(.007625D0+(-.017244D0+.015422D0*Z2)*Z2+.000584D0*Y2) 
      A1=-.000954D0
      A2=.171187D0
      A3=.004567D0
      GO TO 17
 2    IF(Z2.GT.2.20D0) GO TO 3
      GAS1=2.02636D0+.0584931D0*Y2
      GAS2=.454886D0+.027433D0*Y2
      GAS3=.165265D0+.014275D0*Y2
      GAS4=.136685D0+.010071D0*Y2
      GAS5=.058493D0-.027433D0*Z2
      GAS6=-.014275D0+.010071D0*Z2
      GAS7=DEXP((-10.D0+3.D0*Y2)*(Z2-.023D0*Y2-2.025D0))
      DERE=-30.0D0
      DERR=0.285D0
      A1=.008737D0
      A2=.184842D0
      A3=-.302441D0
      GO TO 15
 3    IF(Z2.GT.3.05D0) GO TO 4
      GAS1=1.60804D0+.034791D0*Y2
      GAS2=.188906D0+.010927D0*Y2
      GAS3=.124117D0+.007277D0*Y2
      GAS4=.069839D0+.003985D0*Y2
      GAS5=.034791D0-.010927D0*Z2
      GAS6=-.07277D0+.03985D0*Z2
      GAS7=DEXP(-30.0D0*(Z2+.007D0*Y2-2.691D0))
      DERE=-30.0D0
      DERR=0.21D0
      A1=.017884D0
      A2=.153672D0
      A3=-.930224D0
      GO TO 15
 4    IF(Z2.GT.3.38D0) GO TO 5
      GAS1=1.25672D0+.007073D0*Y2
      GAS2=.039228D0-.000491D0*Y2
      GAS3=-.721798D0-.073753D0*Y2
      GAS4=-.198942D0-.021539D0*Y2
      GAS5=.007073D0+.000491D0*Z2
      GAS6=.073753D0-.021539D0*Z2
      GAS7=DEXP(0.425D0*Y2-50.0D0*Z2+166.7D0)
      DERE=-50.D0
      DERR=0.325D0
      A1=.002379D0
      A2=.217959D0
      A3=.005943D0
      GO TO 15
 5    GAMM=-84.0327D0+(-.331761D0+.001153D0*Y2)*Y2+(72.2066D0+.491914D0*Y2)*Z2 &
          +(-20.3559D0-.070617D0*Y2+1.90979D0*Z2)*Z2*Z2
      GAME=2.304D0*( 72.2066D0+.491914D0*Y2+(-40.7118D0-.141234D0*Y2+5.72937D0*Z2)*Z2)
      GAMR=2.304D0*(-.831761D0+.002306D0*Y2+(.491914D0-.070617D0*Z2)*Z2)
      A1=.006572D0
      A2=.183396D0
      A3=-.135960D0
      GO TO 17
 6    IF(Z2.GT..65D0) GO TO 7
      GAMM=1.4D0
      SNDSQ=E*.560D0
      GO TO 18
 7    IF(Z2.GT.1.54D0) GO TO 8
      GAS1=1.44813D0+.001292D0*Y2
      GAS2=.073510D0+.001948D0*Y2
      GAS3=-.054745D0+.013705D0*Y2
      GAS4=-.055473D0+.021874D0*Y2
      GAS5=.001292D0-.001948D0*Z2
      GAS6=-.013705D0+.021874D0*Z2
      GAS7=DEXP(-10.0D0*(Z2-1.42D0))
      DERE=-1.D0
      DERR=0.D0
      A1=-.001973D0
      A2=.185233D0
      A3=-.059952D0
      GO TO 15
  8   IF(Z2.GT.2.22D0) GO TO 9
      GAS1=1.73158D0+.003902D0*Y2
      GAS2=.272846D0-.006237D0*Y2
      GAS3=-.041419D0-.037475D0*Y2
      GAS4=.016984D0-.018038D0*Y2
      GAS5=.003902D0+.006237D0*Z2
      GAS6=.037475D0-.018038D0*Z2
      GAS7=DEXP((-10.0D0+3.0D0*Y2)*(Z2-.023D0*Y2-2.025D0))
      DERE=3.D0*Y2-10.D0
      DERR=3.D0*Z2+12.15D0*Y2-20.325D0
      A1=-.013027D0
      A2=.07427D0
      A3=.012889D0
      GO TO 15
 9    IF(Z2.GT.2.90D0) GO TO 10
      GAS1=1.59350D0+.075324D0*Y2
      GAS2=.176186D0+.026072D0*Y2
      GAS3=.200838D0+.058536D0*Y2
      GAS6=.099687D0+.025287D0*Y2 !VIENE POR ACA
      GAS5=.075324D0-.026072D0*Z2
      GAS6=-.058536D0+.025287D0*Z2
      GAS7=DEXP((-10.D0+5.D0*Y2)*(Z2-2.7D0))
      DERE=5.D0*Y2-10.D0
      DERR=5.D0*Z2-13.5D0
      A1=.004342D0
      A2=.212192D0
      A3=-.001293D0
      GO TO 15
 10   GAS1=1.12688D0-.025957D0*Y2
      GAS2=-.013602D0-.013772D0*Y2
      GAS3=.127737D0+.087942D0*Y2
      GAS4=.043104D0+.023547D0*Y2
      GAS5=-.025957D0+.013772D0*Z2
      GAS6=-.087942D0+.023547D0*Z2
!C      GAS7=EXP(-20.0D0*Z2+(4.0D0*Z2-13.2D0)*Y2+66.D0)
      GAS7=DEXP((-20.D0+4.D0*Y2)*(Z2-3.3D0))
      DERE=-20.D0+4.D0*Y2
      DERR=4.D0*Z2-13.2D0
      A1=.006348D0
      A2=.209716D0
      A3=-.006001D0
      GO TO 15
 11   IF(Z2.GT..65D0) GO TO 12
      GAMM=1.4D0
      SNDSQ=E*.560D0
      GO TO 18
 12   IF(Z2.GT.1.68D0) GO TO 13
      GAS1=1.45510D0-.000102D0*Y2
      GAS2=.081537D0-.000166D0*Y2
      GAS3=-.128647D0+.049454D0*Y2
      GAS4=-.101036D0+.033518D0*Y2
      GAS5=-.000102D0+.000166D0*Z2
      GAS6=-.049454D0+.033518D0*Z2
      GAS7=DEXP(-15.D0*(Z2-1.420D0))
      DERE=-15.D0
      DERR=0.D0
      A1=.00045D0
      A2=.203892D0
      A3=.101797D0
      GO TO 15
 13   IF(Z2.GT.2.46D0) GO TO 14
      GAS1=1.59608D0-.042426D0*Y2
      GAS2=.192840D0-.029353D0*Y2
      GAS3=.019430D0-.005954D0*Y2
      GAS4=.026097D0-.006164D0*Y2
      GAS5=-.042426D0+.029353D0*Z2
      GAS6=.005954D0-.006164D0*Z2
      GAS7=DEXP(-15.D0*(Z2-2.050D0))
      DERE=-15.D0
      DERR=0.D0
      A1=-.006609D0
      A2=.127637D0
      A3=.297037D0
      GO TO 15
 14   GAS1=1.54363D0-.049071D0*Y2
      GAS2=.153562D0-.029209D0*Y2
      GAS3=.324907D0+.077599D0*Y2
      GAS4=.142408D0+.022071D0*Y2
      GAS5=-.049071D0+.029209D0*Z2
      GAS6=-.077599D0+.022071D0*Z2
      GAS7=DEXP(-10.0D0*(Z2-2.708D0))
      DERE=-10.D0
      DERR=0.D0
      A1=-.000081D0
      A2=.226601D0
      A3=.170922D0
 15   GAS10=1.D0/(1.+GAS7)
 16   GAS8=GAS3-GAS4*Z2
      GAS8=GAS3-GAS4*Z2
      GAS9=GAS8*GAS7*GAS10*GAS10
      GAMM=GAS1-GAS2*Z2-GAS8*GAS10
      GAME=2.304D0*(-GAS2+GAS4*GAS10+GAS9*DERE)
      GAMR=2.304D0*(GAS5+GAS6*GAS10+GAS9*DERR)
 17   SNDSQ=E*(A1+(GAMM-1.D0)*(GAMM+A2*GAME)+A3*GAMR)
 18   A=DSQRT(SNDSQ)
      P=RHO*E*(GAMM-1.D0)
      X2=DLOG10(P/1.0134D5)
      Y2=Y2+.0231264D0
      Z3=X2-Y2
      IF(Y2.GT.-.50D0) GO TO 29
      IF(Y2.GT.-4.50D0)GO TO 24
      IF(Z3.GT..30D0)  GO TO 19
      T=P/(287.D0*RHO)
      RETURN
 19   IF(Z3.GT.1.0D0) GO TO 20
      T=10.D0**(.2718D0+.00074D0*Y2+(.990136D0-.004947D0*Y2)*Z3+(.990717D0 &
          +.175194D0*Y2-(.982407D0+.159233D0*Y2)*Z3)/(1.D0+DEXP(-20.D0*(Z3-.88D0))))
      GO TO 32
 20   IF(Z3.GT.1.35D0) GO TO 21
      T=10.D0**(1.39925D0+.167780D0*Y2+(-.143168D0-.159234D0*Y2)*Z3+(-.027614D0 &
          -.090761D0*Y2+(.307036D0+.121621D0*Y2)*Z3)/(1.D0+DEXP(-20.D0*(Z3-1.17D0))))
      GO TO 32
 21   IF(Z3.GT.1.79D0) GO TO 22
      T=10.D0**(1.11401D0+.002221D0*Y2+(.351875D0+.017246D0*Y2)*Z3+(-1.15099D0 &
          -.173555D0*Y2+(.673342D0-.088399D0*Y2)*Z3)/(1.D0+DEXP(-20.D0*(Z3-1.56D0))))
      GO TO 32

 22   IF(Z3.GT.2.47D0) GO TO 23
      T=10.D0**(1.01722D0-.017918D0*Y2+(.473523D0+.025456D0*Y2)*Z3+(-2.17978D0 &
          -.334716D0*Y2+(.898619D0+.127386D0*Y2)*Z3)/(1.D0+DEXP(-20.D0*(Z3-2.22D0))))
      GO TO 32
 23   T=10.D0**(-45.0871D0-9.00504D0*Y2+(35.8685D0+6.79222D0*Y2)*Z3-(6.77699D0 &
          +1.2737D0*Y2)*Z3*Z3+(-.064705D0+.025325D0*Z3)*Y2*Y2)
      GO TO 32
 24   IF(Z3.GT..48D0) GO TO 25
      T=P/(287.D0*RHO)
      RETURN
 25   IF(Z3.GT..9165D0) GO TO 26
      T=10.D0**(.284312D0+.987912D0*Z3+.001644D0*Y2)
      GO TO 32
 26   IF(Z3.GT.1.478D0) GO TO 27
      T=10.D0**(.502071D0-.01299D0*Y2+(.774818D0+.025397D0*Y2)*Z3+(.009912D0 &
          -.150527D0*Y2+(-.000385D0+.105734D0*Y2)*Z3)/(1.D0+DEXP(-15.D0*(Z3-1.28D0))))
      GO TO 32
 27   IF(Z3.GT.2.176D0) GO TO 28
      T=10.D0**(1.02294D0+.021535D0*Y2+(.427212D0+.0069D0*Y2)*Z3+(-.427823D0 &
          -.211991D0*Y2+(.257096D0+.101192D0*Y2)*Z3)/(1.D0+DEXP(-12.D0*(Z3-1.778D0))))
      GO TO 32
 28   T=10.D0**(1.47540D0+.12962D0*Y2+(.254154D0-.046411D0*Y2)*Z3+(-.221229D0 &
          -.057077D0*Y2+(.158116D0+.03043D0*Y2)*Z3)/(1.D0+DEXP(5.D0*Y2*(Z3-2.40D0))))
      GO TO 32
 29   IF(Z3.GT..48D0) GO TO 30
      T=P/(287.D0*RHO)
      RETURN
 30   IF(Z3.GT.1.07D0) GO TO 31
      T=10.D0**(.279268D0+.992172D0*Z3)
      GO TO 32
 31   T=10.D0**(.233261D0-.056383D0*Y2+(1.19783D0+.063121D0*Y2-.165985D0*Z3)*Z3+ &
          (-.814535D0+.099233D0*Y2+(.602385D0-.067428D0*Y2-.093991D0*Z3)*Z3)/(1.D0+DEXP(( &
          5.D0*Y2-20.D0)*(Z3-1.78D0))))
 32   T=T*151.777778D0
      RETURN
      END
    
!Hacer nulo la velocidad normal
SUBROUTINE NORMALVEL(NNOD,VEL_X,VEL_Y,W_X,W_Y,RHO,INORMV_NODE &
     ,RNORMV_VALUEX,RNORMV_VALUEY,NNORMV)
  
  IMPLICIT REAL(8) (A-H,O-Z)
  
  INTEGER NNOD,NNORMV
  INTEGER INORMV_NODE(NNORMV)
  REAL(8) VEL_X(NNOD),VEL_Y(NNOD)
  REAL(8) RNORMV_VALUEX(NNORMV),RNORMV_VALUEY(NNORMV)
  REAL(8) W_X(NNOD),W_Y(NNOD),RHO(NNOD)

  DO INORMV=1,NNORMV
     INOD=INORMV_NODE(INORMV)

     VX=VEL_X(INOD) ; VY=VEL_Y(INOD)
     WX=W_X(INOD)   ; WY=W_Y(INOD)

     RNX=RNORMV_VALUEX(INORMV) ; RNY=RNORMV_VALUEY(INORMV)
     
     UTP=-RNY*(VX-WX) + RNX*(VY-WY)
     
     VEL_X(INOD)=-RNY*UTP+WX
     VEL_Y(INOD)=RNX*UTP+WY

  END DO
  
  RETURN
END SUBROUTINE NORMALVEL

!Fijar velocidad
SUBROUTINE FIXVEL(NNOD,VEL_X,VEL_Y &
     ,NFIXV,IFIXV_NODE,RFIXV_VALUEX,RFIXV_VALUEY)
  
  IMPLICIT REAL(8) (A-H,O-Z)
  
  INTEGER IFIXV_NODE(NFIXV)
  
  REAL(8) VEL_X(NNOD),VEL_Y(NNOD)
  REAL(8) RFIXV_VALUEX(NFIXV),RFIXV_VALUEY(NFIXV)
  
  DO IFIXV=1,NFIXV
     N1=IFIXV_NODE(IFIXV)
     VEL_X(N1)=RFIXV_VALUEX(IFIXV)
     VEL_Y(N1)=RFIXV_VALUEY(IFIXV)
  END DO
  
  RETURN
END SUBROUTINE FIXVEL

!Fijar valores de T y RHO en contornos
SUBROUTINE FIX(NNOD,RHO,VEL_X,VEL_Y,T,E &
     ,NFIXRHO,IFIXRHO_NODE,RFIXRHO_VALUE &
     ,NFIXT,IFIXT_NODE,RFIXT_VALUE &
     ,FR,GAMM,GAMA)
  
  IMPLICIT REAL(8) (A-H,O-Z)
  
  INTEGER IFIXRHO_NODE(NFIXRHO),IFIXT_NODE(NFIXT)
  
  REAL(8) VEL_X(NNOD),VEL_Y(NNOD),T(NNOD),GAMM(NNOD)
  REAL(8) RHO(NNOD),E(NNOD)
  REAL(8) RFIXRHO_VALUE(NFIXRHO),RFIXT_VALUE(NFIXT)
  
  DO IFIXRHO=1,NFIXRHO
     INOD=IFIXRHO_NODE(IFIXRHO)
     RHO(INOD)=RFIXRHO_VALUE(IFIXRHO)
  END DO
  
  DO IFIXT=1,NFIXT
     INOD=IFIXT_NODE(IFIXT)
     GM=GAMM(INOD)-1.D0
     T(INOD)=RFIXT_VALUE(IFIXT)
     E(INOD)=T(INOD)*FR/GM+.5D0*(VEL_X(INOD)**2.D0+VEL_Y(INOD)**2.D0)
  END DO
  
  RETURN
END SUBROUTINE FIX

!Imprimir resultados en formato .flavia
SUBROUTINE PRINTFLAVIA(FR,GAMM,RHO,VEL_X,VEL_Y,P,T,E,RMACH,XPOS,YPOS &
     ,NNOD,ITER,MOVIE,FILE,ILONG)
  IMPLICIT REAL(8) (A-H,O-Z)
  
  INTEGER NNOD
  
  REAL(8) RHO(NNOD),VEL_X(NNOD),VEL_Y(NNOD)
  REAL(8) P(NNOD),T(NNOD),E(NNOD),GAMM(NNOD),RMACH(NNOD)
  REAL(8) XPOS(NNOD),YPOS(NNOD),U(4,NNOD)
  CHARACTER FILE*80
  
  !CCCC----> IMPRESION DE RESULTADOS
  !CCCC-----------------------------
  IF (MOVIE.EQ.0)THEN
     OPEN(2,FILE=FILE(1:ILONG)//'.flavia.res',STATUS='UNKNOWN')
  END IF
  !CCCC----> ESCRITURA DE VELOCIDADES
  !CCCC------------------------------      
  WRITE(2,'(A15,5I6)') 'VELOCITY',2,ITER,2,1,1
  WRITE(2,'(A)') 'VEL_X'
  WRITE(2,'(A)') 'VEL_Y'
  
  DO INOD=1,NNOD
     WRITE(2,'(I6,3F16.6)') INOD,VEL_X(INOD),VEL_Y(INOD)
  END DO

  !CCCC----> ESCRITURA DE POSICIONES
  !CCCC------------------------------      
  WRITE(2,'(A15,5I6)') 'POSITION',2,ITER,2,1,1
  WRITE(2,'(A)') 'X'
  WRITE(2,'(A)') 'Y'
  
  DO INOD=1,NNOD     
     WRITE(2,'(I6,3F16.6)') INOD,XPOS(INOD),YPOS(INOD)
  END DO
  
  !CCCC----> ESCRITURA DE LAS DENSIDADES
  !CCCC---------------------------------      
  WRITE(2,'(A15,5I6)') 'DENSITY',2,ITER,1,1,1
  WRITE(2,'(A)') 'DENSITY'
  
  DO INOD=1,NNOD
     WRITE(2,'(I6,1E16.8)') INOD,RHO(INOD)
  END DO
  
  !CCCC----> ESCRITURA DE LAS PRESIONES
  !CCCC--------------------------------
  WRITE(2,'(A15,5I6)') 'PRESSURE',2,ITER,1,1,1
  WRITE(2,'(A)') 'PRESSURE'
  
  DO INOD=1,NNOD
     WRITE(2,'(I6,1E16.8)') INOD,P(INOD)
  END DO
  
  !CCCC  -----> ESCRITURA DE LAS TEMPERATURAS
  !CCCC  ------------------------------------
  WRITE(2,'(A15,5I6)') 'TEMPERATURE',2,ITER,1,1,1
  WRITE(2,'(A)') 'TEMPERATURE'
  
  DO INOD=1,NNOD
     WRITE(2,'(I6,1E16.8)') INOD,T(INOD)   
  END DO

  !CCCC----> ESCRITURA DEL NUMERO DE MACH
  !CCCC-----------------------------------
  WRITE(2,'(A15,5I6)') 'Mach_Number',2,ITER,2,1,1
  WRITE(2,'(A)') 'Mach_Number'
  WRITE(2,'(A)') 'Mach_Number_Tgas'
  
  DO INOD=1,NNOD
     VEL=DSQRT(VEL_X(INOD)*VEL_X(INOD)+VEL_Y(INOD)*VEL_Y(INOD))
     VC=DSQRT(GAMM(INOD)*FR*T(INOD))
     WRITE(2,'(I6,2F16.8)') INOD,VEL/VC,RMACH(INOD)
     
  END DO
      
  !CCCC----> ESCRITURA DE LA ENERGIA INTERNA
  !CCCC-----------------------------------
  WRITE(2,'(A15,5I6)') 'Internal_Energy',2,ITER,1,1,1
  WRITE(2,'(A)') 'Internal_Energy'
  
  DO INOD=1,NNOD
     WRITE(2,'(I6,1E16.8)') INOD,E(INOD)
  END DO

  !CCCC----> ESCRITURA DE GAMA
  !CCCC-----------------------------------    
  WRITE(2,'(A15,5I6)') 'GAMA',2,ITER,1,1,1
  WRITE(2,'(A)') 'GAMA'
  
  DO INOD=1,NNOD
     WRITE(2,'(I6,1E16.8)') INOD,GAMM(INOD)       
  END DO

  IF (MOVIE.EQ.0)THEN
     CLOSE(2)
  END IF
  
  RETURN
END SUBROUTINE PRINTFLAVIA

!Calcular fuerzas para cada uno de los sets
SUBROUTINE FORCES(NSETS,NSET_NUMB,IELEM_SETS,ISET,NNOD,IFF &
     ,X,Y,P &
     ,FX,FY,RM)
  
  IMPLICIT REAL(8) (A-H,O-Z)
  
  INTEGER IELEM_SETS(10),ISET(2,10,IFF)
  
  REAL(8) X(NNOD),Y(NNOD)
  REAL(8) P(NNOD)
  REAL(8) FX(10),FY(10),RM(10)
  
  DO I=1,NSET_NUMB
     FX(I)=0.D0
     FY(I)=0.D0
     RM(I)=0.D0
  END DO
  
  X_ROT=0.D0
  Y_ROT=0.D0
  
  DO ISET_NUMB=1,NSET_NUMB
     DO II=1,IELEM_SETS(ISET_NUMB)
        N1=ISET(1,ISET_NUMB,II)
        N2=ISET(2,ISET_NUMB,II)
        
        D_PRESS=(P(N1)+P(N2))/2.D0
        RLX=X(N1)-X(N2)
        RLY=Y(N2)-Y(N1)
        
        DFX=D_PRESS*RLY
        DFY=D_PRESS*RLX
        
        FX(ISET_NUMB)=FX(ISET_NUMB)+DFX
        FY(ISET_NUMB)=FY(ISET_NUMB)+DFY
        
        XC=(X(N1)+X(N2))/2.D0
        YC=(Y(N2)+Y(N1))/2.D0
        
        RM(ISET_NUMB)=RM(ISET_NUMB)-DFY*(XC-X_ROT)+DFX*(YC-Y_ROT)
        
     END DO
  END DO
  RETURN
END SUBROUTINE FORCES
      
!Calcular fuerzas viscosas
      SUBROUTINE FORCE_VISC(NELEM,NNOD,IFF &
          ,IELEM_SETS,ISET,N,NSET_NUMB,UINF,VINF,RHOINF,TINF &
          ,X,Y,P,T,VEL_X,VEL_Y,DNX,DNY,RMU,RHO &
          ,F_VX,F_VY)
      

      IMPLICIT REAL(8) (A-H,O-Z)

      INTEGER IELEM_SETS(10),ISET(3,10,IFF),N(3,NELEM)

      REAL(8) X(NNOD),Y(NNOD),P(NNOD),VEL_X(NNOD),VEL_Y(NNOD),RHO(NNOD)
      REAL(8) DNX(3,NELEM),DNY(3,NELEM),T(NNOD)
      REAL(8) F_VX(10),F_VY(10)
      
      F_VX=0.D0 ; F_VY=0.D0
      OPEN(1,FILE='SKIN.DAT',STATUS='UNKNOWN')
      DO ISET_NUMB=1,NSET_NUMB

         DO II=1,IELEM_SETS(ISET_NUMB)
            
            NN1=ISET(1,ISET_NUMB,II)
            NN2=ISET(2,ISET_NUMB,II)
            IELEM=ISET(3,ISET_NUMB,II)

            TEMP=(T(NN1)+T(NN2))/2.D0
	    FMU= 1.716d-5*162.6/(TEMP-110.55)*(TEMP/273.15)**.75D0 !SUTHERLAND
            RLY=-(X(NN2)-X(NN1))
            RLX=Y(NN2)-Y(NN1)
            RMOD=DSQRT(RLX*RLX+RLY*RLY)

            DX=-RLY
            DY=RLX

            RLX=RLX/RMOD
            RLY=RLY/RMOD
            
            DUX=0.D0 ; DUY=0.D0 ; DVX=0.D0 ; DVY=0.D0
            PRESS=0.D0
            DO JJ=1,3
               NN=N(JJ,IELEM)
!CCCC  ----> DERIVADA DE VEL_X               
               DUX=DUX+DNX(JJ,IELEM)*VEL_X(NN)
               DUY=DUY+DNY(JJ,IELEM)*VEL_X(NN)
!CCCC  ----> DERIVADA DE VEL_Y 
               DVX=DVX+DNX(JJ,IELEM)*VEL_Y(NN)
               DVY=DVY+DNY(JJ,IELEM)*VEL_Y(NN)
               PRESS=PRESS+P(NN)
            END DO

            PRESS=PRESS/3.D0
!CCCC  ----> TENSOR 2D
            TXX=-FMU*2.D0/3.D0*(DUX+DVY)+2.D0*FMU*DUX
            TXY=FMU*(DUY+DVX)
            
            TYX=TXY
            TYY=-FMU*2.D0/3.D0*(DUX+DVY)+2.D0*FMU*DVY

!CCCC  ----> CALCULO PARA SACAR EL SKIN FRICTION
            TTX=(-PRESS+TXX)*RLX+TXY*RLY
            TTY=TYX*RLX+(-PRESS+TYY)*RLY
            TMOD=-TTX*RLY+TTY*RLX
            UU=UINF*UINF+VINF*VINF
            SKIN=TMOD/(.5d0*RHOINF*UU)

            F_VX(ISET_NUMB)=F_VX(ISET_NUMB)+((-PRESS+TXX)*RLX+TXY*RLY)*DX
            F_VY(ISET_NUMB)=F_VY(ISET_NUMB)+(TYX*RLX+(-PRESS+TYY)*RLY)*DY
            WRITE(1,*)SKIN,(X(NN2)+X(NN1))/2.D0
         END DO
      END DO
      CLOSE(1)
      RETURN
      END

!Imprimir resultados para Restart
      SUBROUTINE PRINTREST(ITER,NNOD,U,T,GAMM,TIME,FILE,ILONG)

      IMPLICIT REAL(8) (A-H,O-Z)
      
      REAL(8) U(4,NNOD),T(NNOD),GAMM(NNOD)

      CHARACTER FILE*80
      
      OPEN(1,FILE=FILE(1:ILONG)//'.RST',FORM='UNFORMATTED',STATUS='UNKNOWN')

      WRITE(1) ITER,TIME
      
      DO INOD=1,NNOD
         WRITE(1) (U(J,INOD),J=1,4),T(INOD),GAMM(INOD)
      END DO
      
      CLOSE(1) 
      
      RETURN
      END

!Verificar errores de lectura
      SUBROUTINE READ_ERROR(IERROR,IEL) 

      IMPLICIT REAL(8) (A-H,O-Z)

      CHARACTER(10) ERTYPE(5) 
      
      INTEGER IERROR,IEL 
      
      DATA ERTYPE /'NODOS','ELEMENTOS','NO_TRAC','FIX_VEL','NORM_VEL'/ 
      
      IF (IERROR.NE.0) THEN 
         WRITE(*,'(A)') 'ERROR EN LA LECTURA DE ',ERTYPE(IERROR) 
         STOP  
      END IF 
      
      RETURN 
      END  

!Calcular terminos fuente
SUBROUTINE FUENTE(NNOD,NELEM,N,AREA,DNX,DNY,W_X,W_Y,DTL &
     ,U,RHS,FGX,FGY,QH)
  
  IMPLICIT REAL(8) (A-H,O-Z)
  
  INTEGER N(3,NELEM)
  
  REAL(8) U(4,NNOD),DNX(3,NELEM),DNY(3,NELEM)
  REAL(8) AREA(NELEM),RHS(4,NNOD)
  REAL(8) W_X(NNOD),W_Y(NNOD),DTL(NELEM)
  REAL(8) ALF(3),BET(3)
  
  DATA ALF/.5D0,.5D0,0.D0/
  DATA BET/0.D0,.5D0,.5D0/
  
  NGAUSS=3    !PTOS DE GAUSS DONDE VOY A INTERGRAR
  
  DO IELEM=1,NELEM
     
     N1=N(1,IELEM) ; N2=N(2,IELEM) ; N3=N(3,IELEM)
 
     RNX1=DNX(1,IELEM) ; RNX2=DNX(2,IELEM) ; RNX3=DNX(3,IELEM)
     RNY1=DNY(1,IELEM) ; RNY2=DNY(2,IELEM) ; RNY3=DNY(3,IELEM)
     
     !CCCC  ----> DERIVADA DE RHO           
     DRX= RNX1*U(1,N1)+RNX2*U(1,N2)+RNX3*U(1,N3)
     DRY= RNY1*U(1,N1)+RNY2*U(1,N2)+RNY3*U(1,N3)
     !CCCC  ----> DERIVADA DE RHO.VEL_X           
     DRUX= RNX1*U(2,N1)+RNX2*U(2,N2)+RNX3*U(2,N3)
     DRUY= RNY1*U(2,N1)+RNY2*U(2,N2)+RNY3*U(2,N3)
     !CCCC  ----> DERIVADA DE RHO.VEL_Y           
     DRVX= RNX1*U(3,N1)+RNX2*U(3,N2)+RNX3*U(3,N3)
     DRVY= RNY1*U(3,N1)+RNY2*U(3,N2)+RNY3*U(3,N3)
     !CCCC  ----> DERIVADA DE RHO.ET           
     DREX= RNX1*U(4,N1)+RNX2*U(4,N2)+RNX3*U(4,N3)
     DREY= RNY1*U(4,N1)+RNY2*U(4,N2)+RNY3*U(4,N3)
     !CCCC  ----> DERIVADA DE LA VELOCIDAD DE LA MALLA

     AR=AREA(IELEM)*DTL(IELEM)/3.D0
     
     DO J=1,NGAUSS
        
        RN1=1.D0-ALF(J)-BET(J)
        RN2=ALF(J)
        RN3=BET(J)
        WX=RN1*W_X(N1)+RN2*W_X(N2)+RN3*W_X(N3)
        WY=RN1*W_Y(N1)+RN2*W_Y(N2)+RN3*W_Y(N3)

        !CCCC  ----> ENSAMBLE DEL RIGHT HAND SIDE DEL NODO N1
        RHS(1,N1)=RHS(1,N1)-RN1*(DRX *WX+DRY *WY)*AR   
        RHS(2,N1)=RHS(2,N1)-RN1*(DRUX*WX+DRUY*WY)*AR 
        RHS(3,N1)=RHS(3,N1)-RN1*(DRVX*WX+DRVY*WY)*AR 
        RHS(4,N1)=RHS(4,N1)-RN1*(DREX*WX+DREY*WY)*AR   
        
        !CCCC  ----> ENSAMBLE DEL RIGHT HAND SIDE DEL
        RHS(1,N2)=RHS(1,N2)-RN2*(DRX *WX+DRY *WY)*AR   
        RHS(2,N2)=RHS(2,N2)-RN2*(DRUX*WX+DRUY*WY)*AR 
        RHS(3,N2)=RHS(3,N2)-RN2*(DRVX*WX+DRVY*WY)*AR 
        RHS(4,N2)=RHS(4,N2)-RN2*(DREX*WX+DREY*WY)*AR   
        
        !CCCC  ----> ENSAMBLE DEL RIGHT HAND SIDE DEL
        RHS(1,N3)=RHS(1,N3)-RN3*(DRX *WX+DRY *WY)*AR   
        RHS(2,N3)=RHS(2,N3)-RN3*(DRUX*WX+DRUY*WY)*AR 
        RHS(3,N3)=RHS(3,N3)-RN3*(DRVX*WX+DRVY*WY)*AR 
        RHS(4,N3)=RHS(4,N3)-RN3*(DREX*WX+DREY*WY)*AR   
     END DO
  END DO
  
  RETURN
END SUBROUTINE FUENTE

!Obtener cantidad de letras del nombre del archivo de entrada
      FUNCTION LONG_FILE(FILE)

      IMPLICIT REAL(8) (A-H,O-Z)

      CHARACTER FILE*80

      INTEGER LONG_FILE
      
      DO I=1,80
        IF (FILE(I:I).EQ.' ') THEN
          LONG_FILE=I-1
          RETURN
        END IF
      END DO

      WRITE(10,'(A,A)')'ERROR.... ( hay un error con el nombre del file)'
      STOP
      
      END

!Calcular nodos vecidos y ensamblar laplaciano
 SUBROUTINE LAPLACE(NNOD,NELEM,N,AR,DNX,DNY,S,NN1,NN2,IAUX &
      ,NPOS,IND,INDEL,ADIAG,RAUX)
      
   IMPLICIT REAL(8) (A-H,O-Z)
   
   INTEGER N(3,NELEM),NN1(NNOD*15),NN2(NNOD*15) 
   INTEGER IND(NNOD),INDEL(10,NNOD),IAUX(50) 
   
   REAL(8) S(NNOD*15)
   REAL(8) RAUX(NNOD),ADIAG(NNOD)
   REAL(8) DNX(3,NELEM),DNY(3,NELEM),AR(NELEM)
   
   IPOS=0
   DO INOD=1,NNOD 
      IND(INOD)=0 
   END DO
   
   DO IELEM=1,NELEM 
      DO J=1,3
         N1=N(J,IELEM) 
         IND(N1)=IND(N1)+1 
         INDEL(IND(N1),N1)=IELEM 
      END DO
   END DO
   
   DO INOD=1,NNOD 
      
      NCON=1
      IAUX(1)=INOD
      DO INDICE=1,IND(INOD)
         IELEM=INDEL(INDICE,INOD)
         DO J=1,3
            DO ICON=1,NCON
               IF (IAUX(ICON).EQ.N(J,IELEM)) THEN
                  GOTO 100
               END IF
            END DO
            NCON=NCON+1
            IAUX(NCON)=N(J,IELEM)
100         CONTINUE
         END DO
      END DO
      
      DO INDICE=1,IND(INOD) 
         IELEM=INDEL(INDICE,INOD)
         AREA=AR(IELEM)
         DO I=1,3
            N1=N(I,IELEM) 
            
            IF (N1.EQ.INOD) THEN 
               DO J=1,3
                  N2=N(J,IELEM) 
                  RAUX(N2)=RAUX(N2)+(DNX(I,IELEM)*DNX(J,IELEM)+DNY(I,IELEM)*DNY(J,IELEM))/AREA**2.                    
               END DO
            END IF
         END DO
      END DO
      
      ADIAG(INOD)=RAUX(INOD)        
      DO ICON=1,NCON 
         IPOS=IPOS+1
         S(IPOS)=RAUX(IAUX(ICON))
         RAUX(IAUX(ICON))=0.D0
         NN1(IPOS)=INOD
         NN2(IPOS)=IAUX(ICON)	      
      END DO
   END DO
   
   NPOS=IPOS 
   
   RETURN
 END SUBROUTINE LAPLACE

!Pressure switch
SUBROUTINE NEWPRES(NNOD,NPOS,NN1,NN2,P,PS)

  IMPLICIT REAL(8) (A-H,O-Z)
  INTEGER(4) NN1(NPOS),NN2(NPOS)
  REAL(8)S1(NNOD),S2(NNOD),PS(NNOD),P(NNOD)

  S1=0.D0
  S2=0.D0
  
  DO IPOS=1,NPOS
     N1=NN1(IPOS)
     N2=NN2(IPOS)
     AUX=P(N1)-P(N2)
     S1(N1)=S1(N1)+AUX
     S2(N1)=S2(N1)+DABS(AUX)
  END DO
  
  DO I=1,NNOD
     PS(I)=DABS(S1(I))/(S2(I))
     IF(S2(I).LT.5.D-2)PS(I)=0.D0
     IF(PS(I).LT..2D0)PS(I)=0.D0 !MODIFICAR SI ES NECESARIO
  END DO 
  RETURN
END SUBROUTINE NEWPRES

!Calcular tamano optimo para refinamiento adaptativo
SUBROUTINE NEW_SIZE(NNOD,NELEM,N,DNX,DNY,AREA,M,HH,U &
     ,AUX)

  USE DATOS_REFINAMIENTO
  IMPLICIT REAL*8 (A-H,O-Z)
  
  INTEGER N(3,NELEM)
  REAL(8) DNX(3,NELEM),DNY(3,NELEM),AREA(NELEM),HH(NELEM),U(4,NNOD)
  REAL(8) HH_NEW(NELEM)
  REAL(8) SXX(NNOD),SXY(NNOD),SYY(NNOD),AUX(NNOD),M(NNOD)
  
  SXX=0.D0;  SYY=0.D0;  SXY=0.D0;   AUX=0.D0

  DO IELEM=1,NELEM
     N1=N(1,IELEM) ; N2=N(2,IELEM) ; N3=N(3,IELEM)
     
     RNX1=DNX(1,IELEM) ; RNX2=DNX(2,IELEM) ; RNX3=DNX(3,IELEM)
     RNY1=DNY(1,IELEM) ; RNY2=DNY(2,IELEM) ; RNY3=DNY(3,IELEM)
     
     !CCCC  ----> DERIVADA DE RHO.VEL_X           
     DUX= RNX1*U(2,N1)+RNX2*U(2,N2)+RNX3*U(2,N3)
     DUY= RNY1*U(2,N1)+RNY2*U(2,N2)+RNY3*U(2,N3)
     !CCCC  ----> DERIVADA DE RHO.VEL_Y           
     DVX= RNX1*U(3,N1)+RNX2*U(3,N2)+RNX3*U(3,N3)
     DVY= RNY1*U(3,N1)+RNY2*U(3,N2)+RNY3*U(3,N3)
     
     SIGMXX=2*DUX;   SIGMYY=2*DVY;    SIGMXY=DUY+DVX         
     
     DO I=1,3
        N1=N(I,IELEM)
        SXX(N1)=SXX(N1)+SIGMXX*AREA(IELEM)
        SYY(N1)=SYY(N1)+SIGMYY*AREA(IELEM)
        SXY(N1)=SXY(N1)+SIGMXY*AREA(IELEM)
        AUX(N1)=AUX(N1)+AREA(IELEM)
     END DO
  END DO
  
  DO INOD=1,NNOD
     SXX(INOD)=SXX(INOD)/AUX(INOD)
     SYY(INOD)=SYY(INOD)/AUX(INOD)
     SXY(INOD)=SXY(INOD)/AUX(INOD)
  END DO
  
  
  ESIG=0.D0;   SIG=0.D0     
  DO IELEM=1,NELEM
     N1=N(1,IELEM) ; N2=N(2,IELEM) ; N3=N(3,IELEM)
     
     RNX1=DNX(1,IELEM) ; RNX2=DNX(2,IELEM) ; RNX3=DNX(3,IELEM)
     RNY1=DNY(1,IELEM) ; RNY2=DNY(2,IELEM) ; RNY3=DNY(3,IELEM)
     
     !CCCC  ----> DERIVADA DE RHO.VEL_X           
     DUX= RNX1*U(2,N1)+RNX2*U(2,N2)+RNX3*U(2,N3)
     DUY= RNY1*U(2,N1)+RNY2*U(2,N2)+RNY3*U(2,N3)
     !CCCC  ----> DERIVADA DE RHO.VEL_Y           
     DVX= RNX1*U(3,N1)+RNX2*U(3,N2)+RNX3*U(3,N3)
     DVY= RNY1*U(3,N1)+RNY2*U(3,N2)+RNY3*U(3,N3)
     
     SIGMXX=2*DUX;   SIGMYY=2*DVY;    SIGMXY=DUY+DVX         
     
     ESIGXX=SIGMXX-(SXX(N1)+SXX(N2)+SXX(N3))/3.D0
     SIGXX=(SXX(N1)+SXX(N2)+SXX(N3))/3.D0
     
     ESIGYY=SIGMYY-(SYY(N1)+SYY(N2)+SYY(N3))/3.D0
     SIGYY=(SYY(N1)+SYY(N2)+SYY(N3))/3.D0
     
     ESIGXY=SIGMXY-(SXY(N1)+SXY(N2)+SXY(N3))/3.D0
     SIGXY=(SXY(N1)+SXY(N2)+SXY(N3))/3.D0
         
     ESIG=ESIG+(ESIGXX*ESIGXX+ESIGYY*ESIGYY+ESIGXY*ESIGXY)*AREA(IELEM)
     SIG=SIG+(SIGXX*SIGXX+SIGYY*SIGYY+SIGXY*SIGXY)*AREA(IELEM)
  END DO
  
  ERR_AVERAGE=ETA_REFIN*DSQRT((SIG+ESIG)/DFLOAT(NELEM))

  DO IELEM=1,NELEM
     
     N1=N(1,IELEM) ; N2=N(2,IELEM) ; N3=N(3,IELEM)
     
     RNX1=DNX(1,IELEM) ; RNX2=DNX(2,IELEM) ; RNX3=DNX(3,IELEM)
     RNY1=DNY(1,IELEM) ; RNY2=DNY(2,IELEM) ; RNY3=DNY(3,IELEM)
     
     !CCCC  ----> DERIVADA DE RHO.VEL_X           
     DUX= RNX1*U(2,N1)+RNX2*U(2,N2)+RNX3*U(2,N3)
     DUY= RNY1*U(2,N1)+RNY2*U(2,N2)+RNY3*U(2,N3)
     !CCCC  ----> DERIVADA DE RHO.VEL_Y           
     DVX= RNX1*U(3,N1)+RNX2*U(3,N2)+RNX3*U(3,N3)
     DVY= RNY1*U(3,N1)+RNY2*U(3,N2)+RNY3*U(3,N3)
     
     SIGMXX=2*DUX;   SIGMYY=2*DVY;    SIGMXY=DUY+DVX         
     
     ESIGXX=SIGMXX-(SXX(N1)+SXX(N2)+SXX(N3))/3.D0
     
     ESIGYY=SIGMYY-(SYY(N1)+SYY(N2)+SYY(N3))/3.D0
         
     ESIGXY=SIGMXY-(SXY(N1)+SXY(N2)+SXY(N3))/3.D0
         
     ESIG=(ESIGXX*ESIGXX+ESIGYY*ESIGYY+ESIGXY*ESIGXY)*AREA(IELEM)
     
     EPS_I=ESIG/ERR_AVERAGE
     HH_NEW(IELEM)=HH(IELEM)/EPS_I
     
  END DO  

  AUX=0.D0
  DO I=1,NELEM
     AR=AREA(I)
     DO J=1,3
        N1=N(J,I)
        AUX(N1)=AUX(N1)+HH_NEW(I)*AR
     END DO
  END DO
  DO I=1,NNOD
     AUX(I)=AUX(I)/M(I)
     IF(AUX(I).GT.HHMAX_REFIN) AUX(I)=HHMAX_REFIN
     IF(AUX(I).LT.HHMIN_REFIN) AUX(I)=HHMIN_REFIN
  END DO
  
  RETURN
END SUBROUTINE NEW_SIZE

!Estimador de error
SUBROUTINE ESTIMADOR_ERR(NNOD,NELEM,N,DNX,DNY,HH,M,AREA,BETA,VAR,AUX)
  USE DATOS_REFINAMIENTO
  IMPLICIT REAL(8) (A-H,O-Z)
  INTEGER N(3,NELEM)
  REAL(8)DNX(3,NELEM),DNY(3,NELEM),HH(NELEM),VAR(NNOD)
  REAL(8) TITA(NELEM),AUX(NNOD),AREA(NNOD),M(NNOD)
  TITA=0.D0 ; TITA_PROM=0.D0
  DO IELEM=1,NELEM
     N1=N(1,IELEM) ; N2=N(2,IELEM) ; N3=N(3,IELEM)
     
     RNX1=DNX(1,IELEM) ; RNX2=DNX(2,IELEM) ; RNX3=DNX(3,IELEM)
     RNY1=DNY(1,IELEM) ; RNY2=DNY(2,IELEM) ; RNY3=DNY(3,IELEM)
     
     !CCCC  ----> DERIVADA DE LA VARIABLE           
     DUX= RNX1*VAR(N1)+RNX2*VAR(N2)+RNX3*VAR(N3)
     DUY= RNY1*VAR(N1)+RNY2*VAR(N2)+RNY3*VAR(N3)

     TITA(IELEM)=DSQRT((DABS(DUX)+DABS(DUY))*.5D0)*HH(IELEM)
     TITA_PROM=TITA_PROM+TITA(IELEM)
  END DO
  TITA_PROM=TITA_PROM/NELEM

  TITA_SD=0.D0
  DO I=1,NELEM
     TITA_SD=TITA_SD+(TITA(I)-TITA_PROM)**2
  END DO
  TITA_SD=DSQRT(TITA_SD/NELEM)

  !CCCC----> SUAVIZO EN A LOS NODOS
  AUX=0.D0
  DO I=1,NELEM
     AR=AREA(I)
     DO J=1,3
        N1=N(J,I)
        AUX(N1)=AUX(N1)+TITA(I)*AR
     END DO
  END DO
  DO I=1,NNOD
     AUX(I)=AUX(I)/M(I)
  END DO

  !CCCC----> APLICO FILTRO PASA ALTO
  !BETA=3.5D0 ! MODIFICO A OJO
  COMPARADOR=TITA_PROM+BETA*TITA_SD
  DO I=1,NNOD
     IF(AUX(I).GT.COMPARADOR)THEN
        AUX(I)=HHMIN_REFIN
     ELSE
        AUX(I)=HHMAX_REFIN
     END IF
  END DO

RETURN
END SUBROUTINE ESTIMADOR_ERR

!Solver (gradientes conjugados)    
SUBROUTINE GRADCONJ2(S,PRESS,B,NN1,NN2,NNOD,NPOS &
     ,IFIXPRES,RFIXP_VALUE,IFIXP,RES,ADIAG &
     ,PK,APK,Z) 
  
  IMPLICIT REAL(8) (A-H,O-Z)
  
  INTEGER NN1(NPOS),NN2(NPOS),IFIXPRES(IFIXP) 
  REAL(8) RES(NNOD),PRESS(NNOD),S(NPOS) 
  REAL(8) B(NNOD),ADIAG(NNOD),RFIXP_VALUE(IFIXP)  
  !CCCC ---> AUXILIARES
  REAL(8)PK(NNOD),APK(NNOD),Z(NNOD)
  
  CONJERR=1.D-10

  DO IN=1,IFIXP 
     PRESS(IFIXPRES(IN))=RFIXP_VALUE(IN) 
  END DO
  
  CALL RESIDUO2(RES,S,PRESS,NN1,NN2,NPOS,NNOD,IFIXPRES,IFIXP) 
  
  DO INOD=1,NNOD  
     RES(INOD)=B(INOD)-RES(INOD)  
  END DO
 
  DO IN=1,IFIXP 
     RES(IFIXPRES(IN))=0.D0 
  END DO
  
  RR_12=PKAPK2(RES,RES,NNOD)
  K=0 
  DO WHILE (DABS(RR_12).GT.CONJERR.AND.K.LT.1000) 
     
     K=K+1 
     DO INOD=1,NNOD 
        Z(INOD)=RES(INOD)/ADIAG(INOD) 
     END DO
     
     RR_12=PKAPK2(RES,Z,NNOD) 
     
     IF (K.EQ.1) THEN 
        
            DO IK=1,NNOD 
               PK(IK)=Z(IK) 
            END DO 
            
         ELSE 
            
            BET=RR_12/RR_22 
            
            DO IK=1,NNOD 
               PK(IK)=Z(IK)+BET*PK(IK) 
            END DO 
            
         END IF 
         
         CALL RESIDUO2(APK,S,PK,NN1,NN2,NPOS,NNOD,IFIXPRES,IFIXP) 
         ALF=RR_12/PKAPK2(APK,PK,NNOD) 
         
         DO IK=1,NNOD 
            
            PRESS(IK)=PRESS(IK)+ALF*PK(IK) 
            RES(IK)=RES(IK)-ALF*APK(IK) 
            
         END DO 
         
         RR_22=RR_12 
   
      END DO 
 
      !WRITE(*,'(A,I5)')'  ITERACIONES DE GC....',K
      !WRITE(*,'(A,D14.6)')'      RESIDUO FINAL....',RR_22
      
      RETURN 
    END SUBROUTINE GRADCONJ2
      
SUBROUTINE RESIDUO2(RES,S,PRESS,NN1,NN2,NPOS,NNOD &
     ,IFIXPRES,IFIXP) 
      
  IMPLICIT REAL(8) (A-H,O-Z)
  
  INTEGER NN1(NPOS),NN2(NPOS),IFIXPRES(IFIXP) 
  REAL(8) S(NPOS),RES(NNOD),PRESS(NNOD) 
  
  RES=0.D0 
    
  DO INOD=1,NPOS 
     RES(NN1(INOD))=RES(NN1(INOD))+S(INOD)*PRESS(NN2(INOD))      
  END DO
  
  DO IN=1,IFIXP 
     RES(IFIXPRES(IN))=PRESS(IFIXPRES(IN))*1.D30 
  END DO
  
  RETURN 
END SUBROUTINE RESIDUO2

FUNCTION FRR_12(RES,NNOD) 
  
  IMPLICIT REAL(8) (A-H,O-Z)
  
  REAL(8) RES(NNOD)
  
  FRR_12=0.D0 
  DO INOD=1,NNOD 
     FRR_12=FRR_12+RES(INOD)*RES(INOD) 
  END DO
  
  RETURN 
END FUNCTION FRR_12

FUNCTION PKAPK2(APK,PK,NNOD) 
  
  IMPLICIT REAL(8) (A-H,O-Z)
  REAL(8) APK(NNOD),PK(NNOD) 
  
  PKAPK2=0.D0 
  DO INOD=1,NNOD 
     PKAPK2=PKAPK2+PK(INOD)*APK(INOD) 
  END DO
  
  RETURN 
END FUNCTION PKAPK2

!Smoothing de elementos
SUBROUTINE SMOOTH_MESH(NNOD,NELEM,N,X,Y &
     ,SMOOTH_FIX,SMOOTH_SIM)
  
  IMPLICIT REAL(8) (A-H,O-Z)
  
  INTEGER N(3,NELEM),SMOOTH_SIM(2,NNOD)
  INTEGER NELINOD(NNOD),NE(20,NNOD)
  
  REAL(8) X(NNOD),Y(NNOD)
  REAL(8) MU_0,MU_X,MU_Y,MU_Z,RMU1,RMUMIN,GG,GGI,GAM1,GAMMIN
  REAL(8) MU(100),GX(100),GY(100),XX(4),YY(4)
  
  LOGICAL SMOOTH_FIX(NNOD)
  
  write(*,*) 'suavizado inicial  1'
  
  DO INOD=1,NNOD
     NELINOD(INOD)=0  
  END DO
  
  DO IELEM=1,NELEM
     DO IN=1,3
        INOD=N(IN,IELEM)
        NELINOD(INOD)=NELINOD(INOD)+1
        NE(NELINOD(INOD),INOD)=IELEM
     END DO
  END DO
  
  DELTA=1.D-5
  !CCCC   COMIENZO DE ITERACIONES DE OPTIMIZACION
  
  DO IT=1,1000
     
     
     RMM=1.D10
     ERR=0.D0; xmax=0
     NUMEL=0
     DO INOD=1,NNOD
        
        IF (SMOOTH_FIX(INOD)) THEN               
           
           RMUMIN=1.D10
           DO IEL=1,NELINOD(INOD)
              IELEM=NE(IEL,INOD)
              
              DO IN=1,3
                 N1=N(IN,IELEM)
                 XX(IN)=X(N1)
                 YY(IN)=Y(N1)
                 IF (N1.EQ.INOD) IN_P=IN
              END DO
              
              MU_0=RMU1(XX,YY)
              
              IF (RMM.GT.MU_0) THEN
                 RMM=MU_0
                 IELMIN=IELEM
              END IF
              
              MU(IEL)=MU_0
              IF (MU_0.LT.RMUMIN) THEN
                 RMUMIN=MU_0
                 IELM=IEL
              END IF
              
              XX(IN_P)=X(INOD)+DELTA
              MU_X=RMU1(XX,YY)
              XX(IN_P)=X(INOD)
              GX(IEL)=(MU_X-MU_0)/DELTA
              
              YY(IN_P)=Y(INOD)+DELTA
              MU_Y=RMU1(XX,YY)
              YY(IN_P)=Y(INOD)
              GY(IEL)=(MU_Y-MU_0)/DELTA
           END DO
           
           IF (RMUMIN.GT.0.8D0) CYCLE
           NUMEL=NUMEL+1
           
           GMX=GX(IELM)
           GMY=GY(IELM)
           
           GAM1=0.D0
           GAMMIN=1.D10
           DO IEL=1,NELINOD(INOD)
              GG=GMX*GMX+GMY*GMY
              IF (IEL.NE.IELM) THEN
                 GGI=GX(IEL)*GMX+GY(IEL)*GMY
                 IF (GGI.LT.0) THEN
                    GAM1=(MU(IEL)-RMUMIN)/(GG-GGI)
                    IF (GAMMIN.GT.GAM1.AND.GAM1.NE.0.) GAMMIN=GAM1
                 END IF
              END IF
           END DO
           
           GAMMIN=GAMMIN/8
           IF (GAMMIN.GT.100) GAMMIN=1.D-6
           
           X1=X(INOD)+GMX*GAMMIN*SMOOTH_SIM(1,INOD)
           Y1=Y(INOD)+GMY*GAMMIN*SMOOTH_SIM(2,INOD)               
           
           ERR=ERR+(X(INOD)-X1)**2+(Y(INOD)-Y1)**2
           
           X(INOD)=X1
           Y(INOD)=Y1
           
        END IF
     END DO
     
     IF (NUMEL.EQ.0) EXIT
     
  END DO
  
  WRITE(*,*) 'ERROR FINAL:', ERR, '  PASOS:', IT-1,NUMEL
  
  RETURN 
END SUBROUTINE SMOOTH_MESH

FUNCTION RMU1(XX,YY)
  IMPLICIT REAL(8) (A-H,O-Z)
  REAL(8) XX(4),YY(4)
  
  AREA=XX(2)*YY(3)+XX(3)*YY(1)+XX(1)*YY(2) &
       -(XX(2)*YY(1)+XX(3)*YY(2)+XX(1)*YY(3))
  
  RL21=(XX(2)-XX(1))**2.D0+(YY(2)-YY(1))**2.D0
  RL32=(XX(3)-XX(2))**2.D0+(YY(3)-YY(2))**2.D0
  RL13=(XX(1)-XX(3))**2.D0+(YY(1)-YY(3))**2.D0
  
  RL=RL21+RL32+RL13
  
  RMU1=2.D0*DSQRT(3.D0)*AREA/RL
  
  RETURN
END FUNCTION RMU1

!Multiples formas de entrada de datos (A. Figueroa)
SUBROUTINE VARIABLES
  USE DATOS_ENTRADA

  IMPLICIT REAL(8)(A-H,O-Z)
  
  IF(TINF.EQ.0.D0)TINF=PINF/(FR*RHOINF)
  IF(PINF.EQ.0.)PINF=RHOINF*FR*TINF
  IF(RHOINF.EQ.0.D0) RHOINF=PINF/(FR*TINF)
  
  AUX=DSQRT(UINF**2.D0+VINF**2.D0)
  CINF=DSQRT(GAMA*FR*TINF)
  IF(AUX.EQ.0.D0) UINF=CINF*MACHINF

END SUBROUTINE VARIABLES
