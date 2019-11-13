# if defined (SEDIMENT)
! --------------------------------------------------
!    This is subroutine to update the sediment load
!    called by
!         MAIN
!    Update: 03/14/2016 Babak Tehranirad, University of Delaware
!    Update: 06/09/2017 Fengyan Shi, merge into 3.1beta
! --------------------------------------------------
SUBROUTINE eval_sed_load
	USE GLOBAL
    IMPLICIT NONE
	INTEGER::ISTEP,ISTAGE

    REAL(SP)::Dstar,viscosity
	H=Eta*Gamma3+Depth+Zb
    	viscosity=0.00000100
    	Dstar = D50*((Sdensity-1.)*grav/viscosity**2.)**(1./3.)
    	Tau_cr=(Sdensity-1.)*grav*D50*Shields_cr	

		k_s = 2.5*D50 
		DO J=Jbeg,Jend+1
    		DO I=Ibeg,Iend+1
    			H(I,J)=MAX(H(I,J),MinDepth)
    		ENDDO
    	ENDDO
    IF(Sed_Scheme(1:3)=='Upw')THEN                            ! SED LOAD IF

		DO J=Jbeg,Jend+1
    		DO I=Ibeg,Iend+1
				IF(MASK(I,J).GT.0)THEN

					u_c =SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J))
					u_c1=SQRT(U(I+1,J)*U(I+1,J)+V(I+1,J)*V(I+1,J))
					u_c2=SQRT(U(I-1,J)*U(I-1,J)+V(I-1,J)*V(I-1,J))
					u_c3=SQRT(U(I,J+1)*U(I,J+1)+V(I,J+1)*V(I,J+1))
					u_c4=SQRT(U(I,J-1)*U(I,J-1)+V(I,J-1)*V(I,J-1))					

					ustar_c = 0.4*u_c/(-1.+log(30.*MAX(H(I+1,J),MinDepth)/k_s))					
					ustar_c1= 0.4*u_c1/(-1.+log(30.*MAX(H(I+1,J),MinDepth)/k_s))					
					ustar_c2= 0.4*u_c2/(-1.+log(30.*MAX(H(I+1,J),MinDepth)/k_s))					
					ustar_c3= 0.4*u_c3/(-1.+log(30.*MAX(H(I+1,J),MinDepth)/k_s))					
					ustar_c4= 0.4*u_c4/(-1.+log(30.*MAX(H(I+1,J),MinDepth)/k_s))


					k1=5.93*(ustar_c1+ustar_c)*(MAX(H(I+1,J),MinDepth)+MAX(H(I,J),MinDepth))/4.0
					k2=5.93*(ustar_c2+ustar_c)*(MAX(H(I-1,J),MinDepth)+MAX(H(I,J),MinDepth))/4.0
					k3=5.93*(ustar_c3+ustar_c)*(MAX(H(I,J),MinDepth)+MAX(H(I,J+1),MinDepth))/4.0
					k4=5.93*(ustar_c4+ustar_c)*(MAX(H(I,J),MinDepth)+MAX(H(I,J-1),MinDepth))/4.0
					IF (U(I+1,J).GE.0.0) THEN
						F1=CH(I,J)*U(I+1,J)*(H(I+1,J)+H(I,J))/2.0 &
							-k1*(H(I+1,J)+H(I,J))*(CH(I+1,J)-CH(I,J))/2.0/DX
					ELSE
						F1=CH(I+1,J)*U(I+1,J)*(H(I+1,J)+H(I,J))/2.0 &
							-k1*(H(I+1,J)+H(I,J))*(CH(I+1,J)-CH(I,J))/2.0/DX
					ENDIF
					IF (U(I,J).GE.0.0) THEN
						F2=CH(I-1,J)*U(I,J)*(H(I-1,J)+H(I,J))/2.0 &
							-k2*(H(I-1,J)+H(I,J))*(CH(I,J)-CH(I-1,J))/2.0/DX
					ELSE
						F2=CH(I,J)*U(I,J)*(H(I-1,J)+H(I,J))/2.0 &
							-k2*(H(I-1,J)+H(I,J))*(CH(I,J)-CH(I-1,J))/2.0/DX
					ENDIF

					IF (V(I,J+1).GE.0.0) THEN
						F3=CH(I,J)*V(I,J+1)*(H(I,J)+H(I,J+1))/2.0 &
							-k3*(H(I,J+1)+H(I,J))*(CH(I,J+1)-CH(I,J))/2.0/DY
					ELSE
						F3=CH(I,J+1)*V(I,J+1)*(H(I,J)+H(I,J+1))/2.0 &
							-k3*(H(I,J+1)+H(I,J))*(CH(I,J+1)-CH(I,J))/2.0/DY
					ENDIF
					IF (V(I,J).GE.0.0) THEN
						F4=CH(I,J-1)*V(I,J)*(H(I,J)+H(I,J-1))/2.0 &
							-k4*(H(I,J+1)+H(I,J))*(CH(I,J)-CH(I,J-1))/2.0/DY
					ELSE
						F4=CH(I,J)*V(I,J)*(H(I,J)+H(I,J-1))/2.0 &
							-k4*(H(I,J-1)+H(I,J))*(CH(I,J)-CH(I,J-1))/2.0/DY
					ENDIF
! BC (Shoreline)
						IF (MASK(I+1,J).LT.1) THEN
							F1=0.0
						ENDIF
						IF (MASK(I-1,J).LT.1) THEN
							F2=0.0
						ENDIF
						IF (MASK(I,J-1).LT.1) THEN
							F4=0.0
						ENDIF
						IF (MASK(I,J+1).LT.1) THEN
							F3=0.0
						ENDIF
					
! BC (Wall Condition) X-direction
					IF (npx.EQ.0) THEN
						IF (I.EQ.Ibeg) F2=0.0
					ELSEIF (npx.EQ.px-1) THEN
						IF (I.EQ.Iend+1) F1=0.0
					ENDIF
				
! BC (Wall Condition) y-direction
					IF (py.EQ.1) THEN
						IF (J.EQ.Jbeg) F4=0.0
						IF (J.EQ.Jend) F3=0.0
					ELSE
					IF (npy.EQ.0) THEN
						IF (J.EQ.Jbeg) F4=0.0
					ELSEIF (npy.EQ.py-1) THEN    
						IF (J.EQ.Jend) F3=0.0
					ENDIF
					ENDIF
					


					Delta_c(I,J)=((Pickup(I,J)-D(I,J))*DX*DY-(F1-F2)*DY-(F3-F4)*DX) &
								*DT/DX/DY/H(I,J)      !MAX(H(I,J),MinDepth)
					
					CH(I,J)=CH(I,J)+Delta_c(I,J)
					
				ENDIF
			ENDDO
    	ENDDO
    	
    	
    	
		DO J=Jbeg,Jend+1
    		DO I=Ibeg,Iend+1
				IF(MASK(I,J).GT.0)THEN
					u_c=SQRT(((U(I,J)+U(I+1,J))/2.)**2.+((V(I,J)+V(I,J+1))/2.)**2.)
					tau_xy= 0.16/(1.+log(k_s/(30.*MAX(H(I+1,J),MinDepth))))**2*(u_c**2.)	
					
					IF (tau_xy.GT.Tau_cr) THEN
						c_b=0.015*(((tau_xy-Tau_cr)/Tau_cr)**1.50)*Dstar**(-0.3)
						reduction=MIN(1.0,0.65/c_b)
						c_a=reduction*c_b*D50/(0.01*H(I+1,J))
						Pickup(I,J)=MAX(0.0,c_a*WS)
					ELSE
						Pickup(I,J)=0.0
					ENDIF
				ELSE
					Pickup(I,J)=0.0	
				ENDIF
				IF(IN_Mask_s) THEN
					IF(Zb(I,J).GE.(Zs(I,J)-0.001)) Pickup(I,J)=0.0
				ENDIF
			ENDDO
    	ENDDO

! Calculate the Deposition Rate D(I,J)
    	DO J=Jbeg,Jend+1
    		DO I=Ibeg,Iend+1
				IF(MASK(I,J).GT.0)THEN
!  Cao(2004)
			  D(I,J)=MIN(2.0,(1.-n_porosity)/CH(I,J))*CH(I,J)*WS* &
			  				(1-MIN(2.0,(1.-n_porosity)/CH(I,J))*CH(I,J))**2.

			  ELSE
			  			D(I,J)=0.0
			  ENDIF
			ENDDO
    	ENDDO
    ELSE
    	CH0=CH
   	DO ISTAGE=1,3
	    	CALL sediment_flux
    		CALL sediment_source
    		CALL sediment_solver(ISTAGE)
		    CALL EXCHANGE
    	ENDDO
    	
	ENDIF   			! SED LOAD IF    

! Hard Bottom
!		DO J=Jbeg,Jend+1
 !   		DO I=Ibeg,Iend+1
  !  			IF (IN_Mask_s) THEN
   ! 				IF(Mask_s(I,J).LT.0.5) THEN
    !					IF (Zb(I,J).LT.0.02) Pickup(I,J)=0.0
	!				ENDIF
!				ENDIF
!			ENDDO
!		ENDDO

 	H=Eta*Gamma3+Depth+Zb	

    
END SUBROUTINE eval_sed_load

SUBROUTINE sediment_solver(ISTEP)
!-----------------------------------------------------------
!   Solves Sediment Transport Equation with Runge-Kutta time stepping 
!   Called by 
!      eval_sed_load
!    update: 03/10/2015, Babak Tehranirad
!    Update: 06/09/2017 Fengyan Shi, merge into 3.1beta
!-----------------------------------------------------------
	USE GLOBAL
	IMPLICIT NONE
	INTEGER,INTENT(IN)::ISTEP
	REAL(SP),PARAMETER::n_left=-1.0_SP,n_right=1.0_SP,n_bottom=-1.0_SP,n_top=1.0_SP
	REAL(SP)::F_left,F_right,F_bottom,F_top,SED_Source
	REAL(SP),DIMENSION(Ibeg:Iend,Jbeg:Jend)::RS
    
	DO J=Jbeg,Jend
		DO I=Ibeg,Iend
			F_left=Fs(I,J)
			F_right=Fs(I+1,J)
			F_bottom=Gs(I,J)
			F_top=Gs(I,J+1)
			Sed_Source=Pickup(I,J)-D(I,J)
			RS(I,J)=(-1.0_SP/DX*(F_right*n_right+F_left*n_left) &
                -1.0_SP/DY*(F_top*n_top+F_bottom*n_bottom) &
                +Sed_Source)/H(I,J)
            CH(I,J)=MAX(0.0_SP,ALPHA(ISTEP)*CH0(I,J)+BETA(ISTEP)*(CH(I,J)+DT*RS(I,J))) &
            		*FLOAT(MASK(I,J))
            IF (CH(I,J).GT.(.750_SP/13.0_SP)) THEN
!					Zb(I,J)=Zb(I,J)+(CH(I,J)-(1.0_SP/13.0_SP))*H(I,J)/(1.-n_porosity)
            		 CH(I,J)=(.750_SP/13.0_SP)
            ENDIF
        ENDDO
    ENDDO      
    
    
    
END SUBROUTINE sediment_solver

SUBROUTINE sediment_source
!-----------------------------------------------------------
!   Calculates Pickup and Deposition Rates 
!   Called by 
!      eval_sed_load
!    update: 03/10/2015, Babak Tehranirad
!    Update: 06/09/2017 Fengyan Shi, merge into 3.1beta
!-----------------------------------------------------------
    USE GLOBAL
    IMPLICIT NONE
    REAL(SP)::Dstar,viscosity

    viscosity=0.00000100_SP
    Dstar = D50*((Sdensity-1._SP)*grav/viscosity**2._SP)**(1._SP/3._SP)
    Tau_cr=(Sdensity-1._SP)*grav*D50*Shields_cr	
	k_s = 2.5_SP*D50     	



! Calculate the Deposition Rate D(I,J)
    DO J=Jbeg,Jend
    	DO I=Ibeg,Iend
!			IF(MASK(I,J).GT.0)THEN
!  Cao(2004)
				D(I,J)=MIN(2.0_SP,(1._SP-n_porosity)/CH(I,J))*CH(I,J)*WS* &
			  				(1._SP-MIN(2.0,(1._SP-n_porosity)/CH(I,J))*CH(I,J))**2._SP
!			ELSE
!				D(I,J)=0.0_SP
!			ENDIF
!			IF (CH(I,J).GE.(1.0_SP/13.0_SP)) THEN
!				D(I,J)=ZERO
!			ENDIF
		ENDDO
    ENDDO
    
     ! Calculate the Pickup Rate Pickup(I,J)   	
	DO J=Jbeg,Jend
    	DO I=Ibeg,Iend
			IF(MASK(I,J).GT.0)THEN
				u_c=SQRT(((U(I,J)+U(I+1,J))/2.)**2._SP+((V(I,J) &
								+V(I,J+1))/2._SP)**2._SP)
				tau_xy= 0.16_SP/(1._SP+log(k_s/(30._SP*MAX(H(I+1,J),MinDepth)))) &
							**2._SP*(u_c**2._SP)	
					
				IF (tau_xy.GT.Tau_cr) THEN
					c_b=0.015_SP*(((tau_xy-Tau_cr)/Tau_cr)**2.40_SP)*Dstar**(-0.6_SP)
					reduction=MIN(1.0_SP,0.65_SP/c_b)
					c_a=reduction*c_b*D50/(0.01_SP*H(I+1,J))
					Pickup(I,J)=MAX(0.0_SP,c_a*WS)
				ELSE
					Pickup(I,J)=0.0_SP
				ENDIF
			ELSE
				Pickup(I,J)=0.0_SP	
			ENDIF
			IF (CH(I,J).GE.(.50_SP/13.0_SP)) THEN
				Pickup(I,J)=0.0_SP
			ENDIF
			IF(H(I,J).LE.Mindepth*1.2_SP) Pickup(I,J)=0.0_SP
			IF(IN_Mask_s) THEN
					IF(Zb(I,J).GE.(Zs(I,J)-0.001)) Pickup(I,J)=0.0
			ENDIF
		ENDDO
    ENDDO

END SUBROUTINE sediment_source


SUBROUTINE sediment_flux
!-----------------------------------------------------------
!   Update bed elevation 
!   Called by 
!      eval_sed_load
!   update: 03/10/2015, Babak Tehranirad
!    Update: 06/09/2017 Fengyan Shi, merge into 3.1beta
!-----------------------------------------------------------
    USE GLOBAL
    IMPLICIT NONE

! calculate dc/dx and dc/dy for the diffusion term    
    CALL DERIVATIVE_X_High(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK,DX,CH,CHX)
	CALL DERIVATIVE_Y_High(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK,DY,CH,CHY)
	CHH=CH*H	
	
! construct in x-direction	

 !   CALL CONSTRUCT_X(Mloc,Mloc1,Nloc,DX,U,DelxU,UxL,UxR,Kappa)
 !   CALL DelxFun(DX,Mloc,Nloc,U,DelxU)

	CALL CONSTRUCT_HO_X(Mloc,Nloc,Mloc1,Ibeg,Iend,Jbeg,Jend,DX,MASK,U,UxL,UxR)
    CALL CONSTRUCT_HO_X(Mloc,Nloc,Mloc1,Ibeg,Iend,Jbeg,Jend,DX,MASK,V,VxL,VxR)
    CALL CONSTRUCT_HO_X(Mloc,Nloc,Mloc1,Ibeg,Iend,Jbeg,Jend,DX,MASK,CH,CHHxL,CHHxR)
    CALL CONSTRUCT_HO_X(Mloc,Nloc,Mloc1,Ibeg,Iend,Jbeg,Jend,DX,MASK,CHX,CHXxL,CHXxR)
    CALL CONSTRUCT_HO_X(Mloc,Nloc,Mloc1,Ibeg,Iend,Jbeg,Jend,DX,MASK,Eta,EtaRxL,EtaRxR)
     HxL=EtaRxL+Depthx+Zb
     HxR=EtaRxR+Depthx+Zb

	FsL(1:Mloc1,1:Nloc)=CHHxL(1:Mloc1,1:Nloc)*UxL(1:Mloc1,1:Nloc)*HxL(1:Mloc1,1:Nloc) &
						-HxL(1:Mloc1,1:Nloc)*CHXxL(1:Mloc1,1:Nloc)*5.93_SP* &
						HxL(1:Mloc1,1:Nloc)*SQRT(UxL(1:Mloc1,1:Nloc)*UxL(1:Mloc1,1:Nloc) &
						+VxL(1:Mloc1,1:Nloc)*VxL(1:Mloc1,1:Nloc))* &
						0.4_SP/(-1+log(30._SP*MAX(HxL(1:Mloc1,1:Nloc),MinDepth)/k_s))
    	
    FsR(1:Mloc1,1:Nloc)=CHHxR(1:Mloc1,1:Nloc)*UxR(1:Mloc1,1:Nloc)*HxR(1:Mloc1,1:Nloc) &
						-HxR(1:Mloc1,1:Nloc)*CHXxR(1:Mloc1,1:Nloc)*5.93_SP* &
						HxR(1:Mloc1,1:Nloc)*SQRT(UxR(1:Mloc1,1:Nloc)*UxR(1:Mloc1,1:Nloc) &
						+VxR(1:Mloc1,1:Nloc)*VxR(1:Mloc1,1:Nloc))* &
						0.4_SP/(-1+log(30._SP*MAX(HxR(1:Mloc1,1:Nloc),MinDepth)/k_s))
 
 ! construct in y-direction
    CALL CONSTRUCT_HO_Y(Mloc,Nloc,Nloc1,Ibeg,Iend,Jbeg,Jend,DY,MASK,U,UyL,UyR)
    CALL CONSTRUCT_HO_Y(Mloc,Nloc,Nloc1,Ibeg,Iend,Jbeg,Jend,DY,MASK,V,VyL,VyR)
    CALL CONSTRUCT_HO_Y(Mloc,Nloc,Nloc1,Ibeg,Iend,Jbeg,Jend,DY,MASK,CH,CHHyL,CHHyR)
    CALL CONSTRUCT_HO_Y(Mloc,Nloc,Nloc1,Ibeg,Iend,Jbeg,Jend,DY,MASK,CHY,CHYyL,CHYyR)
    CALL CONSTRUCT_HO_Y(Mloc,Nloc,Nloc1,Ibeg,Iend,Jbeg,Jend,DY,MASK,Eta,EtaRyL,EtaRyR)
    HyL=EtaRyL+Depthy+Zb
    HyR=EtaRyR+Depthy+Zb
    
    GsL(1:Mloc,1:Nloc1)=CHHyL(1:Mloc,1:Nloc1)*VyL(1:Mloc,1:Nloc1)*HyL(1:Mloc,1:Nloc1) &
    					-HyL(1:Mloc,1:Nloc1)* &
						CHYyL(1:Mloc,1:Nloc1)*5.93_SP*HyL(1:Mloc,1:Nloc1)*SQRT( &
						UyL(1:Mloc,1:Nloc1)*UyL(1:Mloc,1:Nloc1)+VyL(1:Mloc,1:Nloc1) &
						*VyL(1:Mloc,1:Nloc1))*0.4_SP/(-1+log(30._SP*MAX( &
						HyL(1:Mloc,1:Nloc1),MinDepth)/k_s)) 
	
	GsR(1:Mloc,1:Nloc1)=CHHyR(1:Mloc,1:Nloc1)*VyR(1:Mloc,1:Nloc1)*HyR(1:Mloc,1:Nloc1) &
						-HyR(1:Mloc,1:Nloc1)* &
						CHYyR(1:Mloc,1:Nloc1)*5.93_SP*HyR(1:Mloc,1:Nloc1)*SQRT( &
						UyR(1:Mloc,1:Nloc1)*UyR(1:Mloc,1:Nloc1)+VyR(1:Mloc,1:Nloc1) &
						*VyR(1:Mloc,1:Nloc1))*0.4_SP/(-1+log(30._SP*MAX( &
						HyR(1:Mloc,1:Nloc1),MinDepth)/k_s)) 
	
! calculate wave speed for Rieman solver
	CALL WAVE_SPEED(Mloc,Nloc,Mloc1,Nloc1,UxL,UxR,VyL,VyR,HxL,HxR,HyL,HyR, &
            SxL,SxR,SyL,SyR)  
            
! HLLC Rieman Solver
	CALL HLLC(Mloc1,Nloc,SxL,SxR,FsL,FsR,CHHxL*HxL,CHHxR*HxR,Fs)
    CALL HLLC(Mloc,Nloc1,SyL,SyR,GsL,GsR,CHHyL*HyL,CHHyR*HyR,Gs) 

!Boundary Conditions

! four sides of computational domain
!WEST BC
# if defined (PARALLEL)
        if ( n_west .eq. MPI_PROC_NULL ) then
# endif

# if defined (COUPLING)
	IF(IN_DOMAIN_WEST)THEN
		DO J=Jbeg,Kstart_WEST-1
      		Fs(Ibeg,J)=ZERO
      	ENDDO
     	DO J=Kend_WEST+1,Jend
      		Fs(Ibeg,J)=ZERO
      	ENDDO
	ENDIF
# else
	DO J=Jbeg,Jend
    	Fs(Ibeg,J)=ZERO
	ENDDO
# endif 

# if defined (PARALLEL)
      endif
# endif

!EAST BC
# if defined (PARALLEL)
        if ( n_east .eq. MPI_PROC_NULL ) then
# endif

# if defined (COUPLING)
	IF(IN_DOMAIN_EAST)THEN
		DO J=Jbeg,Kstart_EAST-1
			Fs(Iend1,J)=ZERO
		ENDDO
		DO J=Kend_EAST+1,Jend
    		Fs(Iend1,J)=ZERO
    	ENDDO
	ENDIF
# else
	DO J=Jbeg,Jend
    	Fs(Iend1,J)=ZERO
    ENDDO
# endif 

# if defined (PARALLEL)
      endif
# endif

!South BC
# if defined (PARALLEL)
      if ( n_suth .eq. MPI_PROC_NULL ) then
# endif

# if defined (COUPLING)
	IF(IN_DOMAIN_SOUTH)THEN
    	DO I=Ibeg,Kstart_SOUTH-1
      		Gs(I,Jbeg)=ZERO
     	ENDDO
     	DO I=Kend_SOUTH+1,Iend
      		Gs(I,Jbeg)=ZERO
     	ENDDO
	ENDIF
# else
    DO I=Ibeg,Iend
      Gs(I,Jbeg)=ZERO
    ENDDO
# endif  

# if defined (PARALLEL)
      endif
# endif

!North BC
# if defined (PARALLEL)
      if ( n_nrth .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
	IF(IN_DOMAIN_NORTH)THEN
		DO I=Ibeg,Kstart_NORTH-1
      		Gs(I,Jend1)=ZERO
      	ENDDO
     	DO I=Kend_NORTH+1,Iend
			Gs(I,Jend1)=ZERO
    	ENDDO
	ENDIF
# else
    DO I=Ibeg,Iend
    	Gs(I,Jend1)=ZERO
    ENDDO
# endif 
# if defined (PARALLEL)
     endif
# endif

! mask points (Dry Points)
! Jeff pointed out the loop should be Jbeg-1, Jend+1
! The problem is that the fluxes on the inter-processor boundaries may be
!modified if the point next to the boundary (e.g., in the ghost cells,
!managed by a different processor) is land, but as is the routine doesn't
!check for this. 

	DO j=Jbeg-1,Jend+1
		DO i=Ibeg-1,Iend+1
			IF(MASK(I,J)<1)THEN
        		Gs(I,J)=ZERO
        		Fs(I,J)=ZERO
			ENDIF
		ENDDO
	ENDDO
 ! finisher       
		
END SUBROUTINE sediment_flux

SUBROUTINE update_bed
!-----------------------------------------------------------
!   Update bed elevation 
!   Called by 
!      main
!    update: 08/29/2015, Babak Tehranirad
!    Update: 06/09/2017 Fengyan Shi, merge into 3.1beta
!-----------------------------------------------------------
    USE GLOBAL
    IMPLICIT NONE
	INTEGER::ISTEP
	Counter_s=Counter_s+1
		IF (Counter_s.NE.Morph_step) THEN
			DT1=DT+DT1
			DO J=Jbeg,Jend
				DO I=Ibeg,Iend
					P_ave(I,J)=P_ave(I,J)+Pickup(I,J)
				  IF(IN_Mask_s) THEN
				  IF (Zb(I,J).GE.(Zs(I,J)-0.001)) P_ave(I,J)=0.0_SP
				  ENDIF
				  	D_ave(I,J)=D_ave(I,J)+D(I,J)
			   	ENDDO
			ENDDO
		ELSE

			DT1=DT+DT1

			DO J=Jbeg,Jend
				DO I=Ibeg,Iend
				  P_ave(I,J)= (P_ave(I,J)+Pickup(I,J))/FLOAT(Morph_step)
!				  IF(H(I,J).LE.(MinDepth)*1.4) P_ave(I,J)=0.0_SP
				  IF(IN_Mask_s) THEN
				  IF (Zb(I,J).GE.(Zs(I,J)-0.001)) P_ave(I,J)=0.0_SP
				  ENDIF
				  D_ave(I,J)= ((D_ave(I,J)+D(I,J))/FLOAT(Morph_step))
				ENDDO
			ENDDO
			ZbOld=Zb
			counter_ava=counter_ava+1
			
			IF (MOD(counter_ava,2).EQ.1) THEN
				DO J=Jbeg,Jend
					DO I=Ibeg,Iend	
						DO ISTEP=1,3
							IF (H(I,J).GT.(1.2_SP*(MinDepth))) THEN
								htt= (P_ave(I,J)-D_ave(I,J))/(1.-n_porosity)
							ELSE
								htt=(-D_ave(I,J))/(1.-n_porosity)
							ENDIF
				        	Zb(I,J)=ALPHA(ISTEP)*ZbOld(I,J)+BETA(ISTEP)*(Zb(I,J)+DT1*htt) 							
						ENDDO
						P_ave(I,J)= 0.0
						D_ave(I,J)= 0.0
!						CH(I,J)=MAX(0.0,CH(I,J)*FLOAT(MASK(I,J)))						
						
						IF(IN_Mask_s) THEN
							IF(Zb(I,J).GE.(Zs(I,J)-0.001)) Zb(I,J)=Zs(I,J)   
							IF(Avalanche)THEN
!x-direction
								IF(ABS(Depth(I,J)+Zb(I,J)-Depth(I-1,J)-Zb(I-1,J))/DX.GT.tan_phi) THEN 
									IF(ABS(Depth(I,J)-Depth(I-1,J))/DX.LT.tan_phi) THEN
										IF(Zb(I,J).LT.Zs(I,J)) THEN
											IF(Zb(I-1,J).LT.Zs(I-1,J)) THEN
					u_c1=0.5_sp*(-Depth(I,J)-Zb(I,J)+Depth(I-1,J)+Zb(I-1,J)+DX*tan_phi)
												Zb(I-1,J)=MIN(Zb(I-1,J)-u_c1,Zs(I-1,J))
												Zb(I,J)=MIN(Zb(I,J)+u_c1,Zs(I,J))
												ava(I,J)=1.0_SP
											ELSEIF(Depth(I-1,J)+Zb(I-1,J).GT.Depth(I,J)+Zb(I,J)) THEN
					u_c1=0.5_sp*(-Depth(I,J)-Zb(I,J)+Depth(I-1,J)+Zb(I-1,J)+DX*tan_phi)
												Zb(I-1,J)=MIN(Zb(I-1,J)-u_c1,Zs(I-1,J))
												Zb(I,J)=MIN(Zb(I,J)+u_c1,Zs(I,J))
												ava(I,J)=1.0_SP
											ENDIF			
										ELSEIF(Depth(I-1,J)+Zb(I-1,J).LT.Depth(I,J)+Zb(I,J)) THEN
											IF(Zb(I-1,J).LT.Zs(I-1,J)) THEN	
					u_c1=0.5_sp*(-Depth(I,J)-Zb(I,J)+Depth(I-1,J)+Zb(I-1,J)+DX*tan_phi)
												Zb(I-1,J)=MIN(Zb(I-1,J)-u_c1,Zs(I-1,J))
												Zb(I,J)=MIN(Zb(I,J)+u_c1,Zs(I,J))
												ava(I,J)=1.0_SP	
											ENDIF			
										ENDIF
									ENDIF
								ENDIF				
!y-direction
								IF(ABS(Depth(I,J)+Zb(I,J)-Depth(I,J-1)-Zb(I,J-1))/DY.GT.tan_phi) THEN 
									IF(ABS(Depth(I,J)-Depth(I,J-1))/DY.LT.tan_phi) THEN
										IF(Zb(I,J).LT.Zs(I,J)) THEN
											IF(Zb(I,J-1).LT.Zs(I,J-1)) THEN
					u_c1=0.5_sp*(-Depth(I,J)-Zb(I,J)+Depth(I,J-1)+Zb(I,J-1)+DY*tan_phi)
												Zb(I,J-1)=MIN(Zb(I,J-1)-u_c1,Zs(I,J-1))
												Zb(I,J)=MIN(Zb(I,J)+u_c1,Zs(I,J))
												ava(I,J)=1.0_SP
											ELSEIF(Depth(I,J-1)+Zb(I,J-1).GT.Depth(I,J)+Zb(I,J)) THEN
					u_c1=0.5_sp*(-Depth(I,J)-Zb(I,J)+Depth(I,J-1)+Zb(I,J-1)+DY*tan_phi)
												Zb(I,J-1)=MIN(Zb(I,J-1)-u_c1,Zs(I,J-1))
												Zb(I,J)=MIN(Zb(I,J)+u_c1,Zs(I,J))
												ava(I,J)=1.0_SP
											ENDIF			
										ELSEIF(Depth(I,J-1)+Zb(I,J-1).LT.Depth(I,J)+Zb(I,J)) THEN
											IF(Zb(I,J-1).LT.Zs(I,J-1)) THEN	
					u_c1=0.5_sp*(-Depth(I,J)-Zb(I,J)+Depth(I,J-1)+Zb(I,J-1)+DY*tan_phi)
												Zb(I,J-1)=MIN(Zb(I,J-1)-u_c1,Zs(I,J-1))
												Zb(I,J)=MIN(Zb(I,J)+u_c1,Zs(I,J))
												ava(I,J)=1.0_SP
											ENDIF				
										ENDIF
									ENDIF
								ENDIF												
							ENDIF	
						ELSE
							IF(Avalanche)THEN
								IF(ABS(Depth(I,J)+Zb(I,J)-Depth(I-1,J)-Zb(I-1,J))/DX.GT.tan_phi) THEN 
									IF(ABS(Depth(I,J)-Depth(I-1,J))/DX.LT.tan_phi) THEN
					u_c1=0.5_sp*(-Depth(I,J)-Zb(I,J)+Depth(I-1,J)+Zb(I-1,J)+DX*tan_phi)
										Zb(I-1,J)=Zb(I-1,J)-u_c1
										Zb(I,J)=Zb(I,J)+u_c1
										ava(I,J)=1.0_SP						
									ENDIF
								ENDIF
								IF(ABS(Depth(I,J)+Zb(I,J)-Depth(I,J-1)-Zb(I,J-1))/DY.GT.tan_phi) THEN 
									IF(ABS(Depth(I,J)-Depth(I,J-1))/DY.LT.tan_phi) THEN
					u_c1=0.5_sp*(-Depth(I,J)-Zb(I,J)+Depth(I,J-1)+Zb(I,J-1)+DY*tan_phi)
										Zb(I,J-1)=Zb(I,J-1)-u_c1
										Zb(I,J)=Zb(I,J)+u_c1
										ava(I,J)=1.0_SP						
									ENDIF
								ENDIF	
							ENDIF
						ENDIF						
					ENDDO
				ENDDO
			ELSE

				DO I=Ibeg,Iend	
					DO J=Jbeg,Jend
						DO ISTEP=1,3
							IF (H(I,J).GT.(1.2_SP*(MinDepth))) THEN
								htt= (P_ave(I,J)-D_ave(I,J))/(1.-n_porosity)
							ELSE
								htt=(-D_ave(I,J))/(1.-n_porosity)
							ENDIF
				        	Zb(I,J)=ALPHA(ISTEP)*ZbOld(I,J)+BETA(ISTEP)*(Zb(I,J)+DT1*htt) 							
						ENDDO
						P_ave(I,J)= 0.0
						D_ave(I,J)= 0.0
!						CH(I,J)=MAX(0.0,CH(I,J)*FLOAT(MASK(I,J)))						
						
						IF(IN_Mask_s) THEN
							IF(Zb(I,J).GE.(Zs(I,J)-0.001)) Zb(I,J)=Zs(I,J)   
								IF(Avalanche)THEN
!x-direction
									IF(ABS(Depth(I,J)+Zb(I,J)-Depth(I-1,J)-Zb(I-1,J))/DX.GT.tan_phi) THEN 
										IF(ABS(Depth(I,J)-Depth(I-1,J))/DX.LT.tan_phi) THEN
											IF(Zb(I,J).LT.Zs(I,J)) THEN
												IF(Zb(I-1,J).LT.Zs(I-1,J)) THEN
					u_c1=0.5_sp*(-Depth(I,J)-Zb(I,J)+Depth(I-1,J)+Zb(I-1,J)+DX*tan_phi)
													Zb(I-1,J)=MIN(Zb(I-1,J)-u_c1,Zs(I-1,J))
													Zb(I,J)=MIN(Zb(I,J)+u_c1,Zs(I,J))
													ava(I,J)=1.0_SP
												ELSEIF(Depth(I-1,J)+Zb(I-1,J).GT.Depth(I,J)+Zb(I,J)) THEN
					u_c1=0.5_sp*(-Depth(I,J)-Zb(I,J)+Depth(I-1,J)+Zb(I-1,J)+DX*tan_phi)
													Zb(I-1,J)=MIN(Zb(I-1,J)-u_c1,Zs(I-1,J))
													Zb(I,J)=MIN(Zb(I,J)+u_c1,Zs(I,J))
													ava(I,J)=1.0_SP
												ENDIF			
											ELSEIF(Depth(I-1,J)+Zb(I-1,J).LT.Depth(I,J)+Zb(I,J)) THEN
												IF(Zb(I-1,J).LT.Zs(I-1,J)) THEN	
					u_c1=0.5_sp*(-Depth(I,J)-Zb(I,J)+Depth(I-1,J)+Zb(I-1,J)+DX*tan_phi)
													Zb(I-1,J)=MIN(Zb(I-1,J)-u_c1,Zs(I-1,J))
													Zb(I,J)=MIN(Zb(I,J)+u_c1,Zs(I,J))
													ava(I,J)=1.0_SP	
												ENDIF			
											ENDIF
										ENDIF
									ENDIF				
!y-direction
									IF(ABS(Depth(I,J)+Zb(I,J)-Depth(I,J-1)-Zb(I,J-1))/DY.GT.tan_phi) THEN 
										IF(ABS(Depth(I,J)-Depth(I,J-1))/DY.LT.tan_phi) THEN
											IF(Zb(I,J).LT.Zs(I,J)) THEN
												IF(Zb(I,J-1).LT.Zs(I,J-1)) THEN
					u_c1=0.5_sp*(-Depth(I,J)-Zb(I,J)+Depth(I,J-1)+Zb(I,J-1)+DY*tan_phi)
													Zb(I,J-1)=MIN(Zb(I,J-1)-u_c1,Zs(I,J-1))
													Zb(I,J)=MIN(Zb(I,J)+u_c1,Zs(I,J))
													ava(I,J)=1.0_SP
												ELSEIF(Depth(I,J-1)+Zb(I,J-1).GT.Depth(I,J)+Zb(I,J)) THEN
					u_c1=0.5_sp*(-Depth(I,J)-Zb(I,J)+Depth(I,J-1)+Zb(I,J-1)+DY*tan_phi)
													Zb(I,J-1)=MIN(Zb(I,J-1)-u_c1,Zs(I,J-1))
													Zb(I,J)=MIN(Zb(I,J)+u_c1,Zs(I,J))
													ava(I,J)=1.0_SP
												ENDIF			
											ELSEIF(Depth(I,J-1)+Zb(I,J-1).LT.Depth(I,J)+Zb(I,J)) THEN
												IF(Zb(I,J-1).LT.Zs(I,J-1)) THEN	
					u_c1=0.5_sp*(-Depth(I,J)-Zb(I,J)+Depth(I,J-1)+Zb(I,J-1)+DY*tan_phi)
													Zb(I,J-1)=MIN(Zb(I,J-1)-u_c1,Zs(I,J-1))
													Zb(I,J)=MIN(Zb(I,J)+u_c1,Zs(I,J))
													ava(I,J)=1.0_SP
												ENDIF				
											ENDIF
										ENDIF
									ENDIF												
								ENDIF	
						ELSE
							IF(Avalanche)THEN
								IF(ABS(Depth(I,J)+Zb(I,J)-Depth(I-1,J)-Zb(I-1,J))/DX.GT.tan_phi) THEN 
									IF(ABS(Depth(I,J)-Depth(I-1,J))/DX.LT.tan_phi) THEN
				u_c1=0.5_sp*(-Depth(I,J)-Zb(I,J)+Depth(I-1,J)+Zb(I-1,J)+DX*tan_phi)
										Zb(I-1,J)=Zb(I-1,J)-u_c1
										Zb(I,J)=Zb(I,J)+u_c1
										ava(I,J)=1.0_SP						
									ENDIF
								ENDIF
							IF(ABS(Depth(I,J)+Zb(I,J)-Depth(I,J-1)-Zb(I,J-1))/DY.GT.tan_phi) THEN 
								IF(ABS(Depth(I,J)-Depth(I,J-1))/DY.LT.tan_phi) THEN
				u_c1=0.5_sp*(-Depth(I,J)-Zb(I,J)+Depth(I,J-1)+Zb(I,J-1)+DY*tan_phi)
									Zb(I,J-1)=Zb(I,J-1)-u_c1
									Zb(I,J)=Zb(I,J)+u_c1
									ava(I,J)=1.0_SP						
								ENDIF
							ENDIF	
						ENDIF
						ENDIF						
					
					ENDDO
				ENDDO
		ENDIF
		    Counter_s=0
			DT1=0.0_SP
		ENDIF
	
! Avalanche
	IF(Avalanche)THEN
		c_dum=0
96 continue
		c_dum=c_dum+1
! Add limit here
		IF(c_dum.GT.1000) GO TO 86
		IF (MOD(counter_ava,2).EQ.1) THEN
			DO J=Jbeg,Jend
				DO I=Ibeg,Iend
					IF(IN_Mask_s) THEN
!x-direction
						IF(ABS(Depth(I,J)+Zb(I,J)-Depth(I-1,J)-Zb(I-1,J))/DX.GT.tan_phi) THEN 
							IF(ABS(Depth(I,J)-Depth(I-1,J))/DX.LT.tan_phi) THEN
								IF(Zb(I,J).LT.Zs(I,J)) THEN
									IF(Zb(I-1,J).LT.Zs(I-1,J)) THEN
					u_c1=0.5_sp*(-Depth(I,J)-Zb(I,J)+Depth(I-1,J)+Zb(I-1,J)+DX*tan_phi)
									Zb(I-1,J)=MIN(Zb(I-1,J)-u_c1,Zs(I-1,J))
									Zb(I,J)=MIN(Zb(I,J)+u_c1,Zs(I,J))
									ava(I,J)=1.0_SP
									GO TO 96
									ELSEIF(Depth(I-1,J)+Zb(I-1,J).GT.Depth(I,J)+Zb(I,J)) THEN
					u_c1=0.5_sp*(-Depth(I,J)-Zb(I,J)+Depth(I-1,J)+Zb(I-1,J)+DX*tan_phi)
									Zb(I-1,J)=MIN(Zb(I-1,J)-u_c1,Zs(I-1,J))
									Zb(I,J)=MIN(Zb(I,J)+u_c1,Zs(I,J))
									ava(I,J)=1.0_SP
									GO TO 96
									ENDIF			
								ELSEIF(Depth(I-1,J)+Zb(I-1,J).LT.Depth(I,J)+Zb(I,J)) THEN
									IF(Zb(I-1,J).LT.Zs(I-1,J)) THEN	
					u_c1=0.5_sp*(-Depth(I,J)-Zb(I,J)+Depth(I-1,J)+Zb(I-1,J)+DX*tan_phi)
									Zb(I-1,J)=MIN(Zb(I-1,J)-u_c1,Zs(I-1,J))
									Zb(I,J)=MIN(Zb(I,J)+u_c1,Zs(I,J))
									ava(I,J)=1.0_SP
									GO TO 96	
									ENDIF			
								ENDIF
							ENDIF
						ENDIF				
!y-direction
						IF(ABS(Depth(I,J)+Zb(I,J)-Depth(I,J-1)-Zb(I,J-1))/DY.GT.tan_phi) THEN 
							IF(ABS(Depth(I,J)-Depth(I,J-1))/DY.LT.tan_phi) THEN
								IF(Zb(I,J).LT.Zs(I,J)) THEN
									IF(Zb(I,J-1).LT.Zs(I,J-1)) THEN
					u_c1=0.5_sp*(-Depth(I,J)-Zb(I,J)+Depth(I,J-1)+Zb(I,J-1)+DY*tan_phi)
									Zb(I,J-1)=MIN(Zb(I,J-1)-u_c1,Zs(I,J-1))
									Zb(I,J)=MIN(Zb(I,J)+u_c1,Zs(I,J))
									ava(I,J)=1.0_SP
									GO TO 96
									ELSEIF(Depth(I,J-1)+Zb(I,J-1).GT.Depth(I,J)+Zb(I,J)) THEN
					u_c1=0.5_sp*(-Depth(I,J)-Zb(I,J)+Depth(I,J-1)+Zb(I,J-1)+DY*tan_phi)
									Zb(I,J-1)=MIN(Zb(I,J-1)-u_c1,Zs(I,J-1))
									Zb(I,J)=MIN(Zb(I,J)+u_c1,Zs(I,J))
									ava(I,J)=1.0_SP
									GO TO 96
									ENDIF			
								ELSEIF(Depth(I,J-1)+Zb(I,J-1).LT.Depth(I,J)+Zb(I,J)) THEN
									IF(Zb(I,J-1).LT.Zs(I,J-1)) THEN	
					u_c1=0.5_sp*(-Depth(I,J)-Zb(I,J)+Depth(I,J-1)+Zb(I,J-1)+DY*tan_phi)
									Zb(I,J-1)=MIN(Zb(I,J-1)-u_c1,Zs(I,J-1))
									Zb(I,J)=MIN(Zb(I,J)+u_c1,Zs(I,J))
									ava(I,J)=1.0_SP
									GO TO 96
									ENDIF				
								ENDIF
							ENDIF
						ENDIF	
					ELSE
						IF(ABS(Depth(I,J)+Zb(I,J)-Depth(I-1,J)-Zb(I-1,J))/DX.GT.tan_phi) THEN 
							IF(ABS(Depth(I,J)-Depth(I-1,J))/DX.LT.tan_phi) THEN
					u_c1=0.5_sp*(-Depth(I,J)-Zb(I,J)+Depth(I-1,J)+Zb(I-1,J)+DX*tan_phi)
							Zb(I-1,J)=Zb(I-1,J)-u_c1
							Zb(I,J)=Zb(I,J)+u_c1
							ava(I,J)=1.0_SP
							GO TO 96													
							ENDIF
						ENDIF
						IF(ABS(Depth(I,J)+Zb(I,J)-Depth(I,J-1)-Zb(I,J-1))/DY.GT.tan_phi) THEN 
							IF(ABS(Depth(I,J)-Depth(I,J-1))/DY.LT.tan_phi) THEN
					u_c1=0.5_sp*(-Depth(I,J)-Zb(I,J)+Depth(I,J-1)+Zb(I,J-1)+DY*tan_phi)
							Zb(I,J-1)=Zb(I,J-1)-u_c1
							Zb(I,J)=Zb(I,J)+u_c1
							ava(I,J)=1.0_SP	
							GO TO 96												
							ENDIF
						ENDIF	
					ENDIF
				ENDDO
			ENDDO
 		ELSE
 			DO I=Ibeg,Iend
 				DO J=Jbeg,Jend
					IF(IN_Mask_s) THEN
!x-direction
						IF(ABS(Depth(I,J)+Zb(I,J)-Depth(I-1,J)-Zb(I-1,J))/DX.GT.tan_phi) THEN 
							IF(ABS(Depth(I,J)-Depth(I-1,J))/DX.LT.tan_phi) THEN
								IF(Zb(I,J).LT.Zs(I,J)) THEN
									IF(Zb(I-1,J).LT.Zs(I-1,J)) THEN
					u_c1=0.5_sp*(-Depth(I,J)-Zb(I,J)+Depth(I-1,J)+Zb(I-1,J)+DX*tan_phi)
									Zb(I-1,J)=MIN(Zb(I-1,J)-u_c1,Zs(I-1,J))
									Zb(I,J)=MIN(Zb(I,J)+u_c1,Zs(I,J))
									ava(I,J)=1.0_SP
									GO TO 96
									ELSEIF(Depth(I-1,J)+Zb(I-1,J).GT.Depth(I,J)+Zb(I,J)) THEN
					u_c1=0.5_sp*(-Depth(I,J)-Zb(I,J)+Depth(I-1,J)+Zb(I-1,J)+DX*tan_phi)
									Zb(I-1,J)=MIN(Zb(I-1,J)-u_c1,Zs(I-1,J))
									Zb(I,J)=MIN(Zb(I,J)+u_c1,Zs(I,J))
									ava(I,J)=1.0_SP
									GO TO 96
									ENDIF			
								ELSEIF(Depth(I-1,J)+Zb(I-1,J).LT.Depth(I,J)+Zb(I,J)) THEN
									IF(Zb(I-1,J).LT.Zs(I-1,J)) THEN	
					u_c1=0.5_sp*(-Depth(I,J)-Zb(I,J)+Depth(I-1,J)+Zb(I-1,J)+DX*tan_phi)
									Zb(I-1,J)=MIN(Zb(I-1,J)-u_c1,Zs(I-1,J))
									Zb(I,J)=MIN(Zb(I,J)+u_c1,Zs(I,J))
									ava(I,J)=1.0_SP
									GO TO 96	
									ENDIF			
								ENDIF
							ENDIF
						ENDIF				
!y-direction
						IF(ABS(Depth(I,J)+Zb(I,J)-Depth(I,J-1)-Zb(I,J-1))/DY.GT.tan_phi) THEN 
							IF(ABS(Depth(I,J)-Depth(I,J-1))/DY.LT.tan_phi) THEN
								IF(Zb(I,J).LT.Zs(I,J)) THEN
									IF(Zb(I,J-1).LT.Zs(I,J-1)) THEN
					u_c1=0.5_sp*(-Depth(I,J)-Zb(I,J)+Depth(I,J-1)+Zb(I,J-1)+DY*tan_phi)
									Zb(I,J-1)=MIN(Zb(I,J-1)-u_c1,Zs(I,J-1))
									Zb(I,J)=MIN(Zb(I,J)+u_c1,Zs(I,J))
									ava(I,J)=1.0_SP
									GO TO 96
									ELSEIF(Depth(I,J-1)+Zb(I,J-1).GT.Depth(I,J)+Zb(I,J)) THEN
					u_c1=0.5_sp*(-Depth(I,J)-Zb(I,J)+Depth(I,J-1)+Zb(I,J-1)+DY*tan_phi)
									Zb(I,J-1)=MIN(Zb(I,J-1)-u_c1,Zs(I,J-1))
									Zb(I,J)=MIN(Zb(I,J)+u_c1,Zs(I,J))
									ava(I,J)=1.0_SP
									GO TO 96
									ENDIF			
								ELSEIF(Depth(I,J-1)+Zb(I,J-1).LT.Depth(I,J)+Zb(I,J)) THEN
									IF(Zb(I,J-1).LT.Zs(I,J-1)) THEN	
					u_c1=0.5_sp*(-Depth(I,J)-Zb(I,J)+Depth(I,J-1)+Zb(I,J-1)+DY*tan_phi)
									Zb(I,J-1)=MIN(Zb(I,J-1)-u_c1,Zs(I,J-1))
									Zb(I,J)=MIN(Zb(I,J)+u_c1,Zs(I,J))
									ava(I,J)=1.0_SP
									GO TO 96
									ENDIF				
								ENDIF
							ENDIF
						ENDIF	
					ELSE
						IF(ABS(Depth(I,J)+Zb(I,J)-Depth(I-1,J)-Zb(I-1,J))/DX.GT.tan_phi) THEN 
							IF(ABS(Depth(I,J)-Depth(I-1,J))/DX.LT.tan_phi) THEN
					u_c1=0.5_sp*(-Depth(I,J)-Zb(I,J)+Depth(I-1,J)+Zb(I-1,J)+DX*tan_phi)
							Zb(I-1,J)=Zb(I-1,J)-u_c1
							Zb(I,J)=Zb(I,J)+u_c1
							ava(I,J)=1.0_SP
							GO TO 96													
							ENDIF
						ENDIF
						IF(ABS(Depth(I,J)+Zb(I,J)-Depth(I,J-1)-Zb(I,J-1))/DY.GT.tan_phi) THEN 
							IF(ABS(Depth(I,J)-Depth(I,J-1))/DY.LT.tan_phi) THEN
					u_c1=0.5_sp*(-Depth(I,J)-Zb(I,J)+Depth(I,J-1)+Zb(I,J-1)+DY*tan_phi)
							Zb(I,J-1)=Zb(I,J-1)-u_c1
							Zb(I,J)=Zb(I,J)+u_c1
							ava(I,J)=1.0_SP	
							GO TO 96												
							ENDIF
						ENDIF	
					ENDIF
				ENDDO
			ENDDO
 		ENDIF
 	ENDIF
	86 continue
	END SUBROUTINE  update_bed
# endif

 
