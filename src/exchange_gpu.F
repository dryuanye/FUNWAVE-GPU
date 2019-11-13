module boundary_condition_module
    use cudafor
# if defined (MGPU)
    use PARAM
# else
    use PARAM,only: GRAV,ZERO,SP
# endif
    use GLOBAL,only: Ibeg,Iend,Iend1,Jbeg,Jend,Jend1,&
        Mloc,Nloc,Mloc1,Nloc1,NGhost,&
        PERIODIC,VISCOSITY_BREAKING,WAVEMAKER_VIS,DIFFUSION_SPONGE,&
        DISP_TIME_LEFT,Gamma2,Gamma3,&
# if defined (COUPLING)
        Kstart_WEST,Kend_WEST,Kstart_EAST,Kend_EAST,&
        Kstart_SOUTH,Kend_SOUTH,Kstart_NORTH,Kend_NORTH,&
        IN_DOMAIN_SOUTH,IN_DOMAIN_NORTH,IN_DOMAIN_WEST,IN_DOMAIN_EAST,&
# endif
        WaveMaker
# if defined (MGPU)
    use GLOBAL, only : n_east,n_west,n_suth,n_nrth,&
                       comm2d, ier,myid,PX,PY,&
                       NumberProcessor,ProcessorID
    use mgpu_utilities
# endif
    use mod_cuda,only: MASK_d,P_d,Q_d,Fx_d,Fy_d,Gx_d,Gy_d,&
        Eta_d,EtaRxL_d,EtaRxR_d,EtaRyL_d,EtaRyR_d,&
        Depthx_d,Depthy_d,AGE_BREAKING_d, nu_break_d,&
        U_d,HU_d,V_d,HV_d,Uxx_d,DUxx_d,Vyy_d,DVyy_d,Uxy_d,DUxy_d,Vxy_d,DVxy_d,&
        ETAT_d,ETATx_d,ETATy_d,Ut_d,Vt_d,Utx_d,Vty_d,Utxx_d,Vtyy_d,&
        Utxy_d,Vtxy_d,DUtxx_d,DVtyy_d,DUtxy_d,DVtxy_d,&
        Ux_d,DUx_d,Vy_d,DVy_d,Uy_d,DUy_d,Vx_d,DVx_d,&
        ETAx_d,ETAy_d,&
        n_left=>n_left_d,n_right=>n_right_d,n_bottom=>n_bottom_d,n_top=>n_top_d,&
        BlockDimX_2D,BlockDimY_2D,grid,tBlock,&
        streamID,istat

    implicit none
    interface PHI_COLL_GPU
        module procedure PHI_COLL_GPU
        module procedure PHI_COLL_INT_GPU
    end interface PHI_COLL_GPU

        


contains
!-------------------------------------------------------------------------------------
!
!    BOUNDARY_CONDITION is subroutine to provide boundary conditions at edges of domain
!
!    HISTORY: 
!      05/06/2010 Fengyan Shi
!      04/05/2011  Jeff Harris, corrected bugs in do-loop
!
! -----------------------------------------------------------------------------------
# if defined (MGPU)
attributes(global) subroutine boundary_condition_X_kernel(Ibeg,Iend,Jbeg,Jend,BlockDimX_2D,n_west,n_east,Gamma3)
# else
attributes(global) subroutine boundary_condition_X_kernel(Ibeg,Iend,Jbeg,Jend,BlockDimX_2D,Gamma3)
# endif
    implicit none
    integer,value :: Ibeg,Iend,Jbeg,Jend,BlockDimX_2D
# if defined (MGPU)
    integer,value :: n_west, n_east
# endif
    real(SP),value :: Gamma3
    integer :: I,J
! set local indexes
    j = threadIdx%x + (blockIdx%x-1)*BlockDimX_2D
!
    if (j>=Jbeg-1 .and. j<=Jend+1) then
        do i = Ibeg-1,Iend+1
            IF(MASK_d(I,J)<1)THEN
                P_d(I,J)=ZERO
# if defined (MGPU)
                IF((I/=Ibeg).or.(n_west.ne.MPI_PROC_NULL))THEN
# else
                IF(I/=Ibeg)THEN
# endif
                    Fx_d(I,J)=0.5_SP*GRAV* &
                        (EtaRxL_d(I,J)*EtaRxL_d(I,J)*Gamma3+ &
                        2.0_SP*EtaRxL_d(I,J)*Depthx_d(I,J))*MASK_d(I-1,J)
                ELSE
                    Fx_d(I,J)=ZERO
                ENDIF
                Gx_d(I,J)=ZERO
                P_d(I+1,J)=ZERO
# if defined (MGPU)
                IF((I/=Iend).or.(n_east.ne.MPI_PROC_NULL))THEN
# else
                IF(I/=Iend)THEN
# endif
                    Fx_d(I+1,J)=0.5_SP*GRAV* &
                        (EtaRxR_d(I+1,J)*EtaRxR_d(I+1,J)*Gamma3+ &
                      2.0_SP*EtaRxR_d(I+1,J)*Depthx_d(I+1,J))*MASK_d(I+1,J)
                ELSE
                    Fx_d(I+1,J)=ZERO
                ENDIF
                Gx_d(I+1,J)=ZERO
                Q_d(I,J)=ZERO
                Fy_d(I,J)=ZERO
            ENDIF
        enddo
    endif
end subroutine boundary_condition_X_kernel

# if defined (MGPU)
attributes(global) subroutine boundary_condition_Y_kernel(Ibeg,Iend,Jbeg,Jend,BlockDimX_2D,n_suth,n_nrth,Gamma3)
# else
attributes(global) subroutine boundary_condition_Y_kernel(Ibeg,Iend,Jbeg,Jend,BlockDimX_2D,Gamma3)
# endif
    implicit none
    integer,value :: Ibeg,Iend,Jbeg,Jend,BlockDimX_2D
    integer,value :: n_suth,n_nrth
    real(SP),value :: Gamma3

    integer :: I,J
    i = threadIdx%x + (blockIdx%x-1)*BlockDimX_2D
!
    if (i>=Ibeg-1 .and. i<=Iend+1) then
        do j = Jbeg-1,Jend+1
            IF(MASK_d(I,J)<1)THEN
# if defined (MGPU)
                IF((J/=Jbeg).or.(n_suth.ne.MPI_PROC_NULL))THEN
# else
                IF(J/=Jbeg)THEN
# endif
                    Gy_d(I,J)=0.5_SP*GRAV* &
                        (EtaRyL_d(I,J)*EtaRyL_d(I,J)*Gamma3+ &
                        2.0_SP*EtaRyL_d(I,J)*Depthy_d(I,J))*MASK_d(I,J-1)
                ELSE
                    Gy_d(I,J)=ZERO
                ENDIF
                Q_d(I,J+1)=ZERO
                Fy_d(I,J+1)=ZERO

# if defined (MGPU)
                IF((J/=Jend).or.(n_nrth.ne.MPI_PROC_NULL))THEN
# else
                IF(J/=Jend)THEN
# endif
                    Gy_d(I,J+1)=0.5_SP*GRAV* &
                      (EtaRyR_d(I,J+1)*EtaRyR_d(I,J+1)*Gamma3+ &
                      2.0_SP*EtaRyR_d(I,J+1)*Depthy_d(I,J+1))*MASK_d(I,J+1)
                ELSE
                    Gy_d(I,J+1)=ZERO
                ENDIF
            ENDIF
        enddo
    endif
end subroutine boundary_condition_Y_kernel

SUBROUTINE BOUNDARY_CONDITION_FLUXES_GPU
     IMPLICIT NONE
     INTEGER :: I,J

! four sides of computational domain

# if defined (MGPU)
        if ( n_west .eq. MPI_PROC_NULL ) then
# endif

# if defined (COUPLING)
   IF(IN_DOMAIN_WEST)THEN
!$cuf kernel do(1) <<<*,*>>>
     DO J=Jbeg,Kstart_WEST-1
         P_d(Ibeg,J)=ZERO
         Fx_d(Ibeg,J)=0.5_SP*GRAV*(EtaRxR_d(Ibeg,J)*EtaRxR_d(Ibeg,J)*Gamma3+ &
                 2.0_SP*EtaRxR_d(Ibeg,J)*Depthx_d(Ibeg,J))
         Gx_d(Ibeg,J)=ZERO
     ENDDO
!$cuf kernel do(1) <<<*,*>>>
     DO J=Kend_WEST+1,Jend
         P_d(Ibeg,J)=ZERO
         Fx_d(Ibeg,J)=0.5_SP*GRAV*(EtaRxR_d(Ibeg,J)*EtaRxR_d(Ibeg,J)*Gamma3+ &
              2.0_SP*EtaRxR_d(Ibeg,J)*Depthx_d(Ibeg,J))
         Gx_d(Ibeg,J)=ZERO
      ENDDO
   ENDIF
# else
   IF(WaveMaker(1:11)=='LEFT_BC_IRR')THEN
     ! do nothing
   ELSE
!$cuf kernel do(1) <<<*,*>>>
     DO J=Jbeg,Jend
         P_d(Ibeg,J)=ZERO
         Fx_d(Ibeg,J)=0.5_SP*GRAV*(EtaRxR_d(Ibeg,J)*EtaRxR_d(Ibeg,J)*Gamma3+ &
              2.0_SP*EtaRxR_d(Ibeg,J)*Depthx_d(Ibeg,J))
         Gx_d(Ibeg,J)=ZERO
      ENDDO
   ENDIF ! left bc wavemaker
# endif 

# if defined (MGPU)
      endif
# endif

# if defined (MGPU)
        if ( n_east .eq. MPI_PROC_NULL ) then
# endif

# if defined (COUPLING)
   IF(IN_DOMAIN_EAST)THEN
!$cuf kernel do(1) <<<*,*>>>
     DO J=Jbeg,Kstart_EAST-1
         P_d(Iend1,J)=ZERO
         Fx_d(Iend1,J)=0.5_SP*GRAV*(EtaRxL_d(Iend1,J)*EtaRxL_d(Iend1,J)*Gamma3 &
              +2.0_SP*EtaRxL_d(Iend1,J)*Depthx_d(Iend1,J))
         Gx_d(Iend1,J)=ZERO
     ENDDO
!$cuf kernel do(1) <<<*,*>>>
     DO J=Kend_EAST+1,Jend
         P_d(Iend1,J)=ZERO
         Fx_d(Iend1,J)=0.5_SP*GRAV*(EtaRxL_d(Iend1,J)*EtaRxL_d(Iend1,J)*Gamma3_d+ &
              2.0_SP*EtaRxL_d(Iend1,J)*Depthx_d(Iend1,J))
         Gx_d(Iend1,J)=ZERO
     ENDDO
   ENDIF
# else
!$cuf kernel do(1) <<<*,*>>>
     DO J=Jbeg,Jend
         P_d(Iend1,J)=ZERO
         Fx_d(Iend1,J)=0.5_SP*GRAV*(EtaRxL_d(Iend1,J)*EtaRxL_d(Iend1,J)*Gamma3+ &
              2.0_SP*EtaRxL_d(Iend1,J)*Depthx_d(Iend1,J))
         Gx_d(Iend1,J)=ZERO
     ENDDO
# endif 

# if defined (MGPU)
      endif
# endif

# if defined (CARTESIAN)
! y direction
   IF(PERIODIC)THEN
!   do nothing
   ELSE
# endif

# if defined (MGPU)
      if ( n_suth .eq. MPI_PROC_NULL ) then
# endif

# if defined (COUPLING)
   IF(IN_DOMAIN_SOUTH)THEN
!$cuf kernel do(1) <<<*,*>>>
     DO I=Ibeg,Kstart_SOUTH-1
         Q_d(I,Jbeg)=ZERO
         Fy_d(I,Jbeg)=ZERO
         Gy_d(I,Jbeg)=0.5_SP*GRAV*(EtaRyR_d(I,Jbeg)*EtaRyR_d(I,Jbeg)*Gamma3+ &
              2.0_SP*EtaRyR_d(I,Jbeg)*Depthy_d(I,Jbeg))
     ENDDO
!$cuf kernel do(1) <<<*,*>>>
     DO I=Kend_SOUTH+1,Iend
         Q_d(I,Jbeg)=ZERO
         Fy_d(I,Jbeg)=ZERO
         Gy_d(I,Jbeg)=0.5_SP*GRAV*(EtaRyR_d(I,Jbeg)*EtaRyR_d(I,Jbeg)*Gamma3+ &
              2.0_SP*EtaRyR_d(I,Jbeg)*Depthy_d(I,Jbeg))
     ENDDO
   ENDIF
# else
!$cuf kernel do(1) <<<*,*>>>
     DO I=Ibeg,Iend
         Q_d(I,Jbeg)=ZERO
         Fy_d(I,Jbeg)=ZERO
         Gy_d(I,Jbeg)=0.5_SP*GRAV*(EtaRyR_d(I,Jbeg)*EtaRyR_d(I,Jbeg)*Gamma3+ &
              2.0_SP*EtaRyR_d(I,Jbeg)*Depthy_d(I,Jbeg))
      ENDDO
# endif  

# if defined (MGPU)
      endif
# endif
# if defined (MGPU)
      if ( n_nrth .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
   IF(IN_DOMAIN_NORTH)THEN
!$cuf kernel do(1) <<<*,*>>>
     DO I=Ibeg,Kstart_NORTH-1
         Q_d(I,Jend1)=ZERO
         Fy_d(I,Jend1)=ZERO
         Gy_d(I,Jend1)=0.5_SP*GRAV*(EtaRyL_d(I,Jend1)*EtaRyL_d(I,Jend1)*Gamma3+ &
              2.0_SP*EtaRyL_d(I,Jend1)*Depthy_d(I,Jend1))
     ENDDO
!$cuf kernel do(1) <<<*,*>>>
     DO I=Kend_NORTH+1,Iend
         Q_d(I,Jend1)=ZERO
         Fy_d(I,Jend1)=ZERO
         Gy_d(I,Jend1)=0.5_SP*GRAV*(EtaRyL_d(I,Jend1)*EtaRyL_d(I,Jend1)*Gamma3+ &
              2.0_SP*EtaRyL_d(I,Jend1)*Depthy_d(I,Jend1))
     ENDDO
   ENDIF
# else
!$cuf kernel do(1) <<<*,*>>>
     DO I=Ibeg,Iend
         Q_d(I,Jend1)=ZERO
         Fy_d(I,Jend1)=ZERO
         Gy_d(I,Jend1)=0.5_SP*GRAV*(EtaRyL_d(I,Jend1)*EtaRyL_d(I,Jend1)*Gamma3+ &
              2.0_SP*EtaRyL_d(I,Jend1)*Depthy_d(I,Jend1))
     ENDDO
# endif 
# if defined (MGPU)
     endif
# endif

# if defined (CARTESIAN)
    ENDIF
# endif


    tBlock = dim3 (BlockDimX_2D,1,1)
    grid = dim3 ( ceiling ( real ( Nloc ) / BlockDimX_2D ) ,1,1)
# if defined (MGPU)
    call boundary_condition_X_kernel<<<grid, tBlock,0,streamID(2)>>>&
        (Ibeg,Iend,Jbeg,Jend,BlockDimX_2D,n_west,n_east,Gamma3)
# else
    call boundary_condition_X_kernel<<<grid, tBlock,0,streamID(2)>>>&
        (Ibeg,Iend,Jbeg,Jend,BlockDimX_2D,Gamma3)
# endif
    grid = dim3 ( ceiling ( real ( Mloc ) / BlockDimX_2D ) ,1,1)
# if defined (MGPU)
    call boundary_condition_Y_kernel<<<grid,tBlock,0,streamID(3)>>>&
            (Ibeg,Iend,Jbeg,Jend,BlockDimX_2D,n_suth,n_nrth,Gamma3)
# else
    call boundary_condition_Y_kernel<<<grid,tBlock,0,streamID(3)>>>&
            (Ibeg,Iend,Jbeg,Jend,BlockDimX_2D,Gamma3)
# endif
    istat = cudaDeviceSynchronize()

END SUBROUTINE BOUNDARY_CONDITION_FLUXES_GPU

!-------------------------------------------------------------------
!
!   This subroutine is used to collect data into ghost cells                                                         
!
!   HISTORY:
!   07/09/2010 Fengyan Shi, use dummy variables 2) add vtype=3
!
!-------------------------------------------------------------------
SUBROUTINE EXCHANGE_DISPERSION_GPU
    IMPLICIT NONE
    INTEGER :: VTYPE
    integer :: i,j,k

# if defined (CARTESIAN)
! first period
    if (PERIODIC) then
        call periodic_NSexchange_gpu(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Uxx_d)
        call periodic_NSexchange_gpu(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DUxx_d)
        call periodic_NSexchange_gpu(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Vyy_d)
        call periodic_NSexchange_gpu(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DVyy_d)
        call periodic_NSexchange_gpu(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Uxy_d)
        call periodic_NSexchange_gpu(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DUxy_d)
        call periodic_NSexchange_gpu(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Vxy_d)
        call periodic_NSexchange_gpu(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DVxy_d)
        IF(Gamma2>ZERO)THEN
            IF(DISP_TIME_LEFT)THEN
                call periodic_NSexchange_gpu(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,ETAT_d)
            ELSE
                call periodic_NSexchange_gpu(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Ut_d)
                call periodic_NSexchange_gpu(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Vt_d)
                call periodic_NSexchange_gpu(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Utx_d)
                call periodic_NSexchange_gpu(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Vty_d)
                call periodic_NSexchange_gpu(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Utxx_d)
                call periodic_NSexchange_gpu(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Vtyy_d)
                call periodic_NSexchange_gpu(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Utxy_d)
                call periodic_NSexchange_gpu(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Vtxy_d)
                call periodic_NSexchange_gpu(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DUtxx_d)
                call periodic_NSexchange_gpu(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DVtyy_d)
                call periodic_NSexchange_gpu(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DUtxy_d)
                call periodic_NSexchange_gpu(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DVtxy_d)
            ENDIF
            call periodic_NSexchange_gpu(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Ux_d)
            call periodic_NSexchange_gpu(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DUx_d)
            call periodic_NSexchange_gpu(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Vy_d)
            call periodic_NSexchange_gpu(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DVy_d)
            call periodic_NSexchange_gpu(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Uy_d)
            call periodic_NSexchange_gpu(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DUy_d)
            call periodic_NSexchange_gpu(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Vx_d)
            call periodic_NSexchange_gpu(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DVx_d)
            call periodic_NSexchange_gpu(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,ETAx_d)
            call periodic_NSexchange_gpu(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,ETAy_d)
        ENDIF

    endif
# endif
!
      ! x-direction
# if defined (MGPU)
    if ( n_west .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_WEST)THEN
# endif

    IF(WaveMaker(1:11)=='LEFT_BC_IRR')THEN
      ! do nothing
    ELSE
!$cuf kernel do(2) <<<*,(4,64),stream=streamID(2)>>>
      DO J=Jbeg,Jend
      DO K=1,Nghost
        !if (VTYPE==2) then
# if defined (CARTESIAN)
          Uxx_d(K,J)=-Uxx_d(Ibeg+Nghost-K,J)
          DUxx_d(K,J)=-DUxx_d(Ibeg+Nghost-K,J)
          if (Gamma2>ZERO) then
              IF(DISP_TIME_LEFT)THEN
              
              ELSE
                  Ut_d(K,J)=-Ut_d(Ibeg+Nghost-K,J)
                  Utxx_d(K,J)=-Utxx_d(Ibeg+Nghost-K,J)
                  DUtxx_d(K,J)=-DUtxx_d(Ibeg+Nghost-K,J)
              ENDIF
              Vx_d(K,J)=-Vx_d(Ibeg+Nghost-K,J)
              DVx_d(K,J)=-DVx_d(Ibeg+Nghost-K,J)
              ETAx_d(K,J)=-ETAx_d(Ibeg+Nghost-K,J)
          endif
# else
          Uxx_d(K,J)=-Uxx_d(Ibeg+Nghost-K,J)
          DUxx_d(K,J)=-DUxx_d(Ibeg+Nghost-K,J)
          Vx_d(K,J)=-Vx_d(Ibeg+Nghost-K,J)
          DVx_d(K,J)=-DVx_d(Ibeg+Nghost-K,J)
          ETAx_d(K,J)=-ETAx_d(Ibeg+Nghost-K,J)
# endif
        !else if (VTYPE==3) then
# if defined (CARTESIAN)
          Vyy_d(K,J)=Vyy_d(Ibeg+Nghost-K,J)
          DVyy_d(K,J)=DVyy_d(Ibeg+Nghost-K,J)
          if (Gamma2>ZERO) then
              IF(DISP_TIME_LEFT)THEN
              
              ELSE
                  Vt_d(K,J)=Vt_d(Ibeg+Nghost-K,J)
                  Vtyy_d(K,J)=Vtyy_d(Ibeg+Nghost-K,J)
                  DVtyy_d(K,J)=DVtyy_d(Ibeg+Nghost-K,J)
              ENDIF
              Uy_d(K,J)=Uy_d(Ibeg+Nghost-K,J)
              DUy_d(K,J)=DUy_d(Ibeg+Nghost-K,J)
              ETAy_d(K,J)=ETAy_d(Ibeg+Nghost-K,J)
          endif
# else
          Vyy_d(K,J)=Vyy_d(Ibeg+Nghost-K,J)
          DVyy_d(K,J)=DVyy_d(Ibeg+Nghost-K,J)
          Uy_d(K,J)=Uy_d(Ibeg+Nghost-K,J)
          DUy_d(K,J)=DUy_d(Ibeg+Nghost-K,J)
          ETAy_d(K,J)=ETAy_d(Ibeg+Nghost-K,J)
# endif
        !else if (VTYPE==1) then
# if defined (CARTESIAN)
          Uxy_d(K,J)=Uxy_d(Ibeg+Nghost-K,J)
          DUxy_d(K,J)=DUxy_d(Ibeg+Nghost-K,J)
          Vxy_d(K,J)=Vxy_d(Ibeg+Nghost-K,J)
          DVxy_d(K,J)=DVxy_d(Ibeg+Nghost-K,J)
          if (Gamma2>ZERO) then
              IF(DISP_TIME_LEFT)THEN
                  ETAT_d(K,J)=ETAT_d(Ibeg+Nghost-K,J)
              ELSE
                  Utx_d(K,J)=Utx_d(Ibeg+Nghost-K,J)
                  Vty_d(K,J)=Vty_d(Ibeg+Nghost-K,J)
                  Utxy_d(K,J)=Utxy_d(Ibeg+Nghost-K,J)
                  Vtxy_d(K,J)=Vtxy_d(Ibeg+Nghost-K,J)
                  DUtxy_d(K,J)=DUtxy_d(Ibeg+Nghost-K,J)
                  DVtxy_d(K,J)=DVtxy_d(Ibeg+Nghost-K,J)
              ENDIF
              Ux_d(K,J)=Ux_d(Ibeg+Nghost-K,J)
              DUx_d(K,J)=DUx_d(Ibeg+Nghost-K,J)
              Vy_d(K,J)=Vy_d(Ibeg+Nghost-K,J)
              DVy_d(K,J)=DVy_d(Ibeg+Nghost-K,J)
          endif
# else
          Uxy_d(K,J)=Uxy_d(Ibeg+Nghost-K,J)
          DUxy_d(K,J)=DUxy_d(Ibeg+Nghost-K,J)
          Vxy_d(K,J)=Vxy_d(Ibeg+Nghost-K,J)
          DVxy_d(K,J)=DVxy_d(Ibeg+Nghost-K,J)
          Ux_d(K,J)=Ux_d(Ibeg+Nghost-K,J)
          DUx_d(K,J)=DUx_d(Ibeg+Nghost-K,J)
          Vy_d(K,J)=Vy_d(Ibeg+Nghost-K,J)
          DVy_d(K,J)=DVy_d(Ibeg+Nghost-K,J)
# endif
      ENDDO
      ENDDO
    ENDIF ! end left_bc wavemaker
# if defined (COUPLING)
    ENDIF
# endif
# if defined (MGPU)
    endif
# endif

# if defined (MGPU)
    if ( n_east .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_EAST)THEN
# endif
!$cuf kernel do(2) <<<*,(4,64),stream=streamID(3)>>>
      DO J=Jbeg,Jend
      DO K=1,Nghost
        !if (VTYPE==2) then
# if defined (CARTESIAN)
          Uxx_d(Iend+K,J)=-Uxx_d(Iend-K+1,J)
          DUxx_d(Iend+K,J)=-DUxx_d(Iend-K+1,J)
          if (Gamma2>ZERO) then
              IF(DISP_TIME_LEFT)THEN
              
              ELSE
                  Ut_d(Iend+K,J)=-Ut_d(Iend-K+1,J)
                  Utxx_d(Iend+K,J)=-Utxx_d(Iend-K+1,J)
                  DUtxx_d(Iend+K,J)=-DUtxx_d(Iend-K+1,J)
              ENDIF
              Vx_d(Iend+K,J)=-Vx_d(Iend-K+1,J)
              DVx_d(Iend+K,J)=-DVx_d(Iend-K+1,J)
              ETAx_d(Iend+K,J)=-ETAx_d(Iend-K+1,J)
          endif
# else
          Uxx_d(Iend+K,J)=-Uxx_d(Iend-K+1,J)
          DUxx_d(Iend+K,J)=-DUxx_d(Iend-K+1,J)
          Vx_d(Iend+K,J)=-Vx_d(Iend-K+1,J)
          DVx_d(Iend+K,J)=-DVx_d(Iend-K+1,J)
          ETAx_d(Iend+K,J)=-ETAx_d(Iend-K+1,J)
# endif
        !else if (VTYPE==3) then
# if defined (CARTESIAN)
          Vyy_d(Iend+K,J)=Vyy_d(Iend-K+1,J)
          DVyy_d(Iend+K,J)=DVyy_d(Iend-K+1,J)
          if (Gamma2>ZERO) then
              IF(DISP_TIME_LEFT)THEN
              
              ELSE
                  Vt_d(Iend+K,J)=Vt_d(Iend-K+1,J)
                  Vtyy_d(Iend+K,J)=Vtyy_d(Iend-K+1,J)
                  DVtyy_d(Iend+K,J)=DVtyy_d(Iend-K+1,J)
              ENDIF
              Uy_d(Iend+K,J)=Uy_d(Iend-K+1,J)
              DUy_d(Iend+K,J)=DUy_d(Iend-K+1,J)
              ETAy_d(Iend+K,J)=ETAy_d(Iend-K+1,J)
          endif
# else
          Vyy_d(Iend+K,J)=Vyy_d(Iend-K+1,J)
          DVyy_d(Iend+K,J)=DVyy_d(Iend-K+1,J)
          Uy_d(Iend+K,J)=Uy_d(Iend-K+1,J)
          DUy_d(Iend+K,J)=DUy_d(Iend-K+1,J)
          ETAy_d(Iend+K,J)=ETAy_d(Iend-K+1,J)
# endif
        !else if (VTYPE==1) then
# if defined (CARTESIAN)
          Uxy_d(Iend+K,J)=Uxy_d(Iend-K+1,J)
          DUxy_d(Iend+K,J)=DUxy_d(Iend-K+1,J)
          Vxy_d(Iend+K,J)=Vxy_d(Iend-K+1,J)
          DVxy_d(Iend+K,J)=DVxy_d(Iend-K+1,J)
          if (Gamma2>ZERO) then
              IF(DISP_TIME_LEFT)THEN
                  ETAT_d(Iend+K,J)=ETAT_d(Iend-K+1,J)
              ELSE
                  Utx_d(Iend+K,J)=Utx_d(Iend-K+1,J)
                  Vty_d(Iend+K,J)=Vty_d(Iend-K+1,J)
                  Utxy_d(Iend+K,J)=Utxy_d(Iend-K+1,J)
                  Vtxy_d(Iend+K,J)=Vtxy_d(Iend-K+1,J)
                  DUtxy_d(Iend+K,J)=DUtxy_d(Iend-K+1,J)
                  DVtxy_d(Iend+K,J)=DVtxy_d(Iend-K+1,J)
              ENDIF
              Ux_d(Iend+K,J)=Ux_d(Iend-K+1,J)
              DUx_d(Iend+K,J)=DUx_d(Iend-K+1,J)
              Vy_d(Iend+K,J)=Vy_d(Iend-K+1,J)
              DVy_d(Iend+K,J)=DVy_d(Iend-K+1,J)
          endif
# else
          Uxy_d(Iend+K,J)=Uxy_d(Iend-K+1,J)
          DUxy_d(Iend+K,J)=DUxy_d(Iend-K+1,J)
          Vxy_d(Iend+K,J)=Vxy_d(Iend-K+1,J)
          DVxy_d(Iend+K,J)=DVxy_d(Iend-K+1,J)
          Ux_d(Iend+K,J)=Ux_d(Iend-K+1,J)
          DUx_d(Iend+K,J)=DUx_d(Iend-K+1,J)
          Vy_d(Iend+K,J)=Vy_d(Iend-K+1,J)
          DVy_d(Iend+K,J)=DVy_d(Iend-K+1,J)
# endif
      ENDDO
      ENDDO
# if defined (COUPLING)
    ENDIF
# endif
# if defined (MGPU)
    endif
# endif

      ! y-direction and corners
# if defined (CARTESIAN)
      IF(.NOT.PERIODIC)THEN
# endif 

# if defined (MGPU)
    if ( n_suth .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_SOUTH)THEN
# endif
!$cuf kernel do(2) <<<*,(64,4),stream=streamID(4)>>>
      DO K=1,Nghost
      DO I=1,Mloc
        !if (VTYPE==2) then
# if defined (CARTESIAN)
          Uxx_d(I,K)=Uxx_d(I,Jbeg+Nghost-K)
          DUxx_d(I,K)=DUxx_d(I,Jbeg+Nghost-K)
          if (Gamma2>ZERO) then
              IF(DISP_TIME_LEFT)THEN
              
              ELSE
                  Ut_d(I,K)=Ut_d(I,Jbeg+Nghost-K)
                  Utxx_d(I,K)=Utxx_d(I,Jbeg+Nghost-K)
                  DUtxx_d(I,K)=DUtxx_d(I,Jbeg+Nghost-K)
              ENDIF
              Vx_d(I,K)=Vx_d(I,Jbeg+Nghost-K)
              DVx_d(I,K)=DVx_d(I,Jbeg+Nghost-K)
              ETAx_d(I,K)=ETAx_d(I,Jbeg+Nghost-K)
          endif
# else
          Uxx_d(I,K)=Uxx_d(I,Jbeg+Nghost-K)
          DUxx_d(I,K)=DUxx_d(I,Jbeg+Nghost-K)
          Vx_d(I,K)=Vx_d(I,Jbeg+Nghost-K)
          DVx_d(I,K)=DVx_d(I,Jbeg+Nghost-K)
          ETAx_d(I,K)=ETAx_d(I,Jbeg+Nghost-K)
# endif
        !else if (VTYPE==3) then
# if defined (CARTESIAN)
          Vyy_d(I,K)=-Vyy_d(I,Jbeg+Nghost-K)
          DVyy_d(I,K)=-DVyy_d(I,Jbeg+Nghost-K)
          if (Gamma2>ZERO) then
              IF(DISP_TIME_LEFT)THEN
              
              ELSE
                  Vt_d(I,K)=-Vt_d(I,Jbeg+Nghost-K)
                  Vtyy_d(I,K)=-Vtyy_d(I,Jbeg+Nghost-K)
                  DVtyy_d(I,K)=-DVtyy_d(I,Jbeg+Nghost-K)
              ENDIF
              Uy_d(I,K)=-Uy_d(I,Jbeg+Nghost-K)
              DUy_d(I,K)=-DUy_d(I,Jbeg+Nghost-K)
              ETAy_d(I,K)=-ETAy_d(I,Jbeg+Nghost-K)
          endif
# else
          Vyy_d(I,K)=-Vyy_d(I,Jbeg+Nghost-K)
          DVyy_d(I,K)=-DVyy_d(I,Jbeg+Nghost-K)
          Uy_d(I,K)=-Uy_d(I,Jbeg+Nghost-K)
          DUy_d(I,K)=-DUy_d(I,Jbeg+Nghost-K)
          ETAy_d(I,K)=-ETAy_d(I,Jbeg+Nghost-K)
# endif
        !else if (VTYPE==1) then
# if defined (CARTESIAN)
          Uxy_d(I,K)=Uxy_d(I,Jbeg+Nghost-K)
          DUxy_d(I,K)=DUxy_d(I,Jbeg+Nghost-K)
          Vxy_d(I,K)=Vxy_d(I,Jbeg+Nghost-K)
          DVxy_d(I,K)=DVxy_d(I,Jbeg+Nghost-K)
          if (Gamma2>ZERO) then
              IF(DISP_TIME_LEFT)THEN
                  ETAT_d(I,K)=ETAT_d(I,Jbeg+Nghost-K)
              ELSE
                  Utx_d(I,K)=Utx_d(I,Jbeg+Nghost-K)
                  Vty_d(I,K)=Vty_d(I,Jbeg+Nghost-K)
                  Utxy_d(I,K)=Utxy_d(I,Jbeg+Nghost-K)
                  Vtxy_d(I,K)=Vtxy_d(I,Jbeg+Nghost-K)
                  DUtxy_d(I,K)=DUtxy_d(I,Jbeg+Nghost-K)
                  DVtxy_d(I,K)=DVtxy_d(I,Jbeg+Nghost-K)
              ENDIF
              Ux_d(I,K)=Ux_d(I,Jbeg+Nghost-K)
              DUx_d(I,K)=DUx_d(I,Jbeg+Nghost-K)
              Vy_d(I,K)=Vy_d(I,Jbeg+Nghost-K)
              DVy_d(I,K)=DVy_d(I,Jbeg+Nghost-K)
          endif
# else
          Uxy_d(I,K)=Uxy_d(I,Jbeg+Nghost-K)
          DUxy_d(I,K)=DUxy_d(I,Jbeg+Nghost-K)
          Vxy_d(I,K)=Vxy_d(I,Jbeg+Nghost-K)
          DVxy_d(I,K)=DVxy_d(I,Jbeg+Nghost-K)
          Ux_d(I,K)=Ux_d(I,Jbeg+Nghost-K)
          DUx_d(I,K)=DUx_d(I,Jbeg+Nghost-K)
          Vy_d(I,K)=Vy_d(I,Jbeg+Nghost-K)
          DVy_d(I,K)=DVy_d(I,Jbeg+Nghost-K)
# endif
      ENDDO
      ENDDO
# if defined (COUPLING)
    ENDIF
# endif
# if defined (MGPU)
    endif
# endif

# if defined (MGPU)
    if ( n_nrth .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_NORTH)THEN
# endif
!$cuf kernel do(2) <<<*,(64,4),stream=streamID(5)>>>
      DO K=1,Nghost
      DO I=1,Mloc
        !if (VTYPE==2) then
# if defined (CARTESIAN)
          Uxx_d(I,Jend+K)=Uxx_d(I,Jend-K+1)
          DUxx_d(I,Jend+K)=DUxx_d(I,Jend-K+1)
          if (Gamma2>ZERO) then
              IF(DISP_TIME_LEFT)THEN
              
              ELSE
                  Ut_d(I,Jend+K)=Ut_d(I,Jend-K+1)
                  Utxx_d(I,Jend+K)=Utxx_d(I,Jend-K+1)
                  DUtxx_d(I,Jend+K)=DUtxx_d(I,Jend-K+1)
              ENDIF
              Vx_d(I,Jend+K)=Vx_d(I,Jend-K+1)
              DVx_d(I,Jend+K)=DVx_d(I,Jend-K+1)
              ETAx_d(I,Jend+K)=ETAx_d(I,Jend-K+1)
          endif
# else
          Uxx_d(I,Jend+K)=Uxx_d(I,Jend-K+1)
          DUxx_d(I,Jend+K)=DUxx_d(I,Jend-K+1)
          Vx_d(I,Jend+K)=Vx_d(I,Jend-K+1)
          DVx_d(I,Jend+K)=DVx_d(I,Jend-K+1)
          ETAx_d(I,Jend+K)=ETAx_d(I,Jend-K+1)
# endif
        !else if (VTYPE==3) then
# if defined (CARTESIAN)
          Vyy_d(I,Jend+K)=-Vyy_d(I,Jend-K+1)
          DVyy_d(I,Jend+K)=-DVyy_d(I,Jend-K+1)
          if (Gamma2>ZERO) then
              IF(DISP_TIME_LEFT)THEN
              
              ELSE
                  Vt_d(I,Jend+K)=-Vt_d(I,Jend-K+1)
                  Vtyy_d(I,Jend+K)=-Vtyy_d(I,Jend-K+1)
                  DVtyy_d(I,Jend+K)=-DVtyy_d(I,Jend-K+1)
              ENDIF
              Uy_d(I,Jend+K)=-Uy_d(I,Jend-K+1)
              DUy_d(I,Jend+K)=-DUy_d(I,Jend-K+1)
              ETAy_d(I,Jend+K)=-ETAy_d(I,Jend-K+1)
          endif
# else
          Vyy_d(I,Jend+K)=-Vyy_d(I,Jend-K+1)
          DVyy_d(I,Jend+K)=-DVyy_d(I,Jend-K+1)
          Uy_d(I,Jend+K)=-Uy_d(I,Jend-K+1)
          DUy_d(I,Jend+K)=-DUy_d(I,Jend-K+1)
          ETAy_d(I,Jend+K)=-ETAy_d(I,Jend-K+1)
# endif
        !else if (VTYPE==1) then
# if defined (CARTESIAN)
          Uxy_d(I,Jend+K)=Uxy_d(I,Jend-K+1)
          DUxy_d(I,Jend+K)=DUxy_d(I,Jend-K+1)
          Vxy_d(I,Jend+K)=Vxy_d(I,Jend-K+1)
          DVxy_d(I,Jend+K)=DVxy_d(I,Jend-K+1)
          if (Gamma2>ZERO) then
              IF(DISP_TIME_LEFT)THEN
                  ETAT_d(I,Jend+K)=ETAT_d(I,Jend-K+1)
              ELSE
                  Utx_d(I,Jend+K)=Utx_d(I,Jend-K+1)
                  Vty_d(I,Jend+K)=Vty_d(I,Jend-K+1)
                  Utxy_d(I,Jend+K)=Utxy_d(I,Jend-K+1)
                  Vtxy_d(I,Jend+K)=Vtxy_d(I,Jend-K+1)
                  DUtxy_d(I,Jend+K)=DUtxy_d(I,Jend-K+1)
                  DVtxy_d(I,Jend+K)=DVtxy_d(I,Jend-K+1)
              ENDIF
              Ux_d(I,Jend+K)=Ux_d(I,Jend-K+1)
              DUx_d(I,Jend+K)=DUx_d(I,Jend-K+1)
              Vy_d(I,Jend+K)=Vy_d(I,Jend-K+1)
              DVy_d(I,Jend+K)=DVy_d(I,Jend-K+1)
          endif
# else
          Uxy_d(I,Jend+K)=Uxy_d(I,Jend-K+1)
          DUxy_d(I,Jend+K)=DUxy_d(I,Jend-K+1)
          Vxy_d(I,Jend+K)=Vxy_d(I,Jend-K+1)
          DVxy_d(I,Jend+K)=DVxy_d(I,Jend-K+1)
          Ux_d(I,Jend+K)=Ux_d(I,Jend-K+1)
          DUx_d(I,Jend+K)=DUx_d(I,Jend-K+1)
          Vy_d(I,Jend+K)=Vy_d(I,Jend-K+1)
          DVy_d(I,Jend+K)=DVy_d(I,Jend-K+1)
# endif
      ENDDO
      ENDDO
# if defined (COUPLING)
    ENDIF
# endif
# if defined (MGPU)
    endif
# endif

# if defined (CARTESIAN)
      ENDIF   ! end NOT PERIODIC
# endif


END SUBROUTINE EXCHANGE_DISPERSION_GPU

!---------------------------------------------------------------------------------------
!
!   EXCHANGE subroutine is used to collect data into ghost cells                                                         
!
!   HISTORY:
!   07/09/2010 Fengyan Shi 
!     1) use dummy variables 2) add vtype=3
!   08/19/2015 Choi, corrected segmentation fault, maybe memory leaking                                       
!
!---------------------------------------------------------------------------------------
SUBROUTINE EXCHANGE_GPU
    IMPLICIT NONE
    INTEGER :: VTYPE
    integer :: i,j
!TO DO by YUAN: switch back if necessary
!# if defined (CARTESIAN)
!    REAL(SP),DIMENSION(Mloc,Nloc),DEVICE :: rMASK
!# endif

# if defined (CARTESIAN)
    VTYPE=1
    CALL PHI_COLL_GPU(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Eta_d,VTYPE,PERIODIC)

    IF(VISCOSITY_BREAKING) THEN
        CALL PHI_COLL_GPU(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,AGE_BREAKING_d,VTYPE,PERIODIC)
    ENDIF
    
    !ykchoi (08.19.2015) :: new variable (WAVEMAKER_VIS) 	
    IF(VISCOSITY_BREAKING .OR. WAVEMAKER_VIS) THEN
        CALL PHI_COLL_GPU(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,nu_break_d,VTYPE,PERIODIC)
    ENDIF

    IF(VISCOSITY_BREAKING .OR. DIFFUSION_SPONGE .OR. WAVEMAKER_VIS)THEN
      CALL PHI_COLL_VARIABLE_LENGTH_GPU(Mloc1,Nloc,Ibeg,Iend1,Jbeg,Jend,Nghost,P_d,VTYPE)
      CALL PHI_COLL_VARIABLE_LENGTH_GPU(Mloc,Nloc1,Ibeg,Iend,Jbeg,Jend1,Nghost,Q_d,VTYPE)
    ENDIF
    !TO DO by YUAN, if necessary switch back
    ! add by YUAN. I dont know why there is a need for rmask
!    rMASK = MASK_d ! for periodic boundary condition
!    CALL PHI_COLL_GPU(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,rMASK,VTYPE,PERIODIC)  
!    MASK_d = rMASK 
    !CALL PHI_COLL_INT_GPU(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,MASK_d,VTYPE,PERIODIC)  
    CALL PHI_COLL_GPU(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,MASK_d,VTYPE,PERIODIC)  

    VTYPE=2
    CALL PHI_COLL_GPU(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,U_d,VTYPE,PERIODIC)
    CALL PHI_COLL_GPU(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,HU_d,VTYPE,PERIODIC)
    VTYPE=3
    CALL PHI_COLL_GPU(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,V_d,VTYPE,PERIODIC)
    CALL PHI_COLL_GPU(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,HV_d,VTYPE,PERIODIC)
# else
    VTYPE=1
    CALL PHI_COLL_GPU(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Eta_d,VTYPE)

    VTYPE=2
    CALL PHI_COLL_GPU(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,U_d,VTYPE)
    CALL PHI_COLL_GPU(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,HU_d,VTYPE)
    VTYPE=3
    CALL PHI_COLL_GPU(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,V_d,VTYPE)
    CALL PHI_COLL_GPU(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,HV_d,VTYPE)
# endif

! etaR x mask is a wrong idea
!    Eta=Eta*MASK
!$cuf kernel do(2) <<<*,*>>>
    do j = 1,Nloc
        do i = 1,Mloc
            U_d(i,j)=U_d(i,j)*MASK_d(i,j)
            V_d(i,j)=V_d(i,j)*MASK_d(i,j)
            HU_d(i,j)=HU_d(i,j)*MASK_d(i,j)
            HV_d(i,j)=HV_d(i,j)*MASK_d(i,j)
        enddo
    enddo

! TO DO by YUAN:  Refer to CPU version Subroutine EXCHANGE, the Pre-processing switch FILTERING part is not
! included here. 
    
END SUBROUTINE EXCHANGE_GPU

!-----------------------------------------------------------------------------------
!
!   PHI_COLL_VARIABLE_LENGTH subroutine is used to collect data into ghost cells                                                         
!
!   HISTORY:
!   07/09/2010 Fengyan Shi                                      
!
!-----------------------------------------------------------------------------------
SUBROUTINE PHI_COLL_VARIABLE_LENGTH_GPU(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,PHI,VTYPE)
    implicit none
    INTEGER,INTENT(IN) :: VTYPE
    INTEGER,INTENT(IN) :: Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost
    REAL(SP),DEVICE,INTENT(INOUT) :: PHI(Mloc,Nloc)
    INTEGER :: I,J,K

      ! x-direction
# if defined (MGPU)
    if ( n_west .eq. MPI_PROC_NULL ) then
# endif
    IF(WaveMaker(1:11)=='LEFT_BC_IRR')THEN
      ! do nothing
    ELSE
!$cuf kernel do(2) <<<*,(4,64)>>>
      DO J=Jbeg,Jend  
      DO K=1,Nghost
        PHI(K,J)=PHI(Ibeg+Nghost-K,J)
      ENDDO
      ENDDO
      !istat=cudaMemcpy2d(PHI(1,Jbeg),Mloc,PHI(Ibeg+2,Jbeg),Mloc,1,Nloc-2*Nghost)
      !istat=cudaMemcpy2d(PHI(2,Jbeg),Mloc,PHI(Ibeg+1,Jbeg),Mloc,1,Nloc-2*Nghost)
      !istat=cudaMemcpy2d(PHI(3,Jbeg),Mloc,PHI(Ibeg,Jbeg),Mloc,1,Nloc-2*Nghost)
    ENDIF
# if defined (MGPU)
    endif
# endif

# if defined (MGPU)
    if ( n_east .eq. MPI_PROC_NULL ) then
# endif
!$cuf kernel do(2) <<<*,(4,64)>>>
      DO J=Jbeg,Jend  
      DO K=1,Nghost
        PHI(Iend+K,J)=PHI(Iend-K+1,J)
      ENDDO
      ENDDO
      !istat=cudaMemcpy2d(PHI(Iend+1,Jbeg),Mloc,PHI(Iend,Jbeg),Mloc,1,Nloc-2*Nghost)
      !istat=cudaMemcpy2d(PHI(Iend+2,Jbeg),Mloc,PHI(Iend-1,Jbeg),Mloc,1,Nloc-2*Nghost)
      !istat=cudaMemcpy2d(PHI(Iend+3,Jbeg),Mloc,PHI(Iend-2,Jbeg),Mloc,1,Nloc-2*Nghost)
# if defined (MGPU)
    endif
# endif

      ! y-direction and corners
# if defined (MGPU)
    if ( n_suth .eq. MPI_PROC_NULL ) then
# endif
!$cuf kernel do(2) <<<*,(4,64)>>>
      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,K)=PHI(I,Jbeg+Nghost-K)
      ENDDO
      ENDDO
      !istat=cudaMemcpy2d(PHI(1,1),Mloc,PHI(1,Jbeg+2),Mloc,Mloc,1)
      !istat=cudaMemcpy2d(PHI(1,2),Mloc,PHI(1,Jbeg+1),Mloc,Mloc,1)
      !istat=cudaMemcpy2d(PHI(1,3),Mloc,PHI(1,Jbeg),Mloc,Mloc,1)
# if defined (MGPU)
    endif
# endif

# if defined (MGPU)
    if ( n_nrth .eq. MPI_PROC_NULL ) then
# endif
!$cuf kernel do(2) <<<*,(4,64)>>>
      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,Jend+K)=PHI(I,Jend-K+1)
      ENDDO
      ENDDO
      !istat=cudaMemcpy2d(PHI(1,Jend+1),Mloc,PHI(1,Jend),Mloc,Mloc,1)
      !istat=cudaMemcpy2d(PHI(1,Jend+2),Mloc,PHI(1,Jend-1),Mloc,Mloc,1)
      !istat=cudaMemcpy2d(PHI(1,Jend+3),Mloc,PHI(1,Jend-2),Mloc,Mloc,1)
# if defined (MGPU)
    endif
# endif

END SUBROUTINE PHI_COLL_VARIABLE_LENGTH_GPU

!-------------------------------------------------------------------------------------
!
!   PHI_COLL subroutine is used to collect data into ghost cells
!
!   HISTORY:
!     05/01/2010  Fengyan Shi
!     09/07/2010 Fengyan Shi, fix:
!       1) u v symmetric problem, 2) remove use global 3) fix bug
!     05/27/2010 Gangfeng Ma, corrected some bugs
!
!-------------------------------------------------------------------------------------
# if defined (CARTESIAN)
SUBROUTINE PHI_COLL_GPU(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,PHI,VTYPE,PERIODIC)
# else
SUBROUTINE PHI_COLL_GPU(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,PHI,VTYPE)
# endif
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: VTYPE
    INTEGER,INTENT(IN) :: Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost
    ! tested : the device variable can be used as dummy and actual variables
    REAL(SP),DEVICE,INTENT(INOUT) :: PHI(Mloc,Nloc)
    integer :: I,J,K
# if defined (CARTESIAN)
    LOGICAL :: PERIODIC
# endif
# if defined (MGPU)
    INTEGER :: len,II,JJ
    integer :: status(MPI_STATUS_SIZE)
    REAL(SP),DIMENSION(Mloc,Nghost) :: xx,send2d
# endif

! periodic first because it is not related to VTYPE
# if defined (CARTESIAN)
    IF(PERIODIC)THEN
# if defined (MGPU)
!   _____________________________ exchange
    IF (PY>1)THEN

      len=Mloc*Nghost
! south
    DO II = 1,PX
    ! by YUAN, EW-----1-------2------PX
    if(myid==ProcessorID(II,1)) then
! send from master
! By YUAN, for periodic bc, send from bottom to upper
        istat = cudaMemcpy2D(send2d(1,1),Mloc,PHI(1,1+Nghost),Mloc,Mloc,Nghost)
        call MPI_SENDRECV(send2d,len,MPI_SP,ProcessorID(II,PY),0,&
                         xx,len,MPI_SP,ProcessorID(II,PY),1,MPI_COMM_WORLD,status,ier)
        istat = cudaMemcpy2D(PHI(1,1),Mloc,xx(1,1),Mloc,Mloc,Nghost)
    endif ! end myid

    if(myid==ProcessorID(II,PY))then
        istat = cudaMemcpy2D(send2d(1,1),Mloc,PHI(1,Jend-Nghost+1),Mloc,Mloc,Nghost)
        call MPI_SENDRECV(send2d,len,MPI_SP,ProcessorID(II,1),1,&
                          xx,len,MPI_SP,ProcessorID(II,1),0,MPI_COMM_WORLD,status,ier)
        istat = cudaMemcpy2D(PHI(1,Jend+1),Mloc,xx(1,1),Mloc,Mloc,Nghost)
    endif

 
    ENDDO  ! end PX
! Add by YUAN for illustration
!-----------------------------
!ghost -- recv south exchange row
!----------------------------
!exchange rows
!----------------------------
!----------------------------
!----------------------------
!exchange rows
!----------------------------
!ghost -- recv north exchange rows
!----------------------------
    ELSE  ! PY = 1
!!$cuf kernel do(2) <<<*,(64,4)>>>
       istat = cudaMemcpy2D(PHI(Ibeg,1),Mloc,PHI(Ibeg,Jend-Nghost+1),Mloc,&
           Mloc-2*Nghost,Nghost)
       istat = cudaMemcpy2D(PHI(Ibeg,Jend+1),Mloc,PHI(Ibeg,Nghost+1),Mloc,&
           Mloc-2*Nghost,Nghost)
    ENDIF

!   ----------------------------- end exchange
# else 
        ! Serial
!$cuf kernel do(2) <<<*,(64,4)>>>
        DO K=1,Nghost
        DO I=1,Mloc
          PHI(I,K)=PHI(I,Jend-Nghost+K)
          PHI(I,Jend+K)=PHI(I,Nghost+K)
        ENDDO
        ENDDO        
# endif 
! end parallel
    ENDIF ! end periodic 

# endif 
! end cartesian

! I added coupling condition 10/14/2012
!   add left_bc wavemaker 09/13/2017

! for Eta
    IF(VTYPE==1) THEN  ! for eta
      ! x-direction
# if defined (MGPU)
    if ( n_west .eq. MPI_PROC_NULL ) then
# endif

# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_WEST)THEN
# endif

    IF(WaveMaker(1:11)=='LEFT_BC_IRR')THEN
      ! do nothing
    ELSE
!$cuf kernel do(2) <<<*,(4,64)>>>
      DO J=Jbeg,Jend  
      DO K=1,Nghost
        PHI(K,J)=PHI(Ibeg+Nghost-K,J)
      ENDDO
      ENDDO
    ENDIF  ! end left bc wavemaker

# if defined (COUPLING)
    ENDIF
# endif

# if defined (MGPU)
    endif
# endif

# if defined (MGPU)
    if ( n_east .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_EAST)THEN
# endif
!$cuf kernel do(2) <<<*,(4,64)>>>
      DO J=Jbeg,Jend  
      DO K=1,Nghost
        PHI(Iend+K,J)=PHI(Iend-K+1,J)
      ENDDO
      ENDDO
# if defined (COUPLING)
    ENDIF
# endif
# if defined (MGPU)
    endif
# endif

      ! y-direction and corners
# if defined (CARTESIAN)
      IF(.NOT.PERIODIC)THEN
# endif 
# if defined (MGPU)
    if ( n_suth .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_SOUTH)THEN
# endif
!$cuf kernel do(2) <<<*,(4,64)>>>
      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,K)=PHI(I,Jbeg+Nghost-K)
      ENDDO
      ENDDO
# if defined (COUPLING)
    ENDIF
# endif
# if defined (MGPU)
    endif
# endif

# if defined (MGPU)
    if ( n_nrth .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_NORTH)THEN
# endif
!$cuf kernel do(2) <<<*,(4,64)>>>
      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,Jend+K)=PHI(I,Jend-K+1)
      ENDDO
      ENDDO
# if defined (COUPLING)
    ENDIF
# endif
# if defined (MGPU)
    endif
# endif
# if defined (CARTESIAN)
      ENDIF
# endif

     ENDIF ! end vtype=1

! for u
    IF(VTYPE==2) THEN  ! for u (x-mirror condition)
      ! x-direction
# if defined (MGPU)
    if ( n_west .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_WEST)THEN
# endif

    IF(WaveMaker(1:11)=='LEFT_BC_IRR')THEN
      ! do nothing
    ELSE
!$cuf kernel do(2) <<<*,(4,64)>>>
      DO J=Jbeg,Jend
      DO K=1,Nghost
        PHI(K,J)=-PHI(Ibeg+Nghost-K,J)
      ENDDO
      ENDDO
    ENDIF ! end left_bc wavemaker
# if defined (COUPLING)
    ENDIF
# endif
# if defined (MGPU)
    endif
# endif

# if defined (MGPU)
    if ( n_east .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_EAST)THEN
# endif
!$cuf kernel do(2) <<<*,(4,64)>>>
      DO J=Jbeg,Jend
      DO K=1,Nghost
        PHI(Iend+K,J)=-PHI(Iend-K+1,J)
      ENDDO
      ENDDO
# if defined (COUPLING)
    ENDIF
# endif
# if defined (MGPU)
    endif
# endif

      ! y-direction and corners
# if defined (CARTESIAN)
      IF(.NOT.PERIODIC)THEN
# endif 

# if defined (MGPU)
    if ( n_suth .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_SOUTH)THEN
# endif
!$cuf kernel do(2) <<<*,(4,64)>>>
      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,K)=PHI(I,Jbeg+Nghost-K)
      ENDDO
      ENDDO
# if defined (COUPLING)
    ENDIF
# endif
# if defined (MGPU)
    endif
# endif

# if defined (MGPU)
    if ( n_nrth .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_NORTH)THEN
# endif
!$cuf kernel do(2) <<<*,(4,64)>>>
      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,Jend+K)=PHI(I,Jend-K+1)
      ENDDO
      ENDDO
# if defined (COUPLING)
    ENDIF
# endif
# if defined (MGPU)
    endif
# endif

# if defined (CARTESIAN)
      ENDIF
# endif

     ENDIF ! end vtype=2

    IF(VTYPE==3) THEN ! for v (y-mirror condition)
! for v
      ! x-direction

# if defined (MGPU)
    if ( n_west .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_WEST)THEN
# endif
    IF(WaveMaker(1:11)=='LEFT_BC_IRR')THEN
      ! do nothing
    ELSE
!$cuf kernel do(2) <<<*,(4,64)>>>
      DO J=Jbeg,Jend
      DO K=1,Nghost
        PHI(K,J)=PHI(Ibeg+Nghost-K,J)
      ENDDO
      ENDDO
    ENDIF ! end left_bc wavemaker
# if defined (COUPLING)
    ENDIF
# endif
# if defined (MGPU)
    endif
# endif

# if defined (MGPU)
    if ( n_east .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_EAST)THEN
# endif

!$cuf kernel do(2) <<<*,(4,64)>>>
      DO J=Jbeg,Jend
      DO K=1,Nghost
        PHI(Iend+K,J)=PHI(Iend-K+1,J)
      ENDDO
      ENDDO
# if defined (COUPLING)
    ENDIF
# endif
# if defined (MGPU)
    endif
# endif

      ! y-direction and corners
# if defined (CARTESIAN)
      IF(.NOT.PERIODIC)THEN
# endif 
  ! end cartesian
# if defined (MGPU)
    if ( n_suth .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_SOUTH)THEN
# endif
!$cuf kernel do(2) <<<*,(4,64)>>>
      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,K)=-PHI(I,Jbeg+Nghost-K)
      ENDDO
      ENDDO
# if defined (COUPLING)
    ENDIF
# endif
# if defined (MGPU)
    endif
# endif

# if defined (MGPU)
    if ( n_nrth .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_NORTH)THEN
# endif
!$cuf kernel do(2) <<<*,(4,64)>>>
      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,Jend+K)=-PHI(I,Jend-K+1)
      ENDDO
      ENDDO
# if defined (COUPLING)
    ENDIF
# endif
# if defined (MGPU)
    endif
# endif
# if defined (CARTESIAN)
      ENDIF
# endif

     ENDIF ! end vtype=3

! for cross-derivatives
    IF(VTYPE==4) THEN ! VTYPE==4 for u and v cross-mirror
     ! x-direction
# if defined (MGPU)
    if ( n_west .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_WEST)THEN
# endif
!$cuf kernel do(2) <<<*,(4,64)>>>
      DO J=Jbeg,Jend
      DO K=1,Nghost
        PHI(K,J)=0.0_SP
      ENDDO
      ENDDO
# if defined (COUPLING)
    ENDIF
# endif
# if defined (MGPU)
    endif
# endif

# if defined (MGPU)
    if ( n_east .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_EAST)THEN
# endif
!$cuf kernel do(2) <<<*,(4,64)>>>
      DO J=Jbeg,Jend
      DO K=1,Nghost
        PHI(Iend+K,J)=0.0_SP
      ENDDO
      ENDDO
# if defined (COUPLING)
    ENDIF
# endif
# if defined (MGPU)
    endif
# endif

      ! y-direction and corners, this one is not an exact solution
# if defined (MGPU)
    if ( n_suth .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_SOUTH)THEN
# endif
!$cuf kernel do(2) <<<*,(4,64)>>>
      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,K)=0.0_SP
      ENDDO
      ENDDO
# if defined (COUPLING)
    ENDIF
# endif
# if defined (MGPU)
    endif
# endif

# if defined (MGPU)
    if ( n_nrth .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_NORTH)THEN
# endif
!$cuf kernel do(2) <<<*,(4,64)>>>
      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,Jend+K)=0.0_SP
      ENDDO
      ENDDO
# if defined (COUPLING)
    ENDIF
# endif
# if defined (MGPU)
    endif
# endif

     ENDIF ! end vtype=4

! for symmetric
    IF(VTYPE==5)THEN
      ! x-direction
# if defined (MGPU)
    if ( n_west .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_WEST)THEN
# endif
!$cuf kernel do(2) <<<*,(4,64)>>>
      DO J=Jbeg,Jend
      DO K=1,Nghost
        PHI(K,J)=PHI(Ibeg+Nghost-K,J)
       ENDDO
      ENDDO
# if defined (COUPLING)
    ENDIF
# endif
# if defined (MGPU)
    endif
# endif

# if defined (MGPU)
    if ( n_east .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_EAST)THEN
# endif
!$cuf kernel do(2) <<<*,(4,64)>>>
      DO J=Jbeg,Jend
      DO K=1,Nghost
        PHI(Iend+K,J)=PHI(Iend-K+1,J)
      ENDDO
      ENDDO
# if defined (COUPLING)
    ENDIF
# endif
# if defined (MGPU)
    endif
# endif

      ! y-direction and corners

# if defined (MGPU)
    if ( n_suth .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_SOUTH)THEN
# endif
!$cuf kernel do(2) <<<*,(4,64)>>>
      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,K)=PHI(I,Jbeg+Nghost-K)
      ENDDO
      ENDDO
# if defined (COUPLING)
    ENDIF
# endif
# if defined (MGPU)
    endif
# endif

# if defined (MGPU)
    if ( n_nrth .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_NORTH)THEN
# endif
!$cuf kernel do(2) <<<*,(4,64)>>>
      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,Jend+K)=PHI(I,Jend-K+1)
      ENDDO
      ENDDO
# if defined (COUPLING)
    ENDIF
# endif
# if defined (MGPU)
    endif
# endif

     ENDIF ! end vtype=5

! for anti-symmetric
     IF(VTYPE==6)THEN
      ! x-direction
# if defined (MGPU)
    if ( n_west .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_WEST)THEN
# endif
!$cuf kernel do(2) <<<*,(4,64)>>>
      DO J=Jbeg,Jend
      DO K=1,Nghost
        PHI(K,J)=-PHI(Ibeg+Nghost-K,J)
      ENDDO
      ENDDO 
# if defined (COUPLING)
    ENDIF
# endif
# if defined (MGPU)
    endif
# endif

# if defined (MGPU)
    if ( n_east .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_EAST)THEN
# endif
!$cuf kernel do(2) <<<*,(4,64)>>>
      DO J=Jbeg,Jend
      DO K=1,Nghost
        PHI(Iend+K,J)=-PHI(Iend-K+1,J)
      ENDDO
      ENDDO 
# if defined (COUPLING)
    ENDIF
# endif
# if defined (MGPU)
    endif
# endif

      ! y-direction and corners
# if defined (MGPU)
    if ( n_suth .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_SOUTH)THEN
# endif
!$cuf kernel do(2) <<<*,(4,64)>>>
      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,K)=-PHI(I,Jbeg+Nghost-K)
      ENDDO
      ENDDO   
# if defined (COUPLING)
    ENDIF
# endif
# if defined (MGPU)
    endif
# endif

# if defined (MGPU)
    if ( n_nrth .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_NORTH)THEN
# endif
!$cuf kernel do(2) <<<*,(4,64)>>>
      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,Jend+K)=-PHI(I,Jend-K+1)
      ENDDO
      ENDDO     
# if defined (COUPLING)
    ENDIF
# endif
# if defined (MGPU)
    endif
# endif
    ENDIF ! end vtype=6

# if defined (MGPU)
    call phi_exch_cuda (PHI)
# endif

END SUBROUTINE PHI_COLL_GPU
!=================================================


# if defined (CARTESIAN)
SUBROUTINE periodic_NSexchange_gpu(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,PHI)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost
    REAL(SP),DEVICE,INTENT(INOUT) :: PHI(Mloc,Nloc)
    integer :: I,J,K
# if defined (MGPU)
    INTEGER :: len,II,JJ
    integer :: status(MPI_STATUS_SIZE)
    REAL(SP),DIMENSION(Mloc,Nghost) :: xx,send2d
# endif

! periodic first because it is not related to VTYPE
# if defined (MGPU)
!   _____________________________ exchange
    IF (PY>1)THEN

      len=Mloc*Nghost
! south
    DO II = 1,PX
    ! by YUAN, EW-----1-------2------PX
    if(myid==ProcessorID(II,1)) then
! send from master
! By YUAN, for periodic bc, send from bottom to upper
        istat = cudaMemcpy2D(send2d(1,1),Mloc,PHI(1,1+Nghost),Mloc,Mloc,Nghost)
        call MPI_SENDRECV(send2d,len,MPI_SP,ProcessorID(II,PY),0,&
                         xx,len,MPI_SP,ProcessorID(II,PY),1,MPI_COMM_WORLD,status,ier)
        istat = cudaMemcpy2D(PHI(1,1),Mloc,xx(1,1),Mloc,Mloc,Nghost)
    endif ! end myid

    if(myid==ProcessorID(II,PY))then
        istat = cudaMemcpy2D(send2d(1,1),Mloc,PHI(1,Jend-Nghost+1),Mloc,Mloc,Nghost)
        call MPI_SENDRECV(send2d,len,MPI_SP,ProcessorID(II,1),1,&
                          xx,len,MPI_SP,ProcessorID(II,1),0,MPI_COMM_WORLD,status,ier)
        istat = cudaMemcpy2D(PHI(1,Jend+1),Mloc,xx(1,1),Mloc,Mloc,Nghost)
    endif

 
    ENDDO  ! end PX
! Add by YUAN for illustration
!-----------------------------
!ghost -- recv south exchange row
!----------------------------
!exchange rows
!----------------------------
!----------------------------
!----------------------------
!exchange rows
!----------------------------
!ghost -- recv north exchange rows
!----------------------------
    ELSE  ! PY = 1
       istat = cudaMemcpy2D(PHI(Ibeg,1),Mloc,PHI(Ibeg,Jend-Nghost+1),Mloc,&
           Mloc-2*Nghost,Nghost)
       istat = cudaMemcpy2D(PHI(Ibeg,Jend+1),Mloc,PHI(Ibeg,Nghost+1),Mloc,&
           Mloc-2*Nghost,Nghost)
    ENDIF

!   ----------------------------- end exchange
# else 
       ! Serial
       istat = cudaMemcpy2D(PHI(1,1),Mloc,PHI(1,Jend-Nghost+1),Mloc,&
           Mloc,Nghost)
       istat = cudaMemcpy2D(PHI(1,Jend+1),Mloc,PHI(1,Nghost+1),Mloc,&
           Mloc,Nghost)
# endif 

end subroutine periodic_NSexchange_gpu
# endif


# if defined (CARTESIAN)
SUBROUTINE PHI_COLL_INT_GPU(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,PHI,VTYPE,PERIODIC)
# else
SUBROUTINE PHI_COLL_INT_GPU(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,PHI,VTYPE)
# endif
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: VTYPE
    INTEGER,INTENT(IN) :: Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost
    ! tested : the device variable can be used as dummy and actual variables
    INTEGER,DEVICE,INTENT(INOUT) :: PHI(Mloc,Nloc)
    integer :: I,J,K
# if defined (CARTESIAN)
    LOGICAL :: PERIODIC
# endif
# if defined (MGPU)
    INTEGER,DIMENSION(1) :: req
    INTEGER :: len,II,JJ
    integer :: status(MPI_STATUS_SIZE)
    INTEGER,DIMENSION(Mloc,Nghost) :: xx,send2d
# endif

! periodic first because it is not related to VTYPE
# if defined (CARTESIAN)
    IF(PERIODIC)THEN
# if defined (MGPU)
!   _____________________________ exchange
    IF (PY>1)THEN

      len=Mloc*Nghost
! south
    DO II = 1,PX

    ! by YUAN, EW-----1-------2------PX
    if(myid==ProcessorID(II,1)) then
! send from master
! By YUAN, for periodic bc, send from bottom to upper
        istat = cudaMemcpy2D(send2d(1,1),Mloc,PHI(1,1+Nghost),Mloc,Mloc,Nghost)
        call MPI_SENDRECV(send2d,len,MPI_INTEGER,ProcessorID(II,PY),0,&
                         xx,len,MPI_INTEGER,ProcessorID(II,PY),1,MPI_COMM_WORLD,status,ier)
        istat = cudaMemcpy2D(PHI(1,1),Mloc,xx(1,1),Mloc,Mloc,Nghost)
    endif ! end myid

    if(myid==ProcessorID(II,PY))then
        istat = cudaMemcpy2D(send2d(1,1),Mloc,PHI(1,Jend-Nghost+1),Mloc,Mloc,Nghost)
        call MPI_SENDRECV(send2d,len,MPI_INTEGER,ProcessorID(II,1),1,&
                          xx,len,MPI_INTEGER,ProcessorID(II,1),0,MPI_COMM_WORLD,status,ier)
        istat = cudaMemcpy2D(PHI(1,Jend+1),Mloc,xx(1,1),Mloc,Mloc,Nghost)
    endif

    ENDDO  ! end PX
! Add by YUAN for illustration
!-----------------------------
!ghost -- recv south exchange row
!----------------------------
!exchange rows
!----------------------------
!----------------------------
!----------------------------
!exchange rows
!----------------------------
!ghost -- recv north exchange rows
!----------------------------
    ELSE  ! PY = 1
       istat = cudaMemcpy2D(PHI(Ibeg,1),Mloc,PHI(Ibeg,Jend-Nghost+1),Mloc,&
           Mloc-2*Nghost,Nghost,cudaMemcpyDeviceToDevice)
       istat = cudaMemcpy2D(PHI(Ibeg,Jend+1),Mloc,PHI(Ibeg,Nghost+1),Mloc,&
           Mloc-2*Nghost,Nghost,cudaMemcpyDeviceToDevice)
    ENDIF

!   ----------------------------- end exchange
# else 
        ! Serial
!$cuf kernel do(2) <<<*,(64,4)>>>
        DO K=1,Nghost
        DO I=1,Mloc
          PHI(I,K)=PHI(I,Jend-Nghost+K)
          PHI(I,Jend+K)=PHI(I,Nghost+K)
        ENDDO
        ENDDO        
# endif 
! end parallel
    ENDIF ! end periodic 

# endif 
! end cartesian

    IF(VTYPE==1) THEN  ! for mask, which is integer
      ! x-direction
# if defined (MGPU)
    if ( n_west .eq. MPI_PROC_NULL ) then
# endif

# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_WEST)THEN
# endif

    IF(WaveMaker(1:11)=='LEFT_BC_IRR')THEN
      ! do nothing
    ELSE
!$cuf kernel do(2) <<<*,(4,64)>>>
      DO J=Jbeg,Jend  
      DO K=1,Nghost
        PHI(K,J)=PHI(Ibeg+Nghost-K,J)
      ENDDO
      ENDDO
    ENDIF  ! end left bc wavemaker

# if defined (COUPLING)
    ENDIF
# endif

# if defined (MGPU)
    endif
# endif

# if defined (MGPU)
    if ( n_east .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_EAST)THEN
# endif
!$cuf kernel do(2) <<<*,(4,64)>>>
      DO J=Jbeg,Jend  
      DO K=1,Nghost
        PHI(Iend+K,J)=PHI(Iend-K+1,J)
      ENDDO
      ENDDO
# if defined (COUPLING)
    ENDIF
# endif
# if defined (MGPU)
    endif
# endif

      ! y-direction and corners
# if defined (CARTESIAN)
      IF(.NOT.PERIODIC)THEN
# endif 
# if defined (MGPU)
    if ( n_suth .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_SOUTH)THEN
# endif
!$cuf kernel do(2) <<<*,(4,64)>>>
      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,K)=PHI(I,Jbeg+Nghost-K)
      ENDDO
      ENDDO
# if defined (COUPLING)
    ENDIF
# endif
# if defined (MGPU)
    endif
# endif

# if defined (MGPU)
    if ( n_nrth .eq. MPI_PROC_NULL ) then
# endif
# if defined (COUPLING)
    IF(.NOT.IN_DOMAIN_NORTH)THEN
# endif
!$cuf kernel do(2) <<<*,(4,64)>>>
      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,Jend+K)=PHI(I,Jend-K+1)
      ENDDO
      ENDDO
# if defined (COUPLING)
    ENDIF
# endif
# if defined (MGPU)
    endif
# endif
# if defined (CARTESIAN)
      ENDIF
# endif

     ENDIF ! end vtype=1

# if defined (MGPU)
    call phi_exch_cuda (PHI)
# endif

END SUBROUTINE PHI_COLL_INT_GPU
end module boundary_condition_module
