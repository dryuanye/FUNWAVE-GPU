!-------------------------------------------------------------------------------------
!
!    Mod_cuda is module to declare, initialize the device variables for GPU version
!
!    HISTORY:
!      04/01/2019  created
!
!-----------------------------------------------------------------------------------
module mod_cuda
    use cudafor
    use cusparse
    use PARAM
    use GLOBAL
    IMPLICIT NONE
    SAVE
!
    integer,parameter :: TriSolver = 0    ! 0, 1-cusparse
    type ( cudaEvent ) :: startEvent , stopEvent
    !pinned variables
    integer(kind=cuda_stream_kind),dimension(8), save :: streamID
    integer(kind=cuda_stream_kind),save :: stream_d2h
    integer :: istat
    real(SP),save :: TotalKernelTime
    real(SP) :: KernelTime
    real(SP) :: t1, t2
    type(cusparseHandle),save :: cusparseh
    integer, parameter :: BlockDimX_2D = 16, BlockDimY_2D = 16
    integer, parameter :: BlockDimX_Inner_2D = 14, BlockDimY_Inner_2D = 14
!    integer, parameter :: BlockDimX_2D = 8, BlockDimY_2D = 8
!    integer, parameter :: BlockDimX_Inner_2D = 6, BlockDimY_Inner_2D = 6
    type(dim3) :: grid, tBlock,grid_step1,grid_step2,tBlock_tran,&
        grid_tran,grid_tran_back
    integer(cuda_count_kind):: freeMem,totalMem
    real :: rFreeMem,rUsedMem,rTotalMem,maxUsedMem,minFreeMem
! constant variables ( Refer to mod_param )
    REAL(SP), constant :: n_left_d=-1.0_SP,n_right_d=1.0_SP,n_bottom_d=-1.0_SP,n_top_d=1.0_SP
    real(SP), constant :: Beta_ref_d
    !parameters for MUSCL numerical schemes
    real(SP), dimension(3), constant :: alpha_d=(/0.0_SP,3.0_SP/4.0_SP,1.0_SP/3.0_SP/)
    real(SP), dimension(3), constant :: beta_d=(/1.0_SP,1.0_SP/4.0_SP,2.0_SP/3.0_SP/)
! station data
! timer
# if defined (CARTESIAN)
      real(SP), device :: DX_d, DY_d
!# else
!      REAL(SP), device :: Lon_West_d,Lat_South_d,Dphi_d,Dtheta_d
# endif
!    real(SP), DT_d, DT_fixed_d
    real(SP), dimension(:,:), allocatable, device :: &
            BreakWaterWidth_d, CD_breakwater_d
!    real(SP), device :: BreakWaterAbsorbCoef_d
! coordinates for Cartesian only
    real(SP), dimension(:), allocatable, device :: Xco_d, Yco_d
! device variables for Riemann Problems
    real(SP),dimension(:,:), allocatable, device :: &
# if defined (CARTESIAN)
        U4xL_d,U4xR_d,V4yL_d,V4yR_d, DelxU4_d,DelxV4_d,DelyV4_d,DelyU4_d,&
	V4xL_d,V4xR_d,U4yL_d,U4yR_d, &   
        U4_d,V4_d,U1p_d,V1p_d, &
        U1pp_d,V1pp_d, &
        U2_d,V2_d,U3_d,V3_d, &
# else
# if defined(ZALPHA)
        U4xL_d,U4xR_d,V4yL_d,V4yR_d, DelxU4_d,DelxV4_d,DelyV4_d,DelyU4_d,&
	V4xL_d,V4xR_d,U4yL_d,U4yR_d, &
        U4_d,V4_d,&
# endif
        Dx_d,Dy_d, Coriolis_d,Lat_theta_d, &
        U1p_d,V1p_d, SlopeX_d,SlopeY_d,&
# endif
        DelxU_d,DelxHU_d,DelxV_d,DelxEtar_d,&
        DelxHV_d, DelyHU_d, &
        DelyU_d,DelyHV_d,DelyV_d,DelyEtar_d,&
        UxL_d,UxR_d,VxL_d,VxR_d,&
        HUxL_d,HUxR_d,HUyL_d,HUyR_d,HxL_d,HxR_d, &
        EtaRxL_d,EtaRxR_d,&
        UyL_d,UyR_d,VyL_d,VyR_d,&
        HVxL_d,HVxR_d,HVyL_d,HVyR_d,HyL_d,HyR_d, &
        EtaRyL_d,EtaRyR_d, &
        PL_d,PR_d,QL_d,QR_d, &
        FxL_d,FxR_d,FyL_d,FyR_d, &
        GxL_d,GxR_d,GyL_d,GyR_d, &
        SxL_d,SxR_d,SyL_d,SyR_d, &
! cross-derivatives 
        Vxy_d,DVxy_d,Uxy_d,DUxy_d, &
! second-derivatives
        Uxx_d,DUxx_d,Vyy_d,DVyy_d, &
! first-derivatives
        Ux_d,Vx_d,Uy_d,Vy_d,DUx_d,DUy_d,DVx_d,DVy_d, &
        ETAx_d,ETAy_d, ETAT_d, ETATx_d,ETATy_d, &
! time-derivative
        U0_d,V0_d,Ut_d,Vt_d,Utx_d,Vty_d,Utxx_d,Utxy_d,Vtxy_d,Vtyy_d,&
        DUtxx_d,DUtxy_d,DVtxy_d,DVtyy_d,DUtx_d,DVty_d,&
! original variables
        Fx_d,Fy_d,U_d,V_d,HU_d,HV_d,&
        Gx_d,Gy_d,P_d,Q_d,SourceX_d,SourceY_d,Int2Flo_d, &
        tmp4preview_d,HeightMax_d,HeightMin_d,VelocityMax_d,&
        MomentumFluxMax_d,VorticityMax_d
# if defined (ITERATION)
    real(SP),dimension(:,:),allocatable, device :: Ui_d,Vi_d,ETAi_d,&
                        UbarOld_d,VbarOld_d,EtaOld_d
!       INTEGER, device :: Ki_d
!       LOGICAL, device :: ACCURATE_d = .FALSE.
# endif

! wetting and drying
        INTEGER,DIMENSION(:,:),ALLOCATABLE,device :: &
                MASK_d,MASK_STRUC_d,MASK9_d, MASKtmp_d
        integer,ALLOCATABLE,DEVICE :: MASK9u_d(:,:), MASK9v_d(:,:)
!        REAL(SP), device :: Dmass_d,WetArea_d,DwetEta_d
! wave maker
! See init codes for WaveMakerCode
        integer :: WaveMakerCode
        REAL(SP),DIMENSION(:),ALLOCATABLE,device :: omgn_ir_d
        REAL(SP),DIMENSION(:),ALLOCATABLE,device:: omgn2D_d
        REAL(SP),DIMENSION(:,:),ALLOCATABLE,device :: Wavemaker_Mass_d
        REAL(SP),DIMENSION(:,:,:),ALLOCATABLE,device ::Cm_d,Sm_d
        INTEGER, device :: Nfreq_d,Ntheta_d 
!      for measure time series
       REAL(SP),DIMENSION(:,:),ALLOCATABLE,device :: WAVE_COMP_d
       REAL(SP),DIMENSION(:),ALLOCATABLE,device :: Beta_genS_d,D_genS_d
        REAL(SP),DIMENSION(:,:),ALLOCATABLE,device :: Cd_d

! sponge
        REAL(SP),DIMENSION(:,:),ALLOCATABLE,device :: SPONGE_d,SpongeMaker_d

! smagorinsky and wave height
      REAL(SP),DIMENSION(:,:),ALLOCATABLE,device :: Umean_d,Vmean_d,&
                  ETAmean_d,Usum_d,Vsum_d,ETAsum_d, nu_smg_d
!      REAL(SP),constant::T_INTV_mean_d = 20.0,T_sum_d=0.0,C_smg_d=0.25
      REAL(SP),DIMENSION(:,:),ALLOCATABLE, device :: &
                WaveHeightRMS_d,WaveHeightAve_d,Emax_d,Emin_d,& 
                HrmsSum_d,HavgSum_d
      INTEGER, DIMENSION(:,:),ALLOCATABLE,device :: Num_Zero_Up_d

      !ykchoi
      REAL(SP),DIMENSION(:,:),ALLOCATABLE,device ::&
              ETA2sum_d,ETA2mean_d,SigWaveHeight_d

! depth H=Eta+Depth, 
      REAL(SP),DIMENSION(:,:),ALLOCATABLE,device :: Depth_d,H_d,&
            DepthNode_d,Depthx_d,Depthy_d
!      REAL(SP),device::Depth_Flat_d, SLP_d,Xslp_d

! updating variables
      REAL(SP),DIMENSION(:,:),ALLOCATABLE,device:: &
              Ubar0_d,Vbar0_d,Eta0_d,Ubar_d,Vbar_d,Eta_d,ETAtmp_d


       REAL(SP),DIMENSION(:,:),ALLOCATABLE,device :: &
               Amp_Ser_d,Phase_LEFT_d
       REAL(SP),DIMENSION(:),ALLOCATABLE,device :: &
               Per_Ser_d,Phase_Ser_d,Theta_Ser_d
       REAL(SP),DIMENSION(:),ALLOCATABLE,device :: &
               Segma_Ser_d,Wave_Number_Ser_d
       REAL(SP),DIMENSION(:,:,:),ALLOCATABLE,device :: &
               Cm_eta_d,Sm_eta_d, Cm_u_d,Sm_u_d,Cm_v_d,Sm_v_d

! eddy viscosity breaking
        REAL(SP),DIMENSION(:,:),ALLOCATABLE,device :: AGE_BREAKING_d, AGE_BREAKING0_d

       REAL(SP),DIMENSION(:,:),ALLOCATABLE,device :: nu_break_d,nu_sponge_d
       REAL(SP),DIMENSION(:,:),ALLOCATABLE,device :: CD_4_SPONGE_d

contains



subroutine var_d2h_output
    implicit none
    ! first download, then use the following MAXVAL...
    ! otherwise, need parallel reduction to yield Max Vals.
    ! Cd_breakwater and Depth should be invariable after init.F. check later
    if (icount == 1) then
        if (out_depth .or. breakwater) then
            DEPTH = DEPTH_d
            !CD_breakwater = CD_breakwater_d
        endif
    endif

    IF(OUT_HMAX) HeightMax = HeightMax_d 

!    istat = cudaMemcpy2dAsync(Eta(1,1),Mloc,Eta_d(1,1),Mloc,Mloc,Nloc,stream=stream_d2h)
!    istat = cudaMemcpy2dAsync(U(1,1),Mloc,U_d(1,1),Mloc,Mloc,Nloc,stream=stream_d2h)
!    istat = cudaMemcpy2dAsync(V(1,1),Mloc,V_d(1,1),Mloc,Mloc,Nloc,stream=stream_d2h)
!    istat = cudaMemcpy2dAsync(H(1,1),Mloc,H_d(1,1),Mloc,Mloc,Nloc,stream=stream_d2h)
!    istat = cudaMemcpy2dAsync(MASK(1,1),Mloc,MASK_d(1,1),Mloc,Mloc,Nloc,stream=stream_d2h)
    Eta = Eta_d
    U = U_d
    V = V_d
    H = H_d  
    MASK = MASK_d
    IF(OUT_MASK9)THEN
       MASK9 = MASK9_d
    ENDIF
    IF(OUT_P)THEN
        P = P_d
    ENDIF
    IF(OUT_Q)THEN
        Q = Q_d
    ENDIF
    IF(OUT_VORmax)THEN
        VorticityMax=VorticityMax_d
    ENDIF
    IF(OUT_AGE)THEN
      IF(SHOW_BREAKING)THEN
         AGE_BREAKING = AGE_BREAKING_d
      ENDIF
    ENDIF
    IF(VISCOSITY_BREAKING)THEN
       IF(OUT_NU)THEN
          nu_break = nu_break_d
       ENDIF
    ENDIF
    IF(OUT_TMP)THEN
       tmp4preview = tmp4preview_d
    ENDIF
end subroutine var_d2h_output




SUBROUTINE ESTIMATE_DT_GPU(DT,TIME)
    IMPLICIT NONE
# if defined (MGPU)
    REAL(SP) :: myvar
# endif
    REAL(SP),INTENT(OUT)::DT
    REAL(SP),INTENT(INOUT)::TIME
    real(SP),device :: CelerityTmp,DxDcele,DyDcele
    real(SP) :: tMinX, tMinY, DT_tmp
    real(SP) :: TMP3

    tMinX=LARGE
    tMinY=LARGE
!$cuf kernel do(2) <<<*,*>>>
    DO J=1,Nloc
    DO I=1,Mloc
! x direction
!        CelerityTmp=ABS(U_d(I,J))+SQRT(GRAV*MAX(H_d(I,J),MinDepthFrc))
# if defined (CARTESIAN)
        DxDcele=DX_d/max(SMALL,ABS(U_d(I,J))+SQRT(GRAV*MAX(H_d(I,J),MinDepthFrc)))
# else
        DxDcele=DX_d(I,J)/max(SMALL,ABS(U_d(I,J))+SQRT(GRAV*MAX(H_d(I,J),MinDepthFrc)))
# endif
        tMinX = min(tMinX, DxDcele)
    ENDDO
    ENDDO
!$cuf kernel do(2) <<<*,*>>>
    DO J=1,Nloc
    DO I=1,Mloc
! y direction
# if defined (CARTESIAN)
        DyDcele=DY_d/max(SMALL,ABS(V_d(I,J))+SQRT(GRAV*MAX(H_d(I,J),MinDepthFrc)))
# else
        DyDcele=DY_d(I,J)/max(SMALL,ABS(V_d(I,J))+SQRT(GRAV*MAX(H_d(I,J),MinDepthFrc)))
# endif
        tMinY = min(tMinY, DyDcele)
    ENDDO
    ENDDO
    TMP3 = min(tMinX,tMinY)
# if defined (MGPU)
     call MPI_ALLREDUCE (TMP3,myvar,1,MPI_SP,MPI_MIN,&
          MPI_COMM_WORLD,ier)
     TMP3 = myvar
# endif
    DT_tmp=CFL*TMP3

    IF(FIXED_DT)THEN
        DT = DT_fixed
        DO WHILE (DT > DT_tmp)
            DT=DT/2.0_SP
        ENDDO
    ELSE
        DT = DT_tmp
    ENDIF

     TIME=TIME+DT

END SUBROUTINE ESTIMATE_DT_GPU

! type ( cudaEvent ) :: startEvent , stopEvent
! real :: time
! integer :: istat
!
! istat = cudaEventCreate ( startEvent )
! istat = cudaEventCreate ( stopEvent )
!
! istat = cudaEventRecord ( startEvent , 0)
! call increment <<<grid ,tBlock >>>( a_d , b )
! istat = cudaEventRecord ( stopEvent , 0)
! istat = cudaEventSynchronize ( stopEvent )
! istat = cudaEventElapsedTime ( time , startEvent , stopEvent )
!
!     write (* ,*) ' Time for kernel execution ( ms ): ' , time
! istat = cudaEventDestroy ( startEvent )
! istat = cudaEventDestroy ( stopEvent )

subroutine cudadeviceQuery
    implicit none
    type ( cudaDeviceProp ) :: prop
    integer :: nDevices , i , ierr
    character(len=40) :: screenout

! Number of CUDA - capable devices
      print *, "GPU Query---GPU Query---GPU Query---GPU Query---GPU Query---"
      ierr = cudaGetDeviceCount ( nDevices )
      if ( nDevices == 0) then
          print *, "/ No CUDA devices found /"
          stop
      else if ( nDevices == 1) then
          print *, "/  One CUDA device found /"
      else
          write(*,99) nDevices
99        format( '/ ', i0, ' CUDA devices found /')
      end if
! Loop over devices
      do i = 0 , nDevices -1
          write (* ,100) i
100       format( ' Device Number : ', i0 )
          ierr = cudaGetDeviceProperties ( prop , i )
! General device info
          write (* , 101)  trim ( prop % name )
101       format( ' Device Name : ',a )
          write (* , 102) prop % major , prop % minor
102       format( ' Compute Capability : ',i0 , '. ' , i0 )
          write (*, 103) prop % multiProcessorCount
103       format( ' Number of Multiprocessors : ', i0 )
          write (* , 104) prop % clockRate
104       format( ' Max Clock Rate ( kHz ): ',i0 )

! Execution Configuration
          write ( *, *) 'Execution Configureation Limits'
          write (* , 105) prop % maxGridSize
105       format( ' Max Grid Dims : ' ,2( i0 , ' x ') , i0 )
          write (* , 106) prop % maxThreadsDim
106       format(' Max Block Dims : ' ,2( i0 , ' x ') , i0 )
          write (* , 107)  prop % maxThreadsPerBlock
107       format( ' Max Threads per Block : ',i0 )
    print *, "GPU Query---GPU Query---GPU Query---GPU Query---GPU Query---"

      enddo

end subroutine cudadeviceQuery

subroutine cusparse_init
    implicit none
! TRIDIAGONAL SOLVER CUSPARSE INIT
            if (TriSolver==1) then
                istat = cusparseCreate(cusparseh)
                if(istat /= CUSPARSE_STATUS_SUCCESS) &
                    write(*,*) 'Cusparse Create failed', istat      
            endif
end subroutine cusparse_init

subroutine cusparse_destroy
    implicit none
            if (TriSolver==1) then
                istat = cusparseDestroy(cusparseh)
                if(istat /= CUSPARSE_STATUS_SUCCESS) &
                    write(*,*) 'Cusparse Destroy failed', istat        
            endif
end subroutine cusparse_destroy

subroutine KernelTimePrint
    implicit none

    print *, 'Tridiagonal solver time consumption : ',TotalKernelTime, ' s.'
end subroutine KernelTimePrint


subroutine check_print(varH,varD)
    implicit none
    real(SP), dimension(Mloc,Nloc) :: varH,varD2H,varDiff
    real(SP),dimension(Mloc,Nloc),device :: varD
    real(SP) :: varMaxDiff
    integer :: i,j
    varD2H = varD
    varDiff = varH-varD2H
    varMaxDiff = maxval(abs(varDiff))
    if (varMaxDiff>0.000001) then
        do i = 1,Mloc
        print *, varH(i, 100), varD2H(i,100),varH(i,200),varD2H(i,200)
        enddo
!        open(1224, file='check_print.txt', status='new')
!        do j = 1,Nloc
!            write(1224, *) (varDiff(i,j), i=1,Mloc)
!        enddo
!        close(1224)
        print *, 'Max Difference between H and D is :', varMaxDiff
        stop
    else
        print *, 'Max Difference between H and D is smaller:', varMaxDiff
    endif
end subroutine check_print

subroutine check_print_int(varH,varD)
    implicit none
    integer, dimension(Mloc,Nloc) :: varH,varD2H,varDiff
    integer,dimension(Mloc,Nloc),device :: varD
    integer :: varMaxDiff
    integer :: i,j
    varD2H = varD
    varDiff = varH-varD2H
    varMaxDiff = maxval(abs(varDiff))
    if (real(varMaxDiff)>0.000001) then
        do i = 1,Mloc
        print *, varH(i, 100), varD2H(i,100),varH(i,200),varD2H(i,200)
        enddo
!        open(1224, file='check_print.txt', status='new')
!        do j = 1,Nloc
!            write(1224, *) (varDiff(i,j), i=1,Mloc)
!        enddo
!        close(1224)
        print *, 'Max Difference between H and D is :', varMaxDiff
        stop
    else
        print *, 'Max Difference between H and D is smaller:', varMaxDiff
    endif
end subroutine check_print_int

subroutine check_outputD(varD)
    implicit none
    real(SP), dimension(Mloc,Nloc) :: varD2H
    real(SP),dimension(Mloc,Nloc),device :: varD
    integer :: i,j
    varD2H = varD
        open(1224, file='check_outputD.txt')
        do j = 1,Nloc
            write(1224, *) (varD2H(i,j), i=1,Mloc)
        enddo
        close(1224)
        print *, 'Output Device var for check'
end subroutine check_outputD

subroutine check_outputH(varH)
    implicit none
    real(SP), dimension(Mloc,Nloc) :: varH
    integer :: i,j
        open(1225, file='check_outputH.txt')
        do j = 1,Nloc
            write(1225, *) (varH(i,j), i=1,Mloc)
        enddo
        close(1225)
        print *, 'Output Host var for check'
end subroutine check_outputH


end module mod_cuda
