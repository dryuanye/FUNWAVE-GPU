!-----------------------------------------------------------------------------------
module init_gpu
    use cudafor
    use PARAM
    use GLOBAL
    use mod_cuda
# if defined (MGPU)
    use mgpu_utilities
# endif
    use DispersionFluxesSources
    IMPLICIT NONE
    SAVE
!
contains

!-------------------------------------------------------------------------------------
!
!    AssignDevice is subroutine to assign GPU devices to MPI ranks.
!
!  HISTORY:
!  07/15/2019  Created
!
!-------------------------------------------------------------------------------------
# if defined (MGPU)
subroutine AssignDevice
    use GLOBAL, only: ier,myid,nprocs 
    implicit none
    integer ::  istat, newComm, newRank, dev
    character (len=MPI_MAX_PROCESSOR_NAME) :: hostname
    integer :: namelength
    integer :: i
!
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ier)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ier)
    !For OPENMPI, version 3.1.3 above have MPI_COMM_SPLIT_TYPE
    call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, &
                         MPI_INFO_NULL, newComm, ier)
    call MPI_COMM_RANK(newComm, newRank, ier)
    call MPI_GET_PROCESSOR_NAME(hostname, namelength, ier)
    dev = newRank
    istat = cudaSetDevice(dev)
    !
    OPEN(1227, FILE='GPULOG.txt')
    if (myid.eq.0) WRITE(1227,*)'--------------GPU LOG FILE -----------------'
    if (myid.eq.0) WRITE(1227,*)'--------------- GPU DEVICE --------------'
    do i=0,nprocs-1
      if(myid == i) then
          write(*,"(A8,I4,A13,I2,A8,A20)") "RANK: ",myid," USING GPU: ",dev, " AT: ", hostname(1:namelength)
          write(1227,"(A8,I4,A13,I2,A8,A20)") "RANK: ",myid," USING GPU: ",dev, " AT: ", hostname(1:namelength)
      endif
    end do
end subroutine AssignDevice
#endif

subroutine CudaMemUse
    USE GLOBAL, ONLY: IER,MYID
    USE MOD_CUDA, ONLY: FREEMEM,TOTALMEM,RFREEMEM,RUSEDMEM,RTOTALMEM,MAXUSEDMEM,MINFREEMEM,ISTAT
    implicit none
! Memory usage in GPU device
    ISTAT=CUDAMEMGETINFO(FREEMEM,TOTALMEM)
    RTOTALMEM = TOTALMEM/(1024.**2)
    RFREEMEM =   FREEMEM/(1024.**2)
    RUSEDMEM = (TOTALMEM-FREEMEM)/(1024**2)
# if defined (MGPU)
    CALL MPI_ALLREDUCE(RUSEDMEM,MAXUSEDMEM,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,IER)
    CALL MPI_ALLREDUCE(RFREEMEM,MINFREEMEM,1,MPI_REAL,MPI_MIN,MPI_COMM_WORLD,IER)
    if (myid==0) then
        WRITE(*,"(A20,F7.1,A3,F7.1,A3,F7.1,A8)") " MAX GPU MEMORY USED: ",MAXUSEDMEM," / ",RTOTALMEM," MBYTES"
        WRITE(*,"(A20,F7.1,A3,F7.1,A3,F7.1,A8)") " MIN GPU MEMORY FREE: ",MINFREEMEM," / ",RTOTALMEM," MBYTES"
        WRITE(1227,"(A20,F7.1,A3,F7.1,A3,F7.1,A8)") " MAX GPU MEMORY USED: ",MAXUSEDMEM," / ",RTOTALMEM," MBYTES"
        WRITE(1227,"(A20,F7.1,A3,F7.1,A3,F7.1,A8)") " MIN GPU MEMORY FREE: ",MINFREEMEM," / ",RTOTALMEM," MBYTES"
        IF (MAXUSEDMEM > 0.8*RTOTALMEM) THEN
            PRINT *, 'WARNING: GPU GLOBAL MEMORY NOT ENOUGH'
            WRITE(1227,"(A30)") 'WARNING: GPU GLOBAL MEMORY NOT ENOUGH'
            STOP
        ENDIF
    endif
# else
    WRITE(*,"(A20,F7.1,A3,F7.1,A3,F7.1,A8)") " MAX GPU MEMORY USED: ",RUSEDMEM," / ",RTOTALMEM," MBYTES"
    WRITE(*,"(A20,F7.1,A3,F7.1,A3,F7.1,A8)") " MIN GPU MEMORY FREE: ",RFREEMEM," / ",RTOTALMEM," MBYTES"
    WRITE(1227,"(A20,F7.1,A3,F7.1,A3,F7.1,A8)") " MAX GPU MEMORY USED: ",RUSEDMEM," / ",RTOTALMEM," MBYTES"
    WRITE(1227,"(A20,F7.1,A3,F7.1,A3,F7.1,A8)") " MIN GPU MEMORY FREE: ",RFREEMEM," / ",RTOTALMEM," MBYTES"
    IF (RUSEDMEM > 0.8*RTOTALMEM) THEN
            PRINT *, 'WARNING: GPU GLOBAL MEMORY NOT ENOUGH'
            WRITE(1227,"(A30)") 'WARNING: GPU GLOBAL MEMORY NOT ENOUGH'
            STOP
    ENDIF
# endif

end subroutine CudaMemUse

subroutine stream_init
    implicit none
    integer :: i
    istat = cudaStreamCreate( stream_d2h  )
    do i = 1,8
        istat = cudaStreamCreate( streamID(i)  )
    enddo

end subroutine stream_init

subroutine stream_destroy
    implicit none
    integer :: i
    istat = cudaStreamDestroy(stream_d2h)
    do i = 1,8
        istat = cudaStreamDestroy( streamID(i)  )
    enddo
end subroutine stream_destroy


subroutine ALLOCATE_VARIABLES_GPU
    implicit none
! coordinate for cartesian only
    allocate (Xco_d(Mloc),Yco_d(Nloc))
! allocate variables
    allocate (DelxU_d(Mloc,Nloc),DelxHU_d(Mloc,Nloc),DelxV_d(Mloc,Nloc),DelxEtar_d(Mloc,Nloc),&
              DelyU_d(Mloc,Nloc),DelyHV_d(Mloc,Nloc),DelyV_d(Mloc,Nloc),DelyEtar_d(Mloc,Nloc),&
              DelxHV_d(Mloc,Nloc),DelyHU_d(Mloc,Nloc), &
! U V HU H in x-direction
              UxL_d(Mloc1,Nloc),UxR_d(Mloc1,Nloc),VxL_d(Mloc1,Nloc),VxR_d(Mloc1,Nloc),&
              HUxL_d(Mloc1,Nloc),HUxR_d(Mloc1,Nloc),HVxL_d(Mloc1,Nloc),HVxR_d(Mloc1,Nloc), &
              HxL_d(Mloc1,Nloc),HxR_d(Mloc1,Nloc), &
! U V HV H in y-direction
              UyL_d(Mloc,Nloc1),UyR_d(Mloc,Nloc1),VyL_d(Mloc,Nloc1),VyR_d(Mloc,Nloc1),&
              HVyL_d(Mloc,Nloc1),HVyR_d(Mloc,Nloc1),HUyL_d(Mloc,Nloc1),HUyR_d(Mloc,Nloc1), &
              HyL_d(Mloc,Nloc1),HyR_d(Mloc,Nloc1), &
! cross-derivatives
              Uxy_d(Mloc,Nloc),Vxy_d(Mloc,Nloc),DUxy_d(Mloc,Nloc),DVxy_d(Mloc,Nloc), &
! second-derivatives
              Uxx_d(Mloc,Nloc),Vyy_d(Mloc,Nloc),DUxx_d(Mloc,Nloc),DVyy_d(Mloc,Nloc), &
! 1st-derivatives
              Ux_d(Mloc,Nloc),Uy_d(Mloc,Nloc),Vx_d(Mloc,Nloc),Vy_d(Mloc,Nloc), &
              DUx_d(Mloc,Nloc),DUy_d(Mloc,Nloc),DVx_d(Mloc,Nloc),DVy_d(Mloc,Nloc), &
              ETAT_d(Mloc,Nloc),ETAx_d(Mloc,Nloc),ETAy_d(Mloc,Nloc), &
              ETATx_d(Mloc,Nloc),ETATy_d(Mloc,Nloc), &
! time-derivatives
              U0_d(Mloc,Nloc),V0_d(Mloc,Nloc),Ut_d(Mloc,Nloc),Vt_d(Mloc,Nloc),&
              Utx_d(Mloc,Nloc),Vty_d(Mloc,Nloc),Utxx_d(Mloc,Nloc),Utxy_d(Mloc,Nloc),&
              Vtxy_d(Mloc,Nloc),Vtyy_d(Mloc,Nloc),&
              DUtxx_d(Mloc,Nloc),DUtxy_d(Mloc,Nloc),&
              DVtxy_d(Mloc,Nloc),DVtyy_d(Mloc,Nloc),DUtx_d(Mloc,Nloc),DVty_d(Mloc,Nloc),&
! P Q Eta, Fx, Fy
              PL_d(Mloc1,Nloc),PR_d(Mloc1,Nloc),QL_d(Mloc,Nloc1),QR_d(Mloc,Nloc1), &
              FxL_d(Mloc1,Nloc),FxR_d(Mloc1,Nloc),FyL_d(Mloc,Nloc1),FyR_d(Mloc,Nloc1), &
              GxL_d(Mloc1,Nloc),GxR_d(Mloc1,Nloc),GyL_d(Mloc,Nloc1),GyR_d(Mloc,Nloc1), &
              EtaRxL_d(Mloc1,Nloc),EtaRxR_d(Mloc1,Nloc), &
              EtaRyL_d(Mloc,Nloc1),EtaRyR_d(Mloc,Nloc1), &
! sponge
              SPONGE_d(Mloc,Nloc), SpongeMaker_d(Mloc,Nloc), &
! original variables at notes
              Fx_d(Mloc1,Nloc),Fy_d(Mloc,Nloc1),&
              U_d(Mloc,Nloc),V_d(Mloc,Nloc), HU_d(Mloc,Nloc),HV_d(Mloc,Nloc),&
              Gx_d(Mloc1,Nloc),Gy_d(Mloc,Nloc1), &
              P_d(Mloc1,Nloc),Q_d(Mloc,Nloc1), &
              SxL_d(Mloc1,Nloc),SxR_d(Mloc1,Nloc), &
              SyL_d(Mloc,Nloc1),SyR_d(Mloc,Nloc1),SourceX_d(Mloc,Nloc), &
              SourceY_d(Mloc,Nloc), &
! others
              Umean_d(Mloc,Nloc),Vmean_d(Mloc,Nloc),ETAmean_d(Mloc,Nloc),&
              Usum_d(Mloc,Nloc),Vsum_d(Mloc,Nloc),ETAsum_d(Mloc,Nloc), &
              nu_smg_d(Mloc,Nloc), &
              Num_Zero_Up_d(Mloc,Nloc), &
              WaveHeightRMS_d(Mloc,Nloc),  &
              WaveHeightAve_d(Mloc,Nloc),  &
              Emax_d(Mloc,Nloc),  &
              Emin_d(Mloc,Nloc), &
              HrmsSum_d(Mloc,Nloc), &
              HavgSum_d(Mloc,Nloc), &
	      !ykchoi
	      ETA2sum_d(Mloc,Nloc), ETA2mean_d(Mloc,Nloc), &
                 SigWaveHeight_d(Mloc,Nloc),  &
# if defined (CARTESIAN)
              U4xL_d(Mloc1,Nloc),U4xR_d(Mloc1,Nloc),&
              V4yL_d(Mloc,Nloc1),V4yR_d(Mloc,Nloc1), &
	! ykchoi added V4xL and U4yL (08/06/15)
			V4xL_d(Mloc1,Nloc),V4xR_d(Mloc1,Nloc),&  
			U4yL_d(Mloc,Nloc1),U4yR_d(Mloc,Nloc1),& 
            DelxU4_d(Mloc,Nloc),DelxV4_d(Mloc,Nloc),&
            DelyV4_d(Mloc,Nloc),DelyU4_d(Mloc,Nloc) &
# else
# if defined (ZALPHA)
              U4xL_d(Mloc1,Nloc),U4xR_d(Mloc1,Nloc),&
              V4yL_d(Mloc,Nloc1),V4yR_d(Mloc,Nloc1), &
        ! ykchoi added V4xL and U4yL (08/06/15)
	      V4xL_d(Mloc1,Nloc),V4xR_d(Mloc1,Nloc),&  
	      U4yL_d(Mloc,Nloc1),U4yR_d(Mloc,Nloc1), &  
            DelxU4_d(Mloc,Nloc),DelxV4_d(Mloc,Nloc),&
            DelyV4_d(Mloc,Nloc),DelyU4_d(Mloc,Nloc),&
# endif
              Dx_d(Mloc,Nloc),Dy_d(Mloc,Nloc), Coriolis_d(Mloc,Nloc), &
              SlopeX_d(Mloc,Nloc),SlopeY_d(Mloc,Nloc), &
              Lat_theta_d(Mloc,Nloc) &
# endif
              )
      ALLOCATE(Depth_d(Mloc,Nloc),H_d(Mloc,Nloc),&
               Depthx_d(Mloc1,Nloc),Depthy_d(Mloc,Nloc1), &
               MASK_d(Mloc,Nloc),DepthNode_d(Mloc1,Nloc1), MASKtmp_d(Mloc,Nloc), &
               MASK_STRUC_d(Mloc,Nloc),MASK9_d(Mloc,Nloc), &
               tmp4preview_d(Mloc,Nloc),Int2Flo_d(Mloc,Nloc),&
               Cd_d(Mloc,Nloc),CD_breakwater_d(Mloc,Nloc), &
               MASK9u_d(1:Mloc1,1:Nloc), MASK9v_d(1:Mloc,1:Nloc1)&
              )
! updating variables
      ALLOCATE(Eta_d(Mloc,Nloc),Eta0_d(Mloc,Nloc),ETAtmp_d(Mloc,Nloc), &
               Ubar0_d(Mloc,Nloc),Vbar0_d(Mloc,Nloc),&
               Ubar_d(Mloc,Nloc),Vbar_d(Mloc,Nloc))
! dispersion updating variables
# if defined (CARTESIAN)
      ALLOCATE(U4_d(Mloc,Nloc),V4_d(Mloc,Nloc),U1p_d(Mloc,Nloc), & 
               V1p_d(Mloc,Nloc),U1pp_d(Mloc,Nloc),V1pp_d(Mloc,Nloc),&
               U2_d(Mloc,Nloc),V2_d(Mloc,Nloc),U3_d(Mloc,Nloc),V3_d(Mloc,Nloc))
# else
# if defined (ZALPHA)
      ALLOCATE(U4_d(Mloc,Nloc),V4_d(Mloc,Nloc) )
# endif
      ALLOCATE(U1p_d(Mloc,Nloc),V1p_d(Mloc,Nloc) )
# endif
# if defined (ITERATION)
      ALLOCATE(Ui_d(Mloc,Nloc), Vi_d(Mloc,Nloc),ETAi_d(Mloc,Nloc),&
        UbarOld_d(Mloc,Nloc),VbarOld_d(Mloc,Nloc),EtaOld_d(Mloc,Nloc) )
# endif

  ! HeightMax will be used not only in output but also meteo module
      ALLOCATE(HeightMax_d(Mloc,Nloc))

      ALLOCATE(WaveMaker_Mass_d(Mloc,Nloc))

      IF(VISCOSITY_BREAKING.OR.SHOW_BREAKING)THEN
       ALLOCATE(AGE_BREAKING_d(Mloc,Nloc))
       !YUAN, new var, details refer to breaker_gpu.F
       ALLOCATE(AGE_BREAKING0_d(Mloc,Nloc))
       ALLOCATE(nu_break_d(Mloc,Nloc))
       nu_break_d=nu_bkg
      ENDIF
	
      IF(WAVEMAKER_VIS)THEN
       ALLOCATE(nu_break_d(Mloc,Nloc))
       nu_break_d=ZERO
      ENDIF

      IF(DIFFUSION_SPONGE)THEN
       ALLOCATE(nu_sponge_d(Mloc,Nloc))
       nu_sponge_d=ZERO
      ENDIF

      IF(OUT_Hmin)THEN
        ALLOCATE(HeightMin_d(Mloc,Nloc))
        HeightMin_d=ZERO
      ENDIF
      IF(OUT_Umax)THEN
        ALLOCATE(VelocityMax_d(Mloc,Nloc))
        VelocityMax_d=ZERO
      ENDIF
      IF(OUT_VORmax)THEN
        ALLOCATE(VorticityMax_d(Mloc,Nloc))
        VorticityMax_d=ZERO
      ENDIF
      IF(OUT_MFmax)THEN
        ALLOCATE(MomentumFluxMax_d(Mloc,Nloc))
        MomentumFluxMax_d=ZERO
      ENDIF
end subroutine ALLOCATE_VARIABLES_GPU



SUBROUTINE INITIALIZATION_GPU
     IMPLICIT NONE

! set zeros
!    T_sum_d = 0.0_SP
! Kernel execution time
    TotalKernelTime = 0.0_SP
!
    Umean_d = 0.0_SP
    Vmean_d = 0.0_SP
    ETAmean_d = 0.0_SP
    nu_smg_d = 0.0_SP
    Num_Zero_Up_d = 0.0_SP
    WaveHeightRMS_d = 0.0_SP
    WaveHeightAve_d = 0.0_SP 
    Emax_d  = 0.0_SP
    Emin_d  = 0.0_SP
    HrmsSum_d  = 0.0_SP
    HavgSum_d  = 0.0_SP
    DelxU_d = 0.0_SP
    DelxHU_d =0.0_SP
    DelxV_d =0.0_SP
    DelxEtar_d =0.0_SP
    DelyU_d =0.0_SP
    DelyHV_d =0.0_SP
    DelyV_d =0.0_SP
    DelyEtar_d =0.0_SP
    DelxHV_d =0.0_SP
    DelyHU_d =0.0_SP
    UxL_d =0.0_SP
    UxR_d =0.0_SP
    VxL_d =0.0_SP
    VxR_d =0.0_SP
    HUxL_d =0.0_SP
    HUxR_d =0.0_SP
    HVxL_d =0.0_SP
    HVxR_d =0.0_SP
    HxL_d =0.0_SP
    HxR_d =0.0_SP
    UyL_d =0.0_SP
    UyR_d =0.0_SP
    VyL_d =0.0_SP
    VyR_d =0.0_SP
    HVyL_d =0.0_SP
    HVyR_d =0.0_SP
    HUyL_d =0.0_SP
    HUyR_d =0.0_SP
    HyL_d =0.0_SP
    HyR_d =0.0_SP
# if defined (CARTESIAN)
    U4xL_d =0.0_SP
    U4xR_d =0.0_SP
    V4yL_d =0.0_SP
    V4yR_d =0.0_SP
# endif
    PL_d=0.0_SP
    PR_d=0.0_SP
    QL_d=0.0_SP
    QR_d=0.0_SP
    FxL_d=0.0_SP
    FxR_d=0.0_SP
    FyL_d=0.0_SP
    FyR_d=0.0_SP
    GxL_d=0.0_SP
    GxR_d=0.0_SP
    GyL_d=0.0_SP
    GyR_d=0.0_SP
    SxL_d=0.0_SP
    SxR_d=0.0_SP
    SyL_d=0.0_SP
    SyR_d=0.0_SP
    Ubar0_d=0.0_SP
    Vbar0_d=0.0_SP
    Fx_d=0.0_SP
    Fy_d=0.0_SP
    Gx_d=0.0_SP
    Gy_d=0.0_SP
!    DepthNode_d=10.0_SP
    SourceX_d=0.0_SP
    SourceY_d=0.0_SP
    MASK_STRUC_d=1
    HeightMax_d=0.0_SP
    WaveMaker_Mass_d = 0.0_SP
!
    Beta_ref_d = -0.531_SP
!allocate and upload device variables involved in the 
!different types of wavemakers
    call allocate_wavemaker_variables_gpu
! coordinate for cartesian only
    Xco_d = Xco
    Yco_d = Yco
! for spherical coordinates and Coriolis parameter
    Dx_d = Dx
    Dy_d = Dy
# if !defined (CARTESIAN)
    Coriolis_d = Coriolis
    Lat_theta_d = Lat_theta
# endif
! bathymetry
    DEPTH_d = DEPTH
    DepthX_d = DepthX
    DepthY_d = DepthY
! get Eta and H
     ETA_d = ETA
     H_d = H
     HU_d = HU
     HV_d = HV
     Ubar_d = Ubar
     Vbar_d = Vbar
     U_d = U           
     V_d = V           
     MASK_d = MASK
     MASK9_d = MASK9
!friction and breaking
    Cd_d = Cd

! corresponding to cal_dispersion in init.F    
     P_d            = P           
     Q_d            = Q           
     U0_d           = U0           
     V0_d           = V0         
     Uxx_d          = Uxx         
     DUxx_d         = DUxx        
     Vyy_d          = Vyy         
     DVyy_d         = DVyy        
     Uxy_d          = Uxy         
     DUxy_d         = DUxy        
     Vxy_d          = Vxy         
     DVxy_d         = DVxy        
     ETAT_d         = ETAT        
     ETATx_d        = ETATx       
     ETATy_d        = ETATy       
     Ut_d           = Ut          
     Vt_d           = Vt          
     Utx_d          = Utx         
     Vty_d          = Vty         
     Utxx_d         = Utxx        
     Vtyy_d         = Vtyy        
     Utxy_d         = Utxy        
     Vtxy_d         = Vtxy        
     DUtxx_d        = DUtxx       
     DVtyy_d        = DVtyy       
     DUtxy_d        = DUtxy       
     DVtxy_d        = DVtxy       
     Ux_d           = Ux          
     DUx_d          = DUx         
     Vy_d           = Vy          
     DVy_d          = DVy         
     Uy_d           = Uy          
     DUy_d          = DUy         
     Vx_d           = Vx          
     DVx_d          = DVx         
     ETAx_d         = ETAx        
     ETAy_d         = ETAy        
     DUtx_d         = DUtx
     DVty_d         = DVty
# if defined (CARTESIAN)
     U4_d    = U4
     V4_d    = V4
     U1p_d          = U1p
     V1p_d          = V1p
     U1pp_d         = U1pp
     V1pp_d         = V1pp
     U2_d           = U2
     V2_d           = V2
     U3_d           = U3
     V3_d           = V3
# else
# if defined (ZALPHA)
     U4_d    = U4
     V4_d    = V4
# endif
     U1p_d          = U1p
     V1p_d          = V1p
# endif

     IF(OUT_VORmax) THEN
         VorticityMax_d = VorticityMax
     ENDIF
! sponge layer 
    SPONGE_d=1.0_SP
    SpongeMaker_d=1.0_SP
     IF(DIRECT_SPONGE)THEN
         SPONGE_d = SPONGE
     ENDIF

     IF(DIFFUSION_SPONGE)THEN
         nu_sponge_d = nu_sponge
     ENDIF

     IF(FRICTION_SPONGE)THEN
         IF(.NOT.ALLOCATED(CD_4_SPONGE_d)) ALLOCATE(CD_4_SPONGE_d(Mloc,Nloc))
         CD_4_SPONGE_d = CD_4_SPONGE
     ENDIF

     IF(WaveMaker(1:3)=='ABS')THEN
         SpongeMaker_d = SpongeMaker
     ENDIF
    IF(SHOW_BREAKING)THEN
        AGE_BREAKING_d = 0.0_SP
    ENDIF
! read breakwater 
     IF(BREAKWATER)THEN
         CD_breakwater_d = CD_breakwater
          BreakWaterWidth_d = BreakWaterWidth
     ENDIF  ! endif breakwater
! end breakwater
! read obstacle structures 
     IF(OBSTACLE)THEN
        MASK_STRUC_d = MASK_STRUC
     ENDIF

# if defined (COUPLING)
! TO DO BY YUAN
     CALL COUPLING_INITIAL
# endif
 
# if defined(CARTESIAN)
# else
    SlopeX_d = SlopeX
    SlopeY_d = SlopeY
# endif
!!!!!WaveMaker Codes: GPU kernels dont allow character strings
    SELECT CASE (trim(adjustl(WAVEMAKER)))
        CASE ("INI_REC")
           WaveMakerCode = 0
        CASE ("LEF_SOL")
           WaveMakerCode = 1
        CASE ("INI_SOL")
           WaveMakerCode = 2
        CASE ("INI_OTH")
           WaveMakerCode = 3
        CASE ("WK_REG")
           WaveMakerCode = 4
        CASE ("WK_IRR")
           WaveMakerCode = 5
        CASE ("JON_2D")
           WaveMakerCode = 6
        CASE ("JON_1D")
           WaveMakerCode = 7
        CASE ("TMA_1D")
           WaveMakerCode = 8
        CASE ("WK_TIME_SERIES")
           WaveMakerCode = 9
        CASE ("WK_TIME")
           WaveMakerCode = 9
        CASE ("WAVE_DATA2D")
           WaveMakerCode = 10
        CASE ("WAVE_DATA")
           WaveMakerCode = 10
        CASE ("INI_GAU")
           WaveMakerCode = 11
        CASE ("GAUSIAN")
           WaveMakerCode = 11
        CASE ("LEFT_BC_IRR")
           WaveMakerCode = 12
        CASE ("N_WAVE")
           WaveMakerCode = 13
        CASE ("INI_DIP")
           WaveMakerCode = 14
        CASE ("ABS")
           WaveMakerCode = 15
    END SELECT
! Tridiagonal solver cusparse init
    if (TriSolver == 1) then
        call cusparse_init
    endif
# if defined (MGPU)
    call mpi_datatype
# endif

END SUBROUTINE INITIALIZATION_GPU

subroutine allocate_wavemaker_variables_gpu
    implicit none 

    IF(WaveMaker(1:9)=='WK_DATA2D')THEN
        ALLOCATE (Cm_d(Mloc,Nloc,NumFreq),Sm_d(Mloc,Nloc,NumFreq))
        ALLOCATE (OMGN2D_d(NumFreq))
        Cm_d = Cm
        Sm_d = Sm
        OMGN2D_d = OMGN2D
    ENDIF
    IF(WaveMaker(1:7)=='WK_TIME')THEN
         ALLOCATE(Beta_genS_d(NumWaveComp),D_genS_d(NumWaveComp) )
         Beta_genS_d = Beta_genS
         D_genS_d = D_genS
    ENDIF
    IF(WaveMaker(1:6)=='WK_IRR'.OR.WaveMaker(1:6)=='TMA_1D'  &
       .OR.WaveMaker(1:6)=='JON_1D'.OR.WaveMaker(1:6)=='JON_2D')THEN
         ALLOCATE(omgn_ir_d(Nfreq),Cm_d(Mloc,Nloc,Nfreq),Sm_d(Mloc,Nloc,Nfreq))
         omgn_ir_d = omgn_ir
         Cm_d = Cm
         Sm_d = Sm
     ENDIF
     IF(WAVEMAKER(1:3)=='ABS'.OR.WaveMaker(1:11)=='LEFT_BC_IRR')THEN
         IF(WAVE_DATA_TYPE(1:4)=='DATA')THEN
             ALLOCATE(Cm_eta_d(Mloc,Nloc,Numfreq),Sm_eta_d(Mloc,Nloc,Numfreq), &
                 Cm_u_d(Mloc,Nloc,Numfreq),Sm_u_d(Mloc,Nloc,Numfreq),&
                 Cm_v_d(Mloc,Nloc,Numfreq),Sm_v_d(Mloc,Nloc,Numfreq) )
             Cm_eta_d = Cm_eta
             Sm_eta_d = Sm_eta
             Cm_u_d = Cm_u
             Sm_u_d = Sm_u
             Cm_v_d = Cm_v
             Sm_v_d = Sm_v
         ELSE
             ALLOCATE (Amp_Ser_d(NumFreq,NumDir), Wave_Number_Ser_d(NumFreq), &
                 Per_Ser_d(NumFreq),Theta_Ser_d(NumDir),Segma_Ser_d(NumFreq), &
                 Phase_Ser_d(NumFreq))
             ALLOCATE(Cm_eta_d(Mloc,Nloc,Numfreq),Sm_eta_d(Mloc,Nloc,Numfreq), &
                 Cm_u_d(Mloc,Nloc,Numfreq),Sm_u_d(Mloc,Nloc,Numfreq),&
                 Cm_v_d(Mloc,Nloc,Numfreq),Sm_v_d(Mloc,Nloc,Numfreq) )
             Amp_Ser_d = Amp_Ser
             Wave_Number_Ser_d = Wave_Number_Ser
             Per_Ser_d = Per_Ser
             Theta_Ser_d = Theta_Ser
             Segma_Ser_d = Segma_Ser
             Phase_Ser_d = Phase_Ser
             Cm_eta_d = Cm_eta
             Sm_eta_d = Sm_eta
             Cm_u_d = Cm_u
             Sm_u_d = Sm_u
             Cm_v_d = Cm_v
             Sm_v_d = Sm_v
         ENDIF
     ENDIF
end subroutine allocate_wavemaker_variables_gpu

end module init_gpu
