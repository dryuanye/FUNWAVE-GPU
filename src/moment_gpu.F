!-------------------------------------------------------------------------------------
!
!    DispersionFluxesSources is the module to calculate Dispersion, Fluxes,
!        and source terms at GPU device
!    
!    HISTORY: 
!
!-------------------------------------------------------------------------------------
module DispersionFluxesSources
    use cudafor
    use PARAM
    use GLOBAL
    use mod_cuda
    use boundary_condition_module
    implicit none
    real(SP),dimension(:,:),allocatable,device ::DU_d,DV_d,DUt_d,DVt_d
    real(SP),dimension(:,:),allocatable,device ::omega
    !real(SP),dimension(:,:),allocatable,device :: Din_d
    real(SP),dimension(:,:),allocatable,device :: D1in_d,D2in_d,D3in_d,D4in_d,D5in_d,D6in_d,D7in_d

contains
subroutine cal_dispersion_gpu
    implicit none
    integer :: istat
    integer ::ierrSync, ierrAsync
    if (.NOT. allocated(omega)) allocate( omega(Mloc,Nloc))
    if (.NOT. allocated(DU_d)) allocate( DU_d(Mloc,Nloc))
    if (.NOT. allocated(DV_d)) allocate( DV_d(Mloc,Nloc))
    if (.NOT. allocated(DUt_d)) allocate( DUt_d(Mloc,Nloc))
    if (.NOT. allocated(DVt_d)) allocate( DVt_d(Mloc,Nloc))
! variables for GPU kernels
    tBlock = dim3 (BlockDimX_2D, BlockDimY_2D ,1)
    grid = dim3 ( ceiling ( real ( Mloc ) / BlockDimX_Inner_2D ) , &
                  ceiling ( real ( Nloc ) / BlockDimY_Inner_2D ) , 1) 
# if defined (CARTESIAN)
    call cal_dispersion_step1_kernel<<<grid,tBlock>>>&
        (Ibeg,Iend,Jbeg,Jend,Mloc,Nloc,DT,&
        MinDepthFrc,Gamma2,SHOW_BREAKING)
# else
    call cal_dispersion_step1_kernel<<<grid,tBlock>>>&
        (Ibeg,Iend,Jbeg,Jend,Mloc,Nloc,DT,&
        MinDepthFrc,SHOW_BREAKING)
# endif
    !istat = cudaDeviceSynchronize()
!    ierrSync = cudaGetLastError()
!    ierrAsync = cudaDeviceSynchronize()
!    if (ierrSync /= cudaSuccess) print *, &
!       'Sync kernel error:', cudaGetErrorString(ierrSync)
!    if (ierrAsync /= cudaSuccess) print *, &
!        'Async kernel error:', cudaGetErrorString(ierrAsync)

    tBlock = dim3 (BlockDimX_2D, BlockDimY_2D ,1)
    grid = dim3 ( ceiling ( real ( Mloc ) / BlockDimX_2D ) , &
                  ceiling ( real ( Nloc ) / BlockDimY_2D ) , 1) 
# if defined (CARTESIAN)
    call cal_dispersion_step1p5_kernel<<<grid,tBlock>>>&    
    (Mloc,Nloc,DT,MinDepthFrc,Gamma2)
# else
    call cal_dispersion_step1p5_kernel<<<grid,tBlock>>>&    
    (Mloc,Nloc,DT,MinDepthFrc)
# endif

    tBlock = dim3 (BlockDimX_2D, BlockDimY_2D ,1)
    grid = dim3 ( ceiling ( real ( Mloc ) / BlockDimX_Inner_2D ) , &
                  ceiling ( real ( Nloc ) / BlockDimY_Inner_2D ) , 1) 
# if defined (CARTESIAN)
    call cal_dispersion_step2_kernel<<<grid,tBlock>>>&
        (Ibeg,Iend,Jbeg,Jend,Gamma2)
# else
    call cal_dispersion_step2_kernel<<<grid,tBlock>>>&
        (Ibeg,Iend,Jbeg,Jend)
# endif
    !istat = cudaDeviceSynchronize()
!    ierrSync = cudaGetLastError()
!    ierrAsync = cudaDeviceSynchronize()
!    if (ierrSync /= cudaSuccess) print *, &
!       'Sync kernel error:', cudaGetErrorString(ierrSync)
!    if (ierrAsync /= cudaSuccess) print *, &
!        'Async kernel error:', cudaGetErrorString(ierrAsync)

    tBlock = dim3 (BlockDimX_2D,1,1)
    grid = dim3 ( ceiling ( real ( max(Mloc,Nloc) ) / BlockDimX_2D ),1,1) 
    call cal_dispersion_step3_kernel<<<grid,tBlock>>>&
        (Ibeg,Iend,Jbeg,Jend,Mloc,Nloc,n_west,n_east,n_suth,n_nrth,WaveMakerCode)
    !istat = cudaDeviceSynchronize()
!    ierrSync = cudaGetLastError()
!    ierrAsync = cudaDeviceSynchronize()
!    if (ierrSync /= cudaSuccess) print *, &
!       'Sync kernel error:', cudaGetErrorString(ierrSync)
!    if (ierrAsync /= cudaSuccess) print *, &
!        'Async kernel error:', cudaGetErrorString(ierrAsync)
    
    call EXCHANGE_DISPERSION_GPU
    
    tBlock = dim3 (BlockDimX_2D, BlockDimY_2D ,1)
    grid = dim3 ( ceiling ( real ( Mloc ) / BlockDimX_2D ) , &
                  ceiling ( real ( Nloc ) / BlockDimY_2D ) , 1) 
# if defined (CARTESIAN)
    call cal_dispersion_step4_kernel<<<grid,tBlock>>>&
        (Ibeg,Iend,Jbeg,Jend,Mloc,Nloc,Gamma2,Beta_1,Beta_2)
# else
    call cal_dispersion_step4_kernel<<<grid,tBlock>>>&
        (Ibeg,Iend,Jbeg,Jend,Mloc,Nloc,Beta_1,Beta_2)
# endif
!    ierrSync = cudaGetLastError()
!    ierrAsync = cudaDeviceSynchronize()
!    if (ierrSync /= cudaSuccess) print *, &
!       'Sync kernel error cal_dispersion step4:', cudaGetErrorString(ierrSync)
!    if (ierrAsync /= cudaSuccess) print *, &
!        'Async kernel error cal_dispersion step4:', cudaGetErrorString(ierrAsync)
!YUAN MGPU
    !call PHI_COLL_GPU(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,U4_d,2,PERIODIC)
    !call PHI_COLL_GPU(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,V4_d,3,PERIODIC)
# if defined (MGPU)
# if defined (CARTESIAN) || defined (ZALPHA)
    call phi_exch_cuda(U4_d)
    call phi_exch_cuda(V4_d)
# endif
# endif

# if defined (CARTESIAN)
    IF(gamma2>0.0_SP)THEN
        call cal_dispersion_step5_kernel<<<grid,tBlock>>>&
            (Ibeg,Iend,Jbeg,Jend,Beta_1,Beta_2)
!        ierrSync = cudaGetLastError()
!        ierrAsync = cudaDeviceSynchronize()
!        if (ierrSync /= cudaSuccess) print *, &
!           'Sync kernel error cal_dispersion step5:', cudaGetErrorString(ierrSync)
!        if (ierrAsync /= cudaSuccess) print *, &
!            'Async kernel error cal_dispersion step5:', cudaGetErrorString(ierrAsync)

        call cal_dispersion_step6_kernel<<<grid,tBlock>>>&
            (Ibeg,Iend,Jbeg,Jend,Gamma2,Beta_1,Beta_2,b2,OUT_VORmax)
!        ierrSync = cudaGetLastError()
!        ierrAsync = cudaDeviceSynchronize()
!        if (ierrSync /= cudaSuccess) print *, &
!           'Sync kernel error cal_dispersion step6:', cudaGetErrorString(ierrSync)
!        if (ierrAsync /= cudaSuccess) print *, &
!            'Async kernel error cal_dispersion step6:', cudaGetErrorString(ierrAsync)
    ENDIF
# endif

end subroutine cal_dispersion_gpu

attributes(device) subroutine derivative_x_kernel&
        (MASK,DX,Uin_L,Uin_R,Uout)
    implicit none
    real(SP),intent(in) :: Uin_L,Uin_R,DX
    real(SP),intent(out) :: Uout
    integer,intent(in) :: MASK
    Uout= (Uin_R-Uin_L)/DX/2.0_SP*MASK
end subroutine derivative_x_kernel

attributes(device) subroutine derivative_y_kernel&
        (MASK,DY,Uin_L,Uin_U,Uout)
    implicit none
    real(SP),intent(in) :: Uin_L,Uin_U,DY
    real(SP),intent(out) :: Uout
    integer,intent(in):: MASK
    Uout= (Uin_U-Uin_L)/DY/2.0_SP*MASK
end subroutine derivative_y_kernel

attributes(device) subroutine derivative_xx_kernel&
        (MASK,DX,Uin_L,Uin_R,Uin,Uout)
    implicit none
    real(SP),intent(in) :: Uin,Uin_L,Uin_R,DX
    real(SP),intent(out) :: Uout
    integer,intent(in):: MASK
    Uout= (Uin_R-2.0_SP*Uin+Uin_L)/DX/DX*MASK
end subroutine derivative_xx_kernel

attributes(device) subroutine derivative_yy_kernel&
        (MASK,DY,Uin_L,Uin_U,Uin,Uout)
    implicit none
    real(SP),intent(in) :: Uin,Uin_L,Uin_U,DY
    real(SP),intent(out) :: Uout
    integer,intent(in) :: MASK
    Uout= (Uin_U-2.0_SP*Uin+Uin_L)/DY/DY*MASK
end subroutine derivative_yy_kernel

attributes(device) subroutine derivative_xy_kernel&
        (MASK,DX,DY,Uin_LL,Uin_LU,Uin_RL,Uin_RU,Uout)
    implicit none
    real(SP),intent(in) :: Uin_LL,Uin_LU,Uin_RL,Uin_RU,DX,DY
    real(SP),intent(out) :: Uout
    integer,intent(in) :: MASK
    real(SP) :: tmp1,tmp2
    tmp1=(Uin_RU-Uin_RL)/2.0_SP/DY
    tmp2=(Uin_LU-Uin_LL)/2.0_SP/DY
    Uout= (tmp1-tmp2)/2.0_SP/DX*MASK
end subroutine derivative_xy_kernel

# if defined (CARTESIAN)
attributes(global) subroutine cal_dispersion_step1_kernel&
    (Ibeg,Iend,Jbeg,Jend,Mloc,Nloc,DT,&
    MinDepthFrc,Gamma2,SHOW_BREAKING)
# else
attributes(global) subroutine cal_dispersion_step1_kernel&
    (Ibeg,Iend,Jbeg,Jend,Mloc,Nloc,DT,&
    MinDepthFrc,SHOW_BREAKING)
# endif
    implicit none 
    integer,value :: Ibeg,Iend,Jbeg,Jend,Mloc,Nloc
    real(SP),value :: DT,MinDepthFrc
# if defined (CARTESIAN)
    real(SP),value :: Gamma2
# endif
    logical,value :: SHOW_BREAKING
    real(SP),dimension(BlockDimX_2D,BlockDimY_2D),shared :: &
            U_sh,V_sh,Eta_sh
    real(SP) :: DXg,DYg
    integer :: tMASK
    integer :: i,j,tx,ty

! set local indexes
    tx = threadIdx%x
    ty = threadIdx%y
! Global thread and  block ID
    i = tx + (blockIdx%x-1)*BlockDimX_Inner_2D    ! Inner must be smaller than BlockDim
    j = ty + (blockIdx%y-1)*BlockDimY_Inner_2D
! init shared memory
    if (i>=Ibeg-3 .AND. i<=Iend+3 .AND. j>=Jbeg-3 .AND. j<=Jend+3) then
        U_sh(tx,ty) = U_d(i,j)
        V_sh(tx,ty) = V_d(i,j)
        Eta_sh(tx,ty) = Eta_d(i,j)
    endif
    call syncthreads()

    if (i>=Ibeg-2 .AND. i<=Iend+2 .AND. j>=Jbeg-2 .AND. j<=Jend+2 .and. &
        tx>1 .and. ty>1 .and. tx<blockDim%x .and. ty<blockDim%y) then
# if defined (CARTESIAN)
        DXg=DX_d
        DYg=DY_d
# else
        DXg=DX_d(i,j)
        DYg=DY_d(i,j)
# endif
        tMASK   = MASK9_d(i,j)
!uxx
        call derivative_xx_kernel(tMASK,DXg,U_sh(tx-1,ty),U_sh(tx+1,ty),U_sh(tx,ty),Uxx_d(i,j))
! uxy
        call derivative_xy_kernel(tMASK,DXg,DYg,U_sh(tx-1,ty-1),U_sh(tx-1,ty+1),U_sh(tx+1,ty-1),U_sh(tx+1,ty+1),Uxy_d(i,j))
! vxy
        call derivative_xy_kernel(tMASK,DXg,DYg,V_sh(tx-1,ty-1),V_sh(tx-1,ty+1),V_sh(tx+1,ty-1),V_sh(tx+1,ty+1), Vxy_d(i,j))
! vyy
        call derivative_yy_kernel(tMASK,DYg,V_sh(tx,ty-1),V_sh(tx,ty+1),V_sh(tx,ty),Vyy_d(i,j))
# if defined (CARTESIAN)
! gamma2.ne.0
        IF(Gamma2>0.0_SP)THEN
            call derivative_x_kernel(tMASK,DXg,U_sh(tx-1,ty),U_sh(tx+1,ty),Ux_d(i,j))
            call derivative_x_kernel(tMASK,DXg,V_sh(tx-1,ty),V_sh(tx+1,ty),Vx_d(i,j))
            call derivative_y_kernel(tMASK,DYg,U_sh(tx,ty-1),U_sh(tx,ty+1),Uy_d(i,j))
            call derivative_y_kernel(tMASK,DYg,V_sh(tx,ty-1),V_sh(tx,ty+1),Vy_d(i,j))
            call derivative_x_kernel(tMASK,DXg,Eta_sh(tx-1,ty),Eta_sh(tx+1,ty),ETAx_d(i,j))
            call derivative_y_kernel(tMASK,DYg,Eta_sh(tx,ty-1),Eta_sh(tx,ty+1),ETAy_d(i,j))
        ELSEIF(SHOW_BREAKING)THEN
# else
        IF(SHOW_BREAKING)THEN
# endif
            call derivative_x_kernel(tMASK,DXg,Eta_sh(tx-1,ty),Eta_sh(tx+1,ty),ETAx_d(i,j))
            call derivative_y_kernel(tMASK,DYg,Eta_sh(tx,ty-1),Eta_sh(tx,ty+1),ETAy_d(i,j))
        ENDIF
    endif

end subroutine cal_dispersion_step1_kernel

# if defined (CARTESIAN)
attributes(global) subroutine cal_dispersion_step1p5_kernel&    
    (Mloc,Nloc,DT,MinDepthFrc,Gamma2)
# else
attributes(global) subroutine cal_dispersion_step1p5_kernel&    
    (Mloc,Nloc,DT,MinDepthFrc)
# endif
    implicit none 
    integer,value :: Mloc,Nloc
    real(SP),value :: DT,MinDepthFrc
# if defined (CARTESIAN)
    real(SP),value :: Gamma2
# endif
    real(SP) :: DXg,DYg
    integer :: i,j,tx,ty

! set local indexes
    tx = threadIdx%x
    ty = threadIdx%y
! Global thread and  block ID
    i = tx + (blockIdx%x-1)*BlockDimX_2D 
    j = ty + (blockIdx%y-1)*BlockDimY_2D
    if (i>=1 .and. i<Mloc .and. j>=1 .and. j<Nloc) then
# if defined (CARTESIAN)
        DXg=DX_d
        DYg=DY_d
# else
        DXg=DX_d(i,j)
        DYg=DY_d(i,j)
# endif
! DU DV
        DU_d(I,J)=Max(Depth_d(I,J),MinDepthFrc)*U_d(I,J)
        DV_d(I,J)=Max(Depth_d(I,J),MinDepthFrc)*V_d(I,J)
        ETAT_d(I,J)=-(P_d(I+1,J)-P_d(I,J))/DXg-(Q_d(I,J+1)-Q_d(I,J))/DYg
    endif

! ETAT
# if defined (CARTESIAN)
    IF(Gamma2>0.0_SP)THEN
        if (i>=1 .and. i<=Mloc .and. j>=1 .and. j<=Nloc) then
            Ut_d(I,J) = (U_d(I,J)-U0_d(I,J)) / DT   !ykchoi
            Vt_d(I,J) = (V_d(I,J)-V0_d(I,J)) / DT   !ykchoi
            DUt_d(I,J)=Max(Depth_d(I,J),MinDepthFrc)*Ut_d(I,J)
            DVt_d(I,J)=Max(Depth_d(I,J),MinDepthFrc)*Vt_d(I,J)
        endif

!    ELSEIF(SHOW_BREAKING .OR. WAVEMAKER_VIS)THEN
!# else
!    IF(SHOW_BREAKING .OR. WAVEMAKER_VIS)THEN
    ENDIF
!add by YUAN: why need ELSEIF in the above line, no codes within the elseif branch
# endif

end subroutine cal_dispersion_step1p5_kernel

# if defined (CARTESIAN)
attributes(global) subroutine cal_dispersion_step2_kernel&
    (Ibeg,Iend,Jbeg,Jend,Gamma2)
# else
attributes(global) subroutine cal_dispersion_step2_kernel&
    (Ibeg,Iend,Jbeg,Jend)
# endif
    implicit none 
    integer,value :: Ibeg,Iend,Jbeg,Jend
# if defined (CARTESIAN)
    real(SP),value :: Gamma2
# endif
    real(SP),dimension(BlockDimX_2D,BlockDimY_2D),shared :: &
        DU_sh,DV_sh,Ut_sh,Vt_sh,DUt_sh,DVt_sh
    real(SP) :: DXg,DYg
    integer :: tMASK
    integer :: i,j,tx,ty

! set local indexes
    tx = threadIdx%x
    ty = threadIdx%y
! Global thread and  block ID
    i = tx + (blockIdx%x-1)*BlockDimX_Inner_2D    ! Inner must be smaller than BlockDim
    j = ty + (blockIdx%y-1)*BlockDimY_Inner_2D
! init shared memory
    if (i>=Ibeg-3 .AND. i<=Iend+3 .AND. j>=Jbeg-3 .AND. j<=Jend+3) then
        DU_sh(tx,ty) = DU_d(i,j)
        DV_sh(tx,ty) = DV_d(i,j)
# if defined (CARTESIAN)
        if(Gamma2>0.0_SP)then
            Ut_sh(tx,ty) = Ut_d(i,j)
            Vt_sh(tx,ty) = Vt_d(i,j)
            DUt_sh(tx,ty) = DUt_d(i,j)
            DVt_sh(tx,ty) = DVt_d(i,j)
        endif
# endif
    endif
    call syncthreads()
    if (i>=Ibeg-2 .AND. i<=Iend+2 .AND. j>=Jbeg-2 .AND. j<=Jend+2 .and. &
        tx>1 .and. ty>1 .and. tx<blockDim%x .and. ty<blockDim%y) then
# if defined (CARTESIAN)
        DXg=DX_d
        DYg=DY_d
# else
        DXg=DX_d(i,j)
        DYg=DY_d(i,j)
# endif
        tMASK   = MASK9_d(i,j)
! DUxx
        call derivative_xx_kernel(tMASK,DXg,DU_sh(tx-1,ty),DU_sh(tx+1,ty),DU_sh(tx,ty),DUxx_d(i,j))
! DUxy
        call derivative_xy_kernel(tMASK,DXg,DYg,DU_sh(tx-1,ty-1),DU_sh(tx-1,ty+1),DU_sh(tx+1,ty-1),DU_sh(tx+1,ty+1),DUxy_d(i,j))
! DVxy
        call derivative_xy_kernel(tMASK,DXg,DYg,DV_sh(tx-1,ty-1),DV_sh(tx-1,ty+1),DV_sh(tx+1,ty-1),DV_sh(tx+1,ty+1),DVxy_d(i,j))
! DVyy
        call derivative_yy_kernel(tMASK,DYg,DV_sh(tx,ty-1),DV_sh(tx,ty+1),DV_sh(tx,ty),DVyy_d(i,j))
# if defined (CARTESIAN)
        IF(Gamma2>0.0_SP)THEN
            call derivative_x_kernel(tMASK,DXg,DU_sh(tx-1,ty),DU_sh(tx+1,ty),DUx_d(i,j))
            call derivative_x_kernel(tMASK,DXg,DV_sh(tx-1,ty),DV_sh(tx+1,ty),DVx_d(i,j))
            call derivative_y_kernel(tMASK,DYg,DU_sh(tx,ty-1),DU_sh(tx,ty+1),DUy_d(i,j))
            call derivative_y_kernel(tMASK,DYg,DV_sh(tx,ty-1),DV_sh(tx,ty+1),DVy_d(i,j))
            call derivative_x_kernel(tMASK,DXg,Ut_sh(tx-1,ty),Ut_sh(tx+1,ty),Utx_d(i,j))
            call derivative_y_kernel(tMASK,DYg,Vt_sh(tx,ty-1),Vt_sh(tx,ty+1),Vty_d(i,j))
            call derivative_xx_kernel(tMASK,DXg,Ut_sh(tx-1,ty),Ut_sh(tx+1,ty),Ut_sh(tx,ty),Utxx_d(i,j))
            call derivative_yy_kernel(tMASK,DYg,Vt_sh(tx,ty-1),Vt_sh(tx,ty+1),Vt_sh(tx,ty),Vtyy_d(i,j))
            call derivative_xy_kernel(tMASK,DXg,DYg,Ut_sh(tx-1,ty-1),Ut_sh(tx-1,ty+1),Ut_sh(tx+1,ty-1),Ut_sh(tx+1,ty+1),Utxy_d(i,j))
            call derivative_xy_kernel(tMASK,DXg,DYg,Vt_sh(tx-1,ty-1),Vt_sh(tx-1,ty+1),Vt_sh(tx+1,ty-1),Vt_sh(tx+1,ty+1),Vtxy_d(i,j))
            call derivative_x_kernel(tMASK,DXg,DUt_sh(tx-1,ty),DUt_sh(tx+1,ty),DUtx_d(i,j))
            call derivative_y_kernel(tMASK,DYg,DVt_sh(tx,ty-1),DVt_sh(tx,ty+1),DVty_d(i,j))
            call derivative_xx_kernel(tMASK,DXg,DUt_sh(tx-1,ty),DUt_sh(tx+1,ty),DUt_sh(tx,ty),DUtxx_d(i,j))
            call derivative_yy_kernel(tMASK,DYg,DVt_sh(tx,ty-1),DVt_sh(tx,ty+1),DVt_sh(tx,ty),DVtyy_d(i,j))
            call derivative_xy_kernel(tMASK,DXg,DYg,DUt_sh(tx-1,ty-1),DUt_sh(tx-1,ty+1),DUt_sh(tx+1,ty-1),DUt_sh(tx+1,ty+1),DUtxy_d(i,j))
            call derivative_xy_kernel(tMASK,DXg,DYg,DVt_sh(tx-1,ty-1),DVt_sh(tx-1,ty+1),DVt_sh(tx+1,ty-1),DVt_sh(tx+1,ty+1),DVtxy_d(i,j))
        ENDIF
# endif
    endif
end subroutine cal_dispersion_step2_kernel


attributes(global) subroutine cal_dispersion_step3_kernel&
    (Ibeg,Iend,Jbeg,Jend,Mloc,Nloc,n_west,n_east,n_suth,n_nrth,WaveMakerCode)
    implicit none 
    integer,value :: Ibeg,Iend,Jbeg,Jend,Mloc,Nloc
    integer,value :: n_west,n_east,n_suth,n_nrth
    integer,value :: WaveMakerCode
    integer :: i,tx

! set local indexes
    tx = threadIdx%x
! Global thread and  block ID
    i = tx + (blockIdx%x-1)*BlockDimX_2D 

!  boundary conditions
# if defined (MGPU)
    if(n_west.eq.MPI_PROC_NULL) then
# endif
# if defined(COUPLING)
    IF(.NOT.IN_DOMAIN_WEST)THEN
# endif
!LEFT_BC_IRR
    IF (WaveMakerCode==12)THEN
     ! do nothing
    ELSE
        if (i>=1 .and. i<=Nloc) then
            Uxy_d(Ibeg,i)=0.0_SP
            DUxy_d(Ibeg,i)=0.0_SP
            Vxy_d(Ibeg,i)=0.0_SP
            DVxy_d(Ibeg,i)=0.0_SP
            Utxy_d(Ibeg,i)=0.0_SP
            DUtxy_d(Ibeg,i)=0.0_SP
            Vtxy_d(Ibeg,i)=0.0_SP
            DVtxy_d(Ibeg,i)=0.0_SP
        endif
    ENDIF ! left_bc wavemaker
# if defined(COUPLING)  
   ENDIF 
# endif
  
# if defined (MGPU)
    endif  
# endif

# if defined (MGPU)
    if(n_east.eq.MPI_PROC_NULL) then
# endif
# if defined(COUPLING)
    IF(.NOT.IN_DOMAIN_EAST)THEN
# endif
    if (i>=1 .and. i<=Nloc) then
        Uxy_d(Iend,i)=0.0_SP
        DUxy_d(Iend,i)=0.0_SP
        Vxy_d(Iend,i)=0.0_SP
        DVxy_d(Iend,i)=0.0_SP
        Utxy_d(Iend,i)=0.0_SP
        DUtxy_d(Iend,i)=0.0_SP
        Vtxy_d(Iend,i)=0.0_SP
        DVtxy_d(Iend,i)=0.0_SP
    endif
# if defined(COUPLING)
    ENDIF
# endif
# if defined (MGPU)
    endif  
# endif
  
# if defined (MGPU)
    if(n_suth.eq.MPI_PROC_NULL) then
# endif
# if defined(COUPLING)
    IF(.NOT.IN_DOMAIN_SOUTH)THEN
# endif
    if (i>=1 .and. i<=Mloc) then
        Uxy_d(I,Jbeg)=0.0_SP
        DUxy_d(I,Jbeg)=0.0_SP
        Vxy_d(I,Jbeg)=0.0_SP
        DVxy_d(I,Jbeg)=0.0_SP
        Utxy_d(I,Jbeg)=0.0_SP
        DUtxy_d(I,Jbeg)=0.0_SP
        Vtxy_d(I,Jbeg)=0.0_SP
        DVtxy_d(I,Jbeg)=0.0_SP
    endif
# if defined(COUPLING)
    ENDIF
# endif
# if defined (MGPU)
    endif  
# endif

# if defined (MGPU)
    if(n_nrth.eq.MPI_PROC_NULL) then
# endif
# if defined(COUPLING)
    IF(.NOT.IN_DOMAIN_NORTH)THEN
# endif
    if (i>=1 .and. i<=Mloc) then
        Uxy_d(I,Jend)=0.0_SP
        DUxy_d(I,Jend)=0.0_SP
        Vxy_d(I,Jend)=0.0_SP
        DVxy_d(I,Jend)=0.0_SP
        Utxy_d(I,Jend)=0.0_SP
        DUtxy_d(I,Jend)=0.0_SP
        Vtxy_d(I,Jend)=0.0_SP
        DVtxy_d(I,Jend)=0.0_SP
    endif 
# if defined(COUPLING)
    ENDIF
# endif
# if defined (MGPU)
    endif  
# endif
     
end subroutine cal_dispersion_step3_kernel


# if defined (CARTESIAN)
attributes(global) subroutine cal_dispersion_step4_kernel&    
    (Ibeg,Iend,Jbeg,Jend,Mloc,Nloc,Gamma2,Beta_1,Beta_2)
# else
attributes(global) subroutine cal_dispersion_step4_kernel&    
    (Ibeg,Iend,Jbeg,Jend,Mloc,Nloc,Beta_1,Beta_2)
# endif
    implicit none
    real(SP),value :: Beta_1,Beta_2
# if defined (CARTESIAN)
    real(SP),value :: Gamma2
# endif
    integer,value :: Ibeg,Iend,Jbeg,Jend,Mloc,Nloc
    real(SP) ::  tMASK,tDepth, UxxVxy,UxyVyy,HUxxHVxy,HUxyHVyy,reta,ken1,&
            ken2
    real(SP) :: DXg,DYg
    integer :: i,j,tx,ty

! set local indexes
    tx = threadIdx%x
    ty = threadIdx%y
! Global thread and  block ID
    i = tx + (blockIdx%x-1)*BlockDimX_2D    ! Inner must be smaller than BlockDim
    j = ty + (blockIdx%y-1)*BlockDimY_2D
! calculate V1p  without nonlinear dispersion
    if (i>=1 .AND. i<=Mloc .AND. j>=1 .AND. j<=Nloc) then
# if defined (CARTESIAN)
        DXg=DX_d
        DYg=DY_d
# else
        DXg=DX_d(i,j)
        DYg=DY_d(i,j)
# endif
        tMASK = MASK9_d(i,j)
        tDepth = DEPTH_d(i,j)
# if defined (CARTESIAN)
        U4_d(I,J)=(1.0_SP/3.0_SP-Beta_1+0.5_SP*Beta_1*Beta_1)*tDepth*tDepth*(Uxx_d(I,J)+Vxy_d(I,J)) &
                +(Beta_1-1.0_SP/2.0_SP)*tDepth*(DUxx_d(I,J)+DVxy_d(I,J))
        V4_d(I,J)=(1.0_SP/3.0_SP-Beta_1+0.5_SP*Beta_1*Beta_1)*tDepth*tDepth*(Uxy_d(I,J)+Vyy_d(I,J)) &
                +(Beta_1-1.0_SP/2.0_SP)*tDepth*(DUxy_d(I,J)+DVyy_d(I,J)) 

        !------[ykchoi (04/14/2017)
        IF(gamma2>0.0_SP)THEN
	    UxxVxy = Uxx_d(I,J) + Vxy_d(I,J)
	    UxyVyy = Uxy_d(I,J) + Vyy_d(I,J)

	    HUxxHVxy = DUxx_d(I,J) + DVxy_d(I,J)
	    HUxyHVyy = DUxy_d(I,J) + DVyy_d(I,J)

	    !rh = Depth(I,J)
	    reta = Eta_d(I,J)

  	    ken1 = ( 1.0_SP/6.0_SP - Beta_1 + Beta_1*Beta_1 )*tDepth*reta*Beta_2   &
	          + ( 1.0_SP/2.0_SP*Beta_1*Beta_1 - 1.0_SP/6.0_SP )*reta*reta*Beta_2*Beta_2
	    ken2 = ( Beta_1 - 1.0_SP/2.0_SP )*reta*Beta_2

	    U4_d(I,J) = U4_d(I,J) + gamma2*tMASK*( ken1*UxxVxy + ken2*HUxxHVxy )
	    V4_d(I,J) = V4_d(I,J) + gamma2*tMASK*( ken1*UxyVyy + ken2*HUxyHVyy )
	 ENDIF
	 !------ykchoi (04/14/2017)]
# endif

# if defined (ZALPHA)
         U4_d(I,J)=(1.0_SP/3.0_SP-Beta_1+0.5_SP*Beta_1*Beta_1)  &
               *tDepth*tDepth*(Uxx_d(I,J)+Vxy_d(I,J)) &
                +(Beta_1-1.0_SP/2.0_SP)*tDepth*(DUxx_d(I,J)+DVxy_d(I,J))
         V4_d(I,J)=(1.0_SP/3.0_SP-Beta_1+0.5_SP*Beta_1*Beta_1)  &
               *tDepth*tDepth*(Uxy_d(I,J)+Vyy_d(I,J)) &
                +(Beta_1-1.0_SP/2.0_SP)*tDepth*(DUxy_d(I,J)+DVyy_d(I,J))
# endif

# if defined (CARTESIAN)
	 U1p_d(I,J)=0.5_SP*(1.0_SP-Beta_1)*(1.0_SP-Beta_1)  &
                *tDepth*tDepth*(Uxx_d(I,J)+Vxy_d(I,J)) &
               +(Beta_1-1.0_SP)*tDepth*(DUxx_d(I,J)+DVxy_d(I,J))
       
	 V1p_d(I,J)=0.5_SP*(1.0_SP-Beta_1)*(1.0_SP-Beta_1)  &
                *tDepth*tDepth*(Uxy_d(I,J)+Vyy_d(I,J)) &
               +(Beta_1-1.0_SP)*tDepth*(DUxy_d(I,J)+DVyy_d(I,J))
# else
	 U1p_d(I,J)=0.5_SP*(1.0_SP-Beta_1)*(1.0_SP-Beta_1)  &
                *tDepth*tDepth*(Uxx_d(I,J)+Vxy_d(I,J)) &
               +(Beta_1-1.0_SP)*tDepth*(DUxx_d(I,J)+DVxy_d(I,J))

         V1p_d(I,J)=0.5_SP*(1.0_SP-Beta_1)*(1.0_SP-Beta_1)  &
                *tDepth*tDepth*(Uxy_d(I,J)+Vyy_d(I,J)) &
               +(Beta_1-1.0_SP)*tDepth*(DUxy_d(I,J)+DVyy_d(I,J))
# if defined(SPH_EXTRA_DISP)
! extra terms for spherical za^2/2(-1/r0 tan theta vx +za(-1/r0 tan theta HVx), if this option if off, it will be included in a1 and b1 …
 
       U1p_d(I,J)=U1p_d(I,J)+0.5_SP*(1.0_SP-Beta_1)  &
                *tDepth*tDepth*(-1.0_SP/R_earth*  &
                TAN(Lat_theta_d(I,J))*0.5_SP*(VxL_d(I,J)+VxR_d(I,J))) &
               +(Beta_1-1.0_SP)*tDepth*(-1.0_SP/R_earth* &
                TAN(Lat_theta_d(I,J))*0.5_SP*(HVxL_d(I,J)+HVxR_d(I,J)))

! extra terms in y: za^2/2(1/r0 tan theta ux -1/r0 tan theta vy -1/r0^2/cos^2theta v)
!                  +za (1/r0 tan theta HUx -1/r0 tan theta HVy -1/r0^2/cos^2 HV)
       V1p_d(I,J)=V1p_d(I,J)+0.5_SP*(1.0_SP-Beta_1) &
                 *tDepth*tDepth &
                /R_earth*( TAN(Lat_theta_d(I,J))*0.5_SP*(UxL_d(I,J)+UxR_d(I,J)) &
                          -TAN(Lat_theta_d(I,J))*0.5_SP*(VyL_d(I,J)+VyR_d(I,J)) &
                          -V_d(I,J)/R_earth/COS(Lat_theta_d(I,J))/COS(Lat_theta_d(I,J)) ) &
               +(Beta_1-1.0_SP)*tDepth &
                /R_earth*( 0.5_SP*(HUxL_d(I,J)+HUxR_d(I,J))*TAN(Lat_theta_d(I,J)) &
                          -0.5_SP*(HVyL_d(I,J)+HVyR_d(I,J))*TAN(Lat_theta_d(I,J)) &
                          -HV_d(I,J)/R_earth/COS(Lat_theta_d(I,J))/COS(Lat_theta_d(I,J)) &        
                 )
# endif
! end spherical extra dispersion terms which can be implemented in TIME_LEFT
# endif

    endif
end subroutine cal_dispersion_step4_kernel


# if defined (CARTESIAN)
attributes(global) subroutine cal_dispersion_step5_kernel&    
    (Ibeg,Iend,Jbeg,Jend,Beta_1,Beta_2)
    implicit none
    real(SP),value :: Beta_1,Beta_2
    integer,value :: Ibeg,Iend,Jbeg,Jend
    real(SP) :: UxxVxy,UxyVyy,HUxxHVxy,HUxyHVyy, &
                reta,ken1,ken2,ken3,ken4,tETAx,tETAy
    real(SP) :: DXg,DYg
    integer :: i,j,tx,ty

! set local indexes
    tx = threadIdx%x
    ty = threadIdx%y
! Global thread and  block ID
    i = tx + (blockIdx%x-1)*BlockDimX_2D
    j = ty + (blockIdx%y-1)*BlockDimY_2D

    if (i>=Ibeg .AND. i<=Iend .AND. j>=Jbeg .AND. j<=Jend ) then
        tETAx = ETAx_d(i,j)  
        tETAy = ETAy_d(i,j) 
! 
        UxxVxy=Uxx_d(i,j)+Vxy_d(i,j) 
        UxyVyy=Uxy_d(i,j)+Vyy_d(i,j)
        HUxxHVxy=DUxx_d(i,j)+DVxy_d(i,j)
        HUxyHVyy=DUxy_d(i,j)+DVyy_d(i,j)
        reta=Eta_d(i,j) 

        U1pp_d(I,J)=-reta*Beta_2*tETAx*Beta_2*(Utx_d(I,J)+Vty_d(I,J)) - 0.5_SP*reta*reta*Beta_2*Beta_2*(Utxx_d(I,J)+Vtxy_d(I,J))&
                  -tETAx*Beta_2*(DUtx_d(I,J)+DVty_d(I,J)) -reta*Beta_2*(DUtxx_d(I,J)+DVtxy_d(I,J))

        V1pp_d(I,J)=-reta*Beta_2*tETAy*Beta_2*(Utx_d(I,J)+Vty_d(I,J)) - 0.5_SP*reta*reta*Beta_2*Beta_2*(Utxy_d(I,J)+Vtyy_d(I,J))&
                  -tETAy*Beta_2*(DUtx_d(I,J)+DVty_d(I,J)) -reta*Beta_2*(DUtxy_d(I,J)+DVtyy_d(I,J))
        
	  !ykchoi (2015. 08.06.)
	  !remaining terms in U1p, V1p in Shi et al. (2012) publication   
	  ken1 = Beta_1*( 1.0_SP - Beta_1 )*Depth_d(i,j)*ETAT_d(I,J)*Beta_2 - Beta_1*Beta_1*reta*Beta_2*ETAT_d(I,J)*Beta_2
	  ken2 = Beta_1*( 1.0_SP - Beta_1 )*Depth_d(i,j)*reta*Beta_2 - 0.5*Beta_1*Beta_1*reta*reta*Beta_2*Beta_2
	  ken3 = Beta_1*ETAT_d(I,J)*Beta_2
	  ken4 = Beta_1*reta*Beta_2

	  U1pp_d(I,J) = U1pp_d(I,J) - ken1*UxxVxy - ken2*( Utxx_d(I,J)+Vtxy_d(I,J) ) + ken3*HUxxHVxy + ken4*(DUtxx_d(I,J)+DVtxy_d(I,J))
	  V1pp_d(I,J) = V1pp_d(I,J) - ken1*UxyVyy - ken2*( Utxy_d(I,J)+Vtyy_d(I,J) ) + ken3*HUxyHVyy + ken4*(DUtxy_d(I,J)+DVtyy_d(I,J))
    endif
end subroutine cal_dispersion_step5_kernel
# endif


# if defined (CARTESIAN)
attributes(global) subroutine cal_dispersion_step6_kernel&    
    (Ibeg,Iend,Jbeg,Jend,Gamma2,Beta_1,Beta_2,b2,OUT_VORmax)
    implicit none
    real(SP),value :: Beta_1,Beta_2,Gamma2,b2
    integer,value :: Ibeg,Iend,Jbeg,Jend
    logical,value :: OUT_VORmax
    real(SP) :: UxxVxy,UxyVyy,HUxxHVxy,HUxyHVyy, &
                UxxVxy_x,UxxVxy_y,UxyVyy_x,UxyVyy_y, &
                HUxxHVxy_x,HUxxHVxy_y,HUxyHVyy_x,HUxyHVyy_y, &
                rh,rhx,rhy,reta,ken1,ken2,ken3,ken4,ken5, &
                tU,tV,tMASK,tDepth,tETAx,tETAy,tUx,tVy
    real(SP) :: omega_0,omega_1
    integer :: i,j,tx,ty

! set local indexes
    tx = threadIdx%x
    ty = threadIdx%y
! Global thread and  block ID
    i = tx + (blockIdx%x-1)*BlockDimX_2D 
    j = ty + (blockIdx%y-1)*BlockDimY_2D

    if (i>=Ibeg .AND. i<=Iend .AND. j>=Jbeg .AND. j<=Jend) then
        tU = U_d(i,j) 
        tV = V_d(i,j)
        tETAx = ETAx_d(i,j)
        tETAy = ETAy_d(i,j)
        tUx = Ux_d(i,j)
        tVy = Vy_d(i,j)
! 
        UxxVxy=Uxx_d(i,j)+Vxy_d(i,j) 
        UxyVyy=Uxy_d(i,j)+Vyy_d(i,j)
        UxxVxy_x=(Uxx_d(i+1,j)+Vxy_d(i+1,j)-Uxx_d(i-1,j)-Vxy_d(i-1,j))/2.0_SP/DX_d
        UxxVxy_y=(Uxx_d(i,j+1)+Vxy_d(i,j+1)-Uxx_d(i,j-1)-Vxy_d(i,j-1))/2.0_SP/DY_d
        UxyVyy_x=(Uxy_d(i+1,j)+Vyy_d(i+1,j)-Uxy_d(i-1,j)-Vyy_d(i-1,j))/2.0_SP/DX_d 
        UxyVyy_y=(Uxy_d(i,j+1)+Vyy_d(i,j+1)-Uxy_d(i,j-1)-Vyy_d(i,j-1))/2.0_SP/DY_d

        HUxxHVxy=DUxx_d(i,j)+DVxy_d(i,j) 
        HUxyHVyy=DUxy_d(i,j)+DVyy_d(i,j) 
        HUxxHVxy_x=(DUxx_d(i+1,j)+DVxy_d(i+1,j)-DUxx_d(i-1,j)-DVxy_d(i-1,j))/2.0_SP/DX_d 
        HUxxHVxy_y=(DUxx_d(i,j+1)+DVxy_d(i,j+1)-DUxx_d(i,j-1)-DVxy_d(i,j-1))/2.0_SP/DY_d 
        HUxyHVyy_x=(DUxy_d(i+1,j)+DVyy_d(i+1,j)-DUxy_d(i-1,j)-DVyy_d(i-1,j))/2.0_SP/DX_d
        HUxyHVyy_y=(DUxy_d(i,j+1)+DVyy_d(i,j+1)-DUxy_d(i,j-1)-DVyy_d(i,j-1))/2.0_SP/DY_d

        rh=Depth_d(i,j) 
        rhx=(Depth_d(i+1,j)-Depth_d(i-1,j))/2.0_SP/DX_d 
        rhy=(Depth_d(i,j+1)-Depth_d(i,j-1))/2.0_SP/DY_d
        reta=Eta_d(i,j) 

          
         ken1=(Beta_1-1.0_SP)*(rhx+tETAx)*Beta_2
         ken2=(Beta_1-1.0_SP)*(rh+reta)*Beta_2
         ken3=( (1.0_SP-Beta_1)*(1.0_SP-Beta_1)*rh*rhx*Beta_2*Beta_2-Beta_1*(1.0_SP-Beta_1)*(rhx*reta*Beta_2+rh*tETAx*Beta_2) &
                    +(Beta_1*Beta_1-1.0_SP)*reta*tETAx*Beta_2*Beta_2 )
         ken4=( 0.5_SP*(1.0_SP-Beta_1)*(1.0_SP-Beta_1)*rh*rh*Beta_2*Beta_2-Beta_1*(1.0_SP-Beta_1)*rh*reta*Beta_2 &
                      +0.5_SP*(Beta_1*Beta_1-1.0_SP)*reta*reta*Beta_2*Beta_2 )
         ken5=( (1.0_SP-Beta_1)*(1.0_SP-Beta_1)*rh*rhy*Beta_2*Beta_2-Beta_1*(1.0_SP-Beta_1)*(rhy*reta*Beta_2+rh*tETAy*Beta_2) &
                    +(Beta_1*Beta_1-1.0_SP)*reta*tETAy*Beta_2*Beta_2 )

        U2_d(I,J)=ken1*(tU*HUxxHVxy+tV*HUxyHVyy) &
                +ken2*(tUx*HUxxHVxy+tU*HUxxHVxy_x &
                    +Vx_d(i,j)*HUxyHVyy+tV*HUxyHVyy_x) &
                +ken3 & 
                   *(tU*UxxVxy+tV*UxyVyy) &
                +ken4  &
                   *(tUx*UxxVxy+tU*UxxVxy_x+Vx_d(i,j)*UxyVyy+tV*UxyVyy_x) &
                +Beta_2*Beta_2*(DUx_d(I,J)+DVy_d(I,J)+reta*Beta_2*(tUx+tVy))  &
                   *(HUxxHVxy+tETAx*Beta_2*(tUx+tVy)+reta*Beta_2*UxxVxy)

        ! ykchoi (15. 08. 06.)
	  ! ken1 should be modified as following term.
        ken1=(Beta_1-1.0_SP)*(rhy+tETAy)*Beta_2

        V2_d(I,J)=ken1*(tV*HUxxHVxy+tV*HUxyHVyy) &
                +ken2*(Uy_d(i,j)*HUxxHVxy+tU*HUxxHVxy_y &
                    +tVy*HUxyHVyy+tV*HUxyHVyy_y) &
                +ken5 & 
                   *(tU*UxxVxy+tV*UxyVyy) &
                +ken4  &
                   *(Uy_d(i,j)*UxxVxy+tU*UxxVxy_y+tVy*UxyVyy+tV*UxyVyy_y) &
                +Beta_2*Beta_2*(DUx_d(I,J)+DVy_d(I,J)+reta*Beta_2*(tUx+tVy))  &
                   *(HUxyHVyy+tETAy*Beta_2*(tUx+tVy)+reta*Beta_2*UxyVyy)

        omega_0=Vx_d(i,j)-Uy_d(i,j)

! ykchoi (15. 08. 06.)
! omega_1 term is modified as Shi et al. (2012) publication
        omega_1=( b2*rhx + Beta_1*tETAx )*Beta_2*( HUxyHVyy + (b2*rh+Beta_1*reta)*Beta_2*UxyVyy )  &
              - ( b2*rhy + Beta_1*tETAy )*Beta_2*( HUxxHVxy + (b2*rh+Beta_1*reta)*Beta_2*UxxVxy )

        omega(I,J)=omega_0+omega_1

        IF(OUT_VORmax) THEN
        IF(abs(omega(I,J)).GT.VorticityMax_d(I,J)) THEN
        VorticityMax_d(I,J)=omega(I,J)
        ENDIF
        ENDIF

       ken1=((Beta_1-1.0_SP/2.0_SP)*(reta+rh)*Beta_2)
       ken2=(1.0_SP/3.0_SP-Beta_1+0.5_SP*Beta_1*Beta_1)*rh*rh*Beta_2*Beta_2  &
               + (1.0_SP/6.0_SP-Beta_1+Beta_1*Beta_1)*rh*reta*Beta_2 &
               +(1.0_SP/2.0_SP*Beta_1*Beta_1-1.0_SP/6.0_SP)*reta*reta*Beta_2*Beta_2

       U3_d(I,J)=-tV*omega_1 - omega_0 &
                 *(ken1*HUxyHVyy &
                   +ken2*UxyVyy)
! ykchoi

	 V3_d(I,J) = tU*omega_1 + omega_0 &
	            *(ken1*HUxxHVxy &
                    +ken2*UxxVxy)
! ykchoi
    endif
end subroutine cal_dispersion_step6_kernel
# endif


subroutine fluxes_gpu
    implicit none
    IF(HIGH_ORDER(1:3)=='FOU') THEN       !ykchoi (08/28/2016)
      ! for the fourth-order, Choi used 4th + combined minmod and van leer
      ! it is more stable than 4th + minmod
        CALL CONSTRUCTION_HO_GPU
        CALL wave_speed_gpu
    ELSE                                     
        CALL DelxyFun_gpu
        CALL CONSTRUCTION_GPU
        CALL wave_speed_gpu
    ENDIF

    IF(CONSTR(1:3)=='HLL')THEN
        CALL FLUX_AT_INTERFACE_HLL_GPU
    ELSE
        CALL FLUX_AT_INTERFACE_GPU
    ENDIF

   CALL BOUNDARY_CONDITION_FLUXES_GPU

end subroutine fluxes_gpu

subroutine FLUX_AT_INTERFACE_GPU

! for averaging approach for predictor
!$cuf kernel do(2) <<<*,*>>>
    do j = 1,Nloc  
        do i = 1,Mloc1
            P_d(i,j)=0.5_SP*(PR_d(i,j)+PL_d(i,j))
            Fx_d(i,j)=0.5_SP*(FxR_d(i,j)+FxL_d(i,j))
            Gx_d(i,j)=0.5_SP*(GxR_d(i,j)+GxL_d(i,j))
        enddo
    enddo
!$cuf kernel do(2) <<<*,*>>>
    do j = 1,Nloc1 
        do i = 1,Mloc
            Q_d(i,j)=0.5_SP*(QR_d(i,j)+QL_d(i,j))
            Fy_d(i,j)=0.5_SP*(FyR_d(i,j)+FyL_d(i,j))
            Gy_d(i,j)=0.5_SP*(GyR_d(i,j)+GyL_d(i,j))
        enddo
    enddo
end subroutine FLUX_AT_INTERFACE_GPU

subroutine FLUX_AT_INTERFACE_HLL_GPU
    implicit none
! variables for GPU kernels
    tBlock = dim3 (BlockDimX_2D, BlockDimY_2D ,1)
    grid = dim3 ( ceiling ( real ( Mloc ) / BlockDimX_2D ) , &
                  ceiling ( real ( Nloc ) / BlockDimY_2D ) , 1)
    CALL hll_kernel<<<grid,tBlock>>>(Mloc1,Nloc,SxL_d,SxR_d,PL_d,PR_d,EtaRxL_d,EtaRxR_d,P_d)
    CALL hll_kernel<<<grid,tBlock>>>(Mloc,Nloc1,SyL_d,SyR_d,QL_d,QR_d,EtaRyL_d,EtaRyR_d,Q_d)
    CALL hll_kernel<<<grid,tBlock>>>(Mloc1,Nloc,SxL_d,SxR_d,FxL_d,FxR_d,HUxL_d,HUxR_d,Fx_d)
    CALL hll_kernel<<<grid,tBlock>>>(Mloc,Nloc1,SyL_d,SyR_d,FyL_d,FyR_d,HUyL_d,HUyR_d,Fy_d)
    CALL hll_kernel<<<grid,tBlock>>>(Mloc1,Nloc,SxL_d,SxR_d,GxL_d,GxR_d,HVxL_d,HVxR_d,Gx_d)
    CALL hll_kernel<<<grid,tBlock>>>(Mloc,Nloc1,SyL_d,SyR_d,GyL_d,GyR_d,HVyL_d,HVyR_d,Gy_d)
end subroutine FLUX_AT_INTERFACE_HLL_GPU
        
attributes(global) subroutine hll_kernel(M,N,SL,SR,FL,FR,UL,UR,FOUT)
    INTEGER,INTENT(IN),VALUE::M,N
    REAL(SP),INTENT(IN),DIMENSION(M,N),DEVICE::SL,SR,FL,FR,UL,UR
    REAL(SP),INTENT(OUT),DIMENSION(M,N),DEVICE::FOUT
    real(SP) :: tSR,tSL
    integer :: i,j,tx,ty

! set local indexes
    tx = threadIdx%x
    ty = threadIdx%y
! Global thread and  block ID
    i = tx + (blockIdx%x-1)*BlockDimX_2D
    j = ty + (blockIdx%y-1)*BlockDimY_2D
!
    IF (I>=1 .AND. I<=M .AND. J>=1 .AND. J<=N) THEN
        tSR = SR(I,J)
        tSL = SL(I,J)
        IF(SL(I,J)>=ZERO) THEN
            FOUT(I,J)=FL(I,J)
        ELSEIF(SR(I,J)<=ZERO) THEN
            FOUT(I,J)=FR(I,J)
        ELSE
            FOUT(I,J)=SR(I,J)*FL(I,J)-SL(I,J)*FR(I,J)+SL(I,J)*SR(I,J)*(UR(I,J)-UL(I,J))
            IF((ABS(SR(I,J)-SL(I,J)))<SMALL)THEN
                FOUT(I,J)=FOUT(I,J)/SMALL
            ELSE
                FOUT(I,J)=FOUT(I,J)/(SR(I,J)-SL(I,J))
            ENDIF
        ENDIF
    ENDIF
end subroutine hll_kernel    
    
!TO DO by YUAN, for HIGH_ORDER not equal to FOURTH
!subroutine construction


attributes(global) subroutine DelyFun_kernel(M,N,DIN,DOUT)
    implicit none
    INTEGER,INTENT(IN),VALUE::M,N
    REAL(SP),INTENT(IN),DIMENSION(M,N),DEVICE::DIN
    REAL(SP),INTENT(OUT),DIMENSION(M,N),DEVICE::DOUT
    real(SP),dimension(BlockDimX_2D,BlockDimY_2D),shared :: &
            DIN_sh
    integer :: i,j,tx,ty
    real(SP) :: TMP1,TMP2

! set local indexes
    tx = threadIdx%x
    ty = threadIdx%y
! Global thread and  block ID
    i = tx + (blockIdx%x-1)*16
    j = ty + (blockIdx%y-1)*BlockDimY_Inner_2D
! init shared memory
    if (i>=1 .AND. i<=M .AND. j>=1 .AND. j<=N) then
        DIN_sh(tx,ty) = DIN(i,j)
    endif
    call syncthreads()

! van Leer Limiter     
    if (i>=1 .AND. i<=M .AND. j>=2 .AND. j<=N-1 .and. &
        tx>=1 .and. ty>1 .and. tx<=blockDim%x .and. ty<blockDim%y) then
# if defined (CARTESIAN)
        TMP1=(DIN_sh(tx,ty+1)-DIN_sh(tx,ty))/DY_d
        TMP2=(DIN_sh(tx,ty)-DIN_sh(tx,ty-1))/DY_d
# else
        TMP1=(DIN_sh(tx,ty+1)-DIN_sh(tx,ty))/DY_d(i,j)
        TMP2=(DIN_sh(tx,ty)-DIN_sh(tx,ty-1))/DY_d(i,j-1)
# endif
        IF((ABS(TMP1)+ABS(TMP2))<SMALL)THEN
            DOUT(I,J)=ZERO
        ELSE
            DOUT(I,J)=(TMP1*ABS(TMP2)+ABS(TMP1)*TMP2)/(ABS(TMP1)+ABS(TMP2))
        ENDIF
    endif
    if (j==1 .AND. i>=1 .AND. i<=M) then
# if defined (CARTESIAN)
      DOUT(I,J)=(DIN(I,J+1)-DIN(I,J))/DY_d
# else
      DOUT(I,J)=(DIN(I,J+1)-DIN(I,J))/DY_d(I,J)
# endif
    endif
    if (j==N .AND. i>=1 .AND. i<=M) then
# if defined (CARTESIAN)
      DOUT(I,J)=(DIN(I,J)-DIN(I,J-1))/DY_d
# else
      DOUT(I,J)=(DIN(I,J)-DIN(I,J-1))/DY_d(I,J)
# endif
    endif

end subroutine DelyFun_kernel

attributes(global) subroutine DelxFun_kernel(M,N,DIN,DOUT)
    implicit none
    INTEGER,INTENT(IN),VALUE::M,N
    REAL(SP),INTENT(IN),DIMENSION(M,N),DEVICE::DIN
    REAL(SP),INTENT(OUT),DIMENSION(M,N),DEVICE::DOUT
    real(SP),dimension(BlockDimX_2D,BlockDimY_2D),shared :: &
            DIN_sh
    integer :: i,j,tx,ty
    real(SP) :: TMP1,TMP2

! set local indexes
    tx = threadIdx%x
    ty = threadIdx%y
! Global thread and  block ID
    i = tx + (blockIdx%x-1)*BlockDimX_Inner_2D
    j = ty + (blockIdx%y-1)*16
! init shared memory
    if (i>=1 .AND. i<=M .AND. j>=1 .AND. j<=N) then
        DIN_sh(tx,ty) = DIN(i,j)
    endif
    call syncthreads()
! van Leer Limiter     
    if (i>=2 .AND. i<=M-1 .AND. j>=1 .AND. j<=N .and. &
        tx>1 .and. ty>=1 .and. tx<blockDim%x .and. ty<=blockDim%y) then
# if defined (CARTESIAN)
        TMP1=(DIN_sh(tx+1,ty)-DIN_sh(tx,ty))/DX_d
        TMP2=(DIN_sh(tx,ty)-DIN_sh(tx-1,ty))/DX_d
# else
        TMP1=(DIN_sh(tx+1,ty)-DIN_sh(tx,ty))/DX_d(I,J)
        TMP2=(DIN_sh(tx,ty)-DIN_sh(tx-1,ty))/DX_d(I-1,J)
# endif
        IF((ABS(TMP1)+ABS(TMP2))<SMALL)THEN
            DOUT(I,J)=0.0_SP
        ELSE
            DOUT(I,J)=(TMP1*ABS(TMP2)+ABS(TMP1)*TMP2)/(ABS(TMP1)+ABS(TMP2))
        ENDIF
    endif
    if (i==1 .AND. j>=1 .AND. j<=N) then
# if defined (CARTESIAN)
      DOUT(I,J)=(DIN(I+1,J)-DIN(I,J))/DX_d
# else
      DOUT(I,J)=(DIN(I+1,J)-DIN(I,J))/DX_d(I,J)
# endif
    endif
    if (i==M .AND. j>=1 .AND. j<=N) then
# if defined (CARTESIAN)
      DOUT(I,J)=(DIN(I,J)-DIN(I-1,J))/DX_d
# else
      DOUT(I,J)=(DIN(I,J)-DIN(I-1,J))/DX_d(I,J)
# endif
    endif
end subroutine DelxFun_kernel

subroutine DelxyFun_gpu
    implicit none
! variables for GPU kernels
    tBlock = dim3 (BlockDimX_2D, BlockDimY_2D ,1)
    grid = dim3 ( ceiling ( real ( Mloc ) / BlockDimX_Inner_2D ) , &
                  ceiling ( real ( Nloc ) / 16 ) , 1) 
! compute DelxFun
    CALL DelxFun_kernel<<<grid,tBlock>>>(Mloc,Nloc,Eta_d,DelxEtar_d)
    CALL DelxFun_kernel<<<grid,tBlock>>>(Mloc,Nloc,U_d,DelxU_d)
    CALL DelxFun_kernel<<<grid,tBlock>>>(Mloc,Nloc,V_d,DelxV_d)
    CALL DelxFun_kernel<<<grid,tBlock>>>(Mloc,Nloc,HU_d,DelxHU_d)
    CALL DelxFun_kernel<<<grid,tBlock>>>(Mloc,Nloc,HV_d,DelxHV_d)

    grid = dim3 ( ceiling ( real ( Mloc ) / 16 ) , &
        ceiling ( real ( Nloc ) / BlockDimY_Inner_2D ) , 1) 
! compute DelyFun_kernel
    CALL DelyFun_kernel<<<grid,tBlock>>>(Mloc,Nloc,Eta_d,DelyEtar_d)
    CALL DelyFun_kernel<<<grid,tBlock>>>(Mloc,Nloc,U_d,DelyU_d)
    CALL DelyFun_kernel<<<grid,tBlock>>>(Mloc,Nloc,V_d,DelyV_d)
    CALL DelyFun_kernel<<<grid,tBlock>>>(Mloc,Nloc,HV_d,DelyHV_d)
    CALL DelyFun_kernel<<<grid,tBlock>>>(Mloc,Nloc,HU_d,DelyHU_d)
end subroutine DelxyFun_gpu

attributes(global) subroutine wave_speed_kernel&
    (Mloc,Nloc,Mloc1,Nloc1,NGhost,&
    UL,UR,VL,VR,HxL,HxR,HyL,HyR,SxL,SxR,SyL,SyR)
    implicit none
    integer,value :: Mloc,Nloc,Mloc1,Nloc1,NGhost
    REAL(SP),INTENT(IN),DIMENSION(Mloc1,Nloc),DEVICE::UL,UR,HxL,HxR
    REAL(SP),INTENT(IN),DIMENSION(Mloc,Nloc1),DEVICE::VL,VR,HyL,HyR
    REAL(SP),INTENT(OUT),DIMENSION(Mloc1,Nloc),DEVICE::SxL,SxR
    REAL(SP),INTENT(OUT),DIMENSION(Mloc,Nloc1),DEVICE::SyL,SyR         
    REAL(SP)::SQR_PHI_L,SQR_PHI_R,SQR_PHI_S,U_S
    integer :: i,j,tx,ty

! set local indexes
    tx = threadIdx%x
    ty = threadIdx%y
! Global thread and  block ID
    i = tx + (blockIdx%x-1)*BlockDimX_2D
    j = ty + (blockIdx%y-1)*BlockDimY_2D
! x interface
    if (i>=1+Nghost .AND. i<=Mloc1-Nghost .AND. &
        j>=1+Nghost .AND. j<=Nloc-Nghost) then
        SQR_PHI_L=SQRT(GRAV*ABS(HxL(I,J)))
        SQR_PHI_R=SQRT(GRAV*ABS(HxR(I,J)))
        SQR_PHI_S=0.5*(SQR_PHI_L+SQR_PHI_R)+0.25*(UL(I,J)-UR(I,J))  
        U_S=0.5*(UL(I,J)+UR(I,J))+SQR_PHI_L-SQR_PHI_R
        SxL(I,J)=MIN(UL(I,J)-SQR_PHI_L,U_S-SQR_PHI_S)
        SxR(I,J)=MAX(UR(I,J)+SQR_PHI_R,U_S+SQR_PHI_S)
    endif

! y interface
    if (i>=1+Nghost .AND. i<=Mloc-Nghost .AND. &
        j>=1+Nghost .AND. j<=Nloc1-Nghost) then
        SQR_PHI_L=SQRT(GRAV*ABS(HyL(I,J)))
        SQR_PHI_R=SQRT(GRAV*ABS(HyR(I,J)))
        SQR_PHI_S=0.5*(SQR_PHI_L+SQR_PHI_R)+0.25*(VL(I,J)-VR(I,J))  
        U_S=0.5*(VL(I,J)+VR(I,J))+SQR_PHI_L-SQR_PHI_R
        SyL(I,J)=MIN(VL(I,J)-SQR_PHI_L,U_S-SQR_PHI_S)
        SyR(I,J)=MAX(VR(I,J)+SQR_PHI_R,U_S+SQR_PHI_S)
    endif
end subroutine wave_speed_kernel    

subroutine wave_speed_gpu
    implicit none
    
! variables for GPU kernels
    tBlock = dim3 (BlockDimX_2D, BlockDimY_2D ,1)
    grid = dim3 ( ceiling ( real ( Mloc ) / BlockDimX_2D ) , &
                  ceiling ( real ( Nloc ) / BlockDimY_2D ) , 1)
    call wave_speed_kernel<<<grid,tBlock>>>&
        (Mloc,Nloc,Mloc1,Nloc1,NGhost,&
        UxL_d,UxR_d,VyL_d,VyR_d,HxL_d,HxR_d,HyL_d,HyR_d,&
        SxL_d,SxR_d,SyL_d,SyR_d)
! ghost cells, this does not really matter
!$cuf kernel do(2) <<<*,(4,64)>>>
    DO J=1+Nghost,Nloc-Nghost
        DO I=1,Nghost
            SxL_d(I,J)=SxL_d(Nghost+1,J)
            SxR_d(I,J)=SxR_d(Nghost+1,J)
        ENDDO
    ENDDO
!$cuf kernel do(2) <<<*,(4,64)>>>
    DO J=1+Nghost,Nloc-Nghost
        DO I=Mloc1-Nghost+1,Mloc
            SxL_d(I,J)=SxL_d(Mloc1-Nghost,J)
            SxR_d(I,J)=SxR_d(Mloc1-Nghost,J)       
        ENDDO
    ENDDO
!$cuf kernel do(2) <<<*,(64,4)>>>
    DO J=1,Nghost
        DO I=1,Mloc1
            SxL_d(I,J)=SxL_d(I,Nghost+1)
            SxR_d(I,J)=SxR_d(I,Nghost+1)
        ENDDO
    ENDDO
!$cuf kernel do(2) <<<*,(64,4)>>>
    DO J=Nloc-Nghost+1,Nloc
        DO I=1,Mloc1
            SxL_d(I,J)=SxL_d(I,Nloc-Nghost)
            SxR_d(I,J)=SxR_d(I,Nloc-Nghost)
        ENDDO
    ENDDO
!$cuf kernel do(2) <<<*,(64,4)>>>
    DO J=1,Nghost
        DO I=1+Nghost,Mloc-Nghost
            SyL_d(I,J)=SyL_d(I,Nghost+1)
            SyR_d(I,J)=SyR_d(I,Nghost+1)
        ENDDO
    ENDDO
!$cuf kernel do(2) <<<*,(64,4)>>>
    DO J=Nloc1-Nghost+1,Nloc1
        DO I=1+Nghost,Mloc-Nghost
            SyL_d(I,J)=SyL_d(I,Nloc1-Nghost)
            SyR_d(I,J)=SyR_d(I,Nloc1-Nghost)       
        ENDDO
    ENDDO
!$cuf kernel do(2) <<<*,(4,64)>>>
    DO J=1,Nloc1
        DO I=1,Nghost
            SyL_d(I,J)=SyL_d(Nghost+1,J)
            SyR_d(I,J)=SyR_d(Nghost+1,J)
        ENDDO
    ENDDO
!$cuf kernel do(2) <<<*,(4,64)>>>
    DO J=1,Nloc1
        DO I=Mloc-Nghost+1,Mloc
            SyL_d(I,J)=SyL_d(Mloc-Nghost,J)
            SyR_d(I,J)=SyR_d(Mloc-Nghost,J)
        ENDDO
    ENDDO
end subroutine wave_speed_gpu

    
subroutine CONSTRUCTION_GPU
    implicit none
    integer :: istat
    
!D2D transfer
    istat = cudaMemcpy2D(MASK9u_d(1,1),Mloc1,MASK9_d(1,1),Mloc,&
        Mloc,Nloc,cudaMemcpyDeviceToDevice)
    istat = cudaMemcpy2D(MASK9v_d(1,1),Mloc,MASK9_d(1,1),Mloc,&
        Mloc,Nloc,cudaMemcpyDeviceToDevice)
    istat = cudaMemcpy2D(MASK9u_d(Mloc1,1),Mloc1,MASK9_d(Mloc,1),Mloc,&
        1,Nloc,cudaMemcpyDeviceToDevice)
    istat = cudaMemcpy2D(MASK9v_d(1,Nloc1),Mloc,MASK9_d(1,Nloc),Mloc,&
        Mloc,1,cudaMemcpyDeviceToDevice)

! construct in x-direction
# if defined (CARTESIAN)
    IF(DISPERSION)THEN
        grid = dim3 ( ceiling ( real ( Mloc ) / BlockDimX_Inner_2D ) , &
                  ceiling ( real ( Nloc ) / 16 ) , 1) 
        CALL DelxFun_kernel<<<grid,tBlock>>>(Mloc,Nloc,U4_d,DelxU4_d)
        CALL DelxFun_kernel<<<grid,tBlock>>>(Mloc,Nloc,V4_d,DelxV4_d)
    ENDIF
# else
# if defined(ZALPHA)
    IF(DISPERSION)THEN
        grid = dim3 ( ceiling ( real ( Mloc ) / BlockDimX_Inner_2D ) , &
                  ceiling ( real ( Nloc ) / 16 ) , 1) 
        CALL DelxFun_kernel<<<grid,tBlock>>>(Mloc,Nloc,U4_d,DelxU4_d)
        CALL DelxFun_kernel<<<grid,tBlock>>>(Mloc,Nloc,V4_d,DelxV4_d)
    ENDIF
# endif
# endif
! variables for GPU kernels
    tBlock = dim3 (BlockDimX_2D, BlockDimY_2D ,1)
    ! Mloc1=Mloc+1
    grid = dim3 ( ceiling ( real ( (Mloc+1) ) / BlockDimX_2D ) , &
                  ceiling ( real ( Nloc ) / BlockDimY_2D ) , 1)
    call construct_x_kernel<<<grid,tBlock,0,streamID(2)>>>&
            (Mloc,Nloc,Mloc1,U_d,DelxU_d,UxL_d,UxR_d)
    call construct_x_kernel<<<grid,tBlock,0,streamID(3)>>>&
            (Mloc,Nloc,Mloc1,V_d,DelxV_d,VxL_d,VxR_d)
    call construct_x_kernel<<<grid,tBlock,0,streamID(4)>>>&
            (Mloc,Nloc,Mloc1,HU_d,DelxHU_d,HUxL_d,HUxR_d)
    call construct_x_kernel<<<grid,tBlock,0,streamID(5)>>>&
            (Mloc,Nloc,Mloc1,HV_d,DelxHV_d,HVxL_d,HVxR_d)
    call construct_x_kernel<<<grid,tBlock,0,streamID(6)>>>&
            (Mloc,Nloc,Mloc1,Eta_d,DelxEtar_d,EtaRxL_d,EtaRxR_d)
# if defined (CARTESIAN)
    IF(DISPERSION)THEN
        call construct_x_kernel<<<grid,tBlock,0,streamID(7)>>>&
            (Mloc,Nloc,Mloc1,U4_d,DelxU4_d,U4xL_d,U4xR_d)
        call construct_x_kernel<<<grid,tBlock,0,streamID(8)>>>&
            (Mloc,Nloc,Mloc1,V4_d,DelxV4_d,V4xL_d,V4xR_d)
    ENDIF
# else
# if defined(ZALPHA)
    IF(DISPERSION)THEN
        call construct_x_kernel<<<grid,tBlock,0,streamID(7)>>>&
            (Mloc,Nloc,Mloc1,U4_d,DelxU4_d,U4xL_d,U4xR_d)
        call construct_x_kernel<<<grid,tBlock,0,streamID(8)>>>&
            (Mloc,Nloc,Mloc1,V4_d,DelxV4_d,V4xL_d,V4xR_d)
    ENDIF
# endif
# endif
   istat = cudaDeviceSynchronize()
!$cuf kernel do(2) <<<*,*>>>
    do j = 1,Nloc
        do i = 1,Mloc1 
            HxL_d(i,j)=EtaRxL_d(i,j)+Depthx_d(i,j)
            HxR_d(i,j)=EtaRxR_d(i,j)+Depthx_d(i,j)
        enddo
    enddo

# if defined(CARTESIAN)
!$cuf kernel do(2) <<<*,*>>>
    do j = 1,Nloc
        do i = 1,Mloc1 
            PL_d(i,j)=HUxL_d(i,j) &
            + Gamma1*MASK9u_d(i,j)*HxL_d(i,j)*U4xL_d(i,j)

            PR_d(i,j)=HUxR_d(i,j) &
            + Gamma1*MASK9u_d(i,j)*HxR_d(i,j)*U4xR_d(i,j)

            FxL_d(i,j)=Gamma3*PL_d(i,j)*(UxL_d(i,j)  &
            +Gamma1*MASK9u_d(i,j)*U4xL_d(i,j)) &
            +0.5*GRAV*((EtaRxL_d(i,j))*(EtaRxL_d(i,j))*Gamma3  &
            +2.0_SP*(EtaRxL_d(i,j))*(Depthx_d(i,j)))

            FxR_d(i,j)=Gamma3*PR_d(i,j)*(UxR_d(i,j)  &
            +Gamma1*MASK9u_d(i,j)*U4xR_d(i,j)) &
            +0.5*GRAV*((EtaRxR_d(i,j))*(EtaRxR_d(i,j))*Gamma3  &
            +2.0_SP*(EtaRxR_d(i,j))*(Depthx_d(i,j)))
        enddo
    enddo
# else
# if defined(ZALPHA)
!$cuf kernel do(2) <<<*,*>>>
    do j = 1,Nloc
        do i = 1,Mloc1 
            PL_d(i,j)=HUxL_d(i,j) &
            + Gamma1*MASK9u_d(i,j)*HxL_d(i,j)*U4xL_d(i,j)

            PR_d(i,j)=HUxR_d(i,j) &
            + Gamma1*MASK9u_d(i,j)*HxR_d(i,j)*U4xR_d(i,j)

            FxL_d(i,j)=Gamma3*PL_d(i,j)*(UxL_d(i,j)  &
               +Gamma1*MASK9u_d(i,j)*U4xL_d(i,j)) &
               +0.5*GRAV*((EtaRxL_d(i,j))*(EtaRxL_d(i,j))*Gamma3  &
               +2.0_SP*(EtaRxL_d(i,j))*(Depthx_d(i,j)))

            FxR_d(i,j)=Gamma3*PR_d(i,j)*(UxR_d(i,j)  &
               +Gamma1*MASK9u_d(i,j)*U4xR_d(i,j)) &
               +0.5*GRAV*((EtaRxR_d(i,j))*(EtaRxR_d(i,j))*Gamma3  &
               +2.0_SP*(EtaRxR_d(i,j))*(Depthx_d(i,j)))
        enddo
    enddo

# else
!$cuf kernel do(2) <<<*,*>>>
    do j = 1,Nloc
        do i = 1,Mloc1 
            PL_d(i,j)=HUxL_d(i,j)
            PR_d(i,j)=HUxR_d(i,j)
            FxL_d(i,j)=Gamma3*PL_d(i,j)*UxL_d(i,j) &
               +0.5*GRAV*((EtaRxL_d(i,j))*(EtaRxL_d(i,j))*Gamma3 &
               +2.0_SP*(EtaRxL_d(i,j))*(Depthx_d(i,j)))
            FxR_d(i,j)=Gamma3*PR_d(i,j)*UxR_d(i,j) &
               +0.5*GRAV*((EtaRxR_d(i,j))*(EtaRxR_d(i,j))*Gamma3 &
               +2.0_SP*(EtaRxR_d(i,j))*(Depthx_d(i,j)))
        enddo
    enddo
# endif
# endif

# if defined (CARTESIAN)

!$cuf kernel do(2) <<<*,*>>>
    do j = 1,Nloc
        do i = 1,Mloc1 
            GxL_d(i,j)=Gamma3*HxL_d(i,j)*    &
              ( UxL_d(i,j) + Gamma1*MASK9u_d(i,j)*U4xL_d(i,j) )*   &
              ( VxL_d(i,j) + Gamma1*MASK9u_d(i,j)*V4xL_d(i,j) )
               
            GxR_d(i,j)=Gamma3*HxR_d(i,j)*    &
              ( UxR_d(i,j) + Gamma1*MASK9u_d(i,j)*U4xR_d(i,j) )*   &
              ( VxR_d(i,j) + Gamma1*MASK9u_d(i,j)*V4xR_d(i,j) )
        enddo
    enddo
# else
!$cuf kernel do(2) <<<*,*>>>
    do j = 1,Nloc
        do i = 1,Mloc1 
            GxL_d(i,j)=Gamma3*HxL_d(i,j)*UxL_d(i,j)*VxL_d(i,j)
            GxR_d(i,j)=Gamma3*HxR_d(i,j)*UxR_d(i,j)*VxR_d(i,j)
        enddo
    enddo
# endif



! construct in y-direction
# if defined (CARTESIAN)
    IF(DISPERSION)THEN
        grid = dim3 ( ceiling ( real ( Mloc ) / 16 ) , &
            ceiling ( real ( Nloc ) / BlockDimY_Inner_2D ) , 1) 
        CALL DelyFun_kernel<<<grid,tBlock>>>(Mloc,Nloc,V4_d,DelyV4_d)
        CALL DelyFun_kernel<<<grid,tBlock>>>(Mloc,Nloc,U4_d,DelyU4_d)
    ENDIF
# else
# if defined(ZALPHA)
    IF(DISPERSION)THEN
        grid = dim3 ( ceiling ( real ( Mloc ) / 16 ) , &
            ceiling ( real ( Nloc ) / BlockDimY_Inner_2D ) , 1) 
        CALL DelyFun_kernel<<<grid,tBlock>>>(Mloc,Nloc,V4_d,DelyV4_d)
        CALL DelyFun_kernel<<<grid,tBlock>>>(Mloc,Nloc,U4_d,DelyU4_d)
    ENDIF
# endif
# endif
    ! Nloc1=Nloc+1
    tBlock = dim3 (BlockDimX_2D, BlockDimY_2D ,1)
    grid = dim3 ( ceiling ( real ( Mloc ) / BlockDimX_2D ) , &
                  ceiling ( real ( (Nloc+1) ) / BlockDimY_2D ) , 1)
    call construct_y_kernel<<<grid,tBlock,0,streamID(2)>>>&
            (Mloc,Nloc,Nloc1,U_d,DelyU_d,UyL_d,UyR_d)
    call construct_y_kernel<<<grid,tBlock,0,streamID(3)>>>&
            (Mloc,Nloc,Nloc1,V_d,DelyV_d,VyL_d,VyR_d)
    call construct_y_kernel<<<grid,tBlock,0,streamID(4)>>>&
            (Mloc,Nloc,Nloc1,HU_d,DelyHU_d,HUyL_d,HUyR_d)
    call construct_y_kernel<<<grid,tBlock,0,streamID(5)>>>&
            (Mloc,Nloc,Nloc1,HV_d,DelyHV_d,HVyL_d,HVyR_d)
    call construct_y_kernel<<<grid,tBlock,0,streamID(6)>>>&
            (Mloc,Nloc,Nloc1,Eta_d,DelyEtar_d,EtaRyL_d,EtaRyR_d)
# if defined (CARTESIAN)
    IF(DISPERSION)THEN
        call construct_y_kernel<<<grid,tBlock,0,streamID(7)>>>&
            (Mloc,Nloc,Nloc1,V4_d,DelyV4_d,V4yL_d,V4yR_d)
        call construct_y_kernel<<<grid,tBlock,0,streamID(8)>>>&
            (Mloc,Nloc,Nloc1,U4_d,DelyU4_d,U4yL_d,U4yR_d)
    ENDIF
# else
# if defined(ZALPHA)
    IF(DISPERSION)THEN
        call construct_y_kernel<<<grid,tBlock,0,streamID(7)>>>&
            (Mloc,Nloc,Nloc1,V4_d,DelyV4_d,V4yL_d,V4yR_d)
        call construct_y_kernel<<<grid,tBlock,0,streamID(8)>>>&
            (Mloc,Nloc,Nloc1,U4_d,DelyU4_d,U4yL_d,U4yR_d)
    ENDIF
# endif
# endif
   istat = cudaDeviceSynchronize()

!$cuf kernel do(2) <<<*,*>>>
    do j = 1,Nloc1
        do i = 1,Mloc 
            HyL_d(i,j)=EtaRyL_d(i,j)+Depthy_d(i,j)
            HyR_d(i,j)=EtaRyR_d(i,j)+Depthy_d(i,j)
        enddo
    enddo

# if defined (CARTESIAN)
!$cuf kernel do(2) <<<*,*>>>
    do j = 1,Nloc1
        do i = 1,Mloc 
            QL_d(i,j)=HVyL_d(i,j)   &
              + Gamma1*MASK9v_d(i,j)*HyL_d(i,j)*V4yL_d(i,j)

            QR_d(i,j)=HVyR_d(i,j)    &
              + Gamma1*MASK9v_d(i,j)*HyR_d(i,j)*V4yR_d(i,j)

            GyL_d(i,j)=Gamma3*QL_d(i,j)*(VyL_d(i,j)   &
               +Gamma1*MASK9v_d(i,j)*V4yL_d(i,j)) &
               +0.5*GRAV*((EtaRyL_d(i,j))*(EtaRyL_d(i,j))*Gamma3   &
               +2.0_SP*(EtaRyL_d(i,j))*(Depthy_d(i,j)))

            GyR_d(i,j)=Gamma3*QR_d(i,j)*(VyR_d(i,j)    &
               +Gamma1*MASK9v_d(i,j)*V4yR_d(i,j)) &
               +0.5*GRAV*((EtaRyR_d(i,j))*(EtaRyR_d(i,j))*Gamma3    &
               +2.0_SP*(EtaRyR_d(i,j))*(Depthy_d(i,j)))
        enddo
    enddo
# else
# if defined(ZALPHA)
!$cuf kernel do(2) <<<*,*>>>
    do j = 1,Nloc1
        do i = 1,Mloc 
            QL_d(i,j)=HVyL_d(i,j)   &
              + Gamma1*MASK9v_d(i,j)*HyL_d(i,j)*V4yL_d(i,j)

            QR_d(i,j)=HVyR_d(i,j)    &
              + Gamma1*MASK9v_d(i,j)*HyR_d(i,j)*V4yR_d(i,j)

            GyL_d(i,j)=Gamma3*QL_d(i,j)*(VyL_d(i,j)   &
               +Gamma1*MASK9v_d(i,j)*V4yL_d(i,j)) &
               +0.5*GRAV*((EtaRyL_d(i,j))*(EtaRyL_d(i,j))*Gamma3   &
               +2.0_SP*(EtaRyL_d(i,j))*(Depthy_d(i,j)))

            GyR_d(i,j)=Gamma3*QR_d(i,j)*(VyR_d(i,j)    &
               +Gamma1*MASK9v_d(i,j)*V4yR_d(i,j)) &
               +0.5*GRAV*((EtaRyR_d(i,j))*(EtaRyR_d(i,j))*Gamma3    &
               +2.0_SP*(EtaRyR_d(i,j))*(Depthy_d(i,j)))
        enddo
    enddo

# else
!$cuf kernel do(2) <<<*,*>>>
    do j = 1,Nloc1
        do i = 1,Mloc 
            QL_d(i,j)=HVyL_d(i,j)
            QR_d(i,j)=HVyR_d(i,j)
            GyL_d(i,j)=Gamma3*QL_d(i,j)*VyL_d(i,j) &
               +0.5*GRAV*((EtaRyL_d(i,j))*(EtaRyL_d(i,j))*Gamma3 & 
               +2.0_SP*(EtaRyL_d(i,j))*(Depthy_d(i,j)))
            GyR_d(i,j)=Gamma3*QR_d(i,j)*VyR_d(i,j) &
               +0.5*GRAV*((EtaRyR_d(i,j))*(EtaRyR_d(i,j))*Gamma3 &
               +2.0_SP*(EtaRyR_d(i,j))*(Depthy_d(i,j)))
        enddo
    enddo
# endif
# endif

# if defined (CARTESIAN)

!$cuf kernel do(2) <<<*,*>>>
    do j = 1,Nloc1
        do i = 1,Mloc 
            FyL_d(i,j)=Gamma3*HyL_d(i,j)*   &
             ( UyL_d(i,j) + Gamma1*MASK9v_d(i,j)*U4yL_d(i,j) )*   &
             ( VyL_d(i,j) + Gamma1*MASK9v_d(i,j)*V4yL_d(i,j) )
            
            FyR_d(i,j)=Gamma3*HyR_d(i,j)*   &
             ( UyR_d(i,j) + Gamma1*MASK9v_d(i,j)*U4yR_d(i,j) )*   &
             ( VyR_d(i,j) + Gamma1*MASK9v_d(i,j)*V4yR_d(i,j) )
        enddo
    enddo
# else
!$cuf kernel do(2) <<<*,*>>>
    do j = 1,Nloc1
        do i = 1,Mloc 
            FyL_d(i,j)=Gamma3*HyL_d(i,j)*UyL_d(i,j)*VyL_d(i,j)
            FyR_d(i,j)=Gamma3*HyR_d(i,j)*UyR_d(i,j)*VyR_d(i,j)
        enddo
    enddo
# endif

end subroutine CONSTRUCTION_GPU



subroutine CONSTRUCTION_HO_GPU
    implicit none
    integer :: istat
    !real(SP),dimension(Mloc,Nloc),device :: Din_d
    
!    if (.NOT. allocated(Din_d)) allocate(Din_d(Mloc,Nloc))
    if (.NOT. allocated(D1in_d)) allocate(D1in_d(Mloc,Nloc))
    if (.NOT. allocated(D2in_d)) allocate(D2in_d(Mloc,Nloc))
    if (.NOT. allocated(D3in_d)) allocate(D3in_d(Mloc,Nloc))
    if (.NOT. allocated(D4in_d)) allocate(D4in_d(Mloc,Nloc))
    if (.NOT. allocated(D5in_d)) allocate(D5in_d(Mloc,Nloc))
    if (.NOT. allocated(D6in_d)) allocate(D6in_d(Mloc,Nloc))
    if (.NOT. allocated(D7in_d)) allocate(D7in_d(Mloc,Nloc))
	
!D2D transfer
    istat = cudaMemcpy2D(MASK9u_d(1,1),Mloc1,MASK9_d(1,1),Mloc,&
        Mloc,Nloc,cudaMemcpyDeviceToDevice)
    istat = cudaMemcpy2D(MASK9v_d(1,1),Mloc,MASK9_d(1,1),Mloc,&
        Mloc,Nloc,cudaMemcpyDeviceToDevice)
    istat = cudaMemcpy2D(MASK9u_d(Mloc1,1),Mloc1,MASK9_d(Mloc,1),Mloc,&
        1,Nloc,cudaMemcpyDeviceToDevice)
    istat = cudaMemcpy2D(MASK9v_d(1,Nloc1),Mloc,MASK9_d(1,Nloc),Mloc,&
        Mloc,1,cudaMemcpyDeviceToDevice)

! construct in x-direction
! variables for GPU kernels
    tBlock = dim3 (BlockDimX_2D, BlockDimY_2D ,1)
    ! Mloc1=Mloc+1
    grid_step1 = dim3 ( ceiling ( real ( (Mloc+1) ) / 12 ) , &
                  ceiling ( real ( Nloc ) / 16 ) , 1)
    grid_step2 = dim3 ( ceiling ( real ( (Mloc+1) ) / BlockDimX_2D ) , &
                  ceiling ( real ( Nloc ) / BlockDimY_2D ) , 1)
    call construct_HO_x_step1_kernel<<<grid_step1,tBlock,0,streamID(2)>>>&
            (Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,U_d,D1in_d)
    call construct_HO_x_step2_kernel<<<grid_step2,tBlock,0,streamID(2)>>>&
            (Mloc,Nloc,Mloc1,Ibeg,Iend,Jbeg,Jend,U_d,D1in_d,UxL_d,UxR_d)
    call construct_HO_x_step1_kernel<<<grid_step1,tBlock,0,streamID(3)>>>&
            (Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,V_d,D2in_d)
    call construct_HO_x_step2_kernel<<<grid_step2,tBlock,0,streamID(3)>>>&
            (Mloc,Nloc,Mloc1,Ibeg,Iend,Jbeg,Jend,V_d,D2in_d,VxL_d,VxR_d)
    call construct_HO_x_step1_kernel<<<grid_step1,tBlock,0,streamID(4)>>>&
            (Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,HU_d,D3in_d)
    call construct_HO_x_step2_kernel<<<grid_step2,tBlock,0,streamID(4)>>>&
            (Mloc,Nloc,Mloc1,Ibeg,Iend,Jbeg,Jend,HU_d,D3in_d,HUxL_d,HUxR_d)
    call construct_HO_x_step1_kernel<<<grid_step1,tBlock,0,streamID(5)>>>&
            (Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,HV_d,D4in_d)
    call construct_HO_x_step2_kernel<<<grid_step2,tBlock,0,streamID(5)>>>&
            (Mloc,Nloc,Mloc1,Ibeg,Iend,Jbeg,Jend,HV_d,D4in_d,HVxL_d,HVxR_d)
    call construct_HO_x_step1_kernel<<<grid_step1,tBlock,0,streamID(6)>>>&
            (Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Eta_d,D5in_d)
    call construct_HO_x_step2_kernel<<<grid_step2,tBlock,0,streamID(6)>>>&
            (Mloc,Nloc,Mloc1,Ibeg,Iend,Jbeg,Jend,Eta_d,D5in_d,EtaRxL_d,EtaRxR_d)
# if defined (CARTESIAN)
    IF(DISPERSION)THEN
        call construct_HO_x_step1_kernel<<<grid_step1,tBlock,0,streamID(7)>>>&
                (Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,U4_d,D6in_d)
        call construct_HO_x_step2_kernel<<<grid_step2,tBlock,0,streamID(7)>>>&
                (Mloc,Nloc,Mloc1,Ibeg,Iend,Jbeg,Jend,U4_d,D6in_d,U4xL_d,U4xR_d)
        call construct_HO_x_step1_kernel<<<grid_step1,tBlock,0,streamID(8)>>>&
                (Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,V4_d,D7in_d)
        call construct_HO_x_step2_kernel<<<grid_step2,tBlock,0,streamID(8)>>>&
                (Mloc,Nloc,Mloc1,Ibeg,Iend,Jbeg,Jend,V4_d,D7in_d,V4xL_d,V4xR_d)
    ENDIF
# else
# if defined(ZALPHA)
    IF(DISPERSION)THEN
        call construct_HO_x_step1_kernel<<<grid_step1,tBlock,0,streamID(7)>>>&
                (Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,U4_d,D6in_d)
        call construct_HO_x_step2_kernel<<<grid_step2,tBlock,0,streamID(7)>>>&
                (Mloc,Nloc,Mloc1,Ibeg,Iend,Jbeg,Jend,U4_d,D6in_d,U4xL_d,U4xR_d)
        call construct_HO_x_step1_kernel<<<grid_step1,tBlock,0,streamID(8)>>>&
                (Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,V4_d,D7in_d)
        call construct_HO_x_step2_kernel<<<grid_step2,tBlock,0,streamID(8)>>>&
                (Mloc,Nloc,Mloc1,Ibeg,Iend,Jbeg,Jend,V4_d,D7in_d,V4xL_d,V4xR_d)
    ENDIF
# endif
# endif
   istat = cudaDeviceSynchronize()
   
!$cuf kernel do(2) <<<*,*>>>
    do j = 1,Nloc
        do i = 1,Mloc1 
            HxL_d(i,j)=EtaRxL_d(i,j)+Depthx_d(i,j)
            HxR_d(i,j)=EtaRxR_d(i,j)+Depthx_d(i,j)
        enddo
    enddo

# if defined (CARTESIAN)
!$cuf kernel do(2) <<<*,*>>>
    do j = 1,Nloc
        do i = 1,Mloc1 
            PL_d(i,j)=HUxL_d(i,j) &
            + Gamma1*MASK9u_d(i,j)*HxL_d(i,j)*U4xL_d(i,j)

            PR_d(i,j)=HUxR_d(i,j) &
            + Gamma1*MASK9u_d(i,j)*HxR_d(i,j)*U4xR_d(i,j)

            FxL_d(i,j)=Gamma3*PL_d(i,j)*(UxL_d(i,j)  &
            +Gamma1*MASK9u_d(i,j)*U4xL_d(i,j)) &
            +0.5*GRAV*((EtaRxL_d(i,j))*(EtaRxL_d(i,j))*Gamma3  &
            +2.0_SP*(EtaRxL_d(i,j))*(Depthx_d(i,j)))

            FxR_d(i,j)=Gamma3*PR_d(i,j)*(UxR_d(i,j)  &
            +Gamma1*MASK9u_d(i,j)*U4xR_d(i,j)) &
            +0.5*GRAV*((EtaRxR_d(i,j))*(EtaRxR_d(i,j))*Gamma3  &
            +2.0_SP*(EtaRxR_d(i,j))*(Depthx_d(i,j)))
        enddo
    enddo
# else
# if defined(ZALPHA)
!$cuf kernel do(2) <<<*,*>>>
    do j = 1,Nloc
        do i = 1,Mloc1 
            PL_d(i,j)=HUxL_d(i,j) &
            + Gamma1*MASK9u_d(i,j)*HxL_d(i,j)*U4xL_d(i,j)

            PR_d(i,j)=HUxR_d(i,j) &
            + Gamma1*MASK9u_d(i,j)*HxR_d(i,j)*U4xR_d(i,j)

            FxL_d(i,j)=Gamma3*PL_d(i,j)*(UxL_d(i,j)  &
               +Gamma1*MASK9u_d(i,j)*U4xL_d(i,j)) &
               +0.5*GRAV*((EtaRxL_d(i,j))*(EtaRxL_d(i,j))*Gamma3  &
               +2.0_SP*(EtaRxL_d(i,j))*(Depthx_d(i,j)))

            FxR_d(i,j)=Gamma3*PR_d(i,j)*(UxR_d(i,j)  &
               +Gamma1*MASK9u_d(i,j)*U4xR_d(i,j)) &
               +0.5*GRAV*((EtaRxR_d(i,j))*(EtaRxR_d(i,j))*Gamma3  &
               +2.0_SP*(EtaRxR_d(i,j))*(Depthx_d(i,j)))
        enddo
    enddo

# else
!$cuf kernel do(2) <<<*,*>>>
    do j = 1,Nloc
        do i = 1,Mloc1 
            PL_d(i,j)=HUxL_d(i,j)
            PR_d(i,j)=HUxR_d(i,j)
            FxL_d(i,j)=Gamma3*PL_d(i,j)*UxL_d(i,j) &
               +0.5*GRAV*((EtaRxL_d(i,j))*(EtaRxL_d(i,j))*Gamma3 &
               +2.0_SP*(EtaRxL_d(i,j))*(Depthx_d(i,j)))
            FxR_d(i,j)=Gamma3*PR_d(i,j)*UxR_d(i,j) &
               +0.5*GRAV*((EtaRxR_d(i,j))*(EtaRxR_d(i,j))*Gamma3 &
               +2.0_SP*(EtaRxR_d(i,j))*(Depthx_d(i,j)))
        enddo
    enddo
# endif
# endif

# if defined (CARTESIAN)

!$cuf kernel do(2) <<<*,*>>>
    do j = 1,Nloc
        do i = 1,Mloc1 
            GxL_d(i,j)=Gamma3*HxL_d(i,j)*    &
              ( UxL_d(i,j) + Gamma1*MASK9u_d(i,j)*U4xL_d(i,j) )*   &
              ( VxL_d(i,j) + Gamma1*MASK9u_d(i,j)*V4xL_d(i,j) )
               
            GxR_d(i,j)=Gamma3*HxR_d(i,j)*    &
              ( UxR_d(i,j) + Gamma1*MASK9u_d(i,j)*U4xR_d(i,j) )*   &
              ( VxR_d(i,j) + Gamma1*MASK9u_d(i,j)*V4xR_d(i,j) )
        enddo
    enddo
# else
!$cuf kernel do(2) <<<*,*>>>
    do j = 1,Nloc
        do i = 1,Mloc1 
            GxL_d(i,j)=Gamma3*HxL_d(i,j)*UxL_d(i,j)*VxL_d(i,j)
            GxR_d(i,j)=Gamma3*HxR_d(i,j)*UxR_d(i,j)*VxR_d(i,j)
        enddo
    enddo
# endif



! construct in y-direction
    ! Nloc1=Nloc+1
    grid_step1 = dim3 ( ceiling ( real ( Mloc ) / 16 ) , &
                  ceiling ( real ( (Nloc+1) ) / 12 ) , 1)
    grid_step2 = dim3 ( ceiling ( real ( Mloc ) / BlockDimX_2D ) , &
                  ceiling ( real ( (Nloc+1) ) / BlockDimY_2D ) , 1)
    call construct_HO_y_step1_kernel<<<grid_step1,tBlock,0,streamID(2)>>>&
            (Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,U_d,D1in_d)
    call construct_HO_y_step2_kernel<<<grid_step2,tBlock,0,streamID(2)>>>&
            (Mloc,Nloc,Nloc1,Ibeg,Iend,Jbeg,Jend,U_d,D1in_d,UyL_d,UyR_d)
    call construct_HO_y_step1_kernel<<<grid_step1,tBlock,0,streamID(3)>>>&
            (Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,V_d,D2in_d)
    call construct_HO_y_step2_kernel<<<grid_step2,tBlock,0,streamID(3)>>>&
            (Mloc,Nloc,Nloc1,Ibeg,Iend,Jbeg,Jend,V_d,D2in_d,VyL_d,VyR_d)
    call construct_HO_y_step1_kernel<<<grid_step1,tBlock,0,streamID(4)>>>&
            (Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,HU_d,D3in_d)
    call construct_HO_y_step2_kernel<<<grid_step2,tBlock,0,streamID(4)>>>&
            (Mloc,Nloc,Nloc1,Ibeg,Iend,Jbeg,Jend,HU_d,D3in_d,HUyL_d,HUyR_d)
    call construct_HO_y_step1_kernel<<<grid_step1,tBlock,0,streamID(5)>>>&
            (Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,HV_d,D4in_d)
    call construct_HO_y_step2_kernel<<<grid_step2,tBlock,0,streamID(5)>>>&
            (Mloc,Nloc,Nloc1,Ibeg,Iend,Jbeg,Jend,HV_d,D4in_d,HVyL_d,HVyR_d)
    call construct_HO_y_step1_kernel<<<grid_step1,tBlock,0,streamID(6)>>>&
            (Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Eta_d,D5in_d)
    call construct_HO_y_step2_kernel<<<grid_step2,tBlock,0,streamID(6)>>>&
            (Mloc,Nloc,Nloc1,Ibeg,Iend,Jbeg,Jend,Eta_d,D5in_d,EtaRyL_d,EtaRyR_d)
# if defined (CARTESIAN)
    IF(DISPERSION)THEN
        call construct_HO_y_step1_kernel<<<grid_step1,tBlock,0,streamID(7)>>>&
                (Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,U4_d,D6in_d)
        call construct_HO_y_step2_kernel<<<grid_step2,tBlock,0,streamID(7)>>>&
                (Mloc,Nloc,Nloc1,Ibeg,Iend,Jbeg,Jend,U4_d,D6in_d,U4yL_d,U4yR_d)
        call construct_HO_y_step1_kernel<<<grid_step1,tBlock,0,streamID(8)>>>&
                (Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,V4_d,D7in_d)
        call construct_HO_y_step2_kernel<<<grid_step2,tBlock,0,streamID(8)>>>&
                (Mloc,Nloc,Nloc1,Ibeg,Iend,Jbeg,Jend,V4_d,D7in_d,V4yL_d,V4yR_d)
    ENDIF
# else
# if defined(ZALPHA)
    IF(DISPERSION)THEN
        call construct_HO_y_step1_kernel<<<grid_step1,tBlock,0,streamID(7)>>>&
                (Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,U4_d,D6in_d)
        call construct_HO_y_step2_kernel<<<grid_step2,tBlock,0,streamID(7)>>>&
                (Mloc,Nloc,Nloc1,Ibeg,Iend,Jbeg,Jend,U4_d,D6in_d,U4yL_d,U4yR_d)
        call construct_HO_y_step1_kernel<<<grid_step1,tBlock,0,streamID(8)>>>&
                (Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,V4_d,D7in_d)
        call construct_HO_y_step2_kernel<<<grid_step2,tBlock,0,streamID(8)>>>&
                (Mloc,Nloc,Nloc1,Ibeg,Iend,Jbeg,Jend,V4_d,D7in_d,V4yL_d,V4yR_d)
    ENDIF
# endif 
# endif 
    istat = cudaDeviceSynchronize()
    
!$cuf kernel do(2) <<<*,*>>>
    do j = 1,Nloc1
        do i = 1,Mloc 
            HyL_d(i,j)=EtaRyL_d(i,j)+Depthy_d(i,j)
            HyR_d(i,j)=EtaRyR_d(i,j)+Depthy_d(i,j)
        enddo
    enddo

# if defined (CARTESIAN)
!$cuf kernel do(2) <<<*,*>>>
    do j = 1,Nloc1
        do i = 1,Mloc 
            QL_d(i,j)=HVyL_d(i,j)   &
              + Gamma1*MASK9v_d(i,j)*HyL_d(i,j)*V4yL_d(i,j)

            QR_d(i,j)=HVyR_d(i,j)    &
              + Gamma1*MASK9v_d(i,j)*HyR_d(i,j)*V4yR_d(i,j)

            GyL_d(i,j)=Gamma3*QL_d(i,j)*(VyL_d(i,j)   &
               +Gamma1*MASK9v_d(i,j)*V4yL_d(i,j)) &
               +0.5*GRAV*((EtaRyL_d(i,j))*(EtaRyL_d(i,j))*Gamma3   &
               +2.0_SP*(EtaRyL_d(i,j))*(Depthy_d(i,j)))

            GyR_d(i,j)=Gamma3*QR_d(i,j)*(VyR_d(i,j)    &
               +Gamma1*MASK9v_d(i,j)*V4yR_d(i,j)) &
               +0.5*GRAV*((EtaRyR_d(i,j))*(EtaRyR_d(i,j))*Gamma3    &
               +2.0_SP*(EtaRyR_d(i,j))*(Depthy_d(i,j)))
        enddo
    enddo
# else
# if defined(ZALPHA)
!$cuf kernel do(2) <<<*,*>>>
    do j = 1,Nloc1
        do i = 1,Mloc 
            QL_d(i,j)=HVyL_d(i,j)   &
              + Gamma1*MASK9v_d(i,j)*HyL_d(i,j)*V4yL_d(i,j)

            QR_d(i,j)=HVyR_d(i,j)    &
              + Gamma1*MASK9v_d(i,j)*HyR_d(i,j)*V4yR_d(i,j)

            GyL_d(i,j)=Gamma3*QL_d(i,j)*(VyL_d(i,j)   &
               +Gamma1*MASK9v_d(i,j)*V4yL_d(i,j)) &
               +0.5*GRAV*((EtaRyL_d(i,j))*(EtaRyL_d(i,j))*Gamma3   &
               +2.0_SP*(EtaRyL_d(i,j))*(Depthy_d(i,j)))

            GyR_d(i,j)=Gamma3*QR_d(i,j)*(VyR_d(i,j)    &
               +Gamma1*MASK9v_d(i,j)*V4yR_d(i,j)) &
               +0.5*GRAV*((EtaRyR_d(i,j))*(EtaRyR_d(i,j))*Gamma3    &
               +2.0_SP*(EtaRyR_d(i,j))*(Depthy_d(i,j)))
        enddo
    enddo

# else
!$cuf kernel do(2) <<<*,*>>>
    do j = 1,Nloc1
        do i = 1,Mloc 
            QL_d(i,j)=HVyL_d(i,j)
            QR_d(i,j)=HVyR_d(i,j)
            GyL_d(i,j)=Gamma3*QL_d(i,j)*VyL_d(i,j) &
               +0.5*GRAV*((EtaRyL_d(i,j))*(EtaRyL_d(i,j))*Gamma3 & 
               +2.0_SP*(EtaRyL_d(i,j))*(Depthy_d(i,j)))
            GyR_d(i,j)=Gamma3*QR_d(i,j)*VyR_d(i,j) &
               +0.5*GRAV*((EtaRyR_d(i,j))*(EtaRyR_d(i,j))*Gamma3 &
               +2.0_SP*(EtaRyR_d(i,j))*(Depthy_d(i,j)))
        enddo
    enddo
# endif
# endif

# if defined (CARTESIAN)

!$cuf kernel do(2) <<<*,*>>>
    do j = 1,Nloc1
        do i = 1,Mloc 
            FyL_d(i,j)=Gamma3*HyL_d(i,j)*   &
             ( UyL_d(i,j) + Gamma1*MASK9v_d(i,j)*U4yL_d(i,j) )*   &
             ( VyL_d(i,j) + Gamma1*MASK9v_d(i,j)*V4yL_d(i,j) )
            
            FyR_d(i,j)=Gamma3*HyR_d(i,j)*   &
             ( UyR_d(i,j) + Gamma1*MASK9v_d(i,j)*U4yR_d(i,j) )*   &
             ( VyR_d(i,j) + Gamma1*MASK9v_d(i,j)*V4yR_d(i,j) )
        enddo
    enddo
# else
!$cuf kernel do(2) <<<*,*>>>
    do j = 1,Nloc1
        do i = 1,Mloc 
            FyL_d(i,j)=Gamma3*HyL_d(i,j)*UyL_d(i,j)*VyL_d(i,j)
            FyR_d(i,j)=Gamma3*HyR_d(i,j)*UyR_d(i,j)*VyR_d(i,j)
        enddo
    enddo
# endif

end subroutine CONSTRUCTION_HO_GPU


attributes(global) subroutine construct_x_step1_kernel(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Vin,Din)
    implicit none
    integer,value :: Mloc,Nloc,Ibeg,Iend,Jbeg,Jend
    REAL(SP) :: TXP1,TXP2,TXP3,DVP1,DVP2,DVP3
    real(SP),dimension(Mloc,Nloc),device,intent(in) :: Vin
    real(SP),dimension(Mloc,Nloc),device,intent(out) :: Din
    real(SP),dimension(BlockDimX_2D,BlockDimY_2D),shared :: &
            Vin_sh
    integer :: i,j,tx,ty

! set local indexes
    tx = threadIdx%x
    ty = threadIdx%y
! Global thread and  block ID
! BlockDimX_Inner_2D should be 12?, while Y direction should be 16
    !i = tx + (blockIdx%x-1)*BlockDimX_Inner_2D  
    !j = ty + (blockIdx%y-1)*BlockDimY_Inner_2D
    i = tx + (blockIdx%x-1)*12  
    j = ty + (blockIdx%y-1)*16
! init shared memory
    if (i>=Ibeg-3 .AND. i<=Iend+3 .AND. j>=Jbeg .AND. j<=Jend) then
        Vin_sh(tx,ty) = Vin(i,j)
    endif
    call syncthreads()

    if (i>=Ibeg-1 .AND. i<=Iend+2 .AND. j>=Jbeg .AND. j<=Jend .and. &
      tx>2 .and. ty>=1 .and. tx<BLockDimX_2D .and. ty<=BlockDimY_2D)then
        TXP1=Vin_sh(tx-1,ty)-Vin_sh(tx-2,ty)
        TXP2=Vin_sh(tx,ty)-Vin_sh(tx-1,ty)
        TXP3=Vin_sh(tx+1,ty)-Vin_sh(tx,ty)
        if (TXP1.ge.0.0_SP) then
           DVP1=MAX(0.0_SP,MIN(TXP1,2.0*TXP2,2.0*TXP3))          
        else
           DVP1=MIN(0.0_SP,MAX(TXP1,2.0*TXP2,2.0*TXP3))          
        endif
        if (TXP2.ge.0.0_SP) then
           DVP2=MAX(0.0_SP,MIN(TXP2,2.0*TXP3,2.0*TXP1))          
        else
           DVP2=MIN(0.0_SP,MAX(TXP2,2.0*TXP3,2.0*TXP1))          
        endif
        if (TXP3.ge.0.0_SP) then
           DVP3=MAX(0.0_SP,MIN(TXP3,2.0*TXP1,2.0*TXP2))          
        else
           DVP3=MIN(0.0_SP,MAX(TXP3,2.0*TXP1,2.0*TXP2))          
        endif

        ! dry D-2 I-1 I D+1, lower-order
        IF(MASK_d(I-2,J)==0.OR.MASK_d(I+1,J)==0)THEN
          TXP2=Vin_sh(tx,ty)-Vin_sh(tx-1,ty)
          TXP1=TXP2
          TXP3=TXP2
          if (TXP1.ge.0.0_SP) then
             DVP1=MAX(0.0_SP,MIN(TXP1,2.0*TXP2,2.0*TXP3))
          else
             DVP1=MIN(0.0_SP,MAX(TXP1,2.0*TXP2,2.0*TXP3))
          endif
          if (TXP2.ge.0.0_SP) then
             DVP2=MAX(0.0_SP,MIN(TXP2,2.0*TXP3,2.0*TXP1))
          else
             DVP2=MIN(0.0_SP,MAX(TXP2,2.0*TXP3,2.0*TXP1))
          endif
          if (TXP3.ge.0.0_SP) then
             DVP3=MAX(0.0_SP,MIN(TXP3,2.0*TXP1,2.0*TXP2))          
          else
             DVP3=MIN(0.0_SP,MAX(TXP3,2.0*TXP1,2.0*TXP2))          
          endif
        ! here actually DVP1=DVP2=DVP3=TXP2
        ENDIF
        ! dry I-2 D-1 D I+1, zero gradient
        IF(MASK_d(I-1,J)==0.OR.MASK_d(I,J)==0)THEN
          DVP1=ZERO
          DVP2=ZERO
          DVP3=ZERO
        ENDIF
        Din(I,J)=TXP2-1.0_SP/6.0_SP*(DVP3-2.0_SP*DVP2+DVP1)
    endif
end subroutine construct_x_step1_kernel

attributes(global) subroutine construct_x_kernel(M,N,M1,Vin,Din,OutL,OutR)
    implicit none
    integer,value :: M,N,M1
    real(SP),dimension(M,N),device,intent(in) :: Din,Vin
    real(SP),dimension(M1,N),device,intent(out) :: OutL,OutR
    integer :: i,j,tx,ty

! set local indexes
    tx = threadIdx%x
    ty = threadIdx%y
! Global thread and  block ID
    i = tx + (blockIdx%x-1)*BlockDimX_2D
    j = ty + (blockIdx%y-1)*BlockDimY_2D
!
    if (i>=2 .AND. i<=M .AND. j>=1 .AND. j<=N) then

# if defined (CARTESIAN)
          OutL(I,J) = Vin(I-1,J) + 0.5_SP*DX_d*Din(I-1,J)
          OutR(I,J) = Vin(I,J) - 0.5_SP*DX_d*Din(I,J)
# else
          OutL(I,J) = Vin(I-1,J) + 0.5_SP*DX_d(I-1,J)*Din(I-1,J)
          OutR(I,J) = Vin(I,J) - 0.5_SP*DX_d(I,J)*Din(I,J)
# endif
    endif
    if (i==M1 .AND. j>=1 .AND. j<=N) then
# if defined (CARTESIAN)
        OutL(I,J) = Vin(M,J) + 0.5_SP*DX_d*Din(M,J)
# else
        OutL(I,J) = Vin(M,J) + 0.5_SP*DX_d(M,J)*Din(M,J)
# endif
        OutR(I,J) = OutL(I,J)
    endif
    if (i==1 .AND. j>=1 .AND. j<=N) then
# if defined (CARTESIAN)
        OutR(I,J) = Vin(1,J) - 0.5_SP*DX_d*Din(1,J)
# else
        OutR(I,J) = Vin(1,J) - 0.5_SP*DX_d(1,J)*Din(1,J)
# endif
        OutL(I,J) = OutR(I,J)
    endif

end subroutine construct_x_kernel

attributes(global) subroutine construct_y_kernel(M,N,N1,Vin,Din,OutL,OutR)
    implicit none
    integer,value :: M,N,N1
    real(SP),dimension(M,N),device,intent(in) :: Din,Vin
    real(SP),dimension(M,N1),device,intent(out) :: OutL,OutR
    integer :: i,j,tx,ty

! set local indexes
    tx = threadIdx%x
    ty = threadIdx%y
! Global thread and  block ID
    i = tx + (blockIdx%x-1)*BlockDimX_2D
    j = ty + (blockIdx%y-1)*BlockDimY_2D
!
    if (i>=1 .AND. i<=M .AND. j>=2 .AND. j<=N) then

# if defined (CARTESIAN)
           OutL(I,J) = Vin(I,J-1) + 0.5_SP*DY_d*Din(I,J-1)
           OutR(I,J) = Vin(I,J) - 0.5_SP*DY_d*Din(I,J)
# else
           OutL(I,J) = Vin(I,J-1) + 0.5_SP*DY_d(I,J-1)*Din(I,J-1)
           OutR(I,J) = Vin(I,J) - 0.5_SP*DY_d(I,J)*Din(I,J)
# endif
    endif
    if (j==N1 .AND. i>=1 .AND. i<=M) then
# if defined (CARTESIAN)
        OutL(I,J) = Vin(I,N) + 0.5_SP*DY_d*Din(I,N)
# else
        OutL(I,J) = Vin(I,N) + 0.5_SP*DY_d(I,N)*Din(I,N)
# endif
        OutR(I,J)=OutL(I,J)
    endif
    if (j==1 .AND. i>=1 .AND. i<=M) then
# if defined (CARTESIAN)
        OutR(I,J) = Vin(I,1) - 0.5_SP*DY_d*Din(I,1)
# else
        OutR(I,J) = Vin(I,1) - 0.5_SP*DY_d(I,1)*Din(I,1)
# endif
        OutL(I,J) = OutR(I,J)
    endif

end subroutine construct_y_kernel

attributes(global) subroutine construct_HO_x_step1_kernel(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Vin,Din)
    implicit none
    integer,value :: Mloc,Nloc,Ibeg,Iend,Jbeg,Jend
    REAL(SP) :: TXP1,TXP2,TXP3,DVP1,DVP2,DVP3
    real(SP),dimension(Mloc,Nloc),device,intent(in) :: Vin
    real(SP),dimension(Mloc,Nloc),device,intent(out) :: Din
    real(SP),dimension(BlockDimX_2D,BlockDimY_2D),shared :: &
            Vin_sh
    integer :: i,j,tx,ty

! set local indexes
    tx = threadIdx%x
    ty = threadIdx%y
! Global thread and  block ID
! BlockDimX_Inner_2D should be 12?, while Y direction should be 16
    !i = tx + (blockIdx%x-1)*BlockDimX_Inner_2D  
    !j = ty + (blockIdx%y-1)*BlockDimY_Inner_2D
    i = tx + (blockIdx%x-1)*12  
    j = ty + (blockIdx%y-1)*16
! init shared memory
    if (i>=Ibeg-3 .AND. i<=Iend+3 .AND. j>=Jbeg .AND. j<=Jend) then
        Vin_sh(tx,ty) = Vin(i,j)
    endif
    call syncthreads()

    if (i>=Ibeg-1 .AND. i<=Iend+2 .AND. j>=Jbeg .AND. j<=Jend .and. &
      tx>2 .and. ty>=1 .and. tx<BLockDimX_2D .and. ty<=BlockDimY_2D)then
        TXP1=Vin_sh(tx-1,ty)-Vin_sh(tx-2,ty)
        TXP2=Vin_sh(tx,ty)-Vin_sh(tx-1,ty)
        TXP3=Vin_sh(tx+1,ty)-Vin_sh(tx,ty)
        if (TXP1.ge.0.0_SP) then
           DVP1=MAX(0.0_SP,MIN(TXP1,2.0*TXP2,2.0*TXP3))          
        else
           DVP1=MIN(0.0_SP,MAX(TXP1,2.0*TXP2,2.0*TXP3))          
        endif
        if (TXP2.ge.0.0_SP) then
           DVP2=MAX(0.0_SP,MIN(TXP2,2.0*TXP3,2.0*TXP1))          
        else
           DVP2=MIN(0.0_SP,MAX(TXP2,2.0*TXP3,2.0*TXP1))          
        endif
        if (TXP3.ge.0.0_SP) then
           DVP3=MAX(0.0_SP,MIN(TXP3,2.0*TXP1,2.0*TXP2))          
        else
           DVP3=MIN(0.0_SP,MAX(TXP3,2.0*TXP1,2.0*TXP2))          
        endif

        ! dry D-2 I-1 I D+1, lower-order
        IF(MASK_d(I-2,J)==0.OR.MASK_d(I+1,J)==0)THEN
          TXP2=Vin_sh(tx,ty)-Vin_sh(tx-1,ty)
          TXP1=TXP2
          TXP3=TXP2
          if (TXP1.ge.0.0_SP) then
             DVP1=MAX(0.0_SP,MIN(TXP1,2.0*TXP2,2.0*TXP3))
          else
             DVP1=MIN(0.0_SP,MAX(TXP1,2.0*TXP2,2.0*TXP3))
          endif
          if (TXP2.ge.0.0_SP) then
             DVP2=MAX(0.0_SP,MIN(TXP2,2.0*TXP3,2.0*TXP1))
          else
             DVP2=MIN(0.0_SP,MAX(TXP2,2.0*TXP3,2.0*TXP1))
          endif
          if (TXP3.ge.0.0_SP) then
             DVP3=MAX(0.0_SP,MIN(TXP3,2.0*TXP1,2.0*TXP2))          
          else
             DVP3=MIN(0.0_SP,MAX(TXP3,2.0*TXP1,2.0*TXP2))          
          endif
        ! here actually DVP1=DVP2=DVP3=TXP2
        ENDIF
        ! dry I-2 D-1 D I+1, zero gradient
        IF(MASK_d(I-1,J)==0.OR.MASK_d(I,J)==0)THEN
          DVP1=ZERO
          DVP2=ZERO
          DVP3=ZERO
        ENDIF
        Din(I,J)=TXP2-1.0_SP/6.0_SP*(DVP3-2.0_SP*DVP2+DVP1)
    endif
end subroutine construct_HO_x_step1_kernel
        
attributes(global) subroutine construct_HO_x_step2_kernel(Mloc,Nloc,Mloc1,Ibeg,Iend,Jbeg,Jend,Vin,Din,OutL,OutR)
    implicit none
    integer,value :: Mloc,Nloc,Mloc1,Ibeg,Iend,Jbeg,Jend
    REAL(SP) :: VAN1, VAN2, RAT
    REAL(SP) :: TMP1,TMP2
    real(SP),dimension(Mloc,Nloc),device,intent(in) :: Din,Vin
    real(SP),dimension(Mloc1,Nloc),device,intent(out) :: OutL,OutR
    integer :: i,j,tx,ty

! set local indexes
    tx = threadIdx%x
    ty = threadIdx%y
! Global thread and  block ID
    i = tx + (blockIdx%x-1)*BlockDimX_2D
    j = ty + (blockIdx%y-1)*BlockDimY_2D
!
    if (i>=Ibeg .AND. i<=Iend+1 .AND. j>=Jbeg .AND. j<=Jend) then
        TMP1 = Din(I-1,J);   TMP2 = Din(I,J);
        IF( ABS(TMP1).le.SMALL ) TMP1 = SMALL*SIGN( 1.0_SP, TMP1 )
        IF( ABS(TMP2).le.SMALL ) TMP2 = SMALL*SIGN( 1.0_SP, TMP2 )
        RAT = TMP2/TMP1
        VAN1 = 0.0_SP;
        IF( abs(1.0 + RAT) > SMALL ) VAN1 = ( RAT + ABS(RAT) )/( 1.0_SP + RAT)

        RAT = TMP1/TMP2
        VAN2 = 0.0_SP;
        IF( abs(1.0 + RAT) > SMALL ) VAN2 = ( RAT + ABS(RAT) )/( 1.0_SP + RAT)

        OutL(I,J) = Vin(I-1,J) + ( 1.0_SP/6.0_SP )*( VAN1*TMP1 + 2.0_SP*VAN2*TMP2 )
!!!!!!!!!!!!
        TMP1 = Din(I,J);   TMP2 = Din(I+1,J);

        IF( ABS(TMP1).le.SMALL ) TMP1 = SMALL*SIGN( 1.0_SP, TMP1 )
        IF( ABS(TMP2).le.SMALL ) TMP2 = SMALL*SIGN( 1.0_SP, TMP2 )
        RAT = TMP2/TMP1
        VAN1 = 0.0_SP;
        IF( abs(1.0 + RAT) > SMALL ) VAN1 = ( RAT + ABS(RAT) )/( 1.0_SP + RAT)

        RAT = TMP1/TMP2
        VAN2 = 0.0_SP;
        IF( abs(1.0 + RAT) > SMALL ) VAN2 = ( RAT + ABS(RAT) )/( 1.0_SP + RAT)

        OutR(I,J) = Vin(I,J) - ( 1.0_SP/6.0_SP )*( 2.0_SP*VAN1*TMP1 + VAN2*TMP2 )
    endif

end subroutine construct_HO_x_step2_kernel

!subroutine CONSTRUCT_HO_Y_GPU
!    implicit none
!    real(SP),dimension(Mloc,Nloc),intent(in),device :: Vin_d
!    real(SP),dimension(Mloc1,Nloc),intent(out),device :: OutL_d,OutR_d
!    real(SP),dimension(Mloc,Nloc),device :: Din_d
!    
!    Din_d = 0.0_SP
!! variables for GPU kernels
!    tBlock = dim3 (BlockDimX_2D, BlockDimY_2D ,1)
!    grid = dim3 ( ceiling ( real ( Mloc ) / BlockDimX_Inner_2D ) , &
!                  ceiling ( real ( Nloc ) / BlockDimY_Inner_2D ) , 1)
!    call construct_HO_x_step1_kernel<<<grid,tBlock>>>&
!            (Vin_d,Din_d)
!    call cudaDeviceSynchronize()
!    call construct_HO_x_step2_kernel<<<grid,tBlock>>>&
!            ()
!end subroutine CONSTRUCT_HO_X_GPU

attributes(global) subroutine construct_HO_y_step1_kernel(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Vin,Din)
    implicit none
    integer,value :: Mloc,Nloc,Ibeg,Iend,Jbeg,Jend
    REAL(SP) :: TYP1,TYP2,TYP3,DVP1,DVP2,DVP3
    real(SP),dimension(Mloc,Nloc),device,intent(in) :: Vin
    real(SP),dimension(Mloc,Nloc),device,intent(out) :: Din
    real(SP),dimension(BlockDimX_2D,BlockDimY_2D),shared :: &
            Vin_sh
    integer :: i,j,tx,ty

! set local indexes
    tx = threadIdx%x
    ty = threadIdx%y
! Global thread and  block ID
    i = tx + (blockIdx%x-1)*16
    j = ty + (blockIdx%y-1)*12
! init shared memory
    if (i>=Ibeg .AND. i<=Iend .AND. j>=Jbeg-3 .AND. j<=Jend+3) then
        Vin_sh(tx,ty) = Vin(i,j)
    endif
    call syncthreads()

    if (i>=Ibeg .AND. i<=Iend .AND. j>=Jbeg-1 .AND. j<=Jend+2 .and. &
     tx>=1 .and. ty>2 .and. tx<=BlockDimX_2D .and. ty<BlockDimY_2D) then
        TYP1=Vin_sh(tx,ty-1)-Vin_sh(tx,ty-2)
        TYP2=Vin_sh(tx,ty)-Vin_sh(tx,ty-1)
        TYP3=Vin_sh(tx,ty+1)-Vin_sh(tx,ty)
        if (TYP1.ge.0.0_SP) then
           DVP1=MAX(0.0_SP,MIN(TYP1,2.0*TYP2,2.0*TYP3))          
        else
           DVP1=MIN(0.0_SP,MAX(TYP1,2.0*TYP2,2.0*TYP3))          
        endif
        if (TYP2.ge.0.0_SP) then
           DVP2=MAX(0.0_SP,MIN(TYP2,2.0*TYP3,2.0*TYP1))          
        else
           DVP2=MIN(0.0_SP,MAX(TYP2,2.0*TYP3,2.0*TYP1))          
        endif
        if (TYP3.ge.0.0_SP) then
           DVP3=MAX(0.0_SP,MIN(TYP3,2.0*TYP1,2.0*TYP2))          
        else
           DVP3=MIN(0.0_SP,MAX(TYP3,2.0*TYP1,2.0*TYP2))          
        endif

        ! dry D-2 J-1 J D+1, lower-order
        IF(MASK_d(I,J-2)==0.OR.MASK_d(I,J+1)==0)THEN
          TYP2=Vin(I,J)-Vin(I,J-1)
          TYP1=TYP2
          TYP3=TYP2
          if (TYP1.ge.0.0_SP) then
             DVP1=MAX(0.0_SP,MIN(TYP1,2.0*TYP2,2.0*TYP3))
          else
             DVP1=MIN(0.0_SP,MAX(TYP1,2.0*TYP2,2.0*TYP3))
          endif
          if (TYP2.ge.0.0_SP) then
             DVP2=MAX(0.0_SP,MIN(TYP2,2.0*TYP3,2.0*TYP1))
          else
             DVP2=MIN(0.0_SP,MAX(TYP2,2.0*TYP3,2.0*TYP1))
          endif
          if (TYP3.ge.0.0_SP) then
             DVP3=MAX(0.0_SP,MIN(TYP3,2.0*TYP1,2.0*TYP2))
          else
             DVP3=MIN(0.0_SP,MAX(TYP3,2.0*TYP1,2.0*TYP2))
          endif
        ! here actually DVP1=DVP2=DVP3=TYP2
        ENDIF    
        ! dry J-2 D-1 D J+1, zero gradient
        IF(MASK_d(I,J-1)==0.OR.MASK_d(I,J)==0)THEN
          DVP1=ZERO
          DVP2=ZERO
          DVP3=ZERO
        ENDIF
        Din(I,J)=TYP2-1.0_SP/6.0_SP*(DVP3-2.0_SP*DVP2+DVP1)
    endif
end subroutine construct_HO_y_step1_kernel
        
attributes(global) subroutine construct_HO_y_step2_kernel(Mloc,Nloc,Nloc1,Ibeg,Iend,Jbeg,Jend,Vin,Din,OutL,OutR)
    implicit none
    integer,value :: Mloc,Nloc,Nloc1,Ibeg,Iend,Jbeg,Jend
    REAL(SP) :: VAN1, VAN2, RAT
    REAL(SP) :: TMP1,TMP2
    real(SP),dimension(Mloc,Nloc),device,intent(in) :: Din,Vin
    real(SP),dimension(Mloc,Nloc1),device,intent(out) :: OutL,OutR
    integer :: i,j,tx,ty

! set local indexes
    tx = threadIdx%x
    ty = threadIdx%y
! Global thread and  block ID
    i = tx + (blockIdx%x-1)*BlockDimX_2D
    j = ty + (blockIdx%y-1)*BlockDimY_2D
!
    if (i>=Ibeg .AND. i<=Iend .AND. j>=Jbeg .AND. j<=Jend+1) then
        TMP1 = Din(I,J-1);   TMP2 = Din(I,J);
        IF( ABS(TMP1).le.SMALL ) TMP1 = SMALL*SIGN( 1.0_SP, TMP1 )
        IF( ABS(TMP2).le.SMALL ) TMP2 = SMALL*SIGN( 1.0_SP, TMP2 )
        RAT = TMP2/TMP1
        VAN1 = 0.0_SP;
        IF( abs(1.0 + RAT) > SMALL ) VAN1 = ( RAT + ABS(RAT) )/( 1.0_SP + RAT)

        RAT = TMP1/TMP2
        VAN2 = 0.0_SP;
        IF( abs(1.0 + RAT) > SMALL ) VAN2 = ( RAT + ABS(RAT) )/( 1.0_SP + RAT)

        OutL(I,J) = Vin(I,J-1) + ( 1.0_SP/6.0_SP )*( VAN1*TMP1 + 2.0_SP*VAN2*TMP2 )

	!!!!!!!!!!!!!!!

        TMP1 = Din(I,J);   TMP2 = Din(I,J+1);

        IF( ABS(TMP1).le.SMALL ) TMP1 = SMALL*SIGN( 1.0_SP, TMP1 )
        IF( ABS(TMP2).le.SMALL ) TMP2 = SMALL*SIGN( 1.0_SP, TMP2 )
        RAT = TMP2/TMP1
        VAN1 = 0.0_SP;
        IF( abs(1.0 + RAT) > SMALL ) VAN1 = ( RAT + ABS(RAT) )/( 1.0_SP + RAT)

        RAT = TMP1/TMP2
        VAN2 = 0.0_SP;
        IF( abs(1.0 + RAT) > SMALL ) VAN2 = ( RAT + ABS(RAT) )/( 1.0_SP + RAT)

        OutR(I,J) = Vin(I,J) - ( 1.0_SP/6.0_SP )*( 2.0_SP*VAN1*TMP1 + VAN2*TMP2 )
    endif


end subroutine construct_HO_y_step2_kernel



attributes(global) subroutine SourceTerms_WaveMaker_kernel&
        (Ibeg,Iend,Jbeg,Jend,Mloc,Nloc,iista,jjsta,WaveMakerCode,&
        Xc_WK,Width_WK,Yc_WK,Ywidth_WK,Time_ramp, &
        Tperiod,TIME,D_gen,Beta_gen,rlamda,Nfreq, &
        FreqPeak,NumWaveComp,PeakPeriod,NumFreq, &
        VISCOSITY_BREAKING,WAVEMAKER_VIS,DIFFUSION_SPONGE,&
        nu_vis)
    implicit none
    integer,value :: Ibeg,Iend,Jbeg,Jend,Mloc,Nloc,iista,jjsta
    integer,value :: WaveMakerCode,Nfreq,NumWaveComp,NumFreq
    real(SP),value :: Xc_WK,Width_WK,Yc_WK,Ywidth_WK,Time_ramp, &
        Tperiod,TIME,D_gen,Beta_gen,rlamda,FreqPeak, &
        PeakPeriod
    logical,value :: VISCOSITY_BREAKING,WAVEMAKER_VIS,DIFFUSION_SPONGE
    real(SP) :: xmk,ymk,Dxg,Dyg,WK_Source
    integer :: kd,kf
    integer :: i,j
    real,dimension(Mloc,Nloc),intent(out),device :: nu_vis

! Global thread and  block ID
    i = threadIdx%x + (blockIdx%x-1)*BlockDimX_2D 
    j = threadIdx%y + (blockIdx%y-1)*BlockDimY_2D
!
# if defined (CARTESIAN)
        DXg=DX_d
        DYg=DY_d
# else
! only for wavemaker
        DXg=DX_d(1,1)
        DYg=DY_d(1,1)
# endif
!
    if (i>=1 .AND. i<=Mloc .AND. j>=1 .AND. j<=Nloc) then

# if defined (MGPU)
        xmk=(I-Ibeg)*DXg + (iista-1)*DXg
        ymk=(J-Jbeg)*DYg + (jjsta-1)*DYg
# else
        xmk=(I-Ibeg)*DXg
        ymk=(J-Jbeg)*DYg
# endif
        IF(ABS(xmk-Xc_WK)<Width_WK.AND. &
            ABS(ymk-Yc_WK)<Ywidth_WK/2.0_SP)THEN
! WK_REG
            IF(WaveMakerCode==4)THEN          
                WaveMaker_Mass_d(I,J)=TANH(PI/(Time_ramp*Tperiod)*TIME)*D_gen &
                 *EXP(-Beta_gen*(xmk-Xc_WK)**2)&
                 *SIN(rlamda*(ymk-ZERO)-2.0_SP*PI/Tperiod*TIME)    
            ENDIF

!  added TMA_1D and JON_1D 02/24/2017
!WK_IRR,TMA_1D,JON_1D,JON_2D
            IF(WaveMakerCode==5.OR.WaveMakerCode==8  &
                .OR.WaveMakerCode==7.OR.WaveMakerCode==6)THEN
                 WK_Source=ZERO
                 DO kf=1,Nfreq
                     WK_Source=WK_Source+TANH(PI/(Time_ramp/FreqPeak)*TIME)*(Cm_d(I,J,kf) &
                            *COS(OMGN_IR_d(KF)*TIME) &
                            +Sm_d(I,J,kf)*SIN(OMGN_IR_d(KF)*TIME))
                 ENDDO
                 WaveMaker_Mass_d(I,J)=WK_Source

            ENDIF
!WK_TIME
            IF(WaveMakerCode==9)THEN
                WK_Source=ZERO
                DO kf=1,NumWaveComp
                  WK_Source=WK_Source &
                    +TANH(PI/(Time_ramp*PeakPeriod)*TIME)*D_genS_d(kf) &
                      *EXP(-Beta_genS_d(kf)*(xmk-Xc_WK)**2)&
                      *COS(2.0_SP*PI/WAVE_COMP_d(kf,1)*TIME-WAVE_COMP_d(kf,3)) 
                ENDDO
                WaveMaker_Mass_d(I,J)=WK_Source
            ENDIF
!WK_DATA2D
            IF(WaveMakerCode==10)THEN
                WK_Source=ZERO
                DO kf=1,NumFreq
                    WK_Source=WK_Source+TANH(PI/(Time_ramp/FreqPeak)*TIME)*(Cm_d(I,J,kf) &
                            *COS(OMGN2D_d(KF)*TIME) &
                            +Sm_d(I,J,kf)*SIN(OMGN2D_d(KF)*TIME))
                ENDDO
                WaveMaker_Mass_d(I,J)=WK_Source
            ENDIF

        ENDIF ! end maker domain

        nu_vis(i,j)=0.0_SP
        IF(VISCOSITY_BREAKING .OR. WAVEMAKER_VIS)THEN
            nu_vis(i,j) = nu_break_d(i,j)
        ENDIF

        IF(DIFFUSION_SPONGE)THEN
            nu_vis(i,j) = nu_vis(i,j) + nu_sponge_d(i,j)
        ENDIF
    endif
end subroutine SourceTerms_WaveMaker_kernel


# if defined (CARTESIAN)
attributes(global) subroutine SourceTerms_step1_kernel&
        (Ibeg,Iend,Jbeg,Jend,Mloc,Nloc,iista,jjsta,&
        MinDepthFrc,Gamma1,Gamma2,WaveMakerCd,&
        Xc_WK,Width_WK,Yc_WK,Ywidth_WK,Time_ramp, &
        Tperiod,TIME,D_gen,Beta_gen,rlamda,Nfreq, &
        FreqPeak,NumWaveComp,PeakPeriod,NumFreq, &
        BREAKWATER,WaveMakerCurrentBalance,FRICTION_SPONGE)
# else
attributes(global) subroutine SourceTerms_step1_kernel&
        (Ibeg,Iend,Jbeg,Jend,Mloc,Nloc,iista,jjsta,&
        MinDepthFrc,Gamma1,WaveMakerCd,&
        Xc_WK,Width_WK,Yc_WK,Ywidth_WK,Time_ramp, &
        Tperiod,TIME,D_gen,Beta_gen,rlamda,Nfreq, &
        FreqPeak,NumWaveComp,PeakPeriod,NumFreq, &
        BREAKWATER,WaveMakerCurrentBalance,FRICTION_SPONGE)
# endif
    implicit none
    integer,value :: Ibeg,Iend,Jbeg,Jend,Mloc,Nloc,iista,jjsta
    integer,value :: NumWaveComp,Nfreq,NumFreq
    real(SP),value :: MinDepthFrc,Gamma1,WaveMakerCd
# if defined (CARTESIAN)
    real(SP),value :: Gamma2
# endif
    real(SP),value :: Xc_WK,Width_WK,Yc_WK,Ywidth_WK,Time_ramp, &
        Tperiod,TIME,D_gen,Beta_gen,rlamda,&
        FreqPeak,PeakPeriod
    logical,value :: BREAKWATER,WaveMakerCurrentBalance,FRICTION_SPONGE
    real(SP) :: xmk,ymk,Dxg,Dyg,WK_Source
    integer :: i,j,tx,ty
!    real(SP),dimension(8,8),shared :: &
    real(SP),dimension(BlockDimX_2D,BlockDimY_2D),shared :: &
            U_sh,V_sh,P_sh,Q_sh
# if defined (CARTESIAN) || defined (ZALPHA)
    real(SP),dimension(BlockDimX_2D,BlockDimY_2D),shared :: U4_sh,V4_sh
# endif
#if defined (MIXING)
    !real(SP),dimension(8,8),shared :: nu_smg_sh
    real(SP),dimension(BlockDimX_2D,BlockDimY_2D),shared :: nu_smg_sh
# endif
    integer :: tMASK
# if defined (CARTESIAN)
        DXg=DX_d
        DYg=DY_d
# else
! only for wavemaker
        DXg=DX_d(1,1)
        DYg=DY_d(1,1)
# endif
! set local indexes
    tx = threadIdx%x
    ty = threadIdx%y
! Global thread and  block ID
    i = tx + (blockIdx%x-1)*BlockDimX_Inner_2D    ! Inner must be smaller than BlockDim
    j = ty + (blockIdx%y-1)*BlockDimY_Inner_2D
    !i = tx + (blockIdx%x-1)*6    ! Inner must be smaller than BlockDim
    !j = ty + (blockIdx%y-1)*6
! init shared memory
    if (i>=Ibeg-1 .AND. i<=Iend+1 .AND. j>=Jbeg-1 .AND. j<=Jend+1) then
        U_sh(tx,ty) = U_d(i,j)
        V_sh(tx,ty) = V_d(i,j)
        P_sh(tx,ty) = P_d(i,j)
        Q_sh(tx,ty) = Q_d(i,j)
# if defined (CARTESIAN) || defined (ZALPHA)
        U4_sh(tx,ty) = U4_d(i,j)
        V4_sh(tx,ty) = V4_d(i,j)
# endif
# if defined (MIXING)
        nu_smg_sh(tx,ty) = nu_smg_d(i,j)
# endif
    endif
    call syncthreads()

    if (i>=Ibeg .AND. i<=Iend .AND. j>=Jbeg .AND. j<=Jend .and. &
        tx>1 .and. ty>1 .and. tx<blockDim%x .and. ty<blockDim%y) then
        tMASK   = MASK9_d(i,j)
!
! depth gradient term
# if defined (CARTESIAN)
! second order, move the second term to left-hand side
        SourceX_d(I,J)=GRAV*(Eta_d(I,J))/DXg*(Depthx_d(I+1,J)-Depthx_d(I,J))*MASK_d(I,J) &
 ! friction
# if defined (MANNING)
                   -(GRAV*Cd_d(I,J)**2/  &
                   Max(H_d(I,J),MinDepthFrc)**(0.333333_SP))  &
                   *U_sh(tx,ty)*SQRT(U_sh(tx,ty)*U_sh(tx,ty)+V_sh(tx,ty)*V_sh(tx,ty)) &
# else
                   -Cd_d(I,J)*U_sh(tx,ty)*SQRT(U_sh(tx,ty)*U_sh(tx,ty)+V_sh(tx,ty)*V_sh(tx,ty)) &
# endif

                       ! dispersion
                       ! (h+eta)(u*\nabla V4 + V4 * \nabla u - v1pp-v2-v3)
                   + Gamma1*tMASK*Max(H_d(I,J),MinDepthFrc)*(         & 
                     U_sh(tx,ty)*0.5_SP*(U4_sh(tx+1,ty)-U4_sh(tx-1,ty))/DXg+V_sh(tx,ty)*0.5_SP*(U4_sh(tx,ty+1)-U4_sh(tx,ty-1))/DYg &
                     +U4_sh(tx,ty)*0.5_SP*(U_sh(tx+1,ty)-U_sh(tx-1,ty))/DXg+V4_sh(tx,ty)*0.5_SP*(U_sh(tx,ty+1)-U_sh(tx,ty-1))/DYg  &
                     -Gamma2*tMASK*(U1pp_d(I,J)+U2_d(I,J)+U3_d(I,J)) &
                     )    &
                        ! Ht(-V4+V1p) = div(M)*(U4-U1p)
                    +Gamma1*tMASK*((P_sh(tx+1,ty)-P_sh(tx,ty))/DXg+(Q_sh(tx,ty+1)-Q_sh(tx,ty))/DYg) &
                      *(U4_sh(tx,ty)-U1p_d(I,J)) &
! wavemaker
                   +WaveMaker_Mass_d(I,J)*U_sh(tx,ty)

        IF(FRICTION_SPONGE) THEN
          ! note that, compared to wei et al, we used flux. so need multiply D
            SourceX_d(I,J) = SourceX_d(I,J) &
                 - CD_4_SPONGE_d(I,J)*U_sh(tx,ty)*SQRT(U_sh(tx,ty)*U_sh(tx,ty)+V_sh(tx,ty)*V_sh(tx,ty)) &
                   *Depth_d(I,J)
        ENDIF

        IF(WaveMakerCurrentBalance)THEN
# if defined (MGPU)
	    xmk=(I-Ibeg)*DXg + (iista-1)*DXg
	    ymk=(J-Jbeg)*DYg + (jjsta-1)*DYg
# else
            xmk=(I-Ibeg)*DXg
            ymk=(J-Jbeg)*DYg
# endif
            IF(ABS(xmk-Xc_WK)<Width_WK.AND. &
                ABS(ymk-Yc_WK)<Ywidth_WK/2.0_SP)THEN
                SourceX_d(I,J) = SourceX_d(I,J) &
                  -WaveMakerCd*U_sh(tx,ty)*SQRT(U_sh(tx,ty)*U_sh(tx,ty)+V_sh(tx,ty)*V_sh(tx,ty))
            ENDIF
        ENDIF ! current balance

        IF(BREAKWATER) THEN
          ! note that, compared to wei et al, we used flux. so need multiply D
            SourceX_d(I,J) = SourceX_d(I,J) &
                 - CD_breakwater_d(I,J)*U_sh(tx,ty)*SQRT(U_sh(tx,ty)*U_sh(tx,ty)+V_sh(tx,ty)*V_sh(tx,ty)) &
                   *Depth_d(I,J)
        ENDIF
          
        SourceY_d(I,J)=GRAV*(Eta_d(I,J))/DYg*(Depthy_d(I,J+1)-Depthy_d(I,J))*MASK_d(I,J) &
                          ! friction
# if defined (MANNING)
                   -(GRAV*Cd_d(I,J)**2/  &
                   Max(H_d(I,J),MinDepthFrc)**(0.333333_SP))  &
                   *V_sh(tx,ty)*SQRT(U_sh(tx,ty)*U_sh(tx,ty)+V_sh(tx,ty)*V_sh(tx,ty)) &
# else
                   -Cd_d(I,J)*V_sh(tx,ty)*SQRT(U_sh(tx,ty)*U_sh(tx,ty)+V_sh(tx,ty)*V_sh(tx,ty)) &
# endif                          
		          ! dispersion
                          ! (h+eta)(u*\nabla V4 + V4 * \nabla u -v1pp-v2-v3)
                   + Gamma1*tMASK*Max(H_d(I,J),MinDepthFrc)*(         & 
                     U_sh(tx,ty)*0.5_SP*(V4_sh(tx+1,ty)-V4_sh(tx-1,ty))/DXg+V_sh(tx,ty)*0.5_SP*(V4_sh(tx,ty+1)-V4_sh(tx,ty-1))/DYg &
                     +U4_sh(tx,ty)*0.5_SP*(V_sh(tx+1,ty)-V_sh(tx-1,ty))/DXg+V4_sh(tx,ty)*0.5_SP*(V_sh(tx,ty+1)-V_sh(tx,ty-1))/DYg  &
                     -Gamma2*tMASK*(V1pp_d(I,J)+V2_d(I,J)+V3_d(I,J)) &
                     )    &
                          ! Ht(-V4+V1p) = div(Q)*(V4-V1p)
                    +Gamma1*tMASK*((P_sh(tx+1,ty)-P_sh(tx,ty))/DXg+(Q_sh(tx,ty+1)-Q_sh(tx,ty))/DYg) &
                      *(V4_sh(tx,ty)-V1p_d(I,J))  &
! wavemaker
                   +WaveMaker_Mass_d(I,J)*V_sh(tx,ty)
       IF(FRICTION_SPONGE) THEN
          ! note that, compared to wei et al, we used flux. so need multiply D
          SourceY_d(I,J) = SourceY_d(I,J) &
                   -CD_4_SPONGE_d(I,J)*V_sh(tx,ty)*SQRT(U_sh(tx,ty)*U_sh(tx,ty)+V_sh(tx,ty)*V_sh(tx,ty)) &
                    *Depth_d(I,J)             
       ENDIF

       IF(WaveMakerCurrentBalance)THEN
# if defined (MGPU)
	      xmk=(I-Ibeg)*DXg + (iista-1)*DXg
	      ymk=(J-Jbeg)*DYg + (jjsta-1)*DYg
# else
            xmk=(I-Ibeg)*DXg
            ymk=(J-Jbeg)*DYg
# endif
         IF(ABS(xmk-Xc_WK)<Width_WK.AND. &
            ABS(ymk-Yc_WK)<Ywidth_WK/2.0_SP)THEN
            SourceY_d(I,J) = SourceY_d(I,J) &
                  -WaveMakerCd*V_sh(tx,ty)*SQRT(U_sh(tx,ty)*U_sh(tx,ty)+V_sh(tx,ty)*V_sh(tx,ty))
         ENDIF
       ENDIF ! current balance


       IF(BREAKWATER) THEN
          ! note that, compared to wei et al, we used flux. so need multiply D
          SourceY_d(I,J) = SourceY_d(I,J) &
                   -CD_BREAKWATER_d(I,J)*V_sh(tx,ty)*SQRT(U_sh(tx,ty)*U_sh(tx,ty)+V_sh(tx,ty)*V_sh(tx,ty)) &
                    *Depth_d(I,J)             
       ENDIF

!# if defined (METEO)
!       IF(WindForce)THEN
!          SourceX(I,J)=SourceX(I,J)+MASK_WIND(I,J)*RHO_AW*Cdw*WindU2D(I,J) &
!                *SQRT(WindU2D(I,J)*WindU2D(I,J) &
!                      + WindV2D(I,J)*WindV2D(I,J))
!          SourceY(I,J)=SourceY(I,J)+MASK_WIND(I,J)*RHO_AW*Cdw*WindV2D(I,J) &
!                *SQRT(WindU2D(I,J)*WindU2D(I,J) &
!                      + WindV2D(I,J)*WindV2D(I,J))
!       ENDIF   ! end wind
!# endif 


# if defined (MIXING)
       SourceX_d(I,J) = SourceX_d(I,J) + 1.0_SP/DXg*( &
                   0.5_SP*(nu_smg_sh(tx,ty)+nu_smg_sh(tx+1,ty))*1.0_SP/DXg*(P_sh(tx+1,ty)-P_sh(tx,ty))  &
                 - 0.5_SP*(nu_smg_sh(tx,ty)+nu_smg_sh(tx-1,ty))*1.0_SP/DXg*(P_sh(tx,ty)-P_sh(tx-1,ty)) ) &
               + 0.5_SP/DYg*( &
                   0.5_SP*(nu_smg_sh(tx,ty+1)+nu_smg_sh(tx,ty))*1.0_SP/DYg*(P_sh(tx,ty+1)-P_sh(tx,ty)) &
                  -0.5_SP*(nu_smg_sh(tx,ty-1)+nu_smg_sh(tx,ty))*1.0_SP/DY*(P_sh(tx,ty)-P_sh(tx,ty-1))  ) &   
               + 1.0_SP/DYg*( &
                   nu_smg_sh(tx,ty+1)*0.5_SP/DXg*(Q_sh(tx+1,ty)sh(tx,ty+1)-Q_sh(tx-1,ty)sh(tx,ty+1)) &
                  -nu_smg_sh(tx,ty-1)*0.5_SP/DXg*(Q_sh(tx+1,ty)sh(tx,ty-1)-Q_sh(tx-1,ty)sh(tx,ty-1))  )

       SourceY_d(I,J) = SourceY_d(I,J) + 1.0_SP/DYg*( &
                   0.5_SP*(nu_smg_sh(tx,ty)+nu_smg_sh(tx,ty+1))*1.0_SP/DYg*(Q_sh(tx,ty+1)-Q_sh(tx,ty))  &
                 - 0.5_SP*(nu_smg_sh(tx,ty)+nu_smg_sh(tx,ty-1))*1.0_SP/DYg*(Q_sh(tx,ty)-Q_sh(tx,ty-1)) ) &
               + 0.5_SP/DXg*( &
                   0.5_SP*(nu_smg_sh(tx+1,ty)+nu_smg_sh(tx,ty))*1.0_SP/DXg*(Q_sh(tx+1,ty)-Q_sh(tx,ty)) &
                  -0.5_SP*(nu_smg_sh(tx-1,ty)+nu_smg_sh(tx,ty))*1.0_SP/DXg*(Q_sh(tx,ty)-Q_sh(tx-1,ty))  ) &   
               + 1.0_SP/DX*( &
                   nu_smg_sh(tx+1,ty)*0.5_SP/DYg*(P_sh(tx+1,ty)sh(tx,ty+1)-P_sh(tx+1,ty)sh(tx,ty-1)) &
                  -nu_smg_sh(tx-1,ty)*0.5_SP/DYg*(P_sh(tx-1,ty)sh(tx,ty+1)-P_sh(tx-1,ty)sh(tx,ty-1))  )

# endif

# else

! second order, move the second term to left-hand side

# if defined (ZALPHA)
       SourceX_d(I,J)=GRAV*(Eta_d(I,J))*SlopeX_d(I,J)*MASK_d(I,J) &
                       ! friction
# if defined (MANNING)
                   -(GRAV*Cd_d(I,J)**2/  &
                   Max(H_d(I,J),MinDepthFrc)**(0.333333_SP))  &
                   *U_sh(tx,ty)*SQRT(U_sh(tx,ty)*U_sh(tx,ty)+V_sh(tx,ty)*V_sh(tx,ty)) &
# else
                   -Cd_d(I,J)*U_sh(tx,ty)*SQRT(U_sh(tx,ty)*U_sh(tx,ty)+V_sh(tx,ty)*V_sh(tx,ty)) &
# endif
                       ! dispersion
                        ! Ht(+V1p) = div(M)*(-U1p)
                    +Gamma1*tMASK*((P_sh(tx+1,ty)-P_sh(tx,ty))/DX_d(I,J)+(Q_sh(tx,ty+1)-Q_sh(tx,ty))/DY_d(I,J)) &
                      *(U4_sh(tx,ty)-U1p_d(I,J)) &
                        ! Coriolis
                    +Coriolis_d(I,J)*0.5_SP*(Q_sh(tx,ty)+Q_sh(tx,ty+1))
          

       SourceY_d(I,J)=GRAV*(Eta_d(I,J))*SlopeY_d(I,J)*MASK_d(I,J) &
                          ! friction
# if defined (MANNING)
                   -(GRAV*Cd_d(I,J)**2/  &
                   Max(H_d(I,J),MinDepthFrc)**(0.333333_SP))  &
                   *V_sh(tx,ty)*SQRT(U_sh(tx,ty)*U_sh(tx,ty)+V_sh(tx,ty)*V_sh(tx,ty)) &
# else
                   -Cd_d(I,J)*V_sh(tx,ty)*SQRT(U_sh(tx,ty)*U_sh(tx,ty)+V_sh(tx,ty)*V_sh(tx,ty)) &
# endif
                          ! dispersion
                          ! Ht(+V1p) = div(Q)*(-V1p)
                    +Gamma1*tMASK*((P_sh(tx+1,ty)-P_sh(tx,ty))/DX_d(I,J)+(Q_sh(tx,ty+1)-Q_sh(tx,ty))/DY_d(I,J)) &
                      *(V4_sh(tx,ty)-V1p_d(I,J)) &
                        ! Coriolis
                    -Coriolis_d(I,J)*0.5_SP*(P_sh(tx,ty)+P_sh(tx+1,ty))
# else
       SourceX_d(I,J)=GRAV*(Eta_d(I,J))*SlopeX_d(I,J)*MASK_d(I,J) &
                       ! friction
# if defined (MANNING)
                   -(GRAV*Cd_d(I,J)**2/  &
                   Max(H_d(I,J),MinDepthFrc)**(0.333333_SP))  &
                   *U_sh(tx,ty)*SQRT(U_sh(tx,ty)*U_sh(tx,ty)+V_sh(tx,ty)*V_sh(tx,ty)) &
# else
                   -Cd_d(I,J)*U_sh(tx,ty)*SQRT(U_sh(tx,ty)*U_sh(tx,ty)+V_sh(tx,ty)*V_sh(tx,ty)) &
# endif
                       ! dispersion
                        ! Ht(+V1p) = div(M)*(-U1p)
                    +Gamma1*tMASK*((P_sh(tx+1,ty)-P_sh(tx,ty))/DX_d(I,J)+(Q_sh(tx,ty+1)-Q_sh(tx,ty))/DY_d(I,J)) &
                      *(-U1p_d(I,J)) &
                        ! Coriolis
                    +Coriolis_d(I,J)*0.5_SP*(Q_sh(tx,ty)+Q_sh(tx,ty+1))
          

       SourceY_d(I,J)=GRAV*(Eta_d(I,J))*SlopeY_d(I,J)*MASK_d(I,J) &
                          ! friction
# if defined (MANNING)
                   -(GRAV*Cd_d(I,J)**2/  &
                   Max(H_d(I,J),MinDepthFrc)**(0.333333_SP))  &
                   *V_sh(tx,ty)*SQRT(U_sh(tx,ty)*U_sh(tx,ty)+V_sh(tx,ty)*V_sh(tx,ty)) &
# else
                   -Cd_d(I,J)*V_sh(tx,ty)*SQRT(U_sh(tx,ty)*U_sh(tx,ty)+V_sh(tx,ty)*V_sh(tx,ty)) &
# endif
                          ! dispersion
                          ! Ht(+V1p) = div(Q)*(-V1p)
                    +Gamma1*tMASK*((P_sh(tx+1,ty)-P_sh(tx,ty))/DX_d(I,J)+(Q_sh(tx,ty+1)-Q_sh(tx,ty))/DY_d(I,J)) &
                      *(-V1p_d(I,J)) &
                        ! Coriolis
                    -Coriolis_d(I,J)*0.5_SP*(P_sh(tx,ty)+P_sh(tx+1,ty))
# endif

# endif

    endif
!# if defined (VESSEL)
!# if defined (VESSEL_PANEL_SOURCE)
!   ! no pressure term for panel source 
!# else
!      SourceX = SourceX + VesselPressureX
!      SourceY = SourceY + VesselPressureY
!# endif
!# endif
!
!# if defined (METEO)
!    IF(AirPressure)THEN
!      SourceX = SourceX + StormPressureX
!      SourceY = SourceY + StormPressureY   
!    ENDIF  
!# endif

end subroutine SourceTerms_step1_kernel


# if defined (CARTESIAN)
attributes(global) subroutine SourceTerms_step2_kernel&
        (Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,nu_vis)
    implicit none
    integer,value :: Mloc,Nloc,Ibeg,Iend,Jbeg,Jend
    real,dimension(Mloc,Nloc),device :: nu_vis
    integer :: i,j,tx,ty
    real(SP),dimension(BlockDimX_2D,BlockDimY_2D),shared :: &
            nu_vis_sh,HU_sh,HV_sh
    real(SP) :: Dxg,Dyg
    DXg=DX_d
    DYg=DY_d
! set local indexes
    tx = threadIdx%x
    ty = threadIdx%y
! Global thread and  block ID
    i = tx + (blockIdx%x-1)*BlockDimX_Inner_2D    ! Inner must be smaller than BlockDim
    j = ty + (blockIdx%y-1)*BlockDimY_Inner_2D
! init shared memory
    if (i>=Ibeg-1 .AND. i<=Iend+1 .AND. j>=Jbeg-1 .AND. j<=Jend+1) then
        nu_vis_sh(tx,ty) = nu_vis(i,j)
        HU_sh(tx,ty) = HU_d(i,j)
        HV_sh(tx,ty) = HV_d(i,j)
    endif
    call syncthreads()

! eddy viscosity breaking
    if (i>=Ibeg .AND. i<=Iend .AND. j>=Jbeg .AND. j<=Jend .and. &
        tx>1 .and. ty>1 .and. tx<blockDim%x .and. ty<blockDim%y) then
!to do by YUAN,default PQ_scheme is false, P and Q need halo 2. 
!     IF(PQ_scheme)THEN

!      it turns out P and Q are not exchanged at processor interface
!      it affects edges, make PQ_scheme=false

!       SourceX(I,J) = SourceX(I,J) + 0.5_SP/DX*( &
!                       nu_vis(I+1,J)* &
!                       1.0_SP/DX*(P(I+2,J)-P(I+1,J)) &
!                      -nu_vis(I-1,J)* &
!                       1.0_SP/DX*(P(I,J)-P(I-1,J)) ) &
!!
!                                   + 1.0_SP/DY*( &
!                       0.5_SP*(nu_vis(I,J+1)+nu_vis(I,J))* &                 
!                       0.5_SP/DY*(P(I,J+1)+P(I+1,J+1)-P(I,J)-P(I+1,J)) &
!                      -0.5_SP*(nu_vis(I,J-1)+nu_vis(I,J))* &
!                       0.5_SP/DY*(P(I,J)+P(I+1,J)-P(I,J-1)-P(I+1,J-1)) )
!
!       SourceY(I,J) = SourceY(I,J) + 0.5_SP/DY*( &
!                       nu_vis(I,J+1)* &
!                       1.0_SP/DY*(Q(I,J+2)-Q(I,J+1)) &
!                      -nu_vis(I,J-1)* &
!                       1.0_SP/DY*(Q(I,J)-Q(I,J-1)) ) &
!!
!                                   + 1.0_SP/DX*( &
!                       0.5_SP*(nu_vis(I+1,J)+nu_vis(I,J))* &                 
!                       0.5_SP/DX*(Q(I+1,J)+Q(I+1,J+1)-Q(I,J)-Q(I,J+1)) &
!                      -0.5_SP*(nu_vis(I-1,J)+nu_vis(I,J))* &
!                       0.5_SP/DX*(Q(I,J)+Q(I,J+1)-Q(I-1,J)-Q(I-1,J+1)) )
!
!     ELSE
       SourceX_d(I,J) = SourceX_d(I,J) + 0.5_SP/DXg*( &
                      (nu_vis_sh(tx+1,ty)+nu_vis_sh(tx,ty)) &
                      *1.0_SP/DXg*(HU_sh(tx+1,ty)-HU_sh(tx,ty)) &
                     -(nu_vis_sh(tx-1,ty)+nu_vis_sh(tx,ty)) &
                      *1.0_SP/DXg*(HU_sh(tx,ty)-HU_sh(tx-1,ty)) ) &
                                   + 0.5_SP/DYg*( &
                      (nu_vis_sh(tx,ty+1)+nu_vis_sh(tx,ty)) &
                      *1.0_SP/DYg*(HU_sh(tx,ty+1)-HU_sh(tx,ty)) &
                     -(nu_vis_sh(tx,ty-1)+nu_vis_sh(tx,ty)) &
                      *1.0_SP/DYg*(HU_sh(tx,ty)-HU_sh(tx,ty-1)) )

       SourceY_d(I,J) = SourceY_d(I,J) + 0.5_SP/DXg*( &
                      (nu_vis_sh(tx+1,ty)+nu_vis_sh(tx,ty)) &
                      *1.0_SP/DXg*(HV_sh(tx+1,ty)-HV_sh(tx,ty)) &
                     -(nu_vis_sh(tx-1,ty)+nu_vis_sh(tx,ty)) &
                      *1.0_SP/DXg*(HV_sh(tx,ty)-HV_sh(tx-1,ty)) ) &
                                   + 0.5_SP/DYg*( &
                      (nu_vis_sh(tx,ty+1)+nu_vis_sh(tx,ty)) &
                      *1.0_SP/DYg*(HV_sh(tx,ty+1)-HV_sh(tx,ty)) &
                     -(nu_vis_sh(tx,ty-1)+nu_vis_sh(tx,ty)) &
                      *1.0_SP/DYg*(HV_sh(tx,ty)-HV_sh(tx,ty-1)) )

!     ENDIF ! end pq_scheme
    endif

end subroutine SourceTerms_step2_kernel
# endif


subroutine SourceTerms_GPU
# if defined (VESSEL)
    USE VESSEL_MODULE
# endif
# if defined (METEO)
    USE METEO_MODULE
# endif

    implicit none
    integer :: istat
    real,dimension(Mloc,Nloc),device :: nu_vis
    real(SP) :: xmk,ymk,WK_Source

! variables for GPU kernels
    tBlock = dim3 (BlockDimX_2D, BlockDimY_2D ,1)
    grid = dim3 ( ceiling ( real ( Mloc ) / BlockDimX_2D ) , &
                  ceiling ( real ( Nloc ) / BlockDimY_2D ) , 1)
    call SourceTerms_WaveMaker_kernel<<<grid, tBlock>>>&
        (Ibeg,Iend,Jbeg,Jend,Mloc,Nloc,iista,jjsta,WaveMakerCode,&
        Xc_WK,Width_WK,Yc_WK,Ywidth_WK,Time_ramp, &
        Tperiod,TIME,D_gen,Beta_gen,rlamda,Nfreq, &
        FreqPeak,NumWaveComp,PeakPeriod,NumFreq, &
        VISCOSITY_BREAKING,WAVEMAKER_VIS,DIFFUSION_SPONGE,&
        nu_vis)
    grid = dim3 ( ceiling ( real ( Mloc ) / BlockDimX_Inner_2D ) , &
                  ceiling ( real ( Nloc ) / BlockDimY_Inner_2D ) , 1)
    ! YUAN: test smaller block size, sourceterm_step1 is too massive for a cuda
    ! core.
# if defined (CARTESIAN)
    call SourceTerms_step1_kernel<<<grid, tBlock>>>&
        (Ibeg,Iend,Jbeg,Jend,Mloc,Nloc,iista,jjsta,&
        MinDepthFrc,Gamma1,Gamma2,WaveMakerCd,&
        Xc_WK,Width_WK,Yc_WK,Ywidth_WK,Time_ramp, &
        Tperiod,TIME,D_gen,Beta_gen,rlamda,Nfreq, &
        FreqPeak,NumWaveComp,PeakPeriod,NumFreq, &
        BREAKWATER,WaveMakerCurrentBalance,FRICTION_SPONGE)
# else
    call SourceTerms_step1_kernel<<<grid, tBlock>>>&
        (Ibeg,Iend,Jbeg,Jend,Mloc,Nloc,iista,jjsta,&
        MinDepthFrc,Gamma1,WaveMakerCd,&
        Xc_WK,Width_WK,Yc_WK,Ywidth_WK,Time_ramp, &
        Tperiod,TIME,D_gen,Beta_gen,rlamda,Nfreq, &
        FreqPeak,NumWaveComp,PeakPeriod,NumFreq, &
        BREAKWATER,WaveMakerCurrentBalance,FRICTION_SPONGE)
# endif
# if defined (CARTESIAN)
    IF(VISCOSITY_BREAKING.OR.DIFFUSION_SPONGE.OR.WAVEMAKER_VIS)THEN
        call SourceTerms_step2_kernel<<<grid, tBlock>>>&
            (Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,nu_vis)
    ENDIF
# endif
! ----
end subroutine SourceTerms_GPU

end module DispersionFluxesSources
