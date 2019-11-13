!------------------------------------------------------------------------------------
!
!      FILE etauv_solver_gpu.F
!
!-------------------------------------------------------------------------------------
!
!    ESTIMATE_HUV is subroutine to calculate eta, ubar and vbar
!      using 3rd-order LK scheme 
!
!    HISTORY: 
!      05/12/2011 Fengyan Shi
!
!-------------------------------------------------------------------------------------
module etauv_solver_module
    use cudafor
    use cusparse
    use PARAM
    use GLOBAL,only:ISTAGE,MGlob,NGlob,NGhost,Ibeg,Iend,Jbeg,Jend,iista,jjsta,&
        Mloc,Nloc,DISPERSION,PERIODIC,ETA_LIMITER,DISP_TIME_LEFT,&
        Gamma1,Gamma2,Gamma3,MinDepth,MinDepthFrc,FroudeCap,&
        TroughLimit,CrestLimit,Xc_WK,Yc_WK,Width_WK,Ywidth_WK,&
        DT,b1,b2,&
        n_east, n_west, n_suth, n_nrth,comm2d, ier,myid,PX,PY,&
        NumberProcessor,ProcessorID
# if defined (MGPU)
    use mgpu_utilities
# endif
    use mod_cuda, only:alpha_d, beta_d, &
        DX_d, DY_d,WaveMakerCode,&
        n_left=>n_left_d,n_right=>n_right_d,n_bottom=>n_bottom_d,n_top=>n_top_d,&
        U_d,V_d,HU_d,HV_d,Ubar_d,Vbar_d,MASK_d,MASK9_d,H_d,Depth_d,&
        EtaOld_d,UbarOld_d,VbarOld_d,Eta_d,Eta0_d,Ubar0_d,Vbar0_d,WaveMaker_Mass_d,&
        Lat_theta_d,U4_d,V4_d,P_d,Q_d,Fx_d,Fy_d,Gx_d,Gy_d,Vxy_d,DVxy_d,&
        Uxy_d,DUxy_d,Ux_d,DUx_d,ETAy_d,SourceX_d,SourceY_d,&
        grid,tBlock,tBlock_tran, grid_tran,grid_tran_back,&
        BlockDimX_2D,BlockDimY_2D,BlockDimX_Inner_2D,BlockDimY_Inner_2D,&
        istat,startEvent,stopEvent,KernelTime,TotalKernelTime,t1,t2,cusparseh,TriSolver,&
        streamID
# if defined (VESSEL_PANEL_SOURCE)
        USE VESSEL_MODULE
# endif
    real(SP), dimension(:,:),allocatable, device :: Y1_d, Y2_d
    real(SP), dimension(:), allocatable, device :: A_beg_d, C_end_d,Y1_end_d,Y2_end_d
    real(SP), dimension(:), allocatable :: A_beg, C_end,Y1_end,Y2_end
# if defined (MGPU)
    real(SP), dimension(:,:),allocatable,pinned :: rmsgx, smsgx
    real(SP), dimension(:,:),allocatable,pinned :: rmsgy, smsgy
    real(SP), dimension(:,:),allocatable,pinned :: rmsgy2, smsgy2
# endif
!
! device variables for tridiagonal solvers
    real(SP), dimension(:,:),ALLOCATABLE, device :: myAx_d, myCx_d, myDx_d, myFx_d
    real(SP), dimension(:,:),ALLOCATABLE, device :: myAy_d, myCy_d, myC1y_d, myDy_d,myD2y_d, myFy_d
    real(SP), dimension(:,:),ALLOCATABLE, device :: myAx_T, myCx_T, myDx_T, myFx_T
    real(SP), dimension(:),ALLOCATABLE, device :: Arow,Brow,Crow,Drow,D1row

        Logical :: Allocated_1 = .False., Allocated_2 = .False.
!interface 
!    integer(4) function cusparseSgtsv2StridedBatch(handle, m, dl, d, du, x, batchCount, batchStride,pBuffer) bind(C,name='cusparseSgtsv2StridedBatch')
!        use iso_c_binding
!        use cusparse
!        implicit none
!        !!pgi$ ignore_tkr handle, m, dl, d, du, x, batchCount, batchStride,pBuffer
!        type(cusparseHandle),value :: handle
!        integer(4),value :: m, batchCount, batchStride
!        real(4),dimension(*), device :: dl, d, du, x 
!        type(c_ptr), value :: pBuffer
!        !character(c_char),device ::pBuffer(*)
!    end function cusparseSgtsv2StridedBatch
!    integer(4) function cusparseSgtsv2StridedBatch_bufferSizeExt(handle, m, dl, d, du, x, batchCount, batchStride,pBufferSize) bind(C,name='cusparseSgtsv2StridedBatch_bufferSizeExt')
!        use iso_c_binding
!        use cusparse
!        implicit none
!        !!pgi$ ignore_tkr handle, m, dl, d, du, x, batchCount, batchStride
!        type(cusparseHandle),value :: handle
!        integer(4),value :: m, batchCount, batchStride
!        real(4),dimension(*), device :: dl, d, du, x 
!        !integer(c_size_t),value :: pBufferSize
!        integer(8),value :: pBufferSize
!    end function cusparseSgtsv2StridedBatch_bufferSizeExt
!end interface



contains
    subroutine ESTIMATE_HUV_GPU(ISTEP)
        implicit none

        integer :: ISTEP
        tBlock = dim3 (BlockDimX_2D, BlockDimY_2D ,1)
        grid = dim3 ( ceiling ( real (Mloc) / BlockDimX_2D ) , &
                ceiling ( real (Nloc) / BlockDimY_2D ) , 1)
        
        call HUV_BAR_KERNEL <<<grid ,tBlock>>> &
            (Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,iista,jjsta,ISTEP,WaveMakerCode,DT,&
            Xc_WK, Yc_WK,Width_WK,Ywidth_WK,Gamma1,Gamma3 &
# if !defined (ITERATION)
            ,ETA_LIMITER, TroughLimit, CrestLimit &
# endif
            )
!
        IF(DISPERSION)THEN
! device variables for tridiagonal solvers
            if (.NOT. Allocated_1)then
            if (TriSolver==0) then
                Allocated_1 = .True.
                allocate( myAx_d(Mloc,Nloc),myCx_d(Mloc,Nloc),&
                          myDx_d(Mloc,Nloc),myFx_d(Mloc,Nloc),myAx_T(Nloc,Mloc),&
                          myCx_T(Nloc,Mloc),myDx_T(Nloc,Mloc),myFx_T(Nloc,Mloc),&
                          myAy_d(Mloc,Nloc),myCy_d(Mloc,Nloc),myC1y_d(Mloc,Nloc),&
                          myDy_d(Mloc,Nloc),myD2y_d(Mloc,Nloc),myFy_d(Mloc,Nloc) )
            else if (TriSolver==1) then
                allocate( Arow((Mloc-2*NGhost)*(Nloc-2*NGhost)))
                allocate( Brow((Mloc-2*NGhost)*(Nloc-2*NGhost)))
                Brow = 1.0_SP
!C1row is not neritten, but D1row is necessary
                allocate( Crow((Mloc-2*NGhost)*(Nloc-2*NGhost)),&
                 D1row((Mloc-2*NGhost)*(Nloc-2*NGhost)),&
                 Drow((Mloc-2*NGhost)*(Nloc-2*NGhost)))
            endif
            endif
        ! prepare myA-D for tridiagonal solver
        ! variables for GPU kernels
            tBlock = dim3 (BlockDimX_2D, BlockDimY_2D ,1)
            grid = dim3 ( ceiling ( real (Mloc) / BlockDimX_2D ) , &
                  ceiling ( real (Nloc) / BlockDimY_2D ) , 1)
            if (TriSolver==0) then
                call triDx_init_kernel<<<grid,tBlock>>>&
                    (Ibeg,Iend,Jbeg,Jend,n_west,b1,b2,Gamma1,MinDepthFrc,WaveMakerCode)
            else if (TriSolver==1) then
                call triDx_cusparse_init_kernel<<<grid,tBlock>>>&
                    (Ibeg,Iend,Jbeg,Jend,Mloc,n_west,NGhost,b1,b2,Gamma1,MinDepthFrc,WaveMakerCode)
            endif
        ! prepare myA-D for tridiagonal solver in y direction
            tBlock = dim3 (BlockDimX_2D, BlockDimY_2D ,1)
            grid = dim3 ( ceiling ( real (Mloc) / BlockDimX_2D ) , &
                  ceiling ( real (Nloc) / BlockDimY_2D ) , 1)
            if (TriSolver==0) then
                call triDy_init_kernel<<<grid,tBlock>>>&
                    (Ibeg,Iend,Jbeg,Jend,b1,b2,Gamma1,Gamma2,MinDepthFrc,WaveMakerCode,DISP_TIME_LEFT)
            else if (TriSolver==1) then
                call triDy_cusparse_init_kernel<<<grid,tBlock>>>&
                    (Ibeg,Iend,Jbeg,Jend,Nloc,NGhost,b1,b2,Gamma1,Gamma2,MinDepthFrc,WaveMakerCode,DISP_TIME_LEFT,PERIODIC)
            endif
            
            if (TriSolver==0) then
                if(PERIODIC) then
                    call periodic_triDx_triDy_cuda
                else
                    call triDx_triDy_cuda
                endif
            else if (TriSolver==1) then
! TO DO BY YUAN: rewrite the subroutines of cusparse solvers for non- and periodic condition
! However, since cudagtsv is not fast as thought, and can not be accelerated by
! Multiple GPUs, so leave it to the future update.
! So, stick to Trisolver = 0, donot try 1 at current version
                call triDx_cusparse
                call triDy_cusparse
            endif
            
            
        ELSE  ! if no dispersion
!$cuf kernel do(2) <<<*, *>>>
            DO J=Jbeg,Jend
            DO I=Ibeg,Iend  
                U_d(I,J)=Ubar_d(I,J)/Max(H_d(I,J),MinDepthFrc)
                V_d(I,J)=Vbar_d(I,J)/Max(H_d(I,J),MinDepthFrc)
            ENDDO
            ENDDO   

        ENDIF  ! end dispersion
            
        tBlock = dim3 (BlockDimX_2D, BlockDimY_2D ,1)
        grid = dim3 ( ceiling ( real (Mloc) / BlockDimX_2D ) , &
              ceiling ( real (Nloc) / BlockDimY_2D ) , 1)
        call post_tridiagonal_kernel<<<grid,tBlock>>>&
            (Ibeg,Iend,Jbeg,Jend,Mloc,Nloc,MinDepthFrc,FroudeCap,DISPERSION) 
     
END SUBROUTINE ESTIMATE_HUV_GPU


# if defined (MGPU)
subroutine triDxDy_mgpu_cuda_v1
!Multiple GPU tridiagonal solver for non-periodic boundary condition
    implicit none
    integer  IbegTridx, IendTridx, JbegTridy,JendTridy
    INTEGER   len
    integer status(MPI_STATUS_SIZE), req


    tBlock = dim3(32,1,1)
    if (.NOT. allocated(rmsgx)) then
        allocate( rmsgx(Nloc,2),smsgx(Nloc,2),rmsgy(Mloc,2),smsgy(Mloc,2))
    endif
!================================
!X
! forward sweep

    if ( n_west .ne. MPI_PROC_NULL ) then
       call MPI_IRECV( rmsgx, 2*Nloc, MPI_SP,&
            n_west, 0, comm2d, req, ier )
       call MPI_WAIT( req, status, ier )
       IbegTridx = Ibeg
       istat = cudaMemcpy2dAsync(myDx_T(1,Ibeg-1),Nloc,rmsgx(1,1),Nloc,Nloc,1,cudaMemcpyHostToDevice,streamID(2))
       istat = cudaMemcpy2dAsync(myCx_T(1,Ibeg-1),Nloc,rmsgx(1,2),Nloc,Nloc,1,cudaMemcpyHostToDevice,streamID(2))
    else
       IbegTridx = Ibeg+1
    endif   
!==================================
!Y
! forward sweep
    if ( n_suth .ne. MPI_PROC_NULL ) then
       call MPI_IRECV( rmsgy, 2*Mloc, MPI_SP,&
            n_suth, 1, comm2d, req, ier )
       call MPI_WAIT(req, status, ier )
       JbegTridy = Jbeg
       istat =  cudaMemcpy2dAsync(myDy_d(1,Jbeg-1),Mloc,rmsgy(1,1),Mloc,Mloc,1,cudaMemcpyHostToDevice,streamID(3))
       istat =  cudaMemcpy2dAsync(myCy_d(1,Jbeg-1),Mloc,rmsgy(1,2),Mloc,Mloc,1,cudaMemcpyHostToDevice,streamID(3))
    else
       JbegTridy = Jbeg+1
    endif   

!=================================
!X
    call triDx_forwardsweep_cuda_kernel<<<dim3(ceiling(real(Nloc)/tBlock%x), 1, 1),tBlock,0,streamID(2)>>>&
        ( myAx_T,myCx_T,myDx_T,Mloc, Nloc,IbegTridx,Iend,Jbeg,Jend)
!==================================
!Y
    call triDy_forwardsweep_cuda_kernel<<<dim3(ceiling(real(Mloc)/tBlock%x), 1, 1),tBlock,0,streamID(3)>>>&
        ( myAy_d,myCy_d,myDy_d,Mloc, Nloc,Ibeg,Iend,JbegTridy,Jend)
!
        !istat = cudaDeviceSynchronize()
!=================================
!X
    if ( n_east .ne. MPI_PROC_NULL ) then
        istat = cudaMemcpy2DAsync(smsgx(1,1),Nloc,myDx_T(1,Iend),Nloc,&
            Nloc,1,cudaMemcpyDeviceToHost,streamID(2))
        istat = cudaMemcpy2DAsync(smsgx(1,2),Nloc,myCx_T(1,Iend),Nloc,&
            Nloc,1,cudaMemcpyDeviceToHost,streamID(2))
        istat = cudaStreamSynchronize(streamID(2))
        call MPI_ISEND( smsgx, 2*Nloc, MPI_SP,&
            n_east, 0, comm2d, req, ier )
        call MPI_WAIT( req, status, ier )
    endif
!==================================
!Y
    if ( n_nrth .ne. MPI_PROC_NULL ) then
        istat = cudaMemcpy2dAsync(smsgy(1,1),Mloc,myDy_d(1,Jend),Mloc,&
            Mloc,1,cudaMemcpyDeviceToHost,streamID(3))
        istat = cudaMemcpy2dAsync(smsgy(1,2),Mloc,myCy_d(1,Jend),Mloc,&
            Mloc,1,cudaMemcpyDeviceToHost,streamID(3))
        istat = cudaStreamSynchronize(streamID(3))
       call MPI_ISEND( smsgy, 2*Mloc, MPI_SP,&
            n_nrth, 1, comm2d, req, ier )
       call MPI_WAIT(req, status, ier )
    endif

!=================================
!X
! back substitution
    if ( n_east .ne. MPI_PROC_NULL ) then
       call MPI_IRECV( rmsgx, 2*Nloc, MPI_SP,&
            n_east, 3, comm2d, req, ier )
       call MPI_WAIT( req, status, ier )
       IendTridx = Iend
       istat = cudaMemcpy2dAsync(myFx_T(1,Iend+1),Nloc,rmsgx(1,1),Nloc,Nloc,1,cudaMemcpyHostToDevice,streamID(2))
    else
       IendTridx = Iend-1
       istat=cudaMemcpy2dAsync(myFx_T(1,Iend),Nloc,myDx_T(1,Iend),Nloc,Nloc,1,cudaMemcpyHostToDevice,streamID(2))
    endif   
!==================================
!Y
    if ( n_nrth .ne. MPI_PROC_NULL ) then
       call MPI_IRECV( rmsgy, 2*Mloc, MPI_SP,&
            n_nrth, 4, comm2d, req, ier )
       call MPI_WAIT(req, status, ier )
       JendTridy = Jend
       istat = cudaMemcpy2dAsync(myFy_d(1,Jend+1),Mloc,rmsgy(1,1),Mloc,Mloc,1,cudaMemcpyHostToDevice,streamID(3))
    else
       JendTridy = Jend-1
       istat=cudaMemcpy2dAsync(myFy_d(1,Jend),Mloc,myDy_d(1,Jend),Mloc,Mloc,1,cudaMemcpyHostToDevice,streamID(3))
    endif   
!
!=================================
!X
    call triDx_backsubstitution_cuda_kernel<<<dim3(ceiling(real(Nloc)/tBlock%x), 1, 1),tBlock,0,streamID(2)>>>&
        ( myCx_T,myDx_T,myFx_T,Mloc, Nloc,Ibeg,IendTridx,Jbeg,Jend)
!==================================
!Y
    call triDy_backsubstitution_cuda_kernel<<<dim3(ceiling(real(Mloc)/tBlock%x), 1, 1),tBlock,0,streamID(3)>>>&
        ( myCy_d,myDy_d,myFy_d,Mloc, Nloc,Ibeg,Iend,Jbeg,JendTridy)
!
    !    istat = cudaDeviceSynchronize()
!=================================
!X
    if ( n_west .ne. MPI_PROC_NULL ) then
        istat = cudaMemcpy2DAsync(smsgx(1,1),Nloc,myFx_T(1,Ibeg),Nloc,&
            Nloc,1,cudaMemcpyDeviceToHost,streamID(2))
        istat = cudaStreamSynchronize(streamID(2))
       call MPI_ISEND( smsgx, 2*Nloc, MPI_SP,&
            n_west, 3, comm2d, req, ier )
       call MPI_WAIT( req, status, ier )
    endif
!==================================
!Y
    if ( n_suth .ne. MPI_PROC_NULL ) then
        istat = cudaMemcpy2dAsync(smsgy(1,1),Mloc,Y1_d(1,Jbeg),Mloc,&
            Mloc,1,cudaMemcpyDeviceToHost,streamID(3))
        istat = cudaStreamSynchronize(streamID(3))
       call MPI_ISEND( smsgy, 2*Mloc, MPI_SP,&
            n_suth, 4, comm2d, req, ier )
       call MPI_WAIT( req, status, ier )
    endif
!=========================

end subroutine triDxDy_mgpu_cuda_v1
# endif


# if defined (MGPU)
subroutine triDxDy_mgpu_cuda_v2
!Version 2
    implicit none
    !real(SP), dimension(Nloc,2),device :: rmsg_d, smsg_d
    integer  IbegTridx, IendTridx, JbegTridy,JendTridy
    INTEGER   len
    integer status(MPI_STATUS_SIZE), req
    integer,dimension(2) :: req2
    integer,dimension(MPI_STATUS_SIZE,2) :: status2


    tBlock = dim3(32,1,1)
    if (.NOT. allocated(rmsgx)) then
        allocate( rmsgx(Nloc,2),smsgx(Nloc,2),rmsgy(Mloc,2),smsgy(Mloc,2),&
        rmsgy2(Mloc,2),smsgy2(Mloc,2))
    endif
!================================
!X
! forward sweep

    if ( n_west .ne. MPI_PROC_NULL ) then
       call MPI_IRECV( rmsgx, 2*Nloc, MPI_SP,&
            n_west, 0, comm2d, req, ier )
       call MPI_WAIT( req, status, ier )
       IbegTridx = Ibeg
       istat = cudaMemcpy2dAsync(myDx_T(1,Ibeg-1),Nloc,rmsgx(1,1),Nloc,Nloc,1,cudaMemcpyHostToDevice,streamID(2))
       istat = cudaMemcpy2dAsync(myCx_T(1,Ibeg-1),Nloc,rmsgx(1,2),Nloc,Nloc,1,cudaMemcpyHostToDevice,streamID(2))
    else
       IbegTridx = Ibeg+1
    endif   
!==================================
!Y1&Y2
! forward sweep
    if ( n_suth .ne. MPI_PROC_NULL ) then
       call MPI_IRECV( rmsgy, 2*Mloc, MPI_SP,&
            n_suth, 1, comm2d, req2(1), ier )
       call MPI_IRECV( rmsgy2, 2*Mloc, MPI_SP,&
            n_suth, 2, comm2d, req2(2), ier )
       call MPI_WAITALL(2, req2, status2, ier )
       JbegTridy = Jbeg
       istat =  cudaMemcpy2dAsync(myDy_d(1,Jbeg-1),Mloc,rmsgy(1,1),Mloc,Mloc,1,cudaMemcpyHostToDevice,streamID(3))
       istat = cudaMemcpy2dAsync(myC1y_d(1,Jbeg-1),Mloc,rmsgy(1,2),Mloc,Mloc,1,cudaMemcpyHostToDevice,streamID(3))
       istat =  cudaMemcpy2dAsync(myD2y_d(1,Jbeg-1),Mloc,rmsgy2(1,1),Mloc,Mloc,1,cudaMemcpyHostToDevice,streamID(4))
       istat = cudaMemcpy2dAsync(myCy_d(1,Jbeg-1),Mloc,rmsgy2(1,2),Mloc,Mloc,1,cudaMemcpyHostToDevice,streamID(4))
    else
       JbegTridy = Jbeg+1
    endif   

!=================================
!X
    call triDx_forwardsweep_cuda_kernel<<<dim3(ceiling(real(Nloc)/tBlock%x), 1, 1),tBlock,0,streamID(2)>>>&
        ( myAx_T,myCx_T,myDx_T,Mloc, Nloc,IbegTridx,Iend,Jbeg,Jend)
!==================================
!Y1
    call triDy_forwardsweep_cuda_kernel<<<dim3(ceiling(real(Mloc)/tBlock%x), 1, 1),tBlock,0,streamID(3)>>>&
        ( myAy_d,myC1y_d,myDy_d,Mloc, Nloc,Ibeg,Iend,JbegTridy,Jend)
!==================================
!Y2
    call triDy_forwardsweep_cuda_kernel<<<dim3(ceiling(real(Mloc)/tBlock%x), 1, 1),tBlock,0,streamID(4)>>>&
         ( myAy_d,myCy_d,myD2y_d,Mloc, Nloc,Ibeg,Iend,JbegTridy,Jend)
!
        !istat = cudaDeviceSynchronize()
!=================================
!X
    if ( n_east .ne. MPI_PROC_NULL ) then
        istat = cudaMemcpy2DAsync(smsgx(1,1),Nloc,myDx_T(1,Iend),Nloc,&
            Nloc,1,cudaMemcpyDeviceToHost,streamID(2))
        istat = cudaMemcpy2DAsync(smsgx(1,2),Nloc,myCx_T(1,Iend),Nloc,&
            Nloc,1,cudaMemcpyDeviceToHost,streamID(2))
        istat = cudaStreamSynchronize(streamID(2))
        call MPI_ISEND( smsgx, 2*Nloc, MPI_SP,&
            n_east, 0, comm2d, req, ier )
        call MPI_WAIT( req, status, ier )
    endif
!==================================
!Y1&Y2
    if ( n_nrth .ne. MPI_PROC_NULL ) then
        istat = cudaMemcpy2dAsync(smsgy(1,1),Mloc,myDy_d(1,Jend),Mloc,&
            Mloc,1,cudaMemcpyDeviceToHost,streamID(3))
        istat = cudaMemcpy2dAsync(smsgy(1,2),Mloc,myC1y_d(1,Jend),Mloc,&
            Mloc,1,cudaMemcpyDeviceToHost,streamID(3))
        istat = cudaMemcpy2dAsync(smsgy2(1,1),Mloc,myD2y_d(1,Jend),Mloc,&
            Mloc,1,cudaMemcpyDeviceToHost,streamID(4))
        istat = cudaMemcpy2dAsync(smsgy2(1,2),Mloc,myCy_d(1,Jend),Mloc,&
            Mloc,1,cudaMemcpyDeviceToHost,streamID(4))
        istat = cudaStreamSynchronize(streamID(3))
        istat = cudaStreamSynchronize(streamID(4))
       call MPI_ISEND( smsgy, 2*Mloc, MPI_SP,&
            n_nrth, 1, comm2d, req2(1), ier )
       call MPI_ISEND( smsgy2, 2*Mloc, MPI_SP,&
            n_nrth, 2, comm2d, req2(2), ier )
       call MPI_WAITALL(2, req2, status2, ier )
    endif

!=================================
!X
! back substitution
    if ( n_east .ne. MPI_PROC_NULL ) then
       call MPI_IRECV( rmsgx, 2*Nloc, MPI_SP,&
            n_east, 3, comm2d, req, ier )
       call MPI_WAIT( req, status, ier )
       IendTridx = Iend
       istat = cudaMemcpy2dAsync(myFx_T(1,Iend+1),Nloc,rmsgx(1,1),Nloc,Nloc,1,cudaMemcpyHostToDevice,streamID(2))
    else
       IendTridx = Iend-1
       istat=cudaMemcpy2dAsync(myFx_T(1,Iend),Nloc,myDx_T(1,Iend),Nloc,Nloc,1,cudaMemcpyHostToDevice,streamID(2))
    endif   
!==================================
!Y1&Y2

    if ( n_nrth .ne. MPI_PROC_NULL ) then
       call MPI_IRECV( rmsgy, 2*Mloc, MPI_SP,&
            n_nrth, 4, comm2d, req2(1), ier )
       call MPI_IRECV( rmsgy2, 2*Mloc, MPI_SP,&
            n_nrth, 5, comm2d, req2(2), ier )
       call MPI_WAITALL(2, req2, status2, ier )
       JendTridy = Jend
       istat = cudaMemcpy2dAsync(Y1_d(1,Jend+1),Mloc,rmsgy(1,1),Mloc,Mloc,1,cudaMemcpyHostToDevice,streamID(3))
       istat = cudaMemcpy2dAsync(Y2_d(1,Jend+1),Mloc,rmsgy2(1,1),Mloc,Mloc,1,cudaMemcpyHostToDevice,streamID(4))
    else
       JendTridy = Jend-1
       istat=cudaMemcpy2dAsync(Y1_d(1,Jend),Mloc,myDy_d(1,Jend),Mloc,Mloc,1,cudaMemcpyHostToDevice,streamID(3))
       istat=cudaMemcpy2dAsync(Y2_d(1,Jend),Mloc,myD2y_d(1,Jend),Mloc,Mloc,1,cudaMemcpyHostToDevice,streamID(4))
    endif   
!
!=================================
!X
    call triDx_backsubstitution_cuda_kernel<<<dim3(ceiling(real(Nloc)/tBlock%x), 1, 1),tBlock,0,streamID(2)>>>&
        ( myCx_T,myDx_T,myFx_T,Mloc, Nloc,Ibeg,IendTridx,Jbeg,Jend)
!==================================
!Y1
    call triDy_backsubstitution_cuda_kernel<<<dim3(ceiling(real(Mloc)/tBlock%x), 1, 1),tBlock,0,streamID(3)>>>&
        ( myC1y_d,myDy_d,Y1_d,Mloc, Nloc,Ibeg,Iend,Jbeg,JendTridy)
!==================================
!Y2
    call triDy_backsubstitution_cuda_kernel<<<dim3(ceiling(real(Mloc)/tBlock%x), 1, 1),tBlock,0,streamID(4)>>>&
        ( myCy_d,myD2y_d,Y2_d,Mloc, Nloc,Ibeg,Iend,Jbeg,JendTridy)
!
    !    istat = cudaDeviceSynchronize()
!=================================
!X
    if ( n_west .ne. MPI_PROC_NULL ) then
        istat = cudaMemcpy2DAsync(smsgx(1,1),Nloc,myFx_T(1,Ibeg),Nloc,&
            Nloc,1,cudaMemcpyDeviceToHost,streamID(2))
        istat = cudaStreamSynchronize(streamID(2))
       call MPI_ISEND( smsgx, 2*Nloc, MPI_SP,&
            n_west, 3, comm2d, req, ier )
       call MPI_WAIT( req, status, ier )
    endif
!==================================
!Y1&Y2
    if ( n_suth .ne. MPI_PROC_NULL ) then
        istat = cudaMemcpy2dAsync(smsgy(1,1),Mloc,Y1_d(1,Jbeg),Mloc,&
            Mloc,1,cudaMemcpyDeviceToHost,streamID(3))
        istat = cudaMemcpy2dAsync(smsgy2(1,1),Mloc,Y2_d(1,Jbeg),Mloc,&
            Mloc,1,cudaMemcpyDeviceToHost,streamID(4))
        istat = cudaStreamSynchronize(streamID(3))
        istat = cudaStreamSynchronize(streamID(4))
       call MPI_ISEND( smsgy, 2*Mloc, MPI_SP,&
            n_suth, 4, comm2d, req2(1), ier )
       call MPI_ISEND( smsgy2, 2*Mloc, MPI_SP,&
            n_suth, 5, comm2d, req2(2), ier )
       call MPI_WAITALL(2, req2, status2, ier )
    endif
!=========================

end subroutine triDxDy_mgpu_cuda_v2
# endif

!---------------------------------------------------------------------------------------!
!                                                                                       !
!   Subroutine triDx_mgpu_cuda is to call tridiagonal solver kernel in x direction           !
!                                                                                       !
!         History : created on 04/04 2019                                               !
!                                                                                       !
!---------------------------------------------------------------------------------------!
# if defined (MGPU)
subroutine triDx_mgpu_cuda
    implicit none
    real(SP), dimension(:,:),allocatable,pinned :: rmsg, smsg
    !real(SP), dimension(Nloc,2),device :: rmsg_d, smsg_d
    integer :: IbegTridx, IendTridx
    INTEGER :: len
    integer status(MPI_STATUS_SIZE), req

    allocate( rmsg(Nloc,2), smsg(Nloc,2) )

! forward sweep

        if ( n_west .ne. MPI_PROC_NULL ) then
           call MPI_IRECV( rmsg, 2*Nloc, MPI_SP,&
                n_west, 0, comm2d, req, ier )
           call MPI_WAIT( req, status, ier )
           IbegTridx = Ibeg
           !istat = cudaMemcpy2d(rmsg_d,Nloc,rmsg,Nloc,Nloc,2)
           istat = cudaMemcpy2dAsync(myDx_T(1,Ibeg-1),Nloc,rmsg(1,1),Nloc,Nloc,1,cudaMemcpyHostToDevice,streamID(2))
           istat = cudaMemcpy2dAsync(myCx_T(1,Ibeg-1),Nloc,rmsg(1,2),Nloc,Nloc,1,cudaMemcpyHostToDevice,streamID(2))
       else
           IbegTridx = Ibeg+1
        endif   

!
     call triDx_forwardsweep_cuda_kernel<<<dim3(ceiling(real(Nloc)/tBlock%x), 1, 1),tBlock,0,streamID(2)>>>&
         ( myAx_T,myCx_T,myDx_T,Mloc, Nloc,IbegTridx,Iend,Jbeg,Jend)
!
        if ( n_east .ne. MPI_PROC_NULL ) then
            istat = cudaMemcpy2DAsync(smsg(1,1),Nloc,myDx_T(1,Iend),Nloc,&
                Nloc,1,cudaMemcpyDeviceToHost,streamID(2))
            istat = cudaMemcpy2DAsync(smsg(1,2),Nloc,myCx_T(1,Iend),Nloc,&
                Nloc,1,cudaMemcpyDeviceToHost,streamID(2))
            istat = cudaStreamSynchronize(streamID(2))
            call MPI_ISEND( smsg, 2*Nloc, MPI_SP,&
                n_east, 0, comm2d, req, ier )
            call MPI_WAIT( req, status, ier )
        endif

! back substitution
        if ( n_east .ne. MPI_PROC_NULL ) then
           call MPI_IRECV( rmsg, 2*Nloc, MPI_SP,&
                n_east, 1, comm2d, req, ier )
           call MPI_WAIT( req, status, ier )
           IendTridx = Iend
           !istat = cudaMemcpy2d(rmsg_d(1,1),Nloc,rmsg(1,1),Nloc,Nloc,2)
           istat = cudaMemcpy2dAsync(myFx_T(1,Iend+1),Nloc,rmsg(1,1),Nloc,Nloc,1,cudaMemcpyHostToDevice,streamID(2))
        else
           IendTridx = Iend-1
           istat=cudaMemcpy2dAsync(myFx_T(1,Iend),Nloc,myDx_T(1,Iend),Nloc,Nloc,1,cudaMemcpyHostToDevice,streamID(2))
        endif   
     call triDx_backsubstitution_cuda_kernel<<<dim3(ceiling(real(Nloc)/tBlock%x), 1, 1),tBlock,0,streamID(2)>>>&
         ( myCx_T,myDx_T,myFx_T,Mloc, Nloc,Ibeg,IendTridx,Jbeg,Jend)
!
        if ( n_west .ne. MPI_PROC_NULL ) then
            istat = cudaMemcpy2DAsync(smsg(1,1),Nloc,myFx_T(1,Ibeg),Nloc,&
                Nloc,1,cudaMemcpyDeviceToHost,streamID(2))
            istat = cudaStreamSynchronize(streamID(2))
           call MPI_ISEND( smsg, 2*Nloc, MPI_SP,&
                n_west, 1, comm2d, req, ier )
           call MPI_WAIT( req, status, ier )
        endif

    deallocate( rmsg, smsg )
end subroutine triDx_mgpu_cuda
# endif

subroutine triDx_triDy_cuda
    implicit none
!
    tBlock_tran = dim3(BlockDimX_2D,BlockDimY_2D,1)
    grid_tran = dim3(ceiling(real(Nloc)/tBlock_tran%x), ceiling(real(Mloc)/tBlock_tran%y), 1)
    grid_tran_back = dim3( ceiling(real(Mloc)/tBlock_tran%x), ceiling(real(Nloc)/tBlock_tran%y), 1 )

    call transpose_kernel<<<grid_tran,tBlock_tran>>>(myAx_d, myAx_T, Nloc, Mloc)
    call transpose_kernel<<<grid_tran,tBlock_tran>>>(myCx_d, myCx_T, Nloc, Mloc)
    call transpose_kernel<<<grid_tran,tBlock_tran>>>(myDx_d, myDx_T, Nloc, Mloc)
# if defined (MGPU)
    call triDxDy_mgpu_cuda_v1
# else
     tBlock = dim3(32,1,1)
     call triDx_cuda_kernel<<<dim3(ceiling(real(Nloc)/tBlock%x),1,1),tBlock, 0,streamID(2)>>>&
         ( myAx_T,myCx_T,myDx_T,myFx_T,Mloc, Nloc,Ibeg,Iend,Jbeg,Jend)
     call triDy_cuda_kernel<<<dim3(ceiling(real(Mloc)/tBlock%x),1,1),tBlock, 0,streamID(3)>>>&
         (myAy_d, myCy_d, myDy_d, myFy_d,&
                                        Mloc, Nloc,&
                                        Ibeg,Iend,&
                                        Jbeg,Jend)
# endif

!Make sure all previous kernels and copys completed
    istat = cudaDeviceSynchronize()
! Be aware of that the myF_d was substituted with U_d
    call transpose_kernel<<<grid_tran_back, tBlock_tran>>>(myFx_T, U_d,Mloc,Nloc)
!
    V_d=myFy_d

end subroutine triDx_triDy_cuda

subroutine periodic_triDx_triDy_cuda
    implicit none
    integer :: i,j
# if defined (MGPU)
    INTEGER :: req
    INTEGER :: l,II
    integer :: status(MPI_STATUS_SIZE)
    REAL(SP),DIMENSION(Mloc) :: xx
# endif
    real(SP),dimension(Mloc),device :: sbeta_d
    real(SP),dimension(Mloc) :: sbeta
            
    if (.NOT. Allocated_2) then
        Allocated_2 = .TRUE.
        allocate( Y1_d(Mloc,Nloc),Y2_d(Mloc,Nloc),&
                  A_beg_d(Mloc),C_end_d(Mloc),&
                  Y1_end_d(Mloc),Y2_end_d(Mloc))
        allocate(A_beg(Mloc),C_end(Mloc),&
                 Y1_end(Mloc),Y2_end(Mloc))
    endif
# if defined (MGPU)
    IF (PY>1)THEN
! south
        DO II = 1,PX
            if(myid==ProcessorID(II,1)) then
              istat = cudaMemcpy2D(A_beg_d(1),Mloc,myAy_d(1,Jbeg),Mloc,Mloc,1)
              istat = cudaMemcpy(A_beg,A_beg_d,Mloc)
!!$cuf kernel do(1) <<<*,*>>>
!                do i = 1,Mloc
!                  A_beg_d(i)=myA_d(i,Jbeg)
!                enddo
! send from master
                call MPI_SEND(A_beg,Mloc,MPI_SP,ProcessorID(II,PY),0,MPI_COMM_WORLD,ier)

            endif ! end myid

            if(myid==ProcessorID(II,PY))then
              call MPI_IRECV(xx,Mloc,MPI_SP,ProcessorID(II,1),0,MPI_COMM_WORLD,req,ier)
              call MPI_WAIT(req,status,ier)
              istat = cudaMemcpy(A_beg_d,xx,Mloc)
            endif
! north
            if(myid==ProcessorID(II,PY)) then
              istat = cudaMemcpy2D(C_end_d(1),Mloc,myCy_d(1,Jend),Mloc,Mloc,1)
              istat = cudaMemcpy(C_end,C_end_d,Mloc)
! send from master
              call MPI_SEND(C_end,Mloc,MPI_SP,ProcessorID(II,1),0,MPI_COMM_WORLD,ier)
            endif ! end myid

            if(myid==ProcessorID(II,1))then
              call MPI_IRECV(xx,Mloc,MPI_SP,ProcessorID(II,PY),0,MPI_COMM_WORLD,req,ier)
              call MPI_WAIT(req,status,ier)
              istat = cudaMemcpy(C_end_d,xx,Mloc)
            endif
 
        ENDDO  ! end PX

    ELSE  ! PY = 1
        istat = cudaMemcpy2D(A_beg_d(1),Mloc,myAy_d(1,Jbeg),Mloc,Mloc,1)
        istat = cudaMemcpy2D(C_end_d(1),Mloc,myCy_d(1,Jend),Mloc,Mloc,1)
    ENDIF
# else
    istat = cudaMemcpy2D(A_beg_d(1),Mloc,myAy_d(1,Jbeg),Mloc,Mloc,1)
    istat = cudaMemcpy2D(C_end_d(1),Mloc,myCy_d(1,Jend),Mloc,Mloc,1)
# endif
!  initialize d
    myD2y_d = 0.0_SP
# if defined (MGPU)
    DO II=1,PX
      if(myid == ProcessorID(II,1)) then
!$cuf kernel do(1) <<<*,*>>>
        do i= Ibeg, Iend
            myCy_d(i,Jbeg)=myCy_d(i,Jbeg)/(1.0_SP+C_end_d(i)) 
            myDy_d(i,Jbeg)=myDy_d(i,Jbeg)/(1.0_SP+C_end_d(i))
            myD2y_d(i,Jbeg)=1.0_SP/(1.0_SP+C_end_d(i))
        enddo
      endif
      if(myid == ProcessorID(II,PY)) then
!$cuf kernel do(1) <<<*,*>>>
        do i= Ibeg, Iend
            myAy_d(i,Jend)=myAy_d(i,Jend)/(1.0_SP+A_beg_d(i))
            myDy_d(i,Jend)=myDy_d(i,Jend)/(1.0_SP+A_beg_d(i))
            myD2y_d(i,Jend)=-1.0_SP/(1.0_SP+A_beg_d(i))
        enddo
      endif
   ENDDO  ! end PX
# else
!$cuf kernel do(1) <<<*,*>>>
    do i= Ibeg, Iend
        myCy_d(i,Jbeg)=myCy_d(i,Jbeg)/(1.0_SP+C_end_d(i)) 
        myDy_d(i,Jbeg)=myDy_d(i,Jbeg)/(1.0_SP+C_end_d(i))
        myAy_d(i,Jend)=myAy_d(i,Jend)/(1.0_SP+A_beg_d(i))
        myDy_d(i,Jend)=myDy_d(i,Jend)/(1.0_SP+A_beg_d(i))
        myD2y_d(i,Jbeg)=1.0_SP/(1.0_SP+C_end_d(i))
        myD2y_d(i,Jend)=-1.0_SP/(1.0_SP+A_beg_d(i))
    enddo
# endif
!
    myC1y_d=myCy_d
!
    tBlock_tran = dim3(BlockDimX_2D,BlockDimY_2D,1)
    grid_tran = dim3(ceiling(real(Nloc)/tBlock_tran%x), ceiling(real(Mloc)/tBlock_tran%y), 1)
    grid_tran_back = dim3( ceiling(real(Mloc)/tBlock_tran%x), ceiling(real(Nloc)/tBlock_tran%y), 1 )

    call transpose_kernel<<<grid_tran,tBlock_tran>>>(myAx_d, myAx_T, Nloc, Mloc)
    call transpose_kernel<<<grid_tran,tBlock_tran>>>(myCx_d, myCx_T, Nloc, Mloc)
    call transpose_kernel<<<grid_tran,tBlock_tran>>>(myDx_d, myDx_T, Nloc, Mloc)
# if defined (MGPU)
    call triDxDy_mgpu_cuda_v2
# else
     tBlock = dim3(32,1,1)
     call triDx_cuda_kernel<<<dim3(ceiling(real(Nloc)/tBlock%x),1,1),tBlock, 0,streamID(2)>>>&
         ( myAx_T,myCx_T,myDx_T,myFx_T,Mloc, Nloc,Ibeg,Iend,Jbeg,Jend)
     call triDy_cuda_kernel<<<dim3(ceiling(real(Mloc)/tBlock%x),1,1),tBlock, 0,streamID(3)>>>&
         (myAy_d, myC1y_d, myDy_d, Y1_d,&
                                        Mloc, Nloc,&
                                        Ibeg,Iend,&
                                        Jbeg,Jend)
    call triDy_cuda_kernel<<<dim3(ceiling(real(Mloc)/tBlock%x), 1, 1),tBlock, 0,streamID(4)>>>&
        (myAy_d, myCy_d, myD2y_d, Y2_d,&
                                        Mloc, Nloc,&
                                        Ibeg,Iend,&
                                        Jbeg,Jend)
# endif

!Make sure all previous kernels and copys completed
    istat = cudaDeviceSynchronize()
! Be aware of that the myF_d was substituted with U_d
     call transpose_kernel<<<grid_tran_back, tBlock_tran>>>(myFx_T, U_d,Mloc,Nloc)

!  Compute scale factor
# if defined (MGPU)
   IF(PY>1)THEN
! transfer y2 y1 end
     DO II=1,PX
! y2
       if(myid==ProcessorID(II,PY)) then
         istat = cudaMemcpy2D(Y2_end_d(1),Mloc,Y2_d(1,Jend),Mloc,&
             Mloc,1,cudaMemcpyDeviceToDevice)
         istat = cudaMemcpy(Y2_end,Y2_end_d,Mloc)
   ! send from master
         call MPI_SEND(Y2_end,Mloc,MPI_SP,ProcessorID(II,1),0,MPI_COMM_WORLD,ier)
!         call MPI_SEND(Y2_d(1,Jend),1,YType,ProcessorID(II,1),0,MPI_COMM_WORLD,ier)
       endif ! end myid
   
       if(myid==ProcessorID(II,1))then
         call MPI_IRECV(xx,Mloc,MPI_SP,ProcessorID(II,PY),0,MPI_COMM_WORLD,req,ier)
         call MPI_WAIT(req,status,ier)
         istat = cudaMemcpy(Y2_end_d,xx,Mloc)
       endif
   ! y1
       if(myid==ProcessorID(II,PY)) then
         istat = cudaMemcpy2D(Y1_end_d(1),Mloc,Y1_d(1,Jend),Mloc,&
             Mloc,1,cudaMemcpyDeviceToDevice)
         istat = cudaMemcpy(Y1_end,Y1_end_d,Mloc)
   ! send from master
!           call MPI_SEND(Y1_d(1,Jend),1,YType,ProcessorID(II,1),0,MPI_COMM_WORLD,ier)
           call MPI_SEND(Y1_end,Mloc,MPI_SP,ProcessorID(II,1),0,MPI_COMM_WORLD,ier)
       endif ! end myid
   
       if(myid==ProcessorID(II,1))then
         !call MPI_IRECV(Y1_end_d(1),1,YType,ProcessorID(II,PY),0,MPI_COMM_WORLD,req(1),ier)
         call MPI_IRECV(xx,Mloc,MPI_SP,ProcessorID(II,PY),0,MPI_COMM_WORLD,req,ier)
         call MPI_WAIT(req,status,ier)
         istat = cudaMemcpy(Y1_end_d,xx,Mloc)
       endif
   
   ! calculate beta
       if(myid==ProcessorID(II,1))then
!$cuf kernel do(1) <<<*,*>>>
          DO I=Ibeg,Iend
            sbeta_d(i) = ( C_end_d(i)*Y1_d(i,Jbeg) - A_beg_d(i)*Y1_end_d(i)) / &
                ( 1.0_SP - ( C_end_d(i)*Y2_d(i,Jbeg) - A_beg_d(i)*Y2_end_d(i) ) )
         ENDDO
       endif ! end myid
   ! transfer to other II
       if(myid==ProcessorID(II,1))then
           istat = cudaMemcpy(sbeta,sbeta_d, Mloc)
           DO l = 2,PY
             call MPI_SEND(sbeta,Mloc,MPI_SP,ProcessorID(II,l),0,MPI_COMM_WORLD,ier)
           ENDDO
       endif  
       DO l = 2,PY
         if(myid==ProcessorID(II,l))then
             call MPI_IRECV(xx,Mloc,MPI_SP,ProcessorID(II,1),0,MPI_COMM_WORLD,req,ier)
             call MPI_WAIT(req,status,ier)
!             beta_d = xx
            istat = cudaMemcpy(sbeta_d, xx, Mloc)
         endif
       ENDDO      ! end l

    ENDDO ! end PX
     
  ELSE  ! PY=1
!$cuf kernel do(1) <<<*,*>>>
    do i=Ibeg, Iend
        sbeta_d(i) = ( C_end_d(i)*Y1_d(i,Jbeg) - A_beg_d(i)*Y1_d(i,Jend)) / &
            ( 1.0_SP - ( C_end_d(i)*Y2_d(i,Jbeg) - A_beg_d(i)*Y2_d(i,Jend) ) )
    enddo

  ENDIF ! if PY>1 
# else
!$cuf kernel do(1) <<<*,*>>>
    do i=Ibeg, Iend
        sbeta_d(i) = ( C_end_d(i)*Y1_d(i,Jbeg) - A_beg_d(i)*Y1_d(i,Jend)) / &
            ( 1.0_SP - ( C_end_d(i)*Y2_d(i,Jbeg) - A_beg_d(i)*Y2_d(i,Jend) ) )
    enddo
# endif
!$cuf kernel do(2) <<<*,*>>>
    do j = Jbeg, Jend
        do i=Ibeg, Iend        
            myFy_d(i,j)=Y1_d(i,j)+sbeta_d(i)*Y2_d(i,j)
            V_d(i,j)=myFy_d(i,j)
        enddo
    enddo
end subroutine periodic_triDx_triDy_cuda

# if defined (MGPU)
subroutine triDy_mgpu_cuda(myA_D,myC_D,myD_D,myF_D,stream_dy)
    implicit none
    integer(kind=cuda_stream_kind) :: stream_dy
    real(SP), dimension(Mloc,Nloc),device :: myA_d, myC_d, myD_d
    real(SP), dimension(Mloc,Nloc),device :: myF_d
    real(SP), dimension(:,:),allocatable,pinned :: rmsg, smsg
    !real(SP), dimension(Mloc,2),device :: rmsg_d, smsg_d
    integer JbegTridy,JendTridy
    INTEGER :: req
    INTEGER :: len
    integer :: status(MPI_STATUS_SIZE)

    allocate( rmsg(Mloc,2),smsg(Mloc,2))
        
! forward sweep

     tBlock = dim3(32,1,1)
        if ( n_suth .ne. MPI_PROC_NULL ) then
           call MPI_IRECV( rmsg, 2*Mloc, MPI_SP,&
                n_suth, 0, comm2d, req, ier )
           call MPI_WAIT( req, status, ier )
           JbegTridy = Jbeg
           !istat = cudaMemcpy2d(rmsg_d(1,1),Mloc,rmsg(1,1),Mloc,Mloc,2)
           istat =  cudaMemcpy2dAsync(myD_d(1,Jbeg-1),Mloc,rmsg(1,1),Mloc,Mloc,1,cudaMemcpyHostToDevice,stream_dy)
           istat = cudaMemcpy2dAsync(myC_d(1,Jbeg-1),Mloc,rmsg(1,2),Mloc,Mloc,1,cudaMemcpyHostToDevice,stream_dy)
       else
           JbegTridy = Jbeg+1
        endif   
!
     call triDy_forwardsweep_cuda_kernel<<<dim3(ceiling(real(Mloc)/tBlock%x), 1, 1),tBlock,0,stream_dy>>>&
         ( myA_d,myC_d,myD_d,Mloc, Nloc,Ibeg,Iend,JbegTridy,Jend)

        if ( n_nrth .ne. MPI_PROC_NULL ) then
            istat = cudaMemcpy2dAsync(smsg(1,1),Mloc,myD_d(1,Jend),Mloc,&
                Mloc,1,cudaMemcpyDeviceToHost,stream_dy)
            istat = cudaMemcpy2dAsync(smsg(1,2),Mloc,myC_d(1,Jend),Mloc,&
                Mloc,1,cudaMemcpyDeviceToHost,stream_dy)
            istat = cudaStreamSynchronize(stream_dy)
           call MPI_ISEND( smsg, 2*Mloc, MPI_SP,&
                n_nrth, 0, comm2d, req, ier )
           call MPI_WAIT( req, status, ier )
        endif

! back substitution

        if ( n_nrth .ne. MPI_PROC_NULL ) then
           call MPI_IRECV( rmsg, 2*Mloc, MPI_SP,&
                n_nrth, 1, comm2d, req, ier )
           call MPI_WAIT( req, status, ier )
           JendTridy = Jend
           !istat = cudaMemcpy2d(rmsg_d(1,1),Mloc,rmsg(1,1),Mloc,Mloc,2)
           istat = cudaMemcpy2dAsync(myF_d(1,Jend+1),Mloc,rmsg(1,1),Mloc,Mloc,1,cudaMemcpyHostToDevice,stream_dy)
        else
           JendTridy = Jend-1
           istat=cudaMemcpy2dAsync(myF_d(1,Jend),Mloc,myD_d(1,Jend),Mloc,Mloc,1,cudaMemcpyHostToDevice,stream_dy)
        endif   

     call triDy_backsubstitution_cuda_kernel<<<dim3(ceiling(real(Mloc)/tBlock%x), 1, 1),tBlock,0,stream_dy>>>&
         ( myC_d,myD_d,myF_d,Mloc, Nloc,Ibeg,Iend,Jbeg,JendTridy)

        if ( n_suth .ne. MPI_PROC_NULL ) then
            istat = cudaMemcpy2dAsync(smsg(1,1),Mloc,myF_d(1,Jbeg),Mloc,&
                Mloc,1,cudaMemcpyDeviceToHost,stream_dy)
            istat = cudaStreamSynchronize(stream_dy)
           call MPI_ISEND( smsg, 2*Mloc, MPI_SP,&
                n_suth, 1, comm2d, req, ier )
           call MPI_WAIT( req, status, ier )
        endif

    deallocate( rmsg,smsg)
end subroutine triDy_mgpu_cuda
# endif


!---------------------------------------------------------------------------------------!
!                                                                                       !
!   Subroutine triDx_cusparse is to use CuSparse toolbox to batch-solve tridiagnal      !
!                                                                                       !
!         History : created on 04/04 2019                                               !
!                                                                                       !
!---------------------------------------------------------------------------------------!
!!version 1 
subroutine triDx_cusparse
!    use iso_c_binding
    implicit none

# if defined (TIMING)
        call CPU_TIME(t1)
# endif

!        integer(8) :: bufsize
!        character(kind=c_char),device,allocatable,target :: buf(:)
!        !istat = cusparseSgtsv2StridedBatch_bufferSizeExt(cusparseh, Mloc-2*NGhost, Arow, Brow, Crow, Drow, Nloc-2*NGhost, Mloc-2*NGhost,bufsize)
!        allocate(buf(1000000))
!        istat = cusparseSgtsv2StridedBatch(cusparseh, Mloc-2*NGhost, Arow, Brow, Crow, Drow, Nloc-2*NGhost, Mloc-2*NGhost,c_loc(buf))
        istat = cusparseSgtsvStridedBatch(cusparseh, Mloc-2*NGhost, Arow, Brow, Crow, Drow, Nloc-2*NGhost, Mloc-2*NGhost)
        if(istat /= CUSPARSE_STATUS_SUCCESS) then
                write(*,*) 'Cusparse tridiagonal solver failed', istat
                stop
        endif
        
!$cuf kernel do(2) <<<*,*>>>
        do j = Jbeg,Jend
        do i = Ibeg,Iend
            U_d(i,j) = Drow((i-NGhost)+(j-NGhost-1)*(Mloc-2*NGhost))
        enddo
        enddo
# if defined (TIMING)
        istat = cudaDeviceSynchronize()
        call CPU_TIME(t2)
        KernelTime = t2-t1
        TotalKernelTime = TotalKernelTime+KernelTime
# endif
end subroutine triDx_cusparse

subroutine triDy_cusparse
    implicit none

# if defined (TIMING)
        call CPU_TIME(t1)
# endif
        istat = cusparseSgtsvStridedBatch(cusparseh, Nloc-2*NGhost, Arow, Brow, Crow, Drow, Mloc-2*NGhost, Nloc-2*NGhost)
        if(istat /= CUSPARSE_STATUS_SUCCESS) then
                write(*,*) 'Cusparse tridiagonal solver failed', istat
                stop
        endif
        
!$cuf kernel do(2) <<<*,*>>>
        do i = Ibeg,Iend
        do j = Jbeg,Jend
            V_d(i,j) = Drow((j-NGhost)+(i-NGhost-1)*(Jend-Jbeg+1))
        enddo
        enddo
# if defined (TIMING)
        istat = cudaDeviceSynchronize()
        call CPU_TIME(t2)
        KernelTime = t2-t1
        TotalKernelTime = TotalKernelTime+KernelTime
# endif
end subroutine triDy_cusparse


subroutine periodic_triDy_cusparse
    implicit none
    !real(SP), dimension((Mloc-2*NGhost)*(Nloc-2*NGhost)), device :: Y1row, Y2row
    real(SP), dimension(Mloc), device :: sbeta_d
    integer :: i,j
# if defined (MGPU)
    INTEGER:: req
    INTEGER :: len,l,II
    integer :: status(MPI_STATUS_SIZE)
    REAL(SP),DIMENSION(Mloc) :: xx
# endif
            
# if defined (MGPU)
    IF (PY>1)THEN
        len=Mloc
! south
        DO II = 1,PX
            if(myid==ProcessorID(II,1)) then
!$cuf kernel do(1) <<<*,*>>>
              DO i = Ibeg,Iend
                  A_beg_d(I) = Arow(1+(i-NGhost-1)*(Nloc-2*NGhost))
              ENDDO
! send from master
              call MPI_SEND(A_beg_d,len,MPI_SP,ProcessorID(II,PY),0,MPI_COMM_WORLD,ier)

            endif ! end myid

            if(myid==ProcessorID(II,PY))then
              call MPI_IRECV(xx,len,MPI_SP,ProcessorID(II,1),0,MPI_COMM_WORLD,req,ier)
              call MPI_WAIT(req,status,ier)
              A_beg_d=xx
            endif
! north
            if(myid==ProcessorID(II,PY)) then
!$cuf kernel do(1) <<<*,*>>>
              DO i = Ibeg,Iend
                  C_end_d(I) = Crow(Nloc-2*NGhost+(i-NGhost-1)*(Nloc-2*NGhost))
              ENDDO
! send from master
              call MPI_SEND(C_end_d,len,MPI_SP,ProcessorID(II,1),0,MPI_COMM_WORLD,ier)
            endif ! end myid

            if(myid==ProcessorID(II,1))then
              call MPI_IRECV(xx,len,MPI_SP,ProcessorID(II,PY),0,MPI_COMM_WORLD,req,ier)
              call MPI_WAIT(req,status,ier)
              C_end_d = xx
            endif
 
        ENDDO  ! end PX
    ELSE  ! PY = 1
!$cuf kernel do(1) <<<*,*>>>
        DO i = Ibeg,Iend
            A_beg_d(I) = Arow(1+(i-NGhost-1)*(Nloc-2*NGhost))
            C_end_d(I) = Crow(Nloc-2*NGhost+(i-NGhost-1)*(Nloc-2*NGhost))
        ENDDO
    ENDIF
# else
!$cuf kernel do(1) <<<*,*>>>
    DO i = Ibeg,Iend
        A_beg_d(I) = Arow(1+(i-NGhost-1)*NGlob)
        C_end_d(I) = Crow(NGlob+(i-NGhost-1)*NGlob)
    ENDDO
# endif


# if defined (MGPU)
    DO II=1,PX
      if(myid == ProcessorID(II,1)) then
!$cuf kernel do(1) <<<*,*>>>
        do i= Ibeg, Iend
            Crow(1+(i-NGhost-1)*(Nloc-2*NGhost))=Crow(1+(i-NGhost-1)*(Nloc-2*NGhost))/(1+C_end_d(i)) 
            Drow(1+(i-NGhost-1)*(Nloc-2*NGhost))=Drow(1+(i-NGhost-1)*(Nloc-2*NGhost))/(1+C_end_d(i))
        enddo
      endif
      if(myid == ProcessorID(II,PY)) then
!$cuf kernel do(1) <<<*,*>>>
        do i= Ibeg, Iend
            Arow(Nloc-2*NGhost+(i-NGhost-1)*(Nloc-2*NGhost))=Arow(Nloc-2*NGhost+(i-NGhost-1)*(Nloc-2*NGhost))/(1+A_beg_d(i))
            Drow(Nloc-2*NGhost+(i-NGhost-1)*(Nloc-2*NGhost))=Drow(Nloc-2*NGhost+(i-NGhost-1)*(Nloc-2*NGhost))/(1+A_beg_d(i))
        enddo
      endif
   ENDDO  ! end PX
# else
!$cuf kernel do(1) <<<*,*>>>
    do i= Ibeg, Iend
        Crow(1+(i-NGhost-1)*(Nloc-2*NGhost))=Crow(1+(i-NGhost-1)*(Nloc-2*NGhost))/(1+C_end_d(i)) 
        Drow(1+(i-NGhost-1)*(Nloc-2*NGhost))=Drow(1+(i-NGhost-1)*(Nloc-2*NGhost))/(1+C_end_d(i))
        Arow(Nloc-2*NGhost+(i-NGhost-1)*(Nloc-2*NGhost))=Arow(Nloc-2*NGhost+(i-NGhost-1)*(Nloc-2*NGhost))/(1+A_beg_d(i))
        Drow(Nloc-2*NGhost+(i-NGhost-1)*(Nloc-2*NGhost))=Drow(Nloc-2*NGhost+(i-NGhost-1)*(Nloc-2*NGhost))/(1+A_beg_d(i))
    enddo
# endif
!!C1row is not necessary since it is not rewritten,but D1row is necessary
    D1row=Drow   
    Arow(1+(i-NGhost-1)*(Nloc-2*NGhost))=0.0_SP
    Crow(Nloc-2*NGhost+(i-NGhost-1)*(Nloc-2*NGhost))=0.0_SP
    istat = cusparseSgtsvStridedBatch(cusparseh, Nloc-2*NGhost, Arow, Brow, Crow, D1row, Mloc-2*NGhost, Nloc-2*NGhost)
    if(istat /= CUSPARSE_STATUS_SUCCESS) then
        write(*,*) 'Cusparse tridiagonal solver failed', istat
            stop
    endif
    !Y1row = D1row

# if defined (MGPU)
!  initialize d
    Drow = 0.0_SP
    DO II=1,PX
        if(myid==ProcessorID(II,1))then
!$cuf kernel do(1) <<<*,*>>>
          do i=Ibeg, Iend
              Drow(1+(i-NGhost-1)*(Nloc-2*NGhost))=1.0_SP/(1.0_SP+C_end_d(i))
          enddo
        endif
        if(myid==ProcessorID(II,PY))then
!$cuf kernel do(1) <<<*,*>>>
          do i=Ibeg, Iend
              Drow((Nloc-2*NGhost)+(i-NGhost-1)*(Nloc-2*NGhost))=-1.0_SP/(1.0_SP+A_beg_d(i))
          enddo
        endif
    ENDDO ! end PX
# else 
! if serial, NGlob = Nloc-2*NGhost
!$cuf kernel do(1) <<<*,*>>>
    do i=Ibeg, Iend
        Drow(1+(i-NGhost-1)*NGlob)=1.0_SP/(1.0_SP+C_end_d(i))
        Drow(NGlob+(i-NGhost-1)*NGlob)=-1.0_SP/(1.0_SP+A_beg_d(i))
    enddo
!$cuf kernel do(2) <<<*,*>>>
    do i=Ibeg, Iend
        do j = Jbeg+1, Jend-1
            Drow((j-NGhost)+(i-NGhost-1)*NGlob)=ZERO
        enddo
    enddo
# endif

!  Solution	pass two:  Solve By2=w
    istat = cusparseSgtsvStridedBatch(cusparseh,Nloc-2*NGhost, Arow, Brow, Crow, Drow, Mloc-2*NGhost, Nloc-2*NGhost)
    if(istat /= CUSPARSE_STATUS_SUCCESS) then
        write(*,*) 'Cusparse tridiagonal solver failed', istat
            stop
    endif
    !Y2row = Drow

!  Compute scale factor
# if defined (MGPU)
   IF(PY>1)THEN
! transfer y2 y1 end
     DO II=1,PX
! y2
       if(myid==ProcessorID(II,PY)) then
!$cuf kernel do(1) <<<*,*>>>
         do i = Ibeg,Iend
             Y2_end_d(i)=Drow((Nloc-2*NGhost)+(i-NGhost-1)*(Nloc-2*NGhost))
         enddo
   ! send from master
         call MPI_SEND(Y2_end_d,len,MPI_SP,ProcessorID(II,1),0,MPI_COMM_WORLD,ier)
       endif ! end myid
   
       if(myid==ProcessorID(II,1))then
         call MPI_IRECV(xx,len,MPI_SP,ProcessorID(II,PY),0,MPI_COMM_WORLD,req,ier)
         call MPI_WAIT(req,status,ier)
         Y2_end_d = xx
       endif
   ! y1
       if(myid==ProcessorID(II,PY)) then
!$cuf kernel do(1) <<<*,*>>>
         do i = Ibeg,Iend
             Y1_end_d(i)=D1row((Nloc-2*NGhost)+(i-NGhost-1)*(Nloc-2*NGhost))
         enddo
   ! send from master
           call MPI_SEND(Y1_end_d,len,MPI_SP,ProcessorID(II,1),0,MPI_COMM_WORLD,ier)
       endif ! end myid
   
       if(myid==ProcessorID(II,1))then
         call MPI_IRECV(xx,len,MPI_SP,ProcessorID(II,PY),0,MPI_COMM_WORLD,req,ier)
         call MPI_WAIT(req,status,ier)
         Y1_end_d = xx
       endif

   ! calculate beta
       if(myid==ProcessorID(II,1))then
!$cuf kernel do(1) <<<*,*>>>
          DO I=Ibeg,Iend
            sbeta_d(i) = ( C_end_d(i)*D1row(1+(i-NGhost-1)*(Nloc-2*NGhost)) - A_beg_d(i)*Y1_end_d(i) ) / &
                ( 1.0_SP - ( C_end_d(i)*Drow(1+(i-NGhost-1)*(Nloc-2*NGhost)) - A_beg_d(i)*Y2_end_d(i) ) )
         ENDDO
       endif ! end myid
   ! transfer to other II
       if(myid==ProcessorID(II,1))then
           DO l = 2,PY
             call MPI_SEND(sbeta_d,len,MPI_SP,ProcessorID(II,l),0,MPI_COMM_WORLD,ier)
           ENDDO
       endif  
       DO l = 2,PY
         if(myid==ProcessorID(II,l))then
             call MPI_IRECV(xx,len,MPI_SP,ProcessorID(II,1),0,MPI_COMM_WORLD,req,ier)
             call MPI_WAIT(req,status,ier)
             sbeta_d = xx
         endif
       ENDDO      ! end l

    ENDDO ! end PX
     
  ELSE  ! PY=1
!$cuf kernel do(1) <<<*,*>>>
    do i=Ibeg, Iend
        sbeta_d(i) = ( C_end_d(i)*D1row(1+(i-NGhost-1)*(Nloc-2*NGhost)) - A_beg_d(i)*D1row((Nloc-2*NGhost)+(i-NGhost-1)*(Nloc-2*NGhost)) ) / &
            ( 1.0_SP - ( C_end_d(i)*Drow(1+(i-NGhost-1)*(Nloc-2*NGhost)) - A_beg_d(i)*Drow((Nloc-2*NGhost)+(i-NGhost-1)*(Nloc-2*NGhost)) ) )
    enddo

  ENDIF ! if PY>1 

# else

!$cuf kernel do(1) <<<*,*>>>
    do i=Ibeg, Iend
        sbeta_d(i) = ( C_end_d(i)*D1row(1+(i-NGhost-1)*NGlob) - A_beg_d(i)*D1row(NGlob+(i-NGhost-1)*NGlob) ) / &
            ( 1.0_SP - ( C_end_d(i)*Drow(1+(i-NGhost-1)*NGlob) - A_beg_d(i)*Drow(NGlob+(i-NGhost-1)*NGlob) ) )
    enddo
# endif
!$cuf kernel do(2) <<<*,*>>>
    do i=Ibeg, Iend        
        do j = Jbeg, Jend
            V_d(i,j)=D1row((j-NGhost)+(i-NGhost-1)*(Nloc-2*NGhost))+sbeta_d(i)*Drow((j-NGhost)+(i-NGhost-1)*(Nloc-2*NGhost))
        enddo
    enddo
end subroutine periodic_triDy_cusparse

!---------------------------------------------------------------------------------------!
!                                                                                       !
!   Cuda Kernel HUV_BAR_KERNEL is to compute Eta_bar, Ubar,Vbar matrix               !
!                                                                                       !
!         History : created on 04/04 2019                                               !
!                                                                                       !
!---------------------------------------------------------------------------------------!

ATTRIBUTES(GLOBAL) SUBROUTINE HUV_BAR_KERNEL&
            (Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,iista,jjsta,ISTEP,WaveMakerCode,DT,&
            Xc_WK, Yc_WK,Width_WK,Ywidth_WK,Gamma1,Gamma3 &
# if !defined (ITERATION)
            ,ETA_LIMITER, TroughLimit, CrestLimit &
# endif
            )
     IMPLICIT NONE
! Dummy variables     
     INTEGER,VALUE,INTENT(IN)::ISTEP,WaveMakerCode,Ibeg,Iend,Jbeg,Jend,Mloc,Nloc,iista,jjsta
# if !defined (ITERATION)     
     REAL(SP),VALUE,INTENT(IN) :: TroughLimit, CrestLimit
     LOGICAL,VALUE,INTENT(IN) :: ETA_LIMITER
# endif
     REAL(SP),VALUE,INTENT(in) ::Gamma1,Gamma3, Xc_WK, Yc_WK, Width_WK,Ywidth_WK,DT
! set local variables
     REAL(SP) :: F_left,F_right,F_bottom,F_top
     REAL(SP) :: R1,R2,R3
     REAL(SP) :: xmk,ymk
     REAL(SP) :: DXg,DYg
    integer :: I,J,tx,ty
! set local indexes
     tx = threadIdx%x
     ty = threadIdx%y
! Global thread and  block ID
     i = tx + (blockIdx%x-1)*BlockDimX_2D
     j = ty + (blockIdx%y-1)*BlockDimY_2D
! MUSCL-Hancock, Zhou et al., p. 7

! solve eta
     if (i>=Ibeg .AND. i<=Iend .AND. j>=Jbeg .AND. j<=Jend )then
# if defined (CARTESIAN)
         DXg=DX_d
         DYg=DY_d
# else
! only for wavemaker
         DXg=DX_d(1,1)
         DYg=DY_d(1,1)
# endif
         F_left=P_d(i,j)
         F_right=P_d(i+1,J)
         F_bottom=Q_d(i,j)
         F_top=Q_d(i,j+1)
! now work for spherical # if defined (CARTESIAN)
         IF(WaveMakerCode==5.OR.WaveMakerCode==8  &
            .OR.WaveMakerCode==7.OR.WaveMakerCode==6)THEN
# if defined (MGPU)
	     xmk=(I-Ibeg)*DXg + (iista-1)*DXg
	     ymk=(J-Jbeg)*DYg + (jjsta-1)*DYg
# else
            xmk=(I-Ibeg)*DXg
            ymk=(J-Jbeg)*DYg
# endif
         IF(ABS(xmk-Xc_WK)<Width_WK.AND. &
            ABS(ymk-Yc_WK)<Ywidth_WK/2.0_SP)THEN

          R1=-1.0_SP/DXg*(F_right*n_right+F_left*n_left) &
                -1.0_SP/DYg*(F_top*n_top+F_bottom*n_bottom) &
        ! wavemaker
                 +WaveMaker_Mass_d(I,J)
!                +WK_Source      
         ELSE
         R1=-1.0_SP/DXg*(F_right*n_right+F_left*n_left) &
                   -1.0_SP/DYg*(F_top*n_top+F_bottom*n_bottom)
         ENDIF
       ELSEIF(WaveMakerCode==4)THEN
# if defined (MGPU)
	      xmk=(I-Ibeg)*DXg + (iista-1)*DXg
	      ymk=(J-Jbeg)*DYg + (jjsta-1)*DYg
# else
            xmk=(I-Ibeg)*DXg
            ymk=(J-Jbeg)*DYg
# endif
         IF(ABS(xmk-Xc_WK)<Width_WK.AND. &
            ABS(ymk-Yc_WK)<Ywidth_WK/2.0_SP)THEN
          
          R1=-1.0_SP/DXg*(F_right*n_right+F_left*n_left) &
                -1.0_SP/DYg*(F_top*n_top+F_bottom*n_bottom) &
        ! wavemaker 
                 +WaveMaker_Mass_d(I,J)    
         ELSE
         R1=-1.0_SP/DXg*(F_right*n_right+F_left*n_left) &
                   -1.0_SP/DYg*(F_top*n_top+F_bottom*n_bottom)
         ENDIF
       ELSEIF(WaveMakerCode==9)THEN
# if defined (MGPU)
	      xmk=(I-Ibeg)*DXg + (iista-1)*DXg
	      ymk=(J-Jbeg)*DYg + (jjsta-1)*DYg
# else
            xmk=(I-Ibeg)*DXg
            ymk=(J-Jbeg)*DYg
# endif
         IF(ABS(xmk-Xc_WK)<Width_WK.AND. &
            ABS(ymk-Yc_WK)<Ywidth_WK/2.0_SP)THEN
          
          R1=-1.0_SP/DXg*(F_right*n_right+F_left*n_left) &
                -1.0_SP/DYg*(F_top*n_top+F_bottom*n_bottom) &
                +WaveMaker_Mass_d(I,J)      
         ELSE
         R1=-1.0_SP/DXg*(F_right*n_right+F_left*n_left) &
                   -1.0_SP/DYg*(F_top*n_top+F_bottom*n_bottom)
         ENDIF   

       ELSEIF(WaveMakerCode==10)THEN
# if defined (MGPU)
	      xmk=(I-Ibeg)*DXg + (iista-1)*DXg
	      ymk=(J-Jbeg)*DYg + (jjsta-1)*DYg
# else
            xmk=(I-Ibeg)*DXg
            ymk=(J-Jbeg)*DYg
# endif
         IF(ABS(xmk-Xc_WK)<Width_WK.AND. &
            ABS(ymk-Yc_WK)<Ywidth_WK/2.0_SP)THEN
          
          R1=-1.0_SP/DXg*(F_right*n_right+F_left*n_left) &
                -1.0_SP/DYg*(F_top*n_top+F_bottom*n_bottom) &
                +WaveMaker_Mass_d(I,J)      
         ELSE
         R1=-1.0_SP/DXg*(F_right*n_right+F_left*n_left) &
                   -1.0_SP/DYg*(F_top*n_top+F_bottom*n_bottom)
         ENDIF 
   
      ELSE ! no wk_wavemaker, theres bug in version 1.1 Dxg,Dyg should be 
           ! replaced by Dxg() and Dy()
# if defined (CARTESIAN)
        R1=-1.0_SP/DXg*(F_right*n_right+F_left*n_left) &
                   -1.0_SP/DYg*(F_top*n_top+F_bottom*n_bottom)
# if defined (VESSEL_PANEL_SOURCE)
        R1=R1 +VesselFluxGradient_d(I,J)
# endif

# else
        R1=-1.0_SP/DX_d(I,J)*(F_right*n_right+F_left*n_left) &
                   -1.0_SP/DY_d(I,J)*(F_top*n_top+F_bottom*n_bottom)
# endif
      ENDIF
# if defined (CARTESIAN)
! do nothing
# else

# if defined (ZALPHA)
! extra terms for z_alpha: 1/R tan theta (HV+hV4)
        R1=R1+1.0_SP/R_earth*TAN(Lat_theta_d(I,J))*(HV_d(I,J)+ &
                Depth_d(I,J)*V4_d(I,J)*Gamma1)
# else
! extra terms for depth-averaged u: 1/R tan theta HV
        R1=R1+1.0_SP/R_earth*TAN(Lat_theta_d(I,J))*HV_d(I,J)
# endif

# endif

# if defined (ITERATION)
      Eta_d(I,J)=ALPHA_d(ISTEP)*Eta0_d(I,J)+BETA_d(ISTEP)*(EtaOld_d(I,J)+DT*R1)
# else
      Eta_d(I,J)=ALPHA_d(ISTEP)*Eta0_d(I,J)+BETA_d(ISTEP)*(Eta_d(I,J)+DT*R1)

! eta_limiter is used for the case as wave touches the seabed within wavemaker
      IF(ETA_LIMITER)THEN
        IF(Eta_d(I,J)<TroughLimit)Eta_d(I,J)=TroughLimit
        IF(Eta_d(I,J)>CrestLimit)Eta_d(I,J)=CrestLimit
      ENDIF ! end eta_limiter

# endif
    endif !end eta device kernel

! solve ubar
     if (i>=Ibeg .AND. i<=Iend .AND. j>=Jbeg .AND. j<=Jend)then
         F_left=Fx_d(i,j)
         F_right=Fx_d(i+1,J)
         F_bottom=Fy_d(i,j)
         F_top=Fy_d(i,j+1)
# if defined (CARTESIAN)
         R2=-1.0_SP/DX_d*(F_right*n_right+F_left*n_left) &
                       -1.0_SP/DY_d*(F_top*n_top+F_bottom*n_bottom) &
                        +SourceX_d(I,J)
# else
         R2=-1.0_SP/DX_d(i,j)*(F_right*n_right+F_left*n_left) &
                       -1.0_SP/DY_d(i,j)*(F_top*n_top+F_bottom*n_bottom) &
                        +SourceX_d(I,J)
# endif
# if defined (ITERATION)
         Ubar_d(I,J)=ALPHA_d(ISTEP)*Ubar0_d(I,J)+BETA_d(ISTEP)*(UbarOld_d(I,J)+DT*R2)
# else
         Ubar_d(I,J)=ALPHA_d(ISTEP)*Ubar0_d(I,J)+BETA_d(ISTEP)*(Ubar_d(I,J)+DT*R2)
# endif
     endif

! solve vbar
     if (i>=Ibeg .AND. i<=Iend .AND. j>=Jbeg .AND. j<=Jend )then
         F_left=Gx_d(i,j)
         F_right=Gx_d(i+1,J)
         F_bottom=Gy_d(i,j)
         F_top=Gy_d(i,j+1)
# if defined (CARTESIAN)
         R3=-1.0_SP/DX_d*(F_right*n_right+F_left*n_left) &
                       -1.0_SP/DY_d*(F_top*n_top+F_bottom*n_bottom) &
                       +SourceY_d(I,J)
# else
         R3=-1.0_SP/DX_d(i,j)*(F_right*n_right+F_left*n_left) &
                       -1.0_SP/DY_d(i,j)*(F_top*n_top+F_bottom*n_bottom) &
                       +SourceY_d(I,J)
# endif
# if defined (ITERATION)
         Vbar_d(I,J)=ALPHA_d(ISTEP)*Vbar0_d(I,J)+BETA_d(ISTEP)*(VbarOld_d(I,J)+DT*R3)
# else
         Vbar_d(I,J)=ALPHA_d(ISTEP)*Vbar0_d(I,J)+BETA_d(ISTEP)*(Vbar_d(I,J)+DT*R3)
# endif
     endif

!Move from Get_Eta_U_V.... from here
     if (i<=Mloc .AND. j<=Nloc) then
         H_d(i,j)=Eta_d(i,j)*Gamma3+Depth_d(i,j)
     endif

END SUBROUTINE HUV_BAR_KERNEL

!---------------------------------------------------------------------------------------!
!                                                                                       !
!   Cuda Kernel triDx_init_kernel is to prepare tridiagnals for tridiagnal solver in    !
!                                 x direction.                                          !
!         History : created on 04/14 2019                                               !
!                                                                                       !
!---------------------------------------------------------------------------------------!
attributes(global) subroutine triDx_init_kernel&
        (Ibeg,Iend,Jbeg,Jend,n_west,b1,b2,Gamma1,MinDepthFrc,WaveMakerCode)
    implicit none
!dummy variables
    integer,value :: Ibeg,Iend,Jbeg,Jend
    integer,value :: n_west
    real(SP),value :: b1,b2,Gamma1,MinDepthFrc
    integer,value :: WaveMakerCode
!local variables
    real(SP) :: dep,depl,depr
    real(SP) :: tmp1, tmp2, tmp3, tmp4
    real(SP) :: DXg
    integer::  I,J,tx,ty

! set local indexes
     tx = threadIdx%x
     ty = threadIdx%y
! Global thread and  block ID
     i = tx + (blockIdx%x-1)*BlockDimX_2D
     j = ty + (blockIdx%y-1)*BlockDimY_2D

    if (I>=Ibeg .AND. I<=Iend .AND. J>=Jbeg .AND. J<=Jend) then
# if defined (CARTESIAN)
        DXg=DX_d
# else
        DXg=DX_d(i,j)
# endif
        dep=Max(Depth_d(I,J),MinDepthFrc)
        depl=Max(Depth_d(I-1,J),MinDepthFrc)
        depr=Max(Depth_d(I+1,J),MinDepthFrc)
        !No difference between Cartesian and Spherical version if DX(I,J) intead of DXg
        tmp1=Gamma1*MASK9_d(I,J)*(b1/2.0_SP/DXg/DXg*dep*dep + b2/DXg/DXg*depl*dep)
        tmp2=1.0_SP+Gamma1*MASK9_d(I,J)*(-b1/DXg/DXg*dep*dep-2.0_SP*b2/DXg/DXg*dep*dep)
        tmp3=Gamma1*MASK9_d(I,J)*(b1/2.0_SP/DXg/DXg*dep*dep + b2/DXg/DXg*dep*depr)
        tmp4=Ubar_d(I,J)*MASK_d(I,J)/Max(H_d(I,J),MinDepthFrc)  &
             + Gamma1*MASK9_d(I,J)*( -b1/2.0_SP*dep*dep*Vxy_d(I,J)-b2*dep*DVxy_d(I,J))

# if defined (COUPLING)
        IF(IN_DOMAIN_WEST)THEN
            IF(I.eq.Ibeg)THEN
                tmp4=tmp4-tmp1*U_d(I-1,J)
            ENDIF
        ELSEIF(IN_DOMAIN_EAST)THEN
            IF(I.eq.Iend)THEN
                tmp4=tmp4-tmp3*U_d(I+1,J)
            ENDIF
        ENDIF
# endif

!  left_bc wavemaker
# if defined (MGPU)
        if(n_west.eq.MPI_PROC_NULL) then
# endif

        IF (WaveMakerCode==12)THEN
           IF(I.eq.Ibeg)THEN
             tmp4=tmp4-tmp1*U_d(I-1,J)
           ENDIF
        ENDIF

# if defined (MGPU)
        endif  
# endif
        IF(tmp2.NE.0.0_SP.OR.MASK_d(I,J).GT.0)THEN
           myAx_d(I,J)=tmp1/tmp2
           myCx_d(I,J)=tmp3/tmp2
           myDx_d(I,J)=tmp4/tmp2
        ELSE
           myAx_d(I,J)=ZERO
           myCx_d(I,J)=ZERO
           myDx_d(I,J)=ZERO
        ENDIF
    endif
end subroutine triDx_init_kernel

attributes(global) subroutine triDx_cusparse_init_kernel&
        (Ibeg,Iend,Jbeg,Jend,Mloc,n_west,NGhost,b1,b2,Gamma1,MinDepthFrc,WaveMakerCode)
    implicit none
!dummy variables
    integer,value :: Ibeg,Iend,Jbeg,Jend,Mloc,NGhost
    integer,value :: n_west
    real(SP),value :: b1,b2,Gamma1,MinDepthFrc
    integer,value :: WaveMakerCode
!local variables
    real(SP) :: dep,depl,depr
    real(SP) :: tmp1, tmp2, tmp3, tmp4
    real(SP) :: DXg
    integer::  I,J,tx,ty,RowIndex

! set local indexes
     tx = threadIdx%x
     ty = threadIdx%y
! Global thread and  block ID
     i = tx + (blockIdx%x-1)*BlockDimX_2D
     j = ty + (blockIdx%y-1)*BlockDimY_2D

    if (I>=Ibeg .AND. I<=Iend .AND. J>=Jbeg .AND. J<=Jend) then
# if defined (CARTESIAN)
        DXg=DX_d
# else
        DXg=DX_d(i,j)
# endif
        dep=Max(Depth_d(I,J),MinDepthFrc)
        depl=Max(Depth_d(I-1,J),MinDepthFrc)
        depr=Max(Depth_d(I+1,J),MinDepthFrc)
        !No difference between Cartesian and Spherical version if DX(I,J) intead of DXg
        tmp1=Gamma1*MASK9_d(I,J)*(b1/2.0_SP/DXg/DXg*dep*dep + b2/DXg/DXg*depl*dep)
        tmp2=1.0_SP+Gamma1*MASK9_d(I,J)*(-b1/DXg/DXg*dep*dep-2.0_SP*b2/DXg/DXg*dep*dep)
        tmp3=Gamma1*MASK9_d(I,J)*(b1/2.0_SP/DXg/DXg*dep*dep + b2/DXg/DXg*dep*depr)
        tmp4=Ubar_d(I,J)*MASK_d(I,J)/Max(H_d(I,J),MinDepthFrc)  &
             + Gamma1*MASK9_d(I,J)*( -b1/2.0_SP*dep*dep*Vxy_d(I,J)-b2*dep*DVxy_d(I,J))

# if defined (COUPLING)
        IF(IN_DOMAIN_WEST)THEN
            IF(I.eq.Ibeg)THEN
                tmp4=tmp4-tmp1*U_d(I-1,J)
            ENDIF
        ELSEIF(IN_DOMAIN_EAST)THEN
            IF(I.eq.Iend)THEN
                tmp4=tmp4-tmp3*U_d(I+1,J)
            ENDIF
        ENDIF
# endif

!  left_bc wavemaker
# if defined (MGPU)
        if(n_west.eq.MPI_PROC_NULL) then
# endif

        IF (WaveMakerCode==12)THEN
           IF(I.eq.Ibeg)THEN
             tmp4=tmp4-tmp1*U_d(I-1,J)
           ENDIF
        ENDIF

# if defined (MGPU)
        endif  
# endif
        RowIndex = (i-NGhost)+(j-NGhost-1)*(Mloc-2*NGhost)
        IF(tmp2.NE.0.0_SP.OR.MASK_d(I,J).GT.0)THEN
           Arow(RowIndex)=tmp1/tmp2
           Crow(RowIndex)=tmp3/tmp2
           Drow(RowIndex)=tmp4/tmp2
        ELSE
           Arow(RowIndex)=ZERO
           Crow(RowIndex)=ZERO
           Drow(RowIndex)=ZERO
        ENDIF
        if(i==Ibeg) Arow(RowIndex) = 0.0_SP
        if(i==Iend) Crow(RowIndex) = 0.0_SP
        
    endif
end subroutine triDx_cusparse_init_kernel

!---------------------------------------------------------------------------------------!
!                                                                                       !
!   Cuda Kernel triDy_init_kernel is to prepare tridiagnals for tridiagnal solver in    !
!                                 y direction.                                          !
!         History : created on 04/14 2019                                               !
!                                                                                       !
!---------------------------------------------------------------------------------------!

attributes(global) subroutine triDy_init_kernel&
        (Ibeg,Iend,Jbeg,Jend,b1,b2,Gamma1,Gamma2,MinDepthFrc,WaveMakerCode,DISP_TIME_LEFT)
    implicit none
!dummy variables
    integer,value :: Ibeg,Iend,Jbeg,Jend
    real(SP),value :: b1,b2,Gamma1,Gamma2,MinDepthFrc
    integer,value :: WaveMakerCode
    logical,value :: DISP_TIME_LEFT
!local variables
    real(SP) :: dep,depl,depr,reta,retal,retar
    real(SP) :: tmp1, tmp2, tmp3, tmp4
    real(SP) :: DYg
    integer :: I,J,tx,ty

! set local indexes
     tx = threadIdx%x
     ty = threadIdx%y
! Global thread and  block ID
     i = tx + (blockIdx%x-1)*BlockDimX_2D
     j = ty + (blockIdx%y-1)*BlockDimY_2D

    if (I>=Ibeg .AND. I<=Iend .AND. J>=Jbeg .AND. J<=Jend) then
# if defined (CARTESIAN)
        DYg=DY_d
# else
        DYg=DY_d(i,j)
# endif
! y direction
# if defined (MGPU)
        myAy_d(I,J)=ZERO
        myCy_d(I,J)=ZERO
        myDy_d(I,J)=ZERO
# endif
        dep=Max(Depth_d(I,J),MinDepthFrc)
        depl=Max(Depth_d(I,J-1),MinDepthFrc)
        depr=Max(Depth_d(I,J+1),MinDepthFrc)

# if defined (CARTESIAN)
        IF(DISP_TIME_LEFT)THEN
            reta=Eta_d(I,J)
            retal=Eta_d(I,J-1)
            retar=Eta_d(I,J+1)
            tmp1=Gamma1*MASK9_d(I,J)*(b1/2.0_SP/DYg/DYg*dep*dep + b2/DYg/DYg*depl*dep) &
                  -Gamma2*MASK9_d(I,J)*((reta+retal)*depl/2.0_SP/DYg/DYg+(retal+reta)*(retal+reta)/8.0_SP/DYg/DYg)
            tmp2=1.0_SP+Gamma1*MASK9_d(I,J)*(-b1/DYg/DYg*dep*dep-2.0_SP*b2/DYg/DYg*dep*dep) &
                  +Gamma2*MASK9_d(I,J)*((retar+retal+2.0_SP*reta)/2.0_SP/DYg/DYg &
                             +(retar*retar+2.0_SP*reta*reta &
                               +2.0_SP*reta*retar+2.0_SP*retal*reta+retal*retal)/8.0_SP/DYg/DYg)
            tmp3=Gamma1*MASK9_d(I,J)*(b1/2.0_SP/DYg/DYg*dep*dep + b2/DYg/DYg*dep*depr) &
                  -Gamma2*MASK9_d(I,J)*((reta+retar)*depr/2.0_SP/DYg/DYg+(retar+reta)*(retar+reta)/8.0_SP/DYg/DYg)
            tmp4=Vbar_d(I,J)*MASK_d(I,J)/Max(H_d(I,J),MinDepthFrc)  &
                  + Gamma1*MASK9_d(I,J)*(-b1/2.0_SP*dep*dep*Uxy_d(I,J)-b2*dep*DUxy_d(I,J)) &
                 + Gamma2*MASK9_d(I,J)*(reta*reta/2.0_SP*Uxy_d(I,J)+reta*DUxy_d(I,J) &
                                  + reta*ETAy_d(I,J)*Ux_d(I,J) + ETAy_d(I,J)*DUx_d(I,J) )
        ELSE
            tmp1=Gamma1*MASK9_d(I,J)*(b1/2.0_SP/DYg/DYg*dep*dep + b2/DYg/DYg*depl*dep) 
            tmp2=1.0_SP+Gamma1*MASK9_d(I,J)*(-b1/DYg/DYg*dep*dep-2.0_SP*b2/DYg/DYg*dep*dep) 
            tmp3=Gamma1*MASK9_d(I,J)*(b1/2.0_SP/DYg/DYg*dep*dep + b2/DYg/DYg*dep*depr)
            tmp4=Vbar_d(I,J)*MASK_d(I,J)/Max(H_d(I,J),MinDepthFrc)  &
                  + Gamma1*MASK9_d(I,J)*(-b1/2.0_SP*dep*dep*Uxy_d(I,J)-b2*dep*DUxy_d(I,J))
        ENDIF  
# else
        tmp1=Gamma1*MASK9_d(I,J)*(b1/2.0_SP/DYg/DYg*dep*dep + b2/DYg/DYg*depl*dep) 
        tmp2=1.0_SP+Gamma1*MASK9_d(I,J)*(-b1/DYg/DYg*dep*dep-2.0_SP*b2/DYg/DYg*dep*dep) 
        tmp3=Gamma1*MASK9_d(I,J)*(b1/2.0_SP/DYg/DYg*dep*dep + b2/DYg/DYg*dep*depr)
        tmp4=Vbar_d(I,J)*MASK_d(I,J)/Max(H_d(I,J),MinDepthFrc)  &
              + Gamma1*MASK9_d(I,J)*(-b1/2.0_SP*dep*dep*Uxy_d(I,J)-b2*dep*DUxy_d(I,J)) 
! remember to document this part for spherical version
# endif

# if defined (COUPLING)
        IF(IN_DOMAIN_SOUTH)THEN
          IF(J.eq.Jbeg)THEN
            tmp4=tmp4-tmp1*V_d(I,J-1)
          ENDIF
        ELSEIF(IN_DOMAIN_SOUTH)THEN
          IF(J.eq.Jend)THEN
            tmp4=tmp4-tmp3*V_d(I,J+1)
          ENDIF
        ENDIF
# endif
  
        IF(tmp2.NE.0.0_SP.OR.MASK_d(I,J).GT.0)THEN
            myAy_d(I,J)=tmp1/tmp2
            myCy_d(I,J)=tmp3/tmp2
            myDy_d(I,J)=tmp4/tmp2
        ELSE
            myAy_d(I,J)=ZERO
            myCy_d(I,J)=ZERO
            myDy_d(I,J)=ZERO
        ENDIF
    endif
end subroutine triDy_init_kernel

attributes(global) subroutine triDy_cusparse_init_kernel&
        (Ibeg,Iend,Jbeg,Jend,Nloc,NGhost,b1,b2,Gamma1,Gamma2,MinDepthFrc,WaveMakerCode,DISP_TIME_LEFT,PERIODIC)
    implicit none
!dummy variables
    integer,value :: Ibeg,Iend,Jbeg,Jend,NGhost,Nloc
    real(SP),value :: b1,b2,Gamma1,Gamma2,MinDepthFrc
    integer,value :: WaveMakerCode
    logical,value :: PERIODIC,DISP_TIME_LEFT
!local variables
    real(SP) :: dep,depl,depr,reta,retal,retar
    real(SP) :: tmp1, tmp2, tmp3, tmp4
    real(SP) :: DYg
    integer :: I,J,tx,ty,RowIndex

! set local indexes
     tx = threadIdx%x
     ty = threadIdx%y
! Global thread and  block ID
     i = tx + (blockIdx%x-1)*BlockDimX_2D
     j = ty + (blockIdx%y-1)*BlockDimY_2D

    if (I>=Ibeg .AND. I<=Iend .AND. J>=Jbeg .AND. J<=Jend) then
# if defined (CARTESIAN)
        DYg=DY_d
# else
        DYg=DY_d(i,j)
# endif
        RowIndex = (j-NGhost)+(i-NGhost-1)*(Nloc-2*NGhost)
! y direction
# if defined (MGPU)
        Arow(RowIndex)=ZERO
        Crow(RowIndex)=ZERO
        Drow(RowIndex)=ZERO
# endif
        dep=Max(Depth_d(I,J),MinDepthFrc)
        depl=Max(Depth_d(I,J-1),MinDepthFrc)
        depr=Max(Depth_d(I,J+1),MinDepthFrc)

# if defined (CARTESIAN)
        IF(DISP_TIME_LEFT)THEN
            reta=Eta_d(I,J)
            retal=Eta_d(I,J-1)
            retar=Eta_d(I,J+1)
            tmp1=Gamma1*MASK9_d(I,J)*(b1/2.0_SP/DYg/DYg*dep*dep + b2/DYg/DYg*depl*dep) &
                  -Gamma2*MASK9_d(I,J)*((reta+retal)*depl/2.0_SP/DYg/DYg+(retal+reta)*(retal+reta)/8.0_SP/DYg/DYg)
            tmp2=1.0_SP+Gamma1*MASK9_d(I,J)*(-b1/DYg/DYg*dep*dep-2.0_SP*b2/DYg/DYg*dep*dep) &
                  +Gamma2*MASK9_d(I,J)*((retar+retal+2.0_SP*reta)/2.0_SP/DYg/DYg &
                             +(retar*retar+2.0_SP*reta*reta &
                               +2.0_SP*reta*retar+2.0_SP*retal*reta+retal*retal)/8.0_SP/DYg/DYg)
            tmp3=Gamma1*MASK9_d(I,J)*(b1/2.0_SP/DYg/DYg*dep*dep + b2/DYg/DYg*dep*depr) &
                  -Gamma2*MASK9_d(I,J)*((reta+retar)*depr/2.0_SP/DYg/DYg+(retar+reta)*(retar+reta)/8.0_SP/DYg/DYg)
            tmp4=Vbar_d(I,J)*MASK_d(I,J)/Max(H_d(I,J),MinDepthFrc)  &
                  + Gamma1*MASK9_d(I,J)*(-b1/2.0_SP*dep*dep*Uxy_d(I,J)-b2*dep*DUxy_d(I,J)) &
                 + Gamma2*MASK9_d(I,J)*(reta*reta/2.0_SP*Uxy_d(I,J)+reta*DUxy_d(I,J) &
                                  + reta*ETAy_d(I,J)*Ux_d(I,J) + ETAy_d(I,J)*DUx_d(I,J) )
        ELSE
            tmp1=Gamma1*MASK9_d(I,J)*(b1/2.0_SP/DYg/DYg*dep*dep + b2/DYg/DYg*depl*dep) 
            tmp2=1.0_SP+Gamma1*MASK9_d(I,J)*(-b1/DYg/DYg*dep*dep-2.0_SP*b2/DYg/DYg*dep*dep) 
            tmp3=Gamma1*MASK9_d(I,J)*(b1/2.0_SP/DYg/DYg*dep*dep + b2/DYg/DYg*dep*depr)
            tmp4=Vbar_d(I,J)*MASK_d(I,J)/Max(H_d(I,J),MinDepthFrc)  &
                  + Gamma1*MASK9_d(I,J)*(-b1/2.0_SP*dep*dep*Uxy_d(I,J)-b2*dep*DUxy_d(I,J))
        ENDIF  
# else
        tmp1=Gamma1*MASK9_d(I,J)*(b1/2.0_SP/DYg/DYg*dep*dep + b2/DYg/DYg*depl*dep) 
        tmp2=1.0_SP+Gamma1*MASK9_d(I,J)*(-b1/DYg/DYg*dep*dep-2.0_SP*b2/DYg/DYg*dep*dep) 
        tmp3=Gamma1*MASK9_d(I,J)*(b1/2.0_SP/DYg/DYg*dep*dep + b2/DYg/DYg*dep*depr)
        tmp4=Vbar_d(I,J)*MASK_d(I,J)/Max(H_d(I,J),MinDepthFrc)  &
              + Gamma1*MASK9_d(I,J)*(-b1/2.0_SP*dep*dep*Uxy_d(I,J)-b2*dep*DUxy_d(I,J)) 
! remember to document this part for spherical version
# endif

# if defined (COUPLING)
        IF(IN_DOMAIN_SOUTH)THEN
          IF(J.eq.Jbeg)THEN
            tmp4=tmp4-tmp1*V_d(I,J-1)
          ENDIF
        ELSEIF(IN_DOMAIN_SOUTH)THEN
          IF(J.eq.Jend)THEN
            tmp4=tmp4-tmp3*V_d(I,J+1)
          ENDIF
        ENDIF
# endif
        IF(tmp2.NE.0.0_SP.OR.MASK_d(I,J).GT.0)THEN
            Arow(RowIndex)=tmp1/tmp2
            Crow(RowIndex)=tmp3/tmp2
            Drow(RowIndex)=tmp4/tmp2
        ELSE
            Arow(RowIndex)=ZERO
            Crow(RowIndex)=ZERO
            Drow(RowIndex)=ZERO
        ENDIF
        if (PERIODIC) then
!
        else
            if(j==Jbeg) Arow(RowIndex) = 0.0_SP
            if(j==Jend) Crow(RowIndex) = 0.0_SP
        endif
    endif
end subroutine triDy_cusparse_init_kernel


!---------------------------------------------------------------------------------------!
!                                                                                       !
!   Cuda Kernel triDx_forwardsweep_cuda_kernel is tridiagnal solver in x direction                   !
!                                                                                       !
!         History : created on 04/14 2019                                               !
!                                                                                       !
!---------------------------------------------------------------------------------------!

attributes(global) subroutine triDx_forwardsweep_cuda_kernel(A,C,D,&
                                        Mloc, Nloc,&
                                        Ibeg,Iend,&
                                        Jbeg,Jend)
    implicit none
    integer, value, intent(in) :: Mloc, Nloc, Ibeg, Iend, Jbeg, Jend
    real(SP), dimension(Nloc, Mloc),device :: A, C, D
    integer :: i,j

    j = threadIdx%x + (blockIdx%x-1)*blockDim%x   !Number of tridiagonal systems
    if ( j >= Jbeg .and. j <= Jend) then
        do i =Ibeg, Iend     ! Remember the ACDZ has been transposed (Ibeg-Iend now are the 2nd rank)
            if ( A(j,i) .NE. ZERO ) then
                C(j,i)=C(j,i)/A(j,i)/( 1.0_SP/A(j,i)-C(j,i-1) )
                D(j,i)=( D(j,i)/A(j,i)-D(j,i-1) ) / ( 1.0_SP/A(j,i)-C(j,i-1) )
            endif
        enddo
    endif
end subroutine triDx_forwardsweep_cuda_kernel

!---------------------------------------------------------------------------------------!
!                                                                                       !
!   Cuda Kernel triDx__cuda_kernel is tridiagnal solver in x direction                   !
!                                                                                       !
!         History : created on 04/14 2019                                               !
!                                                                                       !
!---------------------------------------------------------------------------------------!

attributes(global) subroutine triDx_backsubstitution_cuda_kernel(C,D,Z,&
                                        Mloc, Nloc,&
                                        Ibeg,Iend,&
                                        Jbeg,Jend)
    implicit none
    integer, value, intent(in) :: Mloc, Nloc, Ibeg, Iend, Jbeg, Jend
    real(SP), dimension(Nloc, Mloc),device :: C,D,Z
    integer :: i,j

    j = threadIdx%x + (blockIdx%x-1)*blockDim%x   !Number of tridiagonal systems
    if ( j >= Jbeg .and. j <= Jend) then
        !YUAN: Note that here Iend is actually IendTridx 
        do i = Iend, Ibeg, -1
            Z(j,i) = D(j,i)-C(j,i)*Z(j,i+1)
        enddo
    endif
end subroutine triDx_backsubstitution_cuda_kernel

!---------------------------------------------------------------------------------------!
!                                                                                       !
!   Cuda Kernel triDx_cuda_kernel is tridiagnal solver in x direction                   !
!                                                                                       !
!         History : created on 04/14 2019                                               !
!                                                                                       !
!---------------------------------------------------------------------------------------!

attributes(global) subroutine triDx_cuda_kernel(A,C,D,Z,&
                                        Mloc, Nloc,&
                                        Ibeg,Iend,&
                                        Jbeg,Jend)
    implicit none
    integer, value, intent(in) :: Mloc, Nloc, Ibeg, Iend, Jbeg, Jend
    real(SP), dimension(Nloc, Mloc), intent(inout),device :: A, C, D
    real(SP), dimension(Nloc, Mloc), intent(out),device :: Z
    integer :: i,j
    j = threadIdx%x + (blockIdx%x-1)*blockDim%x   !Number of tridiagonal systems
    if ( j >= Jbeg .and. j <= Jend) then
        do i =Ibeg+1, Iend     ! Remember the ACDZ has been transposed (Ibeg-Iend now are the 2nd rank)
            if ( A(j,i) .NE. ZERO ) then
                C(j,i)=C(j,i)/A(j,i)/( 1.0_SP/A(j,i)-C(j,i-1) )
                D(j,i)=( D(j,i)/A(j,i)-D(j,i-1) ) / ( 1.0_SP/A(j,i)-C(j,i-1) )
            endif
        enddo
        Z(j, Iend) = D(j,Iend)
        do i = Iend-1, Ibeg, -1
            Z(j,i) = D(j,i)-C(j,i)*Z(j,i+1)
        enddo
    endif
end subroutine triDx_cuda_kernel

!---------------------------------------------------------------------------------------!
!                                                                                       !
!   Cuda Kernel triDy_forwardsweep_cuda_kernel is tridiagnal solver in y direction                   !
!                                                                                       !
!         History : created on 04/04 2019                                               !
!                                                                                       !
!---------------------------------------------------------------------------------------!
! for triDy, not necessary to transpose the ACDZ.
attributes(global) subroutine triDy_forwardsweep_cuda_kernel(A,C,D,&
                                        Mloc, Nloc,&
                                        Ibeg,Iend,&
                                        Jbeg,Jend)
    implicit none
    integer, value, intent(in) :: Mloc, Nloc, Ibeg, Iend, Jbeg, Jend
    real(SP), dimension(Mloc, Nloc),device :: A, C, D
    integer :: i,j
    i = threadIdx%x + (blockIdx%x-1)*blockDim%x   !Number of tridiagonal systems
    if ( i >= Ibeg .and. i <= Iend) then
        do j =Jbeg, Jend
            if ( A(i,j) .NE. ZERO ) then
                C(i,j)=C(i,j)/A(i,j)/( 1.0_SP/A(i,j)-C(i,j-1) )
                D(i,j)=( D(i,j)/A(i,j)-D(i,j-1) ) / ( 1.0_SP/A(i,j)-C(i,j-1) )
            endif
        enddo
    endif
end subroutine triDy_forwardsweep_cuda_kernel

!---------------------------------------------------------------------------------------!
!                                                                                       !
!   Cuda Kernel triDy_backsubstitution__cuda_kernel is tridiagnal solver in y direction                   !
!                                                                                       !
!         History : created on 04/04 2019                                               !
!                                                                                       !
!---------------------------------------------------------------------------------------!
! for triDy, not necessary to transpose the ACDZ.
attributes(global) subroutine triDy_backsubstitution_cuda_kernel(C,D,Z,&
                                        Mloc, Nloc,&
                                        Ibeg,Iend,&
                                        Jbeg,Jend)
    implicit none
    integer, value, intent(in) :: Mloc, Nloc, Ibeg, Iend, Jbeg, Jend
    real(SP), dimension(Mloc, Nloc),device :: C,D,Z
    integer :: i,j
    i = threadIdx%x + (blockIdx%x-1)*blockDim%x   !Number of tridiagonal systems
    if ( i >= Ibeg .and. i <= Iend) then
        !YUAN: Jend-1 is substituted with Jend (JendTridy)
        do j = Jend, Jbeg, -1
            Z(i,j) = D(i,j)-C(i,j)*Z(i,j+1)
        enddo
    endif
end subroutine triDy_backsubstitution_cuda_kernel

!---------------------------------------------------------------------------------------!
!                                                                                       !
!   Cuda Kernel triDy_cuda_kernel is tridiagnal solver in y direction                   !
!                                                                                       !
!         History : created on 04/04 2019                                               !
!                                                                                       !
!---------------------------------------------------------------------------------------!
! for triDy, not necessary to transpose the ACDZ.
attributes(global) subroutine triDy_cuda_kernel(A,C,D,Z,&
                                        Mloc, Nloc,&
                                        Ibeg,Iend,&
                                        Jbeg,Jend)
    implicit none
    integer, value, intent(in) :: Mloc, Nloc, Ibeg, Iend, Jbeg, Jend
    real(SP), dimension(Mloc, Nloc), intent(inout),device :: A, C, D
    real(SP), dimension(Mloc, Nloc), intent(out),device   :: Z
    integer :: i,j
    i = threadIdx%x + (blockIdx%x-1)*blockDim%x   !Number of tridiagonal systems
    if ( i >= Ibeg .and. i <= Iend) then
        do j =Jbeg+1, Jend
            if ( A(i,j) .NE. ZERO ) then
                C(i,j)=C(i,j)/A(i,j)/( 1.0_SP/A(i,j)-C(i,j-1) )
                D(i,j)=( D(i,j)/A(i,j)-D(i,j-1) ) / ( 1.0_SP/A(i,j)-C(i,j-1) )
            endif
        enddo
        Z(i, Jend) = D(i,Jend)
        do j = Jend-1, Jbeg, -1
            Z(i,j) = D(i,j)-C(i,j)*Z(i,j+1)
        enddo
    endif
end subroutine triDy_cuda_kernel


attributes(global) subroutine transpose_kernel(idata,odata,n1,n2)
    implicit none
    integer, value, intent(IN) :: n1,n2
    real(SP), dimension(n2,n1), intent(IN)  :: idata
    real(SP), dimension(n1,n2), intent(OUT) :: odata
    integer :: i,j
    i = threadIdx%x + (blockIdx%x-1)*blockDim%x
    j = threadIdx%y + (blockIdx%y-1)*blockDim%y
    if(i<=n1 .and. j<=n2) then
        odata(i,j) = idata(j,i)
    end if
end subroutine transpose_kernel

!---------------------------------------------------------------------------------------!
!                                                                                       !
!   Cuda Kernel post_tridiagonal_kernel is to do some chore after tri-solver            !
!                                                                                       !
!         History : created on 04/04 2019                                               !
!                                                                                       !
!---------------------------------------------------------------------------------------!
attributes(global) subroutine post_tridiagonal_kernel(Ibeg,Iend,Jbeg,Jend,Mloc,Nloc,MinDepthFrc,FroudeCap,DISPERSION)
    implicit none
! Dummy variables     
    integer, value, intent(in) :: Ibeg,Iend,Jbeg,Jend,Mloc,Nloc
    real(SP), value :: MinDepthFrc,FroudeCap
    logical, value :: DISPERSION
! set local variables
    real(SP) :: Fr,Utotal,Utheta
    integer :: i,j
! Global thread and  block ID
    i = threadIdx%x + (BlockIdx%x-1)*BlockDim%x
    j = threadIdx%y + (BlockIdx%y-1)*BlockDim%y
     if (i>=Ibeg .AND. i<=Iend .AND. j>=Jbeg .AND. j<=Jend)then
         IF(MASK_d(I,J)<1)THEN
             Ubar_d(I,J)=ZERO
             Vbar_d(I,J)=ZERO
             U_d(I,J)=ZERO
             V_d(I,J)=ZERO
             HU_d(I,J)=ZERO
             HV_d(I,J)=ZERO
         ELSE
             HU_d(I,J)=Max(H_d(I,J),MinDepthFrc)*U_d(I,J)
             HV_d(I,J)=Max(H_d(I,J),MinDepthFrc)*V_d(I,J)
! apply Froude cap
             Utotal=SQRT(U_d(I,J)*U_d(I,J)+V_d(I,J)*V_d(I,J))
             Fr=SQRT(GRAV*Max(H_d(I,J),MinDepthFrc))
             IF(Utotal/Fr.gt.FroudeCap)THEN
                 Utheta=ATAN2(V_d(I,J),U_d(I,J))
                 U_d(I,J)=FroudeCap*Fr*COS(Utheta)
                 V_d(I,J)=FroudeCap*Fr*SIN(Utheta)
                 HU_d(I,J)=U_d(I,J)*Max(H_d(I,J),MinDepthFrc)
                 HV_d(I,J)=V_d(I,J)*Max(H_d(I,J),MinDepthFrc)
             ENDIF
! end Froude cap
         ENDIF
     endif
     if(.NOT. DISPERSION) then
        if (i<=Mloc .and. j<=Nloc) then
             Ubar_d(i,j)=HU_d(i,j)
             Vbar_d(i,j)=HV_d(i,j)
         endif
     endif
end subroutine post_tridiagonal_kernel

end module etauv_solver_module
