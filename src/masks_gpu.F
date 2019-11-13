!-------------------------------------------------------------------------------------
!
!    UPDATE_MASK_GPU is subroutine to update mask at GPU device
!       note that mask also be updated in fluxes subroutine
!    HISTORY: 
!       04/28/2019 created
!-------------------------------------------------------------------------------------
module update_mask_module
    use cudafor
    use PARAM, only : SP
    use GLOBAL,only :Ibeg,Iend,Jbeg,Jend,Mloc,Nloc,MGlob,NGlob,NGhost, &
        MinDepth,MinDepthFrc,SWE_ETA_DEP,VISCOSITY_BREAKING
    use mod_cuda, only: MASK_d,MASKtmp_d,MASK_STRUC_d,MASK9_d,Eta_d,ETAtmp_d, Depth_d,DepthX_d,DepthY_d,&
        BlockDimX_2D,BlockDimY_2D,BlockDimX_Inner_2D,BlockDimY_Inner_2D,&
        grid,tBlock
# if defined (MGPU)
    use mgpu_utilities
# endif

contains
attributes(global) subroutine update_mask_kernel&
        (MinDepth,Ibeg,Iend,Jbeg,Jend)
    implicit none
    real(SP),value :: MinDepth
    integer, value :: Ibeg,Iend,Jbeg,Jend
    integer :: i,j
! Global thread and  block ID
    I = threadIdx%x + (blockIdx%x-1)*BlockDimX_2D
    J = threadIdx%y + (blockIdx%y-1)*BlockDimY_2D
    if (i>=Ibeg-2 .AND. i<=Iend+2 .AND. j>=Jbeg-2 .AND. j<=Jend+2) then 
! flood
        IF(MASK_STRUC_d(I,J)==1)THEN
            IF(MASK_d(I,J)<1)THEN
                ! left
                IF(MASK_d(I-1,J)==1.AND.Eta_d(I-1,J)>Eta_d(I,J))THEN
                    MASKtmp_d(I,J)=1
                ENDIF
                ! right
                IF(MASK_d(I+1,J)==1.AND.Eta_d(I+1,J)>Eta_d(I,J))THEN
                    MASKtmp_d(I,J)=1
                ENDIF
                ! bottom
                IF(MASK_d(I,J-1)==1.AND.Eta_d(I,J-1)>Eta_d(I,j))THEN
                    MASKtmp_d(I,J)=1
                ENDIF
                ! top
                IF(MASK_d(I,J+1)==1.AND.Eta_d(I,J+1)>Eta_d(I,j))THEN
                    MASKtmp_d(I,J)=1
                ENDIF
! drying
            ELSE
                IF(Eta_d(I,J)<-Depth_d(I,J))THEN
                    MASKtmp_d(I,J)=0
                    ETAtmp_d(I,J)=MinDepth-Depth_d(I,J)
                ENDIF
            ENDIF    
        ENDIF
!TO DO by YUAN:Now I can not find suitable solution for IGNORE_SLOPE (in nature serially)
! try to solve by 1-d do loop in each gpu kernel
!# if defined (IGNORE_BIG_SLOPE)
!    ! do not truncate depthx and depthy
!# else
!        IF(MASK_d(I,J)<1)THEN
!            DepthX_d(I,J)=Depth_d(I-1,J)
!            DepthX_d(I+1,J)=Depth_d(I+1,J)
!            DepthY_d(I,J)=Depth_d(I,J-1)
!            DepthY_d(I,J+1)=Depth_d(I,J+1)
!        ENDIF  
!# endif  
    endif
!
end subroutine update_mask_kernel

attributes(global) subroutine ignore_big_slope_kernel(Ibeg,Iend,Jbeg,Jend)
    implicit none
    integer :: i,j
    integer,value :: Ibeg,Iend,Jbeg,Jend
! Global thread and  block ID
     i = threadIdx%x + (blockIdx%x-1)*BlockDimX_2D
    if (i>=Ibeg-2 .AND. i<=Iend+2) then 
        do j = Jbeg-2,Jend+2
            IF(MASK_d(I,J)<1)THEN
                DepthY_d(I,J)=Depth_d(I,J-1)
                DepthY_d(I,J+1)=Depth_d(I,J+1)
            ENDIF
        enddo
    endif
    if (i>=Jbeg-2 .AND. i<=Jend+2) then 
        do j = Ibeg-2,Iend+2
            IF(MASK_d(J,I)<1)THEN
                DepthX_d(j,i)=Depth_d(j-1,i)
                DepthX_d(j+1,i)=Depth_d(j+1,i)
            ENDIF
        enddo
    endif
end subroutine ignore_big_slope_kernel

attributes(global) subroutine update_mask9_kernel(MinDepthFrc,SWE_ETA_DEP,Ibeg,Iend,Jbeg,Jend)
    implicit none
    real(SP),value :: MinDepthFrc,SWE_ETA_DEP
    integer, value :: Ibeg,Iend,Jbeg,Jend
    integer, dimension(BlockDimX_2D,BlockDimY_2D), shared :: MASK_sh
    integer :: i,j,tx,ty
! set local indexes
     tx = threadIdx%x
     ty = threadIdx%y
! Global thread and  block ID
     i = tx + (blockIdx%x-1)*BlockDimX_Inner_2D
     j = ty + (blockIdx%y-1)*BlockDimY_Inner_2D
! init shared memory
    if (i>=Ibeg-2 .AND. i<=Iend+2.AND. j>=Jbeg-2 .AND. j<=Jend+2) then 
        MASK_sh(tx,ty) = MASK_d(i,j)
    endif
    call syncthreads()
    if (i>=Ibeg-1 .AND. i<=Iend+1 .AND. j>=Jbeg-1 .AND. j<=Jend+1 .AND.&
        tx>1 .and. ty>1 .and. tx<BlockDimX_2D .and. ty<BlockDimY_2D) then
            MASK9_d(I,J)=MASK_sh(tx,ty)*MASK_sh(tx-1,ty)*MASK_sh(tx+1,ty)  &
                *MASK_sh(tx+1,ty+1)*MASK_sh(tx,ty+1)*MASK_sh(tx-1,ty+1) &
                *MASK_sh(tx+1,ty-1)*MASK_sh(tx,ty-1)*MASK_sh(tx-1,ty-1) 
            IF(ABS(Eta_d(I,J))/MAX(DEPTH_d(I,J),MinDepthFrc)>SWE_ETA_DEP)THEN
                MASK9_d(I,J)=0
            ENDIF
    endif
end subroutine update_mask9_kernel


SUBROUTINE UPDATE_MASK_GPU
     IMPLICIT NONE
     integer :: istat
     integer ::ierrSync, ierrAsync
!
# if defined (MGPU)
    call phi_int_exch_cuda(MASK_d)
    call phi_int_exch_cuda(MASK9_d)
# endif
!
     MASKtmp_d=MASK_d
     ETAtmp_d=ETA_d

     tBlock = dim3 (BlockDimX_2D, BlockDimY_2D ,1)
     grid = dim3 ( ceiling ( real ( Mloc ) / tBlock%x ) , &
            ceiling ( real ( Nloc ) / tBlock%y ) , 1)
     call update_mask_kernel<<<grid,tBlock>>>(MinDepth,Ibeg,Iend,Jbeg,Jend)
     !barrier to make sure MASKtmp_d is ready to d2d copy
     istat = cudaDeviceSynchronize() 
# if defined (IGNORE_BIG_SLOPE)
    ! do not truncate depthx and depthy
# else
     tBlock = dim3 (BlockDimX_2D,1 ,1)
     grid = dim3 ( ceiling ( real ( max(Mloc,Nloc) ) / tBlock%x ) , &
            1 , 1)
    call ignore_big_slope_kernel<<<grid,tBlock>>>(Ibeg,Iend,Jbeg,Jend)
    istat = cudaDeviceSynchronize() 
# endif  
     MASK_d = MASKtmp_d
     ETA_d = ETAtmp_d
     IF(VISCOSITY_BREAKING)THEN
         ! dont use mask9 for viscosity breaking
     ELSE
         tBlock = dim3 (BlockDimX_2D, BlockDimY_2D ,1)
         grid = dim3 ( ceiling ( real ( Mloc ) / BlockDimX_Inner_2D ) , &
             ceiling ( real ( Nloc ) / BlockDimY_Inner_2D ) , 1)
         call update_mask9_kernel<<<grid,tBlock>>>(MinDepthFrc,SWE_ETA_DEP,Ibeg,Iend,Jbeg,Jend) 
!        ierrSync = cudaGetLastError()
!        ierrAsync = cudaDeviceSynchronize()
!        if (ierrSync /= cudaSuccess) print *, &
!           'Sync kernel error update mask9:', cudaGetErrorString(ierrSync)
!        if (ierrAsync /= cudaSuccess) print *, &
!            'Async kernel error update mask9:', cudaGetErrorString(ierrAsync)
    ENDIF
!
# if defined (MGPU)
    call phi_int_exch_cuda(MASK_d)
    call phi_int_exch_cuda(MASK9_d)
# endif

END SUBROUTINE UPDATE_MASK_GPU

end module update_mask_module
