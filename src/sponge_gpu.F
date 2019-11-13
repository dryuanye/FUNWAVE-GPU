!-------------------------------------------------------------------------------------
!
!    SPONGE_DAMPING is subroutine for dampping waves using Larsen-Dancy(1983)
!    type sponge layer 
!    
!    HISTORY: 10/27/2010 Fengyan Shi
!
!-------------------------------------------------------------------------------------
module sponge_module
    use cudafor
    use PARAM, only: ZERO
    use GLOBAL, only: MGlob,NGlob,NGhost, &
                Mloc,Nloc
    use mod_cuda, ONLY: tBlock,grid,BlockDimX_2D,BlockDimY_2D,ETA_d,U_d,V_d,SPONGE_d,MASK_d

contains

attributes(global) subroutine sponge_damping_kernel(Mloc,Nloc)
    implicit none
    integer, value :: Mloc, Nloc
    integer :: i,j
! Global thread and  block ID
    i = threadIdx%x + (blockIdx%x-1)*BlockDimX_2D
    j = threadIdx%y + (blockIdx%y-1)*BlockDimY_2D
    if (i>=1 .AND. i<=Mloc .AND. j>=1 .AND. j<=Nloc) then
        IF(MASK_d(I,J)>0)THEN
            ETA_d(I,J)=ETA_d(I,J)/SPONGE_d(I,J)
        ENDIF
        U_d(I,J)=U_d(I,J)/SPONGE_d(I,J)
        V_d(I,J)=V_d(I,J)/SPONGE_d(I,J)
    endif
end subroutine sponge_damping_kernel

SUBROUTINE SPONGE_DAMPING_GPU
    implicit none
! variables for GPU kernels
    tBlock = dim3 (BlockDimX_2D, BlockDimY_2D ,1)
    grid = dim3 ( ceiling ( real ( Mloc ) / BlockDimX_2D ) , &
                  ceiling ( real ( Nloc ) / BlockDimY_2D ) , 1)
    call sponge_damping_kernel<<<grid, tBlock>>>(Mloc,Nloc)

END SUBROUTINE SPONGE_DAMPING_GPU

end module sponge_module
