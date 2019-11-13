!-------------------------------------------------------------------------------------
!
!    WAVE_BREAKING is subroutine to do wave breaking 
!      using the artificial eddy viscosity scheme 
!      For the shock-capturing based breaking scheme 
!      this subroutine is only used for
!      demonstration or calculating bubbles or foam as 
!      SHOW_BREAKING = T.
!    
!    HISTORY: 
!      11/22/2010 Fengyan Shi
!      04/15/2015 Fengyan Shi, added viscosity breaking
!      09/19/2015 YoungKwang Chio, added viscosity distribution within WaveMaker region
!
!-------------------------------------------------------------------------------------
module wave_breaking_module
    use cudafor
    use PARAM
    use GLOBAL,only:MGlob,NGlob,NGhost,Ibeg,Iend,Jbeg,Jend,iista,jjsta,&
        Mloc,Nloc,MinDepthFrc,SHOW_BREAKING,&
        Xc_WK,Yc_WK,Width_WK,Ywidth_WK,WAVEMAKER_Cbrk,nu_bkg,T_brk,&
        Cbrk1,Cbrk2,WAVEMAKER_VIS,WAVEMAKER_visbrk,visbrk,DT
        
    use mod_cuda, only:DX_d, DY_d,H_d,Depth_d,ETAt_d,ETAx_d,ETAy_d,Eta_d,&
        AGE_BREAKING0_d,AGE_BREAKING_d,nu_break_d,grid,tBlock,&
        BlockDimX_2D,BlockDimY_2D,BlockDimX_Inner_2D,BlockDimY_Inner_2D

contains

SUBROUTINE WAVE_BREAKING_GPU
    IMPLICIT NONE

! variables for GPU kernels
    tBlock = dim3 (BlockDimX_2D, BlockDimY_2D ,1)
    grid = dim3 ( ceiling ( real ( Mloc ) / BlockDimX_Inner_2D ) , &
                  ceiling ( real ( Nloc ) / BlockDimY_Inner_2D ) , 1)
    IF(SHOW_BREAKING)THEN
        AGE_BREAKING0_d=AGE_BREAKING_d
        call breaking_kernel<<<grid, tBlock>>>&
            (Ibeg,Iend,Jbeg,Jend,iista,jjsta,MinDepthFrc,WAVEMAKER_Cbrk,&
            Xc_WK,Width_WK,Yc_WK,Ywidth_WK,nu_bkg,T_brk,Cbrk1,Cbrk2,DT)
    ELSEIF( WAVEMAKER_VIS )THEN
        grid = dim3 ( ceiling ( real ( Mloc ) / BlockDimX_2D ) , &
                  ceiling ( real ( Nloc ) / BlockDimY_2D ) , 1)
        call viscosity_wmaker_kernel<<<grid, tBlock>>>&
            (Mloc,Nloc,NGhost,Ibeg,Iend,Jbeg,Jend,iista,jjsta,&
            MinDepthFrc,WAVEMAKER_Cbrk,WAVEMAKER_visbrk,Xc_WK,Width_WK,Yc_WK,Ywidth_WK,&
            nu_bkg,visbrk)
    ENDIF

END SUBROUTINE WAVE_BREAKING_GPU

!-------------------------------------------------------------------------------------
!
!    VISCOSITY_WMAKER is subroutine 
!      to calculate viscosity inside wavemaker
!
!    HISTORY: 
!      08/19/2015 YoungKwang Choi
!
!-------------------------------------------------------------------------------------
attributes(global) subroutine viscosity_wmaker_kernel&
    (Mloc,Nloc,NGhost,Ibeg,Iend,Jbeg,Jend,iista,jjsta,&
    MinDepthFrc,WAVEMAKER_Cbrk,WAVEMAKER_visbrk,Xc_WK,Width_WK,Yc_WK,Ywidth_WK,&
    nu_bkg,visbrk)
    implicit none
    integer,value :: Mloc,Nloc,NGhost,Ibeg,Iend,Jbeg,Jend,iista,jjsta
    real(SP),value :: MinDepthFrc,WAVEMAKER_Cbrk,Xc_WK,Width_WK,Yc_WK,Ywidth_WK
    real(SP),value :: nu_bkg,visbrk,WAVEMAKER_visbrk
    REAL(SP) :: cap1
    real(SP) ::  tmp2,tmp3
    REAL(SP) :: xmk,ymk,DXg,DYg
    integer :: i,j,tx,ty
! set local indexes
    tx = threadIdx%x
    ty = threadIdx%y
! Global thread and  block ID
    i = tx + (blockIdx%x-1)*BlockDimX_2D    ! Inner must be smaller than BlockDim
    j = ty + (blockIdx%y-1)*BlockDimY_2D
    if (i>=Nghost+1 .AND. i<=Mloc-Nghost .AND. j>=Nghost+1 .AND. j<=Nloc-Nghost) then
        tmp3=SQRT(GRAV*MAX(MinDepthFrc,H_d(I,J)))
        tmp2=visbrk*tmp3

# if defined (CARTESIAN)
        DXg=DX_d
        DYg=DY_d
# else
        DXg=DX_d(I,J)
        DYg=DY_d(I,J)
# endif

! set viscosity
! wavemaker
# if defined (MGPU)
        xmk=(I-Ibeg)*DXg + (iista-1)*DXg
        ymk=(J-Jbeg)*DYg + (jjsta-1)*DYg
# else
        xmk=(I-Ibeg)*DXg
        ymk=(J-Jbeg)*DYg
# endif

! wavemaker doesnt use breaker age

        IF(ABS(xmk-Xc_WK)<Width_WK.AND. &
            ABS(ymk-Yc_WK)<Ywidth_WK/2.0_SP)THEN

            IF(ETAt_d(I,J)>MIN(tmp2,WAVEMAKER_visbrk*tmp3))THEN
                cap1=1.0*(MAX(Depth_d(I,J),MinDepthFrc)+ETA_d(I,J))
                nu_break_d(I,J)=cap1*WAVEMAKER_visbrk*tmp3+nu_bkg
            ELSE
                nu_break_d(I,J)=ZERO+nu_bkg
            ENDIF
          
        ENDIF ! end wavemaker
    endif
end subroutine viscosity_wmaker_kernel

!-------------------------------------------------------------------------------------
!
!    BREAKING is subroutine to calculate viscosity and breaking age
!    for the artificial eddy viscosity breaking scheme
!    
!    HISTORY: 
!     11/22/2010 Fengyan Shi
!     05/08/2019 AGE_BREAKING0 is added to store AGE_BREAKING at previous time step
!
!-------------------------------------------------------------------------------------
attributes(global) subroutine breaking_kernel&
    (Ibeg,Iend,Jbeg,Jend,iista,jjsta,MinDepthFrc,WAVEMAKER_Cbrk,&
    Xc_WK,Width_WK,Yc_WK,Ywidth_WK,nu_bkg,T_brk,Cbrk1,Cbrk2,DT)
    implicit none
    integer,value :: Ibeg,Iend,Jbeg,Jend,iista,jjsta
    real(SP),value :: MinDepthFrc,WAVEMAKER_Cbrk,Xc_WK,Width_WK,Yc_WK,Ywidth_WK
    real(SP),value :: nu_bkg,T_brk,Cbrk1,Cbrk2,DT
    REAL(SP) :: cap1,cap2
    REAL(SP) :: tmp1,tmp2,tmp3
    REAL(SP) :: xmk,ymk,DXg,DYg
    REAL(SP) :: etat_star,B,T_star
    REAL(SP)::C,Angle,AGE1,AGE2,AGE3,propx,propy,propxy
    real(SP),dimension(BlockDimX_2D,BlockDimY_2D),shared :: Age_sh
    integer :: i,j,tx,ty
! set local indexes
    tx = threadIdx%x
    ty = threadIdx%y
! Global thread and  block ID
    i = tx + (blockIdx%x-1)*BlockDimX_Inner_2D    ! Inner must be smaller than BlockDim
    j = ty + (blockIdx%y-1)*BlockDimY_Inner_2D
! init shared memory
! Important,add by YUAN : new var age_breaking0 to store previous values, 
! CPU version rewrite the age_breaking directly, which make parallel impossible
    if (i>=Ibeg-1 .AND. i<=Iend+1 .AND. j>=Jbeg-1 .AND. j<=Jend+1) then
        Age_sh(tx,ty) = AGE_BREAKING0_d(i,j)
    endif
    call syncthreads()
    !
    if (i>=Ibeg-1 .AND. i<=Iend+1 .AND. j>=Jbeg-1 .AND. j<=Jend+1 .and. &
        tx>1 .and. ty>1 .and. tx<blockDim%x .and. ty<blockDim%y) then
!
        tmp3=SQRT(GRAV*MAX(MinDepthFrc,H_d(I,J)))
        tmp1=Cbrk1*tmp3
        tmp2=Cbrk2*tmp3
        IF(ETAt_d(I,J).GE.tmp1.AND.(  &
            Age_sh(tx,ty).EQ.ZERO.OR.Age_sh(tx,ty).GT.T_brk))THEN
            AGE_BREAKING_d(I,J)=DT
        ELSE
            IF(Age_sh(tx,ty).GT.ZERO)THEN
                AGE_BREAKING_d(I,J)=Age_sh(tx,ty)+DT
            ELSE
                tmp1=MAX(SQRT(ETAx_d(I,J)*ETAx_d(I,J)+ETAy_d(I,J)*ETAy_d(I,J)),SMALL)
                C=MIN(ABS(ETAt_d(I,J))/tmp1,SQRT(GRAV*ABS(H_d(I,J))))
! propagation time between a dx, dy and ds
# if defined (CARTESIAN)
                DXg=DX_d
                DYg=DY_d
                propxy=SQRT(DXg*DXg+DYg*DYg)/MAX(C,SMALL)
                propx=SQRT(DXg*DXg)/MAX(C,SMALL)
                propy=SQRT(DYg*DYg)/MAX(C,SMALL)
# else
                DXg=DX_d(I,J)
                DYg=DY_d(I,J)
                propxy=SQRT(DXg*DXg+DYg*DYg)/MAX(C,SMALL)
                propx=SQRT(DXg*DXg)/MAX(C,SMALL)
                propy=SQRT(DYg*DYg)/MAX(C,SMALL)
# endif

!  I did two changes
!  1) conversion 180/pi
!  2) add - to get sh(tx+1,ty) quadrants
!  
                ANGLE=ATAN2(-ETAy_d(I,J),-ETAx_d(I,J))*180.0_SP/PI

                IF(ETAt_d(I,J).GE.tmp2)THEN
! 4 quadrants 
! quadrant 1
                    IF(ANGLE.GE.ZERO.AND.ANGLE.LT.90.0_SP)THEN
                       IF((Age_sh(tx-1,ty)>=DT.AND.Age_sh(tx-1,ty)>propx).OR.&
                          (Age_sh(tx-1,ty-1)>=DT.AND.Age_sh(tx-1,ty-1)>propxy).OR.&
                          (Age_sh(tx,ty-1)>=DT.AND.Age_sh(tx,ty-1)>propy))THEN
                           AGE_BREAKING_d(I,J)=DT
                       ENDIF         
                    ENDIF
! quadrant 2
                    IF(ANGLE.GE.90.0_SP.AND.ANGLE.LT.180.0_SP)THEN
                      IF((Age_sh(tx+1,ty)>=DT.AND.Age_sh(tx+1,ty)>propx).OR.&
                         (Age_sh(tx+1,ty-1)>=DT.AND.Age_sh(tx+1,ty-1)>propxy).OR.&
                         (Age_sh(tx,ty-1)>=DT.AND.Age_sh(tx,ty-1)>propy))THEN
                       AGE_BREAKING_d(I,J)=DT
                      ENDIF         
                    ENDIF
! quadrant 3
                    IF(ANGLE.GE.-180.0_SP.AND.ANGLE.LT.-90.0_SP)THEN
                      IF((Age_sh(tx+1,ty)>=DT.AND.Age_sh(tx+1,ty)>propx).OR.&
                         (Age_sh(tx+1,ty+1)>=DT.AND.Age_sh(tx+1,ty+1)>propxy).OR.&
                         (Age_sh(tx,ty+1)>=DT.AND.Age_sh(tx,ty+1)>propy))THEN
                          AGE_BREAKING_d(I,J)=DT
                      ENDIF         
                    ENDIF
! quadrant 4
                    IF(ANGLE.GE.-90.0_SP.AND.ANGLE.LT.0.0_SP)THEN
                      IF((Age_sh(tx,ty+1)>=DT.AND.Age_sh(tx,ty+1)>propy).OR.&
                         (Age_sh(tx-1,ty+1)>=DT.AND.Age_sh(tx-1,ty+1)>propxy).OR.&
                         (Age_sh(tx-1,ty)>=DT.AND.Age_sh(tx-1,ty)>propx))THEN
                          AGE_BREAKING_d(I,J)=DT
                      ENDIF         
                    ENDIF
                ENDIF
            ENDIF
        ENDIF
! set viscosity

! wavemaker
# if defined (MGPU)
![---ykchoi Jan/23/2018
        xmk=(I-Ibeg)*DXg + (iista-1)*DXg
        ymk=(J-Jbeg)*DYg + (jjsta-1)*DYg
# else
        xmk=(I-Ibeg)*DXg
        ymk=(J-Jbeg)*DYg
# endif

! wavemaker doesnt use breaker age

        IF(ABS(xmk-Xc_WK)<Width_WK.AND. &
            ABS(ymk-Yc_WK)<Ywidth_WK/2.0_SP)THEN

! remove wavebreaking inside wavemaker

            IF(ETAt_d(I,J)>MIN(tmp2,WAVEMAKER_Cbrk*tmp3))THEN
                cap1=1.0*(MAX(Depth_d(I,J),MinDepthFrc)+ETA_d(I,J))
                nu_break_d(I,J)=cap1*WAVEMAKER_Cbrk*tmp3+nu_bkg
            ELSE
                nu_break_d(I,J)=ZERO+nu_bkg
            ENDIF

        ELSE ! outside wavemaker    

            IF(Age_sh(tx,ty)>ZERO.AND.Age_sh(tx,ty)<T_brk.AND.ETAt_d(I,J)>tmp2)THEN
                cap1=1.0*(MAX(Depth_d(I,J),MinDepthFrc)+ETA_d(I,J))
!   
# if defined (VIS_KENNEDY_ORIGINAL)

! etat_star
                etat_star=tmp2
                T_star=5.0_SP*SQRT(MAX(Depth_d(I,J),MinDepthFrc)/GRAV)
                IF(Age_sh(tx,ty)>=0.0_SP.AND.Age_sh(tx,ty)<T_star)THEN
                    etat_star=tmp1+Age_sh(tx,ty)/T_star*(tmp2-tmp1)
                ENDIF
! B
                B = 0.0_SP
                IF(ETAt_d(I,J)>etat_star.AND.ETAt_d(I,J)<=2.0_SP*etat_star)THEN
                    B = ETAt_d(I,J)/etat_star-1.0_SP
                ELSEIF(ETAt_d(I,J)>2.0_SP*etat_star)THEN
                    B = 1.0_SP
                ENDIF
        
                nu_break_d(I,J)=cap1*ABS(ETAt_d(I,J))*B+nu_bkg

# elif defined (VIS_KENNEDY)

! this is a modified Kennedy breaker, I removed etat_star because the TVD is more stable
! I also modified the transition, the original of which caused instability

                B = 0.0_SP
                IF(ETAt_d(I,J)>tmp1.AND.ETAt_d(I,J)<=2.0_SP*tmp1)THEN
                    B = ETAt_d(I,J)/tmp1-1.0_SP
                ELSEIF(ETAt_d(I,J)>2.0_SP*tmp1)THEN
                    B = 1.0_SP
                ENDIF

                nu_break_d(I,J) = cap1*tmp2*(1.0_SP+B) + nu_bkg
       

# elif defined (VIS_STATIC_TRANS)
                nu_break_d(I,J)=cap1*tmp3*(Cbrk2+(Cbrk1-Cbrk2)*(T_brk-Age_sh(tx,ty))/T_brk)+nu_bkg
# else
                nu_break_d(I,J) = cap1*tmp2 + nu_bkg
# endif


            ELSE
                nu_break_d(I,J)=ZERO+nu_bkg
            ENDIF
          
        ENDIF ! end wavemaker

   endif
end subroutine breaking_kernel

end module wave_breaking_module
