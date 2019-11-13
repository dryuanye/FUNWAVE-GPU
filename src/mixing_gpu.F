!-------------------------------------------------------------------------------------
!
!    MIXING_STUFF is subroutine to calculate mixing related, time-averaged properties
!    mean eta is also calculated.
!    
!    HISTORY: 05/02/2011 Fengyan Shi
!
!-------------------------------------------------------------------------------------
module calculate_mean_module
    use cudafor
    use PARAM, only: ZERO,SP
    use GLOBAL, only: Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,T_INTV_mean,&
        T_sum,time,STEADY_TIME,DT,OUT_Umean,OUT_Vmean,OUT_ETAmean,OUT_WaveHeight,&
        Umean,Vmean,ETAmean,WaveHeightRMS,WaveHeightAve,SigWaveHeight
    use mod_cuda, only: U_d,V_d,ETA_d,ETA0_d,Umean_d,Vmean_d, &
        WaveHeightRMS_d,WaveHeightAve_d,ETAmean_d,Usum_d,Vsum_d,ETAsum_d, &
        ETA2sum_d,ETA2mean_d,SigWaveHeight_d,&
        HrmsSum_d,HavgSum_d,Emax_d,Emin_d,Num_Zero_Up_d,&
        BlockDimX_2D,BlockDimY_2D,grid,tBlock

contains
!-------------------------------------------------------------------------------------
!
!    CALCULATE_MEAN is subroutine to calculate mean u v required by 
!      smagorinsky mixing and wave height
!      mean eta is also calculated.
!    
!    HISTORY: 
!      05/02/2011 Fengyan Shi
!                 Young-Kwang Choi added some time-averaging stuff
!
!-------------------------------------------------------------------------------------
attributes(global) subroutine calc_mean_kernel(Mloc,Nloc,T_sum,DT)
    implicit none
    real(SP), value :: T_sum,DT
    integer, value :: Mloc, Nloc
    integer :: i,j
    real(SP) :: DTg
! Global thread and  block ID
    i = threadIdx%x + (blockIdx%x-1)*BlockDimX_2D
    j = threadIdx%y + (blockIdx%y-1)*BlockDimY_2D
    DTg = DT
    if (i>=1 .AND. i<=Mloc .AND. j>=1 .AND. j<=Nloc) then
        ETA2sum_d(i,j) = (Eta_d(i,j)-ETAmean_d(i,j))*(Eta_d(i,j)-ETAmean_d(i,j))*DTg + ETA2sum_d(i,j)
        ETA2mean_d(i,j) = ETA2sum_d(i,j)/T_sum
        Usum_d(i,j)=U_d(i,j)*DTg+Usum_d(i,j)
        Vsum_d(i,j)=V_d(i,j)*DTg+Vsum_d(i,j)
        ETAsum_d(i,j)=ETA_d(i,j)*DTg+ETAsum_d(i,j)
        Umean_d(i,j)=Usum_d(i,j)/T_sum
        Vmean_d(i,j)=Vsum_d(i,j)/T_sum
        ETAmean_d(i,j)=ETAsum_d(i,j)/T_sum
        SigWaveHeight_d(i,j) = 4.004*SQRT( ETA2mean_d(i,j) )  
        IF(Num_Zero_Up_d(I,J)>=2)THEN
            WaveHeightAve_d(I,J)=HavgSum_d(I,J)/Num_Zero_Up_d(I,J)
            WaveHeightRMS_d(I,J)=SQRT(HrmsSum_d(I,J)/Num_Zero_Up_d(I,J))
        ENDIF
    endif
end subroutine calc_mean_kernel


attributes(global) subroutine calc_sum_kernel(Mloc,Nloc,DT)
    implicit none
    REAL(SP)::Tmpe,Tmp_0
    integer, value :: Mloc, Nloc
    real(SP),value :: DT
    integer :: i,j
    real(SP) :: DTg
! Global thread and  block ID
    i = threadIdx%x + (blockIdx%x-1)*BlockDimX_2D
    j = threadIdx%y + (blockIdx%y-1)*BlockDimY_2D
    DTg = DT
    if (i>=1 .AND. i<=Mloc .AND. j>=1 .AND. j<=Nloc) then
        Usum_d(i,j)=U_d(i,j)*DTg+Usum_d(i,j)
        Vsum_d(i,j)=V_d(i,j)*DTg+Vsum_d(i,j)
        ETAsum_d(i,j)=ETA_d(i,j)*DTg+ETAsum_d(i,j)
        ETA2sum_d(i,j) = (Eta_d(i,j)-ETAmean_d(i,j))*(Eta_d(i,j) &
                -ETAmean_d(i,j))*DTg + ETA2sum_d(i,j)
! wave height
        if(Eta_d(i,j)>Emax_d(i,j)) Emax_d(i,j) = Eta_d(i,j)
        if(Eta_d(i,j)<Emin_d(i,j)) Emin_d(i,j) = Eta_d(i,j)
        Tmpe = Eta_d(i,j)-ETAmean_d(i,j)
        Tmp_0 = Eta0_d(i,j)-ETAmean_d(i,j)
        if(Tmpe>Tmp_0.and.Tmpe*Tmp_0<=Zero) then
            Num_Zero_Up_d(i,j) = Num_Zero_Up_d(i,j)+1
            if(Num_Zero_Up_d(i,j)>=2) then
                HavgSum_d(i,j) = HavgSum_d(i,j)+Emax_d(i,j)-Emin_d(i,j)
                HrmsSum_d(i,j) = HrmsSum_d(i,j)+(Emax_d(i,j)-Emin_d(i,j))**2
            endif
            ! reset Emax and Emin to find next wave
            Emax_d(i,j) = -1000.
            Emin_d(i,j) = 1000.
        endif  
    endif
end subroutine calc_sum_kernel

subroutine var_d2h_mean
    implicit none
    IF(OUT_Umean)THEN
      Umean = Umean_d
    ENDIF
    IF(OUT_Vmean)THEN
      Vmean = Vmean_d
    ENDIF
    IF(OUT_ETAmean)THEN
      ETAmean = ETAmean_d
    ENDIF
    IF(OUT_WaveHeight)THEN
      WaveHeightRMS=WaveHeightRMS_d
      WaveHeightAve=WaveHeightAve_d
      SigWaveHeight=SigWaveHeight_d
    ENDIF
end subroutine var_d2h_mean


SUBROUTINE CALCULATE_MEAN_GPU
! calculate mean for smagorinsky s mixing and wave height
    IMPLICIT NONE
! variables for GPU kernels
    tBlock = dim3 (BlockDimX_2D, BlockDimY_2D ,1)
    grid = dim3 ( ceiling ( real ( Mloc ) / BlockDimX_2D ) , &
                  ceiling ( real ( Nloc ) / BlockDimY_2D ) , 1)
!
    T_sum=T_sum+DT

    IF(T_sum.GE.T_INTV_mean)THEN
        call calc_mean_kernel<<<grid, tBlock>>>(Mloc,Nloc,T_sum,DT)
!
        T_sum=T_sum-T_INTV_mean 
        Usum_d=ZERO
        Vsum_d=ZERO
        ETAsum_d=ZERO
        ETA2sum_d=ZERO 
!
        !first download from device to host, then call preview_mean at host
        call var_d2h_mean
        CALL PREVIEW_MEAN
    ELSE
        call calc_sum_kernel<<<grid,tBlock>>>(Mloc,Nloc,DT)
    ENDIF  ! end average time

END SUBROUTINE CALCULATE_MEAN_GPU

SUBROUTINE MIXING_STUFF_GPU
     IMPLICIT NONE
! calculate mean for smagorinsky s mixing and wave height
     IF( time >= STEADY_TIME )THEN    
         CALL CALCULATE_MEAN_GPU
     ENDIF    

END SUBROUTINE MIXING_STUFF_GPU


end module calculate_mean_module
