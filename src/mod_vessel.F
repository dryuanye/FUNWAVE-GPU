!------------------------------------------------------------------------------------
!
!      FILE vessle.F
!
!      This file is part of the FUNWAVE-TVD program under the Simplified BSD license
!
!-------------------------------------------------------------------------------------
! 
!    Copyright (c) 2016, FUNWAVE Development Team
!
!    (See http://www.udel.edu/kirby/programs/funwave/funwave.html
!     for Development Team membership)
!
!    All rights reserved.
!
!    FUNWAVE_TVD is free software: you can redistribute it and/or modify
!    it under the terms of the Simplified BSD License as released by
!    the Berkeley Software Distribution (BSD).
!
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions are met:
!
!    1. Redistributions of source code must retain the above copyright notice, this
!       list of conditions and the following disclaimer.
!    2. Redistributions in binary form must reproduce the above copyright notice,
!    this list of conditions and the following disclaimer in the documentation
!    and/or other materials provided with the distribution.
!
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
!    ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
!    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
!    ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
!    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
!    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!  
!    The views and conclusions contained in the software and documentation are those
!    of the authors and should not be interpreted as representing official policies,
!    either expressed or implied, of the FreeBSD Project.
!  
!-------------------------------------------------------------------------------------
!
!  VESSEL is a module to model ship-wakes    
!
!  HISTORY :
!    10/07/2016  Fengyan Shi
!
!-------------------------------------------------------------------------------------

# if defined (VESSEL)

MODULE VESSEL_MODULE
  USE PARAM
  USE GLOBAL,ONLY : Mloc,Nloc,Nghost,Ibeg,Iend,Jbeg,Jend,DX,DY, &
                    H,ETA,ETA0,Xco,Yco
  USE INPUT_READ
#if defined (PARALLEL)
  USE GLOBAL,ONLY : myid,ier, npx,npy,PX,PY
  USE MPI
# endif
  IMPLICIT NONE
  SAVE

    INTEGER :: NumVessel,Kves
    REAL(SP),DIMENSION(:),ALLOCATABLE :: Xvessel1,Yvessel1,Xvessel2,Yvessel2, &
                                       LengthVessel,WidthVessel, &
                                       AlphaVessel,BetaVessel,Pvessel, &
                                       TimeVessel1,TimeVessel2,ThetaVessel
    REAL(SP),DIMENSION(:,:),ALLOCATABLE :: VesselPressureTotal,VesselPressureEach, &
                                       VesselPressureX,VesselPressureY
    LOGICAL :: OUT_VESSEL = .TRUE.
    REAL(SP),DIMENSION(:),ALLOCATABLE :: ResistanceX,ResistanceY,ResPosX,ResNegX, &
                                         ResPosY,ResNegY
    REAL(SP):: PLOT_INTV_VESSEL,PLOT_COUNT_VESSEL
# if defined (VESSEL_PANEL_SOURCE)
    REAL(SP),DIMENSION(:,:),ALLOCATABLE :: VesselFluxGradient,VesselFluxGradientEach, &
                                           Gz
    REAL(SP),DIMENSION(:),ALLOCATABLE :: Uvel,Vvel
# endif

#if defined (PARALLEL)
    REAL(SP) :: myvar
# endif    
    LOGICAL :: MakeVesselDraft = .TRUE.


!INTERFACE READ_FOUR_TYPE_VALUES
!  Module Procedure VESSEL_INITIAL
!  Module Procedure VESSEL_FORCING
!END INTERFACE

CONTAINS
  
! READ VESSEL

SUBROUTINE VESSEL_INITIAL
  USE GLOBAL,ONLY : itmp1,itmp2,itmp3,itmp4,itmp5,SMALL
# if defined (PARALLEL)
  USE GLOBAL,ONLY : iista,jjsta   !ykchoi Jan/23/2018
# endif
                    
  USE INPUT_READ
  IMPLICIT NONE
  CHARACTER(LEN=80)::FILE_NAME=' '
  CHARACTER(LEN=80)::FILE_VESSEL=' '
  CHARACTER(LEN=80) :: VESSEL_FOLDER=' '
  CHARACTER(LEN=80)::TMP_NAME=' '
  INTEGER :: Ifile,ierr
  CHARACTER(LEN=80):: VesselName

! read vessel number and folder from input.txt
      FILE_NAME='input.txt'

! vessel folder
      CALL READ_STRING(VESSEL_FOLDER,FILE_NAME,'VESSEL_FOLDER',ierr)
# if defined (PARALLEL)
      if (myid.eq.0) WRITE(3,'(A15,A50)')'VESSEL_FOLDER:', VESSEL_FOLDER
# else
      WRITE(3,'(A15,A50)')'VESSEL_FOLDER:', VESSEL_FOLDER
# endif

      CALL READ_INTEGER(NumVessel,FILE_NAME,'NumVessel',ierr)
# if defined (PARALLEL)
      if (myid.eq.0) WRITE(3,'(A12,I3)') 'NumVessel = ',NumVessel
# else
      WRITE(3,'(A12,I3)') 'NumVessel = ',NumVessel
# endif

      CALL READ_LOGICAL(OUT_VESSEL,FILE_NAME,'OUT_VESSEL',ierr)

      ALLOCATE (Xvessel1(NumVessel),Yvessel1(NumVessel),  &
                Xvessel2(NumVessel),Yvessel2(NumVessel),  &
                TimeVessel1(NumVessel),TimeVessel2(NumVessel), &
                LengthVessel(NumVessel),WidthVessel(NumVessel), &
                AlphaVessel(NumVessel),BetaVessel(NumVessel),  &
                Pvessel(NumVessel),ThetaVessel(NumVessel),&
                ResistanceX(NumVessel),ResistanceY(NumVessel), &
                ResPosX(NumVessel),ResNegX(NumVessel), &
                ResPosY(NumVessel),ResNegY(NumVessel) )


# if defined (VESSEL_PANEL_SOURCE)

      ALLOCATE (VesselFluxGradient(Mloc,Nloc),VesselFluxGradientEach(Mloc,Nloc), &
                Gz(Mloc,Nloc))
      ALLOCATE (Uvel(NumVessel),Vvel(NumVessel))
      
      Gz = ZERO
      Uvel = ZERO
      Vvel = ZERO

# else

      ALLOCATE (VesselPressureTotal(Mloc,Nloc), VesselPressureEach(Mloc,Nloc),&
                VesselPressureX(Mloc,Nloc), &
                 VesselPressureY(Mloc,Nloc) )

# endif


! plot vessel intitial
     PLOT_COUNT_VESSEL = 0
     CALL READ_FLOAT(PLOT_INTV_VESSEL,FILE_NAME,'PLOT_INTV_VESSEL',ierr)
     IF(ierr==1)THEN
# if defined (PARALLEL)
      if (myid.eq.0) WRITE(3,'(A50)')'PLOT_INTV_VESSEL not specified, use SMALL'
# else
      WRITE(3,'(A50)')'PLOT_INTV_VESSEL not specified, use SMALL'
# endif
       PLOT_INTV_VESSEL = SMALL
     ENDIF
      
  DO Kves = 1, NumVessel

!  file name
    itmp1=mod(Kves/10000,10)
    itmp2=mod(Kves/1000,10)
    itmp3=mod(Kves/100,10)
    itmp4=mod(Kves/10,10)
    itmp5=mod(Kves,10)
    write(FILE_VESSEL(1:1),'(I1)')itmp1
    write(FILE_VESSEL(2:2),'(I1)')itmp2
    write(FILE_VESSEL(3:3),'(I1)')itmp3
    write(FILE_VESSEL(4:4),'(I1)')itmp4
    write(FILE_VESSEL(5:5),'(I1)')itmp5

    TMP_NAME = TRIM(VESSEL_FOLDER)//'vessel_'//TRIM(FILE_VESSEL)

! check existing

 INQUIRE(FILE=TRIM(TMP_NAME),EXIST=FILE_EXIST)
  IF(.NOT.FILE_EXIST)THEN
# if defined (PARALLEL)
   IF(MYID==0)  &
   WRITE(*,*) TRIM(TMP_NAME), ' specified in ', TRIM(VESSEL_FOLDER), ' but CANNOT BE FOUND. STOP'
   CALL MPI_FINALIZE (ier)
   STOP
# else
    WRITE(*,*) TRIM(TMP_NAME), ' specified in ', TRIM(VESSEL_FOLDER), ' but CANNOT BE FOUND. STOP'
    STOP
# endif
  ENDIF

! open file
  Ifile=Kves+200
  OPEN(Ifile,FILE=TRIM(TMP_NAME))

! read file
         READ(Ifile,*)  ! title
         READ(Ifile,*)  VesselName  ! vessel name
         READ(Ifile,*)  ! length and width
         READ(Ifile,*)  LengthVessel(Kves), WidthVessel(Kves), &
                      AlphaVessel(Kves),BetaVessel(Kves),Pvessel(Kves)
         READ(Ifile,*)  ! t, x, y
         READ(Ifile,*)  TimeVessel2(Kves),Xvessel2(Kves),Yvessel2(Kves)

         TimeVessel1(Kves) = TimeVessel2(Kves)
         Xvessel1(Kves) = Xvessel2(Kves)
         Yvessel1(Kves) = Yvessel2(Kves)

         AlphaVessel(Kves) = Max(SMALL, AlphaVessel(Kves))
         BetaVessel(Kves) = Max(SMALL, BetaVessel(Kves))
         AlphaVessel(Kves) = Min(1.0_SP, AlphaVessel(Kves))
         BetaVessel(Kves) = Min(1.0_SP, BetaVessel(Kves))

# if defined (PARALLEL)
   IF(MYID==0)THEN
   WRITE(3,*) 'Vessel Name: ',  TRIM(VesselName)
   WRITE(3,*) 'Vessel Length', LengthVessel(Kves)
   WRITE(3,*) 'Vessel Width', WidthVessel(Kves)
   WRITE(3,*) 'Vessel Alpha', AlphaVessel(Kves)
   WRITE(3,*) 'Vessel Beta', BetaVessel(Kves)
   WRITE(3,*) 'Vessel P', PVessel(Kves)
   WRITE(3,*) 'Initial Time, X, Y', TimeVessel2(Kves),Xvessel2(Kves),Yvessel2(Kves)
   ENDIF
# else
   WRITE(3,*) 'Vessel Name: ',  TRIM(VesselName)
   WRITE(3,*) 'Vessel Length', LengthVessel(Kves)
   WRITE(3,*) 'Vessel Width', WidthVessel(Kves)
   WRITE(3,*) 'Vessel Alpha', AlphaVessel(Kves)
   WRITE(3,*) 'Vessel Beta', BetaVessel(Kves)
   WRITE(3,*) 'Vessel P', PVessel(Kves)
   WRITE(3,*) 'Initial Time, X, Y', TimeVessel2(Kves),Xvessel2(Kves),Yvessel2(Kves)
# endif

  ENDDO  ! end Kves

End SUBROUTINE VESSEL_INITIAL

SUBROUTINE VESSEL_FORCING
  USE GLOBAL,ONLY : Mloc,Nloc,tmp1,tmp2,SMALL,TIME,ZERO
  USE INPUT_READ
  IMPLICIT NONE
  INTEGER :: Ifile,ierr,I,J
  REAL(SP) :: Xves,Yves,Lves,Wves,Px,Py,DetaX,DetaY

# if defined (VESSEL_PANEL_SOURCE)
  VesselFluxGradient = ZERO
# else
  VesselPressureTotal = ZERO
# endif

  DO Kves = 1,NumVessel

    IF(TIME>TimeVessel1(Kves).AND.TIME>TimeVessel2(Kves)) THEN

         TimeVessel1(Kves)=TimeVessel2(Kves)
         Xvessel1(Kves) = Xvessel2(Kves)
         Yvessel1(Kves) = Yvessel2(Kves)

    Ifile = 200 + Kves

    READ(Ifile,*,END=120)  TimeVessel2(Kves),Xvessel2(Kves),Yvessel2(Kves)

# if defined (PARALLEL)
   IF(MYID==0)THEN
     WRITE(3,*)'Read Vessel # ', Kves
     WRITE(3,*)'T,X,Y = ', TimeVessel2(Kves),Xvessel2(Kves),Yvessel2(Kves)
   ENDIF
# else
     WRITE(3,*)'Read Vessel # ', Kves
     WRITE(3,*)'T,X,Y = ', TimeVessel2(Kves),Xvessel2(Kves),Yvessel2(Kves)
# endif

    ThetaVessel(Kves) = ATAN2(Yvessel2(Kves)-Yvessel1(Kves),  &
                              Xvessel2(Kves)-Xvessel1(Kves))


# if defined (VESSEL_PANEL_SOURCE)
    IF((TimeVessel2(Kves)-TimeVessel1(Kves))> ZERO)THEN
     Uvel(Kves) = (Xvessel2(Kves)-Xvessel1(Kves))/(TimeVessel2(Kves)-TimeVessel1(Kves))
     Vvel(Kves) = (Yvessel2(Kves)-Yvessel1(Kves))/(TimeVessel2(Kves)-TimeVessel1(Kves))
    ENDIF
# endif

    ENDIF ! end time > timevessel2

! calculate force
    tmp2=ZERO
    tmp1=ZERO

    IF(TIME>TimeVessel1(Kves))THEN
      IF(TimeVessel1(Kves).EQ.TimeVessel2(Kves))THEN
        ! no more data
        tmp2=ZERO
        tmp1=ZERO
      ELSE
      tmp2=(TimeVessel2(Kves)-TIME) &
            /MAX(SMALL, ABS(TimeVessel2(Kves)-TimeVessel1(Kves)))
      tmp1=1.0_SP - tmp2;
      ENDIF  ! no more data?
    ENDIF ! time>time_1

    Xves = Xvessel2(Kves)*tmp1 +Xvessel1(Kves)*tmp2
    Yves = Yvessel2(Kves)*tmp1 +Yvessel1(Kves)*tmp2

# if defined (VESSEL_PANEL_SOURCE)

     CALL GREEN_FUNCTION_SOURCE (Xves,Yves)

# else  
   ! pressure source

! rectangular
    VesselPressureEach = ZERO
    ResistanceX(Kves)=ZERO
    ResistanceY(Kves)=ZERO
    ResPosX(Kves) = ZERO
    ResNegX(Kves) = ZERO
    ResPosY(Kves) = ZERO
    ResNegY(Kves) = ZERO
    DO J=1,Nloc
    DO I=1,Mloc
      Lves=(Xco(I)-Xves)*COS(ThetaVessel(Kves)) + (Yco(J)-Yves)*SIN(ThetaVessel(Kves))
      Wves=-(Xco(I)-Xves)*SIN(ThetaVessel(Kves)) + (Yco(J)-Yves)*COS(ThetaVessel(Kves))
# if defined (TORSVIK)
      
      IF(ABS(Lves)<=0.5_SP*LengthVessel(Kves).AND. &
         ABS(Wves)<=0.5_SP*WidthVessel(Kves)) THEN
         VesselPressureEach(I,J) = Pvessel(Kves)  &
                  *COS(PI*Lves/(LengthVessel(Kves)))**2 &
                  *COS(PI*Wves/(WidthVessel(Kves)))**2                  
      ENDIF
# elif defined (TMP)

    IF(ABS(Lves)<=0.5_SP*LengthVessel(Kves).AND. &
         ABS(Wves)<=0.5_SP*WidthVessel(Kves)) THEN

      Px = ZERO
      Py = ZERO

! I multiplide 0.01 here 
     IF(Lves>ZERO)THEN
      IF(ABS(Lves)>0.5_SP*LengthVessel(Kves)*AlphaVessel(Kves)*0.01.AND. &
         ABS(Lves)<0.5_SP*LengthVessel(Kves))THEN
          Px = COS(PI*(Lves-0.5_SP*AlphaVessel(Kves)*0.01*LengthVessel(Kves))  &
                   /((1.0_SP-AlphaVessel(Kves)*0.01)*LengthVessel(Kves)))**2
      ELSEIF(ABS(Lves)<=0.5_SP*LengthVessel(Kves)*AlphaVessel(Kves)*0.1)THEN
          Px = 1.0_SP
      ENDIF
     ENDIF

     IF(Lves<=ZERO)THEN
      IF(ABS(Lves)>0.5_SP*LengthVessel(Kves)*AlphaVessel(Kves).AND. &
         ABS(Lves)<0.5_SP*LengthVessel(Kves))THEN
          Px = COS(PI*(Lves-0.5_SP*AlphaVessel(Kves)*LengthVessel(Kves))  &
                   /((1.0_SP-AlphaVessel(Kves))*LengthVessel(Kves)))**2
      ELSEIF(ABS(Lves)<=0.5_SP*LengthVessel(Kves)*AlphaVessel(Kves))THEN
          Px = 1.0_SP
      ENDIF
     ENDIF

      IF(ABS(Wves)>0.5_SP*WidthVessel(Kves)*BetaVessel(Kves).AND. &
         ABS(Wves)<0.5_SP*WidthVessel(Kves)) THEN
          Py = COS(PI*(Wves-0.5_SP*BetaVessel(Kves)*WidthVessel(Kves))  &
                  /((1.0_SP-BetaVessel(Kves))*WidthVessel(Kves)))**2
      ELSEIF(ABS(Wves)<=0.5_SP*WidthVessel(Kves)*BetaVessel(Kves)) THEN
          Py = 1.0_SP
      ENDIF

         VesselPressureEach(I,J) = Pvessel(Kves)*Px*Py

    ENDIF  ! end inside ship rectangule

# else
   ! Ertekin et al. JFM 1986

    IF(ABS(Lves)<=0.5_SP*LengthVessel(Kves).AND. &
         ABS(Wves)<=0.5_SP*WidthVessel(Kves)) THEN

      Px = ZERO
      Py = ZERO

      IF(ABS(Lves)>0.5_SP*LengthVessel(Kves)*AlphaVessel(Kves).AND. &
         ABS(Lves)<0.5_SP*LengthVessel(Kves))THEN
          
          Px = COS(PI*(Lves-0.5_SP*AlphaVessel(Kves)*LengthVessel(Kves))  &
                   /((1.0_SP-AlphaVessel(Kves))*LengthVessel(Kves)))**2
      ELSEIF(ABS(Lves)<=0.5_SP*LengthVessel(Kves)*AlphaVessel(Kves))THEN
          Px = 1.0_SP
      ENDIF
      IF(ABS(Wves)>0.5_SP*WidthVessel(Kves)*BetaVessel(Kves).AND. &
         ABS(Wves)<0.5_SP*WidthVessel(Kves)) THEN
          Py = COS(PI*(Wves-0.5_SP*BetaVessel(Kves)*WidthVessel(Kves))  &
                  /((1.0_SP-BetaVessel(Kves))*WidthVessel(Kves)))**2
      ELSEIF(ABS(Wves)<=0.5_SP*WidthVessel(Kves)*BetaVessel(Kves)) THEN
          Py = 1.0_SP
      ENDIF

         VesselPressureEach(I,J) = Pvessel(Kves)*Px*Py

    ENDIF  ! end inside ship rectangule

# endif
     
    ENDDO
    ENDDO  ! end grid

! calculate resistance. In the initial code, the resistance was calculated
! inside the last loop. It turns out that resistance was re-counted 
! because of including ghost cells (thanks to Jeff Harris)    

    DO J=Jbeg,Jend
    DO I=Ibeg,Iend
# if defined (TORSVIK)
    ! do nothing
# else

         DetaX = (ETA(I+1,J)-ETA(I-1,J))/2.0_SP
         DetaY = (ETA(I,J+1)-ETA(I,J-1))/2.0_SP

!         ResistanceX(Kves)=ResistanceX(Kves)  &
!                 +VesselPressureEach(I,J)*RHO_WATER*GRAV  &
!                  *DetaX*DY

!         ResistanceY(Kves)=ResistanceY(Kves)  &
!                 +VesselPressureEach(I,J)*RHO_WATER*GRAV  &
!                  *DetaY*DX

         IF(DetaX>=0)THEN
           ResPosX(Kves) = ResPosX(Kves) &
                  +VesselPressureEach(I,J)*RHO_WATER*GRAV  &
                  *DetaX*DY
         ELSE
           ResNegX(Kves) = ResNegX(Kves) &
                  +VesselPressureEach(I,J)*RHO_WATER*GRAV  &
                  *DetaX*DY           
         ENDIF

         IF(DetaY>=0)THEN
           ResPosY(Kves) = ResPosY(Kves) &
                  +VesselPressureEach(I,J)*RHO_WATER*GRAV  &
                  *DetaY*DX
         ELSE
           ResNegY(Kves) = ResNegY(Kves) &
                  +VesselPressureEach(I,J)*RHO_WATER*GRAV  &
                  *DetaY*DX          
         ENDIF

         ResistanceX(Kves) = ResPosX(Kves) + ResNegX(Kves)
         ResistanceY(Kves) = ResPosY(Kves) + ResNegY(Kves)
      
# endif
 ! end torsvik
    ENDDO
    ENDDO


# endif
  ! end difference sources

120 CONTINUE  ! no more data for vessel Kves


# if defined (VESSEL_PANEL_SOURCE)
     VesselFluxGradient = VesselFluxGradient + VesselFluxGradientEach
# else
  ! pressure source

# if defined (PARALLEL)
     call MPI_ALLREDUCE(ResistanceX(Kves),myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier)
     ResistanceX(Kves) = myvar
     call MPI_ALLREDUCE(ResistanceY(Kves),myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier)
     ResistanceY(Kves) = myvar
     call MPI_ALLREDUCE(ResPosX(Kves),myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier)
     ResPosX(Kves) = myvar
     call MPI_ALLREDUCE(ResPosY(Kves),myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier)
     ResPosY(Kves) = myvar
     call MPI_ALLREDUCE(ResNegX(Kves),myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier)
     ResNegX(Kves) = myvar
     call MPI_ALLREDUCE(ResNegY(Kves),myvar,1,MPI_SP,MPI_SUM,MPI_COMM_WORLD,ier)
     ResNegY(Kves) = myvar
# endif

    VesselPressureTotal = VesselPressureTotal+VesselPressureEach

# endif
  ! end different sources

  ENDDO  ! end Kves

! sourceX and sourceY


# if defined (VESSEL_PANEL_SOURCE)
   ! do nothing about pressure gradient for panel source
# else
    DO J=Jbeg,Jend
    DO I=Ibeg,Iend

!   I modified the term to negative 11/22/2016 fyshi

       VesselPressureX(I,J) = -Grav*H(I,J)*  &
               (VesselPressureTotal(I+1,J)-VesselPressureTotal(I-1,J))/2.0_SP  &
               /DX
       VesselPressureY(I,J) = -Grav*H(I,J)*  &
               (VesselPressureTotal(I,J+1)-VesselPressureTotal(I,J-1))/2.0_SP  &
               /DY

    ENDDO
    ENDDO

# endif
  ! end panel or pressure sources


!  make initial draft 09/12/2018
    IF(MakeVesselDraft)THEN
       MakeVesselDraft = .FALSE.
        DO J=1,Nloc
        DO I=1,Mloc
          IF(ABS(VesselPressureTotal(I,J)).GT.SMALL)THEN      
            Eta(I,J) = - VesselPressureTotal(I,J)
            ETA0(I,J) = ETA(I,J)
          ENDIF
        ENDDO
        ENDDO
    ENDIF

END SUBROUTINE VESSEL_FORCING

# if defined (VESSEL_PANEL_SOURCE)
SUBROUTINE GREEN_FUNCTION_SOURCE (Xves,Yves)
  USE GLOBAL,ONLY : Mloc,Nloc,tmp1,tmp2,SMALL,TIME,ZERO, H,U,V,  &
                    DX,DY,Ibeg,Iend,Jbeg,Jend,DEPTH
  IMPLICIT NONE
  INTEGER :: ierr,I,J
  REAL(SP) :: Xves,Yves,Lves,Wves


  VesselFluxGradientEach = ZERO

# if defined (REALISTIC_VESSEL_BODY)
  REAL, DIMENSION(Mloc,Nloc) :: W_farfield,Z_prime
                                ! Z_prime defined 0 at surface and positive below
  INTEGER :: m_trunc, km
  REAL :: r3,r3_2m,r3_4m

  m_trunc = 5
  Gz = ZERO
     
  DO J=Jbeg,Jend
  DO I=Ibeg,Iend 

      Lves=(Xco(I)-Xves)*COS(ThetaVessel(Kves)) + (Yco(J)-Yves)*SIN(ThetaVessel(Kves))
      Wves=-(Xco(I)-Xves)*SIN(ThetaVessel(Kves)) + (Yco(J)-Yves)*COS(ThetaVessel(Kves))

!     the 3D nearfield is solved inside -L to L and -W to W
!     note that this length doubles the vessel length
!     the far field is outside this rectangule 
      
      IF(ABS(Lves)>=LengthVessel(Kves).AND. &
         ABS(Wves)>=WidthVessel(Kves)) THEN

        Z_prime(I,J) = Pvessel(Kves)
        W_farfield = - U(I,J)*(Depth(I+1,J)-Depth(I-1,J))/DX*0.5_SP &
                     - V(I,J)*(Depth(I,J+1)-Depth(I,J-1))/DY*0.5_SP &
                     - (Depth(I,J)-Z_prime(I,J)) &
                     * (U(I+1,J)-U(I-1,J))/DX*0.5_SP &
                      + V(I,J+1)-V(I,J-1))/DY*0.5_SP)

          r3=((Xco(I,J)-Xves)**2+(Yco(I,J)-Yves)**2+Z_prime**2)**(3.0_SP/2.0_SP)
          Gz(I,J) = -2.0_SP*Z_prime(I,J)/MAX(SMALL,r3)
                               ! actually r3 should larger than zero
        DO km = 1, m_trunc

          r3_2m=((Xco(I,J)-Xves)**2+(Yco(I,J)-Yves)**2  &
                 +(Z_prime+2.0_SP*km*H(I,J))**2)**(3.0_SP/2.0_SP)
          r3_4m=((Xco(I,J)-Xves)**2+(Yco(I,J)-Yves)**2  &
                 +(Z_prime-2.0_SP*km*H(I,J))**2)**(3.0_SP/2.0_SP)
          Gz(I,J) = Gz(I,J) -2.0_SP * (-1.0_SP)**km &
                    * ( (Z_prime(I,J) + 2.0_SP*km*Depth(I,J))/r3_2m &
                       +(Z_prime(I,J) + 2.0_SP*km*Depth(I,J))/r3_4m )
        ENDDO

         VesselFluxGradientEach(I,J) = -1.0_SP /2.0_SP /PI *( (Uvel(Kves)-U(I,J))/DX &
                               (Uvel(Kves)-U(I,J))/DY -W_farfield/Depth(I,J) )
                               
     ENDIF  ! end inside vessel body

  ENDDO
  ENDDO

# else
  

! assume single value of submerged vessel body
! we only record Gz at z_prime

    DO J=1,Nloc
    DO I=1,Mloc
      Lves=(Xco(I)-Xves)*COS(ThetaVessel(Kves)) + (Yco(J)-Yves)*SIN(ThetaVessel(Kves))
      Wves=-(Xco(I)-Xves)*SIN(ThetaVessel(Kves)) + (Yco(J)-Yves)*COS(ThetaVessel(Kves))

      
      IF(ABS(Lves)<=0.5_SP*LengthVessel(Kves).AND. &
         ABS(Wves)<=0.5_SP*WidthVessel(Kves)) THEN
         VesselFluxGradientEach(I,J) = Pvessel(Kves)  &
                  *SIN(2.0*PI*Lves/(LengthVessel(Kves))) &
                  *COS(PI*Wves/(WidthVessel(Kves)))**2   
!   print*,i,j,VesselFluxGradientEach(I,J), SIN(2.0*PI*Lves/(LengthVessel(Kves))), &
!          COS(PI*Wves/(WidthVessel(Kves)))**2
      ENDIF
    ENDDO
    ENDDO

# endif
  !  end realistic vessel body, otherwise slender

!open(99,file='tmp.txt')
!do j=1,Nloc
! write(99,108)(VesselFluxGradientEach(i,j),i=1,Mloc)
!enddo
!close(99)
!108 format(5000E12.3)
!stop

END SUBROUTINE GREEN_FUNCTION_SOURCE
# endif 
  ! end panel source


END MODULE VESSEL_MODULE

# endif 
! end vessel
