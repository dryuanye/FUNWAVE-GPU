# if defined (MGPU)
module mgpu_utilities
    use cudafor
    USE PARAM
    USE GLOBAL
    use mod_cuda, only: istat
    integer,save :: HaloWE,HaloSN,HaloWE_int,HaloSN_int,YType
    interface phi_exch_cuda
        module procedure phi_exch_cuda
        module procedure phi_int_exch_cuda
    end interface phi_exch_cuda

contains
subroutine mpi_datatype
    implicit none
    call MPI_TYPE_VECTOR(Nloc,Nghost,Mloc,MPI_SP,HaloWE,ier)
    call MPI_TYPE_COMMIT(HaloWE,ier)
    call MPI_TYPE_VECTOR(Nghost,Mloc,Mloc,MPI_SP,HaloSN,ier)
    call MPI_TYPE_COMMIT(HaloSN,ier)
    call MPI_TYPE_VECTOR(Nloc,Nghost,Mloc,MPI_INTEGER,HaloWE_int,ier)
    call MPI_TYPE_COMMIT(HaloWE_int,ier)
    call MPI_TYPE_VECTOR(Nghost,Mloc,Mloc,MPI_INTEGER,HaloSN_int,ier)
    call MPI_TYPE_COMMIT(HaloSN_int,ier)
    call MPI_TYPE_VECTOR(1,Mloc,Mloc,MPI_SP,YType,ier)
    call MPI_TYPE_COMMIT(YType,ier)
end subroutine mpi_datatype

!SUBROUTINE phi_exch_cuda (PHI)
!    IMPLICIT NONE
!    REAL(SP),INTENT(INOUT),DEVICE :: PHI(Mloc,Nloc)
!
!    !REAL(SP),DIMENSION(Mloc,Nghost),DEVICE :: rNmsg, sNmsg,rSmsg,sSmsg
!    !REAL(SP),DIMENSION(Nghost,Nloc),DEVICE :: rWmsg, sWmsg,rEmsg,sEmsg
!
!! for east-west
!
!    if ( n_west .ne. MPI_PROC_NULL ) then
!       call MPI_SENDRECV( PHI(Ibeg,1),1,HaloWE,n_west, 1,& 
!                       PHI(Ibeg-Nghost,1),1,HaloWE,n_west,0,comm2d,status,ier)
!    endif
!
!    if ( n_east .ne. MPI_PROC_NULL ) then
!       call MPI_SENDRECV( PHI(Iend-Nghost+1,1), 1, HaloWE,n_east, 0,&
!                          PHI(Iend+1,1),1,HaloWE,n_east,1,comm2d,status,ier )
!    endif
!
!! for nrth-suth
!
!    if ( n_suth .ne. MPI_PROC_NULL ) then
!       call MPI_SENDRECV( PHI(1,Jbeg),1,HaloSN, n_suth, 1,&
!                          PHI(1,Jbeg-Nghost),1,HaloSN,n_suth, 0, comm2d,status, ier )
!    endif
!
!    if ( n_nrth .ne. MPI_PROC_NULL ) then
!       call MPI_SENDRECV( PHI(1,Jend-Nghost+1),1,HaloSN,n_nrth, 0,&
!                       PHI(1,Jend+1),1,HaloSN, n_nrth, 1, comm2d,status, ier )
!    endif
!
!END SUBROUTINE phi_exch_cuda
!!-------------------------------------------------------------------------------------
!!
!!    phi_exch_cuda is the subroutine to exchange GPU device variable at processor
!!          interface
!!    
!!    HISTORY: 02/14/2011 Jeff Harris
!!             05/01/2011 Fengyan Shi, implemented into the TVD code
!!
!!-------------------------------------------------------------------------------------
!
!SUBROUTINE phi_int_exch_cuda (PHI)
!    IMPLICIT NONE
!    INTEGER,INTENT(INOUT),DEVICE :: PHI(Mloc,Nloc)
!
!    !INTEGER,DIMENSION(Mloc,Nghost), DEVICE :: rNmsg, sNmsg,rSmsg,sSmsg
!    !INTEGER,DIMENSION(Nghost,Nloc), DEVICE :: rWmsg, sWmsg,rEmsg,sEmsg
!
!! for east-west
!
!    if ( n_west .ne. MPI_PROC_NULL ) then
!       call MPI_SENDRECV( PHI(Ibeg,1),1,HaloWE_int,n_west, 1,& 
!                       PHI(Ibeg-Nghost,1),1,HaloWE_int,n_west,0,comm2d,status,ier)
!    endif
!
!    if ( n_east .ne. MPI_PROC_NULL ) then
!       call MPI_SENDRECV( PHI(Iend-Nghost+1,1), 1, HaloWE_int,n_east, 0,&
!                          PHI(Iend+1,1),1,HaloWE_int,n_east,1,comm2d,status,ier )
!    endif
!
!! for nrth-suth
!
!    if ( n_suth .ne. MPI_PROC_NULL ) then
!       call MPI_SENDRECV( PHI(1,Jbeg),1,HaloSN_int, n_suth, 1,&
!                          PHI(1,Jbeg-Nghost),1,HaloSN_int,n_suth, 0, comm2d,status, ier )
!    endif
!
!    if ( n_nrth .ne. MPI_PROC_NULL ) then
!       call MPI_SENDRECV( PHI(1,Jend-Nghost+1),1,HaloSN_int,n_nrth, 0,&
!                       PHI(1,Jend+1),1,HaloSN_int, n_nrth, 1, comm2d,status, ier )
!    endif
!END SUBROUTINE phi_int_exch_cuda

!==================================================================================

SUBROUTINE phi_exch_cuda (PHI)
    IMPLICIT NONE
    REAL(SP),INTENT(INOUT),DEVICE :: PHI(Mloc,Nloc)

    !REAL(SP)  :: PHI_check(Mloc,Nloc)
    REAL(SP),DIMENSION(Mloc,Nghost) :: rNmsg, sNmsg,rSmsg,sSmsg
    REAL(SP),DIMENSION(Nghost,Nloc) :: rWmsg, sWmsg,rEmsg,sEmsg
    INTEGER :: status(MPI_STATUS_SIZE)
    INTEGER :: len

! for east-west

    len = Nloc * Nghost
    if ( n_west .ne. MPI_PROC_NULL ) then
       istat = cudaMemcpy2D(sWmsg(1,1),Nghost,PHI(Ibeg,1),Mloc,Nghost,Nloc)
       call MPI_SENDRECV(sWmsg ,len,MPI_SP,n_west, 1,& 
                       rWmsg,len,MPI_SP,n_west,0,comm2d,status,ier)
       istat = cudaMemcpy2D(PHI(Ibeg-Nghost,1),Mloc,rWmsg(1,1),Nghost,Nghost,Nloc)
    endif

    if ( n_east .ne. MPI_PROC_NULL ) then
        istat = cudaMemcpy2D(sEmsg(1,1),Nghost,PHI(Iend-Nghost+1,1),Mloc,Nghost,Nloc)
       call MPI_SENDRECV( sEmsg, len, MPI_SP,n_east, 0,&
                          rEmsg,len,MPI_SP,n_east,1,comm2d,status,ier )
       istat = cudaMemcpy2D(PHI(Iend+1,1),Mloc,rEmsg(1,1),Nghost,Nghost,Nloc)
    endif

! for nrth-suth

    len = Mloc * Nghost
    if ( n_suth .ne. MPI_PROC_NULL ) then
       istat = cudaMemcpy2D(sSmsg(1,1),Mloc,PHI(1,Jbeg),Mloc,Mloc,Nghost)
       call MPI_SENDRECV( sSmsg,len,MPI_SP, n_suth, 1,&
                          rSmsg,len,MPI_SP,n_suth, 0, comm2d,status, ier )
        istat = cudaMemcpy2D(PHI(1,Jbeg-Nghost),Mloc,rSmsg(1,1),Mloc,Mloc,Nghost)
        !PHI_check = PHI
        !print * , 'Send to S', PHI_check(900,4), sSmsg(900,1), 'Rec from S',&
        !PHI_check(900,1),rSmsg(900,1)
    endif

    if ( n_nrth .ne. MPI_PROC_NULL ) then
        istat = cudaMemcpy2D(sNmsg(1,1),Mloc,PHI(1,Jend-Nghost+1),Mloc,Mloc,Nghost)
       call MPI_SENDRECV( sNmsg,len,MPI_SP,n_nrth, 0,&
                          rNmsg,len,MPI_SP, n_nrth, 1, comm2d,status, ier )
         istat = cudaMemcpy2D(PHI(1,Jend+1),Mloc,rNmsg(1,1),Mloc,Mloc,Nghost)
        !PHI_check = PHI
        !print * , 'Rec from N', PHI_check(900,504), rNmsg(900,1), 'Send to N',&
        !PHI_check(900,501),sNmsg(900,1)
    endif

END SUBROUTINE phi_exch_cuda
!-------------------------------------------------------------------------------------
!
!    phi_exch_cuda is the subroutine to exchange GPU device variable at processor
!          interface
!    
!    HISTORY: 02/14/2011 Jeff Harris
!             05/01/2011 Fengyan Shi, implemented into the TVD code
!
!-------------------------------------------------------------------------------------

SUBROUTINE phi_int_exch_cuda (PHI)
    IMPLICIT NONE
    INTEGER,INTENT(INOUT),DEVICE :: PHI(Mloc,Nloc)

    INTEGER,DIMENSION(Mloc,Nghost) :: rNmsg, sNmsg,rSmsg,sSmsg
    INTEGER,DIMENSION(Nghost,Nloc) :: rWmsg, sWmsg,rEmsg,sEmsg
    INTEGER :: status(MPI_STATUS_SIZE)
    INTEGER :: len

! for east-west

    len = Nloc * Nghost
    if ( n_west .ne. MPI_PROC_NULL ) then
       istat = cudaMemcpy2D(sWmsg(1,1),Nghost,PHI(Ibeg,1),Mloc,Nghost,Nloc)
       call MPI_SENDRECV(sWmsg ,len,MPI_INTEGER,n_west, 1,& 
                       rWmsg,len,MPI_INTEGER,n_west,0,comm2d,status,ier)
       istat = cudaMemcpy2D(PHI(Ibeg-Nghost,1),Mloc,rWmsg(1,1),Nghost,Nghost,Nloc)
    endif

    if ( n_east .ne. MPI_PROC_NULL ) then
        istat = cudaMemcpy2D(sEmsg(1,1),Nghost,PHI(Iend-Nghost+1,1),Mloc,Nghost,Nloc)
       call MPI_SENDRECV( sEmsg, len, MPI_INTEGER,n_east, 0,&
                          rEmsg,len,MPI_INTEGER,n_east,1,comm2d,status,ier )
       istat = cudaMemcpy2D(PHI(Iend+1,1),Mloc,rEmsg(1,1),Nghost,Nghost,Nloc)
    endif

! for nrth-suth

    len = Mloc * Nghost
    if ( n_suth .ne. MPI_PROC_NULL ) then
       istat = cudaMemcpy2D(sSmsg(1,1),Mloc,PHI(1,Jbeg),Mloc,Mloc,Nghost)
       call MPI_SENDRECV( sSmsg,len,MPI_INTEGER, n_suth, 1,&
                          rSmsg,len,MPI_INTEGER,n_suth, 0, comm2d,status, ier )
        istat = cudaMemcpy2D(PHI(1,Jbeg-Nghost),Mloc,rSmsg(1,1),Mloc,Mloc,Nghost)
    endif

    if ( n_nrth .ne. MPI_PROC_NULL ) then
        istat = cudaMemcpy2D(sNmsg(1,1),Mloc,PHI(1,Jend-Nghost+1),Mloc,Mloc,Nghost)
       call MPI_SENDRECV( sNmsg,len,MPI_INTEGER,n_nrth, 0,&
                          rNmsg,len,MPI_INTEGER, n_nrth, 1, comm2d,status, ier )
         istat = cudaMemcpy2D(PHI(1,Jend+1),Mloc,rNmsg(1,1),Mloc,Mloc,Nghost)
    endif
END SUBROUTINE phi_int_exch_cuda

end module mgpu_utilities
# endif
