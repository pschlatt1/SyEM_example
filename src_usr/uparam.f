!=======================================================================
!     Adam Peplinski; 2015.10.20
!     Set of subroutines to read user and module parameters using 
!     namelists.
!     
!=======================================================================
!***********************************************************************
!     read parameters from the file
      subroutine uprm_read
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            ! NID
      include 'INPUT_DEF'
      include 'INPUT'           ! REAFLE

!     local variables
      integer len, ierr
      integer iunit

      character*132 fname 

!     functions
      integer ltrunc
!-----------------------------------------------------------------------
!     Open parameter file and read contents
      ierr=0
      if (NID.eq.0) then
!     find free unit
         call IO_freeid(iunit, ierr)
!     get file name and open file
         if(ierr.eq.0) then
            call blank(fname,132)
            len = ltrunc(REAFLE,132) - 4
            call chcopy(fname,REAFLE,len)
            fname(len+1:len+6)='.upar'
            write(6,*) 'Openning parameter file: ',trim(fname)
            open (unit=iunit,file=fname,status='old',action='read',
     $           iostat=ierr)
         endif
      endif
      call err_chk(ierr,'Error opening .upar file.$')

!     place to call module _param_in routines
      call uprm_in(iunit)

!     close the file
      ierr=0
      if (NID.eq.0) close(unit=iunit,iostat=ierr)
      call err_chk(ierr,'Error closing .upar file.$')

!     stamp logs
      if (NIO.eq.0) write(*,*) 'User parameter list'
      call uprm_out(6)
      if (NIO.eq.0) write(*,*) 

      return
      end
!***********************************************************************
!     read parameters
      subroutine uprm_in(iunit)
      implicit none

!     argument list
      integer iunit
!-----------------------------------------------------------------------
!     place to call module _param_in routines

      rewind(iunit)
!     user parameters
      call user_param_in(iunit)

      rewind(iunit)
!     Statistics
      call avg_param_in(iunit)

      rewind(iunit)
!     restart
      call chkpt_param_in(iunit)

      rewind(iunit)
!     SEM
      call sem_param_in(iunit)

      return
      end
!***********************************************************************
!     output parameters
      subroutine uprm_out(iunit)
      implicit none

!     argument list
      integer iunit
!-----------------------------------------------------------------------
!     place to call module _param_in routines

!     user parameters
      call user_param_out(iunit)

!     Statistics
      call avg_param_out(iunit)

!     restart
      call chkpt_param_out(iunit)

!     SEM
      call sem_param_out(iunit)

      return
      end
!***********************************************************************
!     read user parameters
      subroutine user_param_in(fid)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            !
      include 'PARALLEL_DEF' 
      include 'PARALLEL'        ! ISIZE, WDSIZE, LSIZE,CSIZE
      include 'USERPAR'         !

      real, parameter :: pi = 3.1415926535897932384626433832795028

!     argument list
      integer fid

!     local variables
      integer ierr

!     namelists; cannot be empty
      namelist /USERPAR/ bent_phi,bent_radius,z_offset
!-----------------------------------------------------------------------
!     default values
      bent_phi = 90.
      bent_radius = 3.1415/4
      z_offset = 0
!     read the file
      ierr=0
      if (NID.eq.0) then
         read(unit=fid,nml=USERPAR,iostat=ierr)
!     add exceptions in the future and check
!     if it is parameter and compiler independent
!     iostat
!     84 - NAMELIST group header not found in external file
!     85 - NAMELIST group header not found in internal file

      endif
      call err_chk(ierr,'Error reading USERPAR parameters.$')

      ! convert to rad
      bent_phi = bent_phi/180.*pi

!     broadcast data
      call bcast(bent_phi,WDSIZE)
      call bcast(bent_radius,WDSIZE)
      call bcast(z_offset,WDSIZE)

      return
      end
!***********************************************************************
!     write user parameters
      subroutine user_param_out(fid)
      implicit none

      include 'SIZE_DEF'
      include 'SIZE'            !
      include 'USERPAR'         !

!     argument list
      integer fid               ! file id

!     local variables
      integer ierr

!     namelists; cannot be empty
      namelist /USERPAR/ bent_phi,bent_radius,z_offset
!-----------------------------------------------------------------------
      ierr=0
      if (NID.eq.0) then
         write(unit=fid,nml=USERPAR,iostat=ierr)
      endif
      call err_chk(ierr,'Error writing USERPAR parameters.$')

      return
      end
!***********************************************************************
