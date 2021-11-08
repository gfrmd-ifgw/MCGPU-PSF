

!* Original penEasy disclaimer
!*********************************************************************************************
!*penEasy copyright and disclaimer
!*penEasy
!*Copyright (c) 2003-2019
!*Universitat Politecnica de Catalunya
!*Permission to use, copy, modify and re-distribute copies of this software, or parts of it,
!*and its documentation for any purpose is hereby granted without fee, provided that the
!*above copyright notice appears in all copies and that both that copyright notice and this
!*permission notice appear in all supporting documentation. The Universitat Politecnica
!*de Catalunya makes no representations about the suitability of this software for any
!*purpose. It is provided \as is" without express or implied warranty.
!**********************************************************************************************




!*******************************************************************
!*                      PHASE SPACE FILE                           *
!*                      Converter                                  *
!*                                                                 *
!* Short description:                                              *
!*   Writes to a file the dynamical state of all particles that    *
!*   reach a given detection material in IAEA format.              *
!* Rev:                                                            *
!*   2019-12-23                                                    *
!*******************************************************************



!*******************************************************************
!*******************************************************************
!*    Modules section                                              *
!*******************************************************************
!*******************************************************************

!*Dummy PENELOPE variables for compatibility

C  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      MODULE TRACK_mod

      SAVE  
      DOUBLE PRECISION :: E,X,Y,Z,U,V,W,WGHT
      INTEGER*4 :: KPAR,IBODY,MAT
      INTEGER*4, DIMENSION (5) :: ILB
      INTEGER*4 :: IPOL=0
      DOUBLE PRECISION :: SP1,SP2,SP3
      LOGICAL :: LAGE =.FALSE.
      DOUBLE PRECISION :: PAGE=0.0D0

      INTEGER*4 :: MHINGE
      DOUBLE PRECISION :: E0SEGM,DESOFT,SSOFT
C
      END MODULE TRACK_mod
C  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

C  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      MODULE PENELOPE_mod
      SAVE  
      INTEGER*4, PARAMETER :: MAXMAT=10
      DOUBLE PRECISION, DIMENSION (3,MAXMAT) :: EABS=50.0D0
      DOUBLE PRECISION, DIMENSION (MAXMAT) :: C1=0.01D0,C2=0.01D0
      DOUBLE PRECISION, DIMENSION (MAXMAT) :: WCC=1.0D2,WCR=1.0D2
      INTEGER*4, PARAMETER :: NMS=1000
      INTEGER*4, PARAMETER :: NEGP=200
      INTEGER*4 :: NMAT
      DOUBLE PRECISION, DIMENSION (MAXMAT) :: DEN=1.0D0, RDEN=1.0D0
      END MODULE PENELOPE_mod
C  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      module RSEEDcommonMod
!*******************************************************************
      implicit none
      save
      integer*4 seed1,seed2
      common/rseed/seed1,seed2
      end
	  
!*End of Dummy PENELOPE variables for compatibility
      module PSFmod
!*******************************************************************
!*    Vars for the PSF tally.                                      *
!*******************************************************************
      implicit none
      save
      integer detmat,psfunit,formatiaea
      integer*8 nhisti,nhlast,npar(4)
	  integer*4 dn, pcount
      end


      module PSFendinoMod
!*******************************************************************
!*    Vars for the endino particle---used for finishing a PSF, see *
!*    tally documentation.                                         *
!*******************************************************************
      implicit none
      save
      integer*4 gkpar,gilb(5)
      real*8 gen,gx,gy,gz,gu,gv,gw
      end


      module PSFiaeaMod
!*******************************************************************
!*    Vars for IAEA-formatted PSFs                                 *
!*******************************************************************
      implicit none
      save
      integer*4 sourceWrite
      end
	  

!*******************************************************************
!*******************************************************************
!*    Routines section                                             *
!*******************************************************************
!*******************************************************************


!*******************************************************************
!*    IAEA specific routines follow.                               *
!*    Not adapted to protons: IAEAwrite will issue an error message*
!*    and abort if called.                                         *
!*******************************************************************
#ifdef IAEAPSF

      subroutine IAEAiniwrite(buffer)
!*******************************************************************
!*    Initializes.                                                 *
!*                                                                 *
!*    Input:                                                       *
!*      buffer -> PSF filename                                     *
!*******************************************************************
      use PSFiaeaMod
      implicit none
      character*80 buffer

      integer i
      integer*4 accessWrite,accessAppend,nXtraReal,nXtraInt
      integer*4 ind,xtraIntType,errorcode
      parameter (accessWrite=2,accessAppend=3)
!      parameter (nXtraReal=0,nXtraInt=6)
	  parameter (nXtraReal=0,nXtraInt=0)

      i = len_trim(buffer)+1
      buffer(i:i) = char(0)       ! C convention for end-of-string

        call iaea_new_source(sourceWrite,buffer,accessWrite,
     &                       errorcode)
        if (errorcode.lt.0) then
          write(*,'(a)') 'IAEAiniwrite:ERROR: Unable to open PSF.'
          stop
        endif


      call iaea_set_extra_numbers(sourceWrite,nXtraReal,nXtraInt)
      errorcode = 0    ! Set to no-error

      ind = 0          ! First extra long-int number
      xtraIntType = 1  ! This IAEA code identifies DeltaN
      call iaea_set_type_extralong_variable(sourceWrite,ind,
     &  xtraIntType)
      errorcode = min(errorcode,xtraIntType)

      xtraIntType = 0  ! IAEA code identifies user-defined extra int, default if IAEAPSFold is active
      ind = 1
#ifndef IAEAPSFold
      ! Using IAEA routines that recognize ILB(1:4) codes
      xtraIntType = 7  ! This IAEA code identifies ilb(1)
#endif
      call iaea_set_type_extralong_variable(sourceWrite,ind,
     &  xtraIntType)
      errorcode = min(errorcode,xtraIntType)

      xtraIntType = 0
      ind = 2
#ifndef IAEAPSFold
      xtraIntType = 6  ! This IAEA code identifies ilb(2)
#endif
      call iaea_set_type_extralong_variable(sourceWrite,ind,
     &  xtraIntType)
      errorcode = min(errorcode,xtraIntType)

      xtraIntType = 0
      ind = 3
#ifndef IAEAPSFold
      xtraIntType = 5  ! This IAEA code identifies ilb(3)
#endif
      call iaea_set_type_extralong_variable(sourceWrite,ind,
     &  xtraIntType)
      errorcode = min(errorcode,xtraIntType)

      xtraIntType = 0
      ind = 4
#ifndef IAEAPSFold
      xtraIntType = 4  ! This IAEA code identifies ilb(4)
#endif
      call iaea_set_type_extralong_variable(sourceWrite,ind,
     &  xtraIntType)
      errorcode = min(errorcode,xtraIntType)

      ind = 5
      xtraIntType = 3  ! This IAEA code identifies ilb(5), already defined in old IAEA codes
      call iaea_set_type_extralong_variable(sourceWrite,ind,
     &  xtraIntType)
      errorcode = min(errorcode,xtraIntType)

      ! Note: xtraIntType=2 identifies LATCH variable in EGS codes

      select case (errorcode)
      case (-1)
        write(*,'(a)') 'IAEAiniwrite:ERROR: '//
     &    'Unable to set an extra variable.'
        stop
      case (-2)
        write(*,'(a)') 'IAEAiniwrite:internalERROR: '//
     &    'Extra variable index is out of range.'
        stop
      case (-3)
        write(*,'(a)') 'IAEAiniwrite:internalERROR: '//
     &    'Extra variable type is out of range.'
        stop
      case default
        if (errorcode.lt.0) then
          write(*,'(a)') 'IAEAiniwrite:internalERROR: '//
     &    'Undefined error while setting an extra variable.'
          stop
        endif
      end select
      end


      subroutine IAEAwrite(dn)
!*******************************************************************
!*    Write a particle to the PSF in IAEA format.                  *
!*                                                                 *
!*    Input:                                                       *
!*      dn -> Incremental history num.                             *
!*******************************************************************
      use TRACK_mod
      use PSFiaeaMod
      implicit none
      integer*4 dn

      integer*4 n_stat,kpariaea,xtraInt(6)
      real*4 xtraFloat(1)
      real*4 eiaea,xiaea,yiaea,ziaea,uiaea,viaea,wiaea,wghtiaea

      ! Convert particle index from PENELOPE to IAEA:
      select case (kpar)
      case (1)
        kpariaea = 2
      case (2)
        kpariaea = 1
      case (3)
        kpariaea = 3
      case (4)
        write(*,'(a)') 'IAEAwrite:ERROR: protons not implemented.'
        stop
      case default
        write(*,'(a)') 'IAEAwrite:ERROR: Invalid KPAR.'
        stop
      end select

      eiaea = real(e)*1.0e-6  ! Energy in MeV
      xiaea = real(x)         ! Vars type should match C++ definition for 'float'
      yiaea = real(y)
      ziaea = real(z)
      uiaea = real(u)
      viaea = real(v)
      wiaea = real(w)
      wghtiaea = real(wght)

      xtraInt(1) = dn
      xtraInt(2) = ilb(1)
      xtraInt(3) = ilb(2)
      xtraInt(4) = ilb(3)
      xtraInt(5) = ilb(4)
      xtraInt(6) = ilb(5)
      xtraFloat(1) = 0.0      ! Dummy value, not written to PSF

      n_stat = dn             ! Should be 0 if same history, >0 else

      call iaea_write_particle(sourceWrite,n_stat,kpariaea,eiaea,
     &  wghtiaea,xiaea,yiaea,ziaea,uiaea,viaea,wiaea,xtraFloat,
     &  xtraInt)

      if (n_stat.eq.-1) then
        write(*,'(a)') 'IAEAwrite:ERROR: Unable to write particle.'
        stop
      endif
      end


      subroutine IAEAwriteReport(mode)
!*******************************************************************
!*    Uptates IAEA header file.                                    *
!*                                                                 *
!*    Input:                                                       *
!*      mode -> -1 if called at end-of-simulation, 0 for dump only,*
!*               1 for tally report only, 2 for both.              *
!*******************************************************************
      use PSFmod
      use PSFiaeaMod
      implicit none
      integer mode

      integer*4 errorcode
      integer(selected_int_kind(R=18)) :: nTopiaea

      call flush(6)  ! Prevents messy output when using Fortran+C

      nTopiaea = nhisti
      call iaea_set_total_original_particles(sourceWrite,nTopiaea)
      if (nTopiaea.lt.0) then
        write(*,'(a)') 'IAEAwriteReport:ERROR: '//
     &    'unable to write no. of histories to header file.'
        stop
      endif

      call iaea_update_header(sourceWrite,errorcode)
      if (errorcode.lt.0) then
        write(*,'(a)') 'IAEAwriteReport:ERROR: '//
     &    'unable to update header file.'
        stop
      endif

      if (mode.eq.-1) then  ! Close the PSF
        call iaea_destroy_source(sourceWrite,errorcode)
        if (errorcode.lt.0) then
          write(*,'(a)') 'IAEAwriteReport:ERROR: '//
     &      'unable to close PSF.'
          ! Since already at end-of-sim, do not stop here
        endif
      endif
      end

#endif
!*******************************************************************
!*    End of IAEA block.                                           *
!*******************************************************************

      integer function finduf()
!*******************************************************************
!*    Finds a valid (i.e. unused) file unit.                       *
!*******************************************************************
      implicit none
      logical used
      integer maxuf
      parameter (maxuf=17)

      finduf = 6
      do
        finduf = finduf+1
        if (finduf.gt.maxuf) then
          write(*,*) ''
          write(*,'(a)') 'finduf:ERROR: Unable to find a valid unit'
          stop
        endif
        inquire(finduf,opened=used)
        if (.not.used) return
      enddo
      end


      subroutine getline(buffer,ufile, errorcode)
!*******************************************************************
!*    Reads a line of text from stdin or file. The line is returned*
!*    only if it is not blank nor a comment line.                  *
!*                                                                 *
!*    Input:                                                       *
!*      ufile: unit file, must be already open; if =0, reads stdin *
!*    Output:                                                      *
!*      buffer: line read, 80 chars max                            *
!*    Comments:                                                    *
!*      Lines starting with a '#' are taken as comments and skipped*
!*******************************************************************
      implicit none
      character*80 buffer
      integer ufile

      integer errorcode

      do                                           ! Read while empty or comment lines
        if (ufile.eq.0) then                       ! Read from std input
          read(*,'(a80)',iostat=errorcode) buffer
        else                                       ! Read from unit file
          read(ufile,'(a80)',iostat=errorcode) buffer
        endif
        if (errorcode.ne.0) then
          !write(*,'(a)') 'getline:ERROR: unable to read line.'
          EXIT
        endif
        if (buffer(1:1).eq.'#') cycle              ! A comment line
        if (verify(buffer,' ').ne.0) return        ! Non-blank line, done
      enddo                                        ! Blank line, read another
      end
	  
      function upcase(string) result(upper)
!*******************************************************************
!*    Converts a string to uppercase.                              *
!*    Adapted from http://www.star.le.ac.uk/~cgp/fortran.html      *
!*    No explicit interface needed for this version, 80 characters *
!*    assumed.                                                     *
!*******************************************************************
      implicit none
      character*80 string,upper
      integer j

      do j=1,len(string)
        if(string(j:j) >= "a" .and. string(j:j) <= "z") then
             upper(j:j) = achar(iachar(string(j:j))-32)
        else
             upper(j:j) = string(j:j)
        end if
      enddo
      end
	  
!******************************************************************
!         MAIN PROGRAM 
      program main
	  
	  use TRACK_mod
	  use PSFmod
	  
	  character*80 fname, line, foutput
	  integer uf, error, finduf, endoffile, nparticle
	  real*8 emax
	  character*80 fmtstr
      parameter (fmtstr='((es12.5),6(1x,es12.5))')
!	  real*4 psfmatrix(10019,7) 
      real*4, allocatable :: psfmatrix(:,:)
	  integer totpar, j
	  
	  write (*,'(A)',advance='no') 'Enter PSF filename'
	  read (*,*) fname
	  uf = finduf()                   ! Find a valid unit for the file
      !open(uf,file=fname,status='old', iostat=error)
	  
	  write (*,'(A)',advance='no') 'Enter maximum energy'
	  read (*,*) emax
	  
	  write (*,'(A)',advance='no') 'Enter output file name'
	  read (*,*) foutput
	  
	  write (*,'(A)',advance='no') 'Enter number of particles'
	  read (*,*) totpar
	  
	  write (*,'(A)',advance='no') 'Enter total of histories'
	  read (*,*) nhisti
	  
	  allocate(psfmatrix(totpar, 7))
	  
	  
	  call IAEAiniwrite(foutput)
	  dn = 0
	  pcount = 0
	  
	  open(uf, file=fname, form='unformatted')
	  read(uf) psfmatrix
	  close(uf)
	  
	  do j=1, totpar
!	  call getline(line, uf, endoffile)
!	  if (endoffile.eq.0) then 
	       dn = 1
		   pcount = pcount + 1
	       !set particle state!
		   X = psfmatrix(j, 1)
		   Y = psfmatrix(j, 2)
		   Z = psfmatrix(j, 3)
		   U = psfmatrix(j, 4)
		   V = psfmatrix(j, 5)
		   W = psfmatrix(j, 6)
		   E = psfmatrix(j, 7)
!		   read(line, *, iostat=endoffile) X,Y,Z,U,V,W,E
! DEBUG		   
!		   write(*,*) E,U,V,W,X,Y,Z
		   WGHT = 1
		   KPAR = 2
		   ilb(1) = 1
		   ilb(2) = 0
		   ilb(3) = 0
		   ilb(4) = 0
		   ilb(5) = 0
		   call IAEAwrite(dn)
!	    else 
!	       write(*,*) 'done'
!		   EXIT
!	   endif 
	  enddo
	  
#ifdef IAEAPSF
      call IAEAwriteReport(-1)
#endif
	  end
!******************************************************************


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
