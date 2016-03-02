!
!                              D   I   S   C   L   A   I   M   E   R
!
!       You are using a BETA version of the program dkes_data_prep2.f, which is currently
!       under development by D. A. Spong at the Fusion Energy Division,
!       Oak Ridge National Laboratory.  Please report any problems or comments
!       to him.  As a BETA version, this program is subject to change
!       and improvement without notice.
!
!
!      This code reads the visc.ext files produced by dkes_data_prep1.f
!       and processes their data to produce
!       tables of transport and viscosity coefficients
!       as a function of collisionality (visc.ext files). The
!       viscosities are calculated using the analysis of
!       Sugama, et al., Phys. Plasmas vol. 9 (2002), pg. 4637]
!
!       Command line arguments:
!        xprep2 temp_file.ext surf_inner surf_outer delta_surf
!         where temp_file.ext is just the filename extension of the temp.*
!         file produced by xrun that you want to process.
!
!        surf_inner = innermost surface used
!        delta_surf = increment in surface number
!        surf_outer = outermost surface used
!        note: these should be the same surfaces as previously run with dkes_data_prep1.f
!
        implicit none
        INTEGER, PARAMETER :: mxpts=200
	REAL*8, DIMENSION(mxpts):: cmul, efield, L11, L31, L33,
     >    mstar, lstar, nstar, mpp, mpt, mtt
        REAL*8, DIMENSION(:), ALLOCATABLE :: efield_all
        REAL*8, DIMENSION(:,:), ALLOCATABLE :: L11_all, L31_all,
     >    L33_all, mstar_all, lstar_all, nstar_all, mpp_all,
     >    mpt_all, mtt_all
        REAL*8 dsdr, wtovth
        CHARACTER*60, DIMENSION(15) :: visc_files
        CHARACTER*60 mstar_file, nstar_file, lstar_file,
     >              mpp_file, mpt_file, mtt_file,
     >              l11_file, l31_file, l33_file, file_nm
	CHARACTER*1 tb
	CHARACTER*20 :: arg1, arg2, arg3, arg4
	CHARACTER*20 :: dr, srf
	CHARACTER*10 :: dr_num
        CHARACTER*60 :: file_out
	INTEGER i, iunit, num_files, j, ipts, ierr, jj, j_inv,
     >   isurf, surf0, outer_surf, delta_surf
	dsdr = 1.      !default value
	wtovth = 0.2   !default value
	num_files = 15
        CALL getarg(1, arg1)
        CALL getarg(2, arg2)
        CALL getarg(3, arg3)
        CALL getarg(4, arg4)
	read(arg2,'(i4)') surf0
	read(arg3,'(i4)') outer_surf
	read(arg4,'(i4)') delta_surf
c
c      Start loop over flux surfaces:
c
       do isurf=surf0,outer_surf,delta_surf
       write(dr_num, '(i4)') isurf
       dr = "s" // trim(adjustl(dr_num)) // "/"
       srf = "s" // trim(adjustl(dr_num))
       file_nm = "visc." // trim(adjustl(arg1)) // "_" //
     >   trim(adjustl(srf)) // "_E"
c
c      Identify and open the visc files that are to be
c       combined for the data input file for lijs:
c
	visc_files(1)= trim(adjustl(dr)) //
     >    trim(adjustl(file_nm))//"0.00000"
	visc_files(2)= trim(adjustl(dr)) //
     >    trim(adjustl(file_nm))//"1m10"
	visc_files(3)= trim(adjustl(dr)) //
     >    trim(adjustl(file_nm))//"1m9"
	visc_files(4)= trim(adjustl(dr)) //
     >    trim(adjustl(file_nm))//"1m8"
	visc_files(5)= trim(adjustl(dr)) //
     >    trim(adjustl(file_nm))//"1m7"
	visc_files(6)= trim(adjustl(dr)) //
     >    trim(adjustl(file_nm))//"1m6"
	visc_files(7)= trim(adjustl(dr)) //
     >    trim(adjustl(file_nm))//"1m5"
	visc_files(8)= trim(adjustl(dr)) //
     >    trim(adjustl(file_nm))//"0.00010"
	visc_files(9)= trim(adjustl(dr)) //
     >    trim(adjustl(file_nm))//"0.00030"
	visc_files(10)= trim(adjustl(dr)) //
     >    trim(adjustl(file_nm))//"0.00100"
	visc_files(11)= trim(adjustl(dr)) //
     >    trim(adjustl(file_nm))//"0.00300"
	visc_files(12)= trim(adjustl(dr)) //
     >    trim(adjustl(file_nm))//"0.01000"
	visc_files(13)= trim(adjustl(dr)) //
     >    trim(adjustl(file_nm))//"0.03000"
	visc_files(14)= trim(adjustl(dr)) //
     >    trim(adjustl(file_nm))//"0.10000"
	visc_files(15)= trim(adjustl(dr)) //
     >    trim(adjustl(file_nm))//"0.30000"
c	write(*,*) visc_files(1)
c	write(*,*) visc_files(14)
	do i=1,num_files
	  iunit = i + 6
c	  write(*,*) i, visc_files(i)
	  open (unit=iunit,file=visc_files(i),status='old')
	end do
c
c     Initialize arrays:
c
c      do i=1,mxpts
c       cmul(i) = 0.; efield(i) = 0.; L11(i) = 0.; L31(i) = 0.
c       L33(i) = 0.; mstar(i) = 0.; lstar(i) = 0.; nstar(i) = 0.
c       mpp(i) = 0.; mpt(i) = 0.; mtt(i) = 0.
c      end do
c
c      Read in one of the visc files and count the number of
c       collisionalities that are present (it is then assumed
c       that all subsequent visc files have the same number of
c       collisionalities.
c
	iunit = 7
        read(iunit,98) tb	
        read(iunit,98) tb
	do j=1,mxpts	
         READ(iunit,97, END=100) cmul(j),tb,efield(j),tb,mpp(j),tb,
     >    mpt(j),tb,mtt(j),tb,mstar(j),tb,
     >    lstar(j),tb,nstar(j),tb,L11(j),tb,L31(j),tb,L33(j),tb,wtovth
         ipts = j
	end do
  100   continue
        rewind(unit=7)
c	write(*,*) ipts
c
c      Allocate arrays for the electric field (1D), and transport
c       coefficients (2D: collisionality and efield).  Open each of
c       the visc files and place the transport coefficient data into
c       the 2D arrays:
c
        allocate(efield_all(num_files), stat=ierr)	
	allocate(L11_all(ipts,num_files),stat=ierr)
	allocate(L31_all(ipts,num_files),stat=ierr)
	allocate(L33_all(ipts,num_files),stat=ierr)
	allocate(mpp_all(ipts,num_files),stat=ierr)
	allocate(mpt_all(ipts,num_files),stat=ierr)
	allocate(mtt_all(ipts,num_files),stat=ierr)
	allocate(mstar_all(ipts,num_files),stat=ierr)
	allocate(nstar_all(ipts,num_files),stat=ierr)
	allocate(lstar_all(ipts,num_files),stat=ierr)
c
c      Initialize arrays to zero:
c
c        do i=1,num_files
c	 efield_all(i) = 0.
c	 do j=1,ipts
c	  L11_all(j,i) = 0.; L31_all(j,i) = 0.; L33_all(j,i) = 0.
c	  mpp_all(j,i) = 0.; mpt_all(j,i) = 0.; mtt_all(j,i) = 0.
c	  mstar_all(j,i) = 0.; nstar_all(j,i) = 0.; lstar_all(j,i) = 0.
c	 end do
c	end do
	do i=1,num_files
	 iunit = i + 6
         read(iunit,98) tb	
         read(iunit,98) tb
	 do j=1,ipts	
          READ(iunit,97) cmul(j),tb,efield(j),tb,mpp(j),tb,
     >    mpt(j),tb,mtt(j),tb,mstar(j),tb,
     >    lstar(j),tb,nstar(j),tb,L11(j),tb,L31(j),tb,L33(j),tb,wtovth
          efield_all(i) = efield(1)
          L11_all(j,i) = L11(j)
          L31_all(j,i) = L31(j)
          L33_all(j,i) = L33(j)
          mstar_all(j,i) = mstar(j)
          nstar_all(j,i) = nstar(j)
          lstar_all(j,i) = lstar(j)
          mpp_all(j,i) = mpp(j)
          mpt_all(j,i) = mpt(j)
          mtt_all(j,i) = mtt(j)
	 end do
c         write(*,*) L11(1),L31(1),L33(1)
	end do
c
c      Write out the transport coefficient data into the file format
c        needed by lijs (separate files used for each transport coefficient):
c
        mstar_file = trim(adjustl(dr)) // "mstar_lijs_" //
     >    trim(adjustl(arg1)) // "_s" // trim(adjustl(dr_num))
        nstar_file = trim(adjustl(dr)) // "nstar_lijs_" //
     >    trim(adjustl(arg1)) // "_s" // trim(adjustl(dr_num))
        lstar_file = trim(adjustl(dr)) // "lstar_lijs_" //
     >    trim(adjustl(arg1)) // "_s" // trim(adjustl(dr_num))
        mpp_file = trim(adjustl(dr)) // "mpp_lijs_" //
     >    trim(adjustl(arg1)) // "_s" // trim(adjustl(dr_num))
        mpt_file = trim(adjustl(dr)) // "mpt_lijs_" //
     >    trim(adjustl(arg1)) // "_s" // trim(adjustl(dr_num))
        mtt_file = trim(adjustl(dr)) // "mtt_lijs_" //
     >    trim(adjustl(arg1)) // "_s" // trim(adjustl(dr_num))
        l11_file = trim(adjustl(dr)) // "l11_lijs_" //
     >    trim(adjustl(arg1)) // "_s" // trim(adjustl(dr_num))
        l31_file = trim(adjustl(dr)) // "l31_lijs_" //
     >    trim(adjustl(arg1)) // "_s" // trim(adjustl(dr_num))
        l33_file = trim(adjustl(dr)) // "l33_lijs_" //
     >    trim(adjustl(arg1)) // "_s" // trim(adjustl(dr_num))
	open(unit=30, file=mstar_file, status="unknown")
	open(unit=31, file=nstar_file, status="unknown")
	open(unit=32, file=lstar_file, status="unknown")
	open(unit=33, file=mpp_file, status="unknown")
	open(unit=34, file=mpt_file, status="unknown")
	open(unit=35, file=mtt_file, status="unknown")
	open(unit=36, file=l11_file, status="unknown")
	open(unit=37, file=l31_file, status="unknown")
	open(unit=38, file=l33_file, status="unknown")
	do i=30,38
	  write(i,*) dsdr, wtovth
	  write(i,*) ipts, num_files
	  do j=1,ipts
	    j_inv = ipts - j + 1      !must reverse cmul order to be increasing
	    write(i,'(1x,e15.7)') cmul(j_inv)
	  end do
	  do j=1,num_files
	    write(i,'(1x,e15.7)') efield_all(j)
	  end do
	  do j=1,num_files
	   do jj=1,ipts
	     j_inv = ipts - jj + 1      !must reverse cmul order to be increasing
	     if(i .eq. 30) write(i,'(1x,e15.7)') mstar_all(j_inv,j)
	     if(i .eq. 31) write(i,'(1x,e15.7)') nstar_all(j_inv,j)
	     if(i .eq. 32) write(i,'(1x,e15.7)') lstar_all(j_inv,j)
	     if(i .eq. 33) write(i,'(1x,e15.7)') mpp_all(j_inv,j)
	     if(i .eq. 34) write(i,'(1x,e15.7)') mpt_all(j_inv,j)
	     if(i .eq. 35) write(i,'(1x,e15.7)') mtt_all(j_inv,j)
	     if(i .eq. 36) write(i,'(1x,e15.7)') l11_all(j_inv,j)
	     if(i .eq. 37) write(i,'(1x,e15.7)') l31_all(j_inv,j)
	     if(i .eq. 38) write(i,'(1x,e15.7)') l33_all(j_inv,j)
	   end do
	  end do
	end do
c
c      Close the visc files and terminate program:
c
	do i=1,num_files
	 iunit = i + 6
	 close(unit=iunit)
	end do
	do i=30,38
	  close(unit=i)
	end do
        deallocate(efield_all)	
	deallocate(L11_all)
	deallocate(L31_all)
	deallocate(L33_all)
	deallocate(mpp_all)
	deallocate(mpt_all)
	deallocate(mtt_all)
	deallocate(mstar_all)
	deallocate(nstar_all)
	deallocate(lstar_all)
	end do      !end of loop over surfaces: isurf=surf0,outer_surf,delta_surf
   98   FORMAT(a1)
   97   FORMAT(12(e12.5,a1),e12.5)
	end
