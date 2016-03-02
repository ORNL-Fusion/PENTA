!  Compile with (Carver, 6/2011):
!  mpif90 -o xrun run_dkes_multi_proc_ver4.f
!  Run on Carver (6/2011) with:
!  mpirun -np 16 ./xrun 280 350 t 1 16 > cms
!  Or submit to batch system with:
!   #PBS -q regular
!   #PBS -A mp202
!   #PBS -l nodes=2:ppn=8
!   #PBS -l walltime=160:00:00
!   #PBS -N dkes
!   #PBS -e my_job.$PBS_JOBID.err
!   #PBS -o my_job.$PBS_JOBID.out
!   #PBS -V
!   #PBS -mbe
!
!   cd $PBS_O_WORKDIR
!   mpirun -np 16 ./xrun 280 350 t 1 16 > cms
!
!    Command line arguments to xrun:
!     1st and 2nd args - these are the starting and ending flux
!      surfaces to be analyzed.  Since the parallelization is
!      over flux surface, the number of flux surfaces should
!      equal the number of processors available
!      (i.e., in this case
!        last surface - first surface + 1 = 17 - 2 + 1 = 16)
!     3rd arg - This is the dir_create logical variable.  This
!      should be t if it is the first run and f if it is a subsequnet
!      run.  For the first run, subdirectories s2, s3, s4, ..., s17
!      are set up to keep the files for the different processors
!      separate.  Also, links are set up in each subdirectory to the
!      two required equilibrium input files: boozmn_ncsx.nc and wout_ncsx.nc
!     4th and 5th args - These determine the range of electric fields
!      (out of the efields array) that will be cycled through.  Time
!      constraints may require only one efield to be done in a run.  In
!      that case, just make these two args the same number.
!
      implicit none
      INCLUDE 'mpif.h'
!     This is a shell program which allows multiple copies
!      of DKES to be run in parallel
      INTEGER, PARAMETER :: itest = 0
      CHARACTER*50, PARAMETER :: boozer_ext = 'rwm'
      CHARACTER*10, PARAMETER :: dkes_exec = 'xdkes'
      CHARACTER*10, PARAMETER :: dkes_prep = 'xdkes_prep'
      CHARACTER*50 :: booz_file, wout_file
      CHARACTER*50 :: arg1, arg2, arg3, arg4, arg5
      CHARACTER*45 :: extension
      CHARACTER*12 :: ch_cmul
      CHARACTER*12 :: line1_end
      CHARACTER*3 :: ch_npe
      CHARACTER*4 :: ch_surf, ch_i
      CHARACTER*80 :: cmd_line1, cmd_line2, cmd_line3, cmd_line4,
     1   cmd_line5, cmd_line6, log_file, cmd, cmd_line0
      INTEGER, PARAMETER :: num_surfs = 40, num_cmul = 32
      INTEGER, PARAMETER :: num_efields = 15
      REAL, DIMENSION(num_cmul):: cmul
      CHARACTER*12, DIMENSION(num_efields) :: ch_efield
      CHARACTER*20 :: ch_num_efield
      REAL, DIMENSION(num_efields):: efields
      INTEGER i, isurf, icount
      INTEGER ier, npes, mype, ie
      INTEGER :: js, inner, outer, ie_upper, ie_lower
      LOGICAL :: dir_create
      CHARACTER*20, DIMENSION(800) :: dirname, dnum
      CHARACTER*20 :: my_dir, my_dirname
c
c   Tabulate the collisionality values and electric fields over
c    which DKES will be run.
c
      cmul(1:num_cmul) = (/40.0,20.0, 10.0, 6.2, 3.84, 0.238E+01,
     2  0.148E+01, 0.916E+00, 0.568E+00, 0.352E+00,
     3  0.218E+00, 0.135E+00, 0.839E-01, 0.520E-01, 0.323E-01,
     4  0.200E-01, 0.124E-01, 0.769E-02, 0.477E-02, 0.296E-02,
     5  0.183E-02, 0.114E-02, 0.704E-03, 0.437E-03, 0.271E-03,
     6  0.168E-03, 0.104E-03, 0.704E-04, 0.437E-04, 0.271E-04,
     6  0.168E-04, 0.10E-04/)
c      cmul(1:num_cmul) = (/1280.0,640.0,320.0,160.0,80.0,40.0,20.0/)
      efields(1:num_efields) = (/0.0000, 0.0001,0.0003,0.001,0.003,
     1  0.01,0.03,0.1,0.3,
     2  1.e-10,1.e-9,1.e-8,1.e-7,1.e-6,1.e-5/)
       ch_efield((num_efields-5):num_efields) = (/"  1m10 ",
     1  "  1m9  ","  1m8  ","  1m7  ","  1m6  ","  1m5  "/)
       do i=1,num_efields-6
         WRITE(ch_efield(i),'(f7.5)') efields(i)
       end do
c
c    Get the command line arguments.  arg1 is the starting flux surface
c     value (integer), arg2 is the ending flux surface value (integer),
c     arg3 is a logical that the user should use to indicate if this is
c     the first run, implying that directories need to be created by xrun
c     to contain the parallel job output (arg3=.true.), or if this is a
c     continuation run and directories have already been created and
c     don't need to be recreated (arg3=.false.), arg4 = index of lower
c     value of electric field, arg5 = index of upper value of electric field.
c
      call getarg(1, arg1)
      call getarg(2, arg2)
      call getarg(3, arg3)
      call getarg(4, arg4)
      call getarg(5, arg5)
      read(arg1, '(i3)') inner
      read(arg2, '(i3)') outer
      read(arg3, '(L1)') dir_create
      read(arg4, '(i2)') ie_lower
      read(arg5, '(i2)') ie_upper
c
c    Form the names of the two equilibrium input files:
c
      wout_file = "wout_" // trim(adjustl(boozer_ext)) // ".nc"
      booz_file = "boozmn_" // trim(adjustl(boozer_ext)) // ".nc"
c
! Get NPES and MYPE.  Requires initialization of MPI.
      CALL MPI_INIT(IER)
      IF (IER .NE. 0) THEN
        WRITE(*,'("MPI_INIT returned IER =")') IER
        STOP
        ENDIF
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPES, IER)
      IF (IER .NE. 0) THEN
        WRITE(*,'("MPI_COMM_SIZE returned IER =")') IER
        STOP
        ENDIF
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYPE, IER)
      IF (IER .NE. 0) THEN
        WRITE(*,'("MPI_COMM_RANK returned IER =")') IER
        STOP
        ENDIF
        
      if(mype .eq. 0) write(*,'("Running with",i3," processors")') npes
      if(mype .eq. 0) write(*,*) wout_file, booz_file
c
c    Set up subdirectories to contain output from parallel jobs.
c     Also, set up links to the wout and boozmn files. 
       i=4*mype+inner
       write(dnum(i),'(i4)') (4*mype+inner)
       dirname(i) = "s" // trim(adjustl(dnum(i)))
       cmd = "mkdir " // dirname(i)
       if(dir_create) write(*,*) trim(adjustl(cmd))
       if(dir_create) call system(trim(adjustl(cmd)))
       cmd_line0 = "cd " // trim(adjustl(dirname(i))) //
     1   "; ln -s ../" //trim(adjustl(wout_file)) // " " //
     2   trim(adjustl(wout_file))
        if(itest .eq. 0 .and. dir_create)
     1    call system(trim(adjustl(cmd_line0)))
       cmd_line0 =  "cd " // trim(adjustl(dirname(i))) //
     1  ";ln -s ../" //trim(adjustl(booz_file)) // " " //
     2  trim(adjustl(booz_file))
        if(itest .eq. 0 .and. dir_create)
     1     call system(trim(adjustl(cmd_line0)))
c
c    Assign different flux surfaces and subdirectories to each
c    parallel job (mype)
c      
      icount = 1
       isurf = 4*mype + inner
       write(my_dir,'(i4)') isurf
       my_dirname = "s" // trim(adjustl(my_dir))
       log_file = trim(adjustl(my_dirname)) // "/log_mpp_cmds" //
     1   trim(adjustl(my_dir))
       if(itest .eq. 1)write(*,*) log_file
       if(itest .eq. 1)
     1    open(unit=10,file=trim(adjustl(log_file)),status="unknown")
       if(itest .eq. 1) write(10,*) cmd_line0
       WRITE(ch_npe,'(i3)') mype
       if(isurf .gt. 1) WRITE(ch_surf,'(i4)') isurf
c        Exclude js = 1 VMEC surface since there seem to be problems
c           running DKES there:
       if(isurf .eq. 1) WRITE(ch_surf,'(i4)') isurf + 1
c
c      Start loop over collisionlaties and electric fields for each
c       processor.  Resolution parameters are defined for DKES depending
c       on collisionality based on previous experience.
c
        DO ie=ie_lower, ie_upper
         WRITE(ch_num_efield,'(e18.12)') efields(ie)
         DO i=1,num_cmul
          WRITE(ch_i,'(i4)') (i + 12)
	  
          if(cmul(i) .ge. 10. .and. cmul(i) .lt. 10000.) then
           WRITE(ch_cmul,'(f7.1)') cmul(i)
           line1_end = "f 3 50"
	   
          else if(cmul(i) .ge. 1. .and. cmul(i) .lt. 10.) then
           WRITE(ch_cmul,'(f4.1)') cmul(i)
           line1_end = "f 3 50"
	   
          else if(cmul(i) .ge. 0.1 .and. cmul(i) .lt. 1.) then
           WRITE(ch_cmul,'(f4.2)') cmul(i)
           line1_end = "f 4 50"
	   
          else if(cmul(i) .ge. 0.01 .and. cmul(i) .lt. 0.1) then
           WRITE(ch_cmul,'(f5.3)') cmul(i)
           line1_end = "f 4 60"
	   
          else if(cmul(i) .ge. 0.001 .and. cmul(i) .lt. 0.01) then
           WRITE(ch_cmul,'(f6.4)') cmul(i)
           line1_end = "f 5 150"
	   
          else if(cmul(i) .ge. 0.0001 .and. cmul(i) .lt. 0.001) then
           WRITE(ch_cmul,'(f7.5)') cmul(i)
           line1_end = "f 6 200"
	   
          else if(cmul(i) .ge. 0.00001 .and. cmul(i) .lt. 0.0001) then
           WRITE(ch_cmul,'(f8.6)') cmul(i)
           line1_end = "f 7 200"
	   
          else if(cmul(i) .ge. 0.000001 .and. cmul(i)
     1      .lt. 0.00001) then
           WRITE(ch_cmul,'(f9.7)') cmul(i)
           line1_end = "f 7 250"
	   
          else if(cmul(i) .ge. 0.0000001 .and. cmul(i)
     1     .lt. 0.000001) then
           WRITE(ch_cmul,'(f10.8)') cmul(i)
           line1_end = "f 8 250"
	   
          endif
c
c     Form text for subsequent system calls.  For each case, first
c      cd into local processor's subdirectory, run DKES, rename its
c      output file (to a less generic name that won't be overwritten
c      by the next run), then append that run's output on to the end
c      of the more global results file.  Note:  somewhat different
c      commands must be used, depending whether this is the first time
c      through or a subsequent run.
c
       cmd_line1 = "cd " // trim(adjustl(my_dirname)) //
     1      "; xdkes" // " " // trim(boozer_ext) // " " //
     2     trim(adjustl(ch_surf)) // " " // trim(adjustl(ch_cmul)) //
     3      " " // trim(adjustl(ch_num_efield)) // " " //
     4      trim(adjustl(line1_end))
       cmd_line2 = "cd " // trim(adjustl(my_dirname)) //
     1   "; mv results." // trim(boozer_ext) //
     2   " " // "results." // trim(boozer_ext) //
     3   "_s" // trim(adjustl(ch_surf)) // "_nu" //
     4   trim(adjustl(ch_cmul)) // "_E" // trim(adjustl(ch_efield(ie)))
      if(i .gt. 1) then
       cmd_line5 = "cd " // trim(adjustl(my_dirname)) //
     1   "; tail -1 results." // trim(boozer_ext) //
     2   "_s" // trim(adjustl(ch_surf)) // "_nu" //
     3   trim(adjustl(ch_cmul)) // "_E" //
     4   trim(adjustl(ch_efield(ie))) //
     4   " >> temp." // trim(boozer_ext) //
     5   "_s" // trim(adjustl(ch_surf)) //
     6  "_E" // trim(adjustl(ch_efield(ie)))
       else if(i .eq. 1) then
       cmd_line5 = "cd " // trim(adjustl(my_dirname)) //
     1   "; cp results." // trim(boozer_ext) //
     2   "_s" // trim(adjustl(ch_surf)) // "_nu" //
     3   trim(adjustl(ch_cmul)) // "_E" //
     4   trim(adjustl(ch_efield(ie))) //
     4   " temp." // trim(boozer_ext) //
     5   "_s" // trim(adjustl(ch_surf)) //
     6  "_E" // trim(adjustl(ch_efield(ie)))
       endif
       if(icount .gt. 1) then
       cmd_line3 = "cd " // trim(adjustl(my_dirname)) //
     1   "; tail -1 results." // trim(boozer_ext) //
     2   "_s" // trim(adjustl(ch_surf)) // "_nu" //
     3   trim(adjustl(ch_cmul)) // "_E" //
     4   trim(adjustl(ch_efield(ie))) //
     4   " >> temp"
       else if(icount .eq. 1) then
       cmd_line3 = "cd " // trim(adjustl(my_dirname)) //
     1   "; cp results." // trim(boozer_ext) //
     2   "_s" // trim(adjustl(ch_surf)) // "_nu" //
     3   trim(adjustl(ch_cmul)) // "_E" //
     4   trim(adjustl(ch_efield(ie))) //
     4   " temp"       
       
       endif
c
c     For debugging purposes, setting itest = 1 writes out the text
c      of the system calls into logfiles (rather than executing them)
c
       if (itest .eq. 1) then
       WRITE(10,*) cmd_line1
       WRITE(10,*) cmd_line2
       WRITE(10,*) cmd_line3
       WRITE(10,*) cmd_line5
c
c     Setting itest = 0 sends the text to the system calls and executes
c      the production run normally.
c
      else if(itest .eq. 0) then
       call system(trim(adjustl(cmd_line1)))
       call system(trim(adjustl(cmd_line2)))
       call system(trim(adjustl(cmd_line3)))
       call system(trim(adjustl(cmd_line5)))
      endif
      icount = icount + 1
      END DO           !i=1,num_cmul
      END DO           !i=1,num_efields
c
c     Create a final_results file and clean up:
c
      cmd_line4 ="cd " // trim(adjustl(my_dirname)) //
     1  "; mv temp final_results." // trim(boozer_ext)
      if(itest .eq. 1) then
       WRITE(10,*) cmd_line4
      else if(itest .eq. 0) then
       call system(trim(adjustl(cmd_line4)))
      endif
      if(itest .eq. 1) close(unit=10)
      call MPI_BARRIER(MPI_COMM_WORLD,IER)
      CALL MPI_FINALIZE(IER)
      END
