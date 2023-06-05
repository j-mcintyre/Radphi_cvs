*---------------------------------------------
* mctuplet.F - a stand-alone host for mctuple.c
* ---------
* Richard Jones
* March 28, 2000
*---------------------------------------------

#define COUNT_IS_CUMULATIVE 1
#define DEFAULT_EVENT_COUNT 10000
#define DEFAULT_INIT_STRING 'init'
#define OUTPUT_FILE_FULL 95
#define OUTPUT_FILE_OVERRUN 96

      program mctuplet
      implicit none
      integer hspace
      parameter (hspace=60000000)
      common /pawc/ hq(hspace)
      common /quest/ iquest(100)
      real hq
      integer iquest
      integer info(13)
      integer istat,icycle
      integer count,start,nleft
      integer ret
      integer id
      integer noutfile
      character*1000 cmd
      character*1000 outfile1,outfile2
      character*1000 outfile
      character*1000 initstring
      character*8000 argv
      integer mctuple
      external mctuple
      integer renamef,readlnf
      external renamef,readlnf
      integer lenocc
      external lenocc
      logical firstfile
      integer i,n
      integer bsize
      data bsize/65536/

C---- Initialization of HBOOK
      call HLIMIT(hspace)
      call HBSET('BSIZE',bsize,istat)

C---- Read and process input data files
      initstring = DEFAULT_INIT_STRING
      n = mctuple(initstring,0,0,0)
      if (n.ne.0) then
        call exit(9)
      endif
      outfile='mctuple.hbook'
      noutfile=0
      firstfile=.TRUE.
      count=DEFAULT_EVENT_COUNT
      nleft=count
      start=0
      id=1
      do i=1,iargc()
        call getarg(i,argv)
        if (index(argv,'--count=').eq.1) then
          read(argv(9:),*) count
          nleft = count
        elseif (index(argv,'--start=').eq.1) then
          read(argv(9:),*) start
        elseif (index(argv,'--id=').eq.1) then
          read(argv(6:),*) id
        elseif (index(argv,'--output=').eq.1) then
          read(argv(10:),'(a1000)') outfile
          noutfile=0
        elseif (argv(1:2).eq.'--') then
          write(6,*) 'Usage: mctuplet [--count=nnnn] [--start=nnn] '//
     +          '[--id=nn] [--output=fn] [-guclvbxtpmoe] [file1 ...'
          call exit(9)
        elseif (index(argv,'-').eq.1) then
          initstring = 'init ' // argv(2:)
          n = mctuple(initstring,0,0,0)
          if (n.ne.0) then
            call exit(9)
          endif
        else
          if (firstfile) then
            open(unit=50,file=outfile,status='old',err=9)
            close(unit=50)
            call HROPEN(50,'RZfile',outfile,'U',bsize,istat)
            call HSCR(id,0,'')
            go to 10
    9       iquest(10) = 256000   ! extend RZ quota to 2^32 bits
            call HROPEN(50,'RZfile',outfile,'NQE',bsize,istat)
   10       firstfile=.FALSE.
          endif
   11     continue
#if !COUNT_IS_CUMULATIVE
          nleft = count
#endif
          ret = mctuple(argv,id,nleft,start)
          if (ret.lt.0) then
            call exit(9)
#if COUNT_IS_CUMULATIVE
          elseif (nleft.le.0) then
            go to 99
#endif
          elseif (ret.eq.OUTPUT_FILE_FULL) then
            call HROUT(0,icycle,' ')
            call HREND('RZfile')
            call HDELET(0)
            close(unit=50)
            if (noutfile.lt.10) then
              write(outfile2,'(2a,i1,1a)') outfile(:lenocc(outfile)),
     +                                     ',',noutfile
            elseif (noutfile.lt.100) then
              write(outfile2,'(2a,i2,1a)') outfile(:lenocc(outfile)),
     +                                     ',',noutfile
            elseif (noutfile.lt.1000) then
              write(outfile2,'(2a,i3,1a)') outfile(:lenocc(outfile)),
     +                                     ',',noutfile
            else
              stop 'OK, this is getting a little extreme!'
            endif
            noutfile=noutfile+1
            write(6,*) 'mctuplet warning- ',
     +                 'out of space on output file:'
            write(6,*) 'saving ',outfile(:lenocc(outfile)),
     +                 ' as ',outfile2(:lenocc(outfile2)),
     +                 ' and going on.'
            n=readlnf(outfile,outfile1)
            if (n.le.0) then
              outfile1=outfile
            else
              outfile1=outfile1(:n)
            endif
            write(cmd,'(4a)') 'mv ',outfile1(:lenocc(outfile1)),
     +                        ' ',outfile2(:lenocc(outfile2))
            call system(cmd,istat)
            if (istat.ne.0) then
              write(7,*) 'I/O Error in command: ',cmd(:lenocc(cmd))
              call exit(8)
            endif
            iquest(10) = 256000   ! extend RZ quota to 2^32 bits
            call HROPEN(50,'RZfile',outfile,'NQE',bsize,istat)
            goto 11
          elseif (ret.eq.OUTPUT_FILE_OVERRUN) then
            call exit(5)
          endif
        endif
      enddo

C---- Termination of HBOOK
   99 call HROUT(0,icycle,' ')
      call HREND('RZfile')
      call exit(0)
      END
