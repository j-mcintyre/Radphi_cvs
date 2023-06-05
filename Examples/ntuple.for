C This file serves as a "prototype" for paw
C since Fortran does not support prototypes
      subroutine ntuple(infile,count,id)
      character*100 infile
      integer count,id
      write(6,*) 'This should never be called!!'
      end
