

c     --------------------------------------------------------------
c              STRIP THE COMMENT LINES FROM AN INPUT FILE
c     --------------------------------------------------------------
      subroutine strip_comments(filein,fileout)
      character*60 dd
      integer filein,fileout
 10   continue
      read(filein,50,end=20) dd
      if((dd(1:1).ne.'#').and.(dd(1:1).ne.';').and.
     %   (dd(1:1).ne.'%').and.(dd(1:1).ne.'=')) then
          write(fileout,50) dd
      endif
 50   format(a16)
      goto 10
 20   continue
      end




c     =============================================================
c          SOME PLATFORM-DEPENDENT FUNCTIONS AND SUBROUTINES
c     =============================================================


#ifdef COMPILER_XLAIX
c     --------------------------------------------------------------
c                         a new sign function
c     --------------------------------------------------------------
      function sgnf(a,b)
      implicit none
      doubleprecision a,b,sgnf
c
      if(a*b.lt.0.d0) then
         sgnf = -1.d0*a
      else
         sgnf = 1.d0*a
      endif
c     
      return
      end
c
#endif


c     --------------------------------------------------------------
c                   FUNCTION: TEST VALIDITY OF NUMBER
c
c     0 = Number is okay
c     1 = Number is INF
c     2 = Number is NAN
c     --------------------------------------------------------------
      function number_invalid(a)
      implicit none
      doubleprecision a,b,c
      integer number_invalid
      logical div,sub
c
      b=a*2.d0
      b=b/2.d0
      c=a-1.d100
c
      div = (b.eq.a)
      sub = (c.lt.a)
c
      if(div.and.sub) then
          number_invalid = 0
      elseif(div) then
          number_invalid = 1
      else
          number_invalid = 2
      endif
c
      return
      end
