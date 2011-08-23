c     =============================================================
c                    NUMERICAL RECIPES ROUTINES
c     =============================================================

c     --------------------------------------------------------------
c                Brent's algorithm for root finding
c     --------------------------------------------------------------
      FUNCTION zbrent(func,x1,x2,tol)
      INTEGER ITMAX
      DOUBLEPRECISION zbrent,tol,x1,x2,func,EPSS
      EXTERNAL func
      PARAMETER (ITMAX=100,EPSS=3.d-8)
      INTEGER iter
      DOUBLEPRECISION a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
      a=x1
      b=x2
      fa=func(a)
      fb=func(b)
      if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.))stop
     *'root must be bracketed for zbrent'
      c=b
      fc=fb
      do 11 iter=1,ITMAX
        if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
          c=a
          fc=fa
          d=b-a
          e=d
        endif
        if(abs(fc).lt.abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
        endif
        tol1=2.*EPSS*abs(b)+0.5*tol
        xm=.5*(c-b)
        if(abs(xm).le.tol1 .or. fb.eq.0.)then
          zbrent=b
          return
        endif
        if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
          s=fb/fa
          if(a.eq.c) then
            p=2.*xm*s
            q=1.-s
          else
            q=fa/fc
            r=fb/fc
            p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
            q=(q-1.)*(r-1.)*(s-1.)
          endif
          if(p.gt.0.) q=-q
          p=abs(p)
          if(2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) then
            e=d
            d=p/q
          else
            d=xm
            e=d
          endif
        else
          d=xm
          e=d
        endif
        a=b
        fa=fb
        if(abs(d) .gt. tol1) then
          b=b+d
        else
          b=b+sign(tol1,xm)
        endif
        fb=func(b)
11    continue
      stop 'zbrent exceeding maximum iterations'
      zbrent=b
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..


c     -------------------------------------------------------------
c                          MAKE A SPLINE
c     -------------------------------------------------------------
      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      DOUBLEPRECISION yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=500)
      INTEGER i,k
      DOUBLEPRECISION p,qn,sig,un,u(NMAX)
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+
     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
     *u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..


c     -------------------------------------------------------------
c                     SUBSTITUTE A SPLINE
c     -------------------------------------------------------------
      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER n
      DOUBLEPRECISION x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      DOUBLEPRECISION a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) stop 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
     *2)/6.
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..


c     -------------------------------------------------------------
c               HUNT FOR LOCATION OF A NUMBER IN A TABLE
c     -------------------------------------------------------------
      subroutine hunt(xx,n,x,jlo)
      integer jlo,n
      doubleprecision x,xx(n)
      integer inc,jhi,jm
      logical ascnd
      ascnd=xx(n).gt.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        goto 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
1       jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          goto 1
        endif
      else
        jhi=jlo
2       jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo
          inc=inc+inc
          goto 2
        endif
      endif
3     if(jhi-jlo.eq.1)return
      jm=(jhi+jlo)/2
      if(x.gt.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      goto 3
      end
C  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..



c     -------------------------------------------------------------
c               HUNT FOR LOCATION OF A NUMBER IN A TABLE
c                  (same as hunt, but other version)
c     -------------------------------------------------------------
      subroutine hunt0(xx,n,x,jlo)
      implicit none
      doubleprecision xx,x
      integer jhi,jlo,jm,inc,n
      dimension xx(0:n)
      logical ascnd
      ascnd=xx(n).gt.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        go to 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
1       jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          go to 1
        endif
      else
        jhi=jlo
2       jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo
          inc=inc+inc
          go to 2
        endif
      endif
3     if(jhi-jlo.eq.1)return
      jm=(jhi+jlo)/2
      if(x.gt.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      go to 3
      end


c     -------------------------------------------------------------
c                         BESSEL FUNCTIONS
c                     (Num Recip 2nd edition)
c     -------------------------------------------------------------
      FUNCTION BESSI0(X)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 Y,P1,P2,P3,P4,P5,P6,P7,
     *    Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9
      DATA P1,P2,P3,P4,P5,P6,P7/1.0D0,3.5156229D0,3.0899424D0,1.2067492D
     *0,
     *    0.2659732D0,0.360768D-1,0.45813D-2/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1,
     *    0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,
     *    0.2635537D-1,-0.1647633D-1,0.392377D-2/
      IF (ABS(X).LT.3.75D0) THEN
        Y=(X/3.75D0)**2
        BESSI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
      ELSE
        AX=ABS(X)
        Y=3.75/AX
        BESSI0=(EXP(AX)/SQRT(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4
     *      +Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
      ENDIF
      RETURN
      END

      FUNCTION BESSI1(X)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 Y,P1,P2,P3,P4,P5,P6,P7,
     *    Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9
      DATA P1,P2,P3,P4,P5,P6,P7/0.5D0,0.87890594D0,0.51498869D0,
     *    0.15084934D0,0.2658733D-1,0.301532D-2,0.32411D-3/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,-0.3988024D-1,
     *    -0.362018D-2,0.163801D-2,-0.1031555D-1,0.2282967D-1,
     *    -0.2895312D-1,0.1787654D-1,-0.420059D-2/
      IF (ABS(X).LT.3.75D0) THEN
        Y=(X/3.75D0)**2
        BESSI1=X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
        AX=ABS(X)
        Y=3.75D0/AX
        BESSI1=(EXP(AX)/SQRT(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+
     *      Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
      ENDIF
      RETURN
      END

      FUNCTION bessj1(x)
      doubleprecision bessj1,x
      doubleprecision ax,xx,z
      DOUBLE PRECISION p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,
     *s1,s2,s3,s4,s5,s6,y
      SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,
     *s5,s6
      DATA r1,r2,r3,r4,r5,r6/72362614232.d0,-7895059235.d0,
     *242396853.1d0,-2972611.439d0,15704.48260d0,-30.16036606d0/,s1,s2,
     *s3,s4,s5,s6/144725228442.d0,2300535178.d0,18583304.74d0,
     *99447.43394d0,376.9991397d0,1.d0/
      DATA p1,p2,p3,p4,p5/1.d0,.183105d-2,-.3516396496d-4,
     *.2457520174d-5,-.240337019d-6/, q1,q2,q3,q4,q5/.04687499995d0,
     *-.2002690873d-3,.8449199096d-5,-.88228987d-6,.105787412d-6/
      if(abs(x).lt.8.)then
        y=x**2
        bessj1=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+
     *y*(s4+y*(s5+y*s6)))))
      else
        ax=abs(x)
        z=8./ax
        y=z**2
        xx=ax-2.356194491
        bessj1=sqrt(.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*
     *p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))*sign(1.d0,x)
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..



c     -------------------------------------------------------------
c                      MODIFIED BESSEL FUNCTIONS
c                      (Num Recip 2nd edition)
c     But these are the K_n(x) Bessel functions divided by exp(-x), 
c     because of numerical practicality.
c     -------------------------------------------------------------
      FUNCTION BESSK0EXP(X)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 Y,P1,P2,P3,P4,P5,P6,P7,
     *    Q1,Q2,Q3,Q4,Q5,Q6,Q7
      DATA P1,P2,P3,P4,P5,P6,P7/-0.57721566D0,0.42278420D0,0.23069756D0,
     *    0.3488590D-1,0.262698D-2,0.10750D-3,0.74D-5/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,-0.7832358D-1,0.2189568D-1,
     *    -0.1062446D-1,0.587872D-2,-0.251540D-2,0.53208D-3/
      IF (X.LE.2.0D0) THEN
        Y=X*X/4.0D0
        BESSK0EXP=((-LOG(X/2.0D0)*BESSI0(X))+(P1+Y*(P2+Y*(P3+
     *        Y*(P4+Y*(P5+Y*(P6+Y*P7)))))))*EXP(X)
      ELSE
        Y=(2.0D0/X)
        BESSK0EXP=(1.D0/SQRT(X))*(Q1+Y*(Q2+Y*(Q3+
     *        Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
      ENDIF
      RETURN
      END

      FUNCTION BESSK1EXP(X)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 Y,P1,P2,P3,P4,P5,P6,P7,
     *    Q1,Q2,Q3,Q4,Q5,Q6,Q7
      DATA P1,P2,P3,P4,P5,P6,P7/1.0D0,0.15443144D0,-0.67278579D0,
     *    -0.18156897D0,-0.1919402D-1,-0.110404D-2,-0.4686D-4/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,0.23498619D0,-0.3655620D-1,
     *    0.1504268D-1,-0.780353D-2,0.325614D-2,-0.68245D-3/
      IF (X.LE.2.0D0) THEN
        Y=X*X/4.0D0
        BESSK1EXP=((LOG(X/2.0D0)*BESSI1(X))+(1.0D0/X)*(P1+Y*(P2+
     *      Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))))*EXP(X)
      ELSE
        Y=2.0D0/X
        BESSK1EXP=(1.D0/SQRT(X))*(Q1+Y*(Q2+Y*(Q3+
     *      Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
      ENDIF
      RETURN
      END

      FUNCTION BESSKEXP(N,X)
      IMPLICIT REAL*8(A-H,O-Z)
      IF (N.LT.2) STOP 'bad argument N in BESSK'
      TOX=2.0D0/X
      BKM=BESSK0EXP(X)
      BK=BESSK1EXP(X)
      DO 11 J=1,N-1
        BKP=BKM+J*TOX*BK
        BKM=BK
        BK=BKP
11    CONTINUE
      BESSKEXP=BK
      RETURN
      END


c     -------------------------------------------------------------
c               THE ORIGINAL BESSEL FUNCTIONS (UNRENORMALIZED)
c     -------------------------------------------------------------
      FUNCTION bessk(n,x)
      INTEGER n
      doubleprecision  bessk,x
CU    USES bessk0,bessk1
      INTEGER j
      doubleprecision  bk,bkm,bkp,tox,bessk0,bessk1
      if (n.lt.2) stop 'bad argument n in bessk'
      tox=2.0/x
      bkm=bessk0(x)
      bk=bessk1(x)
      do 11 j=1,n-1
        bkp=bkm+j*tox*bk
        bkm=bk
        bk=bkp
11    continue
      bessk=bk
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..

      FUNCTION bessk0(x)
      DOUBLEPRECISION bessk0,x
CU    USES bessi0
      DOUBLEPRECISION bessi0
      DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,y
      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7
      DATA p1,p2,p3,p4,p5,p6,p7/-0.57721566d0,0.42278420d0,0.23069756d0,
     *0.3488590d-1,0.262698d-2,0.10750d-3,0.74d-5/
      DATA q1,q2,q3,q4,q5,q6,q7/1.25331414d0,-0.7832358d-1,0.2189568d-1,
     *-0.1062446d-1,0.587872d-2,-0.251540d-2,0.53208d-3/
      if (x.le.2.0) then
        y=x*x/4.0
        bessk0=(-log(x/2.0)*bessi0(x))+(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*
     *(p6+y*p7))))))
      else
        y=(2.0/x)
        bessk0=(exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*
     *q7))))))
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..

      FUNCTION bessk1(x)
      DOUBLEPRECISION bessk1,x
CU    USES bessi1
      DOUBLEPRECISION bessi1
      DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,y
      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7
      DATA p1,p2,p3,p4,p5,p6,p7/1.0d0,0.15443144d0,-0.67278579d0,
     *-0.18156897d0,-0.1919402d-1,-0.110404d-2,-0.4686d-4/
      DATA q1,q2,q3,q4,q5,q6,q7/1.25331414d0,0.23498619d0,-0.3655620d-1,
     *0.1504268d-1,-0.780353d-2,0.325614d-2,-0.68245d-3/
      if (x.le.2.0) then
        y=x*x/4.0
        bessk1=(log(x/2.0)*bessi1(x))+(1.0/x)*(p1+y*(p2+y*(p3+y*(p4+y*
     *(p5+y*(p6+y*p7))))))
      else
        y=2.0/x
        bessk1=(exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*
     *q7))))))
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..

      FUNCTION gammln(xx)
      DOUBLEPRECISION gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..


c     -------------------------------------------------------------
c             BESSEL FUNCTIONS USED BY THE COMPTON MODULE
c     These have been messed up in a way that is useful for the
c     Compton CSK module...
c     -------------------------------------------------------------
      FUNCTION BESKCSK(N,X)
      IMPLICIT REAL*8(A-H,O-Z)
      DATA PI2/1.570796327D0/
      TOX=2.0/X
      BKM=BESSKCSK0(X)
      BK=BESSKCSK1(X)
      DO 11 J=1,N-1
        BKP=BKM+J*TOX*BK
        BKM=BK
        BK=BKP
11    CONTINUE
      BESKCSK=BK
      RETURN
      END
   

      FUNCTION BESSKCSK0(X)
      IMPLICIT REAL*8(A-H,O-Z)
      DATA P1,P2,P3,P4,P5,P6,P7/-0.57721566D0,0.42278420D0,0.23069756D0,
     *    0.3488590D-1,0.262698D-2,0.10750D-3,0.74D-5/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,-0.7832358D-1,0.2189568D-1,
     *    -0.1062446D-1,0.587872D-2,-0.251540D-2,0.53208D-3/
      IF (X.LE.2.0) THEN
        Y=X*X/4.0
        BESSKCSK0=(-LOG(X/2.0)*BESSI0(X))+(P1+Y*(P2+Y*(P3+
     *        Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
        Y=(2.0/X)
C        BESSKCSK0=(EXP(-X)/SQRT(X))*(Q1+Y*(Q2+Y*(Q3+
C     *        Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
        BESSKCSK0=(1.D0/SQRT(X))*(Q1+Y*(Q2+Y*(Q3+
     *        Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
      ENDIF
      RETURN
      END
 
      FUNCTION BESSKCSK1(X)
      IMPLICIT REAL*8(A-H,O-Z)
      DATA P1,P2,P3,P4,P5,P6,P7/1.0D0,0.15443144D0,-0.67278579D0,
     *    -0.18156897D0,-0.1919402D-1,-0.110404D-2,-0.4686D-4/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,0.23498619D0,-0.3655620D-1,
     *    0.1504268D-1,-0.780353D-2,0.325614D-2,-0.68245D-3/
      IF (X.LE.2.0) THEN
        Y=X*X/4.0
        BESSKCSK1=(LOG(X/2.0)*BESSI1(X))+(1.0/X)*(P1+Y*(P2+
     *      Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
         Y=2.0/X
C        BESSKCSK1=(EXP(-X)/SQRT(X))*(Q1+Y*(Q2+Y*(Q3+
C     *      Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
        BESSKCSK1=(1.D0/SQRT(X))*(Q1+Y*(Q2+Y*(Q3+
     *      Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
      ENDIF
      RETURN
      END


      FUNCTION BESKCSK2(X)
      IMPLICIT REAL*8(A-H,O-Z)
c          Returns the value of  K2(x)*exp(x)  <==============
c             (Abramowitz & Stegun p.379)
c
      BESKCSK2 = BESSKCSK0(X) + 2.D0*BESSKCSK1(X)/X
      RETURN
      END
 

c     --------------------------------------------------------------
c              POLINOMIAL INTERPOLATION FOR THREE POINTS
c     --------------------------------------------------------------
      subroutine polint3(xa,ya,x,y)
      doubleprecision x,y,xa(3),ya(3)

      y= (((x-xa(2))*(x-xa(3)))/((xa(1)-xa(2))*(xa(1)-xa(3))))*ya(1)
     % + (((x-xa(1))*(x-xa(3)))/((xa(2)-xa(1))*(xa(2)-xa(3))))*ya(2)
     % + (((x-xa(1))*(x-xa(2)))/((xa(3)-xa(1))*(xa(3)-xa(2))))*ya(3) 

      return
      end


c     --------------------------------------------------------------
c             POLINOMIAL INTERPOLATION FOR FOUR POINTS
c     --------------------------------------------------------------
      subroutine polint4(xa,ya,x,y)
      doubleprecision x,y,xa(4),ya(4)
      
      y = (((x-xa(2))*(x-xa(3))*(x-xa(4)))
     % /((xa(1)-xa(2))*(xa(1)-xa(3))*(xa(1)-xa(4))))*ya(1)  
 
     % + (((x-xa(1))*(x-xa(3))*(x-xa(4)))
     % /((xa(2)-xa(1))*(xa(2)-xa(3))*(xa(2)-xa(4))))*ya(2) 

     % + (((x-xa(1))*(x-xa(2))*(x-xa(4)))
     % /((xa(3)-xa(1))*(xa(3)-xa(2))*(xa(3)-xa(4))))*ya(3) 

     % + (((x-xa(1))*(x-xa(2))*(x-xa(3)))
     % /((xa(4)-xa(1))*(xa(4)-xa(2))*(xa(4)-xa(3))))*ya(4)

      return
      end


c     --------------------------------------------------------------
c          NUMERICAL RECIPES 2*N-POINT GAUSS-LEGENDRE INTEGRATION
c     --------------------------------------------------------------
      SUBROUTINE qgausn(func,a,b,n,ss)
      DOUBLEPRECISION a,b,ss,func
      INTEGER n
      EXTERNAL func
      INTEGER j
      DOUBLEPRECISION dx,xm,xr,w(50),x(50)
      if(2*n.gt.50) then
          write(*,*) 'ERROR in qgausn: 2*n>50'
          stop 13
      endif
      call gauleg(a,b,x,w,2*n)
      xm=0.5*(b+a)
      xr=0.5*(b-a)/(1.d0*n)
      ss=0
      do 11 j=1,2*n
        ss=ss+w(j)*func(x(j))
11    continue
c      stop 13
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..

c     --------------------------------------------------------------
c          NUMERICAL RECIPES GAUSS-LEGENDRE WEIGHTS AND ABCISSAS
c     --------------------------------------------------------------
      SUBROUTINE gauleg(x1,x2,x,w,n)
      INTEGER n
      DOUBLEPRECISION x1,x2,x(n),w(n)
      DOUBLE PRECISION EPS
      PARAMETER (EPS=3.d-14)
      INTEGER i,j,m
      DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
        z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        continue
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
        if(abs(z-z1).gt.EPS)goto 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..

c     --------------------------------------------------------
c      LU Decomposition
c     --------------------------------------------------------
      SUBROUTINE ludcmp(a,n,np,indx,d)
      INTEGER n,np,indx(n),NMAX
      DOUBLEPRECISION d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0d-20)
      INTEGER i,imax,j,k
      DOUBLEPRECISION aamax,dum,sum,vv(NMAX)
      d=1.
      do 12 i=1,n
        aamax=0.d0
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.d0) stop 'singular matrix in ludcmp'
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.d0
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.d0)a(j,j)=TINY
        if(j.ne.n)then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..


c     --------------------------------------------------------
c      LU Back substitution
c     --------------------------------------------------------
      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      DOUBLEPRECISION a(np,np),b(n)
      INTEGER i,ii,j,ll
      DOUBLEPRECISION sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.d0) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..



c     -------------------------------------------------------------
c                        SORTING ROUTINES
c     -------------------------------------------------------------
      SUBROUTINE indexx(n,arr,indx)
      INTEGER n,indx(n),M,NSTACK
      DOUBLEPRECISION arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      DOUBLEPRECISION a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK)stop 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
C  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..


      SUBROUTINE ray_sort(n,ra,rb,rc,ird,ire,wksp,iwksp,iiwksp)
      INTEGER n,iwksp(n),iiwksp(n),ird(n),ire(n)
      DOUBLEPRECISION ra(n),rb(n),rc(n),wksp(n)
CU    USES indexx
      INTEGER j
      call indexx(n,ra,iwksp)
      do 11 j=1,n
        wksp(j)=ra(j)
11    continue
      do 12 j=1,n
        ra(j)=wksp(iwksp(j))
12    continue
      do 13 j=1,n
        wksp(j)=rb(j)
13    continue
      do 14 j=1,n
        rb(j)=wksp(iwksp(j))
14    continue
      do 15 j=1,n
        wksp(j)=rc(j)
15    continue
      do 16 j=1,n
        rc(j)=wksp(iwksp(j))
16    continue
      do 17 j=1,n
        iiwksp(j)=ird(j)
17    continue
      do 18 j=1,n
        ird(j)=iiwksp(iwksp(j))
18    continue
      do 19 j=1,n
        iiwksp(j)=ire(j)
19    continue
      do 20 j=1,n
        ire(j)=iiwksp(iwksp(j))
20    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..


      SUBROUTINE sort(n,arr)
      INTEGER n,M,NSTACK
      DOUBLEPRECISION arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      DOUBLEPRECISION a,temp
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          do 11 i=j-1,1,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
11        continue
          i=0
2         arr(i+1)=a
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l+1).gt.arr(l))then
          temp=arr(l+1)
          arr(l+1)=arr(l)
          arr(l)=temp
        endif
        i=l+1
        j=ir
        a=arr(l)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        goto 3
5       arr(l)=arr(j)
        arr(j)=a
        jstack=jstack+2
        if(jstack.gt.NSTACK)stop 'NSTACK too small in sort'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
C  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..


c     --------------------------------------------------------------
c                    SOME MORE SPECIAL FUNCTIONS
c     --------------------------------------------------------------
      FUNCTION EI(X)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION AS(6),AL(4)
      DATA AS/-0.57721566D0,0.99999193D0,-0.24991055D0,
     %           0.05519968D0,-0.00976004D0,0.00107857D0/
      DATA AL/2.334733D0,0.250621D0,3.330657D0,1.681534D0/
c
c     this is the exponential-integral function Ei(x)
c
      IF(X.LT.1.) THEN
       EI = -LOG(X)
       DO I = 1,6
        EI = EI + AS(I)*X**(I-1) 
       ENDDO
       EI = EI*EXP(X)
      ELSE
       EI = (X*X+AL(1)*X+AL(2))/(X*X+AL(3)*X+AL(4))
       EI = EI/X
      ENDIF
      EI = -EI
      RETURN
      END



      FUNCTION GAMMP(A,X)
      IMPLICIT REAL*8(A-H,O-Z)
      IF(X.LT.0..OR.A.LE.0.)STOP
      IF(X.LT.A+1.d0)THEN
        CALL GSER(GAMMP,A,X,GLN)
      ELSE
        CALL GCF(GAMMCF,A,X,GLN)
        GAMMP=1.d0-GAMMCF
      ENDIF
      RETURN
      END
   

      FUNCTION GAMMQ(A,X)
      IMPLICIT REAL*8(A-H,O-Z)
      IF(X.LT.0..OR.A.LE.0.)STOP
      IF(X.LT.A+1.d0)THEN
        CALL GSER(GAMSER,A,X,GLN)
        GAMMQ=1.d0-GAMSER
      ELSE
        CALL GCF(GAMMQ,A,X,GLN)
      ENDIF
      RETURN
      END
   

      SUBROUTINE GSER(GAMSER,A,X,GLN)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (ITMAX=100,EPS=3.D-7)
      GLN=GAMMLN(A)
      IF(X.LE.0.)THEN
        IF(X.LT.0.)STOP
        GAMSER=0.
        RETURN
      ENDIF
      AP=A
      SUM=1./A
      DEL=SUM
      DO 11 N=1,ITMAX
        AP=AP+1.d0
        DEL=DEL*X/AP
        SUM=SUM+DEL
        IF(ABS(DEL).LT.ABS(SUM)*EPS)GO TO 1
11    CONTINUE
      STOP 'A too large, ITMAX too small'
1     GAMSER=SUM*EXP(-X+A*LOG(X)-GLN)
      RETURN


      END

      SUBROUTINE GCF(GAMMCF,A,X,GLN)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (ITMAX=100,EPS=3.D-7)
      GLN=GAMMLN(A)
      GOLD=0.
      A0=1.
      A1=X
      B0=0.
      B1=1.
      FAC=1.
      DO 11 N=1,ITMAX
        AN=FLOAT(N)
        ANA=AN-A
        A0=(A1+A0*ANA)*FAC
        B0=(B1+B0*ANA)*FAC
        ANF=AN*FAC
        A1=X*A0+ANF*A1
        B1=X*B0+ANF*B1
        IF(A1.NE.0.)THEN
          FAC=1./A1
          G=B1*FAC
          IF(ABS((G-GOLD)/G).LT.EPS)GO TO 1
          GOLD=G
        ENDIF
11    CONTINUE
      STOP 'A too large, ITMAX too small'
1     GAMMCF=EXP(-X+A*LOG(X)-GLN)*G
      RETURN
      END
      

      FUNCTION ERFC(X)
      IMPLICIT REAL*8(A-H,O-Z)
      IF(X.LT.0.)THEN
        ERFC=1.+GAMMP(.5d0,X**2)
      ELSE
        ERFC=GAMMQ(.5d0,X**2)
      ENDIF
      RETURN
      END


c     --------------------------------------------------------------
c                      MAKE INDEXED FILENAME
c
c     This routine makes filenames of the type dusttemp_4.dat etc.
c     --------------------------------------------------------------
      subroutine make_indexed_filename(base,index,ext,filename)
      implicit none
      character*80 base
      integer index
      character*80 ext
      character*80 filename
      character*12 ch
c     
      if((index.lt.0).or.(index.ge.1000)) then
          write(*,*) 'ERROR in make_indexed_filename()'
          stop 729
      endif
      if(index.lt.10) then
          write(ch,11) index
 11       format(I1)
      elseif(index.lt.100) then
          write(ch,12) index
 12       format(I2)
      elseif(index.lt.1000) then
          write(ch,13) index
 13       format(I3)
      endif
      filename = base(1:len_trim(base))//
     %           ch(1:len_trim(ch))//
     %           ext(1:len_trim(ext))
c
      end


