      program main
      implicit none
      real*8  :: x(1000), y(1000), y2(1000)
      real*8  :: xn, yn, q, xmin, xmax
      integer :: i, n

      x  = 0.d0
      y  = 0.d0
      y2 = 0.d0

      open(1,file="PEC_shifted.dat", status="old")

      read(1,*)
      read(1,*)
      n=0
      do
        read(1,*, end=1)x(n+1),q,y(n+1)
        n=n+1
      enddo
 1    continue
      close(1)

      print*,"No. of original data points:",n

      call spline(x(:n),y(:n),n,1.d30,1.d30,y2(:n))

      print*,"Provide new point:"
      read(*,*)xn

      xmin=x(1)
      xmax=x(n)
      if(xn<xmin .or. xn>xmax) &
     & stop "STOP: new point outside the range of original points"

      call splint(x(:n),y(:n),y2(:n),n,xn,yn)
      print*,"Value at new point:",yn

      end program main
! C------------------------------ spline --------------------------------
      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION x(n),y(n),y2(n)
      DIMENSION u(N)
      if (yp1.gt..99e30) then
      y2(1)=0.
      u(1)=0.
      else
      y2(1)=-0.5
      u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do i=2,n-1
      sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
      p=sig*y2(i-1)+2.
      y2(i)=(sig-1.)/p
      u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
      &/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      enddo
      if (ypn.gt..99e30) then
      qn=0.
      un=0.
      else
      qn=0.5
      un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do k=n-1,1,-1
      y2(k)=y2(k)*y2(k+1)+u(k)
      enddo
      return
      END
!!!
      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION xa(n),y2a(n),ya(n)
      klo=1
      khi=n
    1 if (khi-klo.gt.1) then
      k=(khi+klo)/2
      if(xa(k).gt.x)then
      khi=k
      else
      klo=k
      endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) STOP "bad xa input in splint"
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      return
      END        