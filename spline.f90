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
