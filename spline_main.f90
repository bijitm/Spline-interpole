      program spine_main
      implicit none
      integer, parameter :: mxgrid = 10000, mxcol = 10
      character(100) :: fbasename, fileinput, fileoutput, line
      character(3)   :: fextnsion
      integer :: nignore, nycol, ngrid, ngrid_new, igrid
      integer :: iline, icol
      real*8  :: x(mxgrid), y(mxcol,mxgrid), y2(mxcol,mxgrid)
      real*8  :: xfirst, xlast
      real*8  :: xmin, xmax, dx_new, x_new, y_new(mxcol)

      namelist/params/fbasename,fextnsion,nignore,nycol,dx_new, &
                      xfirst,xlast,ngrid_new

      ! Default values
      nignore = 0       ! no of lines containing text at the beginning of the file
      nycol   = 1       ! no of data columns
      dx_new  = -9999d0 ! new dx
      ngrid_new = -9999 ! new ngrid
      xfirst  = -9999d0 ! first data point
      xlast   = -9999d0 ! last data point

      read(5,params)
      write(*,*)"This is program for spline interpolation"
      write(*,*)

      if (nycol>mxcol) stop "nycol>mxcol. Increase mxcol."

      if (dx_new==-9999d0.and.ngrid_new==-9999) &
              stop "Provide either dx_new or ngrid_new"

      if (dx_new/=-9999d0.and.ngrid_new/=-9999) then
         print*, "Both dx_new and ngrid_new are provided: ", &
         "ngrid_new will be ignored"
         ngrid_new=-9999
      endif

      if (ngrid_new==0) stop "Cannot have ngrid_new = 0"

      fileinput=trim(fbasename)//'.'//trim(fextnsion)
      fileoutput=trim(fbasename)//'_interpolated.'//trim(fextnsion)

      open(1,file=fileinput,status='old')
      open(10,file=fileoutput,status='unknown')

      do iline=1,nignore
        read(1,'(A)') line
        write(10,'(A)') line
      enddo

      x  = 0.d0
      y  = 0.d0
      y2 = 0.d0

      ngrid=0
      do
        read(1,*, end=1)x(ngrid+1),y(:nycol,ngrid+1)
        ngrid=ngrid+1
      enddo
 1    continue
      close(1)

      xmin = min(x(1),x(ngrid))
      xmax = max(x(1),x(ngrid))

      if (ngrid==0) stop "No data point provided"

      if (ngrid>mxgrid) stop "ngrid larger than mxgrid"

      if (xfirst/=-9999d0.and.(xfirst<xmin.or.xfirst>xmax)) &
              stop "xfirst outside range"
      if (xlast/=-9999d0.and.(xlast<xmin.or.xlast>xmax)) &
              stop "xlast outside range"

      if (xfirst==-9999d0) xfirst = xmin
      if (xlast ==-9999d0) xlast  = xmax

      write(*,*) "No. of input data points:",ngrid

      do icol=1,nycol
        call spline(x(:ngrid),y(icol,:ngrid),ngrid,1.d30,1.d30,y2(icol,:ngrid))
      enddo

      if (dx_new==-9999d0.and.ngrid_new>1) dx_new = (xlast-xfirst)/dble(ngrid_new-1)
      if (ngrid_new==1) dx_new=0d0
      if (ngrid_new==-9999) ngrid_new = int((xlast-xfirst)/dx_new)+1

      write(*,*) "No. of new data points:",ngrid_new

      do igrid = 1, ngrid_new
        x_new = xfirst + dble(igrid-1)*dx_new
        do icol=1,nycol
          call splint(x(:ngrid),y(icol,:ngrid),y2(icol,:ngrid),ngrid,x_new,y_new(icol))
        enddo
        write(10,'(f18.8,10es15.6)')x_new, y_new(:nycol)
      enddo
      close(10)

      write(*,*)"Interpolation successfully done and written to file: ",trim(fileoutput)

      end program spine_main
