# A simple makefile to compile the code
main = spline_main.f90
subroutine = spline.f90
all = $(subroutine) $(main)
spline: $(all)
	  gfortran $(all) -o interpole

