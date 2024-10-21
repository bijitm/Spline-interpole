# A simple makefile to compile the code
main = spline-interpol.f90

spline: $(main) 
	  gfortran $(main) -o spline

