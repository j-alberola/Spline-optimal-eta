

GF=	gfortran

Spline: Program_roots
	$(GF) -o Spline.x Spline.f90 Roots.o

Program_roots:
	$(GF) -c -o Roots.o Program_third_party.f90
