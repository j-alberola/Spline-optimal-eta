 module spline_module

use PolynomialRoots

implicit none


integer, parameter ::dp = kind(0.d0)

!
! CLOSE THE FILES AFTER THE SUBROUTINES
!

contains

!Read the points in file and storing data in a matrix (supossing file contains:
! eta_value real_energy imag_ener -eta*real_derivatvie -eta*imag_deriv
!
!
! MAYBE STORE EACH VARIABLE SEPARETLY IN A VECTOR
!
subroutine Reading (M,num_points)

real*8, allocatable, intent(out) :: M (:,:)
integer, intent(out) :: num_points
integer :: i, j, ios, stat
character(len=1000) :: file_name
character(len=1) :: foo


write (6,*) 'Insert file name'
read(5,*) file_name

open (UNIT=10, FILE=file_name,ACTION='read',IOSTAT=stat,STATUS='OLD')
if (stat .ne. 0) then
        write (*,*) 'File not found. Introduce a valid name'
        stop
end if 

num_points = 0
do
        read (10,*,iostat=ios) foo
        if (ios /= 0) exit
        num_points = num_points + 1
end do


Rewind 10
write (*,*) NUM_POINTS
allocate (M(num_points, 7))
do i = 1, num_points
        read (10, *) (M (i, j), j = 1,5)
!
!  In the quantum package output, the derivatives are stored as -eta*(dE/deta)
!
        M (i,4) = -M(i,4)/M(i,1)
        M (i,5) = -M(i,5)/M(i,1)
end do


!
! Calculate second derivatives
! PUT INTO PROPER LOOPS 

M(1,6) = (M(2,4)-M(1,4))/(M(2,1)-M(1,1))
M(1,7) = (M(2,5)-M(1,5))/(M(2,1)-M(1,1))


do i = 2, num_points-1

M(i,6) = (M(i+1,4)-M(i-1,4))/(M(i+1,1)-M(i-1,1))
M(i,7) = (M(i+1,5)-M(i-1,5))/(M(i+1,1)-M(i-1,1))

end do


M(num_points,6) = (M(num_points,4)-M(num_points-1,4))/(M(num_points,1)-M(num_points-1,1))
M(num_points,7) = (M(num_points,5)-M(num_points-1,5))/(M(num_points,1)-M(num_points-1,1))

open (UNIT=20, FILE='second_derivative_debug')


do i = 1, num_points

write(20,*) (M(i,j), j=1,7)

end do


close(10)
close(20)
end subroutine

!
! SUBROTUINES TO GENERATE COEFFICIENTS FOR THE POLYNOMYAL
!
subroutine polynomials_coefficients(M,num_points,h,z,poly_coeff,poly_coeff_deriv,poly_coeff_deriv2)

real*8, intent(in) :: M(num_points,7)
real*8, intent(in) :: h
integer, intent(in) :: z
integer, intent(in) :: num_points
real*8, intent(out) :: poly_coeff(6,2)
real*8, intent(out) :: poly_coeff_deriv(5,2)
real*8, intent(out) :: poly_coeff_deriv2(4,2)
integer :: i
real*8, dimension(2) :: p0, p1, dp0, dp1, d2p0, d2p1

!! NO NEED TO PASS THE WHOLE MATRIX

!!!CHECK MATRIX THING. MAYBE CHANGE THE STORAGE OF THE INITIAL VALUES. INSTEAD 
! OF A MATRIX USE DIFERENT VECTORS
!
! JUST RENAMING VARIABLES FOR BETTER UNDERSTANDING AND CODING
!

do i = 1, 2

p0(i) = M(z,i+1)
p1(i) = M(z+1,i+1)
dp0(i) = M(z,i+3)
dp1(i) = M(z+1,i+3)
d2p0(i) = M(z,i+5)
d2p1(i) = M(z+1,i+5)

end do



do i = 1,2
!!REORDE 1 AND 2
   Poly_coeff(1,i) = -3.d0*dp0(i)*h-3.d0*dp1(i)*h-(d2p0(i)*h**2)/2.d0+(d2p1(i)*h**2)/2.d0-6.d0*p0(i)+6.d0*p1(i)
   Poly_coeff(2,i) = 8*dp0(i)*h+7*dp1(i)*h+(3.d0*d2p0(i)*h**2)/2.d0-d2p1(i)*h**2+15.d0*p0(i)-15.d0*p1(i)
   Poly_coeff(3,i) = -6.d0*dp0(i)*h-4*dp1(i)*h-(3*d2p0(i)*h**2)/2.d0+(d2p1(i)*h**2)/2.d0-10*p0(i)+10*p1(i)
   Poly_coeff(4,i) = (d2p0(i)*h**2)/2.d0
   Poly_coeff(5,i) = dp0(i)*h
   Poly_coeff(6,i) = p0(i)

   Poly_coeff_deriv (1,i) = 5.d0 * Poly_coeff(1,i)/(M(z+1,1)-M(z,1))
   Poly_coeff_deriv (2,i) = 4.d0 * Poly_coeff(2,i)/(M(z+1,1)-M(z,1))
   Poly_coeff_deriv (3,i) = 3.d0 * Poly_coeff(3,i)/(M(z+1,1)-M(z,1))
   Poly_coeff_deriv (4,i) = 2.d0 * Poly_coeff(4,i)/(M(z+1,1)-M(z,1))
   Poly_coeff_deriv (5,i) = Poly_coeff(5,i)/(M(z+1,1)-M(z,1))

   Poly_coeff_deriv2 (1,i) = 4.d0 * Poly_coeff_deriv(1,i)/(M(z+1,1)-M(z,1))
   Poly_coeff_deriv2 (2,i) = 3.d0 * Poly_coeff_deriv(2,i)/(M(z+1,1)-M(z,1))
   Poly_coeff_deriv2 (3,i) = 2.d0 * Poly_coeff_deriv(3,i)/(M(z+1,1)-M(z,1))
   Poly_coeff_deriv2 (4,i) = Poly_coeff_deriv(4,i)/(M(z+1,1)-M(z,1))

end do
end subroutine

!
! SUBROUTINES TO GENERATE ADDITIONAL POINTS (EQUALLY SPACED) FOR THE POLYNOMYAL
!

subroutine polynomial_evaluation (M1,M2,num_points,Poly_coeff,Poly_coeff_deriv,Poly_coeff_deriv2,np_poly)

real*8, intent(in) :: M1, M2
real*8, intent(in) :: Poly_coeff(6,2), Poly_coeff_deriv(5,2), Poly_coeff_deriv2(4,2) 
integer, intent(in) :: num_points
real*8 :: eta_inter, t
real*8 :: poly_real, poly_imag, poly_real_deriv, poly_imag_deriv, poly_real_deriv2, poly_imag_deriv2
integer :: j, np_poly
real*8 :: coeffs(9)


do j = 1,np_poly
       eta_inter = M1*REAL(np_poly-j,8)/(REAL(np_poly-1,8)) + M2*(REAL(j-1,8)/(REAL(np_poly-1,8)))
       t = (eta_inter-M1)/(M2-M1)
       poly_real = interpol(t,Poly_coeff(1,1), Poly_coeff(2,1), Poly_coeff(3,1), Poly_coeff(4,1), Poly_coeff(5,1), Poly_coeff(6,1))
       poly_imag = interpol(t,Poly_coeff(1,2), Poly_coeff(2,2), Poly_coeff(3,2), Poly_coeff(4,2), Poly_coeff(5,2), Poly_coeff(6,2))

       poly_real_deriv = interpol_deriv(t,Poly_coeff_deriv(1,1), Poly_coeff_deriv(2,1), Poly_coeff_deriv(3,1),&
       Poly_coeff_deriv(4,1), Poly_coeff_deriv(5,1))
       poly_imag_deriv = interpol_deriv(t,Poly_coeff_deriv(1,2), Poly_coeff_deriv(2,2), Poly_coeff_deriv(3,2),&
       Poly_coeff_deriv(4,2), Poly_coeff_deriv(5,2))

       poly_real_deriv2 = interpol_deriv2(t,Poly_coeff_deriv2(1,1), Poly_coeff_deriv2(2,1), Poly_coeff_deriv2(3,1),&
       Poly_coeff_deriv2(4,1))
       poly_imag_deriv2 = interpol_deriv2(t,Poly_coeff_deriv2(1,2), Poly_coeff_deriv2(2,2), Poly_coeff_deriv2(3,2),&
       Poly_coeff_deriv2(4,2))

       call velocity_deriv(eta_inter,t,M1,M2,Poly_coeff_deriv(1,1), Poly_coeff_deriv(2,1), Poly_coeff_deriv(3,1),&
       Poly_coeff_deriv(4,1), Poly_coeff_deriv(5,1),Poly_coeff_deriv(1,2), Poly_coeff_deriv(2,2), Poly_coeff_deriv(3,2),&
       Poly_coeff_deriv(4,2), Poly_coeff_deriv(5,2),coeffs)
       write (30,*) eta_inter,poly_real,poly_imag,poly_real_deriv,poly_imag_deriv,&
                    poly_real_deriv2, poly_imag_deriv2,&
                    eta_inter*SQRT((poly_real_deriv)**2+(poly_imag_deriv)**2)
end do
end subroutine

!subroutine polynomial_evaluation_U (M1,M2,num_points,Poly_coeff,Poly_coeff_deriv,np_poly)
!
!real*8, intent(in) :: M1, M2
!real*8, intent(in) :: Poly_coeff(4,2), Poly_coeff_deriv (3,2)
!integer, intent(in) :: num_points
!real*8 :: eta_inter
!real*8 :: poly_real, poly_imag, poly_real_deriv, poly_imag_deriv
!real*8 :: poly_real_U, poly_imag_U, poly_real_deriv_U, poly_imag_deriv_U
!integer :: j, np_poly
!
!
!do j = 1,np_poly
!       eta_inter = M1*REAL(np_poly-j,8)/(REAL(np_poly-1,8)) + M2*(REAL(j-1,8)/(REAL(np_poly-1,8)))
!
!       poly_real = interpol(eta_inter-M1,Poly_coeff(1,1), Poly_coeff(2,1), Poly_coeff(3,1), Poly_coeff(4,1))
!       poly_imag = interpol(eta_inter-M1,Poly_coeff(1,2), Poly_coeff(2,2), Poly_coeff(3,2), Poly_coeff(4,2))
!
!       poly_real_deriv = interpol_deriv(eta_inter-M1,Poly_coeff_deriv(1,1), Poly_coeff_deriv(2,1), Poly_coeff_deriv(3,1))
!       poly_imag_deriv = interpol_deriv(eta_inter-M1,Poly_coeff_deriv(1,2), Poly_coeff_deriv(2,2), Poly_coeff_deriv(3,2))
!
!       poly_real_U = interpol_U(eta_inter, poly_real, poly_real_deriv) 
!       poly_imag_U = interpol_U(eta_inter, poly_imag, poly_imag_deriv)
!!
!       poly_real_deriv_U = interpol_deriv_U(eta_inter, Poly_coeff(1,1), Poly_coeff(2,1), M1) 
!       poly_imag_deriv_U = interpol_deriv_U(eta_inter, Poly_coeff(1,2), Poly_coeff(2,2), M1)
!
!       write (30,*) eta_inter,poly_real_U,poly_imag_U,poly_real_deriv_U,poly_imag_deriv_U,&
!                   eta_inter*SQRT((poly_real_deriv_U)**2+(poly_imag_deriv_U)**2)
!end do
!end subroutine
!
!
subroutine fitting (M,num_points)

real*8, intent(in) :: M(num_points, 7)
integer, intent(in) :: num_points
real*8 :: Poly_coeff(6,2),Poly_coeff_deriv(5,2), Poly_coeff_deriv2(4,2)
integer :: z, np_poly
real*8 :: eta_step


write (6,*) 'Number of points of the polynomyal'
read (5,*) np_poly
open (unit=30, file='Polynomial_points')
do z = 1, num_points-1
       eta_step = M(z+1,1)-M(z,1)
       call polynomials_coefficients (M,num_points,eta_step,z,Poly_coeff,poly_coeff_deriv,poly_coeff_deriv2)
       call polynomial_evaluation (M(z,1),M(z+1,1),num_points,Poly_coeff,poly_coeff_deriv,poly_coeff_deriv2,np_poly)
end do
end subroutine

!subroutine fitting_U (M,num_points)
!
!real*8, intent(in) :: M(num_points, 5)
!integer, intent(in) :: num_points
!real*8 :: Poly_coeff(6,2),Poly_coeff_deriv(5,2)
!integer :: z, np_poly
!real*8 :: eta_step
!
!write (6,*) 'Number of points of the polynomyal'
!read (5,*) np_poly
!open (unit=30, file='Polynomial_points')
!do z = 1, num_points-1
!       eta_step = M(z+1,1)-M(z,1)
!       call polynomials_coefficients (M,num_points,eta_step,z,Poly_coeff,poly_coeff_deriv)
!       call polynomial_evaluation_U (M(z,1),M(z+1,1),num_points,Poly_coeff,poly_coeff_deriv,np_poly)
!end do
!end subroutine
!
real*8 function interpol (eta, c1, c2, c3 ,c4, c5, c6) result (poly)

real*8, intent(in) :: c1, c2, c3, c4, c5, c6
real*8, intent(in) :: eta

poly = c1*eta**5 + c2*eta**4 + c3*eta**3 + c4*eta**2 + c5*eta + c6
end function


real*8 function interpol_deriv (eta, c1, c2, c3, c4, c5) result (poly_deriv)

real*8, intent(in) :: c1, c2, c3, c4, c5
real*8, intent(in) :: eta

poly_deriv = c1*eta**4 + c2*eta**3 + c3*eta**2 + c4*eta + c5
end function


real*8 function interpol_deriv2 (eta, c1, c2, c3 ,c4) result (poly)

real*8, intent(in) :: c1, c2, c3, c4
real*8, intent(in) :: eta

poly = c1*eta**3 + c2*eta**2 + c3*eta + c4
end function


!
!  MAYBE WRITE INTERPO_U FULL EXRPESION oe WRITE INTERPOL_DERIV_U IN FUNCTION OF ETA, POLY POLY_DERIV, POLY_DERIV_2 
!

!real*8 function interpol_U (eta, poly, poly_deriv) result (poly_U)
!
!real*8, intent(in) :: poly, poly_deriv
!real*8, intent(in) :: eta
!
!poly_U = poly-eta*poly_deriv
!end function
!
!real*8 function interpol_deriv_U (eta, c1, c2, eta1) result (poly_U_deriv)
!
!real*8, intent(in) :: c1, c2
!real*8, intent(in) :: eta, eta1
!
!poly_U_deriv = -(6.d0*c1*eta**2+(2.d0*c2-6.d0*eta1*c1)*eta)
!end function
!

!
!SUBROUTINES USED FOR THE CALCUALTION OF THE MINIMA
!
!
!subroutine minimum_find (M,num_points)
!
!real*8, intent(in) :: M(num_points, 5)
!integer, intent(in) :: num_points
!complex(dp) :: zeros(4)
!real*8 :: coeffs(5)
!real*8 :: Poly_coeff(4,2),Poly_coeff_deriv(3,2)
!real*8 :: eta_step
!real*8 :: minimum, initial_point
!integer :: z, i
!
!
!do z = 1, num_points-1
!       eta_step = M(z+1,1)-M(z,1)
!       call polynomials_coefficients (M,num_points,eta_step,z,Poly_coeff,poly_coeff_deriv)
!       call velocity_deriv (M(z,1),Poly_coeff_deriv(1,1), Poly_coeff_deriv(2,1), Poly_coeff_deriv(3,1),&
!                            Poly_coeff_deriv(1,2), Poly_coeff_deriv(2,2), Poly_coeff_deriv(3,2),coeffs)
!       call QuarticRoots(Coeffs,zeros)
!       do i =1,4
!          if (AIMAG(zeros(i)) .eq. 0.d0 .and. DBLE(zeros(i)) .gt. M(z,1) .and. DBLE(zeros(i)) .lt. M(z+1,1)) then
!             minimum = DBLE(zeros(i))
!             if (velocity_deriv2(coeffs,DBLE(zeros(i))) .gt. 0.d0) then
!                initial_point = M(z,1)
!                call polynomial_evaluation_minima (minimum,initial_point,Poly_coeff,poly_coeff_deriv)
!             end if
!          end if
!       end do
!end do
!end subroutine
!
!subroutine minimum_find_U (M,num_points)
!
!real*8, intent(in) :: M(num_points, 5)
!integer, intent(in) :: num_points
!complex(dp) :: zeros(4)
!real*8 :: coeffs(5)
!real*8 :: Poly_coeff(4,2),Poly_coeff_deriv(3,2)
!real*8 :: eta_step
!real*8 :: minimum, initial_point
!integer :: z, i
!
!
!do z = 1, num_points-1
!       eta_step = M(z+1,1)-M(z,1)
!       call polynomials_coefficients (M,num_points,eta_step,z,Poly_coeff,poly_coeff_deriv)
!       call velocity_deriv_U (Poly_coeff(1,1),Poly_coeff(2,1),Poly_coeff(1,2),Poly_coeff(2,2),M(z,1),coeffs)
!       call QuarticRoots(Coeffs,zeros)
!       do i =1,4
!          if (AIMAG(zeros(i)) .eq. 0.d0 .and. DBLE(zeros(i)) .gt. M(z,1) .and. DBLE(zeros(i)) .lt. M(z+1,1)) then
!             minimum = DBLE(zeros(i))
!             if (velocity_deriv2_U(coeffs,DBLE(zeros(i))) .gt. 0.d0) then
!                initial_point = M(z,1)
!                call polynomial_evaluation_minima_U (minimum,initial_point,Poly_coeff,poly_coeff_deriv)
!             end if
!          end if
!       end do
!end do
!end subroutine
!
!
!!
!! a,b,c and d,e,f are the coefficents of the derivative of the real and the complex polynomial
!!
subroutine velocity_deriv (eta,t,eta1,eta2,a,b,c,d,e,f,g,h,i,j,coeffs)
!COEFFICIENTS FOR THE DERIVATIVE OF THE VELOCITY OF TH 0th ORDER ENERGY
!real*8, intent(in) :: eta
real*8, intent(in) :: eta1, eta2
real*8, intent(in) :: eta, t
real*8, intent(in) :: a,b,c,d,e
real*8, intent(in) :: f,g,h,i,j
real*8, intent(out) :: coeffs(9)
real*8 :: step, poly


!
! NEEDS EXTENSIVE CLEANING + d0 TO ALL NUMBERS MULTYPLING/DIVIDING
!
!

step = eta2-eta1

coeffs(1) = e**2 + (d*e*eta1)/(-eta1 + eta2) + (eta1*i*j)/(-eta1 + eta2) + j**2 
coeffs(2) = (2*d*e + (d**2*eta1)/(-eta1 + eta2) + (2*c*e*eta1)/(-eta1 + eta2) - (d*e*eta1)/(-eta1 + eta2) +&
            (d*e*eta2)/(-eta1 + eta2) + (eta1*i**2)/(-eta1 + eta2) + (2*eta1*h*j)/(-eta1 + eta2) +&
            2*i*j -    (eta1*i*j)/(-eta1 + eta2) + (eta2*i*j)/(-eta1 + eta2))
coeffs(3) = (d**2 + 2*c*e + (3*c*d*eta1)/(-eta1 + eta2) - (d**2*eta1)/(-eta1 + eta2) +&
            (3*b*e*eta1)/(-eta1 + eta2) - (2*c*e*eta1)/(-eta1 + eta2) + (d**2*eta2)/(-eta1 + eta2) +&
            (2*c*e*eta2)/(-eta1 + eta2) + (3*eta1*h*i)/(-eta1 + eta2) + i**2 - (eta1*i**2)/(-eta1 + eta2) +&
            (eta2*i**2)/(-eta1 + eta2) + (3*eta1*g*j)/(-eta1 + eta2) + 2*h*j - (2*eta1*h*j)/(-eta1 + eta2) +&
            (2*eta2*h*j)/(-eta1 + eta2))
coeffs(4) = (2*c*d + 2*b*e + (2*c**2*eta1)/(-eta1 + eta2) + (4*b*d*eta1)/(-eta1 + eta2) - (3*c*d*eta1)/(-eta1 + eta2) +&
            (4*a*e*eta1)/(-eta1 + eta2) - (3*b*e*eta1)/(-eta1 + eta2) + (3*c*d*eta2)/(-eta1 + eta2) +&
            (3*b*e*eta2)/(-eta1 + eta2) + (2*eta1*h**2)/(-eta1 + eta2) + (4*eta1*g*i)/(-eta1 + eta2) +&
            2*h*i - (3*eta1*h*i)/(-eta1 + eta2) + (3*eta2*h*i)/(-eta1 + eta2) + (4*eta1*f*j)/(-eta1 + eta2) +&
            2*g*j - (3*eta1*g*j)/(-eta1 + eta2) + (3*eta2*g*j)/(-eta1 + eta2))
coeffs(5) = (c**2 + 2*b*d + 2*a*e + (5*b*c*eta1)/(-eta1 + eta2) - (2*c**2*eta1)/(-eta1 + eta2) +&
            (5*a*d*eta1)/(-eta1 + eta2) - (4*b*d*eta1)/(-eta1 + eta2) - (4*a*e*eta1)/(-eta1 + eta2) +&
            (2*c**2*eta2)/(-eta1 + eta2) + (4*b*d*eta2)/(-eta1 + eta2) + (4*a*e*eta2)/(-eta1 + eta2) +&
            (5*eta1*g*h)/(-eta1 + eta2) + h**2 - (2*eta1*h**2)/(-eta1 + eta2) + (2*eta2*h**2)/(-eta1 + eta2) +&
            (5*eta1*f*i)/(-eta1 + eta2) + 2*g*i - (4*eta1*g*i)/(-eta1 + eta2) + (4*eta2*g*i)/(-eta1 + eta2) +&
            2*f*j - (4*eta1*f*j)/(-eta1 + eta2) + (4*eta2*f*j)/(-eta1 + eta2))
coeffs(6) = (2*b*c + 2*a*d + (3*b**2*eta1)/(-eta1 + eta2) + (6*a*c*eta1)/(-eta1 + eta2) -&
            (5*b*c*eta1)/(-eta1 + eta2) -  (5*a*d*eta1)/(-eta1 + eta2) + (5*b*c*eta2)/(-eta1 + eta2) +&
            (5*a*d*eta2)/(-eta1 + eta2) + (3*eta1*g**2)/(-eta1 + eta2) + (6*eta1*f*h)/(-eta1 + eta2) +&
            2*g*h - (5*eta1*g*h)/(-eta1 + eta2) + (5*eta2*g*h)/(-eta1 + eta2) + 2*f*i - (5*eta1*f*i)/(-eta1 + eta2) +&
            (5*eta2*f*i)/(-eta1 + eta2))
coeffs(7) = (b**2 + 2*a*c + (7*a*b*eta1)/(-eta1 + eta2) - (3*b**2*eta1)/(-eta1 + eta2) - (6*a*c*eta1)/(-eta1 + eta2) +&
            (3*b**2*eta2)/(-eta1 + eta2) + (6*a*c*eta2)/(-eta1 + eta2) + (7*eta1*f*g)/(-eta1 + eta2) + g**2 -&
            (3*eta1*g**2)/(-eta1 + eta2) + (3*eta2*g**2)/(-eta1 + eta2) + 2*f*h - (6*eta1*f*h)/(-eta1 + eta2) +&
            (6*eta2*f*h)/(-eta1 + eta2))
coeffs(8) = (2*a*b + (4*a**2*eta1)/(-eta1 + eta2) - (7*a*b*eta1)/(-eta1 + eta2) + (7*a*b*eta2)/(-eta1 + eta2) +&
            (4*eta1*f**2)/(-eta1 + eta2) + 2*f*g - (7*eta1*f*g)/(-eta1 + eta2) + (7*eta2*f*g)/(-eta1 + eta2))
coeffs(9) = (a**2 - (4*a**2*eta1)/(-eta1 + eta2) + (4*a**2*eta2)/(-eta1 + eta2) + f**2 -&
            (4*eta1*f**2)/(-eta1 + eta2) + (4*eta2*f**2)/(-eta1 + eta2))



!coeffs = coeffs/eta1

poly = coeffs(9)*t**8 + coeffs(8)*t**7 + coeffs(7)*t**6 + coeffs(6)*t**5 + coeffs(5)*t**4 + coeffs(4)*t**3 +&
       coeffs(3)*t**2 + coeffs(2)*t + coeffs(1)

write (50,*) eta,t,poly
end subroutine

!real*8 function velocity_deriv2 (coeffs, eta) result (res)
!!EVALUATION OF THE 2nd DERIVATIVE OF THE VELOCITY OF TH 1st ORDER ENERGY
!real*8, intent(in) :: coeffs(5), eta
!
!res = 4.d0*coeffs(5)*eta**3+3.d0*coeffs(4)*eta**2+2.d0*coeffs(3)*eta+coeffs(2)
!end function
!
!!!
!!! REVISE THIS SUBROUTINE I WROTE 2+ AND IT SHOULD BE 2*
!!!
!subroutine velocity_deriv_U (a,b,c,d,eta_1,coeffs)
!!COEFFICIENTS FOR THE DERIVATIVE OF THE VELOCITY OF TH 1st ORDER ENERGY
!real*8, intent(in) :: a, b, c, d, eta_1
!real*8, intent(out) :: coeffs(5)
!
!! Compute the terms step by step
!coeffs (:) = 0.d0
!
!coeffs(3) = 2.0D0 * b * (2.0D0 * b - 6.0D0 * a * eta_1) - 6.0D0 * a * eta_1 *(2.0D0 * b - 6.0D0*a * eta_1) + &
!           (2.0D0 * b - 6.0D0 * a * eta_1)**2 + 2.0D0 * d * (2.0D0 * d - 6.0D0 * c * eta_1) -&
!            6.0D0 * c * eta_1 * (2.0D0 * d - 6.0D0 * c * eta_1) + (2.0D0 * d - 6.0D0 * c * eta_1)**2
!coeffs(4) = 12.0D0 * a * b + 12.0D0 * c * d - 36.0D0 * a**2 * eta_1 - 36.0D0 * c**2 * eta_1 + &
!            24.0D0 * a * (2.0D0 * b - 6.0D0 * a * eta_1) + 24.0D0 * c * (2.0D0 * d - 6.0D0 * c * eta_1)
!coeffs(5) = (108.0D0 * a**2 + 108.0D0 * c**2)
!
!! Final result
!!poly = coeffs(5)*eta**4 + coeffs(4)*eta**3 + coeffs(3)*eta**2
!end subroutine
!
!real*8 function velocity_deriv2_U (coeffs, eta) result (res)
!!EVALUATION OF THE 2nd DERIVATIVE OF THE VELOCITY OF THE 1st ORDER ENERGY
!real*8, intent(in) :: coeffs(5), eta
!
!res = 4.d0*coeffs(5)*eta**3+3.d0*coeffs(4)*eta**2+2.d0*coeffs(3)*eta
!end function
!
!subroutine polynomial_evaluation_minima (eta_inter,initial_point,Poly_coeff,Poly_coeff_deriv)
!
!real*8, intent(in) :: eta_inter, initial_point
!real*8, intent(in) :: Poly_coeff(4,2), Poly_coeff_deriv (3,2)
!real*8 :: poly_real, poly_imag, poly_real_deriv, poly_imag_deriv
!
!
!       open (UNIT=20, FILE='Minima.dat')
!
!       poly_real = interpol(eta_inter-initial_point,Poly_coeff(1,1), Poly_coeff(2,1), Poly_coeff(3,1), Poly_coeff(4,1))
!       poly_imag = interpol(eta_inter-initial_point,Poly_coeff(1,2), Poly_coeff(2,2), Poly_coeff(3,2), Poly_coeff(4,2))
!
!       poly_real_deriv = interpol_deriv(eta_inter-initial_point,Poly_coeff_deriv(1,1), Poly_coeff_deriv(2,1), Poly_coeff_deriv(3,1))
!       poly_imag_deriv = interpol_deriv(eta_inter-initial_point,Poly_coeff_deriv(1,2), Poly_coeff_deriv(2,2), Poly_coeff_deriv(3,2))
!
!       write (20,*) eta_inter, poly_real, poly_imag, poly_real_deriv, poly_imag_deriv,&
!                    eta_inter*SQRT(poly_real_deriv**2+poly_imag_deriv**2),&
!                    poly_real-eta_inter*poly_real_deriv, poly_imag-eta_inter*poly_imag_deriv
!
!end subroutine
!
!subroutine polynomial_evaluation_minima_U (eta_inter,initial_point,Poly_coeff,Poly_coeff_deriv)
!
!real*8, intent(in) :: eta_inter, initial_point
!real*8, intent(in) :: Poly_coeff(4,2), Poly_coeff_deriv (3,2)
!real*8 :: poly_real, poly_imag, poly_real_deriv, poly_imag_deriv
!real*8 :: poly_real_U, poly_imag_U, poly_real_deriv_U, poly_imag_deriv_U
!
!       open (UNIT=20, FILE='Minima.dat')
!
!       poly_real = interpol(eta_inter-initial_point,Poly_coeff(1,1), Poly_coeff(2,1), Poly_coeff(3,1), Poly_coeff(4,1))
!       poly_imag = interpol(eta_inter-initial_point,Poly_coeff(1,2), Poly_coeff(2,2), Poly_coeff(3,2), Poly_coeff(4,2))
!
!       poly_real_deriv = interpol_deriv(eta_inter-initial_point,Poly_coeff_deriv(1,1), Poly_coeff_deriv(2,1), Poly_coeff_deriv(3,1))
!       poly_imag_deriv = interpol_deriv(eta_inter-initial_point,Poly_coeff_deriv(1,2), Poly_coeff_deriv(2,2), Poly_coeff_deriv(3,2))
!
!       poly_real_U = interpol_U(eta_inter, poly_real, poly_real_deriv) 
!       poly_imag_U = interpol_U(eta_inter, poly_imag, poly_imag_deriv)
!!
!       poly_real_deriv_U = interpol_deriv_U(eta_inter, Poly_coeff(1,1), Poly_coeff(2,1), initial_point) 
!       poly_imag_deriv_U = interpol_deriv_U(eta_inter, Poly_coeff(1,2), Poly_coeff(2,2), initial_point)
!
!       write (20,*) eta_inter,poly_real_U,poly_imag_U,poly_real_deriv_U,poly_imag_deriv_U,&
!                   eta_inter*SQRT((poly_real_deriv_U)**2+(poly_imag_deriv_U)**2)
!
!end subroutine
!
!
!
!subroutine evaluation_at_eta_value(M,num_points)
!
!
!real*8, intent(in) :: M(num_points,5)
!integer, intent(in) :: num_points
!real*8 :: evaluated_eta, eta_step
!real*8 :: Poly_coeff(4,2), poly_coeff_deriv(3,2)
!integer :: z, segment_saved
!
!
!write(*,*) 'Specify the eta value'
!read(*,*) evaluated_eta
!
!
!!
!! IN CASE YOU WANT TO READ THE FILE CONTAINING INFORMATION FOR N2
!!
!
!!
!!  ADD A WARNING IN CASE THE VALUE OF EVALUATED_ETA IS OUTSIDE THE INTERVAL OF GIVEN POINTS
!!
!do z = 1, num_points-1
!       if (evaluated_eta .ge. M(z,1) .and. evaluated_eta .le. M(z+1,1)) then
!          eta_step = M(z+1,1)-M(z,1)
!          segment_saved = z
!          exit 
!       end if
!end do
!
!call polynomials_coefficients (M,num_points,eta_step,segment_saved,Poly_coeff,poly_coeff_deriv)
!call polynomial_evaluation_minima (evaluated_eta,M(segment_saved,1),Poly_coeff,Poly_coeff_deriv)
!end subroutine
!
end module

!
!
!

program Spline

Use spline_module

implicit none


real*8, allocatable :: M(:,:)
integer :: num_points
integer :: option2
character(len=50) :: option


call Reading (M,num_points)


open (UNIT=50, FILE='VELOCITY_DERIV')
write (6,*) 'Write Polynomial_fitting or Minima_find or Neutral_Interpolation in function of option desired'
read(5,*) option
if (option .eq. 'Polynomial_fitting') then
    write (6,*) 'Write 0 or 1 in function if the value searched corresponds to the 0th or 1st order energy'
    read(5,*) option2
    if (option2 .eq. 0) then
       call fitting (M,num_points)
    else if (option2 .eq. 1) then
!       call fitting_U (M,num_points)
    else
       write (*,*) 'Bad option, try again'
    end if
!else if (option .eq. 'Minima_find') then
!    write (6,*) 'Write 0 or 1 in function if the value searched corresponds to the 0th or 1st order energy'
!    read(5,*) option2
!    if (option2 .eq. 0) then    
!       call minimum_find (M,num_points)
!    else if (option2 .eq. 1) then
!       call minimum_find_U (M,num_points)
!    else
!       write (*,*) 'Bad option, try again'
!    end if
!else if (option .eq. 'Neutral_Interpolation') then
!    call evaluation_at_eta_value(M,num_points)
!else
    write (*,*) 'Bad option chosen, try again'
end if



!


end program

