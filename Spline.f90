 module spline_module

use PolynomialRoots

implicit none


integer, parameter ::dp = kind(0.d0)


!
! ADD OTHER FILES FOR POLYNOMYAL :FOR EXAMPLE
!

contains

!Read the points in file and storing data in a matrix (supossing file contains:
! eta_value real_energy imag_ener real_derivatvie imag_deriv
!
!
! MAYBE STORE EACH VARIABLE SEPARETLY IN A VECTOR
!
subroutine Reading (M,num_points)

real*8, allocatable, intent(out) :: M (:,:)
integer, intent(out) :: num_points
integer :: i, ios
character(len=1000) :: file_name
character(len=1) :: foo


write (6,*) 'Insert file name'
read(5,*) file_name

open (UNIT=10, FILE=file_name)

num_points = 0
do
        read (10,*,iostat=ios) foo
        if (ios /= 0) exit
        num_points = num_points + 1
end do


Rewind 10

allocate (M(num_points, 5))
do i = 1, num_points
        read (10, *) M (i, :)
!
!  In the quantum package output, the derivatives are stored as -eta*(dE/deta)
!
        M (i,4) = -M(i,4)/M(i,1)
        M (i,5) = -M(i,5)/M(i,1)
end do




end subroutine

!
! SUBROTUINES TO GENERATE COEFFICIENTS FOR THE POLYNOMYAL
!
subroutine polynomials_coefficients(M,num_points,h,z,poly_coeff,poly_coeff_deriv)

real*8, intent(in) :: M(num_points,5)
real*8, intent(in) :: h
integer, intent(in) :: z
integer, intent(in) :: num_points
real*8, intent(out) :: poly_coeff(4,2)
real*8, intent(out) :: poly_coeff_deriv(3,2)
integer :: i

!! NO NEED TO PASS THE WHOLE MATRIX

!!!CHECK MATRIX THING
do i = 1,2
!!REORDE 1 AND 2
   Poly_coeff(1,i) = (2.d0*M(z,i+1)+h*M(z,i+3)-2.d0*M(z+1,i+1)+h*M(z+1,i+3))/h**3
   Poly_coeff(2,i) =(-3.d0*M(z,i+1)-2.d0*h*M(z,i+3)+3.d0*M(z+1,i+1)-h*M(z+1,i+3))/h**2
   Poly_coeff(3,i) = M(z,i+3)
   Poly_coeff(4,i) = M(z,i+1)

   Poly_coeff_deriv (1,i) = 3.d0 * Poly_coeff(1,i)
   Poly_coeff_deriv (2,i) = 2.d0 * Poly_coeff(2,i)
   Poly_coeff_deriv (3,i) = Poly_coeff(3,i)

end do
end subroutine

!
! SUBROTUINES TO GENERATE ADDITIONAL POINTS FOR THE POLYNOMYAL
!

subroutine polynomial_evaluation (M,num_points,z,Poly_coeff,Poly_coeff_deriv,np_poly)

real*8, intent(in) :: M(num_points,5)
real*8, intent(in) :: Poly_coeff(4,2), Poly_coeff_deriv (3,2)
integer, intent(in) :: z, num_points
real*8 :: eta_inter
real*8 :: poly_real, poly_imag, poly_real_deriv, poly_imag_deriv
integer :: j, np_poly
real*8 :: coeffs(5),poly

do j = 1,np_poly
       eta_inter = M(z,1)*REAL(np_poly-1-(j-1),8)/(REAL(np_poly-1,8)) + M(z+1,1)*(REAL(j-1,8)/(REAL(np_poly-1,8)))

       poly_real = interpol(eta_inter-M(z,1),Poly_coeff(1,1), Poly_coeff(2,1), Poly_coeff(3,1), Poly_coeff(4,1))
       poly_imag = interpol(eta_inter-M(z,1),Poly_coeff(1,2), Poly_coeff(2,2), Poly_coeff(3,2), Poly_coeff(4,2))

       poly_real_deriv = interpol_deriv(eta_inter-M(z,1),Poly_coeff_deriv(1,1), Poly_coeff_deriv(2,1), Poly_coeff_deriv(3,1))
       poly_imag_deriv = interpol_deriv(eta_inter-M(z,1),Poly_coeff_deriv(1,2), Poly_coeff_deriv(2,2), Poly_coeff_deriv(3,2))


       write (*,*) eta_inter,poly_real,poly_imag,poly_real_deriv,poly_imag_deriv,&
                   eta_inter*SQRT((poly_real_deriv)**2+(poly_imag_deriv)**2)
end do
end subroutine

subroutine polynomial_evaluation_U (M,num_points,z,Poly_coeff,Poly_coeff_deriv,np_poly)

real*8, intent(in) :: M(num_points,5)
real*8, intent(in) :: Poly_coeff(4,2), Poly_coeff_deriv (3,2)
integer, intent(in) :: z, num_points
real*8 :: eta_inter
real*8 :: poly_real, poly_imag, poly_real_deriv, poly_imag_deriv
real*8 :: poly_real_U, poly_imag_U, poly_real_deriv_U, poly_imag_deriv_U
integer :: j, np_poly
real*8 :: coeffs(5),poly

do j = 1,np_poly
       eta_inter = M(z,1)*REAL(np_poly-1-(j-1),8)/(REAL(np_poly-1,8)) + M(z+1,1)*(REAL(j-1,8)/(REAL(np_poly-1,8)))

       poly_real = interpol(eta_inter-M(z,1),Poly_coeff(1,1), Poly_coeff(2,1), Poly_coeff(3,1), Poly_coeff(4,1))
       poly_imag = interpol(eta_inter-M(z,1),Poly_coeff(1,2), Poly_coeff(2,2), Poly_coeff(3,2), Poly_coeff(4,2))

       poly_real_deriv = interpol_deriv(eta_inter-M(z,1),Poly_coeff_deriv(1,1), Poly_coeff_deriv(2,1), Poly_coeff_deriv(3,1))
       poly_imag_deriv = interpol_deriv(eta_inter-M(z,1),Poly_coeff_deriv(1,2), Poly_coeff_deriv(2,2), Poly_coeff_deriv(3,2))

       poly_real_U = interpol_U(eta_inter, poly_real, poly_real_deriv) 
       poly_imag_U = interpol_U(eta_inter, poly_imag, poly_imag_deriv)
!
       poly_real_deriv_U = interpol_deriv_U(eta_inter, Poly_coeff(1,1), Poly_coeff(2,1), M(z,1)) 
       poly_imag_deriv_U = interpol_deriv_U(eta_inter, Poly_coeff(1,2), Poly_coeff(2,2), M(z,1))

       write (*,*) eta_inter,poly_real_U,poly_imag_U,poly_real_deriv_U,poly_imag_deriv_U,&
                   eta_inter*SQRT((poly_real_deriv_U)**2+(poly_imag_deriv_U)**2)
end do
end subroutine


subroutine fitting (M,num_points)

real*8, intent(in) :: M(num_points, 5)
integer, intent(in) :: num_points
real*8 :: Poly_coeff(4,2),Poly_coeff_deriv(3,2)
integer :: z, np_poly
real*8 :: eta_step

write (6,*) 'Number of points of the polynomyal'
read (5,*) np_poly
open (unit=30, file='coeffs')
do z = 1, num_points-1
       eta_step = M(z+1,1)-M(z,1)
       call polynomials_coefficients (M,num_points,eta_step,z,Poly_coeff,poly_coeff_deriv)
       call polynomial_evaluation (M,num_points,z,Poly_coeff,poly_coeff_deriv,np_poly)
end do
end subroutine

subroutine fitting_U (M,num_points)

real*8, intent(in) :: M(num_points, 5)
integer, intent(in) :: num_points
real*8 :: Poly_coeff(4,2),Poly_coeff_deriv(3,2)
integer :: z, np_poly
real*8 :: eta_step

write (6,*) 'Number of points of the polynomyal'
read (5,*) np_poly
open (unit=30, file='coeffs')
do z = 1, num_points-1
       eta_step = M(z+1,1)-M(z,1)
       call polynomials_coefficients (M,num_points,eta_step,z,Poly_coeff,poly_coeff_deriv)
       call polynomial_evaluation_U (M,num_points,z,Poly_coeff,poly_coeff_deriv,np_poly)
end do
end subroutine

real*8 function interpol (eta, c1, c2, c3 ,c4) result (poly)

real*8, intent(in) :: c1, c2, c3, c4
real*8, intent(in) :: eta

poly = c1*eta**3 + c2*eta **2 + c3*eta + c4
end function


real*8 function interpol_deriv (eta, c1, c2, c3) result (poly_deriv)

real*8, intent(in) :: c1, c2, c3
real*8, intent(in) :: eta

poly_deriv = c1*eta**2 + c2*eta + c3
end function

real*8 function interpol_U (eta, poly, poly_deriv) result (poly_U)

real*8, intent(in) :: poly, poly_deriv
real*8, intent(in) :: eta

poly_U = poly-eta*poly_deriv
end function

real*8 function interpol_deriv_U (eta, c1, c2, eta1) result (poly_U_deriv)

real*8, intent(in) :: c1, c2
real*8, intent(in) :: eta, eta1

poly_U_deriv = -(6.d0*c1*eta**2+(2.d0*c2-6.d0*eta1*c1)*eta)
end function


!
!SUBROUTINES USED FOR THE CALCUALTION OF THE MINIMA
!



!
! a,b,c and d,e,f are the coefficents of the derivative of the real and the complex polynomial
!
subroutine velocity_deriv (eta1,a,b,c,d,e,f,coeffs)
!COEFFICIENTS FOR THE DERIVATIVE OF THE VELOCITY OF TH 0th ORDER ENERGY
!real*8, intent(in) :: eta
real*8, intent(in) :: eta1
real*8, intent(in) :: a,b,c
real*8, intent(in) :: d,e,f
real*8, intent(out) :: coeffs(5)



coeffs(1) = c**2+f**2-2.d0*b*c*eta1-2.d0*e*f*eta1+b**2*eta1**2+&
         2.d0*a*c*eta1**2+e**2*eta1**2+2.d0*d*f*eta1**2-2.d0*a*b*eta1**3-&
         2.d0*d*e*eta1**3+a**2*eta1**4+d**2*eta1**4
coeffs(2) = 3.d0*b*c+3.d0*e*f-3.d0*b**2*eta1-6.d0*a*c*eta1-3.d0*e**2*eta1-&
         6.d0*d*f*eta1+9.d0*a*b*eta1**2+9.d0*d*e*eta1**2-&
         6.d0*a**2*eta1**3-6.d0*d**2*eta1**3
coeffs(3) = 2.d0*b**2+4.d0*a*c+2.d0*e**2+4.d0*d*f-12.d0*a*b*eta1-&
         12.d0*d*e*eta1+12.d0*a**2*eta1**2+12.d0*d**2*eta1**2
coeffs(4) = 5.d0*a*b+5.d0*d*e-10.d0*a**2*eta1-10.d0*d**2*eta1
coeffs(5) = 3.d0*a**2+3.d0*d**2
end subroutine

real*8 function velocity_deriv2 (coeffs, eta) result (res)
!EVALUATION OF THE 2nd DERIVATIVE OF THE VELOCITY OF TH 0th ORDER ENERGY
real*8, intent(in) :: coeffs(5), eta

res = 4.d0*coeffs(5)*eta**3+3.d0*coeffs(4)*eta**2+2.d0*coeffs(3)*eta+coeffs(2)

end function

!!
!! REVISE THIS SUBROUTINE I WROTE 2+ AND IT SHOULD BE 2*
!!
subroutine velocity_deriv_U (eta,a,b,c,d,eta_1,coeffs,poly)
!COEFFICIENTS FOR THE DERIVATIVE OF THE VELOCITY OF TH 1st ORDER ENERGY
real*8, intent(in) :: eta, a, b, c, d, eta_1
real*8, intent(out) :: coeffs(5), poly

! Compute the terms step by step
coeffs (:) = 0.d0

coeffs(3) = 2.0D0 * (2.0D0 * b - 6.0D0 * a * eta_1) + b * (2.0D0 * b - 6.0D0 * a * eta_1) - 6.0D0 * a * eta_1 *&
        (2.0D0 * b - 6.0D0*a * eta_1) + &
    (2.0D0 * b - 6.0D0 * a * eta_1)**2 + 2.0D0 * (2.0D0 * d - 6.0D0 * c * eta_1) + d * (2.0D0 * d - 6.0D0 * c * eta_1) - &
    6.0D0 * c * eta_1 * (2.0D0 * d - 6.0D0 * c * eta_1) + (2.0D0 * d - 6.0D0 * c * eta_1)**2
coeffs(4) = 12.0D0 * a + 6.0D0 * a * b + 12.0D0 * c + 6.0D0 * c * d - 36.0D0 * a**2 * eta_1 - 36.0D0 * c**2 * eta_1 + &
    24.0D0 * a * (2.0D0 * b - 6.0D0 * a * eta_1) + 24.0D0 * c * (2.0D0 * d - 6.0D0 * c * eta_1)
coeffs(5) = (108.0D0 * a**2 + 108.0D0 * c**2)

! Final result
poly = coeffs(5)*eta**4 + coeffs(4)*eta**3 + coeffs(3)*eta**2

!print *, "The result is:", result


end subroutine


subroutine minimum_find (M,num_points)


real*8, intent(in) :: M(num_points, 5)
integer, intent(in) :: num_points
complex(dp) :: zeros(4)
real*8 :: coeffs(5)
real*8 :: Poly_coeff(4,2),Poly_coeff_deriv(3,2)
real*8 :: eta_step
real*8 :: minimum, initial_point
integer :: z, i,j
real*8 :: eta, poly

!
do z = 1, num_points-1
       eta_step = M(z+1,1)-M(z,1)
       call polynomials_coefficients (M,num_points,eta_step,z,Poly_coeff,poly_coeff_deriv)

!
!    EVALUATION OF THE DERIVATIVE IN ANOTHER SUBROTUINE
       call velocity_deriv (M(z,1),Poly_coeff_deriv(1,1), Poly_coeff_deriv(2,1), Poly_coeff_deriv(3,1),&
                          Poly_coeff_deriv(1,2), Poly_coeff_deriv(2,2), Poly_coeff_deriv(3,2),coeffs)

!       call velocity_deriv_U (1.d0,Poly_coeff(1,1),Poly_coeff(2,1),Poly_coeff(1,2),Poly_coeff(2,2),M(z,1),coeffs,poly)
 
        call QuarticRoots(Coeffs,zeros)
!       write (*,*) zeros
       do i =1,4
          if (AIMAG(zeros(i)) .eq. 0.d0 .and. DBLE(zeros(i)) .gt. M(z,1) .and. DBLE(zeros(i)) .lt. M(z+1,1)) then
!             WRITe (*,*) " MINIMA FOUND"
             minimum = DBLE(zeros(i))
             if (velocity_deriv2(coeffs,DBLE(zeros(i))) .gt. 0.d0) then
                initial_point = M(z,1)
                call polynomial_evaluation_minima (minimum,initial_point,Poly_coeff,poly_coeff_deriv)
             end if
          end if
       end do

end do

end subroutine

subroutine polynomial_evaluation_minima (eta_inter,initial_point,Poly_coeff,Poly_coeff_deriv)

real*8, intent(in) :: eta_inter, initial_point
real*8, intent(in) :: Poly_coeff(4,2), Poly_coeff_deriv (3,2)
real*8 :: poly_real, poly_imag, poly_real_deriv, poly_imag_deriv


       open (UNIT=20, FILE='Minima.dat')

       poly_real = interpol(eta_inter-initial_point,Poly_coeff(1,1), Poly_coeff(2,1), Poly_coeff(3,1), Poly_coeff(4,1))
       poly_imag = interpol(eta_inter-initial_point,Poly_coeff(1,2), Poly_coeff(2,2), Poly_coeff(3,2), Poly_coeff(4,2))

       poly_real_deriv = interpol_deriv(eta_inter-initial_point,Poly_coeff_deriv(1,1), Poly_coeff_deriv(2,1), Poly_coeff_deriv(3,1))
       poly_imag_deriv = interpol_deriv(eta_inter-initial_point,Poly_coeff_deriv(1,2), Poly_coeff_deriv(2,2), Poly_coeff_deriv(3,2))

       write (20,*) eta_inter, poly_real, poly_imag, poly_real_deriv, poly_imag_deriv,&
                    eta_inter*SQRT(poly_real_deriv**2+poly_imag_deriv**2),&
                    poly_real-eta_inter*poly_real_deriv, poly_imag-eta_inter*poly_imag_deriv

end subroutine



subroutine evaluation_at_eta_value(M,num_points)


real*8, intent(in) :: M(num_points,5)
integer, intent(in) :: num_points
real*8 :: evaluated_eta, eta_step
real*8 :: Poly_coeff(4,2), poly_coeff_deriv(3,2)
integer :: z, segment_saved


write(*,*) 'Specify the eta value'
read(*,*) evaluated_eta


!
! IN CASE YOU WANT TO READ THE FILE CONTAINING INFORMATION FOR N2
!

!
! CHANGE M TO ONLY PASS  GIVEM_VALUES?


!
!  ADD A WARNING IN CASE THE VALUE OF EVALUATED_ETA IS OUTSIDE THE INTERVAL OF GIVEN POINTS
!
do z = 1, num_points-1
       if (evaluated_eta .ge. M(z,1) .and. evaluated_eta .le. M(z+1,1)) then
          eta_step = M(z+1,1)-M(z,1)
          segment_saved = z
          exit 
       end if
end do

call polynomials_coefficients (M,num_points,eta_step,segment_saved,Poly_coeff,poly_coeff_deriv)
call polynomial_evaluation_minima (evaluated_eta,M(segment_saved,1),Poly_coeff,Poly_coeff_deriv)



end subroutine


end module


program Spline

Use spline_module

implicit none


real*8, allocatable :: M(:,:)
integer :: num_points
integer :: i
integer :: option2
character(len=50) :: option


call Reading (M,num_points)


write (6,*) 'Write Polynomial_fitting or Minima_find or Neutral_Interpolation in function of option desired'
read(5,*) option
if (option .eq. 'Polynomial_fitting') then
    write (6,*) 'Write 0 or 1 in function if the value searched corresponds to th 0th or 1st order energy'
    read(5,*) option2
    if (option2 .eq. 0) then
       call fitting (M,num_points)
    else if (option2 .eq. 1) then
       call fitting_U (M,num_points)
    else
       write (*,*) 'Bad option, try again'
    end if
else if (option .eq. 'Minima_find') then
    call minimum_find (M,num_points)
else if (option .eq. 'Neutral_Interpolation') then
    call evaluation_at_eta_value(M,num_points)
else
    write (*,*) 'Bad option chosen, try again'
end if



!


end program

