module spline_module

use PolynomialRoots 

implicit none 


integer, parameter ::dp = kind(0.d0)


contains

!Read the points in file and storing data in a matrix (supossing file contains:
! eta_value real_energy imag_ener real_derivatvie imag_deriv
!
subroutine Reading (M,num_points)

real*8, allocatable, intent(out) :: M (:,:) 
integer, intent(out) :: num_points
integer :: i, ios
character(len=30) :: file_name
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
!!!
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


subroutine polynomial_evaluation (M,num_points,z,Poly_coeff,Poly_coeff_deriv,np_poly)

real*8, intent(in) :: M(num_points,5)
real*8, intent(in) :: Poly_coeff(4,2), Poly_coeff_deriv (3,2)
integer, intent(in) :: z, num_points
real*8 :: eta_inter
real*8 :: poly_real, poly_imag, poly_real_deriv, poly_imag_deriv, poly_real_deriv2, poly_imag_deriv2
integer :: j, np_poly
real*8 :: coeffs(5),poly


!!DREAL or dp as in the complex CHANGE
do j = 1,np_poly
       eta_inter = M(z,1)*(REAL(np_poly-1,8)-REAL(j-1,8))/(REAL(np_poly-1,8)) + M(z+1,1)*(REAL(j-1,8)/(REAL(np_poly-1,8)))
       
       poly_real = interpol(eta_inter-M(z,1),Poly_coeff(1,1), Poly_coeff(2,1), Poly_coeff(3,1), Poly_coeff(4,1))
       poly_imag = interpol(eta_inter-M(z,1),Poly_coeff(1,2), Poly_coeff(2,2), Poly_coeff(3,2), Poly_coeff(4,2)) 
         
       poly_real_deriv = interpol_deriv(eta_inter-M(z,1),Poly_coeff_deriv(1,1), Poly_coeff_deriv(2,1), Poly_coeff_deriv(3,1))
       poly_imag_deriv = interpol_deriv(eta_inter-M(z,1),Poly_coeff_deriv(1,2), Poly_coeff_deriv(2,2), Poly_coeff_deriv(3,2))

       call minimum_find_U (eta_inter,Poly_coeff(1,1),Poly_coeff(2,1),Poly_coeff(1,2),Poly_coeff(2,2),M(z,1),coeffs,poly)
!       poly_real_deriv2 = interpol_deriv2(eta_inter-M(z,1),2.d0*Poly_coeff_deriv(1,1),& 
        !                  Poly_coeff_deriv(2,1), eta_inter)
!       poly_imag_deriv2 = interpol_deriv2(eta_inter-M(z,1),2.d0*Poly_coeff_deriv(1,2),& 
        !                  Poly_coeff_deriv(2,2), eta_inter)

        write (*,*) eta_inter,interpol_deriv_U(eta_inter,Poly_coeff(1,1),Poly_coeff(2,1),M(z,1)),&
                    interpol_deriv_U(eta_inter,Poly_coeff(1,2),Poly_coeff(2,2),M(z,1)),&
                    eta_inter*SQRT(interpol_deriv_U(eta_inter,Poly_coeff(1,1),Poly_coeff(2,1),M(z,1))**2+&
                                   interpol_deriv_U(eta_inter,Poly_coeff(1,2),Poly_coeff(2,2),M(z,1))**2),&
                    poly 
    !   write (*,*) eta_inter, poly_real, poly_imag, poly_real_deriv, poly_imag_deriv
!                   poly_real-eta_inter*poly_real_deriv, poly_imag-eta_inter*poly_imag_deriv,&
!                   eta_inter*dsqrt(poly_real_deriv2**2+poly_imag_deriv2**2)
                   !eta_inter*dsqrt(poly_real_deriv**2+poly_imag_deriv**2)
end do

end subroutine



subroutine minimum_find_U (eta,a,b,c,d,eta_1,coeffs,poly)

real*8, intent(in) :: eta, a, b, c, d, eta_1
real*8, intent(out) :: coeffs(5), poly

! Compute the terms step by step
coeffs (:) = 0.d0
coeffs(5) = (108.0D0 * a**2 + 108.0D0 * c**2)

coeffs(4) = 12.0D0 * a + 6.0D0 * a * b + 12.0D0 * c + 6.0D0 * c * d - 36.0D0 * a**2 * eta_1 - 36.0D0 * c**2 * eta_1 + &
    24.0D0 * a * (2.0D0 * b - 6.0D0 * a * eta_1) + 24.0D0 * c * (2.0D0 * d - 6.0D0 * c * eta_1)

coeffs(3) = 2.0D0 * (2.0D0 * b - 6.0D0 * a * eta_1) + b * (2.0D0 * b - 6.0D0 * a * eta_1) - 6.0D0 * a * eta_1 *&
        (2.0D0 * b - 6.0D0*a * eta_1) + &
    (2.0D0 * b - 6.0D0 * a * eta_1)**2 + 2.0D0 * (2.0D0 * d - 6.0D0 * c * eta_1) + d * (2.0D0 * d - 6.0D0 * c * eta_1) - &
    6.0D0 * c * eta_1 * (2.0D0 * d - 6.0D0 * c * eta_1) + (2.0D0 * d - 6.0D0 * c * eta_1)**2

! Final result
poly = coeffs(5)*eta**4 + coeffs(4)*eta**3 + coeffs(3)*eta**2

!print *, "The result is:", result


end subroutine





real*8 function interpol_deriv_U (eta, c1, c2, eta1) result (poly)

real*8, intent(in) :: c1, c2
real*8, intent(in) :: eta, eta1

poly = 6.d0*c1*eta**2+(2*c2-6.d0*eta1*c1)*eta

end function


subroutine polynomial_evaluation_minima (point_min,initial_point,Poly_coeff,Poly_coeff_deriv)

real*8, intent(in) :: point_min, initial_point
real*8, intent(in) :: Poly_coeff(4,2), Poly_coeff_deriv (3,2)
real*8 :: eta_inter
real*8 :: poly_real, poly_imag, poly_real_deriv, poly_imag_deriv



       eta_inter = point_min

       poly_real = interpol(eta_inter-initial_point,Poly_coeff(1,1), Poly_coeff(2,1), Poly_coeff(3,1), Poly_coeff(4,1))
       poly_imag = interpol(eta_inter-initial_point,Poly_coeff(1,2), Poly_coeff(2,2), Poly_coeff(3,2), Poly_coeff(4,2))

       poly_real_deriv = interpol_deriv(eta_inter-initial_point,Poly_coeff_deriv(1,1), Poly_coeff_deriv(2,1), Poly_coeff_deriv(3,1))
       poly_imag_deriv = interpol_deriv(eta_inter-initial_point,Poly_coeff_deriv(1,2), Poly_coeff_deriv(2,2), Poly_coeff_deriv(3,2))

       write (*,*) eta_inter, poly_real, poly_imag, poly_real_deriv, poly_imag_deriv,&
                   poly_real-eta_inter*poly_real_deriv, poly_imag-eta_inter*poly_imag_deriv
                   !eta_inter*dsqrt(poly_real_deriv**2+poly_imag_deriv**2)

end subroutine





subroutine fitting (M,num_points)

real*8, intent(in) :: M(num_points, 5)
integer, intent(in) :: num_points
real*8 :: Poly_coeff(4,2),Poly_coeff_deriv(3,2)
integer :: z, np_poly
real*8 :: eta_step

write (6,*) 'Number of points of the polynomyal'
read (5,*) np_poly
do z = 1, num_points-1
       eta_step = M(z+1,1)-M(z,1)
       call polynomials_coefficients (M,num_points,eta_step,z,Poly_coeff,poly_coeff_deriv)

       call polynomial_evaluation (M,num_points,z,Poly_coeff,poly_coeff_deriv,np_poly)
end do

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

do z = 1, num_points-1
       eta_step = M(z+1,1)-M(z,1)
       call polynomials_coefficients (M,num_points,eta_step,z,Poly_coeff,poly_coeff_deriv)
!do j = 0,40
 !      eta = M(z,1)+eta_step*0.025*j
 !      call velocity_deriv (eta,M(z,1),Poly_coeff_deriv(1,1), Poly_coeff_deriv(2,1), Poly_coeff_deriv(3,1),&
!                           Poly_coeff_deriv(1,2), Poly_coeff_deriv(2,2), Poly_coeff_deriv(3,2),coeffs,poly)
  !     write (*,*) 'Poly', eta, poly
!end do

       call minimum_find_U (1.d0,Poly_coeff(1,1),Poly_coeff(2,1),Poly_coeff(1,2),Poly_coeff(2,2),M(z,1),coeffs,poly)

       call QuarticRoots(Coeffs,zeros)
       do i =1,4
          if (AIMAG(zeros(i)) .eq. 0.d0 .and. DBLE(zeros(i)) .gt. M(z,1) .and. DBLE(zeros(i)) .lt. M(z+1,1)) then
             WRITe (*,*) " MINIMA FOUND"
             minimum = DBLE(zeros(i))
             initial_point = M(z,1)
             call polynomial_evaluation_minima (minimum,initial_point,Poly_coeff,poly_coeff_deriv)
          end if
       end do

end do

end subroutine



!!PUT INTENTS
real*8 function interpol (t, c1, c2, c3 ,c4) result (poly)

real*8, intent(in) :: c1, c2, c3, c4
real*8, intent(in) :: t

poly = c1*t**3 + c2*t**2 + c3*t + c4

end function


real*8 function interpol_deriv (t, c1, c2, c3) result (poly)

real*8, intent(in) :: c1, c2, c3
real*8, intent(in) :: t

poly = c1*t**2 + c2*t + c3

end function


!
! a,b,c and d,e,f are the coefficents of the derivative of the real and the complex polynomial
!
subroutine velocity_deriv (eta,eta1,a,b,c,d,e,f,coeffs,poly)

real*8, intent(in) :: eta
real*8, intent(in) :: eta1
real*8, intent(in) :: a,b,c
real*8, intent(in) :: d,e,f
real*8, intent(out) :: coeffs(5)
real*8, intent(out) :: poly

coeffs(5) = 3.d0*a**2+3.d0*d**2
coeffs(4) = 5.d0*a*b+5.d0*d*e-10.d0*a**2*eta1-10.d0*d**2*eta1
coeffs(3) = 2.d0*b**2+4.d0*a*c+2.d0*e**2+4.d0*d*f-12.d0*a*b*eta1-&
         12.d0*d*e*eta1+12.d0*a**2*eta1**2+12.d0*d**2*eta1**2
coeffs(2) = 3.d0*b*c+3.d0*e*f-3.d0*b**2*eta1-6.d0*a*c*eta1-3.d0*e**2*eta1-&
         6.d0*d*f*eta1+9.d0*a*b*eta1**2+9.d0*d*e*eta1**2-&
         6.d0*a**2*eta1**3-6.d0*d**2*eta1**3
coeffs(1) = c**2+f**2-2.d0*b*c*eta1-2.d0*e*f*eta1+b**2*eta1**2+&
         2.d0*a*c*eta1**2+e**2*eta1**2+2.d0*d*f*eta1**2-2.d0*a*b*eta1**3-&
         2.d0*d*e*eta1**3+a**2*eta1**4+d**2*eta1**4


poly = coeffs(5)*eta**4.d0+coeffs(4)*eta**3.d0+coeffs(3)*eta**2.d0+coeffs(2)*eta+coeffs(1)

end subroutine

end module


program Spline

Use spline_module  

implicit none


real*8, allocatable :: M(:,:)
integer :: num_points
integer :: i 
character(len=50) :: option


call Reading (M,num_points)

!do i = 1, num_points
!       write (*,*) M(i,:)
!end do

write (6,*) 'Write Polynomial_fitting or Minima_find in function of option desired'
read(5,*) option
if (option .eq. 'Polynomial_fitting') then
    call fitting (M,num_points)
else if (option .eq. 'Minima_find') then
    call minimum_find (M,num_points) 
else 
    write (*,*) 'Bad option chosen, try again'
end if




end program


