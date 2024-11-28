module spline_module

implicit none 

contains

!Read the points in file and storing data in a matrix (supossing file contains:
! eta_value real_energy imag_ener real_derivatvie iamg_deriv
!Imaginary part should be stored as complex?
!
subroutine Reading (M,num_points)

real*8, allocatable, intent(out) :: M (:,:) 
integer, intent(out) :: num_points
integer :: i, ios
character(len=30) :: file_name
character(len=1) :: foo


!write (6,*) 'Insert file name'
!read(5,*) file_name
        
open (UNIT=10, FILE='datafile.dat')

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


subroutine var_change (M,num_points,M_coeff)

real*8, intent(in) :: M(num_points, 5)
real*8 :: eta_inter
real*8 :: M_coeff_deriv(3,2), M_coeff2(4,2)
real*8, intent(out) :: M_coeff(4, 2)
integer, intent(in) :: num_points
real*8 :: t(11)
real*8 :: h, num_diff
integer :: i,j,z


! In a mtrix the terms of the polynomial for both the real and the imaginary part are stored


do z = 1, num_points-1
       h = (M(z+1,1)-M(z,1))
!       num_diff = (M(z+1,2)-M(z,2))/(h)
!              write (*,*) 'Numerical derivative', num_diff
!              write (*,*) 'Analytical derivative', M(z,4)

!
!  Change j loop in function fo the number of points desired
!
       do j = 1,4
              eta_inter = M(z,1)*((3.d0-(j-1))/3.d0) + M(z+1,1)*((j-1)/3.d0)
              t(j) = (eta_inter - M(z,1)) / (M(z+1,1)-M(z,1))
              do i = 1,2

                     M_coeff2(1,i) = (2.d0*M(z,i+1)+h*M(z,i+3) - 2*M(z+1,i+1) + h * M(z+1,i+3))
                     M_coeff2(2,i) = (-3.d0*M(z,i+1)+3.d0*M(z+1,i+1)-2.d0*h*M(z,i+3)-h* M(z+1,i+3))
                     M_coeff2(3,i) = h*M(z,i+3)  
                     M_coeff2(4,i) = M(z,i+1)
                      
                     M_coeff(1,i) = M_coeff2(1,i)*t(j)**3.d0
                     M_coeff(2,i) = M_coeff2(2,i)*t(j)**2.d0 
                     M_coeff(3,i) = M_coeff2(3,i)*t(j)
                     M_coeff(4,i) = M_coeff2(4,i)

                     M_coeff_deriv(1,i) = 3.d0*M_coeff2(1,i)*t(j)**2.d0
                     M_coeff_deriv(2,i) = 2.d0*M_coeff2(2,i)*t(j)
                     M_coeff_deriv(3,i) = M_coeff2(3,i) 
              end do                 
              M_coeff_deriv (:,:) = M_coeff_deriv (:,:) / (M(z+1,1)-M(z,1))
              write (*,*) t(j), eta_inter, sum(M_coeff(:,1)), sum(M_coeff(:,2)) , -sum(M_coeff_deriv(:,1))*eta_inter,&
                          -sum(M_coeff_deriv(:,2))*eta_inter
       end do
end do

end subroutine

end module


program Spline

Use spline_module  

implicit none



real*8, allocatable :: M(:,:)
real*8 :: M_coeff (4,2)
integer :: num_points
integer :: i 

call Reading (M,num_points)

!do i = 1, num_points
!       write (*,*) M(i,:)
!end do

call var_change(M,num_points,M_coeff)





end program


