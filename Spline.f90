module spline_module

implicit none 

contains

!Read the points in file and storing data in a matrix (put into subroutine)
!Imaginary part should be stored as complex?
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
end do

end subroutine

!Variable change eta to t??? (Must be generalized for a given 2 points, first try with the first 2 points)

subroutine var_change (M,num_points,M_coeff)

real*8, intent(in) :: M(num_points, 5)
real*8:: eta_inter
real*8, intent(out) :: M_coeff(4, 2)
integer, intent(in) :: num_points
real*8 :: t(11)
integer :: i,j,z


!Just to have a value for the interpolation between the two points
!eta_inter = (M(1,1) + M(2,1))/2.d0
!eta_inter = M(2,1)

do i = 0,10
       eta_inter = M(1,1)*((10.d0-i)/10.d0) + M(2,1)*(i/10.d0)
       t(i+1) = (eta_inter - M(1,1)) / (M(2,1)-M(1,1))
end do

! In a mtrix the terms of the polynomial for both the real and the imaginary part are stored


do z = 1, num_points-1
       do j = 1,11
              eta_inter = M(z,1)*((10.d0-(j-1))/10.d0) + M(z+1,1)*((j-1)/10.d0)
              t(i+1) = (eta_inter - M(z,1)) / (M(z+1,1)-M(z,1))
              do i = 1,2
                     M_coeff(1,i) = (2.d0*t(j)**3.d0 - 3.d0*t(j)**2.d0 + 1.d0) * M(z,i+1)
                     M_coeff(2,i) = (t(j)**3.d0 - 2.d0*t(j)**2.d0 + t(j)) * (M(z+1,1)-M(z,1)) * M(z,i+3)
                     M_coeff(3,i) = (-2.d0*t(j)**3.d0 + 3.d0*t(j)**2.d0) * M(z+1,i+1)
                     M_coeff(4,i) = (t(j)**3.d0 - t(j)**2.d0) * (M(z+1,1)-M(z,1)) * M(z+1,i+3)
!                     write (*,*) 'Real part of the two inital points and the interpolated value'
!                     write (*,*) M(z,2), M(z+1,2), sum(M_coeff(:,1))
!                     write (*,*) 'Imaginary part of the two inital points and the interpolated value'
!                     write (*,*) M(z,3), M(z+1,3), sum(M_coeff(:,2))
              end do                     
              write (*,*) eta_inter, sum(M_coeff(:,1)), sum(M_coeff(:,2)) 
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

!write (*,*) 'Real part of the two inital points and the interpolated value'
!write (*,*) M(1,2), M(2,2), sum(M_coeff(:,1))


!write (*,*) 'Imaginary part of the two inital points and the interpolated value'
!write (*,*) M(1,3), M(2,3), sum(M_coeff(:,2))



end program


