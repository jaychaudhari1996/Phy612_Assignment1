!module variables
implicit none

!constants
integer, parameter :: nrand = 1000	!number of random numbers
integer, parameter :: n = 50		!number of particles
real, parameter :: sbox = 10.0		!size of the box
real, parameter :: dt = 0.001		!time step
real, parameter :: t = 1000.0		!total time
real, parameter :: finaltemp = 0.5	!temperature

real, parameter :: s1,s2,s3,a1,a2,m1,m2	!random number generators
s1=315.
s2=578.
s3=962.
a1=6984815.
a2=9865883.
m1=81916153.
m2=98498151.

!variables
real, dimension (:), allocatable :: rx1,rx2,rvx1,rvx2,rvy1,rvy2,rvz1,rvz2	!random number dump
real, dimension (:), allocatable :: x,y,z,vx,vy,vz				!allocated random numbers
real, dimension (:), allocatable :: xp,yp,zp,v,v2				!calculations
real, dimension (:,:), allocatable :: xr,yr,zr,r				!vectors
integer :: i,j,k								!loopers

allocate (rx1(nrand),rx2(nrand),rvx1(nrand),rvx2(nrand),rvy1(nrand),rvy2(nrand),rvz1(nrand),rvz2(nrand))
allocate (x(n),y(n),z(n),vx(n),vy(n),vz(n))
allocate (xp(n),yp(n),zp(n),v(n),v2(n))
allocate (xr(n,n),yr(n,n),zr(n,n),r(n,n))

end module variables

!**************************************************************************************************
!calculating previous position
subroutine prev_pos()
use variables
do i = 1,n                                            
	xp(i) = x(i) - vx(i)*dt
	yp(i) = y(i) - vy(i)*dt
	zp(i) = z(i) - vz(i)*dt

	if (xp(i) .lt. 0.0) xp(i) = xp(i) + boxl                       !periodic boundary conditions
	if (xp(i) .gt. boxl) xp(i) = xp(i) - boxl
	if (yp(i) .lt. 0.0) yp(i) = yp(i) + boxl
	if (yp(i) .gt. boxl) yp(i) = yp(i) - boxl
	if (zp(i) .lt. 0.0) zp(i) = zp(i) + boxl
	if (zp(i) .gt. boxl) zp(i) = zp(i) - boxl
end do
return
end subroutine
 
program MD
end program MD
