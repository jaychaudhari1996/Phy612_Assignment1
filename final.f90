module variables
implicit none
!constants
integer, parameter :: nrand = 1000	!number of random numbers
integer, parameter :: n = 50		!number of particles
real, parameter :: sbox = 10.0		!size of the box
real, parameter :: dp = 2.0		!diameter of the particle
real, parameter :: vn = 0.5		!velocity normalization factor
real, parameter :: dt = 0.001		!time step
integer, parameter :: t = 1000		!timesteps
real, parameter :: finaltemp = 0.5	!temperature
real,parameter :: eps = 0.01		!effective force distance
real,parameter :: mp = 1.0		!mass of particle
real,parameter :: s1=315.0,s2=578.0,s3=962.0,a1=6984815.0,a2=9865883.0,m1=81916153.0,m2=98498151.0	!random number generators


!variables
real, dimension (nrand) :: rx1,rx2,ry1,ry2,rz1,rz2,rvx1,rvx2,rvy1,rvy2,rvz1,rvz2	!random number dump
real, dimension (n) :: x,y,z,vx,vy,vz				!allocated random numbers
real, dimension (n) :: xp,yp,zp,v,v2,fx,fy,fz,xx,yy,zz		!previous_position,velocity,forces,dump_position calculations
real, dimension (n,n) :: xr,yr,zr,r,ffx,ffy,ffz			!distance, position, forces vectors
integer :: i,j,k,cnt,time					!loopers

real :: ljcut,sfx,sfy,sfz,svx,svy,svz,sv2,avgv,tempn,temp

end module variables


!**************************************************************************************************
!calculating previous position
subroutine prev_pos()
use variables
implicit none
do i = 1,n                                            
	xp(i) = x(i) - vx(i)*dt
	yp(i) = y(i) - vy(i)*dt
	zp(i) = z(i) - vz(i)*dt

	if (xp(i) .lt. 0.0) xp(i) = xp(i) + sbox                     !periodic boundary conditions
	if (xp(i) .gt. sbox) xp(i) = xp(i) - sbox
	if (yp(i) .lt. 0.0) yp(i) = yp(i) + sbox
	if (yp(i) .gt. sbox) yp(i) = yp(i) - sbox
	if (zp(i) .lt. 0.0) zp(i) = zp(i) + sbox
	if (zp(i) .gt. sbox) zp(i) = zp(i) - sbox
end do
return
end subroutine


!**************************************************************************************************
!initializing force
subroutine force()
use variables
implicit none

ljcut = (2**(1/6.0))*dp

do i = 1,n
	fx(i) = 0.0
	fy(i) = 0.0
	fz(i) = 0.0
end do

do i = 1,n-1
	do j = i+1,n

		xr(i,j) = x(i) - x(j)
		yr(i,j) = y(i) - y(j)
		zr(i,j) = z(i) - z(j)

		xr(i,j) = xr(i,j) - sbox*nint(xr(i,j)/sbox)              !minimum image convention
		yr(i,j) = yr(i,j) - sbox*nint(yr(i,j)/sbox)
		zr(i,j) = zr(i,j) - sbox*nint(zr(i,j)/sbox)
		r(i,j) = sqrt(xr(i,j)**2 + yr(i,j)**2 + zr(i,j)**2)

		if (r(i,j) .lt. ljcut) then
        	ffx(i,j) = (48.*eps*xr(i,j)*((dp/r(i,j))**12 - 0.5*(dp/r(i,j))**6))/(r(i,j)**2)
        	ffy(i,j) = (48.*eps*yr(i,j)*((dp/r(i,j))**12 - 0.5*(dp/r(i,j))**6))/(r(i,j)**2)
        	ffz(i,j) = (48.*eps*zr(i,j)*((dp/r(i,j))**12 - 0.5*(dp/r(i,j))**6))/(r(i,j)**2)

		fx(i) = fx(i) + ffx(i,j)
		fy(i) = fy(i) + ffy(i,j)
		fz(i) = fz(i) + ffz(i,j)
		fx(j) = fx(j) - ffx(i,j)
		fy(j) = fy(j) - ffy(i,j)
		fz(j) = fz(j) - ffz(i,j)
		end if        
	end do
end do

sfx = 0.0
sfy = 0.0
sfz = 0.0

do i = 1,n
	sfx = sfx + fx(i)
	sfy = sfy + fy(i)
	sfz = sfz + fz(i)
end do

return
end subroutine  


!**************************************************************************************************
!updating position
subroutine update_pos()
use variables
implicit none

do i = 1,n
        xx(i) = x(i)
        yy(i) = y(i)
        zz(i) = z(i)
	x(i) = 2.*x(i) - xp(i) + ((dt**2)*fx(i))/mp
	y(i) = 2.*y(i) - yp(i) + ((dt**2)*fy(i))/mp
	z(i) = 2.*z(i) - zp(i) + ((dt**2)*fz(i))/mp
	
	if (x(i) .lt. 0.0) x(i) = x(i) + sbox                       !periodic boundary conditions
	if (x(i) .gt. sbox) x(i) = x(i) - sbox
	if (y(i) .lt. 0.0) y(i) = y(i) + sbox
	if (y(i) .gt. sbox) y(i) = y(i) - sbox
	if (z(i) .lt. 0.0) z(i) = z(i) + sbox
	if (z(i) .gt. sbox) z(i) = z(i) - sbox

        xp(i) = xx(i)
        yp(i) = yy(i)
        zp(i) = zz(i)
end do
return
end subroutine 


!**************************************************************************************************
!main program
program MD
use variables

open(20,file="force_time.dat")
open(30,file="momentum_time.dat")
open(40,file="temperature_time.dat")

x(1) = 1.0
y(1) = 1.0
z(1) = 1.0

rx1(1) = s1
ry1(1) = s2
rz1(1) = s3


!**************************************************************************************************
!Assigning positions
do i = 2,nrand
	rx1(i) = mod((rx1(i-1)*a1),m1)
	rx2(i) = sbox*(rx1(i)/m1)
	ry1(i) = mod((ry1(i-1)*a1),m1)
	ry2(i) = sbox*(ry1(i)/m1)
	rz1(i) = mod((rz1(i-1)*a1),m1)
	rz2(i) = sbox*(rz1(i)/m1)
end do

do i = 2,n
	j=1
	k=1
	do while (j.lt.i)
		if (((rx2(k)-x(j))**2 + (ry2(k)-y(j))**2 + (rz2(k)-z(j))**2) .gt. (4.0)) then
		j = j+1
		else
		k = k + 1
		j = 1
		end if
	end do
	x(i) = rx2(k)
	y(i) = ry2(k)
	z(i) = rz2(k)

end do

do i = 2,n
	do j = i+1,n
		xr(i,j) = x(i) - x(j)
		yr(i,j) = y(i) - y(j)
		zr(i,j) = z(i) - z(j)

		xr(i,j) = xr(i,j) - boxl*nint(xr(i,j)/boxl)		!minimum image convention
		yr(i,j) = yr(i,j) - boxl*nint(yr(i,j)/boxl)
		zr(i,j) = zr(i,j) - boxl*nint(zr(i,j)/boxl)

		r(i,j) = sqrt(xr(i,j)**2 + yr(i,j)**2 + zr(i,j)**2)	!distance between every pairs
	
		if (r(i,j) .lt. sig) then                                !check whether there is any overlap
		cnt = cnt+1
		end if
	end do
	
	if ((x(i).gt.sbox) .or. (y(i).gt.sbox) .or. (z(i).gt.sbox)) then    !check whether all particles are inside box
	print*, "Particle outside the box",i
	end if

end do

print*,"Number of overlaps =",cnt


!**************************************************************************************************
!Assigning velocites

rvx1(1) = s1
rvy1(1) = s2
rvz1(1) = s3

avgv = 0.0
sv2 = 0.0
svx = 0.0
svy = 0.0
svz = 0.0

do i = 2,n
	rvx1(i) = mod((rvx1(i-1)*a2),m2)
	vx(i) = (2.*(rvx1(i)/m2) - 1.)*vn
	rvy1(i) = mod((rvy1(i-1)*a2),m2)
	vy(i) = (2.*(rvy1(i)/m2) - 1.)*vn
	rvz1(i) = mod((rvz1(i-1)*a2),m2)
	vz(i) = (2.*(rvz1(i)/m2) - 1.)*vn
end do

do i = 1,n
	svx = svx + vx(i)
	svy = svy + vy(i)
	svz = svz + vz(i)
	avgv = avgv + sqrt(vx(i)**2 + vy(i)**2 + vz(i)**2)
end do

svx = svx/real(n)
svy = svy/real(n)
svz = svz/real(n)
avgv = avgv/real(n) 

do i = 1,n                              !assign centre of mass velocity zero
	vx(i) = vx(i) - svx
	vy(i) = vy(i) - svy
	vz(i) = vz(i) - svz

	v2(i) = vx(i)**2 + vy(i)**2 + vz(i)**2                     
	sv2 = sv2 + v2(i)
end do

sv2 = sv2/real(n)                                         !average velocity
temp = sv2/3.0

print*,"Initial temperature=",temp
tempn = (finaltemp/temp)**0.5                                 !scaling factor

sv2 = 0.0

do i = 1,n                                              !velocity rescaling
	vx(i) = vx(i)*tempn
	vy(i) = vy(i)*tempn
	vz(i) = vz(i)*tempn

	v2(i) = vx(i)**2 + vy(i)**2 + vz(i)**2
	v(i) = sqrt(v2(i))
	sv2 = sv2 + v2(i)
end do

sv2 = sv2/real(n)
temp = sv2/3.0
print*,"Final temperature =",temp
	print*,sum(vx),sum(vy),sum(vz)

	call prev_pos()

do time = 1,t
	call force()
	write (20,*) sfx,sfy,sfz,time
	call update_pos()
end do

end program MD
