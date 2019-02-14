program theateq
use derivatemp
use constants
implicit none
real*8, allocatable, dimension(:,:,:,:) :: T
real*8, allocatable, dimension(:,:,:) :: dT
real*8 :: tstep = 1E-3, dx
integer, parameter :: timesteps = 1000000
integer :: ilow, ihigh, i, j, k, l, n = 2, resolution = 100
integer :: speed = 9001

n = 10*n + 1
ilow = int(0.3*n+0.7)
ihigh = int(0.7*n+0.3)
dx = boxwidth/float(n-1)
print*,'How much computing POWER?'
read*, speed
if (speed .le. 9000) then
allocate(T(n,n,n,2),dT(n,n,n))
T(:,:,:,1) = 30
T(ilow:ihigh,ilow:ihigh,ilow:ihigh,1) = 100
do i = 1,timesteps-1
	!K1
	call tderiv(T(:,:,:,1), dT)
	T(:,:,:,2) = T(:,:,:,1) + (tstep/6.0)*dT
	!K2
	call tderiv(T(:,:,:,1) + (tstep/2.0)*dT, dT)
	T(:,:,:,2) = T(:,:,:,2) + (tstep/3.0)*dT
	!K3
	call tderiv(T(:,:,:,1) + (tstep/2.0)*dT, dT)
	T(:,:,:,2) = T(:,:,:,2) + (tstep/3.0)*dT
	!K4
	call tderiv(T(:,:,:,1) + (tstep)*dT, dT)
	T(:,:,:,2) = T(:,:,:,2) + (tstep/6.0)*dT
	T(:,:,:,1) = T(:,:,:,2)
end do
open(unit=100, file='ans.txt')
print*,'Writing data...'
do i = 1,n
	do j = 1,n
		do k = 1,n
		write(100,*) i, j, k, 1, T(i,j,k,2)
		end do
	end do
end do
close(100)


else
allocate(T(n,n,n,timesteps), dT(n,n,n))
T(:,:,:,1) = 30
T(ilow:ihigh,ilow:ihigh,ilow:ihigh,1) = 100
do i = 1,timesteps-1
	!K1
	call tderiv(T(:,:,:,i), dT)
	T(:,:,:,i+1) = T(:,:,:,i) + (tstep/6.0)*dT
	!K2
	call tderiv(T(:,:,:,i) + (tstep/2.0)*dT, dT)
	T(:,:,:,i+1) = T(:,:,:,i+1) + (tstep/3.0)*dT
	!K3
	call tderiv(T(:,:,:,i) + (tstep/2.0)*dT, dT)
	T(:,:,:,i+1) = T(:,:,:,i+1) + (tstep/3.0)*dT
	!K4
	call tderiv(T(:,:,:,i) + (tstep)*dT, dT)
	T(:,:,:,i+1) = T(:,:,:,i+1) + (tstep/6.0)*dT
end do
print*,'Writing data...'

open(unit=100, file = 'ans.txt')
do l = 1, int(timesteps/resolution)
	do i = 1,n
		do j = 1,n
			do k = 1,n
			write(100,*) i, j, k, l,  T(i,j,k,resolution*l)
			end do
		end do
	end do
end do
close(100)
end if

deallocate(dT,T)

open(unit=200,file = 'info.txt')
write(200,*) n
write(200,*) timesteps
write(200,*) boxwidth/float(n-1)
write(200,*) resolution
close(200)

end program theateq
