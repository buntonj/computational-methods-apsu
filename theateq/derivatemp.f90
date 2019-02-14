module derivatemp
use constants
implicit none
contains

subroutine tderiv(T, dT)
real*8, dimension(:,:,:) :: T, dT
real*8 :: dx
integer :: n, i, j, k, ilow, ihigh
dT = 0
n = size(T,1)
dx = boxwidth/dble(n-1)
ilow = int(0.3*n+0.7)
ihigh = int(0.7*n+0.3)
!Corners
!1st index = 1 edges absorb heat from sun
dT(1,1,1) = sigma*((Tsun**4)/(csil*denssil*dx))*(Rsun/Rorbit)**2 - sigma*(T(1,1,1)**4)/(csil*denssil*dx) &
			+ alphasil/(dx**2)*(T(1,1,2)+T(1,2,1)+T(2,1,1)-3.0*T(1,1,1))

dT(1,1,n) = sigma*((Tsun**4)/(csil*denssil*dx))*(rsun/rorbit)**2 - sigma*(T(1,1,n)**4)/(csil*denssil*dx) &
			+ alphasil/(dx**2)*(T(1,1,n-1)+T(1,2,n)+T(2,1,n)-3.0*T(1,1,n))

dT(1,n,n) = sigma*((Tsun**4)/(csil*denssil*dx))*(rsun/rorbit)**2 - sigma*(T(1,n,n)**4)/(csil*denssil*dx) & 
			+ alphasil/(dx**2)*(T(1,n,n-1)+T(1,n-1,n)+T(2,1,n)-3.0*T(1,n,n))

dT(1,n,1) = sigma*((Tsun**4)/(csil*denssil*dx))*(rsun/rorbit)**2 - sigma*(T(1,n,1)**4)/(csil*denssil*dx) &
			+ alphasil/(dx**2)*(T(1,n,2)+T(1,n-1,1)+T(2,n,1)-3.0*T(1,n,1))


dT(n,n,n) = -sigma*(T(n,n,n)**4)/(csil*denssil*dx)*(dx**2) &
		+ alphasil/(dx**2)*(T(n,n,n-1)+T(n,n-1,n)+T(n-1,n,n)-3.0*T(n,n,n))
dT(n,n,1) = -sigma*(T(n,n,1)**4)/(csil*denssil*dx)*(dx**2) &
		+ alphasil/(dx**2)*(T(n,n,2)+T(n,n-1,1)+T(n-1,n,1)-3.0*T(n,n,1))
dT(n,1,1) = -sigma*(T(n,1,1)**4)/(csil*denssil*dx)*(dx**2) &
		+ alphasil/(dx**2)*(T(n,1,2)+T(n,2,1)+T(n-1,1,1)-3.0*T(n,1,1))
dT(n,1,n) = -sigma*(T(n,1,n)**4)/(csil*denssil*dx)*(dx**2) &
		+ alphasil/(dx**2)*(T(n,1,n-1)+T(n,2,n)+T(n-1,1,n)-3.0*T(n,1,n))


!Edges
do i = 2, n-1
	!1st index = 1 edges absorb heat from sun
	dT(1,1,i) =  sigma*((Tsun**4)/(csil*denssil*dx))*(rsun/rorbit)**2 - sigma*(T(1,1,i)**4)/(csil*denssil*dx) &
			+ alphasil/(dx**2)*(T(1,1,i+1)+T(1,1,i-1)+T(1,2,i)+T(2,1,i)-4.0*T(1,1,i))

	dT(1,n,i) = sigma*((Tsun**4)/(csil*denssil*dx))*(rsun/rorbit)**2 - sigma*(T(1,n,i)**4)/(csil*denssil*dx) &
			+ alphasil/(dx**2)*(T(1,n,i+1)+T(1,n,i-1)+T(1,n-1,i)+T(2,n,i)-4.0*T(1,n,i))

	dT(1,i,n) = sigma*((Tsun**4)/(csil*denssil*dx))*(rsun/rorbit)**2 - sigma*(T(1,i,n)**4)/(csil*denssil*dx) &
			+ alphasil/(dx**2)*(T(1,i+1,n)+T(1,i-1,n)+T(1,i,n-1)+T(2,i,n)-4.0*T(1,i,n))

	dT(1,i,1) = sigma*((Tsun**4)/(csil*denssil*dx))*(rsun/rorbit)**2 - sigma*(T(1,i,1)**4)/(csil*denssil*dx) &
			+ alphasil/(dx**2)*(T(1,i+1,1)+T(1,i-1,1)+T(2,i,1)+T(1,i,2)-4.0*T(1,i,1))


	dT(n,1,i) = -sigma*(T(n,1,i)**4)/(csil*denssil*dx) &
			+ alphasil/(dx**2)*(T(n,1,i-1)+T(n,1,i+1)+T(n,2,i)+T(n-1,1,i)-4.0*T(n,1,i))
	dT(n,n,i) = -sigma*(T(n,n,i)**4)/(csil*denssil*dx) &
			+ alphasil/(dx**2)*(T(n,n,i-1)+T(n,n,i+1)+T(n,n-1,i)+T(n-1,n,i)-4.0*T(n,n,i))
	dT(n,i,1) = -sigma*(T(n,i,1)**4)/(csil*denssil*dx) &
			+ alphasil/(dx**2)*(T(n,i,2)+T(n,i+1,1)+T(n,i-1,1)+T(n-1,i,1)-4.0*T(n,i,1))
	dT(n,i,n) = -sigma*(T(n,i,n)**4)/(csil*denssil*dx) &
			+ alphasil/(dx**2)*(T(n,i,n-1)+T(n,i+1,n)+T(n,i-1,n)+T(n-1,i,n)-4.0*T(n,i,n))
	dT(i,1,1) = -sigma*(T(i,1,1)**4)/(csil*denssil*dx) &
			+ alphasil/(dx**2)*(T(i,1,2)+T(i,2,1)+T(i+1,1,1)+T(i-1,1,1)-4.0*T(i,1,1))
	dT(i,n,1) = -sigma*(T(i,n,1)**4)/(csil*denssil*dx) &
			+ alphasil/(dx**2)*(T(i,n,2)+T(i,n-1,1)+T(i+1,n,1)+T(i-1,n,1)-4.0*T(i,n,1))
	dT(i,1,n) = -sigma*(T(i,1,n)**4)/(csil*denssil*dx) &
			+ alphasil/(dx**2)*(T(i,1,n-1)+T(i,2,n)+T(i+1,1,n)+T(i-1,1,n)-4.0*T(i,1,n))
	dT(i,n,n) = -sigma*(T(i,n,n)**4)/(csil*denssil*dx) &
			+ alphasil/(dx**2)*(T(i,n,n-1)+T(i,n-1,n)+T(i+1,n,n)+T(i-1,n,n)-4.0*T(i,n,n))

end do


! Faces
do i = 2, n-1
	do j = 2, n-1
	dT(1,i,j) = sigma*((Tsun**4)/(csil*denssil*dx))*(rsun/rorbit)**2 - sigma*(T(1,i,j)**4)/(denssil*dx) &
			+ alphasil/(dx**2)*(T(1,i,j+1)+T(1,i,j-1)+T(1,i+1,j)+T(1,i-1,j)+T(2,i,j)-5.0*T(1,i,j))
			
	dT(i,1,j) = -sigma*(T(i,1,j)**4)/(csil*denssil*dx) &
			+ alphasil/(dx**2)*(T(i,1,j+1)+T(i,1,j-1)+T(i+1,1,j)+T(i-1,1,j)+T(i,2,j)-5.0*T(i,1,j))
	dT(i,j,1) = -sigma*(T(i,j,1)**4)/(csil*denssil*dx) &
			+ alphasil/(dx**2)*(T(i,j+1,1)+T(i,j-1,1)+T(i+1,j,1)+T(i-1,j,1)+T(i,j,2)-5.0*T(i,j,1))
	dT(n,i,j) = -sigma*(T(n,i,j)**4)/(csil*denssil*dx) &
			+ alphasil/(dx**2)*(T(n,i,j+1)+T(n,i,j-1)+T(n,i+1,j)+T(n,i-1,j)+T(n-1,i,j)-5.0*T(n,i,j))
	dT(i,n,j) = -sigma*(T(i,n,j)**4)/(csil*denssil*dx) &
			+ alphasil/(dx**2)*(T(i,n,j+1)+T(i,n,j-1)+T(i+1,n,j)+T(i-1,n,j)+T(i,n-1,j)-5.0*T(i,n,j))
	dT(i,j,n) = -sigma*(T(i,j,n)**4)/(csil*denssil*dx) &
			+ alphasil/(dx**2)*(T(i,j-1,n)+T(i,j+1,n)+T(i+1,j,n)+T(i-1,j,n)+T(i,j,n-1)-5.0*T(i,j,n))
	end do
end do

do i = 2,n-1
	do j = 2,n-1
		do k = 2,n-1
		dT(i,j,k) = alphasty/(dx**2)*(T(i,j,k+1)+T(i,j,k-1)+T(i,j+1,k)+T(i,j-1,k)+T(i-1,j,k)+T(i+1,j,k)-6.0*T(i,j,k))
		end do
	end do
end do

dT(ilow:ihigh,ilow:ihigh,ilow:ihigh) = 0

end subroutine
end module
