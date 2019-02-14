module integrateme
use constants
implicit none
contains

function integral(n, m)
real*8 :: integral
integer, intent(in) :: n, m
real*8 :: x, dx
real*8, dimension(5) :: fx
integer :: nsteps, i, j

if ((n .gt. 0) .and. (m .le. 0)) then
	integral = 0.0d0
	return
else if ((m .gt. 0)  .and. (n .le. 0)) then
	integral = 0.0d0
	return
end if

nsteps = 4*2*nmax
dx = l/(4*dble(nsteps))
integral = 0.0
x = 0.0

if (n .gt. 0) then
	do i = 1,nsteps
		x = dble(i-1)*4.0*dx
		do j = 1,5
		fx(j) = depth*dsin(dble(n)*pi*x/l)*dsin(dble(m)*pi*x/l)*dexp(-sd*(x/l)**2)
		x = x + dx
		end do
	integral = integral + (4.0/45.0)*dx*(7.0*fx(1) + 32.0*fx(2) + 12.0*fx(3) + 32.0*fx(4) + 7.0*fx(5))
	end do
else 
	do i = 1,nsteps
		x = dble(i-1)*4.0*dx
		do j = 1,5
		fx(j) = depth*dcos(dble(n)*pi*x/l)*dcos(dble(m)*pi*x/l)*dexp(-sd*(x/l)**2)
		x = x + dx
		end do
	integral = integral + (4.0/45.0)*dx*(7.0*fx(1) + 32.0*fx(2) + 12.0*fx(3) + 32.0*fx(4) + 7.0*fx(5))
	end do
end if
return
end function

end module integrateme
