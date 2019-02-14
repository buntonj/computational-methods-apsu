program quantumwell
use constants
use integrateme
implicit none
real*8, allocatable, dimension(:,:,:,:) :: wfns
real*8, allocatable, dimension(:,:) :: H 
real*8, allocatable, dimension(:,:) :: oneDeewfns
real*8, allocatable, dimension(:) :: x
real*8 :: dx
integer, parameter :: resolution = 100
integer :: k1x, k1y, k1z, k2x, k2y, k2z
integer :: m1, m2, m3, n, a, b, i
character:: oned

! LAPACK DSYEV ARGUMENTS
integer :: numel, info
real*8, allocatable, dimension(:) :: eigenvalues, work
integer :: lwork
character*1 :: jobz = 'V', uplo = 'U'

print*, '1D?'
read(*,*) oned

if(oned == 'n') then

print*,'Starting sin/sin/sin terms...'

!----------------3D ATTEMPT--------------!
! x: sin
! y: sin
! z: sin

allocate(H(nmax**3,nmax**3), wfns(elevels,0:resolution,0:resolution,0:resolution))

do k1x = 1,nmax
	do k1y = 1, nmax
		do k1z = 1, nmax
			do k2x = k1x, nmax
				do k2y = k1y, nmax
					do k2z = k1z, nmax
					a = (k1x-1)*nmax**2 + (k1y-1)*nmax + k1z
					b = (k2x-1)*nmax**2 + (k2y-1)*nmax + k2z

					H(a,b) = integral(k1x,k2x)*integral(k1y,k2y)*integral(k1z,k2z)*(2/L)**3
					H(b,a) = H(a,b)
					end do
				end do
			end do
		H(a,a) = ((pi*hbar/l)**2)/m *(k1x**2 + k1y**2 + k1z**2) + H(a,a)
		end do
	end do
end do

numel = size(H,1)
lwork = 3*numel
allocate(eigenvalues(numel), work(lwork),x(0:resolution))
call dsyev(jobz,uplo,numel,H,numel,eigenvalues,work,lwork,info)
wfns = 0.0d0

dx = 2*L/(size(x)-1)
do i = 0,resolution
x(i) = -L + float(i)*dx
end do
print*, dx
do n = 1,elevels
	do k1x = 1,nmax
		do k1y = 1,nmax
			do k1z = 1,nmax
			a = (k1x-1)*nmax**2 + (k1y-1)*nmax + k1z
				do m1 = 0, resolution
					do m2 = 0, resolution
						 do m3 = 0, resolution
						 wfns(n,m1,m2,m3) = wfns(n,m1,m2,m3) + (2/l)**(3/2)*H(n,a)*sin(k1x*pi/l*x(m1))&
						 	*sin(k1y*pi*x(m2)/l)*sin(k1z*pi*x(m3)/l)
						 end do
					end do
				end do
			end do
		end do
	end do
end do

open(unit=100,file='sinsinsinvectors.txt')
print*,"Writing sin/sin/sin data..."
do m1 = 1,elevels
	do k1x = 0,resolution
		do k1y = 0,resolution
			do k1z = 0, resolution
			!Columns: energy eigenvalue, i, j, k, psi(x,y,z)
			write(100,*) eigenvalues(m1), k1x, k1y, k1z, wfns(m1,k1x,k1y,k1z)
			end do
		end do
	end do
end do
deallocate(H,eigenvalues,work,wfns)



!----------------3D ATTEMPT--------------!
! x: sin
! y: sin
! z: cos
print*,'Starting sin/sin/cos terms...'
allocate(H((nmax+1)*nmax**2,(nmax+1)*nmax**2), wfns(elevels,0:resolution,0:resolution,0:resolution))

do k1x = 1,nmax
	do k1y = 1, nmax
		do k1z = 0, -nmax, -1
			do k2x = k1x, nmax
				do k2y = k1y, nmax
					do k2z = k1z, -nmax,-1
					a = (k1x-1)*nmax*(nmax+1) + (k1y-1)*(nmax+1) - k1z + 1
					b = (k2x-1)*nmax*(nmax+1) + (k2y-1)*(nmax+1) - k2z + 1

					H(a,b) = integral(k1x,k2x)*integral(k1y,k2y)*integral(k1z,k2z)*(2/L)**3
					H(b,a) = H(a,b)
					end do
				end do
			end do
		H(a,a) = ((pi*hbar/l)**2)/m *(k1x**2 + k1y**2 + k1z**2) + H(a,a)
		end do
	end do
end do

numel = size(H,1)
lwork = 3*numel
allocate(eigenvalues(numel), work(lwork))
call dsyev(jobz,uplo,numel,H,numel,eigenvalues,work,lwork,info)
wfns = 0.0d0


do n = 1,elevels
	do k1x = 1,nmax
		do k1y = 1,nmax
			do k1z = 0,-nmax,-1
			a = (k1x-1)*nmax*(nmax+1) + (k1y-1)*(nmax+1) - k1z + 1
				do m1 = 0, resolution
					do m2 = 0, resolution
						 do m3 = 0, resolution
						 wfns(n,m1,m2,m3) = wfns(n,m1,m2,m3) + (2/l)**(3/2)*H(n,a)*sin(k1x*pi/l*x(m1))&
						 	*sin(k1y*pi*x(m2)/l)*cos(k1z*pi*x(m3)/l)
						 end do
					end do
				end do
			end do
		end do
	end do
end do

open(unit=101,file='sinsincosvectors.txt')
print*,"Writing sin/sin/cos data..."
do m1 = 1,elevels
	do k1x = 0,resolution
		do k1y = 0,resolution
			do k1z = 0, resolution
			!Columns: energy eigenvalue, i, j, k, psi(x,y,z)
			write(101,*) eigenvalues(m1), k1x, k1y, k1z, wfns(m1,k1x,k1y,k1z)
			end do
		end do
	end do
end do
close(101)
deallocate(H,eigenvalues,work,wfns)


!----------------3D ATTEMPT--------------!
! x: sin
! y: cos
! z: cos
print*,'Starting sin/cos/cos terms...'
allocate(H(((nmax+1)**2)*nmax,((nmax+1)**2)*nmax), wfns(elevels,0:resolution,0:resolution,0:resolution))

do k1x = 1,nmax
	do k1y = 0, -nmax, -1
		do k1z = 0, -nmax, -1
			do k2x = k1x, nmax
				do k2y = k1y, -nmax,-1
					do k2z = k1z, -nmax,-1
					a = (k1x-1)*(nmax+1)**2 + (-k1y)*(nmax) - k1z + 1
					b = (k2x-1)*(nmax+1)**2 + (-k2y)*(nmax) - k2z + 1

					H(a,b) = integral(k1x,k2x)*integral(k1y,k2y)*integral(k1z,k2z)*(2/L)**3
					H(b,a) = H(a,b)
					end do
				end do
			end do
		H(a,a) = ((pi*hbar/l)**2)/m *(k1x**2 + k1y**2 + k1z**2) + H(a,a)
		end do
	end do
end do

numel = size(H,1)
lwork = 3*numel
allocate(eigenvalues(numel), work(lwork))
call dsyev(jobz,uplo,numel,H,numel,eigenvalues,work,lwork,info)
wfns = 0.0d0


do n = 1,elevels
	do k1x = 1,nmax
		do k1y = 0,-nmax,-1
			do k1z = 0,-nmax,-1
			a = (k1x-1)*(nmax+1)**2 + (-k1y)*(nmax) - k1z + 1
				do m1 = 0, resolution
					do m2 = 0, resolution
						 do m3 = 0, resolution
						 wfns(n,m1,m2,m3) = wfns(n,m1,m2,m3) + (2/l)**(3/2)*H(n,a)*sin(k1x*pi/l*x(m1))&
						 	*cos(k1y*pi*x(m2)/l)*cos(k1z*pi*x(m3)/l)
						 end do
					end do
				end do
			end do
		end do
	end do
end do

open(unit=102,file='sincoscosvectors.txt')
print*,"Writing sin/cos/cos data..."
do m1 = 1,elevels
	do k1x = 0,resolution
		do k1y = 0,resolution
			do k1z = 0, resolution
			!Columns: energy eigenvalue, i, j, k, psi(x,y,z)
			write(102,*) eigenvalues(m1), k1x, k1y, k1z, wfns(m1,k1x,k1y,k1z)
			end do
		end do
	end do
end do
close(102)
deallocate(H,eigenvalues,work,wfns)



!----------------3D ATTEMPT--------------!
! x: cos
! y: cos
! z: cos
print*,'Starting cos/cos/cos terms...'
allocate(H((nmax+1)**3,(nmax+1)**3), wfns(elevels,0:resolution,0:resolution,0:resolution))

do k1x = 0,-nmax,-1
	do k1y = 0, -nmax, -1
		do k1z = 0, -nmax, -1
			do k2x = k1x, -nmax,-1
				do k2y = k1y, -nmax,-1
					do k2z = k1z, -nmax,-1
					a = (-k1x)*(nmax+1)**2 + (-k1y)*(nmax+1) - k1z + 1
					b = (-k2x)*(nmax+1)**2 + (-k2y)*(nmax+1) - k2z + 1

					H(a,b) = integral(k1x,k2x)*integral(k1y,k2y)*integral(k1z,k2z)*(2/L)**3
					H(b,a) = H(a,b)
					end do
				end do
			end do
		H(a,a) = ((pi*hbar/l)**2)/m *(k1x**2 + k1y**2 + k1z**2) + H(a,a)
		end do
	end do
end do

numel = size(H,1)
lwork = 3*numel
allocate(eigenvalues(numel), work(lwork))
call dsyev(jobz,uplo,numel,H,numel,eigenvalues,work,lwork,info)
wfns = 0.0d0


do n = 1,elevels
	do k1x = 0,-nmax,-1
		do k1y = 0,-nmax,-1
			do k1z = 0,-nmax,-1
			a = (-k1x)*(nmax+1)**2 + (-k1y)*(nmax+1) - k1z + 1
				do m1 = 0, resolution
					do m2 = 0, resolution
						 do m3 = 0, resolution
						 wfns(n,m1,m2,m3) = wfns(n,m1,m2,m3) + (2/l)**(3/2)*H(n,a)*cos(k1x*pi/l*x(m1))&
						 	*cos(k1y*pi*x(m2)/l)*cos(k1z*pi*x(m3)/l)
						 end do
					end do
				end do
			end do
		end do
	end do
end do

open(unit=104,file='coscoscosvectors.txt')
print*,"Writing cos/cos/cos data..."
do m1 = 1,elevels
	do k1x = 0,resolution
		do k1y = 0,resolution
			do k1z = 0, resolution
			!Columns: energy eigenvalue, i, j, k, psi(x,y,z)
			write(104,*) eigenvalues(m1), k1x, k1y, k1z, wfns(m1,k1x,k1y,k1z)
			end do
		end do
	end do
end do
close(104)
deallocate(H,eigenvalues,work,wfns)

deallocate(x)


else

!------------------1 - DIMENSIONAL CASE------------------!
!1-D TRY

!Sin Terms
allocate(H(nmax,nmax), oneDeewfns(elevels,resolution))
H = 0.0D0
oneDeewfns = 0.0D0

do k1z = 1, nmax
	do k2z = k1z, nmax
	a = k1z
	b = k2z
	H(a,b) = integral(k1z,k2z)*(2.0/L)
	H(b,a) = H(a,b)
	end do
	H(a,a) = ((pi*hbar*k1z/l)**2)/(2*m) + H(a,a)
end do

numel = size(H,1)
lwork = 3*numel

open(unit=800, file='before.txt')
do a = 1,nmax
	write(800,*) (H(a,b), b = 1,nmax)
end do
close(800)

allocate(eigenvalues(numel), work(lwork),x(resolution))
call dsyev(jobz,uplo,numel,H,numel,eigenvalues,work,lwork,info)
onedeewfns = 0.0d0
if(info .ne. 0) then
print*,'LAPACK error:', info
endif

dx = 2*L/resolution
do i = 1,resolution
x(i) = -L + float(i)*dx
end do

open(unit=200, file='eigendata.txt')
do a = 1,nmax
	write(200,*) (H(a,b), b = 1,nmax)
end do
close(200)

do n = 1,elevels
	do k1z = 1,nmax
		 do m3 = 1,resolution
		 oneDeewfns(n,m3) = oneDeewfns(n,m3) + sqrt(2.0/l)*H(n,k1z)*sin(dble(k1z)*pi*x(m3)/l)
		 end do
	end do
	oneDeewfns(n,:) = oneDeewfns(n,:)/(norm2(oneDeewfns(n,:)))
end do

open(unit=100,file='sinvectors.txt')
print*,"Writing sin data..."
do m1 = 1,elevels
	do k1z = 1, resolution
	!Columns: energy eigenvalue, i, j, k, psi(x,y,z)
	write(100,*) m1, eigenvalues(m1),  k1z, x(k1z), oneDeewfns(m1,k1z)
	end do
end do
deallocate(H,eigenvalues,work,onedeewfns,x)
close(100)



!Cos terms
allocate(H(0:nmax,0:nmax), oneDeewfns(elevels,resolution))

do k1z = 0, -nmax, -1
	do k2z = k1z, -nmax,-1
	a = abs(k1z)
	b = abs(k2z)
	H(a,b) = integral(k1z,k2z)*(2.0/L)
	H(b,a) = H(a,b)
	end do
end do

do k1z = 0,-nmax, -1
	a =  abs(k1z)
	H(a,a) = ((pi*hbar*k1z/l)**2.0)/m + H(a,a)
end do

numel = size(H,1)
lwork = 3*numel
open(unit=200, file='eigendata.txt')
do a = 0,nmax
	write(200,*) (H(a,b), b = 0,nmax)
end do
close(200)

allocate(eigenvalues(numel), work(lwork),x(resolution))
call dsyev(jobz,uplo,numel,H,numel,eigenvalues,work,lwork,info)

if(info .ne. 0) then
print*,'LAPACK error:', info
endif

onedeewfns = 0.0d0
dx = 2*L/resolution

do i = 1,resolution
x(i) = -L + float(i)*dx
end do


do n = 1,elevels
	do k1z = 0,-nmax,-1
		 do m3 = 1,resolution
		 oneDeewfns(n,m3) = oneDeewfns(n,m3) + (2.0/l)**(1/2)*H(n,abs(k1z))*cos(dble(k1z)*pi*x(m3)/l)
		 end do
	end do
	oneDeewfns(n,:) = oneDeewfns(n,:)/(norm2(oneDeewfns(n,:)))
end do

open(unit=400,file='cosvectors.txt')
print*,"Writing cos data..."
do m1 = 1,elevels
	do k1z = 1, resolution
	!Columns: energy eigenvalue, i, j, k, psi(x,y,z)
	write(400,*) m1, eigenvalues(m1),  k1z, x(k1z), oneDeewfns(m1,k1z)
	end do
end do

open(unit=500, file='potential.txt')
do i = 1,resolution
	write(500,*) x(i), depth*dexp(-sd*(x(i)/l)**2)
end do
close(500)

deallocate(H,eigenvalues,work,onedeewfns,x)


end if
close(400)
open(unit=300, file='info.txt')
write(300,*) nmax
write(300,*) depth
write(300,*) sd
write(300,*) l
write(300,*) resolution
write(300,*) dx
close(300)
end program quantumwell
