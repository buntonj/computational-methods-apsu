program rocketship
use constants
use gravipot
implicit none

real*8, allocatable, dimension(:,:) :: r, v, f, kr, kv
real*8, allocatable, dimension(:) :: E, p, m, thrustdirection
integer, parameter :: tstage0index = tstage0/dt, tstage1index = tstage1/dt
integer :: i, j, stage = 0, alert = 0, crashalert = 0
integer :: tsteps = (tstage0index + tstage1index) + 200000

allocate(r(2,0:tsteps),v(2,0:tsteps),f(2,0:tsteps),p(0:tsteps),E(0:tsteps),m(0:tsteps),kr(2,4),kv(2,4),thrustdirection(2))

m(0) = totalmass

r(1,0) = 0.0D0
r(2,0) = Re

v(1,0) = 0! Re*pi/(12*60*60)
v(2,0) = 0
thrustdirection(1) = 0.05
thrustdirection(2) = 1
p(0) = m(0)*norm2(v(:,0))

! DO NOT TOUCH ANYMORE YOU FUCKING JACKASS
do i = 1,tstage0index
	!K1 CALCULATION
	!call gforce(r(:,i-1), v(:,i-1), F(:,i-1), m(i-1), stage)
	call gforce(r(:,i-1), thrustdirection, F(:,i-1), m(i-1), stage)
	kv(:,1) = (dt/m(i-1))*f(:,i-1)
	kr(:,1) = dt*(v(:,i-1))
	m(i) =  m(i-1) - stage0fuelmass/(2.0*tstage0)
	
	!K2 CALCULATION
	!call gforce(r(:,i-1) + kr(:,1)/2.0, v(:,i-1) + kv(:,1)/2.0, kv(:,2), m(i), stage)
	call gforce(r(:,i-1) + kr(:,1)/2.0, thrustdirection + kv(:,1)/2.0, kv(:,2), m(i), stage)
	kv(:,2) = (dt/m(i))*kv(:,2)
	kr(:,2) = dt*(v(:,i-1)+kv(:,2)/2.0)

	!K3 CALCULATION
	!call gforce(r(:,i-1) + kr(:,2)/2.0, v(:,i-1) + kv(:,2)/2.0, kv(:,3), m(i), stage)
	call gforce(r(:,i-1) + kr(:,2)/2.0, thrustdirection + kv(:,2)/2.0, kv(:,3), m(i), stage)
	kv(:,3) = (dt/m(i))*kv(:,3)
	kr(:,3) = dt*(v(:,i-1) + kv(:,3)/2.0)

	m(i) =  m(i) - stage0fuelmass/(2.0*tstage0)
	!K4 CALCULATION
	!call gforce(r(:,i-1) + kr(:,3), v(:,i-1) + kv(:,3), kv(:,4), m(i), stage)
	call gforce(r(:,i-1) + kr(:,3), thrustdirection + kv(:,3), kv(:,4), m(i), stage)
	kv(:,4) = (dt/m(i))*kv(:,4)
	kr(:,4) = dt*(v(:,i-1) + kv(:,4))

	!AVERAGE
	v(:,i) = v(:,i-1) + (1.0/6.0)*(kv(:,1) + 2.0*kv(:,2) + 2.0*kv(:,3) + kv(:,4))
	r(:,i) = r(:,i-1) + (1.0/6.0)*(kr(:,1) + 2.0*kr(:,2) + 2.0*kr(:,3) + kr(:,4))
	call energy(r(:,i), v(:,i), E(i), m(i))
	thrustdirection = v(:,i)
	kv = 0.0d0
	kr = 0.0d0
	p(i) = m(i) * norm2(v(:,i))
	if (norm2(r(:,i)) .lt. Re) then
		print*, 'At timestep :', i
		print*,'Rocket crashed into earth!'
		goto 100
	else if (norm2(v(:,i)) .gt. 0.99*sqrt(GG*Me/norm2(r(:,i))) .and. norm2(v(:,i)) .le. 1.01*sqrt(GG*Me/norm2(r(:,i)))) then
		print*, 'At timestep :', i
		print*, 'CIRCULAR ORBIT CONDITIONS STAGE 0, SWITCHING to STAGE 1!'
		stage = 2
		j = i
		goto 200
	end if
end do

200 continue
m(i + 1) = stage1mass
stage = 1

do i = j,tsteps
	!K1 CALCULATION
	call gforce(r(:,i-1), v(:,i-1), F(:,i-1), m(i-1), stage)
	kv(:,1) = (dt/m(i-1))*f(:,i-1)
	kr(:,1) = dt*(v(:,i-1)+kv(:,1))
	m(i) =  m(i-1) - stage0fuelmass/(2.0*tstage0)*(2-stage)

	!K2 CALCULATION
	call gforce(r(:,i-1) + kr(:,1)/2.0, v(:,i-1) + kv(:,1)/2.0, kv(:,2), m(i-1), stage)
	kv(:,2) = (dt/m(i))*kv(:,2)
	kr(:,2) = dt*(v(:,i-1)+kv(:,2)/2.0)

	!K3 CALCULATION
	call gforce(r(:,i-1) + kr(:,2)/2.0, v(:,i-1) + kv(:,2)/2.0, kv(:,3), m(i-1), stage)
	kv(:,3) = (dt/m(i))*kv(:,3)
	kr(:,3) = dt*(v(:,i-1) + kv(:,3)/2.0)

	m(i) =  m(i) - stage0fuelmass/(2.0*tstage0)*(2-stage)
	!K4 CALCULATION
	call gforce(r(:,i-1) + kr(:,3), v(:,i-1) + kv(:,3), kv(:,4), m(i-1), stage)
	kv(:,4) = (dt/m(i))*kv(:,4)
	kr(:,4) = dt*(v(:,i-1) + kv(:,4))

	!AVERAGE
	v(:,i) = v(:,i-1) + (1.0/6.0)*(kv(:,1) + 2.0*kv(:,2) + 2.0*kv(:,3) + kv(:,4))
	r(:,i) = r(:,i-1) + (1.0/6.0)*(kr(:,1) + 2.0*kr(:,2) + 2.0*kr(:,3) + kr(:,4))
	call energy(r(:,i), v(:,i), E(i), m(i))
	
	p(i) = m(i) * norm2(v(:,i))

	kv = 0.0d0
	kr = 0.0d0
	if (norm2(r(:,i)) .lt. Re) then
		print*, 'Rocket crashed into earth!'
	else if (norm2(v(:,i)) .ge. sqrt(GG*Me/norm2(r(:,i)))*0.99&
		.and. norm2(v(:,i)) .le. sqrt(GG*Me/norm2(r(:,i)))*1.01) then
		stage = 2
		print*, 'CIRCULAR ORBIT IN STAGE 1!'
		goto 300
	end if
end do

300 continue
do i = j,tsteps
	!K1 CALCULATION
	call gforce(r(:,i-1), v(:,i-1), F(:,i-1), m(j), stage)
	kv(:,1) = (dt/m(j))*f(:,i-1)
	kr(:,1) = dt*(v(:,i-1)+kv(:,1))

	!K2 CALCULATION
	call gforce(r(:,i-1) + kr(:,1)/2.0, v(:,i-1) + kv(:,1)/2.0, kv(:,2), m(j), stage)
	kv(:,2) = (dt/m(j))*kv(:,2)
	kr(:,2) = dt*(v(:,i-1)+kv(:,2)/2.0)

	!K3 CALCULATION
	call gforce(r(:,i-1) + kr(:,2)/2.0, v(:,i-1) + kv(:,2)/2.0, kv(:,3), m(j), stage)
	kv(:,3) = (dt/m(j))*kv(:,3)
	kr(:,3) = dt*(v(:,i-1) + kv(:,3)/2.0)

	m(i) =  m(i) - stage0fuelmass/(2.0*tstage0)*(2-stage)
	!K4 CALCULATION
	call gforce(r(:,i-1) + kr(:,3), v(:,i-1) + kv(:,3), kv(:,4), m(j), stage)
	kv(:,4) = (dt/m(j))*kv(:,4)
	kr(:,4) = dt*(v(:,i-1) + kv(:,4))

	!AVERAGE
	v(:,i) = v(:,i-1) + (1.0/6.0)*(kv(:,1) + 2.0*kv(:,2) + 2.0*kv(:,3) + kv(:,4))
	r(:,i) = r(:,i-1) + (1.0/6.0)*(kr(:,1) + 2.0*kr(:,2) + 2.0*kr(:,3) + kr(:,4))
	call energy(r(:,i), v(:,i), E(i), m(i))
	p(i) = m(i) * norm2(v(:,i))
	
	kv = 0.0d0
	kr = 0.0d0
	if (norm2(r(:,i)) .lt. Re .and. crashalert == 0) then
		print*, 'Rocket crashed into earth!'
		crashalert = 1
	else if (norm2(v(:,i)) .ge. sqrt(GG*Me/norm2(r(:,i)))*0.9&
		.and. norm2(v(:,i)) .le. sqrt(GG*Me/norm2(r(:,i)))*1.1 .and. alert == 0) then
		stage = 2
		print*, 'CIRCULAR ORBIT IN FREE FLOATING!'
		alert = 1
	end if
end do


100 continue
open(unit=100, file='pos.txt')
open(unit=200, file='mass.txt')
open(unit=300, file='vel.txt')
open(unit=400, file='E.txt')
open(unit=500, file='P.txt')
do i = 0,tsteps
	write(100,*) r(:,i)
	write(200,*) m(i)
	write(300,*) v(:,i)
	write(400,*) E(i)
	write(500,*) p(i)
end do
close(100)
close(200)
close(300)
close(400)
close(500)

deallocate(r,v,F,E,p,m)
end program rocketship
