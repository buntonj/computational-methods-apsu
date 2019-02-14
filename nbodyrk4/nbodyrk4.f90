program nbodyrk4
use constants
use morse
use lennardjones
implicit none

real*8, parameter :: pi = acos(-1.0)
real*8, parameter :: kb = 1.38065D-23
real*8, parameter :: T = 100, dt = 1D-14
integer, parameter :: edgex = 5, edgey = 5, edgez = 5, tsteps = 10000
integer, parameter :: nbodies = edgex*edgey*edgez
real*8, dimension(1:3,1:nbodies,0:tsteps+1) :: r = 0, v = 0, F = 0
real*8, dimension(3,0:tsteps) :: p = 0
real*8, dimension(0:tsteps) :: E = 0
real*8, dimension(3,1:nbodies,4) :: KR = 0, KV = 0
real*8 :: average
integer :: i, j, k, it, counter, start, finish, persecond, time
real, dimension(3,2) :: rand
character :: pottype

! -O3 for third level optimization
call system_clock(count_rate=persecond)
print*,'Lennard-Jones or Morse?'
print*,'(L/M)'
read(*,*) pottype

! PERFORM LENNARD-JONES POTENTIAL SIMULATION
if(pottype == 'L') then
print*,'Performing Lennard-Jones potential simulation...'
!do time = 1,10
! INITIALIZE r(:,:,0)
counter = 1
do i = 1, edgex
	do j = 1, edgey
		do k = 1, edgez
		r(1,counter,0) = (i-1)*rm
		r(2,counter,0) = (j-1)*rm
		r(3,counter,0) = (k-1)*rm
		counter = counter + 1
		end do
	end do
end do
print*,'Position initialized...'

! INITIALIZE v(:,:,0)
do i = 1, nbodies
	call random_number(rand)
	v(:,i,0) = sqrt(-2.0*Kb*T/m*log(rand(:,1)))*cos(2*pi*rand(:,2))*0.95
end do
print*,'Velocity initialized...'

!call system_clock(COUNT=start)
! CALCULATE INITIAL ENERGY and MOMENTUM
call ljenergy(r(:,:,0), v(:,:,0), E(0))
call ljmomentum(v(:,:,0),p(:,0))
print*,'Force, energy, momentum initialized...'

! TIME LOOP
do it = 1,tsteps
	call ljforce(r(:,:,it-1),F(:,:,it-1))
	!Calculate K1s
	KV(:,:,1) = dt*F(:,:,it-1)/m
	KR(:,:,1) = dt*V(:,:,it-1)
	!Calculate K2s
	call ljforce((r(:,:,it-1)+KR(:,:,1)/2),KV(:,:,2))
	KV(:,:,2) = dt*KV(:,:,2)/m
	KR(:,:,2) = dt*(V(:,:,it-1)+KV(:,:,1)/2)
	!Calculate K3s
	call ljforce((r(:,:,it-1)+KR(:,:,2)/2),KV(:,:,3))
	KV(:,:,3) = dt*KV(:,:,3)/m
	KR(:,:,3) = dt*(V(:,:,it-1)+KV(:,:,2)/2)
	!Calculate K4s
	call ljforce((r(:,:,it-1)+KR(:,:,3)),KV(:,:,4))
	KV(:,:,4) = dt*KV(:,:,4)/m
	KR(:,:,4) = dt*(V(:,:,it-1)+KV(:,:,3))
	!Average all four Ks
	v(:,:,it) = v(:,:,it-1) + (1.0/6.0)*(KV(:,:,1)+2*KV(:,:,2)+2*KV(:,:,3)+KV(:,:,4))
	r(:,:,it) = r(:,:,it-1) + (1.0/6.0)*(KR(:,:,1)+2*KR(:,:,2)+2*KR(:,:,3)+KR(:,:,4))
	KV(:,:,:) = 0
	KR(:,:,:) = 0
	call ljenergy(r(:,:,it), v(:,:,it), E(it))
	call ljmomentum(v(:,:,it),p(:,it))
end do ! End timesteps loop
!call system_clock(COUNT=finish)
!average = average + (float(finish)-float(start))/float(persecond)
!end do
!print*,'Average computation time (seconds):', average/10


!PERFORM MORSE POTENTIAL SIMULATION
else if (pottype == 'M') then
print*,'Performing Morse potential simulation...'

!do time = 1,10
! INITIALIZE r(:,:,0)
counter = 1
do i = 1, edgex
	do j = 1, edgey
		do k = 1, edgez
		r(1,counter,0) = (i-1)*re
		r(2,counter,0) = (j-1)*re
		r(3,counter,0) = (k-1)*re
		counter = counter + 1
		end do
	end do
end do
print*,'Position initialized...'

! INITIALIZE v(:,:,0)
do i = 1, nbodies
	call random_number(rand)
	v(:,i,0) = sqrt(-2.0*Kb*T/m*log(rand(:,1)))*cos(2*pi*rand(:,2))*0.95
end do
print*,'Velocity initialized...'

!call system_clock(count=start)
! CALCULATE INITIAL ENERGY and MOMENTUM
call morseenergy(r(:,:,0), v(:,:,0), E(0))
call morsemomentum(v(:,:,0),p(:,0))
print*,'Energy, momentum initialized...'

! TIME LOOP
do it = 1,tsteps
	call morseforce(r(:,:,it-1),F(:,:,it-1))
	!Calculate K1s
	KV(:,:,1) = dt*F(:,:,it-1)/m
	KR(:,:,1) = dt*V(:,:,it-1)
	!Calculate K2s
	call morseforce((r(:,:,it-1)+KR(:,:,1)/2),KV(:,:,2))
	KV(:,:,2) = dt*KV(:,:,2)/m
	KR(:,:,2) = dt*(V(:,:,it-1)+KV(:,:,1)/2)
	!Calculate K3s
	call morseforce((r(:,:,it-1)+KR(:,:,2)/2),KV(:,:,3))
	KV(:,:,3) = dt*KV(:,:,3)/m
	KR(:,:,3) = dt*(V(:,:,it-1)+KV(:,:,2)/2)
	!Calculate K4s
	call morseforce((r(:,:,it-1)+KR(:,:,3)),KV(:,:,4))
	KV(:,:,4) = dt*KV(:,:,4)/m
	KR(:,:,4) = dt*(V(:,:,it-1)+KV(:,:,3))
	!Average all four Ks
	v(:,:,it) = v(:,:,it-1) + (1.0/6.0)*(KV(:,:,1)+2*KV(:,:,2)+2*KV(:,:,3)+KV(:,:,4))
	r(:,:,it) = r(:,:,it-1) + (1.0/6.0)*(KR(:,:,1)+2*KR(:,:,2)+2*KR(:,:,3)+KR(:,:,4))
	KV(:,:,:) = 0
	KR(:,:,:) = 0
	call morseenergy(r(:,:,it), v(:,:,it), E(it))
	call morsemomentum(v(:,:,it),p(:,it))
end do ! End timesteps loop
!call system_clock(COUNT=finish)
!average = average + float(finish-start)/float(persecond)
!end do
!print*,'Average computation time (seconds):', average/10
end if


! OUTPUT DATA
print*,'Writing data...'
call outputdata(r, v, F, E, P, nbodies, tsteps, dt, pottype)

end program nbodyrk4


!DATA OUTPUT SUBROUTINE
subroutine outputdata(r, v, F, E, P, nbodies, tsteps, dt, pottype)
integer :: nbodies, tsteps, it
real*8, dimension(3,1:nbodies,0:tsteps+1) :: r, v, F
real*8, dimension(3,0:tsteps) :: P
real*8, dimension(0:tsteps) :: E
real*8 :: dt
character :: pottype

open(unit = 100, file = 'pos.txt')
open(unit = 200, file = 'vel.txt')
open(unit = 300, file = 'force.txt')
open(unit = 400, file = 'E.txt')
open(unit = 500, file = 'P.txt')
open(unit = 600, file = 'lasttype.txt')

do it = 0, tsteps
do i = 1, nbodies
write(100,*) i, r(:,i,it), norm2(r(:,i,it)), it*dt !particle #, x, y, z, time
write(200,*) i, v(:,i,it), norm2(v(:,i,it)), it*dt !VEL
write(300,*) i, F(:,i,it), norm2(F(:,i,it)), it*dt !F
end do
write(400,*) it*dt, E(it) ! total E
write(500,*) it*dt, p(:,it), norm2(p(:,it)) !total P in each direction, and mag
end do
write(600,*) pottype
write(600,*) nbodies
write(600,*) tsteps

close(100)
close(200)
close(300)
close(400)
close(500)
close(600)

end subroutine

!RANDOM SEED SUBROUTINE
subroutine init_random_seed()
     integer :: i, n, clock
     integer, dimension(:), allocatable :: seed

     CALL RANDOM_SEED(size = n)
     allocate(seed(n))

     CALL SYSTEM_CLOCK(COUNT=clock)

     seed = clock + 37 * (/ (i - 1, i = 1, n) /)
     CALL RANDOM_SEED(PUT = seed)

     deallocate(seed)
end subroutine

!RANDOM INDEX FUNCTION
integer function randomindex(length)
     real :: random
	call random_number(random)
     randomindex = floor(random*length + 1.0)
return
end
