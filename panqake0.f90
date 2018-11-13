program panqake
implicit none

! Lump element model

! TODO - use script for parameters
real :: a, b
real, parameter :: mu0 = 1.26e-6, pi = 3.141593653
integer, parameter :: nt = 100 , ndp = 1
integer, parameter :: np = 1
integer :: i, j, k, t, tt
real :: Ip = 0, Ip0 = 0, Ip1 = 100, dIp = 10
real :: dt = 1e-3, t0 = 0, t1 = 10, th = 100
real :: ra

real, dimension(np) :: Iz, Lz, Rz,dIz
real, dimension(np) :: Rc


print *, "INIT PANQAKE"

open ( unit = 1, file = "currIz.txt")
open ( unit = 2, file = "currIp.txt")
open ( unit = 3, file = "currIr.txt")
open ( unit = 4, file = "currAll.txt")

! init I array to Ip
do i = 1, np, 1
	Iz(i) = 0
	Lz(i) = mu0*nt*nt/4e-3*100
	Rz(i) = 1e-6
	Rc(i) = 100
	dIz(i) = 0
end do

! differential equation - ramp coil
t1 = Ip1/dIp
tt = (t1 - t0 + th)/dt
Ip = Ip0
dIp = dIp * dt
!Ip1 = tt * dIp

!print *, "Time		Current"
do t = 0, tt, 1
	do i = 1, np, 1
			dIz(i) = (1/Lz(i))*(Rc(i)*Ip - Iz(i)*(Rz(i) + Rc(i)))*dt

			!print *, t*dt, Ip, Iz(i)
			write(1,*) t*dt, Ip, Iz(i)
			write(2,*) t*dt, Ip
			write(3,*) t*dt, Ip - Iz(i)
			write(4,*) t*dt, Ip, Iz(i)

			Iz(i) = Iz(i) + dIz(i)
	end do
	if (Ip < Ip1) then
		Ip = Ip + dIp
	else 
		Ip = Ip1
	end if
end do

end program panqake
