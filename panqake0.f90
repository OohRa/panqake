program panqake
implicit none

! Lump element model

! TODO - use script for parameters
real :: a, b
real, parameter :: mu0 = 1.26e-6, pi = 3.141593653
real, parameter :: sct = 4e-3
integer, parameter :: nt = 100 , ndp = 1
integer, parameter :: np = 1
!integer, parameter :: dp = selected_real_kind(18,300)
integer :: i, j, k, t, tt,ostep
real :: Ip = 0, Ip0 = 0, Ip1 = 10, dIp = 10
real*8 :: dt = 1e-6, otime = 1e-3, time, t0 = 0, t1 = 10, th = 0, timefloor, teps = 1e-4
real :: r0 = 7e-3, zt

real, dimension(np) :: z0
real, dimension(np) :: Iz, Lz, Rz, dIz
real, dimension(np) :: Rc


print *, "INIT PANQAKE"

open ( unit = 1, file = "currIz.txt")
open ( unit = 2, file = "currIp.txt")
open ( unit = 3, file = "currIr.txt")
open ( unit = 4, file = "currAll.txt")

! init I array to Ip
zt = ndp * sct
do i = 1, np, 1
	Iz(i) = 0
    z0(i) = -zt/2+i*sct
	Lz(i) = mu0/2 * r0**2*nt / (r0**2)**(3.0/2)
	Rz(i) = 1e-6
	Rc(i) = 1
	dIz(i) = 0
end do

! differential equation - ramp coil
t1 = Ip1/dIp
tt = (t1 - t0 + th)/dt
Ip = Ip0
dIp = dIp * dt
teps = dt*10
ostep = otime/dt
!Ip1 = tt * dIp

! Data output for post processing
!print *, "Time		Current"
print *, "Total time steps: ", tt
print*, "Time step: ", dt
do t = 0, tt, 1
    time = t*dt
	do i = 1, np, 1
        dIz(i) = (1/Lz(i))*(Rc(i)*Ip - Iz(i)*(Rz(i) + Rc(i)))*dt
        ! output at specific time intervals 
        !print *, time, t, (t/ostep)
        if (modulo(t,ostep) == 0) then
            print 75, t*dt, Ip, Iz(i), Ip-Iz(i)
            75 format( 1x, F6.4, 1x, 3F6.2, 1x)
            !print *, z0(i), Lz(i)
            write(1,*) t*dt, Ip, Iz(i)
            write(2,*) t*dt, Ip
            write(3,*) t*dt, Ip - Iz(i)
            write(4,*) t*dt, Ip, Iz(i)
        end if

        Iz(i) = Iz(i) + dIz(i)
	end do
	if (Ip < Ip1) then
		Ip = Ip + dIp
	else 
		Ip = Ip1
	end if
end do

end program panqake
