program panqake
implicit none

! Lump element model

! TODO - use script for parameters
real :: a, b
real, parameter :: mu0 = 1.26e-6, pi = 3.141593653
real, parameter :: sct = 4e-3
integer, parameter :: nt = 100 , ndp = 3
integer, parameter :: np = ndp*2
!integer, parameter :: dp = selected_real_kind(18,300)
integer :: i, j, k, t, tt,ostep
real*8 :: Ip = 0, Ip0 = 0, Ip1 = 0.1, dIp = 0.1
real*8 :: dt = 1e-6, otime = 1e-3, time, t0 = 0, t1 = 0, th = 0, timefloor, teps = 1e-4
real :: r0 = 7e-3, zt, Aphi, dz

real, dimension(np) :: z0
real*8, dimension(np) :: Iz, Rz, dIz, Itemp, Ic		! Ic needs to be function of field and temp
real, dimension(np,np) :: Lz, iLz
real, dimension(np) :: Rc




print *, "INIT PANQAKE"

open ( unit = 1, file = "currIz.txt")
open ( unit = 2, file = "currIp.txt")
open ( unit = 3, file = "currIr.txt")
open ( unit = 4, file = "currAll.txt")

! init I array to Ip
zt = 2*ndp*sct
do i = 1, np
	Iz(i) = 0
	Itemp(i) = 0
    z0(i) = -zt/2 + sct/2 + (i-1)*sct
	Rz(i) = 1
	Rc(i) = 1
	dIz(i) = 0
	Ic(i) = 450
end do

! Print out general parameters 
print *, "Pancake heights"
print *, z0
print *, "Superconducting matrix resistance"
print *, Rz
print *, "Characteristic coil resistance"
print *, Rc


! Inductance matrix
do i = 1, np
	do j = 1, np
		Aphi = pi*r0**2
		dz = z0(j) - z0(i)
		if (i /= j ) then
			Lz(i,j) = 0
		else
			Lz(i,j) = nt*nt*Aphi*mu0/2*r0**2/((dz**2+r0**2)**(3.0/2))/1e3
		end if
	end do
end do

! invert inductance matrix		why?
call inverse(Lz,iLz,np)
print *, "Inductance matrix L"
do i=1,np
	print *, (Lz(i,j),j=1,np)
end do
print *, "Inverse inductance matrix iL"
do i=1,np
	print *, (iLz(i,j),j=1,np)
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
print*, "Time	Ip	Iz(i...n)"

! main loop
do t = 0, tt
    time = t*dt
	do i = 1, np
		Itemp(i) = Iz(i)
		!Iz(i) = Rc(i)*Ip
		do j = 1, np
				!dIz(i) = dIz(i) + iLz(i,j)*(Rc(j)*Ip - Iz(j)*(Rz(j) + Rc(j)))*dt
			!Iz(i) = Iz(i) + Lz(i,j)*dIz(j)
		end do
		!Iz(i) = Iz(i)/(Rz(i)+Rc(i))
		!Iz(i) = 1/(Rz(i)+Rc(i))*(Rc(i)*Ip-Lz(i,i)*dIz(i))
		dIz(i) = 1/Lz(i,i)*(Rc(j)*Ip - Iz(j)*(Rz(j) + Rc(j)))*dt
		if (Iz(i) > Ic(i)) then
			!Rz(i) = Rz(i)*(1+Ic(i)/Iz(i))
		end if
	end do

	do i = 1, np
		!dIz(i) = (Iz(i) - Itemp(i))/dt
		Iz(i) = Iz(i) + dIz(i)
	end do
	!print *, t, dIz(1:)
		! output at specific time intervals 
		!print *, time, t, (t/ostep), modulo(t,ostep)
        if (modulo(t,ostep) < 1e-6) then
			print 74, t*dt, Ip, Iz(1:)
			74 format( 1x, F6.3, 1x, F10.6, 1x, 99F13.9, 1x)
			write(4,74) t*dt, Ip, Iz(1:)
			do, i = 1, 0, 1
					!print *, Iz
					!print *, z0(i), Lz(i)
					!write(1,*) t*dt, Ip, Iz(i)
					!write(2,*) t*dt, Ip
					!write(3,*) t*dt, Ip - Iz(i)
					!write(4,*) t*dt, Ip, Iz(i)
			end do
        end if

	! Ramp Power supply current
	if (Ip < Ip1) then
		Ip = Ip + dIp
	else 
		Ip = Ip1
	end if
end do

end

! Inverse matrix routine
! credited to Alex G.
subroutine inverse(a,c,n)
implicit none
integer :: n
real, dimension(n,n) :: a, c 
real, dimension(n,n) :: L, U
real, dimension(n) :: b, d, x
real :: coeff
integer :: i, j, k

! initialization for matrices L, U, and b
L=0.0
U=0.0
b=0.0

! forward elimination
do k=1, n-1
	do i=k+1,n
		coeff=a(i,k)/a(k,k)
		L(i,k) = coeff
		do j=k+1,n
			a(i,j) = a(i,j)-coeff*a(k,j)
		end do
	end do
end do

! L and U matrices
! L matrix is a matrix of the elimination coefficient
! diagonal elements are 1.0
do i=1,n
	L(i,i) = 1.0
end do
! U matrix  is the upper triangular part of A
do j=1,n
	do i=1,j
		U(i,j) = a(i,j)
	end do
end do

! compute columns of the inverse matrix C
do k=1,n
	b(k)=1.0
	d(1) = b(1)
! solve Ld=b using forward substitution
		do i=2,n
			d(i)=b(i)
			do j=1,i-1
				d(i) = d(i) - L(i,j)*d(j)
			end do
		end do

		! solve Ux=d using back substitution
		x(n) = d(n)/U(n,n)
		do i = n-1,1,-1
			x(i) = d(i)
			do j=n,i+1,-1
				x(i)=x(i)-U(i,j)*x(j)
			end do
			x(i) = x(i)/u(i,i)
		end do

		! fill solutions x(n) into column k of C
		do i=1,n
			c(i,k) = x(i)
		end do
		b(k) = 0.0
end do
end subroutine inverse


