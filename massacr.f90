! ----------------------------------------------------------------------------------%%
!
! MASSACR
!
! SUMMARY: main method runs fluid dynamic simulation coupled to geochemical
!          simulation and writes selected model output to file
!
! TO RUN: make -f theMakeFile
!		  mpirun -np 4 ./massacr
!
! ----------------------------------------------------------------------------------%%

PROGRAM main

use globals
use initialize
!use alteration
!use netcdf

implicit none

include 'mpif.h'
! INCLUDE "IPhreeqc.f90.inc"
save

! functions within massacr.f90
interface

	! solves thermal energy equation
	function h_next (h, psi, rho_in, phi_in, u_in, v_in, frac6_in, temp6_in, dt_in)
		use globals
		use initialize
		implicit none
		! integers
		integer :: i, j, n, ii, m=3
		! inputs
		real(4) :: sx, sy, qx, qy, rho_in(xn,yn), phi_in(xn,yn)
		! velocity stuff
		real(4) :: uf(xn,yn), vf(xn,yn), u_in(xn,yn), v_in(xn,yn)
		real(4) :: u(xn,yn), v(xn,yn), uLong((xn-2)*(yn-2)), vLong((xn-2)*(yn-2))
		real(4) ::  velocities0(xn,2*yn)
		! matrix stuff
		real(4) :: h(xn,yn), h_next(xn,yn), psi(xn,yn)
		! real(4) :: aa((xn-2)*(yn-2),(xn-2)*(yn-2)), a((xn-2)*(yn-2),(xn-2)*(yn-2)+1)
		real(4) :: aBand((xn-2)*(yn-2),5), bBand((xn-2)*(yn-2),5)
		! real(4) :: bb((xn-2)*(yn-2),(xn-2)*(yn-2)), b((xn-2)*(yn-2),(xn-2)*(yn-2)+1)
		real(4) :: h0(xn,yn), uVec((xn-2)*(yn-2)), h_nextRow((xn-2)*(yn-2))
		real(4) :: kMatLong((xn-2)*(yn-2))
		real(4) :: mn(xn,yn)
		real(4) :: sxMat(xn,yn), syMat(xn,yn), sxLong((xn-2)*(yn-2)), syLong((xn-2)*(yn-2))
		real(4) :: qxMat(xn,yn), qyMat(xn,yn), qxLong((xn-2)*(yn-2)), qyLong((xn-2)*(yn-2))
		real(4) :: frac6_in(yn,2), temp6_in(yn,2), dt_in
	end function h_next

	! solves streamfunction vorticity equation
	function psi_next (h, rhs0, psi, rho_in, phi_in, perm_in, band_in, permx, permy, stage,frac6_in)
		use globals
		use initialize
		implicit none
		! integers
		integer :: i, j, ii, n, m, stage
		! inputs
		real(4) :: rhs0(xn,yn), rhs1(xn,yn), rhsLong(longP)
		real(4) :: h(xn,yn), psi(xn,yn), rho_in(xn,yn), phi_in(xn,yn), perm_in(xn,yn)
		! matrix stuff
		real(4) :: uVec(longP), psiLong((xn)*(yn)), psi_nextRow(longP)
		real(4) :: psi_next(xn,yn)
		real(4) :: mn(xn,yn)
		! back to band
		real(4) :: aBand0(longP,4*((yn/2)-2) + 3), band_in(longP,4*((yn/2)-2) + 3)
		real(4) :: rhoLong(longP)
		real(4) :: permx(xn,yn), permy(xn,yn), frac6_in(yn,2)
	end function psi_next



	function make_band(perm_in,phi_in,permx,permy,rho_in)
		use globals
		use initialize
		implicit none
		integer :: i, j, ii, n, m
		real(4) :: perm_in(xn,yn), phi_in(xn,yn), rho_in(xn,yn)
		real(4) :: permx(xn,yn), permy(xn,yn), permLong(longP)
		real(4) :: permxLong(longP), permyLong(longP)
		real(4) :: innerBand(longP,2*((yn/2)-2) + 1), make_band(longP,2*((yn/2)-2) + 1)
		real(4) :: permx_left(xn,yn), permx_right(xn,yn), permy_bottom(xn,yn), permy_top(xn,yn)
		real(4) :: permx_left_long(longP), permx_right_long(longP), permy_bottom_long(longP), permy_top_long(longP)
		real(4) :: perm_long(longP)
	end function make_band




	! calculates fluid density
	function rho_next (h_in)
		use globals
		use initialize
		implicit none
		integer :: i,j
		real(4) :: h_in(xn,yn), rho_next(xn,yn)
	end function rho_next

	! calculates fluid density
	function rho_one(h_in)
		use globals
		use initialize
		implicit none
		integer :: i,j
		real(4) :: h_in, rho_one
	end function rho_one

	! calculates viscosity
	function visc_next(h_in)
		use globals
		use initialize
		implicit none
		integer :: i,j
		real(4) :: h_in(xn,yn), visc_next(xn,yn)
	end function visc_next

	! calculates velocities from streamfunction values
	function velocities(psi)
		use globals
		use initialize
		implicit none
		real(4) :: velocities(xn,2*yn), psi(xn,yn)
		real(4) :: u0(xn,yn), v0(xn,yn)
	end function velocities

	! calculates partial derivative of any 1D or 2D array
	function partial(array,rows,cols,d1,d2,dim)
		use globals
		use initialize
		implicit none
		integer :: rows, cols, dim, i, j, ii, jj
		real(4) :: array(rows,cols), d1, d2, d
		real(4) :: partial(rows,cols)
	end function partial

	! calculates partial derivative of any 1D or 2D array
	function partial_edge(array,rows,cols,d1,d2,dim)
		use globals
		use initialize
		implicit none
		integer :: rows, cols, dim, i, j, ii, jj
		real(4) :: array(rows,cols), d1, d2, d
		real(4) :: partial_edge(rows,cols)
	end function partial_edge

	! calculates partial derivative of any 1D or 2D array
	function partial_edge_p(array,rows,cols,d1,d2,dim)
		use globals
		use initialize
		implicit none
		integer :: rows, cols, dim, i, j, ii, jj
		real(4) :: array(rows,cols), d1, d2, d
		real(4) :: partial_edge_p(rows,cols)
	end function partial_edge_p

	! writes 2D array to file
	function write_matrix ( m, n, table, filename )
		use globals
		implicit none
		integer :: m, n, j, output_status, unit0, reclen
		character ( len = * ) filename
		character ( len = 30 ) string
		real(4)  :: table(m,n) , write_matrix
	end function write_matrix

	! writes 1D array to file
	function write_vec ( n, vector, filename )
		use globals
		implicit none
		integer :: n, j, output_status, unit0
		character ( len = * ) filename
		real(4)  :: vector(n), write_vec
	end function write_vec

end interface

!-DECLARE EVERYTHING

! dependent variable arrays
real(4) :: h(xn,yn), psi(xn,yn), pside(xn,yn) ! xn rows deep & yn columns wide
real(4) :: hmat((xn*tn/(mstep*ar)),yn), psimat((xn*tn/(mstep*ar)),yn), rhomat((xn*tn/(mstep*ar)),yn)
real(4) :: velocities0(xn,2*yn)
real(4) :: phimat((xn*tn/(mstep*ar)),yn), umat((xn*tn/(mstep*ar)),yn), vmat((xn*tn/(mstep*ar)),yn), rhsmat((xn*tn/(mstep*ar)),yn)
real(4) :: permxMat((xn*tn/(mstep*ar)),yn), permyMat((xn*tn/(mstep*ar)),yn)
real(4) :: psiCoarseMat((xn*tn/(mstep*ar)),yn), uCoarseMat((xn*tn/(mstep*ar)),yn), vCoarseMat((xn*tn/(mstep*ar)),yn)
real(4) :: permmat((xn*tn/(mstep*ar)),yn)
real(4) :: u(xn,yn), v(xn,yn),  permx(xn,yn), permy(xn,yn), uTest(xn,yn), vTest(xn,yn), psiTest(xn,yn), hTest(xn,yn)


! autumn performance review
integer :: counti, countf, count_rate
real :: timeBit

! material properties
real(4) :: rho(xn,yn), visc(xn,yn)
real(4) :: rhs0(xn,yn)
integer :: unit
real(4) :: phiCoarse(xn/cell,yn/cell)
real(4) :: phi0(xn,yn), phi(xn,yn)

! netCDF & output stuff
integer :: xInt, yInt, tInt, hInt, uInt, vInt
integer :: ncid
integer :: x_dimid, y_dimid, t_dimid, h_dimid, u_dimid, v_dimid
integer :: x_varid, y_varid, t_varid, h_varid, u_varid, v_varid
integer :: i, j, ii, m, n, jj
real(4) :: yep



! benchmark stuff
real(4) :: nusseltLocalv(xn,1), nuBar




! solute transport stuff
real(4) :: uTransport(xn/cell,yn/cell), vTransport(xn/cell,yn/cell)
real(4) :: uCoarse(xn/cell,yn/cell), vCoarse(xn/cell,yn/cell)
real(4) :: psiCoarse(xn/cell,yn/cell), velocitiesCoarse0(xn/cell,2*yn/cell)
real(4) :: permeability0(xn,yn)

! message passing stuff
integer, parameter :: max_rows = 10000000
integer, parameter :: send_data_tag = 2001, return_data_tag = 2002
integer :: my_id, root_process, ierr, status(MPI_STATUS_SIZE)
integer :: num_procs, an_id, num_rows_to_receive
integer :: avg_rows_per_process, num_rows, num_rows_to_send
integer :: end_row, sender, start_row, num_rows_received
real(4) :: vector(max_rows), vector2(max_rows), partial_sum, sum
real(4) :: local_mean, global_mean
real(4) :: hLocal(xn*(yn/2)), dt_local
integer :: order



INTEGER(KIND=4) :: id, all=190


! MIGRATING ALTERATION STUFF





!path = "/home/navah/input/"
path2 = ""














write(*,*) "testing..."

!-INITIALIZE ALL PROCESSORS

! process #0 is the root process
root_process = 0

! initialize a process
call MPI_INIT ( ierr )

! find out the process ID and how many processes were started so far
call MPI_COMM_RANK (MPI_COMM_WORLD, my_id, ierr)
call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)

! print out current processor id
write(*,*) "my_id:", my_id
write(*,*) " "

! what to do if you are the master processor
if (my_id .eq. root_process) then

!-DO STUFF WITH THE MASTER PROCESSOR

! initialize domain geometry
call init()




permeability0 = permeability
phiCoarse = .1
phi = .1
phi = 0.1
phi0 = phi



! boundary & initial condtions for flow equations
psi=0.0
frac6 = 0.0
! frac6(:,1) = 1.0e-3
! frac6(:,2) = -1.0e-3
temp6 = 280.0
psi(1,1:yn) = bcyPsi(1,1:yn)
psi(xn,1:yn) = bcyPsi(2,1:yn)
psi(1:xn,1) = bcxPsi(1:xn,1)
psi(1:xn,yn) = bcxPsi(1:xn,2)
psi = 0.0!dx*scope*1000.0 - dx*scope*0.0
rho = rho_fluid
visc = viscosity

h = param_tsw

!no slant IC, 11/15/15
do i=1,xn
	do ii=1,yn

		 ! (.68-(0.00133*i*1.0))

		h(i,1) = param_tsw!480.0 + ( first-(0.0016*i*factor)) * dy/1.8

		h(i,ii) = param_tsw!480.0*(h(i,1)/h(1,1)) + (param_tsw-(480.0*(h(i,1)/h(1,1)) ))*((-y_min)+y(ii))/((-y_min))
		!h(i,ii) = h(i,ii) - (400.0 + (param_tsw-400.0)*((-y_min)-max(param_o_rhs,param_o))/((max(param_o_rhs,param_o))))
	end do
end do

! do i=1,xn
! 	do ii=1,yn
! 		if (y(ii) .ge. -max(param_o_rhs,param_o)) then
! 			h(i,ii) = 275.0
! 		end if
! 	end do
! end do


! put initial values into array for file
!hmat(1:xn,1:yn) = h
!psimat(1:xn,1:yn) = psi
!umat(1:xn,1:yn) = u
!vmat(1:xn,1:yn) = v

uTransport = 0.0
vTransport = 0.0
uCoarse = 0.0
vCoarse = 0.0







permx = partial((phi/(permeability)),xn,yn,dx,dy,1)
permy = partial((phi/(permeability)),xn,yn,dx,dy,2)

outerBand = make_band(permeability,phi,permx,permy,rho)
outerBand = band(outerBand,2*((yn/2)-2) + 1,longP)

! permx = 0.0
! permy = 0.0






!-DYNAMICS LOOP
	! this is the main loop that does all the solving for tn timesteps
	do j = crashstep, tn
	!do j = 2, 50


! 			! JDF
! 			! move altered cells
! 			if ( mod(j,nint(tn*cell/(xn*2.3))) .eq. 0) then
! 				write(*,*) "shift"
! 				do i =2,(xn/cell)
! 	 				primaryShift(i,:,:) = primary(i-1,:,:)
! 	 				secondaryShift(i,:,:) = secondary(i-1,:,:)
! 				end do
! 				primaryShift(1,:,5) = 9.67700
! 				secondaryShift(1,:,:) = 0.0
! 				primary = primaryShift
! 				secondary = secondaryShift
! 			end if

	! 		! solve thermal energy equation

if (j .eq. crashstep) then
	write(*,*) "STARTING STEP:" , j
	write(*,*) " "
	write(*,*) " "
	write(*,*) " "
	OPEN(UNIT=8, status = 'replace', FILE=trim(path_final) // 'dynamicStep.txt')
	write(*,*) "opened"
	write(8,*) 0
	close ( 8 )
	!velocitiesCoarse0 = 0.0
end if


if (mod(j,mstep) .eq. 0) then
	write(*,*) "STARTING STEP:" , j
	write(*,*) " "
	write(*,*) " "
	write(*,*) " "
	OPEN(UNIT=8, status = 'replace', FILE=trim(path_final) // 'dynamicStep.txt')
	write(*,*) "opened"
	write(8,*) j/mstep
	close ( 8 )
	!velocitiesCoarse0 = 0.0
end if




if (mod(j,mstep/10) .eq. 0) then
	write(*,*) "WRITING TO DYNAMIC SUB STEP" , j
	write(*,*) " "
	write(*,*) " "
	write(*,*) " "
		! OPEN(UNIT=88, status = 'replace', FILE=trim(path_final) // 'dynamicSubStep.txt')
		! write(88,*) mod(j,mstep)
		! close ( 88 )
end if

! write(*,*) " "
! write(*,*) " "
! write(*,*) " "
! write(*,*) " "
! write(*,*) "j step:" , j

if (restart .ne. 1) then



		dt_bit = dt


		h = h_next(h, psi,rho,phi,u,v,frac6,temp6,dt_bit)
		h = h_bc(h)


		! short outcrop outflow condition
		do jj=2,yn-1
		do i=1,xn

			if ((v(i,jj) .ge. 0.0) .and. (mask(i,jj) .eq. 25.0)  ) then
				h(i,jj+1) = (4.0/3.0)*h(i,jj) - (1.0/3.0)*h(i,jj-1)
			end if

			if ((v(i,jj) .ge. 0.0) .and. (mask(i,jj) .eq. 12.5) ) then
				h(i,jj+1) = (4.0/3.0)*h(i,jj) - (1.0/3.0)*h(i,jj-1)
			end if
			if ((v(i,jj) .ge. 0.0) .and. (mask(i,jj) .eq. 17.5) ) then
				h(i,jj+1) = (4.0/3.0)*h(i,jj) - (1.0/3.0)*h(i,jj-1)
			end if

		end do
		end do

		! inner boundaries outflow condtion
		do jj=2,yn-1
		do i=1,xn

				if ((mask(i,jj) .eq. 5.0) .and. (mask(i,jj) .eq. 5.0) .and. (u(i,jj) .gt. 0.0)) then
					h(i+1,jj) = (4.0/3.0)*h(i,jj) - (1.0/3.0)*h(i-1,jj)
				end if
				if ((mask(i,jj) .eq. 12.5) .and. (u(i,jj) .gt. 0.0)) then
					h(i+1,jj) = (4.0/3.0)*h(i,jj) - (1.0/3.0)*h(i-1,jj)
				end if

				! since these are also a 50.0 bc, they maybe aren't wilcock's but isothermal??
				if ((mask(i,jj) .eq. 10.0) .and. (mask(i,jj) .eq. 10.0).and. (u(i,jj) .lt. 0.0)) then
					h(i-1,jj) = (4.0/3.0)*h(i,jj) - (1.0/3.0)*h(i+1,jj)
				end if
				if ((mask(i,jj) .eq. 17.5) .and. (u(i,jj) .lt. 0.0)) then
					h(i-1,jj) = (4.0/3.0)*h(i,jj) - (1.0/3.0)*h(i+1,jj)
				end if

		end do
		end do

		! insulating boundaries during spinup phase
		if (iter .eq. 0) then
			do ii=2,yn-1
					if ((mask(f_index1,ii) .eq. 6.0) .or. (mask(f_index1,ii) .eq. 6.5) .or. (mask(f_index1,ii) .eq. 6.1) .or. (mask(f_index1,ii) .eq. 6.05)) then
						temp6(ii,2) = (4.0/3.0)*h(f_index1,ii) - (1.0/3.0)*h(f_index1+1,ii)
					end if
					if ((mask(f_index1-1,ii) .eq. 3.0) .or. (mask(f_index1-1,ii) .eq. 3.5) .or. (mask(f_index1-1,ii) .eq. 3.1) .or. (mask(f_index1-1,ii) .eq. 3.05)) then
						temp6(ii,1) = (4.0/3.0)*h(f_index1-1,ii) - (1.0/3.0)*h(f_index1-1-1,ii)
					end if

					if ((mask(f_index1,ii) .eq. 6.5)) then
						temp6(ii+1,2) = (4.0/3.0)*temp6(ii,2) - (1.0/3.0)*temp6(ii-1,2)
					end if

					if ((mask(f_index1-1,ii) .eq. 3.5)) then
						temp6(ii+1,1) = (4.0/3.0)*temp6(ii,1) - (1.0/3.0)*temp6(ii-1,1)
					end if
			end do
		end if

		! fracture temperature set by vertical heat advection after spinup phase
		if ((iter .eq. 1)) then

! 		temp6_mid = temp6
!
				do ii=2,yn-1
					if ((mask(f_index1-1,ii) .eq. 3.05)) then
						temp6(ii,1) = (4.0/3.0)*h(f_index1-1,ii) - (1.0/3.0)*h(f_index1-1-1,ii)
					end if
				end do

				do ii=2,yn-1
					if ((mask(f_index1-1,ii) .eq. 3.1) .or. (mask(f_index1-1,ii) .eq. 3.0) .or. (mask(f_index1-1,ii) .eq. 3.5)) then
						temp6(ii,1) = temp6(ii-1,1) - (dy*lambdaMat(f_index1-1,ii)/(dx*4179.0*frac6(ii,1))) * (h(f_index1-1,ii) - h(f_index1-2,ii))
					end if
				end do

! 				temp6 = temp6_mid
!
! 				do ii=2,yn-1
! 					if ((mask(f_index1-1,ii) .eq. 3.1)) then
! 						temp6_a(ii) = 0.0
! 						temp6_b(ii) = 1.0 + frac6(ii,1)
! 					end if
! 				end do
!

! 				if (mod(j,50) .eq. 0) then
!
! 				do jj=1,10000000
!
! 					temp6 = temp6_mid
! 					do ii=2,yn-1
!
! 						if ((mask(f_index1-1,ii) .eq. 3.1) .or. (mask(f_index1-1,ii) .eq. 3.0) .or. (mask(f_index1-1,ii) .eq. 3.5)) then
! 							temp6_mid(ii,1) = temp6(ii,1) - frac6(ii,1)*(dt/(rho_fluid*param_f_dx*10000000.0))*(temp6(ii,1) - temp6(ii-1,1))/dy
! 						end if
! ! 						if (mask(f_index1-1,ii) .eq. 3.0) then
! ! 							temp6_mid(ii,1) = temp6(ii,1) - frac6(ii,1)*(dt/(rho_fluid*param_f_dx*10000000.0))*(temp6(ii,1) - temp6(ii-1,1))/dy
! ! 						end if
! ! 						if (mask(f_index1-1,ii) .eq. 3.5) then
! ! 							temp6_mid(ii,1) = temp6(ii,1) - frac6(ii,1)*(dt/(rho_fluid*param_f_dx*10000000.0))*(temp6(ii,1) - temp6(ii-1,1))/dy
! ! 						end if
! 					end do
!
! 				end do
! ! 				write(*,*) "mod 50"
! ! 				write(*,*) temp6(:,1)
! ! 				write(*,*) " "
!
! ! ! TESTING SOMETHING
! ! do ii=2,yn-1
! ! 	if ((mask(f_index1-1,ii) .eq. 3.1) .or. (mask(f_index1-1,ii) .eq. 3.0) .or. (mask(f_index1-1,ii) .eq. 3.5)) then
! ! 		temp6_mid(ii,1) = 328.0
! ! 	end if
! ! end do
!
! 				end if



				! top of fracture is no heat flow out? wilcox condition?
				do ii=yn/2,yn-1
						if ((mask(f_index1,ii) .eq. 6.5)) then
							temp6(ii+1,2) = (4.0/3.0)*temp6(ii,2) - (1.0/3.0)*temp6(ii-1,2)
						end if

						if ((mask(f_index1-1,ii) .eq. 3.5)) then
							temp6(ii+1,1) = (4.0/3.0)*temp6(ii,1) - (1.0/3.0)*temp6(ii-1,1)
						end if
				end do

		end if


		do ii=1,yn-1
			temp6(ii,2) = temp6(ii,1)
		end do



		rho = rho_next(h)
		visc = visc_next(h)




		! solve streamfunction-vorticity equation

		rhs0 = -1.0*partial((rho/rho_fluid)-1.0,xn,yn,dx,dy,1)

		frac6_last = frac6

		! find properties at top and base of fracture
		do jj=yn/2,yn-1
			if (maskP(f_index1-1,jj) .eq. 3.1) then
				h_base = temp6(jj,1) ! h(f_index1-1,jj)!
				y_base = y(jj)
				jj_base = jj
			end if
			if (maskP(f_index1-1,jj) .eq. 3.5) then
				h_top = param_tsw ! temp6(jj,1)!
				y_top = y(jj)
				jj_top = jj
			end if
		end do

		h_adjacent = sum(temp6(jj_base:jj_top,1))/(jj_top-jj_base+1)

! 		write(*,*) "tube temps"
! 		write(*,*) temp6(jj_base:jj_top,1)
!
! 		write(*,*) "h_adjacent"
! 		write(*,*) h_adjacent
!
! 		write(*,*) "jj_top"
! 		write(*,*) jj_top
!
! 		write(*,*) "jj_base"
! 		write(*,*) jj_base
!
! 		write(*,*) "rho_adjacent"
! 		write(*,*) rho_one(h_adjacent)
!

		if ((j .gt. spinup-1)) then
			iter = 1
		end if

! 		if ((j .eq. spinup-1)) then
! 			h_adjacent = 273.0
! 		end if


		if (iter .eq. 1) then
			!write(*,*) "loop of fracj 1"
			do jj=yn/2,yn-1
				if ((maskP(f_index1,jj) .eq. 6.0) .or. (maskP(f_index1,jj) .eq. 6.5) .or. (maskP(f_index1,jj) .eq. 6.1)) then
					!frac6(jj,1) = -param_f_por*param_f_dx*(param_f_dx*param_f_dx*rho_fluid*grav/((y_top-y_base)*viscosity*12.0))*(rho_one(h_top)*y_top - rho_one(h_base)*y_base - rho_fluid*(y_top-y_base))
					frac6(jj,1) = -param_f_por*param_f_dx*(param_f_dx*param_f_dx*rho_one(h_base)*grav/((y_top-y_base)*viscosity*12.0))*(rho_one(h_top)*(y_top-y_top) - rho_one(h_base)*(y_base-y_top) - rho_one(param_tsw)*(y_top-y_base))
					!frac6(jj,1) = param_f_por*psi(f_index1-1,jj) + (dx*(12.0*(1.0e-16)*1.0)/(param_f_dx*param_f_dx*param_f_dx))

! 					write(*,*) "full frac6 jj 1"
! 					write(*,*) frac6(jj,1)
!
! 					write(*,*) "extra bit"
! 					write(*,*) (dx*(12.0*(1.0e-16)*1.0)/(param_f_dx*param_f_dx*param_f_dx))
! 					write(*,*) " "

					!# temporary frac6 prescription
					frac6_fix = param_myr_fix * param_h * rho_fluid / ( (3.14e7) * phi(1,1) )
					if (j .gt. tn/10) then
						frac6(jj,1) = frac6_fix !0.30
					end if

				end if
			end do
			!write(*,*) " "




			psi = psi_next(h, rhs0, psi, rho, phi, permeability, outerBand, permx, permy, j/mstep,frac6)
			psi = psi_bc(psi)

			do jj=yn/2,yn-1
				do i=1,xn
					if ((maskP(i,jj) .eq. 50.0) .and. (i .lt. f_index1)) then
						psi(i,jj+1) = maxval(frac6(:,1))
					end if
! 					if ((maskP(i,jj) .eq. 3.5)) then
! 						psi(i,jj+1) = maxval(frac6(:,1))
! 					end if
				end do
			end do


		end if


		! run with frac6 = 0 during spinup phase
		if (iter .eq. 0) then

			psi = psi_next(h, rhs0, psi, rho, phi, permeability, outerBand, permx, permy, j/mstep,frac6)
			psi = psi_bc(psi)

		end if






end if ! end if restart .ne. 1

			! get velocities from streamfunction
			velocities0 = velocities(psi)
			u = phi*velocities0(1:xn,1:yn)/(rho_fluid)!*phi
			v = phi*velocities0(1:xn,yn+1:2*yn)/(rho_fluid)!*phi

! 			do ii=1,yn
! 				do i=2,xn-2
! 					if ((maskP(i,ii) .eq. 6.0)) then
! 						v(i,ii) = -1.0*(frac6(ii,2) - frac6(ii,1))/(rho(i,ii)*dx)
! 					end if
! 				end do
! 			end do

! 			do jj=2,yn-1
! 				do i=1,xn
! 					if ((maskP(i,jj) .eq. 2.5)) then
! 						u(i,jj) = 0.0
! 						v(i,jj) = 0.0
! 					end if
! 					if ((maskP(i,jj) .eq. 7.5)) then
! 						u(i,jj) = 0.0
! 						v(i,jj) = 0.0
! 					end if
! 				end do
! 			end do



!if ((j .ge. thresh) .and. (mod(j,mstep) .eq. 0)) then

!do i = 1,cstep
! 				iso(:,:,1) = solute_next(iso(:,:,1),u,v,1.0)
! 				iso(:,:,1) = iso(:,:,1)*exp(-(3.949e-12)*dt)
!
! 				iso(:,:,2) = solute_next(iso(:,:,2),u,v,1.0)
!end do
!
! if ((param_trace .eq. 1)) then
!
! 				isoTrace(:,:,1) = particles_next(isoTrace(:,:,1),u,v,1.0,ison,particle_sat)
!
! 				do jj = 1,cstep
! 				isoTrace(3,:,1) = isoTrace(3,:,1)*exp(-(3.949e-12)*dt/cstep)
! 				end do
!
! ! 				inertTrace(:,:) = particles_next(inertTrace,u,v,10.0,inertn,inert_sat)
! end if





! 	write(*,*) "about to start mstep loop"

	! things only done every mth timestep go here
	if (mod(j,mstep) .eq. 0) then

! 		! OUTER BAND THING
! 		outerBand = make_band(permeability,phi)
! 		outerBand = band(outerBand,2*(yn-2) + 1,(xn-2)*(yn-2))
! 		permx = partial_edge((phi/(permeability)),xn,yn,dx,dy,1)
! 		permy = partial_edge((phi/(permeability)),xn,yn,dx,dy,2)

		! make coarse grid average velocities
! 		uTransport = (uTransport/(mstep*wscale))
! 		vTransport = (vTransport/(mstep*wscale))
		!uTransport(1,:) = 0.0
		!uTransport(xn/cell,:) = 0.0

! 		write(*,*) "SOLUTE COURANT NUMBERS"
! 		!write(*,*) (dt*mstep/(cstep*dx*cell))*maxval(abs(uTransport))
! 		!write(*,*) (dt*mstep/(cstep*dy*cell))*maxval(abs(vTransport))
! 		write(*,*) dt*mstep/(cstep*dy*dy*cell*cell)

!NOT RUNNING TRANSPORT RIGHT NOW, JUST CHEMISTRY



!-ISO
! u_1d = param_f_dx*param_f_dx*param_f_dx*param_f_por*grav*22.0/(viscosity*12.0*param_h)
! ! cstep_int = dx*cstep/(maxval(u/phi(1,1))*mstep)
! ! cstep_num = 6.28e10/(cstep_int*mstep)
!
! cstep_int = 6.28e11/tn
! cstep_num = 10.0*cstep_int*mstep*maxval(abs(v/phi(1,1)))/dy
!
! if (cstep_num .le. 1) then
! 	cstep_num = 1
! end if


! 		write(*,*) "starting cstep advection"
!
! !
! 		!do i = 1,cstep_num
! 		do i = 1,cstep
!
! ! 			iso(:,:,1) = solute_next(iso(:,:,1),u/phi(1,1),v/phi(1,1),1.0)
! ! 			iso(:,:,1) = iso(:,:,1)*exp(-(3.85e-12)*mstep*cstep_int/(cstep_num))
! !
! ! 			iso(:,:,2) = solute_next(iso(:,:,2),u/phi(1,1),v/phi(1,1),1.0)
!
! ! 			n=1 ! pH
! ! 			solute(:,:,n) = solute_next(solute(:,:,n),u,v,sea(n))
!   			n=2 ! alk
!   			solute(:,:,n) = solute_next(solute(:,:,n),u/phi,v/phi,sea(n))
!  			!n=3 ! water?
!  	 		!solute(:,:,n) = solute_next(solute(:,:,n),u,v,sea(n))
!  			n=4 ! c
!  	 		solute(:,:,n) = solute_next(solute(:,:,n),u/phi,v/phi,sea(n))
! 			n=5 ! ca
! 			solute(:,:,n) = solute_next(solute(:,:,n),u/phi,v/phi,sea(n))
! 			n=6 ! mg
! 	 		solute(:,:,n) = solute_next(solute(:,:,n),u/phi,v/phi,sea(n))
!  			n=7 ! na
!  	 		solute(:,:,n) = solute_next(solute(:,:,n),u/phi,v/phi,sea(n))
!  			n=8 ! k
!  	 		solute(:,:,n) = solute_next(solute(:,:,n),u/phi,v/phi,sea(n))
!  			n=9 ! fe
!  	 		solute(:,:,n) = solute_next(solute(:,:,n),u/phi,v/phi,sea(n))
!  			n=10 ! s
!  	 		solute(:,:,n) = solute_next(solute(:,:,n),u/phi,v/phi,sea(n))
!  			n=11 ! si
!  	 		solute(:,:,n) = solute_next(solute(:,:,n),u/phi,v/phi,sea(n))
!  			n=12 ! cl
!  	 		solute(:,:,n) = solute_next(solute(:,:,n),u/phi,v/phi,sea(n))
!  			n=13 ! al
!  	 		solute(:,:,n) = solute_next(solute(:,:,n),u/phi,v/phi,sea(n))
! !  			n=14 ! inert
! !  	 		solute(:,:,n) = solute_next(solute(:,:,n),u/phi,v/phi,0.1)
! ! 			do jj = 1,yn/cell
! ! 			do ii = 1,xn/cell
! ! 				solute(ii,jj,13) = max(solute(ii,jj,13),1.0e-8)
! ! 			end do
! ! 			end do
! 		end do
! 		write(*,*) "done with cstep advection"
!
!
! 		write(*,*) "about to stretch everything for transport to nodes"


		! stretch everything out






		! add timestep's output to output arrays
		if (mod(j,mstep*ar) .eq. 0) then
			 rhsmat(1+xn*(j/(mstep*ar)-1):xn*(j/(mstep*ar)),1:yn) = rhs0
			 rhomat(1+xn*(j/(mstep*ar)-1):xn*(j/(mstep*ar)),1:yn) = rho
			 hmat(1+xn*(j/(mstep*ar)-1):xn*(j/(mstep*ar)),1:yn) = h
			 psimat(1+xn*(j/(mstep*ar)-1):xn*(j/(mstep*ar)),1:yn) = psi
			 umat(1+xn*(j/(mstep*ar)-1):xn*(j/(mstep*ar)),1:yn) = u
			 vmat(1+xn*(j/(mstep*ar)-1):xn*(j/(mstep*ar)),1:yn) = v
			 permmat(1+xn*(j/(mstep*ar)-1):xn*(j/(mstep*ar)),1:yn) = permeability
			 permxMat(1+xn*(j/(mstep*ar)-1):xn*(j/(mstep*ar)),1:yn) = permx
			 permyMat(1+xn*(j/(mstep*ar)-1):xn*(j/(mstep*ar)),1:yn) = permy


			!  psiCoarseMat(1+(xn/cell)*(j/(mstep*ar)-1):(xn/cell)*(j/(mstep*ar)),1:yn/cell) = psiCoarse
			!  uCoarseMat(1+(xn/cell)*(j/(mstep*ar)-1):(xn/cell)*(j/(mstep*ar)),1:yn/cell) = uTransport
			!  vCoarseMat(1+(xn/cell)*(j/(mstep*ar)-1):(xn/cell)*(j/(mstep*ar)),1:yn/cell) = vTransport
	 		! ! reset coarse grid velocities for next timestep
	 		! uTransport = 0.0
	 		! vTransport = 0.0

! 			 primaryMat(1+(xn/cell)*(j/(mstep*ar)-1):(xn/cell)*(j/(mstep*ar)),1:yn/cell,:) = primary
! 			 secondaryMat(1+(xn/cell)*(j/(mstep*ar)-1):(xn/cell)*(j/(mstep*ar)),1:yn/cell,:)= secondary
! 			 soluteMat(1+(xn/cell)*(j/(mstep*ar)-1):(xn/cell)*(j/(mstep*ar)),1:yn/cell,:) = solute
! 			 mediumMat(1+(xn/cell)*(j/(mstep*ar)-1):(xn/cell)*(j/(mstep*ar)),1:yn/cell,:) = medium
! ! 			 isoMat(1+(xn/cell)*(j/(mstep*ar)-1):(xn/cell)*(j/(mstep*ar)),1:yn/cell,:) = iso
! 			 saturationMat(1+(xn/cell)*(j/(mstep*ar)-1):(xn/cell)*(j/(mstep*ar)),1:yn/cell,:) = saturation

! 			 isoTraceMat(:,1+ison*(j/(mstep*ar)-1):ison*(j/(mstep*ar)),:) = isoTrace
! 			 inertTraceMat(:,1+inertn*(j/(mstep*ar)-1):inertn*(j/(mstep*ar))) = inertTrace




!
! !-WRITE CHECKPOINT TO FILE
!
! !yep = write_vec ( 1, real(crashstep,kind=4), trim(path) // 'checkpoint/crashstep.txt' )
!
! yep = write_matrix ( xn, yn, real(permeability,kind=4), trim(path) // 'checkpoint/permeability.txt' )



if (restart .ne. 1) then



! yep = write_matrix ( xn, yn, real(rhs0,kind=4), trim(path) // 'rhs.txt' )
! yep = write_matrix ( xn, yn, real(psi,kind=4), trim(path) // 'psi.txt' )
yep = write_matrix ( xn, yn, real(h,kind=4), trim(path) // 'h.txt' )
!
! yep = write_matrix ( xn, yn, real(u,kind=4), trim(path) // 'u.txt' )
! yep = write_matrix ( xn, yn, real(v,kind=4), trim(path) // 'v.txt' )








!-WRITE EVERYTHING TO FILE



yep = write_matrix ( yn, 2, real(frac6, kind = 4), trim(path) // 'frac6.txt' )
yep = write_matrix ( yn, 2, real(temp6, kind = 4), trim(path) // 'temp6.txt' )
! yep = write_matrix ( xn, yn/2, real(mask(:,(yn/2)+1:), kind = 4), trim(path) // 'mask.txt' )
yep = write_matrix ( xn, yn/2, real(maskP(:,(yn/2)+1:), kind = 4), trim(path) // 'maskP.txt' )
yep = write_vec ( xn, real(x,kind=4), trim(path) // 'x.txt' )
yep = write_vec ( yn/2, real(y(yn/2:),kind=4), trim(path) // 'y.txt' )
yep = write_matrix (xn*tn/(mstep*ar), yn/2, real(hmat(:,(yn/2)+1:), kind = 4), trim(path) // 'hMat.txt' )
yep = write_matrix (xn*tn/(mstep*ar), yn/2, real(psimat(:,(yn/2)+1:),kind=4), trim(path) // 'psiMat.txt' )
! yep = write_matrix (xn*tn/(mstep*ar), yn/2, real(rhomat(:,(yn/2)+1:),kind=4), trim(path) // 'rhoMat.txt' )
! yep = write_matrix (xn*tn/(mstep*ar), yn/2, real(umat(:,(yn/2)+1:), kind = 4), trim(path) // 'uMat.txt' )
! yep = write_matrix (xn*tn/(mstep*ar), yn/2, real(vmat(:,(yn/2)+1:),kind=4), trim(path) // 'vMat.txt' )
! yep = write_matrix ( xn, yn/2,real(lambdaMat(:,(yn/2)+1:),kind=4), trim(path) // 'lambdaMat.txt' )
yep = write_matrix ( xn, yn/2,real(permeability(:,(yn/2)+1:),kind=4), trim(path) // 'permeability.txt' )
! yep = write_matrix ( xn, yn/2,real(phi(:,(yn/2)+1:),kind=4), trim(path) // 'phi.txt' )
!yep = write_matrix ( xn/cell, yn/cell, real(phiCoarse,kind=4), trim(path) // 'phiCoarse.txt' )

!if (maxval(medium(:,:,5)) .eq. 1.0) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SAVING TIME !!!!!!!!!!!!!!!!!!!!!!!!

end if ! end write only if restart ne 1




	 	end if
		! end mstep*ar loop




end if
! end mstep timestep loop, finally


end do
! end all timestep loop




write(*,*) " "
write(*,*) "ALL DONE!"
write(*,*) tn
write(*,*) "steps"
write(*,*) tn/mstep
write(*,*) "msteps"




! what to do if you are a slave processor
else









! end loop through processors
end if




! close up shop
call MPI_FINALIZE ( ierr )


END PROGRAM main



! ----------------------------------------------------------------------------------%%
!
! H_NEXT
!
! SUMMARY: computes the 2D temperature profile for the current timestep
!
! INPUTS: h(xn,yn) : temperature profile of previous timestep
!         psi(xn,yn) : 2D streamfunction array
!         rho_in(xn,yn) : 2D density array
!         flux(xn,2) : top and bottom heat flux boundary conditions
!
! RETURNS: h_next(xn,yn) : temperature profile of current timestep
!
! ----------------------------------------------------------------------------------%%

function h_next (h, psi, rho_in, phi_in, u_in, v_in, frac6_in, temp6_in, dt_in)

use globals
use initialize
implicit none

! interface
!
! 	function velocities(psi)
! 		use globals
! 		use initialize
! 		implicit none
! 		real(4) :: velocities(xn,2*yn), psi(xn,yn)
! 		real(4) :: u0(xn,yn), v0(xn,yn)
! 	end function velocities
!
! end interface

! declare errthing

! integers
integer :: i, j, n, ii, m=3
! inputs
real(4) :: sx, sy, qx, qy, rho_in(xn,yn), flux(xn,2), phi_in(xn,yn)
! velocity stuff
real(4) :: uf(xn,yn), vf(xn,yn), u_in(xn,yn), v_in(xn,yn)
real(4) :: u(xn,yn), v(xn,yn), uLong((xn-2)*(yn-2)), vLong((xn-2)*(yn-2))
real(4) ::  velocities0(xn,2*yn)
! matrix stuff
real(4) :: h(xn,yn), h_next(xn,yn), psi(xn,yn)
! real(4) :: aa((xn-2)*(yn-2),(xn-2)*(yn-2)), a((xn-2)*(yn-2),(xn-2)*(yn-2)+1)
real(4) :: aBand((xn-2)*(yn-2),5), bBand((xn-2)*(yn-2),5)
! real(4) :: bb((xn-2)*(yn-2),(xn-2)*(yn-2)), b((xn-2)*(yn-2),(xn-2)*(yn-2)+1)
real(4) :: h0(xn,yn), hMid(xn,yn), uVec((xn-2)*(yn-2)), h_nextRow((xn-2)*(yn-2))
real(4) :: kMatLong((xn-2)*(yn-2))
real(4) :: mn(xn,yn)
real(4) :: sxMat(xn,yn), syMat(xn,yn), sxLong((xn-2)*(yn-2)), syLong((xn-2)*(yn-2))
real(4) :: qxMat(xn,yn), qyMat(xn,yn), qxLong((xn-2)*(yn-2)), qyLong((xn-2)*(yn-2))
real(4) :: frac6_in(yn,2), temp6_in(yn,2), dt_in


! ! calculate velocities from streamfunction values
! velocities0 = velocities(psi)
! uf = velocities0(1:xn,1:yn)
! vf = velocities0(1:xn,yn+1:2*yn)
u = -1.0*u_in
v = -1.0*v_in
! uf = -1.0*uf!*rho_in
! vf = -1.0*vf!*rho_in


! do ii=2,yn-1
! 	do i=1,xn
! 		if ((maskP(i,ii) .eq. 2.5)) then
! 			u(i,ii) = 0.0
! 			v(i,ii) = 0.0
! 		end if
! 		if ((maskP(i,ii) .eq. 7.5)) then
! 			u(i,ii) = 0.0
! 			v(i,ii) = 0.0
! 		end if
! 	end do
! end do

! do ii=2,yn
! 	do i=1,xn
! 		if ((maskP(i,ii) .ne. 0.0) .and. (maskP(i,ii-1) .ne. 0.0)) then
! 			u(i,ii) = (uf(i,ii) + uf(i,ii-1))/(2.0)
! 		end if
! 	end do
! end do
!
! do ii=1,yn
! 	do i=2,xn-2
! 		if ((maskP(i,ii) .ne. 0.0) .and. (maskP(i-1,ii) .ne. 0.0)) then
! 			v(i,ii) = (vf(i,ii) + vf(i-1,ii))/(2.0)
! 		end if
! 	end do
! end do
!
!
!
! do ii=1,yn
! 	do i=2,xn-2
! 		if ((maskP(i,ii) .eq. 9.0)) then
! 			v(i,ii) = -1.0*phi_in(i,ii)*(psi(i+1,ii) - psi(i,ii))/dx
! 		end if
! 		if ((maskP(i,ii) .eq. 3.0)) then
! 			v(i,ii) = -1.0*phi_in(i,ii)*(psi(i,ii) - psi(i-1,ii))/dx
! 		end if
! 	end do
! end do
!
! do ii=1,yn
! 	do i=2,xn-2
! 		if ((maskP(i,ii) .eq. 6.0)) then
! 			v(i,ii) = (frac6_in(ii,2) - frac6_in(ii,1))/(2.0*param_f_dx)
! 		end if
! 	end do
! end do

! do ii=1,yn
! 	do i=2,xn-2
! 		if ((maskP(i,ii) .eq. 6.0)) then
! 			v(i,ii) = -1.0*(frac6_in(ii,2) - frac6_in(ii,1))/(rho_in(i,ii)*2.0*param_f_dx)
! 		end if
! 	end do
! end do

! ! fracture vertical velocity
! v(xn-1,:) = v(xn-1,:)*param_f_dx/dx






uLong = -1.0*reshape(u(2:xn-1,2:yn-1), (/(xn-2)*(yn-2)/))
vLong = -1.0*reshape(transpose(v(2:xn-1,2:yn-1)), (/(xn-2)*(yn-2)/))


mn = h



h0 = h
hMid = h

qx = dt/(dx)
qy = dt/(dy)
sx = (2.0*dt_in*lambdaMat(1,1))/(dx*dx*rho_fluid*cp)
sy = (2.0*dt_in*lambdaMat(1,1))/(dy*dy*rho_fluid*cp)

! qxMat = 1.0*dt*4179.0*phi_in/(dx*(rho_in*(phi_in) + (2200.0*(1.0-phi_in)))*cp)
! qyMat = 1.0*dt*4179.0*phi_in/(dy*(rho_in*(phi_in) + (2200.0*(1.0-phi_in)))*cp)
! ! sxMat = (2.0*dt*(lambdaMat*(1.0-phi_in) + 0.6*phi_in))/(dx*dx*((rho_in*phi_in) + (2200.0*(1.0-phi_in)))*cp)
! ! syMat = (2.0*dt*(lambdaMat*(1.0-phi_in) + 0.6*phi_in))/(dy*dy*((rho_in*phi_in) + (2200.0*(1.0-phi_in)))*cp)
! sxMat = 2.0*dt*(lambdaMat)/(dx*dx*((rho_in*phi_in) + (2200.0*(1.0-phi_in)))*cp)
! syMat = 2.0*dt*(lambdaMat)/(dy*dy*((rho_in*phi_in) + (2200.0*(1.0-phi_in)))*cp)


qxMat = dt_in*4179.0/(dx*cp)
qyMat = dt_in*4179.0/(dy*cp)
sxMat = (2.0*dt_in*(lambdaMat*(1.0-phi_in) + 0.6*phi_in))/(dx*dx*rho_in*cp)
syMat = (2.0*dt_in*(lambdaMat*(1.0-phi_in) + 0.6*phi_in))/(dy*dy*rho_in*cp)


! ! print stability conditions at each timestep
! write(*,*) " "
! write(*,*) "velocity check"
! write(*,"(F10.5)") maxval(abs(u))*maxval(abs(qxMat))
! write(*,"(F10.5)") maxval(abs(v))*maxval(abs(qyMat))
! write(*,*) "conduction check"
! write(*,"(F10.5)") maxval(abs(syMat))
! write(*,*) " "
!
! write(*,*) "velocity check"
! write(*,*) maxval(abs(u))
! write(*,*) maxval(abs(v))

! vertical boundary conditions

! stretch
stretch = h0
stretch(2,:) = stretch(1,:)
stretch(xn-1,:) = stretch(xn,:)
do ii=2,yn-2
	do i=2,xn-1
		if ((mask(i,ii) .eq. 5.0) .or. (mask(i,ii) .eq. 12.5) ) then
			stretch(i,ii) = stretch(i+1,ii)
		end if

		if ((mask(i,ii) .eq. 3.0) .or. (mask(i,ii) .eq. 3.5) .or. (mask(i,ii) .eq. 3.1) .or. (mask(i,ii) .eq. 3.05)) then
			stretch(i,ii) = temp6_in(ii,1)
		end if

		if ((mask(i,ii) .eq. 10.0) .or. (mask(i,ii) .eq. 17.5)) then
			stretch(i,ii) = stretch(i-1,ii)
		end if

		if ((mask(i,ii) .eq. 6.0) .or. (mask(i,ii) .eq. 6.5) .or. (mask(i,ii) .eq. 6.1) .or. (mask(i,ii) .eq. 6.05)) then
			stretch(i,ii) = temp6_in(ii,2)
		end if
	end do
end do

stretchLong = reshape(stretch(2:xn-1,2:yn-1), (/(xn-2)*(yn-2)/))

do ii=2,yn-1
do i=2,xn-1
	! right outcrop (left boundary)
	if ((mask(i,ii) .eq. 17.5)) then
		h(i,ii) = h(i,ii) + h0(i-1,ii)*sxMat(i,ii)/2.0
	end if
	if ((mask(i,ii) .eq. 10.0)) then
		h(i,ii) = h(i,ii) + h0(i-1,ii)*sxMat(i,ii)/2.0
	end if
	if ((mask(i,ii) .eq. 6.0) .or. (mask(i,ii) .eq. 6.5) .or. (mask(i,ii) .eq. 6.1) .or. (mask(i,ii) .eq. 6.05)) then
		h(i,ii) = h(i,ii) + temp6_in(ii,2)*sxMat(i,ii)/2.0
	end if

	! left outcrop (right boundary)
	if ((mask(i,ii) .eq. 12.5)) then
		h(i,ii) = h(i,ii) + h0(i+1,ii)*sxMat(i,ii)/2.0
	end if
	if ((mask(i,ii) .eq. 5.0)) then
		h(i,ii) = h(i,ii) + h0(i+1,ii)*sxMat(i,ii)/2.0
	end if
	if ((mask(i,ii) .eq. 3.0) .or. (mask(i,ii) .eq. 3.5) .or. (mask(i,ii) .eq. 3.1) .or. (mask(i,ii) .eq. 3.05)) then
		h(i,ii) = h(i,ii) + temp6_in(ii,1)*sxMat(i,ii)/2.0
	end if



end do

if (mask(1,ii) .eq. 1.0) then
	h(2,ii) = h(2,ii) + h0(1,ii)*sxMat(2,ii)/2.0  ! left
end if
if (mask(1,ii) .eq. 25.0) then
	h(2,ii) = h(2,ii) + h0(1,ii)*sxMat(2,ii)/2.0  ! top left corner
end if

if (mask(xn,ii) .eq. 1.0) then
	h(xn-1,ii) = h(xn-1,ii) + h0(xn,ii)*sxMat(xn-1,ii)/2.0  ! right
end if
if (mask(xn,ii) .eq. 50.0) then
	h(xn-1,ii) = h(xn-1,ii) + h0(xn,ii)*sxMat(xn-1,ii)/2.0  ! top right corner
end if
end do



!h(2,2) = h(2,2) + h0(1,2)*sxMat(2,2)/2.0  ! bottom left corner
!h(xn-1,2) = h(xn-1,2) + h0(xn,2)*sxMat(xn-1,2)/2.0  ! bottom right corner

uVec = reshape(h(2:xn-1,2:yn-1), (/(xn-2)*(yn-2)/))
h0Long = reshape(h0(2:xn-1,2:yn-1), (/(xn-2)*(yn-2)/))
sxLong = reshape(sxMat(2:xn-1,2:yn-1), (/(xn-2)*(yn-2)/))
syLong = reshape(syMat(2:xn-1,2:yn-1), (/(xn-2)*(yn-2)/))
qxLong = reshape(qxMat(2:xn-1,2:yn-1), (/(xn-2)*(yn-2)/))
qyLong = reshape(qyMat(2:xn-1,2:yn-1), (/(xn-2)*(yn-2)/))

! make the band
aBand = 0.0
aBand(:,1) = -sxLong/2.0 - uLong*qxLong/2.0
aBand(:,2) = 1.0+sxLong
aBand(:,3) = -sxLong/2.0 + uLong*qxLong/2.0

! bottom left corner, 2 and 3
aBand(1,1) =  0.0
aBand(1,2) = 1.0 + sxLong(1)/1.0 - uLong(1)*qxLong(1)!/2.0
aBand(1,3) = -sxLong(1)/2.0 + uLong(1)*qxLong(1)!/2.0

! top right corner, 1 and 2
aBand((xn-2)*(yn-2),1) = -sxLong((xn-2)*(yn-2))/2.0 - uLong((xn-2)*(yn-2))*qxLong((xn-2)*(yn-2))!/2.0
aBand((xn-2)*(yn-2),2) = 1.0 + sxLong((xn-2)*(yn-2))/1.0 + uLong((xn-2)*(yn-2))*qxLong((xn-2)*(yn-2))!/2.0
aBand((xn-2)*(yn-2),3) =  0.0



do i = 2,(xn-2)*(yn-2)-1

	! flow left anywhere, 2 and 3
	if (uLong(i) .lt. 0.0) then

		aBand(i,1) = -sxLong(i)/2.0
		aBand(i,2) = 1.0+sxLong(i) - uLong(i)*qxLong(i)
		aBand(i,3) = -sxLong(i)/2.0 + uLong(i)*qxLong(i)

	end if


	! flow right anywhere, 1 and 2
	if (uLong(i) .gt. 0.0) then

		aBand(i,1) = -sxLong(i)/2.0 - uLong(i)*qxLong(i)
		aBand(i,2) = 1.0+sxLong(i) + uLong(i)*qxLong(i)
		aBand(i,3) = -sxLong(i)/2.0

	end if



	! left edges, default 2 and 3
	if (((mod(i-1,xn-2).eq.0)) .or. (maskLong(i).eq.10.0) .or. (maskLong(i).eq.17.5) .or. (maskLong(i).eq.6.0) .or. (maskLong(i).eq.6.5) .or. (maskLong(i).eq.6.1) .or. (maskLong(i).eq.6.05)) then
			aBand(i,1) =  0.0
			aBand(i,2) = 1.0 + sxLong(i) - uLong(i)*qxLong(i)
			aBand(i,3) = -sxLong(i)/2.0 + uLong(i)*qxLong(i)
	end if



	! left edge flowing to the right, 1 & 2
	if (uLong(i) .gt. 0.0) then

		! left edge of right outcrop
		if ((maskLong(i).eq.10.0) .or. (maskLong(i).eq.17.5) .or. (maskLong(i).eq.6.0) .or. (maskLong(i).eq.6.5) .or. (maskLong(i).eq.6.1) .or. (maskLong(i).eq.6.05)) then
				aBand(i,1) =  0.0
				aBand(i,2) = 1.0 + sxLong(i) + uLong(i)*qxLong(i)
				aBand(i,3) = -sxLong(i)/2.0
				uVec(i) = uVec(i) + uLong(i)*qxLong(i)*stretchLong(i)
		end if

		! left edge but not uppper left corner
		if ((mod(i-1,xn-2).eq.0) .and. (maskLong(i) .ne. 25.0)) then
				aBand(i,1) =  0.0
				aBand(i,2) = 1.0 + sxLong(i) + uLong(i)*qxLong(i)
				aBand(i,3) = -sxLong(i)/2.0
				uVec(i) = uVec(i) + uLong(i)*qxLong(i)*stretchLong(i)
		end if


	end if


	! right edge, 1 and 2 by default
	if ((mod(i,xn-2) .eq. 0) .or. (maskLong(i).eq.5.0) .or. (maskLong(i).eq.12.5) .or. (maskLong(i).eq.3.0) .or. (maskLong(i).eq.3.5) .or. (maskLong(i).eq.3.1) .or. (maskLong(i).eq.3.05)) then
			aBand(i,1) = -sxLong(i)/2.0 - uLong(i)*qxLong(i)
			aBand(i,2) = 1.0 + sxLong(i) + uLong(i)*qxLong(i)
			aBand(i,3) =  0.0
	end if


	! right edge flowing to the left, 2 and 3
	if (uLong(i) .lt. 0.0) then

		! right edge of left outcrop
		if ((maskLong(i).eq.5.0) .or. (maskLong(i).eq.12.5) .or. (maskLong(i).eq.3.0) .or. (maskLong(i).eq.3.5) .or. (maskLong(i).eq.3.1) .or. (maskLong(i).eq.3.05)) then
				aBand(i,1) = -sxLong(i)/2.0
				aBand(i,2) = 1.0 + sxLong(i) - uLong(i)*qxLong(i)
				aBand(i,3) =  0.0
				uVec(i) = uVec(i) - uLong(i)*qxLong(i)*stretchLong(i)
		end if

		! right edge but not upper right corner
		if ((mod(i,xn-2) .eq. 0) .and. (maskLong(i) .ne. 25.0)) then
				aBand(i,1) = -sxLong(i)/2.0
				aBand(i,2) = 1.0 + sxLong(i) - uLong(i)*qxLong(i)
				aBand(i,3) =  0.0
				uVec(i) = uVec(i) - uLong(i)*qxLong(i)*stretchLong(i)
		end if


	end if


	! upper left corner, default 2 and 3
	if ((mod(i-1,xn-2).eq.0) .and. (maskLong(i) .eq. 25.0)) then
			aBand(i,1) =  0.0
			aBand(i,2) = 1.0 + sxLong(i)/1.0 - uLong(i)*qxLong(i)!/2.0
			aBand(i,3) = -sxLong(i)/2.0 + uLong(i)*qxLong(i)!/2.0
	end if

	! upper right corner, default 1 and 2
	if ((mod(i,xn-2).eq.0) .and. (maskLong(i) .eq. 25.0)) then
			aBand(i,1) = -sxLong(i)/2.0 - uLong(i)*qxLong(i)!/2.0
			aBand(i,2) = 1.0 + sxLong(i)/1.0 + uLong(i)*qxLong(i)!/2.0
			aBand(i,3) =  0.0
	end if

	! bottom right corner, 1 and 2
	if (i.eq.xn-2) then
			aBand(i,1) = -sxLong(i)/2.0 - uLong(i)*qxLong(i)!/2.0
			aBand(i,2) = 1.0 + sxLong(i)/1.0 + uLong(i)*qxLong(i)!/2.0
			aBand(i,3) =  0.0
	end if






!

end do

do i=1,(xn-2)*(yn-2)
	! mask
	if ((maskLong(i) .eq. 0.0) .or. (maskLong(i) .eq. 600.0)) then
		aBand(i,2) = 1.0
		aBand(i,1) = 0.0
		aBand(i,3) = 0.0
	end if

end do

! make sure solver doesn't go out of bounds
do i=1,((yn-2)-1)
	ii = i*(xn-2)
	aBand(ii,3) = 0.0
	aBand(ii+1,1) = 0.0
end do

!!!!!!!!!!!! THIS !!!!!!!!!!!
h_nextRow = tridiag(aBand(:,1),aBand(:,2),aBand(:,3),uVec,(xn-2)*(yn-2))
h(2:xn-1,2:yn-1) = reshape(h_nextRow, (/xn-2, yn-2/))




stretch = h0
stretch(:,2) = stretch(:,1)
do ii=2,yn-1
	do i=2,xn-1
		if ((mask(i,ii) .eq. 25.0) .or. (mask(i,ii) .eq. 50.0)  .or. (mask(i,ii) .eq. 17.5)  .or. (mask(i,ii) .eq. 12.5)  .or. (mask(i,ii) .eq. 3.5)  .or. (mask(i,ii) .eq. 6.5)) then
			stretch(i,ii) = stretch(i,ii+1)
		end if

	end do
end do

stretchLongT = reshape(transpose(stretch(2:xn-1,2:yn-1)), (/(xn-2)*(yn-2)/))




!h0LongT = reshape(transpose(h0(2:xn-1,2:yn-1)), (/(xn-2)*(yn-2)/))
!h0T = transpose(h0)
sxLong = reshape(transpose(sxMat(2:xn-1,2:yn-1)), (/(xn-2)*(yn-2)/))
syLong = reshape(transpose(syMat(2:xn-1,2:yn-1)), (/(xn-2)*(yn-2)/))
qxLong = reshape(transpose(qxMat(2:xn-1,2:yn-1)), (/(xn-2)*(yn-2)/))
qyLong = reshape(transpose(qyMat(2:xn-1,2:yn-1)), (/(xn-2)*(yn-2)/))

! horizontal boundary conditions
h(3:xn-2,2) = h(3:xn-2,2) + h0(3:xn-2,1)*syMat(3:xn-2,2)/2.0 ! bottom


! top of sediment
do ii=2,yn-1
do i=2,xn-1
	! top of sediment and any short outcrops
	if ((mask(i,ii) .eq. 50.0) .or. (mask(i,ii) .eq. 2.0) .or. (mask(i,ii) .eq. 7.0)) then
		h(i,ii) = h(i,ii) + h0(i,ii+1)*(syMat(i,ii)/2.0) ! top
	end if

	if ((mask(i,ii) .eq. 25.0) .and. (i .ge. 3) .and. (i .le. xn-2)) then
		h(i,ii) = h(i,ii) + h0(i,ii+1)*(syMat(i,ii)/2.0) ! top
	end if
	if ((mask(i,ii) .eq. 25.0) .and. (i .eq. 2)) then
		h(i,ii) = h(i,ii) + h0(i,ii+1)*(syMat(i,ii)/2.0) ! top left corner
	end if
	if ((mask(i,ii) .eq. 25.0) .and. (i .eq. xn-1)) then
		h(i,ii) = h(i,ii) + h0(i,ii+1)*(syMat(i,ii)/2.0) ! top right corner
	end if

	if ((mask(i,ii) .eq. 1.0) .and. (ii .eq. 2) .and. (i .eq. 2)) then
		h(i,ii) = h(i,ii) + h0(i,ii-1)*(syMat(i,ii)/2.0) ! bottom left corner
	end if
	if ((mask(i,ii) .eq. 1.0) .and. (ii .eq. 2) .and. (i .eq. xn-1)) then
		h(i,ii) = h(i,ii) + h0(i,ii-1)*(syMat(i,ii)/2.0) ! bottom right corner
	end if


	if ((mask(i,ii) .eq. 12.5) .or. (mask(i,ii) .eq. 17.5) .or. (mask(i,ii) .eq. 3.5) .or. (mask(i,ii) .eq. 6.5)) then
		h(i,ii) = h(i,ii) + h0(i,ii+1)*(syMat(i,ii)/2.0) ! top
	end if


end do
end do

!genTrans = transpose(h)
!h_nextRow = reshape(genTrans(2:yn-1,2:xn-1), (/(xn-2)*(yn-2)/))
h_nextRow = reshape(transpose(h(2:xn-1,2:yn-1)), (/(xn-2)*(yn-2)/))

! make the band
bBand = 0.0
bBand(:,1) = -syLong/2.0 - vLong(:)*qyLong/2.0
bBand(:,2) = 1.0+syLong
bBand(:,3) = -syLong/2.0 + vLong*qyLong/2.0

! bottom left corner
bBand(1,1) =  0.0
bBand(1,2) = 1.0 + syLong(1)/1.0 - vLong(1)*qyLong(1)!/2.0
bBand(1,3) = -syLong(1)/2.0 + vLong(1)*qyLong(1)!/2.0

! top right corner
bBand((xn-2)*(yn-2),1) = -syLong((xn-2)*(yn-2))/2.0 - vLong((xn-2)*(yn-2))*qyLong((xn-2)*(yn-2))!/2.0
bBand((xn-2)*(yn-2),2) = 1.0 + syLong((xn-2)*(yn-2))/1.0 + vLong((xn-2)*(yn-2))*qyLong((xn-2)*(yn-2))!/2.0
bBand((xn-2)*(yn-2),3) =  0.0
do i = 2,(xn-2)*(yn-2)-1


	! flow going down anywhere, 2 and 3
	if (vLong(i) .lt. 0.0) then

		bBand(i,1) = -syLong(i)/2.0
		bBand(i,2) = 1.0+syLong(i) - vLong(i)*qyLong(i)
		bBand(i,3) = -syLong(i)/2.0 + vLong(i)*qyLong(i)

	end if

	! flow coming up anywhere, 1 and 2
	if (vLong(i) .gt. 0.0) then

		bBand(i,1) = -syLong(i)/2.0 - vLong(i)*qyLong(i)
		bBand(i,2) = 1.0+syLong(i) + vLong(i)*qyLong(i)
		bBand(i,3) = -syLong(i)/2.0

	end if


	!!!!! TOP EDGES


	! bottom rrow, default 2 and 3 !!
	if (mod(i-1,yn-2) .eq. 0) then
			bBand(i,1) =  0.0
			bBand(i,2) = 1.0 + syLong(i) - vLong(i)*qyLong(i)
			bBand(i,3) = -syLong(i)/2.0 + vLong(i)*qyLong(i)
	end if

	! bottom rrow, if flow is coming up (1 and 2)
	if (vLong(i) .gt. 0.0) then

		! bottom row but not bottom right corner
		if ((mod(i-1,yn-2) .eq. 0) .and. (i .ne. (xn-2)*(yn-2)-(yn-2)+1)) then
				bBand(i,1) =  0.0
				bBand(i,2) = 1.0 + syLong(i) + vLong(i)*qyLong(i)
				bBand(i,3) = -syLong(i)/2.0
				h_nextRow(i) = h_nextRow(i) + vLong(i)*qyLong(i)*stretchLongT(i)
		end if

	end if


	! last/top edge, default 1 and 2
	if ((maskLongT(i).eq.25.0) .or. (maskLongT(i).eq.50.0) .or. (maskLongT(i).eq.2.0) .or. (maskLongT(i).eq.7.0)) then
			bBand(i,1) = -syLong(i)/2.0 - vLong(i)*qyLong(i)
			bBand(i,2) = 1.0 + syLong(i) + vLong(i)*qyLong(i)
			bBand(i,3) =  0.0

	end if

	if ((maskLongT(i).eq.12.5).or.(maskLongT(i).eq.17.5).or.(maskLongT(i).eq.3.5).or.(maskLongT(i).eq.6.5)) then
			bBand(i,1) = -syLong(i)/2.0 - vLong(i)*qyLong(i)!/2.0
			bBand(i,2) = 1.0 + syLong(i)/1.0 + vLong(i)*qyLong(i)!/2.0
			bBand(i,3) =  0.0
	end if


	! top edges if flow coming down (2, 3)
	if (vLong(i) .lt. 0.0) then

		if ((maskLongT(i).eq.50.0) .or. (maskLongT(i).eq.2.0) .or. (maskLongT(i).eq.7.0)) then
				bBand(i,1) = -syLong(i)/2.0
				bBand(i,2) = 1.0 + syLong(i) - vLong(i)*qyLong(i)
				bBand(i,3) =  0.0
				h_nextRow(i) = h_nextRow(i) - vLong(i)*qyLong(i)*stretchLongT(i)
		end if

		! top but not top left or top right corner
		if ((maskLongT(i).eq.25.0) .and. (i .gt. yn-2) .and. (i .le. (xn-2)*(yn-2)-(yn-2))) then
				bBand(i,1) = -syLong(i)/2.0
				bBand(i,2) = 1.0 + syLong(i) - vLong(i)*qyLong(i)
				bBand(i,3) =  0.0
				h_nextRow(i) = h_nextRow(i) - vLong(i)*qyLong(i)*stretchLongT(i)

		end if

		if ((maskLongT(i).eq.12.5).or.(maskLongT(i).eq.17.5).or.(maskLongT(i).eq.3.5).or.(maskLongT(i).eq.6.5)) then
				bBand(i,1) = -syLong(i)/2.0
				bBand(i,2) = 1.0 + syLong(i)/1.0 - vLong(i)*qyLong(i)!/2.0
				bBand(i,3) =  0.0
				h_nextRow(i) = h_nextRow(i) - vLong(i)*qyLong(i)*stretchLongT(i)!/2.0
		end if

	end if




	! upper left corner
	if ((i .le. (yn-2)) .and. (maskLongT(i) .eq. 25.0)) then
			bBand(i,1) = -syLong(i)/2.0 - vLong(i)*qyLong(i)!/2.0
			bBand(i,2) = 1.0 + syLong(i)/1.0 + vLong(i)*qyLong(i)!/2.0
			bBand(i,3) =  0.0
	end if

	! upper right corner
	if ((i .gt. (yn-2)*(xn-2)-(yn-2)) .and. (maskLongT(i) .eq. 25.0)) then
			bBand(i,1) = -syLong(i)/2.0 - vLong(i)*qyLong(i)!/2.0
			bBand(i,2) = 1.0 + syLong(i)/1.0 + vLong(i)*qyLong(i)!/2.0
			bBand(i,3) =  0.0

	end if



	! bottom right corner
	if (i .eq. (yn-2)*(xn-2)-(yn-2)+1) then
			bBand(i,1) =  0.0
			bBand(i,2) = 1.0 + syLong(i)/1.0 - vLong(i)*qyLong(i)!/2.0
			bBand(i,3) = -syLong(i)/2.0 + vLong(i)*qyLong(i)!/2.0
	end if




end do

do i=1,(xn-2)*(yn-2)
	! mask
	if ((maskLongT(i) .eq. 0.0) .or. (maskLongT(i) .eq. 600.0)) then
		bBand(i,2) = 1.0
		bBand(i,1) = 0.0
		bBand(i,3) = 0.0
	end if

end do

! make sure solver doesn't go out of bounds
do i=1,((xn-2)-1)
	ii = i*(yn-2)
	bBand(ii,3) = 0.0
	bBand(ii+1,1) = 0.0
end do

h_nextRow = tridiag(bBand(:,1),bBand(:,2),bBand(:,3),h_nextRow,(xn-2)*(yn-2))
h_next(2:xn-1,2:yn-1) = transpose(reshape(h_nextRow, (/yn-2, xn-2/)))
hMid = h_next
!
!
!

! ! something awful right here
! ! possible not-doing-upwind correction?
! do ii=2,yn-1
! do i = 2,xn-1
! 	if ( (hMid(i,ii).gt.h0(i-1,ii)).and.(hMid(i,ii).gt.h0(i+1,ii)).and. &
! 	(hMid(i,ii).gt.h0(i,ii-1)).and.(hMid(i,ii).gt.h0(i,ii+1))) then
! 		h_next(i,ii) = (h0(i-1,ii) + h0(i+1,ii) + h0(i,ii-1) + h0(i,ii+1))/4.0
! 	end if
!
! 	if ( (hMid(i,ii).lt.h0(i-1,ii)).and.(hMid(i,ii).lt.h0(i+1,ii)).and. &
! 	(hMid(i,ii).lt.h0(i,ii-1)).and.(hMid(i,ii).lt.h0(i,ii+1))) then
! 		h_next(i,ii) = (h0(i-1,ii) + h0(i+1,ii) + h0(i,ii-1) + h0(i,ii+1))/4.0
! 	end if
!
!
! end do
! end do




! check out how this equation is converging to steady-state
!write(*,*) "deltaT"
!write(*,*) maxval(abs((mn(2:xn-1,2:yn-1)-h_next(2:xn-1,2:yn-1))/h_next(2:xn-1,2:yn-1)))

return

end function h_next


! ----------------------------------------------------------------------------------%%
!
! PSI_NEXT
!
! SUMMARY: computes the 2D streamfunction array of the current timestep
!
! INPUTS: h(xn,yn) : temperature profile
!         rhs0(xn,yn) : right hand side of streamfunction-vorticity equation
!         psi(xn,yn) : 2D streamfunction array of previous timestep
!         top_in(xn,1) : permeable upper boundary
!         rho_in(xn,yn) : 2D density array
!
! RETURNS: psi_next(xn,yn): 2D streamfunction array for current timestep
!
! ----------------------------------------------------------------------------------%%

function psi_next (h, rhs0, psi, rho_in, phi_in, perm_in, band_in, permx, permy, stage, frac6_in)

use globals
use initialize
implicit none

! interface
!
! 	function partial(array,rows,cols,d1,d2,dim)
! 		use globals
! 		use initialize
! 		implicit none
! 		integer :: rows, cols, dim, i, j, ii, jj
! 		real(4) :: array(rows,cols), d1, d2, d
! 		real(4) :: partial(rows,cols)
! 	end function partial
!
! end interface

! declare errthing

! integers
integer :: i, j, ii, n, m, stage
! inputs
real(4) :: rhs0(xn,yn), rhs1(xn,yn), rhsLong(longP)
real(4) :: h(xn,yn), psi(xn,yn), rho_in(xn,yn), phi_in(xn,yn), perm_in(xn,yn)
! matrix stuff
real(4) :: uVec(longP), psiLong((xn)*(yn)), psi_nextRow(longP)
real(4) :: psi_next(xn,yn)
real(4) :: mn(xn,yn)
! back to band
real(4) :: aBand0(longP,2*((yn/2)-2) + 1), band_in(longP,2*((yn/2)-2) + 1)
real(4) :: rhoLong(longP)
real(4) :: permx(xn,yn), permy(xn,yn), frac6_in(yn,2)
real(4) :: permx_left(xn,yn), permx_right(xn,yn), permy_bottom(xn,yn), permy_top(xn,yn)
real(4) :: psi_f

psi_next = 0.0

mn = psi

m = 2*((yn/2)-2) + 1

rhs1 = rhs0

rho_in = rho_fluid

!phi_in = 1.0

permx_left = phi_in / ((grav*rho_fluid/viscosity)*perm_in*rho_fluid)
permx_right = phi_in / ((grav*rho_fluid/viscosity)*perm_in*rho_fluid)
permy_bottom = phi_in / ((grav*rho_fluid/viscosity)*perm_in*rho_fluid)
permy_top = phi_in / ((grav*rho_fluid/viscosity)*perm_in*rho_fluid)
do ii=2,yn-1
do i=2,xn-1
	permx_left(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_fluid + perm_in(i-1,ii)*rho_fluid) / 2.0)
	permx_right(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_in(i,ii) + perm_in(i+1,ii)*rho_fluid) / 2.0)

	if (maskP(i,ii) .eq. 5.0) then
		permx_right(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_in(i,ii)))
	end if

	permy_bottom(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_fluid + perm_in(i,ii-1)*rho_fluid) / 2.0)
	permy_top(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_fluid + perm_in(i,ii+1)*rho_fluid) / 2.0)

	if ((maskP(i,ii) .eq. 6.0) .or. (maskP(i,ii) .eq. 6.5) .or. (maskP(i,ii) .eq. 6.1) .or. (maskP(i,ii) .eq. 6.05) .or. (maskP(i,ii) .eq. 6.01)) then
		permx_left(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*perm_in(i,ii)*rho_fluid)!phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_in(i,ii) + (param_f_dx*param_f_dx/3.0)*rho_in(i,ii)) / 2.0)!
	end if
	if ((maskP(i,ii) .eq. 3.0) .or. (maskP(i,ii) .eq. 3.5) .or. (maskP(i,ii) .eq. 3.1) .or. (maskP(i,ii) .eq. 3.05) .or. (maskP(i,ii) .eq. 3.01)) then
		permx_right(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*perm_in(i,ii)*rho_fluid)!phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_in(i,ii) + (param_f_dx*param_f_dx/3.0)*rho_in(i,ii)) / 2.0)!
	end if
end do
end do

do ii=2,yn-1
do i=2,xn-1
	if ((maskP(i,ii) .eq. 50.0) .or. (maskP(i,ii) .eq. 3.5) .or. (maskP(i,ii) .eq. 6.5)) then
		permy_top(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_fluid) / 1.0)
	end if
end do
end do


! do ii=yn/2,yn-1
! do i=2,xn-1
! 	if ((y(ii) .ge. sed(i)-dy) .and. (any(maskP(i,:) .eq. 50.0))) then
! 		permy_top(i,ii) = phi_in(i,ii) / ((1e-16 + 1e-16) / 2.0)
! 	end if
! end do
! end do

do ii=yn/2,yn
	if (maskP(1,ii) .eq. 1.0) then
	 ! left
	 rhs1(2,ii) = rhs1(2,ii) + psi(1,ii)*permx_left(2,ii)/(dx*dx)
	end if

	if (maskP(xn,ii) .eq. 1.0) then
	! right
	 rhs1(xn-1,ii) = rhs1(xn-1,ii) + psi(xn,ii)*permx_right(xn-1,ii)/(dx*dx)
	end if

	! bottom left corner
	if (maskP(1,ii) .eq. 100.0) then
	 ! left
		rhs1(2,ii) = rhs1(2,ii) + psi(2,ii-1)*permy_bottom(2,ii)/(dy*dy)
		rhs1(2,ii) = rhs1(2,ii) + psi(1,ii)*permx_left(2,ii)/(dx*dx)
	end if

	! bottom right corner
	if (maskP(xn,ii) .eq. 100.0) then
	 ! left
		rhs1(xn-1,ii) = rhs1(xn-1,ii) + psi(xn-1,ii-1)*permy_bottom(xn-1,ii)/(dy*dy)
		rhs1(xn-1,ii) = rhs1(xn-1,ii) + psi(xn,ii)*permx_right(xn-1,ii)/(dx*dx)
	end if


end do




 ! mask boundary conditions
 ! vertical boundary conditions
do ii=yn/2,yn-1
do i=2,xn-1

	if ((maskP(i,ii) .eq. 17.5)) then
	    rhs1(i,ii) = rhs1(i,ii) + psi(i-1,ii)*permx_left(i,ii)/(dx*dx)
	end if

	if ((maskP(i,ii) .eq. 10.0)) then
	    !rhs1(i,ii) = rhs1(i,ii) + psi(i-1,ii)*permx_left(i,ii)/(dx*dx)
		rhs1(i,ii) = rhs1(i,ii) + ((4.0/3.0)*psi(i,ii) - (1.0/3.0)*psi(i+1,ii))*permx_left(i,ii)/(dx*dx)
	end if

	if ((maskP(i,ii) .eq. 6.0) .or. (maskP(i,ii) .eq. 6.5) .or. (maskP(i,ii) .eq. 6.1) .or. (maskP(i,ii) .eq. 6.05) .or. (maskP(i,ii) .eq. 6.01)) then
	    rhs1(i,ii) = rhs1(i,ii) + frac6_in(ii,2)*permx_left(i,ii)/(dx*dx)
	end if

	if ((maskP(i,ii) .eq. 12.5)) then
	    rhs1(i,ii) = rhs1(i,ii) + psi(i+1,ii)*permx_right(i,ii)/(dx*dx)
	end if

	if ((maskP(i,ii) .eq. 5.0)) then
		!rhs1(i,ii) = rhs1(i,ii) + psi(i+1,ii)*permx_right(i,ii)/(dx*dx)
		rhs1(i,ii) = rhs1(i,ii) + ((4.0/3.0)*psi(i,ii) - (1.0/3.0)*psi(i-1,ii))*permx_right(i,ii)/(dx*dx)
	end if

	if ((maskP(i,ii) .eq. 3.05) .or. (maskP(i,ii) .eq. 3.01)) then
		rhs1(i,ii) = rhs1(i,ii) + frac6_in(ii,1)*permx_right(i,ii)/(dx*dx)
	end if

	if ((maskP(i,ii) .eq. 3.0) .or. (maskP(i,ii) .eq. 3.5) .or. (maskP(i,ii) .eq. 3.1)) then
		rhs1(i,ii) = rhs1(i,ii) + frac6_in(ii,1)*permx_right(i,ii)/(dx*dx)
		!rhs1(i,ii) = rhs1(i,ii) + (frac6_in(ii,1)/(dx*dx))*phi_in(i,ii)*2.0*viscosity/(grav*rho_fluid*(rho_fluid*perm_in(i,ii) + (rho_fluid*param_f_dx*param_f_dx/12.0) ))
		!rhs1(i,ii) = rhs1(i,ii) + (frac6_in(ii,1)/(dx*dx))*(phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_in(i,ii) + (param_f_dx*param_f_dx/12.0)*rho_fluid) / 2.0))
	end if


end do
end do




!
! ! at the bottom of everything
do ii=yn/2,yn-1
do i=3,xn-2
	if ((maskP(i,ii) .eq. 100.0) .or. (maskP(i,ii) .eq. 3.01) .or. (maskP(i,ii) .eq. 6.01)) then
		rhs1(i,ii) = rhs1(i,ii) + psi(i,ii-1)*permy_bottom(i,ii)/(dy*dy)
	end if
end do
end do

! do ii=yn/2,yn
! do i=2,xn-1
!
! if ((maskP(i,ii) .eq. 50.0) .or. (maskP(i,ii) .eq. 25.0) .or. (maskP(i,ii) .eq. 12.5) .or. (maskP(i,ii) .eq. 17.5)) then
! 	rhs1(i,ii) = rhs1(i,ii) + psi(i,ii+1)*permy_top(i,ii)/(dy*dy)
! end if
!
! end do
! end do

do ii=yn/2,yn
do i=2,xn-1

if ((maskP(i,ii) .eq. 25.0) .or. (maskP(i,ii) .eq. 12.5) .or. (maskP(i,ii) .eq. 17.5) .or. (maskP(i,ii) .eq. 3.5) .or. (maskP(i,ii) .eq. 6.5)) then
	rhs1(i,ii) = rhs1(i,ii) + psi(i,ii+1)*permy_top(i,ii)/(dy*dy)
end if

end do
end do

do ii=yn/2,yn
do i=2,xn-1

if ((maskP(i,ii) .eq. 50.0)) then
	rhs1(i,ii) = rhs1(i,ii) + psi(i,ii+1)*permy_top(i,ii)/(dy*dy)
end if

end do
end do

! ! CORNER 2
! do ii=yn/2,yn
! do i=2,xn-1
!
! if ((maskP(i,ii) .eq. 50.0) .and. (maskP(i-1,ii) .ne. 50.0)) then
! 	rhs1(i,ii) = rhs1(i,ii) + ((4.0/3.0)*psi(i,ii) - (1.0/3.0)*psi(i,ii-1)) *permy_top(i,ii)/(dy*dy)
! end if
!
! if ((maskP(i,ii) .eq. 50.0) .and. (maskP(i+1,ii) .ne. 50.0)) then
! 	rhs1(i,ii) = rhs1(i,ii) + ((4.0/3.0)*psi(i,ii) - (1.0/3.0)*psi(i,ii-1)) *permy_top(i,ii)/(dy*dy)
! end if
!
! end do
! end do


do ii=1,yn
do i=1,xn
	if ((maskP(i,ii) .eq. 0.0) .or. (maskP(i,ii) .eq. 600.0)) then
		rhs1(i,ii) = 0.0
	end if
end do
end do


uVec = reshape(transpose(rhs1( 2:xn-1 , (yn/2)+2:yn-1 )),(/longP/))

psi_next = 0.0





! THIS IS WHERE THE BAND IS MADE
aband0 = band_in

! use the banded solver here
psi_nextRow = solve(aBand0,uVec,2*((yn/2)-2) + 1,longP)
psi_next(2:xn-1,(yn/2)+2:yn-1) = transpose(reshape(psi_nextRow, (/(yn/2)-2, xn-2/)))




!write(*,*) "deltaPSI"
!write(*,*) maxval(abs((mn(2:xn-1,2:yn-1)-psi_next(2:xn-1,2:yn-1))/psi_next(2:xn-1,2:yn-1)))


return

end function psi_next




! ----------------------------------------------------------------------------------%%
!
! MAKE BAND
!
! SUMMARY: only make the big matrix every once in a while
!
!
! ----------------------------------------------------------------------------------%%





function make_band(perm_in,phi_in,permx,permy,rho_in)


	use globals
	use initialize
	implicit none


	interface

	function partial(array,rows,cols,d1,d2,dim)
		use globals
		use initialize
		implicit none
		integer :: rows, cols, dim, i, j, ii, jj
		real(4) :: array(rows,cols), d1, d2, d
		real(4) :: partial(rows,cols)
	end function partial

	function partial_edge(array,rows,cols,d1,d2,dim)
		use globals
		use initialize
		implicit none
		integer :: rows, cols, dim, i, j, ii, jj
		real(4) :: array(rows,cols), d1, d2, d
		real(4) :: partial_edge(rows,cols)
	end function partial_edge

	function partial_edge_p(array,rows,cols,d1,d2,dim)
		use globals
		use initialize
		implicit none
		integer :: rows, cols, dim, i, j, ii, jj
		real(4) :: array(rows,cols), d1, d2, d
		real(4) :: partial_edge_p(rows,cols)
	end function partial_edge_p

	end interface


	integer :: i, j, ii, n, m
	real(4) :: perm_in(xn,yn), phi_in(xn,yn), rho_in(xn,yn)
	real(4) :: permx(xn,yn), permy(xn,yn), permLong(longP)
	real(4) :: permxLong(longP), permyLong(longP)
	real(4) :: innerBand(longP,2*((yn/2)-2) + 1), make_band(longP,2*((yn/2)-2) + 1)
	real(4) :: permx_left(xn,yn), permx_right(xn,yn), permy_bottom(xn,yn), permy_top(xn,yn)
	real(4) :: permx_left_long(longP), permx_right_long(longP), permy_bottom_long(longP), permy_top_long(longP)
	real(4) :: perm_long(longP)

!	phi_in = 1.0

	rho_in = rho_fluid

	permx_left = phi_in / ((grav*rho_fluid/viscosity)*perm_in*rho_fluid)
	permx_right = phi_in / ((grav*rho_fluid/viscosity)*perm_in*rho_fluid)
	permy_bottom = phi_in / ((grav*rho_fluid/viscosity)*perm_in*rho_fluid)
	permy_top = phi_in / ((grav*rho_fluid/viscosity)*perm_in*rho_fluid)
	do ii=2,yn-1
	do i=2,xn-1
		permx_left(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_fluid + perm_in(i-1,ii)*rho_fluid) / 2.0)
		permx_right(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_in(i,ii) + perm_in(i+1,ii)*rho_fluid) / 2.0)

		if (maskP(i,ii) .eq. 5.0) then
			permx_right(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_in(i,ii)))
		end if

		permy_bottom(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_fluid + perm_in(i,ii-1)*rho_fluid) / 2.0)
		permy_top(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_fluid + perm_in(i,ii+1)*rho_fluid) / 2.0)

		if ((maskP(i,ii) .eq. 6.0) .or. (maskP(i,ii) .eq. 6.5) .or. (maskP(i,ii) .eq. 6.1) .or. (maskP(i,ii) .eq. 6.05) .or. (maskP(i,ii) .eq. 6.01)) then
			permx_left(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*perm_in(i,ii)*rho_fluid)!phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_in(i,ii) + (param_f_dx*param_f_dx/3.0)*rho_in(i,ii)) / 2.0)!
		end if
		if ((maskP(i,ii) .eq. 3.0) .or. (maskP(i,ii) .eq. 3.5) .or. (maskP(i,ii) .eq. 3.1) .or. (maskP(i,ii) .eq. 3.05) .or. (maskP(i,ii) .eq. 3.01)) then
			permx_right(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*perm_in(i,ii)*rho_fluid)!phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_in(i,ii) + (param_f_dx*param_f_dx/3.0)*rho_in(i,ii)) / 2.0)!
		end if
	end do
	end do

	do ii=2,yn-1
	do i=2,xn-1
		if ((maskP(i,ii) .eq. 50.0) .or. (maskP(i,ii) .eq. 3.5) .or. (maskP(i,ii) .eq. 6.5)) then
			permy_top(i,ii) = phi_in(i,ii) / ((grav*rho_fluid/viscosity)*(perm_in(i,ii)*rho_fluid) / 1.0)
		end if
	end do
	end do


	permx_left_long = reshape(transpose(permx_left(2:xn-1,(yn/2)+2:yn-1)),(/longP/))
	permx_right_long = reshape(transpose(permx_right(2:xn-1,(yn/2)+2:yn-1)),(/longP/))
	permy_bottom_long = reshape(transpose(permy_bottom(2:xn-1,(yn/2)+2:yn-1)),(/longP/))
	permy_top_long = reshape(transpose(permy_top(2:xn-1,(yn/2)+2:yn-1)),(/longP/))


! 	do ii=2,yn-1
! 	do i=2,xn-1
! 		if (maskP(i,ii) .eq. 50.0) then
! 			permy_top(:,ii) = perm_in(:,ii)! phi_in(i,ii) / ((perm_in(i,ii) + perm_in(i,ii+1)) / 2.0)
! 		end if
! 	end do
! 	end do

	! make the band
	innerBand = 0.0
	make_band = 0.0
	m = 2*((yn/2)-2) + 1


	! trial (7 days) inner band
	innerBand = 0.0
		do i = 1,longP
			! diagonal
			innerBand(i,(m+1)/2) = (permx_right_long(i) + permx_left_long(i))/(dx*dx) + (permy_top_long(i) + permy_bottom_long(i))/(dy*dy)


			! off-diagonals
			! all but left edge
			if ((i .gt. ((yn/2)-2)) .and. (maskPLongT(i) .ne. 10.0) .and. (maskPLongT(i) .ne. 17.5) .and. (maskPLongT(i) .ne. 6.0) .and. (maskPLongT(i) .ne. 6.5) .and. (maskPLongT(i) .ne. 6.1) .and. (maskPLongT(i) .ne. 6.05) .and. (maskPLongT(i) .ne. 6.01)) then
				innerBand(i,1) = -permx_left_long(i)/(dx*dx)
			end if
			! all but right edge
			if ((i .le. (longP)-((yn/2)-2)) .and. (maskPLongT(i) .ne. 5.0) .and. (maskPLongT(i) .ne. 12.5) .and. (maskPLongT(i) .ne. 3.0) .and. (maskPLongT(i) .ne. 3.5) .and. (maskPLongT(i) .ne. 3.1) .and. (maskPLongT(i) .ne. 3.05) .and. (maskPLongT(i) .ne. 3.01)) then
				innerBand(i,m) = -permx_right_long(i)/(dx*dx)
			end if


			! all but bottom
			if ((maskPLongT(i) .ne. 100.0) .or. (maskPLongT(i) .ne. 3.01) .or. (maskPLongT(i) .ne. 6.01)) then
				innerBand(i,(m+1)/2-1) = -permy_bottom_long(i)/(dy*dy)
			end if

			! all but top
			if ((maskPLongT(i) .ne. 50.0) .and. (maskPLongT(i) .ne. 25.0) .and. (maskPLongT(i) .ne. 12.5) .and. (maskPLongT(i) .ne. 17.5) .and. (maskPLongT(i) .ne. 3.5) .and. (maskPLongT(i) .ne. 6.5)) then
				innerBand(i,(m+1)/2+1) = -permy_top_long(i)/(dy*dy)
			end if

! 			if (maskPLongT(i) .eq. 2.0) then
! 				innerBand(i,m) = 0.0 ! skip right
! 				innerBand(i,(m+1)/2+1) = 0.0 ! skip top
! 			end if
!
! 			if (maskPLongT(i) .eq. 7.0) then
! 				innerBand(i,1) = 0.0 ! skip left
! 				innerBand(i,(m+1)/2+1) = 0.0 ! skip top
! 			end if

		end do

	do i = 1,longP

		! mask
		if ((maskPLongT(i) .eq. 0.0) .or. (maskPLongT(i) .eq. 600.0)) then
			innerBand(i,:) = 0.0
			innerBand(i,(m+1)/2) = 1.0!(2.0)/(permLong(i)*dx*dx) + (2.0)/(permLong(i)*dy*dy)
		end if

	end do

	!write(*,*) innerBand(1,:)

	make_band = innerBand

return
end function make_band














! ! ----------------------------------------------------------------------------------%%
! !
! ! ALT_NEXT
! !
! ! SUMMARY: solves for equilibrium at a single grid cell using PHREEQC
! !
! ! INPUTS: temp : temperature of grid cell
! !         timestep : time elapsed
! !         primaryList(5) : amounts of primary minerals
! !         secondaryList(16) : amounts of secondary minerals
! !         soluteList(11) : concentrations of solutes
! !
! ! RETURNS: alt_next(1,altnum): returns everything from PHREEQC in a big pile
! !          and it gets parsed in the main method's geochem loop
! !
! ! ----------------------------------------------------------------------------------%%
!
! function alt_next (temp, timestep, primaryList, secondaryList, soluteList, mediumList)
! use globals
! use initialize
! use alteration
! implicit none
!
! interface
!
! end interface
!
! ! declare errthing
!
! integer :: order
! real(4) :: temp, timestep
! real(4) :: alt_next(1,167)
! real(4) :: alter0(1,167)
! real(4) :: primaryList(g_pri), secondaryList(g_sec), soluteList(g_sol), mediumList(g_med)
!
! ! use the alteration module
! alter0 = alter(temp-272.9, timestep, primaryList, secondaryList, soluteList, mediumList)
!
! ! rename it for a reason that i now forget
! alt_next = alter0
!
! end function alt_next


! ----------------------------------------------------------------------------------%%
!
! RHO_NEXT
!
! SUMMARY : solves for density using linear thermally expansive equation of state
!
! INPUTS : h_in(xn,yn) : 2D temperature array of current timestep
!
! RETURNS : rho_next(xn,yn) : 2D density array of current timestep
!
! ----------------------------------------------------------------------------------%%

function rho_next(h_in)

use globals
use initialize
implicit none

! declare errthing
integer :: i,j
real(4) :: h_in(xn,yn), rho_next(xn,yn)


do i=1,xn
	do j = 1,yn
		rho_next(i,j) = rho_fluid*(1.0 - alpha*(h_in(i,j)-273.0))
	end do
end do

return

end function rho_next



! ----------------------------------------------------------------------------------%%
!
! RHO_ONE
!
!
! ----------------------------------------------------------------------------------%%

function rho_one(h_in)

use globals
use initialize
implicit none

! declare errthing
integer :: i,j
real(4) :: h_in, rho_one


		rho_one = rho_fluid*(1.0 - alpha*((h_in-273.0)))


return

end function rho_one


! ----------------------------------------------------------------------------------%%
!
! VISC_NEXT
!
! SUMMARY : solves for density using linear thermally expansive equation of state
!
! INPUTS : h_in(xn,yn) : 2D temperature array of current timestep
!
! RETURNS : rho_next(xn,yn) : 2D density array of current timestep
!
! ----------------------------------------------------------------------------------%%

function visc_next(h_in)

use globals
use initialize
implicit none

! declare errthing
integer :: i,j
real(4) :: h_in(xn,yn), visc_next(xn,yn)


do i=1,xn
	do j = 1,yn
		!visc_next(i,j) = (2.44e-5)*10.0**(247.8/(h_in(i,j)-140.0))
		visc_next(i,j) = .001!.001!.0002 !.00350*exp(-(h_in(i,j)-273.16)/35.0) + .0002
		! from nist.gov
	end do
end do

return

end function visc_next




! ----------------------------------------------------------------------------------%%
!
! VELOCITIES
!
! SUMMARY : computes the darcy velocity (specific discharge) from the streamfunction
!           using finite difference partial derivatives
!
! INPUTS : psi(xn,yn) : 2D streamfunction array of current timestep
!
! RETURNS : velocities(xn,2*yn) : both u and v velocities in one matrix
!
! ----------------------------------------------------------------------------------%%


function velocities(psi)

use globals
use initialize
implicit none

interface

	function partial_edge(array,rows,cols,d1,d2,dim)
		use globals
		use initialize
		implicit none
		integer :: rows, cols, dim, i, j, ii, jj
		real(4) :: array(rows,cols), d1, d2, d
		real(4) :: partial_edge(rows,cols)
	end function partial_edge

	function partial_edge_p(array,rows,cols,d1,d2,dim)
		use globals
		use initialize
		implicit none
		integer :: rows, cols, dim, i, j, ii, jj
		real(4) :: array(rows,cols), d1, d2, d
		real(4) :: partial_edge_p(rows,cols)
	end function partial_edge_p

end interface

! declare errthing
integer :: i,ii
real(4) :: velocities(xn,2*yn), psi(xn,yn)
real(4) :: u0(xn,yn), v0(xn,yn)

u0 = partial_edge_p(psi,xn,yn,dx,dy,2)
v0 = -partial_edge_p(psi,xn,yn,dx,dy,1)

velocities(1:xn,1:yn) = u0
velocities(1:xn,yn+1:2*yn) = v0

return
end function velocities






! ----------------------------------------------------------------------------------%%
!
! PARTIAL
!
! SUMMARY : versatile function for solving for second-order accurate partial
!           derivatives of 1D or 2D arrays with respect to specified dimension
!
! INPUTS : array(rows,cols) : array to be partially differentiated
!          rows : number of rows
!          cols : number of columns
!          d1 : grid spacing in first dimension
!          d2 : grid spacing in second dimension
!          dim : dimension you differentiate w.r.t.
!
! ----------------------------------------------------------------------------------%%

function partial(array,rows,cols,d1,d2,dim)

use globals
use initialize
implicit none

! declare errthing
integer :: rows, cols, dim, i, j, ii, jj
real(4) :: array(rows,cols), d1, d2, d
real(4) :: partial(rows,cols)

! write(*,*) "dim"
! write(*,*) rows
! write(*,*) cols

! figure out which direction derivative goes (dx or dy)

partial = 0.0

if (dim .eq. 1) then
	ii = 1
	jj = 0
	d = d1

	! compute edges beforehand
	partial(1,:) = ( -3.0*array(1,:) + 4.0*array(2,:) -array(3,:)) / (2.0*d)
	partial(rows,:) = ( 3.0*array(rows,:) - 4.0*array(rows-1,:) + array(rows-2,:) ) / (2.0*d)

	do i = 2,rows-1
		do j = 1,cols
			partial(i,j) = (array(i+1,j) - array(i-1,j))/(2.0*d)
		end do
	end do

	do i = 2,rows-1
		do j = 1,cols
			if ((maskP(i,j) .eq. 3.0) .or. (maskP(i,j) .eq. 3.5) .or. (maskP(i,j) .eq. 3.1) .or. (maskP(i,j) .eq. 3.05) .or. (maskP(i,j) .eq. 3.01) .or. (mask(i,j) .eq. 3.05)) then
				partial(i,j) = (array(i,j) - array(i-1,j))/d
				!partial(i,j) = (3.0*array(i,j) - 4.0*array(i-1,j) + array(i-2,j))/(2.0*d)
			end if
			if ((maskP(i,j) .eq. 6.0) .or. (maskP(i,j) .eq. 6.5) .or. (maskP(i,j) .eq. 6.1) .or. (maskP(i,j) .eq. 6.05) .or. (maskP(i,j) .eq. 6.01) .or. (mask(i,j) .eq. 6.05)) then
				partial(i,j) = (array(i+1,j) - array(i,j))/d
				!partial(i,j) = (-3.0*array(i,j) + 4.0*array(i+1,j) - array(i+2,j))/(2.0*d)
			end if
		end do
	end do

end if


if (dim .eq. 2) then
	ii = 0
	jj = 1
	d = d2

	! compute edges beforehand
	partial(:,1) = ( -3.0*array(:,1) + 4.0*array(:,2) -array(:,3)) / (2.0*d)
	partial(:,cols) = ( 3.0*array(:,cols) - 4.0*array(:,cols-1) + array(:,cols-2) ) / (2.0*d)
end if

! ! use central difference method ignoring edges (already done)
! do i = 2-jj,rows-1+jj
!     do j = 2-ii,cols-1+ii
!     	partial(i,j) = (array(i+ii,j+jj)-array(i-ii,j-jj))/(2.0*d)
! 	end do
! end do


! if (dim .eq. 1) then
! 	do i=2,rows-1
! 		do j=2,cols-1
!
! 			if (maskP(i,j) .eq. 5.0) then
! 				partial(i,j) = (array(i,j) - array(i-1,j))/d
! 				partial(i+1,j) = (array(i+1,j) - array(i,j))/d
! 			end if
!
! 			if (maskP(i,j) .eq. 10.0) then
! 				partial(i,j) = (array(i+1,j) - array(i,j))/d
! 				partial(i-1,j) = (array(i,j) - array(i-1,j))/d
! 			end if
!
!
! 		end do
! 	end do
!
!
! end if

! ! if jj=1, dim=2, cols=yn, rows=xn

! ! compute edges beforehand
! partial(:,1) = ( -3.0*array(:,1) + 4.0*array(:,2) -array(:,3)) / (2.0*d)
! partial(:,yn) = ( 3.0*array(:,yn) - 4.0*array(:,yn-1) + array(:,yn-2) ) / (2.0*d)


! do i = 1,xn
!     do j = 2,yn-1
!     	partial(i,j) = (array(i,j+1)-array(i,j-1))/(2.0*d)
! 	end do
! end do






! ! if ii=1, dim=1, cols=yn, rows=xn

! ! compute edges beforehand
! partial(1,:) = ( -3.0*array(1,:) + 4.0*array(2,:) -array(3,:)) / (2.0*d)
! partial(xn,:) = ( 3.0*array(xn,:) - 4.0*array(xn-1,:) + array(xn-2,:) ) / (2.0*d)


! do i = 2,xn-1
!     do j = 1,yn
!     	partial(i,j) = (array(i+1,j)-array(i-1,j))/(2.0*d)
! 	end do
! end do

return
end function partial







! ----------------------------------------------------------------------------------%%
!
! PARTIAL_EDGE
!
! SUMMARY : versatile function for solving for second-order accurate partial
!           derivatives of 1D or 2D arrays with respect to specified dimension
!
! INPUTS : array(rows,cols) : array to be partially differentiated
!          rows : number of rows
!          cols : number of columns
!          d1 : grid spacing in first dimension
!          d2 : grid spacing in second dimension
!          dim : dimension you differentiate w.r.t.
!
! ----------------------------------------------------------------------------------%%

function partial_edge(array,rows,cols,d1,d2,dim)

use globals
use initialize
implicit none

! declare errthing
integer :: rows, cols, dim, i, j, ii, jj
real(4) :: array(rows,cols), d1, d2, d
real(4) :: partial_edge(rows,cols)

! write(*,*) "dim"
! write(*,*) rows
! write(*,*) cols


! figure out which direction derivative goes (dx or dy)
if (dim .eq. 1) then
	ii = 1
	jj = 0
	d = d1
end if

if (dim .eq. 2) then
	ii = 0
	jj = 1
	d = d2
end if
partial_edge = 0.0
! use central difference method ignoring edges (already done)
do i = 2-jj,rows-1+jj
    do j = 2-ii,cols-1+ii
		if ((mask(i,j) .ne. 0.0)) then
    		partial_edge(i,j) = (array(i+ii,j+jj)-array(i-ii,j-jj))/(2.0*d)
		end if
	end do
end do

! do i = 2,rows-1
!     do j = 2,cols-1
! 		if ((mask(i,j) .ne. 0.0) .or. (mask(i+jj,j) .ne. 0.0) .or. (mask(i-jj,j-ii) .ne. 0.0)) then
!     		partial_edge(i,j) = (array(i+ii,j+jj)-array(i-ii,j-jj))/(2.0*d)
! 		end if
! 	end do
! end do

do i = 2,rows-1
    do j = 2,cols-1
		if ((mask(i,j) .ne. 0.0) .or. (mask(i-1,j) .ne. 0.0) .or. (mask(i+1,j) .ne. 0.0) .or. (mask(i,j-1) .ne. 0.0) .or. (mask(i,j+1) .ne. 0.0) ) then
    		partial_edge(i,j) = (array(i+ii,j+jj)-array(i-ii,j-jj))/(2.0*d)
		end if
	end do
end do


! figure out which direction derivative goes (dx or dy)
if (dim .eq. 1) then
	! compute edges beforehand
	partial_edge(1,:) = ( -3.0*array(1,:) + 4.0*array(2,:) -array(3,:)) / (2.0*d)
	partial_edge(rows,:) = ( 3.0*array(rows,:) - 4.0*array(rows-1,:) + array(rows-2,:) ) / (2.0*d)

	! inner edges
	do j=2,cols-1
		do i=2,rows-1
			if ((mask(i-1,j) .eq. 12.5)) then
				partial_edge(i,j) =( 3.0*array(i,j) - 4.0*array(i-1,j) + array(i-2,j) ) / (2.0*d)!(array(i+1,j) - array(i,j))/(d)
			end if
			if ((mask(i-1,j) .eq. 5.0) .and. (mask(i-1,j) .eq. 5.0)) then
				partial_edge(i,j) = ( 3.0*array(i,j) - 4.0*array(i-1,j) + array(i-2,j) ) / (2.0*d)!(array(i+1,j) - array(i,j))/(d)
			end if

			if ((mask(i+1,j) .eq. 17.5)) then
				partial_edge(i,j) = ( -3.0*array(i,j) + 4.0*array(i+1,j) - array(i+2,j) ) / (2.0*d)
			end if
			if ((mask(i+1,j) .eq. 10.0) .and. (mask(i+1,j) .eq. 10.0)) then
				partial_edge(i,j) = ( -3.0*array(i,j) + 4.0*array(i+1,j) - array(i+2,j) ) / (2.0*d)!(array(i,j) - array(i-1,j))/(d)
			end if


			if (mask(i,j-1) .eq. 12.5) then
				partial_edge(i,j) = ( 3.0*array(i,j) - 4.0*array(i-1,j) + array(i-2,j) ) / (2.0*d)!(array(i+1,j) - array(i,j))/(d)
			end if

			if (mask(i,j-1) .eq. 17.5) then
				partial_edge(i,j) =( -3.0*array(i,j) + 4.0*array(i+1,j) - array(i+2,j) ) / (2.0*d)!(array(i,j) - array(i-1,j))/(d)
			end if


		end do
	end do


	do j=1,cols
		do i=2,rows-1

			if (mask(i-1,j-1) .eq. 12.5) then
				partial_edge(i,j) = partial_edge(i-1,j)
			end if

			if (mask(i+1,j-1) .eq. 17.5) then
				partial_edge(i,j) = partial_edge(i+1,j)
			end if

		end do
	end do


end if




if (dim .eq. 2) then
	! compute edges beforehand
	partial_edge(:,1) = ( -3.0*array(:,1) + 4.0*array(:,2) -array(:,3)) / (2.0*d)
	partial_edge(:,cols) = ( 3.0*array(:,cols) - 4.0*array(:,cols-1) + array(:,cols-2) ) / (2.0*d)


	! inner edges
	do j=2,cols
		do i=1,rows
			if ((mask(i,j-1).eq.25.0) .or. (mask(i,j-1).eq.12.5)  .or. (mask(i,j-1).eq.17.5)) then
				partial_edge(i,j) = ( 3.0*array(i,j) - 4.0*array(i,j-1) + array(i,j-2) ) / (2.0*d)
			end if
		end do
	end do

	do j=2,cols-1
		do i=2,rows-1
			if ((mask(i,j-1).eq.50.0) .and. (mask(i,j-1).eq.50.0) .and. (mask(i,j-1).eq.50.0)) then
				partial_edge(i,j) =  ( 3.0*array(i,j) - 4.0*array(i,j-1) + array(i,j-2) ) / (2.0*d)!
			end if
		end do
	end do



	! JANUARY DELETION...
	do j=2,cols
		do i=2,rows-1

			if (mask(i-1,j) .eq. 12.5) then
				partial_edge(i,j) = ( 3.0*array(i,j) - 4.0*array(i,j-1) + array(i,j-2) ) / (2.0*d)
			end if

			if (mask(i+1,j) .eq. 17.5) then
				partial_edge(i,j) = ( 3.0*array(i,j) - 4.0*array(i,j-1) + array(i,j-2) ) / (2.0*d)
			end if

		end do
	end do

	do j=2,cols
		do i=2,rows-1

			if ((mask(i-1,j-1) .eq. 12.5)) then
				partial_edge(i,j) = partial_edge(i,j-1)
			end if

			if ((mask(i+1,j-1) .eq. 17.5)) then
				partial_edge(i,j) = partial_edge(i,j-1)
			end if


		end do
	end do



end if





return
end function partial_edge










! ----------------------------------------------------------------------------------%%
!
! partial_edge_p
!
! SUMMARY : versatile function for solving for second-order accurate partial
!           derivatives of 1D or 2D arrays with respect to specified dimension
!
! INPUTS : array(rows,cols) : array to be partially differentiated
!          rows : number of rows
!          cols : number of columns
!          d1 : grid spacing in first dimension
!          d2 : grid spacing in second dimension
!          dim : dimension you differentiate w.r.t.
!
! ----------------------------------------------------------------------------------%%


function partial_edge_p(array,rows,cols,d1,d2,dim)

use globals
use initialize
implicit none

! declare errthing
integer :: rows, cols, dim, i, j, ii, jj
real(4) :: array(rows,cols), d1, d2, d
real(4) :: partial_edge_p(rows,cols)

! write(*,*) "dim"
! write(*,*) rows
! write(*,*) cols

partial_edge_p = 0.0

! figure out which direction derivative goes (dx or dy)
if (dim .eq. 1) then
	ii = 1
	jj = 0
	d = d1

	do i=2,rows-1
		do j=1,cols
			if (maskP(i,j) .ne. 0.0) then
				partial_edge_p(i,j) = (array(i+1,j) - array(i-1,j)) / (2.0*d)
			end if
		end do
	end do

	do i=2,rows-1
		do j=1,cols
			if ((maskP(i,j) .eq. 5.0) .or. (maskP(i,j) .eq. 2.5)) then
				partial_edge_p(i,j) = (array(i,j) - array(i-1,j)) / (1.0*d)
			end if
		end do
	end do

	do i=1,rows-1
		do j=1,cols
! 			if (maskP(i,j) .eq. 6.0) then
! 				partial_edge_p(i,j) = 0.0
! 			end if
			if ((maskP(i,j) .eq. 6.0) .or. (maskP(i,j) .eq. 6.5) .or. (maskP(i,j) .eq. 6.1) .or. (maskP(i,j) .eq. 6.05) .or. (maskP(i,j) .eq. 6.01)) then
				partial_edge_p(i,j) = (array(i+1,j) - array(i,j)) / d
			end if
			if ((maskP(i,j) .eq. 3.0) .or. (maskP(i,j) .eq. 3.5) .or. (maskP(i,j) .eq. 3.1) .or. (maskP(i,j) .eq. 3.05) .or. (maskP(i,j) .eq. 3.01)) then
				partial_edge_p(i,j) = (array(i,j) - array(i-1,j)) / d
			end if
		end do
	end do
!
	! fracture vertical velocity
	!partial_edge_p(xn-1,:) = (array(xn,:) - array(xn-1,:)) / param_f_dx


end if

if (dim .eq. 2) then
	ii = 0
	jj = 1
	d = d2

	do i=1,rows
		do j=2,cols-1
			if ((maskP(i,j) .ne. 0.0)) then
				partial_edge_p(i,j) = (array(i,j+1) - array(i,j-1)) / (2.0*d)
			end if
		end do
	end do

	do i=1,rows
		do j=2,cols-1
			if ((maskP(i,j) .eq. 50.0) .or. (maskP(i,j) .eq. 3.5) .or. (maskP(i,j) .eq. 6.5)) then
				partial_edge_p(i,j) = (array(i,j) - array(i,j-1)) / (1.0*d)
			end if
		end do
	end do



end if




return
end function partial_edge_p




! ----------------------------------------------------------------------------------%%
!
! WRITE_VEC
!
! SUMMARY: Write linspace-style vector to file
!
! INPUTS: n : dimension
!         vector : vector with data
!         filename : file name
!
! RETURNS: write_vec
!
! ----------------------------------------------------------------------------------%%

function write_vec ( n, vector, filename )
use globals
  implicit none
  integer :: n, j, output_status, unit0
  character ( len = * ) filename
  real(4)  :: vector(n), write_vec



  unit0 = get_unit ()
  open ( unit = unit0, file = filename, status = 'replace', iostat = output_status )
  if ( output_status /= 0 ) then
    write ( *, '(a,i8)' ) 'COULD NOT OPEN OUTPUT FILE "' // &
      trim ( filename ) // '" USING UNIT ', unit0
    unit0 = -1
    stop
  end if


  if ( 0 < n ) then
    do j = 1, n
      write ( unit0, '(2x,g24.16)' ) vector(j)
    end do

  end if


  close ( unit = unit0 )
  write_vec = 1.0
  return
end function write_vec




! ----------------------------------------------------------------------------------%%
!
! WRITE_MATRIX
!
! SUMMARY: Write 2d array to file
!
! INPUTS: m,n : 2d dimensinons
!         table : 2d array with data
!         filename : file name
!
! RETURNS: write_matrix
!
! ----------------------------------------------------------------------------------%%

function write_matrix ( m, n, table, filename )
use globals
  implicit none
  integer :: m, n, j, output_status, unit0, reclen
  character ( len = * ) filename
  character ( len = 30 ) string
  real(4)  :: table(m,n) , write_matrix



  INQUIRE(iolength=reclen)table
  unit0 = get_unit ()
  open ( unit = unit0, file = filename, &
    status = 'replace', iostat = output_status, buffered='YES', buffercount=500)

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) 'Could not open the output file "' // &
      trim ( filename ) // '" on unit ', unit0
    unit0 = -1
    stop
  end if



!
 	write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 14, '.', 6, ')'
! 	!write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 14, '.', 6, ')'
!
!     do j = 1, n
!       write ( unit0, string) table(1:m,j)
!     end do

    do j = 1, n
      write ( unit0, 400) table(1:m,j)
    end do
400 FORMAT(<m>g14.6)


  close ( unit = unit0 )
  write_matrix = 2.0
  return
end function write_matrix
