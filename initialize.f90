module initialize

use globals
implicit none

save

integer :: g, gg, f, ff, long, longP, maskBit
real(4) :: x(xn), y(yn), t(tn)
real(4) :: rho0(xn,yn)
real(4) :: bcx0(xn,2), bcy0(2,yn), bcxPsi(xn,2), bcyPsi(2,yn), ic0(xn,yn)
real(4) :: kMat(xn,yn), lambdaMat(xn,yn), lambdaMatW(xn,yn), porosity(xn,yn), reactive(xn,yn)
real(4) :: sedx
real(4) :: sed(xn), sed1(xn), sed2(xn), sed3(xn)

! EXTRA BULLSHIT
real(4) :: permeable(xn,4), permeability(xn,yn)
real(4) :: mask(xn,yn), maskLong((xn-2)*(yn-2)), maskLongT((xn-2)*(yn-2)), h0Long((xn-2)*(yn-2)), h0LongT((xn-2)*(yn-2)), h0T(yn,xn)
real(4) :: maskPLongFull(xn*yn)
real(4) :: maskLongU((xn-2)*(yn-0)), maskLongTV((xn-0)*(yn-2))
real(4) :: maskP(xn,yn), maskPLong((xn-2)*((yn/2)-2)), maskPLongT((xn-2)*((yn/2)-2))
real(4) :: outerBand((xn-2)*((yn/2)-2),2*((yn/2)-2) + 1), bigBand((xn-2)*((yn/2)-2),4*((yn/2)-2) + 3)
real(4) :: stretch(xn,yn), stretchLong((xn-2)*(yn-2)), stretchT(xn,yn), stretchLongT((xn-2)*(yn-2))
real(4) :: inter, slope, inter2, slope2, boundary, buffer, edge, edge2




! coordinates for optimization
real(4) :: fives(2,yn), tens(2,yn), twentyfives(2,xn), fifties(2,xn)

! shell file parameters
character(len=300) :: path, path2, path_final, crashstring, restartstring, iso_path
character(len=300) :: param_o_string, param_w_string, param_w_rhs_string, param_h_string, param_o_rhs_string, param_tsw_string
character(len=300) :: param_dic_string, param_scope_string, param_trace_string, param_ch_string, param_f_dx_string, param_f_k_string
character(len=300) :: param_paq_string, param_ch_rhs_string, param_f_freq_string, param_f_por_string, param_h_s_string
integer :: in, crashstep, restart, param_trace
real(4):: param_o, param_w, param_w_rhs, param_h, param_o_rhs, param_tsw, param_dic, param_scope, param_ch
real(4) :: param_paq, param_ch_rhs, param_f_dx, param_f_k, param_f_freq, param_f_por, param_h_s


! TRANSPOSED
real(4) :: hTrans(yn,xn), psiTrans(yn,xn), permeabilityTrans(yn,xn), phiTrans(yn,xn)


! D2 STUFF
integer :: xn2, yn2, long2, m2
real(4) :: perm2 = 1e-12


real(4) :: frac6(yn,2), frac6_last(yn,2), temp6(yn,2), temp6_last(yn,2), temp6_mid(yn,2)
integer :: f_index1 = xn-5, iter = 0, spinup = 20!50000


real(4) :: temp6_a(yn), temp6_b(yn), temp6_c(yn), temp6_rhs(yn)


contains

! ----------------------------------------------------------------------------------%%
!
! SUBROUTINE TO INITIALIZE, JUST CALL IT
!
! ----------------------------------------------------------------------------------%%

subroutine init ()
use globals
integer :: m,n
!integer :: i,ii


in = iargc()
call getarg(1,restartstring)
call getarg(2,path)
call getarg(3,crashstring)
call getarg(4,param_o_string)
call getarg(5,param_w_string)
call getarg(6,param_w_rhs_string)
call getarg(7,param_h_string)
call getarg(8,param_o_rhs_string)
call getarg(9,param_paq_string)
call getarg(10,param_f_dx_string)
call getarg(11,param_f_por_string)
call getarg(12,param_h_s_string)


! parameter 1
read (restartstring, *) restart

! parameter 2 : path
path_final = path

! parameter 3
read (crashstring, *) crashstep

! parameter 4
read (param_o_string, *) param_o
param_o = param_o

! parameter 5
read (param_w_string, *) param_w

! parameter 6
read (param_w_rhs_string, *) param_w_rhs

! parameter 7
read (param_h_string, *) param_h
param_h = param_h

! parameter 8
read (param_o_rhs_string, *) param_o_rhs
param_o_rhs = param_o_rhs

! parameter 9
read (param_paq_string, *) param_paq

! parameter 10
read (param_f_dx_string, *) param_f_dx
param_f_dx = 10.0**(param_f_dx)

! parameter 11
read (param_f_por_string, *) param_f_por

! parameter 12
read (param_h_s_string, *) param_h_s



!permf = rho_fluid*grav*4.0*param_f_dx*param_f_dx/(12.0*viscosity)
permf = param_f_dx*param_f_dx/3.0

! SET UP THINGS THAT CAN'T BE DONE IN THE MODULE FOR WHATEVER REASON
dx = ( x_max - x_min ) / real ( xn - 1, kind = 4 )
x = linspace ( xn, x_min,x_max )
dy = ( y_max - y_min ) / real ( yn - 1, kind = 4 )
y = linspace ( yn, y_min, y_max )
dt = ( t_max - t_min ) / real ( tn - 1, kind = 4 )
t = linspace ( tn, t_min, t_max)


param_tsw = 275.0

! SCALE TO RESOLUTION


! BOUNDARY CONDITIONS
ic0(:,:) = 273.16 ! IC
bcx0(:,1) = 273.16 ! bottom
bcx0(:,2) = 273.16 ! top
bcy0(1,:) = 273.16 ! left
bcy0(2,:) = 273.16 ! right
bcyPsi(1,1:yn) = 0.0 ! left
bcyPsi(2,1:yn) = 0.0 ! right
bcxPsi(1:xn,1) = 0.0 ! bottom
bcxPsi(1:xn,2) = 0.0 ! top



! PERMEABILITY
lambdaMat = 1.8
permeability = 1e-17

! HEAT TRANSFER PARAMETERS
kMat = 2.0/(1000.0*4186.0)
ki=2.0/(1000.0*4186.0)


!reactive = 1.0









!-PERMEABILITY SET UP

slope = param_w
slope2 = x_max-param_w_rhs
buffer = 1000.0
edge = param_w
edge2 = x_max-param_w_rhs

! sed = -150.0
! sed1 = -300.0


sed = 0.0


! sed1 = -250.0 + (2.0*dy)
sed1 = (-1.0*param_h_s) + (2.0*dy)

!# param_o stuff
param_o = param_o - param_h + sed1(1)

!sed1 = sum(sed1)/xn

sed1((param_w/dx)+1:) = sed1(:xn-(param_w/dx))
sed1(1:param_w/dx) = sed1((param_w/dx)+1)

sed((param_w/dx)+1:) = sed(:xn-(param_w/dx))
sed(1:param_w/dx) = sed((param_w/dx)+1)

sed1(f_index1-2:f_index1+2) = sum(sed1(f_index1-3:f_index1+2))/6.0
sed(f_index1-2:f_index1+2) = sum(sed(f_index1-3:f_index1+2))/6.0

sed2 = sed-sed1!-sed1! - sed1


sed2((param_w/dx)+1:) = sed2(:xn-(param_w/dx))
sed2(1:param_w/dx) = sed2((param_w/dx)+1)

!sed2 = (-sum(sed2-sed1)/xn)*(sed2-sed1)

!sed = sum(sed)/xn
!sed = -100.0

! if (param_o_rhs .gt. param_o) then
! 	sed = sed-(param_o_rhs)
! end if

if (param_o .ge. param_o_rhs) then
	sed = sed-(param_o)
end if


! if (param_o_rhs .gt. param_o) then
! 	sed1 = sed1-(param_o_rhs)
! end if

if (param_o .ge. param_o_rhs) then
	sed1 = sed1-(param_o)
end if

do gg=1,yn
do g=1,xn
	if ((y(gg) .ge. sed1(g)) .and. (x(g) .gt. edge) .and. (x(g) .lt. edge2)) then
		 !lambdaMat(g,gg) = 1.2/((sed2(g))/(maxval(sed2)))
		 lambdaMat(g,gg) = 1.2!*(sum(sed2)/(xn*sed2(g)))

		 !lambdaMat(g,gg) = (sed2(g)/minval(sed2))*(sed2(g)/minval(sed2))*100.0*1.2*(sum(sed1)/xn)/(sed2(g)*(sed2(g)-sed1(g)))
		 !lambdaMat(g,gg) = 1.2*sed2(g)*xn/sum(sed2)

		!lambdaMat(g,gg) = 1.2/((sed(g)-sed1(g))/(minval(sed-sed1)))

	end if
end do
end do


! with sediment cap
sed3 = sed1 - (param_h)



	! the mask
	mask = 1.0
	do gg=1,yn
	do g =2,xn
			if ((x(g) .ge. edge) .and. (x(g) .le. edge2+5000.0) .and. (y(gg) .ge. sed(g))) then
				mask(g,gg) = 0.0
			end if
	end do
	end do

	do g =1,xn
			! left outcrop top
			if ((x(g) .lt. edge)) then
				if (param_o_rhs .gt. param_o) then
					mask(g,yn-(abs(param_o-param_o_rhs)/dy)-1) = 25.0
					mask(g,yn-(abs(param_o-param_o_rhs)/dy):yn) = 0.0
				end if
				if (param_o_rhs .lt. param_o) then
					mask(g,yn-1) = 25.0
					mask(g,yn:yn) = 0.0
				end if
				if (param_o_rhs .eq. param_o) then
					mask(g,yn-1) = 25.0
					mask(g,yn-0:yn) = 0.0
				end if
			end if

! 			! right outcrop top
! 			if ((x(g) .gt. edge2)) then
! 				if (param_o_rhs .gt. param_o) then
! 					mask(g,yn-1) = 25.0
! 					mask(g,yn-0:yn) = 0.0
! 				end if
! 				if (param_o_rhs .lt. param_o) then
! 					mask(g,yn-(abs(param_o-param_o_rhs)/dy)-1) = 25.0
! 					mask(g,yn-(abs(param_o-param_o_rhs)/dy):yn) = 0.0
! 				end if
! 				if (param_o_rhs .eq. param_o) then
! 					mask(g,yn-1) = 25.0
! 					mask(g,yn-0:yn) = 0.0
! 				end if
! 			end if
	end do

	! inner vertical edges
	do gg=1,yn
	do g =2,xn-1
		if ((mask(g,gg) - mask(g-1,gg) .eq. -1.0) .and. (x(g) .le. edge+dx)) then
			mask(g-1,gg) = 5.0
		end if
		if ((mask(g,gg) - mask(g-1,gg) .eq. 1.0) .and. (x(g) .ge. edge2-dx)) then
			mask(g,gg) = 10.0
		end if
	end do
	end do

	do gg=2,yn-1
	do g =1,xn
		if ((mask(g,gg) .eq. 1.0) .and. (mask(g,gg+1) .eq. 0.0)) then
			mask(g,gg) = 50.0
		end if
	end do
	end do

	do gg=2,yn-1
	do g =2,xn-1
		! left upper corner
		if ((mask(g,gg) .eq. 25.0) .and. (mask(g+1,gg-1) .eq. 5.0)) then
			mask(g+1,gg) = 12.5
		end if

		! right upper corner
		if ((mask(g,gg) .eq. 25.0) .and. (mask(g-1,gg-1) .eq. 10.0)) then
			mask(g-1,gg) = 17.5
		end if

	end do
	end do


	do gg=2,yn-1
	do g =2,xn-1
		! left bottom corner
		if ((mask(g,gg) .eq. 5.0) .and. (mask(g+1,gg-1) .eq. 50.0)) then
			mask(g,gg-1) = 2.5
! 			mask(g+1,gg) = 2.0
!
! 			mask(g,gg) = 1.0
! 			mask(g+1,gg-1) = 1.0
		end if

		! right bottom corner
		if ((mask(g,gg) .eq. 10.0) .and. (mask(g-1,gg-1) .eq. 50.0)) then
			mask(g,gg-1) = 7.5
! 			mask(g-1,gg) = 7.0
!
! 			mask(g,gg) = 1.0
! 			mask(g-1,gg-1) = 1.0
		end if
	end do
	end do






	maskP = mask

	do gg=1,yn-1
	do g =1,xn
		if (y(gg) .le. sed3(g)) then
			maskP(g,gg) = 0.0
		end if

	end do
	end do

	do gg=2,yn-3
	do g =1,xn

! 		if ((y(gg+1) .le. minval(sed3)) .and. (y(gg+2) .gt. minval(sed3))) then
! 			maskP(:,gg) = 100.0
! 			code = gg
! 		end if

if (gg .eq. yn/2 + 3) then
	maskP(:,gg) = 100.0
	code = gg
end if

	end do
	end do


	do gg=2,yn-3
	do g =1,xn

		if ((y(gg) .le. sed3(g)) .and. (gg .gt. code)) then
			maskP(:,gg) = 1.0
		end if

	end do
	end do









	fives = 1.0
	tens = 1.0
	twentyfives = 1.0
	fifties = 1.0
		do gg=1,yn
		do g=1,xn
			if ((maskP(g,gg) .eq. 5.0) .or. (maskP(g,gg) .eq. 12.5)) then
				fives(1,gg) = g
				fives(2,gg) = gg
			end if

			if ((maskP(g,gg) .eq. 10.0) .or. (maskP(g,gg) .eq. 17.5)) then
				tens(1,gg) = g
				tens(2,gg) = gg
			end if
		end do
		end do


		do g=1,xn
		do gg=1,yn
			if ((maskP(g,gg) .eq. 25.0) .or. (maskP(g,gg) .eq. 12.5) .or. (maskP(g,gg) .eq. 17.5)) then
				twentyfives(1,g) = g
				twentyfives(2,g) = gg
			end if

			if (maskP(g,gg) .eq. 50.0) then
				fifties(1,g) = g
				fifties(2,g) = gg
			end if
		end do
		end do

		write(*,*) "fives"
		write(*,*) fives


	! TOO SIMPLE
	permeability = param_paq
	do gg=1,yn-2
	do g=1,xn
		if ((any(maskP(g,:) .eq. 50.0)) .and. (y(gg) .ge. sed1(g)) .and. (y(gg) .le. sed(g)+dy) .and. (x(g) .le. x_max-param_w_rhs+5000.0)) then
			permeability(g:g,gg:gg) = 1e-17
		end if


		if ((y(gg) .lt. sed3(g))) then
			permeability(g,gg) = 1e-17
		end if

		if ((x(g) .le. param_w) .or. (x(g) .ge. x_max-param_w_rhs)) then
		end if

		if (g .gt. f_index1 - 1) then
		end if

	end do
	end do


	! 369 fracture goes HERE
	do gg=2,yn-3
	do g =1,xn

		if ((maskP(g,gg) .eq. 1.0) .and. (g .eq. f_index1) .and. (y(gg) .ge. sed1(g))) then
			maskP(g,gg) = 6.0
			maskP(g-1,gg) = 3.0
! 			maskP(g+1,gg) = 9.0

			mask(g,gg) = 6.0
			mask(g-1,gg) = 3.0
! 			mask(g+1,gg) = 9.0

! 			lambdaMat(g,gg) = 6.0
! 			lambdaMat(g-1,gg) = 6.0

			!permeability(g-1:g,gg) = param_paq
		end if

		!if ((maskP(g,gg) .eq. 100.0) .and. (g .eq. f_index1) .and. (y(gg) .ge. sed1(g))) then
		if ((maskP(g,gg-1) .eq. 1.0) .and. (maskP(g,gg) .eq. 6.0)) then

			maskP(g,code) = 6.01
			maskP(g,code+1:gg-1) = 6.05
			maskP(g,gg-1) = 6.1
			maskP(g-1,code) = 3.01
			maskP(g-1,code+1:gg-1) = 3.05
			maskP(g-1,gg-1) = 3.1
! 			maskP(g+1,gg) = 9.0

			mask(g,2:gg-1) = 6.05
			!permeability(g,2:gg-1) = 1e-18
			mask(g,gg-1) = 6.1
			mask(g-1,2:gg-1) = 3.05
			!permeability(g-1,2:gg-1) = 1e-18
			mask(g-1,gg-1) = 3.1
! 			mask(g+1,gg) = 9.0

! 			lambdaMat(g,gg) = 6.0
! 			lambdaMat(g-1,gg) = 6.0

			!permeability(g-1:g,gg) = param_paq
		end if

		if ((maskP(g,gg) .eq. 50.0) .and. (g .eq. f_index1) .and. (y(gg) .ge. sed1(g))) then
			maskP(g,gg) = 6.5
			maskP(g-1,gg) = 3.5
! 			maskP(g+1,gg) = 9.0

			mask(g,gg) = 6.5
			mask(g-1,gg) = 3.5
! 			mask(g+1,gg) = 9.0

!  			lambdaMat(g,gg:gg+1) = 6.0
!  			lambdaMat(g-1,gg:gg+1) = 6.0

			!permeability(g-1:g,gg:) = param_paq
		end if

	end do
	end do

! 	do gg=1,yn
! 	do g=2,xn-1
! 		if ((g .eq. 6) .and. (y(gg) .gt. sed3(g)+dy+dy)) then
! 			permeability(g,gg) = 1e-10
! 		end if
! 		if ((g .eq. 6) .and. (y(gg) .gt. sed3(g)+dy+dy) .and. (maskP(g,gg) .eq. 1.0)) then
! 			mask(g,gg) = 6.0
! 			maskP(g,gg) = 6.0
!
! 			mask(g+1,gg) = 9.0
! 			maskP(g+1,gg) = 9.0
! 		end if
! 	end do
! 	end do
!
! 	do gg=2,yn-1
! 	do g=2,xn-1
!
! 		if ((maskP(g,gg) .eq. 6.0) .and. (maskP(g,gg-1) .eq. 1.0)) then
! ! 			mask(g,gg) = 25.0
! ! 			maskP(g,gg) = 25.0
! !
! ! 			mask(g+1,gg) = 1.0
! ! 			maskP(g+1,gg) = 1.0
!
! 		end if
! 	end do
! 	end do


	! high lambda in deep basalt
	do gg=1,yn
	do g=1,xn
		if ((y(gg) .lt. sed1(g)) .and. (permeability(g,gg) .eq. 1e-18)) then
			lambdaMat(g,gg) = 1.8
		end if
	end do
	end do


	active_cells = 0
	do gg=1,yn
	do g=1,xn
		if (maskP(g,gg) .ne. 0.0) then
			active_cells = active_cells + 1
		end if
	end do
	end do
! 	ison = particle_sat*active_cells

	long = (xn-2)*(yn-2)
	longP = (xn-2)*((yn/2)-2)

! 	! SINGLE FRACTURE BOX
! 	mask = 1.0
! 	mask(21,yn/2+5:) = 6.0
! 	!lambdaMat(xn/2,yn/2+1:) = 0.01
! 	mask(21-1,yn/2+5:) = 3.0
! 	mask(21+1,yn/2+5:) = 9.0
! 	mask(:,yn-1) = 25.0
! 	mask(:,yn) = 0.0
! 	maskP = mask
! 	maskP(:,:yn/2+4) = 0.0
! 	maskP(:,yn/2+4) = 100.0


	maskLong = reshape(mask(2:xn-1,2:yn-1), (/long/))
	maskLongT = reshape(transpose(mask(2:xn-1,2:yn-1)), (/long/))

	maskLongU = reshape(mask(2:xn-1,1:yn), (/(xn-2)*(yn-0)/))
	maskLongTV = reshape(transpose(mask(1:xn,2:yn-1)), (/(xn-0)*(yn-2)/))

	maskPLong = reshape(maskP(2:xn-1,(yn/2)+2:yn-1), (/longP/))
	maskPLongT = reshape(transpose(maskP(2:xn-1,(yn/2)+2:yn-1)), (/longP/))



return

end subroutine init






function h_bc(h_in)

	use globals
	real(4) :: h_in(xn,yn), h_bc(xn,yn), rip_lith_y(xn)
	integer :: p, pp

	! do p = 1,xn/4
	! rip_lith_y(p) = .55 - 0.00035*p
	! end do
    !
	! do p = xn/4,xn
	! rip_lith_y(p) = .55 - 0.00035*(xn/4) - 0.000175*(p - xn/4)
	! end do

	do p = 1,xn
		!rip_lith_y(p) = 0.5 / (((0.75+0.035*dx*p/1000.0)**0.5))
		rip_lith_y(p) = 0.55
	end do

!((0.75+0.035*50*x/1000)

! 	rip_lith_y = (/0.602087135445, &
! & 0.60116684559, 0.600246555734, 0.599326265879, 0.598405976024, 0.597485686169, 0.596565396314, &
! & 0.595645106459, 0.594724816604, 0.593804526749, 0.592884236894, 0.591963947039, 0.591043657184, &
! & 0.590223660487, 0.589514072973, 0.588804485459, 0.588094897945, 0.587385310431, 0.586675722917, &
! & 0.585966135403, 0.585256547889, 0.584546960375, 0.583837372862, 0.583127785348, 0.582418197834, &
! & 0.58170861032, 0.580999022806, 0.580289435292, 0.579579847778, 0.578870260264, 0.57816067275, &
! & 0.577451085236, 0.576741497722, 0.576031910208, 0.575322322694, 0.57461273518, 0.573903147666, &
! & 0.573193560152, 0.572483972639, 0.571774385125, 0.571064797611, 0.570355210097, 0.569645622583, &
! & 0.568936035069, 0.568226447555, 0.567516860041, 0.566810372198, 0.566178703989, 0.565547035781, &
! & 0.564915367572, 0.564283699364, 0.563652031155, 0.563020362947, 0.562388694739, 0.56175702653, &
! & 0.561125358322, 0.560493690113, 0.559862021905, 0.559230353696, 0.558598685488, 0.557967017279, &
! & 0.557335349071, 0.556703680862, 0.556072012654, 0.555440344446, 0.554808676237, 0.554177008029, &
! & 0.55354533982, 0.552913671612, 0.552282003403, 0.551650335195, 0.551018666986, 0.550386998778, &
! & 0.54975533057, 0.549123662361, 0.548491994153, 0.547860325944, 0.547228657736, 0.546596989527, &
! & 0.545965321319, 0.54533365311, 0.544701984902, 0.544070316694, 0.543438648485, 0.542806980277, &
! & 0.542175312068, 0.54154364386, 0.540911975651, 0.540280307443, 0.539648639234, 0.539016971026, &
! & 0.538385302817, 0.537753634609, 0.537121966401, 0.536490298192, 0.535858629984, 0.535226961775, &
! & 0.534595293567, 0.533963625358, 0.53333195715, 0.532700288941, 0.532068620733, 0.531436952525, &
! & 0.530805284316, 0.530173616108, 0.529541947899, 0.528910279691, 0.528278611482, 0.527646943274, &
! & 0.527015275065, 0.526383606857, 0.525751938648, 0.52512027044, 0.524488602232, 0.523856934023, &
! & 0.523225265815, 0.522593597606, 0.521961929398, 0.521330261189, 0.520698592981, 0.520066924772, &
! & 0.519435256564, 0.518803588356, 0.518171920147, 0.517540251939, 0.51690858373, 0.516276915522, &
! & 0.515681666368, 0.515130941731, 0.514580217093, 0.514029492455, 0.513478767818, 0.51292804318, &
! & 0.512377318542, 0.511826593905, 0.511275869267, 0.510725144629, 0.510174419992, 0.509623695354, &
! & 0.509072970716, 0.508522246079, 0.507971521441, 0.507420796803, 0.506870072166, 0.506319347528, &
! & 0.50576862289, 0.505217898253, 0.504667173615, 0.504116448977, 0.503565724339, 0.503014999702, &
! & 0.502464275064, 0.501913550426, 0.501362825789, 0.500812101151, 0.500261376513, 0.499710651876, &
! & 0.499159927238, 0.498679752216, 0.498249474112, 0.497819196008, 0.497388917903, 0.496958639799, &
! & 0.496528361695, 0.496098083591, 0.495667805487, 0.495237527382, 0.494807249278, 0.494376971174, &
! & 0.49394669307, 0.493516414966, 0.493086136861, 0.492655858757, 0.492225580653, 0.491795302549, &
! & 0.491365024445, 0.490934746341, 0.490504468236, 0.490074190132, 0.489643912028, 0.489213633924, &
! & 0.48878335582, 0.488353077715, 0.487922799611, 0.487492521507, 0.487062243403, 0.486631965299, &
! & 0.486201687194, 0.48577140909, 0.485341130986, 0.484873859367, 0.48440077448, 0.483927689594, &
! & 0.483454604708, 0.482981519822, 0.482508434936, 0.48203535005, 0.481562265164, 0.481089180278, &
! & 0.480616095391, 0.480143010505, 0.479669925619, 0.479196840733, 0.478723755847, 0.478250670961, &
! & 0.477777586075, 0.477304501189, 0.476831416302, 0.476358331416, 0.47588524653, 0.475412161644, &
! & 0.474939076758, 0.474465991872, 0.473992906986, 0.4735198221, 0.473046737213, 0.472573652327, &
! & 0.472100567441, 0.471627482555, 0.471154397669, 0.470681312783, 0.470208227897, 0.469735143011, &
! & 0.469262058124, 0.468788973238, 0.468315888352, 0.467944418472, 0.46759004727, 0.467235676069, &
! & 0.466881304867, 0.466526933666, 0.466172562464, 0.465818191263, 0.465463820061, 0.46510944886, &
! & 0.464755077658, 0.464400706457, 0.464046335255, 0.463691964054, 0.463337592852, 0.462983221651, &
! & 0.462628850449, 0.462274479248, 0.461920108046, 0.461565736845, 0.461211365643, 0.460856994442, &
! & 0.46050262324, 0.460148252039, 0.459793880837, 0.459439509636, 0.459085138434, 0.458730767233, &
! & 0.458376396031, 0.45802202483, 0.457667653628, 0.457313282427, 0.456958911225, 0.456604540024, &
! & 0.456250168822, 0.455895797621, 0.455541426419, 0.455187055218, 0.454832684016, 0.454478312815, &
! & 0.454123941613, 0.453769570412, 0.45341519921, 0.453060828009, 0.452706456807, 0.452352085606, &
! & 0.451997714404, 0.451643343203, 0.451288972001, 0.4509346008, 0.450580229598, 0.450225858397, &
! & 0.449871487195, 0.449517115994, 0.449162744792, 0.448808373591, 0.448454002389, 0.448099631188, &
! & 0.447745259986, 0.447390888785, 0.447036517584, 0.446682146382, 0.446327775181, 0.445973403979, &
! & 0.445619032778, 0.445264661576, 0.444910290375, 0.444555919173, 0.444201547972, 0.44384717677, &
! & 0.443492805569, 0.443138434367, 0.442784063166, 0.442429691964, 0.442075320763, 0.441720949561, &
! & 0.44136657836, 0.441012207158, 0.440657835957, 0.440303464755, 0.439949093554, 0.439594722352, &
! & 0.439240351151, 0.438885979949, 0.438531608748, 0.438177237546, 0.437822866345, 0.437468495143, &
! & 0.437114123942, 0.43675975274, 0.436405381539, 0.436051010337, 0.435696639136, 0.435342267934, &
! & 0.434987896733, 0.434633525531, 0.43427915433, 0.433924783128, 0.433570411927, 0.433216040725, &
! & 0.432861669524, 0.432507298322, 0.432152927121, 0.431798555919, 0.431444184718, 0.431089813516, &
! & 0.430735442315, 0.430422684669, 0.430144535282, 0.429866385895, 0.429588236509, 0.429310087122, &
! & 0.429031937735, 0.428753788348, 0.428475638961, 0.428197489574, 0.427919340187, 0.427641190801, &
! & 0.427363041414, 0.427084892027, 0.42680674264, 0.426528593253, 0.426250443866, 0.42597229448, &
! & 0.425694145093, 0.425415995706, 0.425137846319, 0.424859696932, 0.424581547545, 0.424303398158, &
! & 0.424025248772, 0.423747099385, 0.423468949998, 0.423190800611, 0.422912651224, 0.422634501837, &
! & 0.422356352451, 0.422078203064, 0.421800053677, 0.42152190429, 0.421243754903, 0.420965605516, &
! & 0.420687456129, 0.420409306743, 0.420131157356, 0.419853007969, 0.419574858582, 0.419296709195, &
! & 0.419018559808, 0.418740410422, 0.418462261035, 0.418184111648, 0.417905962261, 0.417627812874, &
! & 0.417349663487, 0.4170715141, 0.416793364714, 0.416515215327, 0.41623706594, 0.415958916553, &
! & 0.415680767166, 0.415402617779, 0.415124468393, 0.414878071636, 0.414638695561, 0.414399319487, &
! & 0.414159943412, 0.413920567337, 0.413681191262, 0.413441815188, 0.413202439113, 0.412963063038, &
! & 0.412723686964, 0.412484310889, 0.412244934814, 0.41200555874, 0.411766182665, 0.41152680659, &
! & 0.411287430516, 0.411048054441, 0.410808678366, 0.410569302292, 0.410329926217, 0.410090550142, &
! & 0.409851174068, 0.409611797993, 0.409372421918, 0.409133045844, 0.408893669769, 0.408654293694, &
! & 0.40841491762, 0.408175541545, 0.40793616547, 0.407696789396, 0.407457413321, 0.407218037246, &
! & 0.406978661172, 0.406739285097, 0.406499909022, 0.406260532948, 0.406021156873, 0.405781780798, &
! & 0.405542404724, 0.405303028649, 0.405063652574, 0.4048242765, 0.404584900425, 0.40434552435, &
! & 0.404106148276, 0.403866772201, 0.403627396126, 0.403388020051, 0.403148643977, 0.402914577103, &
! & 0.402697703397, 0.402480829691, 0.402263955985, 0.402047082279, 0.401830208573, 0.401613334867, &
! & 0.401396461161, 0.401179587455, 0.400962713749, 0.400745840043, 0.400528966337, 0.400312092631, &
! & 0.400095218925, 0.399878345219, 0.399661471513, 0.399444597807, 0.399227724101, 0.399010850395, &
! & 0.398793976689, 0.398577102983, 0.398360229277, 0.398143355571, 0.397926481865, 0.397709608159, &
! & 0.397492734453, 0.397275860747, 0.397058987041, 0.396842113335, 0.396625239629, 0.396408365923, &
! & 0.396191492217, 0.395974618511, 0.395757744805, 0.395540871099, 0.395323997393, 0.395107123687, &
! & 0.394890249981, 0.394673376275, 0.394456502569, 0.394239628863, 0.394022755157, 0.393805881451, &
! & 0.393589007745, 0.393372134039, 0.393155260333, 0.392938386627, 0.392721512921, 0.3925131885, &
! & 0.392334317409, 0.392155446318, 0.391976575228, 0.391797704137, 0.391618833046, 0.391439961955, &
! & 0.391261090865, 0.391082219774, 0.390903348683, 0.390724477592, 0.390545606501, 0.390366735411, &
! & 0.39018786432, 0.390008993229, 0.389830122138, 0.389651251048, 0.389472379957, 0.389293508866, &
! & 0.389114637775, 0.388935766685, 0.388756895594, 0.388578024503, 0.388399153412, 0.388220282321, &
! & 0.388041411231, 0.38786254014, 0.387683669049, 0.387504797958, 0.387325926868, 0.387147055777, &
! & 0.386968184686, 0.386789313595, 0.386610442505, 0.386431571414, 0.386252700323, 0.386073829232, &
! & 0.385894958141, 0.385716087051, 0.38553721596, 0.385358344869, 0.385179473778, 0.385000602688, &
! & 0.384821731597, 0.384642860506, 0.384463989415, 0.384285118325, 0.384106247234, 0.383927376143, &
! & 0.383763195563, 0.383625514403, 0.383487833244, 0.383350152084, 0.383212470925, 0.383074789766, &
! & 0.382937108606, 0.382799427447, 0.382661746287, 0.382524065128, 0.382386383969, 0.382248702809, &
! & 0.38211102165, 0.38197334049, 0.381835659331, 0.381697978171, 0.381560297012, 0.381422615853, &
! & 0.381284934693, 0.381147253534, 0.381009572374, 0.380871891215, 0.380734210055, 0.380596528896, &
! & 0.380458847737, 0.380321166577, 0.38027111787, 0.38027111787, 0.38027111787, 0.38027111787, &
! & 0.38027111787, 0.38027111787, 0.38027111787, 0.38027111787, 0.38027111787, 0.38027111787, &
! & 0.38027111787, 0.38027111787, 0.38027111787, 0.38027111787, 0.38027111787, 0.38027111787, &
! & 0.38027111787, 0.38027111787, 0.38027111787, 0.38027111787/)

	h_bc = h_in


	! top of outcrops
	do pp=1,yn
	do p=1,xn
		if (mask(p,pp) .eq. 0.0) then
			 h_bc(p,pp) = param_tsw
		end if
	end do
	end do



	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!    VERTICAL OUTER   !!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	do pp=1,yn
		if (mask(xn,pp) .ne. 0.0) then
			h_bc(xn,pp) = (4.0/3.0)*h_in(xn-1,pp) - (1.0/3.0)*h_in(xn-2,pp) ! right
		end if
	end do

	do pp=1,yn
		if (mask(1,pp) .ne. 0.0) then
			h_bc(1,pp) = (4.0/3.0)*h_in(2,pp) - (1.0/3.0)*h_in(3,pp) ! left
		end if
	end do


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!   HORIZONTAL OUTER  !!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! bottom
	do p = 1,xn
		!h_bc(p,1) = h_in(p,2) + ( .48) * dy/2.0
		h_bc(p,1) = h_in(p,2) + ( rip_lith_y(p)) * dy/1.8
	end do

	! two lines recent...
	 	h_bc(xn,1) = h_bc(xn,2)
	 	h_bc(1,1) = h_bc(1,2)


	! top of outcrops
	do p=1,xn
 	do pp=2,yn-1
		if (mask(p,pp) .eq. 25.0) then
			 h_bc(p,pp+1) = param_tsw
		end if

		if (mask(p,pp) .eq. 12.5) then
			 h_bc(p,pp+1) = param_tsw
		end if
		if (mask(p,pp) .eq. 17.5) then
			 h_bc(p,pp+1) = param_tsw
		end if
	end do
	end do

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!   HORIZONTAL INNER  !!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! top of sediment
	do p=2,xn-1
 	do pp=2,yn-1
		if ((mask(p,pp) .eq. 50.0) .or. (mask(p,pp) .eq. 3.5) .or. (mask(p,pp) .eq. 6.5)) then
			 h_bc(p,pp+1) = param_tsw
		end if
	end do
	end do

	do p=1,xn
 	do pp=2,yn-1
		if ((mask(p,pp) .eq. 12.5) .or. (mask(p,pp) .eq. 17.5)) then
			 h_bc(p,pp+1) = param_tsw
		end if
	end do
	end do


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!    VERTICAL INNER   !!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	do p=2,xn-1
 	do pp=2,yn

		! left outcrop
		if ((mask(p,pp) .eq. 12.5)) then
			h_bc(p+1,pp) = param_tsw
		end if

		if ((mask(p,pp) .eq. 5.0)) then
			h_bc(p+1,pp) = param_tsw
		end if

		! right outcrop
		if ((mask(p,pp) .eq. 17.5)) then
			h_bc(p-1,pp) = param_tsw
		end if

		if ((mask(p,pp) .eq. 10.0)) then
			h_bc(p-1,pp) = param_tsw
		end if

	end do
	end do



! 	do p=2,xn-1
!  	do pp=2,yn
!
! 		! left outcrop
! 		if ((mask(p,pp) .ne. 5.0) .and. (maskP(p,pp) .eq. 5.0)) then
! 			h_bc(p+1,pp) = param_tsw + 5.0
! 		end if
!
! 		! right outcrop
! 		if ((mask(p,pp) .ne. 10.0) .and. (maskP(p,pp) .eq. 10.0)) then
! 			h_bc(p-1,pp) = param_tsw + 5.0
! 		end if
!
! 	end do
! 	end do

! 	!optimized
! 	h_bc(fives(:,1)+1,fives(:,2)) = param_tsw
! 	h_bc(tens(:,1)-1,tens(:,2)) = param_tsw
! 	h_bc(twelve(1,1)+1,twelve(1,2)) = param_tsw
! 	h_bc(seventeen(1,1)-1,seventeen(1,2)) = param_tsw

! 	!!!! ACTUAL CORNERS
! 	do i=2,xn-1
! 		do ii=1,yn
!
! 			! right upper corner
! 			if ((mask(i,ii) .eq. 17.5)) then
! 			h_bc(i-1,ii+1) = h_bc(i,ii+1)
! 			end if
!
! 			! left upper corner
! 			if ((mask(i,ii) .eq. 12.5)) then
! 			h_bc(i+1,ii+1) = h_bc(i,ii+1)
! 			end if
! 		end do
! 	end do
!

! 	! SINGLE FRACTURE BOX
! 	h_bc(:,1) = 285.0
! 	! top of sediment
! 	do p=1,xn
! 	do pp=1,yn-1
! 		if ((mask(p,pp) .eq. 25.0)) then
! 			 h_bc(p,pp+1) = 275.0
! 		end if
!
! 		if ((mask(p,pp) .eq. 1.0) .and. (maskP(p,pp) .eq. 0.0)) then
! 			 h_bc(p,pp+1) = 285.0
! 		end if
! 	end do
! 	end do

return
end function h_bc










function psi_bc(psi_in)

	use globals
	real(4) :: psi_in(xn,yn), psi_bc(xn,yn)
	integer :: p,pp

	psi_bc = psi_in



	do pp=1,yn
    do p=1,xn
		if ((maskP(p,pp) .eq. 0.0)) then
			psi_bc(p,pp) = 0.0
		end if
	end do
	end do


    do p=1,xn
    	do pp=2,yn-1

    		! top of outcrops
    		if ((maskP(p,pp) .eq. 25.0) .or. (maskP(p,pp) .eq. 12.5) .or. (maskP(p,pp) .eq. 17.5)) then
   				psi_bc(p,pp+1) = (4.0/3.0)*psi_in(p,pp) - (1.0/3.0)*psi_in(p,pp-1)
    		end if

    	end do
    end do



	do pp=1,yn
		! left
		if (maskP(1,pp) .ne. 0.0) then
			psi_bc(1,pp) = 0.0
		end if
		! right
		if (maskP(xn,pp) .ne. 0.0) then
			psi_bc(xn,pp) =0.0
		end if
	end do

	! 12/09

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!   HORIZONTAL OUTER  !!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	do p=1,xn
		do pp=2,yn-1
			if ((maskP(p,pp) .eq. 100.0) .or. (maskP(p,pp) .eq. 3.01) .or. (maskP(p,pp) .eq. 6.01)) then
				psi_bc(p,pp-1) = 0.0
			end if

		end do
	end do






	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!   HORIZONTAL INNER  !!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




	do p=2,xn-1
		do pp=2,yn-1


				! right outcrop (left boundary)
				if ((maskP(p,pp) .eq. 10.0)) then
					psi_bc(p-1,pp) = (4.0/3.0)*psi_in(p,pp) - (1.0/3.0)*psi_in(p+1,pp)
				end if

				if ((maskP(p,pp) .eq. 17.5)) then
					psi_bc(p-1,pp) = (4.0/3.0)*psi_in(p,pp) - (1.0/3.0)*psi_in(p+1,pp)
				end if

				! left outcrop (right boundary)
				if ((maskP(p,pp) .eq. 5.0)) then
					psi_bc(p+1,pp) = (4.0/3.0)*psi_in(p,pp) - (1.0/3.0)*psi_in(p-1,pp)
				end if

				if ((maskP(p,pp) .eq. 12.5)) then
					psi_bc(p+1,pp) = (4.0/3.0)*psi_in(p,pp) - (1.0/3.0)*psi_in(p-1,pp)
				end if

		end do
	end do

	do p=2,xn
		do pp=2,yn-1

			! top of sediment
			if ((maskP(p,pp) .eq. 50.0) .or. (maskP(p,pp) .eq. 3.5) .or. (maskP(p,pp) .eq. 6.5)) then
				psi_bc(p,pp+1) = (4.0/3.0)*psi_in(p,pp) - (1.0/3.0)*psi_in(p,pp-1)
			end if



		end do
	end do


! 	do p=2,xn
! 		do pp=2,yn-1
!
! 			! top of sediment
! 			if ((maskP(p,pp) .eq. 50.0) .and. (p .lt. 20)) then
! 				psi_bc(p,pp+1) = psi_bc(20,pp+1)
! 			end if
!
!
!
! 		end do
! 	end do

	! SUSTAINED

! 	do p=1,xn
! 		do pp=2,yn-1
!
! ! 			! top of sediment
! ! 			if ((maskP(p,pp) .eq. 50.0) .and. (p .eq. xn)) then
! ! 				psi_bc(p,pp+1) = 0.0
! ! 			end if
!
! ! 			if ((maskP(p,pp) .eq. 50.0) .and. (p .eq. xn-1)) then
! ! 				psi_bc(p,pp+1) = 0.1e-5
! ! 			end if
!
! ! 			if ((maskP(p,pp) .eq. 50.0) .and. (p .eq. xn-2)) then
! ! 				psi_bc(p,pp+1) = 0.2e-5
! ! 			end if
! !
! ! 			if ((maskP(p,pp) .eq. 50.0) .and. (p .eq. xn-3)) then
! ! 				psi_bc(p,pp+1) = 0.3e-5
! ! 			end if
! !
! ! 			if ((maskP(p,pp) .eq. 50.0) .and. (p .le. xn-4)) then
! ! 				psi_bc(p,pp+1) = 0.4e-5
! ! 			end if
!
!
!
!
!
! 		end do
! 	end do
!
	psi_bc(f_index1,:) = 0.0


! ! SINGLE FRACTURE BOX
! psi_bc(1,:) = 0.0
! psi_bc(xn,:) = 0.0
!
! do p=2,xn
! 	do pp=2,yn-1
!
! 		! top of sediment
! 		if ((maskP(p,pp) .eq. 25.0)) then
! 			psi_bc(p,pp+1) = 0.0
! 		end if
!
! 		! top of sediment
! 		if ((maskP(p,pp) .eq. 100.0)) then
! 			psi_bc(p,pp-1) = 0.0
! 		end if
!
! 	end do
! end do



	return
end function psi_bc























function psi_mod(psi_in)

	use globals
	real(4) :: psi_in(xn,yn), psi_mod(xn,yn)
	integer :: p,pp

	psi_mod = psi_in

	do pp=1,yn
    do p=1,xn
		if (mask(p,pp) .eq. 0.0) then
			psi_mod(p,pp) = dx*(5.0e-12)*1000.0 - dx*(5.0e-12)*0.0
		end if
	end do
	end do




	return
end function psi_mod

end module initialize
