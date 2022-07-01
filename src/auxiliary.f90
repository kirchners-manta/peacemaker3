!=========================================================================================
! Peacemaker -- A Quantum Cluster Equilibrium Code.
!
! Copyright 2004-2006 Barbara Kirchner, University of Bonn
! Copyright 2007-2012 Barbara Kirchner, University of Leipzig
! Copyright 2013-2022 Barbara Kirchner, University of Bonn
!
! This file is part of Peacemaker.
!
! Peacemaker is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Peacemaker is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Peacemaker.  If not, see <http://www.gnu.org/licenses/>
!=========================================================================================
! This module provides auxiliary functions that don't really fit anywhere else.
module auxiliary
    use iso_varying_string
    use kinds
    use constants
    implicit none
    private
    !=====================================================================================
    ! Public entities.
    public :: range_t
    public :: set_range
    public :: check_range
    public :: string2int
    public :: string2real
    public :: array_sample
    public :: write_range
    public :: swap
    public :: ln_factorial
    public :: derivative
    public :: progress_bar
    !=====================================================================================
    ! Interfaces
    interface array_sample
        module procedure array_sample_real
        module procedure array_sample_integer
    end interface array_sample
    !=====================================================================================
    ! Data type representing an interval.
    type :: range_t
        real(dp) :: first, last, delta
        integer :: num
    end type range_t
    !=====================================================================================
    contains
        !=================================================================================
        ! Converts a string to real. Returns ios == 0, if the conversion was succesful.
        function string2real(s, ios)
            type(varying_string), intent(in) :: s
            integer, intent(out) :: ios
            real(dp):: string2real

            character(:), allocatable :: tmp

            allocate(tmp, source = char(s))
            read(tmp, *, iostat = ios) string2real
            deallocate(tmp)
        end function string2real
        !=================================================================================
        ! Converts a string to integer. Returns ios == 0, if the conversion was succesful.
        function string2int(s, ios)
            type(varying_string), intent(in) :: s
            integer, intent(out) :: ios
            integer:: string2int

            character(:), allocatable :: tmp

            allocate(tmp, source = char(s))
            read(tmp, *, iostat = ios) string2int
            deallocate(tmp)
        end function string2int
        !=================================================================================
        ! Writes the first and the last value of an array to screen. If the array is
        ! bigger than 2 values, dots are inserted in between. Does not advance the ouput.
        subroutine array_sample_real(array)
            real(dp), dimension(:), intent(in) :: array
            integer:: n

            n = size(array)
            if (n == 0) then
                write(*, '(G0.6)', advance = "no") "-"
            else if (n == 1) then
                write(*, '(G0.6)', advance = "no") array(1)
            else if (n == 2) then
                write(*, '(G0.6, 1X, G0.6)', advance = "no") array(1), array(2)
            else
                write(*, '(G0.6, 1X, A, 1X, G0.6)', advance = "no") &
                    array(1), "...", array(n)
            end if
        end subroutine array_sample_real
        !=================================================================================
        ! Pretty prints a range.
        subroutine write_range(r)
            type(range_t), intent(in) :: r

            if (r%num == 1) then
                write(*, '(G0.6)', advance = "no") r%first
            else if (r%num == 2) then
                write(*, '(G0.6, 1X, G0.6)', advance = "no") r%first, r%last
            else
                write(*, '(G0.6, 1X, A, 1X, G0.6)', advance = "no") r%first, "...", r%last
            end if
        end subroutine write_range
        !=================================================================================
        subroutine set_range(r, first, last, num)
            type(range_t), intent(out) :: r
            real(dp), intent(in) :: first
            real(dp), intent(in) :: last
            integer, intent(in) :: num

            r%first = first
            r%last = last
            r%num = num
            if (num /= 1) then
                r%delta = (last-first)/real(num-1, dp)
            else
                r%delta = 0.0_dp
            end if
        end subroutine set_range
        !=================================================================================
        function check_range(r)
            logical:: check_range
            type(range_t), intent(in) :: r

            check_range = .true.
            if (r%last < r%first) check_range = .false.
            if (r%num < 1) check_range = .false.
        end function check_range
        !=================================================================================
        subroutine array_sample_integer(array)
            integer, dimension(:), intent(in) :: array
            integer:: n

            n = size(array)
            if (n == 1) then
                write(*, '(G0)', advance = "no") array(1)
            else if (n == 2) then
                write(*, '(G0, 1X, G0)', advance = "no") array(1), array(2)
            else
                write(*, '(G0, 1X, A, 1X, G0)', advance = "no") array(1), "...", array(n)
            end if
        end subroutine array_sample_integer
        !=================================================================================
        ! Swaps two reals.
        pure subroutine swap(r1, r2)
            real(dp), intent(inout) :: r1
            real(dp), intent(inout) :: r2
            real(dp):: tmp

            tmp = r1
            r1 = r2
            r2 = tmp
        end subroutine swap
        !=================================================================================
        ! Returns the logarithm of the factorial of x. Uses the identity x! = Gamma(x+1)
        ! Uses the Lanczos approximation with P. Godfreys coefficients for the calculation
        ! of the gamma function. Adapated from Numerical recipies.
        pure function ln_factorial(x)
            real(dp), intent(in) :: x
            real(dp):: ln_factorial

            real(dp), dimension(0:14), parameter :: c = [ &
                0.99999999999999709182_dp, &
                57.156235665862923517_dp, &
                -59.597960355475491248_dp, &
                14.136097974741747174_dp, &
                -0.49191381609762019978_dp, &
                0.33994649984811888699e-4_dp, &
                0.46523628927048575665e-4_dp, &
                -0.98374475304879564677e-4_dp, &
                0.15808870322491248884e-3_dp, &
                -0.21026444172410488319e-3_dp, &
                0.21743961811521264320e-3_dp, &
                -0.16431810653676389022e-3_dp, &
                0.84418223983852743293e-4_dp, &
                -0.26190838401581408670e-4_dp, &
                0.36899182659531622704e-5_dp]
            real(dp), parameter :: sqrt2pi = sqrt(2.0_dp*pi)
            real(dp), parameter :: g = 607.0_dp/128.0_dp

            real(dp):: tmp
            real(dp):: series
            real(dp):: denominator
            integer:: i

            tmp = x + g + 0.5_dp;
            tmp = (x+0.5_dp)*log(tmp) - tmp

            denominator = x
            series = c(0)
            do i = 1, 14
                denominator = denominator + 1.0_dp
                series = series + c(i)/denominator
            end do

            ln_factorial = tmp + log(sqrt2pi*series)
        end function ln_factorial
        !=================================================================================
        ! Numerically approximates the derivative of y with respect to x. Uses central
        ! differentiation.
        function derivative(y, x) result(d)
            real(dp), dimension(:), intent(in) :: y
            real(dp), dimension(:), intent(in) :: x
            real(dp), dimension(size(y)) :: d

            integer:: i

            if (size(y) == 1) then
                d(1) = 0.0_dp
            else
                d(1) = (y(2) - y(1))/(x(2) - x(1))
                do i = 2, size(y)-1
                    d(i) = (y(i+1) - y(i-1))/(x(i+1) - x(i-1))
                end do
                d(size(y)) = (y(size(y)) - y(size(y)-1))/(x(size(y)) - x(size(y)-1))
            end if
        end function derivative
        !=================================================================================
        ! Prints a fancy progress bar. Yay!
        subroutine progress_bar(current, total, active, newline)
            use iso_fortran_env
            integer, intent(in) :: current, total
            integer :: nbars
            integer :: i
            logical, intent(in) :: active
            logical, optional :: newline

            if (.not. active) return

            nbars = int(real(current*100, dp) / total)

            write(error_unit, '(A1, 4X, "[")', advance = 'no') achar(13)
            do i = 1, nbars - 1
                write(error_unit, '("=")', advance = 'no')
            end do
            write(error_unit, '(">")', advance = 'no')
            do i = nbars + 1, 100
                write(error_unit, '(" ")', advance = 'no')
            end do
            write(error_unit, '("]", 2X, I3, "% ")', advance = 'no') nbars
            if (present(newline)) then
                if (newline) then
                    write(error_unit, *)
                end if
            end if
        end subroutine progress_bar
        !=================================================================================
end module auxiliary
!=========================================================================================
