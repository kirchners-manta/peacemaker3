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
! This module provides two very simple routines for error handling.
!=========================================================================================
module error
    use, intrinsic :: iso_fortran_env
    use iso_varying_string
    implicit none
    private
    !=====================================================================================
    ! Public entities.
    public :: pmk_error
    public :: pmk_missing_key_error
    public :: pmk_argument_error
    public :: pmk_argument_count_error
    public :: pmk_unphysical_argument_error
    public :: pmk_illegal_range_error
    public :: pmk_warning
    !=====================================================================================
    ! Interfaces.
    interface pmk_unphysical_argument_error
        module procedure pmk_unphysical_argument_error_char
        module procedure pmk_unphysical_argument_error_string
    end interface pmk_unphysical_argument_error
    !=====================================================================================
    interface pmk_argument_error
        module procedure pmk_argument_error_char
        module procedure pmk_argument_error_string
    end interface pmk_argument_error
    !=====================================================================================
    interface pmk_argument_count_error
        module procedure pmk_argument_count_error_char
        module procedure pmk_argument_count_error_string
    end interface pmk_argument_count_error
    !=====================================================================================
    ! Error unit. Change this, if iso_fortran_env is not available on your system.
    integer, parameter :: my_unit = error_unit
    !=====================================================================================
    contains
        !=================================================================================
        ! Writes an error message and stops execution.
        subroutine pmk_error(what)
            character(*), intent(in) :: what

            write(my_unit, '(A,1X,A)') "error;", what
            stop
        end subroutine pmk_error
        !=================================================================================
        ! Writes an error message indicating that an argument is invalid for the given
        ! keyword/section. Terminates execution.
        subroutine pmk_argument_error_string(key, section)
            character(*), intent(in) :: key
            type(varying_string), intent(in) :: section

            call pmk_argument_error_char(key, char(section))
        end subroutine pmk_argument_error_string
        !=================================================================================
        subroutine pmk_argument_error_char(key, section)
            character(*), intent(in) :: key, section

            call pmk_error("illegal argument in for keyword '" // key // &
                "' in section [" // section // "]")
        end subroutine pmk_argument_error_char
        !=================================================================================
        ! Writes an error message indicating that the argument count for the given
        ! keyword/section is invalid. Terminates execution.
        subroutine pmk_argument_count_error_string(key, section)
            character(*), intent(in) :: key
            type(varying_string), intent(in) :: section

            call pmk_argument_count_error(key, char(section))
        end subroutine pmk_argument_count_error_string
        !=================================================================================
        subroutine pmk_argument_count_error_char(key, section)
            character(*), intent(in) :: key, section

            call pmk_error("illegal number of arguments for keyword '" // key // &
                "' in section [" // section // "]")
        end subroutine pmk_argument_count_error_char
        !=================================================================================
        ! Writes an error message indicating that an argument has an unphsyical value.
        ! Terminates execution.
        subroutine pmk_unphysical_argument_error_char(key, section)
            character(*), intent(in) :: key, section

            call pmk_error("unphsyical value for keyword '" // key // &
                "' in section [" // section // "]")
        end subroutine pmk_unphysical_argument_error_char
        !=================================================================================
        subroutine pmk_unphysical_argument_error_string(key, section)
            character(*), intent(in) :: key
            type(varying_string), intent(in) :: section

            call pmk_unphysical_argument_error_char(key, char(section))
        end subroutine pmk_unphysical_argument_error_string
        !=================================================================================
        ! Writes an error message indicating that the given keyword is missing. Terminates
        ! execution.
        subroutine pmk_missing_key_error(key, section)
            character(*), intent(in) :: key
            type(varying_string), intent(in) :: section

            call pmk_error("missing keyword '" // key // "' in section [" // &
                char(section) // "]")
        end subroutine pmk_missing_key_error
        !=================================================================================
        ! Writes an error message indicating that there was something wrong with the
        ! range specification.
        subroutine pmk_illegal_range_error(key, section)
            character(*), intent(in) :: key, section

            call pmk_error("illegal range specification for keyword '" // key // &
                "' in section [" // section // "]")
        end subroutine pmk_illegal_range_error
        !=================================================================================
        ! Writes a warning.
        subroutine pmk_warning(what)
            character(*), intent(in) :: what

            write(my_unit, '(A,1X,A)') "warning;", what
        end subroutine pmk_warning
        !=================================================================================
end module error
!=========================================================================================
