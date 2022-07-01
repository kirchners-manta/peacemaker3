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
! This module provides global data.
module shared_data
    use kinds
    use auxiliary, only: range_t
    implicit none
    private
    !=====================================================================================
    ! Public entities
    public :: global_data
    !=====================================================================================
    ! Data type storing global data.
    type :: global_data_t
        integer :: qce_iterations, newton_iterations, grid_iterations, optimizer, nconverged
        integer, dimension(:), allocatable :: degree
        real(dp) :: press, mtot, vexcl
        real(dp), dimension(:), allocatable :: ntot
        real(dp) :: max_deviation, vdamp, rotor_cutoff
        type(range_t) :: amf, bxv, temp
        type(range_t) :: amf_temp, bxv_temp
        real(dp), dimension(:), allocatable :: monomer_amounts
        logical :: progress_bar, imode
    end type global_data_t
    !=====================================================================================
    ! Instantiation of the data type specified above. Data in here is shared across
    ! threads and shall not be modified during a QCE calculation.
    type(global_data_t):: global_data
    !=====================================================================================
end module shared_data
!=========================================================================================
