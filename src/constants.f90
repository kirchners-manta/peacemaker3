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
! This module provides access to natural constants.
!=========================================================================================
module constants
    use kinds
    implicit none
    !=====================================================================================
    ! Natural constants.
    real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)
    real(dp), parameter :: planck = 6.62606957e-34_dp ! J*s
    real(dp), parameter :: avogadro = 6.0221413e23_dp ! 1/mol
    real(dp), parameter :: kb = 1.3806488e-23_dp ! J/K
    real(dp), parameter :: speed_of_light = 299792458.0_dp ! m/s
    real(dp), parameter :: amu = 1.660538921e-27_dp ! kg
    !=====================================================================================
    ! Derived constants for convenience.
    real(dp), parameter :: gas_constant = avogadro*kb
    real(dp), parameter :: hbar = planck/(2.0_dp*pi)
    !=====================================================================================
    ! Auxiliary constants.
    real(dp), parameter :: global_eps = 1.0e-10_dp
    !=====================================================================================
end module constants
!=========================================================================================
