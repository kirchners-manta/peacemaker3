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
! This module provides kind parameters.
!=========================================================================================
module kinds
    implicit none
    !=====================================================================================
    ! Real kind parameters.
    integer, parameter :: sp = selected_real_kind(p = 7)
    integer, parameter :: dp = selected_real_kind(p = 15)
    integer, parameter :: qp = selected_real_kind(p = 33)
    !=====================================================================================
end module kinds
!=========================================================================================
