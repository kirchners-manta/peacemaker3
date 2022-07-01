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
! This module provides atomic data.
module atomic_data
    use kinds
    implicit none
    private
    !=====================================================================================
    ! Public entities.
    public :: periodic_table
    !=====================================================================================
    ! Data type storing atomic data.
    type :: element_t
        integer :: atomic_number
        character(2) :: atomic_symbol
        real(dp) :: atomic_mass ! in amu
    end type element_t
    !=====================================================================================
    ! Class that provides access to the periodic table of elements.
    type :: periodic_table_t
        integer :: n
        type(element_t), dimension(:), allocatable :: elements
    contains
        procedure, public :: init => periodic_table_init
        procedure, public :: mass => periodic_table_mass
    end type periodic_table_t
    !=====================================================================================
    ! Instantiation of the periodic_table_t type.
    type(periodic_table_t):: periodic_table
    !=====================================================================================
    contains
        !=================================================================================
        ! Initializes the periodic table.
        subroutine periodic_table_init(table)
            use constants
            class(periodic_table_t):: table

            table%n = 110
            allocate(table%elements(table%n))

            table%elements( 1) = element_t( 1, "H ",  1.008)
            table%elements( 2) = element_t( 2, "He",  4.003)
            table%elements( 3) = element_t( 3, "Li",  6.94 )
            table%elements( 4) = element_t( 4, "Be",  9.012)
            table%elements( 5) = element_t( 5, "B ", 10.81 )
            table%elements( 6) = element_t( 6, "C ", 12.01 )
            table%elements( 7) = element_t( 7, "N ", 14.01 )
            table%elements( 8) = element_t( 8, "O ", 16.00 )
            table%elements( 9) = element_t( 9, "F ", 19.00 )
            table%elements(10) = element_t(10, "Ne", 20.18 )
            table%elements(11) = element_t(11, "Na", 22.99 )
            table%elements(12) = element_t(12, "Mg", 24.31 )
            table%elements(13) = element_t(13, "Al", 26.98 )
            table%elements(14) = element_t(14, "Si", 28.09 )
            table%elements(15) = element_t(15, "P ", 30.97 )
            table%elements(16) = element_t(16, "S ", 32.06 )
            table%elements(17) = element_t(17, "Cl", 35.45 )
            table%elements(18) = element_t(18, "Ar", 39.95 )
            table%elements(19) = element_t(19, "Ac", 227.028)
            table%elements(20) = element_t(20, "Am", 243)
            table%elements(21) = element_t(21, "Sb", 121.757)
            table%elements(22) = element_t(22, "As", 74.92159)
            table%elements(23) = element_t(23, "At", 210)
            table%elements(24) = element_t(24, "Ba", 137.327)
            table%elements(25) = element_t(25, "Bk", 247)
            table%elements(26) = element_t(26, "Bi", 208.98037)
            table%elements(27) = element_t(27, "Bh", 262)
            table%elements(28) = element_t(28, "Br", 79.904)
            table%elements(29) = element_t(29, "Cd", 112.411)
            table%elements(30) = element_t(30, "Ca", 40.078)
            table%elements(31) = element_t(31, "Cf", 251)
            table%elements(32) = element_t(32, "Ce", 140.115)
            table%elements(33) = element_t(33, "Cs", 132.90543)
            table%elements(34) = element_t(34, "Cr", 51.9961)
            table%elements(35) = element_t(35, "Co", 58.93320)
            table%elements(36) = element_t(36, "Cu", 63.546)
            table%elements(37) = element_t(37, "Cm", 247)
            table%elements(38) = element_t(38, "Db", 262)
            table%elements(39) = element_t(39, "Dy", 162.50)
            table%elements(40) = element_t(40, "Es", 252)
            table%elements(41) = element_t(41, "Er", 167.26)
            table%elements(42) = element_t(42, "Eu", 151.965)
            table%elements(43) = element_t(43, "Fm", 257)
            table%elements(44) = element_t(44, "Fr", 223)
            table%elements(45) = element_t(45, "Gd", 157.25)
            table%elements(46) = element_t(46, "Ga", 69.723)
            table%elements(47) = element_t(47, "Ge", 72.61)
            table%elements(48) = element_t(48, "Au", 196.96654)
            table%elements(49) = element_t(49, "Hf", 178.49)
            table%elements(50) = element_t(50, "Hs", 265)
            table%elements(51) = element_t(51, "Ho", 164.93032)
            table%elements(52) = element_t(52, "In", 114.82)
            table%elements(53) = element_t(53, "I" , 126.90447)
            table%elements(54) = element_t(54, "Ir", 192.22)
            table%elements(55) = element_t(55, "Fe", 55.847)
            table%elements(56) = element_t(56, "Kr", 83.80)
            table%elements(57) = element_t(57, "La", 138.9055)
            table%elements(58) = element_t(58, "Lr", 262)
            table%elements(59) = element_t(59, "Pb", 207.2)
            table%elements(60) = element_t(60, "Lu", 174.967)
            table%elements(61) = element_t(61, "Mn", 54.93805)
            table%elements(62) = element_t(62, "Mt", 266)
            table%elements(63) = element_t(63, "Md", 258)
            table%elements(64) = element_t(64, "Hg", 200.59)
            table%elements(65) = element_t(65, "Mo", 95.94)
            table%elements(66) = element_t(66, "Nd", 144.24)
            table%elements(67) = element_t(67, "Np", 237.048)
            table%elements(68) = element_t(68, "Ni", 58.6934)
            table%elements(69) = element_t(69, "Nb", 92.90638)
            table%elements(70) = element_t(70, "No", 259)
            table%elements(71) = element_t(71, "Os", 190.2)
            table%elements(72) = element_t(72, "Pd", 106.42)
            table%elements(73) = element_t(73, "Pt", 195.08)
            table%elements(74) = element_t(74, "Pu", 244)
            table%elements(75) = element_t(75, "Po", 209)
            table%elements(76) = element_t(76, "K" , 39.0983)
            table%elements(77) = element_t(77, "Pr", 140.90765)
            table%elements(78) = element_t(78, "Pm", 145)
            table%elements(79) = element_t(79, "Pa", 231.0359)
            table%elements(80) = element_t(80, "Ra", 226.025)
            table%elements(81) = element_t(81, "Rn", 222)
            table%elements(82) = element_t(82, "Re", 186.207)
            table%elements(83) = element_t(83, "Rh", 102.90550)
            table%elements(84) = element_t(84, "Rb", 85.4678)
            table%elements(85) = element_t(85, "Ru", 101.07)
            table%elements(86) = element_t(86, "Rf", 261)
            table%elements(87) = element_t(87, "Sm", 150.36)
            table%elements(88) = element_t(88, "Sc", 44.955910)
            table%elements(89) = element_t(89, "Sg", 263)
            table%elements(90) = element_t(90, "Se", 78.96)
            table%elements(91) = element_t(91, "Ag", 107.8682)
            table%elements(92) = element_t(92, "Sr", 87.62)
            table%elements(93) = element_t(93, "Ta", 180.9479)
            table%elements(94) = element_t(94, "Tc", 98)
            table%elements(95) = element_t(95, "Te", 127.60)
            table%elements(96) = element_t(96, "Tb", 158.92534)
            table%elements(97) = element_t(97, "Tl", 204.3833)
            table%elements(98) = element_t(98, "Th", 232.0381)
            table%elements(99) = element_t(99, "Tm", 168.93421)
            table%elements(100) = element_t(100, "Sn", 118.710)
            table%elements(101) = element_t(101, "Ti", 47.88)
            table%elements(102) = element_t(102, "W" , 183.85)
            table%elements(103) = element_t(103, "U" , 238.0289)
            table%elements(104) = element_t(104, "V" , 50.9415)
            table%elements(105) = element_t(105, "Xe", 131.29)
            table%elements(106) = element_t(106, "Yb", 173.04)
            table%elements(107) = element_t(107, "Y" , 88.90585)
            table%elements(108) = element_t(108, "Zn", 65.39)
            table%elements(109) = element_t(109, "Zr", 91.224)
            table%elements(110) = element_t(110, "Xx", 0.0)
            
        end subroutine periodic_table_init
        !=================================================================================
        ! Returns the mass of an element.
        function periodic_table_mass(table, atomic_symbol) result(mass)
            use error
            class(periodic_table_t), intent(in) :: table
            character(2), intent(in) :: atomic_symbol
            real(dp):: mass

            integer:: i
            do i = 1, table%n
                if (table%elements(i)%atomic_symbol == atomic_symbol) then
                    mass = table%elements(i)%atomic_mass
                    return
                end if
            end do

            mass = 0.0_dp
            call pmk_error("unknown element '" // atomic_symbol // "'")
        end function periodic_table_mass
        !=================================================================================
end module atomic_data
!=========================================================================================
