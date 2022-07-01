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
! This module collects headers, footers and other information that is written to the
! screen.
!=========================================================================================
module info
    implicit none
    private

    public :: print_usage_info
    public :: print_welcome_info
    public :: print_citation_info

    contains
        !=================================================================================
        ! Prints a message how to run peacemaker.
        subroutine print_usage_info()
            write(*, "(A)") "Usage: peacemaker [input] [clusterset]"
        end subroutine print_usage_info
        !=================================================================================
        ! Prints the header that is written at the start of the program.
        subroutine print_welcome_info()
            write(*, *)
            write(*, "(A)") "    Peacemaker 3.1.0 - a Quantum Cluster Equilibrium (QCE) code:             "
            write(*, *)
            write(*, "(A)") "        Copyright 2004-2006 Barbara Kirchner, University of Bonn"
            write(*, "(A)") "        Copyright 2007-2012 Barbara Kirchner, University of Leipzig          "
            write(*, "(A)") "        Copyright 2013-2022 Barbara Kirchner, University of Bonn             "
            write(*, *)
            write(*, "(A)") "        This program is free software: you can redistribute it and/or modify "
            write(*, "(A)") "        it under the terms of the GNU General Public License as published by "
            write(*, "(A)") "        the Free Software Foundation, either version 3 of the License, or    "
            write(*, "(A)") "        (at your option) any later version.                                  "
            write(*, *)
            write(*, "(A)") "        This program is distributed in the hope that it will be useful,      "
            write(*, "(A)") "        but WITHOUT ANY WARRANTY; without even the implied warranty of       "
            write(*, "(A)") "        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        "
            write(*, "(A)") "        GNU General Public License for more details.                         "
            write(*, *)
            write(*, "(A)") "        You should have received a copy of the GNU General Public License    "
            write(*, "(A)") "        along with this program.  If not, see <http://www.gnu.org/licenses/>."
            write(*, *)
            write(*, "(A)") "    The Peacemaker team (in alphabetical order):                             "
            write(*, *)
            write(*, "(A)") "        * Johannes Ingenmey                                                  "
            write(*, "(A)") "        * Eva Perlt                                                          "
            write(*, "(A)") "        * Michael von Domaros                                                "
            write(*, "(A)")
            write(*, "(A)") "        * Barbara Kirchner                                                   "
            write(*, *)
            write(*, "(A1, A, A)") achar(27), '[1m', &
                            "    Please let us know that you use Peacemaker by sending a short email to   "
            write(*, "(A)") "    qce@thch.uni-bonn.de. This will help us reaching out to you in case of   "
            write(*, "(A, A1, A)") &
                            "    critical bug fixes or new releases.                                      ", achar(27), '[0m'
            write(*, *)
        end subroutine print_welcome_info
        !=================================================================================
        ! Prints citation info.
        subroutine print_citation_info()
            write(*, '(A1, A)', advance='no') achar(27), '[1m'
            write(*, "(A)") '    Please always cite:'
            write(*, *)
            write(*, "(A)") '        * Michael von Domaros, Eva Perlt, Johannes Ingenmey, Gwydyon Marchelli, Barbara Kirchner:'
            write(*, "(A)") '          "Peacemaker 2: Making clusters talk about binary mixtures and neat liquids".'
            write(*, "(A)") '          SoftwareX (2018).'
            write(*, *)
            write(*, "(A)") '        * Barbara Kirchner, Christian Spickermann, Sebastian B. C. Lehmann, Eva Perlt,'
            write(*, "(A)") '          Johanna Langner, Michael von Domaros, Patricia Reuther, Frank Uhlig, Miriam Kohagen,'
            write(*, "(A)") '          Marc Brüssel:'
            write(*, "(A)") '          "What can clusters tell us about the bulk? PEACEMAKER: ' // &
                'Extended quantum cluster equilibrium calculations".'
            write(*, "(A)") '          Comp. Physics Commun. 182 (2011), 1428.'
            write(*, *)
            write(*, "(A)") '        * Barbara Kirchner:'
            write(*, "(A)") '          "Cooperative versus dispersion effects: ' // &
                'What is more important in an associated liquid such as water?".'
            write(*, "(A)") '          J. Chem. Phys. 123 (2005), 204116.'
            write(*, *)
            write(*, "(A)") '    For binary mixtures, please cite:'
            write(*, *)
            write(*, "(A)") '        * Marc Brüssel, Eva Perlt, Sebastian B. C. Lehmann, Michael von Domaros, Barbara Kirchner:'
            write(*, "(A)") '          "Binary systems from quantum cluster equilibrium theory".'
            write(*, "(A)") '          J. Chem. Phys. 135 (2011), 194113.'
            write(*, *)
            write(*, "(A)") '    For anharmonicites, please cite:'
            write(*, *)
            write(*, "(A)") '        * Michael von Domaros, Eva Perlt:'
            write(*, "(A)") '          "Anharmonic effects in the quantum cluster equilibrium method".'
            write(*, "(A)") '          J. Chem. Phys. 146 (2017), 154502.'
            write(*, *)
            write(*, "(A)") '    For acid constants, please cite:'
            write(*, *)
            write(*, "(A)") '        * Eva Perlt, Michael von Domaros, Barbara Kirchner, Ralf Ludwig, Frank Weinhold".'
            write(*, "(A)") '          "Predicting the ionic product of water".'
            write(*, "(A)") '          Sci. Rep. 7 (2017), 10244.'
            write(*, '(A1, A)') achar(27), '[0m'
        end subroutine print_citation_info
        !=================================================================================
end module info
!=========================================================================================
