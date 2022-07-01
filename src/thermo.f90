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
! This module implements all thermodynamic functions and their output.
module thermo
    use kinds
    implicit none
    private
    !=====================================================================================
    ! Public entities.
    public :: calculate_thermo, calculate_cluster_thermo
    !=====================================================================================
    contains
        !=================================================================================
        ! Calculates the system partition function.
        ! Attention: This subroutine is missing the summand arising from the particle
        ! indistinguishability.
        subroutine calculate_lnq_sys(ntemp, nclust, pop, lnq_clust, lnq_sys)
            integer, intent(in) :: ntemp, nclust
            real(dp), dimension(ntemp, nclust), intent(in) :: pop, lnq_clust
            real(dp), dimension(ntemp), intent(out) :: lnq_sys

            integer :: itemp, iclust

            lnq_sys = 0.0_dp
            do itemp = 1, ntemp
                do iclust = 1, nclust
                    lnq_sys(itemp) = lnq_sys(itemp) + &
                        pop(itemp, iclust)*lnq_clust(itemp, iclust)
                end do
            end do
        end subroutine calculate_lnq_sys
        !=================================================================================
        ! Adds the part of the system partition function arising from particle
        ! indistinguishability.
        subroutine add_lnq_indi(ntemp, nclust, pop, lnq_sys)
            use auxiliary, only: ln_factorial
            integer, intent(in) :: ntemp, nclust
            real(dp), dimension(ntemp, nclust), intent(in) :: pop
            real(dp), dimension(ntemp), intent(inout) :: lnq_sys

            integer :: itemp, iclust
            
            do itemp = 1, ntemp
                do iclust = 1, nclust
                    lnq_sys(itemp) = lnq_sys(itemp) - ln_factorial(pop(itemp, iclust))
                end do
            end do
        end subroutine add_lnq_indi
        !=================================================================================
        ! Calculates the Helmholtz free energy.
        subroutine calculate_helmholtz_energy(ntemp, temp, lnq, a)
            use constants, only: kb

            integer, intent(in) :: ntemp
            real(dp), dimension(ntemp), intent(in) :: temp, lnq
            real(dp), dimension(ntemp), intent(out) :: a

            integer :: itemp
            do itemp = 1, ntemp
                a(itemp) = -kb*temp(itemp)*lnq(itemp)
            end do
        end subroutine calculate_helmholtz_energy
        !=================================================================================
        ! Calculates the Gibbs free enthalpy.
        subroutine calculate_gibbs_enthalpy(ntemp, temp, lnq, vol, press, g)
            integer, intent(in) :: ntemp
            real(dp), dimension(ntemp), intent(in) :: temp, lnq, vol
            real(dp), intent(in) :: press
            real(dp), dimension(ntemp), intent(out) :: g

            real(dp), dimension(ntemp) :: a
            integer :: itemp

            call calculate_helmholtz_energy(ntemp, temp, lnq, a)
            do itemp = 1, ntemp
                g(itemp) = a(itemp) + press*vol(itemp)
            end do
        end subroutine calculate_gibbs_enthalpy
        !=================================================================================
        ! Calculates the internal energy.
        subroutine calculate_internal_energy(ntemp, temp, dlnq, u)
            use constants, only: kb

            integer, intent(in) :: ntemp
            real(dp), dimension(ntemp), intent(in) :: temp, dlnq
            real(dp), dimension(ntemp), intent(out) :: u

            integer :: itemp

            do itemp = 1, ntemp
                u(itemp) = kb*temp(itemp)**2*dlnq(itemp)
            end do
        end subroutine calculate_internal_energy
        !=================================================================================
        ! Calculates the enthalpy.
        subroutine calculate_enthalpy(ntemp, temp, dlnq, vol, press, h)
            integer, intent(in) :: ntemp
            real(dp), dimension(ntemp), intent(in) :: temp, dlnq, vol
            real(dp), intent(in) :: press
            real(dp), dimension(ntemp), intent(out) :: h

            integer :: itemp
            real(dp), dimension(ntemp) :: u

            call calculate_internal_energy(ntemp, temp, dlnq, u)
            do itemp = 1, ntemp
                h(itemp) = u(itemp) + press*vol(itemp)
            end do
        end subroutine calculate_enthalpy
        !=================================================================================
        ! Calculates the entropy.
        subroutine calculate_entropy(ntemp, temp, lnq, dlnq, s)
            integer, intent(in) :: ntemp
            real(dp), dimension(ntemp), intent(in) :: temp, lnq, dlnq
            real(dp), dimension(ntemp), intent(out) :: s

            real(dp), dimension(ntemp) :: u, a

            call calculate_internal_energy(ntemp, temp, dlnq, u)
            call calculate_helmholtz_energy(ntemp, temp, lnq, a)
            s = (u-a)/temp
        end subroutine calculate_entropy
        !=================================================================================
        ! Calculates the volumetric expansion coefficient.
        subroutine calculate_expansion_coefficient(ntemp, vol, dvol, alpha)
            integer, intent(in) :: ntemp
            real(dp), dimension(ntemp), intent(in) :: vol, dvol
            real(dp), dimension(ntemp), intent(out) :: alpha

            integer :: itemp

            if (ntemp > 1) then
                do itemp = 1, ntemp
                    alpha(itemp) = 1.0_dp/vol(itemp)*dvol(itemp)
                end do
            else
                alpha = 0.0_dp
            end if
        end subroutine calculate_expansion_coefficient
        !=================================================================================
        ! Calculates the heat capacity at constant volume.
        subroutine calculate_cv(ntemp, temp, dlnq, ddlnq, cv)
            use constants, only: kb

            integer, intent(in) :: ntemp
            real(dp), dimension(ntemp), intent(in) :: temp, dlnq, ddlnq
            real(dp), dimension(ntemp), intent(out) :: cv

            integer :: itemp

            do itemp = 1, ntemp
                cv(itemp) = 2.0_dp*kb*temp(itemp)*dlnq(itemp) + &
                    kb*temp(itemp)**2*ddlnq(itemp)
            end do
        end subroutine calculate_cv
        !=================================================================================
        ! Calculates the heat capacity at constant pressure.
        subroutine calculate_cp(ntemp, temp, dlnq, ddlnq, dvol, press, cp)
            integer, intent(in) :: ntemp
            real(dp), dimension(ntemp), intent(in) :: temp, dlnq, ddlnq, dvol
            real(dp), intent(in) :: press
            real(dp), dimension(ntemp), intent(out) :: cp

            integer :: itemp
            real(dp), dimension(ntemp) :: cv

            call calculate_cv(ntemp, temp, dlnq, ddlnq, cv)
            if (ntemp > 1) then
                do itemp = 1, ntemp
                    cp(itemp) = cv(itemp) + dvol(itemp)*press
                end do
            else
                cp = 0.0_dp
            end if

        end subroutine calculate_cp
        !=================================================================================
        ! Writes volume, exclusion volume and status code to file.
        subroutine write_volume(ntemp, converged, temp, vol, vexcl, alpha, solution)
            integer, intent(in) :: ntemp
            real(dp), dimension(ntemp), intent(in) :: temp, vol, alpha
            real(dp), intent(in) :: vexcl
            integer, dimension(ntemp), intent(in) :: solution
            logical, dimension(ntemp), intent(in) :: converged

            integer :: myunit
            integer :: itemp

            open(newunit = myunit, action = "write", status = "unknown", &
                file = "volume.dat")
            write(myunit, '(A1,A12,4(1X,A13))') &
                "#", &
                "T/K", &
                "V/dm^3", &
                "V(excl)/dm^3", &
                "alpha/K^-1", &
                "status"
            do itemp = 1, ntemp
                if (converged(itemp)) then
                    write(myunit, '(4(ES13.6,1X),I35)') &
                        temp(itemp), &
                        1.0e3_dp*vol(itemp), &
                        1.0e3_dp*vexcl, &
                        alpha(itemp), &
                        solution(itemp)
                end if
            end do
            close(myunit)
        end subroutine write_volume
        !=================================================================================
        ! Writes thermodynamic functions that do not depend on any derivatives.
        subroutine write_thermo0(ntemp, converged, temp, a, g)
            integer, intent(in) :: ntemp
            logical, dimension(ntemp), intent(in) :: converged
            real(dp), dimension(ntemp), intent(in) :: temp, a, g

            integer :: myunit
            integer :: itemp

            open(newunit = myunit, action = "write", status = "unknown", &
                file = "thermo0.dat")
            write(myunit, '(A1,A12,2(1X,A13))') &
                "#", &
                "T/K", &
                "A/kJ", &
                "G/kJ"
            do itemp = 1, ntemp
                if (converged(itemp)) then
                    write(myunit, '(ES20.12,2(1X,ES20.12))') &
                        temp(itemp), &
                        1.0e-3_dp*a(itemp), &
                        1.0e-3_dp*g(itemp)
                end if
            end do
            close(myunit)
        end subroutine write_thermo0
        !=================================================================================
        ! Writes thermodynamic functions that depend on first derivatives.
        subroutine write_thermo1(ntemp, converged, temp, u, h, s)
            integer, intent(in) :: ntemp
            logical, dimension(ntemp), intent(in) :: converged
            real(dp), dimension(ntemp), intent(in) :: temp, u, h, s

            integer :: myunit
            integer :: itemp

            open(newunit = myunit, action = "write", status = "unknown", &
                file = "thermo1.dat")
            write(myunit, '(A1,A12,3(1X,A13))') &
                "#", &
                "T/K", &
                "U/kJ", &
                "H/kJ", &
                "S/J*K^-1"
            do itemp = 1, ntemp
                if (converged(itemp)) then
                    write(myunit, '(ES13.6,3(1X,ES13.6))') &
                        temp(itemp), &
                        1.0e-3_dp*u(itemp), &
                        1.0e-3_dp*h(itemp), &
                        s(itemp)
                end if
            end do
            close(myunit)
        end subroutine write_thermo1
        !=================================================================================
        ! Writes thermodynamic functions that depend on second derivatives.
        subroutine write_thermo2(ntemp, converged, temp, cv, cp)
            integer, intent(in) :: ntemp
            logical, dimension(ntemp), intent(in) :: converged
            real(dp), dimension(ntemp), intent(in) :: temp, cv, cp

            integer :: myunit
            integer :: itemp

            open(newunit = myunit, action = "write", status = "unknown", &
                file = "thermo2.dat")
            write(myunit, '(A1,A12,2(1X,A13))') &
                "#", &
                "T/K", &
                "Cv/J*K^-1", &
                "Cp/J*K^-1"
            do itemp = 1, ntemp
                if (converged(itemp)) then
                    write(myunit, '(ES13.6,2(1X,ES13.6))') &
                        temp(itemp), &
                        cv(itemp), &
                        cp(itemp)
                end if
            end do
            close(myunit)
        end subroutine write_thermo2
        !=================================================================================
        ! Writes partition functions.
        subroutine write_partition_functions(ntemp, converged, temp, lnq, dlnq, ddlnq)
            integer, intent(in) :: ntemp
            logical, dimension(ntemp), intent(in) :: converged
            real(dp), dimension(ntemp), intent(in) :: temp, lnq, dlnq, ddlnq

            integer :: myunit
            integer :: itemp

            open(newunit = myunit, action = "write", status = "unknown", &
                file = "partition_functions.dat")
            write(myunit, '(A1,A12,3(1X,A13))') &
                "#", &
                "T/K", &
                "lnQ", &
                "dlnQ", &
                "ddlnQ"
            do itemp = 1, ntemp
                if (converged(itemp)) then
                    write(myunit, '(ES13.6,3(1X,ES13.6))') &
                        temp(itemp), &
                        lnq(itemp), &
                        dlnq(itemp), &
                        ddlnq(itemp)
                end if
            end do
            close(myunit)
        end subroutine write_partition_functions
        !=================================================================================
        ! Writes thermodynamic functions of clusters.
        subroutine write_cluster_thermo(ntemp, nclust, converged, labels, temp, E, header, &
                factor, filename)
            use iso_varying_string
            integer, intent(in) :: ntemp, nclust
            logical, dimension(ntemp), intent(in) :: converged
            type(varying_string), dimension(nclust), intent(in) :: labels
            real(dp), dimension(ntemp), intent(in) :: temp
            real(dp), dimension(ntemp, nclust), intent(in) :: E
            character(*), intent(in) :: header, filename
            real(dp), intent(in) :: factor

            integer :: myunit
            integer :: itemp, iclust
            character(80) :: fmtspec

            ! Write cluster specific thermodynamics
            write(fmtspec, '(A,G0,A)') '(A1,A12,', nclust, '(1X,A13))' 
            open(newunit = myunit, action = "write", status = "unknown", &
                file = filename)
            write(myunit, '(A)') "# T/K; " // header
            write(myunit, fmtspec) "#", "T/K", (char(labels(iclust)), iclust = 1, nclust)
            write(fmtspec, '(A,G0,A)') '(ES13.6,', nclust, '(1X,ES13.6))' 
            do itemp = 1, ntemp
                if (converged(itemp)) write(myunit, fmtspec) temp(itemp), &
                    (E(itemp, iclust)*factor, iclust = 1, nclust)
            end do
            close(myunit)
        end subroutine write_cluster_thermo
        !=================================================================================
        ! Writes contributions to the a thermodynamic function.
        subroutine write_contributions(ntemp, converged, temp, x, vib, rot, trans, elec, &
                mf, indi, header, factor, filename)
            integer, intent(in) :: ntemp
            logical, dimension(ntemp), intent(in) :: converged
            real(dp), dimension(ntemp), intent(in) :: temp, x, vib, rot, trans, elec, &
                mf, indi
            character(*), intent(in) :: header, filename
            real(dp), intent(in) :: factor

            integer :: myunit
            integer :: itemp

            open(newunit = myunit, action = "write", status = "unknown", &
                file = filename)
            write(myunit, '(A)') "# T/K; " // header
            write(myunit, '(A1,A12,7(1X,A13))') &
                "#", &
                "T", &
                "tot", &
                "vib", &
                "rot", &
                "trans", &
                "elec", &
                "mf", &
                "indi"
            do itemp = 1, ntemp
                if (converged(itemp)) then
                    write(myunit, '(ES13.6,7(1X,ES13.6))') &
                        temp(itemp), &
                        factor*x(itemp), &
                        factor*vib(itemp), &
                        factor*rot(itemp), &
                        factor*trans(itemp), &
                        factor*elec(itemp), &
                        factor*mf(itemp), &
                        factor*indi(itemp)
                end if
            end do
            close(myunit)
        end subroutine write_contributions
        !=================================================================================
        ! Calculates and writes cluster specific thermodynamic functions.
        subroutine calculate_cluster_thermo(ntemp, nclust, pop, press, temp, vol, &
                   lnq_clust, dlnq_clust, ddlnq_clust, dvol, converged, labels)
            use partition_functions, only: pf_t
            use input, only: pmk_input
            use iso_varying_string
            use constants
            use cluster
            use auxiliary, only: ln_factorial
            integer, intent(in) :: ntemp, nclust
            real(dp), dimension(ntemp, nclust), intent(in) :: pop
            real(dp), intent(in) :: press
            real(dp), dimension(ntemp), intent(in) :: temp, vol, dvol
            type(pf_t), dimension(ntemp, nclust), intent(in) :: lnq_clust, dlnq_clust, &
                ddlnq_clust
            logical, dimension(ntemp), intent(in) :: converged
            type(varying_string), dimension(nclust), intent(in) :: labels

            integer :: iclust, itemp
            real(dp), dimension(ntemp, nclust) :: a, g, u, s, h, cv, cp, indi

            indi = 0.0_dp
            do itemp = 1, ntemp
                do iclust = 1, nclust
!                    write(*,*) pop(itemp, iclust)
                    indi(itemp, iclust) = temp(itemp) * kb * ln_factorial(pop(itemp, iclust))
                end do
            end do

            ! Calculate and write thermodynamic functions of clusters.
            do iclust = 1, nclust
              call calculate_helmholtz_energy(ntemp, temp, lnq_clust(:,iclust)%qtot*avogadro, a(:,iclust))
              call calculate_gibbs_enthalpy(ntemp, temp, lnq_clust(:,iclust)%qtot*avogadro, &
                     vol, press, g(:,iclust))
              call calculate_internal_energy(ntemp, temp, dlnq_clust(:,iclust)%qtot*avogadro, u(:,iclust))
              call calculate_enthalpy(ntemp, temp, dlnq_clust(:,iclust)%qtot*avogadro, vol, press, h(:,iclust))
              call calculate_entropy(ntemp, temp, lnq_clust(:,iclust)%qtot*avogadro, &
                     dlnq_clust(:,iclust)%qtot*avogadro, s(:,iclust))
              call calculate_cv(ntemp, temp, dlnq_clust(:,iclust)%qtot, ddlnq_clust(:,iclust)%qtot, cv(:,iclust))
              call calculate_cp(ntemp, temp, dlnq_clust(:,iclust)%qtot, ddlnq_clust(:,iclust)%qtot, dvol, press, cp(:,iclust))
            end do

            call write_cluster_thermo(ntemp, nclust, converged, labels, temp, a, "A/kJ", 1.0e-3_dp, "helmholtz_clusters.dat")
            call write_cluster_thermo(ntemp, nclust, converged, labels, temp, g, "G/kJ", 1.0e-3_dp, "gibbs_clusters.dat")
            call write_cluster_thermo(ntemp, nclust, converged, labels, temp, u, "U/kJ", 1.0e-3_dp, "internal_clusters.dat")
            call write_cluster_thermo(ntemp, nclust, converged, labels, temp, h, "H/kJ", 1.0e-3_dp, "enthalpy_clusters.dat")
            call write_cluster_thermo(ntemp, nclust, converged, labels, temp, s, "S/J*K^-1", 1.0_dp, "entropy_clusters.dat")
            call write_cluster_thermo(ntemp, nclust, converged, labels, temp, indi, "I/kJ", 1.0e-3_dp, "indi_clusters.dat")
            call write_cluster_thermo(ntemp, nclust, converged, labels, temp, cv, "Cv/J*K^-1", 1.0_dp, "cv_clusters.dat")
            call write_cluster_thermo(ntemp, nclust, converged, labels, temp, cp, "Cp/J*K^-1", 1.0_dp, "cp_clusters.dat")
            call write_cluster_thermo(ntemp, nclust, converged, labels, temp, lnq_clust%qtot, &
                  "ln(q)", 1.0_dp, "lnq_clusters.dat")
            call write_cluster_thermo(ntemp, nclust, converged, labels, temp, dlnq_clust%qtot, &
                  "d ln(q) / dT", 1.0_dp, "dlnq_clusters.dat")
            call write_cluster_thermo(ntemp, nclust, converged, labels, temp, ddlnq_clust%qtot, & 
                  "d^2 ln(q) / dT^2", 1.0_dp, "ddlnq_clusters.dat")
        end subroutine calculate_cluster_thermo
        !=================================================================================
        ! Calculates and writes thermodynamic functions.
        subroutine calculate_thermo(ntemp, nclust, pop, press, temp, vol, vexcl, &
                lnq_clust, dlnq_clust, ddlnq_clust, dvol, solution, converged)
            use partition_functions, only: pf_t
            use input, only: pmk_input
            integer, intent(in) :: ntemp, nclust
            real(dp), dimension(ntemp, nclust), intent(in) :: pop
            real(dp), intent(in) :: press, vexcl
            real(dp), dimension(ntemp), intent(in) :: temp, vol, dvol
            type(pf_t), dimension(ntemp, nclust), intent(in) :: lnq_clust, dlnq_clust, &
                ddlnq_clust
            integer, dimension(ntemp), intent(in) :: solution
            logical, dimension(ntemp), intent(in) :: converged

            real(dp), dimension(ntemp) :: lnq, dlnq, ddlnq
            real(dp), dimension(ntemp) :: lnq_vib, dlnq_vib, ddlnq_vib
            real(dp), dimension(ntemp) :: lnq_rot, dlnq_rot, ddlnq_rot
            real(dp), dimension(ntemp) :: lnq_elec, dlnq_elec, ddlnq_elec
            real(dp), dimension(ntemp) :: lnq_int, dlnq_int, ddlnq_int
            real(dp), dimension(ntemp) :: lnq_trans, dlnq_trans, ddlnq_trans
            real(dp), dimension(ntemp) :: lnq_indi, dlnq_indi, ddlnq_indi
            real(dp), dimension(ntemp) :: vib, rot, trans, elec, mf, indi
            real(dp), dimension(ntemp) :: a, g, u, s, h, alpha, cv, cp

            ! Calculate system partition function.
            call calculate_lnq_sys(ntemp, nclust, pop, lnq_clust(:,:)%qtot, lnq)
            call add_lnq_indi(ntemp, nclust, pop, lnq)
            call calculate_lnq_sys(ntemp, nclust, pop, dlnq_clust(:,:)%qtot, dlnq)
            call calculate_lnq_sys(ntemp, nclust, pop, ddlnq_clust(:,:)%qtot, ddlnq)

            ! Calculate contributions to the system partition function.
            call calculate_lnq_sys(ntemp, nclust, pop, lnq_clust(:,:)%qvib, lnq_vib)
            call calculate_lnq_sys(ntemp, nclust, pop, dlnq_clust(:,:)%qvib, dlnq_vib)
            call calculate_lnq_sys(ntemp, nclust, pop, ddlnq_clust(:,:)%qvib, ddlnq_vib)
            call calculate_lnq_sys(ntemp, nclust, pop, lnq_clust(:,:)%qrot, lnq_rot)
            call calculate_lnq_sys(ntemp, nclust, pop, dlnq_clust(:,:)%qrot, dlnq_rot)
            call calculate_lnq_sys(ntemp, nclust, pop, ddlnq_clust(:,:)%qrot, ddlnq_rot)
            call calculate_lnq_sys(ntemp, nclust, pop, lnq_clust(:,:)%qelec, lnq_elec)
            call calculate_lnq_sys(ntemp, nclust, pop, dlnq_clust(:,:)%qelec, dlnq_elec)
            call calculate_lnq_sys(ntemp, nclust, pop, ddlnq_clust(:,:)%qelec, ddlnq_elec)
            call calculate_lnq_sys(ntemp, nclust, pop, lnq_clust(:,:)%qint, lnq_int)
            call calculate_lnq_sys(ntemp, nclust, pop, dlnq_clust(:,:)%qint, dlnq_int)
            call calculate_lnq_sys(ntemp, nclust, pop, ddlnq_clust(:,:)%qint, ddlnq_int)
            call calculate_lnq_sys(ntemp, nclust, pop, lnq_clust(:,:)%qtrans, lnq_trans)
            call calculate_lnq_sys(ntemp, nclust, pop, dlnq_clust(:,:)%qtrans, dlnq_trans)
            call calculate_lnq_sys(ntemp, nclust, pop, ddlnq_clust(:,:)%qtrans, &
                ddlnq_trans)
            lnq_indi = 0.0_dp
            call add_lnq_indi(ntemp, nclust, pop, lnq_indi)
            dlnq_indi = 0.0_dp
            ddlnq_indi = 0.0_dp

            ! Write volume, exclusion volume, and status code.
            call calculate_expansion_coefficient(ntemp, vol, dvol, alpha)
            call write_volume(ntemp, converged, temp, vol, vexcl, alpha, solution)

            ! Calculate and write thermodynamic functions.
            call calculate_helmholtz_energy(ntemp, temp, lnq, a)
            call calculate_gibbs_enthalpy(ntemp, temp, lnq, vol, press, g)
            call calculate_internal_energy(ntemp, temp, dlnq, u)
            call calculate_enthalpy(ntemp, temp, dlnq, vol, press, h)
            call calculate_entropy(ntemp, temp, lnq, dlnq, s)
            call calculate_cv(ntemp, temp, dlnq, ddlnq, cv)
            call calculate_cp(ntemp, temp, dlnq, ddlnq, dvol, press, cp)
            call write_thermo0(ntemp, converged, temp, a, g)
            call write_thermo1(ntemp, converged, temp, u, h, s)
            call write_thermo2(ntemp, converged, temp, cv, cp)

            ! Calculate and write contributions to the Helmholtz free energy.
            if (pmk_input%contrib .or. pmk_input%helmholtz_contrib) then
                call calculate_helmholtz_energy(ntemp, temp, lnq_vib, vib)
                call calculate_helmholtz_energy(ntemp, temp, lnq_rot, rot)
                call calculate_helmholtz_energy(ntemp, temp, lnq_trans, trans)
                call calculate_helmholtz_energy(ntemp, temp, lnq_elec, elec)
                call calculate_helmholtz_energy(ntemp, temp, lnq_int, mf)
                call calculate_helmholtz_energy(ntemp, temp, lnq_indi, indi)
                call write_contributions(ntemp, converged, temp, a, vib, rot, trans, &
                    elec, mf, indi, "A/kJ", 1.0e-3_dp, "helmholtz_contrib.dat")
            end if

            ! Calculate and write contributions to the internal energy.
            if (pmk_input%contrib .or. pmk_input%internal_contrib) then
                call calculate_internal_energy(ntemp, temp, dlnq_vib, vib)
                call calculate_internal_energy(ntemp, temp, dlnq_rot, rot)
                call calculate_internal_energy(ntemp, temp, dlnq_trans, trans)
                call calculate_internal_energy(ntemp, temp, dlnq_elec, elec)
                call calculate_internal_energy(ntemp, temp, dlnq_int, mf)
                call calculate_internal_energy(ntemp, temp, dlnq_indi, indi)
                call write_contributions(ntemp, converged, temp, u, vib, rot, trans, &
                    elec, mf, indi, "U/kJ", 1.0e-3_dp, "internal_contrib.dat")
            end if

            ! Calculate and write contributions to the entropy.
            if (pmk_input%contrib .or. pmk_input%entropy_contrib) then
                call calculate_entropy(ntemp, temp, lnq_vib, dlnq_vib, vib)
                call calculate_entropy(ntemp, temp, lnq_rot, dlnq_rot, rot)
                call calculate_entropy(ntemp, temp, lnq_trans, dlnq_trans, trans)
                ! Electronic contributions to entropy are zero.
                elec = 0.0_dp
                ! Mean field contributions to entropy are zero.
                mf = 0.0_dp
                call calculate_entropy(ntemp, temp, lnq_indi, dlnq_indi, indi)
                call write_contributions(ntemp, converged, temp, s, vib, rot, trans, &
                    elec, mf, indi, "S/J*K^-1", 1.0_dp, "entropy_contrib.dat")
            end if

            ! Calculate and write contributions to the heat capacity at constant volume.
            if (pmk_input%contrib .or. pmk_input%cv_contrib) then
                call calculate_cv(ntemp, temp, dlnq_vib, ddlnq_vib, vib)
                call calculate_cv(ntemp, temp, dlnq_rot, ddlnq_rot, rot)
                call calculate_cv(ntemp, temp, dlnq_trans, ddlnq_trans, trans)
                ! Electronic contributions to Cv are zero.
                elec = 0.0_dp
                ! Mean field contributions to Cv are zero.
                mf = 0.0_dp
                call calculate_cv(ntemp, temp, dlnq_indi, ddlnq_indi, indi)
                call write_contributions(ntemp, converged, temp, cv, vib, rot, trans, &
                    elec, mf, indi, "Cv/J*K^-1", 1.0_dp, "cv_contrib.dat")
            end if

            ! Calculate and write the total partition function and its derivatives.
            if (pmk_input%contrib) then
                call write_partition_functions(ntemp, converged, temp, lnq, dlnq, ddlnq)
            end if
        end subroutine calculate_thermo
        !=================================================================================
end module thermo
!=========================================================================================
