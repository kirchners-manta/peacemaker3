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
! This module provides the cluster_t data type, which represents a cluster, and associated
! procedures. The cluster_t data type contains static information about the clusters,
! that are shared across all QCE calculations. By a QCE calculation we mean a complete
! cycle of iterations for a given thermodynamic state (N, p, T) and parameter set 
! (amf, bxv). Peacemaker usually performs many of these calculations, as it samples
! temperatures, parameters, etc. Thus quantities such as moments of inertia and
! frequencies go in here, but temperature or pressure don't, as these may differ across
! multiple QCE calculations.
module cluster
    use kinds
    use iso_varying_string
    use config
    use error
    use input
    use auxiliary
    use constants
    implicit none
    private
    !=====================================================================================
    ! Public entities.
    public :: clusterset
    public :: monomer
    public :: setup_clusterset
    public :: print_clusterset
    public :: check_clusterset
    !=====================================================================================
    ! The cluster_t data type, which represents a cluster. Reasonable defaults should be
    ! set here, if there are any.
    type :: cluster_t
        ! Total mass of the cluster in amu.
        real(dp) :: mass
        ! Rotational symmetry number sigma.
        integer :: sigma = 1
        ! Vibrational frequencies in 1/cm.
        real(dp), dimension(:), allocatable :: frequencies
        ! Principal moments of inertia in amu*Angstrom^2.
        real(dp), dimension(:), allocatable :: inertia
        ! Composition array (number of monomers).
        integer, dimension(:), allocatable :: composition
        ! Adiabatic interaction energy in kJ/mol.
        real(dp) :: energy
        ! Cluster volume in Angstrom^3
        real(dp) :: volume
        ! Frequency scaling factor.
        real(dp) :: fscale = 1.0_dp
        ! Anharmonicity constant
        real(dp) :: anharmonicity = 0.0_dp
        ! Flags classifying the cluster.
        logical :: linear = .false.
        logical :: atom = .false.
        logical :: monomer = .false.
        ! Cluster label.
        type(varying_string) :: label
    end type cluster_t
    !=====================================================================================
    ! The clusterset array and monomer array. The clusterset array is the central quantity
    ! in Peacemaker. It contains all the static information about all the clusters. The
    ! monomer array contains the indices of monomers in the cluster set array.
    type(cluster_t), dimension(:), allocatable, target :: clusterset
    integer, dimension(:), allocatable :: monomer
    !=====================================================================================
    contains
        !=================================================================================
        ! Given the clusterset configuration, this procedure sets up the clusterset array.
        subroutine setup_clusterset(cfg)
            ! The clusterset configuration.
            type(config_t), intent(inout) :: cfg
    
            ! A pointer to the current record.
            type(record_t), pointer :: p
            ! A pointer to the current cluster.
            type(cluster_t), pointer :: c
            ! The number of clusters
            integer:: nr_clusters
            ! Loop and status stuff.
            integer:: i
            integer:: j
            
            nr_clusters = cfg%nr_sections()
            if (nr_clusters == 0) call pmk_error("empty clusterset")
            allocate(clusterset(nr_clusters))
            do i = 1, nr_clusters
                ! Let c point to the current cluster for convenience.
                c => clusterset(i)
                ! Get the cluster label.
                c%label = cfg%get_section_label(i)
    
                ! Check for the monomer flag.
                p => cfg%get_record(c%label, "monomer")
                if (associated(p)) then
                    call process_monomer_record(c, p%nr_args)
                    ! Get the monomer volume.
                    p => cfg%get_record(c%label, "volume")
                    if (associated(p)) then
                        call process_volume_record(c, p%nr_args, p%args)
                    else
                        call pmk_missing_key_error("volume", c%label)
                    end if
                end if
    
                ! Get moments of inertia and mass from xyz file.
                p => cfg%get_record(c%label, "coordinates")
                if (associated(p)) then
                    call process_coordinates_record(c, p%nr_args, p%args)
                else
                    call pmk_missing_key_error("coordinates", c%label)
                end if
    
                ! Get the cluster composition.
                p => cfg%get_record(c%label, "composition")
                if (associated(p)) then
                    call process_composition_record(c, p%nr_args, p%args)
                else
                    call pmk_missing_key_error("composition", c%label)
                end if
    
                ! Get the frequency scaling factor.
                p => cfg%get_record(c%label, "frequency_scale")
                if (associated(p)) then
                    call process_frequency_scale_record(c, p%nr_args, p%args)
                end if
    
                ! Get the frequencies.
                p => cfg%get_record(c%label, "frequencies")
                if (associated(p)) then
                    call process_frequencies_record(c, p%nr_args, p%args)
                else
                    call pmk_missing_key_error("frequencies", c%label)
                end if
    
                ! Get the anharmonicity constant.
                p => cfg%get_record(c%label, "anharmonicity")
                if (associated(p)) then
                    call process_anharmonicity_record(c, p%nr_args, p%args)
                end if
    
                ! Get the interaction energy.
                p => cfg%get_record(c%label, "energy")
                if (associated(p)) then
                    call process_energy_record(c, p%nr_args, p%args)
                else
                    call pmk_missing_key_error("energy", c%label)
                end if
    
                ! Get the rotational symmetry number sigma.
                p => cfg%get_record(c%label, "sigma")
                if (associated(p)) then
                    call process_sigma_record(c, p%nr_args, p%args)
                else
                    call pmk_missing_key_error("sigma", c%label)
                end if
    
            end do
    
            ! Set up the monomer array, assuming that everything is alright. We will check
            ! for errors later.
            allocate(monomer(pmk_input%components))
            monomer = 0
            do i = 1, nr_clusters
                if (clusterset(i)%monomer) then
                    do j = 1, pmk_input%components
                        if (clusterset(i)%composition(j) == 1) monomer(j) = i
                    end do
                end if
            end do
    
            ! Calculate cluster volumes, assuming that everything is alright. We will
            ! check for errors later.
            do i = 1, nr_clusters
                if (clusterset(i)%monomer) cycle
                clusterset(i)%volume = 0.0_dp
                do j = 1, pmk_input%components
                    clusterset(i)%volume = clusterset(i)%volume + &
                        real(clusterset(i)%composition(j), dp) * &
                        clusterset(monomer(j))%volume
                end do
            end do
        end subroutine setup_clusterset
        !=================================================================================
        ! Prints information about the processed cluster set to the screen.
        subroutine print_clusterset()
            ! A pointer to the current cluster for convenience.
            type(cluster_t), pointer :: c
            integer:: i
    
            write(*,'(4X,A)') "Using the following clusterset:"
            write(*,*)
            do i = 1, size(clusterset)
                c => clusterset(i)
    
                ! Print cluster label.
                if (c%monomer) then
                    write(*,'(8X,3A)') "[", char(c%label), "] (monomer)"
                else
                    write(*,'(8X,3A)') "[", char(c%label), "]"
                end if
    
                ! Print composition.
                write(*, '(12X,A,1X)', advance = "no") "composition:"
                call array_sample(c%composition)
                write(*, *)
    
                ! Print rotational symmetry number, energy, volume, mass.
                write(*, '(12X,A,1X,G0)') "sigma:", c%sigma
                write(*, '(12X,A,1X,G0.6,1X,A)') "energy:", c%energy, "[kJ/mol]"
                write(*, '(12X,A,1X,G0.6,1X,A)') "volume:", c%volume, "[A^3]"
                write(*, '(12X,A,1X,G0.6,1X,A)') "mass:", c%mass, "[amu]"
    
                ! Print intertia array.
                write(*, '(12X,A,1X,3(G0.6,1X))', advance = "no") "inertia:", &
                    c%inertia(:)
                write(*, '(A)') "[amu*Angstrom^2]"
    
                ! Print anhamronicity constant.
                if (c%anharmonicity > 0.0_dp) then
                    write(*, '(12X,A,1X,G0.6)') "anharmonicity constant", c%anharmonicity
                end if
    
                ! Print frequencies and scaling factor.
                if (abs(c%fscale-1.0_dp) > global_eps) then
                    write(*, '(12X,A,G0.6,A)', advance = "no") &
                        "frequencies (scaled by ", c%fscale, "): "
                else
                    write(*, '(12X,A,1X)', advance = "no") "frequencies:"
                end if
                call array_sample(c%frequencies)
                write(*, '(1X,A)') "[1/cm]"
                write(*, *)
            end do
        end subroutine print_clusterset
        !=================================================================================
        ! Processes the monomer record.
        subroutine process_monomer_record(c, nr_args)
            type(cluster_t), pointer, intent(inout) :: c
            integer, intent(in) :: nr_args
    
            if (nr_args == 0) then
                c%monomer = .true.
            else
                call pmk_argument_count_error("monomer", c%label)
            end if
        end subroutine process_monomer_record
        !=================================================================================
        ! Processes a coordinates record.
        subroutine process_coordinates_record(c, nr_args, args)
            use atomic_data, only: periodic_table
            use constants
            type(cluster_t), pointer, intent(inout) :: c
            integer, intent(in) :: nr_args
            type(varying_string), dimension(nr_args), intent(in) :: args
    
            integer:: i
            integer:: ios
            integer:: my_unit
            integer:: nr_atoms
    
            character(2), dimension(:), allocatable :: label
            real(dp), dimension(:), allocatable :: mass
            real(dp), dimension(:, :), allocatable :: xyz
            real(dp), dimension(3) :: com
            real(dp), dimension(3, 3) :: inertia
            real(dp), dimension(3, 3) :: identity
            real(dp), dimension(3, 3) :: B
            real(dp):: p
            real(dp):: p1
            real(dp):: p2
            real(dp):: q
            real(dp):: r
            real(dp):: phi
            real(dp):: tmp
            real(dp), dimension(3) :: eig
            integer:: n
    
            if (nr_args == 1) then
                ! Open unit.
                open(newunit = my_unit, file = char(args(1)), action = 'read', &
                    status = 'old', iostat = ios)
                if (ios /= 0) call pmk_error("could not open '" // char(args(1)) // "'")
    
                ! Read number of atoms.
                read(my_unit, *, iostat = ios) nr_atoms
                if (ios /= 0) call pmk_error("illegal file format in '" // &
                    char(args(1)) // "'")
    
                ! Comment line.
                read(my_unit, *, iostat = ios)
    
                ! Read coordinates.
                allocate(label(nr_atoms))
                allocate(mass(nr_atoms))
                allocate(xyz(nr_atoms, 3))
                do i = 1, nr_atoms
                    read(my_unit, *, iostat = ios) &
                        label(i), xyz(i, 1), xyz(i, 2), xyz(i, 3)
                    if (ios /= 0) call pmk_error("unexpected end of file in '" // &
                        char(args(1)) // "'")
                end do
    
                ! Assign masses and calculate total mass.
                do i = 1, nr_atoms
                    mass(i) = periodic_table%mass(label(i))
                end do
                c%mass = sum(mass)
    
                ! Calculate center of mass and shift to origin.
                com = 0.0_dp
                do i = 1, nr_atoms
                    com(:) = com(:) + mass(i)*xyz(i, :)
                end do
                com(:) = com(:) / c%mass
                do i = 1, nr_atoms
                    xyz(i, :) = xyz(i, :) - com(:)
                end do
    
                ! Calculate inertia tensor.
                inertia = 0.0_dp
                do i = 1, nr_atoms
                    inertia(1,1) = inertia(1,1) + mass(i)*(xyz(i,2)**2 + xyz(i,3)**2)
                    inertia(2,2) = inertia(2,2) + mass(i)*(xyz(i,1)**2 + xyz(i,3)**2)
                    inertia(3,3) = inertia(3,3) + mass(i)*(xyz(i,1)**2 + xyz(i,2)**2)
                    inertia(1,2) = inertia(1,2) - mass(i)*xyz(i,1)*xyz(i,2)
                    inertia(1,3) = inertia(1,3) - mass(i)*xyz(i,1)*xyz(i,3)
                    inertia(2,3) = inertia(2,3) - mass(i)*xyz(i,2)*xyz(i,3)
                end do
                inertia(2, 1) = inertia(1, 2)
                inertia(3, 1) = inertia(1, 3)
                inertia(3, 2) = inertia(2, 3)
    
                ! Diagonalize inertia tensor.
                identity = 0.0_dp
                identity(1,1) = 1.0_dp
                identity(2,2) = 1.0_dp
                identity(3,3) = 1.0_dp
                p1 = inertia(1,2)**2 + inertia(1,3)**2 + inertia(2,3)**2
                if (p1 <= global_eps) then
                    ! Tensor is already diagonal.
                    eig(1) = inertia(1,1)
                    eig(2) = inertia(2,2)
                    eig(3) = inertia(3,3)
                else
                    q = (inertia(1,1) + inertia(2,2) + inertia(3,3))/3.0_dp
                    p2 = (inertia(1,1)-q)**2 + (inertia(2,2)-q)**2 + &
                        (inertia(3,3)-q)**2 + 2.0_dp*p1
                    p = sqrt(p2/6.0_dp)
                    B = (inertia(:, :) - q*identity(:, :))/p
                    r = B(1,1)*B(2,2)*B(3,3) + B(1,2)*B(2,3)*B(3,1) + &
                        B(1,3)*B(2,1)*B(3,2) - B(1,3)*B(2,2)*B(3,1) - &
                        B(1,2)*B(2,1)*B(3,3) - B(1,1)*B(2,3)*B(3,2)
                    r = 0.5_dp * r
    
                    if (r <= -1.0_dp) then
                        phi = pi/3.0_dp
                    else if (r >= 1.0_dp) then
                        phi = 0.0_dp
                    else
                        phi = acos(r)/3.0_dp
                    end if
    
                    eig(1) = q + 2.0_dp*p*cos(phi)
                    eig(3) = q + 2.0_dp*p*cos(phi + (2.0_dp*pi/3.0_dp))
                    eig(2) = 3.0_dp*q - eig(1) - eig(3)
                end if
    
                ! Sort the eigenvalues.
                if (eig(1) < eig(2)) then
                    tmp = eig(2)
                    eig(2) = eig(1)
                    eig(1) = tmp
                end if
                if (eig(1) < eig(3)) then
                    tmp = eig(3)
                    eig(3) = eig(1)
                    eig(1) = tmp
                end if
                if (eig(2) < eig(3)) then
                    tmp = eig(3)
                    eig(3) = eig(2)
                    eig(2) = tmp
                end if
    
                ! Check for atoms, linear molecules, and assign moments of inertia.
                n = 3 - count(eig <= global_eps)
                if (n == 3) then
                    allocate(c%inertia(3))
                    c%inertia = eig
                else if (n == 2) then
                    c%linear = .true.
                    allocate(c%inertia(1))
                    c%inertia(1) = eig(1)
                else if (n == 0) then
                    c%atom = .true.
                end if
    
                ! Clean up.
                close(my_unit)
                deallocate(label)
                deallocate(mass)
                deallocate(xyz)
            else
                call pmk_argument_count_error("coordinates", c%label)
            end if
        end subroutine process_coordinates_record
        !=================================================================================
        ! Processes a composition record.
        subroutine process_composition_record(c, nr_args, args)
            type(cluster_t), pointer, intent(inout) :: c
            integer, intent(in) :: nr_args
            type(varying_string), dimension(nr_args), intent(in) :: args
    
            integer:: i
            integer:: ios
    
            if (nr_args == pmk_input%components) then
                allocate(c%composition(pmk_input%components))
                do i = 1, pmk_input%components
                    c%composition(i) = string2int(args(i), ios)
                    if (ios /= 0) call pmk_argument_error("composition", c%label)
                end do
            else
                call pmk_argument_count_error("composition", c%label)
            end if
        end subroutine process_composition_record
        !=================================================================================
        ! Reads frequencies from a frequency file. A frequency file is similar to an
        ! xyz file. It contains the number of vibrational frequencies in the first line,
        ! followed by an empty/comment line, followed by one line per frequency.
        ! Sorts the frequencies.
        subroutine process_frequencies_record(c, nr_args, args)
            type(cluster_t), pointer, intent(inout) :: c
            integer, intent(in) :: nr_args
            type(varying_string), dimension(nr_args), intent(in) :: args
    
            integer:: nr_frequencies
            integer:: nr_zeros
            integer:: my_unit
            integer:: ios
            integer:: i
            integer:: j
            real(dp):: tmp
    
            if (nr_args == 1) then
                ! Open unit.
                open(newunit = my_unit, file = char(args(1)), action = 'read', &
                    status = 'old', iostat = ios)
                if (ios /= 0) call pmk_error("could not open '" // char(args(1)) // "'")
    
                ! Read number of frequencies.
                read(my_unit, *, iostat = ios) nr_frequencies
                if (ios /= 0) call pmk_error("illegal file format in '" // &
                    char(args(1)) // "'")
    
                ! Comment line.
                read(my_unit, *, iostat = ios)
    
                ! Read number of zeros.
                nr_zeros = 0
                do i = 1, nr_frequencies
                    read(my_unit, *, iostat = ios) tmp
                    if (ios /= 0) call pmk_error("unexpected end of file in '" // &
                        char(args(1)) // "'")
                    if (abs(tmp) <= global_eps) then
                        nr_zeros = nr_zeros + 1
                    end if
                end do
                rewind(my_unit)
                allocate(c%frequencies(nr_frequencies - nr_zeros))
    
                ! Read frequencies.
                nr_zeros = 0
                read(my_unit, *, iostat = ios)
                read(my_unit, *, iostat = ios)
                do i = 1, nr_frequencies
                    read(my_unit, *, iostat = ios) tmp
                    if (abs(tmp) <= global_eps) then
                        nr_zeros = nr_zeros + 1
                    else
                        c%frequencies(i - nr_zeros) = tmp
                    end if
                end do
    
                ! Sort frequencies.
                do i = 1, size(c%frequencies) - 1
                    do j = i + 1, size(c%frequencies)
                        if (c%frequencies(j) < c%frequencies(i)) then
                            tmp = c%frequencies(i)
                            c%frequencies(i) = c%frequencies(j)
                            c%frequencies(j) = tmp
                        end if
                    end do
                end do
    
                c%frequencies = c%fscale*c%frequencies
                close(my_unit)
            else
                call pmk_argument_count_error("frequencies", c%label)
            end if
        end subroutine process_frequencies_record
        !=================================================================================
        ! Processes the adiabatic interaction energy.
        subroutine process_energy_record(c, nr_args, args)
            type(cluster_t), pointer, intent(inout) :: c
            integer, intent(in) :: nr_args
            type(varying_string), dimension(nr_args), intent(in) :: args
    
            integer:: ios
    
            if (nr_args == 1) then
                c%energy = string2real(args(1), ios)
                if (ios /= 0) call pmk_argument_error("energy", c%label)
            else
                call pmk_argument_count_error("energy", c%label)
            end if
        end subroutine process_energy_record
        !=================================================================================
        ! Processes the rotational symmetry number.
        subroutine process_sigma_record(c, nr_args, args)
            type(cluster_t), pointer, intent(inout) :: c
            integer, intent(in) :: nr_args
            type(varying_string), dimension(nr_args), intent(in) :: args
    
            integer:: ios
    
            if (nr_args == 1) then
                c%sigma = string2int(args(1), ios)
                if (ios /= 0) call pmk_argument_error("sigma", c%label)
            else
                call pmk_argument_count_error("sigma", c%label)
            end if
        end subroutine process_sigma_record
        !=================================================================================
        ! Processes the volume record.
        subroutine process_volume_record(c, nr_args, args)
            type(cluster_t), pointer, intent(inout) :: c
            integer, intent(in) :: nr_args
            type(varying_string), dimension(nr_args), intent(in) :: args
    
            integer:: ios
    
            if (nr_args == 1) then
                c%volume = string2real(args(1), ios)
                if (ios /= 0) call pmk_argument_error("volume", c%label)
            else
                call pmk_argument_count_error("volume", c%label)
            end if
        end subroutine process_volume_record
        !=================================================================================
        ! Processes the frequency scaling factor record.
        subroutine process_frequency_scale_record(c, nr_args, args)
            type(cluster_t), pointer, intent(inout) :: c
            integer, intent(in) :: nr_args
            type(varying_string), dimension(nr_args), intent(in) :: args
    
            integer:: ios
    
            if (nr_args == 1) then
                c%fscale = string2real(args(1), ios)
                if (ios /= 0) call pmk_argument_error("frequency_scale", c%label)
            else
                call pmk_argument_count_error("frequency_scale", c%label)
            end if
        end subroutine process_frequency_scale_record
        !=================================================================================
        ! Processes the anharmonicity record.
        subroutine process_anharmonicity_record(c, nr_args, args)
            type(cluster_t), pointer, intent(inout) :: c
            integer, intent(in) :: nr_args
            type(varying_string), dimension(nr_args), intent(in) :: args
    
            integer:: ios
    
            if (nr_args == 1) then
                c%anharmonicity = string2real(args(1), ios)
                if (ios /= 0) call pmk_argument_error("anharmonicity", c%label)
            else
                call pmk_argument_count_error("anharmonicity", c%label)
            end if
        end subroutine process_anharmonicity_record
        !=================================================================================
        ! Performs sanity checks on the clusterset.
        subroutine check_clusterset()
            ! A pointer to the current cluster.
            type(cluster_t), pointer :: c
            ! Loop and status stuff.
            integer:: i
            integer:: j
    
            ! Check monomer count
            if (count(clusterset%monomer) /= pmk_input%components) &
                call pmk_error("invalid number of monomers")
    
            do i = 1, size(clusterset)
                c => clusterset(i)
    
                ! Check compositions
                if (any(c%composition < 0)) &
                    call pmk_unphysical_argument_error("composition", c%label)
    
                ! Check sigma
                if (c%sigma < 0) &
                    call pmk_unphysical_argument_error("sigma", c%label)
    
                ! Check the anharmoncity constant.
                if (c%anharmonicity < 0.0_dp) &
                    call pmk_unphysical_argument_error("anharmonicity", c%label)
    
                ! Check the frequency scaling factor.
                if (c%fscale <= 0.0_dp) &
                    call pmk_unphysical_argument_error("frequency_scale", c%label)
    
                ! Check frequencies
                do j = 1, size(c%frequencies)
                    if (c%frequencies(j) < 0.0_dp) &
                        call pmk_unphysical_argument_error("frequencies", c%label)
                end do
    
                ! Check volume
                if (c%volume <= 0.0_dp) &
                    call pmk_unphysical_argument_error("volume", c%label)
    
                ! Check monomers.
                if (c%monomer) then
    
                    ! Check sum of composition.
                    if (sum(c%composition) /= 1) call pmk_error(&
                        "monomer/composition mismatch in cluster '[" // &
                        char(c%label) // "]'")
    
                    ! Check that this is the only monomer for the given component, by
                    ! comparing it to the monomer array. Assume that there are multiple
                    ! monomers for a particular component. Then because of the way we set
                    ! up the monomer array in setup_clusterset(), the monomer array will
                    ! point to the last 'monomer' for this component. We can use this fact
                    ! to check for multiply defined monomers.
                    do j = 1, pmk_input%components
                        ! Check for missing monomer first
                        if (monomer(j) == 0) call pmk_error("missing monomers")
                        if (c%composition(j) == 1) then
                            if (monomer(j) /= i) call pmk_error( &
                                "more than one monomer per component")
                        end if
                    end do
                end if
            end do
        end subroutine check_clusterset
        !=================================================================================
end module cluster
!=========================================================================================
