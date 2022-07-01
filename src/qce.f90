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
! This module implements the QCE algorithm, the core of Peacemaker.
module qce
    use omp_lib
    use kinds
    use cluster
    use constants
    use auxiliary, only: range_t
    use shared_data
    use partition_functions
    implicit none
    private
    !=====================================================================================
    ! Public entities
    public :: qce_prepare
    public :: qce_start
    public :: qce_finalize
    !=====================================================================================
    ! Data type storing reference_data.
    type :: reference_t
        logical :: compare, compare_isobar, compare_density, compare_phase_transition
        real(dp) :: density_weight, isobar_weight, phase_transition_weight
        real(dp) :: phase_transition, density, density_temperature
        real(dp), dimension(:), allocatable :: isobar_temperature, isobar_volume
    end type reference_t
    !=====================================================================================
    ! Instantiation of the data type above.
    type(reference_t):: reference
    !=====================================================================================
    ! Data type representing an isobar. Instantions of this data type are distributed
    ! over threads. Together with global_data, they fully describe a QCE calculation.
    type :: isobar_t
        real(dp) :: amf, bxv, amf_temp, bxv_temp, error
        integer, dimension(:), allocatable :: solution
        logical, dimension(:), allocatable :: converged
        real(dp), dimension(:), allocatable :: vol, temp, gibbs
        real(dp), dimension(:, :), allocatable :: populations
        type(pf_t), dimension(:, :), allocatable :: lnq
    end type isobar_t
    !=====================================================================================
    ! Instantiation of the data type specified above. This is the data that qce_main()
    ! works with.
    type(isobar_t):: ib
    type(isobar_t):: best_ib
    !=====================================================================================
    contains
        !=================================================================================
        ! Prepares the isobar array. Use of input data is only legit inside this
        ! subroutine. Unit conversion may be performed here.
        subroutine qce_prepare()
            ! Use of input data is only legit inside this subroutine.
            use input
    
            ! Assign global data.
            global_data%press = pmk_input%pressure*1.0e5_dp ! Conversion to Pa
            global_data%max_deviation = pmk_input%max_deviation
            global_data%vdamp = pmk_input%volume_damping_factor
            global_data%rotor_cutoff = pmk_input%rotor_cutoff
            global_data%qce_iterations = pmk_input%qce_iterations
            global_data%newton_iterations = pmk_input%newton_iterations
            global_data%grid_iterations = pmk_input%grid_iterations
            global_data%optimizer = pmk_input%optimizer
            global_data%imode = pmk_input%imode
            global_data%amf = pmk_input%amf
            global_data%bxv = pmk_input%bxv
            global_data%amf_temp = pmk_input%amf_temp
            global_data%bxv_temp = pmk_input%bxv_temp
            global_data%temp = pmk_input%temperature
            global_data%progress_bar = pmk_input%progress_bar
    
            allocate(global_data%monomer_amounts(size(monomer)))
            global_data%monomer_amounts = pmk_input%monomer_amounts
            call initialize_conserved_quantities()

            allocate(global_data%degree(size(monomer)))
            call initialize_degree()
    
            global_data%nconverged = 0
    
            ! Assign reference data.
            reference%compare = pmk_input%compare
    
            reference%compare_isobar = pmk_input%compare_isobar
            if (reference%compare_isobar) then
                reference%isobar_weight = pmk_input%ref_isobar_weight
                allocate(reference%isobar_temperature( &
                    size(pmk_input%ref_isobar_temperature)))
                reference%isobar_temperature = pmk_input%ref_isobar_temperature
                allocate(reference%isobar_volume( &
                    size(pmk_input%ref_isobar_volume)))
                reference%isobar_volume = 1.0e-3_dp*pmk_input%ref_isobar_volume
            end if
    
            reference%compare_density = pmk_input%compare_density
            if (reference%compare_density) then
                reference%density_weight = pmk_input%ref_density_weight
                reference%density = pmk_input%ref_density
                reference%density_temperature = pmk_input%ref_density_temperature
            end if
    
            reference%compare_phase_transition = pmk_input%compare_phase_transition
            if (reference%compare_phase_transition) then
                reference%phase_transition_weight = pmk_input%ref_phase_transition_weight
                reference%phase_transition = pmk_input%ref_phase_transition
            end if
        end subroutine qce_prepare
        !=================================================================================
        ! Initializes the conserved quantities (total number of particles and total mass
        ! of the system).
        subroutine initialize_conserved_quantities()
            integer:: i
    
            ! Particle number
            allocate(global_data%ntot(size(monomer)))
            global_data%ntot = global_data%monomer_amounts*avogadro
    
            ! Mass
            global_data%mtot = 0.0_dp
            do i = 1, size(monomer)
                global_data%mtot = global_data%mtot + &
                    global_data%monomer_amounts(i)*clusterset(monomer(i))%mass
            end do
            global_data%mtot = global_data%mtot / 1000.0_dp
    
            ! Excluded volume
            global_data%vexcl = 0.0_dp
            do i = 1, size(monomer)
                global_data%vexcl = global_data%vexcl + &
                    global_data%monomer_amounts(i)*clusterset(monomer(i))%volume
            end do
            global_data%vexcl = global_data%vexcl*avogadro*1.0e-30_dp
        end subroutine initialize_conserved_quantities
        !=================================================================================
        ! Initializes the degree of the population and mass polynomials.
        subroutine initialize_degree()
            integer:: i
            integer:: j
    
            global_data%degree = 0
            do i = 1, size(clusterset)
                do j = 1, size(monomer)
                    if (clusterset(i)%composition(j) > global_data%degree(j)) &
                        global_data%degree(j) = clusterset(i)%composition(j)
                end do
            end do
        end subroutine initialize_degree
        !=================================================================================
        ! Starts all QCE calculations. The main loop in here is OpenMP parallelized.
        subroutine qce_start()
            use auxiliary, only: progress_bar
            integer:: igrid
            integer:: iamf
            integer:: ibxv
            integer:: iamf_temp
            integer:: ibxv_temp
            integer:: itemp
            integer:: nr_isobars_computed, nr_isobars_total
    
            best_ib%error = huge(0.0)
            nr_isobars_computed = 0
            nr_isobars_total = global_data%amf%num*global_data%bxv%num*global_data%grid_iterations * &
                               global_data%amf_temp%num * global_data%bxv_temp%num
            
            iamf_temp = 1
            ibxv_temp = 1
            
            ! Outer loop that decreases the grid size in each iteration
            do igrid = 1, global_data%grid_iterations

                !$OMP PARALLEL DEFAULT(none), &
                !$OMP& PRIVATE(ib, iamf, ibxv, iamf_temp, ibxv_temp, itemp), &
                !$OMP& SHARED(global_data, clusterset, monomer, best_ib, reference, &
                !$OMP& nr_isobars_computed, nr_isobars_total)
                
                !$OMP DO COLLAPSE(4), SCHEDULE(GUIDED)
                do iamf = 1, global_data%amf%num
                    do ibxv = 1, global_data%bxv%num
                        do iamf_temp = 1, global_data%amf_temp%num
                            do ibxv_temp = 1, global_data%bxv_temp%num
                
                                allocate(ib%temp(global_data%temp%num))
                                allocate(ib%vol(global_data%temp%num))
                                allocate(ib%gibbs(global_data%temp%num))
                                allocate(ib%converged(global_data%temp%num))
                                allocate(ib%solution(global_data%temp%num))
                                allocate(ib%lnq(global_data%temp%num, size(clusterset)))
                                allocate(ib%populations(global_data%temp%num, size(clusterset)))
                        
                                ib%amf = global_data%amf%first + (iamf-1)*global_data%amf%delta
                                ib%amf = ib%amf/avogadro**2
                                ib%bxv = global_data%bxv%first + (ibxv-1)*global_data%bxv%delta
        
                                ib%amf_temp = global_data%amf_temp%first + (iamf_temp-1)*global_data%amf_temp%delta
                                ib%amf_temp = ib%amf_temp/avogadro**2
                                ib%bxv_temp = global_data%bxv_temp%first + (ibxv_temp-1)*global_data%bxv_temp%delta
                                
                                do itemp = 1, global_data%temp%num
                                    ib%temp(itemp) = global_data%temp%first &
                                        + (itemp-1)*global_data%temp%delta
                                end do
                        
                                call qce_main(ib)
                        
                                !$OMP CRITICAL
                                global_data%nconverged = global_data%nconverged + count(ib%converged)
                                if (ib%error < best_ib%error) best_ib = ib
                                !$OMP END CRITICAL
                        
                                deallocate(ib%temp)
                                deallocate(ib%vol)
                                deallocate(ib%gibbs)
                                deallocate(ib%converged)
                                deallocate(ib%solution)
                                deallocate(ib%lnq)
                                deallocate(ib%populations)
                        
                                !$OMP ATOMIC
                                nr_isobars_computed = nr_isobars_computed + 1
                                !$OMP END ATOMIC
                
#ifdef _OPENMP
                                if (omp_get_thread_num() == 0) then
                                    call progress_bar(nr_isobars_computed, nr_isobars_total, global_data%progress_bar)
                                end if
#else           
                                    call progress_bar(nr_isobars_computed, nr_isobars_total, global_data%progress_bar)
#endif          
                            end do
                        end do
                    end do
                end do
                !$OMP END DO
                !$OMP END PARALLEL
            
            !TODO Implement decreasing factor as an input parameter
            ! Decreasing grid size by a factor of 1/5 around the current minimum.
            global_data%amf%delta = global_data%amf%delta/5.0_dp
            global_data%bxv%delta = global_data%bxv%delta/5.0_dp
            
            global_data%amf%first = best_ib%amf*avogadro**2 - global_data%amf%delta*(global_data%amf%num-1)/2.0_dp
            global_data%bxv%first = best_ib%bxv - global_data%bxv%delta*(global_data%bxv%num-1)/2.0_dp
            
            global_data%amf_temp%delta = global_data%amf_temp%delta/5.0_dp
            global_data%bxv_temp%delta = global_data%bxv_temp%delta/5.0_dp
            
            global_data%amf_temp%first = best_ib%amf_temp*avogadro**2 - &
                                         global_data%amf_temp%delta*(global_data%amf_temp%num-1)/2.0_dp
            global_data%bxv_temp%first = best_ib%bxv_temp - &
                                         global_data%bxv_temp%delta*(global_data%bxv_temp%num-1)/2.0_dp
            
            end do
            

            ! Start interface mode if specified in input file.
            if (global_data%imode) then
                call interface()
            ! Starts the optimizer if any parameter was chosen to be optimized.
            else if (global_data%optimizer /= 0) then
                call downhill_simplex()
            end if


            call progress_bar(nr_isobars_computed, nr_isobars_total, global_data%progress_bar, newline=.true.)
        end subroutine qce_start
        !=================================================================================
        ! Starts interface mode.
        subroutine interface()
            real(dp), dimension(:), allocatable :: params
            integer:: myunit, stat, n_params, optimizer
            logical:: do_sp
            logical:: stop_im
            
            write(0,*) "Interface mode"

            do_sp = .false.
            stop_im = .false.

            n_params=1
            optimizer = global_data%optimizer
            if (optimizer >= 1000) then
                optimizer = optimizer - 1000
                n_params = n_params + 1
            end if
            if (optimizer >= 100) then
                optimizer = optimizer - 100
                n_params = n_params + 1
            end if
            if (optimizer >= 10) then
                optimizer = optimizer - 10
                n_params = n_params + 1
            end if
            if (optimizer == 1) then
                n_params = n_params + 1
            end if

            allocate(params(n_params))

            ! Starts the actual interface mode. In each iteration of the loop, Peacemaker
            ! will look for a file named "imode.input". If it exists, Peacemaker will read
            ! a set of parameters from it and perform a single point QCE calculation. The
            ! calculated error will be written into a file "sp.out". If the file "stop_imode"
            ! exists, Peacemaker will end the loop and finalize the results.
            do while (.true.)
                inquire(file="imode.input", exist=do_sp)
                inquire(file="stop_imode", exist=stop_im)
                if (stop_im) exit
                if (do_sp) then
                    open(newunit = myunit, iostat=stat, action = "read", status = "unknown", &
                        file = "imode.input")
                    read(myunit,*) params(1:n_params-1)
                    close(myunit, status = "delete")
                    
                    call single_qce(params)
                    
                    open(newunit = myunit, action = "write", status = "unknown", &
                        file = "sp.out", position = "append")
                    write(myunit, '(G0.10)') params(n_params)
                    close(myunit)
                end if
                do_sp = .false.
            end do
            
            open(newunit = myunit, action = "read", status = "unknown", &
                file = "stop_imode")
            close(myunit, status = "delete")

            deallocate(params)
            
        end subroutine interface
        !=================================================================================
        ! Starts the Downhill-Simplex algorithm to optimize amf and bxv. The simplex consists
        ! of three points x1, x2, and x3, each corresponding to a specific combination of
        ! amf and bxv.
        subroutine downhill_simplex()
            real(dp), dimension(:,:), allocatable :: simplex
            real(dp), dimension(:), allocatable :: m, r, e, c
            real(dp), dimension(4,5) :: simplex_start
            real(dp):: diff, crit
            integer:: i, j, n, optimizer
            logical:: opt_a, opt_b, opt_at, opt_bt
            
            opt_a = .false.
            opt_b = .false.
            opt_at = .false.
            opt_bt = .false.
            
            ! The simplex_start array contains a set of vectors that are used to construct
            ! the n-simplex, that is a triangle, tetrahedron, or pentachoron (4-simplex) around
            ! the current best set of parameters.
            simplex_start(1,:) = (/-1.0_dp,  1.0_dp,  0.0_dp,  0.0_dp, 0.0_dp/)
            simplex_start(2,:) = (/-1.0_dp, -1.0_dp,  1.0_dp,  0.0_dp, 0.0_dp/)
            simplex_start(3,:) = (/-1.0_dp, -1.0_dp, -1.0_dp,  1.0_dp, 0.0_dp/)
            simplex_start(4,:) = (/-1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, 1.0_dp/)
            
            n=1
            optimizer = global_data%optimizer
            if (optimizer >= 1000) then
                opt_bt = .true.
                optimizer = optimizer - 1000
                n = n + 1
            end if
            if (optimizer >= 100) then
                opt_at = .true.
                optimizer = optimizer - 100
                n = n + 1
            end if
            if (optimizer >= 10) then
                opt_b = .true.
                optimizer = optimizer - 10
                n = n + 1
            end if
            if (optimizer == 1) then
                opt_a = .true.
                n = n + 1
            end if
            
            allocate(simplex(n,n))
            allocate(m(n))
            allocate(r(n))
            allocate(e(n))
            allocate(c(n))
            
            do i = 1, n
                j = 1
                if (opt_a) then
                    simplex(i,j) = best_ib%amf*avogadro**2.0_dp+simplex_start(j,i)*0.05_dp
                    j = j + 1
                end if
                if (opt_b) then
                    simplex(i,j) = best_ib%bxv+simplex_start(j,i)*0.05_dp
                    j = j + 1
                end if
                if (opt_at) then
                    simplex(i,j) = best_ib%amf_temp*avogadro**2.0_dp+simplex_start(j,i)*0.00005_dp
                    j = j + 1
                end if
                if (opt_bt) then
                    simplex(i,j) = best_ib%bxv_temp+simplex_start(j,i)*0.0005_dp
                end if
            end do
                
            !TODO: Implement convergence criterion as input parameter.
            diff = 1.0e12_dp
            crit = 1.0e-12_dp
                        
            do i = 1, n
                call single_qce(simplex(i,:))
            end do
            
            do while (abs(diff) > crit)
                call reorder(simplex)
                
                diff = simplex(1,n) - simplex(n,n)
                
                call middle(simplex, m)
                call reflection(simplex, m, r)
                
                if (r(n) < simplex(1,n)) then
                    call expansion(simplex, m, e)
                    
                    if (r(n) < e(n)) then
                        simplex(n,:) = r(:)
                    else
                        simplex(n,:) = e(:)
                    end if
                    
                    cycle
                end if
                
                
                if (r(n) < simplex(n-1,n)) then
                    simplex(n,:) = r(:)
                    cycle
                end if
                
                if (r(n) < simplex(n,n)) then
                    call contraction(r, m, c)
                else
                    call contraction(simplex(n,1:n), m, c)
                end if
                
                if (c(n) < simplex(n,n)) then
                    simplex(n,:) = c(:)
                    cycle
                end if
                
                call compression(simplex)
            
            end do
            
            deallocate(simplex)
            deallocate(m)
            deallocate(r)
            deallocate(e)
            deallocate(c)
            
            
        end subroutine downhill_simplex
        !=================================================================================
        ! Reorders simplex x1, x2, ..., xn from best to worst
        subroutine reorder(simplex)
            real(dp), dimension(:,:), allocatable, intent(inout) :: simplex
            real(dp), dimension(:), allocatable :: swap
            integer:: i, n
            logical:: quit
            
            n = size(simplex, 1)
            
            allocate(swap(n))
            
            quit = .false.
            
            do while (.not. quit)
                quit = .true.
                do i = 1, n-1
                    if (simplex(i,n) > simplex(i+1,n)) then
                        swap(:) = simplex(i+1,:)
                        simplex(i+1,:) = simplex(i,:)
                        simplex(i,:) = swap(:)
                        quit = .false.
                    end if
                end do
            end do
            
            deallocate(swap)
            
        end subroutine reorder
        !=================================================================================
        ! Calculates midpoint m of n-1 best points x1, ..., xn-1
        subroutine middle(simplex, m)
            real(dp), dimension(:,:), intent(in) :: simplex
            real(dp), dimension(:), intent(out) :: m
            integer:: i, n
            
            n = size(simplex, 1)
            
            do i = 1, n-1
                m(i) = SUM(simplex(:n-1,i))/real(n-1, dp)
            end do
            
            call single_qce(m)
            
        end subroutine middle
        !=================================================================================
        ! Reflects the worst point xn at midpoint m to form point r
        subroutine reflection(simplex, m, r)
            real(dp), dimension(:,:), intent(in) :: simplex
            real(dp), dimension(:), intent(in) :: m
            real(dp), dimension(:), intent(out) :: r
            real(dp):: alpha
            integer:: i, n
            
            n = size(simplex, 1)
            
            alpha = 1.0_dp
            
            do i = 1, n-1
                r(i) = (1.0_dp + alpha) * m(i) - alpha * simplex(n,i)
            end do
            
            call single_qce(r)
        end subroutine reflection
        !=================================================================================
        ! Expands the worst point xn at midpoint m to form point e
        subroutine expansion(simplex, m, e)
            real(dp), dimension(:,:), intent(in) :: simplex
            real(dp), dimension(:), intent(in) :: m
            real(dp), dimension(:), intent(out) :: e
            real(dp):: gamma
            integer:: i, n
            
            n = size(simplex, 1)
            
            gamma = 2.0_dp
            
            do i = 1, n-1
                e(i) = (1.0_dp + gamma) * m(i) - gamma * simplex(n,i)
            end do
            
            call single_qce(e)
        end subroutine expansion
        !=================================================================================
        ! Contracts the better point h of r or xn with the midpoint m to form point c
        subroutine contraction(h, m, c)
            real(dp), dimension(:), intent(in) :: h, m
            real(dp), dimension(:), intent(out) :: c
            real(dp):: beta
            integer:: i
            
            beta = 0.5_dp
            
            do i = 1, size(c)
                c(i) = beta * m(i) + (1.0_dp - beta) * h(i)
            end do
            
            call single_qce(c)
        
        end subroutine contraction
        !=================================================================================
        ! Compresses the simplex
        subroutine compression(simplex)
            real(dp), dimension(:,:), intent(inout) :: simplex
            real(dp):: sigma
            integer:: i, j, n
            
            n = size(simplex, 1)

            sigma = 0.5_dp
            
            do i = 1, n
                do j = 1, n
                    simplex(i,j) = sigma * simplex(1,j) + (1.0_dp - sigma) * simplex(i,j)
                end do
                call single_qce(simplex(i,:))
            end do

        end subroutine compression
        !=================================================================================
        ! Performs a single QCE calculation for a given point params.
        subroutine single_qce(params)
            real(dp), dimension(:), intent(inout) :: params
            real(dp):: amf
            real(dp):: bxv
            real(dp):: amf_temp
            real(dp):: bxv_temp
            real(dp):: error
            integer:: itemp

            amf = best_ib%amf*avogadro**2
            bxv = best_ib%bxv
            amf_temp = best_ib%amf_temp*avogadro**2
            bxv_temp = best_ib%bxv_temp
            
            !TODO: Find better solution for this.
            select case (global_data%optimizer)
            case (1)
                amf = params(1)
            case (10)
                bxv = params(1)
            case (11)
                amf = params(1)
                bxv = params(2)
            case (100)
                amf_temp = params(1)
            case (101)
                amf = params(1)
                amf_temp = params(2)
            case (110)
                bxv = params(1)
                amf_temp = params(2)
            case (111)
                amf = params(1)
                bxv = params(2)
                amf_temp = params(3)
            case (1000)
                bxv_temp = params(1)
            case (1001)
                amf = params(1)
                bxv_temp = params(2)
            case (1010)
                bxv = params(1)
                bxv_temp = params(2)
            case (1011)
                amf = params(1)
                bxv = params(2)
                bxv_temp = params(3)
            case (1100)
                amf_temp = params(1)
                bxv_temp = params(2)
            case (1101)
                amf = params(1)
                amf_temp = params(2)
                bxv_temp = params(3)
            case (1110)
                bxv = params(1)
                amf_temp = params(2)
                bxv_temp = params(3)
            case (1111)
                amf = params(1)
                bxv = params(2)
                amf_temp = params(3)
                bxv_temp = params(4)
            case default
                amf = params(1)
                bxv = params(2)
            end select

            allocate(ib%temp(global_data%temp%num))
            allocate(ib%vol(global_data%temp%num))
            allocate(ib%gibbs(global_data%temp%num))
            allocate(ib%converged(global_data%temp%num))
            allocate(ib%solution(global_data%temp%num))
            allocate(ib%lnq(global_data%temp%num, size(clusterset)))
            allocate(ib%populations(global_data%temp%num, size(clusterset)))
    
            ib%amf = amf/avogadro**2
            ib%bxv = bxv
            ib%amf_temp = amf_temp/avogadro**2
            ib%bxv_temp = bxv_temp
            do itemp = 1, global_data%temp%num
                ib%temp(itemp) = global_data%temp%first &
                    + (itemp-1)*global_data%temp%delta
            end do
    
            call qce_main(ib)
    
            params(size(params)) = ib%error


            if (ib%error < best_ib%error) then
                best_ib = ib
            end if
    
            deallocate(ib%temp)
            deallocate(ib%vol)
            deallocate(ib%gibbs)
            deallocate(ib%converged)
            deallocate(ib%solution)
            deallocate(ib%lnq)
            deallocate(ib%populations)
    
        end subroutine single_qce
        !=================================================================================
        ! Calculates a QCE isobar.
        subroutine qce_main(ib)
            type(isobar_t), intent(inout) :: ib
    
            real(dp):: amf
            real(dp):: bxv
            real(dp):: v0
            real(dp):: vdamp
            real(dp):: vol
            real(dp):: gibbs
            real(dp):: error
            logical:: converged
            logical:: copy
            type(pf_t), dimension(size(clusterset)) :: lnq
            real(dp), dimension(size(clusterset)) :: populations
    
            integer:: itemp
    
            ! We perform two full QCE cycles for each temperature. In cycle One, we start
            ! from the ideal gas volume at the highest temperature. We decrease
            ! temperature and use the last converged volume as initial guess.
            ! In cycle Two, we start from a very small volume at the lowest
            ! temperature. We increase the temperature and use the last converged volume
            ! as initial guess.
            ! We choose the solution that led to smaller Gibbs free enthalpy.
    
            ! Cycle One
            vdamp = 1.0_dp - global_data%vdamp
            converged = .false.
            do itemp = size(ib%temp), 1, -1
                ! Use volume from previous temperature as initial guess, if available.
                ! Otherwise use ideal gas volume.
                if (.not. converged) then
                    v0 = ideal_gas_volume(ib%temp(itemp))
                else
                    v0 = vol
                end if
                
                bxv = ib%bxv + ib%temp(itemp) * ib%bxv_temp
                
                ! Perform QCE iteration.
                call qce_iteration(0.0_dp, bxv, ib%temp(itemp), v0, vdamp, vol, &
                    gibbs, populations(:), lnq(:), 1, converged)
                ! Copy results.
                ib%gibbs(itemp) = gibbs
                ib%vol(itemp) = vol
                ib%lnq(itemp, :) = lnq(:)
                ib%converged(itemp) = converged
                ib%populations(itemp, :) = populations(:)
                ib%solution(itemp) = 1 ! Note below
                if(converged) ib%solution(itemp) = ib%solution(itemp) + 100 ! Note below.
            end do
    
            ! Cycle Two
            vdamp = 1.0_dp + global_data%vdamp
            converged = .false.
            do itemp = 1, size(ib%temp)
                ! Use volume from previous temperature as initial guess, if available.
                ! Otherwise use damped ideal gas volume.
                if (.not. converged) then
                    v0 = 1.0e-2_dp*ideal_gas_volume(ib%temp(itemp))
                else
                    !v0 = vol
                    ! Test:
                    v0 = min(vol, 1.0e-1*ideal_gas_volume(ib%temp(itemp)))
                end if
                ! Perform QCE iteration.

                amf = ib%amf + ib%temp(itemp) * ib%amf_temp
                bxv = ib%bxv + ib%temp(itemp) * ib%bxv_temp

                call qce_iteration(amf, bxv, ib%temp(itemp), v0, vdamp, vol, &
                    gibbs, populations(:), lnq(:), 2, converged)
                ! Copy results, if necessary.
                if (converged) then
                    ib%solution(itemp) = ib%solution(itemp) + 10 ! Note below.
                    if (.not. ib%converged(itemp)) then
                        copy = .true.
                    else if (gibbs < ib%gibbs(itemp)) then
                        copy = .true.
                    else
                        copy = .false.
                    end if
                    if (copy) then
                        ib%gibbs(itemp) = gibbs
                        ib%vol(itemp) = vol
                        ib%lnq(itemp, :) = lnq(:)
                        ib%converged(itemp) = converged
                        ib%populations(itemp, :) = populations(:)
                        ib%solution(itemp) = ib%solution(itemp) + 1 ! Note below
                    end if
                end if
            end do
    
            ! Determine isobar quality, if necessary.
            ib%error = 0.0_dp
            if (reference%compare_density) then
                call compare_density(ib, reference%density_temperature, &
                    reference%density, error)
                ib%error = ib%error + reference%density_weight*error
            end if
            if (reference%compare_isobar) then
                call compare_isobar(ib, reference%isobar_temperature, &
                    reference%isobar_volume, error)
                ib%error = ib%error + reference%isobar_weight*error
            end if
            if (reference%compare_phase_transition) then
                call compare_phase_transition(ib, reference%phase_transition, error)
                ib%error = ib%error + reference%phase_transition_weight*error
            end if
    
            ! Note on ib%solution.
            ! This array indicates which solutions have converged and which solution was
            ! chosen. The first digit equals 1, if solution 1 has converged. The second
            ! digit equals 1, if solution 2 has converged. The third digits indicates
            ! which solution was chosen.
        end subroutine qce_main
        !=================================================================================
        ! Calculates the ideal gas volume at the given temperature.
        function ideal_gas_volume(temp)
            real(dp), intent(in) :: temp
            real(dp):: ideal_gas_volume
    
            ! The ideal gas volume.
            ideal_gas_volume = sum(global_data%ntot)*kb*temp/global_data%press
        end function ideal_gas_volume
        !=================================================================================
        ! Performs multiple QCE iterations and checks for convergence.
        subroutine qce_iteration(amf, bxv, temp, v0, vdamp, vol, gibbs, populations, &
            lnq, cyclus, converged)
            real(dp), intent(in) :: amf
            real(dp), intent(in) :: bxv
            real(dp), intent(in) :: v0
            real(dp), intent(in) :: vdamp
            real(dp), intent(in) :: temp
            real(dp), intent(out) :: vol
            real(dp), intent(out) :: gibbs
            type(pf_t), dimension(size(clusterset)), intent(out) :: lnq
            real(dp), dimension(size(clusterset)), intent(out) :: populations
            integer, intent(in) :: cyclus
            logical, intent(out) :: converged
    
            integer:: iteration
            logical:: success
            real(dp):: vdamp_local
            real(dp):: old_vol
    
            ! Initialize the volume, populations, and Gibbs free enthalpy.
            vol = v0
            vdamp_local = vdamp
            gibbs = huge(0.0)
            call initialize_populations(populations)

            converged = .false.
            qce_loop: do iteration = 1, global_data%qce_iterations
                ! Calculate the cluster partition functions. In the first iteration, all
                ! of them need to be calculated. Later, we only need to update those,
                ! that depend on the volume.
                if (iteration == 1) then
                    call calculate_lnq(lnq, amf, bxv, temp, vol)
                else
                    call update_lnq(lnq, amf, bxv, temp, vol)
                end if
    
                ! Calculate new populations.
                call calculate_populations(populations, lnq, cyclus, success)
                if (.not. success) then
                    ! Without new populations, we can't proceed to get a new volume. In
                    ! order to proceed anyway, we will use a slightly damped initial
                    ! volume guess. The degree of damping depends on the number of times
                    ! that something in the QCE iteration has failed already.
                    vol = v0*vdamp_local
                    vdamp_local = vdamp_local*vdamp
                    call initialize_populations(populations)
                    cycle qce_loop
                end if
    
                ! Calculate new volume.
                old_vol = vol
                call calculate_volume(vol, vdamp, amf, bxv, temp, populations, success)
                if (.not. success) then
                    ! We couldn't get any physical volume and can't proceed to
                    ! recalculate the partition functions. In order to proceed anyway, we
                    ! will use a slightly damped initial volume guess. The degree of
                    ! damping depends on the number of times that something in the QCE
                    ! iteration has failed already.
                    vol = v0*vdamp_local
                    vdamp_local = vdamp_local*vdamp
                    call initialize_populations(populations)
                    cycle qce_loop
                end if
    
                ! Check for convergence.
                call check_convergence(gibbs, temp, vol, populations, lnq, converged)
                if (converged) exit qce_loop
            end do qce_loop
        end subroutine qce_iteration
        !=================================================================================
        ! Initializes populations. Assumes monomers, only.
        subroutine initialize_populations(populations)
            real(dp), dimension(size(clusterset)), intent(out) :: populations
    
            integer:: i
    
            ! Assume only monomers.
            populations = 0.0_dp
            do i = 1, size(monomer)
                populations(monomer(i)) = global_data%monomer_amounts(i)*avogadro
            end do
        end subroutine initialize_populations
        !=================================================================================
        ! Calculates populations, by solving the corresponding polynomials.
        ! Note that the coefficients of the polynomials in this subroutine are multiplied
        ! by the term (N_1^tot+N_2^tot)**(i+j) for numerical reasons. Thus the results of
        ! the Newton algorithm are populations divided by the same factor.
        subroutine calculate_populations(populations, lnq, cyclus, success)
            use polynomial
            real(dp), dimension(size(clusterset)), intent(out) :: populations
            type(pf_t), dimension(size(clusterset)), intent(in) :: lnq
            logical, intent(out) :: success
            integer, intent(in) :: cyclus
    
            integer:: iclust
            integer:: i, j
            integer:: n_comp
            integer:: indx_clust
            real(dp):: perm
            real(dp):: coeff
            real(dp), dimension(size(monomer)) :: sum_of_composition
            real(dp), dimension(size(monomer)) :: monomer_populations
            real(dp), dimension(size(monomer), size(clusterset)) :: pop_coeffs
            real(dp), dimension(:), allocatable :: coeffs_all
            integer, dimension(:), allocatable :: degree
    
            allocate(degree(size(monomer)))
            degree(:) = global_data%degree(:) + 1

            ! Calculate coefficients for >> each cluster <<
            do iclust = 1, size(clusterset)
                associate(c => clusterset(iclust))
                    ! Instead of the old mass polynomial, new linearly independent variants of the population
                    ! polynomial are created by counting some of the components multiple times, thus, generating
                    ! enough equations to solve the population polynomial for any number of components.
                    do i = 1, size(monomer)
                        sum_of_composition(i) = real(sum(c%composition), dp) + real((i-1) * c%composition(i), dp)
                    end do
    
                    ! TEST
                    ! Calculate: ln[q(i)/(q(1)^i*q(2)^j)] + (i+j)*ln(N_1^tot+N_2^tot)
                    coeff = lnq(iclust)%qtot
                    do i = 1, size(monomer)
                        coeff = coeff - c%composition(i)*lnq(monomer(i))%qtot
                    end do
                    coeff = coeff + sum_of_composition(1)*log(sum(global_data%ntot))
    
                    ! Calculate coefficients of the population polynomial and its variants.
                    do i = 1, size(monomer)
                        pop_coeffs(i, iclust) = exp(coeff)*sum_of_composition(i)/(sum(global_data%ntot)+(i-1)*global_data%ntot(i))
                    end do

                end associate
            end do
!            write(*,*) "lnq", lnq%qtot

            ! Calculate coefficients for >> each possible composition <<. We simulate a multi-dimensional
            ! array as a linear or one-dimensional array. n_comp is the number of possible cluster compositions
            ! and thus the number of coefficients in the population polynomial (some of which may be 0).
            ! indx_clust is the position of a cluster in the array coeffs_all.
            n_comp = 1
            do i = 1, size(degree)
                n_comp = n_comp * degree(i)
            end do
            allocate(coeffs_all(0:n_comp * size(monomer) -1))
            coeffs_all = 0.0_dp
            do i = 0, size(degree)-1
                coeffs_all(i*n_comp) = -1.0_dp
            end do

            do iclust = 1, size(clusterset)
               indx_clust = clusterset(iclust)%composition(1)
               do j = 2, size(monomer)
                   indx_clust = indx_clust + clusterset(iclust)%composition(j) * product(degree(:j-1))
               end do
               do i = 1, size(monomer)
                   perm = GAMMA(real(sum(clusterset(iclust)%composition)+1, dp))
                   do j = 1, size(monomer)
                       perm = perm/GAMMA(real(clusterset(iclust)%composition(j)+1, dp))
                   end do
                   coeffs_all(indx_clust + (i-1)*n_comp) = coeffs_all(indx_clust + (i-1)*n_comp) + pop_coeffs(i,iclust)
               end do
            end do

    
            ! Initial guess of monomer populations.
            if (cyclus == 1) then
                ! Gas phase cycle.
                monomer_populations(:) = &
                    global_data%monomer_amounts(:)*avogadro/sum(global_data%ntot)
            else
                ! Liquid phase cycle.
                monomer_populations(:) = 1.0e-2_dp * &
                    global_data%monomer_amounts(:)*avogadro/sum(global_data%ntot)
            end if
            
            ! Solve the polynomials.
            call newton(size(degree), global_data%degree, size(coeffs_all), coeffs_all, monomer_populations, &
                global_data%newton_iterations, success)
            
            ! Check for unphysical solution.
            if (any(monomer_populations < 0.0_dp)) success = .false.
    
            ! Calculate the remaining populations.
            if (success) then
                ! Get rid of the scaling factor N^tot and calculate the remaining
                ! populations.
                monomer_populations = monomer_populations*sum(global_data%ntot)/1.0_dp
                call calculate_remaining_populations(populations, monomer_populations, &
                    lnq)
            end if
            

            end subroutine calculate_populations
        !=================================================================================
        ! Given the monomer populations, this calculates all other populations.
        subroutine calculate_remaining_populations(populations, monomer_populations, &
            lnq)
            type(pf_t), dimension(size(clusterset)), intent(in) :: lnq
            real(dp), dimension(size(monomer)), intent(in) :: monomer_populations
            real(dp), dimension(size(clusterset)), intent(out) :: populations
    
            integer:: i
            integer:: j
            real(dp):: tmp
    
            ! Assign monomer populations.
            do i = 1, size(monomer)
                populations(monomer(i)) = monomer_populations(i)
            end do
            ! Calculate the remaining populations.
            do i = 1, size(clusterset)
                associate(c => clusterset(i))
                    if (c%monomer) cycle
                    tmp = 0.0_dp
                    do j = 1, size(monomer)
                        tmp = tmp + real(c%composition(j), dp)* &
                            (log(monomer_populations(j)) - lnq(monomer(j))%qtot)
                    end do
                    tmp = tmp + lnq(i)%qtot
                    populations(i) = exp(tmp) * 1.0_dp
                end associate
            end do
        end subroutine calculate_remaining_populations
        !=================================================================================
        ! This subroutine calculates the new volume.
        subroutine calculate_volume(vol, vdamp, amf, bxv, temp, populations, success)
            use polynomial
            real(dp), intent(in) :: vdamp
            real(dp), intent(in) :: amf
            real(dp), intent(in) :: bxv
            real(dp), intent(in) :: temp
            real(dp), intent(out) :: vol
            real(dp), dimension(size(clusterset)), intent(in) :: populations
            logical, intent(out) :: success
    
            integer :: i, j
            real(dp), dimension(0:3) :: coeffs
            complex(dp), dimension(3) :: roots
            logical, dimension(3) :: valid_roots

            ! Calculate the coefficients.
            coeffs(0) = amf*bxv*global_data%vexcl*sum(global_data%ntot)**2
            coeffs(1) = -amf*sum(global_data%ntot)**2    
            coeffs(2) = kb*temp*sum(populations) + global_data%press*bxv*global_data%vexcl
            coeffs(3) = -global_data%press
    
            ! Solve the volume polynomial.
            call solve_polynomial3(coeffs, roots)
    
            ! Check for physical roots.
            where (abs(aimag(roots)) <= global_eps .and. &
                (real(roots, dp) - bxv*global_data%vexcl) >= global_eps)
                valid_roots = .true.
            else where
                valid_roots = .false.
            end where
    
            ! Use the smallest/largest valid volume, depending on the temperature cycle.
            if (any(valid_roots)) then
                if (vdamp <= 1.0_dp) then
                    vol = maxval(real(roots, dp), mask = valid_roots)
                else
                    vol = minval(real(roots, dp), mask = valid_roots)
                end if
                success = .true.
            else
                success = .false.
            end if
        end subroutine calculate_volume
        !=================================================================================
        ! This subroutine checks whether the QCE iteration has converged. It checks the
        ! deviation of the Gibbs free enthalpy from the previous run. The Gibbs free
        ! enthalpy ensures that both populations and volume have converged.
        subroutine check_convergence(gibbs, temp, vol, populations, lnq, &
            converged)
            use auxiliary, only: ln_factorial
            type(pf_t), dimension(size(clusterset)), intent(in) :: lnq
            real(dp), intent(in) :: vol
            real(dp), intent(in) :: temp
            real(dp), intent(inout) :: gibbs
            real(dp), dimension(size(clusterset)), intent(in) :: populations
            logical, intent(out) :: converged
            real(dp):: new_gibbs
            real(dp):: deviation
            integer:: iclust
    
            converged = .false.
            ! Calculate new Gibbs energy.
            new_gibbs = 0.0_dp
            do iclust = 1, size(clusterset)
                new_gibbs = new_gibbs - ln_factorial(populations(iclust)) + &
                    populations(iclust)*lnq(iclust)%qtot
            end do
            new_gibbs = -kb*temp*new_gibbs + global_data%press*vol
            deviation = new_gibbs/gibbs-1.0_dp
            gibbs = new_gibbs
    
            if (abs(deviation) <= global_data%max_deviation) converged = .true.
        end subroutine check_convergence
        !=================================================================================
        ! Compares to an experimental isobar.
        subroutine compare_isobar(ib, temp, vol, error)
            type(isobar_t), intent(in) :: ib
            real(dp), dimension(:), intent(in) :: temp
            real(dp), dimension(:), intent(in) :: vol
            real(dp), intent(out) :: error
    
            integer:: i
            real(dp):: vcalc
    
            error = 0.0_dp
            do i = 1, size(temp)
                vcalc = determine_volume(ib, temp(i))
                error = error + ((vcalc - vol(i))/vol(i))**2
            end do
            error = error / real(size(temp), dp)
        end subroutine compare_isobar
        !=================================================================================
        ! Compare to an experimental density at a given temperature.
        subroutine compare_density(ib, temp, density, error)
            type(isobar_t), intent(in) :: ib
            real(dp), intent(in) :: temp
            real(dp), intent(in) :: density
            real(dp), intent(out) :: error
    
            real(dp):: vol
            real(dp):: vcalc
    
            vol = global_data%mtot/(1.0e3_dp*density)
            vcalc = determine_volume(ib, temp)
            error = ((vcalc - vol)/vol)**2
        end subroutine compare_density
        !=================================================================================
        ! Compares to an experimental phase transition.
        subroutine compare_phase_transition(ib, pt, error)
            type(isobar_t), intent(in) :: ib
            real(dp), intent(in) :: pt
            real(dp), intent(out) :: error
    
            real(dp):: ptcalc
    
            ptcalc = determine_phase_transition(ib)
            error = ((ptcalc - pt)/pt)**2
        end subroutine compare_phase_transition
        !=================================================================================
        ! Determines the volume of a QCE isobar at a given temperature.
        ! Interpolates linearly.
        function determine_volume(ib, temp)
            type(isobar_t), intent(in) :: ib
            real(dp), intent(in) :: temp
            real(dp):: determine_volume
    
            integer:: itemp
            real(dp):: x
    
            ! Narrow down temperature frame.
            search: do itemp = 1, size(ib%temp)
                if (ib%temp(itemp) >= temp) exit search
            end do search
    
            ! Temperature is between [itemp, itemp+1). Calculate interpolation factor.
            x = (temp - ib%temp(itemp))/(ib%temp(itemp+1) - ib%temp(itemp))
    
            ! Interpolate volume.
            determine_volume = (1.0_dp - x)*ib%vol(itemp) + x*ib%vol(itemp+1)
        end function determine_volume
        !=================================================================================
        ! Determines the phase transition of a QCE isobar, empirically by looking for the
        ! largest volume jump. This may not actually be a phase transition.
        function determine_phase_transition(ib)
            type(isobar_t), intent(in) :: ib
            real(dp), dimension(2) :: m, b
            real(dp):: determine_phase_transition
    
            integer:: itemp
            integer:: itemp_max
            real(dp):: dv
            real(dp):: dv_max
    
            dv_max = 0.0_dp
            itemp_max = 2
            do itemp = 2, size(ib%temp)
                if (ib%converged(itemp) .and. ib%converged(itemp-1) &
                    .and. (ib%solution(itemp-1) - ib%solution(itemp) == 1)) then
                    dv = ib%vol(itemp) - ib%vol(itemp-1)
                    if (dv > dv_max) then
                        dv_max = dv
                        itemp_max = itemp
                    end if
                end if
            end do
            
            ! Phase transition point is calculated by approximating the Gibbs energy's slope before
            ! and after the volume jump
            if ((itemp_max+1 <= size(ib%temp)) .and. (itemp_max-2 > 0)) then
                m(1) = (ib%gibbs(itemp_max+1) - ib%gibbs(itemp_max))/(ib%temp(itemp_max+1) - ib%temp(itemp_max))
                m(2) = (ib%gibbs(itemp_max-1) - ib%gibbs(itemp_max-2))/(ib%temp(itemp_max-1) - ib%temp(itemp_max-2))
                
                b(1) = ib%gibbs(itemp_max) - m(1) * ib%temp(itemp_max)
                b(2) = ib%gibbs(itemp_max-1) - m(2) * ib%temp(itemp_max-1)

                determine_phase_transition = &
                    (b(1) - b(2))/(m(2) - m(1))
            else 
                determine_phase_transition = itemp_max
            end if
        end function determine_phase_transition
        !=================================================================================
        ! This subroutine prints info about the QCE calculation, calculates properties
        ! and writes output files.
        subroutine qce_finalize()
            ! Determine number of iterations and number of converged iterations both
            ! total and for the best isobar.
            write(*,'(4X,A,1X,G0,A,G0)') "Number of converged iterations:", &
                count(best_ib%converged), "/", global_data%temp%num
            if (reference%compare) then
                write(*,'(4X,A,1X,G0,A,G0)') &
                    "Number of converged iterations (total):", &
                    global_data%nconverged, "/", &
                    global_data%amf%num*global_data%bxv%num*global_data%temp%num*global_data%grid_iterations
                write(*,*)

                select case (size(monomer))
                case (1)
                    write(*,'(4X,A,G0.10,A,G0.10)') "Best isobar found for: amf = ", &
                        best_ib%amf*avogadro**2, " J*m^3/mol^2, bxv = ", &
                        best_ib%bxv
                    write(*,'(4X,A,G0.10,A,G0.10)') "Best isobar found for: amf_temp = ", &
                        best_ib%amf_temp*avogadro**2, " J*m^3/mol^2, bxv_temp = ", &
                        best_ib%bxv_temp
                case default
                    write(*,'(4X,A,G0.6,A,G0.6)') "Best isobar found for: amf_mix = ", &
                        best_ib%amf*avogadro**2, " J*m^3/mol^2, bxv_mix = ", &
                        best_ib%bxv
                end select
                write(*,'(4X,A,G0.6)') "Error: ", best_ib%error
                if (reference%compare_phase_transition) write(*,'(4X,A,G0.6,A)') &
                    "Calculated phase transition: ", &
                    determine_phase_transition(best_ib), " K"
                if (reference%compare_density) write(*,'(4X,A,G0.6,A)') &
                    "Calculated density: ", &
                    1.0e3_dp*global_data%mtot / &
                    (1.0e6_dp*determine_volume(best_ib, reference%density_temperature)), &
                    " g/cm^3"
            end if
            write(*,*)
    
            ! Perform Post-Processing.
            call post_processing(global_data%temp%num, size(clusterset), clusterset(:)%label, best_ib%temp, &
                best_ib%lnq, best_ib%populations, best_ib%vol, &
                best_ib%bxv*global_data%vexcl, best_ib%amf, best_ib%bxv, best_ib%amf_temp, best_ib%bxv_temp, &
                global_data%press, sum(global_data%ntot), best_ib%solution, best_ib%converged)
        end subroutine qce_finalize
        !=================================================================================
        ! This is the place, where all thermodynamic properties are caclulated using
        ! (hopefully) easy-to understand code. The idea is that this subroutine gets
        ! a couple of arrays containg parition functions, volume, populations, etc.
        ! and calculates everything there is to calculate from these in one place.
        ! Results are also written to file here.
        subroutine post_processing(ntemp, nclust, labels, temp, lnq_clust, pop, vol, vexcl, amf, &
            bxv, amf_temp, bxv_temp, press, ntot, solution, converged)
            use auxiliary, only: ln_factorial, derivative
            use input, only: pmk_input
            use lengths
            use iso_varying_string
            use thermo
            integer, intent(in) :: ntemp
            integer, intent(in) :: nclust
            real(dp), intent(in) :: press
            real(dp), intent(in) :: ntot
            real(dp), intent(in) :: vexcl
            type(pf_t), dimension(ntemp, nclust) :: lnq_clust
            logical, dimension(ntemp), intent(in) :: converged
            integer, dimension(ntemp), intent(in) :: solution
            real(dp), dimension(ntemp), intent(in) :: vol
            real(dp), dimension(ntemp), intent(in) :: temp
            real(dp), intent(in) :: amf
            real(dp), intent(in) :: bxv
            real(dp), intent(in) :: amf_temp
            real(dp), intent(in) :: bxv_temp
            real(dp), dimension(ntemp, nclust), intent(in) :: pop
            type(varying_string), dimension(nclust), intent(in) :: labels

            integer:: iclust
            integer:: itemp
            integer:: myunit
            type(pf_t), dimension(ntemp, nclust) :: dlnq_clust, ddlnq_clust
            real(dp), dimension(ntemp) :: dvol
            real(dp), dimension(ntemp, nclust) :: pop_norm
            real(dp), dimension(ntemp, nclust) :: conc
            type(pf_t), dimension(nclust) :: d, dd
            character(fmt_len) :: fmtspec

            ! Calculate derivatives.
            do itemp = 1, ntemp
                call calculate_dlnq(d, amf, bxv, amf_temp, bxv_temp, temp(itemp), vol(itemp))
                call calculate_ddlnq(dd, amf, bxv, amf_temp, bxv_temp, temp(itemp), vol(itemp))
                dlnq_clust(itemp, :) = d
                ddlnq_clust(itemp, :) = dd
            end do
            dvol = derivative(vol, temp)

            ! Calculate and write thermodynamic quantities.
            call calculate_thermo(ntemp, nclust, pop, press, temp, vol, vexcl, &
                lnq_clust, dlnq_clust, ddlnq_clust, dvol, solution, converged)
            
            if (pmk_input%contrib) then
                call calculate_cluster_thermo(ntemp, nclust, pop, press, temp, vol, &
                    lnq_clust, dlnq_clust, ddlnq_clust, dvol, converged, labels)
            end if


            ! Normalize populations.
            do iclust = 1, nclust
                pop_norm(:, iclust) = &
                    pop(:, iclust)*real(sum(clusterset(iclust)%composition), dp)/ntot
            end do
    
            ! Write populations.
            write(fmtspec, '(A,G0,A)') '(A1,A12,', nclust, '(1X,A13))' 
            open(newunit = myunit, action = "write", status = "unknown", &
                file = "populations.dat")
            write(myunit, fmtspec) "#", "T/K", (char(labels(iclust)), iclust = 1, nclust)
            write(fmtspec, '(A,G0,A)') '(ES13.6,', nclust, '(1X,ES13.6))' 
            do itemp = 1, ntemp
                if (converged(itemp)) write(myunit, fmtspec) temp(itemp), &
                    (pop_norm(itemp, iclust), iclust = 1, nclust)
            end do
            close(myunit)
    
            ! Calculate concentrations.
            do iclust = 1, nclust
                conc(:, iclust) = pop(:, iclust)/(avogadro*vol*1.0e3_dp)
            end do
    
            ! Write concentrations.
            write(fmtspec, '(A,G0,A)') '(A1,A12,', nclust, '(1X,A13))' 
            open(newunit = myunit, action = "write", status = "unknown", &
                file = "concentrations.dat")
            write(myunit, fmtspec) "#", "T/K", (char(labels(iclust)), iclust = 1, nclust)
            write(fmtspec, '(A,G0,A)') '(ES13.6,', nclust, '(1X,ES13.6))' 
            do itemp = 1, ntemp
                if (converged(itemp)) write(myunit, fmtspec) temp(itemp), &
                    (conc(itemp, iclust), iclust = 1, nclust)
            end do
            close(myunit)
        end subroutine post_processing
        !=================================================================================
end module qce
!=========================================================================================
