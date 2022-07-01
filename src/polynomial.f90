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
! This module implements algorithms for solving polynomials.
!=========================================================================================
module polynomial
    use kinds
    use cluster
    use constants
    implicit none
    private
    !=====================================================================================
    ! Public entities.
    public :: newton
    public :: solve_polynomial3
    !=====================================================================================
    ! Convergence criterion in the Newton algorithm.
    real(dp), parameter :: newton_convergence = global_eps
    !=====================================================================================
    contains
        !=================================================================================
        ! The Newton-Raphson algorithm for the simultaneous solution of n n-dimensional 
        ! polynomials. For dimensions 1-3, the inverse of the Jacobi matrix is calculated
        ! directly. For higher dimensions an LU-decomposition is called.
        subroutine newton(n, d, l, coeffs, x, iterations, success)
            integer, intent(in) :: n                    ! Dimension of system
            integer, dimension(n), intent(in) :: d      ! Degree of each dimension
            integer :: l                                 ! Length of coeffs
            real(dp), dimension(0:l), intent(in) :: coeffs   ! Coefficients
            real(dp), dimension(n), intent(inout) :: x
            integer, intent(in) :: iterations
            logical, intent(out) :: success

            integer :: i, j
            real(dp), dimension(n) :: p, dx, xnew
            real(dp), dimension(n, n) :: pdiff, inverse   ! Jacobi matrix and its inverse
            real(dp) :: det

                l = l/n
            newton_loop: do i = 1, iterations
                do j = 1, n
                    call horner(n, d, l, coeffs((j-1)*l:j*l-1), x, p(j), pdiff(:, j))
                end do
                pdiff = transpose(pdiff)
                
                ! Check for convergence.
                if (all(abs(p) <= newton_convergence)) then
                    success = .true.
                    return
                end if
                
                select case (size(monomer))
                    case (1)
                        dx = -p/pdiff(1,1)
                    case (2)
                        det = pdiff(1,1)*pdiff(2,2) - pdiff(1,2)*pdiff(2,1)
                        dx(1) = p(2)*pdiff(1,2) - p(1)*pdiff(2,2)
                        dx(2) = p(1)*pdiff(2,1) - p(2)*pdiff(1,1)
                        dx = dx/det
                    case (3)
                        ! Calculate elements of inverseed 3x3 matrix
                        inverse(1,1) =   pdiff(2,2)*pdiff(3,3) - pdiff(2,3)*pdiff(3,2)
                        inverse(1,2) = -(pdiff(2,1)*pdiff(3,3) - pdiff(2,3)*pdiff(3,1))
                        inverse(1,3) =   pdiff(2,1)*pdiff(3,2) - pdiff(2,2)*pdiff(3,1)
                        inverse(2,1) = -(pdiff(1,2)*pdiff(3,3) - pdiff(1,3)*pdiff(3,2))
                        inverse(2,2) =   pdiff(1,1)*pdiff(3,3) - pdiff(1,3)*pdiff(3,1)
                        inverse(2,3) = -(pdiff(1,1)*pdiff(3,2) - pdiff(1,2)*pdiff(3,1))
                        inverse(3,1) =   pdiff(1,2)*pdiff(2,3) - pdiff(1,3)*pdiff(2,2)
                        inverse(3,2) = -(pdiff(1,1)*pdiff(2,3) - pdiff(1,3)*pdiff(2,1))
                        inverse(3,3) =   pdiff(1,1)*pdiff(2,2) - pdiff(1,2)*pdiff(2,1)
        
                        ! Determine Newton step.
                        det = pdiff(1,1)*inverse(1,1) + pdiff(1,2)*inverse(1,2) + pdiff(1,3)*inverse(1,3)
                        
                        dx(1) = inverse(1,1)*p(1) + inverse(2,1)*p(2) + inverse(3,1)*p(3)
                        dx(2) = inverse(1,2)*p(1) + inverse(2,2)*p(2) + inverse(3,2)*p(3)
                        dx(3) = inverse(1,3)*p(1) + inverse(2,3)*p(2) + inverse(3,3)*p(3)
                        
                        dx = -dx/det
                    case default
                        call lu_dcmp(n, pdiff, p)
                        dx = -p
                end select

                ! Half step size, if this would move x out of physical bounds.
                do j = 1, n
                    xnew(j) = x(j) + dx(j)
                    do while (xnew(j) < 0.0_dp .or. xnew(j) > 1.0_dp)
                        dx = 0.5 * dx
                        xnew(j) = x(j) + dx(j)
                    end do
                end do
                

                ! Move solution closer to zero.
                x = xnew
            end do newton_loop
            success = .false.
        end subroutine newton
        !=================================================================================
        ! LU-decomposition algorithm taken from "Numerical Recipes" by William H. Press et al.
        ! Solves a matrix equation A * x = b and returns x. For simplicity, the solution x is
        ! stored in the same array as b, and the LU-decomposed matrix is stored in the same
        ! array as A.
        subroutine lu_dcmp(n, lu, x)
            integer, intent(in) :: n
            real(dp), dimension(n, n) :: lu
            real(dp), dimension(n), intent(inout) :: x
            
            integer :: i, j, k, imax, ii, ip
            integer, dimension(n) :: indx
            real(dp) :: SMALL, big, temp, d, summ
            real(dp), dimension(n) :: vv

            imax = 0
            indx(:) = 0
            SMALL = 1.0e-30_dp
            
            ! Decompose
            do i = 1, n
                big = 0.0_dp
                do j = 1, n
                    temp = abs(lu(i,j))
                    if (abs(temp) > big) big = temp
                end do
                vv(i) = 1.0_dp/big
            end do
            
            do j = 1, n
                do i = 1, j-1
                    summ = lu(i,j)
                    do k = 1, i-1
                        summ = summ - lu(i,k) * lu(k,j)
                    end do
                    lu(i,j) = summ
                end do
                

                big = 0.0_dp
                
                do i = j, n
                    summ = lu(i,j)
                    do k = 1, j-1
                        summ = summ - lu(i,k)*lu(k,j)
                    end do
                    lu(i,j) = summ
                    
                    temp = vv(i) * abs(summ)
                    if (temp >= big) then
                        big = temp
                        imax = i
                    end if
                end do
                
                if (j /= imax) then
                    do k = 1, n
                        temp = lu(imax,k)
                        lu(imax,k) = lu(j,k)
                        lu(j,k) = temp
                    end do
                    d = -d
                    vv(imax) = vv(j)
                end if
                
                indx(j) = imax
                if (abs(lu(j,j)) < SMALL) lu(j,j) = SMALL
                if (j /= n) then
                    temp = 1.0_dp/lu(j,j)
                    do i = j+1, n
                        lu(i,j) = lu(i,j) * temp
                    end do
                end if
            end do
            
            ! Solve
            ii = 0
            do i = 1, n
                ip = indx(i)
                summ = x(ip)
                x(ip) = x(i)
                if (ii /= 0) then
                    do j = ii, i-1
                        summ = summ - lu(i,j) * x(j)
                    end do
                else if (abs(summ) > SMALL) then
                    ii = i
                end if
                x(i) = summ
            end do
            
            do i = n, 1, -1
                summ = x(i)
                do j = i+1, n
                    summ = summ - lu(i,j) * x(j)
                end do
                x(i) = summ/lu(i,i)
            end do

        end subroutine lu_dcmp
        !=================================================================================
        ! Solves a third order polynomial.
        subroutine solve_polynomial3(coeffs, roots)
            real(dp), dimension(0:3), intent(in) :: coeffs
            complex(dp), dimension(3), intent(out) :: roots

            complex(dp), dimension(6) :: cof
            complex(dp), parameter :: im = (0.0_dp,1.0_dp)
            complex(dp), parameter :: one_third = cmplx(1.0_dp/3.0_dp, 0.0_dp, dp)
            real(dp), parameter :: two_power_one_third = 2.0_dp**(1.0_dp/3.0_dp)
            real(dp), parameter :: sqrt3 = sqrt(3.0_dp)
            
            roots = (0.0_dp, 0.0_dp)
            cof = (0.0_dp, 0.0_dp)
            
            cof(1) = cmplx(-coeffs(2)/3.0_dp/coeffs(3), 0.0_dp, kind = dp)
            cof(2) = cmplx(-coeffs(2)**2 + 3.0_dp*coeffs(3)*coeffs(1), 0.0_dp, kind = dp)
            cof(3) = cmplx(-2.0_dp*coeffs(2)**3 + 9.0_dp*coeffs(3)*coeffs(2)*coeffs(1) &
                -27.0_dp*coeffs(0)*coeffs(3)**2, 0.0_dp, kind = dp)
            cof(4) = 4.0_dp*cof(2)**3 + cof(3)**2
            cof(5) = cmplx(3.0_dp*coeffs(3), 0.0_dp, kind = dp)
            cof(6) = (cof(3) + sqrt(cof(4)))**one_third
            
            roots(1) = cof(1) - two_power_one_third*cof(2)/(cof(5)*cof(6)) + &
                cof(6)/two_power_one_third/cof(5)
            roots(2) = cof(1) + &
                (1.0_dp + im*sqrt3)*cof(2)/(two_power_one_third**2*cof(5)*cof(6)) - &
                (1.0_dp - im*sqrt3)*cof(6)/(2.0_dp*cof(5)*two_power_one_third)
            roots(3) = cof(1) + &
                (1.0_dp - im*sqrt3)*cof(2)/(two_power_one_third**2*cof(5)*cof(6)) - &
                (1.0_dp + im*sqrt3)*cof(6)/(2.0_dp*cof(5)*two_power_one_third)
        end subroutine solve_polynomial3
        !=================================================================================
        ! Uses Horner's scheme to evaluate an n-dimensional polynomial of degree (d(1),d(2),d(3),...)
        ! and its derivative at a point x, given the coefficients c.
        recursive subroutine horner(n, d, l, c, x, p, pdiff)
            integer, intent(in) :: n, l                     ! number of dimensions and length of c
            integer, dimension(n), intent(in) :: d          ! degree of each dimension
            real(dp), dimension(0:l-1), intent(in) :: c     ! coefficients
            real(dp), dimension(n), intent(in) :: x
            real(dp), intent(out) :: p
            real(dp), dimension(n), intent(out) :: pdiff

            integer :: i, j, m, o
            real(dp) :: b
            real(dp), dimension(n-1) :: bdiff

            if (n == 1) then
                p = c(d(1))
                pdiff = 0.0_dp
                do i = d(1)-1, 0, -1
                    pdiff = x(1)*pdiff + p
                    p = p*x(1) + c(i)
                end do
                return
            end if

            m = 1
            do i = 1, n-1
                m = m * (d(i)+1)
            end do
            m = m * d(n)
            

            call horner(n-1, d(1:n-1), size(c(m:)), c(m:), x(1:n-1), b, bdiff)
            p = b
            do i=1, n-1
                pdiff(i) = bdiff(i)
            end do
            pdiff(n) = 0.0_dp
            do i = d(n), 1, -1
                m = 1
                do j = 1, n-1
                    m = m * (d(j)+1)
                end do
                o = m-1
                m = m * (i-1)

                call horner(n-1, d(1:n-1), size(c(m:m+o)), c(m:m+o), x(1:n-1), b, bdiff)
                do j=1, n-1
                    pdiff(j) = x(n)*pdiff(j) + bdiff(j)
                end do
                pdiff(n) = x(n)*pdiff(n) + p
                p = x(n)*p + b
            end do
        end subroutine horner
        !=================================================================================
end module polynomial
!=========================================================================================
