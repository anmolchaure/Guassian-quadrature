program numeric
    use types
    use gauss_quad
    implicit none

    ! --- parameters ---
    real(rdkind), parameter :: a = 0.0_rdkind
    real(rdkind), parameter :: b = 1.5_rdkind

    ! --- Gaussian quadrature ---
    integer      :: order
    real(rdkind) :: gq_result
    integer      :: gq_fcalls

    ! --- trapezoid ---
    integer, parameter :: n_trap_cases = 13
    integer      :: trap_panels(n_trap_cases)
    integer      :: k
    real(rdkind) :: trap_result
    integer      :: trap_fcalls

    ! --- reference (high-order Gauss, order=6 is our best) ---
    real(rdkind) :: ref_result
    integer      :: ref_fcalls

    ! ----------------------------------------------------------------
    print *, repeat('-', 65)
    print *, 'Numerical integration of f(x) = exp(x^2)  on [0.0, 1.5]'
    print *, repeat('-', 65)

    ! ================================================================
    ! Part 1 – Gaussian quadrature, orders 1 to 6
    ! ================================================================
    print *
    print *, '=== GAUSSIAN QUADRATURE (orders 1 to 6) ==='
    print '(a6, 2x, a22, 2x, a10)', 'Order', 'Approximation', 'f-calls'
    print *, repeat('-', 45)

    gauss_orders: do order = 1, 6
        call gauss_integrate(order, a, b, gq_result, gq_fcalls)
        print '(i6, 2x, f22.15, 2x, i6)', order, gq_result, gq_fcalls
    end do gauss_orders

    ! Use order=6 as our reference for comparison
    call gauss_integrate(6, a, b, ref_result, ref_fcalls)

    print *
    print *, 'Reference value (order=6 Gauss): ', ref_result

    ! ================================================================
    ! Part 2 – Trapezoid method, increasing panel counts
    ! ================================================================
    print *
    print *, '=== TRAPEZOID METHOD ==='
    print '(a10, 2x, a22, 2x, a14, 2x, a8)', &
          'Panels', 'Approximation', 'Abs Error', 'f-calls'
    print *, repeat('-', 62)

    trap_panels = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192]

    trap_cases: do k = 1, n_trap_cases
        call trapezoid_integrate(trap_panels(k), a, b, trap_result, trap_fcalls)
        print '(i10, 2x, f22.15, 2x, es14.4, 2x, i8)', &
              trap_panels(k), trap_result, abs(trap_result - ref_result), trap_fcalls
    end do trap_cases

    print *
    print *, repeat('-', 65)
    print *, 'Done.'

end program numeric
