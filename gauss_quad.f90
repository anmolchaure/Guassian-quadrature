module gauss_quad
    use types
    implicit none


    integer, parameter :: max_half = 3

    contains

    ! ================================================================
    ! Integrand: f(x) = exp(x^2)
    ! ================================================================
    pure function f(x) result(val)
        real(rdkind), intent(in) :: x
        real(rdkind)             :: val
        val = exp(x * x)
    end function f



    subroutine load_gl_nodes_weights(order, nodes, weights)
        integer,      intent(in)  :: order
        real(rdkind), intent(out) :: nodes(order), weights(order)

        real(rdkind) :: np(max_half), wp(max_half)

        if (order == 1) then
            ! degree 1: single node at 0, weight 2
            np(1) = 0.0000000000000000_rdkind
            wp(1) = 2.0000000000000000_rdkind

        else if (order == 2) then
            ! degree 2: one positive node (mirror gives negative)
            np(1) = 0.57735026918962573_rdkind
            wp(1) = 1.0000000000000000_rdkind

        else if (order == 3) then
            ! degree 3: one positive node + zero node
            np(1) = 0.77459666924148340_rdkind
            wp(1) = 0.55555555555555547_rdkind
            np(2) = 0.0000000000000000_rdkind
            wp(2) = 0.88888888888888906_rdkind

        else if (order == 4) then
            ! degree 4: two positive nodes
            np(1) = 0.86113631159405246_rdkind
            wp(1) = 0.34785484513745413_rdkind
            np(2) = 0.33998104358485604_rdkind
            wp(2) = 0.65214515486254587_rdkind

        else if (order == 5) then
            ! degree 5: two positive nodes + zero node
            np(1) = 0.90617984593866419_rdkind
            wp(1) = 0.23692688505618875_rdkind
            np(2) = 0.53846931010568377_rdkind
            wp(2) = 0.47862867049936592_rdkind
            np(3) = 0.0000000000000000_rdkind
            wp(3) = 0.56888888888889078_rdkind

        else if (order == 6) then
            ! degree 6: three positive nodes
            np(1) = 0.93246951420315261_rdkind
            wp(1) = 0.17132449237916891_rdkind
            np(2) = 0.66120938646626592_rdkind
            wp(2) = 0.36076157304813894_rdkind
            np(3) = 0.23861918608319749_rdkind
            wp(3) = 0.46791393457269215_rdkind

        else
            error stop 'load_gl_nodes_weights: order must be 1 to 6'
        end if

        ! --- expand to full symmetric node set ---
        if (order == 1) then
            ! single node at zero
            nodes(1)   =  np(1)
            weights(1) =  wp(1)

        else if (order == 2) then
            ! -x1, +x1
            nodes(1)   = -np(1);  weights(1) = wp(1)
            nodes(2)   =  np(1);  weights(2) = wp(1)

        else if (order == 3) then
            ! -x1, 0, +x1
            nodes(1)   = -np(1);  weights(1) = wp(1)
            nodes(2)   =  np(2);  weights(2) = wp(2)   ! zero node
            nodes(3)   =  np(1);  weights(3) = wp(1)

        else if (order == 4) then
            ! -x1, -x2, +x2, +x1  (x1 > x2)
            nodes(1)   = -np(1);  weights(1) = wp(1)
            nodes(2)   = -np(2);  weights(2) = wp(2)
            nodes(3)   =  np(2);  weights(3) = wp(2)
            nodes(4)   =  np(1);  weights(4) = wp(1)

        else if (order == 5) then
            ! -x1, -x2, 0, +x2, +x1
            nodes(1)   = -np(1);  weights(1) = wp(1)
            nodes(2)   = -np(2);  weights(2) = wp(2)
            nodes(3)   =  np(3);  weights(3) = wp(3)   ! zero node
            nodes(4)   =  np(2);  weights(4) = wp(2)
            nodes(5)   =  np(1);  weights(5) = wp(1)

        else if (order == 6) then
            ! -x1, -x2, -x3, +x3, +x2, +x1
            nodes(1)   = -np(1);  weights(1) = wp(1)
            nodes(2)   = -np(2);  weights(2) = wp(2)
            nodes(3)   = -np(3);  weights(3) = wp(3)
            nodes(4)   =  np(3);  weights(4) = wp(3)
            nodes(5)   =  np(2);  weights(5) = wp(2)
            nodes(6)   =  np(1);  weights(6) = wp(1)
        end if

    end subroutine load_gl_nodes_weights

    ! ================================================================
    ! Gaussian quadrature on [a, b].
    !
    ! Change of variable: t in [-1,1] -> x in [a,b]
    !   x = 0.5*(b-a)*t + 0.5*(b+a)
    !   dx = 0.5*(b-a) dt
    !
    ! Result: 0.5*(b-a) * sum_i  w_i * f(x_i)
    ! ================================================================
    subroutine gauss_integrate(order, a, b, result, fcalls)
        integer,      intent(in)  :: order
        real(rdkind), intent(in)  :: a, b
        real(rdkind), intent(out) :: result
        integer,      intent(out) :: fcalls

        real(rdkind) :: nodes(order), weights(order)
        real(rdkind) :: half_len, mid, x
        integer      :: i

        call load_gl_nodes_weights(order, nodes, weights)

        half_len = 0.5_rdkind * (b - a)
        mid      = 0.5_rdkind * (b + a)
        result   = 0.0_rdkind

        sum_nodes: do i = 1, order
            x      = half_len * nodes(i) + mid
            result = result + weights(i) * f(x)
        end do sum_nodes

        result = half_len * result
        fcalls = order

    end subroutine gauss_integrate

    ! ================================================================
    ! Composite trapezoid rule on [a, b] with n_panels panels.
    !
    !   h   = (b - a) / n
    !   I  ~= h * [ f(a)/2 + f(x_1) + ... + f(x_{n-1}) + f(b)/2 ]
    !
    ! f-calls = n_panels + 1  (one per node including endpoints)
    ! ================================================================
    subroutine trapezoid_integrate(n_panels, a, b, result, fcalls)
        integer,      intent(in)  :: n_panels
        real(rdkind), intent(in)  :: a, b
        real(rdkind), intent(out) :: result
        integer,      intent(out) :: fcalls

        real(rdkind) :: h, x
        integer      :: i

        h      = (b - a) / real(n_panels, rdkind)
        result = 0.5_rdkind * (f(a) + f(b))

        interior: do i = 1, n_panels - 1
            x      = a + real(i, rdkind) * h
            result = result + f(x)
        end do interior

        result = h * result
        fcalls = n_panels + 1

    end subroutine trapezoid_integrate

end module gauss_quad
