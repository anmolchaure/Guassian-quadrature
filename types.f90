module types
    implicit none
    integer, parameter :: ikind  = selected_int_kind(9)
    integer, parameter :: rskind = selected_real_kind(p=6,  r=37)
    integer, parameter :: rdkind = selected_real_kind(p=15, r=300)
    integer, parameter :: rqkind = selected_real_kind(p=33, r=4000)
end module types
