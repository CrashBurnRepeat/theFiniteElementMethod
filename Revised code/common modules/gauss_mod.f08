module gauss_mod
  use init_mod, only : ngaus

contains
  subroutine Gauss (pos, w1)
    implicit none

    real, allocatable, intent (in out) :: pos(:)
    real, allocatable, intent (in out) :: w1(:)

    allocate (pos(ngaus))
    allocate (w1(ngaus))

    select case (ngaus)
    case (1)
      pos = 0.0
      w1  = 2.0
    case (2)
      pos = [-0.5773502692, 0.5773502692]
      w1  = [1.0, 1.0]
    case (3)
      pos = [-0.7745966692, 0.0, 0.7745966692]
      w1  = [0.55555555555, 0.88888888888, 0.55555555555]
    case (4)
      pos = [-0.8611363116, -0.3399810436, &
            0.3399810436, 0.8611363116]
      w1  = [0.3478548451, 0.6521451549,&
            0.6521451549, 0.3478548451]
    case default
      print *, ' Invalid value for number of Gauss points'
      stop
    end select
  end subroutine Gauss
end module gauss_mod
