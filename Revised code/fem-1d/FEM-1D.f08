program FEM-1D
  use init_mod
  use gauss_mod
  use assemb_mod
  implicit none

  real, allocatable  :: pos(:)
  real, allocatable  :: w1(:)

  call Init ()
  call Gauss (ngaus, pos, w1)

  cnew = cold ! Already set?

  call Assemb ()
end program FEM-1D
