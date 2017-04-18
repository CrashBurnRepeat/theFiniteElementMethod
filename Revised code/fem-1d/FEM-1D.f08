program FEM_1D
  use init_mod
  use gauss_mod
  use assemb_mod
  implicit none

  integer            :: kount
  integer            :: ntime
  real               :: time
  real, allocatable  :: pos(:)
  real, allocatable  :: w1(:)

  call Init ()
  call Gauss (pos, w1)

  cnew = cold ! Already set?

  call Assemb (w1)

  kount = 0
  ntime = 0
  time  = 0.0

  !call Bndcon ()

end program FEM_1D
