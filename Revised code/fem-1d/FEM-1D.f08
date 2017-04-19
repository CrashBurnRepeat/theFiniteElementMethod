program FEM_1D
  use init_mod
  use gauss_mod
  use assemb_mod
  implicit none

  integer            :: kount
  integer            :: nstep
  real               :: time
  real, allocatable  :: pos(:)
  real, allocatable  :: w1(:)

  call Init ()

  call Gauss (pos, w1)

  cnew = cold ! Already set?

  call Assemb (pos, w1)

  kount = 0
  ntime = 0
  time  = 0.0

  call Bndcon ()

  do nstep = 1, nstop
    call Matrix ()
    if (ntype == 2) then
      if (kount == kprnt) then
        call Print_data (ntime, time)
        kount = 0
      end if
      time = time + dt
      kount = kount + 1
      ntime = ntime + 1
    end if
    call Resid (time)
    cold = cnew
  end do
  write (*, 10)
10 format (/, 'Maximum number of time steps or iterations reached')
end program FEM_1D
