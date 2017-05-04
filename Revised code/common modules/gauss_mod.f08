module gauss_mod
  use init_mod, only : ngaus,&
                       nnode, knode,&
                       iflg,&
                       cnew, cold,&
                       af, afm,&
                       a, f, p, g,&
                       nts, nnst,&
                       ntype, ntime,&
                       fixed,&
                       dt,&
                       Print_data,&
                       Write_data

  real, allocatable, private :: b(:)

  public  :: Gauss
  public  :: Matrix
  private :: Seidel
  private :: Gaussr
  public  :: Resid

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

  subroutine Matrix ()
    implicit none

    integer :: l, jj
    integer :: i, k
    real    :: rhs

    if (.not. allocated (b)) allocate (b(nnode))
    b = f

    if (ntype /= 1) then ! matches original intent?
      do l = 1, nnode
        do jj = 1, nnode
          rhs = (afm * a(l, jj) + p(l, jj) / dt) * cold(jj)
          b(l) = b(l) + rhs
        end do
      end do
    end if

    if (iflg == 1) then
      call Seidel ()
    else
      do k = 1, nnst
        i = nts(k)
        b(i) = 0.0
      end do
      call Gaussr (cnew, nnode)
    end if
  end subroutine Matrix

  subroutine Seidel ()
    implicit none

    real, parameter :: ER = 0.0001
    integer         :: iter
    integer         :: l, jj
    real            :: amax
    real            :: sum
    real            :: oldval
    real            :: err
    real            :: s

    iter = 0
    sum = 0.0
    amax = 1.0 + ER !Guarantees entry into while loops
    select case (ntype)
    case (1)
      do while (amax > ER)
        amax = 0.0
        do l = 1, nnode
          if (knode(l) == 1) cycle
          oldval = cnew(l)
          sum = 0.0
          do jj = 1, nnode
            if (jj == l) cycle
            sum = sum + a(l, jj) * cnew(jj)
          end do
          if (a(l,l) == 0.0) cycle ! real equality?
          cnew(l) = (-sum + f(l))/a(l,l)
          err = abs (cnew(l) - oldval)
          if (err > amax) amax = err
        end do
        iter = iter + 1
      end do
    case (2)
      do while (amax > ER)
        amax = 0.0
        do l = 1, nnode
          if (knode(l) == 1) cycle
          oldval = cnew(l)
          sum = 0.0
          do jj = 1, nnode
            if (jj == l) cycle
            sum = sum + (af * a(l, jj) + p(l, jj) / dt) * cnew(jj)
          end do
          s = af * a(l,l) + p(l,l) / dt
          cnew(l) = (-sum + b(l)) / s
          err = abs (cnew(l) - oldval)
          if (err > amax) amax = err !why this comparison? seems redundant
        end do
      end do
      ntime = iter
    case default
      print *, 'Invalid value for ntype'
      stop
    end select
  end subroutine Seidel

  subroutine Gaussr (d, n)
    implicit none

    integer, intent (in)  :: n
    real,    intent (out) :: d(:)
    integer               :: i, j, k, k1
    real                  :: r(size (b))
    real                  :: s(size (g, 1), size (g, 2))

    r = b
    s = g

!   Set up Dirichlet values
    do k = 1, nnst
      i = nts(k)
      do j = 1, n
        if (j == i) cycle
        r(j) = r(j) - s(j,i) * fixed(i)
        s(j, i) = 0.0
      end do
      do j =1, n
        if (j == i) cycle
        r(i) = r(i) - s(i,j) * fixed(i)
        s(i, j) = 0.0
      end do
      s(i, i) = 1.0
      r(i) = fixed(i)
    end do

!   Elimination routine
    do k = 1, n
      k1 = k + 1
      r(k) = r(k) / s(k,k)

      if (k == n) then
        exit
      else
        do j = k1, n
          if (s(k, j) == 0) cycle
          s(k, j) = s(k, j) / s(k, k)
          do i = k1, n
            s(i, j) = s(i, j) - s(i, k) * s(k, j)
          end do
          r(j) = r(j) - s(j, k) * r(k)
        end do
      end if
    end do

    do !Potential infinite loop
      k1 = k
      k = k - 1
      if (k == 0) exit
      do j = k1, n
        r(k) = r(k) - s(k, j) * r(j)
      end do
    end do

    d = r
  end subroutine Gaussr

  subroutine Resid (time)
    implicit none

    real, intent (in) :: time
    integer           :: maxres
    integer           :: i
    real              :: r
    real, parameter   :: ERR = 0.0001

    maxres = 0
    do i = 1, nnode
      r = abs (cnew(i) - cold(i))
      if (r > ERR) maxres = 1
    end do
    if (maxres == 1) return ! why not return from inside do?
    if (maxres == 0) then
      write (*, 20) ntime
      call Print_data (ntime, time)
      call Write_data ()
      write (*, 10)
    end if
10 format (/, 1x, 'Solution is finished')
20 format (/, 2x, 'Program has converged in ', i3, ' steps')

  stop !only reached if calculations have converged
  end subroutine Resid
end module gauss_mod
