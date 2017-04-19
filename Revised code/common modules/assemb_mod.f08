module assemb_mod
  use init_mod, only : nelem,&
                       numn,&
                       x,&
                       ngaus,&
                       node, nnode, nodel,&
                       af,&
                       rhocp,&
                       p, a, g, f, h,&
                       dt, isi, isih,&
                       lme, lem,&
                       nnhc, nnqs,&
                       vx, dx, qq,&
                       q, tinf

  real, allocatable, private :: ns(:)
  real, allocatable, private :: nx(:)

  public  :: Assemb
  public  :: Bndcon
  private :: Nodset
  private :: Shape
  private :: Bcside

contains
  subroutine Assemb  (pos, w1)
    implicit none

    real, dimension (:), intent (in)  :: w1
    real, dimension (:), intent (in)  :: pos

    integer :: k, kk, kkk, iq
    integer :: i, j, m, n
    integer :: l, ll, kl, kn
    integer :: isi1, nel
    real    :: advec
    real    :: det
    real    :: mass, masst
    real    :: diff
    real    :: side
    real    :: xsi

    allocate (ns(numn))
    allocate (nx(numn))

    do k = 1, nelem
      call Nodset (k, i, j, m, n)
      do kk = 1, numn
        l = node (k, kk)
        do iq = 1, ngaus
          xsi = pos(iq)
          call Shape (xsi, i, j, m, n, det)
          f(l) = f(l) + ns(kk) * qq(k) * det * w1(iq)
          do kkk = 1, numn
            ll = node (k,kkk)
            diff = dx(k) * nx(kkk) * nx(kk)
            advec = vx(l) * nx (kkk) * ns(kk)
            mass = ns(kkk) * ns(kk)
            a(l,ll) = a(l,ll) + (diff + advec) * det * w1(iq)
            p(l,ll) = p(l,ll) + mass * det * w1(iq) * rhocp
          end do
        end do
      end do
    end do

    if (nnhc /= 0) then
      do k = 1, nnhc
        nel = lme(k)
        call Nodset(nel, i, j, m, n)
        isi1 = isih(k)
        call Bcside (isi1, side, i, j, m, n)
        kn = nodel(isi1)
        ll = node(nel,kn)
        kk = kn
        kl = ll
        masst = ns(kn) * ns(kk)
        a(ll,kl) = a(ll, kl) + masst * h(k) * side
      end do
    end if

    do l = 1, nnode
      do k = 1, nnode
        g(l,ll) = g(l,ll) + af * a(l,ll) + p(l,ll)/dt
      end do
    end do

  end subroutine Assemb

  subroutine Nodset (k, i, j, m, n)
    implicit none

    integer, intent (in)  :: k
    integer, intent (out) :: i, j
    integer, intent (out) :: m, n

    select case (numn)
    case (2)
      i = node(k, 1)
      j = node(k, 2)
    case (3)
      i = node(k, 1)
      j = node(k, 2)
      m = node(k, 3)
    case (4)
      i = node(k, 1)
      j = node(k, 2)
      m = node(k, 3)
      n = node(k, 4)
    case default
      print *, 'Invalid value for numn'
      stop
    end select
  end subroutine Nodset

  subroutine Shape (xsi, i, j, m, n, det)
    implicit none

    real,           intent (in)  :: xsi
    integer,        intent (in)  :: i, j, m, n
    real,           intent (out) :: det
    real                         :: xxsi, xlen
    real                         :: nxsi(numn)
    integer                      :: k

    select case (numn)
    case (2)
      xlen = abs (x(j) - x(i))
      ns = 0.5 * [1.0 - xsi, 1.0 + xsi]
      nxsi  = [-0.5, 0.5]
      xxsi = dot_product (nxsi, x([i, j]))
    case (3)
      xlen = abs (x(m) - x(i))
      ns = [0.5 * xsi * (xsi - 1.0),&
            1.0 - xsi ** 2,&
            0.5 * xsi * (xsi + 1.0)]
      nxsi = [xsi - 0.5, -2.0 * xsi, xsi + 0.5]
      xxsi = dot_product (nxsi, x([i, j, m]))
    case (4)
      xlen = abs (x(n) - x(i))
      ns =   [0.0625 * (1.0 - xsi) * (0.9 * xsi ** 2 - 1.0),&
              0.5625 * (1.0 - xsi **2) * (1.0 - 3.0 * xsi), &
              0.5625 * (1.0 - xsi **2) * (1.0 + 3.0 * xsi), &
              0.0625 * (1.0 + xsi) * (0.9 * xsi ** 2 - 1.0)]
      nxsi = [0.0625 * (1.0 + 0.18 * xsi - 27.0 * xsi ** 2), &
              0.5625 * (-3.0 - 2.0 * xsi + 9.0 * xsi ** 2), &
              0.5625 * (3.0 - 2.0 * xsi - 9.0 * xsi ** 2), &
              0.0625 * (-1.0 + 0.18 * xsi + 27.0 * xsi ** 2)]
      xxsi = dot_product (nxsi, x([i, j, m, n]))
    case default
      print *, ' Invalid value for numn '
      stop
    end select

    det = xxsi
    if (det == 0.0) write (*, 100) ! exact float equality? stop condition?
    do k=1, numn
      nx(k) = nxsi(k)/det
    end do
100 format (2x, 'The determinant  = 0.0')
  end subroutine Shape

  subroutine Bcside (isi1, side, i, j, m, n)
    implicit none

    integer, intent (in) :: isi1, i, j, m, n
    real, intent (out) :: side
    real :: xsi, det
    select case (isi1)
    case (1)
      xsi = -1.0
      call Shape (xsi, i, j, m, n, det)
      side = 1.0
    case (2)
      xsi = 1.0
      call Shape (xsi, i, j, m, n, det)
      side  = 1.0
    case default
      print *, 'isi1 has an invalid value'
      stop
    end select
  end subroutine Bcside

  subroutine Bndcon ()
    implicit none

    integer :: k, kn, kk, ll
    integer :: nel, isi1
    integer :: i, j, m, n
    real    :: side

    if (.not. allocated (ns)) then
      print *, 'ns is not allocated; requires an Assemb call'
      stop
    end if
    if (nnqs /= 0) then
      do k = 1, nnqs
        nel = lem(k)
        call Nodset (nel, i, j, m, n)
        isi1 = isi(k)
        call Bcside (isi1, side, i, j, m, n)
        kn = nodel(isi1)
        kk = node(nel, kn)
        f(kk) = f(kk) + ns(kn) * q(k) * side
      end do
    else if (nnhc /= 0) then
      do k = 1, nnhc
        nel = lme(k)
        call Nodset (nel, i, j, m, n)
        isi1 = isih(k)
        call Bcside (isi1, side, i, j, m, n)
        kn = nodel(isi1)
        ll = node(nel, kn)
        f(ll) = f(ll) + ns(kn) * h(k) * tinf(k) * side
      end do
    else
    end if
  end subroutine Bndcon
end module assemb_mod
