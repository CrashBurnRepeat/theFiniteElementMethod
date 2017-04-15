module assemb_mod
  use init_mod, only : nelem,&
                       numn,&
                       x

  public  :: Assemb
  public  :: Bcside
  private :: Nodset
  private :: Shape

contains
  subroutine Assemb  ()
    implicit none

    integer :: k

    do k = 1, nelem
      call Nodset (k, i, j, m, n)
      do kk = 1, numn
        l = node (k, kk)
        do iq = 1, ngaus
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

    integer, intent (in)            :: k
    integer, intent (out)           :: i, j
    integer, intent (out), optional :: m, n

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
    end select case
  end subroutine Nodset

  subroutine Shape ()
    implicit none
    
  end subroutine Shape

  subroutine Bcside ()

  end subroutine Bcside
end module assemb_mod
