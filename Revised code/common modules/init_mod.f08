module init_mod
implicit none

  integer, private              :: infile_id, outfile_id
  integer, private              :: mtype
  integer, private              :: numdim
  integer, private              :: nvel
  integer, private              :: iaxi
  real, private                 :: to
  real, private                 :: rho
  real, private                 :: cp
  private                       :: SetPC

  integer, public               :: iflg
  integer, public               :: nnode
  integer, public               :: nelem
  integer, public               :: numn, ngaus
  integer, public               :: nstop
  integer, public               :: kprnt
  integer, public               :: ntype
  integer, public               :: ntime
  integer, public               :: nnhc = 0 ! Intended default?
  integer, public               :: nnst, nnqs
  real, public                  :: dt
  real, public                  :: af
  real, public                  :: afm
  real, public                  :: rhocp
  integer, public, dimension(2) :: isi = 0
  integer, public, dimension(2) :: lem = 0
  integer, public, dimension(2) :: lme = 0
  integer, public, dimension(2) :: isih = 0
  integer, public, dimension(2) :: nodel ! No default value?
  real, public, dimension(2)    :: h = 0.0
  real, public, dimension(2)    :: tinf = 0.0
  real, public, allocatable     :: cnew(:)
  real, public, allocatable     :: cold(:)
  integer, public, allocatable  :: nts(:)
  integer, public, allocatable  :: nqs(:)
  real, public, allocatable     :: x(:)
  real, public, allocatable     :: qq(:)
  real, public, allocatable     :: q(:)
  real, public, allocatable     :: f(:)
  real, public, allocatable     :: fixed(:)
  real, public, allocatable     :: dx(:)
  real, public, allocatable     :: vx(:)
  integer, public, allocatable  :: knode(:)
  integer, public, allocatable  :: node(:,:)
  real, public, allocatable     :: p(:,:)
  real, public, allocatable     :: a(:,:)
  real, public, allocatable     :: g(:,:)
  public :: Init
  public :: Print_data
  public :: Write_data

contains
  subroutine Init ()
    implicit none

    character (len=18)            :: title
    real                          :: xmin
    real                          :: xmax
    real                          :: cnt
    integer                       :: i, j, mm, k
    integer                       :: nt
    integer                       :: isi1
    logical                       :: exit_flg = .false.
    character (len=4)             :: word
    character (len=4), parameter  :: STOP_KEYS(3) = &
                                    ['DIRC','FLUX','CONV']

    call SetPC ()
    write (*, 100)
100 format (80('*'))
    read (infile_id, 101) title
101 format (3x,a18)
    write (*, 102) title
102 format (5x, a18)

    read (infile_id, *) mtype,&
                        numdim,&
                        nnode,&
                        nelem,&
                        numn,&
                        nstop,&
                        kprnt,&
                        nvel
    read (infile_id, *) ntype,&
                        dt,&
                        af,&
                        to,&
                        rho,&
                        cp,&
                        iaxi
    mtype = mtype + 1
    if (ntype==0) then
      af = 0.0
      afm = 0.0
    end if
    ngaus = numn
    rhocp = rho * cp

    allocate (cnew(nnode))
    allocate (cold(nnode))
    allocate (nts(nnode))
    allocate (nqs(nnode))
    allocate (x(nnode))
    allocate (qq(nnode))
    allocate (q(nnode))
    allocate (f(nnode))
    allocate (fixed(nnode))
    allocate (dx(nnode))
    allocate (vx(nnode))
    allocate (knode(nnode))

    cnew = 0.0
    cold = to
    nts = 0
    nqs = 0
    qq = 0.0
    q = 0.0
    f = 0.0
    fixed = 0.0
    dx = 0.0
    vx = 0.0
    knode = 0

    allocate (node(nnode,4))
    node = 0

    allocate (p(nnode,nnode))
    allocate (a(nnode,nnode))
    allocate (g(nnode,nnode))

    p = 0.0
    a = 0.0
    g = 0.0

    read (infile_id, *) (i,x(i),j=1,nnode) ! index values?

    xmin = minval (x) ! valid simplification?
    xmax = maxval (x) ! valid simplification?

    do i=1, nelem
      read (infile_id, *) j,&
                          qq(j),&
                          dx(j),&
                          (node(j,mm), mm=1, numn)
    end do

    do i=1, nnode
      read (infile_id, 1015) word, nt, cnt
      if (word==STOP_KEYS(1)) exit
      nts(i) = nt
      cold(nt) = cnt
      knode(nt) = 1
      fixed(nt) = cnt
    end do
    nnst = i - 1

    do i = 1, nnode
      read (infile_id, 1025) word, q(i), lem(i), isi(i)
      if (word == STOP_KEYS(2)) exit
    end do
    nnqs = i - 1

    k = 0
    do i = 1, nnode
      do j = 1, nnst
        if (i==nts(j)) then
          exit_flg = .true.
          exit
        end if
      end do
      if (exit_flg) exit
      k = k + 1
      nqs(k) = i
    end do

    !do i = 1, 1
    i = 1
    read (infile_id, *) h(i), tinf(i), lme(i), isih(i)
    !end do

    nnhc = 1
    cnew = cold ! intended assignment?

    if (nvel==1) then
      read (infile_id, *) (i, vx(i), j=1,nnode)
    end if

    nodel(1) = 1
    select case (numn)
    case (2)
      nodel(2) = 2
    case (3)
      nodel(2) = 3
    case (4)
      nodel(2) = 4
    case default
      print *, ' Error in nodel value'
      stop
    end select

    if (ntype == 1) write (*, 14) ntype
    if (ntype == 2) write (*, 15) ntype
    afm = af - 1.0
    if (ntype /= 1) then
      if (af == 1.0) write (*, 10) dt ! real equality?
      if (af == 0.5) write (*, 11) dt ! real equality?
    end if

    write (*, 12) nstop, kprnt, to, rho, cp
    write (*, 1035) nnode, nelem
    write (*, 1040)
    write (*, 1041)
    do i = 1, nnode
      write (*, 1045) i, x(i)
    end do

    write (*, 1050)
    do i = 1, nelem
      write (*, 1055) i, qq(i), dx(i), (node(i,mm), mm=1, numn)
    end do

    write (*, 1060)
    write (*, 1061)
    do i = 1, nnst
      write (*, 1065) i, nts(i), cold(nts(i))
    end do

    write (*, 1080)
    write (*, 1082)
    do i = 1, nnqs
      isi1 = isi(i)
      write (*, 1085) i, q(i), node(lem(i), nodel(isi1))
    end do

    write (*, 1081)
    write (*, 1083)
    do i = 1, nnhc
      isi1=isih(i)
      write (*, 1086) i, h(i), tinf(i), node(lme(i), nodel(isi1))
    end do

    if (nvel /= 0) then
      write (*, 4040)
      write (*, 4041)
      do i = 1, nnode
        write (*, 1045) i, vx(i)
      end do
    end if

  10 format (10x, 'Fully implicit ', 2x, 'dt = ', f8.4)
  11 format (10x, 'Crank Nicolson ', 2x, 'dt = ', f8.4)
  12 format (5x, 'nstop = ', i5, 2x, 'kprnt = ', i5, 2x,&
            'to = ', f8.4, 2x, 'rho = ', f8.4, 2x, 'cp = ', f8.4)
  14 format (/, 10x, 'ntype = ', i3, 2x, 'steady state')
  15 format (/, 10x, 'ntype = ', i3, 2x, 'time dependent')
1015 format (6x, a4, i5, 5x, f10.5)
1025 format (6x, a4, f10.5, 2i5)
1026 format (6x, a4, 2f10.5, 2i5)
1035 format (/, 10x, 'No. of nodes = ', i3,&
            5x, 'No. of elements = ', i3)
1040 format (/, 5x, 'Summary of nodal coordinates')
1041 format (7x, 'I', 12x, 'X(I)')
1045 format (5x, i3, 7x, f10.3)
1050 format (/, 1x, 'Element number', 4x, 'Source',&
            5x, 'dx', 8x, 'Node numbers')
1055 format (5x, i3, 5x, 2f10.3, 4(4x, i3))
1060 format (/, 7x, 'Nodes were cnew is specified')
1061 format (8x, 'I', 4x, 'Node', 6x, 'cnew')
1065 format (2x, 2(4x, i3), 3x, f8.3)
1080 format (/, 7x, 'Nodes were flux is specified')
1081 format (/, 7x, 'Nodes were tinf is specified')
1082 format (8x, 'I', 4x, 'Flux(q)', 3x, 'Node numbers')
1083 format (8x, 'I', 6x, 'h', 6x, 'tinf', 4x, 'Node numbers')
1085 format (6x, i3, 2x, f8.2, 6x, i3)
1086 format (6x, i3, 2x, 2f8.2, 4x, i3)
4040 format (/, 5x, 'Summary of specified nodal velocities')
4041 format (7x, 'I', 12x, 'vx(i)')

  end subroutine Init

  subroutine SetPC ()
    implicit none

    character (len=12)    :: infile, outfile

    print *, ' Enter name for input file: '
    read *, infile
    open (newunit=infile_id, file=infile, status='old')

    print *, ' Enter name for plot output file: '
    read *, outfile
    open (newunit=outfile_id, file=outfile, status='new')
    write (outfile_id, '(a)') outfile

    print *, ' Enter 1 for Gauss-Seidel iteration or&
            & 2 for Gauss eliminiation'
    read (*, '(i1)') iflg

  end subroutine SetPC

  subroutine Print_data (ntime, time)
    implicit none

    integer, intent (in) :: ntime
    real,    intent (in) :: time
    integer              :: i

    write (*, 120)
    if (ntype == 1) then
      write (*, 403)
    else
      write (*, 402) ntime, time
    end if

    write (*, 400)
    write (*, 401) (i, cnew(i), i=1, nnode)

120 format (/, 20x, 'Variable values')
400 format (/, 5x, 'node', 6x, 'cnew')
401 format (5x, i3, 5x, f8.2)
402 format (/, 5x, 'Time steps', i4, 2x, 'Time = ', f8.3)
403 format (5x, 'Steady state')
  end subroutine Print_data

  subroutine Write_data ()
    implicit none

    integer :: i, j, mm

    write (outfile_id, '(8i4)') mtype,&
                                numdim,&
                                nnode,&
                                nelem,&
                                numn,&
                                nstop,&
                                kprnt,&
                                nvel
    write (outfile_id, '(i4, 5(f8.4, 1x), 2x i4)')&
                                ntype,&
                                dt,&
                                af,&
                                to,&
                                rho,&
                                cp,&
                                iaxi
    do i = 1, nnode
      write (outfile_id, 100) i, x(i)
    end do

    do j = 1, nelem
      write (outfile_id, 101) j, qq(j), dx(j),&
                              (node(j, mm), mm=1, numn)
    end do

    do i = 1, nnode
      write (outfile_id, 100) i, cnew(i), vx(i)
    end do

100 format (5x, i3, 2x, 2(f8.3, 1x))
101 format (5x, i3, 2x, 2(f8.3, 1x), 2x, 4i3)
  end subroutine Write_data
end module init_mod
