!> Public module containing utility subroutines
module utils
  use constants
  implicit none
  !private !> makes default of everything defined here private, use word public to counter.
  !if not used, can define a single private variable as well (everything else is public)
  logical,private :: available_outside = .FALSE. !> This variable is not available outside of the module
  logical,private  :: reuse_normal_rand=.FALSE. !> used by normal_rand function
  real(wp),private :: stored_normal_rand !> used by normal_rand function

  contains
   !> TRUE if the real number is indistinguishable from another real number
  pure function are_equal(number1,number2)
  logical :: are_equal
  real(wp),intent(in) :: number1, number2
  are_equal = is_zero(abs(number1-number2))
end function are_equal

!> TRUE if the real number is indistinguishable from 0
pure function is_zero(number)
  logical :: is_zero
  real(wp),intent(in) :: number
  is_zero = number < epsilon(0.0_wp)
end function is_zero

!> If all integers are equal, return one Otherwise, terminate with error string
function assert_eq(n,string)
  character(len=*), intent(in) :: string !> Error message
  integer, intent(in) :: n(:) !> Array of integers
  integer :: assert_eq
  if ( all(n(2:) == n(1)) ) then 
     assert_eq=n(1)
  else
     write (0,*) "Two integers are different @ ", string
     STOP "Terminated by assert_eq"
  end if
end function assert_eq

!> Calculate mean and standard deviation of a flat array
pure function mean_std(array)
  real(wp) :: mean_std(2)
  real(wp),intent(in) :: array(:)

  mean_std(1) = sum(array(:)) / size(array) ! mean
  mean_std(2) = sum( (array(:) - mean_std(1))**2 ) / size(array) ! variance
  mean_std(2) = sqrt( mean_std(2) / (size(array)-1) )
end function mean_std
!> Solves for a vector u of size N the tridiagonal linear set given by equation
!! Input vectors b (diagonal elements) and r (right-hand sides) have size N ,
!! while a and c (off-diagonal elements) are size N − 1
RECURSIVE SUBROUTINE tridag(a,b,c,r,u)
  IMPLICIT NONE
  REAL(wp), DIMENSION(:), INTENT(IN) :: a,b,c,r
  REAL(wp), DIMENSION(:), INTENT(OUT) :: u
  INTEGER, PARAMETER :: NPAR_TRIDAG=4 !< Serial algorithm is invoked below this value
  INTEGER :: n,n2,nm,nx
  REAL(wp), DIMENSION(size(b)/2) :: y,q,piva
  REAL(wp), DIMENSION(size(b)/2-1) :: x,z
  REAL(wp), DIMENSION(size(a)/2) :: pivc
  n=size(b)
  ! n=assert_eq((/size(a)+1,size(b),size(c)+1,size(r),size(u)/),'tridag')
  if (n < NPAR_TRIDAG) then
  call tridag_ser(a,b,c,r,u)
  else
  if ( maxval(abs(b(1:n))) < epsilon(0.0_wp) ) STOP "tridag: possible singular matrix"
  n2=size(y)
  nm=size(pivc)
  nx=size(x)
  piva = a(1:n-1:2)/b(1:n-1:2) !Zero the odd a’s and even c’s, giving x, y, z, q
  pivc = c(2:n-1:2)/b(3:n:2)
  y(1:nm) = b(2:n-1:2)-piva(1:nm)*c(1:n-2:2)-pivc*a(2:n-1:2)
  q(1:nm) = r(2:n-1:2)-piva(1:nm)*r(1:n-2:2)-pivc*r(3:n:2)
  if (nm < n2) then
  y(n2) = b(n)-piva(n2)*c(n-1)
  q(n2) = r(n)-piva(n2)*r(n-1)
  end if
  x = -piva(2:n2)*a(2:n-2:2)
  z = -pivc(1:nx)*c(3:n-1:2)
  call tridag(x,y,z,q,u(2:n:2)) ! Recurse and get even u’s
  u(1) = (r(1)-c(1)*u(2))/b(1) ! Substitute and get odd u’s
  u(3:n-1:2) = (r(3:n-1:2)-a(2:n-2:2)*u(2:n-2:2) &
  -c(3:n-1:2)*u(4:n:2))/b(3:n-1:2)
  if (nm == n2) u(n)=(r(n)-a(n-1)*u(n-1))/b(n)
  end if

  contains
  !> Serial tridiagonal solver used by recursive tridag once the problem becomes small enough 
  SUBROUTINE tridag_ser(aa,bb,cc,rr,uu)
     IMPLICIT NONE
     REAL(wp), DIMENSION(:), INTENT(IN) :: aa,bb,cc,rr
     REAL(wp), DIMENSION(:), INTENT(OUT) :: uu
     REAL(wp), DIMENSION(size(b)) :: gam !! One vector of workspace, gam is needed
     INTEGER :: nn,j
     REAL(wp) :: bet
     nn = size(uu)
     ! nn=assert_eq((/size(aa)+1,size(bb),size(cc)+1,size(rr),size(uu)/),"tridag_ser")
     bet=bb(1)
     if ( bet < epsilon(0.0_wp) ) STOP "tridag_ser: Error at code stage &
     & If this happens then you should rewrite your equations &
     & as a set of order N - 1, with u2 trivially eliminated"
     uu(1)=rr(1)/bet
     do j=2,nn ! Decomposition and forward substitution
     gam(j)=cc(j-1)/bet
     bet=bb(j)-aa(j-1)*gam(j)
     if ( bet < epsilon(0.0_wp) ) STOP "tridag_ser: Error at code stage 2"
     uu(j)=(rr(j)-aa(j-1)*uu(j-1))/bet
     end do
     do j=nn-1,1,-1 ! Backsubstitution
     uu(j)=uu(j)-gam(j+1)*uu(j+1)
     end do
  END SUBROUTINE tridag_ser

END SUBROUTINE tridag

!> Places x into ordered array
!! @ detail Given an array xx(1:N ), and given a value x, returns a value j such that x is between
!!  xx(j) and xx(j + 1) xx must be monotonic, either increasing or decreasing j = 0 or
!!  j = N is returned to indicate that x is out of range
pure FUNCTION locate(xx,x)
  IMPLICIT NONE
  REAL(wp), DIMENSION(:), INTENT(IN) :: xx
  REAL(wp), INTENT(IN) :: x
  INTEGER :: locate
  
  INTEGER :: n,jl,jm,ju
  LOGICAL :: ascnd
  n=size(xx)
  ascnd = (xx(n) >= xx(1)) ! True if ascending order of table, false otherwise
  jl=0   ! Initialize lower
  ju=n+1 ! and upper limits
  do
     if (ju-jl <= 1) exit ! Repeat until this condition is satisfied
     jm=(ju+jl)/2 ! Compute a midpoint,
     if (ascnd .eqv. (x >= xx(jm))) then
        jl=jm ! and replace either the lower limit
     else
        ju=jm ! or the upper limit, as appropriate
     end if
  end do
  if ( abs(x - xx(1)) < epsilon(0.0_wp) ) then ! set the output, being careful with the endpoints
     locate=1
  else if ( abs(x - xx(n)) < epsilon(0.0_wp) ) then
     locate=n-1
  else
     locate=jl
  end if
END FUNCTION locate

!> Given 2 3d-vectors, compute their cross produnct
pure function cross_prod(x,y)
  real(wp) :: cross_prod(3)
  real(wp),intent(in),dimension(3) :: x, y
  cross_prod(1) = x(2)*y(3) - x(3)*y(2)
  cross_prod(2) =-x(1)*y(3) + x(3)*y(1)
  cross_prod(3) = x(1)*y(2) - x(2)*y(1)
end function cross_prod

!> Generate normally-distibuted (0 mean, 1 variance) random number using the Box-Muller Method
function normal_rand()
  real(wp) :: normal_rand
  real(wp) :: rnum1,rnum2 !> Random numbers on [0,1]
  real(wp) :: rsq !> distance from origin squared

  if (reuse_normal_rand) then
     normal_rand = stored_normal_rand
     reuse_normal_rand = .FALSE.
  else
     do
        ! 2 uniform random numbers on the square from -1 to 1
        call random_number(rnum1)
        call random_number(rnum2)
        rnum1 = 2*rnum1 - 1
        rnum2 = 2*rnum2 - 1
        rsq = rnum1**2 + rnum2**2
        if ( rsq < 1.0_wp ) exit ! repeat until they are inside the unit circle
     end do
     rsq = sqrt(-2 * log(rsq)/rsq)
     reuse_normal_rand  = .TRUE.
     stored_normal_rand = rnum1 * rsq
     normal_rand        = rnum2 * rsq
  end if
end function normal_rand

!> Initialize random seed from /dev/urandom. If that fails, use XOR of time and PID
subroutine init_random_seed()
  use iso_fortran_env, only: int64
  integer, allocatable :: seed(:)
  integer :: i, n, un, istat, dt(8), pid,getpid
  integer(int64) :: t

  call random_seed(size = n)
  allocate(seed(n))
  ! First try if the OS provides a random number generator
  open(newunit=un, file="/dev/urandom", access="stream", &
     form="unformatted", action="read", status="old",&
     iostat=istat)
  if (istat == 0) then
     read(un) seed
     close(un)
  else
     ! Fallback to XOR:ing the current time and pid. The PID
     ! is useful in case one launches multiple instances of the
     ! same program in parallel.
     call system_clock(t)
     if (t == 0) then
        call date_and_time(values=dt)
        t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
           + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
           + dt(3) * 24_int64 * 60 * 60 * 1000 &
           + dt(5) * 60 * 60 * 1000 &
           + dt(6) * 60 * 1000 + dt(7) * 1000 &
           + dt(8)
     end if
     pid = getpid()
     t = ieor(t, int(pid, kind(t)))
     do i = 1, n
        seed(i) = lcg(t)
     end do
  end if
  write(6,*)
  write(6,*) "The Seed is ", seed
  write(6,*)
  call random_seed(put=seed)
  deallocate(seed)
contains
  ! This simple PRNG might not be good enough for real work,
  ! but is sufficient for seeding a better PRNG.
  function lcg(s)
     integer :: lcg
     integer(int64) :: s
     if (s == 0) then
        s = 104729
     else
        s = mod(s, 4294967296_int64)
     end if
     s = mod(s * 279470273_int64, 4294967291_int64)
     lcg = int(mod(s, int(huge(0), int64)), kind(0))
  end function lcg
end subroutine init_random_seed

SUBROUTINE spline(x,y,yp1,ypn,y2)
  !> @brief Returns second derivatives for interpolating function @see splint
  !> @details Given arrays x and y of length N containing a tabulated function, i.e., y_i = f(xi), with x1 <
  !! x2 < ... < xN , and given values yp1 and ypn for the first derivative of the interpolating
  !! function at points 1 and N, respectively, this routine returns an array y2 of length N
  !! that contains the second derivatives of the interpolating function at the tabulated points
  !! xi If yp1 and/or ypn are equal to 1 × 1030 or larger, the routine is signaled to set the
  !!corresponding boundary condition for a natural spline, with zero second derivative on that
  !! boundary

  real(wp), intent(in) :: ypn !> first derivative at the last point
  real(wp), intent(in) :: yp1 !> first derivative at the first point
  real(wp), dimension(:), intent(in)  :: x !> array of data for fit y = f(x)
  real(wp), dimension(:), intent(in)  :: y !> array of data for fit y = f(x)
  real(wp), dimension(:), intent(out) :: y2 !> second derivatives of the interpolating function at xi

  integer :: n
  real(wp), dimension(size(x)) :: a,b,c,r
  n = size(x)
!   n = assert_eq((/size(x),size(y),size(y2)/),'spline')
  c(1:n-1)=x(2:n)-x(1:n-1) ! Set up the tridiagonal equations
  r(1:n-1)=6*((y(2:n)-y(1:n-1))/c(1:n-1))
  r(2:n-1)=r(2:n-1)-r(1:n-2)
  a(2:n-1)=c(1:n-2)
  b(2:n-1)=2*(c(2:n-1)+a(2:n-1))
  b(1)=1; b(n)=1
  if (yp1 > 0.99e30_wp) then !The lower boundary condition is set either to be “natural”
     r(1)=0; c(1)=0
  else ! or else to have a specified first derivative
     r(1)=(3/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
     c(1)=0.5_wp
  end if ! The upper boundary condition is set either to be “natural”
  if (ypn > 0.99e30_wp) then
     r(n)=0; a(n)=0
  else ! or else to have a specified first derivative
     r(n)=(-3/(x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)
     a(n)=0.5_wp
  end if
  call tridag(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n))
END SUBROUTINE spline

!> Return cubic spline interpolated value at x
function splint(xa,ya,y2a,x)
  REAL(wp), INTENT(IN) :: xa(:)  !> tabulate x [same as for spline]
  REAL(wp), INTENT(IN) :: ya(:)  !> tabulate y=f(x) [same as for spline]
  REAL(wp), INTENT(IN) :: y2a(:) !> second derivatives [output of spline]
  REAL(wp), INTENT(IN) :: x !> evaluate the function here
  REAL(wp) :: splint
  
  INTEGER :: khi,klo,n
  REAL(wp) :: a,b,h
  n = size(xa)
!   n=assert_eq((/size(xa),size(ya),size(y2a)/),'splint')
  klo=max(min(locate(xa,x),n-1),1) ! bisection (log N)
  khi=klo+1 ! klo and khi now bracket the input value of x
  h=xa(khi)-xa(klo)
  if (is_zero(h)) STOP "ERRROR: Bad xa input in splint The xa's must be distinct."
  a=(xa(khi)-x)/h ! Cubic spline polynomial is now evaluated
  b=(x-xa(klo))/h
  splint=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6
end function splint

!> Computes sqrt(a**2 + b**2) without destructive underflow or overflow
pure FUNCTION pythag(a,b)
  REAL(WP), INTENT(IN) :: a,b
  REAL(WP) :: pythag
  REAL(WP) :: absa,absb
  absa = abs(a); absb = abs(b)
  if (absa > absb) then
     pythag = absa * sqrt(1 + (absb/absa)**2)
  else
     if (is_zero(absb)) then
        pythag = 0
     else
        pythag = absb*sqrt(1 + (absa/absb)**2)
     end if
  end if
END FUNCTION pythag

!> !DESTROYS a! Return sorted eigenvalues of real symmetric a 
function eigvals(a)
  real(wp),allocatable :: eigvals(:) !> on-diagonal of the tridiagonal matrix and eigenvalues
  real(wp),intent(inout) :: a(:,:)
  ! real(wp),allocatable :: aa(:,:) !> Copy of matrix a
  real(wp),allocatable :: e(:) !> off-diagonal of the tridiagonal matrix
  integer :: n
  n = size(a,dim=1)
  allocate(eigvals(n),e(n))!,aa(n,n))
  ! aa(:,:) = a(:,:)
  call tred2(a,eigvals,e,novectors=.true.)
  call tqli(eigvals,e)
  call eigsrt(eigvals)
end function eigvals
!> Sort the eigenvalues into descending order and rearranges the columns of v
!! Given the eigenvalues d and (optionally) eigenvectors v as output tqli
SUBROUTINE eigsrt(d,v)
  ! USE nrtype; USE nrutil, ONLY :assert_eq,imaxloc,swap
  REAL(wp), INTENT(INOUT) :: d(:)
  REAL(WP), INTENT(INOUT),optional :: v(:,:)
  real(wp) :: temp
  INTEGER :: i,j,n,ii
  n=size(d)
  do i=1,n-1
     j = sum(minloc(d(i:n))+i-1)
     if (j /= i) then ! SWAP eigenvalues and eigenvectors
        temp = d(i); d(i) = d(j); d(j) = temp
        if (present(v)) then
           do ii=1,n
              temp = v(ii,i); v(ii,i) = v(ii,j); v(ii,j) = temp
           end do
        end if
     end if
  end do
END SUBROUTINE eigsrt

!> eigenvalues and eigenvectors of a real symmetric tridiagonal matrix 
!! (or real, symmetric matrix previously reduced by tred2)
!! Uses QL algorithm with implicit shifts
SUBROUTINE tqli(d,e,z)
  ! USE nr, ONLY :pythag
  REAL(wp), INTENT(INOUT) :: d(:) !> [N] IN: matrix diagonal. OUT: eigenvalues
  REAL(wp), INTENT(INOUT) :: e(:) !> [N] IN: subdiagonal elements with e(1) arbitrary. OUT: destroyed
  REAL(wp), OPTIONAL, INTENT(INOUT) :: z(:,:) !> [N,N] If the eigenvectors are desired, 
  !! IN: the eigenvectors of a matrix that has been reduced by tred2 OR z=Identity (if matrix already tridiagonal)
  !! OUT: kth column of z returns the normalized eigenvector corresponding to d(k).
  INTEGER  :: i,iter,l,m,n,ndum
  REAL(wp) :: b,c,dd,f,g,p,r,s
  REAL(wp) :: ff(size(e))
  n = size(d)
  if (present(z)) ndum = n
  e(:) = eoshift(e(:),1) ! Convenient to renumber the elements of e
  do l=1,n
     iter=0
     iterate: do
        do m=l,n-1 ! Look for a single small subdiagonal element to split the matrix
           dd=abs(d(m))+abs(d(m+1))
           if (abs(e(m))+dd == dd) exit
        end do
        if (m == l) exit iterate
        if (iter == 300) exit ! too many iterations
        iter = iter+1
        g = (d(l+1)-d(l))/(2*e(l)) ! Form shift
        r = pythag(g,1.0_wp)
        g = d(m)-d(l)+e(l)/(g+sign(r,g)) !This is d_m − k_s
        s=1; c=1; p=0
        do i=m-1,l,-1 ! A plane rotation as in the original QL,
           ! followed by Givens rotations to restore tridiagonal form.
           f=s*e(i); b=c*e(i); r=pythag(f,g); e(i+1)=r
           if (is_zero(r)) then ! Recover from underflow
              d(i+1)=d(i+1)-p; e(m)=0
              cycle iterate
           end if
           s=f/r; c=g/r; g=d(i+1)-p
           r=(d(i)-g)*s + 2*c*b
           p=s*r; d(i+1)=g+p; g=c*r-b
           if (present(z)) then ! Form eigenvectors.
              ff(1:n)=z(1:n,i+1)
              z(1:n,i+1)=s*z(1:n,i)+c*ff(1:n)
              z(1:n,i)=c*z(1:n,i)-s*ff(1:n)
           end if
        end do
        d(l)=d(l)-p; e(l)=g; e(m)=0
     end do iterate
  end do
END SUBROUTINE tqli

!> m x n outer product matrix of vector a(m) * b(n)
pure function outerprod(a,b)
  real(wp),allocatable :: outerprod(:,:)
  real(wp),intent(in),dimension(:) :: a,b !> [N], [M] input vectors
  allocate(outerprod(size(a),size(b)))
  outerprod = spread(a,dim=2,ncopies=size(b)) * spread(b,dim=1,ncopies=size(a))
end function outerprod

!> Reduction of real symmetric matrix a to tridiagonal form d[N], e[N].
SUBROUTINE tred2(a,d,e,novectors)
REAL(wp),  INTENT(INOUT) :: a(:,:) !> [N,N] IN: real, symmetric matrix a; OUT: orthogonal matrix Q effecting the transformation
REAL(wp), INTENT(OUT) :: d(:) !> [N] diagonal elements
real(wp), intent(out) :: e(:) !> [N] off-diagonal elements, with e(1)=0
LOGICAL, OPTIONAL, INTENT(IN) :: novectors !> If true, only eigenvecors are needed next, so a is junk
INTEGER :: i,j,l,n
REAL(wp) :: f,g,h,hh,scale
REAL(wp) :: gg(size(a,1))
LOGICAL, SAVE :: yesvec=.true.
n=size(a,1)
if (present(novectors)) yesvec = .not. novectors
do i=n,2,-1
  l=i-1; h=0
  if (l > 1) then
     scale=sum(abs(a(i,1:l)))
     if (is_zero(scale)) then ! Skip transformation
        e(i)=a(i,l)
     else
        a(i,1:l)=a(i,1:l)/scale ! Use scaled a’s for transformation
        h=sum(a(i,1:l)**2) ! Form σ in h
        f=a(i,l)
        g=-sign(sqrt(h),f)
        e(i)=scale*g
        h=h-f*g ! Now h is equation (11.2.4)
        a(i,l)=f-g ! Store u in the ith row of a
        if (yesvec) a(1:l,i)=a(i,1:l)/h ! Store u/H in ith column of a
        do j=1,l !Store elements of p in temporarily unused elements of e
           e(j) = (dot_product(a(j,1:j),a(i,1:j)) &
                 +dot_product(a(j+1:l,j),a(i,j+1:l))) / h
        end do
        f=dot_product(e(1:l),a(i,1:l))
        hh=f/(h+h) !Form K, equation (11.2.11)
        e(1:l)=e(1:l)-hh*a(i,1:l) !Form q and store in e overwriting p
        do j=1,l !Reduce a, equation (11.2.13)
           a(j,1:j)=a(j,1:j)-a(i,j)*e(1:j)-e(j)*a(i,1:j)
        end do
     end if
  else
     e(i)=a(i,l)
  end if
  d(i)=h
end do
if (yesvec) d(1)=0
e(1)=0
do i=1,n !Begin accumulation of transformation matrices
  if (yesvec) then
  l=i-1
  if (.not. is_zero(d(i))) then !This block skipped when i=1. Use u and u/H stored in a to form P·Q
     gg(1:l)=matmul(a(i,1:l),a(1:l,1:l))
     a(1:l,1:l)=a(1:l,1:l)-outerprod(a(1:l,i),gg(1:l))
  end if
  d(i)=a(i,i) !Reset row and column of a to identity matrix for next iteration
  a(i,i)=1; a(i,1:l)=0; a(1:l,i)=0
else
  d(i)=a(i,i)
end if
end do
END SUBROUTINE tred2

!    !> Return d=Eigenvalues and eigenvectors of a real symmetric a(N,N). On output,
!    elements of a above the diagonal are destroyed. d is a vector of length N that returns the
!    eigenvalues of a. v is an N × N matrix whose columns contain, on output, the . nrot returns.

!    SUBROUTINE jacobi(a,d,v,nrot)
!       ! get_diag,nrerror,unit_matrix,& upper_triangle
!       INTEGER,  INTENT(OUT) :: nrot !> the number of Jacobi rotations that were required
!       REAL(wp), INTENT(OUT) :: d(:) !> [N] Eigenvalues
!       REAL(wp), INTENT(INOUT) :: a(:,:)
!       REAL(wp), INTENT(OUT) :: v(:,:) !> [N,N] normalized eigenvectors of a
!       INTEGER(I4B) :: i,ip,iq,n
!       REAL(wp) :: c,g,h,s,sm,t,tau,theta,tresh
!       REAL(wp), DIMENSION(size(d)) :: b,z
!       n=assert_eq((/size(a,1),size(a,2),size(d),size(v,1),size(v,2)/),’jacobi’)
!       call unit_matrix(v(:,:)) Initialize v to the identity matrix.
!       b(:)=get_diag(a(:,:)) Initialize b and d to the diagonal of
!       a.d(:)=b(:)
!       z(:)=0.0 This vector will accumulate terms of
!       the form tapq as in eq. (11.1.14).nrot=0
!       do i=1,50
!       sm=sum(abs(a),mask=upper_triangle(n,n)) Sum off-diagonal elements.
!       if (sm == 0.0) RETURN
!       The normal return, which relies on quadratic convergence to machine underflow.
!       tresh=merge(0.2_wp*sm/n**2,0.0_wp, i < 4 )
!       On the first three sweeps, we will rotate only if tresh exceeded.
!       do ip=1,n-1
!       do iq=ip+1,n
!       g=100*abs(a(ip,iq))
!       After four sweeps, skip the rotation if the off-diagonal element is small.
!       if ((i > 4) .and. (abs(d(ip))+g == abs(d(ip))) &
!       .and. (abs(d(iq))+g == abs(d(iq)))) then
!       a(ip,iq)=0.0
!       else if (abs(a(ip,iq)) > tresh) then
!       h=d(iq)-d(ip)
!       if (abs(h)+g == abs(h)) then
!       t=a(ip,iq)/h t = 1/(2θ)
!       else
!       theta=h/a(ip,iq)/2 Equation (11.1.10).
!       t=1/(abs(theta)+sqrt(1+theta**2))
!       if (theta < 0.0) t=-t
!       end if
!       c=1/sqrt(1+t**2)
!       s=t*c
!       tau=s/(1+c)
!       h=t*a(ip,iq)
!       z(ip)=z(ip)-h
!       z(iq)=z(iq)+h
!       d(ip)=d(ip)-h
!       d(iq)=d(iq)+h
!       a(ip,iq)=0.0
!       call jrotate(a(1:ip-1,ip),a(1:ip-1,iq))
! Case of rotations 1 ≤ j < p.
! call jrotate(a(ip,ip+1:iq-1),a(ip+1:iq-1,iq))
! Case of rotations p < j < q.
! call jrotate(a(ip,iq+1:n),a(iq,iq+1:n))
! Case of rotations q < j ≤ n.
! call jrotate(v(:,ip),v(:,iq))
! nrot=nrot+1
! end if
! end do
! end do
! b(:)=b(:)+z(:)
! d(:)=b(:) Update d with the sum of tapq ,
! z(:)=0.0 and reinitialize z.
! end do
! call nrerror(’too many iterations in jacobi’)
! CONTAINS
! SUBROUTINE jrotate(a1,a2)
! REAL(wp), DIMENSION(:), INTENT(INOUT) :: a1,a2
! REAL(wp), DIMENSION(size(a1)) :: wk1
! wk1(:)=a1(:)
! a1(:)=a1(:)-s*(a2(:)+a1(:)*tau)
! a2(:)=a2(:)+s*(wk1(:)-a2(:)*tau)
! END SUBROUTINE jrotate
! END SUBROUTINE jacobi
end module utils