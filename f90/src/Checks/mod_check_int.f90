!-- Integrability scan of an integrand (-checkint).
!-- Sends subsets of the Vegas variables to the edge, x = d or 1-d, and prints the median of
!-- |f| d^k (k = number of variables at the edge) for three depths d. The integrand is
!-- integrable if this falls with d; it stays constant if a counterterm is missing or wrong.
!-- Combined limits use shallower depths (the amplitudes lose precision when d^k is tiny).
!-- The test points are the last points with f /= 0 of a Vegas warm-up (rec_integrand).
module mod_check_int
  use mod_types
  use mod_histo, only: nohistos
  use mod_kinematics_nlo, only: NLOdim_vegas, xE, xRHO
  use mod_kinematics_nnlo_tc_is_eu_ac, only: NNLOdim_vegas, xE5, xE6, xRHO5, xRHO6
  use mod_proc_parms, only: override, buff, buff_rr, buff_r, set_buff
  use mod_limvals, only: limval_nnlo, limval_onlo
  use mod_auxfunctions, only: mypause
  implicit none
  private

  public :: check_integrability, start_recording, rec_integrand, check_limits_rr, check_limits_onlo

  abstract interface
     function integrand_if(xx,ff,weight)
       use mod_types
       integer    :: integrand_if
       real(dp15) :: xx(30),ff(1),weight
     end function integrand_if
  end interface
  procedure(integrand_if), pointer, save :: target_f => null()

  integer, parameter  :: npts = 100, ndep = 3
  real(dp15), save :: pts(30,npts) = 0.5_dp15
  integer, save    :: nrec = 0, ndim_rec = 0
  !-- depth(:,k) for k variables at the edge
  real(dp15), parameter :: depth(ndep,4) = reshape([1e-2_dp15,1e-4_dp15,1e-6_dp15, 1e-3_dp15,1e-4_dp15,1e-5_dp15, &
       1e-2_dp15,3e-3_dp15,1e-3_dp15, 3e-2_dp15,1e-2_dp15,3e-3_dp15],[ndep,4])

contains

  subroutine start_recording(integrand,ndim)
    procedure(integrand_if) :: integrand
    integer, intent(in) :: ndim
    target_f => integrand
    ndim_rec = ndim
    nrec = 0
  end subroutine start_recording

  !-- the integrand, keeping the last npts points where it does not vanish
  function rec_integrand(xx,ff,weight)
    integer    :: rec_integrand
    real(dp15) :: xx(30),ff(1),weight
    rec_integrand = target_f(xx,ff,weight)
    if (ff(1) /= 0 .and. ff(1) == ff(1)) then
       nrec = nrec + 1
       pts(1:ndim_rec,mod(nrec-1,npts)+1) = xx(1:ndim_rec)
    endif
  end function rec_integrand

  subroutine check_integrability(ndim,integrand)
    integer, intent(in) :: ndim
    interface
       function integrand(xx,ff,weight)
         use mod_types
         integer    :: integrand
         real(dp15) :: xx(30),ff(1),weight
       end function integrand
    end interface
    integer    :: var(4),state(4),nvar,i,j,nlim
    character(len=3) :: name(4)

    if (nrec < npts) then
       write(6,'(a,i6,a)') '# integrability scan: only',nrec,' points with f /= 0 in the warm-up, increase -vegasNc0'
       return
    endif
    nohistos = .true.

    write(6,'(a)') '# integrability scan: median of |f| d^k, k = number of variables at the edge'
    write(6,'(a,i4,a)') '# depths d: 1e-2,1e-4,1e-6 (k=1), 1e-3..1e-5 (k=2), 1e-2..1e-3 (k=3), 3e-2..3e-3 (k=4);',npts,' points; n = points with f /= 0'
    write(6,'(a)') '# limit                     n    d1        d2        d3      verdict'

    !-- every variable alone, both edges
    do i = 1, ndim
       do j = 1, 2
          var(1) = i; state(1) = j
          call scan(1)
       enddo
    enddo

    !-- combined limits of the unresolved-parton variables
    if (ndim == NNLOdim_vegas) then
       nvar = 4; var = [xE5,xE6,xRHO5,xRHO6]; name = ['E5 ','E6 ','R5 ','R6 ']
    elseif (ndim == NLOdim_vegas .or. ndim == NLOdim_vegas + 1) then !-- one resolved parton (S sectors: + z)
       nvar = 2; var(1:2) = [xE,xRHO]; name(1:2) = ['E  ','R  ']
    else
       return
    endif
    !-- state: 0 = generic, 1 = x -> 0, 2 = x -> 1 (angles only)
    nlim = 3**nvar
    do i = 1, nlim - 1
       do j = 1, nvar
          state(j) = mod(i/3**(j-1),3)
       enddo
       if (count(state(1:nvar) /= 0) < 2) cycle
       if (any(state(1:nvar) == 2 .and. name(1:nvar)(1:1) == 'E')) cycle
       call scan(nvar,name)
    enddo

  contains

    subroutine scan(n,nm)
      integer, intent(in) :: n
      character(len=3), intent(in), optional :: nm(:)
      real(dp15) :: x(30),ff(1),wgt,val(npts),med(ndep)
      integer    :: id,ip,k,l,nz(ndep),ir
      character(len=24) :: label
      character(len=8)  :: verdict

      label = ''
      do l = 1, n
         if (state(l) == 0) cycle
         if (present(nm)) then
            label = trim(label)//' '//trim(nm(l))//merge('>0','>1',state(l) == 1)
         else
            write(label,'(a,i2,a)') 'x',var(l),merge('>0','>1',state(l) == 1)
         endif
      enddo
      k = count(state(1:n) /= 0)

      do id = 1, ndep
         nz(id) = 0
         do ip = 1, npts
            x = pts(:,ip)
            do l = 1, n
               if (state(l) == 1) x(var(l)) = depth(id,k)
               if (state(l) == 2) x(var(l)) = 1 - depth(id,k)
            enddo
            ff = 0; wgt = 1
            ir = integrand(x,ff,wgt)
            if (ff(1) /= 0 .and. ff(1) == ff(1)) then
               nz(id) = nz(id) + 1
               val(nz(id)) = abs(ff(1))*depth(id,k)**k
            endif
         enddo
         med(id) = median(val(1:nz(id)))
      enddo

      !-- compare the first depth with enough points to the last one
      id = 1
      do while (id < ndep .and. nz(id) < npts/10)
         id = id + 1
      enddo
      if (nz(id) < npts/10 .or. id == ndep .or. nz(ndep) < npts/10) then
         verdict = 'few'
      elseif (med(ndep) < merge(0.3_dp15,0.6_dp15,k <= 2)*med(id)) then
         verdict = 'ok'
      else
         verdict = 'FAIL'
      endif
      write(6,'(a2,a24,i4,3es10.2,2x,a)') '  ',adjustl(label),maxval(nz),med,verdict
    end subroutine scan

  end subroutine check_integrability

  !-- -checklim for the RR sectors. For random points that pass the cuts, the variables of each limit are set to
  !-- 1e-1..1e-10 and H and the counterterms are printed side by side from limval_nnlo (slots as in Checks/mod_check_lim:
  !-- 1 H, 2 S5, 3 S6, 4 X, 5 Y, 6 S5+S6, ...), with S5 = xE5, S6 = xE6, X = xRHO5 (TC or C5), Y = xRHO6 (C, C5, C6, C56).
  !-- Sign (-1)^(number of limits), so the terms that must cancel agree. Limits and columns: the slots the sector fills.
  !-- The sectors draw their last two variables themselves: all variables are passed with override.
  subroutine check_limits_rr(integrand)
    procedure(integrand_if) :: integrand
    integer, parameter :: var(4) = [xE5,xE6,xRHO5,xRHO6]
    integer, parameter :: slot(0:15) = [1,2,3,6,4,7,9,12,5,8,10,13,11,14,15,16]
    character(2), parameter :: nm(4) = ['S5','S6','X ','Y ']
    real(dp15) :: y0(30),y(30),ff(1),w
    real(dp)   :: val(16,10,15),sgn(16)
    logical    :: used(16)
    character(11) :: lab(16)
    integer    :: b,k,il,ir,nempty

    do b = 0,15
       sgn(slot(b)) = (-1)**popcnt(b)
       lab(slot(b)) = 'H'
       if (b > 0) lab(slot(b)) = ''
       do k = 1,4
          if (.not. btest(b,k-1)) cycle
          if (len_trim(lab(slot(b))) > 0) lab(slot(b)) = trim(lab(slot(b)))//'+'
          lab(slot(b)) = trim(lab(slot(b)))//trim(nm(k))
       enddo
    enddo
    nohistos = .true.
    buff = 0; buff_rr = 0; call set_buff()
    override = .true.
    nempty = 0
    write(6,'(a)') '# -checklim RR: S5 = x(E5), S6 = x(E6), X = x(rho5) (TC or C5), Y = x(rho6) (C, C5, C6 or C56)'
    do
       call random_number(y0)
       if (eval(y0) == 0) cycle
       used = .false.
       do b = 1,15
          do il = 1,10
             y = y0
             do k = 1,4
                if (btest(b,k-1)) y(var(k)) = 10._dp15**(-il)
             enddo
             w = eval(y)
             val(:,il,b) = sgn*limval_nnlo
             used = used .or. limval_nnlo /= 0
          enddo
       enddo
       if (all(val(1,2:8,:) == 0)) then            !-- H cut at every depth used for the verdicts
          nempty = nempty + 1
          if (nempty == 1000) stop 'check_limits_rr: H vanishes in every limit at 1000 points that pass the cuts'
          cycle
       endif
       write(6,'(/,a,4es12.4)') ' point: x(E5), x(E6), x(rho5), x(rho6) =', y0(var)
       do b = 1,15
          if (.not. used(slot(b))) cycle
          write(6,'(/,a,/,a,*(a16))') ' limit '//trim(lab(slot(b))), '       x(E5)       x(E6)     x(rho5)     x(rho6)', &
               pack(adjustr(lab),used)
          do il = 1,10
             y = y0
             do k = 1,4
                if (btest(b,k-1)) y(var(k)) = 10._dp15**(-il)
             enddo
             write(6,'(4es12.4,*(es16.8))') y(var), pack(val(:,il,b),used)
          enddo
       enddo
       call uncovered_head()
       do b = 1,15
          if (.not. used(slot(b))) call uncovered(lab(slot(b)),popcnt(b),[(sum(sgn*val(:,il,b)), il = 1,10)])
       enddo
       call mypause()
    enddo

  contains

    real(dp) function eval(x)                  !-- runs the sector at x, returns H
      real(dp15), intent(in) :: x(30)
      real(dp15) :: xx(30)
      xx = x
      limval_nnlo = 0
      ff = 0; w = 1
      ir = integrand(xx,ff,w)
      eval = limval_nnlo(1)
    end function eval

  end subroutine check_limits_rr

  !-- -checklim for the sectors with one resolved emission: NLO real, RV, and S with a resolved parton (s_o*, at fixed z).
  !-- Same idea as check_limits_rr, reading
  !-- limval_onlo: kinematic group g (as stored by the sector: boosted leg 1, boosted leg 2, z = 1), slot 6(g-1)+k,
  !-- k = H, S, C1 (C for a final-state emitter), C2, SC1 (SC), SC2. Limits: S x(E) -> 0, C1 x(rho) -> 0, C2 x(rho) -> 1.
  subroutine check_limits_onlo(integrand)
    procedure(integrand_if) :: integrand
    character(3), parameter :: nm(6) = ['H  ','S  ','C1 ','C2 ','SC1','SC2']
    integer, parameter :: sgn(6) = [1,-1,-1,-1,1,1]
    integer, parameter :: soft(6) = [0,1,0,0,1,1], coll(6) = [0,0,1,2,1,2]   !-- limit of each slot
    real(dp15) :: y0(30),y(30),ff(1),w
    real(dp)   :: val(18,10,6)
    logical    :: used(18)
    character(6) :: lab(18)
    integer    :: g,k,il,ir,l,nempty

    do g = 1,3
       do k = 1,6
          write(lab(6*(g-1)+k),'(a,"/",i1)') trim(nm(k)),g
       enddo
    enddo
    nohistos = .true.
    buff = 0; buff_r = 0; call set_buff()
    override = .true.
    nempty = 0
    write(6,'(a)') '# -checklim NLO/RV/ONLO: S = x(E) -> 0, C1 = x(rho) -> 0, C2 = x(rho) -> 1; X/g = term X of kinematic group g'
    do
       call random_number(y0)
       if (eval(y0) == 0) cycle
       used = .false.
       do l = 2,6
          do il = 1,10
             y = point(l,il)
             w = eval(y)
             val(:,il,l) = [(sgn(mod(k-1,6)+1)*limval_onlo(k), k = 1,18)]
             used = used .or. limval_onlo /= 0
          enddo
       enddo
       if (all(val(1::6,2:8,2:6) == 0)) then       !-- H cut at every depth used for the verdicts
          nempty = nempty + 1
          if (nempty == 1000) stop 'check_limits_onlo: H vanishes in every limit at 1000 points that pass the cuts'
          cycle
       endif
       write(6,'(/,a,2es12.4)') ' point: x(E), x(rho) =', y0(xE), y0(xRHO)
       do l = 2,6
          if (.not. any(used(l::6))) cycle
          write(6,'(/,a,/,a,*(a16))') ' limit '//trim(nm(l)), '        x(E)      x(rho)', pack(adjustr(lab),used)
          do il = 1,10
             y = point(l,il)
             write(6,'(2es12.4,*(es16.8))') y(xE), y(xRHO), pack(val(:,il,l),used)
          enddo
       enddo
       call uncovered_head()
       do l = 2,6
          do g = 1,3
             if (.not. used(6*(g-1)+l) .and. any(val(6*(g-1)+1,:,l) /= 0)) &
                  call uncovered(lab(6*(g-1)+l),soft(l)+min(coll(l),1),[(sum(sgn*val(6*g-5:6*g,il,l)), il = 1,10)])
          enddo
       enddo
       call mypause()
    enddo

  contains

    function point(l,il) result(x)             !-- y0 moved into limit l at depth 1e-il
      integer, intent(in) :: l,il
      real(dp15) :: x(30),d
      d = 10._dp15**(-il)
      x = y0
      if (soft(l) == 1) x(xE) = d
      if (coll(l) == 1) x(xRHO) = d
      if (coll(l) == 2) x(xRHO) = 1 - d
    end function point

    real(dp) function eval(x)                  !-- runs the sector at x, returns the largest |H| of the groups
      real(dp15), intent(in) :: x(30)
      real(dp15) :: xx(30)
      xx = x
      limval_onlo = 0
      ff = 0; w = 1
      ir = integrand(xx,ff,w)
      eval = maxval(abs(limval_onlo(1::6)))
    end function eval

  end subroutine check_limits_onlo

  !-- W-awareness (PLAN C1.4): in a limit without its own counterterm (neutrino, gamma gamma TC for W; soft gluon covered
  !-- by C1 + C2) the integrand H + all counterterms must be integrable: x^k |sum| (k = variables at the edge) must fall
  !-- between depths 1e-2 and 1e-8, else a counterterm is missing. Points where the sum is cut or not finite are skipped.
  subroutine uncovered_head()
    write(6,'(/,a,/,a)') ' limits without own counterterm: x^k |H + counterterms| must fall (k = variables at the edge)', &
         ' limit             1e-02       1e-04       1e-06       1e-08  verdict'
  end subroutine uncovered_head

  subroutine uncovered(name,k,h)
    character(*), intent(in) :: name
    integer, intent(in)  :: k
    real(dp), intent(in) :: h(:)
    real(dp) :: f(10)
    integer  :: il,i1,i2
    character(10) :: verdict
    f = [(abs(h(il))*10._dp**(-k*il), il = 1,10)]
    i1 = findloc(f(2:8) /= 0 .and. f(2:8) == f(2:8),.true.,dim=1) + 1
    i2 = findloc(f(2:8) /= 0 .and. f(2:8) == f(2:8),.true.,dim=1,back=.true.) + 1
    if (i1 == 1 .or. i2 - i1 < 2) then
       verdict = 'cut'
    elseif (f(i2) < 0.3_dp*f(i1)) then
       verdict = 'integrable'
    else
       verdict = 'SINGULAR'
    endif
    write(6,'(1x,a14,4es12.4,2x,a)') name, f(2:8:2), trim(verdict)
  end subroutine uncovered

  function median(v) result(m)
    real(dp15), intent(in) :: v(:)
    real(dp15) :: m,w(size(v)),t
    integer :: i,j
    m = 0
    if (size(v) == 0) return
    w = v
    do i = 2, size(w)
       t = w(i); j = i - 1
       do while (j >= 1)
          if (w(j) <= t) exit
          w(j+1) = w(j); j = j - 1
       enddo
       w(j+1) = t
    enddo
    m = w((size(w)+1)/2)
  end function median

end module mod_check_int
