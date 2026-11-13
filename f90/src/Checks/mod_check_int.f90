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
  implicit none
  private

  public :: check_integrability, start_recording, rec_integrand

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
