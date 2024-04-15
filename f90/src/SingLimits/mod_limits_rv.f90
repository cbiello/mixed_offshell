module mod_limits_rv
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_process
  use mod_auxfunctions
  implicit none
  private

  public :: col_rvewk_is_ns, col_rvqcd_is_ns
  public :: col_rvewk_is_qg
  public :: col_rvqcd_is_qa

  !-- be careful to the z definition
  !-- these limits are supposed to be called with 1-z, i.e. with x2/(x2-1) instead of 1/(1-x2)
  !-- this is not to loose precision
  
contains

  !-- 1 -> PLOOP, 2 -> PTREE
  !-- z = eg/ei so that the soft singularity is for z->0 [precision]
  !-- non-abelian part removed -> no mu dependence
  function col_rvqcd_is_ns(proc) result(res)
    type(KinConfig), intent(in)  :: proc
    real(dp) :: res(2)
    real(dp) :: z,s5i,etaij
    real(dp) :: omz,lomz,lz,leta,Li2z
    real(dp) :: pnew,ptree

    z     = proc%Lim_KinInv(1)
    s5i   = proc%Lim_KinInv(2)
    etaij = proc%Lim_KinInv(3)

    omz = one-z
    lz = log(z)
    lomz = log(omz)
    leta = log(etaij)
    Li2z = real(dilog2(z),kind=dp)

    ptree = (one+(one-z)**2)/z
    ptree = ptree/(one-z) !-- ``flux'' 1/z factor, with z <-> 1-z

    pnew  = -one

    res(1) = -pnew + &
         (2*Li2z + 3*lomz + 2*leta*lomz + 2*lomz*lz - lomz**2)*ptree

    res(2) = ptree

    res = res * two/s5i ! * Qq**2
    
  end function col_rvqcd_is_ns

  
  !-- 1 -> PLOOP prop to Qq^2, 2 -> PLOOP prop to Qq Ql, 3 --> PTREE
  !-- z = eg/ei so that the soft singularity is for z->0 [precision]
  !-- non-abelian part removed -> no mu dependence
  !-- part 2 is only there because of a bad definition of I1EW
  function col_rvewk_is_ns(i,proc) result(res)
    integer, intent(in) :: i
    type(KinConfig), intent(in)  :: proc
    real(dp) :: res(3)
    real(dp) :: z,s5i,etaij,si3,si4
    real(dp) :: omz,lomz,lz,leta,Li2z
    real(dp) :: pnew,ptree

    z     = proc%Lim_KinInv(1)
    s5i   = proc%Lim_KinInv(2)
    etaij = proc%Lim_KinInv(3)

    si3 = two*scr(proc%AmpMom(:,i),proc%AmpMom(:,3))
    si4 = two*scr(proc%AmpMom(:,i),proc%AmpMom(:,4))
    
    omz = one-z
    lz = log(z)
    lomz = log(omz)
    leta = log(etaij)
    Li2z = real(dilog2(z),kind=dp)

    ptree = (one+(one-z)**2)/z
    ptree = ptree/(one-z) !-- ``flux'' 1/z factor, with z <-> 1-z

    pnew  = -one

    res(1) = -pnew + &
         (2*Li2z + 3*lomz + 2*leta*lomz + 2*lomz*lz - lomz**2)*ptree

    res(2) = -two*log(si3/si4)*lomz * ptree
    res(3) = ptree

    res = res * two/s5i ! * Cf
    
  end function col_rvewk_is_ns
  
  !-- 1 -> PLOOP prop to Qq^2, 2 -> PLOOP prop to Qq Ql, 3 --> PTREE
  !-- z = eg/ei so that the soft singularity is for z->0 [precision]
  !-- non-abelian part removed -> no mu dependence
  !-- part 2 is only there because of a bad definition of I1EW
  function col_rvewk_is_qg(i,proc) result(res)
    integer, intent(in) :: i
    type(KinConfig), intent(in)  :: proc
    real(dp) :: res(3)
    real(dp) :: z,s5i,etaij,si3,si4
    real(dp) :: omz,lomzOz,lz,leta,Li2z
    real(dp) :: pnew,ptree

    z     = proc%Lim_KinInv(1)
    s5i   = proc%Lim_KinInv(2)
    etaij = proc%Lim_KinInv(3)

    si3 = two*scr(proc%AmpMom(:,i),proc%AmpMom(:,3))
    si4 = two*scr(proc%AmpMom(:,i),proc%AmpMom(:,4))
    
    omz = one-z
    lz = log(z)
    lomzOz = log(omz/z)
    leta = log(etaij)
    Li2z = real(dilog2(z),kind=dp)

    ptree = z**2+omz**2
    ptree = ptree/omz

    pnew  = one

    res(1) = -pnew + &
         (5._dp/three*pisq + two*lomzOz*leta - lz**2 - lomzOz**2 + three*lomzOz-two*Li2z)*ptree

    res(2) = -two*log(si3/si4)*lomzOz * ptree
    res(3) = ptree

    res = res * two/s5i ! * Tr
    
  end function col_rvewk_is_qg

  !-- 1 -> PLOOP, 2 --> PTREE
  !-- z = eg/ei so that the soft singularity is for z->0 [precision]
  !-- non-abelian part removed -> no mu dependence
  function col_rvqcd_is_qa(proc) result(res)
    type(KinConfig), intent(in)  :: proc
    real(dp) :: res(2)
    real(dp) :: z,s5i,etaij
    real(dp) :: omz,lomzOz,lz,leta,Li2z
    real(dp) :: pnew,ptree

    z     = proc%Lim_KinInv(1)
    s5i   = proc%Lim_KinInv(2)
    etaij = proc%Lim_KinInv(3)

    omz = one-z
    lz = log(z)
    lomzOz = log(omz/z)
    leta = log(etaij)
    Li2z = real(dilog2(z),kind=dp)

    ptree = z**2+omz**2
    ptree = ptree/omz

    pnew  = one

    res(1) = -pnew + &
         (5._dp/three*pisq + two*lomzOz*leta - lz**2 - lomzOz**2 + three*lomzOz-two*Li2z)*ptree

    res(2) = ptree

    res = res * two/s5i ! * xn * Qq^2, xn = aveqa/aveqq
    
  end function col_rvqcd_is_qa
  
end module mod_limits_rv

