module mod_auxlumi_wp
  use mod_types; use common_def; use mod_consts_dp
  implicit none
  private

  public :: ns_lumi_lo, ns_lumi_nlo, ns_lumi_nnlo
  public :: qg_lumi_nlo, qg_lumi_lo
  public :: gq_lumi_nlo, gq_lumi_lo
  !
  public :: ns_lumi_nlo_sping, ns_lumi_gq_sping, ns_lumi_qg_sping
  public :: qg_lumi_nnlo, qg_lumi_qg_sping, qg_lumi_ns, gq_lumi_ns, &
            qg_lumi_rv, qg_lumi_rv_lo, gq_lumi_rv, gq_lumi_rv_lo
  public :: gq_lumi_nnlo,gq_lumi_qg_sping
  public :: s_lumi_nnlo, s_lumi_gq_sping, s_lumi_qg_sping
  public :: s_lumi_xj, s_lumi_xj_id
  public :: s_lumi_ix, s_lumi_ix_id
  public :: gg_lumi_nnlo, gg_lumi_qxg, gg_lumi_qg, gg_lumi_gq, gg_lumi_gqx, &
            gg_lumi_qqx, gg_lumi_qxq
  public :: s_lumi_gq,s_lumi_qg,ns_lumi_gq,ns_lumi_qg
! 
  public :: ns_lumi_nnlo_bbbar,s_lumi_nnlo_lbbar,s_lumi_nnlo_bl,s_lumi_nnlo_bbarl

contains

  function ns_lumi_lo(res,f1,f2) result(respdf)
    real(dp), intent(in) :: res(2),f1(-6:6),f2(-6:6)
    real(dp) :: respdf

    !-- u dx > nu ep, c sx > nu ep
    respdf = res(1)*sum(f1(2:nf_pdf:2)*f2(-1:-3:-2))

    !-- dx u > nu ep, sx c > nu ep
    respdf = respdf + res(2)*sum(f1(-1:-3:-2)*f2(2:nf_pdf:2))

!    respdf = res(1)*(f1(2)*f2(-1)+f1(4)*f2(-3))
!    respdf = respdf + res(2)*(f1(-1)*f2(+2)+f1(-3)*f2(+4))

  end function ns_lumi_lo

  function ns_lumi_nlo(res_j,f1,f2) result(respdf)
    real(dp), intent(in) :: res_j(6),f1(-6:6),f2(-6:6)
    real(dp) :: respdf


    !-- u dx > g ve ep, c sx > g ve ep
    respdf = res_j(1)*sum(f1(2:nf_pdf:2)*f2(-1:-3:-2))
    
    !-- dx u > g ve ep, sx c > g =ve ep
    respdf = respdf + res_j(4)*sum(f1(-1:-3:-2)*f2(2:nf_pdf:2))
    
    ! print *, 'respdf', respdf
    ! pause

  end function ns_lumi_nlo

  function ns_lumi_nnlo(res,f1,f2) result(respdf)
    ! including b final states
    real(dp), intent(in) :: res(12),f1(-6:6),f2(-6:6)
    real(dp) :: respdf

    !-- u dx > nu ep, c sx > nu ep
    respdf = (res(1)+res(5))*sum(f1(2:nf_pdf:2)*f2(-1:-3:-2))

    !-- dx u > nu ep, sx c > nu ep
    respdf = respdf + (res(2)+res(8))*sum(f1(-1:-3:-2)*f2(2:nf_pdf:2))

    ! u ux 
    respdf = respdf + res(3)*sum(f1(2:nf_pdf:2)*f2(-2:-nf_pdf:-2))

    ! d dx 
    respdf = respdf + res(4)*sum(f1(1:3:2)*f2(-1:-3:-2))

    ! ux u 
    respdf = respdf + res(6)*sum(f1(-2:-nf_pdf:-2)*f2(2:nf_pdf:2))

    ! dx d 
    respdf = respdf + res(7)*sum(f1(-1:-3:-2)*f2(1:3:2))

    ! b bx
    respdf = respdf + res(9)*f1(5)*f2(-5)

    ! bx b
    respdf = respdf + res(10)*f1(-5)*f2(5)

    !-- u dx > b bx nu ep, c sx > b bx nu ep
    respdf = respdf + res(11)*sum(f1(2:nf_pdf:2)*f2(-1:-3:-2))

    !-- dx u > b bx nu ep, sx c > b bx nu ep
    respdf = respdf + res(12)*sum(f1(-1:-3:-2)*f2(2:nf_pdf:2))



!    respdf = res(1)*(f1(2)*f2(-1)+f1(4)*f2(-3))
!    respdf = respdf + res(2)*(f1(-1)*f2(+2)+f1(-3)*f2(+4))

  end function ns_lumi_nnlo

  function ns_lumi_nnlo_bbbar(res,f1,f2) result(respdf)
    ! separate contributions for b bbar in final state
    real(dp), intent(in) :: res(12),f1(-6:6),f2(-6:6)
    real(dp) :: respdf

    !-- u dx > b bx nu ep, c sx > b bx nu ep
    respdf = res(11)*sum(f1(2:nf_pdf:2)*f2(-1:-3:-2))

    !-- dx u > b bx nu ep, sx c > b bx nu ep
    respdf = respdf + res(12)*sum(f1(-1:-3:-2)*f2(2:nf_pdf:2))

  end function ns_lumi_nnlo_bbbar


  function ns_lumi_nlo_sping(res_j_sping,f1,f2) result(respdf_sping)
    complex(dp), intent(in) :: res_j_sping(6,-1:1,-1:1)
    real(dp), intent(in) :: f1(-6:6),f2(-6:6)
    complex(dp) :: respdf_sping(-1:1,-1:1)


    !-- u dx > g ve ep, c sx > g ve ep
    respdf_sping(:,:) = res_j_sping(1,:,:)*sum(f1(2:nf_pdf:2)*f2(-1:-3:-2))
    
    !-- dx u > g ve ep, cs c > g ve ep
    respdf_sping(:,:) = respdf_sping(:,:) + res_j_sping(4,:,:)*sum(f1(-1:-3:-2)*f2(2:nf_pdf:2))
    
    ! print *, 'respdf_sping(1,1)', respdf_sping(1,1)
    ! pause

  end function ns_lumi_nlo_sping

  !-- ns luminosity multiplying gq reduced matrix element, for collinear limits
  function ns_lumi_gq_sping(res_j_sping,f1,f2) result(respdf_sping)
    complex(dp), intent(in) :: res_j_sping(6,-1:1,-1:1)
    real(dp), intent(in) :: f1(-6:6),f2(-6:6)
    complex(dp) :: respdf_sping(-1:1,-1:1)

    !--

    !-- u dx > jj em ep --> g dx
    respdf_sping(:,:) = res_j_sping(3,:,:)*sum(f1(2:nf_pdf:2)*f2(-1:-3:-2))
    
    !-- dx u > jj em ep --> g u
    respdf_sping(:,:) = respdf_sping(:,:) + res_j_sping(2,:,:)*sum(f1(-1:-3:-2)*f2(2:nf_pdf:2))

    ! print *, 'respdf_sping(1,1)', respdf_sping(1,1)
    ! pause

  end function ns_lumi_gq_sping

  function ns_lumi_qg_sping(res_j_sping,f1,f2) result(respdf_sping)
    complex(dp), intent(in) :: res_j_sping(6,-1:1,-1:1)
    real(dp), intent(in) :: f1(-6:6),f2(-6:6)
    complex(dp) :: respdf_sping(-1:1,-1:1)

    !-- u dx > jj em ep --> u g
    respdf_sping(:,:) = res_j_sping(5,:,:)*sum(f1(2:nf_pdf:2)*f2(-1:-3:-2))
    
    !-- dx u > jj em ep --> dx g
    respdf_sping(:,:) = respdf_sping(:,:) + res_j_sping(6,:,:)*sum(f1(-1:-3:-2)*f2(2:nf_pdf:2))

    ! print *, 'respdf_sping(1,1)', respdf_sping(1,1)
    ! pause
     
  end function ns_lumi_qg_sping

  !-- ns luminosity multiplying gq reduced matrix element, for subtraction terms
  function ns_lumi_gq(res_j,f1,f2) result(respdf)
    real(dp), intent(in) :: res_j(6)
    real(dp), intent(in) :: f1(-6:6),f2(-6:6)
    complex(dp) :: respdf

    !-- u dx > jj em ep --> g dx
    respdf = res_j(3)*sum(f1(2:nf_pdf:2)*f2(-1:-3:-2))
    
    !-- dx u > jj em ep --> g u
    respdf = respdf + res_j(2)*sum(f1(-1:-3:-2)*f2(2:nf_pdf:2))

    ! print *, 'respdf', respdf
    ! pause

  end function ns_lumi_gq

  function ns_lumi_qg(res_j,f1,f2) result(respdf)
    real(dp), intent(in) :: res_j(6)
    real(dp), intent(in) :: f1(-6:6),f2(-6:6)
    complex(dp) :: respdf


    !-- u dx > jj em ep --> u g
    respdf = res_j(5)*sum(f1(2:nf_pdf:2)*f2(-1:-3:-2))
    
    !-- dx u > jj em ep --> dx g
    respdf = respdf + res_j(6)*sum(f1(-1:-3:-2)*f2(2:nf_pdf:2))

    ! print *, 'respdf', respdf
    ! pause
     
  end function ns_lumi_qg

  !-- 

  function qg_lumi_nlo(res_j,f1,f2) result(respdf)
    real(dp), intent(in) :: res_j(6),f1(-6:6),f2(-6:6)
    real(dp) :: respdf


    !-- u g > d ve ep, c g > s ve ep
    respdf = res_j(5)*sum(f1(2:nf_pdf:2))*f2(0)
    
    !-- dx g > ux ve ep, sx g > cx ve ep
    respdf = respdf + res_j(6)*sum(f1(-1:-3:-2))*f2(0)
    

    ! print *, 'respdf', respdf
    ! pause

  end function qg_lumi_nlo

  function qg_lumi_qg_sping(res_j_sping,f1,f2) result(respdf_sping)
    complex(dp), intent(in) :: res_j_sping(6,-1:1,-1:1)
    real(dp), intent(in) :: f1(-6:6),f2(-6:6)
    complex(dp) :: respdf_sping(-1:1,-1:1)


    !-- u g > d ve ep, c g > s ve ep
    respdf_sping(:,:) = res_j_sping(5,:,:)*sum(f1(2:nf_pdf:2))*f2(0)
    
    !-- dx g > ux ve ep, sx g > cx ve ep
    respdf_sping(:,:) = respdf_sping(:,:) + res_j_sping(6,:,:)*sum(f1(-1:-3:-2))*f2(0)

    ! print *, 'respdf_sping(1,1)', respdf_sping(1,1)
    ! pause
    
  end function qg_lumi_qg_sping


  function gq_lumi_qg_sping(res_j_sping,f1,f2) result(respdf_sping)
    complex(dp), intent(in) :: res_j_sping(6,-1:1,-1:1)
    real(dp), intent(in) :: f1(-6:6),f2(-6:6)
    complex(dp) :: respdf_sping(-1:1,-1:1)


    !-- g u > d ve ep, g c > s ve ep
    respdf_sping(:,:) = res_j_sping(2,:,:)*sum(f2(2:nf_pdf:2))*f1(0)
    
    !-- g dx > ux ve ep, g sx > cx ve ep
    respdf_sping(:,:) = respdf_sping(:,:) + res_j_sping(3,:,:)*sum(f2(-1:-3:-2))*f1(0)

    ! print *, 'respdf_sping(1,1)', respdf_sping(1,1)
    ! pause
    
  end function gq_lumi_qg_sping

  function qg_lumi_lo(res,f1,f2) result(respdf)
    real(dp), intent(in) :: res(2),f1(-6:6),f2(-6:6)
    real(dp) :: respdf

    !-- u g > d ve ep, c g > s ve ep
    respdf = res(1)*sum(f1(2:nf_pdf:2))*f2(0)
    
    !-- dx g > ux ve ep, sx g > cx ve ep
    respdf = respdf + res(2)*sum(f1(-1:-3:-2))*f2(0)
    
    ! print *, 'respdf', respdf
    ! pause

  end function qg_lumi_lo

  !--

  function gq_lumi_nlo(res_j,f1,f2) result(respdf)
    real(dp), intent(in) :: res_j(6),f1(-6:6),f2(-6:6)
    real(dp) :: respdf


    !-- g u > d ve ep, g c > s ve ep
    respdf = res_j(2)*f1(0)*sum(f2(2:nf_pdf:2))

    !-- g dx > ux ve ep, g sx > cx ve ep
    respdf = respdf + res_j(3)*f1(0)*sum(f2(-1:-3:-2))
    
  end function gq_lumi_nlo
    
  function gq_lumi_lo(res,f1,f2) result(respdf)
    real(dp), intent(in) :: res(2),f1(-6:6),f2(-6:6)
    real(dp) :: respdf

    !-- g u > d ve ep, g c > s ve ep
    respdf = res(2)*f1(0)*sum(f2(2:nf_pdf:2))

    !-- g dx > ux ve ep, g sx > cx ve ep
    respdf = respdf + res(1)*f1(0)*sum(f2(-1:-3:-2))
    
  end function gq_lumi_lo
  
  !--

  function qg_lumi_nnlo(res_jj,f1,f2) result(respdf)
    real(dp), intent(in) :: res_jj(2),f1(-6:6),f2(-6:6)
    real(dp) :: respdf

    !-- u g > d g em ep, c g > s g em ep
    respdf = res_jj(1)*sum(f1(2:nf_pdf:2))*f2(0)
    
    !-- dx g > ux g em ep, sx g > cx g em ep
    respdf = respdf + res_jj(2)*sum(f1(-1:-3:-2))*f2(0)

    ! print *, 'respdf', respdf
    ! pause

  end function qg_lumi_nnlo


 function gq_lumi_nnlo(res_jj,f1,f2) result(respdf)
    real(dp), intent(in) :: res_jj(2),f1(-6:6),f2(-6:6)
    real(dp) :: respdf

    !-- g u  > d g ve ep, c g > s g ve ep
    respdf = res_jj(1)*sum(f2(2:nf_pdf:2))*f1(0)
    
    !-- g dx > ux g ve ep, sx g > cx g ve ep
    respdf = respdf + res_jj(2)*sum(f2(-1:-3:-2))*f1(0)

  end function gq_lumi_nnlo

  function qg_lumi_ns(res_j,f1,f2) result(respdf)
    real(dp), intent(in) :: res_j(6),f1(-6:6),f2(-6:6)
    real(dp) :: respdf
    !--

    !-- u g > d g --> udx
    respdf = res_j(1)*sum(f1(2:nf_pdf:2))*f2(0)
    
    !-- dx g > ux g --> dx u
    respdf = respdf + res_j(4)*sum(f1(-1:-3:-2))*f2(0)

    ! print *, 'respdf', respdf
    ! pause

  end function qg_lumi_ns

 function gq_lumi_ns(res_j,f1,f2) result(respdf)
    real(dp), intent(in) :: res_j(6),f1(-6:6),f2(-6:6)
    real(dp) :: respdf

    !-- g u  > d g --> dx u
    respdf = res_j(4)*f1(0)*sum(f2(2:nf_pdf:2))
    
    !-- g dx  > ux g --> u dx 
    respdf = respdf + res_j(1)*f1(0)*sum(f2(-1:-3:-2))


  end function gq_lumi_ns

  !--

  function s_lumi_nnlo(res_jj,f1,f2) result(respdf)
    real(dp), intent(in) :: res_jj(26),f1(-6:6),f2(-6:6)
    real(dp) :: respdf
    real(dp) :: restmp
    integer :: i,j
    
    respdf = zero

    !-- d u > d d em vex b bx, s c > s s  em vex b bx   
    respdf = respdf + res_jj(1)*(f1(1)*f2(2) + f1(3)*f2(4))

    !-- d c > d s em vex b bx, s u > s d em vex b bx, u c > u s em vex b bx
    respdf = respdf + res_jj(2)*(f1(1)*f2(4) + f1(3)*f2(2) + f1(2)*f2(4))

    !-- d sx > d cx em vex b bx, s dx > s ux em vex b bx, u sx > u cx em vex b bx, c dx > c ux em vex b bx                                                                                               
    respdf = respdf + res_jj(3)*(f1(1)*f2(-3) + f1(3)*f2(-1) + f1(2)*f2(-3) + f1(4)*f2(-1) + f1(1)*f2(-1) + f1(3)*f2(-3))

    !-- u u > d u em vex b bx, c c > s c em vex b bx                                                                                                                                                     
    respdf = respdf + res_jj(4)*(f1(2)*f2(2) + f1(4)*f2(4))

    !--u c > d c em vex b bx
    respdf = respdf + res_jj(5)*(f1(2)*f2(4))

    !-- u sx > d sx em vex b bx, u cx > d cx em vex b bx, c dx > s dx em vex b bx, c ux > s ux em vex b bx
    respdf = respdf + res_jj(6)*(f1(2)*f2(-3) + f1(2)*f2(-4) + f1(4)*f2(-1) + f1(4)*f2(-2) + f1(2)*f2(-2) + f1(4)*f2(-4))

    !-- dx dx > dx ux em vex b bx, sx sx > sx cx em vex b bx
    respdf = respdf + res_jj(7)*(f1(-1)*f2(-1) + f1(-3)*f2(-3))

    !-- dx sx > dx cx em vex b bx
    respdf = respdf + res_jj(8)*(f1(-1)*f2(-3))

    !-- dx sx > sx ux em vex b bx, dx cx > cx ux em vex b bx, sx ux > ux cx em vex b bx
    respdf = respdf + res_jj(9)*(f1(-1)*f2(-3) + f1(-1)*f2(-4) + f1(-3)*f2(-2))

    !-- dx ux > ux ux em vex b bx, sx cx > cx cx em vex b bx    
    respdf = respdf + res_jj(10)*(f1(-1)*f2(-2) + f1(-3)*f2(-4))

    !-- u d > d d em vex b bx, c s > s s em vex b bx        
    respdf = respdf + res_jj(11)*(f1(2)*f2(1) + f1(4)*f2(3))

    !-- c d > d s em vex b bx, u s > s d em vex b bx, c u > u s em vex b bx
    respdf = respdf + res_jj(12)*(f1(4)*f2(1) + f1(2)*f2(3) + f1(4)*f2(2))

    !-- sx d > d cx em vex b bx, dx s > s ux em vex b bx, sx u > u cx em vex b bx, dx c > c ux em vex b bx
    respdf = respdf + res_jj(13)*(f1(-3)*f2(1) + f1(-1)*f2(3) + f1(-3)*f2(2) + f1(-1)*f2(4) + f1(-1)*f2(1) + f1(-3)*f2(3))

    !-- c u > d c em vex b bx                 
    respdf = respdf + res_jj(14)*(f1(4)*f2(2))
    
    !-- sx u > d sx em vex b bx, cx u > d cx em vex b bx, dx c > s dx em vex b bx, ux c > s ux em vex b bx
    respdf = respdf + res_jj(15)*(f1(-3)*f2(2) + f1(-4)*f2(2) + f1(-1)*f2(4) + f1(-2)*f2(4)+ f1(-2)*f2(2) + f1(-4)*f2(4))

    !-- sx dx > dx cx em vex b bx
    respdf = respdf + res_jj(16)*(f1(-3)*f2(-1))

    !-- sx dx > sx ux em vex b bx, cx dx > cx ux em vex b bx, ux sx > ux cx em vex b bx
    respdf = respdf + res_jj(17)*(f1(-3)*f2(-1) + f1(-4)*f2(-1) + f1(-2)*f2(-3))

    !-- ux dx > ux ux em vex b bx, cx sx > cx cx em vex b bx
    respdf = respdf + res_jj(18)*(f1(-2)*f2(-1) + f1(-4)*f2(-3))

    !--b u > b d em vex b bx, b c > b s em vex b bx
    respdf = respdf + res_jj(19)*f1(5)*(f2(2) + f2(4))

    !--u bx > d bx em vex b bx, c bx > s bx em vex b bx
    respdf = respdf + res_jj(20)*(f1(2)+f1(4))*f2(-5)

    !--b dx > b ux em vex b bx, b sx > b cx em vex b bx
    respdf = respdf + res_jj(21)*f1(5)*(f2(-1)+f2(-3))

    !--bx dx > bx ux em vex b bx, bx sx > bx x em vex b bx
    respdf = respdf + res_jj(22)*f1(-5)*(f2(-1)+f2(-3))

    !- u b > b d em vex b bx, c b > b s em vex b bx
    respdf = respdf + res_jj(23)*(f1(2)+f1(4))*f2(5)

    !--bx u > d bx em vex b bx, bx c > s bx em vex b bx
    respdf = respdf + res_jj(24)*f1(-5)*(f2(2)+f2(4))

    !--dx b > b ux em vex b bx, sx b  > b cx em vex b bx
    respdf = respdf + res_jj(25)*(f1(-1)+f1(-3))*f2(5)

    !--dx bx > bx ux em vex b bx , sx bx > bx cx em vex b bx
    respdf = respdf + res_jj(26)*(f1(-1)+f1(-3))*f2(-5)


  end function s_lumi_nnlo


  function s_lumi_nnlo_bl(res_jj,f1,f2) result(respdf)
    real(dp), intent(in) :: res_jj(26),f1(-6:6),f2(-6:6)
    real(dp) :: respdf
    real(dp) :: restmp
    integer :: i,j
    
    respdf = zero

    !--b u > b d em vex b bx, b c > b s em vex b bx
    respdf = respdf + res_jj(19)*f1(5)*(f2(2) + f2(4))

    !--b dx > b ux em vex b bx, b sx > b cx em vex b bx
    respdf = respdf + res_jj(21)*f1(5)*(f2(-1)+f2(-3))

    !- u b > b d em vex b bx, c b > b s em vex b bx
    respdf = respdf + res_jj(23)*(f1(2)+f1(4))*f2(5)

    !--dx b > b ux em vex b bx, sx b  > b cx em vex b bx
    respdf = respdf + res_jj(25)*(f1(-1)+f1(-3))*f2(5)

  end function s_lumi_nnlo_bl

  function s_lumi_nnlo_lbbar(res_jj,f1,f2) result(respdf)
    real(dp), intent(in) :: res_jj(26),f1(-6:6),f2(-6:6)
    real(dp) :: respdf
    real(dp) :: restmp
    integer :: i,j

    respdf = zero

    !--u bx > d bx em vex b bx, c bx > s bx em vex b bx
    respdf = respdf + res_jj(20)*(f1(2)+f1(4))*f2(-5)

    !--bx u > d bx em vex b bx, bx c > s bx em vex b bx
    respdf = respdf + res_jj(24)*f1(-5)*(f2(2)+f2(4))

  end function s_lumi_nnlo_lbbar

  function s_lumi_nnlo_bbarl(res_jj,f1,f2) result(respdf)
    real(dp), intent(in) :: res_jj(26),f1(-6:6),f2(-6:6)
    real(dp) :: respdf
    real(dp) :: restmp
    integer :: i,j

    respdf = zero

    !--bx dx > bx ux em vex b bx, bx sx > bx x em vex b bx
    respdf = respdf + res_jj(22)*f1(-5)*(f2(-1)+f2(-3))

    !--dx bx > bx ux em vex b bx , sx bx > bx cx em vex b bx
    respdf = respdf + res_jj(26)*(f1(-1)+f1(-3))*f2(-5)


  end function s_lumi_nnlo_bbarl




    

  function s_lumi_gq_sping(res_j_sping,f1,f2) result(respdf_sping)
    complex(dp), intent(in) :: res_j_sping(6,-1:1,-1:1)
    real(dp), intent(in) :: f1(-6:6),f2(-6:6)
    complex(dp) :: respdf_sping(-1:1,-1:1)
    real(dp) :: restmp
    integer :: i,j


    respdf_sping=czero
    respdf_sping(:,:) = respdf_sping(:,:) + res_j_sping(2,:,:)*(f1(-5)+f1(-4)+f1(-3)+f1(-2)+f1(1)+f1(2)+f1(3)+f1(4)+f1(5))*f2(2)    ! g u
    respdf_sping(:,:) = respdf_sping(:,:) + res_j_sping(2,:,:)*(f1(-5)+f1(-4)+f1(-2)+f1(-1)+f1(1)+f1(2)+f1(3)+f1(4)+f1(5))*f2(4)    ! g c
    respdf_sping(:,:) = respdf_sping(:,:) + res_j_sping(3,:,:)*(f1(-5)+f1(-4)+f1(-3)+f1(-2)+f1(-1)+f1(1)+f1(3)+f1(4)+f1(5))*f2(-1)  ! g dx
    respdf_sping(:,:) = respdf_sping(:,:) + res_j_sping(3,:,:)*(f1(-5)+f1(-4)+f1(-3)+f1(-2)+f1(-1)+f1(1)+f1(2)+f1(3)+f1(5))*f2(-3)  ! g sx

  end function s_lumi_gq_sping

  function s_lumi_gq(res_j,f1,f2) result(respdf)
    real(dp), intent(in) :: res_j(6)
    real(dp), intent(in) :: f1(-6:6),f2(-6:6)
     real(dp) :: respdf
    real(dp) :: restmp
    integer :: i,j

!    restmp = 0
!    do i = -nf_pdf, nf_pdf
!       do j = 2, nf_pdf, 2
!          if (i .eq. -j .or. i == 0) cycle
!          restmp = restmp + f1(i)*f2(j)
!       enddo
!    enddo
!    respdf = restmp * res_j(2) !-- g u
!
!    restmp = 0
!    do i = -nf_pdf, nf_pdf
!       do j = 1, 3, 2    ! removed b quark here
!          if (i .eq. j .or. i==0) cycle
!          restmp = restmp + f1(i)*f2(-j)
!       enddo
!    enddo
!    respdf = respdf + restmp * res_j(3) !-- g dx

    respdf = zero
    respdf = respdf + res_j(2)*(f1(-5)+f1(-4)+f1(-3)+f1(-2)+f1(1)+f1(2)+f1(3)+f1(4)+f1(5))*f2(2)    ! g u
    respdf = respdf + res_j(2)*(f1(-5)+f1(-4)+f1(-2)+f1(-1)+f1(1)+f1(2)+f1(3)+f1(4)+f1(5))*f2(4)    ! g c
    respdf = respdf + res_j(3)*(f1(-5)+f1(-4)+f1(-3)+f1(-2)+f1(-1)+f1(1)+f1(3)+f1(4)+f1(5))*f2(-1)  ! g dx
    respdf = respdf + res_j(3)*(f1(-5)+f1(-4)+f1(-3)+f1(-2)+f1(-1)+f1(1)+f1(2)+f1(3)+f1(5))*f2(-3)  ! g sx



  end function s_lumi_gq

  function s_lumi_qg_sping(res_j_sping,f1,f2) result(respdf_sping)
    complex(dp), intent(in) :: res_j_sping(6,-1:1,-1:1)
    real(dp), intent(in) :: f1(-6:6),f2(-6:6)
    complex(dp) :: respdf_sping(-1:1,-1:1)
    real(dp) :: restmp
    integer :: i,j


    respdf_sping=czero
    respdf_sping(:,:) = respdf_sping(:,:) + res_j_sping(5,:,:)*(f2(-5)+f2(-4)+f2(-3)+f2(-2)+f2(1)+f2(2)+f2(3)+f2(4)+f2(5))*f1(2)    ! u g
    respdf_sping(:,:) = respdf_sping(:,:) + res_j_sping(5,:,:)*(f2(-5)+f2(-4)+f2(-2)+f2(-1)+f2(1)+f2(2)+f2(3)+f2(4)+f2(5))*f1(4)    ! c g
    respdf_sping(:,:) = respdf_sping(:,:) + res_j_sping(6,:,:)*(f2(-5)+f2(-4)+f2(-3)+f2(-2)+f2(-1)+f2(1)+f2(3)+f2(4)+f2(5))*f1(-1)  ! dx g
    respdf_sping(:,:) = respdf_sping(:,:) + res_j_sping(6,:,:)*(f2(-5)+f2(-4)+f2(-3)+f2(-2)+f2(-1)+f2(1)+f2(2)+f2(3)+f2(5))*f1(-3)  ! sx g


!    restmp = 0
!    do i = -nf_pdf, nf_pdf
!       do j = 2, nf_pdf, 2
!          if (i .eq. -j .or. i == 0) cycle
!          restmp = restmp + f1(j)*f2(i)
!       enddo
!    enddo
!    respdf_sping(:,:) = restmp * res_j_sping(5,:,:) !-- u g
!
!    restmp = 0
!    do i = -nf_pdf, nf_pdf
!       do j = 1, 3, 2
!          if (i .eq. j .or. i==0) cycle
!          restmp = restmp + f1(-j)*f2(i)
!       enddo
!    enddo
!    respdf_sping(:,:) = respdf_sping(:,:) + restmp * res_j_sping(6,:,:) !-- dx g

    
  end function s_lumi_qg_sping

  function s_lumi_qg(res_j,f1,f2) result(respdf)
    real(dp), intent(in) :: res_j(6)
    real(dp), intent(in) :: f1(-6:6),f2(-6:6)
    real(dp) :: respdf
    real(dp) :: restmp
    integer :: i,j


!    restmp = 0
!    do i = -nf_pdf, nf_pdf
!       do j = 2, nf_pdf, 2
!          if (i .eq. -j .or. i == 0) cycle
!          restmp = restmp + f1(j)*f2(i)
!       enddo
!    enddo
!    respdf = restmp * res_j(5) !-- u g
!
!    restmp = 0
!    do i = -nf_pdf, nf_pdf
!       do j = 1, 3, 2
!          if (i .eq. j .or. i==0) cycle
!          restmp = restmp + f1(-j)*f2(i)
!       enddo
!    enddo
!    respdf = respdf + restmp * res_j(6) !-- dx g

    respdf = zero
    respdf = respdf + res_j(5)*(f2(-5)+f2(-4)+f2(-3)+f2(-2)+f2(1)+f2(2)+f2(3)+f2(4)+f2(5))*f1(2)    ! u g
    respdf = respdf + res_j(5)*(f2(-5)+f2(-4)+f2(-2)+f2(-1)+f2(1)+f2(2)+f2(3)+f2(4)+f2(5))*f1(4)    ! c g
    respdf = respdf + res_j(6)*(f2(-5)+f2(-4)+f2(-3)+f2(-2)+f2(-1)+f2(1)+f2(3)+f2(4)+f2(5))*f1(-1)  ! dx g
    respdf = respdf + res_j(6)*(f2(-5)+f2(-4)+f2(-3)+f2(-2)+f2(-1)+f2(1)+f2(2)+f2(3)+f2(5))*f1(-3)  ! sx g

  end function s_lumi_qg

  !-- singlet, limits on leg 1
  function s_lumi_xj(res,f1,f2) result(respdf)
    real(dp), intent(in) :: res(2),f1(-6:6),f2(-6:6)
    real(dp) :: respdf
    real(dp) :: restmp
    integer :: i,j

    respdf = zero
    respdf = respdf + res(1)*(f1(-5)*f2(-1) + f1(-4)*f2(-1) + f1(-3)*f2(-1) + f1(-2)*f2(-1) + f1(-1)*f2(-1) + f1(1)*f2(-1) + f1(3)*f2(-1) + f1(4)*f2(-1) + f1(5)*f2(-1)+ &
         f1(-5)*f2(-3) + f1(-4)*f2(-3) + f1(-3)*f2(-3) + f1(-2)*f2(-3) + f1(-1)*f2(-3) + f1(1)*f2(-3) + f1(2)*f2(-3) + f1(3)*f2(-3) +f1(5)*f2(-3))
    respdf = respdf + res(2)*(f1(-5)*f2(2) + f1(-4)*f2(2) + f1(-3)*f2(2) +f1(-2)*f2(2) + f1(1)*f2(2) + f1(2)*f2(2) + f1(3)*f2(2) + f1(4)*f2(2) + f1(5)*f2(2)+ &
         f1(-5)*f2(4) + f1(-4)*f2(4)  + f1(-2)*f2(4) + f1(-1)*f2(4) + f1(1)*f2(4) + f1(2)*f2(4) + f1(3)*f2(4) + f1(4)*f2(4) +f1(5)*f2(4))

!    restmp = 0
!    do i = -nf_pdf, nf_pdf
!       do j = 2, nf_pdf, 2
!          if (i .eq. -j .or. i==0) cycle
!          restmp = restmp + f1(i)*f2(j)
!       enddo
!    enddo
!    respdf =  restmp * res(2) !-- ux d
!
!    restmp = 0
!    do i = -nf_pdf, nf_pdf
!       do j = 1, 3, 2
!          if (i .eq. j .or. i==0) cycle
!          restmp = restmp + f1(i)*f2(-j)
!       enddo
!    enddo
!    respdf = respdf + restmp * res(1) !-- u dx
    

    ! print *, 'respdf', respdf
    ! pause

  end function s_lumi_xj

  !-- singlet, limits on leg 1, identical final state
  function s_lumi_xj_id(res,f1,f2) result(respdf)
    real(dp), intent(in) :: res(2),f1(-6:6),f2(-6:6)
    real(dp) :: respdf

    respdf = sum(f1(-2:-nf_pdf:-2)*f2(-1:-3:-2))*res(1) + &
         sum(f1( 1: 3: 2)*f2( 2: nf_pdf: 2))*res(2)

  end function s_lumi_xj_id

  !-- singlet, limits on leg 2
  function s_lumi_ix(res,f1,f2) result(respdf)
    real(dp), intent(in) :: res(2),f1(-6:6),f2(-6:6)
    real(dp) :: respdf
    real(dp) :: restmp
    integer :: i,j


    respdf = zero
    respdf = respdf + res(2)*(f2(-5)*f1(-1) + f2(-4)*f1(-1) + f2(-3)*f1(-1) + f2(-2)*f1(-1) + f2(-1)*f1(-1) + f2(1)*f1(-1) + f2(3)*f1(-1) + f2(4)*f1(-1) + f2(5)*f1(-1) + &
         f2(-5)*f1(-3) + f2(-4)*f1(-3) + f2(-3)*f1(-3) + f2(-2)*f1(-3) + f2(-1)*f1(-3) + f2(1)*f1(-3) + f2(2)*f1(-3) + f2(3)*f1(-3) + f2(5)*f1(-3))
    respdf = respdf + res(1)*(f2(-5)*f1(2) + f2(-4)*f1(2) + f2(-3)*f1(2) +f2(-2)*f1(2) + f2(1)*f1(2) + f2(2)*f1(2) + f2(3)*f1(2) + f2(4)*f1(2) + f2(5)*f1(2)+ &
         f2(-5)*f1(4) + f2(-4)*f1(4)  + f2(-2)*f1(4) + f2(-1)*f1(4) + f2(1)*f1(4) + f2(2)*f1(4) + f2(3)*f1(4) + f2(4)*f1(4) + f2(5)*f1(4))


!    restmp = 0
!    do i = -nf_pdf, nf_pdf
!       do j = 2, nf_pdf, 2
!          if (i .eq. -j .or. i==0) cycle
!          restmp = restmp + f1(j)*f2(i)
!       enddo
!    enddo
!    respdf =  restmp * res(1) !-- u dx
!
!    restmp = 0
!    do i = -nf_pdf, nf_pdf
!       do j = 1, 3, 2
!          if (i .eq. j .or. i==0) cycle
!          restmp = restmp + f1(-j)*f2(i)
!       enddo
!    enddo
!    respdf = respdf + restmp * res(2) !-- dx u
    
  end function s_lumi_ix

  !-- singlet, limits on leg 2, identical final state
  function s_lumi_ix_id(res,f1,f2) result(respdf)
    real(dp), intent(in) :: res(2),f1(-6:6),f2(-6:6)
    real(dp) :: respdf


    respdf = sum(f1( 2: nf_pdf: 2)*f2( 1: 3: 2))*res(1) + &
         sum(f1(-1:-3:-2)*f2(-2:-nf_pdf:-2))*res(2)

    
  end function s_lumi_ix_id

  !--

  function gg_lumi_nnlo(res_jj,f1,f2) result(respdf)
    real(dp), intent(in) :: res_jj(1),f1(-6:6),f2(-6:6)
    real(dp) :: respdf

    !-- g g > d ux ve ep, g g > s cx ve ep
    respdf = res_jj(1)*(nup*f1(0)*f2(0))

  end function gg_lumi_nnlo

  function gg_lumi_qxg(res_j,f1,f2) result(respdf)
    real(dp), intent(in) :: res_j(6),f1(-6:6),f2(-6:6)
    real(dp) :: respdf

    !-- g g > d ux em ep, g g > s cx em ep --> dx g > ux ve ep, sx g > cx ve ep 
    respdf = res_j(6)*(nup*f1(0)*f2(0))

  end function gg_lumi_qxg

  function gg_lumi_qg(res_j,f1,f2) result(respdf)
    real(dp), intent(in) :: res_j(6),f1(-6:6),f2(-6:6)
    real(dp) :: respdf

    !-- g g > d ux ve ep, g g > s cx ve ep -- > u g > d ve ep, c g > s ve ep
    respdf = res_j(5)*(nup*f1(0)*f2(0))

  end function gg_lumi_qg

  function gg_lumi_gqx(res_j,f1,f2) result(respdf)
    real(dp), intent(in) :: res_j(6),f1(-6:6),f2(-6:6)
    real(dp) :: respdf

    !-- g g > d ux ve ep, g g > s cx ve ep -- > g dx > ux ve ep, d sx > cx ve ep
    respdf = res_j(3)*(nup*f1(0)*f2(0))

  end function gg_lumi_gqx

  function gg_lumi_gq(res_j,f1,f2) result(respdf)
    real(dp), intent(in) :: res_j(6),f1(-6:6),f2(-6:6)
    real(dp) :: respdf

    !-- g g > d ux ve ep, g g > s cx ve ep > g u > d ve ep, g c > s ve ep
    respdf = res_j(2)*(nup*f1(0)*f2(0))

  end function gg_lumi_gq

  function gg_lumi_qqx(res,f1,f2) result(respdf)
    real(dp), intent(in) :: res(2),f1(-6:6),f2(-6:6)
    real(dp) :: respdf

    !-- g g > d ux ve ep, g g > s cx ve ep --> u dx > ve ep, c sx --> ve ep
    respdf = res(1)*(nup*f1(0)*f2(0))


  end function gg_lumi_qqx

  function gg_lumi_qxq(res,f1,f2) result(respdf)
    real(dp), intent(in) :: res(2),f1(-6:6),f2(-6:6)
    real(dp) :: respdf

    !-- g g > d ux em ep, g g > s cx em ep -- > dx u > ve ep, sx c > ve ep
    respdf = res(2)*(nup*f1(0)*f2(0))

  end function gg_lumi_qxq
  
  !-- 

  function qg_lumi_rv(res,f1,f2) result(respdf)
    real(dp), intent(in) :: res(2),f1(-6:6),f2(-6:6)
    real(dp) :: respdf


    !--u g > d ve ep, c g > s ve ep
    respdf = res(1)*sum(f1(2:nf_pdf:2))*f2(0)
    
    !-- dx g > ux ve ep, sx g > cx ve ep
    respdf = respdf + res(2)*sum(f1(-1:-3:-2))*f2(0)

  end function qg_lumi_rv

  function qg_lumi_rv_lo(res,f1,f2) result(respdf)
    real(dp), intent(in) :: res(2),f1(-6:6),f2(-6:6)
    real(dp) :: respdf

    !-- u g > d ve ep, c g > s ve ep
    respdf = res(1)*sum(f1(2:nf_pdf:2))*f2(0)
    
    !-- dx g > ux ve ep, sx g > cx ve ep
    respdf = respdf + res(2)*sum(f1(-1:-3:-2))*f2(0)

    ! print *, 'respdf', respdf
    ! pause
       
  end function qg_lumi_rv_lo

  function gq_lumi_rv(res,f1,f2) result(respdf)
    real(dp), intent(in) :: res(2),f1(-6:6),f2(-6:6)
    real(dp) :: respdf


    !-- g u > d ve ep, g c > s ve ep
    respdf = res(1)*f1(0)*sum(f2(2:nf_pdf:2))

    !-- g dx > ux ve ep, g sx > cx ve ep
    respdf = respdf + res(2)*f1(0)*sum(f2(-1:-3:-2))

    ! print *, 'respdf', respdf
    ! pause
       
  end function gq_lumi_rv

  function gq_lumi_rv_lo(res,f1,f2) result(respdf)
    real(dp), intent(in) :: res(2),f1(-6:6),f2(-6:6)
    real(dp) :: respdf


    !-- g u > d ve ep, g c > s ve ep --> ux u
    respdf = res(2)*f1(0)*sum(f2(2:nf_pdf:2))

    !-- g dx > ux ve ep, g sx > cx ve ep
    respdf = respdf + res(1)*f1(0)*sum(f2(-1:-3:-2))

    ! print *, 'respdf', respdf
    ! pause
       
  end function gq_lumi_rv_lo
  
end module mod_auxlumi_wp

