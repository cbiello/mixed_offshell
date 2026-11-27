module mod_kinematicsNLO
  use mod_types; use mod_consts_dp; use common_def
  use mod_process
  use mod_auxfunctions; use mod_procfunctions
  use mod_kinematics_lept
  use mod_kinematics_gen
  implicit none
  private

  public :: kinematics_nlo

contains 

  !-- we factor out an extra 1/(8*pisq) from the phase space
  !-- which multiplied by gs**2 in the ME reconstructs asontwopi
  subroutine kinematics_nlo(tau,ylab,yr,HardProc,SoftLim,Col1Lim,Col2Lim) 
    real(dp), intent(in) :: tau,ylab,yr(6)
    type(KinConfig)   :: HardProc,SoftLim,Col1Lim,Col2Lim
    real(dp) :: xi1,xi2,xi1xi2jac,sqrts,slocal,weight
    real(dp) :: p1(4),p2(4),p4(4),p124(4),p4coll(4)  
    real(dp) :: x1,x2,cos4,sin4,phi4
    real(dp) :: s12,s14,s24,z
    real(dp) :: mV,mV2,kallenF
    real(dp) :: flux

    call initialize_config(HardProc)
    call initialize_config(SoftLim)
    call initialize_config(Col1Lim)
    call initialize_config(Col2Lim)
    
    x1 =   yr(1) 
    x2 =   yr(2)

    cos4 = one - two*x2
    sin4 = sqrt(one-cos4**2)
    phi4 = twopi*yr(3) 

    !   full kinematics 

    xi1 = sqrt(tau/(one-x1)) * exp(+ylab) * sqrt((one-x1*x2)/(one-x1+x1*x2))
    xi2 = sqrt(tau/(one-x1)) * exp(-ylab) / sqrt((one-x1*x2)/(one-x1+x1*x2))
    xi1xi2jac = one/(one-x1)

    if (xi1.gt.one .or. xi2.gt.one) then

       HardProc%flag = .true.

    else

       slocal = sh*xi1*xi2
       sqrts = sqrt(slocal)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       
       p4 = half*sqrts*x1 *(/ one, sin4*cos(phi4), sin4*sin(phi4), cos4 /)

       p124 = p1 + p2 - p4

       MV2 = slocal*(one-x1)
       MV  =  sqrt(MV2)
       call get_lept_mom(MV,p124,yr(4:6),kallenF,HardProc)
       
       HardProc%AmpMom(:,1) = p1
       HardProc%AmpMom(:,2) = p2
       HardProc%AmpMom(:,3) = p4
 
       HardProc%KinInv = (/mv2/)

       HardProc%PartFrac=(/xi1,xi2/)

       weight = one/8.0_dp/pi*kallenF * ( x1 ) * ( slocal/two)
       weight = weight/(two*slocal) ! -- flux
       weight = weight * xi1xi2jac
       HardProc%wgt = weight

    endif

       
    !  soft kinematics 

    xi1 = sqrt(tau) * exp(+ylab)
    xi2 = sqrt(tau) * exp(-ylab)
    xi1xi2jac = one

    if (xi1.gt.one .or. xi2.gt.one) then

       SoftLim%flag = .true.

    else

       slocal = sh*xi1*xi2
       sqrts = sqrt(slocal)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 

       p4 = half*sqrts*x1 *(/ one, sin4*cos(phi4), sin4*sin(phi4), cos4 /)

       p124 = p1 + p2

       MV2 = slocal
       MV  = sqrts
       call get_lept_mom(MV,p124,yr(4:6),kallenF,SoftLim)

       SoftLim%LimMom(:,1) = p1
       SoftLim%LimMom(:,2) = p2
       SoftLim%LimMom(:,3) = p4

       SoftLim%AmpMom(:,1) = p1
       SoftLim%AmpMom(:,2) = p2

       s12 = mv2
       s14 = s12*x1*x2
       s24 = s12*x1*(one-x2)

       SoftLim%KinInv = (/s12,s14,s24/)

       SoftLim%PartFrac=(/xi1,xi2/)

       weight = one/8.0_dp/pi*kallenF * ( x1 ) * ( slocal/two)
       weight = weight/(two*slocal) ! -- flux
       weight = weight * xi1xi2jac
       SoftLim%wgt = weight

    endif

    !-- 1-collinear kinematics 

    xi1 = sqrt(tau/(one-x1)) * exp(+ylab) / sqrt(one-x1)
    xi2 = sqrt(tau/(one-x1)) * exp(-ylab) * sqrt(one-x1)
    xi1xi2jac = one/(one-x1)

    if (xi1.gt.one .or. xi2.gt.one) then

       Col1Lim%flag = .true.

    else

       slocal = sh*xi1*xi2
       sqrts = sqrt(slocal)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 

       p4 = half*sqrts*x1 *(/ one, sin4*cos(phi4), sin4*sin(phi4), cos4 /)

       p4coll = half*sqrts*x1 *(/ one,zero,zero,one /)

       p124 = p1 + p2 - p4coll

       MV2 = slocal*(one-x1)
       MV  =  sqrt(MV2)
       call get_lept_mom(MV,p124,yr(4:6),kallenF,Col1Lim)

       Col1Lim%LimMom(:,1) = p1-p4coll
       Col1Lim%LimMom(:,2) = p2
       Col1Lim%LimMom(:,3) = p4

       Col1Lim%AmpMom(:,1) = p1-p4coll
       Col1Lim%AmpMom(:,2) = p2

       Col1Lim%PartFrac = (/xi1,xi2/)
       
       s14 = slocal*x1*x2
       z = one-x1
       Col1Lim%KinInv = (/z,s14,mv2/)

       weight = one/8.0_dp/pi*kallenF * ( x1 ) * ( slocal/two)
       weight = weight/(two*slocal) ! -- flux
       weight = weight * xi1xi2jac
       Col1Lim%wgt = weight

    endif

       
    !-- 2-collinear kinematics 

    xi1 = sqrt(tau/(one-x1)) * exp(+ylab) * sqrt(one-x1)
    xi2 = sqrt(tau/(one-x1)) * exp(-ylab) / sqrt(one-x1)
    xi1xi2jac = one/(one-x1)
    if (xi1.gt.one .or. xi2.gt.one) then

       Col2Lim%flag = .true.

    else

       slocal = sh*xi1*xi2
       sqrts = sqrt(slocal)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 

       p4 = half*sqrts*x1*(/ one, sin4*cos(phi4), sin4*sin(phi4), cos4 /)

       p4coll = half*sqrts*x1*(/ one,zero,zero,-one /)
       
       p124 = p1 + p2 - p4coll
       
       MV2 = slocal*(one-x1)
       MV  =  sqrt(MV2)
       call get_lept_mom(MV,p124,yr(4:6),kallenF,Col2Lim)

       Col2Lim%LimMom(:,1) = p1
       Col2Lim%LimMom(:,2) = p2-p4coll
       Col2Lim%LimMom(:,3) = p4

       Col2Lim%AmpMom(:,1) = p1
       Col2Lim%AmpMom(:,2) = p2-p4coll

       Col2Lim%PartFrac = (/xi1,xi2/)

       s24 = slocal*x1*(one-x2)
       z = one-x1
       Col2Lim%KinInv = (/z,s24,mv2/)

       weight = one/8.0_dp/pi*kallenF * ( x1 ) * ( slocal/two)
       weight = weight/(two*slocal) ! -- flux
       weight = weight*xi1xi2jac
       Col2Lim%wgt = weight

    endif


      
  end subroutine kinematics_nlo


end module mod_kinematicsNLO


