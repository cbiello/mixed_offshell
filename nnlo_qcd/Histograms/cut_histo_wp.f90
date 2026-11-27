!-- observables and cumulative observables. Format must be !-- xxx = number, with these spaces
!-- Nobs = 6
!-- Nnoncum = 5
!-- cutname = 'inclusive' !-- maximum 50 characters

!!!!
! Just a copy of DY Z 

  subroutine init_histo()
    integer :: nhisto

    nhisto = 1

    Histo(nhisto)%obs_code="rate"
    Histo(nhisto)%LowLimit=zero
    Histo(nhisto)%UpLimit=one
    Histo(nhisto)%Step=one

    nhisto = nhisto + 1
    Histo(nhisto)%obs_code="q2"
    Histo(nhisto)%LowLimit=0.0_dp
    Histo(nhisto)%UpLimit=350.0_dp
    Histo(nhisto)%Step=5.0_dp

    nhisto = nhisto + 1
    Histo(nhisto)%obs_code="pT_Z"
    Histo(nhisto)%LowLimit=zero
    Histo(nhisto)%UpLimit=350.0_dp
    Histo(nhisto)%Step=5.0_dp

    nhisto = nhisto + 1
    Histo(nhisto)%obs_code="y_Z"
    Histo(nhisto)%LowLimit=-5.0_dp
    Histo(nhisto)%UpLimit=5.0_dp
    Histo(nhisto)%Step=0.1_dp

    nhisto = nhisto + 1
    Histo(nhisto)%obs_code="eta_l1"
    Histo(nhisto)%LowLimit=-3.0_dp
    Histo(nhisto)%UpLimit=3.0_dp
    Histo(nhisto)%Step=0.1_dp

    !-- cumulative below

    nhisto = nhisto + 1
    Histo(nhisto)%obs_code="c_pT_Z"
    Histo(nhisto)%LowLimit=0
    Histo(nhisto)%UpLimit=350.0_dp
    Histo(nhisto)%Step=5.0_dp

    if (nhisto.ne.nobs) then
       print *, 'uninitialized histograms'
       stop
    endif

    Histo(:)%Nbins = (Histo(:)%UpLimit - Histo(:)%LowLimit)/Histo(:)%Step

    return

  end subroutine init_histo

  subroutine cut_histo(xa,xb,npart,pobs,makecut,binarray,weightarray,mu)
    integer, intent(in) :: npart
    real(dp), intent(in) :: xa,xb,pobs(4,npart)
    real(dp), intent(out) :: weightarray(nobs),mu
    integer, intent(out) :: binarray(nobs)
    logical, intent(out) :: makecut
    real(dp) :: pobs_lab(4,npart)
    real(dp) :: obs(nobs)
    !-- process-dependent part below
    real(dp) :: pV(4),pl1(4),pl2(4),mV,yV,ptlept1,ptlept2,dphi,ptV,etalept1,etalept2
    real(dp)  :: ptj1,ptj2
    real(dp) :: pjet(4,3),pin(4,2)
    integer  :: njet

    makecut = .true.
    pobs_lab(:,:) = zero
    binarray(:) = 0
    weightarray(:) = one
    
    call boostz(xa,xb,npart,pobs,pobs_lab)
   ! debug
    !!Hjj rate
    !pin(1:4,1:2) = pobs_lab(1:4,3:npart-2)
    !call ktjet(pin,npart-4,pjet,njet)
    !!if (njet  .lt. 2 ) return

    mu = sqrt(musq_in)
    pl1 = pobs_lab(:,npart-1) 
    pl2 = pobs_lab(:,npart)
    pV = pobs_lab(:,npart-1) + pobs_lab(:,npart)

    !---------------------------------------------------
    
    !-- process-dependent part below

    mV = sqrt(scr(pV,pV))
    

    yV = get_y(pV)
    ptV = get_pt(pV)
    ptlept1 = get_pt(pl1)
    ptlept2 = get_pt(pl2)
    etalept1 = get_y(pl1)
    
    etalept2 = get_y(pl2)


    if ( abs(etalept1) .gt. 5.0_dp .or. abs(etalept2) .gt. 5.0_dp ) return
    !-----------------------------------------------------------------
    !-- fill observables

    !-- standard
    Obs(1) = 0.5_dp
    Obs(2) = mV
    Obs(3) = ptV
    Obs(4) = yV
    Obs(5) = abs(etalept1)

    !-- cumulative
    Obs(6) = ptV

    call find_bin(obs,binarray,weightarray)

    makecut = .false.

    return


  end subroutine cut_histo

  
  subroutine ktjet(pin,  nin, pjet, njet)
    integer, intent(in) :: nin
    real(dp), intent(in) :: pin(4,nin)  
    integer, intent(out) :: njet
    real(dp), intent(out) :: pjet(4,3)
    real(dp) :: pti,dbi,ptj,dbj
    real(dp) :: yi,yj,dy,d_cosphi,dphi,r12sq,dij
    real(dp) :: dmin,rsq
    real(dp) :: prehard(4,nin),prejet(4,nin)
    integer :: nh,nj,ijmin(2)
    integer :: i,j,le
    real(dp), parameter :: ktclass = -one   ! antikt
    real(dp), parameter :: ptjetcut = 20.0_dp
    real(dp), parameter :: Rjet = 0.4_dp


    rsq = Rjet**2

    le = 0

    njet = 0
    pjet = zero

    do i=1,nin
       if (pin(1,i).ne.zero.and.(abs(pin(4,i)).ne.abs(pin(1,i)))) then 
          le = le+1
          prehard(:,le) = pin(:,i)
       endif
    enddo

    nj = 0
    nh = le
    prejet = zero

!   recombine

    do while (nh.gt.0) ! recombination loop

       pti = sqrt(prehard(2,nh)**2+prehard(3,nh)**2)
       dbi = pti**(2*ktclass)
       dmin = dbi
       ijmin = (/nh,0/)

       do i = 1, nh-1 ! compute db, dij
          
          pti = sqrt(prehard(2,i)**2+prehard(3,i)**2)
          dbi = pti**(2*ktclass)
       
          if (dbi.lt.dmin) then
             dmin = dbi
             ijmin = (/i,0/)
          endif

          do j = i+1,nh
             ptj = sqrt(prehard(2,j)**2+prehard(3,j)**2)
             dbj = ptj**(2*ktclass)

             yi = half * log((prehard(1,i)+prehard(4,i))/(prehard(1,i)-prehard(4,i)))
             yj = half * log((prehard(1,j)+prehard(4,j))/(prehard(1,j)-prehard(4,j)))

             dy = yi-yj

             ! cos(a-b) = cos(a)cos(b)+sin(a)sin(b)
             d_cosphi = (prehard(2,i)*prehard(2,j)+prehard(3,i)*prehard(3,j))/pti/ptj
             if (d_cosphi.gt.one) then
                dphi = zero
             elseif (d_cosphi.lt.-one) then
                dphi = pi
             else
                dphi = acos(d_cosphi) ! 0<dphi<pi
             endif
             
             r12sq = dy**2 + dphi**2

             dij = min(dbi,dbj) * r12sq/rsq

             if (dij.lt.dmin) then
                dmin = dij
                ijmin = (/min(i,j),max(i,j)/)
             endif

          enddo

       enddo ! compute db, dij

       if (ijmin(2).eq.0) then ! no recombination, 1 new jet candidate
          nj = nj+1
          prejet(:,nj) = prehard(:,ijmin(1))

          do i = ijmin(1),nh-1
             prehard(:,i) = prehard(:,i+1)
          enddo

       else ! recombine, no new jet candidates

          prehard(:,ijmin(1)) = prehard(:,ijmin(1))+prehard(:,ijmin(2))

          do i = ijmin(2),nh-1
             prehard(:,i) = prehard(:,i+1)
          enddo

       endif

       nh = nh - 1

    enddo ! recombination loop
    

    do i = 1, nj
       pti = sqrt(prejet(2,i)**2+prejet(3,i)**2)
       if (pti.gt.ptjetcut) then
          njet = njet + 1
          pjet(:,njet) = prejet(:,i)
       endif

    enddo

    return


  end subroutine ktjet



  
