module mod_my_vegas
  use mod_types
  use mod_vegas_parms
  use mod_histo
  implicit none
  private

  public :: vegas_integrate

contains

  !-- wrapper for go_vegas_gen and grids
  subroutine vegas_integrate(ndim_in,integrand,vg_result,vg_error,vg_chi2)
    integer, intent(in) :: ndim_in
    real(dp15), intent(out) :: vg_result,vg_error,vg_chi2
    interface
       function integrand(xx,ff,weight)
         use mod_types
         integer    :: integrand
         real(dp15) :: xx(30),ff(1),weight
       end function integrand
    end interface

    if (writegrid) then
       write (6,*) '# writing out vegas grid'
       call go_vegas_gen(ndim_in,integrand,vg_result,vg_error,vg_chi2,grid_out=gridfile)
    elseif (readgrid) then
       write (6,*) '# reading in vegas grid'
       call go_vegas_gen(ndim_in,integrand,vg_result,vg_error,vg_chi2,grid_in=gridfile)
    else
       call go_vegas_gen(ndim_in,integrand,vg_result,vg_error,vg_chi2)
    endif

  end subroutine vegas_integrate

  subroutine go_vegas_gen(ndim,integrand,res,err,chi2,grid_in,grid_out,io_grid)
    integer, intent(in)   :: ndim
    real(vp), intent(out) :: res,err,chi2
    integer :: idum,iogrid
    character(*), intent(in), optional :: grid_in,grid_out
    integer, intent(in), optional      :: io_grid
    integer, parameter :: verbose = 0
    !
    interface
       function integrand(xx,ff,weight)
         use mod_types
         integer    :: integrand
         real(dp15) :: xx(30),ff(1),weight
       end function integrand
    end interface
    
    if (present(io_grid)) then
       iogrid = io_grid
    else
       iogrid = 66
    endif

    call init_idum(seed,idum)
    
    if (present(grid_in)) then

       !-- read from the grid
       nohistos = .false.
       call clear_histo()
       !
       call my_vegas_lepage(integrand,idum,1,ndim,VegasNc1,VegasIt1,verbose,res,err,chi2,grid_in=grid_in,io_grid=iogrid)

    else

       !-- warmup
       nohistos = .true.
       call clear_histo()
       !
       if(present(grid_out)) then
          call my_vegas_lepage(integrand,idum,0,ndim,VegasNc0,VegasIt0,verbose,res,err,chi2,grid_out=grid_out,io_grid=iogrid)
       else
          call my_vegas_lepage(integrand,idum,0,ndim,VegasNc0,VegasIt0,verbose,res,err,chi2)
       endif

       !-- full run
       call init_idum(seed,idum) !-- re-initialise random number to have identical results when reading-in grid
       nohistos = .false.
       call clear_histo()
       !
       call my_vegas_lepage(integrand,idum,1,ndim,VegasNc1,VegasIt1,verbose,res,err,chi2)
       
    endif

    return
    
  end subroutine go_vegas_gen

  subroutine init_idum(seed,idum)
    integer, intent(in)  :: seed
    integer, intent(out) :: idum

    if (seed.eq.0) then
       idum = -19
    else
       idum = -seed
    endif

  end subroutine init_idum
    
  
  !----------------------------------------------------------------------------------------------------------
  
  !-- vegas algorithm below
  
  subroutine my_vegas_lepage(integrand,idum,init,ndim,ncall,itmx,verb,res,err,chi2a,grid_in,grid_out,io_grid)
    integer, intent(in)   :: init,ndim,ncall,itmx,verb 
    real(vp), intent(out) :: res,err,chi2a
    character(*), intent(in), optional :: grid_in,grid_out
    integer, intent(in), optional :: io_grid
    integer :: iogrid,idum !-- idum -> random seed, must start negative
    logical :: gridexists
    !
    interface
       function integrand(xx,ff,weight)
         use mod_types
         integer    :: integrand
         real(dp15) :: xx(30),ff(1),weight
       end function integrand
    end interface
    integer :: ffdum
    real(dp15) :: ff(1)
    !
    real(vp), parameter :: alph=1.5_vp,tiny=1.0e-30_vp
    integer,  parameter :: mxdim=25,ndmx=50
    !
    integer, save :: i,it,j,k,mds,nd,ndo,ng,npg
    integer, dimension(mxdim), save :: ia,kg
    real(vp), save :: calls,dv2g,dxg,f,f2,f2b,fb,rc,ti,tsi,wgt,xjac,xn,xnd,xo
    real(vp), dimension(ndmx,mxdim), save :: d,di,xi
    real(vp), dimension(mxdim), save :: dt,x
    real(vp), dimension(ndmx), save :: r,xin
    real(vp), save :: schi,si,swgt
    real(vp) :: wtmax = 1.0_vp !-- to match MCFM printout style

    if (init <= 0) then !-- normal entry
       mds=1 !-- mds=0 -> no stratified sampling
       ndo=1
       xi(1,:)=1.0_vp
    endif

    if (init <= 1) then !-- inherit grid, don't inherit result
       si=0.0_vp
       swgt=0.0_vp
       schi=0.0_vp
       
       !-- readin the grid
       if (present(grid_in)) then
          write(6,*) '************* Reading in the grid from ', adjustl(trim(grid_in)), &
               ' **************'
          
          if(present(io_grid)) then
             iogrid = io_grid
          else
             iogrid = 66
          endif
          
          inquire(file=grid_in,exist=gridexists)
          if (.not. gridexists) then
             write (6,*) 'grid file ',trim(grid_in), ' does not exist'
             stop
          endif

          open (unit=iogrid,file=grid_in,status='old')
          read(iogrid,203) mds,ndo,nd
          do j = 1,ndim
             read(iogrid,203) (xi(i,j),i=1,nd)
          enddo
          close (iogrid)
          
       endif
       !-- end readin the grid
       
    endif

    if (init <= 2) then !-- inherit grid and result
       nd=ndmx
       ng=1

       if (mds /= 0) then !-- stratified sampling
          ng=int((ncall/2.0_vp+0.25_vp)**(1.0_vp/ndim))
          mds=1

          if ((2*ng-ndmx) >=0) then
             mds=-1
             npg=ng/ndmx+1
             nd=ng/npg
             ng=npg*nd
          endif

       endif

       k=ng**ndim
       npg=max(ncall/k,2)
       calls=real(npg,kind=vp)*real(k,kind=vp)
       dxg=1.0_vp/ng
       dv2g=(calls*dxg**ndim)**2/npg/npg/(npg-1.0_vp)
       xnd=nd
       dxg=dxg*xnd
       xjac=1.0_vp/calls

       !-- do rebinning if necessary
       if (nd /= ndo) then
          r(1:max(nd,ndo))=1.0_vp
          do j=1,ndim
             call rebin(ndo/xnd,nd,r,xin,xi(:,j))
          enddo
          ndo=nd
       endif

       it = 0
       if (verb >=0) write(6,200) ndim,calls,it,itmx,verb,&
            alph,mds,nd
       call flush(6)

    endif

    do it=1,itmx !-- iteration loop
       ti=0.0_vp
       tsi=0.0_vp
       kg(:)=1
       d(1:nd,:)=0.0_vp
       di(1:nd,:)=0.0_vp
       iterate: do
          fb=0.0_vp
          f2b=0.0_vp
          do k=1,npg
             wgt=xjac
             do j=1,ndim
                xn=(kg(j)-my_ran1(idum))*dxg+1.0_vp
                ia(j)=max(min(int(xn),ndmx),1)
                if (ia(j) > 1) then
                   xo=xi(ia(j),j)-xi(ia(j)-1,j)
                   rc=xi(ia(j)-1,j)+(xn-ia(j))*xo
                else
                   xo=xi(ia(j),j)
                   rc=(xn-ia(j))*xo
                end if
                x(j)=rc
                wgt=wgt*xo*xnd
             end do
             ffdum=integrand(x(1:ndim),ff,wgt)
             f=wgt*ff(1)
             f2=f*f
             fb=fb+f
             f2b=f2b+f2
             do j=1,ndim
                di(ia(j),j)=di(ia(j),j)+f
                if (mds >=0) d(ia(j),j)=d(ia(j),j)+f2
             end do
          end do
          f2b=sqrt(f2b*npg)
          f2b=(f2b-fb)*(f2b+fb)
          if(f2b <= 0.0_vp) f2b=tiny
          ti=ti+fb
          tsi=tsi+f2b
          if (mds < 0) then !-- stratified sampling
             do j=1,ndim
                d(ia(j),j)=d(ia(j),j)+f2b
             end do
          end if
          do k=ndim,1,-1
             kg(k)=mod(kg(k),ng)+1
             if (kg(k) /= 1) cycle iterate
          end do
          exit iterate
       end do iterate
       
       !-- final results for iteration
       tsi=tsi*dv2g
       wgt=1.0_vp/tsi
       si=si+wgt*ti
       schi=schi+wgt*ti**2
       swgt=swgt+wgt
       res=si/swgt
       chi2a=max((schi-si*res)/(it-0.999_vp),0.0_vp)
       err=sqrt(1.0_vp/swgt)
       tsi=sqrt(tsi)

       if (verb >=0 ) then
          !-- it  -> iteration
          !-- ti  -> iteration integral
          !-- res -> accum. imtegral
          !-- tsi -> iteration std. dev.
          !-- err -> accum. std. dev.
          !-- wtmax -> max weight, fixed to 1, to match MCFM printount
          !-- chi2a -> chi^2
          write(6,201) it,ti,res,tsi,err,wtmax,chi2a
          !-- histogramming
          call prepare_and_write(init,it,ti,res,tsi,err,chi2a)
          !-- end histogramming
          call flush(6)
          if (verb /= 0) then
             do j=1,ndim
                write(6,202) j,(xi(i,j),di(i,j),&
                     i=1,nd) !-- same as MCFM
             enddo
          end if
       end if
       do j=1,ndim !-- refine the grid, damped by alph
          xo=d(1,j)
          xn=d(2,j)
          d(1,j)=(xo+xn)/2.0_vp
          dt(j)=d(1,j)
          do i=2,nd-1
             rc=xo+xn
             xo=xn
             xn=d(i+1,j)
             d(i,j)=(rc+xn)/3.0_vp
             dt(j)=dt(j)+d(i,j)
          end do
          d(nd,j)=(xo+xn)/2.0_vp
          dt(j)=dt(j)+d(nd,j)
       end do
       where (d(1:nd,:) < tiny) d(1:nd,:) = tiny
       do j=1,ndim
          r(1:nd)=((1.0_vp-d(1:nd,j)/dt(j))/(log(dt(j))-log(d(1:nd,j))))**alph
          rc=sum(r(1:nd))
          call rebin(rc/xnd,nd,r,xin,xi(:,j))
       end do
    end do !-- end iteration loop
    
    !-- writeout the grid
    if (present(grid_out)) then
       write(6,*) '************* Writing out the grid to ', adjustl(trim(grid_out)), &
            ' **************'

       if(present(io_grid)) then
          iogrid = io_grid
       else
          iogrid = 66
       endif
       
       open (unit=iogrid,file=grid_out,status='unknown')
       write(iogrid,203) mds,ndo,nd
       do j = 1,ndim
          write(iogrid,203) (xi(i,j),i=1,nd)
       enddo
       close (iogrid)
       
    endif
    !-- end writeout the grid

200 format(/,' Input parameters for vegas: ndim =',i3, &
         '  ncall =',f10.0/28x,' it =',i5,'  itmx =',i5/28x, &
         ' verb =',i3,'  alph =',f5.2/28x,' mds = ',i3,'  nd =',i4/28x)

201 format(/'************* Integration by Vegas (iteration ',i3, &
         ') **************' / '*',63x,'*'/, &
         '*  integral  = ',g14.8,2x, &
         ' accum. integral = ',g14.8,'*'/, &
         '*  std. dev. = ',g14.8,2x, &
         ' accum. std. dev = ',g14.8,'*'/, &
         '*   max. wt. = ',g14.6,35x,'*'/,'*',63x,'*'/, &
         '**************   chi**2/iteration = ', &
         g10.4,'   ****************' /)

202 format(1X,' data for axis',i2,/,' ',6x,'x',7x,'  delt i ', &
         2x,'conv','ce   ',11x,'x',7x,'  delt i ',2x,'conv','ce  ' &
         ,11x,'x',7x,'   delt i ',2x,'conv','CE  ',/, &
         (1X,' ',3g12.4,5x,3g12.4,5x,3g12.4))

203 format(/(5z16))
    
  end subroutine my_vegas_lepage

  subroutine rebin(rc,nd,r,xin,xi)
    real(vp), intent(in) :: rc
    integer, intent(in) :: nd
    real(vp), dimension(:), intent(in) :: r
    real(vp), dimension(:), intent(out) :: xin
    real(vp), dimension(:), intent(inout) :: xi
    integer :: i,k
    real(vp) :: dr,xn,xo

    k=0
    xo=0.0_vp
    dr=0.0_vp
    
    do i=1,nd-1
       do
          if (rc <= dr) exit
          k=k+1
          dr=dr+r(k)
       end do
       if (k > 1) xo=xi(k-1)
       xn=xi(k)
       dr=dr-rc
       xin(i)=xn-(xn-xo)*dr/r(k)
    enddo
    xi(1:nd-1)=xin(1:nd-1)
    xi(nd)=1.0_vp
    
  end subroutine rebin

  function my_ran1(idum)
    real(vp) :: my_ran1
    integer :: idum
    integer, parameter :: ia=16807
    integer, parameter :: im=2147483647
    integer, parameter :: iq=127773
    integer, parameter :: ir=2836
    integer, parameter :: ntab=32
    integer, parameter :: ndiv=int(1+(im-1._vp)/ntab) !-- not checked
    real(vp), parameter :: am=1.0_vp/im
    real(vp), parameter :: eps=3e-16_vp
    real(vp), parameter :: rnmx=1.0_vp-eps
    integer, save :: iv(ntab) = 0
    integer, save :: iy = 0
    integer :: j,k

    if (idum.le.0 .or. iy.eq.0) then
       idum = max(-idum,1)
       do j=ntab+8,1,-1
          k=idum/iq
          idum=ia*(idum-k*iq)-ir*k
          if (idum.lt.0) idum=idum+im
          if (j.le.ntab) iv(j)=idum
       end do
       iy=iv(1)
    end if

    k=idum/iq
    idum=ia*(idum-k*iq)-ir*k
    if (idum.lt.0) idum=idum+im
    j=1+iy/ndiv
    iy=iv(j)
    iv(j)=idum
    my_ran1=min(am*real(iy,kind=vp),rnmx)

    return

  end function my_ran1

end module mod_my_vegas
