module mod_aux_sectors
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_process
  use mod_lumi
  use mod_auxfunctions
  use mod_hoppet_tools
  use mod_histo, only: nweights, nohistos
  implicit none
  integer, public, parameter :: ipdf = nweights !-- how many weights at once
  private

  public :: get_prefactor
  public :: get_respdf,get_respdf_vect, get_respdf_gen
  public :: get_respdf_hoppet
  public :: get_respdf_hoppet_gen

  public :: fill_sv_logs

  public :: ns_lumi,aa_lumi
  public :: gq_lumi,qg_lumi
  public :: aq_lumi,qa_lumi
  public :: qq_lumi,qQp_lumi,qQpb_lumi,qQpb_lumi_wp, qQpb_lumi_wm
  public :: ga_lumi,ag_lumi
  
  public :: ns_lumi_splitb
  public :: gq_lumi_splitb,qg_lumi_splitb

  public :: check_ff

  public :: gen_lumi

  public :: multiply_IS_charges

  public :: transition
  
  interface transition
        module procedure transition_res
        module procedure transition_res_sping
    end interface

contains

  subroutine get_prefactor(muR,as_ord,aem_ord,pref)
    real(dp), intent(in)    :: muR
    integer,intent(in)      :: as_ord,aem_ord
    real(dp), intent(out)   :: pref
    real(dp)                :: as,asontwopi
    real(dp15)              :: alphasPDF

    as = alphasPDF(real(muR,kind=dp15))
    asontwopi = as/twopi

    pref = asontwopi**as_ord * (aem/twopi)**aem_ord
    pref = pref * units

  end subroutine get_prefactor


    subroutine get_respdf_gen(as_ord,aem_ord,proc,amp2,respdf)

    integer, intent(in) :: as_ord,aem_ord
    type(KinConfig), intent(in) :: proc
    real(dp), intent(in)        :: amp2(-5:7,-5:7)
    real(dp), intent(out)       :: respdf(ipdf)


    if (ipdf == 1) then
       call get_respdf_one_gen(as_ord,aem_ord,proc,amp2,respdf)
    elseif (ipdf == 3) then
       call get_respdf_three_gen(as_ord,aem_ord,proc,amp2,respdf)
    endif

    contains

      subroutine get_respdf_one_gen(as_ord,aem_ord,proc,amp2,respdf)
        integer, intent(in) :: as_ord,aem_ord
        type(KinConfig), intent(in) :: proc
        real(dp), intent(in)        :: amp2(-5:7,-5:7)
        real(dp), intent(out)       :: respdf(:)
        real(dp) :: f1(-6:7),f2(-6:7),pref
        
        !-- central scale
        call get_prefactor(proc%mur(1),as_ord,aem_ord,pref)
        call get_pdf_qed(proc%PartFrac(1),proc%PartFrac(2),proc%muf(1),f1,f2)
        respdf(1) = pref*gen_lumi(amp2,f1,f2)
    
      end subroutine get_respdf_one_gen

      subroutine get_respdf_three_gen(as_ord,aem_ord,proc,amp2,respdf)
        integer, intent(in) :: as_ord,aem_ord
        type(KinConfig), intent(in) :: proc
        real(dp), intent(in)        :: amp2(-5:7,-5:7)
        real(dp), intent(out)       :: respdf(:)
        real(dp) :: f1(-6:7),f2(-6:7),pref
        
        !-- central scale
        call get_prefactor(proc%mur(1),as_ord,aem_ord,pref)
        call get_pdf_qed(proc%PartFrac(1),proc%PartFrac(2),proc%muf(1),f1,f2)
        respdf(1) = pref*gen_lumi(amp2,f1,f2)

        !-- do not compute the other scales if we're not doing histograms
        if (nohistos) then

           respdf(2:) = zero

        else
               
           !-- half
           call get_prefactor(proc%mur(1)*half,as_ord,aem_ord,pref)
           call get_pdf_qed(proc%PartFrac(1),proc%PartFrac(2),proc%muf(1)*half,f1,f2)
           respdf(2) = pref*gen_lumi(amp2,f1,f2)
           
           !-- twice
           call get_prefactor(proc%mur(1)*two,as_ord,aem_ord,pref)
           call get_pdf_qed(proc%PartFrac(1),proc%PartFrac(2),proc%muf(1)*two,f1,f2)
           respdf(3) = pref*gen_lumi(amp2,f1,f2)

        endif
    
      end subroutine get_respdf_three_gen

    end subroutine get_respdf_gen

  subroutine get_respdf(lumi,as_ord,aem_ord,proc,amp2,respdf)
    interface
       function lumi(res,f1,f2)
         use mod_types
         real(dp), intent(in) :: res(:,:),f1(-6:),f2(-6:)
         real(dp) :: lumi
       end function lumi
    end interface
    integer, intent(in) :: as_ord,aem_ord
    type(KinConfig), intent(in) :: proc
    real(dp), intent(in)        :: amp2(:,:)
    real(dp), intent(out)       :: respdf(ipdf)

    if (ipdf == 1) then
       call get_respdf_one(lumi,as_ord,aem_ord,proc,amp2,respdf)
    elseif (ipdf == 3) then
       call get_respdf_three(lumi,as_ord,aem_ord,proc,amp2,respdf)
    endif

    contains

      subroutine get_respdf_one(lumi,as_ord,aem_ord,proc,amp2,respdf)
        interface
           function lumi(res,f1,f2)
             use mod_types
             real(dp), intent(in) :: res(:,:),f1(-6:),f2(-6:)
             real(dp) :: lumi
           end function lumi
        end interface
        integer, intent(in) :: as_ord,aem_ord
        type(KinConfig), intent(in) :: proc
        real(dp), intent(in)        :: amp2(:,:)
        real(dp), intent(out)       :: respdf(:)
        real(dp) :: f1(-6:7),f2(-6:7),pref
        
        !-- central scale
        call get_prefactor(proc%mur(1),as_ord,aem_ord,pref)
        call get_pdf_qed(proc%PartFrac(1),proc%PartFrac(2),proc%muf(1),f1,f2)
        respdf(1) = pref*lumi(amp2,f1,f2)
    
      end subroutine get_respdf_one

      subroutine get_respdf_three(lumi,as_ord,aem_ord,proc,amp2,respdf)
        interface
           function lumi(res,f1,f2)
             use mod_types
             real(dp), intent(in) :: res(:,:),f1(-6:),f2(-6:)
             real(dp) :: lumi
           end function lumi
        end interface
        integer, intent(in) :: as_ord,aem_ord
        type(KinConfig), intent(in) :: proc
        real(dp), intent(in)        :: amp2(:,:)
        real(dp), intent(out)       :: respdf(:)
        real(dp) :: f1(-6:7),f2(-6:7),pref
        
        !-- central scale
        call get_prefactor(proc%mur(1),as_ord,aem_ord,pref)
        call get_pdf_qed(proc%PartFrac(1),proc%PartFrac(2),proc%muf(1),f1,f2)
        respdf(1) = pref*lumi(amp2,f1,f2)

        !-- do not compute the other scales if we're not doing histograms
        if (nohistos) then

           respdf(2:) = zero

        else
               
           !-- half
           call get_prefactor(proc%mur(1)*half,as_ord,aem_ord,pref)
           call get_pdf_qed(proc%PartFrac(1),proc%PartFrac(2),proc%muf(1)*half,f1,f2)
           respdf(2) = pref*lumi(amp2,f1,f2)
           
           !-- twice
           call get_prefactor(proc%mur(1)*two,as_ord,aem_ord,pref)
           call get_pdf_qed(proc%PartFrac(1),proc%PartFrac(2),proc%muf(1)*two,f1,f2)
           respdf(3) = pref*lumi(amp2,f1,f2)

        endif
    
      end subroutine get_respdf_three

    end subroutine get_respdf

    
      
  !-- evolve PDFs once if we have more limits with same kinematics
  !-- amp2(q/qb,dn/up,i_limit)
  !-- respdf(i_limit,ipdf), note the different order w.r.t. above
  subroutine get_respdf_vect(lumi,as_ord,aem_ord,proc,amp2,respdf)
    interface
       function lumi(res,f1,f2)
         use mod_types
         real(dp), intent(in) :: res(:,:),f1(-6:),f2(-6:)
         real(dp) :: lumi
       end function lumi
    end interface
    integer, intent(in) :: as_ord,aem_ord
    type(KinConfig), intent(in) :: proc
    real(dp), intent(in)        :: amp2(:,:,:)
    real(dp), intent(out)       :: respdf(:,:)
    
    if (ipdf == 1) then
       call get_respdf_vect_one(lumi,as_ord,aem_ord,proc,amp2,respdf)
    elseif (ipdf == 3) then
       call get_respdf_vect_three(lumi,as_ord,aem_ord,proc,amp2,respdf)
    endif

  contains

    subroutine get_respdf_vect_one(lumi,as_ord,aem_ord,proc,amp2,respdf)
      interface
         function lumi(res,f1,f2)
           use mod_types
           real(dp), intent(in) :: res(:,:),f1(-6:),f2(-6:)
           real(dp) :: lumi
         end function lumi
      end interface
      integer, intent(in) :: as_ord,aem_ord
      type(KinConfig), intent(in) :: proc
      real(dp), intent(in)        :: amp2(:,:,:)
      real(dp), intent(out)       :: respdf(:,:)
      real(dp) :: f1(-6:7),f2(-6:7),pref
      integer  :: i 
            
      call get_prefactor(proc%mur(1),as_ord,aem_ord,pref)
      call get_pdf_qed(proc%PartFrac(1),proc%PartFrac(2),proc%muf(1),f1,f2)
      
      do i = 1,size(amp2,3)
         respdf(i,1) = pref*lumi(amp2(:,:,i),f1,f2)
      enddo
      
    end subroutine get_respdf_vect_one

    subroutine get_respdf_vect_three(lumi,as_ord,aem_ord,proc,amp2,respdf)
      interface
         function lumi(res,f1,f2)
           use mod_types
           real(dp), intent(in) :: res(:,:),f1(-6:),f2(-6:)
           real(dp) :: lumi
         end function lumi
      end interface
      integer, intent(in) :: as_ord,aem_ord
      type(KinConfig), intent(in) :: proc
      real(dp), intent(in)        :: amp2(:,:,:)
      real(dp), intent(out)       :: respdf(:,:)
      real(dp) :: f1(-6:7),f2(-6:7),pref
      integer  :: i 

      !-- central scale
      call get_prefactor(proc%mur(1),as_ord,aem_ord,pref)
      call get_pdf_qed(proc%PartFrac(1),proc%PartFrac(2),proc%muf(1),f1,f2)
      
      do i = 1,size(amp2,3)
         respdf(i,1) = pref*lumi(amp2(:,:,i),f1,f2)
      enddo

      !-- do not compute other scales if we're not doing histograms
      if (nohistos) then

         respdf(:,2:) = zero

      else
         
         !-- half scale
         call get_prefactor(proc%mur(1)*half,as_ord,aem_ord,pref)
         call get_pdf_qed(proc%PartFrac(1),proc%PartFrac(2),proc%muf(1)*half,f1,f2)
         
         do i = 1,size(amp2,3)
            respdf(i,2) = pref*lumi(amp2(:,:,i),f1,f2)
         enddo
         
         !-- twice scale
         call get_prefactor(proc%mur(1)*two,as_ord,aem_ord,pref)
         call get_pdf_qed(proc%PartFrac(1),proc%PartFrac(2),proc%muf(1)*two,f1,f2)
         
         do i = 1,size(amp2,3)
            respdf(i,3) = pref*lumi(amp2(:,:,i),f1,f2)
         enddo

      endif
         
    end subroutine get_respdf_vect_three

  end subroutine get_respdf_vect

  !--

    subroutine get_respdf_hoppet_gen(myPDFs1,myPDFs2,as_ord,aem_ord,proc,amp2,respdf,myPDFs1_Lmu,myPDFs2_Lmu)
    use hoppet_v1, except => dp
    interface
       function lumi(res,f1,f2)
         use mod_types
         real(dp), intent(in) :: res(:,:),f1(-6:),f2(-6:)
         real(dp) :: lumi
       end function lumi
    end interface
    type(pdf_table), intent(in)           :: myPDFs1,myPDFs2
    type(pdf_table), intent(in), optional :: myPDFs1_Lmu(:),myPDFs2_Lmu(:)
    integer, intent(in) :: as_ord,aem_ord
    type(KinConfig), intent(in) :: proc
    real(dp), intent(in)        :: amp2(:,:)
    real(dp), intent(out)       :: respdf(:)

    if (ipdf == 1) then
       call get_respdf_hoppet_one_gen(myPDFs1,myPDFs2,as_ord,aem_ord,proc,amp2,respdf,myPDFs1_Lmu,myPDFs2_Lmu)
    elseif (ipdf == 3) then
       call get_respdf_hoppet_three_gen(myPDFs1,myPDFs2,as_ord,aem_ord,proc,amp2,respdf,myPDFs1_Lmu,myPDFs2_Lmu)
    endif

  contains

  subroutine get_respdf_hoppet_one_gen(myPDFs1,myPDFs2,as_ord,aem_ord,proc,amp2,respdf,myPDFs1_Lmu,myPDFs2_Lmu)
      use hoppet_v1, except => dp
      type(pdf_table), intent(in)           :: myPDFs1,myPDFs2
      type(pdf_table), intent(in), optional :: myPDFs1_Lmu(:),myPDFs2_Lmu(:)
      integer, intent(in) :: as_ord,aem_ord
      type(KinConfig), intent(in) :: proc
      real(dp), intent(in)        :: amp2(-5:7,-5:7)
      real(dp), intent(out)       :: respdf(:)
      real(dp) :: pref,mu2ref
      real(dp15) :: f1_dp15(-6:7),f2_dp15(-6:7)

      call get_prefactor(proc%mur(1),as_ord,aem_ord,pref)

      if (present(myPDFs1_Lmu)) then
         mu2ref = proc%mu2ref
         if (present(myPDFs2_Lmu)) then
            call get_pdf_hoppet_qed(myPDFs1,myPDFs2,real(proc%PartFrac(1),kind=dp15),real(proc%PartFrac(2),kind=dp15),&
                 real(proc%muf(1),kind=dp15),f1_dp15,f2_dp15,&
                 myPDFs1_Lmu=myPDFs1_Lmu,myPDFs2_Lmu=myPDFs2_Lmu,mu2ref=real(mu2ref,kind=dp15))
         else
            call get_pdf_hoppet_qed(myPDFs1,myPDFs2,real(proc%PartFrac(1),kind=dp15),real(proc%PartFrac(2),kind=dp15),&
                 real(proc%muf(1),kind=dp15),f1_dp15,f2_dp15,&
                 myPDFs1_Lmu=myPDFs1_Lmu,mu2ref=real(mu2ref,kind=dp15))
         endif
      elseif (present(myPDFs2_Lmu)) then
         mu2ref = proc%mu2ref
         call get_pdf_hoppet_qed(myPDFs1,myPDFs2,real(proc%PartFrac(1),kind=dp15),real(proc%PartFrac(2),kind=dp15),&
              real(proc%muf(1),kind=dp15),f1_dp15,f2_dp15,&
              myPDFs2_Lmu=myPDFs2_Lmu,mu2ref=real(mu2ref,kind=dp15))
      else
         call get_pdf_hoppet_qed(myPDFs1,myPDFs2,real(proc%PartFrac(1),kind=dp15),real(proc%PartFrac(2),kind=dp15),&
              real(proc%muf(1),kind=dp15),f1_dp15,f2_dp15)
      endif

      respdf(1) = pref*gen_lumi(amp2,real(f1_dp15,kind=dp),real(f2_dp15,kind=dp))

    end subroutine get_respdf_hoppet_one_gen

    subroutine get_respdf_hoppet_three_gen(myPDFs1,myPDFs2,as_ord,aem_ord,proc,amp2,respdf,myPDFs1_Lmu,myPDFs2_Lmu)
      use hoppet_v1, except => dp
      type(pdf_table), intent(in)           :: myPDFs1,myPDFs2
      type(pdf_table), intent(in), optional :: myPDFs1_Lmu(:),myPDFs2_Lmu(:)
      integer, intent(in) :: as_ord,aem_ord
      type(KinConfig), intent(in) :: proc
      real(dp), intent(in)        :: amp2(-5:7,-5:7)
      real(dp), intent(out)       :: respdf(:)
      real(dp) :: pref,mu2ref
      real(dp15) :: f1_dp15(-6:7),f2_dp15(-6:7)

      !-- central scale
      call get_prefactor(proc%mur(1),as_ord,aem_ord,pref)

      if (present(myPDFs1_Lmu)) then
         mu2ref = proc%mu2ref
         if (present(myPDFs2_Lmu)) then
            call get_pdf_hoppet_qed(myPDFs1,myPDFs2,real(proc%PartFrac(1),kind=dp15),real(proc%PartFrac(2),kind=dp15),&
                 real(proc%muf(1),kind=dp15),f1_dp15,f2_dp15,&
                 myPDFs1_Lmu=myPDFs1_Lmu,myPDFs2_Lmu=myPDFs2_Lmu,mu2ref=real(mu2ref,kind=dp15))
         else
            call get_pdf_hoppet_qed(myPDFs1,myPDFs2,real(proc%PartFrac(1),kind=dp15),real(proc%PartFrac(2),kind=dp15),&
                 real(proc%muf(1),kind=dp15),f1_dp15,f2_dp15,&
                 myPDFs1_Lmu=myPDFs1_Lmu,mu2ref=real(mu2ref,kind=dp15))
         endif
      elseif (present(myPDFs2_Lmu)) then
         mu2ref = proc%mu2ref
         call get_pdf_hoppet_qed(myPDFs1,myPDFs2,real(proc%PartFrac(1),kind=dp15),real(proc%PartFrac(2),kind=dp15),&
              real(proc%muf(1),kind=dp15),f1_dp15,f2_dp15,&
              myPDFs2_Lmu=myPDFs2_Lmu,mu2ref=real(mu2ref,kind=dp15))
      else
         call get_pdf_hoppet_qed(myPDFs1,myPDFs2,real(proc%PartFrac(1),kind=dp15),real(proc%PartFrac(2),kind=dp15),&
              real(proc%muf(1),kind=dp15),f1_dp15,f2_dp15)
      endif

      respdf(1) = pref*gen_lumi(amp2,real(f1_dp15,kind=dp),real(f2_dp15,kind=dp))

      !-- do not compute the other scales if we're not doing histograms
      if (nohistos) then

         respdf(2:) = zero

      else
      
         !-- half
         call get_prefactor(proc%mur(1)*half,as_ord,aem_ord,pref)
         
         if (present(myPDFs1_Lmu)) then
            mu2ref = proc%mu2ref
            if (present(myPDFs2_Lmu)) then
               call get_pdf_hoppet_qed(myPDFs1,myPDFs2,real(proc%PartFrac(1),kind=dp15),real(proc%PartFrac(2),kind=dp15),&
                    real(proc%muf(1)*half,kind=dp15),f1_dp15,f2_dp15,&
                    myPDFs1_Lmu=myPDFs1_Lmu,myPDFs2_Lmu=myPDFs2_Lmu,mu2ref=real(mu2ref,kind=dp15))
            else
               call get_pdf_hoppet_qed(myPDFs1,myPDFs2,real(proc%PartFrac(1),kind=dp15),real(proc%PartFrac(2),kind=dp15),&
                    real(proc%muf(1)*half,kind=dp15),f1_dp15,f2_dp15,&
                    myPDFs1_Lmu=myPDFs1_Lmu,mu2ref=real(mu2ref,kind=dp15))
            endif
         elseif (present(myPDFs2_Lmu)) then
            mu2ref = proc%mu2ref
            call get_pdf_hoppet_qed(myPDFs1,myPDFs2,real(proc%PartFrac(1),kind=dp15),real(proc%PartFrac(2),kind=dp15),&
                 real(proc%muf(1)*half,kind=dp15),f1_dp15,f2_dp15,&
                 myPDFs2_Lmu=myPDFs2_Lmu,mu2ref=real(mu2ref,kind=dp15))
         else
            call get_pdf_hoppet_qed(myPDFs1,myPDFs2,real(proc%PartFrac(1),kind=dp15),real(proc%PartFrac(2),kind=dp15),&
                 real(proc%muf(1)*half,kind=dp15),f1_dp15,f2_dp15)
         endif
         
         respdf(2) = pref*gen_lumi(amp2,real(f1_dp15,kind=dp),real(f2_dp15,kind=dp))
         
         !-- two
         call get_prefactor(proc%mur(1)*two,as_ord,aem_ord,pref)
         
         if (present(myPDFs1_Lmu)) then
            mu2ref = proc%mu2ref
            if (present(myPDFs2_Lmu)) then
               call get_pdf_hoppet_qed(myPDFs1,myPDFs2,real(proc%PartFrac(1),kind=dp15),real(proc%PartFrac(2),kind=dp15),&
                    real(proc%muf(1)*two,kind=dp15),f1_dp15,f2_dp15,&
                    myPDFs1_Lmu=myPDFs1_Lmu,myPDFs2_Lmu=myPDFs2_Lmu,mu2ref=real(mu2ref,kind=dp15))
            else
               call get_pdf_hoppet_qed(myPDFs1,myPDFs2,real(proc%PartFrac(1),kind=dp15),real(proc%PartFrac(2),kind=dp15),&
                    real(proc%muf(1)*two,kind=dp15),f1_dp15,f2_dp15,&
                    myPDFs1_Lmu=myPDFs1_Lmu,mu2ref=real(mu2ref,kind=dp15))
            endif
         elseif (present(myPDFs2_Lmu)) then
            mu2ref = proc%mu2ref
            call get_pdf_hoppet_qed(myPDFs1,myPDFs2,real(proc%PartFrac(1),kind=dp15),real(proc%PartFrac(2),kind=dp15),&
                 real(proc%muf(1)*two,kind=dp15),f1_dp15,f2_dp15,&
                 myPDFs2_Lmu=myPDFs2_Lmu,mu2ref=real(mu2ref,kind=dp15))
         else
            call get_pdf_hoppet_qed(myPDFs1,myPDFs2,real(proc%PartFrac(1),kind=dp15),real(proc%PartFrac(2),kind=dp15),&
                 real(proc%muf(1)*two,kind=dp15),f1_dp15,f2_dp15)
         endif
         
         respdf(3) = pref*gen_lumi(amp2,real(f1_dp15,kind=dp),real(f2_dp15,kind=dp))

      endif
         
    end subroutine get_respdf_hoppet_three_gen
 
  end subroutine get_respdf_hoppet_gen


  subroutine get_respdf_hoppet(myPDFs1,myPDFs2,lumi,as_ord,aem_ord,proc,amp2,respdf,myPDFs1_Lmu,myPDFs2_Lmu)
    use hoppet_v1, except => dp
    interface
       function lumi(res,f1,f2)
         use mod_types
         real(dp), intent(in) :: res(:,:),f1(-6:),f2(-6:)
         real(dp) :: lumi
       end function lumi
    end interface
    type(pdf_table), intent(in)           :: myPDFs1,myPDFs2
    type(pdf_table), intent(in), optional :: myPDFs1_Lmu(:),myPDFs2_Lmu(:)
    integer, intent(in) :: as_ord,aem_ord
    type(KinConfig), intent(in) :: proc
    real(dp), intent(in)        :: amp2(:,:)
    real(dp), intent(out)       :: respdf(:)
    
    if (ipdf == 1) then
       call get_respdf_hoppet_one(myPDFs1,myPDFs2,lumi,as_ord,aem_ord,proc,amp2,respdf,myPDFs1_Lmu,myPDFs2_Lmu)
    elseif (ipdf == 3) then
       call get_respdf_hoppet_three(myPDFs1,myPDFs2,lumi,as_ord,aem_ord,proc,amp2,respdf,myPDFs1_Lmu,myPDFs2_Lmu)
    endif

  contains
  
    subroutine get_respdf_hoppet_one(myPDFs1,myPDFs2,lumi,as_ord,aem_ord,proc,amp2,respdf,myPDFs1_Lmu,myPDFs2_Lmu)
      use hoppet_v1, except => dp
      interface
         function lumi(res,f1,f2)
           use mod_types
           real(dp), intent(in) :: res(:,:),f1(-6:),f2(-6:)
           real(dp) :: lumi
         end function lumi
      end interface
      type(pdf_table), intent(in)           :: myPDFs1,myPDFs2
      type(pdf_table), intent(in), optional :: myPDFs1_Lmu(:),myPDFs2_Lmu(:)
      integer, intent(in) :: as_ord,aem_ord
      type(KinConfig), intent(in) :: proc
      real(dp), intent(in)        :: amp2(:,:)
      real(dp), intent(out)       :: respdf(:)
      real(dp) :: pref,mu2ref
      real(dp15) :: f1_dp15(-6:7),f2_dp15(-6:7)

      call get_prefactor(proc%mur(1),as_ord,aem_ord,pref)

      if (present(myPDFs1_Lmu)) then
         mu2ref = proc%mu2ref
         if (present(myPDFs2_Lmu)) then
            call get_pdf_hoppet_qed(myPDFs1,myPDFs2,real(proc%PartFrac(1),kind=dp15),real(proc%PartFrac(2),kind=dp15),&
                 real(proc%muf(1),kind=dp15),f1_dp15,f2_dp15,&
                 myPDFs1_Lmu=myPDFs1_Lmu,myPDFs2_Lmu=myPDFs2_Lmu,mu2ref=real(mu2ref,kind=dp15))
         else
            call get_pdf_hoppet_qed(myPDFs1,myPDFs2,real(proc%PartFrac(1),kind=dp15),real(proc%PartFrac(2),kind=dp15),&
                 real(proc%muf(1),kind=dp15),f1_dp15,f2_dp15,&
                 myPDFs1_Lmu=myPDFs1_Lmu,mu2ref=real(mu2ref,kind=dp15))
         endif
      elseif (present(myPDFs2_Lmu)) then
         mu2ref = proc%mu2ref
         call get_pdf_hoppet_qed(myPDFs1,myPDFs2,real(proc%PartFrac(1),kind=dp15),real(proc%PartFrac(2),kind=dp15),&
              real(proc%muf(1),kind=dp15),f1_dp15,f2_dp15,&
              myPDFs2_Lmu=myPDFs2_Lmu,mu2ref=real(mu2ref,kind=dp15))
      else
         call get_pdf_hoppet_qed(myPDFs1,myPDFs2,real(proc%PartFrac(1),kind=dp15),real(proc%PartFrac(2),kind=dp15),&
              real(proc%muf(1),kind=dp15),f1_dp15,f2_dp15)
      endif

      respdf(1) = pref*lumi(amp2,real(f1_dp15,kind=dp),real(f2_dp15,kind=dp))

    end subroutine get_respdf_hoppet_one

    subroutine get_respdf_hoppet_three(myPDFs1,myPDFs2,lumi,as_ord,aem_ord,proc,amp2,respdf,myPDFs1_Lmu,myPDFs2_Lmu)
      use hoppet_v1, except => dp
      interface
         function lumi(res,f1,f2)
           use mod_types
           real(dp), intent(in) :: res(:,:),f1(-6:),f2(-6:)
           real(dp) :: lumi
         end function lumi
      end interface
      type(pdf_table), intent(in)           :: myPDFs1,myPDFs2
      type(pdf_table), intent(in), optional :: myPDFs1_Lmu(:),myPDFs2_Lmu(:)
      integer, intent(in) :: as_ord,aem_ord
      type(KinConfig), intent(in) :: proc
      real(dp), intent(in)        :: amp2(:,:)
      real(dp), intent(out)       :: respdf(:)
      real(dp) :: pref,mu2ref
      real(dp15) :: f1_dp15(-6:7),f2_dp15(-6:7)

      !-- central scale
      call get_prefactor(proc%mur(1),as_ord,aem_ord,pref)

      if (present(myPDFs1_Lmu)) then
         mu2ref = proc%mu2ref
         if (present(myPDFs2_Lmu)) then
            call get_pdf_hoppet_qed(myPDFs1,myPDFs2,real(proc%PartFrac(1),kind=dp15),real(proc%PartFrac(2),kind=dp15),&
                 real(proc%muf(1),kind=dp15),f1_dp15,f2_dp15,&
                 myPDFs1_Lmu=myPDFs1_Lmu,myPDFs2_Lmu=myPDFs2_Lmu,mu2ref=real(mu2ref,kind=dp15))
         else
            call get_pdf_hoppet_qed(myPDFs1,myPDFs2,real(proc%PartFrac(1),kind=dp15),real(proc%PartFrac(2),kind=dp15),&
                 real(proc%muf(1),kind=dp15),f1_dp15,f2_dp15,&
                 myPDFs1_Lmu=myPDFs1_Lmu,mu2ref=real(mu2ref,kind=dp15))
         endif
      elseif (present(myPDFs2_Lmu)) then
         mu2ref = proc%mu2ref
         call get_pdf_hoppet_qed(myPDFs1,myPDFs2,real(proc%PartFrac(1),kind=dp15),real(proc%PartFrac(2),kind=dp15),&
              real(proc%muf(1),kind=dp15),f1_dp15,f2_dp15,&
              myPDFs2_Lmu=myPDFs2_Lmu,mu2ref=real(mu2ref,kind=dp15))
      else
         call get_pdf_hoppet_qed(myPDFs1,myPDFs2,real(proc%PartFrac(1),kind=dp15),real(proc%PartFrac(2),kind=dp15),&
              real(proc%muf(1),kind=dp15),f1_dp15,f2_dp15)
      endif

      respdf(1) = pref*lumi(amp2,real(f1_dp15,kind=dp),real(f2_dp15,kind=dp))

      !-- do not compute the other scales if we're not doing histograms
      if (nohistos) then

         respdf(2:) = zero

      else
      
         !-- half
         call get_prefactor(proc%mur(1)*half,as_ord,aem_ord,pref)
         
         if (present(myPDFs1_Lmu)) then
            mu2ref = proc%mu2ref
            if (present(myPDFs2_Lmu)) then
               call get_pdf_hoppet_qed(myPDFs1,myPDFs2,real(proc%PartFrac(1),kind=dp15),real(proc%PartFrac(2),kind=dp15),&
                    real(proc%muf(1)*half,kind=dp15),f1_dp15,f2_dp15,&
                    myPDFs1_Lmu=myPDFs1_Lmu,myPDFs2_Lmu=myPDFs2_Lmu,mu2ref=real(mu2ref,kind=dp15))
            else
               call get_pdf_hoppet_qed(myPDFs1,myPDFs2,real(proc%PartFrac(1),kind=dp15),real(proc%PartFrac(2),kind=dp15),&
                    real(proc%muf(1)*half,kind=dp15),f1_dp15,f2_dp15,&
                    myPDFs1_Lmu=myPDFs1_Lmu,mu2ref=real(mu2ref,kind=dp15))
            endif
         elseif (present(myPDFs2_Lmu)) then
            mu2ref = proc%mu2ref
            call get_pdf_hoppet_qed(myPDFs1,myPDFs2,real(proc%PartFrac(1),kind=dp15),real(proc%PartFrac(2),kind=dp15),&
                 real(proc%muf(1)*half,kind=dp15),f1_dp15,f2_dp15,&
                 myPDFs2_Lmu=myPDFs2_Lmu,mu2ref=real(mu2ref,kind=dp15))
         else
            call get_pdf_hoppet_qed(myPDFs1,myPDFs2,real(proc%PartFrac(1),kind=dp15),real(proc%PartFrac(2),kind=dp15),&
                 real(proc%muf(1)*half,kind=dp15),f1_dp15,f2_dp15)
         endif
         
         respdf(2) = pref*lumi(amp2,real(f1_dp15,kind=dp),real(f2_dp15,kind=dp))
         
         !-- two
         call get_prefactor(proc%mur(1)*two,as_ord,aem_ord,pref)
         
         if (present(myPDFs1_Lmu)) then
            mu2ref = proc%mu2ref
            if (present(myPDFs2_Lmu)) then
               call get_pdf_hoppet_qed(myPDFs1,myPDFs2,real(proc%PartFrac(1),kind=dp15),real(proc%PartFrac(2),kind=dp15),&
                    real(proc%muf(1)*two,kind=dp15),f1_dp15,f2_dp15,&
                    myPDFs1_Lmu=myPDFs1_Lmu,myPDFs2_Lmu=myPDFs2_Lmu,mu2ref=real(mu2ref,kind=dp15))
            else
               call get_pdf_hoppet_qed(myPDFs1,myPDFs2,real(proc%PartFrac(1),kind=dp15),real(proc%PartFrac(2),kind=dp15),&
                    real(proc%muf(1)*two,kind=dp15),f1_dp15,f2_dp15,&
                    myPDFs1_Lmu=myPDFs1_Lmu,mu2ref=real(mu2ref,kind=dp15))
            endif
         elseif (present(myPDFs2_Lmu)) then
            mu2ref = proc%mu2ref
            call get_pdf_hoppet_qed(myPDFs1,myPDFs2,real(proc%PartFrac(1),kind=dp15),real(proc%PartFrac(2),kind=dp15),&
                 real(proc%muf(1)*two,kind=dp15),f1_dp15,f2_dp15,&
                 myPDFs2_Lmu=myPDFs2_Lmu,mu2ref=real(mu2ref,kind=dp15))
         else
            call get_pdf_hoppet_qed(myPDFs1,myPDFs2,real(proc%PartFrac(1),kind=dp15),real(proc%PartFrac(2),kind=dp15),&
                 real(proc%muf(1)*two,kind=dp15),f1_dp15,f2_dp15)
         endif
         
         respdf(3) = pref*lumi(amp2,real(f1_dp15,kind=dp),real(f2_dp15,kind=dp))

      endif
         
    end subroutine get_respdf_hoppet_three
    
  end subroutine get_respdf_hoppet

  !-- scale-variation logs

  subroutine fill_sv_logs(muf2,muref2,svlogs)
    implicit none
    real(dp), intent(in)  :: muf2,muref2
    real(dp), intent(out) :: svlogs(ipdf)

    if (ipdf == 1) then
      svlogs = log(muf2/muref2)
    elseif(ipdf == 3) then
      svlogs(1) = log(muf2/muref2)
      svlogs(2) = log(muf2/muref2/4)
      svlogs(3) = log(muf2/muref2*4)
    endif

  end subroutine fill_sv_logs

  !-- luminosities


  function gen_lumi(res,f1,f2) result(respdf)
    ! added by Raoul
    real(dp), intent(in) :: res(-5:7,-5:7),f1(-6:),f2(-6:)
    real(dp) :: respdf
    integer  :: i1,i2

    respdf = zero
    do i1= -5,7
       do i2 = -5,7
          respdf = respdf + f1(i1)*f2(i2)*res(i1,i2)
       enddo
    enddo

  end function gen_lumi

    
  
  function ns_lumi(res,f1,f2) result(respdf)
    real(dp), intent(in) :: res(:,:),f1(-6:),f2(-6:)
    real(dp) :: respdf

    respdf = zero

    ! u(1) ux(2) + c(1) cx(2)
    respdf = respdf + res(1,2)*(f1(2)*f2(-2) + f1(4)*f2(-4))

    ! d(1) dx(2) + s(1) sx(2) + b(1) bx(2)
    respdf = respdf + flag_down*res(1,1)*(f1(1)*f2(-1) + f1(3)*f2(-3) + f1(5)*f2(-5))

    ! ux(1) u(2) + cx(1) + c(2)
    respdf = respdf + res(2,2)*(f1(-2)*f2(2) + f1(-4)*f2(4))

    ! dx(1) d(2) + sx(1) s(2) + bx(1) b(2)
    respdf = respdf + flag_down*res(2,1)*(f1(-1)*f2(1) + f1(-3)*f2(3) + f1(-5)*f2(5))

  end function ns_lumi
  
  function ns_lumi_splitb(res,f1,f2) result(respdf)
    real(dp), intent(in) :: res(:,:),f1(-6:),f2(-6:)
    real(dp) :: respdf

    respdf = zero

    ! u(1) ux(2) + c(1) cx(2)
    respdf = respdf + res(1,2)*(f1(2)*f2(-2) + f1(4)*f2(-4))

    ! d(1) dx(2) + s(1) sx(2) 
    respdf = respdf + flag_down*res(1,1)*(f1(1)*f2(-1) + f1(3)*f2(-3))

    ! b(1) bx(2)
    respdf = respdf + flag_down*res(1,3)*(f1(5)*f2(-5)) 
    
    ! ux(1) u(2) + cx(1) + c(2)
    respdf = respdf + res(2,2)*(f1(-2)*f2(2) + f1(-4)*f2(4))

    ! dx(1) d(2) + sx(1) s(2) 
    respdf = respdf + flag_down*res(2,1)*(f1(-1)*f2(1) + f1(-3)*f2(3))

    ! bx(1) b(2)
    respdf = respdf + flag_down*res(2,3)*(f1(-5)*f2(5))
    
  end function ns_lumi_splitb

  function qq_lumi(res,f1,f2) result(respdf)
    real(dp), intent(in) :: res(:,:),f1(-6:),f2(-6:)
    real(dp) :: respdf

    respdf = zero

    ! u(1) u(2) + c(1) c(2)
    respdf = respdf + res(1,2)*(f1(2)*f2(2) + f1(4)*f2(4))

    ! d(1) d(2) + s(1) s(2) + b(1) b(2)
    respdf = respdf + flag_down*res(1,1)*(f1(1)*f2(1) + f1(3)*f2(3) + f1(5)*f2(5)  ) 

    ! ux(1) ux(2) + cx(1) + cx(2)
    respdf = respdf + res(2,2)*(f1(-2)*f2(-2) + f1(-4)*f2(-4))

    ! dx(1) dx(2) + sx(1) sx(2) + bx(1) bx(2)
    respdf = respdf + flag_down*res(2,1)*(f1(-1)*f2(-1) + f1(-3)*f2(-3) + f1(-5)*f2(-5)  )

  end function qq_lumi

  function qQp_lumi(res,f1,f2) result(respdf)
    real(dp), intent(in) :: res(:,:),f1(-6:),f2(-6:)
    real(dp) :: respdf

    respdf = zero

    ! d(1) u(2) + s(1) c(2)
    respdf = respdf + flag_down*res(1,1)*(f1(1)*f2(2) + f1(3)*f2(4))

    ! dx(1) ux(2) + sx(1) cx(2)
    respdf = respdf + flag_down*res(2,1)*(f1(-1)*f2(-2) + f1(-3)*f2(-4))

    ! u(1) d(2) + c(1) s(2)
    respdf = respdf + flag_down*res(1,2)*(f1(2)*f2(1) + f1(4)*f2(3))

    ! ux(1) dx(2) + cx(1) sx(2)
    respdf = respdf + flag_down*res(2,2)*(f1(-2)*f2(-1) + f1(-4)*f2(-3))

  end function qQp_lumi

  function qQpb_lumi(res,f1,f2) result(respdf)
    real(dp), intent(in) :: res(:,:),f1(-6:),f2(-6:)
    real(dp) :: respdf

    respdf = zero

    ! d(1) ux(2) + s(1) cx(2)
    respdf = respdf + flag_down*res(1,1)*(f1(1)*f2(-2) + f1(3)*f2(-4))

    ! dx(1) u(2) + sx(1) c(2)
    respdf = respdf + flag_down*res(2,1)*(f1(-1)*f2(2) + f1(-3)*f2(4))

    ! u(1) dx(2) + c(1) sx(2)
    respdf = respdf + flag_down*res(1,2)*(f1(2)*f2(-1) + f1(4)*f2(-3))

    ! ux(1) d(2) + cx(1) s(2)
    respdf = respdf + flag_down*res(2,2)*(f1(-2)*f2(1) + f1(-4)*f2(3))

  end function qQpb_lumi

  function qQpb_lumi_wp(res,f1,f2) result(respdf)
    real(dp), intent(in) :: res(:,:),f1(-6:),f2(-6:)
    real(dp) :: respdf

    respdf = zero

    ! d(1) ux(2) + s(1) cx(2)
    ! respdf = respdf + flag_down*res(2,1)*(f1(1)*f2(-2) + f1(3)*f2(-4))

    ! dx(1) u(2) + sx(1) c(2)
    respdf = respdf + flag_down*res(1,2)*(f1(-1)*f2(2) + f1(-3)*f2(4))

    ! u(1) dx(2) + c(1) sx(2)
    respdf = respdf + flag_down*res(1,1)*(f1(2)*f2(-1) + f1(4)*f2(-3))

    ! ux(1) d(2) + cx(1) s(2)
    ! respdf = respdf + flag_down*res(1,2)*(f1(-2)*f2(1) + f1(-4)*f2(3))

  end function qQpb_lumi_wp


  function qQpb_lumi_wm(res,f1,f2) result(respdf)
    real(dp), intent(in) :: res(:,:),f1(-6:),f2(-6:)
    real(dp) :: respdf

    respdf = zero

    ! d(1) ux(2) + s(1) cx(2)
    respdf = respdf + flag_down*res(1,1)*(f1(1)*f2(-2) + f1(3)*f2(-4))

    ! dx(1) u(2) + sx(1) c(2)
    ! respdf = respdf + flag_down*res(2,2)*(f1(-1)*f2(2) + f1(-3)*f2(4))

    ! u(1) dx(2) + c(1) sx(2)
    ! respdf = respdf + flag_down*res(1,1)*(f1(2)*f2(-1) + f1(4)*f2(-3))

    ! ux(1) d(2) + cx(1) s(2)
    respdf = respdf + flag_down*res(1,2)*(f1(-2)*f2(1) + f1(-4)*f2(3))

  end function qQpb_lumi_wm
  
  function aa_lumi(res,f1,f2) result(respdf)
    real(dp), intent(in) :: res(:,:),f1(-6:),f2(-6:)
    real(dp) :: respdf

    respdf = res(1,1) * f1(7) * f2(7)

  end function aa_lumi
  
  function gq_lumi_splitb(res,f1,f2) result(respdf)
    real(dp), intent(in) :: res(:,:),f1(-6:),f2(-6:)
    real(dp) :: respdf

    respdf = zero

    ! g(1) ux(2) + g(1) cx(2)
    respdf = respdf + res(1,2)*(f2(-2) + f2(-4))

    ! g(1) dx(2) + g(1) sx(2)
    respdf = respdf + flag_down*res(1,1)*(f2(-1) + f2(-3)) 

    ! g(1) bx(2)
    respdf = respdf + flag_down*res(1,3)*(f2(-5)) 
    
    ! g(1) u(2) + g(1) + c(2)
    respdf = respdf + res(2,2)*(f2(2) + f2(4))

    ! g(1) d(2) + g(1) s(2) 
    respdf = respdf + flag_down*res(2,1)*(f2(1) + f2(3))

    ! g(1) b(2)
    respdf = respdf + flag_down*res(2,3)*(f2(5))

    respdf = respdf * f1(0)
    
  end function gq_lumi_splitb

  function qg_lumi_splitb(res,f1,f2) result(respdf)
    real(dp), intent(in) :: res(:,:),f1(-6:),f2(-6:)
    real(dp) :: respdf

    respdf = zero

    ! g(1) ux(2) + g(1) cx(2)
    respdf = respdf + res(1,2)*(f1(-2) + f1(-4))

    ! g(1) dx(2) + g(1) sx(2) 
    respdf = respdf + flag_down*res(1,1)*(f1(-1) + f1(-3)) 

    ! g(1) bx(2)
    respdf = respdf + flag_down*res(1,3)*(f1(-5)) 
    
    ! g(1) u(2) + g(1) + c(2)
    respdf = respdf + res(2,2)*(f1(2) + f1(4))

    ! g(1) d(2) + g(1) s(2) 
    respdf = respdf + flag_down*res(2,1)*(f1(1) + f1(3))

    ! g(1) b(2)
    respdf = respdf + flag_down*res(2,3)*(f1(5))

    respdf = respdf * f2(0)
    
  end function qg_lumi_splitb

  function aq_lumi(res,f1,f2) result(respdf)
    real(dp), intent(in) :: res(:,:),f1(-6:),f2(-6:)
    real(dp) :: respdf

    respdf = zero

    ! a(1) ux(2) + a(1) cx(2)
    respdf = respdf + res(1,2)*(f2(-2) + f2(-4))

    ! a(1) dx(2) + a(1) sx(2) + a(1) bx(2)
    respdf = respdf + flag_down*res(1,1)*(f2(-1) + f2(-3) + f2(-5)) 
    
    ! a(1) u(2) + a(1) + c(2)
    respdf = respdf + res(2,2)*(f2(2) + f2(4))

    ! a(1) d(2) + a(1) s(2) + a(1) b(2)
    respdf = respdf + flag_down*res(2,1)*(f2(1) + f2(3) + f2(5))

    respdf = respdf * f1(7)
    
  end function aq_lumi

  function qa_lumi(res,f1,f2) result(respdf)
    real(dp), intent(in) :: res(:,:),f1(-6:),f2(-6:)
    real(dp) :: respdf
    respdf = zero
    ! u(1) a(2) + c(1) a(2)
    respdf = respdf + res(1,2)*(f1(2) + f1( 4))
    ! d(1) a(2) + s(1) a(2) + b(1) a(2)
    respdf = respdf + flag_down*res(1,1)*(f1(1) + f1(3) + f1(5))
    ! ux(1) a(2) + cx(1) a(2)
    respdf = respdf + res(2,2)*(f1(-2) + f1(-4))
    ! dx(1) a(2) + sx(1) a(2) + bx(1) a(2)
    respdf = respdf + flag_down*res(2,1)*(f1(-1) + f1(-3) + f1(-5))

    respdf = respdf * f2(7)
  end function qa_lumi

  function gq_lumi(res,f1,f2) result(respdf)
    real(dp), intent(in) :: res(:,:),f1(-6:),f2(-6:)
    real(dp) :: respdf
    respdf = zero
    ! g(1) ux(2) + g(1) cx(2)
    respdf = respdf + res(1,2)*(f1(0)*f2(-2) + f1(0)*f2(-4))
    ! g(1) dx(2) + g(1) sx(2) + g(1) bx(2)
    respdf = respdf + flag_down*res(1,1)*(f1(0)*f2(-1) + f1(0)*f2(-3) + f1(0)*f2(-5) )
    ! g(1) u(2) + g(1) c(2)
    respdf = respdf + res(2,2)*(f1(0)*f2(2) + f1(0)*f2(4))
    ! g(1) d(2) + g(1) s(2) + g(1) b(2)
    respdf = respdf + flag_down*res(2,1)*(f1(0)*f2(1) + f1(0)*f2(3) + f1(0)*f2(5) )
  end function gq_lumi

  function qg_lumi(res,f1,f2) result(respdf)
    real(dp), intent(in) :: res(:,:),f1(-6:),f2(-6:)
    real(dp) :: respdf
    respdf = zero
    ! u(1) g(2) + c(1) g(2)
    respdf = respdf + res(1,2)*(f2(0)*f1( 2) + f2(0)*f1( 4))
    ! d(1) g(2) + s(1) g(2) + b(1) g(2)
    respdf = respdf + flag_down*res(1,1)*(f2(0)*f1( 1) + f2(0)*f1( 3) + f2(0)*f1( 5))
    ! ux(1) g(2) + cx(1) g(2)
    respdf = respdf + res(2,2)*(f2(0)*f1(-2) + f2(0)*f1(-4))
    ! dx(1) g(2) + sx(1) g(2) + bx(1) g(2)
    respdf = respdf + flag_down*res(2,1)*(f2(0)*f1(-1) + f2(0)*f1(-3) + f2(0)*f1(-5))
  end function qg_lumi

  function ga_lumi(res,f1,f2) result(respdf)
    real(dp), intent(in) :: res(:,:),f1(-6:),f2(-6:)
    real(dp) :: respdf
    respdf = res(1,1)*f1(0)*f2(7)
  end function ga_lumi

  function ag_lumi(res,f1,f2) result(respdf)
    real(dp), intent(in) :: res(:,:),f1(-6:),f2(-6:)
    real(dp) :: respdf
    respdf = res(1,1)*f1(7)*f2(0)
  end function ag_lumi




  function multiply_IS_charges(res_in,leg) result(res_out)
    real(dp), intent(in)  :: res_in(-5:7,-5:7)
    integer, intent(in)   :: leg
    real(dp) :: res_out(-5:7,-5:7)
    integer               :: i,j

    

    do i = -5,7
       do j = -5,7
          if (leg .eq. 1) then
             res_out(i,j) = res_in(i,j)*Qsq_IS(i)
          elseif (leg .eq. 2) then
             res_out(i,j) = res_in(i,j)*Qsq_IS(j)
          endif
       enddo
    enddo

  end function multiply_IS_charges
  !--

  !-- NaN checks
  subroutine check_ff(ff,xx,myFint)
    real(dp15) :: ff(1)
    real(dp), intent(in) :: xx(:),myFint(:)
    integer :: i

    if (is_nan(real(ff(1),dp))) then
       print *, 'NaN, setting to zero'

       do i = 1,size(xx)
          print *, 'xx',i,xx(i)
       enddo

       do i = 1,size(myFint)
          print *, 'Fint',i,myFint(i)
       enddo

       ff(1) = 0

       icount_nan = icount_nan + 1

    else

       icount_good = icount_good + 1
       if (ff(1).ne.zero) icount_acc = icount_acc + 1
           
    end if

  end subroutine check_ff


  ! >--------------------------------------------------------------------------
  ! | module procedures of function transision
  ! | - transition_res for argument res(-6:6,-6:6)
  ! | - transition_res_sping for argument res_sping(-1:1,-1:1,-6:6,-6:6)
  ! >--------------------------------------------------------------------------
  !-- private helper function that generates the transition matrices
  subroutine get_transisition_matrix(trans, trans_matrix)
    implicit none
    
    real(dp), intent(out) :: trans_matrix(-5:7,-5:7)
    character(len=*), intent(in) :: trans
    
    integer :: i1, i2
    
    trans_matrix = zero
    
    select case (trans)
    case ('none')
       do i1 = -5, 7, 1
          trans_matrix(i1,i1) = one
       end do
    case ('g -> q')
       do i1 = -5, 7, 1
          if (i1 .ne. 0 .and. i1 .ne. 7) trans_matrix(0,i1) = one
       end do
    case ('ga -> q')
       do i1 = -5, 7, 1
          if (i1 .ne. 0 .and. i1 .ne. 7) trans_matrix(7,i1) = one
       end do
    case ('q -> q')
       do i1 = -int(Nf), int(Nf), 1
          do i2 = -int(Nf), int(Nf), 1
             if ((i1 .ne. 0) .and. (i2 .ne. 0) .and. (i1 .ne. 7) .and. (i2 .ne. 7)) trans_matrix(i1,i2) = one
          end do
       end do
    case ('q -> qb')
       do i1 = -6, 6, 1
          if (i1 .ne. 0 .and. i1 .ne. 7) trans_matrix(i1,-i1) = one
       end do
    case ('q -> g')
       do i1 = -int(Nf), int(Nf), 1
          if (i1 .ne. 0 .and. i1 .ne. 7) trans_matrix(i1,0) = one
       end do
    case ('q -> ga')
       do i1 = -int(Nf), int(Nf), 1
          if (i1 .ne. 0 .and. i1 .ne. 7) trans_matrix(i1,7) = one
       end do
    case default
       print *, 'error: transition function called with unknown splitting.'
       stop
    end select
    
    return
    

  end subroutine get_transisition_matrix
  
  function transition_res(beam_1, beam_2, res)
    implicit none
    
    real(dp) :: transition_res(-6:6,-6:6)
    character(len=*), intent(in) :: beam_1
    character(len=*), intent(in) :: beam_2
    real(dp), intent(in) :: res(-6:6,-6:6)
    real(dp) :: trans_matrix_b1(-6:6,-6:6), trans_matrix_b2(-6:6,-6:6)
    
    call get_transisition_matrix(beam_1, trans_matrix_b1)
    call get_transisition_matrix(beam_2, trans_matrix_b2)
    
    transition_res = matmul(matmul(trans_matrix_b1, res), transpose(trans_matrix_b2))
    
    return
  end function transition_res
  
  function transition_res_sping(beam_1, beam_2, res_sping)
    implicit none
    complex(dp) :: transition_res_sping(-1:1,-1:1,-6:6,-6:6)
    character(len=*), intent(in) :: beam_1
    character(len=*), intent(in) :: beam_2
    complex(dp), intent(in) :: res_sping(-1:1,-1:1,-6:6,-6:6)
    real(dp) :: trans_matrix_b1(-6:6,-6:6), trans_matrix_b2(-6:6,-6:6)
    integer :: h1, h2
    
    call get_transisition_matrix(beam_1, trans_matrix_b1)
    call get_transisition_matrix(beam_1, trans_matrix_b2)
    
    do h1 = -1, 1, 2
       do h2 = -1, 1, 2
          transition_res_sping(h1,h2,:,:) = matmul(matmul(trans_matrix_b1, res_sping(h1,h2,:,:)), transpose(trans_matrix_b2))
       end do
    end do
    
    return
  end function transition_res_sping
  

 
end module mod_aux_sectors
