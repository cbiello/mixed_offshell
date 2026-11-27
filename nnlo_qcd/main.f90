program color_singlet
  use mod_chooser
  !-- debug
  !use mod_subtrfn_nnlo_h_qg
  !real (dp) :: shat,mu2,z,eta41,eta42,eta41log,e4
  !z=0.123_dp
  !shat = 11.23_dp
  !mu2 = 24.0_dp
  !eta41 = 0.456_dp
  !eta42 = one-eta41
  !eta41log = eta41
  !e4 = 0.321_dp

  !print *, qg_qg_sub_qgz124(shat,mu2,z,eta41,eta42,eta41log)
  !print *, qg_qg_sub_qgz124_pls(shat,mu2,z,eta41,eta42,eta41log)
  !print *, qg_qg_sub_qg1z24(shat,mu2,z,eta41,eta42)
  !print *, qg_qg_sub_qg1z24_pls(shat,mu2,z,eta41,eta42)
  !print *, qg_qg_sub_qg124(shat,mu2,E4,eta41,eta42,eta41log)
  !print *, qg_qg_sub_qqb1z24(shat,mu2,z,eta41,eta42)
  !print *, qg_qg_sub_ggz124(shat,mu2,z,eta41,eta42,eta41log)
  !stop
  
  print *, '******************************'
  print *, '*    Color singlet @ NNLO    *'
  print *, '******************************'

  call chooser()
end  program color_singlet

subroutine write_from_lepage(itnumb)
  use mod_types; use common_def; use mod_consts_dp
  use mod_auxfunctions
  use mod_histograms
  implicit none
  integer :: itnumb,i,j,myiternumber
  real(dp) :: HistoGlobalArray_tmp(nobs*numofblocks,nbinmax),ErrorGlobalArray_tmp(nobs*numofblocks,nbinmax)
  character :: tmpfile*(100),ch_iternumber*(3)
  real(dp) :: HistoGlobalArray_tmp_markus(nobs*numofblocks,nbinmax),mytestM
  character :: outputfileM*(100),tmpfileM*(100)
  real(dp) :: mytest,mytest_err
  logical :: flag_nan

  myiternumber = itnumb

  include './Histograms/prepare_and_write.f90'
  include './Histograms/check_sumbin.f90'

  HistoGlobalArray = HistoGlobalArray_tmp
  ErrorGlobalArray = ErrorGlobalArray_tmp
  evalcount = 1

  return

end subroutine write_from_lepage
