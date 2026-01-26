subroutine G_test (inp) bind(C)
!  use handyG
  print *,inp
end subroutine G_test

subroutine G_weight1_fortran (re1,im1,reX,imX,re_ans,im_ans) bind(C)
  use handyG
  real(kind=prec) :: re1,im1,reX,imX,re_ans,im_ans
  complex(kind=prec) :: a1, x, ans
  call clearcache
  a1 = complex(re1, im1)
  x  = complex(reX, imX)
  ans = G([a1],x)
  re_ans = real(ans)
  im_ans = aimag(ans)
end subroutine G_weight1_fortran

subroutine G_weight2_fortran (re1,im1,re2,im2,reX,imX,re_ans,im_ans) bind(C)
  use handyG
  real(kind=prec) :: re1,im1,re2,im2,reX,imX,re_ans,im_ans
  complex(kind=prec) :: a1, a2, x, ans
  call clearcache
  a1 = complex(re1, im1)
  a2 = complex(re2, im2)
  x  = complex(reX, imX)
  ans = G([a1,a2],x)
  re_ans = real(ans)
  im_ans = aimag(ans)
end subroutine G_weight2_fortran

subroutine G_weight3_fortran (re1,im1,re2,im2,re3,im3,reX,imX,re_ans,im_ans) bind(C)
  use handyG
  real(kind=prec) :: re1,im1,re2,im2,re3,im3,reX,imX,re_ans,im_ans
  complex(kind=prec) :: a1, a2, a3, x, ans
  call clearcache
  a1 = complex(re1, im1)
  a2 = complex(re2, im2)
  a3 = complex(re3, im3)
  x  = complex(reX, imX)
  ans = G([a1,a2,a3],x)
  re_ans = real(ans)
  im_ans = aimag(ans)
end subroutine G_weight3_fortran

subroutine G_weight4_fortran (re1,im1,re2,im2,re3,im3,re4,im4,reX,imX,re_ans,im_ans) bind(C)
  use handyG
  real(kind=prec) :: re1,im1,re2,im2,re3,im3,re4,im4,reX,imX,re_ans,im_ans
  complex(kind=prec) :: a1, a2, a3, a4, x, ans
  call clearcache
  a1 = complex(re1, im1)
  a2 = complex(re2, im2)
  a3 = complex(re3, im3)
  a4 = complex(re4, im4)
  x  = complex(reX, imX)
  ans = G([a1,a2,a3,a4],x)
  re_ans = real(ans)
  im_ans = aimag(ans)
end subroutine G_weight4_fortran
