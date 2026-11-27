  print *, 'allocate NNLO sub124 for gg'

  allocate(MyKinSet%AllKinConfig(4))

  allocate(MyKinSet%AllKinConfig(1)%LimMom(4,5))    ! hard qg z124
  allocate(MyKinSet%AllKinConfig(1)%AmpMom(4,5))    
  allocate(MyKinSet%AllKinConfig(1)%KinInv(2))
  !
  allocate(MyKinSet%AllKinConfig(2)%LimMom(4,5))    ! col 42 qg z124
  allocate(MyKinSet%AllKinConfig(2)%AmpMom(4,4))    
  allocate(MyKinSet%AllKinConfig(2)%KinInv(2))

  allocate(MyKinSet%AllKinConfig(3)%LimMom(4,5))    ! hard gq 1z24
  allocate(MyKinSet%AllKinConfig(3)%AmpMom(4,5))    
  allocate(MyKinSet%AllKinConfig(3)%KinInv(2))
  !
  allocate(MyKinSet%AllKinConfig(4)%LimMom(4,5))    ! col 41 gq 1z24
  allocate(MyKinSet%AllKinConfig(4)%AmpMom(4,4))    
  allocate(MyKinSet%AllKinConfig(4)%KinInv(2))


  do j=1,4
     MyKinSet%AllKinConfig(j)%NumLept=2
  enddo


