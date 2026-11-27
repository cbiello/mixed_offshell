  print *, 'allocate NNLO sub12/VV'

  allocate(MyKinSet%AllKinConfig(1))

  allocate(MyKinSet%AllKinConfig(1)%LimMom(4,4))    
  allocate(MyKinSet%AllKinConfig(1)%AmpMom(4,4))    
  allocate(MyKinSet%AllKinConfig(1)%KinInv(1))

  MyKinSet%AllKinConfig(1)%NumLept=2
