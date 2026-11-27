  print *, 'allocate NNLO sub124 for qg'

  allocate(MyKinSet%AllKinConfig(12))

  allocate(MyKinSet%AllKinConfig(1)%LimMom(4,5))    ! hard qg z124
  allocate(MyKinSet%AllKinConfig(1)%AmpMom(4,5))    
  allocate(MyKinSet%AllKinConfig(1)%KinInv(2))
  !
  allocate(MyKinSet%AllKinConfig(2)%LimMom(4,5))    ! col 41 qg z124
  allocate(MyKinSet%AllKinConfig(2)%AmpMom(4,4))    
  allocate(MyKinSet%AllKinConfig(2)%KinInv(2))

  allocate(MyKinSet%AllKinConfig(3)%LimMom(4,5))    ! hard qg 1z24
  allocate(MyKinSet%AllKinConfig(3)%AmpMom(4,5))    
  allocate(MyKinSet%AllKinConfig(3)%KinInv(2))
  !
  allocate(MyKinSet%AllKinConfig(4)%LimMom(4,5))    ! col 41 qg 1z24
  allocate(MyKinSet%AllKinConfig(4)%AmpMom(4,4))    
  allocate(MyKinSet%AllKinConfig(4)%KinInv(2))


  allocate(MyKinSet%AllKinConfig(5)%LimMom(4,5))    ! hard qg 124
  allocate(MyKinSet%AllKinConfig(5)%AmpMom(4,5))    
  allocate(MyKinSet%AllKinConfig(5)%KinInv(2))
  !
  allocate(MyKinSet%AllKinConfig(6)%LimMom(4,5))    ! col 41 qg 124
  allocate(MyKinSet%AllKinConfig(6)%AmpMom(4,4))    
  allocate(MyKinSet%AllKinConfig(6)%KinInv(2))

  !
  allocate(MyKinSet%AllKinConfig(7)%LimMom(4,5))    ! hard gg z124
  allocate(MyKinSet%AllKinConfig(7)%AmpMom(4,5))    
  allocate(MyKinSet%AllKinConfig(7)%KinInv(2))
  !
  allocate(MyKinSet%AllKinConfig(8)%LimMom(4,5))    ! soft gg z124
  allocate(MyKinSet%AllKinConfig(8)%AmpMom(4,4))  
  allocate(MyKinSet%AllKinConfig(8)%KinInv(2))
  !
  allocate(MyKinSet%AllKinConfig(9)%LimMom(4,5))    ! c41 gg z124
  allocate(MyKinSet%AllKinConfig(9)%AmpMom(4,4))    
  allocate(MyKinSet%AllKinConfig(9)%KinInv(2))
  !
  allocate(MyKinSet%AllKinConfig(10)%LimMom(4,5))   ! c42 gg z124
  allocate(MyKinSet%AllKinConfig(10)%AmpMom(4,4))
  allocate(MyKinSet%AllKinConfig(10)%KinInv(2))
  !
  allocate(MyKinSet%AllKinConfig(11)%LimMom(4,5))   ! soft c41 gg z124  
  allocate(MyKinSet%AllKinConfig(11)%AmpMom(4,4))    
  allocate(MyKinSet%AllKinConfig(11)%KinInv(2))
  !
  allocate(MyKinSet%AllKinConfig(12)%LimMom(4,5))   ! soft c42 gg z124  
  allocate(MyKinSet%AllKinConfig(12)%AmpMom(4,4))
  allocate(MyKinSet%AllKinConfig(12)%KinInv(2))

  do j=1,12
     MyKinSet%AllKinConfig(j)%NumLept=2
  enddo


