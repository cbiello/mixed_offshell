  print *, 'allocate NNLO rr for qq/qr'

  allocate(MyKinSet%AllKinConfig(7))

  allocate(MyKinSet%AllKinConfig(1)%LimMom(4,6)) !-- HardProc
  allocate(MyKinSet%AllKinConfig(1)%AmpMom(4,6)) 
  allocate(MyKinSet%AllKinConfig(1)%KinInv(5))
  !
  allocate(MyKinSet%AllKinConfig(2)%LimMom(4,5)) !-- C41
  allocate(MyKinSet%AllKinConfig(2)%AmpMom(4,5))  
  allocate(MyKinSet%AllKinConfig(2)%KinInv(5))
  !
  allocate(MyKinSet%AllKinConfig(3)%LimMom(4,5)) !-- C51
  allocate(MyKinSet%AllKinConfig(3)%AmpMom(4,5))  
  allocate(MyKinSet%AllKinConfig(3)%KinInv(5))
  !
  allocate(MyKinSet%AllKinConfig(4)%LimMom(4,5)) !-- C41
  allocate(MyKinSet%AllKinConfig(4)%AmpMom(4,5))  
  allocate(MyKinSet%AllKinConfig(4)%KinInv(5))
  !
  allocate(MyKinSet%AllKinConfig(5)%LimMom(4,5)) !-- C52
  allocate(MyKinSet%AllKinConfig(5)%AmpMom(4,5))  
  allocate(MyKinSet%AllKinConfig(5)%KinInv(5))
  !
  allocate(MyKinSet%AllKinConfig(6)%LimMom(4,4)) !-- C4152
  allocate(MyKinSet%AllKinConfig(6)%AmpMom(4,4))  
  allocate(MyKinSet%AllKinConfig(6)%KinInv(5))
  !
  allocate(MyKinSet%AllKinConfig(7)%LimMom(4,4)) !-- C4251
  allocate(MyKinSet%AllKinConfig(7)%AmpMom(4,4))  
  allocate(MyKinSet%AllKinConfig(7)%KinInv(5))

  do j=1,7
     MyKinSet%AllKinConfig(j)%NumLept=2
  enddo
