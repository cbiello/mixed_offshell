  if (&
       sec .eq. 'rr4252a' .or. sec .eq. 'rr4151a' .or. &
       sec .eq. 'rr4252b' .or. sec .eq. 'rr4151b' .or. &
       sec .eq. 'rr4252c' .or. sec .eq. 'rr4151c' .or. &
       sec .eq. 'rr4252d' .or. sec .eq. 'rr4151d' ) then

     print *, 'allocate NNLO rr for gq: RR'

     allocate(MyKinSet%AllKinConfig(8))

     allocate(MyKinSet%AllKinConfig(1)%LimMom(4,6))    
     allocate(MyKinSet%AllKinConfig(1)%AmpMom(4,6))    
     allocate(MyKinSet%AllKinConfig(1)%KinInv(6))
     allocate(MyKinSet%AllKinConfig(1)%EikInv(3))
     !
     allocate(MyKinSet%AllKinConfig(2)%LimMom(4,4))  
     allocate(MyKinSet%AllKinConfig(2)%AmpMom(4,5))  
     allocate(MyKinSet%AllKinConfig(2)%KinInv(6))
     allocate(MyKinSet%AllKinConfig(2)%EikInv(3))
     !
     allocate(MyKinSet%AllKinConfig(3)%LimMom(4,4))  
     allocate(MyKinSet%AllKinConfig(3)%AmpMom(4,4))  
     allocate(MyKinSet%AllKinConfig(3)%KinInv(6))
     allocate(MyKinSet%AllKinConfig(3)%EikInv(3))
     !
     allocate(MyKinSet%AllKinConfig(4)%LimMom(4,4))  
     allocate(MyKinSet%AllKinConfig(4)%AmpMom(4,4))  
     allocate(MyKinSet%AllKinConfig(4)%KinInv(6))
     allocate(MyKinSet%AllKinConfig(4)%EikInv(3))
     !
     allocate(MyKinSet%AllKinConfig(5)%LimMom(4,4))  
     allocate(MyKinSet%AllKinConfig(5)%AmpMom(4,5))  
     allocate(MyKinSet%AllKinConfig(5)%KinInv(6))
     allocate(MyKinSet%AllKinConfig(5)%EikInv(3))
     
     allocate(MyKinSet%AllKinConfig(6)%LimMom(4,4))  
     allocate(MyKinSet%AllKinConfig(6)%AmpMom(4,5))  
     allocate(MyKinSet%AllKinConfig(6)%KinInv(6))
     allocate(MyKinSet%AllKinConfig(6)%EikInv(3))
     !
     allocate(MyKinSet%AllKinConfig(7)%LimMom(4,4))  
     allocate(MyKinSet%AllKinConfig(7)%AmpMom(4,4))  
     allocate(MyKinSet%AllKinConfig(7)%KinInv(6))
     allocate(MyKinSet%AllKinConfig(7)%EikInv(3))
     !
     allocate(MyKinSet%AllKinConfig(8)%LimMom(4,4))  
     allocate(MyKinSet%AllKinConfig(8)%AmpMom(4,5))  
     allocate(MyKinSet%AllKinConfig(8)%KinInv(6))
     allocate(MyKinSet%AllKinConfig(8)%EikInv(3))

     do j=1,8
        MyKinSet%AllKinConfig(j)%NumLept=2
     enddo

  elseif (sec .eq. 'rr4251' .or. sec .eq. 'rr4152' ) then

     print *, 'allocate NNLO rr for gq: DC'

     allocate(MyKinSet%AllKinConfig(8))

     allocate(MyKinSet%AllKinConfig(1)%LimMom(4,6))    ! HARD
     allocate(MyKinSet%AllKinConfig(1)%AmpMom(4,6))    
     allocate(MyKinSet%AllKinConfig(1)%KinInv(6))
     allocate(MyKinSet%AllKinConfig(1)%EikInv(3))
     !
     allocate(MyKinSet%AllKinConfig(2)%LimMom(4,4))     !S5
     allocate(MyKinSet%AllKinConfig(2)%AmpMom(4,5))  
     allocate(MyKinSet%AllKinConfig(2)%KinInv(6))
     allocate(MyKinSet%AllKinConfig(2)%EikInv(3))
     !
     allocate(MyKinSet%AllKinConfig(3)%LimMom(4,4))     ! C5
     allocate(MyKinSet%AllKinConfig(3)%AmpMom(4,5))  
     allocate(MyKinSet%AllKinConfig(3)%KinInv(6))
     allocate(MyKinSet%AllKinConfig(3)%EikInv(3))
     !
     allocate(MyKinSet%AllKinConfig(4)%LimMom(4,4))    !C5S5
     allocate(MyKinSet%AllKinConfig(4)%AmpMom(4,5))  
     allocate(MyKinSet%AllKinConfig(4)%KinInv(6))
     allocate(MyKinSet%AllKinConfig(4)%EikInv(3))
     !
     allocate(MyKinSet%AllKinConfig(5)%LimMom(4,4))     !C4
     allocate(MyKinSet%AllKinConfig(5)%AmpMom(4,5))  
     allocate(MyKinSet%AllKinConfig(5)%KinInv(6))
     allocate(MyKinSet%AllKinConfig(5)%EikInv(3))
     !
     allocate(MyKinSet%AllKinConfig(6)%LimMom(4,4))     !C4S5
     allocate(MyKinSet%AllKinConfig(6)%AmpMom(4,4))  
     allocate(MyKinSet%AllKinConfig(6)%KinInv(6))
     allocate(MyKinSet%AllKinConfig(6)%EikInv(3))
     !
     allocate(MyKinSet%AllKinConfig(7)%LimMom(4,4))     !C4C5
     allocate(MyKinSet%AllKinConfig(7)%AmpMom(4,4))  
     allocate(MyKinSet%AllKinConfig(7)%KinInv(6))
     allocate(MyKinSet%AllKinConfig(7)%EikInv(3))
     !
     allocate(MyKinSet%AllKinConfig(8)%LimMom(4,4))     !C4C5S5
     allocate(MyKinSet%AllKinConfig(8)%AmpMom(4,4))  
     allocate(MyKinSet%AllKinConfig(8)%KinInv(6))
     allocate(MyKinSet%AllKinConfig(8)%EikInv(3))

     do j=1,8
        MyKinSet%AllKinConfig(j)%NumLept=2
     enddo

  else

     print *, 'wrong allocation for rr gq sectors'
     stop

  endif
