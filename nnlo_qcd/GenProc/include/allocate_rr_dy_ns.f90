  if (&
       sec .eq. 'rr4151a' .or. sec .eq. 'rr4252a' .or. &
       sec .eq. 'rr4151b' .or. sec .eq. 'rr4252b' .or. &
       sec .eq. 'rr4151c' .or. sec .eq. 'rr4252c' .or. &
       sec .eq. 'rr4151d' .or. sec .eq. 'rr4252d' ) then

     print *, 'allocate NNLO rr for qqb: TC'

     allocate(MyKinSet%AllKinConfig(16))
     !
     allocate(MyKinSet%AllKinConfig(1)%LimMom(4,6))    
     allocate(MyKinSet%AllKinConfig(1)%AmpMom(4,6))    
     allocate(MyKinSet%AllKinConfig(1)%KinInv(5))
     !
     allocate(MyKinSet%AllKinConfig(2)%LimMom(4,4))  
     allocate(MyKinSet%AllKinConfig(2)%AmpMom(4,4))  
     allocate(MyKinSet%AllKinConfig(2)%KinInv(5))
     !
     allocate(MyKinSet%AllKinConfig(3)%LimMom(4,4))  
     allocate(MyKinSet%AllKinConfig(3)%AmpMom(4,5))  
     allocate(MyKinSet%AllKinConfig(3)%KinInv(5))
     !
     allocate(MyKinSet%AllKinConfig(4)%LimMom(4,4))  
     allocate(MyKinSet%AllKinConfig(4)%AmpMom(4,4))  
     allocate(MyKinSet%AllKinConfig(4)%KinInv(5))
     !
     allocate(MyKinSet%AllKinConfig(5)%LimMom(4,4))  
     allocate(MyKinSet%AllKinConfig(5)%AmpMom(4,4))  
     allocate(MyKinSet%AllKinConfig(5)%KinInv(5))
     !
     allocate(MyKinSet%AllKinConfig(6)%LimMom(4,4))  
     allocate(MyKinSet%AllKinConfig(6)%AmpMom(4,4))  
     allocate(MyKinSet%AllKinConfig(6)%KinInv(5))
     !
     allocate(MyKinSet%AllKinConfig(7)%LimMom(4,4))  
     allocate(MyKinSet%AllKinConfig(7)%AmpMom(4,4))  
     allocate(MyKinSet%AllKinConfig(7)%KinInv(5))
     !
     allocate(MyKinSet%AllKinConfig(8)%LimMom(4,4))  
     allocate(MyKinSet%AllKinConfig(8)%AmpMom(4,4))  
     allocate(MyKinSet%AllKinConfig(8)%KinInv(5))
     !
     allocate(MyKinSet%AllKinConfig(9)%LimMom(4,4))  
     allocate(MyKinSet%AllKinConfig(9)%AmpMom(4,5))  
     allocate(MyKinSet%AllKinConfig(9)%KinInv(5))
     !
     if (&
          sec.eq.'rr4151b'.or.sec.eq.'rr4252b' .or. & !-- C54DS
          sec.eq.'rr4151d'.or.sec.eq.'rr4252d') then
        allocate(MyKinSet%AllKinConfig(10)%LimMom(4,3)) 
        allocate(MyKinSet%AllKinConfig(10)%AmpMom(4,4))  
        allocate(MyKinSet%AllKinConfig(10)%KinInv(6))
     else
        allocate(MyKinSet%AllKinConfig(10)%LimMom(4,3)) 
        allocate(MyKinSet%AllKinConfig(10)%AmpMom(4,4))  
        allocate(MyKinSet%AllKinConfig(10)%KinInv(5))
     endif
     !
     allocate(MyKinSet%AllKinConfig(11)%LimMom(4,4))  
     allocate(MyKinSet%AllKinConfig(11)%AmpMom(4,5))  
     allocate(MyKinSet%AllKinConfig(11)%KinInv(5))
     !
     allocate(MyKinSet%AllKinConfig(12)%LimMom(4,3))  
     allocate(MyKinSet%AllKinConfig(12)%AmpMom(4,4))  
     allocate(MyKinSet%AllKinConfig(12)%KinInv(5))
     !
     allocate(MyKinSet%AllKinConfig(13)%LimMom(4,4))  
     allocate(MyKinSet%AllKinConfig(13)%AmpMom(4,4))  
     allocate(MyKinSet%AllKinConfig(13)%KinInv(5))
     !
     allocate(MyKinSet%AllKinConfig(14)%LimMom(4,3))  
     allocate(MyKinSet%AllKinConfig(14)%AmpMom(4,4))  
     allocate(MyKinSet%AllKinConfig(14)%KinInv(5))
     !
     allocate(MyKinSet%AllKinConfig(15)%LimMom(4,4))  
     allocate(MyKinSet%AllKinConfig(15)%AmpMom(4,4))  
     allocate(MyKinSet%AllKinConfig(15)%KinInv(5))
     !
     allocate(MyKinSet%AllKinConfig(16)%LimMom(4,3))  
     allocate(MyKinSet%AllKinConfig(16)%AmpMom(4,4))  
     allocate(MyKinSet%AllKinConfig(16)%KinInv(5))

     do j=1,16
        MyKinSet%AllKinConfig(j)%NumLept=2
     enddo

  elseif (sec .eq. 'rr4152' .or. sec .eq. 'rr4251') then

     print *, 'allocate NNLO rr for gg: DC'

     allocate(MyKinSet%AllKinConfig(16))

     allocate(MyKinSet%AllKinConfig(1)%LimMom(4,6))    ! HARD
     allocate(MyKinSet%AllKinConfig(1)%AmpMom(4,6))    
     allocate(MyKinSet%AllKinConfig(1)%KinInv(5))
     !
     allocate(MyKinSet%AllKinConfig(2)%LimMom(4,4))    ! DS
     allocate(MyKinSet%AllKinConfig(2)%AmpMom(4,4))  
     allocate(MyKinSet%AllKinConfig(2)%KinInv(5))
     !
     allocate(MyKinSet%AllKinConfig(3)%LimMom(4,4))     !S5
     allocate(MyKinSet%AllKinConfig(3)%AmpMom(4,5))  
     allocate(MyKinSet%AllKinConfig(3)%KinInv(5))
     !
     allocate(MyKinSet%AllKinConfig(4)%LimMom(4,4))     !DSS5
     allocate(MyKinSet%AllKinConfig(4)%AmpMom(4,4))  
     allocate(MyKinSet%AllKinConfig(4)%KinInv(5))
     !
     allocate(MyKinSet%AllKinConfig(5)%LimMom(4,4))     !C4
     allocate(MyKinSet%AllKinConfig(5)%AmpMom(4,5))  
     allocate(MyKinSet%AllKinConfig(5)%KinInv(5))
     !
     allocate(MyKinSet%AllKinConfig(6)%LimMom(4,4))      !C4DS 
     allocate(MyKinSet%AllKinConfig(6)%AmpMom(4,5))  
     allocate(MyKinSet%AllKinConfig(6)%KinInv(5))
     !
     allocate(MyKinSet%AllKinConfig(7)%LimMom(4,4))     !C4S5
     allocate(MyKinSet%AllKinConfig(7)%AmpMom(4,4))  
     allocate(MyKinSet%AllKinConfig(7)%KinInv(5))
     !
     allocate(MyKinSet%AllKinConfig(8)%LimMom(4,4))     ! C4DSS5
     allocate(MyKinSet%AllKinConfig(8)%AmpMom(4,4))  
     allocate(MyKinSet%AllKinConfig(8)%KinInv(5))
     !
     allocate(MyKinSet%AllKinConfig(9)%LimMom(4,4))     ! C5
     allocate(MyKinSet%AllKinConfig(9)%AmpMom(4,5))  
     allocate(MyKinSet%AllKinConfig(9)%KinInv(5))
     !
     allocate(MyKinSet%AllKinConfig(10)%LimMom(4,3))    !C5DS
     allocate(MyKinSet%AllKinConfig(10)%AmpMom(4,4))  
     allocate(MyKinSet%AllKinConfig(10)%KinInv(5))
     !
     allocate(MyKinSet%AllKinConfig(11)%LimMom(4,4))    !C5S5
     allocate(MyKinSet%AllKinConfig(11)%AmpMom(4,5))  
     allocate(MyKinSet%AllKinConfig(11)%KinInv(5))
     !
     allocate(MyKinSet%AllKinConfig(12)%LimMom(4,3))    !C5DSS5
     allocate(MyKinSet%AllKinConfig(12)%AmpMom(4,4))  
     allocate(MyKinSet%AllKinConfig(12)%KinInv(5))
     !
     allocate(MyKinSet%AllKinConfig(13)%LimMom(4,4))     !C4C5
     allocate(MyKinSet%AllKinConfig(13)%AmpMom(4,4))  
     allocate(MyKinSet%AllKinConfig(13)%KinInv(5))
     !
     allocate(MyKinSet%AllKinConfig(14)%LimMom(4,3))     ! C4C5DS
     allocate(MyKinSet%AllKinConfig(14)%AmpMom(4,4))  
     allocate(MyKinSet%AllKinConfig(14)%KinInv(5))
     !
     allocate(MyKinSet%AllKinConfig(15)%LimMom(4,4))     !C4C5S5
     allocate(MyKinSet%AllKinConfig(15)%AmpMom(4,4))  
     allocate(MyKinSet%AllKinConfig(15)%KinInv(5))
     !
     allocate(MyKinSet%AllKinConfig(16)%LimMom(4,3))     !C4C5DSS5
     allocate(MyKinSet%AllKinConfig(16)%AmpMom(4,4))  
     allocate(MyKinSet%AllKinConfig(16)%KinInv(5))

     do j=1,16
        MyKinSet%AllKinConfig(j)%NumLept=2
     enddo

  else

     print *, 'wrong allocation for rr sector in gg'
     stop

  endif
