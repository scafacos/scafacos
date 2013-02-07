      ! caladdr version with different calls for fwd usage (shcoord pass2)
      ! caladdr5: including int3x,y,z with int3exis2
      ! caladdr4: optimize load/store without dependencies
      subroutine caladdr1(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       17287
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2, 0,-2, 0,-2, 0,-2, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2,-2, 0, 0,-2,-2, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2,-2,-2,-2, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1, 0,-1, 0,-1, 0,-1, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1,-1, 0, 0,-1,-1, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1,-1,-1,-1, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz      
      
      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(7)),int3z(fiz+viz(7))))
          case(2)
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(7)),fint3z(fiz+viz(7))))
          end select
          vcboxaddr( 20) = vixyz(8)
          vcboxaddr( 21) = vixyz(8) + 8
          vcboxaddr( 22) = vixyz(8) +16
          vcboxaddr( 23) = vixyz(8) +24
          vcboxaddr( 24) = vixyz(8) + 32
          vcboxaddr( 25) = vixyz(8) + 40
          vcboxaddr( 26) = vixyz(8) + 48
          vcboxaddr( 27) = vixyz(8) + 56
          dist =  0+child
          fwdstart = 21
	  case(2)
          select case(int3exist)
          case(1)
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(7)),int3z(fiz+viz5(7))))
          case(2)
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(7)),fint3z(fiz+viz5(7))))
          end select
          vcboxaddr( 20) = vixyz(8)
          vcboxaddr( 21) = vixyz(8) + 1
          vcboxaddr( 22) = vixyz(8) + 2
          vcboxaddr( 23) = vixyz(8) + 3
          vcboxaddr( 24) = vixyz(8) + 4
          vcboxaddr( 25) = vixyz(8) + 5
          vcboxaddr( 26) = vixyz(8) + 6
          vcboxaddr( 27) = vixyz(8) + 7
          
          fwdstart = 21	  
	  end select
        case(2)
	
	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	    
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(7)),int3z(fiz+pviz(7))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(7)),int3z(fiz+pviz(7))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(7)),fint3z(fiz+pviz(7))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(7)),fint3z(fiz+pviz(7))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(2)
          vcboxaddr(  3) = vixyz(2) + poffset(1)
          vcboxaddr(  4) = vixyz(3)
          vcboxaddr(  5) = vixyz(3) + poffset(2)
          vcboxaddr(  6) = vixyz(4)
          vcboxaddr(  7) = vixyz(4) + poffset(1)
          vcboxaddr(  8) = vixyz(4) + poffset(2)
          vcboxaddr(  9) = vixyz(4) + poffset(3)
          vcboxaddr( 10) = vixyz(5)
          vcboxaddr( 11) = vixyz(5) + poffset(4)
          vcboxaddr( 12) = vixyz(6)
          vcboxaddr( 13) = vixyz(6) + poffset(1)
          vcboxaddr( 14) = vixyz(6) + poffset(4)
          vcboxaddr( 15) = vixyz(6) + poffset(5)
          vcboxaddr( 16) = vixyz(7)
          vcboxaddr( 17) = vixyz(7) + poffset(2)
          vcboxaddr( 18) = vixyz(7) + poffset(4)
          vcboxaddr( 19) = vixyz(7) + poffset(6)
          vcboxaddr( 20) = vixyz(8)
          vcboxaddr( 21) = vixyz(8) + poffset(1)
          vcboxaddr( 22) = vixyz(8) + poffset(2)
          vcboxaddr( 23) = vixyz(8) + poffset(3)
          vcboxaddr( 24) = vixyz(8) + poffset(4)
          vcboxaddr( 25) = vixyz(8) + poffset(5)
          vcboxaddr( 26) = vixyz(8) + poffset(6)
          vcboxaddr( 27) = vixyz(8) + poffset(7)
          dist =  0+child
	  fwdstart = 21	 
      end select
      end subroutine caladdr1
      subroutine caladdr2(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       20956
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2, 2,-2, 2,-2, 2,-2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2,-2, 0, 0,-2,-2, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2,-2,-2,-2, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1, 1,-1, 1,-1, 1,-1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1,-1, 0, 0,-1,-1, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1,-1,-1,-1, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)
              vixyz(7)=ior(int3x(fix+vix(7)),ior(int3y(fiy+viy(7)),int3z(fiz+viz(7))))
              vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(7)),int3z(fiz+viz(7))))
          case(2)
              vixyz(7)=ior(fint3x(fix+vix(7)),ior(fint3y(fiy+viy(7)),fint3z(fiz+viz(7))))
              vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(7)),fint3z(fiz+viz(7))))
          end select          
          vcboxaddr( 17) = vixyz(7) + 8
          vcboxaddr( 18) = vixyz(7) + 16
          vcboxaddr( 19) = vixyz(7) + 24
          vcboxaddr( 20) = vixyz(7) + 32
          vcboxaddr( 21) = vixyz(7) + 40
          vcboxaddr( 22) = vixyz(7) + 48
          vcboxaddr( 23) = vixyz(7) + 56
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) +16
          vcboxaddr( 26) = vixyz(8) + 32
          vcboxaddr( 27) = vixyz(8) + 48
          dist =  8+child
          fwdstart = 18
	  case(2)
          select case(int3exist)
          case(1)
              vixyz(7)=ior(int3x(fix+vix5(7)),ior(int3y(fiy+viy5(7)),int3z(fiz+viz5(7))))
              vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(7)),int3z(fiz+viz5(7))))
          case(2)
              vixyz(7)=ior(fint3x(fix+vix5(7)),ior(fint3y(fiy+viy5(7)),fint3z(fiz+viz5(7))))
              vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(7)),fint3z(fiz+viz5(7))))
          end select
          vcboxaddr( 17) = vixyz(7) + 1
          vcboxaddr( 18) = vixyz(7) + 2
          vcboxaddr( 19) = vixyz(7) + 3
          vcboxaddr( 20) = vixyz(7) + 4
          vcboxaddr( 21) = vixyz(7) + 5
          vcboxaddr( 22) = vixyz(7) + 6
          vcboxaddr( 23) = vixyz(7) + 7
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + 2
          vcboxaddr( 26) = vixyz(8) + 4
          vcboxaddr( 27) = vixyz(8) + 6
          
          fwdstart = 18
	  end select
        case(2)

	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	
          select case(int3exist)
          case(1)
              vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
              vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
              vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
              vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
              vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
              vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
              vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(7)),int3z(fiz+pviz(7))))
              vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(7)),int3z(fiz+pviz(7))))
          case(2)
              vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
              vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
              vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
              vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
              vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
              vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
              vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(7)),fint3z(fiz+pviz(7))))
              vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(7)),fint3z(fiz+pviz(7))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(1) + poffset(1)
          vcboxaddr(  3) = vixyz(2)
          vcboxaddr(  4) = vixyz(3)
          vcboxaddr(  5) = vixyz(3) + poffset(1)
          vcboxaddr(  6) = vixyz(3) + poffset(2)
          vcboxaddr(  7) = vixyz(3) + poffset(3)
          vcboxaddr(  8) = vixyz(4)
          vcboxaddr(  9) = vixyz(4) + poffset(2)
          vcboxaddr( 10) = vixyz(5)
          vcboxaddr( 11) = vixyz(5) + poffset(1)
          vcboxaddr( 12) = vixyz(5) + poffset(4)
          vcboxaddr( 13) = vixyz(5) + poffset(5)
          vcboxaddr( 14) = vixyz(6)
          vcboxaddr( 15) = vixyz(6) + poffset(4)
          vcboxaddr( 16) = vixyz(7)
          vcboxaddr( 17) = vixyz(7) + poffset(1)
          vcboxaddr( 18) = vixyz(7) + poffset(2)
          vcboxaddr( 19) = vixyz(7) + poffset(3)
          vcboxaddr( 20) = vixyz(7) + poffset(4)
          vcboxaddr( 21) = vixyz(7) + poffset(5)
          vcboxaddr( 22) = vixyz(7) + poffset(6)
          vcboxaddr( 23) = vixyz(7) + poffset(7)
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + poffset(2)
          vcboxaddr( 26) = vixyz(8) + poffset(4)
          vcboxaddr( 27) = vixyz(8) + poffset(6)
          dist =  8+child
	  fwdstart = 18
      end select
      end subroutine caladdr2
      subroutine caladdr3(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       20347
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2,-2, 0, 0,-2,-2, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2, 2,-2, 2,-2, 2,-2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2,-2,-2,-2, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1,-1, 0, 0,-1,-1, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1, 1,-1, 1,-1, 1,-1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1,-1,-1,-1, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)
            vixyz(7)=ior(int3x(fix+vix(7)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(6)),int3z(fiz+viz(6))))
          case(2)
            vixyz(7)=ior(fint3x(fix+vix(7)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(6)),fint3z(fiz+viz(6))))
          end select          
          vcboxaddr( 18) = vixyz(7) + 16
          vcboxaddr( 19) = vixyz(7) + 24
          vcboxaddr( 20) = vixyz(7) + 32
          vcboxaddr( 21) = vixyz(7) + 40
          vcboxaddr( 22) = vixyz(7) + 48
          vcboxaddr( 23) = vixyz(7) + 56
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + 8
          vcboxaddr( 26) = vixyz(8) + 32
          vcboxaddr( 27) = vixyz(8) + 40
          dist = 16+child
          fwdstart = 19
          case(2)
          select case(int3exist)
          case(1)
            vixyz(7)=ior(int3x(fix+vix5(7)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(6)),int3z(fiz+viz5(6))))
          case(2)
            vixyz(7)=ior(fint3x(fix+vix5(7)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(6)),fint3z(fiz+viz5(6))))
          end select          
          vcboxaddr( 18) = vixyz(7) + 2
          vcboxaddr( 19) = vixyz(7) + 3
          vcboxaddr( 20) = vixyz(7) + 4
          vcboxaddr( 21) = vixyz(7) + 5
          vcboxaddr( 22) = vixyz(7) + 6
          vcboxaddr( 23) = vixyz(7) + 7
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + 1
          vcboxaddr( 26) = vixyz(8) + 4
          vcboxaddr( 27) = vixyz(8) + 5
	  
          fwdstart = 19
	  end select
        case(2)
	
	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(6)),int3z(fiz+pviz(6))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(6)),int3z(fiz+pviz(6))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(6)),fint3z(fiz+pviz(6))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(6)),fint3z(fiz+pviz(6))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(1) + poffset(2)
          vcboxaddr(  3) = vixyz(2)
          vcboxaddr(  4) = vixyz(3)
          vcboxaddr(  5) = vixyz(3) + poffset(1)
          vcboxaddr(  6) = vixyz(3) + poffset(2)
          vcboxaddr(  7) = vixyz(3) + poffset(3)
          vcboxaddr(  8) = vixyz(4)
          vcboxaddr(  9) = vixyz(4) + poffset(1)
          vcboxaddr( 10) = vixyz(5)
          vcboxaddr( 11) = vixyz(5) + poffset(2)
          vcboxaddr( 12) = vixyz(5) + poffset(4)
          vcboxaddr( 13) = vixyz(5) + poffset(6)
          vcboxaddr( 14) = vixyz(6)
          vcboxaddr( 15) = vixyz(6) + poffset(4)
          vcboxaddr( 16) = vixyz(7)
          vcboxaddr( 17) = vixyz(7) + poffset(1)
          vcboxaddr( 18) = vixyz(7) + poffset(2)
          vcboxaddr( 19) = vixyz(7) + poffset(3)
          vcboxaddr( 20) = vixyz(7) + poffset(4)
          vcboxaddr( 21) = vixyz(7) + poffset(5)
          vcboxaddr( 22) = vixyz(7) + poffset(6)
          vcboxaddr( 23) = vixyz(7) + poffset(7)
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + poffset(1)
          vcboxaddr( 26) = vixyz(8) + poffset(4)
          vcboxaddr( 27) = vixyz(8) + poffset(5)

          dist = 16+child
	  fwdstart = 19
      end select
      end subroutine caladdr3
      subroutine caladdr4(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       27157
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2, 2,-2, 2,-2, 2,-2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2,-2, 2, 2,-2,-2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2,-2,-2,-2, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1, 1,-1, 1,-1, 1,-1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1,-1, 1, 1,-1,-1, 1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1,-1,-1,-1, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)
            vixyz(5)=ior(int3x(fix+vix(5)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(6)=ior(int3x(fix+vix(6)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(7)=ior(int3x(fix+vix(7)),ior(int3y(fiy+viy(7)),int3z(fiz+viz(7))))
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(7)),int3z(fiz+viz(7))))
          case(2)
            vixyz(5)=ior(fint3x(fix+vix(5)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(6)=ior(fint3x(fix+vix(6)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(7)=ior(fint3x(fix+vix(7)),ior(fint3y(fiy+viy(7)),fint3z(fiz+viz(7))))
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(7)),fint3z(fiz+viz(7))))
          end select          
          vcboxaddr( 13) = vixyz(5) +24
          vcboxaddr( 14) = vixyz(5) + 32
          vcboxaddr( 15) = vixyz(5) + 40
          vcboxaddr( 16) = vixyz(5) + 48
          vcboxaddr( 17) = vixyz(5) + 56
          vcboxaddr( 18) = vixyz(6)
          vcboxaddr( 19) = vixyz(6) +16
          vcboxaddr( 20) = vixyz(6) + 32
          vcboxaddr( 21) = vixyz(6) + 48
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + 8
          vcboxaddr( 24) = vixyz(7) + 32
          vcboxaddr( 25) = vixyz(7) + 40
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + 32
          dist = 24+child
          fwdstart = 14
          case(2)
          select case(int3exist)
          case(1)          
            vixyz(5)=ior(int3x(fix+vix5(5)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(6)=ior(int3x(fix+vix5(6)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(7)=ior(int3x(fix+vix5(7)),ior(int3y(fiy+viy5(7)),int3z(fiz+viz5(7))))
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(7)),int3z(fiz+viz5(7))))
          case(2)          
            vixyz(5)=ior(fint3x(fix+vix5(5)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(6)=ior(fint3x(fix+vix5(6)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(7)=ior(fint3x(fix+vix5(7)),ior(fint3y(fiy+viy5(7)),fint3z(fiz+viz5(7))))
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(7)),fint3z(fiz+viz5(7))))
          end select         
	  vcboxaddr( 13) = vixyz(5) + 3
          vcboxaddr( 14) = vixyz(5) + 4
          vcboxaddr( 15) = vixyz(5) + 5
          vcboxaddr( 16) = vixyz(5) + 6
          vcboxaddr( 17) = vixyz(5) + 7
          vcboxaddr( 18) = vixyz(6)
          vcboxaddr( 19) = vixyz(6) + 2
          vcboxaddr( 20) = vixyz(6) + 4
          vcboxaddr( 21) = vixyz(6) + 6
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + 1
          vcboxaddr( 24) = vixyz(7) + 4
          vcboxaddr( 25) = vixyz(7) + 5
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + 4
   
          fwdstart = 14
	  end select
        case(2)

	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(7)),int3z(fiz+pviz(7))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(7)),int3z(fiz+pviz(7))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(7)),fint3z(fiz+pviz(7))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(7)),fint3z(fiz+pviz(7))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(1) + poffset(1)
          vcboxaddr(  3) = vixyz(1) + poffset(2)
          vcboxaddr(  4) = vixyz(1) + poffset(3)
          vcboxaddr(  5) = vixyz(2)
          vcboxaddr(  6) = vixyz(2) + poffset(2)
          vcboxaddr(  7) = vixyz(3)
          vcboxaddr(  8) = vixyz(3) + poffset(1)
          vcboxaddr(  9) = vixyz(4)
          vcboxaddr( 10) = vixyz(5)
          vcboxaddr( 11) = vixyz(5) + poffset(1)
          vcboxaddr( 12) = vixyz(5) + poffset(2)
	  vcboxaddr( 13) = vixyz(5) + poffset(3)
          vcboxaddr( 14) = vixyz(5) + poffset(4)
          vcboxaddr( 15) = vixyz(5) + poffset(5)
          vcboxaddr( 16) = vixyz(5) + poffset(6)
          vcboxaddr( 17) = vixyz(5) + poffset(7)
          vcboxaddr( 18) = vixyz(6)
          vcboxaddr( 19) = vixyz(6) + poffset(2)
          vcboxaddr( 20) = vixyz(6) + poffset(4)
          vcboxaddr( 21) = vixyz(6) + poffset(6)
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + poffset(1)
          vcboxaddr( 24) = vixyz(7) + poffset(4)
          vcboxaddr( 25) = vixyz(7) + poffset(5)
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + poffset(4)
	  
          dist = 24+child
	  fwdstart = 14
      end select
      end subroutine caladdr4
      subroutine caladdr5(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       18625
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2,-2, 0, 0,-2,-2, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2,-2,-2,-2, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2, 2,-2, 2,-2, 2,-2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1,-1, 0, 0,-1,-1, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1,-1,-1,-1, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1, 1,-1, 1,-1, 1,-1, 1/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)          
            vixyz(7)=ior(int3x(fix+vix(7)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(6)),int3z(fiz+viz(6))))
          case(2)          
            vixyz(7)=ior(fint3x(fix+vix(7)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(6)),fint3z(fiz+viz(6))))
          end select          
          vcboxaddr( 20) = vixyz(7) + 32
          vcboxaddr( 21) = vixyz(7) + 40
          vcboxaddr( 22) = vixyz(7) + 48
          vcboxaddr( 23) = vixyz(7) + 56
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + 8
          vcboxaddr( 26) = vixyz(8) + 16
          vcboxaddr( 27) = vixyz(8) + 24
          dist = 32+child
          fwdstart = 21
          case(2)
          select case(int3exist)
          case(1)          
            vixyz(7)=ior(int3x(fix+vix5(7)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(6)),int3z(fiz+viz5(6))))
          case(2)          
            vixyz(7)=ior(fint3x(fix+vix5(7)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(6)),fint3z(fiz+viz5(6))))
          end select          
          vcboxaddr( 20) = vixyz(7) + 4
          vcboxaddr( 21) = vixyz(7) + 5
          vcboxaddr( 22) = vixyz(7) + 6
          vcboxaddr( 23) = vixyz(7) + 7
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + 1
          vcboxaddr( 26) = vixyz(8) + 2
          vcboxaddr( 27) = vixyz(8) + 3
          
          fwdstart = 21	
          end select
        case(2)
	
	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(6)),int3z(fiz+pviz(6))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(6)),int3z(fiz+pviz(6))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(6)),fint3z(fiz+pviz(6))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(6)),fint3z(fiz+pviz(6))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(1) + poffset(4)
          vcboxaddr(  3) = vixyz(2)
          vcboxaddr(  4) = vixyz(3)
          vcboxaddr(  5) = vixyz(3) + poffset(1)
          vcboxaddr(  6) = vixyz(3) + poffset(4)
          vcboxaddr(  7) = vixyz(3) + poffset(5)
          vcboxaddr(  8) = vixyz(4)
          vcboxaddr(  9) = vixyz(4) + poffset(1)
          vcboxaddr( 10) = vixyz(5)
          vcboxaddr( 11) = vixyz(5) + poffset(2)
          vcboxaddr( 12) = vixyz(5) + poffset(4)
          vcboxaddr( 13) = vixyz(5) + poffset(6)
          vcboxaddr( 14) = vixyz(6)
          vcboxaddr( 15) = vixyz(6) + poffset(2)
          vcboxaddr( 16) = vixyz(7)
          vcboxaddr( 17) = vixyz(7) + poffset(1)
          vcboxaddr( 18) = vixyz(7) + poffset(2)
          vcboxaddr( 19) = vixyz(7) + poffset(3)
	  vcboxaddr( 20) = vixyz(7) + poffset(4)
          vcboxaddr( 21) = vixyz(7) + poffset(5)
          vcboxaddr( 22) = vixyz(7) + poffset(6)
          vcboxaddr( 23) = vixyz(7) + poffset(7)
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + poffset(1)
          vcboxaddr( 26) = vixyz(8) + poffset(2)
          vcboxaddr( 27) = vixyz(8) + poffset(3)
	  
          dist = 32+child
	  fwdstart = 21
      end select
      end subroutine caladdr5
      subroutine caladdr6(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       24637
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2, 2,-2, 2,-2, 2,-2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2,-2,-2,-2, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2,-2, 2, 2,-2,-2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1, 1,-1, 1,-1, 1,-1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1,-1,-1,-1, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1,-1, 1, 1,-1,-1, 1, 1/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)
            vixyz(5)=ior(int3x(fix+vix(5)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(6)=ior(int3x(fix+vix(6)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(7)=ior(int3x(fix+vix(7)),ior(int3y(fiy+viy(7)),int3z(fiz+viz(7))))
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(7)),int3z(fiz+viz(7))))
          case(2)
            vixyz(5)=ior(fint3x(fix+vix(5)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(6)=ior(fint3x(fix+vix(6)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(7)=ior(fint3x(fix+vix(7)),ior(fint3y(fiy+viy(7)),fint3z(fiz+viz(7))))
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(7)),fint3z(fiz+viz(7))))
          end select          
          vcboxaddr( 15) = vixyz(5) + 40
          vcboxaddr( 16) = vixyz(5) + 48
          vcboxaddr( 17) = vixyz(5) + 56
          vcboxaddr( 18) = vixyz(6)
          vcboxaddr( 19) = vixyz(6) +16
          vcboxaddr( 20) = vixyz(6) + 32
          vcboxaddr( 21) = vixyz(6) + 48
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + 8
          vcboxaddr( 24) = vixyz(7) +16
          vcboxaddr( 25) = vixyz(7) +24
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) +16
          dist = 40+child
          fwdstart =  16	
          case(2)
          select case(int3exist)
          case(1)
            vixyz(5)=ior(int3x(fix+vix5(5)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(6)=ior(int3x(fix+vix5(6)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(7)=ior(int3x(fix+vix5(7)),ior(int3y(fiy+viy5(7)),int3z(fiz+viz5(7))))
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(7)),int3z(fiz+viz5(7))))
          case(2)
            vixyz(5)=ior(fint3x(fix+vix5(5)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(6)=ior(fint3x(fix+vix5(6)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(7)=ior(fint3x(fix+vix5(7)),ior(fint3y(fiy+viy5(7)),fint3z(fiz+viz5(7))))
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(7)),fint3z(fiz+viz5(7))))
          end select
          vcboxaddr( 15) = vixyz(5) + 5
          vcboxaddr( 16) = vixyz(5) + 6
          vcboxaddr( 17) = vixyz(5) + 7
          vcboxaddr( 18) = vixyz(6)
          vcboxaddr( 19) = vixyz(6) + 2
          vcboxaddr( 20) = vixyz(6) + 4
          vcboxaddr( 21) = vixyz(6) + 6
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + 1
          vcboxaddr( 24) = vixyz(7) + 2
          vcboxaddr( 25) = vixyz(7) + 3
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + 2
  
          fwdstart =  16
	  end select
        case(2)

	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(7)),int3z(fiz+pviz(7))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(7)),int3z(fiz+pviz(7))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(7)),fint3z(fiz+pviz(7))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(7)),fint3z(fiz+pviz(7))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(1) + poffset(1)
          vcboxaddr(  3) = vixyz(1) + poffset(4)
          vcboxaddr(  4) = vixyz(1) + poffset(5)
          vcboxaddr(  5) = vixyz(2)
          vcboxaddr(  6) = vixyz(2) + poffset(4)
          vcboxaddr(  7) = vixyz(3)
          vcboxaddr(  8) = vixyz(3) + poffset(1)
          vcboxaddr(  9) = vixyz(4)
          vcboxaddr( 10) = vixyz(5)
          vcboxaddr( 11) = vixyz(5) + poffset(1)
          vcboxaddr( 12) = vixyz(5) + poffset(2)
          vcboxaddr( 13) = vixyz(5) + poffset(3)
          vcboxaddr( 14) = vixyz(5) + poffset(4)
          vcboxaddr( 15) = vixyz(5) + poffset(5)
          vcboxaddr( 16) = vixyz(5) + poffset(6)
          vcboxaddr( 17) = vixyz(5) + poffset(7)
          vcboxaddr( 18) = vixyz(6)
          vcboxaddr( 19) = vixyz(6) + poffset(2)
          vcboxaddr( 20) = vixyz(6) + poffset(4)
          vcboxaddr( 21) = vixyz(6) + poffset(6)
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + poffset(1)
          vcboxaddr( 24) = vixyz(7) + poffset(2)
          vcboxaddr( 25) = vixyz(7) + poffset(3)
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + poffset(2)
	  
          dist = 40+child
	  fwdstart =  16
      end select
      end subroutine caladdr6
      subroutine caladdr7(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       23581
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2,-2,-2,-2, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2, 2,-2, 2,-2, 2,-2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2,-2, 2, 2,-2,-2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1,-1,-1,-1, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1, 1,-1, 1,-1, 1,-1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1,-1, 1, 1,-1,-1, 1, 1/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)
            vixyz(5)=ior(int3x(fix+vix(5)),ior(int3y(fiy+viy(1)),int3z(fiz+viz(1))))
            vixyz(6)=ior(int3x(fix+vix(6)),ior(int3y(fiy+viy(2)),int3z(fiz+viz(2))))
            vixyz(7)=ior(int3x(fix+vix(7)),ior(int3y(fiy+viy(3)),int3z(fiz+viz(3))))
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(4)),int3z(fiz+viz(4))))
          case(2)
            vixyz(5)=ior(fint3x(fix+vix(5)),ior(fint3y(fiy+viy(1)),fint3z(fiz+viz(1))))
            vixyz(6)=ior(fint3x(fix+vix(6)),ior(fint3y(fiy+viy(2)),fint3z(fiz+viz(2))))
            vixyz(7)=ior(fint3x(fix+vix(7)),ior(fint3y(fiy+viy(3)),fint3z(fiz+viz(3))))
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(4)),fint3z(fiz+viz(4))))
          end select
          vcboxaddr( 16) = vixyz(5) + 48
          vcboxaddr( 17) = vixyz(5) + 56
          vcboxaddr( 18) = vixyz(6)
          vcboxaddr( 19) = vixyz(6) + 8
          vcboxaddr( 20) = vixyz(6) + 32
          vcboxaddr( 21) = vixyz(6) + 40
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + 8
          vcboxaddr( 24) = vixyz(7) +16
          vcboxaddr( 25) = vixyz(7) +24
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + 8
          dist = 48+child
          fwdstart = 17	
          case(2)
          select case(int3exist)
          case(1)
            vixyz(5)=ior(int3x(fix+vix5(5)),ior(int3y(fiy+viy5(1)),int3z(fiz+viz5(1))))
            vixyz(6)=ior(int3x(fix+vix5(6)),ior(int3y(fiy+viy5(2)),int3z(fiz+viz5(2))))
            vixyz(7)=ior(int3x(fix+vix5(7)),ior(int3y(fiy+viy5(3)),int3z(fiz+viz5(3))))
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(4)),int3z(fiz+viz5(4))))
          case(2)
            vixyz(5)=ior(fint3x(fix+vix5(5)),ior(fint3y(fiy+viy5(1)),fint3z(fiz+viz5(1))))
            vixyz(6)=ior(fint3x(fix+vix5(6)),ior(fint3y(fiy+viy5(2)),fint3z(fiz+viz5(2))))
            vixyz(7)=ior(fint3x(fix+vix5(7)),ior(fint3y(fiy+viy5(3)),fint3z(fiz+viz5(3))))
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(4)),fint3z(fiz+viz5(4))))
          end select
	  vcboxaddr( 16) = vixyz(5) + 6
          vcboxaddr( 17) = vixyz(5) + 7
          vcboxaddr( 18) = vixyz(6)
          vcboxaddr( 19) = vixyz(6) + 1
          vcboxaddr( 20) = vixyz(6) + 4
          vcboxaddr( 21) = vixyz(6) + 5
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + 1
          vcboxaddr( 24) = vixyz(7) + 2
          vcboxaddr( 25) = vixyz(7) + 3
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + 1
	  
          fwdstart = 17
	  end select
        case(2)
	
	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(4)),int3z(fiz+pviz(4))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(4)),int3z(fiz+pviz(4))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(4)),fint3z(fiz+pviz(4))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(4)),fint3z(fiz+pviz(4))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(1) + poffset(2)
          vcboxaddr(  3) = vixyz(1) + poffset(4)
          vcboxaddr(  4) = vixyz(1) + poffset(6)
          vcboxaddr(  5) = vixyz(2)
          vcboxaddr(  6) = vixyz(2) + poffset(4)
          vcboxaddr(  7) = vixyz(3)
          vcboxaddr(  8) = vixyz(3) + poffset(2)
          vcboxaddr(  9) = vixyz(4)
          vcboxaddr( 10) = vixyz(5)
          vcboxaddr( 11) = vixyz(5) + poffset(1)
          vcboxaddr( 12) = vixyz(5) + poffset(2)
          vcboxaddr( 13) = vixyz(5) + poffset(3)
          vcboxaddr( 14) = vixyz(5) + poffset(4)
          vcboxaddr( 15) = vixyz(5) + poffset(5)
          vcboxaddr( 16) = vixyz(5) + poffset(6)
          vcboxaddr( 17) = vixyz(5) + poffset(7)
          vcboxaddr( 18) = vixyz(6)
          vcboxaddr( 19) = vixyz(6) + poffset(1)
          vcboxaddr( 20) = vixyz(6) + poffset(4)
          vcboxaddr( 21) = vixyz(6) + poffset(5)
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + poffset(1)
          vcboxaddr( 24) = vixyz(7) + poffset(2)
          vcboxaddr( 25) = vixyz(7) + poffset(3)
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + poffset(1)

          dist = 48+child
          fwdstart = 17
      end select
      end subroutine caladdr7
      subroutine caladdr8(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       34339
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2, 2,-2, 2,-2, 2,-2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2,-2, 2, 2,-2,-2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2,-2,-2,-2, 2, 2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1, 1,-1, 1,-1, 1,-1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1,-1, 1, 1,-1,-1, 1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1,-1,-1,-1, 1, 1, 1, 1/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+vix(1)),ior(int3y(fiy+viy(1)),int3z(fiz+viz(1))))
            vixyz(2)=ior(int3x(fix+vix(2)),ior(int3y(fiy+viy(1)),int3z(fiz+viz(1))))
            vixyz(3)=ior(int3x(fix+vix(3)),ior(int3y(fiy+viy(3)),int3z(fiz+viz(3))))
            vixyz(4)=ior(int3x(fix+vix(4)),ior(int3y(fiy+viy(3)),int3z(fiz+viz(3))))
            vixyz(5)=ior(int3x(fix+vix(5)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(6)=ior(int3x(fix+vix(6)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(7)=ior(int3x(fix+vix(7)),ior(int3y(fiy+viy(7)),int3z(fiz+viz(7))))
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(7)),int3z(fiz+viz(7))))
          case(2)
            vixyz(1)=ior(fint3x(fix+vix(1)),ior(fint3y(fiy+viy(1)),fint3z(fiz+viz(1))))
            vixyz(2)=ior(fint3x(fix+vix(2)),ior(fint3y(fiy+viy(1)),fint3z(fiz+viz(1))))
            vixyz(3)=ior(fint3x(fix+vix(3)),ior(fint3y(fiy+viy(3)),fint3z(fiz+viz(3))))
            vixyz(4)=ior(fint3x(fix+vix(4)),ior(fint3y(fiy+viy(3)),fint3z(fiz+viz(3))))
            vixyz(5)=ior(fint3x(fix+vix(5)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(6)=ior(fint3x(fix+vix(6)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(7)=ior(fint3x(fix+vix(7)),ior(fint3y(fiy+viy(7)),fint3z(fiz+viz(7))))
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(7)),fint3z(fiz+viz(7))))
          end select          
          vcboxaddr(  8) = vixyz(1) + 56
          vcboxaddr(  9) = vixyz(2)
          vcboxaddr( 10) = vixyz(2) +16
          vcboxaddr( 11) = vixyz(2) + 32
          vcboxaddr( 12) = vixyz(2) + 48
          vcboxaddr( 13) = vixyz(3)
          vcboxaddr( 14) = vixyz(3) + 8
          vcboxaddr( 15) = vixyz(3) + 32
          vcboxaddr( 16) = vixyz(3) + 40
          vcboxaddr( 17) = vixyz(4)
          vcboxaddr( 18) = vixyz(4) + 32
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + 8
          vcboxaddr( 21) = vixyz(5) +16
          vcboxaddr( 22) = vixyz(5) +24
          vcboxaddr( 23) = vixyz(6)
          vcboxaddr( 24) = vixyz(6) +16
          vcboxaddr( 25) = vixyz(7)
          vcboxaddr( 26) = vixyz(7) + 8
          vcboxaddr( 27) = vixyz(8)
          dist = 56+child
          fwdstart =  9
          case(2)
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+vix5(1)),ior(int3y(fiy+viy5(1)),int3z(fiz+viz5(1))))
            vixyz(2)=ior(int3x(fix+vix5(2)),ior(int3y(fiy+viy5(1)),int3z(fiz+viz5(1))))
            vixyz(3)=ior(int3x(fix+vix5(3)),ior(int3y(fiy+viy5(3)),int3z(fiz+viz5(3))))
            vixyz(4)=ior(int3x(fix+vix5(4)),ior(int3y(fiy+viy5(3)),int3z(fiz+viz5(3))))
            vixyz(5)=ior(int3x(fix+vix5(5)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(6)=ior(int3x(fix+vix5(6)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(7)=ior(int3x(fix+vix5(7)),ior(int3y(fiy+viy5(7)),int3z(fiz+viz5(7))))
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(7)),int3z(fiz+viz5(7))))
          case(2)
            vixyz(1)=ior(fint3x(fix+vix5(1)),ior(fint3y(fiy+viy5(1)),fint3z(fiz+viz5(1))))
            vixyz(2)=ior(fint3x(fix+vix5(2)),ior(fint3y(fiy+viy5(1)),fint3z(fiz+viz5(1))))
            vixyz(3)=ior(fint3x(fix+vix5(3)),ior(fint3y(fiy+viy5(3)),fint3z(fiz+viz5(3))))
            vixyz(4)=ior(fint3x(fix+vix5(4)),ior(fint3y(fiy+viy5(3)),fint3z(fiz+viz5(3))))
            vixyz(5)=ior(fint3x(fix+vix5(5)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(6)=ior(fint3x(fix+vix5(6)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(7)=ior(fint3x(fix+vix5(7)),ior(fint3y(fiy+viy5(7)),fint3z(fiz+viz5(7))))
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(7)),fint3z(fiz+viz5(7))))
          end select          
	  
          vcboxaddr(  8) = vixyz(1) + 7
          vcboxaddr(  9) = vixyz(2)
          vcboxaddr( 10) = vixyz(2) + 2
          vcboxaddr( 11) = vixyz(2) + 4
          vcboxaddr( 12) = vixyz(2) + 6
          vcboxaddr( 13) = vixyz(3)
          vcboxaddr( 14) = vixyz(3) + 1
          vcboxaddr( 15) = vixyz(3) + 4
          vcboxaddr( 16) = vixyz(3) + 5
          vcboxaddr( 17) = vixyz(4)
          vcboxaddr( 18) = vixyz(4) + 4
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + 1
          vcboxaddr( 21) = vixyz(5) + 2
          vcboxaddr( 22) = vixyz(5) + 3
          vcboxaddr( 23) = vixyz(6)
          vcboxaddr( 24) = vixyz(6) + 2
          vcboxaddr( 25) = vixyz(7)
          vcboxaddr( 26) = vixyz(7) + 1
          vcboxaddr( 27) = vixyz(8)
	  
          fwdstart =  9
          end select
        case(2)

	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(7)),int3z(fiz+pviz(7))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(7)),int3z(fiz+pviz(7))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(7)),fint3z(fiz+pviz(7))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(7)),fint3z(fiz+pviz(7))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(1) + poffset(1)
          vcboxaddr(  3) = vixyz(1) + poffset(2)
          vcboxaddr(  4) = vixyz(1) + poffset(3)
          vcboxaddr(  5) = vixyz(1) + poffset(4)
          vcboxaddr(  6) = vixyz(1) + poffset(5)
          vcboxaddr(  7) = vixyz(1) + poffset(6)
          vcboxaddr(  8) = vixyz(1) + poffset(7)
          vcboxaddr(  9) = vixyz(2)
          vcboxaddr( 10) = vixyz(2) + poffset(2)
          vcboxaddr( 11) = vixyz(2) + poffset(4)
          vcboxaddr( 12) = vixyz(2) + poffset(6)
          vcboxaddr( 13) = vixyz(3)
          vcboxaddr( 14) = vixyz(3) + poffset(1)
          vcboxaddr( 15) = vixyz(3) + poffset(4)
          vcboxaddr( 16) = vixyz(3) + poffset(5)
          vcboxaddr( 17) = vixyz(4)
          vcboxaddr( 18) = vixyz(4) + poffset(4)
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + poffset(1)
          vcboxaddr( 21) = vixyz(5) + poffset(2)
          vcboxaddr( 22) = vixyz(5) + poffset(3)
          vcboxaddr( 23) = vixyz(6)
          vcboxaddr( 24) = vixyz(6) + poffset(2)
          vcboxaddr( 25) = vixyz(7)
          vcboxaddr( 26) = vixyz(7) + poffset(1)
          vcboxaddr( 27) = vixyz(8)
          dist = 56+child
          fwdstart =  9
      end select
      end subroutine caladdr8
      subroutine caladdr9(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       23608
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2, 0,-2, 0,-2, 0,-2, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2,-2, 2, 2,-2,-2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2,-2,-2,-2, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1, 0,-1, 0,-1, 0,-1, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1,-1, 1, 1,-1,-1, 1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1,-1,-1,-1, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)
            vixyz(6)=ior(int3x(fix+vix(6)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(7)=ior(int3x(fix+vix(7)),ior(int3y(fiy+viy(7)),int3z(fiz+viz(7))))
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(7)),int3z(fiz+viz(7))))
          case(2)
            vixyz(6)=ior(fint3x(fix+vix(6)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(7)=ior(fint3x(fix+vix(7)),ior(fint3y(fiy+viy(7)),fint3z(fiz+viz(7))))
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(7)),fint3z(fiz+viz(7))))
          end select          
          vcboxaddr( 16) = vixyz(6) +16
          vcboxaddr( 17) = vixyz(6) +24
          vcboxaddr( 18) = vixyz(6) + 32
          vcboxaddr( 19) = vixyz(6) + 40
          vcboxaddr( 20) = vixyz(6) + 48
          vcboxaddr( 21) = vixyz(6) + 56
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + 32
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + 8
          vcboxaddr( 26) = vixyz(8) + 32
          vcboxaddr( 27) = vixyz(8) + 40
          dist = 64+child
          fwdstart = 17

          case(2)
          select case(int3exist)
          case(1)
            vixyz(6)=ior(int3x(fix+vix5(6)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(7)=ior(int3x(fix+vix5(7)),ior(int3y(fiy+viy5(7)),int3z(fiz+viz5(7))))
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(7)),int3z(fiz+viz5(7))))
          case(2)
            vixyz(6)=ior(fint3x(fix+vix5(6)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(7)=ior(fint3x(fix+vix5(7)),ior(fint3y(fiy+viy5(7)),fint3z(fiz+viz5(7))))
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(7)),fint3z(fiz+viz5(7))))
          end select          
          vcboxaddr( 16) = vixyz(6) + 2
          vcboxaddr( 17) = vixyz(6) + 3
          vcboxaddr( 18) = vixyz(6) + 4
          vcboxaddr( 19) = vixyz(6) + 5
          vcboxaddr( 20) = vixyz(6) + 6
          vcboxaddr( 21) = vixyz(6) + 7
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + 4
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + 1
          vcboxaddr( 26) = vixyz(8) + 4
          vcboxaddr( 27) = vixyz(8) + 5
	  
          fwdstart = 17
	  end select
        case(2)
	
	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(7)),int3z(fiz+pviz(7))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(7)),int3z(fiz+pviz(7))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(7)),fint3z(fiz+pviz(7))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(7)),fint3z(fiz+pviz(7))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(1) + poffset(2)
          vcboxaddr(  3) = vixyz(2)
          vcboxaddr(  4) = vixyz(2) + poffset(1)
          vcboxaddr(  5) = vixyz(2) + poffset(2)
          vcboxaddr(  6) = vixyz(2) + poffset(3)
          vcboxaddr(  7) = vixyz(3)
          vcboxaddr(  8) = vixyz(4)
          vcboxaddr(  9) = vixyz(4) + poffset(1)
          vcboxaddr( 10) = vixyz(5)
          vcboxaddr( 11) = vixyz(5) + poffset(2)
          vcboxaddr( 12) = vixyz(5) + poffset(4)
          vcboxaddr( 13) = vixyz(5) + poffset(6)
          vcboxaddr( 14) = vixyz(6)
          vcboxaddr( 15) = vixyz(6) + poffset(1)
          vcboxaddr( 16) = vixyz(6) + poffset(2)
          vcboxaddr( 17) = vixyz(6) + poffset(3)
          vcboxaddr( 18) = vixyz(6) + poffset(4)
          vcboxaddr( 19) = vixyz(6) + poffset(5)
          vcboxaddr( 20) = vixyz(6) + poffset(6)
          vcboxaddr( 21) = vixyz(6) + poffset(7)
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + poffset(4)
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + poffset(1)
          vcboxaddr( 26) = vixyz(8) + poffset(4)
          vcboxaddr( 27) = vixyz(8) + poffset(5)
	  
          dist = 64+child
	  fwdstart = 17
      end select
      end subroutine caladdr9
      subroutine caladdr10(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       27361
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2,-2, 2, 2,-2,-2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2, 2,-2, 2,-2, 2,-2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2,-2,-2,-2, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1,-1, 1, 1,-1,-1, 1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1, 1,-1, 1,-1, 1,-1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1,-1,-1,-1, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)          
            vixyz(5)=ior(int3x(fix+vix(5)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(6)=ior(int3x(fix+vix(6)),ior(int3y(fiy+viy(6)),int3z(fiz+viz(6))))
            vixyz(7)=ior(int3x(fix+vix(7)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(6)),int3z(fiz+viz(6))))
          case(2)          
            vixyz(5)=ior(fint3x(fix+vix(5)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(6)=ior(fint3x(fix+vix(6)),ior(fint3y(fiy+viy(6)),fint3z(fiz+viz(6))))
            vixyz(7)=ior(fint3x(fix+vix(7)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(6)),fint3z(fiz+viz(6))))
          end select          
          vcboxaddr( 13) = vixyz(5) +24
          vcboxaddr( 14) = vixyz(5) + 32
          vcboxaddr( 15) = vixyz(5) + 40
          vcboxaddr( 16) = vixyz(5) + 48
          vcboxaddr( 17) = vixyz(5) + 56
          vcboxaddr( 18) = vixyz(6)
          vcboxaddr( 19) = vixyz(6) + 8
          vcboxaddr( 20) = vixyz(6) + 32
          vcboxaddr( 21) = vixyz(6) + 40
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) +16
          vcboxaddr( 24) = vixyz(7) + 32
          vcboxaddr( 25) = vixyz(7) + 48
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + 32
          dist = 72+child
          fwdstart = 14
          
	  case(2)
          select case(int3exist)
          case(1)          
            vixyz(5)=ior(int3x(fix+vix5(5)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(6)=ior(int3x(fix+vix5(6)),ior(int3y(fiy+viy5(6)),int3z(fiz+viz5(6))))
            vixyz(7)=ior(int3x(fix+vix5(7)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(6)),int3z(fiz+viz5(6))))
          case(2)          
            vixyz(5)=ior(fint3x(fix+vix5(5)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(6)=ior(fint3x(fix+vix5(6)),ior(fint3y(fiy+viy5(6)),fint3z(fiz+viz5(6))))
            vixyz(7)=ior(fint3x(fix+vix5(7)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(6)),fint3z(fiz+viz5(6))))
          end select
          vcboxaddr( 13) = vixyz(5) + 3
          vcboxaddr( 14) = vixyz(5) + 4
          vcboxaddr( 15) = vixyz(5) + 5
          vcboxaddr( 16) = vixyz(5) + 6
          vcboxaddr( 17) = vixyz(5) + 7
          vcboxaddr( 18) = vixyz(6)
          vcboxaddr( 19) = vixyz(6) + 1
          vcboxaddr( 20) = vixyz(6) + 4
          vcboxaddr( 21) = vixyz(6) + 5
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + 2
          vcboxaddr( 24) = vixyz(7) + 4
          vcboxaddr( 25) = vixyz(7) + 6
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + 4
         
          fwdstart = 14
	  end select
	  case(2)

	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	  
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(6)),int3z(fiz+pviz(6))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(6)),int3z(fiz+pviz(6))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(6)),fint3z(fiz+pviz(6))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(6)),fint3z(fiz+pviz(6))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(1) + poffset(1)
          vcboxaddr(  3) = vixyz(1) + poffset(2)
          vcboxaddr(  4) = vixyz(1) + poffset(3)
          vcboxaddr(  5) = vixyz(2)
          vcboxaddr(  6) = vixyz(2) + poffset(1)
          vcboxaddr(  7) = vixyz(3)
          vcboxaddr(  8) = vixyz(3) + poffset(2)
          vcboxaddr(  9) = vixyz(4)
          vcboxaddr( 10) = vixyz(5)
          vcboxaddr( 11) = vixyz(5) + poffset(1)
          vcboxaddr( 12) = vixyz(5) + poffset(2)
          vcboxaddr( 13) = vixyz(5) + poffset(3)
          vcboxaddr( 14) = vixyz(5) + poffset(4)
          vcboxaddr( 15) = vixyz(5) + poffset(5)
          vcboxaddr( 16) = vixyz(5) + poffset(6)
          vcboxaddr( 17) = vixyz(5) + poffset(7)
          vcboxaddr( 18) = vixyz(6)
          vcboxaddr( 19) = vixyz(6) + poffset(1)
          vcboxaddr( 20) = vixyz(6) + poffset(4)
          vcboxaddr( 21) = vixyz(6) + poffset(5)
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + poffset(2)
          vcboxaddr( 24) = vixyz(7) + poffset(4)
          vcboxaddr( 25) = vixyz(7) + poffset(6)
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + poffset(4)
	  
          dist = 72+child
          fwdstart = 14
      end select
      end subroutine caladdr10
      subroutine caladdr11(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       21256
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2, 0,-2, 0,-2, 0,-2, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2,-2,-2,-2, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2,-2, 2, 2,-2,-2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1, 0,-1, 0,-1, 0,-1, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1,-1,-1,-1, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1,-1, 1, 1,-1,-1, 1, 1/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)          
            vixyz(6)=ior(int3x(fix+vix(6)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(7)=ior(int3x(fix+vix(7)),ior(int3y(fiy+viy(7)),int3z(fiz+viz(7))))
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(7)),int3z(fiz+viz(7))))
          case(2)          
            vixyz(6)=ior(fint3x(fix+vix(6)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(7)=ior(fint3x(fix+vix(7)),ior(fint3y(fiy+viy(7)),fint3z(fiz+viz(7))))
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(7)),fint3z(fiz+viz(7))))          
          end select
          vcboxaddr( 18) = vixyz(6) + 32
          vcboxaddr( 19) = vixyz(6) + 40
          vcboxaddr( 20) = vixyz(6) + 48
          vcboxaddr( 21) = vixyz(6) + 56
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) +16
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + 8
          vcboxaddr( 26) = vixyz(8) +16
          vcboxaddr( 27) = vixyz(8) +24
          dist = 80+child
          fwdstart = 19
          case(2)
          select case(int3exist)
          case(1)         
            vixyz(6)=ior(int3x(fix+vix5(6)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(7)=ior(int3x(fix+vix5(7)),ior(int3y(fiy+viy5(7)),int3z(fiz+viz5(7))))
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(7)),int3z(fiz+viz5(7))))
          case(2)         
            vixyz(6)=ior(fint3x(fix+vix5(6)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(7)=ior(fint3x(fix+vix5(7)),ior(fint3y(fiy+viy5(7)),fint3z(fiz+viz5(7))))
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(7)),fint3z(fiz+viz5(7))))
          end select
          vcboxaddr( 18) = vixyz(6) + 4
          vcboxaddr( 19) = vixyz(6) + 5
          vcboxaddr( 20) = vixyz(6) + 6
          vcboxaddr( 21) = vixyz(6) + 7
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + 2
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + 1
          vcboxaddr( 26) = vixyz(8) + 2
          vcboxaddr( 27) = vixyz(8) + 3
          
          fwdstart = 19
	  end select
        case(2)
	
	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(7)),int3z(fiz+pviz(7))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(7)),int3z(fiz+pviz(7))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(7)),fint3z(fiz+pviz(7))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(7)),fint3z(fiz+pviz(7))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(1) + poffset(4)
          vcboxaddr(  3) = vixyz(2)
          vcboxaddr(  4) = vixyz(2) + poffset(1)
          vcboxaddr(  5) = vixyz(2) + poffset(4)
          vcboxaddr(  6) = vixyz(2) + poffset(5)
          vcboxaddr(  7) = vixyz(3)
          vcboxaddr(  8) = vixyz(4)
          vcboxaddr(  9) = vixyz(4) + poffset(1)
          vcboxaddr( 10) = vixyz(5)
          vcboxaddr( 11) = vixyz(5) + poffset(2)
          vcboxaddr( 12) = vixyz(5) + poffset(4)
          vcboxaddr( 13) = vixyz(5) + poffset(6)
          vcboxaddr( 14) = vixyz(6)
          vcboxaddr( 15) = vixyz(6) + poffset(1)
          vcboxaddr( 16) = vixyz(6) + poffset(2)
          vcboxaddr( 17) = vixyz(6) + poffset(3)
          vcboxaddr( 18) = vixyz(6) + poffset(4)
          vcboxaddr( 19) = vixyz(6) + poffset(5)
          vcboxaddr( 20) = vixyz(6) + poffset(6)
          vcboxaddr( 21) = vixyz(6) + poffset(7)
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + poffset(2)
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + poffset(1)
          vcboxaddr( 26) = vixyz(8) + poffset(2)
          vcboxaddr( 27) = vixyz(8) + poffset(3)
	  
          dist = 80+child
          fwdstart = 19
      end select
      end subroutine caladdr11
      subroutine caladdr12(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       25249
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2,-2, 2, 2,-2,-2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2,-2,-2,-2, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2, 2,-2, 2,-2, 2,-2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1,-1, 1, 1,-1,-1, 1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1,-1,-1,-1, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1, 1,-1, 1,-1, 1,-1, 1/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)
            vixyz(5)=ior(int3x(fix+vix(5)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(6)=ior(int3x(fix+vix(6)),ior(int3y(fiy+viy(6)),int3z(fiz+viz(6))))
            vixyz(7)=ior(int3x(fix+vix(7)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(6)),int3z(fiz+viz(6))))
          case(2)
            vixyz(5)=ior(fint3x(fix+vix(5)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(6)=ior(fint3x(fix+vix(6)),ior(fint3y(fiy+viy(6)),fint3z(fiz+viz(6))))
            vixyz(7)=ior(fint3x(fix+vix(7)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(6)),fint3z(fiz+viz(6))))   
          end select       
          vcboxaddr( 15) = vixyz(5) + 40
          vcboxaddr( 16) = vixyz(5) + 48
          vcboxaddr( 17) = vixyz(5) + 56
          vcboxaddr( 18) = vixyz(6)
          vcboxaddr( 19) = vixyz(6) + 8
          vcboxaddr( 20) = vixyz(6) +16
          vcboxaddr( 21) = vixyz(6) +24
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) +16
          vcboxaddr( 24) = vixyz(7) + 32
          vcboxaddr( 25) = vixyz(7) + 48
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) +16
          dist = 88+child
          fwdstart =  16
          case(2)
          select case(int3exist)
          case(1)
            vixyz(5)=ior(int3x(fix+vix5(5)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(6)=ior(int3x(fix+vix5(6)),ior(int3y(fiy+viy5(6)),int3z(fiz+viz5(6))))
            vixyz(7)=ior(int3x(fix+vix5(7)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(6)),int3z(fiz+viz5(6))))
          case(2)
            vixyz(5)=ior(fint3x(fix+vix5(5)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(6)=ior(fint3x(fix+vix5(6)),ior(fint3y(fiy+viy5(6)),fint3z(fiz+viz5(6))))
            vixyz(7)=ior(fint3x(fix+vix5(7)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(6)),fint3z(fiz+viz5(6))))
          end select
          vcboxaddr( 15) = vixyz(5) + 5
          vcboxaddr( 16) = vixyz(5) + 6
          vcboxaddr( 17) = vixyz(5) + 7
          vcboxaddr( 18) = vixyz(6)
          vcboxaddr( 19) = vixyz(6) + 1
          vcboxaddr( 20) = vixyz(6) + 2
          vcboxaddr( 21) = vixyz(6) + 3
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + 2
          vcboxaddr( 24) = vixyz(7) + 4
          vcboxaddr( 25) = vixyz(7) + 6
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + 2
          
          fwdstart =  16
	  end select
        case(2)

	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(6)),int3z(fiz+pviz(6))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(6)),int3z(fiz+pviz(6))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(6)),fint3z(fiz+pviz(6))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(6)),fint3z(fiz+pviz(6))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(1) + poffset(1)
          vcboxaddr(  3) = vixyz(1) + poffset(4)
          vcboxaddr(  4) = vixyz(1) + poffset(5)
          vcboxaddr(  5) = vixyz(2)
          vcboxaddr(  6) = vixyz(2) + poffset(1)
          vcboxaddr(  7) = vixyz(3)
          vcboxaddr(  8) = vixyz(3) + poffset(4)
          vcboxaddr(  9) = vixyz(4)
          vcboxaddr( 10) = vixyz(5)
          vcboxaddr( 11) = vixyz(5) + poffset(1)
          vcboxaddr( 12) = vixyz(5) + poffset(2)
          vcboxaddr( 13) = vixyz(5) + poffset(3)
          vcboxaddr( 14) = vixyz(5) + poffset(4)
          vcboxaddr( 15) = vixyz(5) + poffset(5)
          vcboxaddr( 16) = vixyz(5) + poffset(6)
          vcboxaddr( 17) = vixyz(5) + poffset(7)
          vcboxaddr( 18) = vixyz(6)
          vcboxaddr( 19) = vixyz(6) + poffset(1)
          vcboxaddr( 20) = vixyz(6) + poffset(2)
          vcboxaddr( 21) = vixyz(6) + poffset(3)
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + poffset(2)
          vcboxaddr( 24) = vixyz(7) + poffset(4)
          vcboxaddr( 25) = vixyz(7) + poffset(6)
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + poffset(2)
	  
          dist = 88+child
          fwdstart =  16
      end select
      end subroutine caladdr12
      subroutine caladdr13(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       31222
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2, 0,-2, 0,-2, 0,-2, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2,-2, 2, 2,-2,-2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2,-2,-2,-2, 2, 2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1, 0,-1, 0,-1, 0,-1, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1,-1, 1, 1,-1,-1, 1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1,-1,-1,-1, 1, 1, 1, 1/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)
            vixyz(2)=ior(int3x(fix+vix(2)),ior(int3y(fiy+viy(1)),int3z(fiz+viz(1))))
            vixyz(3)=ior(int3x(fix+vix(3)),ior(int3y(fiy+viy(3)),int3z(fiz+viz(3))))
            vixyz(4)=ior(int3x(fix+vix(4)),ior(int3y(fiy+viy(3)),int3z(fiz+viz(3))))
            vixyz(5)=ior(int3x(fix+vix(5)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(6)=ior(int3x(fix+vix(6)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(7)=ior(int3x(fix+vix(7)),ior(int3y(fiy+viy(7)),int3z(fiz+viz(7))))
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(7)),int3z(fiz+viz(7))))
          case(2)
            vixyz(2)=ior(fint3x(fix+vix(2)),ior(fint3y(fiy+viy(1)),fint3z(fiz+viz(1))))
            vixyz(3)=ior(fint3x(fix+vix(3)),ior(fint3y(fiy+viy(3)),fint3z(fiz+viz(3))))
            vixyz(4)=ior(fint3x(fix+vix(4)),ior(fint3y(fiy+viy(3)),fint3z(fiz+viz(3))))
            vixyz(5)=ior(fint3x(fix+vix(5)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(6)=ior(fint3x(fix+vix(6)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(7)=ior(fint3x(fix+vix(7)),ior(fint3y(fiy+viy(7)),fint3z(fiz+viz(7))))
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(7)),fint3z(fiz+viz(7))))
          end select          
          vcboxaddr( 11) = vixyz(2) + 48
          vcboxaddr( 12) = vixyz(2) + 56
          vcboxaddr( 13) = vixyz(3)
          vcboxaddr( 14) = vixyz(3) + 32
          vcboxaddr( 15) = vixyz(4)
          vcboxaddr( 16) = vixyz(4) + 8
          vcboxaddr( 17) = vixyz(4) + 32
          vcboxaddr( 18) = vixyz(4) + 40
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) +16
          vcboxaddr( 21) = vixyz(6)
          vcboxaddr( 22) = vixyz(6) + 8
          vcboxaddr( 23) = vixyz(6) +16
          vcboxaddr( 24) = vixyz(6) +24
          vcboxaddr( 25) = vixyz(7)
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + 8
          dist = 96+child
          fwdstart = 12
          case(2)
          select case(int3exist)
          case(1)
            vixyz(2)=ior(int3x(fix+vix5(2)),ior(int3y(fiy+viy5(1)),int3z(fiz+viz5(1))))
            vixyz(3)=ior(int3x(fix+vix5(3)),ior(int3y(fiy+viy5(3)),int3z(fiz+viz5(3))))
            vixyz(4)=ior(int3x(fix+vix5(4)),ior(int3y(fiy+viy5(3)),int3z(fiz+viz5(3))))
            vixyz(5)=ior(int3x(fix+vix5(5)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(6)=ior(int3x(fix+vix5(6)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(7)=ior(int3x(fix+vix5(7)),ior(int3y(fiy+viy5(7)),int3z(fiz+viz5(7))))
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(7)),int3z(fiz+viz5(7))))
          case(2)
            vixyz(2)=ior(fint3x(fix+vix5(2)),ior(fint3y(fiy+viy5(1)),fint3z(fiz+viz5(1))))
            vixyz(3)=ior(fint3x(fix+vix5(3)),ior(fint3y(fiy+viy5(3)),fint3z(fiz+viz5(3))))
            vixyz(4)=ior(fint3x(fix+vix5(4)),ior(fint3y(fiy+viy5(3)),fint3z(fiz+viz5(3))))
            vixyz(5)=ior(fint3x(fix+vix5(5)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(6)=ior(fint3x(fix+vix5(6)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(7)=ior(fint3x(fix+vix5(7)),ior(fint3y(fiy+viy5(7)),fint3z(fiz+viz5(7))))
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(7)),fint3z(fiz+viz5(7))))  
          end select        
          vcboxaddr( 11) = vixyz(2) + 6
          vcboxaddr( 12) = vixyz(2) + 7
          vcboxaddr( 13) = vixyz(3)
          vcboxaddr( 14) = vixyz(3) + 4
          vcboxaddr( 15) = vixyz(4)
          vcboxaddr( 16) = vixyz(4) + 1
          vcboxaddr( 17) = vixyz(4) + 4
          vcboxaddr( 18) = vixyz(4) + 5
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + 2
          vcboxaddr( 21) = vixyz(6)
          vcboxaddr( 22) = vixyz(6) + 1
          vcboxaddr( 23) = vixyz(6) + 2
          vcboxaddr( 24) = vixyz(6) + 3
          vcboxaddr( 25) = vixyz(7)
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + 1
	            
          fwdstart = 12
	  end select
	  case(2)
	  
	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	  
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(7)),int3z(fiz+pviz(7))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(7)),int3z(fiz+pviz(7))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(7)),fint3z(fiz+pviz(7))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(7)),fint3z(fiz+pviz(7))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(1) + poffset(2)
          vcboxaddr(  3) = vixyz(1) + poffset(4)
          vcboxaddr(  4) = vixyz(1) + poffset(6)
          vcboxaddr(  5) = vixyz(2)
          vcboxaddr(  6) = vixyz(2) + poffset(1)
          vcboxaddr(  7) = vixyz(2) + poffset(2)
          vcboxaddr(  8) = vixyz(2) + poffset(3)
          vcboxaddr(  9) = vixyz(2) + poffset(4)
          vcboxaddr( 10) = vixyz(2) + poffset(5)
	  vcboxaddr( 11) = vixyz(2) + poffset(6)
          vcboxaddr( 12) = vixyz(2) + poffset(7)
          vcboxaddr( 13) = vixyz(3)
          vcboxaddr( 14) = vixyz(3) + poffset(4)
          vcboxaddr( 15) = vixyz(4)
          vcboxaddr( 16) = vixyz(4) + poffset(1)
          vcboxaddr( 17) = vixyz(4) + poffset(4)
          vcboxaddr( 18) = vixyz(4) + poffset(5)
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + poffset(2)
          vcboxaddr( 21) = vixyz(6)
          vcboxaddr( 22) = vixyz(6) + poffset(1)
          vcboxaddr( 23) = vixyz(6) + poffset(2)
          vcboxaddr( 24) = vixyz(6) + poffset(3)
          vcboxaddr( 25) = vixyz(7)
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + poffset(1)
	  
          dist = 96+child
          fwdstart = 12
      end select
      end subroutine caladdr13
      subroutine caladdr14(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       35947
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2,-2,-2,-2, 2, 2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2, 2,-2, 2,-2, 2,-2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2,-2, 2, 2,-2,-2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1,-1,-1,-1, 1, 1, 1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1, 1,-1, 1,-1, 1,-1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1,-1, 1, 1,-1,-1, 1, 1/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+vix(1)),ior(int3y(fiy+viy(1)),int3z(fiz+viz(1))))
            vixyz(2)=ior(int3x(fix+vix(2)),ior(int3y(fiy+viy(2)),int3z(fiz+viz(2))))
            vixyz(3)=ior(int3x(fix+vix(3)),ior(int3y(fiy+viy(3)),int3z(fiz+viz(3))))
            vixyz(4)=ior(int3x(fix+vix(4)),ior(int3y(fiy+viy(4)),int3z(fiz+viz(4))))
            vixyz(5)=ior(int3x(fix+vix(5)),ior(int3y(fiy+viy(1)),int3z(fiz+viz(1))))
            vixyz(6)=ior(int3x(fix+vix(6)),ior(int3y(fiy+viy(2)),int3z(fiz+viz(2))))
            vixyz(7)=ior(int3x(fix+vix(7)),ior(int3y(fiy+viy(3)),int3z(fiz+viz(3))))
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(4)),int3z(fiz+viz(4))))
           case(2)
            vixyz(1)=ior(fint3x(fix+vix(1)),ior(fint3y(fiy+viy(1)),fint3z(fiz+viz(1))))
            vixyz(2)=ior(fint3x(fix+vix(2)),ior(fint3y(fiy+viy(2)),fint3z(fiz+viz(2))))
            vixyz(3)=ior(fint3x(fix+vix(3)),ior(fint3y(fiy+viy(3)),fint3z(fiz+viz(3))))
            vixyz(4)=ior(fint3x(fix+vix(4)),ior(fint3y(fiy+viy(4)),fint3z(fiz+viz(4))))
            vixyz(5)=ior(fint3x(fix+vix(5)),ior(fint3y(fiy+viy(1)),fint3z(fiz+viz(1))))
            vixyz(6)=ior(fint3x(fix+vix(6)),ior(fint3y(fiy+viy(2)),fint3z(fiz+viz(2))))
            vixyz(7)=ior(fint3x(fix+vix(7)),ior(fint3y(fiy+viy(3)),fint3z(fiz+viz(3))))
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(4)),fint3z(fiz+viz(4))))     
          end select    
          vcboxaddr(  8) = vixyz(1) + 56
          vcboxaddr(  9) = vixyz(2)
          vcboxaddr( 10) = vixyz(2) + 8
          vcboxaddr( 11) = vixyz(2) + 32
          vcboxaddr( 12) = vixyz(2) + 40
          vcboxaddr( 13) = vixyz(3)
          vcboxaddr( 14) = vixyz(3) + 8
          vcboxaddr( 15) = vixyz(3) +16
          vcboxaddr( 16) = vixyz(3) +24
          vcboxaddr( 17) = vixyz(4)
          vcboxaddr( 18) = vixyz(4) + 8
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) +16
          vcboxaddr( 21) = vixyz(5) + 32
          vcboxaddr( 22) = vixyz(5) + 48
          vcboxaddr( 23) = vixyz(6)
          vcboxaddr( 24) = vixyz(6) + 32
          vcboxaddr( 25) = vixyz(7)
          vcboxaddr( 26) = vixyz(7) +16
          vcboxaddr( 27) = vixyz(8)
          dist =104+child
          fwdstart =  9

          case(2)
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+vix5(1)),ior(int3y(fiy+viy5(1)),int3z(fiz+viz5(1))))
            vixyz(2)=ior(int3x(fix+vix5(2)),ior(int3y(fiy+viy5(2)),int3z(fiz+viz5(2))))
            vixyz(3)=ior(int3x(fix+vix5(3)),ior(int3y(fiy+viy5(3)),int3z(fiz+viz5(3))))
            vixyz(4)=ior(int3x(fix+vix5(4)),ior(int3y(fiy+viy5(4)),int3z(fiz+viz5(4))))
            vixyz(5)=ior(int3x(fix+vix5(5)),ior(int3y(fiy+viy5(1)),int3z(fiz+viz5(1))))
            vixyz(6)=ior(int3x(fix+vix5(6)),ior(int3y(fiy+viy5(2)),int3z(fiz+viz5(2))))
            vixyz(7)=ior(int3x(fix+vix5(7)),ior(int3y(fiy+viy5(3)),int3z(fiz+viz5(3))))
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(4)),int3z(fiz+viz5(4))))
          case(2)
            vixyz(1)=ior(fint3x(fix+vix5(1)),ior(fint3y(fiy+viy5(1)),fint3z(fiz+viz5(1))))
            vixyz(2)=ior(fint3x(fix+vix5(2)),ior(fint3y(fiy+viy5(2)),fint3z(fiz+viz5(2))))
            vixyz(3)=ior(fint3x(fix+vix5(3)),ior(fint3y(fiy+viy5(3)),fint3z(fiz+viz5(3))))
            vixyz(4)=ior(fint3x(fix+vix5(4)),ior(fint3y(fiy+viy5(4)),fint3z(fiz+viz5(4))))
            vixyz(5)=ior(fint3x(fix+vix5(5)),ior(fint3y(fiy+viy5(1)),fint3z(fiz+viz5(1))))
            vixyz(6)=ior(fint3x(fix+vix5(6)),ior(fint3y(fiy+viy5(2)),fint3z(fiz+viz5(2))))
            vixyz(7)=ior(fint3x(fix+vix5(7)),ior(fint3y(fiy+viy5(3)),fint3z(fiz+viz5(3))))
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(4)),fint3z(fiz+viz5(4))))
          end select
          vcboxaddr(  8) = vixyz(1) + 7
          vcboxaddr(  9) = vixyz(2)
          vcboxaddr( 10) = vixyz(2) + 1
          vcboxaddr( 11) = vixyz(2) + 4
          vcboxaddr( 12) = vixyz(2) + 5
          vcboxaddr( 13) = vixyz(3)
          vcboxaddr( 14) = vixyz(3) + 1
          vcboxaddr( 15) = vixyz(3) + 2
          vcboxaddr( 16) = vixyz(3) + 3
          vcboxaddr( 17) = vixyz(4)
          vcboxaddr( 18) = vixyz(4) + 1
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + 2
          vcboxaddr( 21) = vixyz(5) + 4
          vcboxaddr( 22) = vixyz(5) + 6
          vcboxaddr( 23) = vixyz(6)
          vcboxaddr( 24) = vixyz(6) + 4
          vcboxaddr( 25) = vixyz(7)
          vcboxaddr( 26) = vixyz(7) + 2
          vcboxaddr( 27) = vixyz(8)
          
          fwdstart =  9
        end select	  
        case(2)
	
	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(4)),int3z(fiz+pviz(4))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(4)),int3z(fiz+pviz(4))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(4)),fint3z(fiz+pviz(4))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(4)),fint3z(fiz+pviz(4))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(1) + poffset(1)
          vcboxaddr(  3) = vixyz(1) + poffset(2)
          vcboxaddr(  4) = vixyz(1) + poffset(3)
          vcboxaddr(  5) = vixyz(1) + poffset(4)
          vcboxaddr(  6) = vixyz(1) + poffset(5)
          vcboxaddr(  7) = vixyz(1) + poffset(6)
          vcboxaddr(  8) = vixyz(1) + poffset(7)
          vcboxaddr(  9) = vixyz(2)
          vcboxaddr( 10) = vixyz(2) + poffset(1)
          vcboxaddr( 11) = vixyz(2) + poffset(4)
          vcboxaddr( 12) = vixyz(2) + poffset(5)
          vcboxaddr( 13) = vixyz(3)
          vcboxaddr( 14) = vixyz(3) + poffset(1)
          vcboxaddr( 15) = vixyz(3) + poffset(2)
          vcboxaddr( 16) = vixyz(3) + poffset(3)
          vcboxaddr( 17) = vixyz(4)
          vcboxaddr( 18) = vixyz(4) + poffset(1)
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + poffset(2)
          vcboxaddr( 21) = vixyz(5) + poffset(4)
          vcboxaddr( 22) = vixyz(5) + poffset(6)
          vcboxaddr( 23) = vixyz(6)
          vcboxaddr( 24) = vixyz(6) + poffset(4)
          vcboxaddr( 25) = vixyz(7)
          vcboxaddr( 26) = vixyz(7) + poffset(2)
          vcboxaddr( 27) = vixyz(8)

          dist =104+child
          fwdstart =  9
      end select
      end subroutine caladdr14
      subroutine caladdr15(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       17143
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2,-2, 0, 0,-2,-2, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2, 0,-2, 0,-2, 0,-2, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2,-2,-2,-2, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1,-1, 0, 0,-1,-1, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1, 0,-1, 0,-1, 0,-1, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1,-1,-1,-1, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(6)),int3z(fiz+viz(6))))
          case(2)
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(6)),fint3z(fiz+viz(6)))) 
          end select         
          vcboxaddr( 20) = vixyz(8)
          vcboxaddr( 21) = vixyz(8) + 8
          vcboxaddr( 22) = vixyz(8) +16
          vcboxaddr( 23) = vixyz(8) +24
          vcboxaddr( 24) = vixyz(8) + 32
          vcboxaddr( 25) = vixyz(8) + 40
          vcboxaddr( 26) = vixyz(8) + 48
          vcboxaddr( 27) = vixyz(8) + 56
          dist =112+child
          fwdstart = 21

          case(2)
          select case(int3exist)
          case(1)
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(6)),int3z(fiz+viz5(6))))
          case(2)
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(6)),fint3z(fiz+viz5(6))))
          end select
          vcboxaddr( 20) = vixyz(8)
          vcboxaddr( 21) = vixyz(8) + 1
          vcboxaddr( 22) = vixyz(8) + 2
          vcboxaddr( 23) = vixyz(8) + 3
          vcboxaddr( 24) = vixyz(8) + 4
          vcboxaddr( 25) = vixyz(8) + 5
          vcboxaddr( 26) = vixyz(8) + 6
          vcboxaddr( 27) = vixyz(8) + 7          
          
          fwdstart = 21
        end select	  
        case(2)

	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(6)),int3z(fiz+pviz(6))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(6)),int3z(fiz+pviz(6))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(6)),fint3z(fiz+pviz(6))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(6)),fint3z(fiz+pviz(6))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(2)
          vcboxaddr(  3) = vixyz(2) + poffset(2)
          vcboxaddr(  4) = vixyz(3)
          vcboxaddr(  5) = vixyz(3) + poffset(1)
          vcboxaddr(  6) = vixyz(4)
          vcboxaddr(  7) = vixyz(4) + poffset(1)
          vcboxaddr(  8) = vixyz(4) + poffset(2)
          vcboxaddr(  9) = vixyz(4) + poffset(3)
          vcboxaddr( 10) = vixyz(5)
          vcboxaddr( 11) = vixyz(5) + poffset(4)
          vcboxaddr( 12) = vixyz(6)
          vcboxaddr( 13) = vixyz(6) + poffset(2)
          vcboxaddr( 14) = vixyz(6) + poffset(4)
          vcboxaddr( 15) = vixyz(6) + poffset(6)
          vcboxaddr( 16) = vixyz(7)
          vcboxaddr( 17) = vixyz(7) + poffset(1)
          vcboxaddr( 18) = vixyz(7) + poffset(4)
          vcboxaddr( 19) = vixyz(7) + poffset(5)
          vcboxaddr( 20) = vixyz(8)
          vcboxaddr( 21) = vixyz(8) + poffset(1)
          vcboxaddr( 22) = vixyz(8) + poffset(2)
          vcboxaddr( 23) = vixyz(8) + poffset(3)
          vcboxaddr( 24) = vixyz(8) + poffset(4)
          vcboxaddr( 25) = vixyz(8) + poffset(5)
          vcboxaddr( 26) = vixyz(8) + poffset(6)
          vcboxaddr( 27) = vixyz(8) + poffset(7)
	  
          dist =112+child
          fwdstart = 21
      end select
      end subroutine caladdr15
      subroutine caladdr16(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       20689
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2,-2,-2,-2, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2, 0,-2, 0,-2, 0,-2, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2,-2, 2, 2,-2,-2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1,-1,-1,-1, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1, 0,-1, 0,-1, 0,-1, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1,-1, 1, 1,-1,-1, 1, 1/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)          
            vixyz(6)=ior(int3x(fix+vix(6)),ior(int3y(fiy+viy(2)),int3z(fiz+viz(2))))
            vixyz(7)=ior(int3x(fix+vix(7)),ior(int3y(fiy+viy(3)),int3z(fiz+viz(3))))
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(4)),int3z(fiz+viz(4))))
          case(2)          
            vixyz(6)=ior(fint3x(fix+vix(6)),ior(fint3y(fiy+viy(2)),fint3z(fiz+viz(2))))
            vixyz(7)=ior(fint3x(fix+vix(7)),ior(fint3y(fiy+viy(3)),fint3z(fiz+viz(3))))
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(4)),fint3z(fiz+viz(4))))    
          end select      
          vcboxaddr( 18) = vixyz(6) + 32
          vcboxaddr( 19) = vixyz(6) + 40
          vcboxaddr( 20) = vixyz(6) + 48
          vcboxaddr( 21) = vixyz(6) + 56
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + 8
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + 8
          vcboxaddr( 26) = vixyz(8) +16
          vcboxaddr( 27) = vixyz(8) +24
          dist =120+child
          fwdstart = 19	

          case(2)
          select case(int3exist)
          case(1)          
            vixyz(6)=ior(int3x(fix+vix5(6)),ior(int3y(fiy+viy5(2)),int3z(fiz+viz5(2))))
            vixyz(7)=ior(int3x(fix+vix5(7)),ior(int3y(fiy+viy5(3)),int3z(fiz+viz5(3))))
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(4)),int3z(fiz+viz5(4))))
          case(2)          
            vixyz(6)=ior(fint3x(fix+vix5(6)),ior(fint3y(fiy+viy5(2)),fint3z(fiz+viz5(2))))
            vixyz(7)=ior(fint3x(fix+vix5(7)),ior(fint3y(fiy+viy5(3)),fint3z(fiz+viz5(3))))
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(4)),fint3z(fiz+viz5(4))))
          end select
          vcboxaddr( 18) = vixyz(6) + 4
          vcboxaddr( 19) = vixyz(6) + 5
          vcboxaddr( 20) = vixyz(6) + 6
          vcboxaddr( 21) = vixyz(6) + 7
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + 1
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + 1
          vcboxaddr( 26) = vixyz(8) + 2
          vcboxaddr( 27) = vixyz(8) + 3          
          
          fwdstart = 19
        end select	  
        case(2)

	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(4)),int3z(fiz+pviz(4))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(4)),int3z(fiz+pviz(4))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(4)),fint3z(fiz+pviz(4))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(4)),fint3z(fiz+pviz(4))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(1) + poffset(4)
          vcboxaddr(  3) = vixyz(2)
          vcboxaddr(  4) = vixyz(2) + poffset(2)
          vcboxaddr(  5) = vixyz(2) + poffset(4)
          vcboxaddr(  6) = vixyz(2) + poffset(6)
          vcboxaddr(  7) = vixyz(3)
          vcboxaddr(  8) = vixyz(4)
          vcboxaddr(  9) = vixyz(4) + poffset(2)
          vcboxaddr( 10) = vixyz(5)
          vcboxaddr( 11) = vixyz(5) + poffset(1)
          vcboxaddr( 12) = vixyz(5) + poffset(4)
          vcboxaddr( 13) = vixyz(5) + poffset(5)
          vcboxaddr( 14) = vixyz(6)
          vcboxaddr( 15) = vixyz(6) + poffset(1)
          vcboxaddr( 16) = vixyz(6) + poffset(2)
          vcboxaddr( 17) = vixyz(6) + poffset(3)
	  vcboxaddr( 18) = vixyz(6) + poffset(4)
          vcboxaddr( 19) = vixyz(6) + poffset(5)
          vcboxaddr( 20) = vixyz(6) + poffset(6)
          vcboxaddr( 21) = vixyz(6) + poffset(7)
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + poffset(1)
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + poffset(1)
          vcboxaddr( 26) = vixyz(8) + poffset(2)
          vcboxaddr( 27) = vixyz(8) + poffset(3)

          dist =120+child
          fwdstart = 19	
      end select
      end subroutine caladdr16
      subroutine caladdr17(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       29002
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2, 2,-2, 2,-2, 2,-2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2,-2, 0, 0,-2,-2, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2,-2,-2,-2, 2, 2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1, 1,-1, 1,-1, 1,-1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1,-1, 0, 0,-1,-1, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1,-1,-1,-1, 1, 1, 1, 1/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)
            vixyz(3)=ior(int3x(fix+vix(3)),ior(int3y(fiy+viy(3)),int3z(fiz+viz(3))))
            vixyz(4)=ior(int3x(fix+vix(4)),ior(int3y(fiy+viy(3)),int3z(fiz+viz(3))))
            vixyz(5)=ior(int3x(fix+vix(5)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(6)=ior(int3x(fix+vix(6)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(7)=ior(int3x(fix+vix(7)),ior(int3y(fiy+viy(7)),int3z(fiz+viz(7))))
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(7)),int3z(fiz+viz(7))))
          case(2)
            vixyz(3)=ior(fint3x(fix+vix(3)),ior(fint3y(fiy+viy(3)),fint3z(fiz+viz(3))))
            vixyz(4)=ior(fint3x(fix+vix(4)),ior(fint3y(fiy+viy(3)),fint3z(fiz+viz(3))))
            vixyz(5)=ior(fint3x(fix+vix(5)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(6)=ior(fint3x(fix+vix(6)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(7)=ior(fint3x(fix+vix(7)),ior(fint3y(fiy+viy(7)),fint3z(fiz+viz(7))))
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(7)),fint3z(fiz+viz(7))))   
          end select       
          vcboxaddr( 12) = vixyz(3) + 40
          vcboxaddr( 13) = vixyz(3) + 48
          vcboxaddr( 14) = vixyz(3) + 56
          vcboxaddr( 15) = vixyz(4)
          vcboxaddr( 16) = vixyz(4) +16
          vcboxaddr( 17) = vixyz(4) + 32
          vcboxaddr( 18) = vixyz(4) + 48
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + 8
          vcboxaddr( 21) = vixyz(6)
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + 8
          vcboxaddr( 24) = vixyz(7) +16
          vcboxaddr( 25) = vixyz(7) +24
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) +16
          dist =128+child
          fwdstart = 13	

          case(2)
          select case(int3exist)
          case(1)
            vixyz(3)=ior(int3x(fix+vix5(3)),ior(int3y(fiy+viy5(3)),int3z(fiz+viz5(3))))
            vixyz(4)=ior(int3x(fix+vix5(4)),ior(int3y(fiy+viy5(3)),int3z(fiz+viz5(3))))
            vixyz(5)=ior(int3x(fix+vix5(5)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(6)=ior(int3x(fix+vix5(6)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(7)=ior(int3x(fix+vix5(7)),ior(int3y(fiy+viy5(7)),int3z(fiz+viz5(7))))
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(7)),int3z(fiz+viz5(7))))
          case(2)
            vixyz(3)=ior(fint3x(fix+vix5(3)),ior(fint3y(fiy+viy5(3)),fint3z(fiz+viz5(3))))
            vixyz(4)=ior(fint3x(fix+vix5(4)),ior(fint3y(fiy+viy5(3)),fint3z(fiz+viz5(3))))
            vixyz(5)=ior(fint3x(fix+vix5(5)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(6)=ior(fint3x(fix+vix5(6)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(7)=ior(fint3x(fix+vix5(7)),ior(fint3y(fiy+viy5(7)),fint3z(fiz+viz5(7))))
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(7)),fint3z(fiz+viz5(7))))
          end select
          vcboxaddr( 12) = vixyz(3) + 5
          vcboxaddr( 13) = vixyz(3) + 6
          vcboxaddr( 14) = vixyz(3) + 7
          vcboxaddr( 15) = vixyz(4)
          vcboxaddr( 16) = vixyz(4) + 2
          vcboxaddr( 17) = vixyz(4) + 4
          vcboxaddr( 18) = vixyz(4) + 6
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + 1
          vcboxaddr( 21) = vixyz(6)
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + 1
          vcboxaddr( 24) = vixyz(7) + 2
          vcboxaddr( 25) = vixyz(7) + 3
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + 2
          
          fwdstart = 13
        end select	  
        case(2)

	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(7)),int3z(fiz+pviz(7))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(7)),int3z(fiz+pviz(7))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(7)),fint3z(fiz+pviz(7))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(7)),fint3z(fiz+pviz(7))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(1) + poffset(1)
          vcboxaddr(  3) = vixyz(1) + poffset(4)
          vcboxaddr(  4) = vixyz(1) + poffset(5)
          vcboxaddr(  5) = vixyz(2)
          vcboxaddr(  6) = vixyz(2) + poffset(4)
          vcboxaddr(  7) = vixyz(3)
          vcboxaddr(  8) = vixyz(3) + poffset(1)
          vcboxaddr(  9) = vixyz(3) + poffset(2)
          vcboxaddr( 10) = vixyz(3) + poffset(3)
          vcboxaddr( 11) = vixyz(3) + poffset(4)
	  vcboxaddr( 12) = vixyz(3) + poffset(5)
          vcboxaddr( 13) = vixyz(3) + poffset(6)
          vcboxaddr( 14) = vixyz(3) + poffset(7)
          vcboxaddr( 15) = vixyz(4)
          vcboxaddr( 16) = vixyz(4) + poffset(2)
          vcboxaddr( 17) = vixyz(4) + poffset(4)
          vcboxaddr( 18) = vixyz(4) + poffset(6)
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + poffset(1)
          vcboxaddr( 21) = vixyz(6)
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + poffset(1)
          vcboxaddr( 24) = vixyz(7) + poffset(2)
          vcboxaddr( 25) = vixyz(7) + poffset(3)
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + poffset(2)
	  
          dist =128+child
          fwdstart = 13	
      end select
      end subroutine caladdr17
      subroutine caladdr18(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       23989
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2,-2,-2,-2, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2,-2, 2, 2,-2,-2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2, 2,-2, 2,-2, 2,-2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1,-1,-1,-1, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1,-1, 1, 1,-1,-1, 1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1, 1,-1, 1,-1, 1,-1, 1/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)          
            vixyz(5)=ior(int3x(fix+vix(5)),ior(int3y(fiy+viy(1)),int3z(fiz+viz(1))))
            vixyz(6)=ior(int3x(fix+vix(6)),ior(int3y(fiy+viy(2)),int3z(fiz+viz(2))))
            vixyz(7)=ior(int3x(fix+vix(7)),ior(int3y(fiy+viy(3)),int3z(fiz+viz(3))))
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(4)),int3z(fiz+viz(4))))
          case(2)          
            vixyz(5)=ior(fint3x(fix+vix(5)),ior(fint3y(fiy+viy(1)),fint3z(fiz+viz(1))))
            vixyz(6)=ior(fint3x(fix+vix(6)),ior(fint3y(fiy+viy(2)),fint3z(fiz+viz(2))))
            vixyz(7)=ior(fint3x(fix+vix(7)),ior(fint3y(fiy+viy(3)),fint3z(fiz+viz(3))))
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(4)),fint3z(fiz+viz(4))))
          end select
          
          vcboxaddr( 16) = vixyz(5) + 48
          vcboxaddr( 17) = vixyz(5) + 56
          vcboxaddr( 18) = vixyz(6)
          vcboxaddr( 19) = vixyz(6) + 8
          vcboxaddr( 20) = vixyz(6) +16
          vcboxaddr( 21) = vixyz(6) +24
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7)+ 8
          vcboxaddr( 24) = vixyz(7) + 32
          vcboxaddr( 25) = vixyz(7) + 40
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + 8
          dist =136+child
          fwdstart = 17

          case(2)
          select case(int3exist)
          case(1)          
            vixyz(5)=ior(int3x(fix+vix5(5)),ior(int3y(fiy+viy5(1)),int3z(fiz+viz5(1))))
            vixyz(6)=ior(int3x(fix+vix5(6)),ior(int3y(fiy+viy5(2)),int3z(fiz+viz5(2))))
            vixyz(7)=ior(int3x(fix+vix5(7)),ior(int3y(fiy+viy5(3)),int3z(fiz+viz5(3))))
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(4)),int3z(fiz+viz5(4))))
          case(2)          
            vixyz(5)=ior(fint3x(fix+vix5(5)),ior(fint3y(fiy+viy5(1)),fint3z(fiz+viz5(1))))
            vixyz(6)=ior(fint3x(fix+vix5(6)),ior(fint3y(fiy+viy5(2)),fint3z(fiz+viz5(2))))
            vixyz(7)=ior(fint3x(fix+vix5(7)),ior(fint3y(fiy+viy5(3)),fint3z(fiz+viz5(3))))
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(4)),fint3z(fiz+viz5(4))))
          end select

          vcboxaddr( 16) = vixyz(5) + 6
          vcboxaddr( 17) = vixyz(5) + 7
          vcboxaddr( 18) = vixyz(6)
          vcboxaddr( 19) = vixyz(6) + 1
          vcboxaddr( 20) = vixyz(6) + 2
          vcboxaddr( 21) = vixyz(6) + 3
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + 1
          vcboxaddr( 24) = vixyz(7) + 4
          vcboxaddr( 25) = vixyz(7) + 5
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + 1
          
          fwdstart = 17	
	  end select
	  case(2)

	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	  
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(4)),int3z(fiz+pviz(4))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(4)),int3z(fiz+pviz(4))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(4)),fint3z(fiz+pviz(4))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(4)),fint3z(fiz+pviz(4))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(1) + poffset(2)
          vcboxaddr(  3) = vixyz(1) + poffset(4)
          vcboxaddr(  4) = vixyz(1) + poffset(6)
          vcboxaddr(  5) = vixyz(2)
          vcboxaddr(  6) = vixyz(2) + poffset(2)
          vcboxaddr(  7) = vixyz(3)
          vcboxaddr(  8) = vixyz(3) + poffset(4)
          vcboxaddr(  9) = vixyz(4)
          vcboxaddr( 10) = vixyz(5)
          vcboxaddr( 11) = vixyz(5) + poffset(1)
          vcboxaddr( 12) = vixyz(5) + poffset(2)
          vcboxaddr( 13) = vixyz(5) + poffset(3)
          vcboxaddr( 14) = vixyz(5) + poffset(4)
          vcboxaddr( 15) = vixyz(5) + poffset(5)
          vcboxaddr( 16) = vixyz(5) + poffset(6)
          vcboxaddr( 17) = vixyz(5) + poffset(7)
          vcboxaddr( 18) = vixyz(6)
          vcboxaddr( 19) = vixyz(6) + poffset(1)
          vcboxaddr( 20) = vixyz(6) + poffset(2)
          vcboxaddr( 21) = vixyz(6) + poffset(3)
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + poffset(1)
          vcboxaddr( 24) = vixyz(7) + poffset(4)
          vcboxaddr( 25) = vixyz(7) + poffset(5)
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + poffset(1)

          dist =136+child
          fwdstart = 17
      end select
      end subroutine caladdr18
      subroutine caladdr19(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       35275
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2, 2,-2, 2,-2, 2,-2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2,-2,-2,-2, 2, 2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2,-2, 2, 2,-2,-2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1, 1,-1, 1,-1, 1,-1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1,-1,-1,-1, 1, 1, 1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1,-1, 1, 1,-1,-1, 1, 1/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+vix(1)),ior(int3y(fiy+viy(1)),int3z(fiz+viz(1))))
            vixyz(2)=ior(int3x(fix+vix(2)),ior(int3y(fiy+viy(1)),int3z(fiz+viz(1))))
            vixyz(3)=ior(int3x(fix+vix(3)),ior(int3y(fiy+viy(3)),int3z(fiz+viz(3))))
            vixyz(4)=ior(int3x(fix+vix(4)),ior(int3y(fiy+viy(3)),int3z(fiz+viz(3))))
            vixyz(5)=ior(int3x(fix+vix(5)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(6)=ior(int3x(fix+vix(6)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(7)=ior(int3x(fix+vix(7)),ior(int3y(fiy+viy(7)),int3z(fiz+viz(7))))
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(7)),int3z(fiz+viz(7))))
          case(2)
            vixyz(1)=ior(fint3x(fix+vix(1)),ior(fint3y(fiy+viy(1)),fint3z(fiz+viz(1))))
            vixyz(2)=ior(fint3x(fix+vix(2)),ior(fint3y(fiy+viy(1)),fint3z(fiz+viz(1))))
            vixyz(3)=ior(fint3x(fix+vix(3)),ior(fint3y(fiy+viy(3)),fint3z(fiz+viz(3))))
            vixyz(4)=ior(fint3x(fix+vix(4)),ior(fint3y(fiy+viy(3)),fint3z(fiz+viz(3))))
            vixyz(5)=ior(fint3x(fix+vix(5)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(6)=ior(fint3x(fix+vix(6)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(7)=ior(fint3x(fix+vix(7)),ior(fint3y(fiy+viy(7)),fint3z(fiz+viz(7))))
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(7)),fint3z(fiz+viz(7))))
          end select
          
          vcboxaddr(  8) = vixyz(1) + 56
          vcboxaddr(  9) = vixyz(2)
          vcboxaddr( 10) = vixyz(2) +16
          vcboxaddr( 11) = vixyz(2) + 32
          vcboxaddr( 12) = vixyz(2) + 48
          vcboxaddr( 13) = vixyz(3)
          vcboxaddr( 14) = vixyz(3) + 8
          vcboxaddr( 15) = vixyz(3) +16
          vcboxaddr( 16) = vixyz(3) +24
          vcboxaddr( 17) = vixyz(4)
          vcboxaddr( 18) = vixyz(4) +16
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + 8
          vcboxaddr( 21) = vixyz(5) + 32
          vcboxaddr( 22) = vixyz(5) + 40
          vcboxaddr( 23) = vixyz(6)
          vcboxaddr( 24) = vixyz(6) + 32
          vcboxaddr( 25) = vixyz(7)
          vcboxaddr( 26) = vixyz(7) + 8
          vcboxaddr( 27) = vixyz(8)
          dist =144+child
          fwdstart =  9	
          case(2)
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+vix5(1)),ior(int3y(fiy+viy5(1)),int3z(fiz+viz5(1))))
            vixyz(2)=ior(int3x(fix+vix5(2)),ior(int3y(fiy+viy5(1)),int3z(fiz+viz5(1))))
            vixyz(3)=ior(int3x(fix+vix5(3)),ior(int3y(fiy+viy5(3)),int3z(fiz+viz5(3))))
            vixyz(4)=ior(int3x(fix+vix5(4)),ior(int3y(fiy+viy5(3)),int3z(fiz+viz5(3))))
            vixyz(5)=ior(int3x(fix+vix5(5)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(6)=ior(int3x(fix+vix5(6)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(7)=ior(int3x(fix+vix5(7)),ior(int3y(fiy+viy5(7)),int3z(fiz+viz5(7))))
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(7)),int3z(fiz+viz5(7))))
          case(2)
            vixyz(1)=ior(fint3x(fix+vix5(1)),ior(fint3y(fiy+viy5(1)),fint3z(fiz+viz5(1))))
            vixyz(2)=ior(fint3x(fix+vix5(2)),ior(fint3y(fiy+viy5(1)),fint3z(fiz+viz5(1))))
            vixyz(3)=ior(fint3x(fix+vix5(3)),ior(fint3y(fiy+viy5(3)),fint3z(fiz+viz5(3))))
            vixyz(4)=ior(fint3x(fix+vix5(4)),ior(fint3y(fiy+viy5(3)),fint3z(fiz+viz5(3))))
            vixyz(5)=ior(fint3x(fix+vix5(5)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(6)=ior(fint3x(fix+vix5(6)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(7)=ior(fint3x(fix+vix5(7)),ior(fint3y(fiy+viy5(7)),fint3z(fiz+viz5(7))))
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(7)),fint3z(fiz+viz5(7))))
          end select
          vcboxaddr(  8) = vixyz(1) + 7
          vcboxaddr(  9) = vixyz(2)
          vcboxaddr( 10) = vixyz(2) + 2
          vcboxaddr( 11) = vixyz(2) + 4
          vcboxaddr( 12) = vixyz(2) + 6
          vcboxaddr( 13) = vixyz(3)
          vcboxaddr( 14) = vixyz(3) + 1
          vcboxaddr( 15) = vixyz(3) + 2
          vcboxaddr( 16) = vixyz(3) + 3
          vcboxaddr( 17) = vixyz(4)
          vcboxaddr( 18) = vixyz(4) + 2
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + 1
          vcboxaddr( 21) = vixyz(5) + 4
          vcboxaddr( 22) = vixyz(5) + 5
          vcboxaddr( 23) = vixyz(6)
          vcboxaddr( 24) = vixyz(6) + 4
          vcboxaddr( 25) = vixyz(7)
          vcboxaddr( 26) = vixyz(7) + 1
          vcboxaddr( 27) = vixyz(8)
          
          fwdstart =  9
	  end select
	case(2)

	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(7)),int3z(fiz+pviz(7))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(7)),int3z(fiz+pviz(7))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(7)),fint3z(fiz+pviz(7))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(7)),fint3z(fiz+pviz(7))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(1) + poffset(1)
          vcboxaddr(  3) = vixyz(1) + poffset(2)
          vcboxaddr(  4) = vixyz(1) + poffset(3)
          vcboxaddr(  5) = vixyz(1) + poffset(4)
          vcboxaddr(  6) = vixyz(1) + poffset(5)
          vcboxaddr(  7) = vixyz(1) + poffset(6)
	  vcboxaddr(  8) = vixyz(1) + poffset(7)
          vcboxaddr(  9) = vixyz(2)
          vcboxaddr( 10) = vixyz(2) + poffset(2)
          vcboxaddr( 11) = vixyz(2) + poffset(4)
          vcboxaddr( 12) = vixyz(2) + poffset(6)
          vcboxaddr( 13) = vixyz(3)
          vcboxaddr( 14) = vixyz(3) + poffset(1)
          vcboxaddr( 15) = vixyz(3) + poffset(2)
          vcboxaddr( 16) = vixyz(3) + poffset(3)
          vcboxaddr( 17) = vixyz(4)
          vcboxaddr( 18) = vixyz(4) + poffset(2)
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + poffset(1)
          vcboxaddr( 21) = vixyz(5) + poffset(4)
          vcboxaddr( 22) = vixyz(5) + poffset(5)
          vcboxaddr( 23) = vixyz(6)
          vcboxaddr( 24) = vixyz(6) + poffset(4)
          vcboxaddr( 25) = vixyz(7)
          vcboxaddr( 26) = vixyz(7) + poffset(1)
          vcboxaddr( 27) = vixyz(8)
	  
          dist =144+child
          fwdstart =  9	
      end select
      end subroutine caladdr19
      subroutine caladdr20(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       24217
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2,-2, 2, 2,-2,-2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2, 0,-2, 0,-2, 0,-2, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2,-2,-2,-2, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1,-1, 1, 1,-1,-1, 1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1, 0,-1, 0,-1, 0,-1, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1,-1,-1,-1, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)          
            vixyz(6)=ior(int3x(fix+vix(6)),ior(int3y(fiy+viy(6)),int3z(fiz+viz(6))))
            vixyz(7)=ior(int3x(fix+vix(7)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(6)),int3z(fiz+viz(6))))
          case(2)          
            vixyz(6)=ior(fint3x(fix+vix(6)),ior(fint3y(fiy+viy(6)),fint3z(fiz+viz(6))))
            vixyz(7)=ior(fint3x(fix+vix(7)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(6)),fint3z(fiz+viz(6))))  
          end select        
          vcboxaddr( 15) = vixyz(6) + 8
          vcboxaddr( 16) = vixyz(6) + 16
          vcboxaddr( 17) = vixyz(6) + 24
          vcboxaddr( 18) = vixyz(6) + 32
          vcboxaddr( 19) = vixyz(6) + 40
          vcboxaddr( 20) = vixyz(6) + 48
          vcboxaddr( 21) = vixyz(6) + 56
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + 32
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + 16
          vcboxaddr( 26) = vixyz(8) + 32
          vcboxaddr( 27) = vixyz(8) + 48
          dist =152+child
          fwdstart =  16

          case(2)
          select case(int3exist)
          case(1)          
            vixyz(6)=ior(int3x(fix+vix5(6)),ior(int3y(fiy+viy5(6)),int3z(fiz+viz5(6))))
            vixyz(7)=ior(int3x(fix+vix5(7)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(6)),int3z(fiz+viz5(6))))
          case(2)          
            vixyz(6)=ior(fint3x(fix+vix5(6)),ior(fint3y(fiy+viy5(6)),fint3z(fiz+viz5(6))))
            vixyz(7)=ior(fint3x(fix+vix5(7)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(6)),fint3z(fiz+viz5(6))))   
          end select       
          vcboxaddr( 15) = vixyz(6) + 1
          vcboxaddr( 16) = vixyz(6) + 2
          vcboxaddr( 17) = vixyz(6) + 3
          vcboxaddr( 18) = vixyz(6) + 4
          vcboxaddr( 19) = vixyz(6) + 5
          vcboxaddr( 20) = vixyz(6) + 6
          vcboxaddr( 21) = vixyz(6) + 7
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + 4
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + 2
          vcboxaddr( 26) = vixyz(8) + 4
          vcboxaddr( 27) = vixyz(8) + 6
          
          fwdstart =  16
	  end select
	  case(2)

	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	  
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(6)),int3z(fiz+pviz(6))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(6)),int3z(fiz+pviz(6))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(6)),fint3z(fiz+pviz(6))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(6)),fint3z(fiz+pviz(6))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(1) + poffset(1)
          vcboxaddr(  3) = vixyz(2)
          vcboxaddr(  4) = vixyz(2) + poffset(1)
          vcboxaddr(  5) = vixyz(2) + poffset(2)
          vcboxaddr(  6) = vixyz(2) + poffset(3)
          vcboxaddr(  7) = vixyz(3)
          vcboxaddr(  8) = vixyz(4)
          vcboxaddr(  9) = vixyz(4) + poffset(2)
          vcboxaddr( 10) = vixyz(5)
          vcboxaddr( 11) = vixyz(5) + poffset(1)
          vcboxaddr( 12) = vixyz(5) + poffset(4)
          vcboxaddr( 13) = vixyz(5) + poffset(5)
          vcboxaddr( 14) = vixyz(6)
          vcboxaddr( 15) = vixyz(6) + poffset(1)
          vcboxaddr( 16) = vixyz(6) + poffset(2)
          vcboxaddr( 17) = vixyz(6) + poffset(3)
          vcboxaddr( 18) = vixyz(6) + poffset(4)
          vcboxaddr( 19) = vixyz(6) + poffset(5)
          vcboxaddr( 20) = vixyz(6) + poffset(6)
          vcboxaddr( 21) = vixyz(6) + poffset(7)
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + poffset(4)
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + poffset(2)
          vcboxaddr( 26) = vixyz(8) + poffset(4)
          vcboxaddr( 27) = vixyz(8) + poffset(6)
	  
          dist =152+child
          fwdstart =  16
      end select
      end subroutine caladdr20
      subroutine caladdr21(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       25765
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2, 0,-2, 0,-2, 0,-2, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2,-2, 0, 0,-2,-2, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2,-2,-2,-2, 2, 2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1, 0,-1, 0,-1, 0,-1, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1,-1, 0, 0,-1,-1, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1,-1,-1,-1, 1, 1, 1, 1/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)          
            vixyz(4)=ior(int3x(fix+vix(4)),ior(int3y(fiy+viy(3)),int3z(fiz+viz(3))))
            vixyz(5)=ior(int3x(fix+vix(5)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(6)=ior(int3x(fix+vix(6)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(7)=ior(int3x(fix+vix(7)),ior(int3y(fiy+viy(7)),int3z(fiz+viz(7))))
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(7)),int3z(fiz+viz(7))))
          case(2)          
            vixyz(4)=ior(fint3x(fix+vix(4)),ior(fint3y(fiy+viy(3)),fint3z(fiz+viz(3))))
            vixyz(5)=ior(fint3x(fix+vix(5)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(6)=ior(fint3x(fix+vix(6)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(7)=ior(fint3x(fix+vix(7)),ior(fint3y(fiy+viy(7)),fint3z(fiz+viz(7))))
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(7)),fint3z(fiz+viz(7))))   
          end select       
          vcboxaddr( 15) = vixyz(4) + 32
          vcboxaddr( 16) = vixyz(4) + 40
          vcboxaddr( 17) = vixyz(4) + 48
          vcboxaddr( 18) = vixyz(4) + 56
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(6)
          vcboxaddr( 21) = vixyz(6) + 8
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) +16
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + 8
          vcboxaddr( 26) = vixyz(8) + 16
          vcboxaddr( 27) = vixyz(8) + 24
          dist =160+child
          fwdstart =  16	
          case(2)   
          select case(int3exist)
          case(1)       
            vixyz(4)=ior(int3x(fix+vix5(4)),ior(int3y(fiy+viy5(3)),int3z(fiz+viz5(3))))
            vixyz(5)=ior(int3x(fix+vix5(5)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(6)=ior(int3x(fix+vix5(6)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(7)=ior(int3x(fix+vix5(7)),ior(int3y(fiy+viy5(7)),int3z(fiz+viz5(7))))
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(7)),int3z(fiz+viz5(7))))
          case(2)       
            vixyz(4)=ior(fint3x(fix+vix5(4)),ior(fint3y(fiy+viy5(3)),fint3z(fiz+viz5(3))))
            vixyz(5)=ior(fint3x(fix+vix5(5)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(6)=ior(fint3x(fix+vix5(6)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(7)=ior(fint3x(fix+vix5(7)),ior(fint3y(fiy+viy5(7)),fint3z(fiz+viz5(7))))
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(7)),fint3z(fiz+viz5(7))))   
          end select       
          vcboxaddr( 15) = vixyz(4) + 4
          vcboxaddr( 16) = vixyz(4) + 5
          vcboxaddr( 17) = vixyz(4) + 6
          vcboxaddr( 18) = vixyz(4) + 7
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(6)
          vcboxaddr( 21) = vixyz(6) + 1
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + 2
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + 1
          vcboxaddr( 26) = vixyz(8) + 2
          vcboxaddr( 27) = vixyz(8) + 3
	  
          fwdstart =  16	
        end select
	case(2)

	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(7)),int3z(fiz+pviz(7))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(7)),int3z(fiz+pviz(7))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(7)),fint3z(fiz+pviz(7))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(7)),fint3z(fiz+pviz(7))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(1) + poffset(4)
          vcboxaddr(  3) = vixyz(2)
          vcboxaddr(  4) = vixyz(2) + poffset(1)
          vcboxaddr(  5) = vixyz(2) + poffset(4)
          vcboxaddr(  6) = vixyz(2) + poffset(5)
          vcboxaddr(  7) = vixyz(3)
          vcboxaddr(  8) = vixyz(3) + poffset(2)
          vcboxaddr(  9) = vixyz(3) + poffset(4)
          vcboxaddr( 10) = vixyz(3) + poffset(6)
          vcboxaddr( 11) = vixyz(4)
          vcboxaddr( 12) = vixyz(4) + poffset(1)
          vcboxaddr( 13) = vixyz(4) + poffset(2)
          vcboxaddr( 14) = vixyz(4) + poffset(3)
          vcboxaddr( 15) = vixyz(4) + poffset(4)
          vcboxaddr( 16) = vixyz(4) + poffset(5)
          vcboxaddr( 17) = vixyz(4) + poffset(6)
          vcboxaddr( 18) = vixyz(4) + poffset(7)
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(6)
          vcboxaddr( 21) = vixyz(6) + poffset(1)
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + poffset(2)
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + poffset(1)
          vcboxaddr( 26) = vixyz(8) + poffset(2)
          vcboxaddr( 27) = vixyz(8) + poffset(3)
          dist =160+child
          fwdstart =  16
      end select
      end subroutine caladdr21
      subroutine caladdr22(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       33163
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2,-2,-2,-2, 2, 2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2, 0,-2, 0,-2, 0,-2, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2,-2, 2, 2,-2,-2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1,-1,-1,-1, 1, 1, 1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1, 0,-1, 0,-1, 0,-1, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1,-1, 1, 1,-1,-1, 1, 1/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)         
            vixyz(2)=ior(int3x(fix+vix(2)),ior(int3y(fiy+viy(2)),int3z(fiz+viz(2))))
            vixyz(3)=ior(int3x(fix+vix(3)),ior(int3y(fiy+viy(3)),int3z(fiz+viz(3))))
            vixyz(4)=ior(int3x(fix+vix(4)),ior(int3y(fiy+viy(4)),int3z(fiz+viz(4))))
            vixyz(5)=ior(int3x(fix+vix(5)),ior(int3y(fiy+viy(1)),int3z(fiz+viz(1))))
            vixyz(6)=ior(int3x(fix+vix(6)),ior(int3y(fiy+viy(2)),int3z(fiz+viz(2))))
            vixyz(7)=ior(int3x(fix+vix(7)),ior(int3y(fiy+viy(3)),int3z(fiz+viz(3))))
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(4)),int3z(fiz+viz(4))))
          case(2)         
            vixyz(2)=ior(fint3x(fix+vix(2)),ior(fint3y(fiy+viy(2)),fint3z(fiz+viz(2))))
            vixyz(3)=ior(fint3x(fix+vix(3)),ior(fint3y(fiy+viy(3)),fint3z(fiz+viz(3))))
            vixyz(4)=ior(fint3x(fix+vix(4)),ior(fint3y(fiy+viy(4)),fint3z(fiz+viz(4))))
            vixyz(5)=ior(fint3x(fix+vix(5)),ior(fint3y(fiy+viy(1)),fint3z(fiz+viz(1))))
            vixyz(6)=ior(fint3x(fix+vix(6)),ior(fint3y(fiy+viy(2)),fint3z(fiz+viz(2))))
            vixyz(7)=ior(fint3x(fix+vix(7)),ior(fint3y(fiy+viy(3)),fint3z(fiz+viz(3))))
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(4)),fint3z(fiz+viz(4))))  
          end select        
          vcboxaddr( 10) = vixyz(2) + 40
          vcboxaddr( 11) = vixyz(2) + 48
          vcboxaddr( 12) = vixyz(2) + 56
          vcboxaddr( 13) = vixyz(3)
          vcboxaddr( 14) = vixyz(3) + 8
          vcboxaddr( 15) = vixyz(4)
          vcboxaddr( 16) = vixyz(4) + 8
          vcboxaddr( 17) = vixyz(4) +16
          vcboxaddr( 18) = vixyz(4) +24
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + 32
          vcboxaddr( 21) = vixyz(6)
          vcboxaddr( 22) = vixyz(6) +16
          vcboxaddr( 23) = vixyz(6) + 32
          vcboxaddr( 24) = vixyz(6) + 48
          vcboxaddr( 25) = vixyz(7)
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) +16
          dist =168+child
          fwdstart = 11
          case(2)
          select case(int3exist)
          case(1)          
            vixyz(2)=ior(int3x(fix+vix5(2)),ior(int3y(fiy+viy5(2)),int3z(fiz+viz5(2))))
            vixyz(3)=ior(int3x(fix+vix5(3)),ior(int3y(fiy+viy5(3)),int3z(fiz+viz5(3))))
            vixyz(4)=ior(int3x(fix+vix5(4)),ior(int3y(fiy+viy5(4)),int3z(fiz+viz5(4))))
            vixyz(5)=ior(int3x(fix+vix5(5)),ior(int3y(fiy+viy5(1)),int3z(fiz+viz5(1))))
            vixyz(6)=ior(int3x(fix+vix5(6)),ior(int3y(fiy+viy5(2)),int3z(fiz+viz5(2))))
            vixyz(7)=ior(int3x(fix+vix5(7)),ior(int3y(fiy+viy5(3)),int3z(fiz+viz5(3))))
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(4)),int3z(fiz+viz5(4))))
          case(2)          
            vixyz(2)=ior(fint3x(fix+vix5(2)),ior(fint3y(fiy+viy5(2)),fint3z(fiz+viz5(2))))
            vixyz(3)=ior(fint3x(fix+vix5(3)),ior(fint3y(fiy+viy5(3)),fint3z(fiz+viz5(3))))
            vixyz(4)=ior(fint3x(fix+vix5(4)),ior(fint3y(fiy+viy5(4)),fint3z(fiz+viz5(4))))
            vixyz(5)=ior(fint3x(fix+vix5(5)),ior(fint3y(fiy+viy5(1)),fint3z(fiz+viz5(1))))
            vixyz(6)=ior(fint3x(fix+vix5(6)),ior(fint3y(fiy+viy5(2)),fint3z(fiz+viz5(2))))
            vixyz(7)=ior(fint3x(fix+vix5(7)),ior(fint3y(fiy+viy5(3)),fint3z(fiz+viz5(3))))
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(4)),fint3z(fiz+viz5(4))))
          end select
          vcboxaddr( 10) = vixyz(2) + 5
          vcboxaddr( 11) = vixyz(2) + 6
          vcboxaddr( 12) = vixyz(2) + 7
          vcboxaddr( 13) = vixyz(3)
          vcboxaddr( 14) = vixyz(3) + 1
          vcboxaddr( 15) = vixyz(4)
          vcboxaddr( 16) = vixyz(4) + 1
          vcboxaddr( 17) = vixyz(4) + 2
          vcboxaddr( 18) = vixyz(4) + 3
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + 4
          vcboxaddr( 21) = vixyz(6)
          vcboxaddr( 22) = vixyz(6) + 2
          vcboxaddr( 23) = vixyz(6) + 4
          vcboxaddr( 24) = vixyz(6) + 6
          vcboxaddr( 25) = vixyz(7)
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + 2
                
          fwdstart = 11
	  end select
	case(2)

	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(4)),int3z(fiz+pviz(4))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(4)),int3z(fiz+pviz(4))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(4)),fint3z(fiz+pviz(4))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(4)),fint3z(fiz+pviz(4))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(1) + poffset(1)
          vcboxaddr(  3) = vixyz(1) + poffset(4)
          vcboxaddr(  4) = vixyz(1) + poffset(5)
          vcboxaddr(  5) = vixyz(2)
          vcboxaddr(  6) = vixyz(2) + poffset(1)
          vcboxaddr(  7) = vixyz(2) + poffset(2)
          vcboxaddr(  8) = vixyz(2) + poffset(3)
          vcboxaddr(  9) = vixyz(2) + poffset(4)
	  vcboxaddr( 10) = vixyz(2) + poffset(5)
          vcboxaddr( 11) = vixyz(2) + poffset(6)
          vcboxaddr( 12) = vixyz(2) + poffset(7)
          vcboxaddr( 13) = vixyz(3)
          vcboxaddr( 14) = vixyz(3) + poffset(1)
          vcboxaddr( 15) = vixyz(4)
          vcboxaddr( 16) = vixyz(4) + poffset(1)
          vcboxaddr( 17) = vixyz(4) + poffset(2)
          vcboxaddr( 18) = vixyz(4) + poffset(3)
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + poffset(4)
          vcboxaddr( 21) = vixyz(6)
          vcboxaddr( 22) = vixyz(6) + poffset(2)
          vcboxaddr( 23) = vixyz(6) + poffset(4)
          vcboxaddr( 24) = vixyz(6) + poffset(6)
          vcboxaddr( 25) = vixyz(7)
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + poffset(2)

          dist =168+child
          fwdstart = 11
      end select
      end subroutine caladdr22
      subroutine caladdr23(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       32110
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2, 0,-2, 0,-2, 0,-2, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2,-2,-2,-2, 2, 2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2,-2, 2, 2,-2,-2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1, 0,-1, 0,-1, 0,-1, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1,-1,-1,-1, 1, 1, 1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1,-1, 1, 1,-1,-1, 1, 1/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)          
            vixyz(2)=ior(int3x(fix+vix(2)),ior(int3y(fiy+viy(1)),int3z(fiz+viz(1))))
            vixyz(3)=ior(int3x(fix+vix(3)),ior(int3y(fiy+viy(3)),int3z(fiz+viz(3))))
            vixyz(4)=ior(int3x(fix+vix(4)),ior(int3y(fiy+viy(3)),int3z(fiz+viz(3))))
            vixyz(5)=ior(int3x(fix+vix(5)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(6)=ior(int3x(fix+vix(6)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(7)=ior(int3x(fix+vix(7)),ior(int3y(fiy+viy(7)),int3z(fiz+viz(7))))
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(7)),int3z(fiz+viz(7))))
          case(2)          
            vixyz(2)=ior(fint3x(fix+vix(2)),ior(fint3y(fiy+viy(1)),fint3z(fiz+viz(1))))
            vixyz(3)=ior(fint3x(fix+vix(3)),ior(fint3y(fiy+viy(3)),fint3z(fiz+viz(3))))
            vixyz(4)=ior(fint3x(fix+vix(4)),ior(fint3y(fiy+viy(3)),fint3z(fiz+viz(3))))
            vixyz(5)=ior(fint3x(fix+vix(5)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(6)=ior(fint3x(fix+vix(6)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(7)=ior(fint3x(fix+vix(7)),ior(fint3y(fiy+viy(7)),fint3z(fiz+viz(7))))
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(7)),fint3z(fiz+viz(7))))      
          end select    
          vcboxaddr( 11) = vixyz(2) + 48
          vcboxaddr( 12) = vixyz(2) + 56
          vcboxaddr( 13) = vixyz(3)
          vcboxaddr( 14) = vixyz(3) +16
          vcboxaddr( 15) = vixyz(4)
          vcboxaddr( 16) = vixyz(4) + 8
          vcboxaddr( 17) = vixyz(4) +16
          vcboxaddr( 18) = vixyz(4) +24
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + 32
          vcboxaddr( 21) = vixyz(6)
          vcboxaddr( 22) = vixyz(6) + 8
          vcboxaddr( 23) = vixyz(6) + 32
          vcboxaddr( 24) = vixyz(6) + 40
          vcboxaddr( 25) = vixyz(7)
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + 8
          dist =176+child
          fwdstart = 12
          case(2)
          select case(int3exist)
          case(1)          
            vixyz(2)=ior(int3x(fix+vix5(2)),ior(int3y(fiy+viy5(1)),int3z(fiz+viz5(1))))
            vixyz(3)=ior(int3x(fix+vix5(3)),ior(int3y(fiy+viy5(3)),int3z(fiz+viz5(3))))
            vixyz(4)=ior(int3x(fix+vix5(4)),ior(int3y(fiy+viy5(3)),int3z(fiz+viz5(3))))
            vixyz(5)=ior(int3x(fix+vix5(5)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(6)=ior(int3x(fix+vix5(6)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(7)=ior(int3x(fix+vix5(7)),ior(int3y(fiy+viy5(7)),int3z(fiz+viz5(7))))
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(7)),int3z(fiz+viz5(7))))
          case(2)          
            vixyz(2)=ior(fint3x(fix+vix5(2)),ior(fint3y(fiy+viy5(1)),fint3z(fiz+viz5(1))))
            vixyz(3)=ior(fint3x(fix+vix5(3)),ior(fint3y(fiy+viy5(3)),fint3z(fiz+viz5(3))))
            vixyz(4)=ior(fint3x(fix+vix5(4)),ior(fint3y(fiy+viy5(3)),fint3z(fiz+viz5(3))))
            vixyz(5)=ior(fint3x(fix+vix5(5)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(6)=ior(fint3x(fix+vix5(6)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(7)=ior(fint3x(fix+vix5(7)),ior(fint3y(fiy+viy5(7)),fint3z(fiz+viz5(7))))
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(7)),fint3z(fiz+viz5(7))))
          end select
          vcboxaddr( 11) = vixyz(2) + 6
          vcboxaddr( 12) = vixyz(2) + 7
          vcboxaddr( 13) = vixyz(3)
          vcboxaddr( 14) = vixyz(3) + 2
          vcboxaddr( 15) = vixyz(4)
          vcboxaddr( 16) = vixyz(4) + 1
          vcboxaddr( 17) = vixyz(4) + 2
          vcboxaddr( 18) = vixyz(4) + 3
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + 4
          vcboxaddr( 21) = vixyz(6)
          vcboxaddr( 22) = vixyz(6) + 1
          vcboxaddr( 23) = vixyz(6) + 4
          vcboxaddr( 24) = vixyz(6) + 5
          vcboxaddr( 25) = vixyz(7)
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + 1
          
           fwdstart = 12
	  end select
	case(2)

	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(7)),int3z(fiz+pviz(7))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(7)),int3z(fiz+pviz(7))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(7)),fint3z(fiz+pviz(7))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(7)),fint3z(fiz+pviz(7))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(1) + poffset(2)
          vcboxaddr(  3) = vixyz(1) + poffset(4)
          vcboxaddr(  4) = vixyz(1) + poffset(6)
          vcboxaddr(  5) = vixyz(2)
          vcboxaddr(  6) = vixyz(2) + poffset(1)
          vcboxaddr(  7) = vixyz(2) + poffset(2)
          vcboxaddr(  8) = vixyz(2) + poffset(3)
          vcboxaddr(  9) = vixyz(2) + poffset(4)
          vcboxaddr( 10) = vixyz(2) + poffset(5)
          vcboxaddr( 11) = vixyz(2) + poffset(6)
          vcboxaddr( 12) = vixyz(2) + poffset(7)
          vcboxaddr( 13) = vixyz(3)
          vcboxaddr( 14) = vixyz(3) + poffset(2)
          vcboxaddr( 15) = vixyz(4)
          vcboxaddr( 16) = vixyz(4) + poffset(1)
          vcboxaddr( 17) = vixyz(4) + poffset(2)
          vcboxaddr( 18) = vixyz(4) + poffset(3)
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + poffset(4)
          vcboxaddr( 21) = vixyz(6)
          vcboxaddr( 22) = vixyz(6) + poffset(1)
          vcboxaddr( 23) = vixyz(6) + poffset(4)
          vcboxaddr( 24) = vixyz(6) + poffset(5)
          vcboxaddr( 25) = vixyz(7)
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + poffset(1)

          dist =176+child
           fwdstart = 12
      end select
      end subroutine caladdr23
      subroutine caladdr24(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       35887
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2,-2, 2, 2,-2,-2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2,-2,-2,-2, 2, 2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2, 2,-2, 2,-2, 2,-2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1,-1, 1, 1,-1,-1, 1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1,-1,-1,-1, 1, 1, 1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1, 1,-1, 1,-1, 1,-1, 1/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+vix(1)),ior(int3y(fiy+viy(1)),int3z(fiz+viz(1))))
            vixyz(2)=ior(int3x(fix+vix(2)),ior(int3y(fiy+viy(2)),int3z(fiz+viz(2))))
            vixyz(3)=ior(int3x(fix+vix(3)),ior(int3y(fiy+viy(1)),int3z(fiz+viz(1))))
            vixyz(4)=ior(int3x(fix+vix(4)),ior(int3y(fiy+viy(2)),int3z(fiz+viz(2))))
            vixyz(5)=ior(int3x(fix+vix(5)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(6)=ior(int3x(fix+vix(6)),ior(int3y(fiy+viy(6)),int3z(fiz+viz(6))))
            vixyz(7)=ior(int3x(fix+vix(7)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(6)),int3z(fiz+viz(6))))
          case(2)
            vixyz(1)=ior(fint3x(fix+vix(1)),ior(fint3y(fiy+viy(1)),fint3z(fiz+viz(1))))
            vixyz(2)=ior(fint3x(fix+vix(2)),ior(fint3y(fiy+viy(2)),fint3z(fiz+viz(2))))
            vixyz(3)=ior(fint3x(fix+vix(3)),ior(fint3y(fiy+viy(1)),fint3z(fiz+viz(1))))
            vixyz(4)=ior(fint3x(fix+vix(4)),ior(fint3y(fiy+viy(2)),fint3z(fiz+viz(2))))
            vixyz(5)=ior(fint3x(fix+vix(5)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(6)=ior(fint3x(fix+vix(6)),ior(fint3y(fiy+viy(6)),fint3z(fiz+viz(6))))
            vixyz(7)=ior(fint3x(fix+vix(7)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(6)),fint3z(fiz+viz(6))))    
          end select      
          vcboxaddr(  8) = vixyz(1) + 56
          vcboxaddr(  9) = vixyz(2)
          vcboxaddr( 10) = vixyz(2) + 8
          vcboxaddr( 11) = vixyz(2) +16
          vcboxaddr( 12) = vixyz(2) +24
          vcboxaddr( 13) = vixyz(3)
          vcboxaddr( 14) = vixyz(3) +16
          vcboxaddr( 15) = vixyz(3) + 32
          vcboxaddr( 16) = vixyz(3) + 48
          vcboxaddr( 17) = vixyz(4)
          vcboxaddr( 18) = vixyz(4) +16
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + 8
          vcboxaddr( 21) = vixyz(5) + 32
          vcboxaddr( 22) = vixyz(5) + 40
          vcboxaddr( 23) = vixyz(6)
          vcboxaddr( 24) = vixyz(6) + 8
          vcboxaddr( 25) = vixyz(7)
          vcboxaddr( 26) = vixyz(7) + 32
          vcboxaddr( 27) = vixyz(8)
          dist =184+child
          fwdstart = 9
          case(2)
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+vix5(1)),ior(int3y(fiy+viy5(1)),int3z(fiz+viz5(1))))
            vixyz(2)=ior(int3x(fix+vix5(2)),ior(int3y(fiy+viy5(2)),int3z(fiz+viz5(2))))
            vixyz(3)=ior(int3x(fix+vix5(3)),ior(int3y(fiy+viy5(1)),int3z(fiz+viz5(1))))
            vixyz(4)=ior(int3x(fix+vix5(4)),ior(int3y(fiy+viy5(2)),int3z(fiz+viz5(2))))
            vixyz(5)=ior(int3x(fix+vix5(5)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(6)=ior(int3x(fix+vix5(6)),ior(int3y(fiy+viy5(6)),int3z(fiz+viz5(6))))
            vixyz(7)=ior(int3x(fix+vix5(7)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(6)),int3z(fiz+viz5(6))))
          case(2)
            vixyz(1)=ior(fint3x(fix+vix5(1)),ior(fint3y(fiy+viy5(1)),fint3z(fiz+viz5(1))))
            vixyz(2)=ior(fint3x(fix+vix5(2)),ior(fint3y(fiy+viy5(2)),fint3z(fiz+viz5(2))))
            vixyz(3)=ior(fint3x(fix+vix5(3)),ior(fint3y(fiy+viy5(1)),fint3z(fiz+viz5(1))))
            vixyz(4)=ior(fint3x(fix+vix5(4)),ior(fint3y(fiy+viy5(2)),fint3z(fiz+viz5(2))))
            vixyz(5)=ior(fint3x(fix+vix5(5)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(6)=ior(fint3x(fix+vix5(6)),ior(fint3y(fiy+viy5(6)),fint3z(fiz+viz5(6))))
            vixyz(7)=ior(fint3x(fix+vix5(7)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(6)),fint3z(fiz+viz5(6))))
          end select
          vcboxaddr(  8) = vixyz(1) + 7
          vcboxaddr(  9) = vixyz(2)
          vcboxaddr( 10) = vixyz(2) + 1
          vcboxaddr( 11) = vixyz(2) + 2
          vcboxaddr( 12) = vixyz(2) + 3
          vcboxaddr( 13) = vixyz(3)
          vcboxaddr( 14) = vixyz(3) + 2
          vcboxaddr( 15) = vixyz(3) + 4
          vcboxaddr( 16) = vixyz(3) + 6
          vcboxaddr( 17) = vixyz(4)
          vcboxaddr( 18) = vixyz(4) + 2
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + 1
          vcboxaddr( 21) = vixyz(5) + 4
          vcboxaddr( 22) = vixyz(5) + 5
          vcboxaddr( 23) = vixyz(6)
          vcboxaddr( 24) = vixyz(6) + 1
          vcboxaddr( 25) = vixyz(7)
          vcboxaddr( 26) = vixyz(7) + 4
          vcboxaddr( 27) = vixyz(8)
          
           fwdstart = 9
        end select	  
        case(2)

	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(6)),int3z(fiz+pviz(6))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(6)),int3z(fiz+pviz(6))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(6)),fint3z(fiz+pviz(6))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(6)),fint3z(fiz+pviz(6))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(1) + poffset(1)
          vcboxaddr(  3) = vixyz(1) + poffset(2)
          vcboxaddr(  4) = vixyz(1) + poffset(3)
          vcboxaddr(  5) = vixyz(1) + poffset(4)
          vcboxaddr(  6) = vixyz(1) + poffset(5)
          vcboxaddr(  7) = vixyz(1) + poffset(6)
          vcboxaddr(  8) = vixyz(1) + poffset(7)
          vcboxaddr(  9) = vixyz(2)
          vcboxaddr( 10) = vixyz(2) + poffset(1)
          vcboxaddr( 11) = vixyz(2) + poffset(2)
          vcboxaddr( 12) = vixyz(2) + poffset(3)
          vcboxaddr( 13) = vixyz(3)
          vcboxaddr( 14) = vixyz(3) + poffset(2)
          vcboxaddr( 15) = vixyz(3) + poffset(4)
          vcboxaddr( 16) = vixyz(3) + poffset(6)
          vcboxaddr( 17) = vixyz(4)
          vcboxaddr( 18) = vixyz(4) + poffset(2)
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + poffset(1)
          vcboxaddr( 21) = vixyz(5) + poffset(4)
          vcboxaddr( 22) = vixyz(5) + poffset(5)
          vcboxaddr( 23) = vixyz(6)
          vcboxaddr( 24) = vixyz(6) + poffset(1)
          vcboxaddr( 25) = vixyz(7)
          vcboxaddr( 26) = vixyz(7) + poffset(4)
          vcboxaddr( 27) = vixyz(8)
	  
          dist =184+child
          fwdstart = 9
      end select
      end subroutine caladdr24
      subroutine caladdr25(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       16399
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2,-2, 0, 0,-2,-2, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2,-2,-2,-2, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2, 0,-2, 0,-2, 0,-2, 0/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1,-1, 0, 0,-1,-1, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1,-1,-1,-1, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1, 0,-1, 0,-1, 0,-1, 0/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)                    
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(6)),int3z(fiz+viz(6))))
          case(2)                    
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(6)),fint3z(fiz+viz(6))))      
          end select    
          vcboxaddr( 20) = vixyz(8)
          vcboxaddr( 21) = vixyz(8) + 8
          vcboxaddr( 22) = vixyz(8) + 16
          vcboxaddr( 23) = vixyz(8) + 24
          vcboxaddr( 24) = vixyz(8) + 32
          vcboxaddr( 25) = vixyz(8) + 40
          vcboxaddr( 26) = vixyz(8) + 48
          vcboxaddr( 27) = vixyz(8) + 56
          dist =192+child
          fwdstart = 21	
          case(2)
          select case(int3exist)
          case(1)          
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(6)),int3z(fiz+viz5(6))))
          case(2)          
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(6)),fint3z(fiz+viz5(6))))     
          end select     
          vcboxaddr( 20) = vixyz(8)
          vcboxaddr( 21) = vixyz(8) + 1
          vcboxaddr( 22) = vixyz(8) + 2
          vcboxaddr( 23) = vixyz(8) + 3
          vcboxaddr( 24) = vixyz(8) + 4
          vcboxaddr( 25) = vixyz(8) + 5
          vcboxaddr( 26) = vixyz(8) + 6
          vcboxaddr( 27) = vixyz(8) + 7
          
          fwdstart = 21
        end select
        case(2)

	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(6)),int3z(fiz+pviz(6))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(6)),int3z(fiz+pviz(6))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(6)),fint3z(fiz+pviz(6))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(6)),fint3z(fiz+pviz(6))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(2)
          vcboxaddr(  3) = vixyz(2) + poffset(4)
          vcboxaddr(  4) = vixyz(3)
          vcboxaddr(  5) = vixyz(3) + poffset(1)
          vcboxaddr(  6) = vixyz(4)
          vcboxaddr(  7) = vixyz(4) + poffset(1)
          vcboxaddr(  8) = vixyz(4) + poffset(4)
          vcboxaddr(  9) = vixyz(4) + poffset(5)
          vcboxaddr( 10) = vixyz(5)
          vcboxaddr( 11) = vixyz(5) + poffset(2)
          vcboxaddr( 12) = vixyz(6)
          vcboxaddr( 13) = vixyz(6) + poffset(2)
          vcboxaddr( 14) = vixyz(6) + poffset(4)
          vcboxaddr( 15) = vixyz(6) + poffset(6)
          vcboxaddr( 16) = vixyz(7)
          vcboxaddr( 17) = vixyz(7) + poffset(1)
          vcboxaddr( 18) = vixyz(7) + poffset(2)
          vcboxaddr( 19) = vixyz(7) + poffset(3)
          vcboxaddr( 20) = vixyz(8)
          vcboxaddr( 21) = vixyz(8) + poffset(1)
          vcboxaddr( 22) = vixyz(8) + poffset(2)
          vcboxaddr( 23) = vixyz(8) + poffset(3)
          vcboxaddr( 24) = vixyz(8) + poffset(4)
          vcboxaddr( 25) = vixyz(8) + poffset(5)
          vcboxaddr( 26) = vixyz(8) + poffset(6)
          vcboxaddr( 27) = vixyz(8) + poffset(7)
	  
          dist =192+child
          fwdstart = 21
      end select
      end subroutine caladdr25
      subroutine caladdr26(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       20452
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2, 2,-2, 2,-2, 2,-2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2,-2,-2,-2, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2,-2, 0, 0,-2,-2, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1, 1,-1, 1,-1, 1,-1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1,-1,-1,-1, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1,-1, 0, 0,-1,-1, 0, 0/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)         
            vixyz(7)=ior(int3x(fix+vix(7)),ior(int3y(fiy+viy(7)),int3z(fiz+viz(7))))
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(7)),int3z(fiz+viz(7))))
          case(2)         
            vixyz(7)=ior(fint3x(fix+vix(7)),ior(fint3y(fiy+viy(7)),fint3z(fiz+viz(7))))
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(7)),fint3z(fiz+viz(7))))   
          end select       
          vcboxaddr( 17) = vixyz(7) + 8
          vcboxaddr( 18) = vixyz(7) +16
          vcboxaddr( 19) = vixyz(7) +24
          vcboxaddr( 20) = vixyz(7) + 32
          vcboxaddr( 21) = vixyz(7) + 40
          vcboxaddr( 22) = vixyz(7) + 48
          vcboxaddr( 23) = vixyz(7) + 56
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) +16
          vcboxaddr( 26) = vixyz(8) + 32
          vcboxaddr( 27) = vixyz(8) + 48
          dist =200+child
          fwdstart = 18	
         case(2)
          select case(int3exist)
          case(1)          
            vixyz(7)=ior(int3x(fix+vix5(7)),ior(int3y(fiy+viy5(7)),int3z(fiz+viz5(7))))
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(7)),int3z(fiz+viz5(7))))
          case(2)          
            vixyz(7)=ior(fint3x(fix+vix5(7)),ior(fint3y(fiy+viy5(7)),fint3z(fiz+viz5(7))))
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(7)),fint3z(fiz+viz5(7))))
          end select
          vcboxaddr( 17) = vixyz(7) + 1
          vcboxaddr( 18) = vixyz(7) + 2
          vcboxaddr( 19) = vixyz(7) + 3
          vcboxaddr( 20) = vixyz(7) + 4
          vcboxaddr( 21) = vixyz(7) + 5
          vcboxaddr( 22) = vixyz(7) + 6
          vcboxaddr( 23) = vixyz(7) + 7
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + 2
          vcboxaddr( 26) = vixyz(8) + 4
          vcboxaddr( 27) = vixyz(8) + 6
          
          fwdstart = 18
        end select
        case(2)
	
	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(7)),int3z(fiz+pviz(7))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(7)),int3z(fiz+pviz(7))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(7)),fint3z(fiz+pviz(7))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(7)),fint3z(fiz+pviz(7))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(1) + poffset(1)
          vcboxaddr(  3) = vixyz(2)
          vcboxaddr(  4) = vixyz(3)
          vcboxaddr(  5) = vixyz(3) + poffset(1)
          vcboxaddr(  6) = vixyz(3) + poffset(4)
          vcboxaddr(  7) = vixyz(3) + poffset(5)
          vcboxaddr(  8) = vixyz(4)
          vcboxaddr(  9) = vixyz(4) + poffset(4)
          vcboxaddr( 10) = vixyz(5)
          vcboxaddr( 11) = vixyz(5) + poffset(1)
          vcboxaddr( 12) = vixyz(5) + poffset(2)
          vcboxaddr( 13) = vixyz(5) + poffset(3)
          vcboxaddr( 14) = vixyz(6)
          vcboxaddr( 15) = vixyz(6) + poffset(2)
          vcboxaddr( 16) = vixyz(7)
          vcboxaddr( 17) = vixyz(7) + poffset(1)
          vcboxaddr( 18) = vixyz(7) + poffset(2)
          vcboxaddr( 19) = vixyz(7) + poffset(3)
          vcboxaddr( 20) = vixyz(7) + poffset(4)
          vcboxaddr( 21) = vixyz(7) + poffset(5)
          vcboxaddr( 22) = vixyz(7) + poffset(6)
          vcboxaddr( 23) = vixyz(7) + poffset(7)
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + poffset(2)
          vcboxaddr( 26) = vixyz(8) + poffset(4)
          vcboxaddr( 27) = vixyz(8) + poffset(6)

          dist =200+child
          fwdstart = 18	
      end select
      end subroutine caladdr26
      subroutine caladdr27(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       19591
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2,-2,-2,-2, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2, 2,-2, 2,-2, 2,-2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2,-2, 0, 0,-2,-2, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1,-1,-1,-1, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1, 1,-1, 1,-1, 1,-1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1,-1, 0, 0,-1,-1, 0, 0/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)          
            vixyz(7)=ior(int3x(fix+vix(7)),ior(int3y(fiy+viy(3)),int3z(fiz+viz(3))))
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(4)),int3z(fiz+viz(4))))
          case(2)          
            vixyz(7)=ior(fint3x(fix+vix(7)),ior(fint3y(fiy+viy(3)),fint3z(fiz+viz(3))))
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(4)),fint3z(fiz+viz(4))))
          end select          
          vcboxaddr( 18) = vixyz(7) + 16
          vcboxaddr( 19) = vixyz(7) + 24
          vcboxaddr( 20) = vixyz(7) + 32
          vcboxaddr( 21) = vixyz(7) + 40
          vcboxaddr( 22) = vixyz(7) + 48
          vcboxaddr( 23) = vixyz(7) + 56
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + 8
          vcboxaddr( 26) = vixyz(8) + 32
          vcboxaddr( 27) = vixyz(8) + 40
          dist =208+child
          fwdstart = 19
          case(2)
          select case(int3exist)
          case(1)          
            vixyz(7)=ior(int3x(fix+vix5(7)),ior(int3y(fiy+viy5(3)),int3z(fiz+viz5(3))))
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(4)),int3z(fiz+viz5(4))))
          case(2)          
            vixyz(7)=ior(fint3x(fix+vix5(7)),ior(fint3y(fiy+viy5(3)),fint3z(fiz+viz5(3))))
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(4)),fint3z(fiz+viz5(4))))
          end select
          vcboxaddr( 18) = vixyz(7) + 2
          vcboxaddr( 19) = vixyz(7) + 3
          vcboxaddr( 20) = vixyz(7) + 4
          vcboxaddr( 21) = vixyz(7) + 5
          vcboxaddr( 22) = vixyz(7) + 6
          vcboxaddr( 23) = vixyz(7) + 7
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + 1
          vcboxaddr( 26) = vixyz(8) + 4
          vcboxaddr( 27) = vixyz(8) + 5
          
          fwdstart = 19
        end select
        case(2)

	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(4)),int3z(fiz+pviz(4))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(4)),int3z(fiz+pviz(4))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(4)),fint3z(fiz+pviz(4))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(4)),fint3z(fiz+pviz(4))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(1) + poffset(2)
          vcboxaddr(  3) = vixyz(2)
          vcboxaddr(  4) = vixyz(3)
          vcboxaddr(  5) = vixyz(3) + poffset(2)
          vcboxaddr(  6) = vixyz(3) + poffset(4)
          vcboxaddr(  7) = vixyz(3) + poffset(6)
          vcboxaddr(  8) = vixyz(4)
          vcboxaddr(  9) = vixyz(4) + poffset(4)
          vcboxaddr( 10) = vixyz(5)
          vcboxaddr( 11) = vixyz(5) + poffset(1)
          vcboxaddr( 12) = vixyz(5) + poffset(2)
          vcboxaddr( 13) = vixyz(5) + poffset(3)
          vcboxaddr( 14) = vixyz(6)
          vcboxaddr( 15) = vixyz(6) + poffset(1)
          vcboxaddr( 16) = vixyz(7)
          vcboxaddr( 17) = vixyz(7) + poffset(1)
          vcboxaddr( 18) = vixyz(7) + poffset(2)
          vcboxaddr( 19) = vixyz(7) + poffset(3)
          vcboxaddr( 20) = vixyz(7) + poffset(4)
          vcboxaddr( 21) = vixyz(7) + poffset(5)
          vcboxaddr( 22) = vixyz(7) + poffset(6)
          vcboxaddr( 23) = vixyz(7) + poffset(7)
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + poffset(1)
          vcboxaddr( 26) = vixyz(8) + poffset(4)
          vcboxaddr( 27) = vixyz(8) + poffset(5)
	  
          dist =208+child
          fwdstart = 19
      end select
      end subroutine caladdr27
      subroutine caladdr28(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       16831
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2, 0,-2, 0,-2, 0,-2, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2,-2,-2,-2, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2,-2, 0, 0,-2,-2, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1, 0,-1, 0,-1, 0,-1, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1,-1,-1,-1, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1,-1, 0, 0,-1,-1, 0, 0/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)          
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(7)),int3z(fiz+viz(7))))
          case(2)          
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(7)),fint3z(fiz+viz(7))))     
          end select     
          vcboxaddr( 20) = vixyz(8)
          vcboxaddr( 21) = vixyz(8) + 8
          vcboxaddr( 22) = vixyz(8) +16
          vcboxaddr( 23) = vixyz(8) +24
          vcboxaddr( 24) = vixyz(8) + 32
          vcboxaddr( 25) = vixyz(8) + 40
          vcboxaddr( 26) = vixyz(8) + 48
          vcboxaddr( 27) = vixyz(8) + 56
          dist =216+child
          fwdstart = 21
        case(2)
          select case(int3exist)
          case(1)          
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(7)),int3z(fiz+viz5(7))))
          case(2)          
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(7)),fint3z(fiz+viz5(7))))    
          end select      
          vcboxaddr( 20) = vixyz(8)
          vcboxaddr( 21) = vixyz(8) + 1
          vcboxaddr( 22) = vixyz(8) + 2
          vcboxaddr( 23) = vixyz(8) + 3
          vcboxaddr( 24) = vixyz(8) + 4
          vcboxaddr( 25) = vixyz(8) + 5
          vcboxaddr( 26) = vixyz(8) + 6
          vcboxaddr( 27) = vixyz(8) + 7
	  
          fwdstart = 21	
	end select  
        case(2)

	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(7)),int3z(fiz+pviz(7))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(7)),int3z(fiz+pviz(7))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(7)),fint3z(fiz+pviz(7))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(7)),fint3z(fiz+pviz(7))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(2)
          vcboxaddr(  3) = vixyz(2) + poffset(1)
          vcboxaddr(  4) = vixyz(3)
          vcboxaddr(  5) = vixyz(3) + poffset(4)
          vcboxaddr(  6) = vixyz(4)
          vcboxaddr(  7) = vixyz(4) + poffset(1)
          vcboxaddr(  8) = vixyz(4) + poffset(4)
          vcboxaddr(  9) = vixyz(4) + poffset(5)
          vcboxaddr( 10) = vixyz(5)
          vcboxaddr( 11) = vixyz(5) + poffset(2)
          vcboxaddr( 12) = vixyz(6)
          vcboxaddr( 13) = vixyz(6) + poffset(1)
          vcboxaddr( 14) = vixyz(6) + poffset(2)
          vcboxaddr( 15) = vixyz(6) + poffset(3)
          vcboxaddr( 16) = vixyz(7)
          vcboxaddr( 17) = vixyz(7) + poffset(2)
          vcboxaddr( 18) = vixyz(7) + poffset(4)
          vcboxaddr( 19) = vixyz(7) + poffset(6)
	  vcboxaddr( 20) = vixyz(8)
          vcboxaddr( 21) = vixyz(8) + poffset(1)
          vcboxaddr( 22) = vixyz(8) + poffset(2)
          vcboxaddr( 23) = vixyz(8) + poffset(3)
          vcboxaddr( 24) = vixyz(8) + poffset(4)
          vcboxaddr( 25) = vixyz(8) + poffset(5)
          vcboxaddr( 26) = vixyz(8) + poffset(6)
          vcboxaddr( 27) = vixyz(8) + poffset(7)

          dist =216+child
          fwdstart = 21
      end select
      end subroutine caladdr28
      subroutine caladdr29(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       23083
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2,-2, 2, 2,-2,-2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2,-2,-2,-2, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2, 0,-2, 0,-2, 0,-2, 0/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1,-1, 1, 1,-1,-1, 1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1,-1,-1,-1, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1, 0,-1, 0,-1, 0,-1, 0/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)
            vixyz(6)=ior(int3x(fix+vix(6)),ior(int3y(fiy+viy(6)),int3z(fiz+viz(6))))
            vixyz(7)=ior(int3x(fix+vix(7)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(6)),int3z(fiz+viz(6))))
          case(2)
            vixyz(6)=ior(fint3x(fix+vix(6)),ior(fint3y(fiy+viy(6)),fint3z(fiz+viz(6))))
            vixyz(7)=ior(fint3x(fix+vix(7)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(6)),fint3z(fiz+viz(6))))
          end select          
          vcboxaddr( 15) = vixyz(6) + 8
          vcboxaddr( 16) = vixyz(6) +16
          vcboxaddr( 17) = vixyz(6) +24
          vcboxaddr( 18) = vixyz(6) + 32
          vcboxaddr( 19) = vixyz(6) + 40
          vcboxaddr( 20) = vixyz(6) + 48
          vcboxaddr( 21) = vixyz(6) + 56
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) +16
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) +16
          vcboxaddr( 26) = vixyz(8) + 32
          vcboxaddr( 27) = vixyz(8) + 48
          dist =224+child
          fwdstart =  16
	case(2)
          select case(int3exist)
          case(1)
            vixyz(6)=ior(int3x(fix+vix5(6)),ior(int3y(fiy+viy5(6)),int3z(fiz+viz5(6))))
            vixyz(7)=ior(int3x(fix+vix5(7)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(6)),int3z(fiz+viz5(6))))
          case(2)
            vixyz(6)=ior(fint3x(fix+vix5(6)),ior(fint3y(fiy+viy5(6)),fint3z(fiz+viz5(6))))
            vixyz(7)=ior(fint3x(fix+vix5(7)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(6)),fint3z(fiz+viz5(6))))
          end select
          vcboxaddr( 15) = vixyz(6) + 1
          vcboxaddr( 16) = vixyz(6) + 2
          vcboxaddr( 17) = vixyz(6) + 3
          vcboxaddr( 18) = vixyz(6) + 4
          vcboxaddr( 19) = vixyz(6) + 5
          vcboxaddr( 20) = vixyz(6) + 6
          vcboxaddr( 21) = vixyz(6) + 7
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + 2
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + 2
          vcboxaddr( 26) = vixyz(8) + 4
          vcboxaddr( 27) = vixyz(8) + 6
                 
          fwdstart =  16
	end select  
        case(2)

	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(6)),int3z(fiz+pviz(6))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(6)),int3z(fiz+pviz(6))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(6)),fint3z(fiz+pviz(6))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(6)),fint3z(fiz+pviz(6))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(1) + poffset(1)
          vcboxaddr(  3) = vixyz(2)
          vcboxaddr(  4) = vixyz(2) + poffset(1)
          vcboxaddr(  5) = vixyz(2) + poffset(4)
          vcboxaddr(  6) = vixyz(2) + poffset(5)
          vcboxaddr(  7) = vixyz(3)
          vcboxaddr(  8) = vixyz(4)
          vcboxaddr(  9) = vixyz(4) + poffset(4)
          vcboxaddr( 10) = vixyz(5)
          vcboxaddr( 11) = vixyz(5) + poffset(1)
          vcboxaddr( 12) = vixyz(5) + poffset(2)
          vcboxaddr( 13) = vixyz(5) + poffset(3)
          vcboxaddr( 14) = vixyz(6)
          vcboxaddr( 15) = vixyz(6) + poffset(1)
          vcboxaddr( 16) = vixyz(6) + poffset(2)
          vcboxaddr( 17) = vixyz(6) + poffset(3)
          vcboxaddr( 18) = vixyz(6) + poffset(4)
          vcboxaddr( 19) = vixyz(6) + poffset(5)
          vcboxaddr( 20) = vixyz(6) + poffset(6)
          vcboxaddr( 21) = vixyz(6) + poffset(7)
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + poffset(2)
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + poffset(2)
          vcboxaddr( 26) = vixyz(8) + poffset(4)
          vcboxaddr( 27) = vixyz(8) + poffset(6)
	  
          dist =224+child
          fwdstart =  16
      end select
      end subroutine caladdr29
      subroutine caladdr30(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       32173
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2,-2,-2,-2, 2, 2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2, 2,-2, 2,-2, 2,-2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2,-2, 0, 0,-2,-2, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1,-1,-1,-1, 1, 1, 1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1, 1,-1, 1,-1, 1,-1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1,-1, 0, 0,-1,-1, 0, 0/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

    
      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)          
            vixyz(3)=ior(int3x(fix+vix(3)),ior(int3y(fiy+viy(3)),int3z(fiz+viz(3))))
            vixyz(4)=ior(int3x(fix+vix(4)),ior(int3y(fiy+viy(4)),int3z(fiz+viz(4))))
            vixyz(5)=ior(int3x(fix+vix(5)),ior(int3y(fiy+viy(1)),int3z(fiz+viz(1))))
            vixyz(6)=ior(int3x(fix+vix(6)),ior(int3y(fiy+viy(2)),int3z(fiz+viz(2))))
            vixyz(7)=ior(int3x(fix+vix(7)),ior(int3y(fiy+viy(3)),int3z(fiz+viz(3))))
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(4)),int3z(fiz+viz(4))))
          case(2)          
            vixyz(3)=ior(fint3x(fix+vix(3)),ior(fint3y(fiy+viy(3)),fint3z(fiz+viz(3))))
            vixyz(4)=ior(fint3x(fix+vix(4)),ior(fint3y(fiy+viy(4)),fint3z(fiz+viz(4))))
            vixyz(5)=ior(fint3x(fix+vix(5)),ior(fint3y(fiy+viy(1)),fint3z(fiz+viz(1))))
            vixyz(6)=ior(fint3x(fix+vix(6)),ior(fint3y(fiy+viy(2)),fint3z(fiz+viz(2))))
            vixyz(7)=ior(fint3x(fix+vix(7)),ior(fint3y(fiy+viy(3)),fint3z(fiz+viz(3))))
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(4)),fint3z(fiz+viz(4))))   
          end select       
          vcboxaddr( 10) = vixyz(3) + 24
          vcboxaddr( 11) = vixyz(3) + 32
          vcboxaddr( 12) = vixyz(3) + 40
          vcboxaddr( 13) = vixyz(3) + 48
          vcboxaddr( 14) = vixyz(3) + 56
          vcboxaddr( 15) = vixyz(4)
          vcboxaddr( 16) = vixyz(4) + 8
          vcboxaddr( 17) = vixyz(4) + 32
          vcboxaddr( 18) = vixyz(4) + 40
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) +16
          vcboxaddr( 21) = vixyz(6)
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) +16
          vcboxaddr( 24) = vixyz(7) + 32
          vcboxaddr( 25) = vixyz(7) + 48
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + 32
          dist =232+child
          fwdstart = 11
	  case(2)
          select case(int3exist)
          case(1)          
            vixyz(3)=ior(int3x(fix+vix5(3)),ior(int3y(fiy+viy5(3)),int3z(fiz+viz5(3))))
            vixyz(4)=ior(int3x(fix+vix5(4)),ior(int3y(fiy+viy5(4)),int3z(fiz+viz5(4))))
            vixyz(5)=ior(int3x(fix+vix5(5)),ior(int3y(fiy+viy5(1)),int3z(fiz+viz5(1))))
            vixyz(6)=ior(int3x(fix+vix5(6)),ior(int3y(fiy+viy5(2)),int3z(fiz+viz5(2))))
            vixyz(7)=ior(int3x(fix+vix5(7)),ior(int3y(fiy+viy5(3)),int3z(fiz+viz5(3))))
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(4)),int3z(fiz+viz5(4))))
          case(2)          
            vixyz(3)=ior(fint3x(fix+vix5(3)),ior(fint3y(fiy+viy5(3)),fint3z(fiz+viz5(3))))
            vixyz(4)=ior(fint3x(fix+vix5(4)),ior(fint3y(fiy+viy5(4)),fint3z(fiz+viz5(4))))
            vixyz(5)=ior(fint3x(fix+vix5(5)),ior(fint3y(fiy+viy5(1)),fint3z(fiz+viz5(1))))
            vixyz(6)=ior(fint3x(fix+vix5(6)),ior(fint3y(fiy+viy5(2)),fint3z(fiz+viz5(2))))
            vixyz(7)=ior(fint3x(fix+vix5(7)),ior(fint3y(fiy+viy5(3)),fint3z(fiz+viz5(3))))
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(4)),fint3z(fiz+viz5(4))))
          end select
          vcboxaddr( 10) = vixyz(3) + 3 
          vcboxaddr( 11) = vixyz(3) + 4
          vcboxaddr( 12) = vixyz(3) + 5
          vcboxaddr( 13) = vixyz(3) + 6
          vcboxaddr( 14) = vixyz(3) + 7
          vcboxaddr( 15) = vixyz(4)
          vcboxaddr( 16) = vixyz(4) + 1
          vcboxaddr( 17) = vixyz(4) + 4
          vcboxaddr( 18) = vixyz(4) + 5
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + 2
          vcboxaddr( 21) = vixyz(6)
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + 2
          vcboxaddr( 24) = vixyz(7) + 4
          vcboxaddr( 25) = vixyz(7) + 6
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + 4
          
          fwdstart = 11
	end select  
        case(2)

	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(4)),int3z(fiz+pviz(4))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(4)),int3z(fiz+pviz(4))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(4)),fint3z(fiz+pviz(4))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(4)),fint3z(fiz+pviz(4))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(1) + poffset(1)
          vcboxaddr(  3) = vixyz(1) + poffset(2)
          vcboxaddr(  4) = vixyz(1) + poffset(3)
          vcboxaddr(  5) = vixyz(2)
          vcboxaddr(  6) = vixyz(2) + poffset(1)
          vcboxaddr(  7) = vixyz(3)
          vcboxaddr(  8) = vixyz(3) + poffset(1)
          vcboxaddr(  9) = vixyz(3) + poffset(2)
          vcboxaddr( 10) = vixyz(3) + poffset(3)
          vcboxaddr( 11) = vixyz(3) + poffset(4)
          vcboxaddr( 12) = vixyz(3) + poffset(5)
          vcboxaddr( 13) = vixyz(3) + poffset(6)
          vcboxaddr( 14) = vixyz(3) + poffset(7)
          vcboxaddr( 15) = vixyz(4)
          vcboxaddr( 16) = vixyz(4) + poffset(1)
          vcboxaddr( 17) = vixyz(4) + poffset(4)
          vcboxaddr( 18) = vixyz(4) + poffset(5)
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + poffset(2)
          vcboxaddr( 21) = vixyz(6)
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + poffset(2)
          vcboxaddr( 24) = vixyz(7) + poffset(4)
          vcboxaddr( 25) = vixyz(7) + poffset(6)
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + poffset(4)
	  
          dist =232+child          
          fwdstart = 11
      end select
      end subroutine caladdr30
      subroutine caladdr31(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       34543
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2,-2, 2, 2,-2,-2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2, 2,-2, 2,-2, 2,-2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2,-2,-2,-2, 2, 2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1,-1, 1, 1,-1,-1, 1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1, 1,-1, 1,-1, 1,-1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1,-1,-1,-1, 1, 1, 1, 1/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+vix(1)),ior(int3y(fiy+viy(1)),int3z(fiz+viz(1))))
            vixyz(2)=ior(int3x(fix+vix(2)),ior(int3y(fiy+viy(2)),int3z(fiz+viz(2))))
            vixyz(3)=ior(int3x(fix+vix(3)),ior(int3y(fiy+viy(1)),int3z(fiz+viz(1))))
            vixyz(4)=ior(int3x(fix+vix(4)),ior(int3y(fiy+viy(2)),int3z(fiz+viz(2))))
            vixyz(5)=ior(int3x(fix+vix(5)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(6)=ior(int3x(fix+vix(6)),ior(int3y(fiy+viy(6)),int3z(fiz+viz(6))))
            vixyz(7)=ior(int3x(fix+vix(7)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(6)),int3z(fiz+viz(6))))
          case(2)
            vixyz(1)=ior(fint3x(fix+vix(1)),ior(fint3y(fiy+viy(1)),fint3z(fiz+viz(1))))
            vixyz(2)=ior(fint3x(fix+vix(2)),ior(fint3y(fiy+viy(2)),fint3z(fiz+viz(2))))
            vixyz(3)=ior(fint3x(fix+vix(3)),ior(fint3y(fiy+viy(1)),fint3z(fiz+viz(1))))
            vixyz(4)=ior(fint3x(fix+vix(4)),ior(fint3y(fiy+viy(2)),fint3z(fiz+viz(2))))
            vixyz(5)=ior(fint3x(fix+vix(5)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(6)=ior(fint3x(fix+vix(6)),ior(fint3y(fiy+viy(6)),fint3z(fiz+viz(6))))
            vixyz(7)=ior(fint3x(fix+vix(7)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(6)),fint3z(fiz+viz(6))))
          end select          
          vcboxaddr(  8) = vixyz(1) + 56
          vcboxaddr(  9) = vixyz(2)
          vcboxaddr( 10) = vixyz(2) + 8
          vcboxaddr( 11) = vixyz(2) + 32
          vcboxaddr( 12) = vixyz(2) + 40
          vcboxaddr( 13) = vixyz(3)
          vcboxaddr( 14) = vixyz(3) +16
          vcboxaddr( 15) = vixyz(3) + 32
          vcboxaddr( 16) = vixyz(3) + 48
          vcboxaddr( 17) = vixyz(4)
          vcboxaddr( 18) = vixyz(4) + 32
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + 8
          vcboxaddr( 21) = vixyz(5) +16
          vcboxaddr( 22) = vixyz(5) +24
          vcboxaddr( 23) = vixyz(6)
          vcboxaddr( 24) = vixyz(6) + 8
          vcboxaddr( 25) = vixyz(7)
          vcboxaddr( 26) = vixyz(7) +16
          vcboxaddr( 27) = vixyz(8)
          dist =240+child
          fwdstart =  9
	case(2)
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+vix5(1)),ior(int3y(fiy+viy5(1)),int3z(fiz+viz5(1))))
            vixyz(2)=ior(int3x(fix+vix5(2)),ior(int3y(fiy+viy5(2)),int3z(fiz+viz5(2))))
            vixyz(3)=ior(int3x(fix+vix5(3)),ior(int3y(fiy+viy5(1)),int3z(fiz+viz5(1))))
            vixyz(4)=ior(int3x(fix+vix5(4)),ior(int3y(fiy+viy5(2)),int3z(fiz+viz5(2))))
            vixyz(5)=ior(int3x(fix+vix5(5)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(6)=ior(int3x(fix+vix5(6)),ior(int3y(fiy+viy5(6)),int3z(fiz+viz5(6))))
            vixyz(7)=ior(int3x(fix+vix5(7)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(6)),int3z(fiz+viz5(6))))
          case(2)
            vixyz(1)=ior(fint3x(fix+vix5(1)),ior(fint3y(fiy+viy5(1)),fint3z(fiz+viz5(1))))
            vixyz(2)=ior(fint3x(fix+vix5(2)),ior(fint3y(fiy+viy5(2)),fint3z(fiz+viz5(2))))
            vixyz(3)=ior(fint3x(fix+vix5(3)),ior(fint3y(fiy+viy5(1)),fint3z(fiz+viz5(1))))
            vixyz(4)=ior(fint3x(fix+vix5(4)),ior(fint3y(fiy+viy5(2)),fint3z(fiz+viz5(2))))
            vixyz(5)=ior(fint3x(fix+vix5(5)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(6)=ior(fint3x(fix+vix5(6)),ior(fint3y(fiy+viy5(6)),fint3z(fiz+viz5(6))))
            vixyz(7)=ior(fint3x(fix+vix5(7)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(6)),fint3z(fiz+viz5(6))))
          end select
          vcboxaddr(  8) = vixyz(1) + 7
          vcboxaddr(  9) = vixyz(2)
          vcboxaddr( 10) = vixyz(2) + 1
          vcboxaddr( 11) = vixyz(2) + 4
          vcboxaddr( 12) = vixyz(2) + 5
          vcboxaddr( 13) = vixyz(3)
          vcboxaddr( 14) = vixyz(3) + 2
          vcboxaddr( 15) = vixyz(3) + 4
          vcboxaddr( 16) = vixyz(3) + 6
          vcboxaddr( 17) = vixyz(4)
          vcboxaddr( 18) = vixyz(4) + 4
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + 1
          vcboxaddr( 21) = vixyz(5) + 2
          vcboxaddr( 22) = vixyz(5) + 3
          vcboxaddr( 23) = vixyz(6)
          vcboxaddr( 24) = vixyz(6) + 1
          vcboxaddr( 25) = vixyz(7)
          vcboxaddr( 26) = vixyz(7) + 2
          vcboxaddr( 27) = vixyz(8)
          
          fwdstart =  9
	end select  
        case(2)

	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(6)),int3z(fiz+pviz(6))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(6)),int3z(fiz+pviz(6))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(6)),fint3z(fiz+pviz(6))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(6)),fint3z(fiz+pviz(6))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(1) + poffset(1)
          vcboxaddr(  3) = vixyz(1) + poffset(2)
          vcboxaddr(  4) = vixyz(1) + poffset(3)
          vcboxaddr(  5) = vixyz(1) + poffset(4)
          vcboxaddr(  6) = vixyz(1) + poffset(5)
          vcboxaddr(  7) = vixyz(1) + poffset(6)
          vcboxaddr(  8) = vixyz(1) + poffset(7)
          vcboxaddr(  9) = vixyz(2)
          vcboxaddr( 10) = vixyz(2) + poffset(1)
          vcboxaddr( 11) = vixyz(2) + poffset(4)
          vcboxaddr( 12) = vixyz(2) + poffset(5)
          vcboxaddr( 13) = vixyz(3)
          vcboxaddr( 14) = vixyz(3) + poffset(2)
          vcboxaddr( 15) = vixyz(3) + poffset(4)
          vcboxaddr( 16) = vixyz(3) + poffset(6)
          vcboxaddr( 17) = vixyz(4)
          vcboxaddr( 18) = vixyz(4) + poffset(4)
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + poffset(1)
          vcboxaddr( 21) = vixyz(5) + poffset(2)
          vcboxaddr( 22) = vixyz(5) + poffset(3)
          vcboxaddr( 23) = vixyz(6)
          vcboxaddr( 24) = vixyz(6) + poffset(1)
          vcboxaddr( 25) = vixyz(7)
          vcboxaddr( 26) = vixyz(7) + poffset(2)
          vcboxaddr( 27) = vixyz(8)

          dist =240+child
          fwdstart =  9
      end select
      end subroutine caladdr31
      subroutine caladdr32(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       16459
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2,-2,-2,-2, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2, 0,-2, 0,-2, 0,-2, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2,-2, 0, 0,-2,-2, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1,-1,-1,-1, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1, 0,-1, 0,-1, 0,-1, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1,-1, 0, 0,-1,-1, 0, 0/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(4)),int3z(fiz+viz(4))))
          case(2)
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(4)),fint3z(fiz+viz(4))))     
          end select     
          vcboxaddr( 20) = vixyz(8)
          vcboxaddr( 21) = vixyz(8) + 8
          vcboxaddr( 22) = vixyz(8) + 16
          vcboxaddr( 23) = vixyz(8) + 24
          vcboxaddr( 24) = vixyz(8) + 32
          vcboxaddr( 25) = vixyz(8) + 40
          vcboxaddr( 26) = vixyz(8) + 48
          vcboxaddr( 27) = vixyz(8) + 56
          dist =248+child
          fwdstart = 21	
	  case(2)
          select case(int3exist)
          case(1)          
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(4)),int3z(fiz+viz5(4))))
          case(2)          
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(4)),fint3z(fiz+viz5(4))))    
          end select      
          vcboxaddr( 20) = vixyz(8)
          vcboxaddr( 21) = vixyz(8) + 1
          vcboxaddr( 22) = vixyz(8) + 2
          vcboxaddr( 23) = vixyz(8) + 3
          vcboxaddr( 24) = vixyz(8) + 4
          vcboxaddr( 25) = vixyz(8) + 5
          vcboxaddr( 26) = vixyz(8) + 6
          vcboxaddr( 27) = vixyz(8) + 7
          
          fwdstart = 21	
	end select  
        case(2)

	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(4)),int3z(fiz+pviz(4))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(4)),int3z(fiz+pviz(4))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(4)),fint3z(fiz+pviz(4))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(4)),fint3z(fiz+pviz(4))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(2)
          vcboxaddr(  3) = vixyz(2) + poffset(2)
          vcboxaddr(  4) = vixyz(3)
          vcboxaddr(  5) = vixyz(3) + poffset(4)
          vcboxaddr(  6) = vixyz(4)
          vcboxaddr(  7) = vixyz(4) + poffset(2)
          vcboxaddr(  8) = vixyz(4) + poffset(4)
          vcboxaddr(  9) = vixyz(4) + poffset(6)
          vcboxaddr( 10) = vixyz(5)
          vcboxaddr( 11) = vixyz(5) + poffset(1)
          vcboxaddr( 12) = vixyz(6)
          vcboxaddr( 13) = vixyz(6) + poffset(1)
          vcboxaddr( 14) = vixyz(6) + poffset(2)
          vcboxaddr( 15) = vixyz(6) + poffset(3)
          vcboxaddr( 16) = vixyz(7)
          vcboxaddr( 17) = vixyz(7) + poffset(1)
          vcboxaddr( 18) = vixyz(7) + poffset(4)
          vcboxaddr( 19) = vixyz(7) + poffset(5)
          vcboxaddr( 20) = vixyz(8)
          vcboxaddr( 21) = vixyz(8) + poffset(1)
          vcboxaddr( 22) = vixyz(8) + poffset(2)
          vcboxaddr( 23) = vixyz(8) + poffset(3)
          vcboxaddr( 24) = vixyz(8) + poffset(4)
          vcboxaddr( 25) = vixyz(8) + poffset(5)
          vcboxaddr( 26) = vixyz(8) + poffset(6)
          vcboxaddr( 27) = vixyz(8) + poffset(7)
          dist =248+child
          fwdstart = 21	
      end select
      end subroutine caladdr32
      subroutine caladdr33(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       21907
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2,-2,-2,-2, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2,-2, 2, 2,-2,-2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2, 0,-2, 0,-2, 0,-2, 0/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1,-1,-1,-1, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1,-1, 1, 1,-1,-1, 1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1, 0,-1, 0,-1, 0,-1, 0/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)          
            vixyz(6)=ior(int3x(fix+vix(6)),ior(int3y(fiy+viy(2)),int3z(fiz+viz(2))))
            vixyz(7)=ior(int3x(fix+vix(7)),ior(int3y(fiy+viy(3)),int3z(fiz+viz(3))))
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(4)),int3z(fiz+viz(4))))
          case(2)          
            vixyz(6)=ior(fint3x(fix+vix(6)),ior(fint3y(fiy+viy(2)),fint3z(fiz+viz(2))))
            vixyz(7)=ior(fint3x(fix+vix(7)),ior(fint3y(fiy+viy(3)),fint3z(fiz+viz(3))))
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(4)),fint3z(fiz+viz(4)))) 
          end select
          vcboxaddr( 16) = vixyz(6) + 16
          vcboxaddr( 17) = vixyz(6) + 24
          vcboxaddr( 18) = vixyz(6) + 32
          vcboxaddr( 19) = vixyz(6) + 40
          vcboxaddr( 20) = vixyz(6) + 48
          vcboxaddr( 21) = vixyz(6) + 56
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + 8
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + 8
          vcboxaddr( 26) = vixyz(8) + 32
          vcboxaddr( 27) = vixyz(8) + 40
          dist =256+child
          fwdstart = 17
	  case(2)
          select case(int3exist)
          case(1)          
            vixyz(6)=ior(int3x(fix+vix5(6)),ior(int3y(fiy+viy5(2)),int3z(fiz+viz5(2))))
            vixyz(7)=ior(int3x(fix+vix5(7)),ior(int3y(fiy+viy5(3)),int3z(fiz+viz5(3))))
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(4)),int3z(fiz+viz5(4))))
          case(2)          
            vixyz(6)=ior(fint3x(fix+vix5(6)),ior(fint3y(fiy+viy5(2)),fint3z(fiz+viz5(2))))
            vixyz(7)=ior(fint3x(fix+vix5(7)),ior(fint3y(fiy+viy5(3)),fint3z(fiz+viz5(3))))
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(4)),fint3z(fiz+viz5(4))))
          end select
          vcboxaddr( 16) = vixyz(6) + 2
          vcboxaddr( 17) = vixyz(6) + 3
          vcboxaddr( 18) = vixyz(6) + 4
          vcboxaddr( 19) = vixyz(6) + 5
          vcboxaddr( 20) = vixyz(6) + 6
          vcboxaddr( 21) = vixyz(6) + 7
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + 1
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + 1
          vcboxaddr( 26) = vixyz(8) + 4
          vcboxaddr( 27) = vixyz(8) + 5
         
          fwdstart = 17
	end select  
        case(2)

	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(4)),int3z(fiz+pviz(4))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(4)),int3z(fiz+pviz(4))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(4)),fint3z(fiz+pviz(4))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(4)),fint3z(fiz+pviz(4))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(1) + poffset(2)
          vcboxaddr(  3) = vixyz(2)
          vcboxaddr(  4) = vixyz(2) + poffset(2)
          vcboxaddr(  5) = vixyz(2) + poffset(4)
          vcboxaddr(  6) = vixyz(2) + poffset(6)
          vcboxaddr(  7) = vixyz(3)
          vcboxaddr(  8) = vixyz(4)
          vcboxaddr(  9) = vixyz(4) + poffset(4)
          vcboxaddr( 10) = vixyz(5)
          vcboxaddr( 11) = vixyz(5) + poffset(1)
          vcboxaddr( 12) = vixyz(5) + poffset(2)
          vcboxaddr( 13) = vixyz(5) + poffset(3)
          vcboxaddr( 14) = vixyz(6)
          vcboxaddr( 15) = vixyz(6) + poffset(1)
          vcboxaddr( 16) = vixyz(6) + poffset(2)
          vcboxaddr( 17) = vixyz(6) + poffset(3)
          vcboxaddr( 18) = vixyz(6) + poffset(4)
          vcboxaddr( 19) = vixyz(6) + poffset(5)
          vcboxaddr( 20) = vixyz(6) + poffset(6)
          vcboxaddr( 21) = vixyz(6) + poffset(7)
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + poffset(1)
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + poffset(1)
          vcboxaddr( 26) = vixyz(8) + poffset(4)
          vcboxaddr( 27) = vixyz(8) + poffset(5)
	  
          dist =256+child
          fwdstart = 17
      end select
      end subroutine caladdr33
      subroutine caladdr34(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       31522
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2, 2,-2, 2,-2, 2,-2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2,-2,-2,-2, 2, 2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2,-2, 0, 0,-2,-2, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1, 1,-1, 1,-1, 1,-1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1,-1,-1,-1, 1, 1, 1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1,-1, 0, 0,-1,-1, 0, 0/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)
            vixyz(3)=ior(int3x(fix+vix(3)),ior(int3y(fiy+viy(3)),int3z(fiz+viz(3))))
            vixyz(4)=ior(int3x(fix+vix(4)),ior(int3y(fiy+viy(3)),int3z(fiz+viz(3))))
            vixyz(5)=ior(int3x(fix+vix(5)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(6)=ior(int3x(fix+vix(6)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(7)=ior(int3x(fix+vix(7)),ior(int3y(fiy+viy(7)),int3z(fiz+viz(7))))
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(7)),int3z(fiz+viz(7))))
          case(2)
            vixyz(3)=ior(fint3x(fix+vix(3)),ior(fint3y(fiy+viy(3)),fint3z(fiz+viz(3))))
            vixyz(4)=ior(fint3x(fix+vix(4)),ior(fint3y(fiy+viy(3)),fint3z(fiz+viz(3))))
            vixyz(5)=ior(fint3x(fix+vix(5)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(6)=ior(fint3x(fix+vix(6)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(7)=ior(fint3x(fix+vix(7)),ior(fint3y(fiy+viy(7)),fint3z(fiz+viz(7))))
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(7)),fint3z(fiz+viz(7))))      
          end select    
          vcboxaddr( 10) = vixyz(3) + 24
          vcboxaddr( 11) = vixyz(3) + 32
          vcboxaddr( 12) = vixyz(3) + 40
          vcboxaddr( 13) = vixyz(3) + 48
          vcboxaddr( 14) = vixyz(3) + 56
          vcboxaddr( 15) = vixyz(4)
          vcboxaddr( 16) = vixyz(4) +16
          vcboxaddr( 17) = vixyz(4) + 32
          vcboxaddr( 18) = vixyz(4) + 48
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + 8
          vcboxaddr( 21) = vixyz(6)
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + 8
          vcboxaddr( 24) = vixyz(7) + 32
          vcboxaddr( 25) = vixyz(7) + 40
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + 32
          dist =264+child
          fwdstart = 11
	  case(2)
          select case(int3exist)
          case(1)
            vixyz(3)=ior(int3x(fix+vix5(3)),ior(int3y(fiy+viy5(3)),int3z(fiz+viz5(3))))
            vixyz(4)=ior(int3x(fix+vix5(4)),ior(int3y(fiy+viy5(3)),int3z(fiz+viz5(3))))
            vixyz(5)=ior(int3x(fix+vix5(5)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(6)=ior(int3x(fix+vix5(6)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(7)=ior(int3x(fix+vix5(7)),ior(int3y(fiy+viy5(7)),int3z(fiz+viz5(7))))
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(7)),int3z(fiz+viz5(7))))
          case(2)
            vixyz(3)=ior(fint3x(fix+vix5(3)),ior(fint3y(fiy+viy5(3)),fint3z(fiz+viz5(3))))
            vixyz(4)=ior(fint3x(fix+vix5(4)),ior(fint3y(fiy+viy5(3)),fint3z(fiz+viz5(3))))
            vixyz(5)=ior(fint3x(fix+vix5(5)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(6)=ior(fint3x(fix+vix5(6)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(7)=ior(fint3x(fix+vix5(7)),ior(fint3y(fiy+viy5(7)),fint3z(fiz+viz5(7))))
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(7)),fint3z(fiz+viz5(7))))
          end select
          vcboxaddr( 10) = vixyz(3) + 3
          vcboxaddr( 11) = vixyz(3) + 4
          vcboxaddr( 12) = vixyz(3) + 5
          vcboxaddr( 13) = vixyz(3) + 6
          vcboxaddr( 14) = vixyz(3) + 7
          vcboxaddr( 15) = vixyz(4)
          vcboxaddr( 16) = vixyz(4) + 2
          vcboxaddr( 17) = vixyz(4) + 4
          vcboxaddr( 18) = vixyz(4) + 6
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + 1
          vcboxaddr( 21) = vixyz(6)
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + 1
          vcboxaddr( 24) = vixyz(7) + 4
          vcboxaddr( 25) = vixyz(7) + 5
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + 4
          
          fwdstart = 11
	end select  
        case(2)

	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(7)),int3z(fiz+pviz(7))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(7)),int3z(fiz+pviz(7))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(7)),fint3z(fiz+pviz(7))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(7)),fint3z(fiz+pviz(7))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(1) + poffset(1)
          vcboxaddr(  3) = vixyz(1) + poffset(2)
          vcboxaddr(  4) = vixyz(1) + poffset(3)
          vcboxaddr(  5) = vixyz(2)
          vcboxaddr(  6) = vixyz(2) + poffset(2)
          vcboxaddr(  7) = vixyz(3)
          vcboxaddr(  8) = vixyz(3) + poffset(1)
          vcboxaddr(  9) = vixyz(3) + poffset(2)
          vcboxaddr( 10) = vixyz(3) + poffset(3)
          vcboxaddr( 11) = vixyz(3) + poffset(4)
          vcboxaddr( 12) = vixyz(3) + poffset(5)
          vcboxaddr( 13) = vixyz(3) + poffset(6)
          vcboxaddr( 14) = vixyz(3) + poffset(7)
          vcboxaddr( 15) = vixyz(4)
          vcboxaddr( 16) = vixyz(4) + poffset(2)
          vcboxaddr( 17) = vixyz(4) + poffset(4)
          vcboxaddr( 18) = vixyz(4) + poffset(6)
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + poffset(1)
          vcboxaddr( 21) = vixyz(6)
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + poffset(1)
          vcboxaddr( 24) = vixyz(7) + poffset(4)
          vcboxaddr( 25) = vixyz(7) + poffset(5)
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + poffset(4)
	  
          dist =264+child
          fwdstart = 11
      end select
      end subroutine caladdr34
      subroutine caladdr35(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       29149
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2,-2,-2,-2, 2, 2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2, 0,-2, 0,-2, 0,-2, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2,-2, 0, 0,-2,-2, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1,-1,-1,-1, 1, 1, 1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1, 0,-1, 0,-1, 0,-1, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1,-1, 0, 0,-1,-1, 0, 0/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)          
            vixyz(4)=ior(int3x(fix+vix(4)),ior(int3y(fiy+viy(4)),int3z(fiz+viz(4))))
            vixyz(5)=ior(int3x(fix+vix(5)),ior(int3y(fiy+viy(1)),int3z(fiz+viz(1))))
            vixyz(6)=ior(int3x(fix+vix(6)),ior(int3y(fiy+viy(2)),int3z(fiz+viz(2))))
            vixyz(7)=ior(int3x(fix+vix(7)),ior(int3y(fiy+viy(3)),int3z(fiz+viz(3))))
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(4)),int3z(fiz+viz(4))))
          case(2)          
            vixyz(4)=ior(fint3x(fix+vix(4)),ior(fint3y(fiy+viy(4)),fint3z(fiz+viz(4))))
            vixyz(5)=ior(fint3x(fix+vix(5)),ior(fint3y(fiy+viy(1)),fint3z(fiz+viz(1))))
            vixyz(6)=ior(fint3x(fix+vix(6)),ior(fint3y(fiy+viy(2)),fint3z(fiz+viz(2))))
            vixyz(7)=ior(fint3x(fix+vix(7)),ior(fint3y(fiy+viy(3)),fint3z(fiz+viz(3))))
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(4)),fint3z(fiz+viz(4))))
          end select
          vcboxaddr( 12) = vixyz(4) + 8
          vcboxaddr( 13) = vixyz(4) +16
          vcboxaddr( 14) = vixyz(4) +24
          vcboxaddr( 15) = vixyz(4) + 32
          vcboxaddr( 16) = vixyz(4) + 40
          vcboxaddr( 17) = vixyz(4) + 48
          vcboxaddr( 18) = vixyz(4) + 56
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(6)
          vcboxaddr( 21) = vixyz(6) +16
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + 32
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) +16
          vcboxaddr( 26) = vixyz(8) + 32
          vcboxaddr( 27) = vixyz(8) + 48
          dist =272+child
          fwdstart = 13
	  case(2)
          select case(int3exist)
          case(1)
	    vixyz(4)=ior(int3x(fix+vix5(4)),ior(int3y(fiy+viy5(4)),int3z(fiz+viz5(4))))
            vixyz(5)=ior(int3x(fix+vix5(5)),ior(int3y(fiy+viy5(1)),int3z(fiz+viz5(1))))
            vixyz(6)=ior(int3x(fix+vix5(6)),ior(int3y(fiy+viy5(2)),int3z(fiz+viz5(2))))
            vixyz(7)=ior(int3x(fix+vix5(7)),ior(int3y(fiy+viy5(3)),int3z(fiz+viz5(3))))
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(4)),int3z(fiz+viz5(4))))
          case(2)
	    vixyz(4)=ior(fint3x(fix+vix5(4)),ior(fint3y(fiy+viy5(4)),fint3z(fiz+viz5(4))))
            vixyz(5)=ior(fint3x(fix+vix5(5)),ior(fint3y(fiy+viy5(1)),fint3z(fiz+viz5(1))))
            vixyz(6)=ior(fint3x(fix+vix5(6)),ior(fint3y(fiy+viy5(2)),fint3z(fiz+viz5(2))))
            vixyz(7)=ior(fint3x(fix+vix5(7)),ior(fint3y(fiy+viy5(3)),fint3z(fiz+viz5(3))))
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(4)),fint3z(fiz+viz5(4))))
          end select
          vcboxaddr( 12) = vixyz(4) + 1
          vcboxaddr( 13) = vixyz(4) + 2
          vcboxaddr( 14) = vixyz(4) + 3
          vcboxaddr( 15) = vixyz(4) + 4
          vcboxaddr( 16) = vixyz(4) + 5
          vcboxaddr( 17) = vixyz(4) + 6
          vcboxaddr( 18) = vixyz(4) + 7
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(6)
          vcboxaddr( 21) = vixyz(6) + 2
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + 4
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + 2
          vcboxaddr( 26) = vixyz(8) + 4
          vcboxaddr( 27) = vixyz(8) + 6
          
          fwdstart = 13
	end select  
        case(2)
	
	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(4)),int3z(fiz+pviz(4))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(4)),int3z(fiz+pviz(4))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(4)),fint3z(fiz+pviz(4))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(4)),fint3z(fiz+pviz(4))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(1) + poffset(1)
          vcboxaddr(  3) = vixyz(2)
          vcboxaddr(  4) = vixyz(2) + poffset(1)
          vcboxaddr(  5) = vixyz(2) + poffset(2)
          vcboxaddr(  6) = vixyz(2) + poffset(3)
          vcboxaddr(  7) = vixyz(3)
          vcboxaddr(  8) = vixyz(3) + poffset(1)
          vcboxaddr(  9) = vixyz(3) + poffset(4)
          vcboxaddr( 10) = vixyz(3) + poffset(5)
          vcboxaddr( 11) = vixyz(4)
          vcboxaddr( 12) = vixyz(4) + poffset(1)
          vcboxaddr( 13) = vixyz(4) + poffset(2)
          vcboxaddr( 14) = vixyz(4) + poffset(3)
          vcboxaddr( 15) = vixyz(4) + poffset(4)
          vcboxaddr( 16) = vixyz(4) + poffset(5)
          vcboxaddr( 17) = vixyz(4) + poffset(6)
          vcboxaddr( 18) = vixyz(4) + poffset(7)
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(6)
          vcboxaddr( 21) = vixyz(6) + poffset(2)
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + poffset(4)
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + poffset(2)
          vcboxaddr( 26) = vixyz(8) + poffset(4)
          vcboxaddr( 27) = vixyz(8) + poffset(6)

          dist =272+child
          fwdstart = 13
      end select
      end subroutine caladdr35
      subroutine caladdr36(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       28117
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2, 0,-2, 0,-2, 0,-2, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2,-2,-2,-2, 2, 2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2,-2, 0, 0,-2,-2, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1, 0,-1, 0,-1, 0,-1, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1,-1,-1,-1, 1, 1, 1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1,-1, 0, 0,-1,-1, 0, 0/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)
            vixyz(4)=ior(int3x(fix+vix(4)),ior(int3y(fiy+viy(3)),int3z(fiz+viz(3))))
            vixyz(5)=ior(int3x(fix+vix(5)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(6)=ior(int3x(fix+vix(6)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(7)=ior(int3x(fix+vix(7)),ior(int3y(fiy+viy(7)),int3z(fiz+viz(7))))
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(7)),int3z(fiz+viz(7))))
          case(2)
            vixyz(4)=ior(fint3x(fix+vix(4)),ior(fint3y(fiy+viy(3)),fint3z(fiz+viz(3))))
            vixyz(5)=ior(fint3x(fix+vix(5)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(6)=ior(fint3x(fix+vix(6)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(7)=ior(fint3x(fix+vix(7)),ior(fint3y(fiy+viy(7)),fint3z(fiz+viz(7))))
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(7)),fint3z(fiz+viz(7))))
          end select
          vcboxaddr( 13) = vixyz(4) + 16
          vcboxaddr( 14) = vixyz(4) + 24
          vcboxaddr( 15) = vixyz(4) + 32
          vcboxaddr( 16) = vixyz(4) + 40
          vcboxaddr( 17) = vixyz(4) + 48
          vcboxaddr( 18) = vixyz(4) + 56
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(6)
          vcboxaddr( 21) = vixyz(6) + 8
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + 32
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + 8
          vcboxaddr( 26) = vixyz(8) + 32
          vcboxaddr( 27) = vixyz(8) + 40
          dist =280+child
          fwdstart = 14
	 case(2)
          select case(int3exist)
          case(1)
            vixyz(4)=ior(int3x(fix+vix5(4)),ior(int3y(fiy+viy5(3)),int3z(fiz+viz5(3))))
            vixyz(5)=ior(int3x(fix+vix5(5)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(6)=ior(int3x(fix+vix5(6)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(7)=ior(int3x(fix+vix5(7)),ior(int3y(fiy+viy5(7)),int3z(fiz+viz5(7))))
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(7)),int3z(fiz+viz5(7))))
          case(2)
            vixyz(4)=ior(fint3x(fix+vix5(4)),ior(fint3y(fiy+viy5(3)),fint3z(fiz+viz5(3))))
            vixyz(5)=ior(fint3x(fix+vix5(5)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(6)=ior(fint3x(fix+vix5(6)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(7)=ior(fint3x(fix+vix5(7)),ior(fint3y(fiy+viy5(7)),fint3z(fiz+viz5(7))))
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(7)),fint3z(fiz+viz5(7))))
          end select
          vcboxaddr( 13) = vixyz(4) + 2
          vcboxaddr( 14) = vixyz(4) + 3
          vcboxaddr( 15) = vixyz(4) + 4
          vcboxaddr( 16) = vixyz(4) + 5
          vcboxaddr( 17) = vixyz(4) + 6
          vcboxaddr( 18) = vixyz(4) + 7
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(6)
          vcboxaddr( 21) = vixyz(6) + 1
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + 4
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + 1
          vcboxaddr( 26) = vixyz(8) + 4
          vcboxaddr( 27) = vixyz(8) + 5
         
          fwdstart = 14
	end select  
        case(2)

	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(7)),int3z(fiz+pviz(7))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(7)),int3z(fiz+pviz(7))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(7)),fint3z(fiz+pviz(7))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(7)),fint3z(fiz+pviz(7))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(1) + poffset(2)
          vcboxaddr(  3) = vixyz(2)
          vcboxaddr(  4) = vixyz(2) + poffset(1)
          vcboxaddr(  5) = vixyz(2) + poffset(2)
          vcboxaddr(  6) = vixyz(2) + poffset(3)
          vcboxaddr(  7) = vixyz(3)
          vcboxaddr(  8) = vixyz(3) + poffset(2)
          vcboxaddr(  9) = vixyz(3) + poffset(4)
          vcboxaddr( 10) = vixyz(3) + poffset(6)
          vcboxaddr( 11) = vixyz(4)
          vcboxaddr( 12) = vixyz(4) + poffset(1)
          vcboxaddr( 13) = vixyz(4) + poffset(2)
          vcboxaddr( 14) = vixyz(4) + poffset(3)
          vcboxaddr( 15) = vixyz(4) + poffset(4)
          vcboxaddr( 16) = vixyz(4) + poffset(5)
          vcboxaddr( 17) = vixyz(4) + poffset(6)
          vcboxaddr( 18) = vixyz(4) + poffset(7)
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(6)
          vcboxaddr( 21) = vixyz(6) + poffset(1)
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + poffset(4)
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + poffset(1)
          vcboxaddr( 26) = vixyz(8) + poffset(4)
          vcboxaddr( 27) = vixyz(8) + poffset(5)
	  
          dist =280+child
          fwdstart = 14
      end select
      end subroutine caladdr36
      subroutine caladdr37(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       33937
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2,-2, 2, 2,-2,-2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2,-2,-2,-2, 2, 2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2, 0,-2, 0,-2, 0,-2, 0/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1,-1, 1, 1,-1,-1, 1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1,-1,-1,-1, 1, 1, 1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1, 0,-1, 0,-1, 0,-1, 0/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)          
            vixyz(2)=ior(int3x(fix+vix(2)),ior(int3y(fiy+viy(2)),int3z(fiz+viz(2))))
            vixyz(3)=ior(int3x(fix+vix(3)),ior(int3y(fiy+viy(1)),int3z(fiz+viz(1))))
            vixyz(4)=ior(int3x(fix+vix(4)),ior(int3y(fiy+viy(2)),int3z(fiz+viz(2))))
            vixyz(5)=ior(int3x(fix+vix(5)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(6)=ior(int3x(fix+vix(6)),ior(int3y(fiy+viy(6)),int3z(fiz+viz(6))))
            vixyz(7)=ior(int3x(fix+vix(7)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(6)),int3z(fiz+viz(6))))
          case(2)          
            vixyz(2)=ior(fint3x(fix+vix(2)),ior(fint3y(fiy+viy(2)),fint3z(fiz+viz(2))))
            vixyz(3)=ior(fint3x(fix+vix(3)),ior(fint3y(fiy+viy(1)),fint3z(fiz+viz(1))))
            vixyz(4)=ior(fint3x(fix+vix(4)),ior(fint3y(fiy+viy(2)),fint3z(fiz+viz(2))))
            vixyz(5)=ior(fint3x(fix+vix(5)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(6)=ior(fint3x(fix+vix(6)),ior(fint3y(fiy+viy(6)),fint3z(fiz+viz(6))))
            vixyz(7)=ior(fint3x(fix+vix(7)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(6)),fint3z(fiz+viz(6)))) 
          end select         
          vcboxaddr(  8) = vixyz(2) + 24
          vcboxaddr(  9) = vixyz(2) + 32
          vcboxaddr( 10) = vixyz(2) + 40
          vcboxaddr( 11) = vixyz(2) + 48
          vcboxaddr( 12) = vixyz(2) + 56
          vcboxaddr( 13) = vixyz(3)
          vcboxaddr( 14) = vixyz(3) + 16
          vcboxaddr( 15) = vixyz(4)
          vcboxaddr( 16) = vixyz(4) +16
          vcboxaddr( 17) = vixyz(4) + 32
          vcboxaddr( 18) = vixyz(4) + 48
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + 8
          vcboxaddr( 21) = vixyz(6)
          vcboxaddr( 22) = vixyz(6) + 8
          vcboxaddr( 23) = vixyz(6) + 32
          vcboxaddr( 24) = vixyz(6) + 40
          vcboxaddr( 25) = vixyz(7)
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + 32
          dist =288+child
          fwdstart =  9
        case(2)
          select case(int3exist)
          case(1)          
            vixyz(2)=ior(int3x(fix+vix5(2)),ior(int3y(fiy+viy5(2)),int3z(fiz+viz5(2))))
            vixyz(3)=ior(int3x(fix+vix5(3)),ior(int3y(fiy+viy5(1)),int3z(fiz+viz5(1))))
            vixyz(4)=ior(int3x(fix+vix5(4)),ior(int3y(fiy+viy5(2)),int3z(fiz+viz5(2))))
            vixyz(5)=ior(int3x(fix+vix5(5)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(6)=ior(int3x(fix+vix5(6)),ior(int3y(fiy+viy5(6)),int3z(fiz+viz5(6))))
            vixyz(7)=ior(int3x(fix+vix5(7)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(6)),int3z(fiz+viz5(6))))
          case(2)          
            vixyz(2)=ior(fint3x(fix+vix5(2)),ior(fint3y(fiy+viy5(2)),fint3z(fiz+viz5(2))))
            vixyz(3)=ior(fint3x(fix+vix5(3)),ior(fint3y(fiy+viy5(1)),fint3z(fiz+viz5(1))))
            vixyz(4)=ior(fint3x(fix+vix5(4)),ior(fint3y(fiy+viy5(2)),fint3z(fiz+viz5(2))))
            vixyz(5)=ior(fint3x(fix+vix5(5)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(6)=ior(fint3x(fix+vix5(6)),ior(fint3y(fiy+viy5(6)),fint3z(fiz+viz5(6))))
            vixyz(7)=ior(fint3x(fix+vix5(7)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(6)),fint3z(fiz+viz5(6))))  
          end select        
          vcboxaddr(  8) = vixyz(2) + 3
          vcboxaddr(  9) = vixyz(2) + 4
          vcboxaddr( 10) = vixyz(2) + 5
          vcboxaddr( 11) = vixyz(2) + 6
          vcboxaddr( 12) = vixyz(2) + 7
          vcboxaddr( 13) = vixyz(3)
          vcboxaddr( 14) = vixyz(3) + 2
          vcboxaddr( 15) = vixyz(4)
          vcboxaddr( 16) = vixyz(4) + 2
          vcboxaddr( 17) = vixyz(4) + 4
          vcboxaddr( 18) = vixyz(4) + 6
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + 1
          vcboxaddr( 21) = vixyz(6)
          vcboxaddr( 22) = vixyz(6) + 1
          vcboxaddr( 23) = vixyz(6) + 4
          vcboxaddr( 24) = vixyz(6) + 5
          vcboxaddr( 25) = vixyz(7)
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + 4
	            
          fwdstart =  9
	end select  
        case(2)

	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(6)),int3z(fiz+pviz(6))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(6)),int3z(fiz+pviz(6))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(6)),fint3z(fiz+pviz(6))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(6)),fint3z(fiz+pviz(6))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(1) + poffset(1)
          vcboxaddr(  3) = vixyz(1) + poffset(2)
          vcboxaddr(  4) = vixyz(1) + poffset(3)
          vcboxaddr(  5) = vixyz(2)
          vcboxaddr(  6) = vixyz(2) + poffset(1)
          vcboxaddr(  7) = vixyz(2) + poffset(2)
          vcboxaddr(  8) = vixyz(2) + poffset(3)
          vcboxaddr(  9) = vixyz(2) + poffset(4)
          vcboxaddr( 10) = vixyz(2) + poffset(5)
          vcboxaddr( 11) = vixyz(2) + poffset(6)
          vcboxaddr( 12) = vixyz(2) + poffset(7)
          vcboxaddr( 13) = vixyz(3)
          vcboxaddr( 14) = vixyz(3) + poffset(2)
          vcboxaddr( 15) = vixyz(4)
          vcboxaddr( 16) = vixyz(4) + poffset(2)
          vcboxaddr( 17) = vixyz(4) + poffset(4)
          vcboxaddr( 18) = vixyz(4) + poffset(6)
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + poffset(1)
          vcboxaddr( 21) = vixyz(6)
          vcboxaddr( 22) = vixyz(6) + poffset(1)
          vcboxaddr( 23) = vixyz(6) + poffset(4)
          vcboxaddr( 24) = vixyz(6) + poffset(5)
          vcboxaddr( 25) = vixyz(7)
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + poffset(4)
	  
          dist =288+child
          fwdstart =  9
      end select
      end subroutine caladdr37
      subroutine caladdr38(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       31831
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2,-2, 2, 2,-2,-2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2, 0,-2, 0,-2, 0,-2, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2,-2,-2,-2, 2, 2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1,-1, 1, 1,-1,-1, 1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1, 0,-1, 0,-1, 0,-1, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1,-1,-1,-1, 1, 1, 1, 1/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)          
            vixyz(2)=ior(int3x(fix+vix(2)),ior(int3y(fiy+viy(2)),int3z(fiz+viz(2))))
            vixyz(3)=ior(int3x(fix+vix(3)),ior(int3y(fiy+viy(1)),int3z(fiz+viz(1))))
            vixyz(4)=ior(int3x(fix+vix(4)),ior(int3y(fiy+viy(2)),int3z(fiz+viz(2))))
            vixyz(5)=ior(int3x(fix+vix(5)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(6)=ior(int3x(fix+vix(6)),ior(int3y(fiy+viy(6)),int3z(fiz+viz(6))))
            vixyz(7)=ior(int3x(fix+vix(7)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(6)),int3z(fiz+viz(6))))
          case(2)          
            vixyz(2)=ior(fint3x(fix+vix(2)),ior(fint3y(fiy+viy(2)),fint3z(fiz+viz(2))))
            vixyz(3)=ior(fint3x(fix+vix(3)),ior(fint3y(fiy+viy(1)),fint3z(fiz+viz(1))))
            vixyz(4)=ior(fint3x(fix+vix(4)),ior(fint3y(fiy+viy(2)),fint3z(fiz+viz(2))))
            vixyz(5)=ior(fint3x(fix+vix(5)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(6)=ior(fint3x(fix+vix(6)),ior(fint3y(fiy+viy(6)),fint3z(fiz+viz(6))))
            vixyz(7)=ior(fint3x(fix+vix(7)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(6)),fint3z(fiz+viz(6))))    
          end select      
          vcboxaddr( 10) = vixyz(2) + 40
          vcboxaddr( 11) = vixyz(2) + 48
          vcboxaddr( 12) = vixyz(2) + 56
          vcboxaddr( 13) = vixyz(3)
          vcboxaddr( 14) = vixyz(3) + 32
          vcboxaddr( 15) = vixyz(4)
          vcboxaddr( 16) = vixyz(4) +16
          vcboxaddr( 17) = vixyz(4) + 32
          vcboxaddr( 18) = vixyz(4) + 48
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + 8
          vcboxaddr( 21) = vixyz(6)
          vcboxaddr( 22) = vixyz(6) + 8
          vcboxaddr( 23) = vixyz(6) +16
          vcboxaddr( 24) = vixyz(6) +24
          vcboxaddr( 25) = vixyz(7)
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) +16
          dist =296+child
          fwdstart = 11
        case(2)
          select case(int3exist)
          case(1)          
            vixyz(2)=ior(int3x(fix+vix5(2)),ior(int3y(fiy+viy5(2)),int3z(fiz+viz5(2))))
            vixyz(3)=ior(int3x(fix+vix5(3)),ior(int3y(fiy+viy5(1)),int3z(fiz+viz5(1))))
            vixyz(4)=ior(int3x(fix+vix5(4)),ior(int3y(fiy+viy5(2)),int3z(fiz+viz5(2))))
            vixyz(5)=ior(int3x(fix+vix5(5)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(6)=ior(int3x(fix+vix5(6)),ior(int3y(fiy+viy5(6)),int3z(fiz+viz5(6))))
            vixyz(7)=ior(int3x(fix+vix5(7)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(6)),int3z(fiz+viz5(6))))
          case(2)          
            vixyz(2)=ior(fint3x(fix+vix5(2)),ior(fint3y(fiy+viy5(2)),fint3z(fiz+viz5(2))))
            vixyz(3)=ior(fint3x(fix+vix5(3)),ior(fint3y(fiy+viy5(1)),fint3z(fiz+viz5(1))))
            vixyz(4)=ior(fint3x(fix+vix5(4)),ior(fint3y(fiy+viy5(2)),fint3z(fiz+viz5(2))))
            vixyz(5)=ior(fint3x(fix+vix5(5)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(6)=ior(fint3x(fix+vix5(6)),ior(fint3y(fiy+viy5(6)),fint3z(fiz+viz5(6))))
            vixyz(7)=ior(fint3x(fix+vix5(7)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(6)),fint3z(fiz+viz5(6))))
          end select
          vcboxaddr( 10) = vixyz(2) + 5
          vcboxaddr( 11) = vixyz(2) + 6
          vcboxaddr( 12) = vixyz(2) + 7
          vcboxaddr( 13) = vixyz(3)
          vcboxaddr( 14) = vixyz(3) + 4
          vcboxaddr( 15) = vixyz(4)
          vcboxaddr( 16) = vixyz(4) + 2
          vcboxaddr( 17) = vixyz(4) + 4
          vcboxaddr( 18) = vixyz(4) + 6
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + 1
          vcboxaddr( 21) = vixyz(6)
          vcboxaddr( 22) = vixyz(6) + 1
          vcboxaddr( 23) = vixyz(6) + 2
          vcboxaddr( 24) = vixyz(6) + 3
          vcboxaddr( 25) = vixyz(7)
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + 2
         
          fwdstart = 11
	end select  
        case(2)
	
	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(6)),int3z(fiz+pviz(6))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(6)),int3z(fiz+pviz(6))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(6)),fint3z(fiz+pviz(6))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(6)),fint3z(fiz+pviz(6))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(1) + poffset(1)
          vcboxaddr(  3) = vixyz(1) + poffset(4)
          vcboxaddr(  4) = vixyz(1) + poffset(5)
          vcboxaddr(  5) = vixyz(2)
          vcboxaddr(  6) = vixyz(2) + poffset(1)
          vcboxaddr(  7) = vixyz(2) + poffset(2)
          vcboxaddr(  8) = vixyz(2) + poffset(3)
          vcboxaddr(  9) = vixyz(2) + poffset(4)
          vcboxaddr( 10) = vixyz(2) + poffset(5)
          vcboxaddr( 11) = vixyz(2) + poffset(6)
          vcboxaddr( 12) = vixyz(2) + poffset(7)
          vcboxaddr( 13) = vixyz(3)
          vcboxaddr( 14) = vixyz(3) + poffset(4)
          vcboxaddr( 15) = vixyz(4)
          vcboxaddr( 16) = vixyz(4) + poffset(2)
          vcboxaddr( 17) = vixyz(4) + poffset(4)
          vcboxaddr( 18) = vixyz(4) + poffset(6)
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + poffset(1)
          vcboxaddr( 21) = vixyz(6)
          vcboxaddr( 22) = vixyz(6) + poffset(1)
          vcboxaddr( 23) = vixyz(6) + poffset(2)
          vcboxaddr( 24) = vixyz(6) + poffset(3)
          vcboxaddr( 25) = vixyz(7)
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + poffset(2)
	  
          dist =296+child
          fwdstart = 11
      end select
      end subroutine caladdr38
      subroutine caladdr39(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       30955
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2,-2,-2,-2, 2, 2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2,-2, 0, 0,-2,-2, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2, 2,-2, 2,-2, 2,-2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1,-1,-1,-1, 1, 1, 1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1,-1, 0, 0,-1,-1, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1, 1,-1, 1,-1, 1,-1, 1/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)
            vixyz(3)=ior(int3x(fix+vix(3)),ior(int3y(fiy+viy(3)),int3z(fiz+viz(3))))
            vixyz(4)=ior(int3x(fix+vix(4)),ior(int3y(fiy+viy(4)),int3z(fiz+viz(4))))
            vixyz(5)=ior(int3x(fix+vix(5)),ior(int3y(fiy+viy(1)),int3z(fiz+viz(1))))
            vixyz(6)=ior(int3x(fix+vix(6)),ior(int3y(fiy+viy(2)),int3z(fiz+viz(2))))
            vixyz(7)=ior(int3x(fix+vix(7)),ior(int3y(fiy+viy(3)),int3z(fiz+viz(3))))
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(4)),int3z(fiz+viz(4))))
          case(2)
            vixyz(3)=ior(fint3x(fix+vix(3)),ior(fint3y(fiy+viy(3)),fint3z(fiz+viz(3))))
            vixyz(4)=ior(fint3x(fix+vix(4)),ior(fint3y(fiy+viy(4)),fint3z(fiz+viz(4))))
            vixyz(5)=ior(fint3x(fix+vix(5)),ior(fint3y(fiy+viy(1)),fint3z(fiz+viz(1))))
            vixyz(6)=ior(fint3x(fix+vix(6)),ior(fint3y(fiy+viy(2)),fint3z(fiz+viz(2))))
            vixyz(7)=ior(fint3x(fix+vix(7)),ior(fint3y(fiy+viy(3)),fint3z(fiz+viz(3))))
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(4)),fint3z(fiz+viz(4))))    
          end select      
          vcboxaddr( 12) = vixyz(3) + 40
          vcboxaddr( 13) = vixyz(3) + 48
          vcboxaddr( 14) = vixyz(3) + 56
          vcboxaddr( 15) = vixyz(4)
          vcboxaddr( 16) = vixyz(4) + 8
          vcboxaddr( 17) = vixyz(4) +16
          vcboxaddr( 18) = vixyz(4) +24
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + 32
          vcboxaddr( 21) = vixyz(6)
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) +16
          vcboxaddr( 24) = vixyz(7) + 32
          vcboxaddr( 25) = vixyz(7) + 48
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) +16
          dist =304+child
          fwdstart = 13
	 case(2)
          select case(int3exist)
          case(1)
            vixyz(3)=ior(int3x(fix+vix5(3)),ior(int3y(fiy+viy5(3)),int3z(fiz+viz5(3))))
            vixyz(4)=ior(int3x(fix+vix5(4)),ior(int3y(fiy+viy5(4)),int3z(fiz+viz5(4))))
            vixyz(5)=ior(int3x(fix+vix5(5)),ior(int3y(fiy+viy5(1)),int3z(fiz+viz5(1))))
            vixyz(6)=ior(int3x(fix+vix5(6)),ior(int3y(fiy+viy5(2)),int3z(fiz+viz5(2))))
            vixyz(7)=ior(int3x(fix+vix5(7)),ior(int3y(fiy+viy5(3)),int3z(fiz+viz5(3))))
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(4)),int3z(fiz+viz5(4))))
          case(2)
            vixyz(3)=ior(fint3x(fix+vix5(3)),ior(fint3y(fiy+viy5(3)),fint3z(fiz+viz5(3))))
            vixyz(4)=ior(fint3x(fix+vix5(4)),ior(fint3y(fiy+viy5(4)),fint3z(fiz+viz5(4))))
            vixyz(5)=ior(fint3x(fix+vix5(5)),ior(fint3y(fiy+viy5(1)),fint3z(fiz+viz5(1))))
            vixyz(6)=ior(fint3x(fix+vix5(6)),ior(fint3y(fiy+viy5(2)),fint3z(fiz+viz5(2))))
            vixyz(7)=ior(fint3x(fix+vix5(7)),ior(fint3y(fiy+viy5(3)),fint3z(fiz+viz5(3))))
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(4)),fint3z(fiz+viz5(4))))
          end select
          vcboxaddr( 12) = vixyz(3) + 5
          vcboxaddr( 13) = vixyz(3) + 6
          vcboxaddr( 14) = vixyz(3) + 7
          vcboxaddr( 15) = vixyz(4)
          vcboxaddr( 16) = vixyz(4) + 1
          vcboxaddr( 17) = vixyz(4) + 2
          vcboxaddr( 18) = vixyz(4) + 3
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + 4
          vcboxaddr( 21) = vixyz(6)
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + 2
          vcboxaddr( 24) = vixyz(7) + 4
          vcboxaddr( 25) = vixyz(7) + 6
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + 2
          
          fwdstart = 13
	end select  
        case(2)

	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(4)),int3z(fiz+pviz(4))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(4)),int3z(fiz+pviz(4))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(4)),fint3z(fiz+pviz(4))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(4)),fint3z(fiz+pviz(4))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(1) + poffset(1)
          vcboxaddr(  3) = vixyz(1) + poffset(4)
          vcboxaddr(  4) = vixyz(1) + poffset(5)
          vcboxaddr(  5) = vixyz(2)
          vcboxaddr(  6) = vixyz(2) + poffset(1)
          vcboxaddr(  7) = vixyz(3)
          vcboxaddr(  8) = vixyz(3) + poffset(1)
          vcboxaddr(  9) = vixyz(3) + poffset(2)
          vcboxaddr( 10) = vixyz(3) + poffset(3)
          vcboxaddr( 11) = vixyz(3) + poffset(4)
          vcboxaddr( 12) = vixyz(3) + poffset(5)
          vcboxaddr( 13) = vixyz(3) + poffset(6)
          vcboxaddr( 14) = vixyz(3) + poffset(7)
          vcboxaddr( 15) = vixyz(4)
          vcboxaddr( 16) = vixyz(4) + poffset(1)
          vcboxaddr( 17) = vixyz(4) + poffset(2)
          vcboxaddr( 18) = vixyz(4) + poffset(3)
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + poffset(4)
          vcboxaddr( 21) = vixyz(6)
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + poffset(2)
          vcboxaddr( 24) = vixyz(7) + poffset(4)
          vcboxaddr( 25) = vixyz(7) + poffset(6)
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + poffset(2)
	  
          dist =304+child
          fwdstart = 13
      end select
      end subroutine caladdr39
      subroutine caladdr40(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       29695
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2,-2, 0, 0,-2,-2, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2,-2,-2,-2, 2, 2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2, 2,-2, 2,-2, 2,-2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1,-1, 0, 0,-1,-1, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1,-1,-1,-1, 1, 1, 1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1, 1,-1, 1,-1, 1,-1, 1/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)
            vixyz(3)=ior(int3x(fix+vix(3)),ior(int3y(fiy+viy(1)),int3z(fiz+viz(1))))
            vixyz(4)=ior(int3x(fix+vix(4)),ior(int3y(fiy+viy(2)),int3z(fiz+viz(2))))
            vixyz(5)=ior(int3x(fix+vix(5)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(6)=ior(int3x(fix+vix(6)),ior(int3y(fiy+viy(6)),int3z(fiz+viz(6))))
            vixyz(7)=ior(int3x(fix+vix(7)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(6)),int3z(fiz+viz(6))))
          case(2)
            vixyz(3)=ior(fint3x(fix+vix(3)),ior(fint3y(fiy+viy(1)),fint3z(fiz+viz(1))))
            vixyz(4)=ior(fint3x(fix+vix(4)),ior(fint3y(fiy+viy(2)),fint3z(fiz+viz(2))))
            vixyz(5)=ior(fint3x(fix+vix(5)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(6)=ior(fint3x(fix+vix(6)),ior(fint3y(fiy+viy(6)),fint3z(fiz+viz(6))))
            vixyz(7)=ior(fint3x(fix+vix(7)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(6)),fint3z(fiz+viz(6))))
          end select
          vcboxaddr( 13) = vixyz(3) + 48
          vcboxaddr( 14) = vixyz(3) + 56
          vcboxaddr( 15) = vixyz(4)
          vcboxaddr( 16) = vixyz(4) + 8
          vcboxaddr( 17) = vixyz(4) +16
          vcboxaddr( 18) = vixyz(4) +24
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + 32
          vcboxaddr( 21) = vixyz(6)
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + 8
          vcboxaddr( 24) = vixyz(7) + 32
          vcboxaddr( 25) = vixyz(7) + 40
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + 8
          dist =312+child
          fwdstart = 14
	case(2)
          select case(int3exist)
          case(1)
            vixyz(3)=ior(int3x(fix+vix5(3)),ior(int3y(fiy+viy5(1)),int3z(fiz+viz5(1))))
            vixyz(4)=ior(int3x(fix+vix5(4)),ior(int3y(fiy+viy5(2)),int3z(fiz+viz5(2))))
            vixyz(5)=ior(int3x(fix+vix5(5)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(6)=ior(int3x(fix+vix5(6)),ior(int3y(fiy+viy5(6)),int3z(fiz+viz5(6))))
            vixyz(7)=ior(int3x(fix+vix5(7)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(6)),int3z(fiz+viz5(6))))
          case(2)
            vixyz(3)=ior(fint3x(fix+vix5(3)),ior(fint3y(fiy+viy5(1)),fint3z(fiz+viz5(1))))
            vixyz(4)=ior(fint3x(fix+vix5(4)),ior(fint3y(fiy+viy5(2)),fint3z(fiz+viz5(2))))
            vixyz(5)=ior(fint3x(fix+vix5(5)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(6)=ior(fint3x(fix+vix5(6)),ior(fint3y(fiy+viy5(6)),fint3z(fiz+viz5(6))))
            vixyz(7)=ior(fint3x(fix+vix5(7)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(6)),fint3z(fiz+viz5(6))))
          end select
          vcboxaddr( 13) = vixyz(3) + 6
          vcboxaddr( 14) = vixyz(3) + 7
          vcboxaddr( 15) = vixyz(4)
          vcboxaddr( 16) = vixyz(4) + 1
          vcboxaddr( 17) = vixyz(4) + 2
          vcboxaddr( 18) = vixyz(4) + 3
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + 4
          vcboxaddr( 21) = vixyz(6)
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + 1
          vcboxaddr( 24) = vixyz(7) + 4
          vcboxaddr( 25) = vixyz(7) + 5
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + 1
          
          fwdstart = 14
	end select  
        case(2)
	
	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(6)),int3z(fiz+pviz(6))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(6)),int3z(fiz+pviz(6))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(6)),fint3z(fiz+pviz(6))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(6)),fint3z(fiz+pviz(6))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(1) + poffset(2)
          vcboxaddr(  3) = vixyz(1) + poffset(4)
          vcboxaddr(  4) = vixyz(1) + poffset(6)
          vcboxaddr(  5) = vixyz(2)
          vcboxaddr(  6) = vixyz(2) + poffset(2)
          vcboxaddr(  7) = vixyz(3)
          vcboxaddr(  8) = vixyz(3) + poffset(1)
          vcboxaddr(  9) = vixyz(3) + poffset(2)
          vcboxaddr( 10) = vixyz(3) + poffset(3)
          vcboxaddr( 11) = vixyz(3) + poffset(4)
          vcboxaddr( 12) = vixyz(3) + poffset(5)
          vcboxaddr( 13) = vixyz(3) + poffset(6)
          vcboxaddr( 14) = vixyz(3) + poffset(7)
          vcboxaddr( 15) = vixyz(4)
          vcboxaddr( 16) = vixyz(4) + poffset(1)
          vcboxaddr( 17) = vixyz(4) + poffset(2)
          vcboxaddr( 18) = vixyz(4) + poffset(3)
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + poffset(4)
          vcboxaddr( 21) = vixyz(6)
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + poffset(1)
          vcboxaddr( 24) = vixyz(7) + poffset(4)
          vcboxaddr( 25) = vixyz(7) + poffset(5)
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + poffset(1)
	  
          dist =312+child
          fwdstart = 14
      end select
      end subroutine caladdr40
      subroutine caladdr41(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       36355
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2,-2,-2,-2, 2, 2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2,-2, 2, 2,-2,-2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2, 2,-2, 2,-2, 2,-2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1,-1,-1,-1, 1, 1, 1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1,-1, 1, 1,-1,-1, 1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1, 1,-1, 1,-1, 1,-1, 1/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+vix(1)),ior(int3y(fiy+viy(1)),int3z(fiz+viz(1))))
            vixyz(2)=ior(int3x(fix+vix(2)),ior(int3y(fiy+viy(2)),int3z(fiz+viz(2))))
            vixyz(3)=ior(int3x(fix+vix(3)),ior(int3y(fiy+viy(3)),int3z(fiz+viz(3))))
            vixyz(4)=ior(int3x(fix+vix(4)),ior(int3y(fiy+viy(4)),int3z(fiz+viz(4))))
            vixyz(5)=ior(int3x(fix+vix(5)),ior(int3y(fiy+viy(1)),int3z(fiz+viz(1))))
            vixyz(6)=ior(int3x(fix+vix(6)),ior(int3y(fiy+viy(2)),int3z(fiz+viz(2))))
            vixyz(7)=ior(int3x(fix+vix(7)),ior(int3y(fiy+viy(3)),int3z(fiz+viz(3))))
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(4)),int3z(fiz+viz(4))))
          case(2)
            vixyz(1)=ior(fint3x(fix+vix(1)),ior(fint3y(fiy+viy(1)),fint3z(fiz+viz(1))))
            vixyz(2)=ior(fint3x(fix+vix(2)),ior(fint3y(fiy+viy(2)),fint3z(fiz+viz(2))))
            vixyz(3)=ior(fint3x(fix+vix(3)),ior(fint3y(fiy+viy(3)),fint3z(fiz+viz(3))))
            vixyz(4)=ior(fint3x(fix+vix(4)),ior(fint3y(fiy+viy(4)),fint3z(fiz+viz(4))))
            vixyz(5)=ior(fint3x(fix+vix(5)),ior(fint3y(fiy+viy(1)),fint3z(fiz+viz(1))))
            vixyz(6)=ior(fint3x(fix+vix(6)),ior(fint3y(fiy+viy(2)),fint3z(fiz+viz(2))))
            vixyz(7)=ior(fint3x(fix+vix(7)),ior(fint3y(fiy+viy(3)),fint3z(fiz+viz(3))))
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(4)),fint3z(fiz+viz(4))))
          end select
          vcboxaddr(  8) = vixyz(1) + 56
          vcboxaddr(  9) = vixyz(2)
          vcboxaddr( 10) = vixyz(2) + 8
          vcboxaddr( 11) = vixyz(2) +16
          vcboxaddr( 12) = vixyz(2) +24
          vcboxaddr( 13) = vixyz(3)
          vcboxaddr( 14) = vixyz(3) + 8
          vcboxaddr( 15) = vixyz(3) + 32
          vcboxaddr( 16) = vixyz(3) + 40
          vcboxaddr( 17) = vixyz(4)
          vcboxaddr( 18) = vixyz(4) + 8
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) +16
          vcboxaddr( 21) = vixyz(5) + 32
          vcboxaddr( 22) = vixyz(5) + 48
          vcboxaddr( 23) = vixyz(6)
          vcboxaddr( 24) = vixyz(6) +16
          vcboxaddr( 25) = vixyz(7)
          vcboxaddr( 26) = vixyz(7) + 32
          vcboxaddr( 27) = vixyz(8)
          dist =320+child
          fwdstart =  9
	case(2)
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+vix5(1)),ior(int3y(fiy+viy5(1)),int3z(fiz+viz5(1))))
            vixyz(2)=ior(int3x(fix+vix5(2)),ior(int3y(fiy+viy5(2)),int3z(fiz+viz5(2))))
            vixyz(3)=ior(int3x(fix+vix5(3)),ior(int3y(fiy+viy5(3)),int3z(fiz+viz5(3))))
            vixyz(4)=ior(int3x(fix+vix5(4)),ior(int3y(fiy+viy5(4)),int3z(fiz+viz5(4))))
            vixyz(5)=ior(int3x(fix+vix5(5)),ior(int3y(fiy+viy5(1)),int3z(fiz+viz5(1))))
            vixyz(6)=ior(int3x(fix+vix5(6)),ior(int3y(fiy+viy5(2)),int3z(fiz+viz5(2))))
            vixyz(7)=ior(int3x(fix+vix5(7)),ior(int3y(fiy+viy5(3)),int3z(fiz+viz5(3))))
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(4)),int3z(fiz+viz5(4))))
          case(2)
            vixyz(1)=ior(fint3x(fix+vix5(1)),ior(fint3y(fiy+viy5(1)),fint3z(fiz+viz5(1))))
            vixyz(2)=ior(fint3x(fix+vix5(2)),ior(fint3y(fiy+viy5(2)),fint3z(fiz+viz5(2))))
            vixyz(3)=ior(fint3x(fix+vix5(3)),ior(fint3y(fiy+viy5(3)),fint3z(fiz+viz5(3))))
            vixyz(4)=ior(fint3x(fix+vix5(4)),ior(fint3y(fiy+viy5(4)),fint3z(fiz+viz5(4))))
            vixyz(5)=ior(fint3x(fix+vix5(5)),ior(fint3y(fiy+viy5(1)),fint3z(fiz+viz5(1))))
            vixyz(6)=ior(fint3x(fix+vix5(6)),ior(fint3y(fiy+viy5(2)),fint3z(fiz+viz5(2))))
            vixyz(7)=ior(fint3x(fix+vix5(7)),ior(fint3y(fiy+viy5(3)),fint3z(fiz+viz5(3))))
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(4)),fint3z(fiz+viz5(4))))
          end select          
          vcboxaddr(  8) = vixyz(1) + 7
          vcboxaddr(  9) = vixyz(2)
          vcboxaddr( 10) = vixyz(2) + 1
          vcboxaddr( 11) = vixyz(2) + 2
          vcboxaddr( 12) = vixyz(2) + 3
          vcboxaddr( 13) = vixyz(3)
          vcboxaddr( 14) = vixyz(3) + 1
          vcboxaddr( 15) = vixyz(3) + 4
          vcboxaddr( 16) = vixyz(3) + 5
          vcboxaddr( 17) = vixyz(4)
          vcboxaddr( 18) = vixyz(4) + 1
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + 2
          vcboxaddr( 21) = vixyz(5) + 4
          vcboxaddr( 22) = vixyz(5) + 6
          vcboxaddr( 23) = vixyz(6)
          vcboxaddr( 24) = vixyz(6) + 2
          vcboxaddr( 25) = vixyz(7)
          vcboxaddr( 26) = vixyz(7) + 4
          vcboxaddr( 27) = vixyz(8)
          
          fwdstart =  9
	end select  
	case(2)

	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(4)),int3z(fiz+pviz(4))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(4)),int3z(fiz+pviz(4))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(4)),fint3z(fiz+pviz(4))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(4)),fint3z(fiz+pviz(4))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(1) + poffset(1)
          vcboxaddr(  3) = vixyz(1) + poffset(2)
          vcboxaddr(  4) = vixyz(1) + poffset(3)
          vcboxaddr(  5) = vixyz(1) + poffset(4)
          vcboxaddr(  6) = vixyz(1) + poffset(5)
          vcboxaddr(  7) = vixyz(1) + poffset(6)
	  vcboxaddr(  8) = vixyz(1) + poffset(7)
          vcboxaddr(  9) = vixyz(2)
          vcboxaddr( 10) = vixyz(2) + poffset(1)
          vcboxaddr( 11) = vixyz(2) + poffset(2)
          vcboxaddr( 12) = vixyz(2) + poffset(3)
          vcboxaddr( 13) = vixyz(3)
          vcboxaddr( 14) = vixyz(3) + poffset(1)
          vcboxaddr( 15) = vixyz(3) + poffset(4)
          vcboxaddr( 16) = vixyz(3) + poffset(5)
          vcboxaddr( 17) = vixyz(4)
          vcboxaddr( 18) = vixyz(4) + poffset(1)
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + poffset(2)
          vcboxaddr( 21) = vixyz(5) + poffset(4)
          vcboxaddr( 22) = vixyz(5) + poffset(6)
          vcboxaddr( 23) = vixyz(6)
          vcboxaddr( 24) = vixyz(6) + poffset(2)
          vcboxaddr( 25) = vixyz(7)
          vcboxaddr( 26) = vixyz(7) + poffset(4)
          vcboxaddr( 27) = vixyz(8)
	  
          dist =320+child
          fwdstart =  9
      end select
      end subroutine caladdr41
      subroutine caladdr42(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       28393
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2,-2, 0, 0,-2,-2, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2, 2,-2, 2,-2, 2,-2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2,-2,-2,-2, 2, 2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1,-1, 0, 0,-1,-1, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1, 1,-1, 1,-1, 1,-1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1,-1,-1,-1, 1, 1, 1, 1/)      
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)
            vixyz(3)=ior(int3x(fix+vix(3)),ior(int3y(fiy+viy(1)),int3z(fiz+viz(1))))
            vixyz(4)=ior(int3x(fix+vix(4)),ior(int3y(fiy+viy(2)),int3z(fiz+viz(2))))
            vixyz(5)=ior(int3x(fix+vix(5)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(6)=ior(int3x(fix+vix(6)),ior(int3y(fiy+viy(6)),int3z(fiz+viz(6))))
            vixyz(7)=ior(int3x(fix+vix(7)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(6)),int3z(fiz+viz(6))))
          case(2)
            vixyz(3)=ior(fint3x(fix+vix(3)),ior(fint3y(fiy+viy(1)),fint3z(fiz+viz(1))))
            vixyz(4)=ior(fint3x(fix+vix(4)),ior(fint3y(fiy+viy(2)),fint3z(fiz+viz(2))))
            vixyz(5)=ior(fint3x(fix+vix(5)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(6)=ior(fint3x(fix+vix(6)),ior(fint3y(fiy+viy(6)),fint3z(fiz+viz(6))))
            vixyz(7)=ior(fint3x(fix+vix(7)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(6)),fint3z(fiz+viz(6))))
          end select
          vcboxaddr( 13) = vixyz(3) + 48
          vcboxaddr( 14) = vixyz(3) + 56
          vcboxaddr( 15) = vixyz(4)
          vcboxaddr( 16) = vixyz(4) + 8
          vcboxaddr( 17) = vixyz(4) + 32
          vcboxaddr( 18) = vixyz(4) + 40
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) +16
          vcboxaddr( 21) = vixyz(6)
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + 8
          vcboxaddr( 24) = vixyz(7) +16
          vcboxaddr( 25) = vixyz(7) +24
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + 8
          dist =328+child
          fwdstart = 14
	case(2)
          select case(int3exist)
          case(1)
            vixyz(3)=ior(int3x(fix+vix5(3)),ior(int3y(fiy+viy5(1)),int3z(fiz+viz5(1))))
            vixyz(4)=ior(int3x(fix+vix5(4)),ior(int3y(fiy+viy5(2)),int3z(fiz+viz5(2))))
            vixyz(5)=ior(int3x(fix+vix5(5)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(6)=ior(int3x(fix+vix5(6)),ior(int3y(fiy+viy5(6)),int3z(fiz+viz5(6))))
            vixyz(7)=ior(int3x(fix+vix5(7)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(6)),int3z(fiz+viz5(6))))
          case(2)
            vixyz(3)=ior(fint3x(fix+vix5(3)),ior(fint3y(fiy+viy5(1)),fint3z(fiz+viz5(1))))
            vixyz(4)=ior(fint3x(fix+vix5(4)),ior(fint3y(fiy+viy5(2)),fint3z(fiz+viz5(2))))
            vixyz(5)=ior(fint3x(fix+vix5(5)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(6)=ior(fint3x(fix+vix5(6)),ior(fint3y(fiy+viy5(6)),fint3z(fiz+viz5(6))))
            vixyz(7)=ior(fint3x(fix+vix5(7)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(6)),fint3z(fiz+viz5(6))))	
          end select  
          vcboxaddr( 13) = vixyz(3) + 6
          vcboxaddr( 14) = vixyz(3) + 7
          vcboxaddr( 15) = vixyz(4)
          vcboxaddr( 16) = vixyz(4) + 1
          vcboxaddr( 17) = vixyz(4) + 4
          vcboxaddr( 18) = vixyz(4) + 5
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + 2
          vcboxaddr( 21) = vixyz(6)
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + 1
          vcboxaddr( 24) = vixyz(7) + 2
          vcboxaddr( 25) = vixyz(7) + 3
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + 1
          
          fwdstart = 14
	end select  
        case(2)

	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(6)),int3z(fiz+pviz(6))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(6)),int3z(fiz+pviz(6))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(6)),fint3z(fiz+pviz(6))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(6)),fint3z(fiz+pviz(6))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(1) + poffset(2)
          vcboxaddr(  3) = vixyz(1) + poffset(4)
          vcboxaddr(  4) = vixyz(1) + poffset(6)
          vcboxaddr(  5) = vixyz(2)
          vcboxaddr(  6) = vixyz(2) + poffset(4)
          vcboxaddr(  7) = vixyz(3)
          vcboxaddr(  8) = vixyz(3) + poffset(1)
          vcboxaddr(  9) = vixyz(3) + poffset(2)
          vcboxaddr( 10) = vixyz(3) + poffset(3)
          vcboxaddr( 11) = vixyz(3) + poffset(4)
          vcboxaddr( 12) = vixyz(3) + poffset(5)
          vcboxaddr( 13) = vixyz(3) + poffset(6)
          vcboxaddr( 14) = vixyz(3) + poffset(7)
          vcboxaddr( 15) = vixyz(4)
          vcboxaddr( 16) = vixyz(4) + poffset(1)
          vcboxaddr( 17) = vixyz(4) + poffset(4)
          vcboxaddr( 18) = vixyz(4) + poffset(5)
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + poffset(2)
          vcboxaddr( 21) = vixyz(6)
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + poffset(1)
          vcboxaddr( 24) = vixyz(7) + poffset(2)
          vcboxaddr( 25) = vixyz(7) + poffset(3)
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + poffset(1)

          dist =328+child
          fwdstart = 14
      end select
      end subroutine caladdr42
      subroutine caladdr43(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       28861
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2,-2,-2,-2, 2, 2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2,-2, 0, 0,-2,-2, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2, 0,-2, 0,-2, 0,-2, 0/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1,-1,-1,-1, 1, 1, 1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1,-1, 0, 0,-1,-1, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1, 0,-1, 0,-1, 0,-1, 0/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)          
            vixyz(4)=ior(int3x(fix+vix(4)),ior(int3y(fiy+viy(4)),int3z(fiz+viz(4))))
            vixyz(5)=ior(int3x(fix+vix(5)),ior(int3y(fiy+viy(1)),int3z(fiz+viz(1))))
            vixyz(6)=ior(int3x(fix+vix(6)),ior(int3y(fiy+viy(2)),int3z(fiz+viz(2))))
            vixyz(7)=ior(int3x(fix+vix(7)),ior(int3y(fiy+viy(3)),int3z(fiz+viz(3))))
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(4)),int3z(fiz+viz(4))))
          case(2)          
            vixyz(4)=ior(fint3x(fix+vix(4)),ior(fint3y(fiy+viy(4)),fint3z(fiz+viz(4))))
            vixyz(5)=ior(fint3x(fix+vix(5)),ior(fint3y(fiy+viy(1)),fint3z(fiz+viz(1))))
            vixyz(6)=ior(fint3x(fix+vix(6)),ior(fint3y(fiy+viy(2)),fint3z(fiz+viz(2))))
            vixyz(7)=ior(fint3x(fix+vix(7)),ior(fint3y(fiy+viy(3)),fint3z(fiz+viz(3))))
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(4)),fint3z(fiz+viz(4))))  
          end select       
          vcboxaddr( 12) = vixyz(4) + 8
          vcboxaddr( 13) = vixyz(4) +16
          vcboxaddr( 14) = vixyz(4) +24
          vcboxaddr( 15) = vixyz(4) + 32
          vcboxaddr( 16) = vixyz(4) + 40
          vcboxaddr( 17) = vixyz(4) + 48
          vcboxaddr( 18) = vixyz(4) + 56
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(6)
          vcboxaddr( 21) = vixyz(6) + 32
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) +16
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) +16
          vcboxaddr( 26) = vixyz(8) + 32
          vcboxaddr( 27) = vixyz(8) + 48
          dist =336+child
          fwdstart = 13
	case(2)
          select case(int3exist)
          case(1)          
            vixyz(4)=ior(int3x(fix+vix5(4)),ior(int3y(fiy+viy5(4)),int3z(fiz+viz5(4))))
            vixyz(5)=ior(int3x(fix+vix5(5)),ior(int3y(fiy+viy5(1)),int3z(fiz+viz5(1))))
            vixyz(6)=ior(int3x(fix+vix5(6)),ior(int3y(fiy+viy5(2)),int3z(fiz+viz5(2))))
            vixyz(7)=ior(int3x(fix+vix5(7)),ior(int3y(fiy+viy5(3)),int3z(fiz+viz5(3))))
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(4)),int3z(fiz+viz5(4))))
          case(2)          
            vixyz(4)=ior(fint3x(fix+vix5(4)),ior(fint3y(fiy+viy5(4)),fint3z(fiz+viz5(4))))
            vixyz(5)=ior(fint3x(fix+vix5(5)),ior(fint3y(fiy+viy5(1)),fint3z(fiz+viz5(1))))
            vixyz(6)=ior(fint3x(fix+vix5(6)),ior(fint3y(fiy+viy5(2)),fint3z(fiz+viz5(2))))
            vixyz(7)=ior(fint3x(fix+vix5(7)),ior(fint3y(fiy+viy5(3)),fint3z(fiz+viz5(3))))
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(4)),fint3z(fiz+viz5(4))))
          end select
          vcboxaddr( 12) = vixyz(4) + 1
          vcboxaddr( 13) = vixyz(4) + 2
          vcboxaddr( 14) = vixyz(4) + 3
          vcboxaddr( 15) = vixyz(4) + 4
          vcboxaddr( 16) = vixyz(4) + 5
          vcboxaddr( 17) = vixyz(4) + 6
          vcboxaddr( 18) = vixyz(4) + 7
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(6)
          vcboxaddr( 21) = vixyz(6) + 4
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + 2
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + 2
          vcboxaddr( 26) = vixyz(8) + 4
          vcboxaddr( 27) = vixyz(8) + 6
          
          fwdstart = 13
	end select  
        case(2)

	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(4)),int3z(fiz+pviz(4))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(4)),int3z(fiz+pviz(4))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(4)),fint3z(fiz+pviz(4))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(4)),fint3z(fiz+pviz(4))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(1) + poffset(1)
          vcboxaddr(  3) = vixyz(2)
          vcboxaddr(  4) = vixyz(2) + poffset(1)
          vcboxaddr(  5) = vixyz(2) + poffset(4)
          vcboxaddr(  6) = vixyz(2) + poffset(5)
          vcboxaddr(  7) = vixyz(3)
          vcboxaddr(  8) = vixyz(3) + poffset(1)
          vcboxaddr(  9) = vixyz(3) + poffset(2)
          vcboxaddr( 10) = vixyz(3) + poffset(3)
          vcboxaddr( 11) = vixyz(4)
          vcboxaddr( 12) = vixyz(4) + poffset(1)
          vcboxaddr( 13) = vixyz(4) + poffset(2)
          vcboxaddr( 14) = vixyz(4) + poffset(3)
          vcboxaddr( 15) = vixyz(4) + poffset(4)
          vcboxaddr( 16) = vixyz(4) + poffset(5)
          vcboxaddr( 17) = vixyz(4) + poffset(6)
          vcboxaddr( 18) = vixyz(4) + poffset(7)
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(6)
          vcboxaddr( 21) = vixyz(6) + poffset(4)
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + poffset(2)
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + poffset(2)
          vcboxaddr( 26) = vixyz(8) + poffset(4)
          vcboxaddr( 27) = vixyz(8) + poffset(6)
	  
          dist =336+child
          fwdstart = 13
      end select
      end subroutine caladdr43
      subroutine caladdr44(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       27685
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2,-2, 0, 0,-2,-2, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2,-2,-2,-2, 2, 2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2, 0,-2, 0,-2, 0,-2, 0/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1,-1, 0, 0,-1,-1, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1,-1,-1,-1, 1, 1, 1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1, 0,-1, 0,-1, 0,-1, 0/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)           
            vixyz(4)=ior(int3x(fix+vix(4)),ior(int3y(fiy+viy(2)),int3z(fiz+viz(2))))
            vixyz(5)=ior(int3x(fix+vix(5)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(6)=ior(int3x(fix+vix(6)),ior(int3y(fiy+viy(6)),int3z(fiz+viz(6))))
            vixyz(7)=ior(int3x(fix+vix(7)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(6)),int3z(fiz+viz(6))))
          case(2)           
            vixyz(4)=ior(fint3x(fix+vix(4)),ior(fint3y(fiy+viy(2)),fint3z(fiz+viz(2))))
            vixyz(5)=ior(fint3x(fix+vix(5)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(6)=ior(fint3x(fix+vix(6)),ior(fint3y(fiy+viy(6)),fint3z(fiz+viz(6))))
            vixyz(7)=ior(fint3x(fix+vix(7)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(6)),fint3z(fiz+viz(6))))  
          end select        
          vcboxaddr( 13) = vixyz(4) + 16
          vcboxaddr( 14) = vixyz(4) + 24
          vcboxaddr( 15) = vixyz(4) + 32
          vcboxaddr( 16) = vixyz(4) + 40
          vcboxaddr( 17) = vixyz(4) + 48
          vcboxaddr( 18) = vixyz(4) + 56
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(6)
          vcboxaddr( 21) = vixyz(6) + 32
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + 8
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + 8
          vcboxaddr( 26) = vixyz(8) + 32
          vcboxaddr( 27) = vixyz(8) + 40
          dist =344+child
          fwdstart = 14
	case(2)
          select case(int3exist)
          case(1)          
            vixyz(4)=ior(int3x(fix+vix5(4)),ior(int3y(fiy+viy5(2)),int3z(fiz+viz5(2))))
            vixyz(5)=ior(int3x(fix+vix5(5)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(6)=ior(int3x(fix+vix5(6)),ior(int3y(fiy+viy5(6)),int3z(fiz+viz5(6))))
            vixyz(7)=ior(int3x(fix+vix5(7)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(6)),int3z(fiz+viz5(6))))
          case(2)          
            vixyz(4)=ior(fint3x(fix+vix5(4)),ior(fint3y(fiy+viy5(2)),fint3z(fiz+viz5(2))))
            vixyz(5)=ior(fint3x(fix+vix5(5)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(6)=ior(fint3x(fix+vix5(6)),ior(fint3y(fiy+viy5(6)),fint3z(fiz+viz5(6))))
            vixyz(7)=ior(fint3x(fix+vix5(7)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(6)),fint3z(fiz+viz5(6))))
          end select
          vcboxaddr( 13) = vixyz(4) + 2
          vcboxaddr( 14) = vixyz(4) + 3
          vcboxaddr( 15) = vixyz(4) + 4
          vcboxaddr( 16) = vixyz(4) + 5
          vcboxaddr( 17) = vixyz(4) + 6
          vcboxaddr( 18) = vixyz(4) + 7
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(6)
          vcboxaddr( 21) = vixyz(6) + 4
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + 1
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + 1
          vcboxaddr( 26) = vixyz(8) + 4
          vcboxaddr( 27) = vixyz(8) + 5
      
          fwdstart = 14
	end select  
        case(2)

	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(6)),int3z(fiz+pviz(6))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(6)),int3z(fiz+pviz(6))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(6)),fint3z(fiz+pviz(6))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(6)),fint3z(fiz+pviz(6))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(1) + poffset(2)
          vcboxaddr(  3) = vixyz(2)
          vcboxaddr(  4) = vixyz(2) + poffset(2)
          vcboxaddr(  5) = vixyz(2) + poffset(4)
          vcboxaddr(  6) = vixyz(2) + poffset(6)
          vcboxaddr(  7) = vixyz(3)
          vcboxaddr(  8) = vixyz(3) + poffset(1)
          vcboxaddr(  9) = vixyz(3) + poffset(2)
          vcboxaddr( 10) = vixyz(3) + poffset(3)
          vcboxaddr( 11) = vixyz(4)
          vcboxaddr( 12) = vixyz(4) + poffset(1)
          vcboxaddr( 13) = vixyz(4) + poffset(2)
          vcboxaddr( 14) = vixyz(4) + poffset(3)
          vcboxaddr( 15) = vixyz(4) + poffset(4)
          vcboxaddr( 16) = vixyz(4) + poffset(5)
          vcboxaddr( 17) = vixyz(4) + poffset(6)
          vcboxaddr( 18) = vixyz(4) + poffset(7)
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(6)
          vcboxaddr( 21) = vixyz(6) + poffset(4)
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + poffset(1)
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + poffset(1)
          vcboxaddr( 26) = vixyz(8) + poffset(4)
          vcboxaddr( 27) = vixyz(8) + poffset(5)
	  
          dist =344+child
          fwdstart = 14
      end select
      end subroutine caladdr44
      subroutine caladdr45(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       25621
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2,-2, 0, 0,-2,-2, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2, 0,-2, 0,-2, 0,-2, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2,-2,-2,-2, 2, 2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1,-1, 0, 0,-1,-1, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1, 0,-1, 0,-1, 0,-1, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1,-1,-1,-1, 1, 1, 1, 1/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)
            vixyz(4)=ior(int3x(fix+vix(4)),ior(int3y(fiy+viy(2)),int3z(fiz+viz(2))))
            vixyz(5)=ior(int3x(fix+vix(5)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(6)=ior(int3x(fix+vix(6)),ior(int3y(fiy+viy(6)),int3z(fiz+viz(6))))
            vixyz(7)=ior(int3x(fix+vix(7)),ior(int3y(fiy+viy(5)),int3z(fiz+viz(5))))
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(6)),int3z(fiz+viz(6))))
          case(2)
            vixyz(4)=ior(fint3x(fix+vix(4)),ior(fint3y(fiy+viy(2)),fint3z(fiz+viz(2))))
            vixyz(5)=ior(fint3x(fix+vix(5)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(6)=ior(fint3x(fix+vix(6)),ior(fint3y(fiy+viy(6)),fint3z(fiz+viz(6))))
            vixyz(7)=ior(fint3x(fix+vix(7)),ior(fint3y(fiy+viy(5)),fint3z(fiz+viz(5))))
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(6)),fint3z(fiz+viz(6))))
          end select
          vcboxaddr( 15) = vixyz(4) + 32
          vcboxaddr( 16) = vixyz(4) + 40
          vcboxaddr( 17) = vixyz(4) + 48
          vcboxaddr( 18) = vixyz(4) + 56
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(6)
          vcboxaddr( 21) = vixyz(6) + 16
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + 8
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + 8
          vcboxaddr( 26) = vixyz(8) +16
          vcboxaddr( 27) = vixyz(8) +24
          dist =352+child
          fwdstart =  16
	 case(2)
          select case(int3exist)
          case(1)
            vixyz(4)=ior(int3x(fix+vix5(4)),ior(int3y(fiy+viy5(2)),int3z(fiz+viz5(2))))
            vixyz(5)=ior(int3x(fix+vix5(5)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(6)=ior(int3x(fix+vix5(6)),ior(int3y(fiy+viy5(6)),int3z(fiz+viz5(6))))
            vixyz(7)=ior(int3x(fix+vix5(7)),ior(int3y(fiy+viy5(5)),int3z(fiz+viz5(5))))
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(6)),int3z(fiz+viz5(6))))
          case(2)
            vixyz(4)=ior(fint3x(fix+vix5(4)),ior(fint3y(fiy+viy5(2)),fint3z(fiz+viz5(2))))
            vixyz(5)=ior(fint3x(fix+vix5(5)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(6)=ior(fint3x(fix+vix5(6)),ior(fint3y(fiy+viy5(6)),fint3z(fiz+viz5(6))))
            vixyz(7)=ior(fint3x(fix+vix5(7)),ior(fint3y(fiy+viy5(5)),fint3z(fiz+viz5(5))))
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(6)),fint3z(fiz+viz5(6))))   
          end select       
          vcboxaddr( 15) = vixyz(4) + 4
          vcboxaddr( 16) = vixyz(4) + 5
          vcboxaddr( 17) = vixyz(4) + 6
          vcboxaddr( 18) = vixyz(4) + 7
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(6)
          vcboxaddr( 21) = vixyz(6) + 2
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + 1
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + 1
          vcboxaddr( 26) = vixyz(8) + 2
          vcboxaddr( 27) = vixyz(8) + 3

          fwdstart =  16
	end select
        case(2)

	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(6)),int3z(fiz+pviz(6))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(5)),int3z(fiz+pviz(5))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(6)),int3z(fiz+pviz(6))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(6)),fint3z(fiz+pviz(6))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(5)),fint3z(fiz+pviz(5))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(6)),fint3z(fiz+pviz(6))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(1) + poffset(4)
          vcboxaddr(  3) = vixyz(2)
          vcboxaddr(  4) = vixyz(2) + poffset(2)
          vcboxaddr(  5) = vixyz(2) + poffset(4)
          vcboxaddr(  6) = vixyz(2) + poffset(6)
          vcboxaddr(  7) = vixyz(3)
          vcboxaddr(  8) = vixyz(3) + poffset(1)
          vcboxaddr(  9) = vixyz(3) + poffset(4)
          vcboxaddr( 10) = vixyz(3) + poffset(5)
          vcboxaddr( 11) = vixyz(4)
          vcboxaddr( 12) = vixyz(4) + poffset(1)
          vcboxaddr( 13) = vixyz(4) + poffset(2)
          vcboxaddr( 14) = vixyz(4) + poffset(3)
          vcboxaddr( 15) = vixyz(4) + poffset(4)
          vcboxaddr( 16) = vixyz(4) + poffset(5)
          vcboxaddr( 17) = vixyz(4) + poffset(6)
          vcboxaddr( 18) = vixyz(4) + poffset(7)
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(6)
          vcboxaddr( 21) = vixyz(6) + poffset(2)
          vcboxaddr( 22) = vixyz(7)
          vcboxaddr( 23) = vixyz(7) + poffset(1)
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + poffset(1)
          vcboxaddr( 26) = vixyz(8) + poffset(2)
          vcboxaddr( 27) = vixyz(8) + poffset(3)
	  
          dist =352+child
          fwdstart =  16
      end select
      end subroutine caladdr45
      subroutine caladdr46(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       34381
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2,-2,-2,-2, 2, 2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2,-2, 2, 2,-2,-2, 2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2, 0,-2, 0,-2, 0,-2, 0/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1,-1,-1,-1, 1, 1, 1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1,-1, 1, 1,-1,-1, 1, 1/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1, 0,-1, 0,-1, 0,-1, 0/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)
            vixyz(2)=ior(int3x(fix+vix(2)),ior(int3y(fiy+viy(2)),int3z(fiz+viz(2))))
            vixyz(3)=ior(int3x(fix+vix(3)),ior(int3y(fiy+viy(3)),int3z(fiz+viz(3))))
            vixyz(4)=ior(int3x(fix+vix(4)),ior(int3y(fiy+viy(4)),int3z(fiz+viz(4))))
            vixyz(5)=ior(int3x(fix+vix(5)),ior(int3y(fiy+viy(1)),int3z(fiz+viz(1))))
            vixyz(6)=ior(int3x(fix+vix(6)),ior(int3y(fiy+viy(2)),int3z(fiz+viz(2))))
            vixyz(7)=ior(int3x(fix+vix(7)),ior(int3y(fiy+viy(3)),int3z(fiz+viz(3))))
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(4)),int3z(fiz+viz(4))))
          case(2)
            vixyz(2)=ior(fint3x(fix+vix(2)),ior(fint3y(fiy+viy(2)),fint3z(fiz+viz(2))))
            vixyz(3)=ior(fint3x(fix+vix(3)),ior(fint3y(fiy+viy(3)),fint3z(fiz+viz(3))))
            vixyz(4)=ior(fint3x(fix+vix(4)),ior(fint3y(fiy+viy(4)),fint3z(fiz+viz(4))))
            vixyz(5)=ior(fint3x(fix+vix(5)),ior(fint3y(fiy+viy(1)),fint3z(fiz+viz(1))))
            vixyz(6)=ior(fint3x(fix+vix(6)),ior(fint3y(fiy+viy(2)),fint3z(fiz+viz(2))))
            vixyz(7)=ior(fint3x(fix+vix(7)),ior(fint3y(fiy+viy(3)),fint3z(fiz+viz(3))))
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(4)),fint3z(fiz+viz(4))))
          end select
          vcboxaddr(  8) = vixyz(2) + 24
          vcboxaddr(  9) = vixyz(2) + 32
          vcboxaddr( 10) = vixyz(2) + 40
          vcboxaddr( 11) = vixyz(2) + 48
          vcboxaddr( 12) = vixyz(2) + 56
          vcboxaddr( 13) = vixyz(3)
          vcboxaddr( 14) = vixyz(3) + 8
          vcboxaddr( 15) = vixyz(4)
          vcboxaddr( 16) = vixyz(4) + 8
          vcboxaddr( 17) = vixyz(4) + 32
          vcboxaddr( 18) = vixyz(4) + 40
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) +16
          vcboxaddr( 21) = vixyz(6)
          vcboxaddr( 22) = vixyz(6) +16
          vcboxaddr( 23) = vixyz(6) + 32
          vcboxaddr( 24) = vixyz(6) + 48
          vcboxaddr( 25) = vixyz(7)
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + 32
          dist =360+child
          fwdstart =  9
        case(2)
          select case(int3exist)
          case(1)
            vixyz(2)=ior(int3x(fix+vix5(2)),ior(int3y(fiy+viy5(2)),int3z(fiz+viz5(2))))
            vixyz(3)=ior(int3x(fix+vix5(3)),ior(int3y(fiy+viy5(3)),int3z(fiz+viz5(3))))
            vixyz(4)=ior(int3x(fix+vix5(4)),ior(int3y(fiy+viy5(4)),int3z(fiz+viz5(4))))
            vixyz(5)=ior(int3x(fix+vix5(5)),ior(int3y(fiy+viy5(1)),int3z(fiz+viz5(1))))
            vixyz(6)=ior(int3x(fix+vix5(6)),ior(int3y(fiy+viy5(2)),int3z(fiz+viz5(2))))
            vixyz(7)=ior(int3x(fix+vix5(7)),ior(int3y(fiy+viy5(3)),int3z(fiz+viz5(3))))
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(4)),int3z(fiz+viz5(4))))
          case(2)
            vixyz(2)=ior(fint3x(fix+vix5(2)),ior(fint3y(fiy+viy5(2)),fint3z(fiz+viz5(2))))
            vixyz(3)=ior(fint3x(fix+vix5(3)),ior(fint3y(fiy+viy5(3)),fint3z(fiz+viz5(3))))
            vixyz(4)=ior(fint3x(fix+vix5(4)),ior(fint3y(fiy+viy5(4)),fint3z(fiz+viz5(4))))
            vixyz(5)=ior(fint3x(fix+vix5(5)),ior(fint3y(fiy+viy5(1)),fint3z(fiz+viz5(1))))
            vixyz(6)=ior(fint3x(fix+vix5(6)),ior(fint3y(fiy+viy5(2)),fint3z(fiz+viz5(2))))
            vixyz(7)=ior(fint3x(fix+vix5(7)),ior(fint3y(fiy+viy5(3)),fint3z(fiz+viz5(3))))
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(4)),fint3z(fiz+viz5(4))))
          end select
          vcboxaddr(  8) = vixyz(2) + 3
          vcboxaddr(  9) = vixyz(2) + 4
          vcboxaddr( 10) = vixyz(2) + 5
          vcboxaddr( 11) = vixyz(2) + 6
          vcboxaddr( 12) = vixyz(2) + 7
          vcboxaddr( 13) = vixyz(3)
          vcboxaddr( 14) = vixyz(3) + 1
          vcboxaddr( 15) = vixyz(4)
          vcboxaddr( 16) = vixyz(4) + 1
          vcboxaddr( 17) = vixyz(4) + 4
          vcboxaddr( 18) = vixyz(4) + 5
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + 2
          vcboxaddr( 21) = vixyz(6)
          vcboxaddr( 22) = vixyz(6) + 2
          vcboxaddr( 23) = vixyz(6) + 4
          vcboxaddr( 24) = vixyz(6) + 6
          vcboxaddr( 25) = vixyz(7)
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + 4
          
          fwdstart =  9
	end select  
        case(2)

	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(4)),int3z(fiz+pviz(4))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(4)),int3z(fiz+pviz(4))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(4)),fint3z(fiz+pviz(4))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(4)),fint3z(fiz+pviz(4))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(1) + poffset(1)
          vcboxaddr(  3) = vixyz(1) + poffset(2)
          vcboxaddr(  4) = vixyz(1) + poffset(3)
          vcboxaddr(  5) = vixyz(2)
          vcboxaddr(  6) = vixyz(2) + poffset(1)
          vcboxaddr(  7) = vixyz(2) + poffset(2)
          vcboxaddr(  8) = vixyz(2) + poffset(3)
          vcboxaddr(  9) = vixyz(2) + poffset(4)
          vcboxaddr( 10) = vixyz(2) + poffset(5)
          vcboxaddr( 11) = vixyz(2) + poffset(6)
          vcboxaddr( 12) = vixyz(2) + poffset(7)
          vcboxaddr( 13) = vixyz(3)
          vcboxaddr( 14) = vixyz(3) + poffset(1)
          vcboxaddr( 15) = vixyz(4)
          vcboxaddr( 16) = vixyz(4) + poffset(1)
          vcboxaddr( 17) = vixyz(4) + poffset(4)
          vcboxaddr( 18) = vixyz(4) + poffset(5)
          vcboxaddr( 19) = vixyz(5)
          vcboxaddr( 20) = vixyz(5) + poffset(2)
          vcboxaddr( 21) = vixyz(6)
          vcboxaddr( 22) = vixyz(6) + poffset(2)
          vcboxaddr( 23) = vixyz(6) + poffset(4)
          vcboxaddr( 24) = vixyz(6) + poffset(6)
          vcboxaddr( 25) = vixyz(7)
          vcboxaddr( 26) = vixyz(8)
          vcboxaddr( 27) = vixyz(8) + poffset(4)
	  
          dist =360+child
          fwdstart =  9
      end select
      end subroutine caladdr46
      subroutine caladdr47(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       18373
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2,-2,-2,-2, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2,-2, 0, 0,-2,-2, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2, 2,-2, 2,-2, 2,-2, 2/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1,-1,-1,-1, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1,-1, 0, 0,-1,-1, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1, 1,-1, 1,-1, 1,-1, 1/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)
            vixyz(7)=ior(int3x(fix+vix(7)),ior(int3y(fiy+viy(3)),int3z(fiz+viz(3))))
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(4)),int3z(fiz+viz(4))))
          case(2)
            vixyz(7)=ior(fint3x(fix+vix(7)),ior(fint3y(fiy+viy(3)),fint3z(fiz+viz(3))))
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(4)),fint3z(fiz+viz(4))))
          end select
          vcboxaddr( 20) = vixyz(7) + 32
          vcboxaddr( 21) = vixyz(7) + 40
          vcboxaddr( 22) = vixyz(7) + 48
          vcboxaddr( 23) = vixyz(7) + 56
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + 8
          vcboxaddr( 26) = vixyz(8) +16
          vcboxaddr( 27) = vixyz(8) +24
          dist =368+child
          fwdstart = 21
          case(2)
          select case(int3exist)
          case(1)
            vixyz(7)=ior(int3x(fix+vix5(7)),ior(int3y(fiy+viy5(3)),int3z(fiz+viz5(3))))
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(4)),int3z(fiz+viz5(4))))
          case(2)
            vixyz(7)=ior(fint3x(fix+vix5(7)),ior(fint3y(fiy+viy5(3)),fint3z(fiz+viz5(3))))
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(4)),fint3z(fiz+viz5(4))))	  
          end select
	  vcboxaddr( 20) = vixyz(7) + 4
          vcboxaddr( 21) = vixyz(7) + 5
          vcboxaddr( 22) = vixyz(7) + 6
          vcboxaddr( 23) = vixyz(7) + 7
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + 1
          vcboxaddr( 26) = vixyz(8) + 2
          vcboxaddr( 27) = vixyz(8) + 3
          
          fwdstart = 21
	  end select
	  case(2)
	  
	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(4)),int3z(fiz+pviz(4))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(4)),int3z(fiz+pviz(4))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(4)),fint3z(fiz+pviz(4))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(4)),fint3z(fiz+pviz(4))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(1) + poffset(4)
          vcboxaddr(  3) = vixyz(2)
          vcboxaddr(  4) = vixyz(3)
          vcboxaddr(  5) = vixyz(3) + poffset(2)
          vcboxaddr(  6) = vixyz(3) + poffset(4)
          vcboxaddr(  7) = vixyz(3) + poffset(6)
          vcboxaddr(  8) = vixyz(4)
          vcboxaddr(  9) = vixyz(4) + poffset(2)
          vcboxaddr( 10) = vixyz(5)
          vcboxaddr( 11) = vixyz(5) + poffset(1)
          vcboxaddr( 12) = vixyz(5) + poffset(4)
          vcboxaddr( 13) = vixyz(5) + poffset(5)
          vcboxaddr( 14) = vixyz(6)
          vcboxaddr( 15) = vixyz(6) + poffset(1)
          vcboxaddr( 16) = vixyz(7)
          vcboxaddr( 17) = vixyz(7) + poffset(1)
          vcboxaddr( 18) = vixyz(7) + poffset(2)
          vcboxaddr( 19) = vixyz(7) + poffset(3)
          vcboxaddr( 20) = vixyz(7) + poffset(4)
          vcboxaddr( 21) = vixyz(7) + poffset(5)
          vcboxaddr( 22) = vixyz(7) + poffset(6)
          vcboxaddr( 23) = vixyz(7) + poffset(7)
          vcboxaddr( 24) = vixyz(8)
          vcboxaddr( 25) = vixyz(8) + poffset(1)
          vcboxaddr( 26) = vixyz(8) + poffset(2)
          vcboxaddr( 27) = vixyz(8) + poffset(3)
	  
          dist =368+child
          fwdstart = 21
      end select
      end subroutine caladdr47
      subroutine caladdr48(int3x,int3y,int3z,boxaddr,px,py,pz,fix,fiy,fiz,child,fwdstart,vcboxaddr,dist,fwdonly,fmmpass,int3exist)
      ! chksm=       16171
      !fwdonly  1:calc higher indexes+own motherbox only 2:calculate all indexes
      use getneighbors_vars
      implicit none
      integer(kind=fmm_integer) int3x(0:*),int3y(0:*),int3z(0:*)
      integer(kind=fmm_integer) :: child,boxaddr,fix,fiy,fiz,px,py,pz
      integer(kind=fmm_integer) :: fwdstart,fwdonly,fmmpass,int3exist
      integer(kind=fmm_integer),dimension(27) :: vcboxaddr
      integer(kind=fmm_integer) :: dist
      integer(kind=fmm_integer),dimension(8), target :: vix = (/-2,-2,-2,-2, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy = (/-2,-2, 0, 0,-2,-2, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz = (/-2, 0,-2, 0,-2, 0,-2, 0/)
      integer(kind=fmm_integer),dimension(8), target :: vix5 = (/-1,-1,-1,-1, 0, 0, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viy5 = (/-1,-1, 0, 0,-1,-1, 0, 0/)
      integer(kind=fmm_integer),dimension(8), target :: viz5 = (/-1, 0,-1, 0,-1, 0,-1, 0/)
      integer(kind=fmm_integer),dimension(7), target :: offset = (/8,16,24,32,40,48,56/)
      integer(kind=fmm_integer),dimension(7), target :: offset5 = (/1,2,3,4,5,6,7/)
      integer(kind=fmm_integer),dimension(:) , pointer :: pvix,pviy,pviz,poffset
      integer(kind=fmm_integer),dimension(8) :: vixyz  

      select case(fwdonly)
        case(1)
        select case(fmmpass)
          case(1)
          select case(int3exist)
          case(1)
            vixyz(8)=ior(int3x(fix+vix(8)),ior(int3y(fiy+viy(4)),int3z(fiz+viz(4))))
          case(2)
            vixyz(8)=ior(fint3x(fix+vix(8)),ior(fint3y(fiy+viy(4)),fint3z(fiz+viz(4))))
          end select
          vcboxaddr( 20) = vixyz(8)
          vcboxaddr( 21) = vixyz(8) + 8
          vcboxaddr( 22) = vixyz(8) + 16
          vcboxaddr( 23) = vixyz(8) + 24
          vcboxaddr( 24) = vixyz(8) + 32
          vcboxaddr( 25) = vixyz(8) + 40
          vcboxaddr( 26) = vixyz(8) + 48
          vcboxaddr( 27) = vixyz(8) + 56
          dist =376+child
          fwdstart = 21
          case(2)
          select case(int3exist)
          case(1)
            vixyz(8)=ior(int3x(fix+vix5(8)),ior(int3y(fiy+viy5(4)),int3z(fiz+viz5(4))))
          case(2)
            vixyz(8)=ior(fint3x(fix+vix5(8)),ior(fint3y(fiy+viy5(4)),fint3z(fiz+viz5(4))))
          end select
          vcboxaddr( 20) = vixyz(8)
          vcboxaddr( 21) = vixyz(8) + 1
          vcboxaddr( 22) = vixyz(8) + 2
          vcboxaddr( 23) = vixyz(8) + 3
          vcboxaddr( 24) = vixyz(8) + 4
          vcboxaddr( 25) = vixyz(8) + 5
          vcboxaddr( 26) = vixyz(8) + 6
          vcboxaddr( 27) = vixyz(8) + 7
          
          fwdstart = 21
	 end select 
	case(2)

	  select case(fmmpass)
	    case(1)
	    pvix => vix
    	    pviy => viy
	    pviz => viz
	    poffset => offset	    	    
	    case(2)
	    pvix => vix5
    	    pviy => viy5
	    pviz => viz5
	    poffset => offset5    	    
	  end select
	  	
          select case(int3exist)
          case(1)
            vixyz(1)=ior(int3x(fix+pvix(1)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(2)=ior(int3x(fix+pvix(2)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(3)=ior(int3x(fix+pvix(3)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(4)=ior(int3x(fix+pvix(4)),ior(int3y(fiy+pviy(4)),int3z(fiz+pviz(4))))
            vixyz(5)=ior(int3x(fix+pvix(5)),ior(int3y(fiy+pviy(1)),int3z(fiz+pviz(1))))
            vixyz(6)=ior(int3x(fix+pvix(6)),ior(int3y(fiy+pviy(2)),int3z(fiz+pviz(2))))
            vixyz(7)=ior(int3x(fix+pvix(7)),ior(int3y(fiy+pviy(3)),int3z(fiz+pviz(3))))
            vixyz(8)=ior(int3x(fix+pvix(8)),ior(int3y(fiy+pviy(4)),int3z(fiz+pviz(4))))
          case(2)
            vixyz(1)=ior(fint3x(fix+pvix(1)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(2)=ior(fint3x(fix+pvix(2)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(3)=ior(fint3x(fix+pvix(3)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(4)=ior(fint3x(fix+pvix(4)),ior(fint3y(fiy+pviy(4)),fint3z(fiz+pviz(4))))
            vixyz(5)=ior(fint3x(fix+pvix(5)),ior(fint3y(fiy+pviy(1)),fint3z(fiz+pviz(1))))
            vixyz(6)=ior(fint3x(fix+pvix(6)),ior(fint3y(fiy+pviy(2)),fint3z(fiz+pviz(2))))
            vixyz(7)=ior(fint3x(fix+pvix(7)),ior(fint3y(fiy+pviy(3)),fint3z(fiz+pviz(3))))
            vixyz(8)=ior(fint3x(fix+pvix(8)),ior(fint3y(fiy+pviy(4)),fint3z(fiz+pviz(4))))
          end select
          vcboxaddr(  1) = vixyz(1)
          vcboxaddr(  2) = vixyz(2)
          vcboxaddr(  3) = vixyz(2) + poffset(4)
          vcboxaddr(  4) = vixyz(3)
          vcboxaddr(  5) = vixyz(3) + poffset(2)
          vcboxaddr(  6) = vixyz(4)
          vcboxaddr(  7) = vixyz(4) + poffset(2)
          vcboxaddr(  8) = vixyz(4) + poffset(4)
          vcboxaddr(  9) = vixyz(4) + poffset(6)
          vcboxaddr( 10) = vixyz(5)
          vcboxaddr( 11) = vixyz(5) + poffset(1)
          vcboxaddr( 12) = vixyz(6)
          vcboxaddr( 13) = vixyz(6) + poffset(1)
          vcboxaddr( 14) = vixyz(6) + poffset(4)
          vcboxaddr( 15) = vixyz(6) + poffset(5)
          vcboxaddr( 16) = vixyz(7)
          vcboxaddr( 17) = vixyz(7) + poffset(1)
          vcboxaddr( 18) = vixyz(7) + poffset(2)
          vcboxaddr( 19) = vixyz(7) + poffset(3)
          vcboxaddr( 20) = vixyz(8)
          vcboxaddr( 21) = vixyz(8) + poffset(1)
          vcboxaddr( 22) = vixyz(8) + poffset(2)
          vcboxaddr( 23) = vixyz(8) + poffset(3)
          vcboxaddr( 24) = vixyz(8) + poffset(4)
          vcboxaddr( 25) = vixyz(8) + poffset(5)
          vcboxaddr( 26) = vixyz(8) + poffset(6)
          vcboxaddr( 27) = vixyz(8) + poffset(7)
	  
          dist =376+child
          fwdstart = 21
      end select
      end subroutine caladdr48
