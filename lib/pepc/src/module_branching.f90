! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2013 Juelich Supercomputing Centre, 
!                         Forschungszentrum Juelich GmbH,
!                         Germany
! 
! PEPC is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! PEPC is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public License
! along with PEPC.  If not, see <http://www.gnu.org/licenses/>.
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Contains mapper functions for space-filling curves
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_branching
    use treevars
    use module_math_tools
    implicit none
    private

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! physical (real) boundaries of local domain
    integer*8 :: left_limit, left_limit_me, right_limit_me, right_limit

    ! virtual boundaries of local domain
    integer*8 :: left_virt_limit, right_virt_limit

    ! estimation
    integer*8 :: D1, D2                  ! partial domains
    integer*8 :: L                       ! inner limit
    integer*8, public :: branch_max_local        ! estimation for local branches
    integer*8, public :: branch_max_global = -1  ! estimation for global branches
    integer*8, allocatable :: branch_level(:)    ! # branches at every level
    integer*8, allocatable :: branch_level_D1(:) ! partial occupancy of D1
    integer*8, allocatable :: branch_level_D2(:) ! partial occupancy of D2

    integer*8, allocatable :: speedup_potenz(:)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    public find_branches
    public branches_initialize
    public branches_finalize
    public branches_initialize_VLD

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  subroutine-implementation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> initialize Virtual Local Domain data
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine branches_initialize_VLD(particles)
        use module_pepc_types
        implicit none
        type(t_particle), intent(in) :: particles(:)
        ! Initialize VLD-stuff
        call get_virtual_local_domain(particles)
        ! CSBE - Cross Sum Branch Node Estimator
        call get_local_apriori_est()
        call get_global_apriori_est()
    end subroutine


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> initialize module, i.e. allocate data arrays
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine branches_initialize()
        use treevars
        implicit none
        integer :: ilevel

        if (allocated(branch_level))    deallocate(branch_level)
        if (allocated(branch_level_D1)) deallocate(branch_level_D1)
        if (allocated(branch_level_D2)) deallocate(branch_level_D2)
        if (allocated(speedup_potenz))  deallocate(speedup_potenz)

        allocate(branch_level(0:nlev))
        allocate(branch_level_D1(0:nlev))
        allocate(branch_level_D2(0:nlev))
        allocate(speedup_potenz(0:nlev))

        speedup_potenz(:) = (/ (2_8**(idim*ilevel),ilevel=0,nlev) /)

    end subroutine branches_initialize

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> finalize module
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine branches_finalize()
        implicit none

        deallocate(branch_level)
        deallocate(branch_level_D1)
        deallocate(branch_level_D2)
        deallocate(speedup_potenz)

    end subroutine branches_finalize


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> setup virtual local domain (VLD)
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine get_virtual_local_domain(particles)
        use module_pepc_types
        implicit none

        type(t_particle), intent(in) :: particles(1:npp+2) ! TODO: remove this +2, see tree_domains for details
          
        if (num_pe > 1) then
            ! get local key limits
            left_limit_me  = particles(  1)%key
            right_limit_me = particles(npp)%key

            ! get key limits for neighbor tasks
            ! and build virtual limits, so that a minimum set a branch nodes comes arround
            ! boundary tasks can access their boundary space fully only need one virtual limit
            if(me.eq.0)then
                right_limit      = particles(npp+1)%key
                right_virt_limit = bpi(right_limit_me,right_limit)
                left_virt_limit  = 2_8**(nlev * idim)
            else if(me.eq.(num_pe-1))then
                left_limit       = particles(npp+1)%key
                left_virt_limit  = bpi(left_limit,left_limit_me)
                right_virt_limit = 2_8**(idim * nlev+1)-1
            else
                left_limit       = particles(npp+2)%key
                right_limit      = particles(npp+1)%key
                left_virt_limit  = bpi(left_limit,left_limit_me)
                right_virt_limit = bpi(right_limit_me,right_limit)
            end if
        else
            left_virt_limit  = 2_8**(idim * nlev)
            right_virt_limit = 2_8**(idim * nlev+1)-1
        end if

    end subroutine get_virtual_local_domain

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> get apriori estimation of local branch nodes (local CSBE)
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine get_local_apriori_est()
        implicit none
          
        ! First find highest power in the Virtual Domain to ensure a correct branch definition
        L = bpi(left_virt_limit,right_virt_limit)
  
        ! divide in two sub-domains
        ! only the last tasks must get 1 particle more
        ! because it s right limit is possibly not presentable
        D1 = L-left_virt_limit
        if(me.eq.(num_pe-1))then
            D2 = right_virt_limit-L+1
        else
            D2 = right_virt_limit-L
        end if
  
        call get_local_csbe(D1, D2, branch_level, branch_level_D1, branch_level_D2)


        ! estimate local number
        branch_max_local = SUM(branch_level(:))

    end subroutine get_local_apriori_est

        
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> get apriori estimation of global branch nodes (global CSBE)
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine get_global_apriori_est()
        implicit none
           include 'mpif.h'

        integer :: ierr
           
        ! get global estimation
        call MPI_ALLREDUCE(branch_max_local, branch_max_global, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_lpepc, ierr)

    end subroutine get_global_apriori_est

        


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> get keys of all estimated global branches
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine get_global_apriori_branches(resultarray, numbranches)
        implicit none
           include 'mpif.h'

        integer*8 :: allD1(num_pe), allD2(num_pe) , allL(num_pe)
        integer :: ierr
        integer :: i
        integer*8, dimension(1:branch_max_global), intent(out) :: resultarray
        integer, dimension(0:num_pe-1), intent(out) :: numbranches
        integer :: resultidx

        ! get D1 and D2 of all tasks
        call MPI_Allgather( D1, 1, MPI_INTEGER8, allD1, 1, MPI_INTEGER8, MPI_COMM_lpepc, ierr)
        call MPI_Allgather( D2, 1, MPI_INTEGER8, allD2, 1, MPI_INTEGER8, MPI_COMM_lpepc, ierr)
        call MPI_Allgather(  L, 1, MPI_INTEGER8, allL,  1, MPI_INTEGER8, MPI_COMM_lpepc, ierr)

        resultidx = 0
        do i=0,num_pe-1
            numbranches(i) = find_possible_branches(allD1(i), allD2(i), allL(i), resultarray, resultidx)
            resultidx = resultidx + numbranches(i)
        end do

    end subroutine get_global_apriori_branches

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> make cross sum and store
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine get_local_csbe(D1, D2, branch_level, branch_level_D1, branch_level_D2)
        implicit none

        integer*8, intent(in) :: D1, D2
        integer*8 :: ilevel, pos
        integer*8, intent(out) :: branch_level(0:nlev), branch_level_D1(0:nlev), branch_level_D2(0:nlev)
 
        ! get estimation number of branches at all levels
        do ilevel=0,nlev
            pos=idim*(nlev-ilevel)
            branch_level_D1(ilevel)=ibits(D1,pos,idim)
            branch_level_D2(ilevel)=ibits(D2,pos,idim)
            branch_level(ilevel)=branch_level_D1(ilevel) + branch_level_D2(ilevel)
        end do
        
    end subroutine get_local_csbe


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> find branches (based on the CSBE)
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine find_branches()
        use module_htable
        use module_spacefilling
        use module_debug, only : pepc_status
        implicit none

        integer :: ilevel, j
        integer*8 :: possible_branch
        integer*8 :: pos

        call pepc_status('FIND BRANCHES')

        nbranch = 0

        ! add placeholder-bit to inner limit L
        ! first get level of limit L
        ilevel = level_from_key(L)

        ! for D1
        pos=L
        do ilevel=0,nlev
            do j=1,int(branch_level_D1(ilevel))
                pos = pos - speedup_potenz((nlev-ilevel))!2**(idim*(nlev-ilevel))
                possible_branch=shift_key_by_level(pos,-(nlev-ilevel))
             
                ! After local build hashtable should contain branch key
                ! otherwise branch does not exists
                ! if entry exists it is counted as branch
                ! otherwise discarded
                if(testaddr(possible_branch))then ! entry exists
                    nbranch = nbranch + 1
                    pebranch(nbranch) = possible_branch
                end if
            end do
        end do

        ! for D2
        pos=L-1
        do ilevel=0,nlev
            do j=1,int(branch_level_D2(ilevel))
                pos = pos + speedup_potenz((nlev-ilevel))!2**(idim*(nlev-ilevel))
                possible_branch=shift_key_by_level(pos,-(nlev-ilevel))
                
                ! After build hashtable should contain branch key
                ! otherwise branch does not exists
                ! if entry exists it is counted as branch
                ! otherwise discarded
                if(testaddr(possible_branch))then ! entry exists
                    nbranch = nbranch + 1
                    pebranch(nbranch) = possible_branch
                end if
            end do
        end do
       
    end subroutine find_branches


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> calculates all possible branch keys in VLD D1 and D2 with
    !> inner limit L
    !> return keys are put into result array, starting from index resultidx
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function find_possible_branches(D1, D2, L, resultarray, resultidx)
        use module_spacefilling
        implicit none
        integer :: find_possible_branches
        integer*8, intent(in) :: D1, D2, L
        integer*8, intent(inout), dimension(1:branch_max_global) :: resultarray
        integer, intent(in) :: resultidx

        integer :: ilevel, j
        integer*8 :: pbranch_level(0:nlev), pbranch_level_D1(0:nlev), pbranch_level_D2(0:nlev)
        integer*8 :: pos

        find_possible_branches = 0

        ! add placeholder-bit to inner limit L
        ! first get level of limit L
        ilevel = level_from_key(L)

        call get_local_csbe(D1, D2, pbranch_level, pbranch_level_D1, pbranch_level_D2)

        ! for D1
        pos=L
        do ilevel=0,nlev
            do j=1,int(pbranch_level_D1(ilevel))
                pos = pos - speedup_potenz((nlev-ilevel))!2**(idim*(nlev-ilevel))
                resultarray(resultidx+find_possible_branches) = shift_key_by_level(pos,-(nlev-ilevel))
                find_possible_branches = find_possible_branches + 1
            end do
        end do

        ! for D2
        pos=L-1
        do ilevel=0,nlev
            do j=1,int(pbranch_level_D2(ilevel))
                pos = pos + speedup_potenz((nlev-ilevel))!2**(idim*(nlev-ilevel))
                resultarray(resultidx+find_possible_branches) = shift_key_by_level(pos,-(nlev-ilevel))
                find_possible_branches = find_possible_branches + 1
            end do
        end do

    end function find_possible_branches


end module module_branching
