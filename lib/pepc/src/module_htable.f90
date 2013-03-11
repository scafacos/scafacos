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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates functions for accessing, manipulating, and verifying hash table data
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module module_htable
    implicit none

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  type declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Hash table datatype - 36 bytes per entry
    type :: t_hash
        integer*8 :: key           !< Key
        integer   :: node          !< Address of particle/pseudoparticle data
        integer   :: link          !< Pointer to next empty address in table in case of collision
        integer   :: leaves        !< # leaves contained within twig (=1 for leaf, npart for root)
        integer   :: childcode     !< Byte code indicating position of children (twig node); particle label (leaf node)
        integer   :: owner         !< Node owner (for branches)
	integer   :: level         !< level_from_key(key)
    end type t_hash

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    type (t_hash), public, target, allocatable :: htable(:) !< hash table
    integer*8,     public ::  hashconst  !< hashing constants

    integer*8, public, parameter :: KEY_INVALID = -1
    integer*8, public, parameter :: KEY_EMPTY   =  0

    ! bits in childcode to be set when children are requested, the request has been sent, and they have arrived
    integer, public, parameter :: CHILDCODE_BIT_REQUEST_POSTED           =  8 !< this bit is used inside the childcode to denote that a request for children information is already in the request queue
    integer, public, parameter :: CHILDCODE_BIT_CHILDREN_AVAILABLE       =  9 !< this bit is used inside the childcode to denote that children information for the node is available in the local hashtable
    integer, public, parameter :: CHILDCODE_BIT_REQUEST_SENT             = 10 !< this bit is used inside the childcode to denote that children information has already been requested from the owner
    integer, public, parameter :: CHILDCODE_BIT_HAS_LOCAL_CONTRIBUTIONS  = 11 !< this bit is set for all nodes that contain some local nodes beneath them
    integer, public, parameter :: CHILDCODE_BIT_HAS_REMOTE_CONTRIBUTIONS = 12 !< this bit is set for all nodes that contain some remote nodes beneath them
    integer, public, parameter :: CHILDCODE_BIT_IS_BRANCH_NODE           = 13 !< this bit is set for all branch nodes (set in tree_exchange)
    integer, public, parameter :: CHILDCODE_BIT_IS_FILL_NODE             = 14 !< this bit is set for all nodes that are above (towards root) branch nodes
    integer, public, parameter :: CHILDCODE_CHILDBYTE                    = b'11111111' !< bits that contain the children information for this node

    integer, public :: maxaddress                    !< max address allowed in #table

    integer, public, allocatable :: free_addr(:)    !< List of free #table addresses (for HASHENTRY routine)
    integer, public, allocatable :: point_free(:)   !< Pointer to free address index
    integer, parameter :: free_lo = 1024             !< min address allowed for resolving collisions (from 4th level up)
    integer :: iused                                  !< counter for collision resolution array free_addr()
    integer :: sum_unused                             !< # free addresses


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    public children_available
    public get_next_node_key
    public get_childkeys
    public make_hashentry
    public key2addr
    public testaddr
    public htable_clear
    public htable_clear_and_insert_root
    public htable_prepare_address_list
    public check_table
    public diagnose_tree
    public htable_remove_keys
    public htable_entry_is_valid
    public htable_entry_is_leaf

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    integer, private, parameter :: start_child_idx = 0 !< index of first child to be used in traversal - do not change, currently not completely implemented
    type (t_hash), private, parameter :: HASHENTRY_EMPTY = t_hash(KEY_EMPTY,0,-1,0,0,0, 0) !< constant for empty hashentry

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  subroutine-implementation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    contains


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> checks whether the given htable-entry is a leaf or a twig
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! TODO: use this function everywhere consistently
    function htable_entry_is_leaf(hashaddr)
      implicit none
      integer, intent(in) :: hashaddr
      logical:: htable_entry_is_leaf
      ! TODO: this way of identifying leaves/twigs is not good -- use number of leaves or (even better) a flag in the childcode, instead
      htable_entry_is_leaf = (htable(hashaddr)%node > 0)
    end function


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> checks whether the given htable-entry is valid and not empty
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! TODO: use this function everywhere consistently
    function htable_entry_is_valid(hashaddr)
      implicit none
      integer, intent(in) :: hashaddr
      logical:: htable_entry_is_valid
      integer*8 :: key

      key = htable(hashaddr)%key

      htable_entry_is_valid = (key .ne. KEY_EMPTY) .and. (key .ne. KEY_INVALID)
    end function


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> empties the htable
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine htable_clear()
        implicit none

        htable = HASHENTRY_EMPTY ! TODO: need list of 'live' adresses to speed this up
                                 ! possible solution: use a "bitmap" of occupied addresses in htable
                                 ! then, only this bitmap has to be cleared upon startup
                                 ! and every test for occupancy of a htable entry is done in this bitmap

    end subroutine


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> empties the htable, inserts the root node and prepares the
    !> collision resolution list
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine htable_clear_and_insert_root()
        use treevars, only : ntwig, me
	use module_spacefilling, only : level_from_key
        implicit none

        call htable_clear()

        ntwig = 1

        htable(1) = t_hash(1_8, -1, -1, 0, IBSET(0, CHILDCODE_BIT_CHILDREN_AVAILABLE), me, level_from_key(1_8))

        call htable_prepare_address_list()

    end subroutine


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> prepares the collision resolution list by traversing the full htable
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine htable_prepare_address_list()
        implicit none
        integer :: i

        ! build list of free addresses for faster collision resolution on insertion into htable
        sum_unused = 0
        iused      = 1   ! reset used-address counter
        do i=0, maxaddress
           if (htable(i)%node == 0 .and. htable(i)%key /=-1 .and. i > free_lo) then
              sum_unused            = sum_unused + 1
              free_addr(sum_unused) = i            ! Free address list for resolving collisions
              point_free(i)         = sum_unused   ! Index
           else
              point_free(i)         = 0
           endif
        enddo

    end subroutine


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> checks whether children for certain htable-address are locally available
    !> or have to be requested
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    pure function children_available(addr)
      implicit none
      logical :: children_available
      integer, intent(in) :: addr

      children_available = btest(htable( addr )%childcode, CHILDCODE_BIT_CHILDREN_AVAILABLE)

    end function


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Function to return key of next node in local tree-walk,
    !> i.e. search for next sibling, uncle, great-uncle etc
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function get_next_node_key(keyin)

        use treevars
        use module_spacefilling

        implicit none
        integer*8 :: get_next_node_key
        integer*8, intent(in) :: keyin

        integer*8 :: search_key, parent_key
        integer   :: parent_addr
        integer   :: parent_child_byte 
        integer   :: search_child_idx

        search_key = keyin

        ! search for next sibling, uncle, great-uncle etc
        do while (search_key > 1) ! loop through parent nodes up to root
            parent_key        = parent_key_from_key(search_key)
            parent_addr       = key2addr( parent_key ,"next_node(), get parent_addr" )
            parent_child_byte = ibits( htable( parent_addr ) % childcode, 0, 2**idim)

            search_child_idx  = child_number_from_key(search_key) ! lower three bits of key

            do ! loop over all siblings
                search_child_idx   = modulo(search_child_idx + 1, 2**idim) ! get next sibling, possibly starting again from first one

                ! if sibling-loop wrapped and reached starting point again --> go up one level
                if ( search_child_idx == start_child_idx ) then
                    search_key = parent_key      ! go up one level
                    exit
                endif

                ! if sibling exists: next_node has been found
                if ( btest(parent_child_byte, search_child_idx) ) then
                    get_next_node_key = child_key_from_parent_key(parent_key, search_child_idx) ! assemble next_node out of parent-key and new sibling-index
                    return
                endif
            end do
        end do

        get_next_node_key  = 1 ! nothing has been found, i.e. top-right corner reached: set pointer=root

    end function get_next_node_key

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> returns the keys of all children, that are attached to the
    !> node at a particular htable address
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine get_childkeys(addr, childnum, childkeys)
        use treevars, only: idim
        use module_spacefilling
        implicit none
        integer, intent(in) :: addr
        integer, intent(out) :: childnum
        integer*8, dimension(:), intent(out) :: childkeys
        integer :: i

        integer*8 :: keyhead
        integer :: childcode

        keyhead   = shift_key_by_level(htable(addr)%key, 1)
        childcode = htable(addr)%childcode
        childnum = 0

        do i=0,2**idim - 1
          if (btest(childcode, i)) then
            childnum            = childnum + 1
            childkeys(childnum) = ior(keyhead, 1_8*i)
          end if
        end do

    end subroutine


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Make entry in hash-table - returns address 'hashaddr'
    !> Resolve collision if necessary
    !> make_hashentry == .true. if anything went fine
    !> make_hashentry == .false. if key already exists in htable, hashaddr is set to
    !> return current address, but the htable-entry itself is *not* modified
    !>
    !> value of hashentry%link is ignored
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function make_hashentry(hashentry, hashaddr)

        use treevars
        use module_debug
        implicit none

        logical :: make_hashentry
        type(t_hash), intent(in) :: hashentry
        integer, intent(out) :: hashaddr ! address in # table returned to calling routine

        integer :: link

        if (.not. testaddr(hashentry%key, hashaddr)) then
          ! this key does not exist in the htable 
          make_hashentry = .true.

          if (hashaddr .eq. -1) then
            ! the first entry is already empty
            hashaddr = int(IAND( hashentry%key, hashconst))

            if (point_free(hashaddr) /= 0) then     ! Check if new address in collision res. list
                free_addr( point_free(hashaddr) ) = free_addr(sum_unused)  ! Replace free address with last on list
                point_free(free_addr(sum_unused)) = point_free(hashaddr)   ! Reset pointer
                point_free(hashaddr)              = 0
                sum_unused = sum_unused - 1
            endif
          else
            ! we are at the end of a linked list --> create new entry
            htable( hashaddr )%link = free_addr(iused)
            hashaddr                = htable( hashaddr )%link
            iused                   = iused + 1
          endif

          ! check if new entry is really empty
          if ((htable(hashaddr)%node /= 0 ) .or. (htable( hashaddr )%key/=0)) then
            write (*,*) 'Something wrong with address list for collision resolution (free_addr in treebuild)'
            write (*,*) 'PE ',me,' key ',hashentry%key,' entry',hashaddr,' used ',iused,'/',sum_unused
            write (*,*) "htable(hashaddr):  ", htable(hashaddr)
            write (*,*) "desired entry:     ", hashentry
            call debug_mpi_abort()
          endif

          link = htable( hashaddr )%link ! would be overwritten by the next code line
          htable( hashaddr ) = hashentry
          htable( hashaddr )%link = link
        else
          ! this key does already exists in the htable - as 'hashaddr' we return its current address
          make_hashentry = .false.
        endif

    end function

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Invalidates entries in the htable, i.e. sets their key to KEY_INVALID
    !> TODO: this function does not free the nodelist-storage and does
    !> not fix any connections inside the tree, esp. concerning the
    !> children-available-flag (it even does not care for them)
    !> additionally, it does not modify nleaf or ntwig which would be
    !> necessary to survive check_table() if the entry really was removed
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine htable_remove_keys(keys, num_keys)
      implicit none
      integer*8, intent(in) :: keys(num_keys)
      integer, intent(in) :: num_keys

      integer :: i

      do i=1,num_keys
        htable(  key2addr(keys(i), 'htable_remove_keys')  )%key = KEY_INVALID
      end do

    end subroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Calculates inverse mapping from the hash-key to the hash-address.
    !>
    !> @param[in] keyin inverse mapping candidate.
    !> @param[in] cmark a description.
    !> @param[out] key2addr the adress if the key exists
    !> @exception if key does not exist, the whole program is aborted
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function key2addr(keyin,cmark)

        use treevars
        use module_debug
        implicit none

        integer*8, intent(in)  :: keyin
        character(LEN=*) :: cmark
        integer :: key2addr

        if (.not. testaddr(keyin, key2addr)) then
          ! could not find key
          DEBUG_WARNING_ALL('("Key not resolved in KEY2ADDR at ",a)', cmark)
          DEBUG_WARNING_ALL('("Bad address, check #-table and key list for PE", I7)', me)
          DEBUG_WARNING_ALL('("key                  (oct) = ", o22)', keyin)
          DEBUG_WARNING_ALL('("initial address      (dez) = ", i22)', int(IAND( keyin, hashconst)))
          DEBUG_WARNING_ALL('("   last address      (dez) = ", i22)', key2addr)
          if (.not. (key2addr == -1)) then
            DEBUG_WARNING_ALL('("htable(lastaddr)%key (oct) = ", o22)', htable(key2addr)%key)
          end if
          DEBUG_WARNING_ALL('("# const              (dez) = ", i22)', hashconst)
          DEBUG_WARNING_ALL('("     maxaddress      (dez) = ", i22)', maxaddress)
          call diagnose_tree()
          call debug_mpi_abort()
        endif

    end function key2addr


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Calculates inverse mapping from the hash-key to the hash-address.
    !>
    !> @param[in] keyin inverse mapping candidate.
    !> @return .true. if the key has been found, .false. otherwise
    !> @param[out] addr address if candidate exists, address of last entry in linked list otherwise, -1 if already the first lookup failed
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function testaddr(keyin, addr)

        use treevars
        implicit none
        include 'mpif.h'

        integer*8, intent(in)  :: keyin
        integer :: ires
        logical :: testaddr
        integer, intent(out), optional :: addr
        integer :: nextaddr, lastaddr

        nextaddr = int(IAND( keyin, hashconst))     ! cell address hash function

        if (htable( nextaddr )%key == KEY_EMPTY) then ! already the first entry is empty
           testaddr                = .false.  ! key does not exist in htable
           if (present(addr)) addr = -1       ! we return -1
           return
        endif

        ires     =  1 ! counter for number of htable lookups

        do while ( htable( nextaddr )%key .ne. keyin )
          lastaddr = nextaddr
          nextaddr = htable( nextaddr )%link    ! look at next linked entry
          ires     = ires + 1

          if (   (nextaddr == -1) & ! reached end of linked list without finding the key --> node is not in htable or htable is corrupt
            .or. (ires >= maxaddress) ) & ! we probed at as many positions as the htable has entries --> circular linked list or htable corrupt
            then
              testaddr                = .false.  ! key does not exist in htable
              if (present(addr)) addr = lastaddr ! we return last entry in the linked list
              return
          endif

        end do

        testaddr                = .true.   ! key has been found in htable
        if (present(addr)) addr = nextaddr ! we return its address

    end function testaddr


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Do some quick checks on the tree structure
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine check_table(callpoint)
        use treevars
        use module_debug
        implicit none
        character(*) :: callpoint
        integer :: nleaf_check, ntwig_check, nleaf_me_check, ntwig_me_check
        logical :: error

        call pepc_status('CHECK TABLE')

        error = .false.

        ntwig_check = count(mask =  htable%node <0 )
        nleaf_check = count(mask =  htable%node >0 )
        nleaf_me_check = count(mask = htable%owner==me .and. htable%node >0 )
        ntwig_me_check = count(mask = htable%owner==me .and. htable%node <0 )

        if (nleaf /= nleaf_check) then
            DEBUG_WARNING('(3a,i0,/,a,i0,a,i0,a,/,a)', 'Table check called ',callpoint,' by PE',me,
                                                       '# leaves in table = ',nleaf_check,'vs ',nleaf,'accumulated',
                                                       'Fixing and continuing for now..')
        !     nleaf = nleaf_check
            error = .true.
        endif

        if (ntwig /= ntwig_check) then
            DEBUG_WARNING('(3a,i0,/,a,i0,a,i0,a,/,a)', 'Table check called ',callpoint,' by PE',me,
                                                       '# twigs in table = ',ntwig_check,'vs ',ntwig,'accumulated',
                                                       'Fixing and continuing for now..')
        !     ntwig = ntwig_check
            error = .true.
        endif

        if (nleaf_me /= nleaf_me_check) then
            DEBUG_WARNING('(3a,i0,/,a,i0,a,i0,a,/,a)', 'Table check called ',callpoint,' by PE',me,
                                                       '# own leaves in table = ',nleaf_me_check,'vs ',nleaf_me,'accumulated',
                                                       'Fixing and continuing for now..')
            nleaf_me = nleaf_me_check
            error = .true.
        endif
        if (ntwig_me /= ntwig_me_check) then
            DEBUG_WARNING('(3a,i0,/,a,i0,a,i0,a,/,a)', 'Table check called ',callpoint,' by PE',me,
                                                       '# own twigs in table = ',ntwig_me_check,'vs ',ntwig_me,'accumulated',
                                                       'Fixing and continuing for now..')

            ntwig_me = ntwig_me_check
            error = .true.
        endif

        if (error) then
          call diagnose_tree()
        endif

    end subroutine check_table



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Print tree structure from hash table to ipefile
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine diagnose_tree(particles)
        use treevars
        use module_pepc_types
        use module_spacefilling
        use module_utils
        use module_debug, only : debug_ipefile_open, debug_ipefile_close, debug_ipefile, pepc_status
        implicit none

        type(t_particle), optional, intent(in) :: particles(1:npp)
        integer*8, dimension(ntwig) :: node_twig      ! twig-nodes
        integer*8, dimension(nleaf) :: node_leaf      ! leaf-nodes

        character(1) :: collision
        integer :: i

        call pepc_status('DIAGNOSE')

        call debug_ipefile_open()

        ! output hash table

        write(debug_ipefile,'(/a)') 'Hash table'

        write(debug_ipefile,'(164x,a35)') &
                    "IS_FILL_NODE              ", &
                    "|IS_BRANCH_NODE           ", &
                    "||HAS_REMOTE_CONTRIBUTIONS", &
                    "|||HAS_LOCAL_CONTRIBUTIONS", &
                    "||||REQUEST_SENT          ", &
                    "|||||CHILDREN_AVAILABLE   ", &
                    "||||||REQUEST_POSTED      "

        write(debug_ipefile,'(5(x,a10),3(x,a22),x,a14,x,a10,4x,a5,a30,/,173("-"),7("V")," 76543210")') &
                     'entry_10', &
                     'entry_8', &
                     'owner', &
                     'level', &
                     'node',  &
                     'key_8', &
                     'key_10', &
                     'parent_8', &
                     'collision link', &
                     'leaves', &
                     'flags', &
                     '||||||| childcod'

        ! write(debug_ipefile,'(154x,a)') " 3      .   2    .     1  .        "
        ! write(debug_ipefile,'(154x,a)') "10987654.32109876.54321098.76543210"

        do i=0,maxaddress
            if (htable_entry_is_valid(i)) then
              ! flag  collisions
              if (htable(i)%link/= -1 ) then
                collision="C"
              else
                collision=" "
              endif

              write (debug_ipefile,'(x,i10,x,o10,3(x,i10),x,o22,x,i22,x,o22,x,a1,x,i12,x,i10,4x,3(b8.8,"."),b8.8)') &
                      i,                   &
                      i,                   &
                      htable(i)%owner,     &
                      level_from_key(htable(i)%key), &
                      htable(i)%node,      &
                      htable(i)%key,       &
                      htable(i)%key,       &
                      parent_key_from_key(htable(i)%key), &
                      collision,           &
                      htable(i)%link,      &
                      htable(i)%leaves,    &
                      ishft(iand(htable(i)%childcode, Z'FF000000'), -24), &
                      ishft(iand(htable(i)%childcode, Z'00FF0000'), -16), &
                      ishft(iand(htable(i)%childcode, Z'0000FF00'), -08), &
                      ishft(iand(htable(i)%childcode, Z'000000FF'), -00)
           end if
        end do

        write (debug_ipefile,'(///a)') 'Tree structure'

        ! get node indices of twig nodes from hash table
        node_twig(  1:ntwig)   = pack(htable(0:maxaddress)%node,mask=htable(0:maxaddress)%node<0)
        call sort(node_twig(:))

        write(debug_ipefile,'(//a/,x,a10,x,a,/,189("-"))') 'Twigs from hash-table', 'node', 'data (see module_interaction_specific::t_tree_node_interaction_data for meaning of the columns)'

        do i=ntwig,1,-1
          write(debug_ipefile,'(x,i10,x)',advance='no') node_twig(i)
          write(debug_ipefile,*) tree_nodes(node_twig(i))
        end do

        ! get node indices of leaf nodes from hash table
        node_leaf(  1:nleaf)   = pack(htable(0:maxaddress)%node,mask=htable(0:maxaddress)%node>0)
        call sort(node_leaf(:))

        write(debug_ipefile,'(//a/,x,a10,x,a,/,189("-"))') 'Leaves from hash-table', 'node', 'data (see module_interaction_specific::t_tree_node_interaction_data for meaning of the columns)'

        do i=1,nleaf
          write(debug_ipefile,'(x,i10,x)',advance='no') node_leaf(i)
          write(debug_ipefile,*) tree_nodes(node_leaf(i))
        end do

        if (present(particles)) then
          ! local particles
          write(debug_ipefile,'(//a/,x,a10,x,a,/,189("-"))') 'Local particles', 'index', 'data (see module_module_pepc_types::t_particle for meaning of the columns)'

          do i=1,npp
            write(debug_ipefile,'(x,i10,x)',advance='no') i
            write(debug_ipefile,*) particles(i)
          end do
        endif

        call debug_ipefile_close()

    end subroutine diagnose_tree


end module module_htable
