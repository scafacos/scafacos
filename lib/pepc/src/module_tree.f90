! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2014 Juelich Supercomputing Centre, 
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

!>
!> Defines a derived type `t_tree` that represents distributed hashed octrees
!> and associated procedures.
!>
module module_tree
    use module_box, only: t_box
    use module_comm_env, only: t_comm_env
    use module_domains, only: t_decomposition
    use pthreads_stuff, only: t_pthread_with_type
    use module_atomic_ops, only: t_atomic_int
    use module_pepc_types
    use, intrinsic :: iso_c_binding
    implicit none
    private

    integer(kind_key), public, parameter :: TREE_KEY_ROOT    =  1_kind_key
    
    !> data type for communicator request queue
    type, public :: t_request_queue_entry
      type(t_tree_node), pointer :: node
      logical :: eager_request
      type(t_request_eager) :: request
      logical :: entry_valid
    end type    

    integer, public, parameter :: TREE_COMM_ANSWER_BUFF_LENGTH   = 10000 !< amount of possible entries in the BSend buffer for shipping child data
    integer, public, parameter :: TREE_COMM_REQUEST_QUEUE_LENGTH = 400000 !< maximum length of request queue

    integer, public, parameter :: TREE_COMM_THREAD_STATUS_STOPPED  = 1
    integer, public, parameter :: TREE_COMM_THREAD_STATUS_STARTING = 2
    integer, public, parameter :: TREE_COMM_THREAD_STATUS_STARTED  = 3
    integer, public, parameter :: TREE_COMM_THREAD_STATUS_STOPPING = 4
    integer, public, parameter :: TREE_COMM_THREAD_STATUS_WAITING  = 5

    !> data type for tree communicator
    type, public :: t_tree_communicator
      ! request queue
      type(t_request_queue_entry) :: req_queue(TREE_COMM_REQUEST_QUEUE_LENGTH)
      type(t_atomic_int), pointer :: req_queue_bottom !< position of queue bottom in array; pushed away from top by tree users
      type(t_atomic_int), pointer :: req_queue_top !< position of queue top in array; pushed towards bottom by communicator only when sending    
      
      ! counters and timers
      integer*8 :: comm_loop_iterations(3) !< number of comm loop iterations (total, sending, receiving)
      real*8 :: timings_comm(3) !< array for storing internal timing information
      integer :: request_balance !< total (#requests - #answers), should be zero after complete traversal
      integer(kind_node) :: sum_ships   !< total number of node ships
      integer(kind_node) :: sum_fetches !< total number of node fetches

      ! thread data
      type(t_pthread_with_type) :: comm_thread
      type(t_atomic_int), pointer :: thread_status
      integer :: processor_id
    end type t_tree_communicator

    !>
    !> A derived type representing a distributed hashed octree over a collection
    !> of particles.
    !>
    type, public :: t_tree
      integer(kind_particle) :: npart       !< number of particles across all ranks
      integer(kind_particle) :: npart_me    !< number of particles on local rank

      integer(kind_node) :: nleaf       !< number of leaves stored locally
      integer(kind_node) :: ntwig       !< number of twigs stored locally
      integer(kind_node) :: nleaf_me    !< number of leaves that originated on this rank
      integer(kind_node) :: ntwig_me    !< number of twigs that originated on this rank
      
      integer(kind_node)   :: nbranch     !< number of branch nodes in tree
      integer(kind_node)   :: nbranch_me  !< number of branch nodes that originated on this rank
      integer(kind_node)   :: nbranch_max_me !< upper limit estimate for number of local branch nodes

      integer(kind_particle) :: nintmax     !< maximum number of interactions
      
      type(t_tree_node), pointer :: nodes(:) !< array of tree nodes
      integer(kind_node) :: nodes_maxentries !< max number of entries in nodes array
      integer(kind_node) :: nodes_nentries   !< number of entries present in nodes array
      integer(kind_node) :: node_root        !< index of the root node in nodes-array
      
      real*8, allocatable :: boxlength2(:) !< precomputed square of maximum edge length of boxes for different levels - used for MAC evaluation
      
      type(t_box) :: bounding_box               !< bounding box enclosing all particles contained in the tree
      type(t_comm_env) :: comm_env              !< communication environment over which the tree is distributed
      type(t_decomposition) :: decomposition    !< permutation of particles inserted into the tree
      type(t_tree_communicator) :: communicator !< associated communicator structure
    end type t_tree

    public tree_create
    public tree_allocated
    public tree_provision_node
    public tree_insert_node
    public tree_insert_node_at_index
    public tree_traverse_to_key
    public tree_node_connect_children
    public tree_check
    public tree_dump
    public tree_stats
    public tree_destroy

    contains

    !>
    !> Create a tree (allocates memory but does not fill the tree)
    !> 
    !> Uses particle numbers (local and global) to estimate the memory needed
    !> for node storage.
    !> A communication environment over which the tree is distributed can be
    !> supplied as an MPI communicator `comm` or a `t_comm_env` in `comm_env`.
    !> If no environment is supplied, a duplicate of the PEPC global environment
    !> is used.
    !>
    subroutine tree_create(t, nl, n, comm, comm_env)
      use module_tree_node, only: NODE_INVALID
      use treevars, only: interaction_list_length_factor, MPI_COMM_lpepc, np_mult, nlev
      use module_interaction_specific, only: get_number_of_interactions_per_particle
      use module_comm_env, only: comm_env_dup, comm_env_mirror
      use module_timings
      use module_debug
      implicit none

      type(t_tree), intent(inout) :: t !< The tree
      integer(kind_particle), intent(in) :: nl !< Number of local particles to be inserted into the tree
      integer(kind_particle), intent(in) :: n !< Total number of particles across communication ranks
      integer, optional, intent(in) :: comm !< An MPI communicator
      type(t_comm_env), optional, intent(in) :: comm_env !< A communication environment

      integer(kind_node) :: maxaddress
      integer :: i

      call pepc_status('ALLOCATE TREE')
      DEBUG_ASSERT(.not. tree_allocated(t))

      ! initialize tree communication environment
      if (present(comm_env)) then
        call comm_env_mirror(comm_env, t%comm_env)
      else if (present(comm)) then
        call comm_env_mirror(comm, t%comm_env)
      else
        call comm_env_dup(MPI_COMM_lpepc, t%comm_env)
      end if

      t%npart = n
      t%npart_me = nl
      t%nleaf = 0
      t%ntwig = 0
      t%nleaf_me = 0
      t%ntwig_me = 0
      t%nbranch = 0
      t%nbranch_me = 0
      t%nbranch_max_me = 0
      t%nintmax = 0

      call timer_start(t_allocate)

      call get_number_of_interactions_per_particle(t%npart, t%nintmax)
      t%nintmax = interaction_list_length_factor * t%nintmax

      ! Space for hash table
      if (np_mult > 0) then
        maxaddress = max(30_kind_node * t%nintmax + 4_kind_node * t%npart_me, 10000_kind_node)
      else
        maxaddress = int(abs(np_mult) * 10000._8, kind_node)
      end if

      DEBUG_ASSERT(.not. associated(t%nodes))
      t%nodes_maxentries = max(maxaddress, 2_kind_node**15)
      allocate(t%nodes(1:t%nodes_maxentries))
      t%nodes_nentries   = 0_kind_node
      t%node_root        = NODE_INVALID

      if (maxaddress <= t%npart_me ) then
        DEBUG_ERROR('("maxaddress = ", I0, " <= t%npart_me = ", I0, ".", / , "You should increase np_mult.")', maxaddress, t%npart_me)
      end if

      call tree_communicator_create(t%communicator)
      
      ! Preprocessed box sizes for each level
      allocate(t%boxlength2(0:nlev))
      t%boxlength2(0) = maxval(t%bounding_box%boxsize)**2
      do i = 1, nlev
        t%boxlength2(i) =  t%boxlength2(i-1)/4._8
      end do

      call timer_stop(t_allocate)
    end subroutine tree_create


    !>
    !> destroy a tree, freeing all memory used
    !>
    subroutine tree_destroy(t)
      use module_tree_node, only: NODE_INVALID
      use module_domains, only: decomposition_allocated, decomposition_destroy
      use module_comm_env, only: comm_env_destroy
      use module_debug
      implicit none

      type(t_tree), intent(inout) :: t !< the tree

      call pepc_status('DEALLOCATE TREE')
      DEBUG_ASSERT(tree_allocated(t))

      call tree_communicator_destroy(t%communicator)
      call comm_env_destroy(t%comm_env)
      if (decomposition_allocated(t%decomposition)) then
        call decomposition_destroy(t%decomposition)
      end if

      DEBUG_ASSERT(associated(t%nodes))
      t%nodes_nentries   = 0_kind_node
      t%node_root        = NODE_INVALID
      t%nodes_maxentries = 0_kind_node
      deallocate(t%nodes)
      
      DEBUG_ASSERT(allocated(t%boxlength2))
      deallocate(t%boxlength2)
    end subroutine tree_destroy


    !>
    !> returns `.true.` if resources have been allocated for the tree `t`
    !>
    function tree_allocated(t)
      implicit none

      logical :: tree_allocated
      type(t_tree), intent(in) :: t

      tree_allocated = associated(t%nodes) .and. tree_communicator_allocated(t%communicator)
    end function tree_allocated


    !>
    !> allocates resources for the tree communicator `c` and initializes its fields
    !>
    subroutine tree_communicator_create(c)
      use, intrinsic :: iso_c_binding
      use module_atomic_ops, only: atomic_allocate_int, atomic_store_int
      use module_debug
      implicit none

      type(t_tree_communicator), intent(inout) :: c

      DEBUG_ASSERT(.not. tree_communicator_allocated(c))
      c%timings_comm = 0.
      
      call atomic_allocate_int(c%req_queue_bottom)
      call atomic_allocate_int(c%req_queue_top)
      call atomic_allocate_int(c%thread_status)
      if (.not. (associated(c%req_queue_bottom) .and. associated(c%req_queue_top) &
                 .and. associated(c%thread_status))) then
        DEBUG_ERROR(*, "atomic_allocate_int() failed!")
      end if
      call atomic_store_int(c%req_queue_bottom, 0)
      call atomic_store_int(c%req_queue_top, 0)
      call atomic_store_int(c%thread_status, TREE_COMM_THREAD_STATUS_STOPPED)

      c%request_balance =  0
      c%req_queue(:)%entry_valid = .false. ! used in send_requests() to ensure that only completely stored entries are sent form the list
      c%sum_ships = 0
      c%sum_fetches = 0
    end subroutine tree_communicator_create


    !>
    !> deallocates the resources of the tree communicator `c`
    !>
    subroutine tree_communicator_destroy(c)
      use, intrinsic :: iso_c_binding
      use module_atomic_ops, only: atomic_deallocate_int
      use module_atomic_ops, only: atomic_load_int
      use module_debug
      implicit none

      type(t_tree_communicator), intent(inout) :: c

      DEBUG_ASSERT(tree_communicator_allocated(c))
      ! TODO: we could just stop the running thread, if only it was not in another module
      if (atomic_load_int(c%thread_status) == TREE_COMM_THREAD_STATUS_STARTED) then
        DEBUG_ERROR(*, "tree_communicator_destroy() called with comm thread still running!")
      end if

      call atomic_deallocate_int(c%req_queue_bottom)
      call atomic_deallocate_int(c%req_queue_top)
      call atomic_deallocate_int(c%thread_status)
    end subroutine tree_communicator_destroy


    !>
    !> returns `.true.` if resources have been allocated for tree communicator `c`
    !>
    function tree_communicator_allocated(c)
      use, intrinsic :: iso_c_binding
      implicit none

      logical :: tree_communicator_allocated
      type(t_tree_communicator), intent(in) :: c

      tree_communicator_allocated = associated(c%thread_status)
    end function tree_communicator_allocated


    !>
    !> reserves storage for a single tree node.
    !>
    !> the index at which to insert the tree node later on is returned 
    !> in `entry_pointer`.
    !>
    function tree_provision_node(t)
      use module_debug
      implicit none

      integer(kind_node) :: tree_provision_node
      type(t_tree), intent(inout) :: t

      DEBUG_ASSERT(tree_allocated(t))
      if (t%nodes_nentries >= t%nodes_maxentries) then
        DEBUG_ERROR('("Node array full. # Entries: ", I0,"/",I0)', t%nodes_nentries, t%nodes_maxentries)
      end if

      t%nodes_nentries = t%nodes_nentries + 1_kind_node
      tree_provision_node = t%nodes_nentries
    end function tree_provision_node


    !>
    !> inserts the node `n` into the tree `t` at the position `i` that was
    !> previously provisioned using `tree_provision_node`.
    !>
    subroutine tree_insert_node_at_index(t, i, n)
      use module_pepc_types, only: t_tree_node, kind_node
      use module_tree_node, only: tree_node_is_leaf
      use module_debug
      implicit none

      type(t_tree), intent(inout) :: t
      integer(kind_node), intent(in) :: i
      type(t_tree_node), intent(in) :: n

      DEBUG_ASSERT(tree_allocated(t))
      ! keep count of leaves / twigs
      if (tree_node_is_leaf(n)) then
        t%nleaf =  t%nleaf + 1
        if (n%owner == t%comm_env%rank) t%nleaf_me = t%nleaf_me + 1
      else
        t%ntwig =  t%ntwig + 1
        if (n%owner == t%comm_env%rank) t%ntwig_me = t%ntwig_me + 1
      end if

      t%nodes(i) = n
    end subroutine tree_insert_node_at_index


    !>
    !> inserts the tree node `n` into the tree `t`.
    !>
    subroutine tree_insert_node(t, n, entry_pointer)
      use module_pepc_types, only: t_tree_node, kind_node
      use module_debug
      implicit none

      type(t_tree), intent(inout) :: t !< Tree into which to insert the node
      type(t_tree_node), intent(in) :: n !< The tree node to insert
      integer(kind_node), optional, intent(out) :: entry_pointer !< where the node was inserted

      integer(kind_node) :: i

      DEBUG_ASSERT(tree_allocated(t))
      i = tree_provision_node(t)
      call tree_insert_node_at_index(t, i, n)
      if (present(entry_pointer)) entry_pointer = i
    end subroutine tree_insert_node


    !>
    !> checks whether a node of key `k` is contained in tree `t`
    !> and returns the node-index in `n`
    !>
    function tree_traverse_to_key(t, k, n)
      use module_spacefilling, only: level_from_key, is_ancestor_of
      use module_tree_node, only: NODE_INVALID, tree_node_get_first_child, &
        tree_node_get_next_sibling
      use module_debug
      implicit none

      logical :: tree_traverse_to_key
      type(t_tree), intent(in) :: t !< the tree
      integer(kind_key), intent(in) :: k !< key to look up
      integer(kind_node), intent(out) :: n
      
      integer(kind_level) :: l, kl
      
      DEBUG_ASSERT(tree_allocated(t))
      kl = level_from_key(k)
      n = t%node_root
      
      ! for every level `l` from root + 1 to the target level
      do l = level_from_key(TREE_KEY_ROOT) + 1_kind_level, kl
        ! enter the list of nodes on level `l`
        n = tree_node_get_first_child(t%nodes(n))

        ! traverse the list until we reach its end (n == NODE_INVALID) or ...
        do while (n /= NODE_INVALID)
          if (is_ancestor_of(t%nodes(n)%key, l, k, kl)) exit ! ... we find the ancestor at level `l`
          n = tree_node_get_next_sibling(t%nodes(n))
        end do

        if (n == NODE_INVALID) exit ! abort search early if end of list reached
      end do
      
      tree_traverse_to_key = n /= NODE_INVALID
      
      if (tree_traverse_to_key) then
        DEBUG_ASSERT_MSG(k == t%nodes(n)%key, '(" : requested key=",o18, " but found key=", o18)', k, t%nodes(n)%key)
      endif
    end function tree_traverse_to_key


    !>
    !> connects `next_sibling` pointers in list of
    !> given nodes `c` and attaches the first of them to
    !> `n` via `first_child` within tree `t`
    !>
    subroutine tree_node_connect_children(t, n, c)
      use module_pepc_types, only: kind_node
      use module_tree_node, only: NODE_INVALID
      use module_debug
      implicit none

      type(t_tree), intent(inout) :: t
      integer(kind_node), intent(in) :: n
      integer(kind_node), intent(in) :: c(:)
      
      integer :: ic, nc

      nc = size(c, kind=kind(nc)); DEBUG_ASSERT(nc > 0)
      ! // DEBUG_ASSERT_MSG(all((c(2:nc)%p%key - c(1:nc - 1)%p%key) >= 1), *, "children are not arranged as expected.")
      ! // DEBUG_ASSERT_MSG(2**idim >= c(nc)%p%key - c(1)%p%key, '("= ", I3, ". Children do not all belong to the same parent.")', c(nc)%p%key - c(1)%p%key)

      t%nodes(n)%first_child = c(1)
      
      do ic = 1, nc - 1
        t%nodes(c(ic))%parent       = n
        t%nodes(c(ic))%next_sibling = c(ic + 1)
      end do
      
      t%nodes(c(nc))%parent       = n
      t%nodes(c(nc))%next_sibling = NODE_INVALID
    end subroutine tree_node_connect_children


    !>
    !> Do some quick checks on the tree structure
    !>
    function tree_check(t, callpoint)
      use module_debug
      use module_pepc_types, only: t_tree_node, kind_node
      use module_debug
      implicit none

      logical :: tree_check
      type(t_tree), intent(in) :: t !< the tree
      character(*), intent(in) :: callpoint !< caller

      integer :: nleaf_check, ntwig_check, nleaf_me_check, ntwig_me_check

      call pepc_status('CHECK TREE')

      tree_check = .true.
      nleaf_check = 0
      nleaf_me_check = 0
      ntwig_check = 0
      ntwig_me_check = 0

      DEBUG_ASSERT(tree_allocated(t))
      call tree_check_helper(t%node_root)

      if (t%nleaf /= nleaf_check) then
        DEBUG_WARNING('(3a,i0,/,a,i0,a,i0,a,/,a)', 'Table check called ',callpoint,' by PE',t%comm_env%rank, '# leaves in table = ',nleaf_check,' vs ',t%nleaf,' accumulated', 'Fixing and continuing for now..')
        tree_check = .false.
      end if

      if (t%ntwig /= ntwig_check) then
        DEBUG_WARNING('(3a,i0,/,a,i0,a,i0,a,/,a)', 'Table check called ',callpoint,' by PE',t%comm_env%rank, '# twigs in table = ',ntwig_check,' vs ',t%ntwig,' accumulated', 'Fixing and continuing for now..')
        tree_check = .false.
      end if

      if (t%nleaf_me /= nleaf_me_check) then
        DEBUG_WARNING('(3a,i0,/,a,i0,a,i0,a,/,a)', 'Table check called ',callpoint,' by PE',t%comm_env%rank, '# own leaves in table = ',nleaf_me_check,' vs ',t%nleaf_me,' accumulated', 'Fixing and continuing for now..')
        tree_check = .false.
      end if

      if (t%ntwig_me /= ntwig_me_check) then
        DEBUG_WARNING('(3a,i0,/,a,i0,a,i0,a,/,a)', 'Table check called ',callpoint,' by PE',t%comm_env%rank, '# own twigs in table = ',ntwig_me_check,' vs ',t%ntwig_me,' accumulated', 'Fixing and continuing for now..')
        tree_check = .false.
      end if

      contains

      recursive subroutine tree_check_helper(nid)
        use module_tree_node, only: tree_node_is_leaf, tree_node_get_first_child, &
          tree_node_get_next_sibling, NODE_INVALID
        use module_pepc_types, only: t_tree_node
        implicit none

        integer(kind_node), intent(in) :: nid
        type(t_tree_node), pointer :: n

        integer(kind_node) :: s, ns

        s  = NODE_INVALID
        ns = NODE_INVALID

        n => t%nodes(nid)
        if (tree_node_is_leaf(n)) then
          nleaf_check = nleaf_check + 1
          if (n%owner == t%comm_env%rank) then
            nleaf_me_check = nleaf_me_check + 1
          end if
          return
        else
          ntwig_check = ntwig_check + 1
          if (n%owner == t%comm_env%rank) then
            ntwig_me_check = ntwig_me_check + 1
          end if
        end if

        s = tree_node_get_first_child(n)
        if (s /= NODE_INVALID) then
          do
            call tree_check_helper(s)
            ns = tree_node_get_next_sibling(t%nodes(s))
            if (ns == NODE_INVALID) exit
            s = ns
          end do
        end if
      end subroutine tree_check_helper
    end function tree_check


    !>
    !> Print tree structure from hash table to ipefile
    !>
    subroutine tree_dump(t, particles)
      use treevars
      use module_pepc_types
      use module_spacefilling
      use module_utils
      use module_debug
      use module_tree_node
      implicit none

      type(t_tree), intent(in) :: t
      type(t_particle), optional, intent(in) :: particles(:)

      integer(kind_node) :: i

      call pepc_status('DIAGNOSE')
      DEBUG_ASSERT(associated(t%nodes))
      call debug_ipefile_open()

      ! output node storage
      write(debug_ipefile,'(/a)') 'Node Storage'

      write(debug_ipefile,'(106x,a48)') &
                  '       REQUEST_POSTED                  ', &
                  '       |     HAS_REMOTE_CONTRIBUTIONS  ', &
                  '       |     |HAS_LOCAL_CONTRIBUTIONS  ', &
                  '       |     ||REQUEST_SENT            ', &
                  '       |     |||CHILDREN_AVAILABLE     ', &
                  '       |     ||||       IS_FILL_NODE   ', &
                  '       |     ||||       |IS_BRANCH_NODE', &
                  '       |     ||||       ||             '

      write(debug_ipefile,'(3(x,a10),x,a22,3(x,a10),x,a10,x,a10,4x,3(a8,x),/,115("-"),a26)') &
                   'node', &
                   'owner', &
                   'level', &
                   'key_8', &
                   'parent_node', &
                   'first_child', &
                   'next_sibling', &
                   'leaves', &
                   'dscndnts', &
                   'flags  |', &
                   'l1  ||||', &
                   'g     ||', &
                   '-------V-----VVVV-------VV'

      ! loop over valid enries in node storage
      do i = 1,t%nodes_nentries
        write (debug_ipefile,'(3(x,i10),x,o22,3(x,i10),2(x,i10),4x,l8,2(".",b8.8))') &
                i, &
                t%nodes(i)%owner, &
                level_from_key(t%nodes(i)%key), &
                t%nodes(i)%key, &
                t%nodes(i)%parent, &
                t%nodes(i)%first_child, &
                t%nodes(i)%next_sibling, &
                t%nodes(i)%leaves, &
                t%nodes(i)%descendants, &
                t%nodes(i)%request_posted, &
                t%nodes(i)%flags_local, &
                t%nodes(i)%flags_global
      end do

      write (debug_ipefile,'(///a)') 'Tree structure'

      write(debug_ipefile,'(//a/,x,a,/,179("-"))') 'Twigs from node-list', 'data (see module_interaction_specific::t_tree_node_interaction_data for meaning of the columns)'

      do i = 1,t%nodes_nentries
        if (.not. tree_node_is_leaf(t%nodes(i))) then
          write(debug_ipefile,*) t%nodes(i)%interaction_data
        endif
      end do

      write(debug_ipefile,'(//a/,x,a,/,179("-"))') 'Leaves from node-list', 'data (see module_interaction_specific::t_tree_node_interaction_data for meaning of the columns)'

      do i = 1,t%nodes_nentries
        if (tree_node_is_leaf(t%nodes(i))) then
          write(debug_ipefile,*) t%nodes(i)%interaction_data
        endif
      end do

      if (present(particles)) then
        ! local particles
        write(debug_ipefile,'(//a/,x,a10,x,a,/,189("-"))') 'Local particles', 'index', 'data (see module_module_pepc_types::t_particle for meaning of the columns)'

        do i = lbound(particles, dim=1), ubound(particles, dim=1)
          write(debug_ipefile,'(x,i10,x)',advance='no') i
          write(debug_ipefile,*) particles(i)
        end do
      end if

      call debug_ipefile_close()
    end subroutine

    !>
    !> gather statistics on the tree structure and dump them to a file
    !>
    subroutine tree_stats(t, u)
      use treevars, only: np_mult
      use module_debug
      implicit none
      include 'mpif.h'

      type(t_tree), intent(in) :: t
      integer, intent(in) :: u

      integer :: i, s
      integer(kind_default) :: ierr
      integer(kind_particle), allocatable :: nparticles(:)
      integer(kind_node), allocatable :: fetches(:), ships(:)
      integer(kind_node), allocatable :: total_keys(:), tot_nleaf(:), tot_ntwig(:)
      integer(kind_particle) :: total_part
      integer(kind_node) :: max_nbranch, min_nbranch, nbranch, branch_max_global
      integer(kind_node) :: gmax_keys
      real, save :: part_imbal = 0.
      integer(kind_particle) ::  part_imbal_max, part_imbal_min
      integer(kind_node) :: nkeys_total

      call pepc_status('STATISTICS')
      DEBUG_ASSERT(tree_allocated(t))

      s = t%comm_env%size
      allocate(nparticles(s), fetches(s), ships(s), total_keys(s), tot_nleaf(s), tot_ntwig(s))

      ! particle distrib
      call MPI_GATHER(t%npart_me,    1, MPI_KIND_PARTICLE, nparticles, 1, MPI_KIND_PARTICLE, 0,  t%comm_env%comm, ierr )
      call MPI_GATHER(t%ntwig_me,    1, MPI_KIND_NODE,     tot_ntwig,  1, MPI_KIND_NODE,     0,  t%comm_env%comm, ierr )
      call MPI_GATHER(t%nleaf_me,    1, MPI_KIND_NODE,     tot_nleaf,  1, MPI_KIND_NODE,     0,  t%comm_env%comm, ierr )
      nkeys_total = t%nleaf + t%ntwig
      call MPI_GATHER(nkeys_total,                         1, MPI_KIND_NODE, total_keys, 1, MPI_KIND_NODE,  0,  t%comm_env%comm, ierr )
      call MPI_GATHER(t%communicator%sum_fetches,          1, MPI_KIND_NODE, fetches,    1, MPI_KIND_NODE,  0,  t%comm_env%comm, ierr )
      call MPI_GATHER(t%communicator%sum_ships,            1, MPI_KIND_NODE, ships,      1, MPI_KIND_NODE,  0,  t%comm_env%comm, ierr )
      call MPI_REDUCE(t%nbranch_me, max_nbranch,           1, MPI_KIND_NODE, MPI_MAX,    0, t%comm_env%comm, ierr )
      call MPI_REDUCE(t%nbranch_me, min_nbranch,           1, MPI_KIND_NODE, MPI_MIN,    0, t%comm_env%comm, ierr )
      call MPI_REDUCE(t%nbranch_me, nbranch,               1, MPI_KIND_NODE, MPI_SUM,    0, t%comm_env%comm, ierr)
      call MPI_REDUCE(t%nbranch_max_me, branch_max_global, 1, MPI_KIND_NODE, MPI_MAX,    0, t%comm_env%comm, ierr)
      call MPI_REDUCE(t%nodes_nentries, gmax_keys,             1, MPI_KIND_NODE, MPI_MAX,    0, t%comm_env%comm, ierr )

      part_imbal_max = MAXVAL(nparticles)
      part_imbal_min = MINVAL(nparticles)
      part_imbal = (part_imbal_max - part_imbal_min) / 1.0 / t%npart * s
      total_part = sum(nparticles)

      if (t%comm_env%first) then
        write (u,'(a20,i7,a22)') 'Tree stats for CPU ', t%comm_env%rank, ' and global statistics'
        write (u,*) '######## GENERAL DATA #####################################################################'
        write (u,'(a50,1i12)') '# procs', s
        write (u,'(a50,i12,f12.2,i12)') 'nintmax, np_mult, maxaddress: ',t%nintmax, np_mult, t%nodes_maxentries
        write (u,'(a50,2i12)') 'npp, npart: ', t%npart_me, t%npart
        write (u,'(a50,2i12)') 'total # nparticles, N/P: ', total_part, int(t%npart/s)
        write (u,'(a50,f12.3,2i12)')   'Particle imbalance ave,min,max: ',part_imbal,part_imbal_min,part_imbal_max
        write (u,*) '######## TREE STRUCTURES ##################################################################'
        write (u,'(a50,3i12)') 'local # leaves, twigs, keys: ', t%nleaf_me, t%ntwig_me, t%nleaf_me + t%ntwig_me
        write (u,'(a50,3i12)') 'non-local # leaves, twigs, keys: ',t%nleaf - t%nleaf_me, t%ntwig - t%ntwig_me, t%nleaf + t%ntwig - t%nleaf_me - t%ntwig_me
        write (u,'(a50,3i12,f12.1,a6,i12)') 'final # leaves, twigs, keys, (max): ', t%nleaf, t%ntwig, t%nleaf + t%ntwig, &
                  real((t%nleaf + t%ntwig), kind(0._8)) / (.01 * real(t%nodes_maxentries, kind(0._8))), ' % of ', t%nodes_maxentries
        write (u,'(a50,1i12,1f12.1, a6,1i12)') 'Global max # keys: ',gmax_keys, real(gmax_keys, kind(0._8))/(.01 * real(t%nodes_maxentries, kind(0._8))), ' % of  ', t%nodes_maxentries
        write (u,*) '######## BRANCHES #########################################################################'
        write (u,'(a50,3i12)') '#branches local, max_global, min_global: ', t%nbranch_me, max_nbranch, min_nbranch
        write (u,'(a50,2i12)') '#branches global sum estimated, sum actual: ', branch_max_global, nbranch
        write (u,'(a50,2i12)') 'max res.space for local branches, global br.: ', t%nbranch_max_me, branch_max_global
        write (u,*) '######## WALK-COMMUNICATION ###############################################################'
        write (u,'(a50,2i12)') 'Max # multipole fetches/ships per cpu: ',maxval(fetches), maxval(ships)
        write (u,'(a50,2i12)') 'Min # multipole fetches/ships per cpu: ',minval(fetches), minval(ships)
        write (u,'(a50,2i12)') 'Local #  multipole fetches & ships: ', t%communicator%sum_fetches, t%communicator%sum_ships
        write (u,'(a50,3i12)') '# of comm-loop iterations (tot,send,recv): ', t%communicator%comm_loop_iterations(:)
        write (u,*) '######## DETAILED DATA ####################################################################'
        write (u,'(2a/(4i10,F8.4,4i15))') '        PE     parts     nleaf     ntwig   ratio        nl_keys', &
                  '       tot_keys        fetches          ships', &
                  (i-1,nparticles(i),tot_nleaf(i),tot_ntwig(i),1.0*tot_nleaf(i)/(1.0*tot_ntwig(i)), &
                  total_keys(i)-(tot_nleaf(i)+tot_ntwig(i)),total_keys(i),fetches(i),ships(i),i=1,s)
      end if

      deallocate(nparticles, fetches, ships, total_keys, tot_nleaf, tot_ntwig)
    end subroutine tree_stats
end module module_tree
