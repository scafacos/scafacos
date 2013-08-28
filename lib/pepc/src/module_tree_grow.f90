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


!>
!> Contains procedures to construct a tree from a collection of particles.
!>
module module_tree_grow
  use module_pepc_types
  implicit none
  private

  public tree_grow

  contains

  !>
  !> Builds the tree from the given particles, redistributes particles
  !> to other MPI ranks if necessary (i.e. reallocates particles and changes np_local)
  !>
  subroutine tree_grow(t, p)
    use module_pepc_types, only: t_particle, t_tree_node, kind_node
    use module_timings
    use module_tree, only: t_tree, tree_create, tree_check, tree_dump
    use module_domains, only: domain_decompose
    use module_debug
    use module_comm_env, only: t_comm_env, comm_env_dup
    use module_box, only: box_create
    use treevars, only: MPI_COMM_lpepc
    use module_spacefilling, only: compute_particle_keys
    implicit none
    include 'mpif.h'

    type(t_tree), intent(inout) :: t !< the tree
    type(t_particle), allocatable, intent(inout) :: p(:) !< input particle data, initializes %x, %data, %work appropriately (and optionally set %label) before calling this function

    type(t_particle) :: bp(2)
    type(t_tree_node), pointer :: root
    type(t_comm_env) :: tree_comm_env
    integer(kind_node), allocatable :: branch_nodes(:)
    integer(kind_particle) :: nl!< local number of simulation particles
    integer(kind_particle) :: n !< total number of simulation particles (across all MPI ranks)
    integer(kind_default) :: ierr

    call pepc_status('GROW TREE')
    !call MPI_BARRIER( MPI_COMM_lpepc, ierr)  ! Wait for everyone to catch up
    call timer_start(t_all)
    call timer_start(t_fields_tree)

    call comm_env_dup(MPI_COMM_lpepc, tree_comm_env)

    ! determine the bounding box of the particle configuration
    call box_create(t%bounding_box, p, tree_comm_env)
    call timer_start(t_domains_keys)
    ! assign SFC coordinate to each particle
    call compute_particle_keys(t%bounding_box, p)
    call timer_stop(t_domains_keys)
    
    ! determine local number of particles
    nl = size(p, kind=kind(nl))

    ! determine total number of particles
    call MPI_ALLREDUCE(nl, n, 1, MPI_KIND_PARTICLE, MPI_SUM, tree_comm_env%comm, ierr)

    ! Domain decomposition: allocate particle keys to PEs
    call domain_decompose(t%decomposition, t%bounding_box, n, p, bp, tree_comm_env)

    ! the local number of particles has changed
    nl = size(p, kind=kind(nl))

    ! allocate the tree structure
    call tree_create(t, nl, n, comm_env = tree_comm_env)

    ! build local part of tree
    call timer_start(t_local)
    call tree_build_from_particles(t, p, bp)

    root => t%nodes(t%node_root)

    if (root%leaves .ne. nl) then
      call tree_dump(t, p)
      DEBUG_ERROR(*, 'did not find all its particles inside the tree after local tree buildup: root_node%leaves =', root%leaves, ' but size(particles) =', nl)
    end if

    if (root%descendants+1 .ne. t%nleaf+t%ntwig) then
      call tree_dump(t, p)
      DEBUG_ERROR(*, '(root%descendants+1 .ne. t%nleaf+t%ntwig) after local tree buildup - some nodes got lost: root_node%descendants =', root%descendants, ' nleaf =', t%nleaf, ' ntwig =', t%ntwig)
    end if

    call timer_stop(t_local)

    if (.not. tree_check(t, "tree_grow: before exchange")) then
      call tree_dump(t, p)
    end if

    ! Should now have multipole information up to root list level(s) (only up to branch level, the information is correct)
    ! By definition, this is complete: each branch node is self-contained.
    ! This information has to be broadcast to the other PEs so that the top levels can be filled in.

    call timer_start(t_exchange_branches)
    call tree_exchange_branches(t, p, bp, branch_nodes)
    call timer_stop(t_exchange_branches)

    ! build global part of tree
    call timer_start(t_global)
    call tree_build_upwards(t, branch_nodes)
    call timer_stop(t_global)
    deallocate(branch_nodes)

    if (.not. tree_check(t, "tree_grow: after exchange")) then
      call tree_dump(t, p)
    end if

    if (root%leaves .ne. t%npart) then
      call tree_dump(t, p)
      DEBUG_ERROR(*, 'did not find all particles inside the htable after global tree buildup: root_node%leaves =', root%leaves, ' but npart_total =', t%npart)
    endif

    call timer_stop(t_fields_tree)
    call pepc_status('TREE GROWN')
  end subroutine tree_grow


  !>
  !> Exchanges tree nodes that are given in `local_branch_keys` with remote ranks.
  !>
  !> Incoming tree nodes are inserted into tree `t`, but the tree above these nodes is not corrected.
  !> Outputs keys of new (and own) tree nodes in `branch_keys`.
  !>
  subroutine tree_exchange(t, num_local_branch_nodes, local_branch_nodes, branch_nodes)
    use module_tree, only: t_tree, tree_insert_node
    use module_tree_node, only: tree_node_pack, tree_node_unpack, TREE_NODE_FLAG_GLOBAL_IS_BRANCH_NODE
    use module_pepc_types, only: t_tree_node, t_tree_node_package, MPI_TYPE_tree_node_package, kind_node
    use module_debug, only : pepc_status
    use module_timings
    use module_debug
    implicit none
    include 'mpif.h'

    type(t_tree), intent(inout) :: t !< the tree
    integer(kind_node), intent(in) :: num_local_branch_nodes !< number of local branch nodes (valid entries in local_branch_nodes(1:num_local_branch_nodes), everything beyond is garbage
    integer(kind_node), intent(in) :: local_branch_nodes(:) !< all local branch nodes
    integer(kind_node), intent(inout), allocatable :: branch_nodes(:) !< all branch nodes in the tree

    integer(kind_default) :: ierr
    integer(kind_pe) :: ip
    type(t_tree_node) :: unpack_node
    type(t_tree_node_package), allocatable :: pack_mult(:), get_mult(:)
    !> these have to be kind_default as MPI_ALLGATHERV expects default integer kind arguments
    integer(kind_default) :: i, j, nbranch, nbranch_sum
    integer(kind_default), allocatable :: nbranches(:), igap(:)

    call timer_start(t_exchange_branches_pack)
    
    nbranch = int(num_local_branch_nodes, kind=kind(nbranch)) ! we need this cast since MPI-strides have to be of kind_default

    ! Pack local branches for shipping
    allocate(pack_mult(nbranch))
    do i = 1, nbranch
      call tree_node_pack(t%nodes(local_branch_nodes(i)), pack_mult(i))
    end do

    call timer_stop(t_exchange_branches_pack)
    call timer_start(t_exchange_branches_admininstrative)

    ! work out stride lengths so that partial arrays placed sequentially in global array
    allocate (nbranches(t%comm_env%size), igap(t%comm_env%size + 1))
    call MPI_ALLGATHER(nbranch, 1, MPI_KIND_DEFAULT, nbranches, 1, MPI_KIND_DEFAULT, t%comm_env%comm, ierr)

    igap(1) = 0
    do ip = 2, t%comm_env%size + 1_kind_pe
      igap(ip) = igap(ip - 1) + nbranches(ip - 1)
    end do

    nbranch_sum = igap(t%comm_env%size + 1)

    allocate(get_mult(1:nbranch_sum), branch_nodes(1:nbranch_sum))

    call timer_stop(t_exchange_branches_admininstrative)
    call timer_start(t_exchange_branches_allgatherv)

    ! actually exchange the branch nodes
    call MPI_ALLGATHERV(pack_mult, nbranch, MPI_TYPE_tree_node_package, get_mult, nbranches, igap, MPI_TYPE_tree_node_package, &
      t%comm_env%comm, ierr)

    deallocate(pack_mult)
    deallocate(nbranches, igap)

    call timer_stop(t_exchange_branches_allgatherv)
    call timer_start(t_exchange_branches_integrate)

    ! Integrate remote branches into local tree
    j = 0
    do i = 1, nbranch_sum
      ! insert all remote branches into local data structures (this does *not* prepare the internal tree connections, but only copies multipole properties and creates the htable-entries)
      if (get_mult(i)%owner /= t%comm_env%rank) then
        call tree_node_unpack(get_mult(i), unpack_node)
        
        call tree_insert_node(t, unpack_node, branch_nodes(i))
      else
        j = j + 1
        branch_nodes(i) = local_branch_nodes(j)
      end if
    end do

    DEBUG_ASSERT_MSG(j == nbranch,*, j, num_local_branch_nodes)

    deallocate(get_mult)

    call timer_stop(t_exchange_branches_integrate)
  end subroutine tree_exchange


  !>
  !> identifies branch nodes in tree `t` and makes them available on all communication ranks.
  !>
  !> the keys of all boundary nodes are returned in `bk`.
  !>
  !> @todo Combining this routine with `tree_exchange()` would be more efficient.
  !>
  subroutine tree_exchange_branches(t, p, bp, bn)
    use module_tree, only: t_tree
    use module_pepc_types, only: t_particle, kind_node
    use module_debug, only: pepc_status
    use module_tree_node
    use module_timings
    implicit none

    type(t_tree), intent(inout) :: t !< the tree
    type(t_particle), intent(in) :: p(:) !< list of local particles
    type(t_particle), intent(in) :: bp(2) !< boundary particles
    integer(kind_node), allocatable, intent(inout) :: bn(:) !< node-indices of all branch nodes

    integer(kind_node) :: i, num_local_branch_nodes
    integer(kind_node), allocatable :: local_branch_nodes(:)
    type(t_tree_node), pointer :: branch

    ! identification of branch nodes
    call timer_start(t_branches_find)
    call find_local_branches(t, p, bp, num_local_branch_nodes, local_branch_nodes)
    t%nbranch_me = num_local_branch_nodes
    call timer_stop(t_branches_find)

    call pepc_status('EXCHANGE BRANCHES')

    call tree_exchange(t, num_local_branch_nodes, local_branch_nodes, bn)
    t%nbranch = size(bn, kind=kind(t%nbranch))

    deallocate(local_branch_nodes)

    do i = 1, t%nbranch
      branch => t%nodes(bn(i))

      ! flag all branch nodes for later identification
      branch%flags_global = ibset(branch%flags_global, TREE_NODE_FLAG_GLOBAL_IS_BRANCH_NODE)

      if (branch%owner /= t%comm_env%rank) then
        ! additionally, we mark all remote branches as remote nodes (this information is propagated upwards later)
        branch%flags_local = ibset(branch%flags_local, TREE_NODE_FLAG_LOCAL_HAS_REMOTE_CONTRIBUTIONS)
      end if
    end do

    contains

    !>
    !> identifies all `nb` local branch nodes in tree `t` and returns their
    !> node indices in `bn`.
    !>
    subroutine find_local_branches(t, p, bp, nb, bn)
      use module_spacefilling, only: shift_key_by_level
      use module_math_tools, only: bpi
      use treevars, only: idim, nlev
      use module_debug, only: pepc_status
      use module_tree, only: tree_traverse_to_key
      implicit none

      type(t_tree), intent(inout) :: t !< the tree
      type(t_particle), intent(in) :: p(:) !< list of local particles
      type(t_particle), intent(in) :: bp(2) !< boundary particles
      integer(kind_node), intent(out) :: nb !< number of local branch nodes
      integer(kind_node), allocatable, intent(out) :: bn(:) !< list of local branch node indices

      integer(kind_level) :: ilevel
      integer(kind_key) :: vld_llim, vld_rlim, L, D1, D2, pos, j, possible_branch
      integer(kind_level) :: branch_level(0:nlev), branch_level_D1(0:nlev), branch_level_D2(0:nlev)
      integer(kind_node) :: found_branch

      call pepc_status('FIND BRANCHES')
      call find_vld_limits(t, p, bp, vld_llim, vld_rlim)

      ! First find highest power in the Virtual Domain to ensure a correct branch definition
      L = bpi(vld_llim, vld_rlim)
      ! divide in two sub-domains
      ! only the last tasks must get 1 particle more
      ! because it s right limit is possibly not presentable
      D1 = L - vld_llim
      D2 = vld_rlim - L
      if (t%comm_env%last) then
        D2 = D2 + 1
      end if

      ! get estimate for number of local branches per level and total
      do ilevel = 0, nlev
        pos = idim * (nlev - ilevel)
        branch_level_D1(ilevel) = int(ibits(D1, pos, idim), kind(branch_level_D1))
        branch_level_D2(ilevel) = int(ibits(D2, pos, idim), kind(branch_level_D2))
        branch_level(ilevel) = branch_level_D1(ilevel) + branch_level_D2(ilevel)
      end do
      t%nbranch_max_me = sum(branch_level(:))
          
      allocate(bn(1:t%nbranch_max_me))
      
      nb = 0
      ! for D1
      pos = L
      do ilevel = 0, nlev
        do j = 1, branch_level_D1(ilevel)
          pos = pos - 2_kind_key**(idim * (nlev - ilevel))
          possible_branch = shift_key_by_level(pos, -(nlev - ilevel))
        
          ! After local build hashtable should contain branch key
          ! otherwise branch does not exists
          ! if entry exists it is counted as branch
          ! otherwise discarded
          if (tree_traverse_to_key(t, possible_branch, found_branch)) then ! entry exists
            nb = nb + 1
            bn(nb) = found_branch
          end if
        end do
      end do

      ! for D2
      pos = L - 1
      do ilevel = 0, nlev
        do j = 1, int(branch_level_D2(ilevel))
          pos = pos + 2_kind_key**(idim * (nlev - ilevel))
          possible_branch = shift_key_by_level(pos, -(nlev - ilevel))
          
          ! After build hashtable should contain branch key
          ! otherwise branch does not exists
          ! if entry exists it is counted as branch
          ! otherwise discarded
          if (tree_traverse_to_key(t, possible_branch, found_branch)) then ! entry exists
            nb = nb + 1
            bn(nb) = found_branch
          end if
        end do
      end do
    end subroutine find_local_branches


    !>
    !> find the virtual limits of the local key domain
    !>
    subroutine find_vld_limits(t, p, bp, l, r)
      use module_math_tools, only: bpi
      use treevars, only: idim, nlev
      implicit none

      type(t_tree), intent(in) :: t !< the tree
      type(t_particle), intent(in) :: p(:) !< list of local particles
      type(t_particle), intent(in) :: bp(2) !< boundary particles
      integer(kind_key), intent(out) :: l !< left virtual key domain limit
      integer(kind_key), intent(out) :: r !< right virtual key domain limit

      integer(kind_key) :: lme, rme, lb, rb
        
      ! get local key limits
      lme = p(1)%key
      rme = p(ubound(p, 1))%key

      ! get key limits for neighbor tasks
      ! and build virtual limits, so that a minimum set of branch nodes comes arround
      ! boundary tasks can access their boundary space fully only need one virtual limit
      l = 2_kind_key**(nlev * idim)
      r = 2_kind_key**(nlev * idim + 1) - 1

      if (.not. t%comm_env%first) then
        lb = bp(1)%key
        l = bpi(lb, lme)
      end if

      if (.not. t%comm_env%last) then
        rb = bp(2)%key
        r = bpi(rme, rb)
      end if
    end subroutine find_vld_limits
  end subroutine tree_exchange_branches


  !>
  !> Builds up the tree `t` from the given start keys `keys` towards root
  !>  - expects, that the nodes that correspond to `keys` already
  !>    have been inserted into the tree
  !>  - missing nodes on the way towards root are added automatically
  !>  - already existing nodes are updated
  !>
  subroutine tree_build_upwards(t, nodes)
    use module_tree, only: t_tree, tree_insert_node, &
      tree_node_connect_children
    use module_debug, only: pepc_status, DBG_TREE, dbg
    use module_tree_node, only: NODE_INVALID
    use module_timings
    use module_sort, only: sort
    use module_spacefilling, only: level_from_key, parent_key_from_key
    use module_pepc_types, only: t_tree_node, kind_node
    use module_debug
    implicit none

    type(t_tree), intent(inout) :: t !< the tree
    integer(kind_node), intent(in) :: nodes(:) !< nodes at which to start the build upwards

    integer(kind_level), allocatable :: key_level(:)
    integer(kind_particle), allocatable :: sort_map(:)
    integer(kind_key), allocatable :: sub_key(:), parent_key(:)
    integer(kind_node), allocatable :: sub_nodes(:), sorted_sub_nodes(:), parent_nodes(:)
    type(t_tree_node) :: parent_node

    integer(kind_node) :: i, nparent, nuniq, k
    integer(kind_node) ::  numnodes, nsub, groupstart, groupend
    integer(kind_level) :: ilevel, maxlevel
    integer(kind_key) :: current_parent_key

    call pepc_status('BUILD TOWARDS ROOT')

    numnodes = size(nodes, kind=kind(numnodes))
    allocate(key_level(numnodes))
    allocate(sub_key(0:numnodes + 1), sub_nodes(1:numnodes + 1), parent_key(0:numnodes + 1), parent_nodes(1:numnodes + 1))
    allocate(sorted_sub_nodes(1:numnodes + 1), sort_map(1:numnodes + 1))

    ! get levels of branch nodes
    do i = 1, numnodes
      key_level(i) = level_from_key(t%nodes(nodes(i))%key)
    end do
    maxlevel = maxval( key_level(:) ) ! Find maximum level

    nparent = 0
    ! iterate through key levels
    do ilevel = maxlevel,1,-1 ! Start at finest level
      ! Collect all keys at this level
      nsub = 0
      do i = 1, numnodes
        if (key_level(i) == ilevel) then
          nsub            = nsub + 1
          sub_key(nsub)   = t%nodes(nodes(i))%key
          sub_nodes(nsub) = nodes(i)
        end if
      end do

      ! Augment list with parent keys checked at previous level
      do i = 1, nparent
        sub_key(nsub + i)   = parent_key(i)
        sub_nodes(nsub + i) = parent_nodes(i)
      end do
      nsub = nsub + nparent

      call sort(sub_key(1:nsub), sort_map(1:nsub)) ! Sort keys

      do i = 1, nsub
        sorted_sub_nodes(i) = sub_nodes(sort_map(i))
      end do

      sub_key(0) = 0 ! remove all duplicates from the list
      nuniq = 0
      do i = 1, nsub
        if (sub_key(i) .ne. sub_key(i-1)) then
          nuniq          = nuniq + 1
          sub_key(nuniq) = sub_key(i)
          sorted_sub_nodes(nuniq) = sorted_sub_nodes(i)
        end if
      end do

      nsub = nuniq
      sub_key(nsub+1) = 0

      ! now, sub_key(1:nsub) contains a list of all keys (unique) at ilevel that
      ! (1) just have been inserted into the tree
      ! (2) have been modified due to additional child data
      ! tree_nodes() and htable()-entries exist for both cases
      ! --> their parents need to be created and/or updated
      i       = 1
      nparent = 0

      do while (i <= nsub)
        ! group keys with the same parent
        current_parent_key = parent_key_from_key(sub_key(i))

        groupstart = i
        do while ((parent_key_from_key(sub_key(i+1)) .eq. current_parent_key) .and. (i+1 <= nsub))
          i = i + 1
        end do
        groupend   = i

        nparent               = nparent + 1
        parent_key(nparent)   = current_parent_key
        parent_nodes(nparent) = NODE_INVALID
        
        do k=groupstart,groupend
          if (t%nodes(sorted_sub_nodes(k))%parent /= NODE_INVALID) then
            ! a parent node is already existing
            parent_nodes(nparent) = t%nodes(sorted_sub_nodes(k))%parent
            exit
          end if
        end do

        ! a node has to be created from scratch
        if (parent_nodes(nparent) == NODE_INVALID) then
          call tree_node_create_from_children(t, parent_node, sorted_sub_nodes(groupstart:groupend), current_parent_key)
          call tree_insert_node(t, parent_node, parent_nodes(nparent))
        else
          call tree_node_update_from_children(t, t%nodes(parent_nodes(nparent)), sorted_sub_nodes(groupstart:groupend), current_parent_key)
        end if
        call tree_node_connect_children(t, parent_nodes(nparent), sorted_sub_nodes(groupstart:groupend))

        ! go on with next group
        i = i + 1
      end do
    end do ! loop over levels

    deallocate(key_level)
    deallocate(sub_key, sub_nodes, parent_key, parent_nodes)
    deallocate(sorted_sub_nodes, sort_map)
  end subroutine tree_build_upwards


  !>
  !> Inserts all given local particles `p` at the correct leaf level in tree `t`
  !> by recursively subdividing the nodes if necessary.
  !>
  !> After function execution, leaf nodes exist for all particles in `p` and
  !> all twig nodes above those leaves until the root node.
  !> Data for leaves are completely valid while those
  !> for twigs have to be updated via a call to `tree_build_upwards()`.
  !>
  !> Upon exit, the `node_leaf` property of any particle is set appropriately.
  !>
  !> @warning In contrast to `tree_build_upwards()`, this function may *not*
  !> be called several times to add further particles etc.
  !>
  subroutine tree_build_from_particles(t, p, bp)
    use module_tree, only: t_tree, TREE_KEY_ROOT
    use module_tree_node, only: NODE_INVALID
    use treevars, only : idim, nlev
    use module_pepc_types, only: t_particle, t_tree_node
    use module_spacefilling, only: level_from_key
    use module_tree_node
    use module_timings
    use module_debug
    implicit none

    type(t_tree), intent(inout) :: t !< the tree
    type(t_particle), intent(inout) :: p(:) !< list of particles
    type(t_particle), intent(in) :: bp(2) !< boundary particles

    type t_keyidx
      integer(kind_key) :: key
      integer(kind_particle) :: idx
    end type

    type(t_keyidx), allocatable :: kidx(:)
    integer(kind_particle) :: i, j

    call pepc_status('INSERT PARTICLES')

    call timer_start(t_build_pure)

    allocate(kidx(size(p, kind = kind_particle) + 2)) ! might not need those two, but what the heck
    i = 0
    if (.not. t%comm_env%first) then
      i = i + 1
      kidx(i) = t_keyidx( bp(1)%key, 0 )
    end if

    do j = 1, size(p, kind = kind_particle)
      i = i + 1
      kidx(i) = t_keyidx( p(j)%key, j )
    end do

    if (.not. t%comm_env%last) then
      i = i + 1
      kidx(i) = t_keyidx( bp(2)%key, 0 )
    end if

    p(:)%node_leaf = NODE_INVALID

    call timer_reset(t_props_leaves)
    call insert_helper(t, TREE_KEY_ROOT, level_from_key(TREE_KEY_ROOT), kidx(1:i), t%node_root)

    deallocate(kidx)
    call timer_stop(t_build_pure)

    t%nleaf_me = t%nleaf !  Retain leaves and twigs belonging to local PE
    t%ntwig_me = t%ntwig

    ! check if we did not miss any particles
    if (any(p(:)%node_leaf == NODE_INVALID)) then
      DEBUG_WARNING_ALL(*, ' did not incorporate all particles into its tree')
    end if

    contains

    !>
    !> Inserts the keys in `ki` below the key `k`.
    !>
    recursive subroutine insert_helper(t, k, l, ki, pi)
      use module_tree, only: t_tree, tree_insert_node, tree_insert_node_at_index, &
        tree_provision_node, tree_node_connect_children
      use module_pepc_types, only: t_tree_node, kind_node
      use module_tree_node
      use module_spacefilling, only: child_key_from_parent_key, is_ancestor_of_particle
      implicit none

      type(t_tree), intent(inout) :: t !< tree in which to find the nodes
      integer(kind_key), intent(in) :: k !< key blow which the `ki` are inserted
      integer(kind_level), intent(in) :: l !< precomputed `level_from_key(k)`
      type(t_keyidx), intent(in) :: ki(:) !< keys (with particles) to be inserted
      integer(kind_node), optional, intent(out) :: pi

      type(t_tree_node) :: this_node
      integer(kind_node) :: inserted_node_idx
      integer(kind_node) :: child_nodes(8)
      integer(kind_key) :: childkey
      integer(kind_level) :: childlevel
      integer :: ichild, nchild
      integer(kind_node) :: si, ip, pstart, pend

      inserted_node_idx = NODE_INVALID
      
      si = size(ki, kind=kind(si))

      if (si < 1) then ! no particles below this key
        ! do nothing
      else if (si == 1) then ! we have arrived at a leaf
        if (ki(1)%idx /= 0) then ! boundary particles have idx == 0, do not really need to be inserted
          call tree_node_create_from_particle(t, this_node, p(ki(1)%idx), k)
          call tree_insert_node(t, this_node, inserted_node_idx)
          p(ki(1)%idx)%node_leaf = inserted_node_idx
        end if
      else ! more particles left, split twig
        if (l >= nlev) then ! no more levels left, cannot split
          DEBUG_WARNING_ALL('("Problem with tree: No more levels. Remaining particles 1..",I0,"  [i, key, label, x, y, z]:")', si)
          DEBUG_ERROR_NO_HEADER('(I6,x,O22.22,x,I16,g20.12,x,g20.12,x,g20.12,x)', ( ip, p(ki(ip)%idx)%key, p(ki(ip)%idx)%label, p(ki(ip)%idx)%x(1:3), ip=1,si ) )
        end if

        inserted_node_idx = tree_provision_node(t)

        nchild = 0
        pstart = lbound(ki, dim = 1, kind = kind(pstart))
        pend = pstart - 1 ! silence compiler warnings about pend being used without initialising
        childlevel = l + 1_kind_level
        do ichild = 0, 2**idim - 1
          childkey = child_key_from_parent_key(k, ichild)

          pend = pstart - 1
          do while (pend < ubound(ki, dim = 1, kind = kind(pend)))
            if (.not. is_ancestor_of_particle(ki(pend + 1)%key, childkey, childlevel)) then; exit; end if
            pend = pend + 1
          end do
          
          call insert_helper(t, childkey, childlevel, ki(pstart:pend), child_nodes(nchild + 1))
          if (child_nodes(nchild + 1) .ne. NODE_INVALID) then; nchild = nchild + 1; end if
          pstart = pend + 1
        end do

        DEBUG_ASSERT(nchild > 0)
        if (.not. pend == ubound(ki, dim = 1)) then
          DEBUG_WARNING_ALL('("Problem with tree: Could not distribute particles among children of ", I0, ". Remaining particles ",I0,"..",I0,"  [i, key, label, x, y, z]:")', k, pend + 1, ubound(ki, dim = 1) )
          DEBUG_ERROR_NO_HEADER('(I6,x,I22.22,x,I16,g20.12,x,g20.12,x,g20.12,x)', ( ip, p(ki(ip)%idx)%key, p(ki(ip)%idx)%label, p(ki(ip)%idx)%x(1:3), ip=pend + 1, ubound(ki, dim = 1) ) )
        end if
        call tree_node_create_from_children(t, this_node, child_nodes(1:nchild), k)
        call tree_insert_node_at_index(t, inserted_node_idx, this_node)
        ! wire up pointers
        call tree_node_connect_children(t, inserted_node_idx, child_nodes(1:nchild))
      end if

      if (present(pi)) then
        pi = inserted_node_idx
      end if
    end subroutine insert_helper
  end subroutine tree_build_from_particles


  !>
  !> populates a tree node with information from a single particle
  !>
  subroutine tree_node_create_from_particle(t, n, p, k)
    use module_pepc_types, only: t_particle, t_tree_node
    use module_tree, only: t_tree
    use module_tree_node
    use module_interaction_specific, only: multipole_from_particle
    use module_spacefilling, only: level_from_key
    use module_timings
    implicit none

    type(t_tree), intent(inout) :: t
    type(t_tree_node), intent(out) :: n
    type(t_particle), intent(in) :: p
    integer(kind_key), intent(in) :: k

    n%flags_global = 0
    n%flags_local  = ibset(n%flags_local, TREE_NODE_FLAG_LOCAL_HAS_LOCAL_CONTRIBUTIONS)
    n%request_posted = .false.
    n%owner        = t%comm_env%rank
    n%key          = k
    n%level        = level_from_key(k)
    n%leaves       = 1
    n%descendants  = 0
    n%parent       = NODE_INVALID
    n%first_child  = NODE_INVALID
    n%next_sibling = NODE_INVALID

    call timer_resume(t_props_leaves)
    call multipole_from_particle(p%x, p%data, n%interaction_data)
    call timer_stop(t_props_leaves)
  end subroutine tree_node_create_from_particle

  
  !>
  !> accumulates properties of child nodes to parent node.
  !> initialises links to NODE_INVALID.
  !>
  subroutine tree_node_create_from_children(t, parent, children, k)
    use module_tree, only: t_tree
    use module_pepc_types, only: t_tree_node, kind_node
    use module_tree_node, only: NODE_INVALID
    implicit none

    type(t_tree), intent(inout) :: t !< tree in which to find the nodes
    type(t_tree_node), intent(inout) :: parent !< parent node
    integer(kind_node), intent(in) :: children(:) !< child nodes
    integer(kind_key), intent(in) :: k !< node key

    parent%parent       = NODE_INVALID
    parent%next_sibling = NODE_INVALID
    parent%first_child  = NODE_INVALID 

    call tree_node_update_from_children(t, parent, children, k)
  end subroutine tree_node_create_from_children


  !>
  !> accumulates properties of child nodes to parent node.
  !> leaves links to node's relatives untouched
  !>
  subroutine tree_node_update_from_children(t, parent, children, k)
    use module_pepc_types, only: t_tree_node, kind_node
    use module_interaction_specific_types, only: t_tree_node_interaction_data
    use module_tree_node
    use module_interaction_specific, only : shift_multipoles_up
    use module_spacefilling, only: parent_key_from_key, level_from_key, child_number_from_key
    use module_debug
    use module_tree, only: t_tree
    implicit none

    type(t_tree), intent(inout) :: t !< tree in which to find the nodes
    type(t_tree_node), intent(inout) :: parent !< parent node
    integer(kind_node), intent(in) :: children(:) !< child nodes
    integer(kind_key), intent(in) :: k !< node key

    type(t_tree_node_interaction_data) :: interaction_data(8)
    integer(kind_node) :: nchild, i
    type(t_tree_node), pointer :: child

    nchild = size(children, kind=kind(nchild))

    parent%key          = k
    parent%request_posted = .false.
    parent%level        = level_from_key(k)
    parent%owner        = t%comm_env%rank

    parent%flags_global = 0
    parent%flags_local  = 0
    parent%leaves       = 0
    parent%descendants  = nchild

    ! Set children_HERE flag parent since we just built it from its children
    parent%flags_local = ibset(parent%flags_local, TREE_NODE_FLAG_LOCAL_CHILDREN_AVAILABLE)

    do i = 1, nchild
      child => t%nodes(children(i))
      interaction_data(i) = child%interaction_data

      ! parents of nodes with local contributions also contain local contributions
      if (btest(child%flags_local, TREE_NODE_FLAG_LOCAL_HAS_LOCAL_CONTRIBUTIONS)) then
        parent%flags_local = ibset(parent%flags_local, TREE_NODE_FLAG_LOCAL_HAS_LOCAL_CONTRIBUTIONS)
      end if

      ! parents of nodes with remote contributions also contain remote contributions
      if (btest(child%flags_local, TREE_NODE_FLAG_LOCAL_HAS_REMOTE_CONTRIBUTIONS)) then
        parent%flags_local = ibset(parent%flags_local, TREE_NODE_FLAG_LOCAL_HAS_REMOTE_CONTRIBUTIONS)
      end if

      ! parents of branch and fill nodes will also be fill nodes
      if (btest(child%flags_global, TREE_NODE_FLAG_GLOBAL_IS_FILL_NODE) .or. btest(child%flags_global, TREE_NODE_FLAG_GLOBAL_IS_BRANCH_NODE)) then
        parent%flags_global = ibset(parent%flags_global, TREE_NODE_FLAG_GLOBAL_IS_FILL_NODE)
      end if

      parent%leaves      = parent%leaves      + child%leaves
      parent%descendants = parent%descendants + child%descendants
    end do

    call shift_multipoles_up(parent%interaction_data, interaction_data(1:nchild))
  end subroutine tree_node_update_from_children
end module module_tree_grow
