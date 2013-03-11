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
!> Contains all tree specific helper routines and data fields
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_tree
    implicit none
    private

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! TODO: move tree_nodes array here

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    public tree_insert_node
    public tree_exchange
    public tree_build_upwards
    public tree_build_from_particles

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !> type for storing key and nideindex together in tree_build_from_particles
    type t_keyidx
      integer :: idx
      integer*8 :: key
    end type

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  subroutine-implementation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    contains


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> if an entry with tree_node%key already exists in htable, then
    !> updates htable and tree_nodes with new data
    !> otherwise: creates new entries
    !>
    !> this routine cannot be used to change a tree_node from leaf to twig or similar
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine tree_update_or_insert_node(tree_node)
        use treevars
        use module_pepc_types
        use module_htable
        implicit none
        include 'mpif.h'

        type(t_tree_node_transport_package), intent(in) :: tree_node
        integer :: addr

        if (testaddr(tree_node%key, addr)) then
          ! the htable-entry and node already exist --> update

          ! if we change the owner from someting else to 'me', we have to keep track of the leaf/twig counters
          if ((htable(addr)%owner .ne. me) .and. (tree_node%owner .eq. me)) then
            if (htable_entry_is_leaf(addr)) then
              nleaf_me = nleaf_me + 1
            else
              ntwig_me = ntwig_me + 1
            endif
          endif

          htable( addr )%leaves    = tree_node%leaves
          htable( addr )%childcode = tree_node%byte
          htable( addr )%owner     = tree_node%owner

          tree_nodes(htable(addr)%node) = tree_node%m
        else
          ! create new htable and nodelist entry
          call tree_insert_node(tree_node)
        endif

    end subroutine


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Inserts a given tree node into the next free position in the tree ( -(ntwig+1) or (nleaf+1) )
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine tree_insert_node(tree_node)
        use treevars
        use module_pepc_types
        use module_htable
        use module_debug
        implicit none

        type(t_tree_node_transport_package), intent(in) :: tree_node
        integer :: hashaddr, lnode


        if (make_hashentry( t_hash(tree_node%key, 0, -1, tree_node%leaves, tree_node%byte, tree_node%owner, tree_node%level), hashaddr)) then
            ! anything is fine - we will have to assign a node number now
            if ( tree_node%leaves == 1 ) then
               nleaf =  nleaf + 1
               lnode =  nleaf
               if (tree_node%owner == me) nleaf_me = nleaf_me+1
            else if ( tree_node%leaves > 1 ) then
               ! twig
               ntwig =  ntwig + 1
               lnode = -ntwig
               if (tree_node%owner == me) ntwig_me = ntwig_me+1
            else
               DEBUG_ERROR(*, "Found a tree node with less than 1 leaf.")
            endif

            ! check for array bound overrun
            if ((ntwig >= maxtwig) .or. (nleaf >= maxleaf)) then
              DEBUG_ERROR('("Tree arrays full. Twigs: ", I0,"/",I0 ,"; Leaves: ", I0,"/",I0)', ntwig, maxtwig, nleaf, maxleaf)
            end if

             htable( hashaddr )%node = lnode

        else
           ! entry with the same key is already existing, so we just overwrite it
           lnode = htable( hashaddr )%node

           DEBUG_WARNING_ALL(*, "PE", me, "has found an already inserted entry while calling make_hashentry(", tree_node%key, lnode, tree_node%leaves, tree_node%byte, tree_node%owner, tree_node%level, hashaddr,") - overwriting it")
        endif

        !insert multipole data into local tree
        tree_nodes( lnode ) = tree_node%m

    end subroutine tree_insert_node

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Accumulates properties of child nodes (given by keys) to parent node
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine shift_nodes_up_key(parent, childkeys, parent_owner)
      use module_pepc_types
      use treevars, only : tree_nodes
      use module_htable
      use module_spacefilling
      implicit none
      type(t_tree_node_transport_package), intent(inout) :: parent
      integer*8, intent(in) :: childkeys(:)
      integer, intent(in) :: parent_owner

      integer :: nchild, i
      type(t_tree_node_transport_package) :: child_nodes(1:8)
      integer :: child_addr, childnumber(1:8)

      nchild = size(childkeys)

      do i=1,nchild
        child_addr     = key2addr(childkeys(i), 'shift_nodes_up_key')
        childnumber(i) = child_number_from_key(childkeys(i))
        child_nodes(i) = t_tree_node_transport_package(childkeys(i),                   &
                                     htable( child_addr )%childcode, &
                                     htable( child_addr )%leaves,    &
                                     htable( child_addr )%owner,     &
                                     htable( child_addr )%level,     &
                         tree_nodes( htable( child_addr )%node ) )
      end do

      call shift_nodes_up(parent, child_nodes(1:nchild), childnumber(1:nchild), parent_owner)

    end subroutine


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Accumulates properties of child nodes to parent node
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine shift_nodes_up(parent, children, childnumber, parent_owner)
      use module_pepc_types
      use module_htable, only : CHILDCODE_BIT_CHILDREN_AVAILABLE
      use module_interaction_specific, only : shift_multipoles_up
      use module_spacefilling
      use module_debug
      use module_htable
      implicit none
        type(t_tree_node_transport_package), intent(inout) :: parent
        type(t_tree_node_transport_package), intent(in) :: children(:)
        integer, intent(in) :: childnumber(:)
        integer, intent(in) :: parent_owner
        integer*8 :: parent_keys(1:8)

        integer :: nchild, i, byte

        nchild = size(children)

        ! check if all keys fit to the same parent
        parent_keys(1:nchild) = parent_key_from_key(children(1:nchild)%key)

        if ( any(parent_keys(2:nchild) .ne. parent_keys(1))) then
          DEBUG_ERROR(*,"Error in shift nodes up: not all supplied children contribute to the same parent node")
        endif

        byte = 0
        do i = 1,nchild
          ! set bits for available children
          byte = IBSET(byte, childnumber(i))
          ! parents of nodes with local contributions also contain local contributions
          if (btest(children(i)%byte, CHILDCODE_BIT_HAS_LOCAL_CONTRIBUTIONS)) byte = ibset(byte, CHILDCODE_BIT_HAS_LOCAL_CONTRIBUTIONS)
          ! parents of nodes with remote contributions also contain remote contributions
          if (btest(children(i)%byte, CHILDCODE_BIT_HAS_REMOTE_CONTRIBUTIONS)) byte = ibset(byte, CHILDCODE_BIT_HAS_REMOTE_CONTRIBUTIONS)
          ! parents of branch and fill nodes will also be fill nodes
          if (btest(children(i)%byte, CHILDCODE_BIT_IS_FILL_NODE) .or. btest(children(i)%byte, CHILDCODE_BIT_IS_BRANCH_NODE)) byte = ibset(byte, CHILDCODE_BIT_IS_FILL_NODE)
        end do


        ! Set children_HERE flag parent since we just built it from its children
        byte =  IBSET( byte, CHILDCODE_BIT_CHILDREN_AVAILABLE )

        parent%key    = parent_keys(1)
        parent%byte   = byte
        parent%leaves = sum(children(1:nchild)%leaves)
        parent%owner  = parent_owner
	parent%level  = level_from_key( parent_keys(1) )

        call shift_multipoles_up(parent%m, children(1:nchild)%m)

    end subroutine



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Exchanges tree nodes that are given in local_branch_keys with remote PEs
    !> incoming tree nodes are inserted into tree_nodes array and htable, but the
    !> tree above these nodes is not corrected
    !> outputs keys of new(and own) htable/tree_node entries in branch_keys
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine tree_exchange(local_branch_keys, nbranch, branch_keys, nbranch_sum)

        use treevars, only : me, num_pe, tree_nodes, nbranches, MPI_COMM_lpepc
        use module_debug, only : pepc_status
        use module_pepc_types
        use module_timings
        use module_htable
        implicit none
        include 'mpif.h'

        integer*8, intent(in) :: local_branch_keys(1:nbranch)
        integer, intent(in) :: nbranch
        integer*8, intent(inout), allocatable :: branch_keys(:)
        integer, intent(out) :: nbranch_sum

        integer :: i,ierr
        type( t_hash ), pointer :: hbranch
        type (t_tree_node_transport_package),allocatable :: pack_mult(:), get_mult(:)
        integer, allocatable :: igap(:)    !  stride lengths of local branch arrays

        if (allocated(branch_keys)) deallocate(branch_keys)

        call timer_start(t_exchange_branches_pack)

        call pepc_status('EXCHANGE BRANCHES')

        ! Pack local branches for shipping
        allocate(pack_mult(nbranch))
        do i=1,nbranch
            hbranch      => htable( key2addr( local_branch_keys(i),'EXCHANGE: info' ) )

            ! additionally, we mark all local branches as branches since this is only done for remote branches during unpack (is used for fill node identification)
            hbranch%childcode = ibset(hbranch%childcode, CHILDCODE_BIT_IS_BRANCH_NODE)

            pack_mult(i) =  t_tree_node_transport_package( local_branch_keys(i), hbranch%childcode, hbranch%leaves, me, hbranch%level, tree_nodes(hbranch%node) )
        end do

        call timer_stop(t_exchange_branches_pack)
        call timer_start(t_exchange_branches_admininstrative)

        call mpi_allgather( nbranch, 1, MPI_INTEGER, nbranches, 1, MPI_INTEGER, MPI_COMM_lpepc, ierr )

        ! work out stride lengths so that partial arrays placed sequentially in global array
        allocate (igap(num_pe+3))

        igap(1) = 0
        do i=2,num_pe+1
            igap(i) = igap(i-1) + nbranches(i-1)
        end do

        nbranch_sum = igap(num_pe+1)

        allocate(get_mult(1:nbranch_sum), branch_keys(1:nbranch_sum))

        call timer_stop(t_exchange_branches_admininstrative)
        call timer_start(t_exchange_branches_allgatherv)

        ! actually exchange the branch nodes
        call MPI_ALLGATHERV(pack_mult, nbranch, MPI_TYPE_tree_node_transport_package, get_mult, nbranches, igap, MPI_TYPE_tree_node_transport_package, MPI_COMM_lpepc, ierr)

        deallocate(pack_mult)
        deallocate (igap)

        call timer_stop(t_exchange_branches_allgatherv)
        call timer_start(t_exchange_branches_integrate)

        ! Integrate remote branches into local tree
        do i = 1,nbranch_sum

            ! insert all remote branches into local data structures (this does *not* prepare the internal tree connections, but only copies multipole properties and creates the htable-entries)
            if (get_mult(i)%owner /= me) then
              ! delete all custom flags from incoming nodes (e.g. CHILDCODE_BIT_CHILDREN_AVAILABLE)
              get_mult(i)%byte = IAND(get_mult(i)%byte, CHILDCODE_CHILDBYTE)
              ! after clearing all bits we have to set the flag for branches again to propagate this property upwards during global buildup
              get_mult(i)%byte = ibset(get_mult(i)%byte, CHILDCODE_BIT_IS_BRANCH_NODE)
              ! additionally, we mark all remote branches as remote nodes (this information is propagated upwards later)
              get_mult(i)%byte = ibset(get_mult(i)%byte, CHILDCODE_BIT_HAS_REMOTE_CONTRIBUTIONS)

              call tree_insert_node(get_mult(i))
            endif
            ! store branch key for later (global tree buildup)
            branch_keys(i) = get_mult(i)%key
        end do
	
        deallocate(get_mult)

        call timer_stop(t_exchange_branches_integrate)

    end subroutine tree_exchange


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Builds up the tree from the given start keys towards root
    !>  - expects, that the nodes that correspond to start_keys already
    !>    have been inserted into htable and tree_nodes array
    !>  - missing nodes on the way towards root are added automatically
    !>  - already existing nodes are updated
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine tree_build_upwards(start_keys, numkeys)

        use treevars, only : me
        use module_debug, only : pepc_status, DBG_TREE, dbg
        use module_timings
        use module_utils
        use module_htable
        use module_spacefilling
        use module_interaction_specific
        use module_pepc_types
        implicit none

        integer*8, intent(in) :: start_keys(1:numkeys)
        integer, intent(in) :: numkeys
        integer, dimension(1:numkeys) :: branch_level
        integer*8, dimension(0:numkeys+1) :: sub_key, parent_key
        type(t_tree_node_transport_package) :: parent_node

        integer :: ilevel, maxlevel, nsub, groupstart, groupend, i, nparent, nuniq
        integer*8 :: current_parent_key

        call pepc_status('BUILD TOWARDS ROOT')

        if (dbg(DBG_TREE)) call check_table('after make_branches ')

        ! get levels of branch nodes
        branch_level(1:numkeys) = level_from_key(start_keys(1:numkeys))
        maxlevel = maxval( branch_level(1:numkeys) )        ! Find maximum level

        nparent = 0

        ! iterate through branch levels
        do ilevel = maxlevel,1,-1                                            ! Start at finest branch level
            ! Collect all branches at this level
            nsub = 0
            do i=1,numkeys
                if (branch_level(i) == ilevel) then
                    nsub          = nsub + 1
                    sub_key(nsub) = start_keys(i)
                endif
            end do

            ! Augment list with parent keys checked at previous level
            sub_key(nsub+1:nsub+nparent) = parent_key(1:nparent)
            nsub                         = nsub + nparent

            call sort(sub_key(1:nsub))                                        ! Sort keys

            sub_key(0)   = 0                                                  ! remove all duplicates from the list
            nuniq = 0
            do i=1,nsub
                if (sub_key(i) .ne. sub_key(i-1)) then
                    nuniq          = nuniq + 1
                    sub_key(nuniq) = sub_key(i)
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

              call shift_nodes_up_key(parent_node, sub_key(groupstart:groupend), me)
              call tree_update_or_insert_node(parent_node)

              nparent             = nparent + 1
              parent_key(nparent) = current_parent_key

              ! go on with next group
              i = i + 1
            end do

        end do

    end subroutine tree_build_upwards


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> clears the htable and inserts all given particles at the correct level
    !> by recursively subdividing the cells if necessary
    !> after function execution, htable- and tree_node-entries for all twig- and
    !> leaf-keys exist. entries for leaves are completely valid while those
    !> for twigs have to be updated via a call to tree_build_upwards()
    !>
    !> upon exit, the key_leaf property of any particle_list entry
    !> is set appropriately
    !>
    !> warning: in contrast to tree_build_upwards(), this function may *not*
    !> be called several times to add further particles etc.
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine tree_build_from_particles(particle_list, nparticles, leaf_keys)
      use treevars, only : nleaf, ntwig, nlev, me, tree_nodes, nleaf_me, ntwig_me
      use module_pepc_types
      use module_spacefilling
      use module_htable
      use module_interaction_specific, only : multipole_from_particle
      use module_timings
      use module_debug
      implicit none
      type(t_particle), intent(inout) :: particle_list(1:nparticles)
      integer, intent(in) :: nparticles
      integer*8, intent(out) :: leaf_keys(1:nparticles)

      type(t_keyidx) :: particles_left(1:2*nparticles) ! each particle might produce one twig --> we need 2*nparticles as storage space
      integer :: i, k, nremaining, nreinserted, level, ibit, hashaddr, nremoved
      integer*8 :: lvlkey

      call pepc_status('INSERT PARTICLES')

      leaf_keys(1:nparticles) = 0_8

      nremaining = nparticles

      do i=1,nparticles
        particles_left(i)%key = particle_list(i)%key
        particles_left(i)%idx = i
      end do

      ! The following code works as follows:
      ! - starting from coarsest level, each particle's key on
      !   that level is computed and inserted into the htable as a leaf
      ! - in case of a collision, the respective htable-entry is turned
      !   into a twig and the former leaf as well as the additional
      !   particle are put onto a list for later processing at the next level
      !
      ! For simplicity, this code makes heavy use of a correspondence between
      ! the index in the particles(:)-array and the node list. This leads to
      ! the construction, that
      !     htable(key2addr(particles(i)%key))%node == i
      ! which is also desirable for later access
      call timer_start(t_build_pure)

      nleaf = nparticles
      level = 0

      do while (nremaining > 0)

        level = level + 1
        nreinserted = 0
        nremoved    = 0

        if (level>nlev) then
           DEBUG_WARNING_ALL('("Problem with tree: No more levels. Remaining particles 1..",I0,"  [i, local index, key, label, x, y, z]:")', nremaining)
           DEBUG_ERROR_NO_HEADER('(2(I6,x),O22.22,x,I16,g20.12,x,g20.12,x,g20.12,x)' , (i, particles_left(i)%idx, particle_list(particles_left(i)%idx)%key, particle_list(particles_left(i)%idx)%label, particle_list(particles_left(i)%idx)%x(1:3), i=1,nremaining ) )
         endif

         ibit = nlev - level ! bit shift factor (0=highest leaf node, nlev-1 = root)

         ! Determine subcell # from key
         ! At a given level, these will be unique
         do i=1,nremaining
           lvlkey = shift_key_by_level( particles_left(i)%key, - ibit )

                                             ! V nodeindex for leaves is identical to original particle index
           if (make_hashentry( t_hash(lvlkey, particles_left(i)%idx, -1, 1, 0, me, level), hashaddr)) then
             ! this key does not exist until now --> has been inserted as leaf
             leaf_keys(particles_left(i)%idx) = lvlkey
             htable(hashaddr)%childcode = ibset(htable(hashaddr)%childcode, CHILDCODE_BIT_HAS_LOCAL_CONTRIBUTIONS) ! we mark this node as having local contributions since it is a leaf, i.e. a node for some local particle
             particles_left(i) = t_keyidx(0, 0_8)
             nremoved = nremoved + 1
           else
             ! the key already exists
             ! do not remove particle from list - we will retry on next level

             ! if current entry at hashaddr is a leaf
             if (htable(hashaddr)%node > 0) then
               ! put it onto our list of unfinished particles again
               nreinserted                              = nreinserted + 1
               particles_left(nremaining + nreinserted) = t_keyidx(htable(hashaddr)%node, particle_list(htable(hashaddr)%node)%key)
               ! remove this entry from leaf_keys array
               leaf_keys(htable(hashaddr)%node) = 0_8
               ! and turn the current entry into a twig
               ntwig                 =  ntwig + 1
               htable(hashaddr)%node = -ntwig
             else
               ! the entry already was a twig --> nothing to do
             endif
           endif

         end do

         ! compact list of remaining particles
         k = 0
         do i=1,nremaining + nreinserted
           if (particles_left(i)%key .ne. 0) then
             k = k + 1
             particles_left(k) = particles_left(i)
           endif
         end do

         nremaining = nremaining + nreinserted - nremoved

         if (nremaining .ne. k) then
           DEBUG_ERROR(*,"Error in bookkeeping in local tree buildup. nremaining = ", nremaining, " but k = ", k)
         endif
      end do

      call timer_stop(t_build_pure)
      call timer_start(t_props_leafs)

      ! now we can use the correspondence between particle list index and tree_node index for setting the multipole properties
      do i=1,nparticles
        call multipole_from_particle(particle_list(i)%x, particle_list(i)%data, tree_nodes(i) )
      end do

      call timer_stop(t_props_leafs)

      nleaf_me = nleaf       !  Retain leaves and twigs belonging to local PE
      ntwig_me = ntwig

      ! copy leaf keys to particle datafield
      particle_list(1:nparticles)%key_leaf = leaf_keys(1:nparticles)

      ! check if we did not miss any particles
      if (any(leaf_keys(1:nparticles) == 0_8)) then
        DEBUG_WARNING_ALL(*, ' did not incorporate all particles into its leaf_keys array')
      endif

    end subroutine




end module module_tree
