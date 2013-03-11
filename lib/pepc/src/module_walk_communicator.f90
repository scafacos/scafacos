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
!>  Encapsulates helper functions that simplify communication during tree traversal
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_walk_communicator
  use treevars
  implicit none

  private

    integer, public :: walk_status
    integer*8, public :: comm_loop_iterations(3) !< number of comm loop iterations (total, sending, receiving)

    !> debug flags - cannot be modified at runtime due to performance reasons
    logical, public, parameter :: walk_comm_debug = .false.

    ! tags to be used in communication
    integer, public, parameter :: TAG_REQUEST_KEY    = 1257 !< message tag for walk communication: message for requesting child data for a certain key
    integer, public, parameter :: TAG_REQUESTED_DATA = 1258 !< message tag for walk communication: message that contains requested child data
    integer, public, parameter :: TAG_FINISHED_PE    = 1259 !< message tag for walk communication: message to rank 0, that the sender has finished with walking, i.e. its walk status == WALK_IAM_FINISHED
    integer, public, parameter :: TAG_FINISHED_ALL   = 1260 !< message tag for walk communication: message from rank 0, that all PEs are finished with walking, i.e. we can set walk_status := WALK_ALL_FINISHED
    integer, public, parameter :: mintag = TAG_REQUEST_KEY
    integer, public, parameter :: maxtag = TAG_FINISHED_ALL

    ! internal status
    integer, public, parameter :: WALK_STILL_RUNNING = 0 !< value of walk_status: there are still particles left to be processed
    integer, public, parameter :: WALK_IAM_FINISHED  = 1 !< value of walk_status: all particles on this PE have finished their walk
    integer, public, parameter :: WALK_ALL_MSG_DONE  = 2 !< value of walk_status: there are no more pending requests or answers in the message queue
    integer, public, parameter :: WALK_I_NOTIFIED_0  = 3 !< value of walk_status: rank 0 has been informed that this PE is finished
    integer, public, parameter :: WALK_ALL_FINISHED  = 4 !< value of walk_status: rank 0 told this pe, that every PE is finished with walking, i.e. no more communication is necessary

    ! internal communication variables - not to be touched from outside the module
    integer*1, private, allocatable, target :: bsend_buffer(:) !< buffer for bsend-alls
    integer, private :: comm_dummy = 123456 !< dummy variable for sending "empty" messages (those, where we are only interested in the tag)
    real*8, public :: timings_comm(3) !< array for storing internal timing information

    ! IDs for internal timing measurement
    integer, public, parameter :: TIMING_COMMLOOP = 1
    integer, public, parameter :: TIMING_RECEIVE  = 2
    integer, public, parameter :: TIMING_SENDREQS = 3

    !> data type for internal request queue
    type, public :: t_request_queue_entry
      integer*8 :: key
      integer   :: addr
      integer   :: owner
    end type


    public init_comm_data
    public uninit_comm_data
    public notify_walk_finished
    public send_request
    public send_data
    public unpack_data
    public send_walk_finished
    public broadcast_walk_finished

  contains


      ! initializes bsend buffer and rwlock objects
      ! returns size of bsend buffer in buffsize in bytes
      subroutine init_comm_data(REQUEST_QUEUE_LENGTH, ANSWER_BUFF_LENGTH)
        use module_pepc_types
        implicit none
        include 'mpif.h'
        integer, intent(in) :: REQUEST_QUEUE_LENGTH, ANSWER_BUFF_LENGTH
        integer :: msg_size_request, msg_size_data
        integer :: buffsize !< size of bsend buffer in bytes
        integer :: ierr

        ! compute upper bounds for request and data message size
        call MPI_PACK_SIZE(1, MPI_INTEGER8, MPI_COMM_lpepc, msg_size_request, ierr)
        msg_size_request = msg_size_request + MPI_BSEND_OVERHEAD

        call MPI_PACK_SIZE(8, MPI_TYPE_tree_node_transport_package, MPI_COMM_lpepc, msg_size_data, ierr)
        msg_size_data = msg_size_data + MPI_BSEND_OVERHEAD

        buffsize = (REQUEST_QUEUE_LENGTH * msg_size_request + ANSWER_BUFF_LENGTH * msg_size_data)
        ! reserve memory for buffered mpi communication
        allocate(bsend_buffer(buffsize))
        ! and tell mpi , where it can be found
        call MPI_BUFFER_ATTACH(bsend_buffer, buffsize, ierr)

        walk_status = WALK_STILL_RUNNING

        timings_comm = 0.

    end subroutine init_comm_data



      subroutine uninit_comm_data()
        implicit none
        include 'mpif.h'
        integer :: ierr
        integer :: buffsize
        integer :: dummy

        ! free our buffer that was reserved for buffered communication
        call MPI_BUFFER_DETACH(dummy, buffsize, ierr) ! FIXME: what is the dummy thought for?
        deallocate(bsend_buffer)

      end subroutine uninit_comm_data



      subroutine notify_walk_finished()
        use treevars
        implicit none

        walk_status = max(walk_status, WALK_IAM_FINISHED)

      end subroutine notify_walk_finished





      subroutine send_walk_finished()
        use treevars
        implicit none
        include 'mpif.h'
        integer :: ierr

        ! notify rank 0 that we are finished with our walk
        call MPI_BSEND(comm_dummy, 1, MPI_INTEGER, 0, TAG_FINISHED_PE, &
          MPI_COMM_lpepc, ierr)

        walk_status = WALK_I_NOTIFIED_0

      end subroutine send_walk_finished


      subroutine broadcast_walk_finished()
        use treevars, only : num_pe, MPI_COMM_lpepc
        use module_debug
        implicit none
        include 'mpif.h'
        integer :: i, ierr

         ! all PEs have to be informed
         ! TODO: need better idea here...
         if (walk_comm_debug) then
           DEBUG_INFO('("PE", I6, " has found out that all PEs have finished walking - telling them to exit now")', me)
         end if

         do i=0,num_pe-1
           call MPI_BSEND(comm_dummy, 1, MPI_INTEGER, i, TAG_FINISHED_ALL, &
             MPI_COMM_lpepc, ierr)
         end do


      end subroutine


    subroutine send_data(requested_key, ipe_sender)
      use module_pepc_types
      use module_debug
      use module_htable
      implicit none
      include 'mpif.h'
      integer*8, intent(in) :: requested_key
      integer, intent(in) :: ipe_sender
      integer :: process_addr
      type(t_tree_node_transport_package), target :: children_to_send(8)
      type(t_tree_node_transport_package), pointer :: c
      integer*8, dimension(8) :: key_child
      integer, dimension(8) :: addr_child, node_child, byte_child, leaves_child, owner_child, level_child
      integer :: j, ic, ierr, nchild

      if (walk_comm_debug) then
        DEBUG_INFO('("PE", I6, " answering request.                         request_key=", O22, ",        sender=", I6)',
                       me, requested_key, ipe_sender )
      end if

      j = 0
      process_addr = key2addr( requested_key,'WALK:send_data:parentkey')       ! get htable addresses

      call get_childkeys(process_addr, nchild, key_child)

      addr_child(1:nchild)   = (/( key2addr( key_child(j),'WALK:send_data:childkey' ),j=1,nchild)/)  ! Table address of children
      node_child(1:nchild)   = htable( addr_child(1:nchild) )%node                        ! Child node index
      byte_child(1:nchild)   = int(IAND( htable( addr_child(1:nchild) )%childcode, CHILDCODE_CHILDBYTE ))! Catch lowest 8 bits of childbyte - filter off requested and here flags
      leaves_child(1:nchild) = htable( addr_child(1:nchild) )%leaves                      ! # contained leaves
      owner_child(1:nchild)  = htable( addr_child(1:nchild) )%owner                       ! real owner of child (does not necessarily have to be identical to me, at least after futural modifications)
      level_child(1:nchild)  = htable( addr_child(1:nchild) )%level                       ! level of child node
      ! Package children properties into user-defined multipole array for shipping
      do ic = 1,nchild
         c=>children_to_send(ic)
           c%m      = tree_nodes(node_child(ic))
           c%key    = key_child(ic)
           c%byte   = byte_child(ic)
           c%leaves = leaves_child(ic)
           c%owner  = owner_child(ic)
	   c%level  = level_child(ic)
      end do
      
      ! Ship child data back to PE that requested it
      call MPI_BSEND( children_to_send(1:nchild), nchild, MPI_TYPE_tree_node_transport_package, &
        ipe_sender, TAG_REQUESTED_DATA, MPI_COMM_lpepc, ierr)

      ! statistics on number of sent children-packages
      sum_ships = sum_ships + 1

    end subroutine send_data


    function send_request(req)
      use module_htable
      implicit none
      include 'mpif.h'
      logical :: send_request
      type(t_request_queue_entry), intent(in) :: req
      integer :: ierr

      if (.not. BTEST( htable(req%addr)%childcode, CHILDCODE_BIT_REQUEST_SENT ) ) then
        ! send a request to PE req_queue_owners(req_queue_top)
        ! telling, that we need child data for particle request_key(req_queue_top)
        call MPI_BSEND(req%key, 1, MPI_INTEGER8, req%owner, TAG_REQUEST_KEY, &
          MPI_COMM_lpepc, ierr)

        htable(req%addr)%childcode = ibset(htable(req%addr)%childcode, CHILDCODE_BIT_REQUEST_SENT )

        send_request = .true.
      else
        send_request = .false.
      end if

    end function



    subroutine unpack_data(child_data, num_children, ipe_sender)
      use module_htable
      use module_tree
      use module_spacefilling
      use module_debug
      use module_atomic_ops
      implicit none
      include 'mpif.h'
      type (t_tree_node_transport_package) :: child_data(num_children) !< child data that has been received
      integer :: num_children !< actual number of valid children in dataset
      integer, intent(in) :: ipe_sender
      integer*8 :: kparent
      integer :: ownerchild, parent_addr(0:num_children), num_parents
      integer :: ic

      num_parents = 0
      parent_addr(0) = maxaddress + 1

      do ic = 1, num_children
        ! save parent address - after (!) inserting all (!) children we can flag it: it`s children are then accessible
        kparent     = parent_key_from_key( child_data(ic)%key )
        parent_addr(num_parents + 1) = key2addr( kparent, 'WALK:unpack_data() - get parent address' )
        if (parent_addr(num_parents) .ne. parent_addr(num_parents + 1)) then
          num_parents = num_parents + 1
        end if

        if (walk_comm_debug) then
          DEBUG_INFO('("PE", I6, " received answer.                            parent_key=", O22, ",        sender=", I6, ",        owner=", I6, ", kchild=", O22)',
                         me, kparent, ipe_sender, ownerchild, child_data(ic)%key)
        end if

        ! tree nodes coming from remote PEs are flagged for easier identification
        child_data(ic)%byte = ibset(child_data(ic)%byte, CHILDCODE_BIT_HAS_REMOTE_CONTRIBUTIONS)

        call tree_insert_node(child_data(ic))
        ! count number of fetched nodes
        sum_fetches = sum_fetches+1
     end do

     call atomic_read_write_barrier()

     ! set 'children-here'-flag for all parent addresses
     ! may only be done *after inserting all* children, hence not(!) during the loop above
     do ic=1,num_parents
         htable( parent_addr(ic) )%childcode = IBSET(  htable( parent_addr(ic) )%childcode, CHILDCODE_BIT_CHILDREN_AVAILABLE) ! Set children_HERE flag for parent node
     end do


    end subroutine unpack_data

end module module_walk_communicator
