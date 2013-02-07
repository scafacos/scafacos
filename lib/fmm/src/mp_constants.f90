module mympi_constants
use iso_c_binding, only: c_int,c_double
implicit none
include 'mpif.h'

	integer, parameter :: MyMPI_Comm = kind(MPI_COMM_NULL)
	integer, parameter :: MyMPI_Op = kind(MPI_OP_NULL)
	integer, parameter :: MyMPI_Datatype = kind(MPI_DATATYPE_NULL)
	integer, parameter :: MyMPI_Errorcode = kind(MPI_DATATYPE_NULL)
	integer, parameter :: MyMPI_Entries = kind(c_int) 
	integer, parameter :: MyMPI_Initialized = kind(c_int) !8
	integer, parameter :: MyMPI_Rank = kind(c_int) 
	integer, parameter :: MyMPI_Timer = kind(c_double) 
end module mympi_constants

module myarmci_constants
use iso_c_binding
implicit none
	integer, parameter :: MyARMCI_Errorcode = c_int
	integer, parameter :: MyARMCI_Sendsize = c_int
	integer, parameter :: MyARMCI_Proc = c_int
	integer, parameter :: MyARMCI_Waitcount = c_int

	integer(MyARMCI_Errorcode), parameter :: MyARMCI_SUCCESS = 0
end module myarmci_constants
