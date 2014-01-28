cimport numpy as np

# expose internal C-types to Cython
# MPI stuff
cdef extern from "mpi.h":
  cdef struct MPI_Comm:
    pass
  cdef MPI_Comm MPI_COMM_WORLD

ctypedef int fcs_int
ctypedef double fcs_float
  
# FCSResult
cdef extern from "fcs.h":
  cdef int FCS_SUCCESS
  cdef int FCS_RESULT_SUCCESS

  cdef struct FCSResult_t:
    pass
  ctypedef FCSResult_t* FCSResult
  fcs_int fcsResult_getReturnCode(FCSResult err)
  const char* fcsResult_getErrorMessage(FCSResult err)
  const char* fcsResult_getErrorSource(FCSResult err)

cdef handleResult(FCSResult result):
    rc = fcsResult_getReturnCode(result)
    if rc == -1 or rc == FCS_SUCCESS or rc == FCS_RESULT_SUCCESS:
        return
    else:
        source = fcsResult_getErrorSource(result)
        message = fcsResult_getErrorMessage(result)
        raise Exception(message)

# FCS type
cdef extern from "fcs.h":
  cdef struct FCS_t:
    pass
  ctypedef FCS_t* FCS
  FCSResult fcs_init(FCS* handle, const char* method, MPI_Comm communicator)
  FCSResult fcs_destroy(FCS handle)
  FCSResult fcs_tune (FCS handle, fcs_int local_particles, fcs_int local_max_particles, fcs_float *positions,  fcs_float *charges)

cdef class scafacos:
    cdef FCS handle

    def __init__(self, method):
        handleResult(fcs_init(&self.handle, method, MPI_COMM_WORLD))

    def tune(self,
             np.ndarray[fcs_float, ndim=2, mode='c'] positions not None,
             np.ndarray[fcs_float, ndim=1, mode='c'] charges not None):
        N = positions.shape[1]
        if N != charges.shape[0]:
            raise Exception("The number of charges and positions must be equal!")
        handleResult(fcs_tune(self.handle, N, N, &positions[0,0], &charges[0]))
        
    def __del__(self):
        handleResult(fcs_destroy(self.handle))
        
        
        
