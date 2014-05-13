cimport numpy as np
import numpy as np
from libc.string cimport memcpy

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

    # constructor and destructor
    FCSResult fcs_init(FCS* handle, const char* method, MPI_Comm communicator)
    FCSResult fcs_destroy(FCS handle)

    # different properties
    const char* fcs_get_method_name(FCS handle)
    FCSResult fcs_set_near_field_flag(FCS handle, fcs_int near_field_flag)
    fcs_int fcs_get_near_field_flag(FCS handle)
    FCSResult fcs_set_total_particles(FCS handle, fcs_int total_particles)
    fcs_int fcs_get_total_particles(FCS handle)

    # box vectors
    FCSResult fcs_set_box_a(FCS handle, const fcs_float *box_a)
    const fcs_float* fcs_get_box_a(FCS handle)
    FCSResult fcs_set_box_b(FCS handle, const fcs_float *box_b)
    const fcs_float* fcs_get_box_b(FCS handle)
    FCSResult fcs_set_box_c(FCS handle, const fcs_float *box_c)
    const fcs_float* fcs_get_box_c(FCS handle)
    FCSResult fcs_set_periodicity(FCS handle, const fcs_int *periodicity)
    const fcs_int *fcs_get_periodicity(FCS handle)
  
    # tune
    FCSResult fcs_tune (FCS handle, fcs_int local_particles,
                        fcs_float *positions,  fcs_float *charges)

    FCSResult fcs_run(FCS handle, fcs_int local_particles,
                      fcs_float *positions, fcs_float *charges,
                      fcs_float *fields, fcs_float *potentials)

  
cdef class scafacos:
    cdef FCS handle

    def __init__(self, method):
        handleResult(fcs_init(&self.handle, method, MPI_COMM_WORLD))
        self.near_field_flag = False

    property total_particles:
        def __get__(self):
            return fcs_get_total_particles(self.handle)

        def __set__(self, n):
            handleResult(fcs_set_total_particles(self.handle, n))
        
    property method:
        def __get__(self):
            return fcs_get_method_name(self.handle)

    property box:
        def __set__(self, np.ndarray[fcs_float, ndim=2, mode='c'] box):
            handleResult(fcs_set_box_a(self.handle, &box[0,0]))
            handleResult(fcs_set_box_b(self.handle, &box[1,0]))
            handleResult(fcs_set_box_c(self.handle, &box[2,0]))

        def __get__(self):
            cdef np.ndarray[fcs_float, ndim=2, mode='c'] box = np.empty((3,3))
            cdef const fcs_float* a = fcs_get_box_a(self.handle)
            cdef const fcs_float* b = fcs_get_box_b(self.handle)
            cdef const fcs_float* c = fcs_get_box_c(self.handle)
            memcpy(&box[0,0], a, 3*sizeof(fcs_float))
            memcpy(&box[1,0], b, 3*sizeof(fcs_float))
            memcpy(&box[2,0], c, 3*sizeof(fcs_float))
            return box

    property periodicity:
        def __set__(self, p):
            cdef fcs_int cp[3]
            if p[0]: cp[0] = 1
            else: cp[0] = 0
            if p[1]: cp[1] = 1
            else: cp[1] = 0
            if p[2]: cp[2] = 1
            else: cp[2] = 0
            handleResult(fcs_set_periodicity(self.handle, cp))

        def __get__(self):
            cdef const fcs_int* p = fcs_get_periodicity(self.handle)
            return bool(p[0]), bool(p[1]), bool(p[2])

    property near_field_flag:
        def __set__(self, f):
            cdef fcs_int cf = int(f)
            fcs_set_near_field_flag(self.handle, cf)

        def __get__(self):
            return bool(fcs_get_near_field_flag(self.handle))
        
    def tune(self,
             np.ndarray[fcs_float, ndim=2, mode='c'] positions not None,
             np.ndarray[fcs_float, ndim=1, mode='c'] charges not None):
        N = positions.shape[1]
        if N != charges.shape[0]:
            raise Exception("The number of charges and positions must be equal!")
        handleResult(fcs_tune(self.handle, N, &positions[0,0], &charges[0]))

    def __call__(self,
            np.ndarray[fcs_float, ndim=2, mode='c'] positions not None,
            np.ndarray[fcs_float, ndim=1, mode='c'] charges not None):
        N = positions.shape[1]
        if N != charges.shape[0]:
            raise Exception("The number of charges and positions must be equal!")

        cdef np.ndarray[fcs_float, ndim=2, mode='c'] fields = np.empty_like(positions)
        cdef np.ndarray[fcs_float, ndim=1, mode='c'] potentials = np.empty_like(charges)
        
        handleResult(fcs_run(self.handle, N, &positions[0,0], &charges[0],
                             &fields[0,0], &potentials[0]))
        return fields, potentials
        
    def __del__(self):
        handleResult(fcs_destroy(self.handle))
        
        
        
