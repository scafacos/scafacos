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

module module_vtk
      use module_base64
      implicit none

      integer, public, parameter :: VTK_VERTEX               =  1
      integer, public, parameter :: VTK_POLY_VERTEX          =  2
      integer, public, parameter :: VTK_LINE                 =  3
      integer, public, parameter :: VTK_POLY_LINE            =  4
      integer, public, parameter :: VTK_TRIANGLE             =  5
      integer, public, parameter :: VTK_TRIANGLE_STRIP       =  6
      integer, public, parameter :: VTK_POLYGON              =  7
      integer, public, parameter :: VTK_PIXEL                =  8
      integer, public, parameter :: VTK_QUAD                 =  9
      integer, public, parameter :: VTK_TETRA                = 10
      integer, public, parameter :: VTK_VOXEL                = 11
      integer, public, parameter :: VTK_HEXAHEDRON           = 12
      integer, public, parameter :: VTK_WEDGE                = 13
      integer, public, parameter :: VTK_PYRAMID              = 14
      integer, public, parameter :: VTK_QUADRATIC_EDGE       = 21
      integer, public, parameter :: VTK_QUADRATIC_TRIANGLE   = 22
      integer, public, parameter :: VTK_QUADRATIC_QUAD       = 23
      integer, public, parameter :: VTK_QUADRATIC_TETRA      = 24
      integer, public, parameter :: VTK_QUADRATIC_HEXAHEDRON = 25

      integer, public, parameter :: VTK_STEP_FIRST  = -1
      integer, public, parameter :: VTK_STEP_NORMAL =  0
      integer, public, parameter :: VTK_STEP_LAST   =  1

      character(*), parameter :: subfolder_collections = "./"
      character(*), parameter :: subfolder_vtk = "./vtk/"
      character(*), parameter :: visitfilename = "timeseries.visit"
      character(*), parameter :: paraviewfilename = "timeseries.pvd"
#ifndef LITTLEENDIAN
      logical, parameter :: bigendian = .true.
#else
      logical, parameter :: bigendian = .false.
#endif

      type vtkfile
        private
          character(40) :: filename
          integer :: filehandle = 96
          integer :: filehandle_par = 97
          integer :: filehandle_visit = 98
          integer :: filehandle_paraview = 99
          character(12) :: byte_order = "BigEndian"
          character(3) :: version = "0.1"
          integer :: my_rank
          integer :: num_pe
          real*8 :: simtime
          integer :: vtk_step
          character(3) ::filesuffix = 'vtk'
          integer :: communicator

          logical, public :: binary = .true. !< flag for switch binary/ascii output in xml-vtk files. should only be modified before calling vtkfile_create (or if the user knows waht he is doing)

        contains
          procedure :: create => vtkfile_create ! filename
          procedure :: create_parallel => vtkfile_create_parallel ! filename, mpi_comm, rank, num_pe --> rank 0 writes .pvtX-file
          procedure :: close => vtkfile_close
          procedure :: write_data_array_header => vtkfile_write_data_array_header
          procedure :: set_communicator => vtkfile_set_communicator

          procedure :: write_data_array_Real4_1  => vtkfile_write_data_array_Real4_1
          procedure :: write_data_array_Real4_3  => vtkfile_write_data_array_Real4_3
          procedure :: write_data_array_Real8_1  => vtkfile_write_data_array_Real8_1
          procedure :: write_data_array_Real8_2  => vtkfile_write_data_array_Real8_2
          procedure :: write_data_array_Real8_3  => vtkfile_write_data_array_Real8_3

          procedure :: write_data_array_Real8_1_field3  => vtkfile_write_data_array_Real8_1_field3
          procedure :: write_data_array_Real8_3_field3  => vtkfile_write_data_array_Real8_3_field3

          procedure :: write_data_array_Int4_1  => vtkfile_write_data_array_Int4_1
          procedure :: write_data_array_Int4_3  => vtkfile_write_data_array_Int4_3
          procedure :: write_data_array_Int8_1  => vtkfile_write_data_array_Int8_1
          procedure :: write_data_array_Int8_3  => vtkfile_write_data_array_Int8_3

          procedure :: write_data_Int4_1        => vtkfile_write_data_Int4_1
          procedure :: write_data_Int8_1        => vtkfile_write_data_Int8_1

          procedure :: write_data_repeat_Int4_1 => vtkfile_write_data_repeat_Int4_1
          procedure :: write_data_repeat_Int8_1 => vtkfile_write_data_repeat_Int8_1

          generic :: write_data_array => write_data_array_Real4_1, & ! name, one-dim real*4, number of entries
                                            write_data_array_Real4_3,  & ! name, three-dim real*4 as three separate arrays, number of entries
                                            write_data_array_Real8_1,  & ! ...
                                            write_data_array_Real8_2,  &
                                            write_data_array_Real8_3,  &
                                            write_data_array_Real8_1_field3,  &
                                            write_data_array_Real8_3_field3,  &
                                            write_data_array_Int4_1,   &
                                            write_data_array_Int4_3,   &
                                            write_data_array_Int8_1,   &
                                            write_data_array_Int8_3,   &
                                            write_data_Int4_1, &
                                            write_data_Int8_1, &
                                            write_data_repeat_Int4_1, &
                                            write_data_repeat_Int8_1
          procedure :: startpointdata => vtkfile_startpointdata
          procedure :: finishpointdata => vtkfile_finishpointdata
          procedure :: startcelldata => vtkfile_startcelldata
          procedure :: finishcelldata => vtkfile_finishcelldata
      end type vtkfile


      type, extends(vtkfile) :: vtkfile_unstructured_grid
        contains
          procedure :: create => vtkfile_unstructured_grid_create ! filename
          procedure :: create_parallel => vtkfile_unstructured_grid_create_parallel ! filename, mpi_comm, rank, num_pe --> rank 0 writes .pvtX-file
          procedure :: write_headers4 => vtkfile_unstructured_grid_write_headers
          procedure :: write_headers8 => vtkfile_unstructured_grid_write_headers8 
          generic :: write_headers => write_headers4, write_headers8 ! number of particles, writes anything incl. <Piece>
                                        
          procedure :: startpoints => vtkfile_unstructured_grid_startpoints
          procedure :: finishpoints => vtkfile_unstructured_grid_finishpoints
          procedure :: startcells => vtkfile_unstructured_grid_startcells
          procedure :: finishcells => vtkfile_unstructured_grid_finishcells
          procedure :: dont_write_cells => vtkfile_unstructured_grid_dont_write_cells
          procedure :: write_final => vtkfile_unstructured_grid_write_final
      end type vtkfile_unstructured_grid


      type, extends(vtkfile) :: vtkfile_rectilinear_grid
          integer, private, dimension(2,3) :: globaldims, mydims
        contains
          procedure :: create => vtkfile_rectilinear_grid_create ! filename
          procedure :: create_parallel => vtkfile_rectilinear_grid_create_parallel ! filename, mpi_comm, rank, num_pe --> rank 0 writes .pvtX-file
          procedure :: write_headers => vtkfile_rectilinear_grid_write_headers ! number of particles, writes anything incl. <Piece>
          procedure :: startcoordinates => vtkfile_rectilinear_grid_startcoordinates
          procedure :: finishcoordinates => vtkfile_rectilinear_grid_finishcoordinates
          procedure :: write_final => vtkfile_rectilinear_grid_write_final
      end type vtkfile_rectilinear_grid

      character(*), private, parameter :: subfolder = subfolder_collections//subfolder_vtk

      contains


      subroutine vtkfile_create(vtk, filename_, step_, simtime_, vtk_step_)
        implicit none
        class(vtkfile) :: vtk
        character(*) :: filename_
        integer :: step_
        real*8 :: simtime_
        integer :: vtk_step_
        call vtk%create_parallel(filename_, step_, 0, 0, simtime_, vtk_step_)
      end subroutine vtkfile_create


      subroutine vtkfile_create_parallel(vtk, filename_, step_, my_rank_, num_pe_, simtime_, vtk_step_, comm_)
        use module_utils
        implicit none
        include 'mpif.h'
        class(vtkfile) :: vtk
        character(*) :: filename_
        character(50) :: fn
        character(6) :: tmp
        real*8 :: simtime_
        integer :: my_rank_, num_pe_, step_
        integer :: vtk_step_
        integer, optional :: comm_

        vtk%num_pe   = max(num_pe_, 1)
        vtk%my_rank  = my_rank_
        vtk%simtime  = simtime_
        vtk%vtk_step = vtk_step_

        if (present(comm_)) then
            vtk%communicator = comm_
        else
            vtk%communicator = MPI_COMM_WORLD
        end if

        write(vtk%filename, '(a, "_", I6.6)') filename_, step_

        write(tmp,'(I6.6)') vtk%my_rank
        fn = trim(subfolder)//trim(vtk%filename)//"."//tmp//"."//trim(vtk%filesuffix)

        call create_directory(trim(subfolder))
        open(vtk%filehandle, file=fn)

        if (vtk%my_rank == 0) then
          open(vtk%filehandle_par, file=trim(subfolder)//trim(vtk%filename)//".p"//trim(vtk%filesuffix))

          if (vtk%vtk_step .eq. VTK_STEP_FIRST) then
            open(vtk%filehandle_visit, file=trim(subfolder_collections)//trim(filename_)//"."//visitfilename,STATUS='UNKNOWN', POSITION = 'REWIND')
            write(vtk%filehandle_visit, '("!NBLOCKS ", I0)') vtk%num_pe

            open(vtk%filehandle_paraview, file=trim(subfolder_collections)//trim(filename_)//"."//paraviewfilename,STATUS='UNKNOWN', POSITION = 'REWIND')
            write(vtk%filehandle_paraview, '("<VTKFile type=""Collection"">", /, "<Collection>")')
          else
            open(vtk%filehandle_visit, file=trim(subfolder_collections)//trim(filename_)//"."//visitfilename,STATUS='UNKNOWN', POSITION = 'APPEND')
            open(vtk%filehandle_paraview, file=trim(subfolder_collections)//trim(filename_)//"."//paraviewfilename,STATUS='UNKNOWN', POSITION = 'APPEND')
          endif
        endif
      end subroutine vtkfile_create_parallel


      subroutine vtkfile_close(vtk)
        implicit none
        class(vtkfile) :: vtk
        close(vtk%filehandle)

        if (vtk%my_rank == 0) then
          if (vtk%vtk_step .eq. VTK_STEP_LAST) then
            write(vtk%filehandle_paraview, '("</Collection>", / , "</VTKFile>")')
          endif

          close(vtk%filehandle_par)
          close(vtk%filehandle_visit)
          close(vtk%filehandle_paraview)
        endif
     end subroutine vtkfile_close


     subroutine vtkfile_set_communicator(vtk, comm)
       implicit none
        class(vtkfile) :: vtk
       integer, intent(in) :: comm

       vtk%communicator = comm
     end subroutine vtkfile_set_communicator


     subroutine vtkfile_write_data_array_header(vtk, name, number_of_components, type)
        implicit none
        class(vtkfile) :: vtk
        character(*) :: name, type
        integer :: number_of_components
        character(6) :: format

        if (vtk%binary) then
          format = "binary"
        else
          format = "ascii"
        end if

        write(vtk%filehandle, '("<DataArray Name=""",a,""" NumberOfComponents=""", I0, """ type=""", a ,""" format=""", a ,""">")') &
                 name, number_of_components, type, trim(format)
        if (vtk%my_rank == 0)  write(vtk%filehandle_par, '("<DataArray Name=""",a,""" NumberOfComponents=""", I0, """ type=""", a ,""" format=""", a ,"""/>")') &
                 name, number_of_components, type, trim(format)
     end subroutine


     subroutine vtkfile_write_data_array_Real4_1(vtk, name, data)
        implicit none
        class(vtkfile) :: vtk
        character(*) :: name
        integer :: ndata, i
        real*4 :: data(:)
        integer*4 :: numbytes
        type(base64_encoder) :: base64
        call vtk%write_data_array_header(name, 1, "Float32")
        
        ndata = size(data, kind=kind(ndata))

        if (vtk%binary) then
          numbytes = ndata*1*4
          call base64%start(vtk%filehandle, bigendian)
          call base64%encode(numbytes)
          call base64%finish()
          call base64%start(vtk%filehandle, bigendian)
          do i=1,ndata
            call base64%encode(data(i))
          end do
          call base64%finish()
        else
          do i=1,ndata
          write(vtk%filehandle, '(G14.6)') data(i)
          end do
        endif
        write(vtk%filehandle, '("</DataArray>")')
     end subroutine vtkfile_write_data_array_Real4_1


     subroutine vtkfile_write_data_array_Real4_3(vtk, name, data1, data2, data3)
        implicit none
        class(vtkfile) :: vtk
        character(*) :: name
        integer :: ndata, i
        real*4 :: data1(:), data2(:), data3(:)
        integer*4 :: numbytes
        type(base64_encoder) :: base64
        call vtk%write_data_array_header(name, 3, "Float32")

        ndata = size(data1, kind=kind(ndata))

        if (vtk%binary) then
          numbytes = ndata*3*4
          call base64%start(vtk%filehandle, bigendian)
          call base64%encode(numbytes)
          call base64%finish()
          call base64%start(vtk%filehandle, bigendian)
          do i=1,ndata
            call base64%encode(data1(i))
            call base64%encode(data2(i))
            call base64%encode(data3(i))
          end do
          call base64%finish()
        else
          do i=1,ndata
          write(vtk%filehandle, '(3G14.6)') data1(i), data2(i), data3(i)
          end do
        endif
        write(vtk%filehandle, '("</DataArray>")')
     end subroutine vtkfile_write_data_array_Real4_3


     subroutine vtkfile_write_data_array_Real8_1(vtk, name, data, scale)
        implicit none
        class(vtkfile) :: vtk
        character(*) :: name
        integer :: ndata, i
        real*8 :: data(:)
        real*8, optional, intent(in) :: scale
        integer*4 :: numbytes
        type(base64_encoder) :: base64
        call vtk%write_data_array_header(name, 1, "Float64")

        ndata = size(data, kind=kind(ndata))

        if (vtk%binary) then
          numbytes = ndata*1*8
          call base64%start(vtk%filehandle, bigendian)
          call base64%encode(numbytes)
          call base64%finish()
          call base64%start(vtk%filehandle, bigendian)
          if (present(scale)) then
            do i=1,ndata
              call base64%encode(data(i)*scale)
            end do
          else
            do i=1,ndata
              call base64%encode(data(i))
            end do
          end if
          call base64%finish()
        else
          if (present(scale)) then
            do i=1,ndata
              write(vtk%filehandle, '(G14.6)') data(i)*scale
            end do
          else
            do i=1,ndata
              write(vtk%filehandle, '(G14.6)') data(i)
            end do
          endif
        endif
        write(vtk%filehandle, '("</DataArray>")')
     end subroutine vtkfile_write_data_array_Real8_1


     subroutine vtkfile_write_data_array_Real8_2(vtk, name, data1, data2)
        implicit none
        class(vtkfile) :: vtk
        character(*) :: name
        integer :: ndata, i
        real*8 :: data1(:), data2(:)
        integer*4 :: numbytes
        type(base64_encoder) :: base64
        call vtk%write_data_array_header(name, 2, "Float64")

        ndata = size(data1, kind=kind(ndata))

        if (vtk%binary) then
          numbytes = ndata*3*8
          call base64%start(vtk%filehandle, bigendian)
          call base64%encode(numbytes)
          call base64%finish()
          call base64%start(vtk%filehandle, bigendian)
          do i=1,ndata
            call base64%encode(data1(i))
            call base64%encode(data2(i))
          end do
          call base64%finish()
        else
          do i=1,ndata
          write(vtk%filehandle, '(3G14.6)') data1(i), data2(i)
          end do
        endif
        write(vtk%filehandle, '("</DataArray>")')
     end subroutine vtkfile_write_data_array_Real8_2


     subroutine vtkfile_write_data_array_Real8_3(vtk, name, data1, data2, data3, scale)
        implicit none
        class(vtkfile) :: vtk
        character(*) :: name
        integer :: ndata, i
        real*8 :: data1(:), data2(:), data3(:)
        real*8, optional, intent(in) :: scale
        integer*4 :: numbytes
        type(base64_encoder) :: base64
        call vtk%write_data_array_header(name, 3, "Float64")

        ndata = size(data1, kind=kind(ndata))

        if (vtk%binary) then
          numbytes = ndata*3*8
          call base64%start(vtk%filehandle, bigendian)
          call base64%encode(numbytes)
          call base64%finish()
          call base64%start(vtk%filehandle, bigendian)
          if (present(scale)) then
            do i=1,ndata
              call base64%encode(scale*data1(i))
              call base64%encode(scale*data2(i))
              call base64%encode(scale*data3(i))
            end do
          else
            do i=1,ndata
              call base64%encode(data1(i))
              call base64%encode(data2(i))
              call base64%encode(data3(i))
            end do
          endif
          call base64%finish()
        else
          if (present(scale)) then
            do i=1,ndata
              write(vtk%filehandle, '(3G14.6)') scale*data1(i), scale*data2(i), scale*data3(i)
            end do
          else
            do i=1,ndata
              write(vtk%filehandle, '(3G14.6)') data1(i), data2(i), data3(i)
            end do
          endif
        endif
        write(vtk%filehandle, '("</DataArray>")')
     end subroutine vtkfile_write_data_array_Real8_3


     subroutine vtkfile_write_data_array_Real8_1_field3(vtk, name, data)
        implicit none
        class(vtkfile) :: vtk
        character(*) :: name
        integer :: ndata(1:3), i,j,k
        real*8 :: data(:,:,:)
        integer*4 :: numbytes
        type(base64_encoder) :: base64
        call vtk%write_data_array_header(name, 1, "Float64")
        
        ndata = shape(data, kind=kind(ndata))

        if (vtk%binary) then
          numbytes = product(ndata)*1*8
          call base64%start(vtk%filehandle, bigendian)
          call base64%encode(numbytes)
          call base64%finish()
          call base64%start(vtk%filehandle, bigendian)
          do k=1,ndata(3)
            do j=1,ndata(2)
              do i=1,ndata(1)
                call base64%encode(data(i,j,k))
              end do
            end do
          end do
          call base64%finish()
        else
          do k=1,ndata(3)
            do j=1,ndata(2)
              do i=1,ndata(1)
                write(vtk%filehandle, '(G14.6)') data(i,j,k)
              end do
            end do
          end do
        endif
        write(vtk%filehandle, '("</DataArray>")')
      end subroutine vtkfile_write_data_array_Real8_1_field3


     subroutine vtkfile_write_data_array_Real8_3_field3(vtk, name, data1, data2, data3)
        implicit none
        class(vtkfile) :: vtk
        character(*) :: name
        integer :: ndata(1:3), i, j, k
        real*8 :: data1(:,:,:), data2(:,:,:), data3(:,:,:)
        integer*4 :: numbytes
        type(base64_encoder) :: base64
        call vtk%write_data_array_header(name, 3, "Float64")

        ndata = shape(data1, kind=kind(ndata))

        if (vtk%binary) then
          numbytes = product(ndata)*3*8
          call base64%start(vtk%filehandle, bigendian)
          call base64%encode(numbytes)
          call base64%finish()
          call base64%start(vtk%filehandle, bigendian)
          do k=1,ndata(3)
            do j=1,ndata(2)
              do i=1,ndata(1)
                call base64%encode(data1(i,j,k))
                call base64%encode(data2(i,j,k))
                call base64%encode(data3(i,j,k))
              end do
            end do
          end do
          call base64%finish()
        else
          do k=1,ndata(3)
            do j=1,ndata(2)
              do i=1,ndata(1)
                write(vtk%filehandle, '(3G14.6)') data1(i,j,k), data2(i,j,k), data3(i,j,k)
              end do
            end do
          end do
        endif
        write(vtk%filehandle, '("</DataArray>")')
     end subroutine vtkfile_write_data_array_Real8_3_field3


     subroutine vtkfile_write_data_array_Int4_1(vtk, name, data)
        implicit none
        class(vtkfile) :: vtk
        character(*) :: name
        integer :: ndata, i
        integer*4 :: data(:)
        integer*4 :: numbytes
        type(base64_encoder) :: base64
        call vtk%write_data_array_header(name, 1, "Int32")

        ndata = size(data, kind=kind(ndata))

        if (vtk%binary) then
          numbytes = ndata*1*4
          call base64%start(vtk%filehandle, bigendian)
          call base64%encode(numbytes)
          call base64%finish()
          call base64%start(vtk%filehandle, bigendian)
          do i=1,ndata
            call base64%encode(data(i))
          end do
          call base64%finish()
        else
          do i=1,ndata
          write(vtk%filehandle, '(I20)') data(i)
          end do
        endif
        write(vtk%filehandle, '("</DataArray>")')
     end subroutine vtkfile_write_data_array_Int4_1


     subroutine vtkfile_write_data_Int4_1(vtk, name, data)
        implicit none
        class(vtkfile) :: vtk
        character(*) :: name
        integer*4 :: data

        call vtk%write_data_array_Int4_1(name, [ data ])
     end subroutine vtkfile_write_data_Int4_1


     subroutine vtkfile_write_data_Int8_1(vtk, name, data)
        implicit none
        class(vtkfile) :: vtk
        character(*) :: name
        integer*8 :: data

        call vtk%write_data_array_Int8_1(name, [ data ])
     end subroutine vtkfile_write_data_Int8_1


     subroutine vtkfile_write_data_array_Int4_3(vtk, name, data1, data2, data3)
        implicit none
        class(vtkfile) :: vtk
        character(*) :: name
        integer :: ndata, i
        integer*4 :: data1(:), data2(:), data3(:)
        integer*4 :: numbytes
        type(base64_encoder) :: base64
        call vtk%write_data_array_header(name, 3, "Int32")

        ndata = size(data1, kind=kind(ndata))

        if (vtk%binary) then
          numbytes = ndata*3*4
          call base64%start(vtk%filehandle, bigendian)
          call base64%encode(numbytes)
          call base64%finish()
          call base64%start(vtk%filehandle, bigendian)
          do i=1,ndata
            call base64%encode(data1(i))
            call base64%encode(data2(i))
            call base64%encode(data3(i))
          end do
          call base64%finish()
        else
          do i=1,ndata
          write(vtk%filehandle, '(3I20)') data1(i), data2(i), data3(i)
          end do
        endif
        write(vtk%filehandle, '("</DataArray>")')
     end subroutine vtkfile_write_data_array_Int4_3


     subroutine vtkfile_write_data_array_Int8_1(vtk, name, data)
        implicit none
        class(vtkfile) :: vtk
        character(*) :: name
        integer :: ndata, i
        integer*8 :: data(:)
        integer*4 :: numbytes
        type(base64_encoder) :: base64
        call vtk%write_data_array_header(name, 1, "Int64")

        ndata = size(data, kind=kind(ndata))

        if (vtk%binary) then
          numbytes = ndata*1*8
          call base64%start(vtk%filehandle, bigendian)
          call base64%encode(numbytes)
          call base64%finish()
          call base64%start(vtk%filehandle, bigendian)
          do i=1,ndata
            call base64%encode(data(i))
          end do
          call base64%finish()
        else
          do i=1,ndata
          write(vtk%filehandle, '(I20)') data(i)
          end do
        endif
        write(vtk%filehandle, '("</DataArray>")')
     end subroutine vtkfile_write_data_array_Int8_1


     subroutine vtkfile_write_data_array_Int8_3(vtk, name, data1, data2, data3)
        implicit none
        class(vtkfile) :: vtk
        character(*) :: name
        integer :: ndata, i
        integer*8 :: data1(:), data2(:), data3(:)
        integer*4 :: numbytes
        type(base64_encoder) :: base64
        call vtk%write_data_array_header(name, 3, "Int64")

        ndata = size(data1, kind=kind(ndata))

        if (vtk%binary) then
          numbytes = ndata*3*8
          call base64%start(vtk%filehandle, bigendian)
          call base64%encode(numbytes)
          call base64%finish()
          call base64%start(vtk%filehandle, bigendian)
          do i=1,ndata
            call base64%encode(data1(i))
            call base64%encode(data2(i))
            call base64%encode(data3(i))
          end do
          call base64%finish()
        else
          do i=1,ndata
          write(vtk%filehandle, '(3I20)') data1(i), data2(i), data3(i)
          end do
        endif
        write(vtk%filehandle, '("</DataArray>")')
     end subroutine vtkfile_write_data_array_Int8_3


     subroutine vtkfile_write_data_repeat_Int4_1(vtk, name, n, v)
        implicit none
        class(vtkfile) :: vtk
        character(*) :: name
        integer*4 :: n
        integer*4 :: v
        integer*4 :: i, numbytes
        type(base64_encoder) :: base64
        call vtk%write_data_array_header(name, 1, "Int32")
 
        if (vtk%binary) then
          numbytes = n*4
          call base64%start(vtk%filehandle, bigendian)
          call base64%encode(numbytes)
          call base64%finish()
          call base64%start(vtk%filehandle, bigendian)
          do i=1,n
            call base64%encode(v)
          end do
          call base64%finish()
        else
          do i=1,n
          write(vtk%filehandle, '(I20)') v
          end do
        endif
        write(vtk%filehandle, '("</DataArray>")')
     end subroutine vtkfile_write_data_repeat_Int4_1


     subroutine vtkfile_write_data_repeat_Int8_1(vtk, name, n, v)
        implicit none
        class(vtkfile) :: vtk
        character(*) :: name
        integer*4 :: n
        integer*8 :: v
        integer*4 :: i, numbytes
        type(base64_encoder) :: base64
        call vtk%write_data_array_header(name, 1, "Int64")
 
        if (vtk%binary) then
          numbytes = n*8
          call base64%start(vtk%filehandle, bigendian)
          call base64%encode(numbytes)
          call base64%finish()
          call base64%start(vtk%filehandle, bigendian)
          do i=1,n
            call base64%encode(v)
          end do
          call base64%finish()
        else
          do i=1,n
          write(vtk%filehandle, '(I20)') v
          end do
        endif
        write(vtk%filehandle, '("</DataArray>")')
     end subroutine vtkfile_write_data_repeat_Int8_1


     subroutine vtkfile_startpointdata(vtk)
        implicit none
        class(vtkfile) :: vtk

        write(vtk%filehandle, '("<PointData>")')
        if (vtk%my_rank == 0) write(vtk%filehandle_par, '("<PPointData>")')
     end subroutine vtkfile_startpointdata


     subroutine vtkfile_finishpointdata(vtk)
        implicit none
        class(vtkfile) :: vtk

        write(vtk%filehandle, '("</PointData>")')
        if (vtk%my_rank == 0) write(vtk%filehandle_par, '("</PPointData>")')
     end subroutine vtkfile_finishpointdata


     subroutine vtkfile_startcelldata(vtk)
        implicit none
        class(vtkfile) :: vtk

        write(vtk%filehandle, '("<CellData>")')
        if (vtk%my_rank == 0) write(vtk%filehandle_par, '("<PCellData>")')
     end subroutine vtkfile_startcelldata


     subroutine vtkfile_finishcelldata(vtk)
        implicit none
        class(vtkfile) :: vtk

        write(vtk%filehandle, '("</CellData>")')
        if (vtk%my_rank == 0) write(vtk%filehandle_par, '("</PCellData>")')
     end subroutine vtkfile_finishcelldata


    ! ########################### Unstructured Grid ################################################

     subroutine vtkfile_unstructured_grid_create(vtk, filename_, step_, simtime_, vtk_step_)
        implicit none
        class(vtkfile_unstructured_grid) :: vtk
        character(*) :: filename_
        integer :: step_
        real*8 :: simtime_
        integer :: vtk_step_
        call vtk%create_parallel(filename_, step_, 0, 0, simtime_, vtk_step_)
     end subroutine vtkfile_unstructured_grid_create


     subroutine vtkfile_unstructured_grid_create_parallel(vtk, filename_, step_, my_rank_, num_pe_, simtime_, vtk_step_, comm_)
        implicit none
        class(vtkfile_unstructured_grid) :: vtk
        character(*) :: filename_
        real*8 :: simtime_
        integer :: my_rank_, num_pe_, step_
        integer :: vtk_step_
        integer, optional :: comm_

        vtk%filesuffix = 'vtu'
        call vtkfile_create_parallel(vtk, filename_, step_, my_rank_, num_pe_, simtime_, vtk_step_, comm_)

     end subroutine vtkfile_unstructured_grid_create_parallel


     subroutine vtkfile_unstructured_grid_write_headers(vtk, npart, ncell)
        implicit none
        class(vtkfile_unstructured_grid) :: vtk
        integer :: npart, ncell

        write(vtk%filehandle, '("<VTKFile type=""UnstructuredGrid"" version=""", a, """ byte_order=""", a, """>")') vtk%version, trim(vtk%byte_order)
        write(vtk%filehandle, '("<UnstructuredGrid GhostLevel=""0"">")')
        write(vtk%filehandle, '("<Piece NumberOfPoints=""", I0, """ NumberOfCells=""", I0, """>")') npart, ncell

        if (vtk%my_rank == 0) then
          write(vtk%filehandle_par, '("<VTKFile type=""PUnstructuredGrid"" version=""", a, """ byte_order=""", a, """>")') vtk%version, trim(vtk%byte_order)
          write(vtk%filehandle_par, '("<PUnstructuredGrid GhostLevel=""0"">")')
        endif
     end subroutine vtkfile_unstructured_grid_write_headers


     subroutine vtkfile_unstructured_grid_write_headers8(vtk, npart, ncell)
        implicit none
        class(vtkfile_unstructured_grid) :: vtk
        integer*8 :: npart, ncell

        write(vtk%filehandle, '("<VTKFile type=""UnstructuredGrid"" version=""", a, """ byte_order=""", a, """>")') vtk%version, trim(vtk%byte_order)
        write(vtk%filehandle, '("<UnstructuredGrid GhostLevel=""0"">")')
        write(vtk%filehandle, '("<Piece NumberOfPoints=""", I0, """ NumberOfCells=""", I0, """>")') npart, ncell

        if (vtk%my_rank == 0) then
          write(vtk%filehandle_par, '("<VTKFile type=""PUnstructuredGrid"" version=""", a, """ byte_order=""", a, """>")') vtk%version, trim(vtk%byte_order)
          write(vtk%filehandle_par, '("<PUnstructuredGrid GhostLevel=""0"">")')
        endif
     end subroutine vtkfile_unstructured_grid_write_headers8


     subroutine vtkfile_unstructured_grid_dont_write_cells(vtk)
        implicit none
        class(vtkfile_unstructured_grid) :: vtk

          call vtk%startcells()
          write(vtk%filehandle, '("<DataArray type=""Int32"" Name=""connectivity"" />")')
          write(vtk%filehandle, '("<DataArray type=""Int32"" Name=""offsets"" />")')
          write(vtk%filehandle, '("<DataArray type=""UInt8"" Name=""types"" />")')
          call vtk%finishcells()
          call vtk%startcelldata()
          call vtk%finishcelldata()

     end subroutine vtkfile_unstructured_grid_dont_write_cells


     subroutine vtkfile_unstructured_grid_write_final(vtk)
        implicit none
        class(vtkfile_unstructured_grid) :: vtk
        integer :: i
        character(6) :: tmp
        character(50) :: fn

          write(vtk%filehandle, '("</Piece>")')
          write(vtk%filehandle, '("</UnstructuredGrid>")')
          write(vtk%filehandle, '("</VTKFile>")')

          if (vtk%my_rank == 0) then
            write(vtk%filehandle_visit, '(/)')

            do i = 0,vtk%num_pe-1
              write(tmp,'(I6.6)') i
              fn = trim(vtk%filename)//"."//tmp//"."//vtk%filesuffix
              write(vtk%filehandle_par, '("<Piece Source=""", a, """/>")') trim(fn)
              write(vtk%filehandle_visit, '(a)') trim(subfolder)//trim(fn)
            end do

            write(vtk%filehandle_par, '("</PUnstructuredGrid>")')
            write(vtk%filehandle_par, '("</VTKFile>")')

            write(vtk%filehandle_paraview, '("<DataSet timestep=""", f0.5,""" file=""", a, """/>")') vtk%simtime, trim(subfolder)//trim(trim(vtk%filename)//".p"//vtk%filesuffix)
          endif
     end subroutine vtkfile_unstructured_grid_write_final


     subroutine vtkfile_unstructured_grid_startpoints(vtk)
        implicit none
        class(vtkfile_unstructured_grid) :: vtk

        write(vtk%filehandle, '("<Points>")')
        if (vtk%my_rank == 0) write(vtk%filehandle_par, '("<PPoints>")')
     end subroutine vtkfile_unstructured_grid_startpoints


     subroutine vtkfile_unstructured_grid_finishpoints(vtk)
        implicit none
        class(vtkfile_unstructured_grid) :: vtk

        write(vtk%filehandle, '("</Points>")')
        if (vtk%my_rank == 0) write(vtk%filehandle_par, '("</PPoints>")')
     end subroutine vtkfile_unstructured_grid_finishpoints


     subroutine vtkfile_unstructured_grid_startcells(vtk)
        implicit none
        class(vtkfile_unstructured_grid) :: vtk

        write(vtk%filehandle, '("<Cells>")')
        if (vtk%my_rank == 0) write(vtk%filehandle_par, '("<PCells>")')
     end subroutine vtkfile_unstructured_grid_startcells


     subroutine vtkfile_unstructured_grid_finishcells(vtk)
        implicit none
        class(vtkfile_unstructured_grid) :: vtk

        write(vtk%filehandle, '("</Cells>")')
        if (vtk%my_rank == 0) write(vtk%filehandle_par, '("</PCells>")')
     end subroutine vtkfile_unstructured_grid_finishcells


    ! ########################### Rectilinear Grid ################################################

     subroutine vtkfile_rectilinear_grid_create(vtk, filename_, step_, simtime_, vtk_step_)
        implicit none
        class(vtkfile_rectilinear_grid) :: vtk
        character(*) :: filename_
        integer :: step_
        real*8 :: simtime_
        integer :: vtk_step_
        call vtk%create_parallel(filename_, step_, 0, 0, simtime_, vtk_step_)
     end subroutine vtkfile_rectilinear_grid_create


     subroutine vtkfile_rectilinear_grid_create_parallel(vtk, filename_, step_, my_rank_, num_pe_, simtime_, vtk_step_, comm_)
        implicit none
        class(vtkfile_rectilinear_grid) :: vtk
        character(*) :: filename_
        real*8 :: simtime_
        integer :: my_rank_, num_pe_, step_
        integer :: vtk_step_
        integer, optional :: comm_

        vtk%filesuffix = 'vtr'
        call vtkfile_create_parallel(vtk, filename_, step_, my_rank_, num_pe_, simtime_, vtk_step_, comm_)

     end subroutine vtkfile_rectilinear_grid_create_parallel


     subroutine vtkfile_rectilinear_grid_write_headers(vtk, globaldims_, mydims_)
        implicit none
        class(vtkfile_rectilinear_grid) :: vtk
        integer, dimension(2, 3), intent(in) :: globaldims_, mydims_

        vtk%globaldims = globaldims_
        vtk%mydims     = mydims_

        write(vtk%filehandle, '("<VTKFile type=""RectilinearGrid"" version=""", a, """ byte_order=""", a, """>")') vtk%version, trim(vtk%byte_order)
        write(vtk%filehandle, '("<RectilinearGrid WholeExtent=""", 6(" ",I0), """>")') vtk%mydims
        write(vtk%filehandle, '("<Piece Extent=""", 6(" ",I0), """>")') vtk%mydims

        if (vtk%my_rank == 0) then
          write(vtk%filehandle_par, '("<VTKFile type=""PRectilinearGrid"" version=""", a, """ byte_order=""", a, """>")') vtk%version, trim(vtk%byte_order)
          write(vtk%filehandle_par, '("<PRectilinearGrid WholeExtent=""", 6(" ",I0), """ GhostLevel=""1"">")') vtk%globaldims
        endif
     end subroutine vtkfile_rectilinear_grid_write_headers


     subroutine vtkfile_rectilinear_grid_write_final(vtk)
        implicit none
        include 'mpif.h'
        class(vtkfile_rectilinear_grid) :: vtk
        integer :: i
        character(6) :: tmp
        character(50) :: fn

        integer, allocatable :: localdims(:,:, :)
        integer :: ierr

          write(vtk%filehandle, '("</Piece>")')
          write(vtk%filehandle, '("</RectilinearGrid>")')
          write(vtk%filehandle, '("</VTKFile>")')


          allocate(localdims(2, 3, vtk%num_pe))

          call MPI_GATHER(vtk%mydims, 2*3, MPI_INTEGER, localdims, 2*3, MPI_INTEGER, 0, vtk%communicator, ierr)

          if (vtk%my_rank == 0) then
            write(vtk%filehandle_visit, '(/)')

            do i = 0,vtk%num_pe-1
              write(tmp,'(I6.6)') i
              fn = trim(vtk%filename)//"."//tmp//"."//vtk%filesuffix
              write(vtk%filehandle_par, '("<Piece Extent=""", 6(" ",I0), """ Source=""", a, """/>")') localdims(:,:,i+1), trim(fn)
              write(vtk%filehandle_visit, '(a)') trim(subfolder)//trim(fn)
            end do

            write(vtk%filehandle_par, '("</PRectilinearGrid>")')
            write(vtk%filehandle_par, '("</VTKFile>")')

            write(vtk%filehandle_paraview, '("<DataSet timestep=""", f0.5,""" file=""", a, """/>")') vtk%simtime, trim(subfolder)//trim(vtk%filename)//".p"//trim(vtk%filesuffix)
          endif

          deallocate(localdims)

     end subroutine vtkfile_rectilinear_grid_write_final


     subroutine vtkfile_rectilinear_grid_startcoordinates(vtk)
        implicit none
        class(vtkfile_rectilinear_grid) :: vtk

        write(vtk%filehandle, '("<Coordinates>")')
        if (vtk%my_rank == 0) write(vtk%filehandle_par, '("<PCoordinates>")')
     end subroutine vtkfile_rectilinear_grid_startcoordinates


     subroutine vtkfile_rectilinear_grid_finishcoordinates(vtk)
        implicit none
        class(vtkfile_rectilinear_grid) :: vtk

        write(vtk%filehandle, '("</Coordinates>")')
        if (vtk%my_rank == 0) write(vtk%filehandle_par, '("</PCoordinates>")')
     end subroutine vtkfile_rectilinear_grid_finishcoordinates
end module module_vtk
