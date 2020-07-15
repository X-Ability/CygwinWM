! -*- F90 -*-
! ERmod - Eneregy Representation Module
! Copyright (C) 2000-2016 Nobuyuki Matubayasi
! Copyright (C) 2010-2016 Shun Sakuraba
! 
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

! mpi module

module mpiproc
#ifdef MPI
  ! MPI
  use mpi
#endif
  implicit none

#ifndef MPI
  integer, parameter :: mpi_status_size = 1
#endif

  ! mpi common variables
  integer :: ierror, mpistatus(mpi_status_size)
  integer :: myrank, nprocs       ! rank number and total number of processes
  integer :: nactiveproc          ! number of active processes
  integer, parameter :: tag_cell = 11, tag_coord = 12, tag_weight = 13
  integer, parameter :: tag_sltcrd = 22, tag_sltwgt = 23
  integer :: mpi_comm_activeprocs

  ! for performance counter
  integer :: prev_state_num
  real(8) :: prev_time
  real(8) :: times(9)
  character(len=8), parameter :: timer_names(9) = (/&
       "RealBlck", &
       "RealSelf", &
       "RecpChrg", &
       "RecpFFT ", &
       "RecpU-V ", &
       "RecpSelf", &
       "RecpGren", &
       "RecpPrep", &
       "RealPrep" /)

contains
  subroutine mpi_setup(type)
    implicit none
    character(len=4) :: type
    integer :: i
#ifdef MPI
    real(4) :: cputime
    real(8), save :: walltime
    if(type == 'init') then
       call mpi_init(ierror)
#ifdef PERF
       walltime = MPI_WTIME()
#endif
       call mpi_rank_size_info
       times(:) = 0.0
    endif
    if(type == 'stop') then
#ifdef PERF
       call CPU_TIME(cputime)
       print *, "rank = ", myrank, ", CPUtime = ", cputime
       if(myrank == 0) then
          print *, "Wall-clock time: ", MPI_WTIME() - walltime
          print *, "Timing breakdown: "
          do i = 1, size(times, 1)
             print *, timer_names(i), times(i)
          end do
       endif
#endif
       call mpi_finalize(ierror)
    endif
#else
    if(type == 'init') call mpi_rank_size_info
#endif
  end subroutine mpi_setup

  subroutine mpi_rank_size_info
    nprocs=1
    myrank=0
#ifdef MPI
    call mpi_comm_size(mpi_comm_world, nprocs, ierror)
    call mpi_comm_rank(mpi_comm_world, myrank, ierror)
#endif
    return
  end subroutine mpi_rank_size_info

  subroutine mpi_abend()
    integer :: ierror
#ifdef MPI
    call mpi_abort(mpi_comm_world, 1, ierror)
#endif
  end subroutine mpi_abend

  subroutine mpi_init_active_group(nactive)
    implicit none
    integer, intent(in) :: nactive

#ifdef MPI
    if(myrank < nactive) then
       call mpi_comm_split(mpi_comm_world, 1, myrank, mpi_comm_activeprocs, ierror)
    else
       call mpi_comm_split(mpi_comm_world, mpi_undefined, 0, mpi_comm_activeprocs, ierror)
    endif

    if(myrank >= nactive .and. mpi_comm_activeprocs /= mpi_comm_null) then
       stop "failed @ mpi_init_active_group"
    endif
#endif
  end subroutine mpi_init_active_group

  subroutine mpi_finish_active_group()
    implicit none
#ifdef MPI
    if(mpi_comm_activeprocs /= mpi_comm_null) then
       call mpi_comm_free(mpi_comm_activeprocs, ierror)
    end if
#endif
  end subroutine mpi_finish_active_group

  ! helper library for reduce variables
  ! Note: mpi_reduce(mpi_in_place, ...) seems to be allowed only on MPI 2.2+
  subroutine mympi_reduce_real(data, data_size, operation, rootrank)
    implicit none
    integer, intent(in) :: data_size, operation, rootrank
    real, intent(inout) :: data(data_size)
    real, allocatable :: buf(:)
    integer :: mpitype

#ifdef MPI
    allocate(buf(data_size))
    select case(kind(data))
    case(4)
       mpitype = mpi_real
    case(8)
       mpitype = mpi_double_precision
    case default
       stop "invalid kind(real) value"
    end select

    call mpi_reduce(data, buf, data_size, mpitype, operation, rootrank, mpi_comm_world, ierror)
    data(:) = buf(:)
    deallocate(buf)
#endif
  end subroutine mympi_reduce_real


  subroutine perf_time(state)
    implicit none
    character(len=4), intent(in), optional :: state
    real(8) :: wall_time

#ifdef MPI
#ifdef PERF
    wall_time = MPI_WTIME()
    if(prev_state_num /= 0) then
       times(prev_state_num) = times(prev_state_num) + (wall_time - prev_time)
    endif
    prev_time = wall_time

    prev_state_num = 0
    if(present(state)) then
       select case(state)
       case("rblk")
          prev_state_num = 1
       case("rslf")
          prev_state_num = 2
       case("kchg")
          prev_state_num = 3
       case("kfft")
          prev_state_num = 4
       case("kuve")
          prev_state_num = 5
       case("kslf")
          prev_state_num = 6
       case("kgrn")
          prev_state_num = 7
       case("kpre")
          prev_state_num = 8
       case("rpre")
          prev_state_num = 9
       end select
    end if
#endif
#endif
  end subroutine perf_time
  

  ! Stop calculation with error message
  subroutine halt_with_error(errtype)
    use engmain, only: stdout
    implicit none
    character(len=7), intent(in) :: errtype

    if(errtype == 'eng_typ') write(stdout, "(A)") " The number of solute types is incorrectly set"
    if(errtype == 'eng_num') write(stdout, "(A)") " The number of solute molecules is incorrectly set"
    if(errtype == 'eng_ins') write(stdout, "(A)") " The solute numbering is incorrect for insertion"
    if(errtype == 'eng_par') write(stdout, "(A)") " The input parameter is incorrectly set"
    if(errtype == 'eng_siz') write(stdout, "(A)") " The number of energy-coordinate meshes is too large"
    if(errtype == 'eng_min') write(stdout, "(A)") " The minimum of the energy coordinate is too large; " // &
         "the ecdmin parameter needs to be smaller"
    if(errtype == 'eng_sft') write(stdout, "(A)") " The eccore parameter is too small; " // &
         "there should be no distribution at energy coordinate larger than eccore in the solution system"
    if(errtype == 'eng_pcr') write(stdout, "(A)") " The pecore parameter is incorrectly set"
    if(errtype == 'eng_ecd') write(stdout, "(A)") " The energy-coordinate system is inconsistent"
    if(errtype == 'eng_per') write(stdout, "(A)") " parameters_er file does not exist"
    if(errtype == 'eng_eci') write(stdout, "(A)") " EcdInfo file does not exist"
    if(errtype == 'eng_ecm') write(stdout, "(A)") " EcdMesh file does not exist"
    if(errtype == 'eng_emf') write(stdout, "(A)") " EcdMesh file has wrong format"
    if(errtype == 'eng_cns') write(stdout, "(A)") " Inconsistency is present in the engproc program"
    if(errtype == 'eng_slb') write(stdout, "(A)") " Slab condition can only used in periodic system"
    if(errtype == 'eng_bug') write(stdout, "(A)") " Bug in engproc.F90"

    if(errtype == 'rcp_fst') write(stdout, "(A)") " The first particle needs to be the solute"
    if(errtype == 'rcp_cns') write(stdout, "(A)") " Inconsistency is present in the recpcal program"

    if(errtype == 'ins_set') write(stdout, "(A)") " The solute specification is incorrectly set"
    if(errtype == 'ins_siz') write(stdout, "(A)") " Inconsistency is present in the setting of the size of bfcoord"
    if(errtype == 'ins_geo') write(stdout, "(A)") " The system geometry is incorrectly set"
    if(errtype == 'ins_ref') write(stdout, "(A)") " RefInfo file is missing"
    if(errtype == 'ins_str') write(stdout, "(A)") " Incorrect size for RefInfo file"
    if(errtype == 'ins_bug') write(stdout, "(A)") " Bug in insertion.F90"

    if(errtype == 'set_slt') write(stdout, "(A)") " The solute type is incorrectly set"
    if(errtype == 'set_num') write(stdout, "(A)") " The number of molecules or atoms is incorrectly set"
    if(errtype == 'set_prs') write(stdout, "(A)") " The system parameters are incorrectly set"
    if(errtype == 'set_ins') write(stdout, "(A)") " The insertion parameters are incorrectly set"
    if(errtype == 'set_reg') write(stdout, "(A)") " The lwreg and/or upreg parameter is incorrectly set"
    if(errtype == 'set_str') write(stdout, "(A)") " The lwstr and/or upstr parameter is incorrectly set"
    if(errtype == 'set_ewa') write(stdout, "(A)") " The Ewald parameters are incorrectly set"
    if(errtype == 'set_trj') write(stdout, "(A)") " Trajectory is shorter than specified in MDinfo"
    if(errtype == 'set_pmt') write(stdout, "(A)") " Permutation index file is invalid"
    if(errtype == 'set_bug') write(stdout, "(A)") " Bug in setconf.F90"

    if(errtype == 'bst_zrw') write(stdout, "(A)") " Division by zero due to inappropriate setting of mass or weight"

    call mpi_abend()                                                     ! MPI
    stop
  end subroutine halt_with_error

  subroutine warning(typ)
    use engmain, only: stdout, force_calculation
    implicit none
    character(len=4), intent(in) :: typ
    if(typ == 'mbin') write(stdout, '(A)') " Warning: the maximum binning energy is too low for this species"
    if(typ == 'emax') write(stdout, '(A)') " Warning: number of total bins in distribution function is too large",&
         " and requiring too much memory"
    if(force_calculation) return
    write(stdout, '(A)') "The program aborts because there is a warning"
    write(stdout, '(A,A)') "If you wish to force program running, specify 'force_calculation = .true.' in parameters_er, ", &
         "at &ene_param section."

    call mpi_abend()
    stop
  end subroutine warning
end module mpiproc                                                       ! MPI
