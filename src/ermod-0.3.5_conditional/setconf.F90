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
!
!
! renaming outside parameters and parameters to avoid conflict
module OUTname
  use engmain, only: trjfile, inffile, io_MDinfo
  use trajectory, only: handle
!
  implicit none
  integer OUTens, OUTbxs, OUTcmb, OUTclt, OUTspo
  real OUTtemp, OUTelc, OUTlwl, OUTupl, OUTscr
  integer OUTew1, OUTew2, OUTew3, OUTms1, OUTms2, OUTms3
  integer OUTnrun, OUTntype
  integer, dimension(:), allocatable :: OUTnmol, OUTsite

  type(handle) :: history_trajectory
  type(handle) :: solute_trajectory

contains

  subroutine OUTinitial
    real, save :: sgmcnv, chgcnv, engcnv, lencnv ! unit conversion factor
    sgmcnv = 2.0 ** (5.0 / 6.0)                  ! from Rmin/2 to sigma
    chgcnv = 1.0 / 1.60217653e-19                ! from C to elementary
    engcnv = 6.0221415e23 / 4.184e3              ! from J to kcal/mol
    lencnv = 10.0                                ! from nm to Angstrom
  end subroutine OUTinitial

  subroutine opentrj
    use trajectory, only: open_trajectory
    call open_trajectory(history_trajectory, trjfile)
  end subroutine opentrj

  subroutine initconf
    call OUTinitial
  end subroutine initconf

  subroutine closetrj
    use trajectory, only: close_trajectory
    call close_trajectory(history_trajectory)
  end subroutine closetrj

  subroutine finiconf
  end subroutine finiconf
 
  ! Initialization - read MDinfo
  ! when ermod is built into MD program, read topologies from parent MD program
  subroutine OUT_MDinfo
    implicit none
    open(unit = io_MDinfo, file = inffile, status = 'old')
    read(io_MDinfo, *) OUTnrun, OUTntype
    allocate( OUTnmol(OUTntype), OUTsite(OUTntype) )
    read(io_MDinfo, *) OUTnmol(1:OUTntype)
    read(io_MDinfo, *) OUTsite(1:OUTntype)
    close(io_MDinfo)
    return
  end subroutine OUT_MDinfo

  ! read system setup (coulomb rule, LJ, etc..)
  ! used when the ermod program is built into MD program and runs on-the-fly
  subroutine OUTintprm
  end subroutine OUTintprm

  ! Get molecular configuration
  subroutine OUTconfig(OUTpos, OUTcell, OUTatm, OUTbox, particle_type, calc_type)
    use mpiproc, only: halt_with_error
    use trajectory, only: read_trajectory, open_trajectory, close_trajectory
    implicit none
    integer, intent(in) :: OUTatm, OUTbox
    character*6, intent(in) :: particle_type
    character*10, intent(in) :: calc_type
    real, intent(out) :: OUTpos(3, OUTatm), OUTcell(3, 3)
    integer :: status
    !
    select case(calc_type)
    case('trjfl_read')              ! reading from a trajectory file
       ! do nothing
    case('on_the_fly')              ! on-the-fly calculation during MD
       stop "Unavailable in the present version"
    case default
       stop "Unknown calc_type in OUTconfig"
    end select
    !
    OUTcell(:, :) = 0.0
    !
    select case(particle_type)
    case('system')                  ! reading the HISTORY file
       call read_trajectory(history_trajectory, OUTatm, (OUTbox == 1), OUTpos, OUTcell, status)
       if(status /= 0) call halt_with_error("set_trj")
    case('solute')                  ! reading the SltConf file
       call read_trajectory(solute_trajectory, OUTatm, .false., OUTpos, OUTcell, status)
       if(status /= 0) then
          ! wrap around
          call close_trajectory(solute_trajectory)
          call open_trajectory(solute_trajectory, "SltConf")
          call read_trajectory(solute_trajectory, OUTatm, .false., OUTpos, OUTcell, status)
          if(status /= 0) then
             stop "Failed to reload solute trajectory"
          endif
       endif
    case default
       stop "Unknown particle_type"
    end select
    !
    return
  end subroutine
!
end module
!
!
module setconf
  implicit none
contains
!
!  setting molecular simulation parameters
!
  subroutine iniparam
    use engmain, only: init_params, &
         iseed, &
         skpcnf, corrcal, selfcal, &
         slttype, wgtslf, wgtsys, wgtins, boxshp, estype, &
         sltspec, hostspec, refspec, ljformat, ljswitch, &
         insorigin, insposition, insorient, insstructure, &
         sltpick, refpick, inscnd, inscfg, &     ! deprecated   
         lwreg, upreg, lwstr, upstr, &
         inptemp, temp, &
         engdiv, maxins, &
         intprm, elecut, lwljcut, upljcut, &
         cmbrule, cltype, screen, ewtoler, splodr, plmode, &
         ew1max, ew2max, ew3max, ms1max, ms2max, ms3max, &
         block_threshold, &
         force_calculation, &
         NO, YES, &
         SYS_NONPERIODIC, SYS_PERIODIC, &
         ES_NVT, ES_NPT, &
         LJFMT_EPS_cal_SGM_nm, LJFMT_EPS_Rminh, LJFMT_EPS_J_SGM_A, &
         LJFMT_A_C, LJFMT_C12_C6, LJFMT_TABLE, &
         LJSWT_POT_CHM, LJSWT_POT_GMX, LJSWT_FRC_CHM, LJSWT_FRC_GMX, &
         LJCMB_ARITH, LJCMB_GEOM, &
         EL_COULOMB, EL_EWALD, EL_PME, EL_PPPM, &
         SLT_SOLN, SLT_REFS_RIGID, SLT_REFS_FLEX, &
         INSORG_ORIGIN, INSORG_NOCHANGE, INSORG_AGGCEN, INSORG_REFSTR, &
         INSPOS_RANDOM, INSPOS_NOCHANGE, &
         INSPOS_SPHERE, INSPOS_SLAB_GENERIC, INSPOS_SLAB_SYMMETRIC, &
         INSPOS_RMSD, INSPOS_GAUSS, &
         INSROT_RANDOM, INSROT_NOCHANGE, &
         INSSTR_NOREJECT, INSSTR_RMSD
    use OUTname, only: OUTintprm, &              ! from outside
         OUTens, OUTbxs, OUTtemp, &              ! from outside
         OUTelc, OUTlwl, OUTupl, &               ! from outside
         OUTcmb, OUTclt, OUTscr, OUTspo, &       ! from outside
         OUTew1, OUTew2, OUTew3, &               ! from outside
         OUTms1, OUTms2, OUTms3                  ! from outside
    use mpiproc, only: halt_with_error                               ! MPI
    implicit none
    real, parameter :: tiny = 1.0e-20
    real :: real_seed
    character(len=3) :: scrtype

    intprm=1                                     ! trajectory reading
    select case(intprm)
    case(0)  ! on-the-fly reading from parent MD setup (currently ineffective)
      call OUTintprm
      boxshp = OUTbxs  ; estype = OUTens  ; inptemp = OUTtemp
      elecut = OUTelc  ; lwljcut = OUTlwl ; upljcut = OUTupl
      cmbrule = OUTcmb ; cltype = OUTclt  ; screen = OUTscr  ; splodr = OUTspo
      ew1max = OUTew1  ; ew2max = OUTew2  ; ew3max = OUTew3
      ms1max = OUTms1  ; ms2max = OUTms2  ; ms3max = OUTms3
    case(1)  ! default settings for trajectory reading
      boxshp = SYS_PERIODIC                      ! periodic boundary
      estype = ES_NPT                            ! constant pressure
      cltype = EL_PME  ; ms1max = 64             ! PME
      inptemp = 300.0                            ! Kelvin
      engdiv = 1                                 ! number of divisions
      screen = 0.0     ; ewtoler = 0.0  
      ewtoler = 1.0e-6 ; elecut = 12.0
      splodr = 4       ; scrtype = 'dis'
      upljcut = elecut ; lwljcut = upljcut - 2.0
    case default
       stop "Unknown intprm"
    end select
    
    sltpick = 0 ; refpick = 0 ; inscnd = 0 ; inscfg = 0    ! deprecated

    insposition = INSPOS_RANDOM
    insorient = INSROT_RANDOM
    insstructure = INSSTR_NOREJECT

    ! block-wise calculation, corresponds to 13 atoms / box
    block_threshold = 4.0

    force_calculation = .false.

    ! only part of constants set here
    call init_params()

    ! default settings
    skpcnf = 1                     ! no skip for trajectory reading
    
    wgtslf = NO
    wgtsys = NO
    wgtins = NO

    selfcal = NO                   ! no construction of self-energy distribution

    select case(slttype)
    case(SLT_SOLN)
       corrcal = NO                ! no calculation of correlation matrix
    case(SLT_REFS_RIGID, SLT_REFS_FLEX)
       corrcal = YES               ! calculation of correlation matrix
       if((cltype == EL_EWALD) .or. (cltype == EL_PME) &
            .or. (cltype == EL_PPPM)) then  ! Ewald, PME and PPPM
          wgtslf = YES
       endif
    case default
       stop "Unknown slttype"
    end select

    sltspec = 1
    hostspec(:) = 0
    refspec(:) = 0

    ljformat = LJFMT_EPS_Rminh
    ljswitch = LJSWT_POT_CHM

    select case(slttype)
    case(SLT_SOLN)
       maxins = 1                  ! not used in calculation of solution
    case(SLT_REFS_RIGID, SLT_REFS_FLEX)
       maxins = 1000               ! default number of particle insertions
    case default
       stop "Unknown slttype"
    end select

    ! Unphysical initialization of lwreg and upreg
    ! The program will terminate if lwreg or upreg stays unphysical when used
    lwreg = -1.0
    upreg = lwreg - 1.0

    ! Unphysical initialization of lwstr and upstr
    ! The program will terminate if lwstr or upstr stays unphysical when used
    lwstr = -1.0
    upstr = lwstr - 1.0

    if(intprm /= 0) then
      cmbrule = LJCMB_ARITH        ! arithmetic mean of LJ sigma
      ew2max = ew1max ; ew3max = ew1max
      ms2max = ms1max ; ms3max = ms1max
    endif

    if(sltpick > 0) sltspec = sltpick                      ! deprecated
    if(refpick > 0) refspec(1) = refpick                   ! deprecated

    select case(inscnd)                                    ! deprecated
    case(0)    ! random
       insposition = INSPOS_RANDOM
    case(1)    ! spherical
       insposition = INSPOS_SPHERE
    case(2)    ! slab (symmetric bilayer)
       insposition = INSPOS_SLAB_SYMMETRIC
    case(3)    ! reference
       insposition = INSPOS_RMSD
    case default
       stop "Unknown inscnd"
    end select

    select case(inscfg)                                    ! deprecated
    case(0)    ! only the intramolecular configuration is from the file
       insorient = INSROT_RANDOM
    case(1)    ! orientation is fixed from the file with random position
       insorient = INSROT_NOCHANGE
    case(2)    ! position and orientation are both fixed from the file
       insorient = INSROT_NOCHANGE ; insposition = INSPOS_NOCHANGE
    case default
       stop "Unknown inscfg"
    end select

    ! default settings done

    ! read again for non-default constants
    call init_params()

    ! initialize random seed
    if(iseed == 0) then
       CALL RANDOM_SEED
       CALL RANDOM_NUMBER(real_seed)
       iseed = 100000 + int(899999 * real_seed)
    endif
    
    ! temperature converted into the unit of kcal/mol
    temp = inptemp * 8.314510e-3 / 4.184

    ! get the screening parameter in Ewald and PME
    if((cltype == EL_EWALD) .or. (cltype == EL_PME) &
         .or. (cltype == EL_PPPM)) then  ! Ewald, PME and PPPM
       if(screen <= tiny) then
          if(ewtoler <= tiny) call halt_with_error('set_ewa')
          screen = getscrn(ewtoler, elecut, scrtype)
       endif
    endif

    plmode = 2                  ! energy calculation parallelization mode

    ! check cltype and related parameters
    select case(cltype)
    case(EL_COULOMB)
    case(EL_EWALD)   ! Ewald parameters, not effective in the current version
       call halt_with_error('set_prs')
       if(ew1max * ew2max * ew3max == 0) call halt_with_error('set_ewa')
    case(EL_PME)     ! PME parameters
       if(ms1max * ms2max * ms3max == 0) call halt_with_error('set_ewa')
    case(EL_PPPM)     ! PPPM parameters
       if(ms1max * ms2max * ms3max == 0) call halt_with_error('set_ewa')
    case default
       stop "Unknown cltype"
    end select

    ! check ljformat parameter
    select case(ljformat)
    case(LJFMT_EPS_cal_SGM_nm, LJFMT_EPS_Rminh, LJFMT_EPS_J_SGM_A, &
         LJFMT_A_C, LJFMT_C12_C6, LJFMT_TABLE)
    case default
       stop "Unknown ljformat"
    end select

    ! check ljswitch parameter
    select case(ljswitch)
    case(LJSWT_POT_CHM, LJSWT_POT_GMX, LJSWT_FRC_CHM, LJSWT_FRC_GMX)
    case default
       stop "Unknown ljswitch"
    end select

    ! check maxins parameter
    select case(slttype)
    case(SLT_SOLN)
       maxins = 1                  ! not used in calculation of solution
    case(SLT_REFS_RIGID)
       if(maxins <  1) call halt_with_error('set_ins')
       ! maxins > 1 makes no sense if coordinate is used as is for rigid solute
       if((insposition == INSPOS_NOCHANGE) .and. &
          (insorient == INSROT_NOCHANGE) .and.   &
          (maxins /= 1)) call halt_with_error('set_ins')
    case(SLT_REFS_FLEX)
       if(maxins <  1) call halt_with_error('set_ins')
    case default
       stop "Unknown slttype"
    end select

    ! check the consistency in parameters for non-periodic system
    if(boxshp == SYS_NONPERIODIC) then
       if((estype == ES_NPT) .or. &
          (cltype /= EL_COULOMB)) call halt_with_error('set_prs')
    endif

    ! check the consistency of insertion with structure-dependent weight
    if(wgtins == YES) then
       if(slttype /= SLT_REFS_FLEX) call halt_with_error('set_ins')
    endif

    ! set insorigin parameter and check insposition parameter
    select case(insposition)
    case(INSPOS_RANDOM)                               ! random
       insorigin = INSORG_ORIGIN
    case(INSPOS_NOCHANGE)                             ! fixed configuration
       insorigin = INSORG_NOCHANGE
    case(INSPOS_SPHERE)                               ! sphere geometry
       insorigin = INSORG_AGGCEN
       ! check lwreg and upreg parameters
       if((lwreg < 0.0) .or. (upreg < 0.0) .or. &
          (lwreg > upreg)) call halt_with_error('set_reg')
    case(INSPOS_SLAB_GENERIC, INSPOS_SLAB_SYMMETRIC)  ! slab configuration
       insorigin = INSORG_AGGCEN
       ! check lwreg and upreg parameters
       if(lwreg > upreg) call halt_with_error('set_reg')
       if(insposition == INSPOS_SLAB_SYMMETRIC) then  ! symmetric bilayer
          if((lwreg < 0.0) .or. (upreg < 0.0)) call halt_with_error('set_reg')
       endif
    case(INSPOS_RMSD, INSPOS_GAUSS)                   ! comparison to reference
       insorigin = INSORG_REFSTR
       ! check lwreg and upreg parameters
       if((lwreg < 0.0) .or. (upreg < 0.0) .or. &
          (lwreg > upreg)) call halt_with_error('set_reg')
    case default
       stop "Unknown insposition"
    end select

    ! check insorigin parameter
    select case(insorigin)
    case(INSORG_ORIGIN)
       if(insposition /= INSPOS_RANDOM) call halt_with_error('set_bug')
    case(INSORG_NOCHANGE)
       if(insposition /= INSPOS_NOCHANGE) call halt_with_error('set_bug')
    case(INSORG_AGGCEN)
       if((insposition /= INSPOS_SPHERE) .and. &
          (insposition /= INSPOS_SLAB_GENERIC) .and. &
          (insposition /= INSPOS_SLAB_SYMMETRIC)) call halt_with_error('set_bug')
    case(INSORG_REFSTR)
       if((insposition /= INSPOS_RMSD) .and. &
          (insposition /= INSPOS_GAUSS)) call halt_with_error('set_bug')
    case default
       stop "Unknown insorigin"
    end select

    ! check insorient parameter
    select case(insorient)
    case(INSROT_RANDOM, INSROT_NOCHANGE)
       ! do nothing
    case default
       stop "Unknown insorient"
    end select

    ! check insstructure parameter
    select case(insstructure)
    case(INSSTR_NOREJECT)    ! no rejection of solute structure
       ! do nothing
    case(INSSTR_RMSD)        ! solute structure rejection with RMSD
       ! check lwstr and upstr parameters
       if((lwstr < 0.0) .or. (upstr < 0.0) .or. &
          (lwstr > upstr)) call halt_with_error('set_str')
       ! check the solute type
       if(slttype == SLT_REFS_RIGID) call halt_with_error('set_slt')
    case default
       stop "Unknown insstructure"
    end select

    return
  end subroutine iniparam

  subroutine check_param
  ! check the consistency between parameters read in the iniparam subroutine
  !                           and those read in the OUT_MDinfo subroutine
    use engmain, only: numtype, slttype, hostspec, refspec, insposition, &
         SLT_SOLN, SLT_REFS_RIGID, SLT_REFS_FLEX, &
         INSPOS_SPHERE, INSPOS_SLAB_GENERIC, INSPOS_SLAB_SYMMETRIC, &
         INSPOS_RMSD, INSPOS_GAUSS
    use mpiproc, only: halt_with_error
    implicit none
    integer :: phys_numtype

    ! number of physical species in the system
    select case(slttype)
    case(SLT_SOLN)
       phys_numtype = numtype
    case(SLT_REFS_RIGID, SLT_REFS_FLEX)
       phys_numtype = numtype - 1
    end select

    ! when restrained relative to aggregate
    ! hostspec is either 0 (undefined)
    !             or between 1 and (numtype for soln, numtype - 1 for refs)
    ! with the insposition specified here, at least one of hostspec is non-zero
    ! Number of non-zero entries of hostspec <= number of species in the system
    if((insposition == INSPOS_SPHERE) .or. &
       (insposition == INSPOS_SLAB_GENERIC) .or. &
       (insposition == INSPOS_SLAB_SYMMETRIC)) then
       if(count( mask = (hostspec(:) < 0) ) > 0) call halt_with_error('set_ins')
       if(count( mask = (hostspec(:) >= 1) ) == 0) call halt_with_error('set_ins')
       if(any(hostspec(:) > phys_numtype)) call halt_with_error('set_ins')
       if(count( mask = (hostspec(:) >= 1) ) > phys_numtype) call halt_with_error('set_ins')
    else
       hostspec(:) = 0
    endif

    ! when restrained against reference
    ! refspec is either 0 (undefined)
    !            or between 1 and (numtype for soln, numtype - 1 for refs)
    ! with the insposition specified here, at least one of refspec is non-zero
    ! Number of non-zero entries of refspec <= number of species in the system
    if((insposition == INSPOS_RMSD) .or. (insposition == INSPOS_GAUSS)) then
       if(count( mask = (refspec(:) < 0) ) > 0) call halt_with_error('set_ins')
       if(count( mask = (refspec(:) >= 1) ) == 0) call halt_with_error('set_ins')
       if(any(refspec(:) > phys_numtype)) call halt_with_error('set_ins')
       if(count( mask = (refspec(:) >= 1) ) > phys_numtype) call halt_with_error('set_ins')
    else
       refspec(:) = 0
    endif

    return
  end subroutine check_param


  real function getscrn(ewtoler, elecut, scrtype)
    implicit none
    character(len=3), intent(in) :: scrtype
    real, intent(in) :: ewtoler, elecut
    real :: ewasml, ewalrg, scrfac, factor
    real, parameter :: error = 1.0e-20
    factor = error + 1.0 ; ewasml = 0.0 ; ewalrg = 1.0e3
    do while(factor > error)
       scrfac = (ewasml + ewalrg) / 2.0
       factor = erfc(scrfac * elecut)
       if(scrtype == 'dis') factor = factor / elecut
       if(factor > ewtoler) then
          ewasml = scrfac
       else
          ewalrg = scrfac
       endif
       factor = abs(factor - ewtoler)
    end do
    getscrn = scrfac
    return
  end function getscrn

  ! Calls OUTinitial / iniparam / OUT_MDinfo, and sets parameters
  subroutine setparam
    use engmain, only: numtype, nummol, numatm, maxcnf, &
         slttype, sltspec, ljformat, &
         moltype, numsite, sluvid, &
         bfcoord, sitemass, charge, &
         ljene_mat, ljlensq_mat, ljtype, ljtype_max, cmbrule, &
         specatm, sitepos, &
         mol_begin_index, mol_end_index, belong_to, mol_charge, &
         solute_file, solvent_file, mol_io, ljtable_file, ljtable_io, &
         LJFMT_EPS_cal_SGM_nm, LJFMT_EPS_Rminh, LJFMT_EPS_J_SGM_A, &
         LJFMT_A_C, LJFMT_C12_C6, LJFMT_TABLE, &
         LJCMB_ARITH, LJCMB_GEOM, &
         SLT_SOLN, SLT_REFS_RIGID, SLT_REFS_FLEX, &
         PT_SOLVENT, PT_SOLUTE, PT_TEST_RIGID, PT_TEST_FLEX
    use OUTname, only: OUTinitial, OUT_MDinfo, &   ! from outside
         OUTnrun, OUTntype, OUTnmol, OUTsite       ! from outside
    use mpiproc, only: halt_with_error
    use utility, only: itoa
    implicit none
    ! only integer power is allowed as the initialization expression (7.1.6.1)
    real, parameter :: sgmcnv = 1.7817974362806784 ! from Rmin/2 to sigma, 2.0**(5.0/6.0)
    real, parameter :: lencnv = 10.0               ! from nm to Angstrom
    real, parameter :: engcnv = 1.0 / 4.184        ! from kJ/mol to kcal/mol
    integer :: pti, stmax, maxsite, uvtype, cmin, cmax, sid, i, ati, m
    integer :: solute_index, cur_solvent, prev_solvent_type, cur_atom
    real :: factor, xst(3)
    real, dimension(:), allocatable :: sitemass_temp, charge_temp
    integer, allocatable :: ljtype_temp(:)
    real, dimension(:), allocatable :: ljlen_temp, ljene_temp
    real, dimension(:), allocatable :: ljlen_temp_table, ljene_temp_table
    integer :: ljtype_found
    logical :: lj_is_new
    integer, dimension(:), allocatable :: pttype, ptcnt, ptsite
    real, dimension(:,:), allocatable :: psite
    character(len=5) :: atmtype
    character(len=80) :: molfile

    call OUTinitial                ! initialization of OUTname module
    call iniparam                  ! initialization of parameters
    call OUT_MDinfo                ! reading MDinfo

    maxcnf = OUTnrun                                 ! from outside
    select case(slttype)
    case(SLT_SOLN)
       numtype = OUTntype                            ! from outside
    case(SLT_REFS_RIGID, SLT_REFS_FLEX)
       numtype = OUTntype + 1                        ! from outside
    end select

    ! check the consistency between the parameters
    ! read in the iniparam and OUT_MDinfo subroutines
    call check_param

    ! pttype is particle type for each molecule group
    allocate( pttype(numtype), ptcnt(numtype), ptsite(numtype) )
    ! Default is physical particle, coordinate existing in (HISTORY) trajectory
    pttype(:) = PT_SOLVENT
    ! Default: same as written in MDinfo
    ptcnt(1:OUTntype) = OUTnmol(1:OUTntype)
    ptsite(1:OUTntype) = OUTsite(1:OUTntype)

    if(slttype == SLT_SOLN) then
       ! Determine which molecule is solute?
       solute_index = 1            ! default for soln
       if((1 <= sltspec) .and. (sltspec <= numtype)) solute_index = sltspec
       pttype(solute_index) = PT_SOLUTE
    else
       ! Test particle information will be defined outside
       ptcnt(numtype) = 1          ! Test particle can only be one molecule

       open(unit = mol_io, file = solute_file, status = 'old')
       ! here only counts no. of lines
       stmax = 0
       do
          read(mol_io, *, end = 99) m
          stmax = stmax + 1
       end do
99     close(mol_io)
       ptsite(numtype) = stmax
       pttype(numtype) = slttype   ! Test particle is the last (for insertion)
       solute_index = numtype
    endif

    ! count up number of mols, 
    ! set max and total no. of atoms 
    nummol = sum( ptcnt(1:numtype) )
    numatm = sum( ptcnt(1:numtype) * ptsite(1:numtype) )

    allocate( moltype(nummol), numsite(nummol), sluvid(nummol) )

    ! make mapping from molecule no. [1..nummol] to particle type [1..numtype]
    cmin = 1
    cmax = 0
    do pti = 1, numtype
       cmax = cmax + ptcnt(pti)
       moltype(cmin:cmax) = pti ! sequential identification
       cmin = cmax + 1
    end do
    if(cmax /= nummol) call halt_with_error("set_num")

    ! Assign the number of sites within molecule
    numsite(1:nummol) = ptsite( moltype(1:nummol) )

    ! Build the solute/solvent specification
    sluvid(1:nummol) = pttype( moltype(1:nummol) )

    ! check if all the solute molecules have the same number of atoms
    stmax = -1
    do i = 1, nummol
       if(moltype(i) == solute_index) then        ! solute
          if(stmax == -1) then                 ! initialize
             stmax = numsite(i)
          else
             if(stmax /= numsite(i)) call halt_with_error("set_slt")
          endif
       endif
    end do

    if(numatm /= sum( numsite(1:nummol) )) stop "something is wrong in setconf::setparam, numatm"


    allocate( bfcoord(3, stmax), sitemass(numatm) )
    allocate( charge(numatm), ljtype(numatm) )
    allocate( sitepos(3, numatm))
    allocate( mol_begin_index(nummol + 1) )
    allocate( belong_to(numatm) )
    allocate( mol_charge(nummol) )

    ! initial setting to zero
    bfcoord(:,:) = 0.0
    sitemass(:) = 0.0
    charge(:) = 0.0

    ! initialize mol_begin_index
    ! mol_begin_index(i) .. (mol_begin_index(i + 1) - 1)
    !    will be the index range for i-th molecule
    ! mol_end_index is defined in the engmain module as a function through
    !    mol_end_index(i) = mol_begin_index(i + 1) - 1
    mol_begin_index(1) = 1
    do i = 1, nummol
       mol_begin_index(i + 1) = mol_begin_index(i) + numsite(i)
    end do
    if(mol_begin_index(nummol + 1) /= numatm + 1) call halt_with_error("set_bug")
    if(mol_end_index(nummol) /= numatm) call halt_with_error("set_bug")

    ! initialize belong_to(map from atom number to molecule no)
    do i = 1, nummol
       belong_to(mol_begin_index(i):mol_end_index(i)) = i
    end do

    ! large enough LJ table size
    allocate( ljlen_temp_table(1:sum( ptsite(:) )), &
              ljene_temp_table(1:sum( ptsite(:) )) )

    ! temporary set of LJ & coordinates
    maxsite = maxval(ptsite(1:numtype))
    allocate( psite(3,maxsite), sitemass_temp(maxsite), charge_temp(maxsite), &
              ljtype_temp(maxsite), ljlen_temp(maxsite), ljene_temp(maxsite) )

    cur_solvent = 0
    cur_atom = 1
    ljtype_max = 0

    ! read molecules specification
    do pti = 1, numtype
       uvtype = pttype(pti)
       if(uvtype == PT_SOLVENT) then            ! solvent
          cur_solvent = cur_solvent + 1
          molfile = solvent_file//trim(adjustl(itoa(cur_solvent)))
       else
          if(ptcnt(pti) > 1) cur_solvent = cur_solvent + 1
          molfile = solute_file                 ! solute / test particle
       endif
       stmax = ptsite(pti)

       ! This part is a bit complicated due to backward compatibility.
       ! for ljtype /= 5, read the table and make table by program
       open(unit = mol_io, file = molfile, status = 'old')
       do sid = 1, stmax
          if(uvtype == SLT_REFS_RIGID) then
             read(mol_io,*) m, atmtype, xst(1:3), psite(1:3,sid)
          else
             read(mol_io,*) m, atmtype, xst(1:3)
          endif
          call getmass(sitemass_temp(sid), atmtype)

          charge_temp(sid) = xst(1)
          if(ljformat == LJFMT_EPS_Rminh) xst(3) = sgmcnv * xst(3)
          if((ljformat == LJFMT_A_C) .or. (ljformat == LJFMT_C12_C6)) then
             if(xst(3) /= 0.0) then
                factor = (xst(2) / xst(3)) ** (1.0 / 6.0)
                xst(2) = xst(3) / (4.0 * (factor ** 6))
                xst(3) = factor
             else
                xst(2) = 0.0
             endif
          endif
          if((ljformat == LJFMT_EPS_J_SGM_A) .or. (ljformat == LJFMT_C12_C6)) then
             xst(2) = engcnv * xst(2)
             xst(3) = lencnv * xst(3)
          endif
          ljene_temp(sid) = xst(2)
          ljlen_temp(sid) = xst(3)
       end do
       close(mol_io)

       if(ljformat == LJFMT_TABLE) then
          ! use numbers directly
          ! No sane system will have the problem with 
          ! string -> double -> int conversion ...
          ljtype_temp(1:stmax) = ljene_temp(1:stmax)
       else
          ! allocate LJ types
          do sid = 1, stmax
             lj_is_new = .true.
             do i = 1, ljtype_max
                ! linear search LJ table
                if((ljlen_temp_table(i) == ljlen_temp(sid)) .and. &
                   (ljene_temp_table(i) == ljene_temp(sid))) then
                   ljtype_found = i
                   lj_is_new = .false.
                   exit
                endif
             end do
             if(lj_is_new) then
                ! new LJ type
                ljtype_max = ljtype_max + 1
                ljlen_temp_table(ljtype_max) = ljlen_temp(sid)
                ljene_temp_table(ljtype_max) = ljene_temp(sid)
                ljtype_found = ljtype_max
             endif
             ljtype_temp(sid) = ljtype_found
          end do
       endif

       do i = 1, ptcnt(pti)
          ljtype(cur_atom:(cur_atom + stmax - 1)) = ljtype_temp(1:stmax)
          charge(cur_atom:(cur_atom + stmax - 1)) = charge_temp(1:stmax)
          sitemass(cur_atom:(cur_atom + stmax - 1)) = sitemass_temp(1:stmax)
          cur_atom = cur_atom + stmax
       end do

       if(uvtype == SLT_REFS_RIGID) bfcoord(1:3, 1:stmax) = psite(1:3, 1:stmax)
    end do

    deallocate( psite, sitemass_temp, charge_temp, &
                ljlen_temp, ljene_temp, ljtype_temp )

    ! Fill LJ table
    if(ljformat == LJFMT_TABLE) then
       ! From table (directly)
       open(unit = ljtable_io, file = ljtable_file, status = 'old', action = 'read')
       read(ljtable_io, *) ljtype_max
       allocate( ljlensq_mat(ljtype_max, ljtype_max), &
                 ljene_mat(ljtype_max, ljtype_max) )
       do i = 1, ljtype_max
          read (ljtable_io, *) ljlensq_mat(i, 1:ljtype_max)
          ljlensq_mat(i, 1:ljtype_max) = ljlensq_mat(i, 1:ljtype_max) ** 2
       end do
       do i = 1, ljtype_max
          read (ljtable_io, *) ljene_mat(i, 1:ljtype_max)
       end do
       close(ljtable_io)
    else
       ! From LJ data
       allocate( ljlensq_mat(ljtype_max, ljtype_max), &
                 ljene_mat(ljtype_max, ljtype_max) )
       do i = 1, ljtype_max
          select case(cmbrule)
          case(LJCMB_ARITH)    ! arithmetic mean
             ljlensq_mat(1:ljtype_max, i) = (( ljlen_temp_table(1:ljtype_max) &
                                             + ljlen_temp_table(i) ) / 2.0) ** 2
          case(LJCMB_GEOM)     ! geometric mean
             ljlensq_mat(1:ljtype_max, i) = ljlen_temp_table(1:ljtype_max) &
                                          * ljlen_temp_table(i)
          case default
             stop "Incorrect cmbrule"
          end select
          ljene_mat(1:ljtype_max, i) = sqrt( ljene_temp_table(1:ljtype_max)  &
                                           * ljene_temp_table(i) )
       end do
    endif
    deallocate( ljlen_temp_table, ljene_temp_table )


    ! conversion to (kcal/mol angstrom)^(1/2)
    ! == sqrt(e^2 * coulomb const * avogadro / (kcal / mol angstrom))
    charge(1:numatm) = 18.22261721 * charge(1:numatm)

    ! get molecule-wise charges
    do i = 1, nummol
       mol_charge(i) = sum( charge(mol_begin_index(i):mol_end_index(i)) )
    end do

    deallocate( pttype, ptcnt, ptsite )
  end subroutine setparam


! Reads up to maxread coordinates serially and distributes frames
! returns number of frames read (EXCLUDING skipped frames)
  subroutine getconf_parallel(maxread, actual_read)
    use engmain, only: skpcnf, boxshp, numsite, sluvid, stdout, &
    ! start of the extension for computing the conditional distributions
                       do_conditional, OrderPrm_io, OrderPrm_read, &
                       OrderPrm_Values, OrderPrm_ArraySize, &
                       YES, &
    ! end of the extension for computing the conditional distributions
                       sitepos, cell, stat_weight_system, &
                       perm_file, perm_io, &
                       PT_SOLVENT, PT_SOLUTE
    use OUTname, only: OUTconfig                     ! from outside
    use mpiproc
    implicit none
    integer, intent(in) :: maxread
    integer, intent(out) :: actual_read
    
    real, dimension(:,:), allocatable :: OUTpos, OUTcell, readpos
    real :: readcell(3, 3)
    real :: weight, readweight
    integer :: i, OUTatm, iproc, nread

    logical, save :: first_time = .true.
    integer, allocatable, save :: permutation(:)
    integer, allocatable :: count_perm(:)
    integer :: stat
    ! start of the extension for computing the conditional distributions
    integer :: cnt_Order, dummy_i, check_i
    integer, parameter :: tag_order = 17
    real, dimension(:), allocatable :: read_order(:)
    if(myrank == 0) allocate( read_order(OrderPrm_ArraySize) )
    ! end of the extension for computing the conditional distributions

    ! sum over solvent & solute in trajectory file (HISTORY); no test particle
    OUTatm = sum( numsite, &
                  mask = ((sluvid == PT_SOLVENT) .or. (sluvid == PT_SOLUTE)) )
    if(myrank == 0) allocate( readpos(3, OUTatm) )
    allocate( OUTpos(3, OUTatm), OUTcell(3, 3) )

    ! first time setup: read index permutation
    if(first_time) then
       first_time = .false.
       open(file = perm_file, unit = perm_io, status = 'old', action = 'read', iostat = stat)
       if(stat == 0) then ! file exists and is successfully opened
          if(myrank == 0) write(stdout, *) "Reading permutation information"
          allocate( permutation(OUTatm) )
          do i = 1, OUTatm
             permutation(i) = i
          end do
          do
             ! the i-th atom in the trajectory (HISTORY) file is set to
             ! the permutation(i)-th atom in the ermod program
             ! All the other variables such as mass and interaction parameters
             ! and the input files have the order of particles AFTER permutation
             ! N.B. Fortran 95 std. 9.4.4.4 allows the reuse of [i]
             read(perm_io, *, end = 99) i, permutation(i)
             if (i < 1 .or. i > OUTatm) then
                print *, "Incorrect range of particle number at ",i
                call halt_with_error('set_pmt')
             endif
             if (permutation(i) < 1 .or. permutation(i) > OUTatm) then
                print *, i, "-th permutation points to ", permutation(i)
                call halt_with_error('set_pmt')
             endif
          end do
99        close(perm_io)
          ! each particle appears once and only once in permutation
          allocate( count_perm(OUTatm) )
          count_perm(:) = 0
          do i = 1, OUTatm
             count_perm(permutation(i)) = count_perm(permutation(i)) + 1
          end do
          if( any(count_perm /= 1) ) then
             do i = 1, OUTatm
                if(count_perm(i) /= 1) then
                   print *, "Atom ", i, " has ", count_perm(i), " entries"
                end if
             end do
             call halt_with_error('set_pmt')
          end if
          deallocate( count_perm )
       endif
    end if
    
    nread = min(nprocs, maxread)

    if(myrank < nread) then
       if(myrank == 0) then              ! rank-0 to read from file
          do iproc = 1, nread
             ! get configuration and store in OUTpos / OUTcell
             do i = 1, skpcnf
                call OUTconfig(readpos, readcell, OUTatm, boxshp, &
                                                  'system','trjfl_read')
                call read_weight(readweight)
    ! start of the extension for computing the conditional distributions
      ! note the consistency with getsolute subroutine in insertion.F90
      ! caution: program will not work as expected when skpcnf > 1
                if((do_conditional == YES) .and. (OrderPrm_read == YES)) then
                   do cnt_Order = 1, OrderPrm_ArraySize
                      read(OrderPrm_io, *, iostat = stat) dummy_i, check_i, read_order(cnt_Order)
                      if((check_i /= cnt_Order) .or. (stat /= 0)) call halt_with_error('set_bug')
                   enddo
                endif
    ! end of the extension for computing the conditional distributions
             end do
             
             if(iproc /= 1) then         ! send the data to other rank
#ifdef MPI
                ! FIXME: rewrite with grouping and scatter?
                call mpi_send(readpos, 3 * OUTatm, mpi_double_precision, &
                              iproc - 1, tag_coord, mpi_comm_world, ierror)
                call mpi_send(readcell, 3 * 3, mpi_double_precision, &
                              iproc - 1, tag_cell, mpi_comm_world, ierror)
                call mpi_send(readweight, 1, mpi_double_precision, &
                              iproc - 1, tag_weight, mpi_comm_world, ierror)
    ! start of the extension for computing the conditional distributions
                call mpi_send(read_order, OrderPrm_ArraySize, mpi_double_precision, &
                              iproc - 1, tag_order, mpi_comm_world, ierror)
    ! end of the extension for computing the conditional distributions
#endif
             else                        ! rank-0 to use the data as read
                ! if this causes memory bottleneck, rewrite with in-place permutation
                OUTpos(:, :) = readpos(:, :)
                OUTcell(:, :) = readcell(:, :)
                weight = readweight
    ! start of the extension for computing the conditional distributions
                OrderPrm_Values(:) = read_order(:)
    ! end of the extension for computing the conditional distributions
             endif
          end do
       else                              ! non-0 rank to receive the data
#ifdef MPI
          call mpi_recv(OUTpos, 3 * OUTatm, mpi_double_precision, &
                        0, tag_coord, mpi_comm_world, mpistatus, ierror)
          call mpi_recv(OUTcell, 3 * 3, mpi_double_precision, &
                        0, tag_cell, mpi_comm_world, mpistatus, ierror)
          call mpi_recv(weight, 1, mpi_double_precision, &
                        0, tag_weight, mpi_comm_world, mpistatus, ierror)
    ! start of the extension for computing the conditional distributions
          call mpi_recv(OrderPrm_Values, OrderPrm_ArraySize, mpi_double_precision, &
                        0, tag_order, mpi_comm_world, mpistatus, ierror)
    ! end of the extension for computing the conditional distributions
#endif
       endif

       if(allocated(permutation)) then
          sitepos(:, permutation(1:OUTatm)) = OUTpos(:, 1:OUTatm)
       else
          sitepos(:, 1:OUTatm) = OUTpos(:, 1:OUTatm)
       endif
       cell(:, :) = OUTcell(:, :)
       stat_weight_system = weight
    endif

    ! start of the extension for computing the conditional distributions
    if(myrank == 0) deallocate( read_order )
    ! end of the extension for computing the conditional distributions
    if(myrank == 0) deallocate( readpos )
    deallocate( OUTpos, OUTcell )
    actual_read = nread

  end subroutine getconf_parallel

  subroutine read_weight(weight)
    use engmain, only: wgtsys, syswgt_file, syswgt_io, stdout, YES
    use mpiproc, only: mpi_setup
    implicit none
    real, intent(out) :: weight
    integer :: dummy, ioerr    
    logical, save :: file_opened = .false.

    if(wgtsys /= YES) then
       weight = 1.0
       return
    endif

    if(.not. file_opened) then
       open(unit = syswgt_io, file = syswgt_file, action = 'read')
       file_opened = .true.
    endif
    
    read(syswgt_io, *, iostat = ioerr) dummy, weight
    if(ioerr /= 0) then 
       write(stdout, *) " The weight file (", syswgt_file, ") is ill formed or its length does not match with that of HISTORY"
       call mpi_setup('stop')
       stop
    endif

  end subroutine read_weight


  subroutine getmass(stmass,atmtype)
    implicit none
    real, parameter :: massH = 1.00794          ! atomic weight (hydrogen)
    real, parameter :: massC = 12.0107          ! atomic weight (carbon)
    real, parameter :: massO = 15.9994          ! atomic weight (oxygen)
    real, parameter :: massN = 14.00674         ! atomic weight (nitrogen)
    real, parameter :: massS = 32.066           ! atomic weight (sulfur)
    real, parameter :: massP = 30.973761        ! atomic weight (phosphorus)
    real, parameter :: massHe = 4.0026          ! atomic weight (helium)
    real, parameter :: massB = 10.811           ! atomic weight (boron)
    real, parameter :: massSi = 28.0855         ! atomic weight (silicon)
    real, parameter :: massLi = 6.941           ! atomic weight (lithium)
    real, parameter :: massNa = 22.989770       ! atomic weight (sodium)
    real, parameter :: massK = 39.0983          ! atomic weight (potassium)
    real, parameter :: massF = 18.9984032       ! atomic weight (fluorine)
    real, parameter :: massCl = 35.4527         ! atomic weight (chlorine)
    real, parameter :: massBr = 79.904          ! atomic weight (bromine)
    real, parameter :: massI =  126.90447       ! atomic weight (iodine)
    real, parameter :: massCa = 40.078          ! atomic weight (calcium)
    real, parameter :: massZn = 65.409          ! atomic weight (zinc)
    real, parameter :: massFe = 55.845          ! atomic weight (iron)
    real, parameter :: massCu = 63.546          ! atomic weight (copper)

    real, intent(out) :: stmass
    character(len=5), intent(in) :: atmtype
    character(len=1) :: eltp1
    character(len=2) :: eltp2
    character(len=3) :: eltp3

    eltp1 = atmtype(1:1)
    if(eltp1 == 'H') stmass = massH
    if(eltp1 == 'C') stmass = massC
    if(eltp1 == 'O') stmass = massO
    if(eltp1 == 'N') stmass = massN
    if(eltp1 == 'S') stmass = massS
    if(eltp1 == 'P') stmass = massP
    if(eltp1 == 'B') stmass = massB
    if(eltp1 == 'K') stmass = massK
    if(eltp1 == 'F') stmass = massF
    if(eltp1 == 'I') stmass = massI
    eltp2 = atmtype(1:2)
    if(eltp2 == 'He') stmass = massHe
    if(eltp2 == 'Si') stmass = massSi
    if(eltp2 == 'Li') stmass = massLi
    if(eltp2 == 'Na') stmass = massNa
    if(eltp2 == 'Cl') stmass = massCl
    if(eltp2 == 'Br') stmass = massBr
    if(eltp2 == 'Ca') stmass = massCa
    if(eltp2 == 'Zn') stmass = massZn
    if(eltp2 == 'Fe') stmass = massFe
    if(eltp2 == 'Cu') stmass = massCu
    if(eltp2 == 'CH') stmass = massC + massH
    eltp3 = atmtype(1:3)
    if(eltp3 == 'CH2') stmass = massC + 2.0 * massH
    if(eltp3 == 'CH3') stmass = massC + 3.0 * massH
    return
  end subroutine getmass
end module
