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

module ptinsrt
  ! test particle insertion of the solute
  implicit none
  real, save :: unrn
  !
  ! insertion against reference structure
  !   insorigin = INSORG_REFSTR: solvent species as superposition reference
  !      refspec: solvent species used as the reference
  !   insstructure = INSSTR_RMSD: solute itself as superposition reference
  !   refstr_file : filename for storing the reference structure
  !   refstr_io : IO for refstr_file
  ! variables for reference structure
  integer                            :: refhost_natom,   refslt_natom
  integer, dimension(:), allocatable :: refhost_specatm, refslt_specatm
  real, dimension(:,:), allocatable  :: refhost_crd,     refslt_crd
  real, dimension(:),   allocatable  :: refhost_weight,  refslt_weight
  real, dimension(:,:), allocatable  ::                  refslt_bestfit
  !
contains
  subroutine instslt(caltype, cntdst, stat_weight_solute)
    use engmain, only: nummol, slttype, numslt, sltlist, iseed, SLT_REFS_FLEX
    use mpiproc, only: halt_with_error
    implicit none
    character(len=4),  intent(in) :: caltype
    integer, optional, intent(in) :: cntdst
    real,    optional, intent(out) :: stat_weight_solute
    integer, save :: insml
    logical :: reject
    
    if((.not. present(cntdst)) .and. (.not. present(stat_weight_solute))) then
       select case(caltype)
       case('init')
          ! sanity check of solute specification
          reject = .false.
          if(numslt /= 1) reject = .true.
          if((numslt == 1) .and. (sltlist(1) /= nummol)) reject = .true.
          if(reject) call halt_with_error('ins_set')
          ! inserted solute is set to the last molecule in the system
          insml = sltlist(1)
          ! initialize the random number used for solute insertion
          call urand_init(iseed)
          ! opening the file for coordinate of flexible solute
          if(slttype == SLT_REFS_FLEX) call getsolute(caltype, insml)
          return
       case('last')
          ! closing the file for coordinate of flexible solute
          if(slttype == SLT_REFS_FLEX) call getsolute(caltype, insml)
          return
       case('proc')
          call halt_with_error('ins_bug')
       case default
          stop "Incorrect caltype in instslt"
       end select
    endif

    if((.not. present(cntdst)) .or. (.not. present(stat_weight_solute))) then
       call halt_with_error('ins_bug')
    endif

    if(slttype == SLT_REFS_FLEX) then
       ! get the configuration and structure-specific weight of flexible solute
       call getsolute(caltype, insml, cntdst, stat_weight_solute)
    else
       ! no structure-specific weight when the solute is rigid
       stat_weight_solute = 1.0
    endif

    reject = .true.
    do while(reject)
       call set_solute_origin(insml)
       call set_shift_com(insml, stat_weight_solute)
       call apply_orientation(insml)
       reject = .false.
       call check_solute_configuration(insml, reject)
       ! user-defined scheme to apply change / reject the solute configuration
       call insscheme(insml, reject)
    end do
    
    return
  end subroutine instslt

  subroutine set_solute_origin(insml)
    use engmain, only: insorigin, numsite, mol_begin_index, mol_end_index, &
                   bfcoord, sitepos, &
                   INSORG_ORIGIN, INSORG_NOCHANGE, INSORG_AGGCEN, INSORG_REFSTR
    implicit none
    integer, intent(in) :: insml
    integer :: molb, mole, nsite
    real :: syscen(3)

    nsite = numsite(insml)
    molb = mol_begin_index(insml)
    mole = mol_end_index(insml) 
    sitepos(1:3, molb:mole) = bfcoord(1:3, 1:nsite)

    select case(insorigin)
    case(INSORG_ORIGIN)
       syscen(:) = (/ 0., 0., 0. /)
       call set_solute_com(insml, syscen)     ! set solute COM to (0,0,0)
    case(INSORG_NOCHANGE)
       ! do nothing and use the coordinate as read from the file
    case(INSORG_AGGCEN)
       call com_aggregate(syscen)             ! get aggregate center (syscen)
       call set_solute_com(insml, syscen)     ! set solute COM to syscen
    case(INSORG_REFSTR)
       call reffit
    case default
       stop "Unknown insorigin in set_solute_origin"
    end select
  end subroutine set_solute_origin

  subroutine set_shift_com(insml, weight)
    use engmain, only: insposition, lwreg, upreg, &
                       numsite, mol_begin_index, mol_end_index, &
                       sitemass, sitepos, &
                       cell, celllen, boxshp, SYS_NONPERIODIC, &
                       INSPOS_RANDOM, INSPOS_NOCHANGE, &
                       INSPOS_SPHERE, &
                       INSPOS_SLAB_GENERIC, INSPOS_SLAB_SYMMETRIC, &
                       INSPOS_RMSD, INSPOS_GAUSS
    use mpiproc, only: halt_with_error
    use bestfit, only: center_of_mass
    implicit none
    integer, intent(in) :: insml
    real, intent(inout) :: weight
    
    integer :: i, insb, inse
    real :: com(3), syscen(3), r(3), norm, dir, t, maxdis, dst, movmax
    real, parameter :: margin_factor = 3.0   ! see below for explanation

    select case(insposition)
    case(INSPOS_RANDOM)
       ! fully random position within periodic box
       if(boxshp == SYS_NONPERIODIC) then    ! system has to be periodic
          call halt_with_error('ins_geo')
       endif
       do i = 1, 3
          call urand(r(i))
       end do
       com(:) = matmul(cell(:, :), r(:))
       call set_solute_com(insml, com)
       return
    case(INSPOS_NOCHANGE)
       ! fixed position as read from the file
       ! do nothing
       return
    case(INSPOS_SPHERE)
       ! spherical random position
       call urand(t)
       t = (t * (upreg ** 3 - lwreg ** 3) + lwreg ** 3) ** (1.0/3.0)
       do
          do i = 1, 3
             call urand(r(i))
          end do
          r(:) = r(:) * 2.0 - 1.0
          norm = sum(r ** 2)
          if((0 < norm) .and. (norm <= 1)) exit
       end do
       norm = sqrt(norm)
       com(:) = t * r(:) / norm
       call shift_solute_com(insml, com)
       return
    case(INSPOS_SLAB_GENERIC, INSPOS_SLAB_SYMMETRIC)
       ! slab random position
       if(boxshp == SYS_NONPERIODIC) then    ! system has to be periodic
          call halt_with_error('ins_geo')
       endif
       call urand(r(1))
       call urand(r(2))
       
       call urand(t)
       t = t * (upreg - lwreg) + lwreg

       if(insposition == INSPOS_SLAB_SYMMETRIC) then  ! symmetric bilayer
          call urand(dir)
          if(dir > 0.5) t = -t
       endif

       r(3) = t / celllen(3)
       com(:) = matmul(cell(:, :), r(:))
       call shift_solute_com(insml, com)
       return
    case(INSPOS_RMSD)
       ! random translation with upreg and molecular size
       insb = mol_begin_index(insml)
       inse = mol_end_index(insml)
       call center_of_mass(numsite(insml), sitepos(1:3, insb:inse), &
                                           sitemass(insb:inse), com)
       maxdis = 0
       do i = insb, inse
          r(1:3) = sitepos(1:3, i) - com(1:3)
          dst = sum( r(1:3) ** 2 )
          if(dst > maxdis) maxdis = dst
       end do
       maxdis = sqrt(maxdis)    ! maximum distance in solute from COM to atom
       ! margin_factor assures that the random translation within movmax
       ! and the following, random rotation encloses
       ! the total domain of solute RMSD < upreg
       movmax = upreg + margin_factor * maxdis
       do i = 1, 3
          call urand(r(i))
       end do
       com(:) = movmax * (2.0 * r(:) - 1.0)
       call shift_solute_com(insml, com)
       return
    case(INSPOS_GAUSS)
       ! 50% mixture of uniform distribution and weighted distribution
       call uniform_gauss_mixture(com, weight)
       call shift_solute_com(insml, com)
       return
    case default
       stop "Unknown insposition in set_shift_com"
    end select

  contains
    subroutine shift_solute_com(insml, com)
      use engmain, only: mol_begin_index, mol_end_index, sitepos
      implicit none
      integer, intent(in) :: insml
      real, intent(in) :: com(3)
      integer :: insb, inse, m
      insb = mol_begin_index(insml)
      inse = mol_end_index(insml)
      do m = 1, 3
         sitepos(m, insb:inse) = sitepos(m, insb:inse) + com(m)
      end do
    end subroutine shift_solute_com

    ! FIXME: this routine does not work for skewed periodic box
    subroutine uniform_gauss_mixture(com, weight)
      use mpiproc, only: myrank
      use engmain, only: PI, celllen, upreg
      implicit none
      real, intent(out) :: com(3)
      real, intent(inout) :: weight

      real, parameter :: uniform_ratio = 0.5

      integer :: i
      real :: scaled_coord(3), sqsum
      real :: r

      real, save :: l_of_sigma(3), z0, z1
      logical, save :: use_uniform = .false.
      logical, save :: first_time = .true.

      if(first_time) then
         l_of_sigma(:) = (celllen(:) / 2) / upreg
         if(myrank == 0) print *, "Lx/2 / sigma = ", l_of_sigma(:)
         z0 = 8 * l_of_sigma(1) * l_of_sigma(2) * l_of_sigma(3)
         z1 = (sqrt(PI) ** 3) * erf(l_of_sigma(1)) * erf(l_of_sigma(2)) * erf(l_of_sigma(3))
         first_time = .false.
      endif

      use_uniform = .not. use_uniform
      if(use_uniform) then
         ! use random position
         do i = 1, 3
            call urand(r)
            scaled_coord(i) = (2 * r - 1) * l_of_sigma(i)
         end do

      else
         ! for weighted insertion
         ! W = x^2
         ! get three N(0, 1) values
         do i = 1, 3
            do 
               r = nrand() / sqrt(2.0)
               if(abs(r) < l_of_sigma(i)) exit
            end do
            scaled_coord(i) = r
         end do

      endif
      sqsum = sum(scaled_coord(:) ** 2)
      ! from WHAM
      weight = weight / (uniform_ratio * 1 / z0 + (1 - uniform_ratio) * exp(-sqsum) / z1)
      com(:) = matmul(cell(:, :) , scaled_coord(:) / (l_of_sigma(:) * 2))
    end subroutine uniform_gauss_mixture
  end subroutine set_shift_com

  subroutine set_solute_com(insml, com)
    use engmain, only: numsite, mol_begin_index, mol_end_index, &
                       sitepos, sitemass
    use bestfit, only: com_shift, com_unshift
    implicit none
    integer, intent(in) :: insml
    real, intent(in) :: com(3)
    integer :: insb, inse, i, n
    real :: tempcom(3)
    
    n = numsite(insml)
    insb = mol_begin_index(insml)
    inse = mol_end_index(insml)
    
    call com_shift(n, sitepos(1:3, insb:inse), sitemass(insb:inse), tempcom)
    do i = 1, 3
       sitepos(i, insb:inse) = sitepos(i, insb:inse) + com(i)
    end do
  end subroutine set_solute_com

  subroutine apply_orientation(insml)
    use engmain, only: insorient, INSROT_RANDOM, INSROT_NOCHANGE, &
                       numsite, mol_begin_index, mol_end_index, &
                       sitepos, sitemass
    use quaternion, only: rotate_inplace
    use bestfit, only: com_shift, com_unshift
    implicit none
    integer, intent(in) :: insml
    integer :: insb, inse, i, n
    real :: com(3), randq(0:3)
    
    n = numsite(insml)
    insb = mol_begin_index(insml)
    inse = mol_end_index(insml)

    select case(insorient)
    case(INSROT_RANDOM)     ! random orientation
       call com_shift(n, sitepos(1:3, insb:inse), sitemass(insb:inse), com)

       ! get random quaternion (prob. success ~ 0.3)
       ! if this is really bad implement:
       ! Marsaglia, G. "Choosing a Point from the Surface of a Sphere."
       !   Ann. Math. Stat. 43, 645-646, 1972.
       do
          do i = 0, 3
             call urand(randq(i))
          end do
          randq(:) = randq(:) * 2.0 - 1.0
          if(sum(randq ** 2) < 1) exit
       end do

       randq(:) = randq(:) / sqrt(sum(randq ** 2)) ! set on unit hyper-sphere surface
       call rotate_inplace(numsite(insml), sitepos(1:3, insb:inse), randq)

       call com_unshift(n, sitepos(1:3, insb:inse), sitemass(insb:inse), com)
       return
    case(INSROT_NOCHANGE)   ! no orientational change
       return
    case default
       stop "Unknown insorient in apply_orientation"
    end select
  end subroutine apply_orientation
  !
  ! user-defined scheme to specify inserted molecule
  ! user may reject snapshot by specifying out_of_range to .true.
  ! or set coordinate in sitepos manually
  subroutine insscheme(insml, out_of_range)
    use engmain, only: hostspec, refspec, &
                       numsite, mol_begin_index, mol_end_index, &
                       specatm, sitemass, sitepos, cell
    use bestfit, only: center_of_mass
    implicit none
    integer, intent(in) :: insml
    logical, intent(inout) :: out_of_range
    if(out_of_range) return
    return
  end subroutine insscheme


  subroutine getsolute(caltype, insml, cntdst, stat_weight)
    use trajectory, only: open_trajectory, close_trajectory
    use engmain, only: slttype, maxins, numsite, bfcoord, stdout, slttrj, &
                       wgtins, sltwgt_file, sltwgt_io, &
                       insstructure, lwstr, upstr, &
    ! start of the extension for computing the conditional distributions
                       do_conditional, OrderPrm_read, &
    ! end of the extension for computing the conditional distributions
                       SLT_REFS_FLEX, INSSTR_NOREJECT, INSSTR_RMSD, YES
    use OUTname, only: OUTconfig, solute_trajectory
    use mpiproc
    use bestfit, only: rmsd_bestfit
    implicit none
    character(len=4),  intent(in) :: caltype
    integer,           intent(in) :: insml
    integer, optional, intent(in) :: cntdst
    real,    optional, intent(out) :: stat_weight
    logical, save :: read_weight
    logical, save :: consecutive_read = .false.  ! .true. for special purpose
    integer, save :: stmax, readmax
    real, dimension(:,:,:), allocatable, save :: solute_crd
    real, dimension(:),     allocatable, save :: solute_wgt
    real, dimension(:,:,:), allocatable :: read_crd
    real, dimension(:),     allocatable :: read_wgt
    real, dimension(:,:),   allocatable :: psite
    integer :: dumint, readcnt, iproc, ioerr
    real :: dumcl(3, 3), weight, rmsd
    logical :: reject
 
    if(slttype /= SLT_REFS_FLEX) call halt_with_error('ins_bug')
 
    if((.not. present(cntdst)) .and. (.not. present(stat_weight))) then
       select case(caltype)
       case('init')
    ! start of the extension for computing the conditional distributions
      ! note the consistency with getconf_parallel subroutine in setconf.F90
      ! caution: program will not work as expected when skpcnf > 1
          if((do_conditional == YES) .and. (OrderPrm_read == YES)) then
             consecutive_read = .true.
          endif
    ! end of the extension for computing the conditional distributions
          if(wgtins == YES) then
             read_weight = .true.
          else
             read_weight = .false.
          endif
          if(consecutive_read) then
             readmax = maxins
          else
             readmax = 1
          endif
          stmax = numsite(insml)

          allocate( solute_crd(3, stmax, readmax), solute_wgt(readmax) )
          if(myrank == 0) then
             call open_trajectory(solute_trajectory, slttrj)
             if(read_weight) open(unit = sltwgt_io, file = sltwgt_file, status = 'old')
          endif
          return
       case('last')
          deallocate( solute_crd, solute_wgt )
          if(myrank == 0) then
             call close_trajectory(solute_trajectory)
             if(read_weight) close(unit = sltwgt_io)
          endif
          return
       case('proc')
          call halt_with_error('ins_bug')
       case default
          stop "Incorrect caltype in getsolute"
       end select
    endif
 
    if((.not. present(cntdst)) .or. (.not. present(stat_weight)) .or. &
       (myrank >= nactiveproc)) call halt_with_error('ins_bug')
 
    ! get the configuration and structure-specific weight
    if((.not. consecutive_read) .or. (cntdst == 1)) then

       ! The following part has a similar program structure as the
       ! coordinate-reading part of getconf_parallel subroutine in setconf.F90
       if(myrank == 0) then               ! rank-0 to read from file
          allocate( read_crd(3, stmax, readmax), read_wgt(readmax) )
          allocate( psite(3, stmax) )
          if(size(psite) /= size(bfcoord)) call halt_with_error('ins_siz')

          do iproc = 1, nactiveproc
             do readcnt = 1, readmax
                reject = .true.
                do while(reject)
                   call OUTconfig(psite, dumcl, stmax, 0, 'solute', 'trjfl_read')
                   if(read_weight) then   ! weight read from a file
                      read(sltwgt_io, *, iostat = ioerr) dumint, weight
                      if(ioerr /= 0) then ! wrap around
                         rewind(sltwgt_io)
                         read(sltwgt_io, *, iostat = ioerr) dumint, weight
                         if(ioerr /= 0) then
                            write(stdout, *) " The weight file (", sltwgt_file, ") is ill-formed"
                            call mpi_setup('stop')
                            stop
                         endif
                      endif
                   else                   ! no structure-specific weight
                      weight = 1.0
                   endif

                   select case(insstructure)
                   case(INSSTR_NOREJECT)  ! no rejection of solute structure
                      reject = .false.
                   case(INSSTR_RMSD)      ! solute structure rejection with RMSD
                      if(size(psite) /= size(refslt_crd)) call halt_with_error('ins_siz')
                      if(stmax /= refslt_natom) call halt_with_error('ins_siz')
                      rmsd = rmsd_bestfit(refslt_natom, refslt_crd, &
                                          psite, refslt_weight)
                      if((lwstr <= rmsd) .and. (rmsd <= upstr)) reject = .false.
                   case default
                      stop "Unknown insstructure in getsolute"
                   end select
                enddo
                read_crd(1:3, 1:stmax, readcnt) = psite(1:3, 1:stmax)
                read_wgt(readcnt) = weight
             enddo

             if(iproc /= 1) then          ! send the data to other rank
#ifdef MPI
                call mpi_send(read_crd, 3 * stmax * readmax, mpi_double_precision, &
                              iproc - 1, tag_sltcrd, mpi_comm_world, ierror)
                call mpi_send(read_wgt, readmax, mpi_double_precision, &
                              iproc - 1, tag_sltwgt, mpi_comm_world, ierror)
#endif
             else                         ! rank-0 to use the data as read
                solute_crd(:, :, :) = read_crd(:, :, :)
                solute_wgt(:) = read_wgt(:)
             endif
          enddo
          deallocate( psite, read_crd, read_wgt )
       else                               ! non-0 rank to receive the data
#ifdef MPI
          call mpi_recv(solute_crd, 3 * stmax * readmax, mpi_double_precision, &
                        0, tag_sltcrd, mpi_comm_world, mpistatus, ierror)
          call mpi_recv(solute_wgt, readmax, mpi_double_precision, &
                        0, tag_sltwgt, mpi_comm_world, mpistatus, ierror)
#endif
       endif
    endif

    if(consecutive_read) then
       readcnt = cntdst
    else
       readcnt = 1
    endif
    bfcoord(1:3, 1:stmax) = solute_crd(1:3, 1:stmax, readcnt)
    stat_weight = solute_wgt(readcnt)
  end subroutine getsolute


  ! returns random value from [0,1)
  ! Any sane compiler implements random_number
  ! (which is included in fortran 95 standards)
  subroutine urand(rndm)          ! uniform random number generator
    implicit none
    real, intent(out) :: rndm
    call random_number(rndm)
  end subroutine urand

  ! Normal random variable N(0,1)
  ! uses Box-Muller method
  real function nrand()
    use engmain, only: PI
    implicit none
    real :: r1, r2
    call urand(r1)
    call urand(r2)
    ! get (0,1] instead of [0, 1)
    r1 = 1 - r1
    nrand = sqrt(-2.0 * log(r1)) * cos(2 * PI * r2)
  end function nrand

  subroutine urand_init(seed)
    use mpiproc, only: myrank
    implicit none
    integer, intent(in) :: seed
    integer :: seedsize
    integer, allocatable :: seedarray(:)

    call random_seed(size = seedsize)
    allocate(seedarray(seedsize))
    seedarray(:) = 1

    seedarray(1) = myrank + seed
    if(seed == 0) call system_clock(count = seedarray(1))

    call random_seed(put = seedarray)
    deallocate(seedarray)
  end subroutine urand_init


  subroutine com_aggregate(aggregate_center)
    use engmain, only: nummol, numsite, hostspec, moltype, &
                       mol_begin_index, mol_end_index, sitemass, sitepos
    use bestfit, only: center_of_mass
    implicit none
    real, intent(out) :: aggregate_center(3)
    real, dimension(:), allocatable   :: agg_mass
    real, dimension(:,:), allocatable :: agg_site
    integer :: num_aggsite, molb, mole, nsite, cnt, i
    num_aggsite = 0
    do i = 1, nummol
       if(any(moltype(i) == hostspec(:))) num_aggsite = num_aggsite + numsite(i)
    enddo
    allocate( agg_mass(num_aggsite), agg_site(3,num_aggsite) )
    cnt = 0
    do i = 1, nummol
       if(any(moltype(i) == hostspec(:))) then
          nsite = numsite(i)
          molb = mol_begin_index(i)
          mole = mol_end_index(i)
          agg_mass(cnt+1:cnt+nsite) = sitemass(molb:mole)
          agg_site(1:3, cnt+1:cnt+nsite) = sitepos(1:3, molb:mole)
          cnt = cnt + nsite
       endif
    enddo
    call center_of_mass(num_aggsite, agg_site, agg_mass, aggregate_center)
    deallocate( agg_mass, agg_site )
  end subroutine com_aggregate
   

  subroutine check_solute_configuration(insml, out_of_range)
    use engmain, only: insorigin, insposition, lwreg, upreg, &
                       numsite, mol_begin_index, mol_end_index, sitepos, &
                       INSORG_REFSTR, INSPOS_RMSD, INSPOS_GAUSS
    use mpiproc, only: halt_with_error
    use bestfit, only: rmsd_nofit
    implicit none
    integer, intent(in) :: insml
    logical, intent(inout) :: out_of_range
    integer :: ptb, pte
    real :: rmsd
    if(out_of_range) return
    if(insorigin == INSORG_REFSTR) then
       if((insposition /= INSPOS_RMSD) .and. (insposition /= INSPOS_GAUSS)) then
          call halt_with_error('ins_bug')
       endif
    else
       return
    endif
    ptb = mol_begin_index(insml)
    pte = mol_end_index(insml)
    if(numsite(insml) /= refslt_natom) call halt_with_error('ins_bug')
    rmsd = rmsd_nofit(refslt_natom, refslt_bestfit, &
                      sitepos(1:3, ptb:pte), refslt_weight)
    if((lwreg > rmsd) .or. (rmsd > upreg)) out_of_range = .true.
    return
  end subroutine check_solute_configuration
!
! fit to reference structure
  subroutine reffit
    use engmain, only: bfcoord, sitepos
    use bestfit, only: fit_a_rotate_b, fit
    implicit none
    integer :: i
    logical, save :: first_time = .true.
    real, dimension(:,:), allocatable :: hostcrd, fit_sltcrd
    if(first_time) then
       allocate( refslt_bestfit(3, refslt_natom) )
       first_time = .false.
    endif
    allocate( hostcrd(3, refhost_natom), fit_sltcrd(3, refslt_natom) )
    do i = 1, refhost_natom
       hostcrd(1:3, i) = sitepos(1:3, refhost_specatm(i))
    end do
    ! 1) fit the reference solvent (specified by refspec) structure from file
    !    onto the current solvent structure from MD
    !    in order to get the best-fit "reference solute configuration"
    call fit_a_rotate_b(refhost_natom, hostcrd, &
                        refhost_crd, refhost_weight, &
                        refslt_natom, refslt_crd, refslt_bestfit)
    ! 2) fit the solute structure read as the bfcoord variable
    call fit(refslt_natom, refslt_bestfit, bfcoord, refslt_weight, fit_sltcrd)
    do i = 1, refslt_natom
       sitepos(1:3, refslt_specatm(i)) = fit_sltcrd(1:3, i)
    end do
    deallocate( hostcrd, fit_sltcrd )
  end subroutine reffit
!
! loading the reference structure
  subroutine load_refstructure
    use engmain, only: nummol, moltype, numsite, sltlist, mol_begin_index, &
                       refspec, refstr_file, refstr_io
    use mpiproc, only: halt_with_error
    implicit none
    integer :: sltmol, atom_count, i, sid, stat
    real :: crd(3), wgt
    character(len=6) :: header

    refhost_natom = 0
    do i = 1, nummol
       if(any(moltype(i) == refspec(:))) refhost_natom = refhost_natom + numsite(i)
    enddo
    if(refhost_natom > 0) then
       allocate( refhost_crd(3, refhost_natom), refhost_weight(refhost_natom) )
       allocate( refhost_specatm(refhost_natom) )
       atom_count = 0
       do i = 1, nummol
          if(any(moltype(i) == refspec(:))) then
             do sid = 1, numsite(i)
                refhost_specatm(atom_count + sid) = mol_begin_index(i) + sid - 1
             end do
             atom_count = atom_count + numsite(i)
          endif
       end do
       if(atom_count /= refhost_natom) call halt_with_error('ins_bug')
    endif

    sltmol = sltlist(1)
    refslt_natom = numsite(sltmol)
    allocate( refslt_crd(3, refslt_natom), refslt_weight(refslt_natom) )
    allocate( refslt_specatm(refslt_natom) )
    do sid = 1, refslt_natom
       refslt_specatm(sid) = mol_begin_index(sltmol) + sid - 1
    end do

    open(unit = refstr_io, file = refstr_file, status = 'old', iostat = stat)
    if(stat /= 0) call halt_with_error('ins_ref')
    atom_count = 0    ! counts the number of lines with ATOM/HETATM header
    do
       read(refstr_io, '(A6)', end = 99) header
       if(header == 'ATOM  ' .or. header == 'HETATM') atom_count = atom_count + 1
    end do
99  rewind(refstr_io)
    if(atom_count /= refhost_natom + refslt_natom) call halt_with_error('ins_str')

    if(refhost_natom > 0) then
       do i = 1, refhost_natom
          call read_refPDF_weight(refhost_specatm(i), crd, wgt)
          refhost_crd(1:3, i) = crd(1:3)
          refhost_weight(i) = wgt
       end do
    endif
    do i = 1, refslt_natom
       call read_refPDF_weight(refslt_specatm(i), crd, wgt)
       refslt_crd(1:3, i) = crd(1:3)
       refslt_weight(i) = wgt
    end do
    close(refstr_io)
    return

  contains
    subroutine read_refPDF_weight(ati, crd, wgt)
      use engmain, only: sitemass
      implicit none
      integer, intent(in) :: ati
      real, intent(out) :: crd(3), wgt
      real, parameter :: massHe = 4.0026          ! atomic weight (helium)
      real :: refindex
      character(len=6) :: header
      do                       ! skip until ATOM/HETATM lines
         read(refstr_io, '(A6)', advance='no') header
         if(header == 'ATOM  ' .or. header == 'HETATM') exit
         read(refstr_io, *)
      end do
      read(refstr_io, '(24X, 3F8.3, F6.2)') crd(1:3), refindex
      if(refindex == 0.0) then ! not counted as an atom in reference structure
         wgt = 0.0
      else                     ! atom in reference with a weight given below
         ! initialize the weight as the atomic mass
         wgt = sitemass(ati)
         !
         ! user can implement his own special selection rule to mask fitting
         ! (e.g. by using B-factor, etc.)
         ! default: hydrogen is masked and the others have the same weight
         if(wgt > 0.95 * massHe) then    ! non-hydrogen
            ! wgt >= massHe, actually, where massHe is the helium atomic weight
            wgt = 1.0
         else                            ! hydrogen
            wgt = 0.0
         endif
         !
         ! comment out the following line if the mass weight is to be used
!        wgt = sitemass(ati)
         !
      endif
      return
    end subroutine read_refPDF_weight
  end subroutine load_refstructure
end module
