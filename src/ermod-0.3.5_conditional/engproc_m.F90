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

module engproc
  implicit none
  character(len=10), parameter :: numbers='0123456789'
  real, parameter :: tiny = 1.0e-30
  integer :: cntdst, slvmax
  integer :: maxdst
  integer :: tagslt
  integer, allocatable :: tagpt(:)

  ! flceng needs to be output in-order
  logical, allocatable :: flceng_stored(:)
  real, allocatable :: flceng(:, :)

contains
  !
  !  procedure for constructing energy distribution functions
  !
  subroutine enginit
    use engmain, only: numtype, nummol, engdiv, corrcal, selfcal, slttype, &
         moltype, sluvid, &
         ermax, numslv, uvmax, uvsoft, esmax, uvspec, &
         uvcrd, edens, ecorr, escrd, eself, &
         voffset, &
         aveuv, slnuv, avediv, minuv, maxuv, numslt, sltlist, &
         ene_confname, io_paramfile, io_flcuv, &
         ecdinfo_file, ecdinfo_io, ecdmesh_file, ecdmesh_io, &
    ! start of the extension for computing the conditional distributions
         maxins, stdout, &
         do_conditional, &
         OrderCrd_file, OrderCrd_io, OrderPrm_file, OrderPrm_io, &
         OrderPrm_read, OrderPrm_Values, OrderPrm_ArraySize, &
         order_species, order_size, &
         order_min, order_max, order_binwidth, &
         order_crd, edcnd, sluvcnd, crcnd, avuvcnd, cndnorm, &
    ! end of the extension for computing the conditional distributions
         SLT_SOLN, SLT_REFS_RIGID, SLT_REFS_FLEX, PT_SOLVENT, NO, YES
    use mpiproc, only: halt_with_error, warning, myrank
    implicit none
    real :: ecdmin, ecfmns, ecmns0, ecdcen, ecpls0, ecfpls, eccore, ecdmax
    real :: eclbin, ecfbin, ec0bin, finfac, ectmvl
    integer :: pecore
    ! pemax :  number of discretization of the solute-solvent energy
    ! pesoft : number of discretization in the soft interaction region
    !    they are constructed from other parameters (pemax = pesoft + pecore)
    integer pemax, pesoft
    integer :: ecprread, meshread, peread
    !
    real, parameter :: infty = 1.0e50      ! essentially equal to infinity
    integer, parameter :: rglmax = 5, large = 10000, too_large_ermax = 15000
    real :: factor, incre, cdrgvl(0:rglmax+1), ecpmrd(large)
    integer :: solute_moltype
    integer :: iduv, i, q, pti, regn, minrg, maxrg, uprgcd(0:rglmax+1), dummy
    integer, dimension(:), allocatable :: tplst
    real, dimension(:,:), allocatable  :: ercrd
    !
    integer :: param_err
    logical :: check_ok, start_line
    namelist /hist/ ecdmin, ecfmns, ecdcen, eccore, ecdmax, &
                    eclbin, ecfbin, ec0bin, finfac, pecore, &
                    ecprread, meshread, peread
    ! start of the extension for computing the conditional distributions
    namelist /conditional/ do_conditional, &
                           order_species, order_min, order_max, order_binwidth
    ! end of the extension for computing the conditional distributions
    !
    allocate( tplst(nummol) )
    numslt = 0
    do i = 1, nummol
       if(sluvid(i) > 0) then
          numslt = numslt + 1
          tplst(numslt) = i
          solute_moltype = moltype(i)
       endif
    enddo
    ! solute must have moltype value equal to solute_moltype
    if(any(sluvid(:) /= PT_SOLVENT .and. moltype(:) /= solute_moltype)) call halt_with_error('eng_typ')
    ! solvent must have moltype value not equal to solute_moltype
    if(any(sluvid(:) == PT_SOLVENT .and. moltype(:) == solute_moltype)) call halt_with_error('eng_typ')
    !
    ! consistency check between slttype and numslt (number of solute molecules)
    check_ok = .true.
    if(numslt <= 0) check_ok = .false.
    if((slttype == SLT_REFS_RIGID) .or. (slttype == SLT_REFS_FLEX)) then
       if(numslt /= 1) check_ok = .false.
    endif
    if(.not. check_ok) call halt_with_error('eng_num')
    !
    allocate( sltlist(numslt) )
    sltlist(1:numslt) = tplst(1:numslt)   ! list of solute molecules
    deallocate( tplst )
    !
    ! solute needs to be the last particle in reference system
    if((slttype == SLT_REFS_RIGID) .or. (slttype == SLT_REFS_FLEX)) then
       if(sltlist(1) /= nummol) call halt_with_error('eng_ins')
    endif
    !
    ! number of solvent species
    if(numslt == 1) then
       numslv = numtype - 1
    elseif(numslt > 1) then     ! solute can also be a solvent species
       numslv = numtype
    else
       call halt_with_error('eng_bug')
    endif
    !
    allocate( uvspec(nummol) )
    uvspec(1:nummol) = moltype(1:nummol)
    if(numslt == 1) then        ! solute is not present in reference solvent
      where(sluvid(:) /= PT_SOLVENT) uvspec(:) = 0 ! solute
      ! after the solute, slide the value of molecule type
      where(sluvid(:) == PT_SOLVENT .and. moltype(:) > solute_moltype) uvspec(:) = uvspec(:) - 1
    endif
    if(numslv /= maxval(uvspec(:))) call halt_with_error('eng_bug')
    !
    allocate( uvmax(numslv), uvsoft(numslv), ercrd(large, 0:numslv) )
    !
    ecprread = NO
    meshread = NO
    ermax = 0
    do pti = 0, numslv
       open(unit = io_paramfile, file = ene_confname, action = "read", iostat = param_err)
       if(param_err == 0) then
          read(io_paramfile, nml = hist)
          close(io_paramfile)
       else
          call halt_with_error('eng_per')
       end if
       if(peread == YES) ecprread = YES   ! deprecated

       ! read coordinate parameters from EcdInfo
       if(ecprread == YES) then
          open(unit = ecdinfo_io, file = ecdinfo_file, action = 'read', iostat = param_err)
          if(param_err /= 0) call halt_with_error('eng_eci')
          read(ecdinfo_io, *)             ! comment line
          do i = 1, large
             read(ecdinfo_io, *, end = 3109) q
             if(q == pti) then
                backspace(ecdinfo_io)
                if(pti == 0) then         ! solute self-energy
                   read(ecdinfo_io, *) iduv, ecpmrd(1:8)
                   pecore = 0               ! no core region for self-energy
                else                      ! solute-solvent interaction energy
                   read(ecdinfo_io, *) iduv, ecpmrd(1:9), pecore
                   ecdmax = ecpmrd(9)
                   if(pecore <= 1) call halt_with_error('eng_pcr')
                endif
                eclbin = ecpmrd(1)
                ecfbin = ecpmrd(2)
                ec0bin = ecpmrd(3)
                finfac = ecpmrd(4)
                ecdmin = ecpmrd(5)
                ecfmns = ecpmrd(6)
                ecdcen = ecpmrd(7)
                eccore = ecpmrd(8)
                if(eccore < tiny) pecore=0
                if(pecore == 1) call halt_with_error('eng_pcr')
                exit
             endif
          enddo
3109      continue
          close(ecdinfo_io)
       end if

       ectmvl = finfac * ecfbin
       ecdmin = ecdmin - ectmvl
       ecfmns = ecfmns - ectmvl
       ecmns0 = ecdcen - ectmvl
       ecpls0 = 2.0 * ecdcen - ecmns0
       ecfpls = 2.0 * ecdcen - ecfmns
       eccore = eccore + ecfpls - ecdcen

       cdrgvl(0) = ecdmin
       cdrgvl(1) = ecfmns
       cdrgvl(2) = ecmns0
       cdrgvl(3) = ecpls0
       cdrgvl(4) = ecfpls
       cdrgvl(5) = eccore
       cdrgvl(6) = ecdmax

       uprgcd(0) = 1
       do regn = 1, rglmax
          if((regn == 1) .or. (regn == 5)) factor = eclbin
          if((regn == 2) .or. (regn == 4)) factor = ecfbin
          if(regn == 3) factor = ec0bin
          incre = cdrgvl(regn) - cdrgvl(regn - 1)
          if(incre < factor) call halt_with_error('eng_ecd')
          iduv = nint(incre / factor)
          uprgcd(regn) = uprgcd(regn - 1) + iduv
       enddo

       pesoft = uprgcd(rglmax) - uprgcd(0)
       pemax = pesoft + pecore
       uprgcd(rglmax + 1) = pemax
       ercrd(uprgcd(0:(rglmax + 1)), pti) = cdrgvl(0:(rglmax + 1))

       if(pemax > large) call halt_with_error('eng_siz')

       if(pecore == 0) i = rglmax         ! no core region
       if(pecore >  0) i = rglmax + 1     ! explicit treatment of core region
       do regn = 1, i
          minrg = uprgcd(regn - 1)
          maxrg = uprgcd(regn)
          if(regn <= rglmax) then
             incre = ercrd(maxrg, pti) - ercrd(minrg, pti)
          endif
          if(regn == (rglmax + 1)) then   ! effective only when pecore > 0
             incre = log(ercrd(maxrg, pti) / ercrd(minrg, pti))
          endif
          factor = incre / real(maxrg - minrg)
          do iduv = minrg, (maxrg - 1)
             incre = factor * real(iduv - minrg)
             if(regn <= rglmax) then
                ercrd(iduv, pti) = ercrd(minrg, pti) + incre
             endif
             if(regn == (rglmax + 1)) then
                ercrd(iduv, pti) = ercrd(minrg, pti) * exp(incre)
             endif
          enddo
       enddo

       ! read coordinate meshes from EcdMesh
       if(meshread == YES) then
          open(unit = ecdmesh_io, file = ecdmesh_file, action = 'read', iostat = param_err)
          if(param_err /= 0) call halt_with_error('eng_ecm')
          start_line = .true.
          do i = 1, large
             read(ecdmesh_io, *, end = 3209) q
             if(q == pti) then
                backspace(ecdmesh_io)
                if( start_line ) then     ! 0-th line of the q-th species
                   start_line = .false.
                   iduv = 0
                   if(pti == 0) then      ! solute self-energy
                      read(ecdmesh_io, *) dummy, pemax
                      pecore = 0            ! no core region for self-energy
                   else                   ! solute-solvent interaction energy
                      read(ecdmesh_io, *) dummy, pemax, pecore
                   endif
                   pesoft = pemax - pecore
                else
                   iduv = iduv + 1        ! iduv-th line of the q-th species
                   read(ecdmesh_io, *) dummy, dummy, ercrd(iduv, pti)
                endif
             endif
          enddo
3209      continue
          close(ecdmesh_io)
          check_ok = .true.
          if(iduv /= pemax) check_ok = .false.
          if(pemax <= pecore) check_ok = .false.
          do iduv = 1, pemax - 1
             if(ercrd(iduv, pti) >= ercrd(iduv + 1, pti)) check_ok = .false.
          enddo
          if(.not. check_ok) call halt_with_error('eng_emf')
       endif

       if(pti == 0) then                  ! solute self-energy
          esmax = pesoft
       else                               ! solute-solvent interaction energy
          uvmax(pti) = pemax
          uvsoft(pti) = pesoft
          ermax = ermax + pemax
       endif
    enddo
    
    allocate( uvcrd(ermax), edens(ermax) )
    i = 0
    do pti = 1, numslv
       pemax = uvmax(pti)
       do iduv = 1, pemax
          i = i + 1
          uvcrd(i) = ercrd(iduv, pti)
       enddo
    enddo

    if(corrcal == YES) then
       if(ermax > too_large_ermax) call warning('emax')
       allocate( ecorr(ermax, ermax) )
    endif

    if(selfcal == YES)  then
      allocate( escrd(esmax), eself(esmax) )
      escrd(1:esmax) = ercrd(1:esmax,0)
    endif
    deallocate( ercrd )

    if(slttype == SLT_SOLN) allocate( aveuv(engdiv, numslv), slnuv(numslv) )
    allocate( avediv(engdiv, 2) )
    allocate( minuv(0:numslv), maxuv(0:numslv) )
     
    do pti = 0, numslv
       minuv(pti) = infty
       maxuv(pti) = -infty
    enddo
    voffset = -infty

    ! start of the extension for computing the conditional distributions
    do_conditional = NO                  ! default = don't do it
    order_min = 0.0
    order_max = 0.0
    order_binwidth = 0.0
    open(unit = io_paramfile, file = ene_confname, action = "read", iostat = param_err)
    if(param_err == 0) then
       read(io_paramfile, nml = conditional)
       close(io_paramfile)
    endif

    if(do_conditional == YES) then
       ! trajectory of order parameter
       inquire(file = OrderPrm_file, exist = check_ok)
       if( check_ok ) then
          OrderPrm_read = YES
          order_species = 0
          if(myrank == 0) open(unit = OrderPrm_io, file = OrderPrm_file, action = "read")
       endif
       if( (order_species <= 0) .or. (order_species > numslv) ) then
          order_species = 0
       else
          if(order_species == solute_moltype) stop " The solute cannot be the solvent species whose energy with solute is set to the order parameter"
!         i = count( mask = (moltype(1:nummol) == order_species) )
!         if(i /= 1) stop " When the interaction energy of a solvent species with solute is the order parameter, the number of that species needs to be 1"
       endif

       ! the following lines are taken from the engconst subroutine
       ! OrderPrm_ArraySize is equal to maxdst in the engconst subroutine
       select case(slttype)
       case(SLT_SOLN)
          OrderPrm_ArraySize = numslt
       case(SLT_REFS_RIGID, SLT_REFS_FLEX)
          OrderPrm_ArraySize = maxins
       end select
       allocate( OrderPrm_Values(OrderPrm_ArraySize) )

       ! setting the meshes for order parameter
       inquire(file = OrderCrd_file, exist = check_ok)
       if( check_ok ) then
          open(unit = OrderCrd_io, file = OrderCrd_file, action = "read")
          order_size = 0
          do
             read(OrderCrd_io, *, end = 997) q
             order_size = order_size + 1
          enddo
997       rewind(OrderCrd_io)
          allocate( order_crd(order_size) )
          do i = 1, order_size
             read(OrderCrd_io, *) q, order_crd(i)
          enddo
          close(OrderCrd_io)
       else
          if((order_max <= order_min) .or. (order_binwidth <= 0)) then
             write(stdout, '(A)') ' order_min, order_max or order_binwidth is incorrectly set'
             call halt_with_error('eng_ecd')
          endif
          order_size = nint( (order_max - order_min) / order_binwidth )
          allocate( order_crd(order_size) )
          do i = 1, order_size
             order_crd(i) = order_min + real(i - 1) * order_binwidth
          enddo
       endif
       check_ok = .true.
       do i = 1, order_size - 1
          if(order_crd(i) >= order_crd(i + 1)) check_ok = .false.
       enddo
       if(.not. check_ok) then
          write(stdout, '(A)') '  Order parameter is incorrectly meshed'
          call halt_with_error('eng_ecd')
       endif

       ! allocation of distribution functions
       allocate( edcnd(ermax, order_size) )
       if(corrcal == YES) then
          factor = real(ermax) * sqrt( real(order_size) )
          if(nint(factor) > too_large_ermax) call warning('emax')
          allocate( crcnd(ermax, ermax, order_size) )
       endif
       if(slttype == SLT_SOLN) then
          allocate( avuvcnd(engdiv, numslv, order_size) )
          allocate( sluvcnd(numslv, order_size) )
       endif
       allocate( cndnorm(order_size) )
    endif
    ! end of the extension for computing the conditional distributions

    call engclear

    ! Output for energy fluctuation
    if(myrank == 0) then
       if(slttype == SLT_SOLN) then
          ! open flcuv file
          open(unit = io_flcuv, file = 'flcuv.tt', status = 'new', action = 'write')
       else
          ! open progress file
          open(unit = io_flcuv, file = 'progress.tt', status = 'new', action = 'write')
       endif
    endif

    return
  end subroutine enginit

  subroutine engclear
    use engmain, only: corrcal, selfcal, slttype, SLT_SOLN, YES, &
    ! start of the extension for computing the conditional distributions
                       do_conditional, edcnd, sluvcnd, crcnd, cndnorm, &
    ! end of the extension for computing the conditional distributions
                       edens, ecorr, eself, slnuv, avslf, engnorm, engsmpl
    implicit none
    edens(:) = 0.0
    if(corrcal == YES) ecorr(:,:) = 0.0
    if(selfcal == YES) eself(:) = 0.0
    if(slttype == SLT_SOLN) slnuv(:) = 0.0
    avslf = 0.0
    engnorm = 0.0
    engsmpl = 0.0
    ! start of the extension for computing the conditional distributions
    if(do_conditional == YES) then
       edcnd(:,:) = 0.0
       if(corrcal == YES) crcnd(:,:,:) = 0.0
       if(slttype == SLT_SOLN) sluvcnd(:,:) = 0.0
       cndnorm(:) = 0.0
    endif
    ! end of the extension for computing the conditional distributions
    return
  end subroutine engclear

  subroutine engproc_cleanup
    use engmain, only: io_flcuv
    use mpiproc, only: myrank
    implicit none
    if(myrank == 0) then
       endfile(io_flcuv)
       close(io_flcuv)
    endif
  end subroutine engproc_cleanup
  !
  !
  subroutine engconst(stnum)
    use engmain, only: nummol, skpcnf, slttype, sluvid, &
                       maxins, numslv, numslt, cltype, cell, &
                       io_flcuv, &
                       SYS_NONPERIODIC, SYS_PERIODIC, &
                       EL_COULOMB, EL_PME, EL_PPPM, ES_NVT, ES_NPT, &
                       SLT_SOLN, SLT_REFS_RIGID, SLT_REFS_FLEX, &
                       PT_SOLVENT, PT_SOLUTE
    use reciprocal, only: recpcal_init, recpcal_spline_greenfunc, &
                          recpcal_prepare_solvent, recpcal_pppm_greenfunc
    use mpiproc                                                      ! MPI
    implicit none
    integer, intent(in) :: stnum
    integer :: i, k, irank
    real :: stat_weight_solute
    integer, dimension(:), allocatable :: tplst
    real, dimension(:),    allocatable :: uvengy
    logical, allocatable :: flceng_stored_g(:,:)
    real, allocatable :: flceng_g(:,:,:)
    real, save :: prevcl(3,3)
    logical :: skipcond
    logical, save :: pme_initialized = .false. ! NOTE: this variable is also used when PPPM is selected.

    call mpi_init_active_group(nactiveproc)
    call sanity_check_sluvid()

    ! for soln: maxdst is number of solute molecules
    ! for refs: maxdst is number of insertions
    select case(slttype)
    case(SLT_SOLN)
       maxdst = numslt
    case(SLT_REFS_RIGID, SLT_REFS_FLEX)
       maxdst = maxins
    end select

    allocate( tplst(nummol) )
    slvmax = 0
    do i = 1, nummol
       ! particle exists in trajectory (not a test particle)
       if((sluvid(i) == PT_SOLVENT) .or. (sluvid(i) == PT_SOLUTE)) then
          slvmax = slvmax + 1
          tplst(slvmax) = i
       end if
    enddo
    allocate( tagpt(slvmax) )
    allocate( uvengy(0:slvmax) )
    tagpt(1:slvmax) = tplst(1:slvmax)  ! and copied from tplst
    deallocate( tplst )

    allocate( flceng_stored(maxdst) )
    allocate( flceng(numslv, maxdst) )
    flceng_stored(:) = .false.
    flceng(:,:) = 0

    if(myrank < nactiveproc) then
       ! Initialize reciprocal space - grid and charges
       call get_inverted_cell
       if(cltype == EL_PME .or. cltype == EL_PPPM) then
          if(.not. pme_initialized) call recpcal_init(slvmax, tagpt)
          
          ! check whether cell size changes
          ! recpcal is called only when cell size differ
          call perf_time("kgrn")
          if((.not. pme_initialized) .or. &
             (any(prevcl(:,:) /= cell(:,:)))) then
             if (cltype == EL_PME) then
                call recpcal_spline_greenfunc()
             elseif (cltype == EL_PPPM) then
                call recpcal_pppm_greenfunc()
             endif
          end if
          call perf_time()
          prevcl(:,:) = cell(:,:)
          
          pme_initialized = .true.

          call perf_time("kpre")
          !$omp parallel do schedule(dynamic) private(k, i)
          do k = 1, slvmax
             i = tagpt(k)
             call recpcal_prepare_solvent(i)
          enddo
          call perf_time()
       endif

       ! cntdst is the pick-up no. of solute molecule from plural solutes (soln)
       ! cntdst is the iteration no. of insertion (refs)
       do cntdst = 1, maxdst
          call get_uv_energy(stnum, stat_weight_solute, uvengy(0:slvmax), skipcond)
          if(skipcond) cycle
          call update_histogram(stat_weight_solute, uvengy(0:slvmax))
       enddo

       select case(slttype)
       case(SLT_SOLN)                           ! for soln: output flceng
          allocate( flceng_stored_g(maxdst, nactiveproc) )
          allocate( flceng_g(numslv, maxdst, nactiveproc) )
          
#ifdef MPI
          ! gather flceng values to rank 0
          call mpi_gather(flceng_stored, maxdst, mpi_logical, &
               flceng_stored_g, maxdst, mpi_logical, &
               0, mpi_comm_activeprocs, ierror)
          call mpi_gather(flceng, numslv * maxdst, mpi_double_precision, &
               flceng_g, numslv * maxdst, mpi_double_precision, &
               0, mpi_comm_activeprocs, ierror)
#else
          flceng_stored_g(:,1) = flceng_stored(:)
          flceng_g(:,:,1) = flceng(:,:)
#endif
          
          if(myrank == 0) then
             do irank = 1, nactiveproc
                do cntdst = 1, maxdst
                   if(flceng_stored_g(cntdst, irank)) then
                      if(maxdst == 1) then
                         write(io_flcuv, 911) &
                                        (stnum + irank - 1) * skpcnf, &
                                        flceng_g(1:numslv, cntdst, irank)
                      else
                         write(io_flcuv, 912) cntdst, &
                                        (stnum + irank - 1) * skpcnf, &
                                        flceng_g(1:numslv, cntdst, irank)
                      endif
911                   format(i9, 9999f15.5)
912                   format(2i9, 9999f15.5)
                   endif
                enddo
             enddo
          endif
          deallocate( flceng_g, flceng_stored_g )

       case(SLT_REFS_RIGID, SLT_REFS_FLEX)      ! for refs: output progress
          if(myrank == 0) write(io_flcuv, *) ( &
                                        (stnum + irank - 1) * skpcnf, &
                                             irank = 1, nactiveproc)
       end select
#ifdef HAVE_FLUSH
       if (myrank == 0) flush(io_flcuv)
#endif
    endif

    deallocate( tagpt, uvengy )
    deallocate( flceng, flceng_stored )

    call mpi_finish_active_group()
    return
  end subroutine engconst
  !
  !
  subroutine engstore(stnum)
    !
    use engmain, only: maxcnf, skpcnf, engdiv, corrcal, selfcal, &
         slttype, wgtslf, &
         plmode, ermax, numslv, esmax, temp, &
         edens, ecorr, eself, &
         aveuv, slnuv, avediv, avslf, minuv, maxuv, &
         engnorm, engsmpl, voffset, voffset_initialized, &
    ! start of the extension for computing the conditional distributions
         do_conditional, &
         order_size, order_crd, edcnd, sluvcnd, crcnd, avuvcnd, cndnorm, &
    ! end of the extension for computing the conditional distributions
         SLT_SOLN, SLT_REFS_RIGID, SLT_REFS_FLEX, NO, YES
    use mpiproc                                                      ! MPI
    implicit none
    ! start of the extension for computing the conditional distributions
    integer :: order_prmid
    real :: order_param
    ! end of the extension for computing the conditional distributions
    integer :: stnum, pti, j, k, iduv, division
    character(len=9) :: engfile
    character(len=13) :: engfile2                         !masuda
    character(len=3) :: suffeng
    character(len=3) :: order_prmid2                      !masuda
    integer, parameter :: eng_io = 51, cor_io = 52, slf_io = 53
    integer, parameter :: ave_io = 54, wgt_io = 55, uvr_io = 56
    real :: voffset_local, voffset_scale
    real :: factor
    real, dimension(:), allocatable :: sve1, sve2
    real, dimension(:, :), allocatable :: sve3
    call mpi_rank_size_info                                          ! MPI
    !

    ! synchronize voffset
    if(wgtslf == YES) then
       voffset_local = voffset
#ifdef MPI
    ! MPI part starts here
       call mpi_allreduce(voffset_local, voffset, 1, &
            mpi_double_precision, mpi_max, mpi_comm_world, ierror)

       ! if uninitialized, use voffset so as not to pollute results with NaN
       if(.not. voffset_initialized) voffset_local = voffset

       ! scale histograms accoording to the maximum voffset
       select case(slttype)
       case(SLT_SOLN)
          voffset_scale = exp((voffset_local - voffset)/temp)
       case(SLT_REFS_RIGID, SLT_REFS_FLEX)
          voffset_scale = exp(-(voffset_local - voffset)/temp)
       end select

       engnorm = engnorm * voffset_scale
       if(selfcal == YES) eself(:) = eself(:) * voffset_scale
       if(slttype == SLT_SOLN) slnuv(:) = slnuv(:) * voffset_scale
       edens(:) = edens(:) * voffset_scale
       if(corrcal == YES) ecorr(:, :) = ecorr(:, :) * voffset_scale
    ! MPI part ends here
#endif
    endif

    ! Gather all information to Master node
#ifdef MPI
    ! MPI part starts here
    if(plmode == 2) then
       call mpi_reduce(avslf, factor, 1, &
            mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierror)
       avslf = factor
       call mpi_reduce(engnorm, factor, 1, &
            mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierror)
       engnorm = factor
       call mpi_reduce(engsmpl, factor, 1, &
            mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierror)
       engsmpl = factor
       if(selfcal == YES) then
         allocate( sve1(esmax) )
         sve1(1:esmax) = eself(1:esmax)
         call mpi_reduce(sve1, eself, esmax, &
              mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierror)
         deallocate( sve1 )
       endif
       if(slttype == SLT_SOLN) call mympi_reduce_real(slnuv, numslv, mpi_sum, 0)
    endif
    allocate( sve1(0:numslv), sve2(0:numslv) )
    sve1(0:numslv) = minuv(0:numslv)
    sve2(0:numslv) = maxuv(0:numslv)
    call mpi_reduce(sve1, minuv, (numslv + 1), &
         mpi_double_precision, mpi_min, 0, mpi_comm_world, ierror)
    call mpi_reduce(sve2, maxuv, (numslv + 1), &
         mpi_double_precision, mpi_max, 0, mpi_comm_world, ierror)
    deallocate( sve1, sve2 )
    allocate( sve1(ermax) )
    sve1(1:ermax) = edens(1:ermax)
    call mpi_reduce(sve1, edens, ermax, &
         mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierror)
    deallocate( sve1 )
    ! MPI part ends here
#endif
    edens(1:ermax) = edens(1:ermax) / engnorm
    if(corrcal == YES) then
#ifdef MPI
    ! MPI part starts here
       allocate( sve3(ermax, ermax) )
       sve3 = ecorr
       call mpi_reduce(sve3, ecorr, (ermax * ermax), &
            mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierror)
       deallocate( sve3 )
    ! MPI part ends here
#endif
       ecorr(1:ermax, 1:ermax) = ecorr(1:ermax, 1:ermax) / engnorm
    endif
    if(selfcal == YES) eself(1:esmax) = eself(1:esmax) / engnorm
    avslf = avslf / engnorm

    ! start of the extension for computing the conditional distributions
    if(do_conditional == YES) then
#ifdef MPI
       if(wgtslf == YES) then
          cndnorm(:) = cndnorm(:) * voffset_scale
          edcnd(:,:) = edcnd(:,:) * voffset_scale
          if(corrcal == YES) crcnd(:,:,:) = crcnd(:,:,:) * voffset_scale
          if(slttype == SLT_SOLN) sluvcnd(:,:) = sluvcnd(:,:) * voffset_scale
       endif
       if(plmode == 2) then
          do order_prmid = 1, order_size
             call mpi_reduce(cndnorm(order_prmid), factor, 1, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierror)
             cndnorm(order_prmid) = factor

             allocate( sve1(ermax) )
             sve1(1:ermax) = edcnd(1:ermax, order_prmid)
             call mpi_reduce(sve1(:), edcnd(:, order_prmid), ermax, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierror)
             deallocate( sve1 )

             if(corrcal == YES) then
                allocate( sve3(ermax, ermax) )
                sve3(:, :) = crcnd(:, :, order_prmid)
                call mpi_reduce(sve3(:, :), crcnd(:, :, order_prmid), &
                                (ermax * ermax), &
                     mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierror)
                deallocate( sve3 )
             endif

             if(slttype == SLT_SOLN) then
                call mympi_reduce_real(sluvcnd(:, order_prmid), numslv, mpi_sum, 0)
             endif
          enddo
       endif
#endif
       cndnorm(:) = cndnorm(:) / engnorm
       edcnd(:, :) = edcnd(:, :) / engnorm
       if(corrcal == YES) crcnd(:, :, :) = crcnd(:, :, :) / engnorm
       if(slttype == SLT_SOLN) sluvcnd(:, :) = sluvcnd(:, :) / engnorm

       if(myrank == 0) then
          division = stnum / (maxcnf / skpcnf / engdiv)
          if(slttype == SLT_SOLN) then
             do order_prmid = 1, order_size
                factor = cndnorm(order_prmid)
                if(factor > tiny) then
                   avuvcnd(division, :, order_prmid) = sluvcnd(:, order_prmid) / factor
                else
                   avuvcnd(division, :, order_prmid) = 0.0
                endif
             enddo
          endif
          if(engdiv == 1) then
             suffeng = '.tt'
          else
             j = division / 10
             k = mod(division, 10)
             suffeng = '.' // numbers(j+1:j+1) // numbers(k+1:k+1)
          endif

          select case(slttype)
          case(SLT_SOLN)
             engfile = 'nmcdsl' // suffeng
          case(SLT_REFS_RIGID, SLT_REFS_FLEX)
             engfile = 'nmcdrf' // suffeng
          end select
          open(unit = eng_io, file = engfile, form = "FORMATTED", action = 'write')
          do order_prmid = 1, order_size
             if(order_prmid < order_size) then
                order_param = (order_crd(order_prmid) + order_crd(order_prmid + 1) ) / 2.0
             else
                order_param = order_crd(order_size) + (order_crd(order_size) - order_crd(order_size - 1) ) / 2.0
             endif
             write(eng_io, '(f15.7,g25.17)') order_param, cndnorm(order_prmid)
          enddo
          endfile(eng_io)
          close(eng_io)

          select case(slttype)
          case(SLT_SOLN)
             engfile = 'ecndsl' // suffeng
          case(SLT_REFS_RIGID, SLT_REFS_FLEX)
             engfile = 'ecndrf' // suffeng
          end select
          open(unit = eng_io, file = engfile, form = "FORMATTED", action = 'write')
          do order_prmid = 1, order_size
             factor = cndnorm(order_prmid)
             if(factor > tiny) then
                edcnd(:, order_prmid) = edcnd(:, order_prmid) / factor
             else
                edcnd(:, order_prmid) = 0.0
             endif
             if(order_prmid > 0) write(eng_io, *)
             do iduv = 1, ermax
                call repval('intn', iduv, factor, pti)
                write(eng_io, '(g15.7,i5,g25.15)') factor, pti, edcnd(iduv, order_prmid)
             enddo
          enddo
          endfile(eng_io)
          close(eng_io)

          if(corrcal == YES) then
             select case(slttype)
             case(SLT_SOLN)
                engfile = 'crcdsl' // suffeng
             case(SLT_REFS_RIGID, SLT_REFS_FLEX)
                engfile = 'crcdrf' // suffeng
             end select
!             open(unit = cor_io, file = engfile, form = "UNFORMATTED", action = 'write')  !masuda
             do order_prmid = 1, order_size
                write(order_prmid2,'(i3.3)')  order_prmid   !masuda
                engfile2 = engfile // '.' //order_prmid2    !masuda
                open(unit = cor_io, file = engfile2, form ="UNFORMATTED",action= 'write')  !masuda
                factor = cndnorm(order_prmid)
                if(factor > tiny) then
                   crcnd(:, :, order_prmid) = crcnd(:, :, order_prmid) / factor
                else
                   crcnd(:, :, order_prmid) = 0.0
                endif
!                if(order_prmid > 0) write(cor_io)   !masuda
                write(cor_io) crcnd(:, :, order_prmid)
                endfile(cor_io)                      !masuda
                close(cor_io)                        !masuda
             enddo
          endif
!          endfile(cor_io)            !masuda
!          close(cor_io)              !masuda

          if((slttype == SLT_SOLN) .and. (stnum == maxcnf)) then
             open(unit = ave_io, file = 'avcnd.tt', action = 'write')
             do order_prmid = 1, order_size
                if(order_prmid > 0) write(ave_io, *)
                do k = 1, engdiv
                   write(ave_io, 751) k, avuvcnd(k, 1:numslv, order_prmid)
                enddo
             enddo
             endfile(ave_io)
             close(ave_io)
          endif
       endif
    endif
    ! end of the extension for computing the conditional distributions

    if(myrank /= 0) return                                            ! MPI
    !
    division = stnum / (maxcnf / skpcnf / engdiv)
    !
    avediv(division, 1) = engnorm / engsmpl
    select case(slttype)
    case(SLT_SOLN)
       aveuv(division, 1:numslv) = slnuv(1:numslv) / engnorm
       avediv(division, 2) = voffset - temp * log(avslf)
    case(SLT_REFS_RIGID, SLT_REFS_FLEX)
       avediv(division, 2) = voffset + temp * log(avslf)
    end select
    !
    if(division == engdiv) then
       if(slttype == SLT_SOLN) then
          open(unit = ave_io, file = 'aveuv.tt', action = 'write')
          do k = 1, engdiv
             write(ave_io, 751) k, aveuv(k, 1:numslv)
751          format(i5, 9999f15.5)
          enddo
          endfile(ave_io)
          close(ave_io)
       endif

       select case(slttype)
       case(SLT_SOLN)
          open(unit = wgt_io, file = 'weight_soln', action = 'write')
       case(SLT_REFS_RIGID, SLT_REFS_FLEX)
          open(unit = wgt_io, file = 'weight_refs', action = 'write')
       end select
       do k = 1, engdiv
          select case(wgtslf)
          case(NO)
             write(wgt_io, '(i5,e20.8)') k, avediv(k,1)
          case(YES)
             write(wgt_io, '(i5,e20.8,e20.7)') k, avediv(k,1), avediv(k,2)
          end select
       enddo
       endfile(wgt_io)
       close(wgt_io)

       open(uvr_io, file = 'uvrange.tt', action = 'write')
       write(uvr_io,'(A)') ' species     minimum        maximum'
       do pti = 0, numslv
          if(maxuv(pti) < 1.0e5) then
             write(uvr_io, '(i5,2f15.5)') pti, minuv(pti), maxuv(pti)
          else
             write(uvr_io, '(i5,f15.5,g18.5)') pti, minuv(pti), maxuv(pti)
          endif
       enddo
       endfile(uvr_io)
       close(uvr_io)
    endif

    if(engdiv == 1) then
       suffeng = '.tt'
    else
       j = division / 10
       k = mod(division, 10)
       suffeng = '.' // numbers(j+1:j+1) // numbers(k+1:k+1)
    endif
    !
    ! storing the solute-solvent pair-energy distribution function
    select case(slttype)
    case(SLT_SOLN)
       engfile = 'engsln' // suffeng
    case(SLT_REFS_RIGID, SLT_REFS_FLEX)
       engfile = 'engref' // suffeng
    end select
    open(unit = eng_io, file = engfile, form = "FORMATTED", action = 'write')
    do iduv = 1, ermax
       call repval('intn', iduv, factor, pti)
       write(eng_io, '(g15.7,i5,g25.15)') factor, pti, edens(iduv)
    enddo
    endfile(eng_io)
    close(eng_io)
    !
    ! storing the solvent-solvent correlation matrix in energy representation
    if(corrcal == YES) then
       select case(slttype)
       case(SLT_SOLN)
          engfile = 'corsln' // suffeng
       case(SLT_REFS_RIGID, SLT_REFS_FLEX)
          engfile = 'corref' // suffeng
       end select
       open(unit = cor_io, file = engfile, form = "UNFORMATTED", action = 'write')
       write(cor_io) ecorr
       endfile(cor_io)
       close(cor_io)
    endif
    !
    ! storing the distribution function of the self-energy of solute
    if(selfcal == YES) then
       engfile = 'slfeng' // suffeng
       open(unit = slf_io, file = engfile, form = "FORMATTED", action = 'write')
       do iduv = 1, esmax
          call repval('self', iduv, factor)
          write(slf_io, '(g15.7,g25.15)') factor, eself(iduv)
       enddo
       endfile(slf_io)
       close(slf_io)
    endif
    !
    return
  end subroutine engstore

  ! Calculate interaction energy between solute and solvent
  subroutine get_uv_energy(stnum, stat_weight_solute, uvengy, out_of_range)
    use engmain, only: maxcnf, skpcnf, slttype, sltlist, cltype, &
                       EL_COULOMB, EL_PME, EL_PPPM, &
                       SLT_SOLN, SLT_REFS_RIGID, SLT_REFS_FLEX
    use ptinsrt, only: instslt
    use realcal, only: realcal_prepare, realcal_proc, realcal_self, &
                       realcal_bare, realcal_cleanup
    use reciprocal, only: recpcal_prepare_solute, recpcal_self_energy, &
                          recpcal_energy
    use mpiproc                                                      ! MPI
    implicit none
    integer, intent(in) :: stnum
    real, intent(out) :: uvengy(0:slvmax), stat_weight_solute
    logical, intent(out) :: out_of_range

    integer :: i, k
    integer(8) :: current_solute_hash
    integer(8) :: solute_hash = 0
    real :: pairep, residual, factor
    real, save :: usreal
    logical, save :: initialized = .false.

    out_of_range = .false.
    ! determine / pick solute structure
    select case(slttype) 
    case(SLT_SOLN)
       tagslt = sltlist(cntdst)
       stat_weight_solute = 1.0
       call check_mol_configuration(out_of_range)
       if(out_of_range) return
    case(SLT_REFS_RIGID, SLT_REFS_FLEX)
       tagslt = sltlist(1)
       if(.not. initialized) call instslt('init')
       initialized = .true.
       call instslt('proc', cntdst, stat_weight_solute)
       if((stnum == maxcnf/skpcnf) .and. (cntdst == maxdst)) call instslt('last')
    end select

    ! At this moment all coordinate in the system is determined
    call realcal_prepare

    uvengy(:) = 0
    ! Calculate system-wide values
    if(cltype == EL_PME .or. cltype == EL_PPPM) then
       call recpcal_prepare_solute(tagslt)
       call perf_time("rblk")
       call realcal_proc(tagslt, tagpt, slvmax, uvengy)
       call perf_time()
       call perf_time("kslf")
       call recpcal_self_energy(uvengy(0))
       call perf_time()
    endif

    ! solute-solute self energy
    pairep = 0.0
    residual = 0.0
    current_solute_hash = get_solute_hash() ! FIXME: if this tuns into a bottleneck, add conditionals
    if(current_solute_hash == solute_hash .or. &
       (slttype == SLT_REFS_RIGID .and. solute_hash /= 0)) then 
       ! For refs calculation, the configuration of solute may change with
       ! random translation and/or rotation upon insertion, 
       ! though the self energy will not change.
       pairep = usreal ! reuse
    else
       call perf_time("rslf")
       call realcal_self(tagslt, pairep) ! calculate self-interaction
       call perf_time()
       usreal = pairep
    endif
    solute_hash = current_solute_hash
    call residual_ene(tagslt, tagslt, residual)
    uvengy(0) = uvengy(0) + pairep + residual

    ! solute-solvent pair
    !$omp parallel do schedule(dynamic) private(i, k, pairep, factor)
    do k = 1, slvmax
       i = tagpt(k)
       if(i == tagslt) cycle

       pairep = 0
       factor = 0
       if(cltype == EL_PME .or. cltype == EL_PPPM) then
          ! called only when PME or PPPM, non-self interaction
          call residual_ene(tagslt, i, pairep)
          call recpcal_energy(tagslt, i, factor)
          pairep = pairep + factor
       else                      ! Bare coulomb solute-solvent interaction
          call realcal_bare(tagslt, i, pairep)
       endif
       !$omp atomic
       uvengy(k) = uvengy(k) + pairep
    enddo

    call realcal_cleanup
  end subroutine get_uv_energy

  subroutine update_histogram(stat_weight_solute, uvengy)
    use engmain, only: wgtslf, estype, slttype, corrcal, selfcal, ermax, &
    ! start of the extension for computing the conditional distributions
                       numslv, moltype, stdout, &
                       do_conditional, OrderPrm_read, OrderPrm_Values, &
                       order_species, order_size, &
                       order_crd, edcnd, sluvcnd, crcnd, cndnorm, &
    ! end of the extension for computing the conditional distributions
                       volume, temp, uvspec, &
                       slnuv, avslf, minuv, maxuv, &
                       edens, ecorr, eself, &
                       stat_weight_system, engnorm, engsmpl, &
                       voffset, voffset_initialized, &
                       SLT_SOLN, SLT_REFS_RIGID, SLT_REFS_FLEX, &
                       ES_NVT, ES_NPT, NO, YES
    use mpiproc
    implicit none
    ! start of the extension for computing the conditional distributions
    integer :: order_prmid
    real :: order_param
    ! end of the extension for computing the conditional distributions
    real, intent(in) :: uvengy(0:slvmax), stat_weight_solute
    integer, allocatable :: insdst(:), engdst(:)
    integer :: i, k, q, iduv, iduvp, pti
    real :: factor, engnmfc, pairep, total_weight

    allocate( insdst(ermax), engdst(ermax) )

    select case(wgtslf)
    case(NO)
       engnmfc = 1.0
    case(YES)
       if(.not. voffset_initialized) then
          voffset = uvengy(0)
          voffset_initialized = .true.
       endif
       factor = uvengy(0) - voffset  ! shifted by offset
       select case(slttype)
       case(SLT_SOLN)
          engnmfc = exp(factor / temp)
       case(SLT_REFS_RIGID, SLT_REFS_FLEX)
          engnmfc = exp(- factor / temp)
       end select
    case default
       stop "Unknown wgtslf"
    end select
    total_weight = stat_weight_system * stat_weight_solute
    if(estype == ES_NPT) call volcorrect(total_weight)
    engnmfc = engnmfc * total_weight
    !
    engnorm = engnorm + engnmfc      ! normalization factor
    engsmpl = engsmpl + 1.0          ! number of sampling
    avslf = avslf + total_weight     ! normalization without solute self-energy

    ! self energy histogram
    if(selfcal == YES) then
      call getiduv(0, uvengy(0), iduv)
      eself(iduv) = eself(iduv) + engnmfc
    endif
    minuv(0) = min(minuv(0), uvengy(0))
    maxuv(0) = max(maxuv(0), uvengy(0))

    insdst(:) = 0                              ! instantaneous histogram
    flceng(:, cntdst) = 0.0                    ! sum of solute-solvent energy
    flceng_stored(cntdst) = .true.
    do k = 1, slvmax
       i = tagpt(k)
       if(i == tagslt) cycle
       pti = uvspec(i)
       if(pti <= 0) call halt_with_error('eng_cns')

       pairep = uvengy(k)
       call getiduv(pti, pairep, iduv)

       ! instantaneous histogram of solute-solvent energy
       insdst(iduv) = insdst(iduv) + 1
       ! sum of solute-solvent energy
       flceng(pti, cntdst) = flceng(pti, cntdst) + pairep

       ! minimum and maxmium of solute-solvent energy
       minuv(pti) = min(minuv(pti), pairep)
       maxuv(pti) = max(maxuv(pti), pairep)
    enddo

    if(slttype == SLT_SOLN) then
       slnuv(:) = slnuv(:) + flceng(:, cntdst) * engnmfc
    endif

    do iduv = 1, ermax
       k = insdst(iduv)
       if(k == 0) cycle
       edens(iduv) = edens(iduv) + engnmfc * real(k)
    enddo
    if(corrcal == YES) then
       do iduv = 1, ermax
          k = insdst(iduv)
          if(k == 0) cycle

          do iduvp = 1, ermax
             q = insdst(iduvp)
             if(q == 0) cycle

             ecorr(iduvp,iduv) = ecorr(iduvp,iduv) + engnmfc * real(k) * real(q)
          enddo
       enddo
    endif

    ! start of the extension for computing the conditional distributions
    if(do_conditional == YES) then
       if( (1 <= order_species) .and. (order_species <= numslv) ) then
          order_param = sum( uvengy(:), mask = (moltype(:) == order_species) )
       else
          if(OrderPrm_read == YES) then
             order_param = OrderPrm_Values(cntdst)
          else
             ! user-defined setting of order parameter
          endif
       endif
       call binsearch(order_crd(1:order_size), order_size, order_param, order_prmid)
       ! smaller than the minimum mesh of order parameter
       if(order_prmid < 1) then
          write(stdout, '(A,g12.4,A)') '  Value of order parameter is ', order_param, ' and is too small'
          call halt_with_error('eng_bug')
       endif
       ! larger than the maximum mesh of order parameter
       if(order_prmid > order_size) call halt_with_error('eng_bug')
       cndnorm(order_prmid) = cndnorm(order_prmid) + engnmfc
       if(slttype == SLT_SOLN) then
          sluvcnd(:, order_prmid) = sluvcnd(:, order_prmid) + flceng(:, cntdst) * engnmfc
       endif
       do iduv = 1, ermax
          k = insdst(iduv)
          if(k == 0) cycle
          edcnd(iduv, order_prmid) = edcnd(iduv, order_prmid) + engnmfc * real(k)
       enddo
       if(corrcal == YES) then
          do iduv = 1, ermax
             k = insdst(iduv)
             if(k == 0) cycle
             do iduvp = 1, ermax
                q = insdst(iduvp)
                if(q == 0) cycle
                crcnd(iduvp,iduv, order_prmid) = crcnd(iduvp,iduv, order_prmid) + engnmfc * real(k) * real(q)
             enddo
          enddo
       endif
    endif
    ! end of the extension for computing the conditional distributions

    deallocate( insdst, engdst )
  end subroutine update_histogram

  !
  subroutine residual_ene(i, j, pairep)
    use engmain, only: screen, volume, mol_charge, cltype, EL_COULOMB, PI
    implicit none
    integer, intent(in) :: i, j
    real, intent(inout) :: pairep
    real :: epcl
    if(cltype == EL_COULOMB) return
    epcl = PI * mol_charge(i) * mol_charge(j) / screen / screen / volume
    if(i == j) epcl = epcl / 2.0   ! self-interaction
    pairep = pairep - epcl
  end subroutine residual_ene
  !
  subroutine volcorrect(weight)
    use engmain, only:  sluvid, cltype, screen, mol_charge, volume, temp, &
                        EL_EWALD, EL_PME, EL_PPPM, PT_SOLVENT, PT_SOLUTE, PI
    implicit none
    real, intent(inout) :: weight
    real :: total_charge, factor
    weight = weight * volume
    if((cltype == EL_EWALD) .or. (cltype == EL_PME) &
         .or. (cltype == EL_PPPM)) then  ! Ewald, PME and PPPM
       ! sum of charges over physical particles
       total_charge = sum( mol_charge, mask = ((sluvid == PT_SOLVENT) &
                                          .or. (sluvid == PT_SOLUTE)) )
       factor = PI * (total_charge ** 2)  / (screen ** 2 ) / volume / 2.0
       weight = weight * exp(factor / temp)
    endif
    return
  end subroutine volcorrect
  !
  subroutine update_cell_info()
    use engmain, only: cell, celllen
    implicit none
    integer :: i
    do i = 1, 3
       celllen(i) = sqrt(sum(cell(:, i) ** 2))
    enddo
  end subroutine update_cell_info

  subroutine get_inverted_cell
    use engmain, only:  cell, invcl, volume
    implicit none
    volume = cell(1,1) * cell(2,2) * cell(3,3) &
           + cell(1,2) * cell(2,3) * cell(3,1) &
           + cell(1,3) * cell(2,1) * cell(3,2) &
           - cell(1,3) * cell(2,2) * cell(3,1) &
           - cell(1,2) * cell(2,1) * cell(3,3) &
           - cell(1,1) * cell(2,3) * cell(3,2)
    invcl(1,1) = cell(2,2) * cell(3,3) - cell(2,3) * cell(3,2)
    invcl(1,2) = cell(1,3) * cell(3,2) - cell(1,2) * cell(3,3)
    invcl(1,3) = cell(1,2) * cell(2,3) - cell(1,3) * cell(2,2)
    invcl(2,1) = cell(2,3) * cell(3,1) - cell(2,1) * cell(3,3)
    invcl(2,2) = cell(1,1) * cell(3,3) - cell(1,3) * cell(3,1)
    invcl(2,3) = cell(1,3) * cell(2,1) - cell(1,1) * cell(2,3)
    invcl(3,1) = cell(2,1) * cell(3,2) - cell(2,2) * cell(3,1)
    invcl(3,2) = cell(1,2) * cell(3,1) - cell(1,1) * cell(3,2)
    invcl(3,3) = cell(1,1) * cell(2,2) - cell(1,2) * cell(2,1)

    invcl(1:3, 1:3) = invcl(1:3, 1:3) / volume
    call update_cell_info
  end subroutine get_inverted_cell


  ! binsearch returns the smallest index (ret) which satisfies
  ! coord(ret) >= v
  ! Special cases are as follows:
  ! v < coord(1)  ==>  ret = 0 (out of range)
  ! v > coord(n)  ==>  ret = n (infinite)
  subroutine binsearch(coord, n, v, ret)
    implicit none 
    real, intent(in) :: coord(n)
    integer, intent(in) :: n
    real, intent(in) :: v
    integer, intent(out) :: ret
    integer :: rmin, rmax, rmid
    if(v < coord(1)) then
       ret = 0
       return
    endif
    if(v > coord(n)) then
       ret = n
       return
    endif

    rmin = 1
    rmax = n + 1
    do
       if(rmax - rmin <= 1) then
          exit
       endif
       rmid = (rmin + rmax - 1) / 2
       if(v > coord(rmid)) then
          rmin = rmid + 1
       else
          rmax = rmid + 1
       endif
    enddo
    ret = rmin - 1
  end subroutine binsearch

  ! returns the position of the bin corresponding to energy coordinate value
  subroutine getiduv(pti, engval, iduv)
    use engmain, only: slttype, uvmax, uvsoft, uvcrd, esmax, escrd, &
                       stdout, SLT_SOLN
    use mpiproc, only: halt_with_error, warning
    implicit none
    integer, intent(in) :: pti
    real, intent(in) :: engval 
    integer, intent(out) :: iduv
    integer :: idmin, idnum, idpti
    real, parameter :: warn_threshold = 1.0e3
    logical, save :: warn_bin_firsttime = .true.

    if(pti >  0) then            ! solute-solvent interaction
       if(pti == 1) then         ! first solvent species
          idmin = 0
       elseif(pti > 1) then      ! second and later solvent species
          idmin = sum( uvmax(1:(pti - 1)) )
       else                      ! bug --- pti must not be smaller than 1
          call halt_with_error('eng_bug')
       endif
       idnum = uvmax(pti)
       call binsearch(uvcrd((idmin + 1):(idmin + idnum)), idnum, engval, idpti)
    elseif(pti == 0) then        ! solute self-energy
       idmin = 0
       idnum = esmax
       call binsearch(escrd(1:idnum), idnum, engval, idpti)
    else                         ! bug --- pti must be non-negative
       call halt_with_error('eng_bug')
    endif

    ! inappropriate setting of energy coordinate
    ! smaller than the minimum energy mesh
    if(idpti <= 0) then
       write(stdout, '(A,g12.4,A,i3,A)') '  energy of ', engval, ' for ', pti, '-th species'
       call halt_with_error('eng_min')
    endif
    ! larger than the maximum energy mesh
    if(idpti > idnum) call halt_with_error('eng_bug')
    ! only for solute-solvent interaction in solution system
    ! larger than the maximum of soft part (linear-graduation part)
    if((slttype == SLT_SOLN) .and. (pti > 0)) then
       if(idpti > uvsoft(pti)) then
          write(stdout, '(A,g12.4,A,i3,A)') '  energy of ', engval, ' for ', pti, '-th species'
          call halt_with_error('eng_sft')
       endif
    endif

    ! Warning if the energy exceeds the maximum binning region and pecore = 0
    ! Since it is hard to see the pecore value at this point,
    ! the energy value is simply shown as a warning
    if((idpti == idnum) .and. (engval < warn_threshold) &
                        .and. (warn_bin_firsttime)) then
       write(stdout, '(A,g12.4,A,i3,A)') '  energy of ', engval, ' for ', pti, '-th species'
       call warning('mbin')
       warn_bin_firsttime = .false.
    endif

    iduv = idmin + idpti

    return
  end subroutine getiduv

  ! Check system consistency: either test particle or solvent must exist
  subroutine sanity_check_sluvid()
    use engmain, only: slttype, nummol, sluvid, SLT_SOLN, SLT_REFS_RIGID, SLT_REFS_FLEX
    use mpiproc, only: halt_with_error
    implicit none

    ! sanity check
    if(any(sluvid(:) < 0) .or. any(sluvid(:) > 3)) call halt_with_error('eng_bug')
    select case(slttype)
    case(SLT_SOLN)
       ! sluvid should be 0 (solvent) or 1 (solute)
       if(any(sluvid(:) >= 2)) call halt_with_error('eng_par')
    case(SLT_REFS_RIGID, SLT_REFS_FLEX)
       ! sluvid should be 0 (solvent), 2, 3 (test particles)
       if(any(sluvid(:) == 1)) call halt_with_error('eng_par')
       ! solvent must exist
       if(all(sluvid(:) /= 0)) call halt_with_error('eng_par')
    end select

    ! solute / test particle must exist
    if(all(sluvid(:) == 0)) call halt_with_error('eng_par')
  end subroutine sanity_check_sluvid

  ! Check whether molecule is within specified region
  subroutine check_mol_configuration(out_of_range)
    use engmain, only: insposition, insstructure, lwreg, upreg, lwstr, upstr, &
                       numsite, mol_begin_index, mol_end_index, &
                       sitepos, boxshp, invcl, celllen, &
                       SYS_NONPERIODIC, &
                       INSPOS_RANDOM, INSPOS_NOCHANGE, &
                       INSPOS_SPHERE, &
                       INSPOS_SLAB_GENERIC, INSPOS_SLAB_SYMMETRIC, &
                       INSPOS_RMSD, INSPOS_GAUSS, &
                       INSSTR_NOREJECT, INSSTR_RMSD
    use ptinsrt, only: refhost_natom, refhost_specatm, &
                       refhost_crd, refhost_weight, &
                       refslt_natom, refslt_crd, refslt_weight, &
                       insscheme
    use bestfit, only: fit_a_rotate_b, rmsd_nofit, rmsd_bestfit
    use mpiproc, only: halt_with_error
    implicit none
    logical, intent(out) :: out_of_range
    real :: dx(3), distance
    integer :: i, ptb, pte
    real, dimension(:,:), allocatable :: hostcrd, refslt_bestfit

    out_of_range = .false.

    select case(insposition)
    case(INSPOS_RANDOM)                               ! random
       ! do nothing
    case(INSPOS_NOCHANGE)                             ! fixed configuration
       ! do nothing
    case(INSPOS_SPHERE)                               ! sphere geometry
       ! The following has the same structure as the corresponding part
       ! in the set_shift_com subroutine within insertion.F90
       call relative_com(tagslt, dx)
       distance = sqrt(dot_product(dx, dx))
       if((lwreg > distance) .or. (distance > upreg)) then
          out_of_range = .true.     ! configuration is rejected
          return
       endif
    case(INSPOS_SLAB_GENERIC, INSPOS_SLAB_SYMMETRIC)  ! slab configuration
       ! constrained to z-axis
       ! The following has the same structure as the corresponding part
       ! in the set_shift_com subroutine within insertion.F90
       if(boxshp == SYS_NONPERIODIC) call halt_with_error('eng_slb')
       call relative_com(tagslt, dx)
       distance = dot_product(invcl(3,:), dx(:)) * celllen(3)
       if(insposition == INSPOS_SLAB_SYMMETRIC) then  ! symmetric bilayer
          distance = abs(distance)
       endif
       if((lwreg > distance) .or. (distance > upreg)) then
          out_of_range = .true.     ! configuration is rejected
          return
       endif
    case(INSPOS_RMSD)                                 ! comparison to reference
       ptb = mol_begin_index(tagslt)
       pte = mol_end_index(tagslt)
       if(numsite(tagslt) /= refslt_natom) call halt_with_error('eng_bug')
       allocate( hostcrd(3, refhost_natom), refslt_bestfit(3, refslt_natom) )
       ! The following has the same structure as the corresponding part
       ! in the reffit subroutine within insertion.F90
       do i = 1, refhost_natom
          hostcrd(1:3, i) = sitepos(1:3, refhost_specatm(i))
       enddo
       call fit_a_rotate_b(refhost_natom, hostcrd, &
                           refhost_crd, refhost_weight, &
                           refslt_natom, refslt_crd, refslt_bestfit)
       ! The following has the same structure as the corresponding part
       ! in the check_solute_configuration subroutine within insertion.F90
       distance = rmsd_nofit(refslt_natom, refslt_bestfit, &
                             sitepos(1:3, ptb:pte), refslt_weight)
       deallocate( hostcrd, refslt_bestfit )
       if((lwreg > distance) .or. (distance > upreg)) then
          out_of_range = .true.     ! configuration is rejected
          return
       endif
    case(INSPOS_GAUSS)                                ! comparison to reference
       ! Exeperimental and under construction...  What is to be written?
    case default
       stop "Unknown insposition in check_mol_configuration"
    end select

    select case(insstructure)
    case(INSSTR_NOREJECT)    ! no rejection of solute structure
       ! do nothing
    case(INSSTR_RMSD)        ! solute structure rejection with RMSD
       ! The following has the same structure as the corresponding part
       ! in the getsolute subroutine within insertion.F90
       ptb = mol_begin_index(tagslt)
       pte = mol_end_index(tagslt)
       if(numsite(tagslt) /= refslt_natom) call halt_with_error('eng_bug')
       distance = rmsd_bestfit(refslt_natom, refslt_crd, &
                               sitepos(1:3, ptb:pte), refslt_weight)
       if((lwstr > distance) .or. (distance > upstr)) then
          out_of_range = .true.     ! structure is rejected
          return
       endif
    case default
       stop "Unknown insstructure in check_mol_configuration"
    end select

    call insscheme(tagslt, out_of_range)
    return

  contains
    subroutine relative_com(tagpt, dx)
      use engmain, only: numsite, mol_begin_index, mol_end_index, &
                         sitemass, sitepos
      use ptinsrt, only: com_aggregate
      use bestfit, only: center_of_mass
      implicit none
      integer, intent(in) :: tagpt
      real, intent(out) :: dx(3)
      integer ptb, pte, stmax
      real :: solute_com(3), aggregate_com(3)
      real, dimension(:), allocatable   :: ptmass
      real, dimension(:,:), allocatable :: ptsite
      stmax = numsite(tagpt)
      ptb = mol_begin_index(tagpt)
      pte = mol_end_index(tagpt)
      allocate( ptmass(stmax), ptsite(3, stmax) )
      ptmass(1:stmax) = sitemass(ptb:pte)
      ptsite(1:3, 1:stmax) = sitepos(1:3, ptb:pte)
      call center_of_mass(stmax, ptsite, ptmass, solute_com)  ! solute COM
      call com_aggregate(aggregate_com)                       ! aggregate COM
      dx(1:3) = solute_com(1:3) - aggregate_com(1:3)
      deallocate( ptmass, ptsite )
      return
    end subroutine relative_com
  end subroutine check_mol_configuration

  ! get the hashed function of solute coordinate
  integer(8) function get_solute_hash()
    use utility, only: hash
    use engmain, only: sitepos, mol_begin_index, numsite
    implicit none

    get_solute_hash = hash(sitepos(1:3, mol_begin_index(tagslt):(mol_begin_index(tagslt+1) - 1)), numsite(tagslt) * 3)
  end function get_solute_hash

  
  subroutine repval(uvtype, iduv, engcoord, spec)
    use engmain, only: numslv, uvmax, uvsoft, uvcrd, esmax, escrd
    use mpiproc, only: halt_with_error
    implicit none
    character(len=4), intent(in) :: uvtype
    integer, intent(in) :: iduv
    real, intent(out) :: engcoord
    integer, intent(out), optional :: spec
    integer :: idpt, cnt, idmin, idmax, idsoft
    select case(uvtype)
    case('self')
       if(present(spec)) call halt_with_error('eng_bug')
       if(iduv <  esmax) engcoord = (escrd(iduv) + escrd(iduv+1)) / 2.0
       if(iduv == esmax) engcoord = escrd(esmax)
       if(iduv >  esmax) call halt_with_error('eng_ecd')
    case('intn')
       if(.not. present(spec)) call halt_with_error('eng_bug')
       idmin = 0
       do cnt = 1, numslv
          idpt = idmin + uvmax(cnt)
          if(iduv <= idpt) exit
          idmin = idpt
       enddo
       spec = cnt
       idsoft = uvsoft(spec)
       idmax = uvmax(spec)
       idpt = iduv - idmin
       if((idpt < 1) .or. (idpt > idmax)) call halt_with_error('eng_ecd')
       if(idpt <= idsoft) then   ! linear graduation
          engcoord = (uvcrd(iduv) + uvcrd(iduv+1)) / 2.0
       else                      ! logarithmic graduation
          if(idpt <  idmax) engcoord = sqrt(uvcrd(iduv) * uvcrd(iduv+1))
          if(idpt == idmax) engcoord = uvcrd(iduv)
       endif
    case default
       stop "Unknown caltype in repval"
    end select
    return
  end subroutine repval
end module engproc
