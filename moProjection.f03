INCLUDE 'moProjection_mod.f03'

      program moProjection
!
!     USAGE:
!           ./moProjection symmMOs.faf MOs.faf [symmMOsIrreps.faf] [pointGroup]
!
!     This is a program that projects MOs in the file MOs.faf onto the MO basis
!     given in the file symmMOs.faf. These two files must include a common
!     atomic orbital basis set and overlap matrix. Additionally, the molecular
!     orientations must be the same.
!
!     It is assumed that the new basis is symmetrized and that there are irrep
!     labels available in symmMOs.faf or symmMOsIrreps.faf. Those two files must
!     have their MOs in the same order, but they do NOT need to be in the same
!     spatial orientation.. Alternatively, one may omit the third command line
!     argument if no others will be used. If there are other optional command
!     line arguments and symmMOsIrreps.faf isn't necessary, use '-' in its
!     place.
!
!     If the [pointGroup] option is given, then the program will print out
!     orbital populations over irreps using the character labels (provided the
!     user supplied point group is one of the ones built in to the program).
!     Currently, the program knows Gaussian's irrep conventions for C2v and D*h.
!
!
!     -H. P. Hratchian, 2024.
!
!
!     USE Connections
!
      use moProjection_mod
!
!     Variable Declarations
!
      implicit none
      integer(kind=int64)::i,j,nCommands,nFafs,nMOSymms
      integer(kind=int64),dimension(:),allocatable::intVecTmp,irrepMOsAlpha,  &
        irrepMOsBeta
      real(kind=real64),dimension(:),allocatable::PAD_weights
      real(kind=real64),dimension(:),allocatable::moColumn,moIrrepPops
      real(kind=real64)::PSP_PAD
      character(len=512)::fafName1,fafName2,fafName3,pointGroup
      character(len=32)::irrepName
      type(mqc_gaussian_unformatted_matrix_file)::faf1,faf2,faf3
      type(mqc_variable)::mqcTmp,mqcTmp1
      type(mqc_variable)::aoOverlap,moCoefficients1alpha,moCoefficients2alpha,  &
        moCoefficients1beta,moCoefficients2beta
!
!     Format Statements
!
 1000 format(1x,'Enter MO Projection Program.')
 1010 format(3x,'Matrix File 1: ',A,/,3x,'Matrix File 2: ',A)
 1011 format(3x,'Matrix File 3: ',A)
 1012 format(3x,'No third file provided. Symmetry labels will be taken from Matrix File 1.')
 1020 format(3x,'Point Group: ',A)
 2000 format(/,1x,'ALPHA Column ',i5)
 2005 format(/,1x,'BETA  Column ',i5)
 2010 format(/,1x,'Pops over irreps')
 2020 format(1x,A7,2x,f14.6)
 2025 format(1x,A7,2x,f14.6,2x)
 2050 format(1x,'PSP_PAD',1x,F14.7,1x)
 8999 format(/,1x,'END OF MO PROJECTION PROGRAM.')
 9000 format(/,1x,'Expected at least 2 command line arguments, but found ',I2,'.')
 9100 format(/,1x,'Confused by the number of command line arguments.')
!
!
!     Start by loading options from the command line.
!
      write(IOut,1000)
      nCommands = command_argument_count()
      if(nCommands.lt.2) then
        write(iOut,9000) nCommands
        goto 999
      endIf
      nFafs = nCommands
      call get_command_argument(1,fafName1)
      call get_command_argument(2,fafName2)
      if(nCommands.ge.3) then
        call get_command_argument(3,fafName3)
        nFafs = 3
        if(fafName3.eq.'-') nFafs = 2
      else
        nFafs = 2
      endIf
      if(nCommands.ge.4) then
        call get_command_argument(4,pointGroup)
        call String_Change_Case(pointGroup,'U')
      else
        pointGroup = 'Unknown'
      endIf
!
!     Echo the input parameters.
!
      if(nFafs.eq.2) then
        write(iOut,1010) TRIM(fafName1),TRIM(fafName2)
        write(iOut,1012)
      elseIf(nFafs.eq.3) then
        write(iOut,1010) TRIM(fafName1),TRIM(fafName2)
        write(iOut,1011) TRIM(fafName3)
      else
        write(iOut,9100)
        goto 999
      endIf
      write(iOut,1020) TRIM(pointGroup)
!
!     Open and load the fafs.
!
      call faf1%load(fafName1)
      call faf2%load(fafName2)
      if(nFafs.eq.3) call faf3%load(fafName3)
!
!     Load the AO overlap matrix from one of the files (we assume the two files
!     have the same AO overlap matrix). Then load the MO coefficients from the
!     two files.
!
      call faf1%getArray('OVERLAP',mqcVarOut=aoOverlap)
      if(DEBUG) call aoOverlap%print(header='Overlap')
      call faf1%getArray('Alpha MO Coefficients',mqcVarOut=moCoefficients1alpha)
      call faf2%getArray('Alpha MO Coefficients',mqcVarOut=moCoefficients2alpha)
      call faf1%getArray('Beta MO Coefficients',mqcVarOut=moCoefficients1beta)
      call faf2%getArray('Beta MO Coefficients',mqcVarOut=moCoefficients2beta)
!
!     Check the orthonormality of the MOs from the two jobs. Only print these
!     outputs if DEBUG is on.
!
      if(DEBUG) then
        call mqc_print(MatMul(Transpose(moCoefficients1alpha),MatMul(aoOverlap,moCoefficients1alpha)),header='C1.S.C1')
        call mqc_print(MatMul(Transpose(moCoefficients2alpha),MatMul(aoOverlap,moCoefficients2alpha)),header='C2.S.C2')
        call mqc_print(MatMul(Transpose(moCoefficients1alpha),MatMul(aoOverlap,moCoefficients2alpha)),header='C1.S.C2')
      endIf
!
!     Grab symmetry array from the first file.
!
      write(*,*) "Andrew here"
      if(nFafs.eq.3) then
        call faf3%getArray('FILE 563 INTEGERS',mqcVarOut=mqcTmp)
        call faf3%getArray('FILE 564 INTEGERS',mqcVarOut=mqcTmp1)
      else
        call faf1%getArray('FILE 563 INTEGERS',mqcVarOut=mqcTmp)
        call faf1%getArray('FILE 564 INTEGERS',mqcVarOut=mqcTmp1)
      endIf
      if(DEBUG) then
        call mqcTmp%print(header='File 563 Integers')
        call mqcTmp1%print(header='File 564 Integers')
      endIf
      irrepMOsAlpha = mqcTmp
      irrepMOsBeta = mqcTmp
      nMOSymms = maxVal(irrepMOsAlpha)
      if(DEBUG) write(iOut,*)' nMOSymms = ',nMOSymms
!
!     Project the second structure's MOs into the basis of the first structure's
!     MOs. Then evaluate each MO's percentage character by irrep.
!

!
!     Obtaining preset PAD_weights here
!
      Allocate(PAD_weights((nMOSymms+1)))
      call pointGroupPADweight(pointGroup,PAD_weights,nMOSymms)
      write(*,*) "PAD_weights size",Size(PAD_weights)
!     First, run the analysis over alpha MOs.
      mqcTmp = MatMul(Transpose(moCoefficients1alpha),MatMul(aoOverlap,moCoefficients2alpha))
      if(Allocated(moIrrepPops)) deAllocate(moIrrepPops)
      Allocate(moIrrepPops(0:nMOSymms))
      do i = 1,Size(mqcTmp,2)
        mqcTmp1 = mqcTmp%column(i)
        write(*,*)' ALPHA Column ',i
        call mqc_print(sum(mqcTmp1),header='  sum of column: ')
        call mqc_print(dot_product(mqcTmp1,mqcTmp1),header='  mag of column: ')
        mqcTmp1 = mqcTmp1*mqcTmp1
        call mqc_print(sum(mqcTmp1),header='  mag again:     ')
        moIrrepPops = float(0)
        do j = 1,Size(mqcTmp1)
          moIrrepPops(irrepMOsAlpha(j)) = moIrrepPops(irrepMOsAlpha(j)) + float(mqcTmp1%getVal([j]))
        endDo
        select case(pointGroup)
        case('D*H','DINFH')
          write(iOut,2010)
          do j = 1,Size(moIrrepPops(1:))
!           write(iOut,2020) TRIM(pointGroupIrrepNameDinfH(j)),moIrrepPops(j)
          endDo
        case('C2V','C02V')
          write(iOut,2010)
          do j = 1,Size(moIrrepPops(1:))
            write(iOut,2020) TRIM(pointGroupIrrepNameC2v(j)),moIrrepPops(j)
          endDo
          write(iOut,2050) Calc_PSP_PAD(moIrrepPops,PAD_weights,Size(moIrrepPops))
        case('D2H','D02H')
          write(iOut,2010)
          do j = 0,Size(moIrrepPops(1:))
            write(iOut,2020) TRIM(pointGroupIrrepNameD2H(j)),moIrrepPops(j)
          endDo
          write(iOut,2050) Calc_PSP_PAD(moIrrepPops,PAD_weights,Size(moIrrepPops))
        case('CS','Cs')
          write(iOut,2010)
          do j = 0,Size(moIrrepPops(1:))
            write(iOut,2020) TRIM(pointGroupIrrepNameCs(j)),moIrrepPops(j)
          endDo
          write(iOut,2050) Calc_PSP_PAD(moIrrepPops,PAD_weights,Size(moIrrepPops))
        case('C2H','C2h')
          write(iOut,2010)
          do j = 0,Size(moIrrepPops(1:))
            write(iOut,2020) TRIM(pointGroupIrrepNameCs(j)),moIrrepPops(j)
          endDo
          write(iOut,2050) Calc_PSP_PAD(moIrrepPops,PAD_weights,Size(moIrrepPops))
        case default
          call mqc_print(moIrrepPops,iOut=iOut,header='Pops over irreps')
        end select
        if(ABS(moIrrepPops(0)).gt.1d-3) write(iOut,*) 'Unassigned Symmetry Population: ',moIrrepPops(0)
        write(*,*)
        write(*,*)
      endDo
!
!     Now, run the analysis over beta MOs.
      mqcTmp = MatMul(Transpose(moCoefficients1beta),MatMul(aoOverlap,moCoefficients2beta))
      if(Allocated(moIrrepPops)) deAllocate(moIrrepPops)
      Allocate(moIrrepPops(0:nMOSymms))
      do i = 1,Size(mqcTmp,2)
        mqcTmp1 = mqcTmp%column(i)
        write(*,*)' BETA  Column ',i
        call mqc_print(sum(mqcTmp1),header='  sum of column: ')
        call mqc_print(dot_product(mqcTmp1,mqcTmp1),header='  mag of column: ')
        mqcTmp1 = mqcTmp1*mqcTmp1
        call mqc_print(sum(mqcTmp1),header='  mag again:     ')
        moIrrepPops = float(0)
        do j = 1,Size(mqcTmp1)
          moIrrepPops(irrepMOsBeta(j)) = moIrrepPops(irrepMOsBeta(j)) + float(mqcTmp1%getVal([j]))
        endDo
        select case(pointGroup)
        case('D*H','DINFH')
          write(iOut,2010)
          do j = 1,Size(moIrrepPops(1:))
!           write(iOut,2020) TRIM(pointGroupIrrepNameDinfH(j)),moIrrepPops(j)
          endDo
        case('C2V','C02V')
          write(iOut,2010)
          do j = 1,Size(moIrrepPops(1:))
            write(iOut,2020) TRIM(pointGroupIrrepNameC2v(j)),moIrrepPops(j)
          endDo
          write(iOut,2050) Calc_PSP_PAD(moIrrepPops,PAD_weights,Size(moIrrepPops))
        case('D2H','D02H')
          write(iOut,2010)
          do j = 0,Size(moIrrepPops(1:))
            write(iOut,2020) TRIM(pointGroupIrrepNameD2H(j)),moIrrepPops(j)
          endDo
          write(iOut,2050) Calc_PSP_PAD(moIrrepPops,PAD_weights,Size(moIrrepPops))
        case('CS','Cs')
          write(iOut,2010)
          do j = 0,Size(moIrrepPops(1:))
            write(iOut,2020) TRIM(pointGroupIrrepNameCs(j)),moIrrepPops(j)
          endDo
          write(iOut,2050) Calc_PSP_PAD(moIrrepPops,PAD_weights,Size(moIrrepPops))
        case default
          call mqc_print(moIrrepPops,iOut=iOut,header='Pops over irreps')
        end select
        if(ABS(moIrrepPops(0)).gt.1d-3) write(iOut,*) 'Unassigned Symmetry Population: ',moIrrepPops(0)
        write(*,*)
        write(*,*)
      end do
!
  999 Continue
      call faf1%closeFile()
      call faf2%closeFile()
      if(nFafs.eq.3) call faf3%closeFile()
      write(iOut,8999)
      end program moProjection
