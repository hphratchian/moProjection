      program moProjection
!
!     This is a test program.
!
!
!     -H. P. Hratchian, 2024.
!
!
!     USE Connections
!
      use mqc_general
      use mqc_molecule
      use mqc_gaussian
      use mqc_algebra2
      use mqc_algebra
      use iso_fortran_env
!
!     Variable Declarations
!
      implicit none
      integer(kind=int64),parameter::iOut=6
      integer(kind=int64)::i,j,nCommands,nMOSymms
      integer(kind=int64),dimension(:),allocatable::intVecTmp,irrepMOsAlpha,  &
        irrepMOsBeta
      real(kind=real64),dimension(:),allocatable::moColumn,moIrrepPops
      character(len=512)::fafName1,fafName2,fafName3
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
 8999 format(/,1x,'END OF MO PROJECTION PROGRAM.')
 9000 format(/,1x,'Expected 2 or 3 command line arguments, but found ',I2,'.')
 9100 format(/,1x,'Confused by the number of command line arguments.')
!
!
!     Start the job by loading the FAFs.
!
      write(IOut,1000)
      nCommands = command_argument_count()
      if(nCommands.lt.2.or.nCommands.gt.3) then
        write(iOut,9000) nCommands
        goto 999
      endIf
      call get_command_argument(1,fafName1)
      call get_command_argument(2,fafName2)
      if(nCommands.eq.2) then
        write(iOut,1010) TRIM(fafName1),TRIM(fafName2)
        write(iOut,1012)
      elseIf(nCommands.eq.3) then
        call get_command_argument(3,fafName3)
        write(iOut,1010) TRIM(fafName1),TRIM(fafName2)
        write(iOut,1011) TRIM(fafName3)
      else
        write(iOut,9100)
        goto 999
      endIf
      call faf1%load(fafName1)
      call faf2%load(fafName2)
      if(nCommands.eq.3) call faf3%load(fafName3)
!
!     Load the AO overlap matrix from one of the files (we assume the two files
!     have the same AO overlap matrix). Then load the MO coefficients from the
!     two files.
!
      call faf1%getArray('OVERLAP',mqcVarOut=aoOverlap)
      call aoOverlap%print(header='Overlap')
      call faf1%getArray('Alpha MO Coefficients',mqcVarOut=moCoefficients1alpha)
      call faf2%getArray('Alpha MO Coefficients',mqcVarOut=moCoefficients2alpha)
      call faf1%getArray('Beta MO Coefficients',mqcVarOut=moCoefficients1beta)
      call faf2%getArray('Beta MO Coefficients',mqcVarOut=moCoefficients2beta)
!
!     Check the orthonormality of the MOs from the two jobs.
!
      call mqc_print(MatMul(Transpose(moCoefficients1alpha),MatMul(aoOverlap,moCoefficients1alpha)),header='C1.S.C1')
      call mqc_print(MatMul(Transpose(moCoefficients2alpha),MatMul(aoOverlap,moCoefficients2alpha)),header='C2.S.C2')
      call mqc_print(MatMul(Transpose(moCoefficients1alpha),MatMul(aoOverlap,moCoefficients2alpha)),header='C1.S.C2')
!
!     Grab symmetry array from the first file.
!
      if(nCommands.eq.3) then
        call faf3%getArray('FILE 563 INTEGERS',mqcVarOut=mqcTmp)
        call faf3%getArray('FILE 564 INTEGERS',mqcVarOut=mqcTmp1)
      else
        call faf1%getArray('FILE 563 INTEGERS',mqcVarOut=mqcTmp)
        call faf3%getArray('FILE 564 INTEGERS',mqcVarOut=mqcTmp1)
      endIf
      call mqcTmp%print(header='File 563 Integers')
      call mqcTmp1%print(header='File 564 Integers')
      irrepMOsAlpha = mqcTmp
      irrepMOsBeta = mqcTmp
      nMOSymms = maxVal(irrepMOsAlpha)
      write(iOut,*)' nMOSymms = ',nMOSymms
!
!     Project the second structure's MOs into the basis of the first structure's
!     MOs. Then evaluate each MO's percentage character by irrep.
!
!     First, run the analysis over alpha MOs.
      mqcTmp = MatMul(Transpose(moCoefficients1alpha),MatMul(aoOverlap,moCoefficients2alpha))
      if(Allocated(moIrrepPops)) deAllocate(moIrrepPops)
      Allocate(moIrrepPops(nMOSymms))
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
        call mqc_print(moIrrepPops,iOut=iOut,header='Pops over irreps')
        write(*,*)
        write(*,*)
      endDo
!
!     Now, run the analysis over beta MOs.
      mqcTmp = MatMul(Transpose(moCoefficients1beta),MatMul(aoOverlap,moCoefficients2beta))
      if(Allocated(moIrrepPops)) deAllocate(moIrrepPops)
      Allocate(moIrrepPops(nMOSymms))
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
        call mqc_print(moIrrepPops,iOut=iOut,header='Pops over irreps')
        write(*,*)
        write(*,*)
      endDo

!!
!!     Get the irrep names.
!!
!      write(*,*)
!      write(*,*)' Hrant - Trying to get the point group irrep names.'
!      call MQC_Gaussian_SetDEBUG(.true.)
!      call faf%getArray('POINT GROUP IRREP NAMES',mqcVarOut=mqcTmp)
!hph-

!
  999 Continue
      call faf1%closeFile()
      call faf2%closeFile()
      if(nCommands.eq.3) call faf3%closeFile()
      write(iOut,8999)
      end program moProjection
