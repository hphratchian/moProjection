      module moProjection_mod
!
!     This module file is meant to support the program moProjection.
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
      logical,parameter::DEBUG=.false.
!
!
      CONTAINS
!
!
!PROCEDURE pointGroupIrrepNameC2v(irrepVal)
      function pointGroupIrrepNameC2v(irrepVal) result(irrepName)
!
!     This function returns the character string for the irrep name
!     corresponding to the irrep integer value according to Gaussian's internal
!     definitions for the C2v point group.
!
!
      implicit none
      integer(kind=int64),intent(in)::irrepVal
      character(len=32)::irrepName
!
      select case(irrepVal)
      case(0)
        irrepName = '?'
      case(1)
        irrepName = 'A1'
      case(2)
        irrepName = 'A2'
      case(3)
        irrepName = 'B1'
      case(4)
        irrepName = 'B2'
      case default
        call mqc_error('Unknown C2v Irrep Value.')
      end select
!
      return
      end function pointGroupIrrepNameC2v
!
!
!PROCEDURE pointGroupIrrepNameDinfH(irrepVal)
      function pointGroupIrrepNameDinfH(irrepVal) result(irrepName)
!
!     This function returns the character string for the irrep name
!     corresponding to the irrep integer value according to Gaussian's internal
!     definitions for the D_infinity_h point group.
!
!
      implicit none
      integer(kind=int64),intent(in)::irrepVal
      character(len=32)::irrepName
!
      select case(irrepVal)
      case(0)
        irrepName = '?'
      case(1)
        irrepName = 'AG'
      case(2)
        irrepName = 'AU'
      case(3)
        irrepName = 'B1G'
      case(4)
        irrepName = 'B1U'
      case(5)
        irrepName = 'B2G'
      case(6)
        irrepName = 'B2U'
      case(7)
        irrepName = 'B3G'
      case(8)
        irrepName = 'B3U'
      case default
        call mqc_error('Unknown DinfH Irrep Value.')
      end select
!
      return
      end function pointGroupIrrepNameDinfH


      end module moProjection_mod
