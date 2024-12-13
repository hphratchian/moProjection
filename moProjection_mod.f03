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
!PROCEDURE pointGroupPADweight(PAD_weights)
      subroutine pointGroupPADweight(pointGroup,PAD_weight,nMOSymms)
!
!     This routine returns the character string for the irrep name
!     corresponding to the irrep integer value according to Gaussian's internal
!     definitions for the D2H point group in the dummy argument irrepName.
!
!     The optional argument PAD_weights will be filled with the numbers of
!     isotropic, perpendicular, and parallel waves expected in a one-electron
!     photo-detachment experiment where the molecule is rotating freely. The
!     numbers of each type are based on using the three primary axis as discrete
!     representations of the rotating system..
!
!
      implicit none
      integer(kind=int64),intent(in)::nMOSymms
      character(len=512),intent(in)::pointGroup
      real(kind=real64),dimension(nMOSymms+1)::PAD_weight
!
!     Andrew October 30th --- PAD Weights below are written as
!     isotropic,perpendicular,and parallel in that order.
!     Andrew November 30th --- PAD Weights below are written NOT
!     isotropic,perpendicular,and parallel in that order.. BUT
!     The total number of parallel waves in x,y,z or z direction.
!     Averaging should be over all three directions NOT the number of total
!     waves. The old way of doing it doesn't make any sense.
      select case(pointGroup)
      case('D2H')
        PAD_weight = [0.0,2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0] 
      case('C2V')
        PAD_weight = [0.0,1.666,0.0,0.0,0.0]
      case default
        call mqc_error('Unknown pointGroup.')
      end select
!
      return
      end subroutine pointGroupPADweight

!PROCEDURE pointGroupIrrepNameC2v(irrepVal)
      function pointGroupIrrepNameC2v(irrepVal) result(irrepname)
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
!PROCEDURE pointGroupIrrepNameD2h(irrepVal)
      function pointGroupIrrepNameD2h(irrepVal) result(irrepName)
!
!     This function returns the character string for the irrep name
!     corresponding to the irrep integer value according to Gaussian's internal
!     definitions for the D2h point group.
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
        call mqc_error('Unknown D2H Irrep Value.')
      end select
!
      return
      end function pointGroupIrrepNameD2h

      function Calc_PSP_PAD(moIrrepPops,PAD_tot_weights,totirreps) result(PSP_PAD)
!
!     This function calculates the PSP PAD value. One takes all the moIrrepPops
!     for the orbital in question, along with the total PAD weights stored in
!     one array, and then multiples the each moirrep pop against the average of
!     the pad_weights for that irrep, taking into account the number of
!     parallel,perpendicular,or isosymmetric waves for each irrep that our
!     possible.
!

      implicit none
      integer(kind=int64),intent(in)::totirreps
      real(kind=real64),intent(in),dimension(totirreps)::moIrrepPops
      real(kind=real64),intent(in),dimension(totirreps+1)::PAD_tot_weights
      real(kind=real64)::PSP_PAD
      integer(kind=int64)::i
    
      PSP_PAD = 0.0
      do i = 1,totirreps
        PSP_PAD = PSP_PAD + moIrrepPops(i+1)*PAD_tot_weights(i+1)
      end do

      end function Calc_PSP_PAD

      end module moProjection_mod
