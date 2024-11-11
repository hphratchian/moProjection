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
      subroutine pointGroupIrrepNameC2v(irrepVal,irrepName,PAD_weights)
!
!     This routine returns the character string for the irrep name
!     corresponding to the irrep integer value according to Gaussian's internal
!     definitions for the C2v point group in the dummy argument irrepName.
!
!     The optional argument PAD_weights will be filled with the numbers of
!     isotropic, perpendicular, and parallel waves expected in a one-electron
!     photo-detachment experiment where the molecule is rotating freely. The
!     numbers of each type are based on using the three primary axis as discrete
!     representations of the rotating system..
!
!
      implicit none
      integer(kind=int64),intent(in)::irrepVal
      character(len=32),intent(out)::irrepName
      integer(kind=int64),dimension(3),optional::PAD_weights
!
!     Andrew October 30th --- PAD Weights below are written as
!     isotropic,perpendicular,and parallel in that order.
      select case(irrepVal)
      case(0)
        irrepName = '?'
        PAD_weights = [ 0,0,0 ]
      case(1)
        irrepName = 'A1'
        PAD_weights = [ 1,0,3 ]
      case(2)
        irrepName = 'A2'
        PAD_weights = [ 0,1,0 ]
      case(3)
        irrepName = 'B1'
        PAD_weights = [ 1,2,0 ]
      case(4)
        irrepName = 'B2'
        PAD_weights = [ 1,2,0 ]
      case default
        call mqc_error('Unknown C2v Irrep Value.')
      end select
!
      return
      end subroutine pointGroupIrrepNameC2v

!PROCEDURE pointGroupIrrepNameD2H(irrepVal)
      subroutine pointGroupIrrepNameD2H(irrepVal,irrepName,PAD_weights)
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
      integer(kind=int64),intent(in)::irrepVal
      character(len=32),intent(out)::irrepName
      integer(kind=int64),dimension(3),optional::PAD_weights
!
!     Andrew October 30th --- PAD Weights below are written as
!     isotropic,perpendicular,and parallel in that order.
      select case(irrepVal)
      case(0)
        irrepName = '?'
        PAD_weights = [ 0,0,0 ]
      case(1)
        irrepName = 'AG'
        PAD_weights = [ 0,0,3 ]
      case(2)
        irrepName = 'AU'
        PAD_weights = [ 0,2,0 ]
      case(3)
        irrepName = 'B1G'
        PAD_weights = [ 0,2,0 ]
      case(4)
        irrepName = 'B1U'
        PAD_weights = [ 0,0,0 ]
      case(5)
        irrepName = 'B2G'
        PAD_weights = [ 0,2,0 ]
      case(6)
        irrepName = 'B2U'
        PAD_weights = [ 0,0,0 ]
      case(7)
        irrepName = 'B3G'
        PAD_weights = [ 0,2,0 ]
      case(8)
        irrepName = 'B3U'
        PAD_weights = [ 0,0,0 ]
      case default
        call mqc_error('Unknown D2H Irrep Value.')
      end select
!
      return
      end subroutine pointGroupIrrepNameD2H
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
      integer(kind=int64),intent(in),dimension(3*totirreps)::PAD_tot_weights
      real(kind=real64)::PSP_PAD,temp_pad
      integer(kind=int64)::i,j,temp,count_temp
    
      PSP_PAD = 0.0
      write(400,*) "This is the story of a girl....", totirreps
      write(400,*) "Who traveled the whole world...", SIZE(moirrepPops)
      write(400,*) "I forgot the other lyrics......", SIZE(PAD_tot_weights)
      do i = 0,totirreps
        temp = 0
        count_temp = 0
        temp_pad = 0.0
        do j = 1,3 
          if(j.eq.3) then
            temp = temp + 2*PAD_tot_weights((i*3)+j)
          end if
          count_temp = count_temp + PAD_tot_weights((i*3)+j)
        end do
        if(count_temp/=0) then
          temp_pad = float(temp)/float(count_temp)
        else
          temp_pad = 0.0
        end if
        PSP_PAD = PSP_PAD + moIrrepPops(i+1)*temp_pad
      end do

      end function Calc_PSP_PAD

      end module moProjection_mod
