
!
! This module sets some parameters used in the ManBo program.
! 
!

MODULE manbo_parameters

  IMPLICIT NONE

  DOUBLE PRECISION, PARAMETER :: aua_to_manbo = 0.496147527
  ! aua_to_manbo converts the atomic unit of acceleration (hartree/bohr*amu) into ManBo unit (angstrom/fs^2)
  DOUBLE PRECISION, PARAMETER :: aup_to_manbo = 83.32475916
  ! aup_to_manbo converts the ManBo unit of linear momentum (amu*angstrom/fs) into atomic unit ("reduced plank's constant"/bohr)
  DOUBLE PRECISION, PARAMETER :: aul_to_manbo = 157.4609749
  ! aul_to_manbo converts the ManBo unit of angular momentum (amu*angstrom^2/fs) into atomic unit ("reduced plank's constant")
  DOUBLE PRECISION, PARAMETER :: const_boltzmann = 1.3806503E-23
  ! const_boltzmann is the Boltzmann constant, given in Joule per Kelvin
  DOUBLE PRECISION, PARAMETER :: convert_energy = 1.66053886E-17
  ! convert_energy is used to convert to Joule the energy unit used in ManBo
  DOUBLE PRECISION, PARAMETER :: convert_jou_har = 2.293712569E17
  ! convert_jou_har is used to convert the energy unit Joule to Hartree (atomic unit)
  DOUBLE PRECISION, PARAMETER :: convert_ang_bohr = 1.889725989
  ! convert_ang_bohr is used to convert the distance unit Angstrom to Bohr (atomic unit)
  DOUBLE PRECISION, PARAMETER :: convert_mol_energy = 627.494390
  ! convert_energy is used to convert hartree unit to kcal/mol

END MODULE manbo_parameters
