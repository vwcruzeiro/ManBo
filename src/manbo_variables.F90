MODULE manbo_variables

  USE manbo_periodic_table
  ! Here we call manbo_periodic_table module
  USE manbo_parameters
  ! Here we call manbo_parameters module

  IMPLICIT NONE

  TYPE atom_in
  ! atom_in contains the information about an atom inside the box
    DOUBLE PRECISION, DIMENSION(3) :: r, v, f
    ! r contains the position coordinates of a given atom, given in angstrom (1E-10 m)
    ! v contains the velocity coordinates of a given atom, given in angstrom/femtosecond (1E5 m/s)
    ! f contains the coordinates of the force acting on a given atom, given in  hartree/bohr (atomic unit)
    DOUBLE PRECISION :: mass, q
    ! mass is the mass of a given atom, given in atomic mass unit (1.66053886E-27 kg)
    ! q is the atomic charge to be used on the eletrostatic embedding
    INTEGER :: an, mol, mol_orig
    ! an is the atomic number of the atom
    ! mol      is the current  molecule numbering to which the atom belongs
    ! mol_orig is the original molecule numbering to which the atom belongs
  END TYPE atom_in
  TYPE atom_out
  ! atom_out contains the information about a replicated atom outside the box
    DOUBLE PRECISION, DIMENSION(3) :: r
    ! r is the vector that contains the position coordinates of a given atom outside the box, given in angstroms (1E-10 m)
    DOUBLE PRECISION :: mass
    ! mass is the mass of a given atom outside the box, and is given in atomic mass unit (1.66053886E-27 kg)
    INTEGER :: an, mol, mol_orig, atom
    ! an is the atomic number of a given atom outside the box
    ! mol      is the current  molecule numbering to which the atom belongs
    ! mol_orig is the original molecule numbering to which the atom belongs
    ! atom is the number of the atom inside the box in which a given atom outside the box has been replicated from
  END TYPE atom_out
  TYPE center_of_mass
    DOUBLE PRECISION, DIMENSION(3) :: r
    ! r is the vector that contains the position coordinates of the mass center of a given molecule, given in angstroms (1E-10 m)
    DOUBLE PRECISION :: mass
    ! mass is the mass of the atoms of a given molecule, and is given in atomic mass unit (1.66053886E-27 kg)
  END TYPE center_of_mass
  TYPE charge_multiplicity
    INTEGER :: q_mol, mul
    ! q_mol contains the charge of a given molecule, and is given in atomic units
    ! mul contains the multiplicity of a given molecule
  END TYPE charge_multiplicity
  TYPE forces
    DOUBLE PRECISION, DIMENSION(3) :: f
    ! f contains force coordinates, and is given in  hartree/bohr (atomic unit)
  END TYPE forces
  TYPE molecules
    INTEGER, DIMENSION(3) :: m
    ! m(1), m(2) and m(3) are molecule numbers, to be used when the algorithm to group molecules is activated
  END TYPE molecules

  TYPE(atom_in), DIMENSION(:), ALLOCATABLE, SAVE :: data_manbo
  ! data_manbo contains all information about the atoms that are inside the box
  TYPE(atom_out), DIMENSION(:), ALLOCATABLE, SAVE :: data_manbo_out
  ! data_manbo_out contains all information about the replicated atoms that are outside the box
  TYPE(center_of_mass), DIMENSION(:), ALLOCATABLE, SAVE :: mc, mc_out, mc_orig, mc_out_orig
  ! mc contains the mass centers of the molecules that are inside the box
  ! mc_out contains the mass centers of the molecules that are outside the box
  ! mc_orig and mc_out_orig are the equivalent of mc and mc_out of original molecules when the molecules were grouped
  TYPE(charge_multiplicity), DIMENSION(:), ALLOCATABLE, SAVE :: char_mul, char_mul_orig
  !  char_mul      contains the current  charge and the multiplicity of all monomers inside the box
  !  char_mul_orig contains the original charge and the multiplicity of all monomers inside the box
  TYPE(molecules), DIMENSION(:), ALLOCATABLE, SAVE :: orig_mols, orig_mols_out
  ! orig_mols is the vector which contains the list of original molecules when the molecules were grouped
  INTEGER :: n_atoms, n_mols, n_mols_orig, n_atoms_out, n_mols_out, n_mols_out_orig, n_qm_procs, qm_prog_memory
  ! n_atoms is the number of atoms which are inside of the box
  ! n_mols      is the current  number of monomers which are inside of the box
  ! n_mols_orig is the original number of monomers which are inside of the box
  ! n_atoms_out is the number of atoms which are outside of the box
  ! n_mols_out      is the current  number of molecules which are outside of the box
  ! n_mols_out_orig is the original number of molecules which are outside of the box
  ! n_qm_procs is the number of processors to be used on the calculations of the quantum chemistry program
  ! qm_prog_method contains the quantum chemistry method to be used by the quantum chemistry program
  DOUBLE PRECISION, DIMENSION(3) :: box_dim
  ! box_dim is the vector which contains the dimensions of the box, and are given in angstroms (1E-10 m)
  DOUBLE PRECISION :: E1, E2, E3
  ! E1, E2 and E3 contains the one-, two- and three-body terms from Many Body Expansion (MBE), given in  hartree (atomic unit)
  DOUBLE PRECISION :: cutoff, delta_t, temperature, target_therm_temp, t_therm, emb_radius, rand_vel_temp
  ! cutoff is the cutoff radius used in calculating the forces, and are given in angstroms (1E-10 m)
  ! delta_t is the intervals of time between the steps, given in femtoseconds (1E-15 s)
  ! temperature is the system's temperature, given in Kelvin
  ! target_therm_temp is the target temperature in case of application of the Berendsen thermostat, given in Kelvin
  ! t_therm is a parameter used in the Berendsen thermostat, given in femtoseconds (1E-15 s)
  ! emb_radius is the cutoff radius of embedding replication, and are given in angstroms (1E-10 m)
  ! rand_vel_temp is the temperature to used to generate random initial velocities, given in Kelvin
  INTEGER :: n_steps, i_print_prop, i_save, i_rst, md_step, expan_order, rand_vel_seed, init_dt, group_monomers
  ! n_steps is the number of steps to be executed by ManBo
  ! i_print_prop is the intervals of steps which the proprieties of the system will be saved
  ! i_save is the intervals of steps which the positions, velocities and other data of the system will be saved
  ! i_rst is the intervals of steps which the restart file will be saved
  ! md_step is the current step of the simulation
  ! expan_order constains the order in which we trucate the MBE (it can be 1, 2 or 3)
  ! rand_vel_seed contains the seed to generate random numbers
  ! init_dt contains the time in which ManBo has started, in miliseconds, encoded as a integer
  ! group_monomers is the number of molecules to be grouped as a single monomer (it can be 1, 2 or 3)
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: orig_repli_mol
  ! orig_repli_mol contains the number of the original molecule associated to a replicated molecule
  LOGICAL :: apply_pbc, apply_therm, use_embedding, use_emb_radius, use_rand_vel, mbe_corr_only, bs_extrapolation
  ! apply_pbc determines if we are going to applicate Periodic Boundary Conditions
  ! apply_therm determines if we are going to applicate the Berendsen thermostat
  ! use_embedding determines if we are going to use eletrostatic embedding on the calculations
  ! use_emb_radius determines if we are going to use an eletrostatic embedding replication cutoff radius
  ! use_rand_vel determines if we are going to use random initial velocities
  ! mbe_corr_only determines if we are going to use the Many-Body Expansion on the correlation energy only
  ! bs_extrapolation determines if we are going to perform extrapolation to infinity basis set
  CHARACTER(LEN=70) :: name_in, name_out, qm_prog, qm_prog_method, qm_prog_basis
  ! name_in is the name of the input file
  ! name_out is the name to be included in the begining of the output files
  ! qm_prog contains the comand to call the quantum cremistry program that will be used for calculations
  ! qm_prog_memory contains the memory to be used by the quantum cremistry program
  ! qm_prog_basis contains the basis set functions to be used by the quantum cremistry program

END MODULE manbo_variables
