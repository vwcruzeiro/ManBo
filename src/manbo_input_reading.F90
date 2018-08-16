MODULE manbo_input_reading

#ifdef USE_PARALLEL
  USE omp_lib
  USE mpi
#elif USE_OPENMPONLY
  USE omp_lib
#endif

  USE manbo_variables
  ! Here we call manbo_variables module
  USE manbo_subroutines
  ! Here we call manbo_subroutines module
  USE manbo_log
  ! Here we call manbo_log module

  IMPLICIT NONE

  CONTAINS

#ifdef USE_PARALLEL
  SUBROUTINE read_input(mpi_id,mpi_size)
#else
  SUBROUTINE read_input(mpi_id)
#endif
  ! This subroutine reads ManBo input file and its variables

    IMPLICIT NONE

#ifdef USE_PARALLEL
    INTEGER, INTENT(in) :: mpi_id,mpi_size
    INTEGER :: ierr
#else
    INTEGER, INTENT(in) :: mpi_id
#endif
    CHARACTER(LEN=200) :: line, text
    INTEGER :: i, j
    DOUBLE PRECISION :: L, kr, val
    INTEGER :: allocate_status
    
    IF (mpi_id == 0) THEN
      delta_t = 1.0
      n_steps = 1
      i_print_prop = 1
      i_save = 1
      i_rst = 1
      n_qm_procs = 1
      use_rand_vel = .FALSE.
      mbe_corr_only = .FALSE.
      bs_extrapolation = .FALSE.
      rand_vel_seed = 10
      apply_pbc = .TRUE.
      expan_order = 2
      group_monomers = 1
      use_embedding = .FALSE.
      qm_prog = "g03"
      qm_prog_memory = 0
      apply_therm = .FALSE.
      t_therm = 10*delta_t
      ! These are default values of ManBo
      
      n_mols = 0
      box_dim = (/ 0.0, 0.0, 0.0 /)
      cutoff = 0.0
      
      ! Getting the input file name
      CALL getarg(1,name_in)
      ! Defining name_out
      name_out = name_in
      CALL delsubstr(name_out,".in ")
      
      OPEN(UNIT=13,FILE=name_in,STATUS='OLD')
        line = ""
        DO WHILE (line/="END_OF_INPUT_VARIABLES")
          READ(13, '(a200)') text
          READ(text, *) line
          IF(line=="box_dim=") THEN
            READ(text, *) line, box_dim(1), box_dim(2), box_dim(3)
          END IF
          IF(line=="cutoff=") THEN
            READ(text, *) line, cutoff
          END IF
          IF(line=="n_mols=") THEN
            READ(text, *) line, n_mols
          END IF
          IF(line=="delta_t=") THEN
            READ(text, *) line, delta_t
          END IF
          IF(line=="n_steps=") THEN
            READ(text, *) line, n_steps
          END IF
          IF(line=="i_print_prop=") THEN
            READ(text, *) line, i_print_prop
          END IF
          IF(line=="i_save=") THEN
            READ(text, *) line, i_save
          END IF
          IF(line=="i_rst=") THEN
            READ(text, *) line, i_rst
          END IF
          IF(line=="n_qm_procs=") THEN
            READ(text, *) line, n_qm_procs
          END IF
          IF(line=="use_rand_vel=") THEN
            READ(text, *) line, use_rand_vel
          END IF
          IF(line=="rand_vel_temp=") THEN
            READ(text, *) line, rand_vel_temp
          END IF
          IF(line=="rand_vel_seed=") THEN
            READ(text, *) line, rand_vel_seed
          END IF
          IF(line=="apply_pbc=") THEN
            READ(text, *) line, apply_pbc
          END IF
          IF(line=="expan_order=") THEN
            READ(text, *) line, expan_order
          END IF
          IF(line=="group_monomers=") THEN
            READ(text, *) line, group_monomers
          END IF
          IF(line=="use_embedding=") THEN
            READ(text, *) line, use_embedding
          END IF
          IF(line=="use_emb_radius=") THEN
            READ(text, *) line, use_emb_radius
          END IF
          IF(line=="emb_radius=") THEN
            READ(text, *) line, emb_radius
          END IF
          IF(line=="qm_prog=") THEN
            READ(text, *) line, qm_prog
          END IF
          IF(line=="qm_prog_memory=") THEN
            READ(text, *) line, qm_prog_memory
          END IF
          IF(line=="qm_prog_method=") THEN
            READ(text, *) line, qm_prog_method
          END IF
          IF(line=="qm_prog_basis=") THEN
            READ(text, *) line, qm_prog_basis
          END IF
          IF(line=="bs_extrapolation=") THEN
            READ(text, *) line, bs_extrapolation
          END IF
          IF(line=="apply_therm=") THEN
            READ(text, *) line, apply_therm
          END IF
          IF(line=="target_therm_temp=") THEN
            READ(text, *) line, target_therm_temp
          END IF
          IF(line=="t_therm=") THEN
            READ(text, *) line, t_therm
          END IF
          IF(line=="mbe_corr_only=") THEN
            READ(text, *) line, mbe_corr_only
          END IF
        END DO
      ! Here we read the initial variables of the input file
      
      n_mols_orig = n_mols
      
      CALL log_start
      ! Here we create the log file of ManBo
      
      IF(n_mols==0 .OR. ABS(cutoff)<1E-5 .OR. (apply_therm .AND. ABS(target_therm_temp)<1E-5)) THEN
        PRINT *, "Compulsory parameters are missing or are equal to zero in the input file. ManBo cannot run."
        CALL log_write("ERROR: Compulsory parameters are missing or are equal to zero in the input file.")
        CALL log_close(1)
        STOP
      END IF
      
      CALL log_write("")
      CALL log_write("Quantum Chemistry Program (QCP) to be used: " // TRIM(ADJUSTL(qm_prog)))
        IF(qm_prog_memory > 0) THEN
          WRITE(line,*) qm_prog_memory
          CALL log_write("Memory required by the QCP: " // TRIM(ADJUSTL(line)) // "MB")
        END IF
      CALL log_write("Quantum Chemistry Method to be used: " // TRIM(ADJUSTL(qm_prog_method)))
        IF (bs_extrapolation) THEN
    CALL log_write("Extrapolation to infinity basis set will be performed (thus cc-pVDZ and cc-pVTZ calculations will be done)")
        ELSE
          CALL log_write("Basis Set Functions to be used: " // TRIM(ADJUSTL(qm_prog_basis)))
        END IF
        IF(bs_extrapolation .AND. TRIM(ADJUSTL(qm_prog_method))/="MP2") THEN
          PRINT *, "Extrapolation to infinity basis set can only be performed with MP2 method."
          CALL log_write("ERROR: Extrapolation to infinity basis set can only be performed with MP2 method,")
          CALL log_write("       on manbo_input_reading.F90")
          CALL log_close(1)
          STOP
        END IF
        IF (mbe_corr_only) THEN
          CALL log_write("Apply the Many-Body Expansion on the correlation energy only: Yes")
        ELSE
          CALL log_write("Apply the Many-Body Expansion on the correlation energy only: No")
        END IF
      WRITE(line,'(f8.4)') cutoff
      CALL log_write("Cutoff Radius to be applied on the Many Body Expansion (MBE) calculations: " // TRIM(ADJUSTL(line)) // " A")
      WRITE(line,*) n_mols
      CALL log_write("Number of molecules of the system: " // TRIM(ADJUSTL(line)))
      
      READ(13, *) line
      READ(13, *) n_atoms
      IF (vector_module(box_dim)<1E-5) THEN
        L = MINVAL(box_dim)
      ELSE
        L = ((FLOAT(n_mols)/20.0)**(1.0/3.0))*8.0
      END IF
      kr = 2.0*(8*((cutoff/L)**3) + 12*((cutoff/L)**2) + 6*(cutoff/L)) + 1.0/n_mols
      ! This is an estimate for allocate data_manbo_out
    END IF
    
#ifdef USE_PARALLEL
    CALL MPI_BCAST(kr, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST(apply_pbc, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST(n_mols,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST(n_atoms, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST(n_steps, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST(name_out, 70, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
#endif /* USE_PARALLEL */
    
    IF(kr<0.1 .OR. .NOT. apply_pbc) THEN
      ALLOCATE(data_manbo(n_atoms), mc(n_mols), char_mul(n_mols), &
               char_mul_orig(n_mols), stat=allocate_status)
    ELSE
      ALLOCATE(data_manbo(n_atoms), data_manbo_out(NINT(kr*n_atoms)),&
               mc(n_mols), char_mul(n_mols), char_mul_orig(n_mols),&
               orig_repli_mol(NINT(kr*n_mols)), stat=allocate_status)
    END IF
    ! Here we allocate some variables
    IF(allocate_status/=0) THEN
      PRINT *, "No memory enough to allocate. ManBo cannot run."
      CALL log_write("ERROR: Error on allocate data_manbo, data_manbo_out or mc")
      CALL log_write("       on manbo_input_reading.F90")
      CALL log_close(1)
      STOP
    END IF
    
    IF (mpi_id == 0) THEN
        DO i=1,n_atoms
          READ(13,*) line, data_manbo(i)%r(1), data_manbo(i)%r(2), data_manbo(i)%r(3), data_manbo(i)%mol_orig
          data_manbo(i)%mol = data_manbo(i)%mol_orig
            DO j=1,SIZE(atom_symbol, DIM=1)
              IF(atom_symbol(j)==line) THEN
                data_manbo(i)%mass = atom_mass(j)
                data_manbo(i)%an = j
                EXIT
              END IF
            END DO
        END DO
        READ(13,*) line
        DO i=1,n_mols
          READ(13,*) j, char_mul_orig(j)%q_mol, char_mul_orig(j)%mul
        END DO
      char_mul = char_mul_orig  
    END IF
    
    IF (mpi_id == 0) THEN
        IF(use_embedding) THEN
          READ(13,*) line
          DO i=1,n_atoms
            READ(13,*) data_manbo(i)%q
          END DO
        END IF
      ! Here we read some variables of data_manbo
      
      WRITE(line,*) n_atoms
      CALL log_write("Number of atoms of the system: " // TRIM(ADJUSTL(line)))
      WRITE(line,*) n_steps
      CALL log_write("Number of steps to be executed: " // TRIM(ADJUSTL(line)))
      WRITE(line,'(f8.4)') delta_t
      CALL log_write("Time interval between the steps: " // TRIM(ADJUSTL(line)) // " femtoseconds")
        IF (use_rand_vel) THEN
          CALL log_write("Generate random initial velocities: Yes")
          WRITE(line,*) rand_vel_temp
          CALL log_write("Temperature to be used to generate random initial velocities: " // TRIM(ADJUSTL(line)) // " K")
          WRITE(line,*) rand_vel_seed
          CALL log_write("Seed to be used to generate the random numbers, for random initial velocities: " // TRIM(ADJUSTL(line)))
        ELSE
          CALL log_write("Generate random initial velocities: No")
        END IF
        IF (apply_pbc) THEN
          WRITE(line,*) "Yes"
        ELSE
          WRITE(line,*) "No"
        END IF
      CALL log_write("Apply Periodic Boundary Conditions (PBC): " // TRIM(ADJUSTL(line)))
      WRITE(line,*) expan_order
      CALL log_write("Order in which the MBE is going to be truncated: " // TRIM(ADJUSTL(line)))
        IF (expan_order>3) THEN
          PRINT *, "The expan_order value cannot be higher than 3. ManBo cannot run."
          CALL log_write("ERROR: expan_order is higher than 3")
          CALL log_close(1)
          STOP
        END IF
      IF (group_monomers > 1) THEN
        WRITE(line,*) group_monomers
        CALL log_write("Number of molecules to be grouped as a single monomer: " // TRIM(ADJUSTL(line)))
          IF (group_monomers>3) THEN
            PRINT *, "The group_monomers value cannot be higher than 3. ManBo cannot run."
            CALL log_write("ERROR: group_monomers is higher than 3")
            CALL log_close(1)
            STOP
          END IF
      END IF
        IF (use_embedding) THEN
          WRITE(line,*) "Yes"
        ELSE
          WRITE(line,*) "No"
        END IF
      CALL log_write("Use electrostatic embedding: " // TRIM(ADJUSTL(line)))
        IF (use_embedding) THEN
          IF (.NOT. apply_pbc) THEN
            IF (use_emb_radius) THEN
              CALL log_write("Use an Embedding Cutoff Radius: No, because PBC are not being applied ")
            ELSE
              CALL log_write("Use an Embedding Cutoff Radius: No")
            END IF
          ELSE
            IF (use_emb_radius) THEN
              WRITE(line,*) "Yes"
            ELSE
              WRITE(line,*) "No"
            END IF
            CALL log_write("Use an Embedding Cutoff Radius: " // TRIM(ADJUSTL(line)))
              IF (use_emb_radius) THEN
                WRITE(line,*) emb_radius
                CALL log_write("Embedding Cutoff Radius to be used: " // TRIM(ADJUSTL(line)) // " A")
              END IF
          END IF
        END IF
        IF (apply_therm) THEN
          CALL log_write("Apply Berendsen thermostat: Yes")
          WRITE(line,*) target_therm_temp
          CALL log_write("Target temperature of the Berendsen thermostat: " // TRIM(ADJUSTL(line)) // " K")
          WRITE(line,*) t_therm
          CALL log_write("Factor T used in Berendsen thermostat: " // TRIM(ADJUSTL(line)) // " picoseconds")
        ELSE
          CALL log_write("Apply Berendsen thermostat: No")
        END IF
      CALL log_write("")
#ifdef USE_PARALLEL
      CALL log_write("|*** Running the Hybrid MPI-OpenMP version of ManBo ***")
      WRITE(line,*) mpi_size
      CALL log_write("|  Number of MPI threads available:       " // TRIM(ADJUSTL(line)))
      WRITE(line,*) omp_get_num_procs()
      CALL log_write("|  Number of OpenMP processors available: " // TRIM(ADJUSTL(line)))
      WRITE(line,*) omp_get_max_threads()
      CALL log_write("|  Number of OpenMP threads available:    " // TRIM(ADJUSTL(line)))
#elif USE_OPENMPONLY
      CALL log_write("|*** Running the OpenMP version of ManBo ***")
      WRITE(line,*) omp_get_num_procs()
      CALL log_write("|  Number of OpenMP processors available: " // TRIM(ADJUSTL(line)))
      WRITE(line,*) omp_get_max_threads()
      CALL log_write("|  Number of OpenMP threads available:    " // TRIM(ADJUSTL(line)))
#else
      CALL log_write("|*** Running the Serial version of ManBo ***")
#endif /* USE_PARALLEL */
      WRITE(line,*) n_qm_procs
      CALL log_write("|  Number of processors to be used by the QM program on each calculation: " // TRIM(ADJUSTL(line)))
      CALL log_write("")
      
      IF (.NOT. use_rand_vel) THEN
        READ(13, *) line
          DO i=1,n_atoms
            READ(13, *) data_manbo(i)%v(1), data_manbo(i)%v(2), data_manbo(i)%v(3)
          END DO
      ELSE
        CALL generate_rand_vel()
      END IF
      ! Here we set the velocities of data_manbo.
      
      IF(use_embedding) THEN
        CALL log_write("Input orientation values entered for the atoms (Coordinates are given in Angstroms and charges in A.U.):")
        CALL log_write("===============================================================================================")
        CALL log_write("  Atom   Symbol     Position X       Position Y       Position Z     Molecule    Emb. Charge   ")
        CALL log_write("===============================================================================================")
          DO i=1,n_atoms
            WRITE(line,'(1x,i4,6x,a2,4x,3(f14.10,3x),2x,i5,6x,f10.6)') i, atom_symbol(data_manbo(i)%an), data_manbo(i)%r(1),&
                                                                        data_manbo(i)%r(2), data_manbo(i)%r(3),&
                                                                        data_manbo(i)%mol_orig,data_manbo(i)%q
            CALL log_write(TRIM(line))
          END DO
        CALL log_write("-----------------------------------------------------------------------------------------------")
      ELSE
        CALL log_write("Input orientation values entered for the atoms (Coordinates are given in Angstroms):")
        CALL log_write("===============================================================================")
        CALL log_write("  Atom   Symbol     Position X       Position Y       Position Z     Molecule  ")
        CALL log_write("===============================================================================")
          DO i=1,n_atoms
            WRITE(line,'(1x,i4,6x,a2,4x,3(f14.10,3x),2x,i5)') i, atom_symbol(data_manbo(i)%an), data_manbo(i)%r(1),&
                                                           data_manbo(i)%r(2), data_manbo(i)%r(3), data_manbo(i)%mol_orig
            CALL log_write(TRIM(line))
          END DO
        CALL log_write("-------------------------------------------------------------------------------")
      END IF
      IF (vector_module(box_dim)<1E-5) THEN
        CALL reorient_box(1)
        CALL log_write("")
        CALL log_write("WARNING: As the box dimensions were not given in the input file, we reoriented the box (with rotations")
        CALL log_write("         and translations) and got the box dimensions")
        CALL log_write("")
        IF(use_embedding) THEN
          CALL log_write("Reoriented values entered for the atoms (Coordinates are given in Angstroms and charges in A.U.):")
          CALL log_write("===============================================================================================")
          CALL log_write("  Atom   Symbol     Position X       Position Y       Position Z     Molecule    Emb. Charge   ")
          CALL log_write("===============================================================================================")
            DO i=1,n_atoms
              WRITE(line,'(1x,i4,6x,a2,4x,3(f14.10,3x),2x,i5,6x,f10.6)') i, atom_symbol(data_manbo(i)%an), data_manbo(i)%r(1),&
                                                                          data_manbo(i)%r(2), data_manbo(i)%r(3),&
                                                                          data_manbo(i)%mol_orig,data_manbo(i)%q
              CALL log_write(TRIM(line))
            END DO
          CALL log_write("-----------------------------------------------------------------------------------------------")
        ELSE
          CALL log_write("Reoriented values entered for the atoms (Coordinates are given in Angstroms):")
          CALL log_write("===============================================================================")
          CALL log_write("  Atom   Symbol     Position X       Position Y       Position Z     Molecule  ")
          CALL log_write("===============================================================================")
            DO i=1,n_atoms
              WRITE(line,'(1x,i4,6x,a2,4x,3(f14.10,3x),2x,i5)') i, atom_symbol(data_manbo(i)%an), data_manbo(i)%r(1),&
                                                             data_manbo(i)%r(2), data_manbo(i)%r(3), data_manbo(i)%mol_orig
              CALL log_write(TRIM(line))
            END DO
          CALL log_write("-------------------------------------------------------------------------------")
        END IF
      ELSE IF (apply_pbc) THEN
        CALL reorient_box(2)
      END IF
      CALL log_write("Dimensions of the box:")
      WRITE(line,'(f10.6)') box_dim(1)
      CALL log_write("  x: " // TRIM(ADJUSTL(line)) // " A")
      WRITE(line,'(f10.6)') box_dim(2)
      CALL log_write("  y: " // TRIM(ADJUSTL(line)) // " A")
      WRITE(line,'(f10.6)') box_dim(3)
      CALL log_write("  z: " // TRIM(ADJUSTL(line)) // " A")
      WRITE(line,'(f12.6)') box_dim(1)*box_dim(2)*box_dim(3)
      CALL log_write("Volume: " // TRIM(ADJUSTL(line)) // " A^3")
      val = 0.0
        DO i=1,n_atoms
          val = val + data_manbo(i)%mass
        END DO
      WRITE(line,'(f12.4)') val
      CALL log_write("Total mass: " // TRIM(ADJUSTL(line)) // " a.u.")
      WRITE(line,'(f12.4)') val*1.66053886/(box_dim(1)*box_dim(2)*box_dim(3))
      CALL log_write("Density: " // TRIM(ADJUSTL(line)) // " g/cm^3")
      IF(cutoff>(MINVAL(box_dim)/2) .AND. apply_pbc) THEN
        CALL log_write("")
        CALL log_write("  WARNING: Entered value for the Cutoff Radius is bigger than half of the smaller side of the box.")
      END IF
      CALL log_write("")
      CALL log_write("Charge and Multiplicity of the Molecules:")
      CALL log_write("====================================")
      CALL log_write("  Molecule   Charge   Multiplicity  ")
      CALL log_write("====================================")
        DO i=1,n_mols
          WRITE(line,'(4x,i4,5x,i3,10x,i2)') i, char_mul_orig(i)%q_mol, char_mul_orig(i)%mul
          CALL log_write(TRIM(line))
        END DO
      CALL log_write("------------------------------------")
      CALL log_write("")
        IF (.NOT. use_rand_vel) THEN
          CALL log_write("Initial velocities of Atoms, obtained from the input file (given in Angstroms/picosecond):")
        ELSE
          CALL log_write("Initial velocities of Atoms, generated randomly (given in Angstroms/picosecond):")
        END IF
      CALL log_write("===========================================================")
      CALL log_write("  Atom     Velocity X       Velocity Y       Velocity Z    ")
      CALL log_write("===========================================================")
        DO i=1,n_atoms
          WRITE(line,'(1x,i4,3x,3(f14.10,3x))') i, data_manbo(i)%v(1), data_manbo(i)%v(2), data_manbo(i)%v(3)
          CALL log_write(TRIM(line))
        END DO
      CALL log_write("-----------------------------------------------------------")
      CALL log_write("")
      ! Here we sent some initial values to ManBo's log
      CLOSE(13)
      ! Closing input file
    END IF
  END SUBROUTINE read_input

END MODULE manbo_input_reading
