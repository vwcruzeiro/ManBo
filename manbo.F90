! ==============================================================================
!                                    ManBo                              
!       Ab Initio Molecular Dynamics Program that uses Many Body Expansion      
! ==============================================================================
!      Created by Vinícius Wilian Dias Cruzeiro and Herbert de Castro Georg
!
! Email: vwcruzeiro@ufl.edu
!

PROGRAM manbo

#ifdef USE_PARALLEL
  USE mpi
#endif /* USE_PARALLEL */

  USE manbo_variables
  ! Obs.: manbo_periodic_table and manbo_parameters are called on manbo_variables module
  USE manbo_input_reading
  USE manbo_forces
  USE manbo_subroutines
  USE manbo_log
  ! Here we call some ManBo's modules

  IMPLICIT NONE
  CHARACTER(LEN=50) :: line
  INTEGER :: temp_dt, ierr, mpi_id, mpi_size
  ! mpi_id is the id of a MPI thread
  ! mpi_size is the total number of MPI threads
  INTEGER, DIMENSION(5) :: dif_temp_dt
  
#ifdef USE_PARALLEL  
  call MPI_INIT ( ierr )

  call MPI_COMM_RANK (MPI_COMM_WORLD, mpi_id, ierr)
  call MPI_COMM_SIZE (MPI_COMM_WORLD, mpi_size, ierr)
#else
  mpi_id = 0
  mpi_size = 1
#endif /* USE_PARALLEL */
  
  CALL read_input(mpi_id,mpi_size)
  ! Here we read and initiate the main variables of the program.
  ! Initial positions (and possibly velocities) are read.

  IF (mpi_id == 0) THEN
    
    md_step = 0
    ! Initial step is setted
    
    CALL log_write("====================================================================")
    CALL log_write(" Step        0 - Calculating the initial forces                     ")
    CALL log_write("====================================================================")
    CALL log_write("")
    
    IF (group_monomers > 1) CALL group_mons
    ! Grouping molecules single monomers
  END IF
  
  CALL calculate_forces(mpi_id)
  ! Initiallize the accelerations
  
  IF (mpi_id == 0) THEN
    CALL save_data_manbo(1)
    CALL print_properties(1)
    ! Saving data of step 0
    
    CALL log_write("")
    CALL log_time("  This step ended at")
    IF(n_steps/=0) THEN
      CALL SYSTEM_CLOCK(temp_dt)
      temp_dt=NINT(REAL(temp_dt - init_dt)/(md_step + 1))
      CALL dt_diference(0,temp_dt,dif_temp_dt)
      WRITE(100,'("   Average time per step: ",i4," days, ",i2," hours, ",i2," minutes, and ",f6.3," seconds.")')&
                   dif_temp_dt(1), dif_temp_dt(2), dif_temp_dt(3), (REAL(dif_temp_dt(4)) + REAL(dif_temp_dt(5))/1000)
      temp_dt=temp_dt*(n_steps - md_step)
      CALL dt_diference(0,temp_dt,dif_temp_dt)
      WRITE(100,'("   The execution is estimated to end in: ",i4," days, ",i2," hours, ",i2," minutes, and ",f6.3," seconds.")')&
                   dif_temp_dt(1), dif_temp_dt(2), dif_temp_dt(3), (REAL(dif_temp_dt(4)) + REAL(dif_temp_dt(5))/1000)
    END IF
    CALL log_write("--------------------------------------------------------------------")
    CALL log_write("")
  END IF
    
    dynamics: DO md_step = 1, n_steps
      ! Update the time
    
      IF (mpi_id == 0) THEN
        CALL log_write("====================================================================")
        WRITE(line,'(" Step ",i8," of ",i8)') md_step, n_steps
        CALL log_write(TRIM(line))
        CALL log_write("====================================================================")
        CALL log_write("")
        
        CALL propagate_positions
        ! Update positions [ r(t) --> r(t+dt) ]  
        CALL propagate_velocities
        ! Update velocities [ v(t) --> v(t+dt/2) ], according to Velocity Verlet algorithm
        
        IF (group_monomers > 1) CALL group_mons
        ! Grouping molecules single monomers
      END IF
      
      CALL calculate_forces(mpi_id)
      ! Calculate energy & forces and update accelerations [ a(t) --> a(t+dt) ]
      
      IF (mpi_id == 0) THEN
        CALL propagate_velocities
        ! Complete the velocity updating [ v(t+dt/2) --> v(t+dt) ]
        
        IF (MOD(md_step, i_print_prop) == 0 .OR. md_step == n_steps) THEN
          CALL print_properties(2)
          ! Calculate and save all properties of the system using the actual positions and velocities
        END IF
        
        IF(apply_therm) THEN
          CALL apply_thermalization
        END IF
        ! If required, we use the Berendsen thermostat
        
        IF (MOD(md_step, i_save) == 0 .OR. md_step == n_steps) THEN
          CALL save_data_manbo(2)
          ! Save positions, velocities and other data as output for future visualization and analysis, and for restarting the simulation
        END IF
        
        IF ((MOD(md_step, i_rst) == 0 .OR. md_step == n_steps) .AND. i_rst > 0) THEN
          CALL save_rst_file()
          ! Save the restart file
        END IF
        
        CALL log_write("")
        CALL log_time("  This step ended at")
        IF(n_steps/=md_step) THEN
          CALL SYSTEM_CLOCK(temp_dt)
          temp_dt=NINT(REAL(temp_dt - init_dt)/(md_step + 1))
          CALL dt_diference(0,temp_dt,dif_temp_dt)
          WRITE(100,'("   Average time per step: ",i4," days, ",i2," hours, ",i2," minutes, and ",f6.3," seconds.")')&
                       dif_temp_dt(1), dif_temp_dt(2), dif_temp_dt(3), (REAL(dif_temp_dt(4)) + REAL(dif_temp_dt(5))/1000)
          temp_dt=temp_dt*(n_steps - md_step)
          CALL dt_diference(0,temp_dt,dif_temp_dt)
          WRITE(100,'("   The execution is estimated to end in ",i5," days, ",i2," hours, ",i2," minutes, and ",f6.3," seconds.")')&
                       dif_temp_dt(1), dif_temp_dt(2), dif_temp_dt(3), (REAL(dif_temp_dt(4)) + REAL(dif_temp_dt(5))/1000)
        END IF
        CALL log_write("--------------------------------------------------------------------")
        CALL log_write("")
      END IF
    END DO dynamics

  IF (mpi_id == 0) THEN    
    CALL log_close(0)
    ! Close the log file of ManBo
  END IF
  
#ifdef USE_PARALLEL
  call MPI_FINALIZE ( ierr )
#endif /* USE_PARALLEL */

END PROGRAM manbo
