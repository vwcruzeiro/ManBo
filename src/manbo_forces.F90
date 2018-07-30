MODULE manbo_forces

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

  SUBROUTINE box_replication
  ! This subroutine replicates the box
  
    IMPLICIT NONE

    INTEGER :: i, j, k, k2, num, x, n1, n2, n3
    DOUBLE PRECISION :: lx2, ly2, lz2
    ! Variables to store the values of half of the box lenght in each dimension
    DOUBLE PRECISION :: mcx, mcy, mcz
    DOUBLE PRECISION :: mcx2, mcy2, mcz2
    DOUBLE PRECISION :: mcx3, mcy3, mcz3
    ! Coordinates of the center of mass of a molecule
    DOUBLE PRECISION :: dx, dy, dz, dxy, dxz, dyz, dxyz, L
    DOUBLE PRECISION :: dx2, dy2, dz2, dxy2, dxz2, dyz2, dxyz2
    DOUBLE PRECISION :: dx3, dy3, dz3, dxy3, dxz3, dyz3, dxyz3
    INTEGER :: allocate_status
    CHARACTER(LEN=100) :: line
    LOGICAL, DIMENSION(3) :: conditions
    TYPE(atom_out), DIMENSION(:), ALLOCATABLE, SAVE :: data_manbo_temp
    
    CALL calculate_center_of_mass(1)
    ! Here we calculate the mass center's of the inside box's molecules
      
    lx2 = box_dim(1)/2
    ly2 = box_dim(2)/2
    lz2 = box_dim(3)/2
  
    k = n_mols
    IF (group_monomers>1) k2 = n_mols_orig
    num = 0

    L = MAXVAL(box_dim)
    x = NINT(2.0*(1.25*(8*((cutoff/L)**3) + 12*((cutoff/L)**2) + 6*(cutoff/L)) + 1.0/n_mols)*n_atoms)
    ! This is an estimate for allocate data_manbo_out

    DEALLOCATE(data_manbo_out, stat=allocate_status)
      IF(allocate_status/=0) THEN
        PRINT *, "No memory enough to deallocate. ManBo cannot run."
        CALL log_write("ERROR: Error on deallocate data_manbo_out on manbo_forces.F90")
        CALL log_close(1)
        STOP
      END IF
    
    ALLOCATE(data_manbo_out(x), stat=allocate_status)
      IF(allocate_status/=0) THEN
        PRINT *, "No memory enough to allocate. ManBo cannot run."
        CALL log_write("ERROR: Error on allocate data_manbo_out on manbo_forces.F90")
        CALL log_close(1)
        STOP
      END IF
    ! Here we allocate data_manbo_out with no data

    IF (group_monomers>1) THEN
      IF (.not. ALLOCATED(orig_mols_out)) THEN
        ALLOCATE(orig_mols_out(NINT(float(x*n_mols)/float(n_atoms))), stat=allocate_status)
          IF(allocate_status/=0) THEN
            PRINT *, "No memory enough to allocate. ManBo cannot run."
            CALL log_write("ERROR: Error on allocate orig_mols_out on manbo_forces.F90")
            CALL log_close(1)
            STOP
          END IF
      END IF
      ! Here we allocate orig_mols_out with no data
      orig_mols_out%m(1) = 0
      orig_mols_out%m(2) = 0
      orig_mols_out%m(3) = 0
    END IF

    DO i=1,n_mols
      IF (.not. group_monomers>1) THEN
        mcx = mc(i)%r(1)
        mcy = mc(i)%r(2)
        mcz = mc(i)%r(3)
        dx = lx2 - ABS(mcx)
        dy = ly2 - ABS(mcy)
        dz = lz2 - ABS(mcz)
        dxy = SQRT(dx*dx + dy*dy)
        dxz = SQRT(dx*dx + dz*dz)
        dyz = SQRT(dy*dy + dz*dz)
        dxyz = SQRT(dx*dx + dy*dy + dz*dz)
      ELSE
        n1 = orig_mols(i)%m(1)
        n2 = orig_mols(i)%m(2)
        n3 = orig_mols(i)%m(3)
        
        mcx = mc_orig(n1)%r(1)
        mcy = mc_orig(n1)%r(2)
        mcz = mc_orig(n1)%r(3)
        dx = lx2 - ABS(mcx)
        dy = ly2 - ABS(mcy)
        dz = lz2 - ABS(mcz)
        dxy = SQRT(dx*dx + dy*dy)
        dxz = SQRT(dx*dx + dz*dz)
        dyz = SQRT(dy*dy + dz*dz)
        dxyz = SQRT(dx*dx + dy*dy + dz*dz)
        
        IF (n2 .ne. 0) THEN
          mcx2 = mc_orig(n2)%r(1)
          mcy2 = mc_orig(n2)%r(2)
          mcz2 = mc_orig(n2)%r(3)
          dx2 = lx2 - ABS(mcx2)
          dy2 = ly2 - ABS(mcy2)
          dz2 = lz2 - ABS(mcz2)
          dxy2 = SQRT(dx2*dx2 + dy2*dy2)
          dxz2 = SQRT(dx2*dx2 + dz2*dz2)
          dyz2 = SQRT(dy2*dy2 + dz2*dz2)
          dxyz2 = SQRT(dx2*dx2 + dy2*dy2 + dz2*dz2)
        END IF
        
        IF (n3 .ne. 0) THEN
          mcx3 = mc_orig(n3)%r(1)
          mcy3 = mc_orig(n3)%r(2)
          mcz3 = mc_orig(n3)%r(3)
          dx3 = lx2 - ABS(mcx3)
          dy3 = ly2 - ABS(mcy3)
          dz3 = lz2 - ABS(mcz3)
          dxy3 = SQRT(dx3*dx3 + dy3*dy3)
          dxz3 = SQRT(dx3*dx3 + dz3*dz3)
          dyz3 = SQRT(dy3*dy3 + dz3*dz3)
          dxyz3 = SQRT(dx3*dx3 + dy3*dy3 + dz3*dz3)
        END IF
      END IF
  
      IF (.not. group_monomers>1) THEN
        conditions(1) = dx < cutoff
        conditions(2) = .FALSE.
        conditions(3) = .FALSE.
      ELSE
        conditions(1) = dx < cutoff
        IF (n2 .ne. 0) THEN
          conditions(2) = dx2 < cutoff
        ELSE
          conditions(2) = .FALSE.
        END IF
        IF (n3 .ne. 0) THEN
          conditions(3) = dx3 < cutoff
        ELSE
          conditions(3) = .FALSE.
        END IF
      END IF
      IF(conditions(1) .or. conditions(2) .or. conditions(3)) THEN
        k = k + 1
        orig_repli_mol(k - n_mols) = i
          DO j=1,n_atoms
            IF(data_manbo(j)%mol==i) THEN
              num = num + 1
              IF (conditions(1)) THEN
                data_manbo_out(num)%r(1) = data_manbo(j)%r(1) - box_dim(1)*sig(mcx)
              ELSE IF (conditions(2)) THEN
                data_manbo_out(num)%r(1) = data_manbo(j)%r(1) - box_dim(1)*sig(mcx2)
              ELSE
                data_manbo_out(num)%r(1) = data_manbo(j)%r(1) - box_dim(1)*sig(mcx3)
              END IF
              data_manbo_out(num)%r(2) = data_manbo(j)%r(2)
              data_manbo_out(num)%r(3) = data_manbo(j)%r(3)
              data_manbo_out(num)%an   = data_manbo(j)%an
              data_manbo_out(num)%mass = data_manbo(j)%mass
              data_manbo_out(num)%mol  = k
              IF (group_monomers>1) THEN
                IF (data_manbo(j)%mol_orig==n1) THEN
                  data_manbo_out(num)%mol_orig = k2 + 1
                ELSE IF (data_manbo(j)%mol_orig==n2) THEN
                  data_manbo_out(num)%mol_orig = k2 + 2
                ELSE IF (data_manbo(j)%mol_orig==n3) THEN
                  data_manbo_out(num)%mol_orig = k2 + 3              
                END IF
              END IF
              data_manbo_out(num)%atom = j
              IF(j==(n_atoms+1)) EXIT
            END IF
          END DO
        IF (group_monomers>1) THEN
          k2 = k2 + 1
          orig_mols_out(k - n_mols)%m(1) = k2
          IF (.not. n2==0) THEN
            k2 = k2 + 1
            orig_mols_out(k - n_mols)%m(2) = k2
          END IF
          IF (.not. n3==0) THEN
            k2 = k2 + 1
            orig_mols_out(k - n_mols)%m(3) = k2
          END IF
        END IF
      END IF
  
      IF (.not. group_monomers>1) THEN
        conditions(1) = dy < cutoff
        conditions(2) = .FALSE.
        conditions(3) = .FALSE.
      ELSE
        conditions(1) = dy < cutoff
        IF (n2 .ne. 0) THEN
          conditions(2) = dy2 < cutoff
        ELSE
          conditions(2) = .FALSE.
        END IF
        IF (n3 .ne. 0) THEN
          conditions(3) = dy3 < cutoff
        ELSE
          conditions(3) = .FALSE.
        END IF
      END IF
      IF(conditions(1) .or. conditions(2) .or. conditions(3)) THEN
        k = k + 1
        orig_repli_mol(k - n_mols) = i
          DO j=1,n_atoms
            IF(data_manbo(j)%mol==i) THEN
              num = num + 1
              data_manbo_out(num)%r(1) = data_manbo(j)%r(1)
              IF (conditions(1)) THEN
                data_manbo_out(num)%r(2) = data_manbo(j)%r(2) - box_dim(2)*sig(mcy)
              ELSE IF (conditions(2)) THEN
                data_manbo_out(num)%r(2) = data_manbo(j)%r(2) - box_dim(2)*sig(mcy2)
              ELSE
                data_manbo_out(num)%r(2) = data_manbo(j)%r(2) - box_dim(2)*sig(mcy3)
              END IF
              data_manbo_out(num)%r(3) = data_manbo(j)%r(3)
              data_manbo_out(num)%an   = data_manbo(j)%an
              data_manbo_out(num)%mass = data_manbo(j)%mass
              data_manbo_out(num)%mol  = k
              IF (group_monomers>1) THEN
                IF (data_manbo(j)%mol_orig==n1) THEN
                  data_manbo_out(num)%mol_orig  = k2 + 1
                ELSE IF (data_manbo(j)%mol_orig==n2) THEN
                  data_manbo_out(num)%mol_orig  = k2 + 2
                ELSE IF (data_manbo(j)%mol_orig==n3) THEN
                  data_manbo_out(num)%mol_orig  = k2 + 3              
                END IF
              END IF
              data_manbo_out(num)%atom = j
              IF(j==(n_atoms+1)) EXIT
            END IF
          END DO
        IF (group_monomers>1) THEN
          k2 = k2 + 1
          orig_mols_out(k - n_mols)%m(1) = k2
          IF (.not. n2==0) THEN
            k2 = k2 + 1
            orig_mols_out(k - n_mols)%m(2) = k2
          END IF
          IF (.not. n3==0) THEN
            k2 = k2 + 1
            orig_mols_out(k - n_mols)%m(3) = k2
          END IF
        END IF
      END IF
  
      IF (.not. group_monomers>1) THEN
        conditions(1) = dz < cutoff
        conditions(2) = .FALSE.
        conditions(3) = .FALSE.
      ELSE
        conditions(1) = dz < cutoff
        IF (n2 .ne. 0) THEN
          conditions(2) = dz2 < cutoff
        ELSE
          conditions(2) = .FALSE.
        END IF
        IF (n3 .ne. 0) THEN
          conditions(3) = dz3 < cutoff
        ELSE
          conditions(3) = .FALSE.
        END IF
      END IF
      IF(conditions(1) .or. conditions(2) .or. conditions(3)) THEN
        k = k + 1
        orig_repli_mol(k - n_mols) = i
          DO j=1,n_atoms
            IF(data_manbo(j)%mol==i) THEN
              num = num + 1
              data_manbo_out(num)%r(1) = data_manbo(j)%r(1)
              data_manbo_out(num)%r(2) = data_manbo(j)%r(2)
              IF (conditions(1)) THEN
                data_manbo_out(num)%r(3) = data_manbo(j)%r(3) - box_dim(3)*sig(mcz)
              ELSE IF (conditions(2)) THEN
                data_manbo_out(num)%r(3) = data_manbo(j)%r(3) - box_dim(3)*sig(mcz2)
              ELSE
                data_manbo_out(num)%r(3) = data_manbo(j)%r(3) - box_dim(3)*sig(mcz3)
              END IF
              data_manbo_out(num)%an   = data_manbo(j)%an
              data_manbo_out(num)%mass = data_manbo(j)%mass
              data_manbo_out(num)%mol  = k
              IF (group_monomers>1) THEN
                IF (data_manbo(j)%mol_orig==n1) THEN
                  data_manbo_out(num)%mol_orig  = k2 + 1
                ELSE IF (data_manbo(j)%mol_orig==n2) THEN
                  data_manbo_out(num)%mol_orig  = k2 + 2
                ELSE IF (data_manbo(j)%mol_orig==n3) THEN
                  data_manbo_out(num)%mol_orig  = k2 + 3              
                END IF
              END IF
              data_manbo_out(num)%atom = j
              IF(j==(n_atoms+1)) EXIT
            END IF
          END DO
        IF (group_monomers>1) THEN
          k2 = k2 + 1
          orig_mols_out(k - n_mols)%m(1) = k2
          IF (.not. n2==0) THEN
            k2 = k2 + 1
            orig_mols_out(k - n_mols)%m(2) = k2
          END IF
          IF (.not. n3==0) THEN
            k2 = k2 + 1
            orig_mols_out(k - n_mols)%m(3) = k2
          END IF
        END IF
      END IF
  
      IF (.not. group_monomers>1) THEN
        conditions(1) = dxy < cutoff
        conditions(2) = .FALSE.
        conditions(3) = .FALSE.
      ELSE
        conditions(1) = dxy < cutoff
        IF (n2 .ne. 0) THEN
          conditions(2) = dxy2 < cutoff
        ELSE
          conditions(2) = .FALSE.
        END IF
        IF (n3 .ne. 0) THEN
          conditions(3) = dxy3 < cutoff
        ELSE
          conditions(3) = .FALSE.
        END IF
      END IF
      IF(conditions(1) .or. conditions(2) .or. conditions(3)) THEN
        k = k + 1
        orig_repli_mol(k - n_mols) = i
          DO j=1,n_atoms
            IF(data_manbo(j)%mol==i) THEN
              num = num + 1
              IF (conditions(1)) THEN
                data_manbo_out(num)%r(1) = data_manbo(j)%r(1) - box_dim(1)*sig(mcx)
                data_manbo_out(num)%r(2) = data_manbo(j)%r(2) - box_dim(2)*sig(mcy)
              ELSE IF (conditions(2)) THEN
                data_manbo_out(num)%r(1) = data_manbo(j)%r(1) - box_dim(1)*sig(mcx2)
                data_manbo_out(num)%r(2) = data_manbo(j)%r(2) - box_dim(2)*sig(mcy2)
              ELSE
                data_manbo_out(num)%r(1) = data_manbo(j)%r(1) - box_dim(1)*sig(mcx3)
                data_manbo_out(num)%r(2) = data_manbo(j)%r(2) - box_dim(2)*sig(mcy3)
              END IF
              data_manbo_out(num)%r(3) = data_manbo(j)%r(3) 
              data_manbo_out(num)%an   = data_manbo(j)%an
              data_manbo_out(num)%mass = data_manbo(j)%mass
              data_manbo_out(num)%mol  = k
              IF (group_monomers>1) THEN
                IF (data_manbo(j)%mol_orig==n1) THEN
                  data_manbo_out(num)%mol_orig  = k2 + 1
                ELSE IF (data_manbo(j)%mol_orig==n2) THEN
                  data_manbo_out(num)%mol_orig  = k2 + 2
                ELSE IF (data_manbo(j)%mol_orig==n3) THEN
                  data_manbo_out(num)%mol_orig  = k2 + 3              
                END IF
              END IF
              data_manbo_out(num)%atom = j
              IF(j==(n_atoms+1)) EXIT
            END IF
          END DO
        IF (group_monomers>1) THEN
          k2 = k2 + 1
          orig_mols_out(k - n_mols)%m(1) = k2
          IF (.not. n2==0) THEN
            k2 = k2 + 1
            orig_mols_out(k - n_mols)%m(2) = k2
          END IF
          IF (.not. n3==0) THEN
            k2 = k2 + 1
            orig_mols_out(k - n_mols)%m(3) = k2
          END IF
        END IF
      END IF
  
      IF (.not. group_monomers>1) THEN
        conditions(1) = dxz < cutoff
        conditions(2) = .FALSE.
        conditions(3) = .FALSE.
      ELSE
        conditions(1) = dxz < cutoff
        IF (n2 .ne. 0) THEN
          conditions(2) = dxz2 < cutoff
        ELSE
          conditions(2) = .FALSE.
        END IF
        IF (n3 .ne. 0) THEN
          conditions(3) = dxz3 < cutoff
        ELSE
          conditions(3) = .FALSE.
        END IF
      END IF
      IF(conditions(1) .or. conditions(2) .or. conditions(3)) THEN
        k = k + 1
        orig_repli_mol(k - n_mols) = i
          DO j=1,n_atoms
            IF(data_manbo(j)%mol==i) THEN
              num = num + 1
              IF (conditions(1)) THEN
                data_manbo_out(num)%r(1) = data_manbo(j)%r(1) - box_dim(1)*sig(mcx)
                data_manbo_out(num)%r(3) = data_manbo(j)%r(3) - box_dim(3)*sig(mcz)
              ELSE IF (conditions(2)) THEN
                data_manbo_out(num)%r(1) = data_manbo(j)%r(1) - box_dim(1)*sig(mcx2)
                data_manbo_out(num)%r(3) = data_manbo(j)%r(3) - box_dim(3)*sig(mcz2)
              ELSE
                data_manbo_out(num)%r(1) = data_manbo(j)%r(1) - box_dim(1)*sig(mcx3)
                data_manbo_out(num)%r(3) = data_manbo(j)%r(3) - box_dim(3)*sig(mcz3)
              END IF
              data_manbo_out(num)%r(2) = data_manbo(j)%r(2)
              data_manbo_out(num)%an   = data_manbo(j)%an
              data_manbo_out(num)%mass = data_manbo(j)%mass
              data_manbo_out(num)%mol  = k
              IF (group_monomers>1) THEN
                IF (data_manbo(j)%mol_orig==n1) THEN
                  data_manbo_out(num)%mol_orig  = k2 + 1
                ELSE IF (data_manbo(j)%mol_orig==n2) THEN
                  data_manbo_out(num)%mol_orig  = k2 + 2
                ELSE IF (data_manbo(j)%mol_orig==n3) THEN
                  data_manbo_out(num)%mol_orig  = k2 + 3              
                END IF
              END IF
              data_manbo_out(num)%atom = j
              IF(j==(n_atoms+1)) EXIT
            END IF
          END DO
        IF (group_monomers>1) THEN
          k2 = k2 + 1
          orig_mols_out(k - n_mols)%m(1) = k2
          IF (.not. n2==0) THEN
            k2 = k2 + 1
            orig_mols_out(k - n_mols)%m(2) = k2
          END IF
          IF (.not. n3==0) THEN
            k2 = k2 + 1
            orig_mols_out(k - n_mols)%m(3) = k2
          END IF
        END IF
      END IF
  
      IF (.not. group_monomers>1) THEN
        conditions(1) = dyz < cutoff
        conditions(2) = .FALSE.
        conditions(3) = .FALSE.
      ELSE
        conditions(1) = dyz < cutoff
        IF (n2 .ne. 0) THEN
          conditions(2) = dyz2 < cutoff
        ELSE
          conditions(2) = .FALSE.
        END IF
        IF (n3 .ne. 0) THEN
          conditions(3) = dyz3 < cutoff
        ELSE
          conditions(3) = .FALSE.
        END IF
      END IF
      IF(conditions(1) .or. conditions(2) .or. conditions(3)) THEN
        k = k + 1
        orig_repli_mol(k - n_mols) = i
          DO j=1,n_atoms
            IF(data_manbo(j)%mol==i) THEN
              num = num + 1
              data_manbo_out(num)%r(1) = data_manbo(j)%r(1)
              IF (conditions(1)) THEN
                data_manbo_out(num)%r(2) = data_manbo(j)%r(2) - box_dim(2)*sig(mcy)
                data_manbo_out(num)%r(3) = data_manbo(j)%r(3) - box_dim(3)*sig(mcz)
              ELSE IF (conditions(2)) THEN
                data_manbo_out(num)%r(2) = data_manbo(j)%r(2) - box_dim(2)*sig(mcy2)
                data_manbo_out(num)%r(3) = data_manbo(j)%r(3) - box_dim(3)*sig(mcz2)
              ELSE
                data_manbo_out(num)%r(2) = data_manbo(j)%r(2) - box_dim(2)*sig(mcy3)
                data_manbo_out(num)%r(3) = data_manbo(j)%r(3) - box_dim(3)*sig(mcz3)
              END IF
              data_manbo_out(num)%an   = data_manbo(j)%an
              data_manbo_out(num)%mass = data_manbo(j)%mass
              data_manbo_out(num)%mol  = k
              IF (group_monomers>1) THEN
                IF (data_manbo(j)%mol_orig==n1) THEN
                  data_manbo_out(num)%mol_orig  = k2 + 1
                ELSE IF (data_manbo(j)%mol_orig==n2) THEN
                  data_manbo_out(num)%mol_orig  = k2 + 2
                ELSE IF (data_manbo(j)%mol_orig==n3) THEN
                  data_manbo_out(num)%mol_orig  = k2 + 3              
                END IF
              END IF
              data_manbo_out(num)%atom = j
              IF(j==(n_atoms+1)) EXIT
            END IF
          END DO
        IF (group_monomers>1) THEN
          k2 = k2 + 1
          orig_mols_out(k - n_mols)%m(1) = k2
          IF (.not. n2==0) THEN
            k2 = k2 + 1
            orig_mols_out(k - n_mols)%m(2) = k2
          END IF
          IF (.not. n3==0) THEN
            k2 = k2 + 1
            orig_mols_out(k - n_mols)%m(3) = k2
          END IF
        END IF
      END IF
  
      IF (.not. group_monomers>1) THEN
        conditions(1) = dxyz < cutoff
        conditions(2) = .FALSE.
        conditions(3) = .FALSE.
      ELSE
        conditions(1) = dxyz < cutoff
        IF (n2 .ne. 0) THEN
          conditions(2) = dxyz2 < cutoff
        ELSE
          conditions(2) = .FALSE.
        END IF
        IF (n3 .ne. 0) THEN
          conditions(3) = dxyz3 < cutoff
        ELSE
          conditions(3) = .FALSE.
        END IF
      END IF
      IF(conditions(1) .or. conditions(2) .or. conditions(3)) THEN
        k = k + 1
        orig_repli_mol(k - n_mols) = i
          DO j=1,n_atoms
            IF(data_manbo(j)%mol==i) THEN
              num = num + 1
              IF (conditions(1)) THEN
                data_manbo_out(num)%r(1) = data_manbo(j)%r(1) - box_dim(1)*sig(mcx)
                data_manbo_out(num)%r(2) = data_manbo(j)%r(2) - box_dim(2)*sig(mcy)
                data_manbo_out(num)%r(3) = data_manbo(j)%r(3) - box_dim(3)*sig(mcz)
              ELSE IF (conditions(2)) THEN
                data_manbo_out(num)%r(1) = data_manbo(j)%r(1) - box_dim(1)*sig(mcx2)
                data_manbo_out(num)%r(2) = data_manbo(j)%r(2) - box_dim(2)*sig(mcy2)
                data_manbo_out(num)%r(3) = data_manbo(j)%r(3) - box_dim(3)*sig(mcz2)
              ELSE
                data_manbo_out(num)%r(1) = data_manbo(j)%r(1) - box_dim(1)*sig(mcx3)
                data_manbo_out(num)%r(2) = data_manbo(j)%r(2) - box_dim(2)*sig(mcy3)
                data_manbo_out(num)%r(3) = data_manbo(j)%r(3) - box_dim(3)*sig(mcz3)
              END IF
              data_manbo_out(num)%an   = data_manbo(j)%an
              data_manbo_out(num)%mass = data_manbo(j)%mass
              data_manbo_out(num)%mol  = k
              IF (group_monomers>1) THEN
                IF (data_manbo(j)%mol_orig==n1) THEN
                  data_manbo_out(num)%mol_orig  = k2 + 1
                ELSE IF (data_manbo(j)%mol_orig==n2) THEN
                  data_manbo_out(num)%mol_orig  = k2 + 2
                ELSE IF (data_manbo(j)%mol_orig==n3) THEN
                  data_manbo_out(num)%mol_orig  = k2 + 3              
                END IF
              END IF
              data_manbo_out(num)%atom = j
              IF(j==(n_atoms+1)) EXIT
            END IF
          END DO
        IF (group_monomers>1) THEN
          k2 = k2 + 1
          orig_mols_out(k - n_mols)%m(1) = k2
          IF (.not. n2==0) THEN
            k2 = k2 + 1
            orig_mols_out(k - n_mols)%m(2) = k2
          END IF
          IF (.not. n3==0) THEN
            k2 = k2 + 1
            orig_mols_out(k - n_mols)%m(3) = k2
          END IF
        END IF
      END IF
    END DO

    WRITE(line,'("  ",i4," atom(s) of ",i4," molecule(s) were successfully replicated.")') num, (k - n_mols)
    CALL log_write(TRIM(line))
    CALL log_write("")

    ALLOCATE(data_manbo_temp(num), stat=allocate_status)
      IF(allocate_status/=0) THEN
        PRINT *, "No memory enough to allocate. ManBo cannot run."
        CALL log_write("ERROR: Error on allocate data_manbo_temp on manbo_forces.F90")
        CALL log_close(1)
        STOP
      END IF
      DO i=1,num
        data_manbo_temp(i)%r(1) = data_manbo_out(i)%r(1)
        data_manbo_temp(i)%r(2) = data_manbo_out(i)%r(2)
        data_manbo_temp(i)%r(3) = data_manbo_out(i)%r(3)
        data_manbo_temp(i)%mass = data_manbo_out(i)%mass
        data_manbo_temp(i)%an = data_manbo_out(i)%an
        data_manbo_temp(i)%mol = data_manbo_out(i)%mol
        IF (group_monomers>1) data_manbo_temp(i)%mol_orig = data_manbo_out(i)%mol_orig
        data_manbo_temp(i)%atom = data_manbo_out(i)%atom
      END DO

    DEALLOCATE(data_manbo_out, stat=allocate_status)
      IF(allocate_status/=0) THEN
        PRINT *, "No memory enough to deallocate. ManBo cannot run."
        CALL log_write("ERROR: Error on deallocate data_manbo_out on manbo_forces.F90")
        CALL log_close(1)
        STOP
      END IF

    ALLOCATE(data_manbo_out(num), stat=allocate_status)
      IF(allocate_status/=0) THEN
        PRINT *, "No memory enough to allocate. ManBo cannot run."
        CALL log_write("ERROR: Error on allocate data_manbo_out on manbo_forces.F90")
        CALL log_close(1)
        STOP
      END IF
      DO i=1,num
        data_manbo_out(i)%r(1) = data_manbo_temp(i)%r(1)
        data_manbo_out(i)%r(2) = data_manbo_temp(i)%r(2)
        data_manbo_out(i)%r(3) = data_manbo_temp(i)%r(3)
        data_manbo_out(i)%mass = data_manbo_temp(i)%mass
        data_manbo_out(i)%an = data_manbo_temp(i)%an
        data_manbo_out(i)%mol = data_manbo_temp(i)%mol
        IF (group_monomers>1) data_manbo_out(i)%mol_orig = data_manbo_temp(i)%mol_orig
        data_manbo_out(i)%atom = data_manbo_temp(i)%atom
      END DO
    DEALLOCATE(data_manbo_temp, stat=allocate_status)
      IF(allocate_status/=0) THEN
        PRINT *, "No memory enough to deallocate. ManBo cannot run."
        CALL log_write("ERROR: Error on deallocate data_manbo_temp on manbo_forces.F90")
        CALL log_close(1)
        STOP
      END IF
    ! Here we reallocated data_manbo_out with the right size
  
    n_mols_out = MAXVAL(data_manbo_out%mol)
    IF (group_monomers>1) n_mols_out_orig = MAXVAL(data_manbo_out%mol_orig)
    n_atoms_out = SIZE(data_manbo_out, DIM=1)
    k = n_mols_out - n_mols
      IF(ALLOCATED(mc_out)) THEN
        DEALLOCATE(mc_out, stat=allocate_status)
          IF(allocate_status/=0) THEN
            PRINT *, "No memory enough to deallocate. ManBo cannot run."
            CALL log_write("ERROR: Error on deallocate mc_out on manbo_forces.F90")
            CALL log_close(1)
            STOP
          END IF
      END IF
      ALLOCATE(mc_out(k), stat=allocate_status)
        IF(allocate_status/=0) THEN
          PRINT *, "No memory enough to allocate. ManBo cannot run."
          CALL log_write("ERROR: Error on allocate mc_out on manbo_forces.F90")
          CALL log_close(1)
          STOP
        END IF
    IF (group_monomers>1) THEN
      k2 = n_mols_out_orig - n_mols_orig
        IF(ALLOCATED(mc_out_orig)) THEN
          DEALLOCATE(mc_out_orig, stat=allocate_status)
            IF(allocate_status/=0) THEN
              PRINT *, "No memory enough to deallocate. ManBo can not run."
              CALL log_write("ERROR: Error on deallocate mc_out_orig on manbo_forces.F90")
              CALL log_close(1)
              STOP
            END IF
        END IF
        ALLOCATE(mc_out_orig(k2), stat=allocate_status)
          IF(allocate_status/=0) THEN
            PRINT *, "No memory enough to allocate. ManBo can not run."
            CALL log_write("ERROR: Error on allocate mc_out_orig on manbo_forces.F90")
            CALL log_close(1)
            STOP
          END IF
    END IF
    
    CALL calculate_center_of_mass(2)
    ! Here we calculate the mass center's of the outside box's molecules
    IF (group_monomers>1) CALL calculate_center_of_mass(4)
    ! Here we calculate the mass center's of the outside box's molecules using the original numbering
  END SUBROUTINE box_replication

  SUBROUTINE calculate_forces(mpi_id)
  ! This subroutine calculates the energy and the forces acting in each atom of the system
    IMPLICIT NONE
    INTEGER, INTENT(in) :: mpi_id
    INTEGER :: i, num, mpi_size
    INTEGER :: n_qm_procs_old, qm_prog_memory_old
    INTEGER, DIMENSION(2) :: char_mul_old
    DOUBLE PRECISION :: A,B,C,D,E1_mbe2,E2_mbe2,E3_mbe2,E1_hfmbe2,E2_hfmbe2,E3_hfmbe2,&
                 E1_mbe3,E2_mbe3,E3_mbe3,E1_hfmbe3,E2_hfmbe3,E3_hfmbe3,E1_hffull2,&
                 E1_hffull3
    CHARACTER(LEN=100) :: line
    CHARACTER(LEN=70) :: method_old
    TYPE temp
      DOUBLE PRECISION, DIMENSION(3) :: f_mbe2,f_hfmbe2,f_mbe3,f_hfmbe3,f_hffull2,f_hffull3
    END TYPE temp
    TYPE(temp), DIMENSION(n_atoms) :: data_temp
    LOGICAL :: use_embedding_old
#ifdef USE_PARALLEL
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
#endif
    
    IF (mpi_id == 0) THEN
      IF(ALLOCATED(data_manbo_out)) THEN
        CALL box_replication
      ELSE
        CALL calculate_center_of_mass(1)
        CALL log_write("  None of the atoms were replicated.")
        CALL log_write("")
        n_mols_out = n_mols
        n_atoms_out = 0
      END IF
    END IF
    
#ifdef USE_PARALLEL
    IF (mpi_id == 0) THEN
      i = system("mkdir -p tempgaufiles_" // TRIM(ADJUSTL(name_out)))
      IF(i/=0) THEN
        PRINT *, "ManBo cannot run. Error to create a temporary folder using the following command:"
        PRINT *, "  mkdir -p tempgaufiles_" // TRIM(ADJUSTL(name_out))
        CALL log_write("ERROR: Error to create a temporary folder using the following command:")
        CALL log_write("       mkdir -p tempgaufiles_" // TRIM(ADJUSTL(name_out)))
        CALL log_close(1)
        STOP
      END IF
    END IF
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif
  
#ifdef USE_PARALLEL
    CALL MPI_BCAST(bs_extrapolation, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST(mbe_corr_only, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size, ierr)
#else
    mpi_size = 1
#endif /* USE_PARALLEL */
    
    IF(.NOT. bs_extrapolation) THEN
      IF(.NOT. mbe_corr_only) THEN
        CALL forces_eval(mpi_id,0)
      ELSE
        IF (mpi_id == 0) THEN
          data_temp%f_mbe2(1) = 0.0
          data_temp%f_mbe2(2) = 0.0
          data_temp%f_mbe2(3) = 0.0
          data_temp%f_hfmbe2(1) = 0.0
          data_temp%f_hfmbe2(2) = 0.0
          data_temp%f_hfmbe2(3) = 0.0
        END IF
        
        CALL forces_eval(mpi_id,0)
        IF (mpi_id == 0) THEN
          WRITE(line,'(3x,"Full EE-MBE Energy done with",a12,":                  ",f16.8," Hartrees")')&
                TRIM(ADJUSTL(qm_prog_method)), (E1+E2+E3)
          CALL log_write(TRIM(line))
          E1_mbe2 = E1
          E2_mbe2 = E2
          E3_mbe2 = E3
          data_temp%f_mbe2(1) = data_manbo%f(1)
          data_temp%f_mbe2(2) = data_manbo%f(2)
          data_temp%f_mbe2(3) = data_manbo%f(3)
          
          method_old = qm_prog_method
          qm_prog_method = "HF"
        END IF
        CALL forces_eval(mpi_id,1)
        IF (mpi_id == 0) THEN
          WRITE(line,'(3x,"Full EE-MBE Energy done with HF:                           ",f16.8," Hartrees")') (E1+E2+E3)
          CALL log_write(TRIM(line))
          E1_hfmbe2 = E1
          E2_hfmbe2 = E2
          E3_hfmbe2 = E3
          data_temp%f_hfmbe2(1) = data_manbo%f(1)
          data_temp%f_hfmbe2(2) = data_manbo%f(2)
          data_temp%f_hfmbe2(3) = data_manbo%f(3)
          
          char_mul_old(1) = char_mul(1)%q_mol
          char_mul_old(2) = char_mul(1)%mul
          DO i=2,n_mols
            char_mul(1)%q_mol = char_mul(1)%q_mol + char_mul(i)%q_mol
            char_mul(1)%mul = char_mul(1)%mul + char_mul(i)%mul-1
          END DO
          n_qm_procs_old = n_qm_procs
          qm_prog_memory_old = qm_prog_memory
#ifdef USE_PARALLEL
          n_qm_procs = omp_get_max_threads()*n_qm_procs
          qm_prog_memory = omp_get_max_threads()*qm_prog_memory
#elif USE_OPENMPONLY
          n_qm_procs = omp_get_max_threads()*n_qm_procs
          qm_prog_memory = omp_get_max_threads()*qm_prog_memory
#endif
          use_embedding_old = use_embedding
          use_embedding = .FALSE.
          CALL forces_eval_full(0)
          WRITE(line,'(3x,"Full HF Energy (no EE-MBE):                                ",f16.8," Hartrees")') E1
          CALL log_write(TRIM(line))
          
          data_manbo%f(1) = data_manbo%f(1) + data_temp%f_mbe2(1) - data_temp%f_hfmbe2(1)
          data_manbo%f(2) = data_manbo%f(2) + data_temp%f_mbe2(2) - data_temp%f_hfmbe2(2)
          data_manbo%f(3) = data_manbo%f(3) + data_temp%f_mbe2(3) - data_temp%f_hfmbe2(3)
          E1 = E1 + E1_mbe2 - E1_hfmbe2
          E2 = E2 + E2_mbe2 - E2_hfmbe2
          E3 = E3 + E3_mbe2 - E3_hfmbe2
          WRITE(line,'(3x,a," Energy with EE-MBE on the correlation only: ",f16.8," Hartrees")')&
                TRIM(ADJUSTL(method_old)), (E1+E2+E3)
          CALL log_write(TRIM(line))
          
          ! Restoring the original values for the next md_step
          qm_prog_method = method_old
          char_mul(1)%q_mol = char_mul_old(1)
          char_mul(1)%mul   = char_mul_old(2)
          n_qm_procs = n_qm_procs_old
          qm_prog_memory = qm_prog_memory_old
          use_embedding = use_embedding_old
        END IF
      END IF
    ELSE
      IF(.NOT. mbe_corr_only) THEN
        IF (mpi_id == 0) THEN
          data_temp%f_mbe2(1) = 0.0
          data_temp%f_mbe2(2) = 0.0
          data_temp%f_mbe2(3) = 0.0
          data_temp%f_hfmbe2(1) = 0.0
          data_temp%f_hfmbe2(2) = 0.0
          data_temp%f_hfmbe2(3) = 0.0
          data_temp%f_mbe3(1) = 0.0
          data_temp%f_mbe3(2) = 0.0
          data_temp%f_mbe3(3) = 0.0
          data_temp%f_hfmbe3(1) = 0.0
          data_temp%f_hfmbe3(2) = 0.0
          data_temp%f_hfmbe3(3) = 0.0
          
          qm_prog_method = "MP2"
          qm_prog_basis = "cc-pVTZ"
        END IF
        CALL forces_eval(mpi_id,0)
        IF (mpi_id == 0) THEN
          WRITE(line,'(3x,"Full EE-MBE Energy computed at MP2/cc-pVTZ:                ",f16.8," Hartrees")') (E1+E2+E3)
          CALL log_write(TRIM(line))
          E1_mbe3 = E1
          E2_mbe3 = E2
          E3_mbe3 = E3
          data_temp%f_mbe3(1) = data_manbo%f(1)
          data_temp%f_mbe3(2) = data_manbo%f(2)
          data_temp%f_mbe3(3) = data_manbo%f(3)
          
          qm_prog_method = "MP2"
          qm_prog_basis = "cc-pVDZ"
        END IF
        CALL forces_eval(mpi_id,1)
        IF (mpi_id == 0) THEN
          WRITE(line,'(3x,"Full EE-MBE Energy computed at MP2/cc-pVDZ:                ",f16.8," Hartrees")') (E1+E2+E3)
          CALL log_write(TRIM(line))
          E1_mbe2 = E1
          E2_mbe2 = E2
          E3_mbe2 = E3
          data_temp%f_mbe2(1) = data_manbo%f(1)
          data_temp%f_mbe2(2) = data_manbo%f(2)
          data_temp%f_mbe2(3) = data_manbo%f(3)
          
          qm_prog_method = "HF"
          qm_prog_basis = "cc-pVTZ"
        END IF
        CALL forces_eval(mpi_id,1)
        IF (mpi_id == 0) THEN
          WRITE(line,'(3x,"Full EE-MBE Energy computed at HF/cc-pVTZ:                 ",f16.8," Hartrees")') (E1+E2+E3)
          CALL log_write(TRIM(line))
          E1_hfmbe3 = E1
          E2_hfmbe3 = E2
          E3_hfmbe3 = E3
          data_temp%f_hfmbe3(1) = data_manbo%f(1)
          data_temp%f_hfmbe3(2) = data_manbo%f(2)
          data_temp%f_hfmbe3(3) = data_manbo%f(3)
          
          qm_prog_method = "HF"
          qm_prog_basis = "cc-pVDZ"
        END IF
        CALL forces_eval(mpi_id,1)
        IF (mpi_id == 0) THEN
          WRITE(line,'(3x,"Full EE-MBE Energy computed at HF/cc-pVDZ:                 ",f16.8," Hartrees")') (E1+E2+E3)
          CALL log_write(TRIM(line))
          E1_hfmbe2 = E1
          E2_hfmbe2 = E2
          E3_hfmbe2 = E3
          data_temp%f_hfmbe2(1) = data_manbo%f(1)
          data_temp%f_hfmbe2(2) = data_manbo%f(2)
          data_temp%f_hfmbe2(3) = data_manbo%f(3)
          
          A = 3.0**3.4/(3.0**3.4-2.0**3.4)
          B = 2.0**3.4/(3.0**3.4-2.0**3.4)
          C = 3.0**2.2/(3.0**2.2-2.0**2.2)
          D = 2.0**2.2/(3.0**2.2-2.0**2.2)
          data_manbo%f(1) = A*data_temp%f_hfmbe3(1)-B*data_temp%f_hfmbe2(1)+C*(data_temp%f_mbe3(1)-data_temp%f_hfmbe3(1))
          data_manbo%f(2) = A*data_temp%f_hfmbe3(2)-B*data_temp%f_hfmbe2(2)+C*(data_temp%f_mbe3(2)-data_temp%f_hfmbe3(2))
          data_manbo%f(3) = A*data_temp%f_hfmbe3(3)-B*data_temp%f_hfmbe2(3)+C*(data_temp%f_mbe3(3)-data_temp%f_hfmbe3(3))
          data_manbo%f(1) = data_manbo%f(1)-D*(data_temp%f_mbe2(1)-data_temp%f_hfmbe2(1))
          data_manbo%f(2) = data_manbo%f(2)-D*(data_temp%f_mbe2(2)-data_temp%f_hfmbe2(2))
          data_manbo%f(3) = data_manbo%f(3)-D*(data_temp%f_mbe2(3)-data_temp%f_hfmbe2(3))
          E1 = A*E1_hfmbe3-B*E1_hfmbe2+C*(E1_mbe3-E1_hfmbe3)-D*(E1_mbe2-E1_hfmbe2)
          E2 = A*E2_hfmbe3-B*E2_hfmbe2+C*(E2_mbe3-E2_hfmbe3)-D*(E2_mbe2-E2_hfmbe2)
          E3 = A*E3_hfmbe3-B*E3_hfmbe2+C*(E3_mbe3-E3_hfmbe3)-D*(E3_mbe2-E3_hfmbe2)
          WRITE(line,'(3x,"HF  Full EE-MBE Energy extrapolated to infinity basis set: ",f16.8," Hartrees")')&
          (A*(E1_hfmbe3+E2_hfmbe3+E3_hfmbe3)-B*(E1_hfmbe2+E2_hfmbe2+E3_hfmbe2))
          CALL log_write(TRIM(line))
          WRITE(line,'(3x,"MP2 Full EE-MBE Energy extrapolated to infinity basis set: ",f16.8," Hartrees")')&
          (E1+E2+E3)
          CALL log_write(TRIM(line))
        END IF
      ELSE
        IF (mpi_id == 0) THEN
          data_temp%f_mbe2(1) = 0.0
          data_temp%f_mbe2(2) = 0.0
          data_temp%f_mbe2(3) = 0.0
          data_temp%f_hfmbe2(1) = 0.0
          data_temp%f_hfmbe2(2) = 0.0
          data_temp%f_hfmbe2(3) = 0.0
          data_temp%f_mbe3(1) = 0.0
          data_temp%f_mbe3(2) = 0.0
          data_temp%f_mbe3(3) = 0.0
          data_temp%f_hfmbe3(1) = 0.0
          data_temp%f_hfmbe3(2) = 0.0
          data_temp%f_hfmbe3(3) = 0.0
          data_temp%f_hffull2(1) = 0.0
          data_temp%f_hffull2(2) = 0.0
          data_temp%f_hffull2(3) = 0.0
          data_temp%f_hffull3(1) = 0.0
          data_temp%f_hffull3(2) = 0.0
          data_temp%f_hffull3(3) = 0.0
          
          qm_prog_method = "MP2"
          qm_prog_basis = "cc-pVTZ"
        END IF
        CALL forces_eval(mpi_id,0)
        IF (mpi_id == 0) THEN
          WRITE(line,'(3x,"Full EE-MBE Energy computed at MP2/cc-pVTZ:                  ",f16.8," Hartrees")') (E1+E2+E3)
          CALL log_write(TRIM(line))
          E1_mbe3 = E1
          E2_mbe3 = E2
          E3_mbe3 = E3
          data_temp%f_mbe3(1) = data_manbo%f(1)
          data_temp%f_mbe3(2) = data_manbo%f(2)
          data_temp%f_mbe3(3) = data_manbo%f(3)
          
          qm_prog_method = "MP2"
          qm_prog_basis = "cc-pVDZ"
        END IF
        CALL forces_eval(mpi_id,1)
        IF (mpi_id == 0) THEN
          WRITE(line,'(3x,"Full EE-MBE Energy computed at MP2/cc-pVDZ:                  ",f16.8," Hartrees")') (E1+E2+E3)
          CALL log_write(TRIM(line))
          E1_mbe2 = E1
          E2_mbe2 = E2
          E3_mbe2 = E3
          data_temp%f_mbe2(1) = data_manbo%f(1)
          data_temp%f_mbe2(2) = data_manbo%f(2)
          data_temp%f_mbe2(3) = data_manbo%f(3)
          
          qm_prog_method = "HF"
          qm_prog_basis = "cc-pVTZ"
        END IF
        CALL forces_eval(mpi_id,1)
        IF (mpi_id == 0) THEN
          WRITE(line,'(3x,"Full EE-MBE Energy computed at HF/cc-pVTZ:                   ",f16.8," Hartrees")') (E1+E2+E3)
          CALL log_write(TRIM(line))
          E1_hfmbe3 = E1
          E2_hfmbe3 = E2
          E3_hfmbe3 = E3
          data_temp%f_hfmbe3(1) = data_manbo%f(1)
          data_temp%f_hfmbe3(2) = data_manbo%f(2)
          data_temp%f_hfmbe3(3) = data_manbo%f(3)
          
          qm_prog_method = "HF"
          qm_prog_basis = "cc-pVDZ"
        END IF
        CALL forces_eval(mpi_id,1)
        IF (mpi_id == 0) THEN
          WRITE(line,'(3x,"Full EE-MBE Energy computed at HF/cc-pVDZ:                   ",f16.8," Hartrees")') (E1+E2+E3)
          CALL log_write(TRIM(line))
          E1_hfmbe2 = E1
          E2_hfmbe2 = E2
          E3_hfmbe2 = E3
          data_temp%f_hfmbe2(1) = data_manbo%f(1)
          data_temp%f_hfmbe2(2) = data_manbo%f(2)
          data_temp%f_hfmbe2(3) = data_manbo%f(3)
          
          char_mul_old(1) = char_mul(1)%q_mol
          char_mul_old(2) = char_mul(1)%mul
          DO i=2,n_mols
            char_mul(1)%q_mol = char_mul(1)%q_mol + char_mul(i)%q_mol
            char_mul(1)%mul = char_mul(1)%mul + char_mul(i)%mul-1
          END DO
          n_qm_procs_old = n_qm_procs
          qm_prog_memory_old = qm_prog_memory
#ifdef USE_PARALLEL
          n_qm_procs = omp_get_max_threads()*n_qm_procs
          qm_prog_memory = omp_get_max_threads()*qm_prog_memory
#elif USE_OPENMPONLY
          n_qm_procs = omp_get_max_threads()*n_qm_procs
          qm_prog_memory = omp_get_max_threads()*qm_prog_memory
#endif
          use_embedding_old = use_embedding
          use_embedding = .FALSE.
        END IF
        
        IF (mpi_size > 1) THEN
#ifdef USE_PARALLEL
          CALL MPI_BCAST(n_mols,            1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
          CALL MPI_BCAST(char_mul(1)%q_mol, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
          CALL MPI_BCAST(char_mul(1)%mul,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
          CALL MPI_BCAST(n_qm_procs,        1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
          DO i=1,n_atoms
            CALL MPI_BCAST(data_manbo(i)%an,  1, MPI_INTEGER,    0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(data_manbo(i)%r,   3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(data_manbo(i)%q,   1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(data_manbo(i)%mol, 1, MPI_INTEGER,    0, MPI_COMM_WORLD, ierr)
          END DO
#endif
        END IF
        
        IF (mpi_id == 0) THEN
          qm_prog_method = "HF"
          qm_prog_basis = "cc-pVTZ"
          CALL forces_eval_full(0)
          WRITE(line,'(3x,"Full Energy (no EE-MBE) computed at HF/cc-pVTZ:              ",f16.8," Hartrees")') E1
          CALL log_write(TRIM(line))
          E1_hffull3 = E1
          data_temp%f_hffull3(1) = data_manbo%f(1)
          data_temp%f_hffull3(2) = data_manbo%f(2)
          data_temp%f_hffull3(3) = data_manbo%f(3)
        END IF
          
        IF (mpi_size == 1 .OR. mpi_id == 1) THEN
          use_embedding = .FALSE.
          qm_prog_method = "HF"
          qm_prog_basis = "cc-pVDZ"
          CALL forces_eval_full(1)
          E1_hffull2 = E1
          data_temp%f_hffull2(1) = data_manbo%f(1)
          data_temp%f_hffull2(2) = data_manbo%f(2)
          data_temp%f_hffull2(3) = data_manbo%f(3)
          IF (mpi_id == 1) THEN
#ifdef USE_PARALLEL
            CALL MPI_SEND(E1_hffull2, 1, MPI_DOUBLE_PRECISION, 0, 11, MPI_COMM_WORLD, ierr)
            CALL MPI_SEND(data_temp%f_hffull2(1), n_atoms*1, MPI_DOUBLE_PRECISION, 0, 12, MPI_COMM_WORLD, ierr)
            CALL MPI_SEND(data_temp%f_hffull2(2), n_atoms*1, MPI_DOUBLE_PRECISION, 0, 12, MPI_COMM_WORLD, ierr)
            CALL MPI_SEND(data_temp%f_hffull2(3), n_atoms*1, MPI_DOUBLE_PRECISION, 0, 12, MPI_COMM_WORLD, ierr)
#endif
          END IF
        END IF
        IF (mpi_size > 1 .AND. mpi_id == 0) THEN
#ifdef USE_PARALLEL
          CALL MPI_RECV(E1_hffull2, 1, MPI_DOUBLE_PRECISION, 1, 11, MPI_COMM_WORLD, status, ierr)
          CALL MPI_RECV(data_temp%f_hffull2(1), n_atoms*1, MPI_DOUBLE_PRECISION, 1, 12, MPI_COMM_WORLD, status, ierr)
          CALL MPI_RECV(data_temp%f_hffull2(2), n_atoms*1, MPI_DOUBLE_PRECISION, 1, 12, MPI_COMM_WORLD, status, ierr)
          CALL MPI_RECV(data_temp%f_hffull2(3), n_atoms*1, MPI_DOUBLE_PRECISION, 1, 12, MPI_COMM_WORLD, status, ierr)
#endif
        END IF
        IF (mpi_id == 0) THEN
          WRITE(line,'(3x,"Full Energy (no EE-MBE) computed at HF/cc-pVDZ:              ",f16.8," Hartrees")') E1_hffull2
          CALL log_write(TRIM(line))
        END IF
          
        IF (mpi_id == 0) THEN
          A = 3.0**3.4/(3.0**3.4-2.0**3.4)
          B = 2.0**3.4/(3.0**3.4-2.0**3.4)
          C = 3.0**2.2/(3.0**2.2-2.0**2.2)
          D = 2.0**2.2/(3.0**2.2-2.0**2.2)
          WRITE(line,'(3x,"HF Full EE-MBE Energy extrapolated to infinity basis set:    ",f16.8," Hartrees")')&
          (A*(E1_hfmbe3+E2_hfmbe3+E3_hfmbe3)-B*(E1_hfmbe2+E2_hfmbe2+E3_hfmbe2))
          CALL log_write(TRIM(line))
          WRITE(line,'(3x,"HF Full Energy (no EE-MBE) extrapolated to infinity b.s.:    ",f16.8," Hartrees")')&
          (A*(E1_hffull3)-B*(E1_hffull2))
          CALL log_write(TRIM(line))
          data_manbo%f(1) = A*data_temp%f_hffull3(1)-B*data_temp%f_hffull2(1)+C*(data_temp%f_mbe3(1)-data_temp%f_hfmbe3(1))
          data_manbo%f(2) = A*data_temp%f_hffull3(2)-B*data_temp%f_hffull2(2)+C*(data_temp%f_mbe3(2)-data_temp%f_hfmbe3(2))
          data_manbo%f(3) = A*data_temp%f_hffull3(3)-B*data_temp%f_hffull2(3)+C*(data_temp%f_mbe3(3)-data_temp%f_hfmbe3(3))
          data_manbo%f(1) = data_manbo%f(1)-D*(data_temp%f_mbe2(1)-data_temp%f_hfmbe2(1))
          data_manbo%f(2) = data_manbo%f(2)-D*(data_temp%f_mbe2(2)-data_temp%f_hfmbe2(2))
          data_manbo%f(3) = data_manbo%f(3)-D*(data_temp%f_mbe2(3)-data_temp%f_hfmbe2(3))
          E1 = A*E1_hffull3-B*E1_hffull2+C*(E1_mbe3-E1_hfmbe3)-D*(E1_mbe2-E1_hfmbe2)
          E2 = C*(E2_mbe3-E2_hfmbe3)-D*(E2_mbe2-E2_hfmbe2)
          E3 = C*(E3_mbe3-E3_hfmbe3)-D*(E3_mbe2-E3_hfmbe2)
          WRITE(line,'(3x,"MP2 extrapolated Energy with EE-MBE on the correlation only: ",f16.8," Hartrees")')&
          (E1+E2+E3)
          CALL log_write(TRIM(line))
          
          ! Restoring the original values for the next md_step
          char_mul(1)%q_mol = char_mul_old(1)
          char_mul(1)%mul   = char_mul_old(2)
          n_qm_procs = n_qm_procs_old
          qm_prog_memory = qm_prog_memory_old
          use_embedding = use_embedding_old
        END IF
      END IF
    END IF
    
#ifdef USE_PARALLEL
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
    IF (mpi_id == 0) THEN
      num = system("rm -R tempgaufiles_" // TRIM(ADJUSTL(name_out)))
      IF(num/=0) THEN
        PRINT *, "ManBo cannot run. Error to delete a temporary folder using the following command:"
        PRINT *, "  rm -R tempgaufiles_" // TRIM(ADJUSTL(name_out))
        CALL log_write("ERROR: Error to delete a temporary folder using the following command:")
        CALL log_write("       rm -R tempgaufiles_" // TRIM(ADJUSTL(name_out)))
        CALL log_close(1)
        STOP
      END IF
    END IF
#else
    num = system("rm -R tempgaufiles_" // TRIM(ADJUSTL(name_out)))
    IF(num/=0) THEN
      PRINT *, "ManBo cannot run. Error to delete a temporary folder using the following command:"
      PRINT *, "  rm -R tempgaufiles_" // TRIM(ADJUSTL(name_out))
      CALL log_write("ERROR: Error to delete a temporary folder using the following command:")
      CALL log_write("       rm -R tempgaufiles_" // TRIM(ADJUSTL(name_out)))
      CALL log_close(1)
      STOP
    END IF
#endif
  END SUBROUTINE calculate_forces
  
  SUBROUTINE forces_eval(mpi_id,p)
  ! This subroutine evaluates the forces and energies
    IMPLICIT NONE
    INTEGER, INTENT(in) :: mpi_id, p
    INTEGER :: allocate_status
    INTEGER :: t1, t2, mpi_size, cnt
    INTEGER :: i, j, k, num, vp_dim, tnpairs, npairs, ntrimers, local
    INTEGER, DIMENSION(3) :: list
    INTEGER, DIMENSION(:,:), ALLOCATABLE, SAVE :: pairs_list, trimers_list
    CHARACTER(LEN=100) :: line
    LOGICAL :: finish
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE :: U, Up
    ! U contains the energy of a monomer, a dimer, or trimer given in  hartree (atomic unit)
    TYPE(forces), DIMENSION(:,:), ALLOCATABLE, SAVE :: force, forcep
    ! force contains the force of a monomer, a dimer, or trimer
    
    IF (mpi_id == 0) THEN
      IF (n_atoms_out > 0) THEN
        vp_dim = NINT(2.0*((1.0+2.0*cutoff/MINVAL(box_dim))**3)*n_mols*((2*cutoff/MINVAL(box_dim))**3))
      ELSE
        vp_dim = NINT(2.0*n_mols*((2*cutoff/MINVAL(box_dim))**3))
      END IF
      ! This is an estimate for allocate val_pairs
    END IF
    
#ifdef USE_PARALLEL
    CALL MPI_BCAST(vp_dim, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST(n_mols, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif
    ALLOCATE(pairs_list(n_mols*vp_dim,2), trimers_list(NINT(FLOAT(n_mols*vp_dim*(n_mols-2))/3.),3),&
             stat=allocate_status)
    IF(allocate_status/=0) THEN
      PRINT *, "No memory enough to allocate. ManBo cannot run."
      CALL log_write("")
      CALL log_write(" ERROR: Error on allocate pairs_list and trimers_list on manbo_forces.F90")
      CALL log_close(1)
      STOP
    END IF
    
    pairs_list = 0
    trimers_list = 0
    ! Setting all values to zero
    
    IF (mpi_id == 0) THEN
      IF (p==0) CALL log_write("  Number of monomers, dimers and trimers used in the calculation of forces:")
      
      WRITE(line,'(3x,i8," monomer(s)")') n_mols
      IF (p==0) CALL log_write(TRIM(line))
      
      ! Creating list of dimers to be computed
      IF (expan_order==2 .OR. expan_order==3) THEN
        list = 0
        num = 0
        DO i=1,n_mols
          k = 0
          DO j=i+1,n_mols_out
            IF (valid_pair(i,j)) THEN
              k = k + 1
              num = num + 1
              pairs_list(num,1) = i
              pairs_list(num,2) = j
            END IF
          END DO
        END DO
        npairs = num
        DO i=n_mols+1,n_mols_out
          k = 0
          DO j=i+1,n_mols_out
            IF (valid_pair(i,j)) THEN
              k = k + 1
              num = num + 1
              pairs_list(num,1) = i
              pairs_list(num,2) = j
            END IF
          END DO
        END DO
        tnpairs = num
      ELSE
        npairs = 0
      END IF
      WRITE(line,'(3x,i8," dimer(s)")') npairs
      IF (p==0) CALL log_write(TRIM(line))
      
      ! Creating list of trimers to be computed
      IF(expan_order==3) THEN        
        num = 0
        DO i=1,npairs-1
          DO j=i+1,npairs
           IF (pairs_list(i,1)==pairs_list(j,1)) THEN
             DO k=j+1,tnpairs
               IF (pairs_list(k,1)==pairs_list(i,2) .AND. pairs_list(k,2)==pairs_list(j,2)) THEN
                 num = num + 1
                 trimers_list(num,1) = pairs_list(i,1)
                 trimers_list(num,2) = pairs_list(i,2)
                 trimers_list(num,3) = pairs_list(j,2)
                 EXIT
               END IF
             END DO
           END IF
          END DO
        END DO
        ntrimers = num
      ELSE
        ntrimers = 0
      END IF
      WRITE(line,'(3x,i8," trimer(s)")') ntrimers
      IF (p==0) CALL log_write(TRIM(line))
    END IF
    
#ifdef USE_PARALLEL
    CALL MPI_BCAST(npairs, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST(ntrimers, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    IF (npairs > 0) THEN
      DO i=1,npairs
        CALL MPI_BCAST(pairs_list(i,1), 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(pairs_list(i,2), 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      END DO
    END IF
    IF (ntrimers > 0) THEN
      DO i=1,ntrimers
        CALL MPI_BCAST(trimers_list(i,1), 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(trimers_list(i,2), 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(trimers_list(i,3), 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      END DO
    END IF
    CALL MPI_BCAST(n_qm_procs, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST(qm_prog, 70, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST(qm_prog_memory, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST(qm_prog_method, 70, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST(qm_prog_basis, 70, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    DO i=1,n_mols
      CALL MPI_BCAST(char_mul(i)%q_mol, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(char_mul(i)%mul,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    END DO
    CALL MPI_BCAST(use_emb_radius, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST(use_embedding,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    DO k=1,n_atoms
      CALL MPI_BCAST(data_manbo(k)%an,  1, MPI_INTEGER,    0, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(data_manbo(k)%r,   3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(data_manbo(k)%q,   1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(data_manbo(k)%mol, 1, MPI_INTEGER,    0, MPI_COMM_WORLD, ierr)
    END DO
    CALL MPI_BCAST(n_mols_out,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST(n_atoms_out, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    DO k=1,n_atoms_out
      CALL MPI_BCAST(data_manbo_out(k)%an,   1, MPI_INTEGER,    0, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(data_manbo_out(k)%r,    3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(data_manbo_out(k)%atom, 1, MPI_INTEGER,    0, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(data_manbo_out(k)%mol,  1, MPI_INTEGER,    0, MPI_COMM_WORLD, ierr)
    END DO
    CALL MPI_BCAST(emb_radius, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST(box_dim,    3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    IF (ALLOCATED(orig_repli_mol)) THEN
      num = SIZE(orig_repli_mol)
      CALL MPI_BCAST(orig_repli_mol, num, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    END IF
    DO i=1,n_mols
      CALL MPI_BCAST(mc(i)%r,    3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(mc(i)%mass, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    END DO
    IF (n_mols_out>n_mols .AND. mpi_id>0) THEN    
      k = n_mols_out - n_mols
        IF(ALLOCATED(mc_out)) THEN
          DEALLOCATE(mc_out)
          ALLOCATE(mc_out(k))
        ELSE
          ALLOCATE(mc_out(k))
        END IF
    END IF
    DO i=1,n_mols_out-n_mols
      CALL MPI_BCAST(mc_out(i)%r,    3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(mc_out(i)%mass, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    END DO
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size, ierr)
#else
    mpi_size = 1
#endif /* USE_PARALLEL */
    
    ALLOCATE(forcep(n_atoms,1+npairs+ntrimers),Up(n_mols+npairs+ntrimers),&
             stat=allocate_status)
    IF(allocate_status/=0) THEN
      PRINT *, "No memory enough to allocate. ManBo cannot run."
      CALL log_write("")
      CALL log_write(" ERROR: Error on allocate forcep and Up on manbo_forces.F90")
      CALL log_close(1)
      STOP
    END IF

    forcep%f(1) = 0.0
    forcep%f(2) = 0.0
    forcep%f(3) = 0.0
    Up = 0.0
    
    finish = .FALSE.
    cnt = 0
    
    !$omp parallel shared( cnt,finish ) firstprivate ( list,local )
    DO WHILE(.NOT. finish)
      !$omp critical
      local = mpi_id+cnt*mpi_size+1
      cnt = cnt + 1
      IF (local .ge. (n_mols+npairs+ntrimers)) THEN
        finish = .TRUE.
      END IF
      !$omp end critical

      ! Computing monomers
      IF (local .le. n_mols) THEN
        list(1) = local
        CALL call_qm_prog(forcep,Up,n_mols+npairs+ntrimers,1+npairs+ntrimers,local,list,1)
      ! Computing dimers
      ELSE IF (local .gt. n_mols .AND. local .le. (n_mols+npairs)) THEN
        list(1) = pairs_list(local-n_mols,1)
        list(2) = pairs_list(local-n_mols,2)
        CALL call_qm_prog(forcep,Up,n_mols+npairs+ntrimers,1+npairs+ntrimers,local,list,2)
      ! Computing trimers
      ELSE IF (local .gt. (n_mols+npairs) .AND. local .le. (n_mols+npairs+ntrimers)) THEN
        list(1) = trimers_list(local-n_mols-npairs,1)
        list(2) = trimers_list(local-n_mols-npairs,2)
        list(3) = trimers_list(local-n_mols-npairs,3)
        CALL call_qm_prog(forcep,Up,n_mols+npairs+ntrimers,1+npairs+ntrimers,local,list,3)
      END IF
    END DO
    !$omp end parallel
    
#ifdef USE_PARALLEL
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif

    IF (mpi_id == 0) THEN
      ALLOCATE(force(n_atoms,1+npairs+ntrimers),U(n_mols+npairs+ntrimers),&
               stat=allocate_status)
      IF(allocate_status/=0) THEN
        PRINT *, "No memory enough to allocate. ManBo cannot run."
        CALL log_write("")
        CALL log_write(" ERROR: Error on allocate force and U on manbo_forces.F90")
        CALL log_close(1)
        STOP
      END IF
      force%f(1) = 0.0
      force%f(2) = 0.0
      force%f(3) = 0.0
      U = 0.0
    END IF
    
#ifdef USE_PARALLEL
    CALL MPI_REDUCE(Up, U, n_mols+npairs+ntrimers, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    DO k=1,n_atoms
      DO i=1,1+npairs+ntrimers
        CALL MPI_REDUCE(forcep(k,i)%f(1), force(k,i)%f(1), 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_REDUCE(forcep(k,i)%f(2), force(k,i)%f(2), 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_REDUCE(forcep(k,i)%f(3), force(k,i)%f(3), 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      END DO
    END DO
#else
    force%f(1) = forcep%f(1)
    force%f(2) = forcep%f(2)
    force%f(3) = forcep%f(3)
    U = Up
#endif /* USE_PARALLEL */
      
    DEALLOCATE(forcep, Up)
  
    IF (mpi_id == 0) THEN
      ! Computing forces
      data_manbo%f(1) = 0.0
      data_manbo%f(2) = 0.0
      data_manbo%f(3) = 0.0
      DO k=1,n_atoms
        data_manbo(k)%f = data_manbo(k)%f + force(k,1)%f
        ! These terms come from the MBE
        DO i=1,npairs
          IF (pairs_list(i,1)==data_manbo(k)%mol .OR. pairs_list(i,2)==data_manbo(k)%mol) THEN
            data_manbo(k)%f = data_manbo(k)%f + force(k,1+i)%f - force(k,1)%f
          END IF
          ! These terms come from the MBE
        END DO
        DO i=1,ntrimers
          IF (trimers_list(i,1)==data_manbo(k)%mol) THEN
            data_manbo(k)%f = data_manbo(k)%f + force(k,1+npairs+i)%f
            CALL pair_num(pairs_list,SIZE(pairs_list,DIM=1),trimers_list(i,1),trimers_list(i,2),j)
            data_manbo(k)%f = data_manbo(k)%f - force(k,1+j)%f
            CALL pair_num(pairs_list,SIZE(pairs_list,DIM=1),trimers_list(i,1),trimers_list(i,3),j)
            data_manbo(k)%f = data_manbo(k)%f - force(k,1+j)%f + force(k,1)%f
          ELSE IF (trimers_list(i,2)==data_manbo(k)%mol) THEN
            data_manbo(k)%f = data_manbo(k)%f + force(k,1+npairs+i)%f
            CALL pair_num(pairs_list,SIZE(pairs_list,DIM=1),trimers_list(i,1),trimers_list(i,2),j)
            data_manbo(k)%f = data_manbo(k)%f - force(k,1+j)%f
            CALL pair_num(pairs_list,SIZE(pairs_list,DIM=1),trimers_list(i,2),trimers_list(i,3),j)
            data_manbo(k)%f = data_manbo(k)%f - force(k,1+j)%f + force(k,1)%f
          ELSE IF (trimers_list(i,3)==data_manbo(k)%mol) THEN
            data_manbo(k)%f = data_manbo(k)%f + force(k,1+npairs+i)%f
            CALL pair_num(pairs_list,SIZE(pairs_list,DIM=1),trimers_list(i,1),trimers_list(i,3),j)
            data_manbo(k)%f = data_manbo(k)%f - force(k,1+j)%f
            CALL pair_num(pairs_list,SIZE(pairs_list,DIM=1),trimers_list(i,2),trimers_list(i,3),j)
            data_manbo(k)%f = data_manbo(k)%f - force(k,1+j)%f + force(k,1)%f
          END IF
        END DO
      END DO
    
      ! Computing energy
      E1 = 0.0
      DO i=1,n_mols
        E1 = E1 + U(i)
        ! These terms come from the MBE
      END DO
      E2 = 0.0
      DO i=1,npairs
        IF (pairs_list(i,2) <= n_mols) THEN
          E2 = E2 + U(n_mols+i) - U(pairs_list(i,1)) - U(pairs_list(i,2))
        ELSE
          E2 = E2 + U(n_mols+i) - U(pairs_list(i,1))
          E2 = E2 - U(orig_repli_mol(pairs_list(i,2) - n_mols))
        END IF
        ! These terms come from the MBE
      END DO
      E3 = 0.0
      DO i=1,ntrimers
        IF (trimers_list(i,2) <= n_mols .AND. trimers_list(i,3) <= n_mols) THEN
          E3 = E3 + U(n_mols+npairs+i)
          CALL pair_num(pairs_list,SIZE(pairs_list,DIM=1),trimers_list(i,1),trimers_list(i,2),j)
          E3 = E3 - U(n_mols+j)
          CALL pair_num(pairs_list,SIZE(pairs_list,DIM=1),trimers_list(i,1),trimers_list(i,3),j)
          E3 = E3 - U(n_mols+j)
          CALL pair_num(pairs_list,SIZE(pairs_list,DIM=1),trimers_list(i,2),trimers_list(i,3),j)
          E3 = E3 - U(n_mols+j)
          E3 = E3 + U(trimers_list(i,1)) + U(trimers_list(i,2)) + U(trimers_list(i,3))
        ELSE IF (trimers_list(i,2) <= n_mols .AND. trimers_list(i,3) > n_mols) THEN
          E3 = E3 + U(n_mols+npairs+i)
          CALL pair_num(pairs_list,SIZE(pairs_list,DIM=1),trimers_list(i,1),trimers_list(i,2),j)
          E3 = E3 - U(n_mols+j)
          CALL pair_num(pairs_list,SIZE(pairs_list,DIM=1),trimers_list(i,1),trimers_list(i,3),j)
          E3 = E3 - U(n_mols+j)
          CALL pair_num(pairs_list,SIZE(pairs_list,DIM=1),trimers_list(i,2),trimers_list(i,3),j)
          E3 = E3 - U(n_mols+j)
          E3 = E3 + U(trimers_list(i,1)) + U(trimers_list(i,2)) + U(orig_repli_mol(trimers_list(i,3) - n_mols))
        ELSE IF (trimers_list(i,2) > n_mols .AND. trimers_list(i,3) > n_mols) THEN
          CALL find_pair(pairs_list,SIZE(pairs_list,DIM=1),trimers_list(i,2),trimers_list(i,3),t1,t2)
          E3 = E3 + U(n_mols+npairs+i)
          CALL pair_num(pairs_list,SIZE(pairs_list,DIM=1),trimers_list(i,1),trimers_list(i,2),j)
          E3 = E3 - U(n_mols+j)
          CALL pair_num(pairs_list,SIZE(pairs_list,DIM=1),trimers_list(i,1),trimers_list(i,3),j)
          E3 = E3 - U(n_mols+j)
          CALL pair_num(pairs_list,SIZE(pairs_list,DIM=1),t1,t2,j)
          E3 = E3 - U(n_mols+j)
          E3 = E3 + U(trimers_list(i,1)) + U(orig_repli_mol(trimers_list(i,2) - n_mols))
          E3 = E3 + U(orig_repli_mol(trimers_list(i,3) - n_mols))
        END IF
        ! These terms come from the MBE
      END DO
      
      DEALLOCATE(force, U)
    END IF
      
    DEALLOCATE(pairs_list, trimers_list)
    
#ifdef USE_PARALLEL
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif
  END SUBROUTINE forces_eval

  SUBROUTINE forces_eval_full(num)
  ! This subroutine evaluates the forces and energies for full calculations
    IMPLICIT NONE
    INTEGER, INTENT(in) :: num
    INTEGER, DIMENSION(3) :: list
    DOUBLE PRECISION, DIMENSION(1) :: U
    TYPE(forces), DIMENSION(n_atoms,1) :: force

    force%f(1) = 0.0
    force%f(2) = 0.0
    force%f(3) = 0.0
    U = 0.0
    
    list(1) = 1
    IF (num == 0) THEN
      CALL call_qm_prog(force,U,1,1,1,list,4)
    ELSE
      CALL call_qm_prog(force,U,1,1,2,list,5)
    END IF
  
    ! Computing forces
    data_manbo%f(1) = force(:,1)%f(1)
    data_manbo%f(2) = force(:,1)%f(2)
    data_manbo%f(3) = force(:,1)%f(3)
    
    ! Computing energy
    E1 = U(1)
    E2 = 0.0
    E3 = 0.0
  END SUBROUTINE forces_eval_full

  SUBROUTINE call_qm_prog(force,U,val,val2,local,list,y)
  ! This subroutine calls a quantum chemistry program for calculations
    IMPLICIT NONE
    INTEGER, DIMENSION(3), INTENT(in) :: list
    INTEGER, INTENT(in) :: val, val2
    DOUBLE PRECISION, DIMENSION(val), INTENT(out) :: U
    TYPE(forces), DIMENSION(n_atoms,val2), INTENT(out) :: force
    INTEGER, INTENT(in) :: local, y

    IF(qm_prog=="g03" .OR. qm_prog=="g09") THEN
      CALL qm_prog_gaussian(force,U,val,val2,local,list,y)
    ELSE
      PRINT *, "Parameter qm_prog invalid. ManBo cannot run."
      CALL log_write("ERROR: Parameter qm_prog invalid. ManBo cannot run.")
      CALL log_close(1)
      STOP
    END IF
  END SUBROUTINE call_qm_prog

  SUBROUTINE qm_prog_gaussian(force,U,val,val2,local,list,y)
  ! This subroutine uses the quantum chemistry program Gaussian for calculations
    IMPLICIT NONE
    INTEGER, DIMENSION(3), INTENT(in) :: list
    INTEGER, INTENT(in) :: local, y, val, val2
    DOUBLE PRECISION, DIMENSION(val), INTENT(out) :: U
    TYPE(forces), DIMENSION(n_atoms,val2), INTENT(out) :: force
    INTEGER :: j, w, num, j1, j2, j3
    INTEGER :: k, system
    DOUBLE PRECISION, DIMENSION(3) :: pos, pos2, mc2, mc3
    DOUBLE PRECISION :: r12, sener
    CHARACTER(LEN=250) :: arq
    CHARACTER(LEN=250) :: cmd
    LOGICAL :: exists
    
    ! Creating Gaussian scratch folder for this execution
    WRITE(arq,*) local
    num = system("mkdir -p tempgaufiles_" // TRIM(ADJUSTL(name_out)) // "/job" // TRIM(ADJUSTL(arq)))
    IF(num/=0) THEN
      PRINT *, "ManBo cannot run. Error to create a temporary folder using the following command:"
      PRINT *, "  mkdir -p tempgaufiles_" // TRIM(ADJUSTL(name_out)) // "/job" // TRIM(ADJUSTL(arq))
      CALL log_write("ERROR: Error to create a temporary folder using the following command:")
      CALL log_write("       mkdir -p tempgaufiles_" // TRIM(ADJUSTL(name_out)) // "/job" // TRIM(ADJUSTL(arq)))
      CALL log_close(1)
      STOP
    END IF
    
    !Creating the GJF file(s) as the Gaussian input
    IF (y==1) THEN
      WRITE(arq,*) list(1)
      cmd = "inp_" // TRIM(ADJUSTL(name_out)) // "_m" // TRIM(ADJUSTL(arq)) // ".gjf"
    ELSE IF (y==2) THEN
      WRITE(arq,*) list(1)
      cmd = "inp_" // TRIM(ADJUSTL(name_out)) // "_d" // TRIM(ADJUSTL(arq))
      WRITE(arq,*) list(2)
      cmd = TRIM(ADJUSTL(cmd)) // "-" // TRIM(ADJUSTL(arq)) // ".gjf"
    ELSE IF (y==3) THEN
      WRITE(arq,*) list(1)
      cmd = "inp_" // TRIM(ADJUSTL(name_out)) // "_t" // TRIM(ADJUSTL(arq))
      WRITE(arq,*) list(2)
      cmd = TRIM(ADJUSTL(cmd)) // "-" // TRIM(ADJUSTL(arq))
      WRITE(arq,*) list(3)
      cmd = TRIM(ADJUSTL(cmd)) // "-" // TRIM(ADJUSTL(arq)) // ".gjf"
    ELSE IF (y==4) THEN
      cmd = "inp_" // TRIM(ADJUSTL(name_out)) // "_full.gjf"
    ELSE IF (y==5) THEN
      cmd = "inp_" // TRIM(ADJUSTL(name_out)) // "_full2.gjf"
    END IF

    OPEN(UNIT=110+local, FILE=cmd)
        IF(qm_prog_memory > 0) THEN
          WRITE(110+local,'("%Mem=",i8,"MB")') qm_prog_memory
        END IF
      WRITE(110+local,'("%nprocshared=",i4)') n_qm_procs
        IF (use_embedding) THEN
          WRITE(110+local,*) "#", TRIM(ADJUSTL(qm_prog_method)), "/", TRIM(ADJUSTL(qm_prog_basis)), " NoSymm Force Charge"
        ELSE
          WRITE(110+local,*) "#", TRIM(ADJUSTL(qm_prog_method)), "/", TRIM(ADJUSTL(qm_prog_basis)), " NoSymm Force"
        END IF
      WRITE(110+local,*) ""
        IF (y==1) THEN
          WRITE(110+local,'("Monomer ",i8)') list(1)
          WRITE(110+local,*) ""
          j = char_mul(list(1))%q_mol
          num = char_mul(list(1))%mul
          WRITE(110+local,'(i3,i2)') j, num
        ELSE IF (y==2) THEN
          WRITE(110+local,'("Dimer ",2(i8," "))') list(1), list(2)
          WRITE(110+local,*) ""
            IF (list(2) <= n_mols) THEN
              j = char_mul(list(1))%q_mol + char_mul(list(2))%q_mol
              num = char_mul(list(1))%mul-1 + char_mul(list(2))%mul-1 + 1
            ELSE
              k = orig_repli_mol(list(2) - n_mols)
              j = char_mul(list(1))%q_mol + char_mul(k)%q_mol
              num = char_mul(list(1))%mul-1 + char_mul(k)%mul-1 + 1
            END IF
          WRITE(110+local,'(i3,i2)') j, num
        ELSE IF (y==3) THEN
          WRITE(110+local,'("Trimer ",3(i8," "))') list(1), list(2), list(3)
          WRITE(110+local,*) ""
            IF (list(2) <= n_mols .AND. list(3) <= n_mols) THEN
              j = char_mul(list(1))%q_mol + char_mul(list(2))%q_mol + char_mul(list(3))%q_mol
              num = char_mul(list(1))%mul-1 + char_mul(list(2))%mul-1 + char_mul(list(3))%mul-1 + 1
            ELSE IF (list(2) <= n_mols .AND. list(3) > n_mols) THEN
              k = orig_repli_mol(list(3) - n_mols)
              j = char_mul(list(1))%q_mol + char_mul(list(2))%q_mol + char_mul(k)%q_mol
              num = char_mul(list(1))%mul-1 + char_mul(list(2))%mul-1 + char_mul(k)%mul-1 + 1
            ELSE
              k = orig_repli_mol(list(3) - n_mols)
              w = orig_repli_mol(list(2) - n_mols)
              j = char_mul(list(1))%q_mol + char_mul(w)%q_mol + char_mul(k)%q_mol
              num = char_mul(list(1))%mul-1 + char_mul(w)%mul-1 + char_mul(k)%mul-1 + 1
            END IF
          WRITE(110+local,'(i3,i2)') j, num
        ELSE IF (y==4 .OR. y==5) THEN
          WRITE(110+local,'("Full calculation")')
          WRITE(110+local,*) ""
          j = char_mul(list(1))%q_mol
          num = char_mul(list(1))%mul
          WRITE(110+local,'(i3,i2)') j, num
        END IF
        IF (y <= 3) THEN
          DO j=1,y
            IF(list(j) <= n_mols) THEN
              DO k=1,n_atoms
                IF(data_manbo(k)%mol==list(j)) THEN
                  WRITE(110+local,'(1x,a3,4x,3(f14.10,3x))') atom_symbol(data_manbo(k)%an), data_manbo(k)%r(1),&
                                                             data_manbo(k)%r(2), data_manbo(k)%r(3)
                END IF
              END DO
            ELSE
              DO k=1,n_atoms_out
                IF(data_manbo_out(k)%mol==list(j)) THEN
                  WRITE(110+local,'(1x,a3,4x,3(f14.10,3x))') atom_symbol(data_manbo_out(k)%an), data_manbo_out(k)%r(1),&
                              data_manbo_out(k)%r(2), data_manbo_out(k)%r(3)
                END IF
              END DO
            END IF
          END DO
          IF(use_embedding .AND. (.NOT. apply_pbc .OR. (apply_pbc .AND. .NOT. use_emb_radius))) THEN
            WRITE(110+local,*) ""
            IF(y==1) THEN
              DO k=1,n_atoms
                IF(data_manbo(k)%mol/=list(1)) THEN
                  WRITE(110+local,'(1x,3(f14.10,3x),f10.6)') data_manbo(k)%r(1), data_manbo(k)%r(2), data_manbo(k)%r(3),&
                                                             data_manbo(k)%q
                END IF
              END DO
              DO k=1,n_atoms_out
                WRITE(110+local,'(1x,3(f14.10,3x),f10.6)') data_manbo_out(k)%r(1),data_manbo_out(k)%r(2),&
                            data_manbo_out(k)%r(3),data_manbo(data_manbo_out(k)%atom)%q
              END DO
            ELSE IF(y==2) THEN
              DO k=1,n_atoms
                IF(data_manbo(k)%mol/=list(1) .AND. data_manbo(k)%mol/=list(2)) THEN
                  WRITE(110+local,'(1x,3(f14.10,3x),f10.6)') data_manbo(k)%r(1), data_manbo(k)%r(2), data_manbo(k)%r(3),&
                                                             data_manbo(k)%q
                END IF
              END DO
              DO k=1,n_atoms_out
                IF(data_manbo_out(k)%mol/=list(2)) THEN
                  WRITE(110+local,'(1x,3(f14.10,3x),f10.6)') data_manbo_out(k)%r(1),data_manbo_out(k)%r(2),&
                              data_manbo_out(k)%r(3),data_manbo(data_manbo_out(k)%atom)%q
                END IF
              END DO
            ELSE IF(y==3) THEN
              DO k=1,n_atoms
                IF(data_manbo(k)%mol/=list(1) .AND. data_manbo(k)%mol/=list(2) &
                   .AND. data_manbo(k)%mol/=list(3)) THEN
                  WRITE(110+local,'(1x,3(f14.10,3x),f10.6)') data_manbo(k)%r(1), data_manbo(k)%r(2), data_manbo(k)%r(3),&
                                                             data_manbo(k)%q
                END IF
              END DO
              DO k=1,n_atoms_out
                IF(data_manbo_out(k)%mol/=list(2) .AND. data_manbo_out(k)%mol/=list(3)) THEN
                  WRITE(110+local,'(1x,3(f14.10,3x),f10.6)') data_manbo_out(k)%r(1),data_manbo_out(k)%r(2),&
                              data_manbo_out(k)%r(3),data_manbo(data_manbo_out(k)%atom)%q
                END IF
              END DO              
            END IF
          ELSE IF(use_embedding .AND. use_emb_radius .AND. apply_pbc) THEN
            WRITE(110+local,*) ""
            num = INT(emb_radius/MINVAL(box_dim)) + 1
            DO j=1,n_mols
              DO j1=-num,num
              DO j2=-num,num
              DO j3=-num,num
                pos(1) = mc(j)%r(1) + j1*box_dim(1)
                pos(2) = mc(j)%r(2) + j2*box_dim(2)
                pos(3) = mc(j)%r(3) + j3*box_dim(3)
                IF(y==1) THEN
                  pos2 = mc(list(1))%r
                  r12 = vector_module(pos-pos2)
                  IF(r12<=emb_radius .AND. .NOT. (j==list(1) .AND. j1==0 .AND. j2==0 .AND. j3==0)) THEN
                    DO k=1,n_atoms
                      IF (data_manbo(k)%mol==j) THEN
                        WRITE(110+local,'(1x,3(f14.10,3x),f10.6)') data_manbo(k)%r(1) + j1*box_dim(1),&
                                                                   data_manbo(k)%r(2) + j2*box_dim(2),&
                                                                   data_manbo(k)%r(3) + j3*box_dim(3), data_manbo(k)%q
                      END IF
                    END DO
                  END IF
                ELSE IF(y==2) THEN
                  IF(list(2)<=n_mols) THEN
                    pos2 = mc(list(1))%mass*mc(list(1))%r + mc(list(2))%mass*mc(list(2))%r
                    pos2 = pos2/(mc(list(1))%mass + mc(list(2))%mass)
                    r12 = vector_module(pos-pos2)
                    mc2 = mc(list(2))%r
                  ELSE
                    pos2 = mc(list(1))%mass*mc(list(1))%r
                    pos2 = pos2 + mc_out(list(2)-n_mols)%mass*mc_out(list(2)-n_mols)%r
                    pos2 = pos2/(mc(list(1))%mass + mc_out(list(2)-n_mols)%mass)
                    r12 = vector_module(pos-pos2)
                    mc2 = mc_out(list(2)-n_mols)%r
                  END IF
                  IF(r12<=emb_radius .AND. vector_module(pos-mc(list(1))%r)>0.001 &
                     .AND. vector_module(pos-mc2)>0.001) THEN
                    DO k=1,n_atoms
                      IF (data_manbo(k)%mol==j) THEN
                        WRITE(110+local,'(1x,3(f14.10,3x),f10.6)') data_manbo(k)%r(1) + j1*box_dim(1),&
                                                                   data_manbo(k)%r(2) + j2*box_dim(2),&
                                                                   data_manbo(k)%r(3) + j3*box_dim(3), data_manbo(k)%q
                      END IF
                    END DO
                  END IF
                ELSE IF(y==3) THEN
                  IF(list(2)<=n_mols .AND. list(3)<=n_mols) THEN
                    pos2 = mc(list(1))%mass*mc(list(1))%r + mc(list(2))%mass*mc(list(2))%r
                    pos2 = pos2 + mc(list(3))%mass*mc(list(3))%r
                    pos2 = pos2/(mc(list(1))%mass + mc(list(2))%mass + mc(list(3))%mass)
                    r12 = vector_module(pos-pos2)
                    mc2 = mc(list(2))%r
                    mc3 = mc(list(3))%r
                  ELSE IF(list(2)<=n_mols .AND. list(3)>n_mols) THEN
                    pos2 = mc(list(1))%mass*mc(list(1))%r + mc(list(2))%mass*mc(list(2))%r
                    pos2 = pos2 + mc_out(list(3)-n_mols)%mass*mc_out(list(3)-n_mols)%r
                    pos2 = pos2/(mc(list(1))%mass + mc(list(2))%mass + mc_out(list(3)-n_mols)%mass)
                    r12 = vector_module(pos-pos2)
                    mc2 = mc(list(2))%r
                    mc3 = mc_out(list(3)-n_mols)%r
                  ELSE IF(list(2)>n_mols .AND. list(3)>n_mols) THEN
                    pos2 = mc(list(1))%mass*mc(list(1))%r
                    pos2 = pos2 + mc_out(list(2)-n_mols)%mass*mc_out(list(2)-n_mols)%r
                    pos2 = pos2 + mc_out(list(3)-n_mols)%mass*mc_out(list(3)-n_mols)%r
                    pos2 = pos2/(mc(list(1))%mass + mc_out(list(2)-n_mols)%mass + mc_out(list(3)-n_mols)%mass)
                    r12 = vector_module(pos-pos2)
                    mc2 = mc_out(list(2)-n_mols)%r
                    mc3 = mc_out(list(3)-n_mols)%r
                  END IF
                  IF(r12<=emb_radius .AND. vector_module(pos-mc(list(1))%r)>0.001 &
                     .AND. vector_module(pos-mc2)>0.001 &
                     .AND. vector_module(pos-mc3)>0.001) THEN
                    DO k=1,n_atoms
                      IF (data_manbo(k)%mol==j) THEN
                        WRITE(110+local,'(1x,3(f14.10,3x),f10.6)') data_manbo(k)%r(1) + j1*box_dim(1),&
                                                                   data_manbo(k)%r(2) + j2*box_dim(2),&
                                                                   data_manbo(k)%r(3) + j3*box_dim(3), data_manbo(k)%q
                      END IF
                    END DO
                  END IF
                END IF
              END DO
              END DO
              END DO      
            END DO
          END IF
        ELSE IF (y==4 .OR. y==5) THEN
          DO k=1,n_atoms
            WRITE(110+local,'(1x,a3,4x,3(f14.10,3x))') atom_symbol(data_manbo(k)%an), data_manbo(k)%r(1),&
                                                       data_manbo(k)%r(2), data_manbo(k)%r(3)
          END DO
          DO k=1,n_atoms_out
            WRITE(110+local,'(1x,a3,4x,3(f14.10,3x))') atom_symbol(data_manbo_out(k)%an), data_manbo_out(k)%r(1),&
                                                       data_manbo_out(k)%r(2), data_manbo_out(k)%r(3)
          END DO
        END IF
      WRITE(110+local,*) ""
      WRITE(110+local,*) ""
    CLOSE(110+local)

    !Executing Gaussian, reading forces and energies
    WRITE(arq,*) local
    cmd = "export GAUSS_SCRDIR=tempgaufiles_" // TRIM(ADJUSTL(name_out)) // "/job" // TRIM(ADJUSTL(arq)) // ";"
    IF (y==1) THEN
      WRITE(arq,*) list(1)
      cmd = TRIM(ADJUSTL(cmd)) // " " // TRIM(ADJUSTL(qm_prog)) // " inp_" // TRIM(ADJUSTL(name_out)) // "_m" //&
            TRIM(ADJUSTL(arq)) // ".gjf"
    ELSE IF (y==2) THEN
      WRITE(arq,*) list(1)
      cmd = TRIM(ADJUSTL(cmd)) // TRIM(ADJUSTL(qm_prog)) // " inp_" // TRIM(ADJUSTL(name_out)) // "_d" // TRIM(ADJUSTL(arq))
      WRITE(arq,*) list(2)
      cmd = TRIM(ADJUSTL(cmd)) // "-" // TRIM(ADJUSTL(arq)) // ".gjf"
    ELSE IF (y==3) THEN
      WRITE(arq,*) list(1)
      cmd = TRIM(ADJUSTL(cmd)) // TRIM(ADJUSTL(qm_prog)) // " inp_" // TRIM(ADJUSTL(name_out)) // "_t" // TRIM(ADJUSTL(arq))
      WRITE(arq,*) list(2)
      cmd = TRIM(ADJUSTL(cmd)) // "-" // TRIM(ADJUSTL(arq))
      WRITE(arq,*) list(3)
      cmd = TRIM(ADJUSTL(cmd)) // "-" // TRIM(ADJUSTL(arq)) // ".gjf"
    ELSE IF (y==4) THEN
      cmd = TRIM(ADJUSTL(cmd)) // TRIM(ADJUSTL(qm_prog)) // " inp_" // TRIM(ADJUSTL(name_out)) // "_full.gjf"
    ELSE IF (y==5) THEN
      cmd = TRIM(ADJUSTL(cmd)) // TRIM(ADJUSTL(qm_prog)) // " inp_" // TRIM(ADJUSTL(name_out)) // "_full2.gjf"
    END IF
    IF(system(cmd)==0) THEN
      IF (y==1) THEN
        WRITE(arq,*) list(1)
        arq = "inp_" // TRIM(ADJUSTL(name_out)) // "_m" // TRIM(ADJUSTL(arq))
      ELSE IF (y==2) THEN
        WRITE(arq,*) list(1)
        cmd = "inp_" // TRIM(ADJUSTL(name_out)) // "_d" // TRIM(ADJUSTL(arq))
        WRITE(arq,*) list(2)
        arq = TRIM(ADJUSTL(cmd)) // "-" // TRIM(ADJUSTL(arq))
      ELSE IF (y==3) THEN
        WRITE(arq,*) list(1)
        cmd = "inp_" // TRIM(ADJUSTL(name_out)) // "_t" // TRIM(ADJUSTL(arq))
        WRITE(arq,*) list(2)
        cmd = TRIM(ADJUSTL(cmd)) // "-" // TRIM(ADJUSTL(arq))
        WRITE(arq,*) list(3)
        arq = TRIM(ADJUSTL(cmd)) // "-" // TRIM(ADJUSTL(arq))
      ELSE IF (y==4) THEN
        arq = "inp_" // TRIM(ADJUSTL(name_out)) // "_full"
      ELSE IF (y==5) THEN
        arq = "inp_" // TRIM(ADJUSTL(name_out)) // "_full2"
      END IF

      num = 0
      IF (y <= 3) THEN
        DO j=1,y
          IF(list(j) <= n_mols) THEN
            DO k=1,n_atoms
              IF(data_manbo(k)%mol==list(j)) num = num + 1
            END DO
          ELSE
            DO k=1,n_atoms_out
              IF(data_manbo_out(k)%mol==list(j)) num = num + 1
            END DO
          END IF
        END DO
      ELSE
        num = n_atoms
      END IF
      num = num + 2
      WRITE(cmd,*) num
      cmd = "grep -A" // TRIM(ADJUSTL(cmd)) // " 'Forces (Hartrees/Bohr)' " // TRIM(ADJUSTL(arq))
      cmd = TRIM(ADJUSTL(cmd)) // ".log | awk 'NR>3{print $3,$4,$5}' > " // TRIM(ADJUSTL(arq)) // ".tmp"
      num = system(cmd)
      
      cmd = TRIM(ADJUSTL(arq)) // ".tmp"
      OPEN(UNIT=110+local, FILE=cmd, STATUS="old")
        IF (y <= 3) THEN
        DO j=1,y
          IF(list(j) <= n_mols) THEN
            DO k=1,n_atoms
              IF (data_manbo(k)%mol==list(j)) THEN
                IF (y==1) THEN
                  READ(110+local,*) force(k,1)%f(1), force(k,1)%f(2), force(k,1)%f(3)
                ELSE
                  READ(110+local,*) force(k,local-n_mols+1)%f(1), force(k,local-n_mols+1)%f(2), force(k,local-n_mols+1)%f(3)
                END IF
              END IF
            END DO
          END IF
        END DO
        ELSE IF (y==4 .OR. y==5) THEN
          DO k=1,n_atoms
            READ(110+local,*) force(k,1)%f(1), force(k,1)%f(2), force(k,1)%f(3)
          END DO        
        END IF
      CLOSE(110+local, status="delete")

        IF (TRIM(ADJUSTL(qm_prog_method))=="MP2") THEN
          cmd = "grep 'EUMP2 =' " // TRIM(ADJUSTL(arq))
          cmd = TRIM(ADJUSTL(cmd)) // ".log | awk '{print $6}' > " // TRIM(ADJUSTL(arq)) // "_2.tmp"
          num = system(cmd) 
        ELSE
          cmd = "grep 'SCF Done:' " // TRIM(ADJUSTL(arq))
          cmd = TRIM(ADJUSTL(cmd)) // ".log | awk '{print $5}' > " // TRIM(ADJUSTL(arq)) // "_2.tmp"
          num = system(cmd)
        END IF
      
      cmd = TRIM(ADJUSTL(arq)) // "_2.tmp"
      OPEN(UNIT=110+local, FILE=cmd, STATUS="old")
        IF (use_embedding .AND. n_mols /= y .AND. y <= 3) THEN
          cmd = "grep 'Self energy of the charges' " // TRIM(ADJUSTL(arq))
          cmd = TRIM(ADJUSTL(cmd)) // ".log | awk '{print $7}' > " // TRIM(ADJUSTL(arq)) // "_3.tmp"
          num = system(cmd)
          cmd = TRIM(ADJUSTL(arq)) // "_3.tmp"
          OPEN(UNIT=1110+local, FILE=cmd, STATUS="old")
          READ(1110+local,*) sener
          CLOSE(1110+local, status="delete")
        END IF
        IF (y<=3) THEN
          READ(110+local,*) U(local)
            IF (use_embedding .AND. n_mols /= y) THEN
              U(local) = U(local) - sener
            END IF
            ! Here we subtract the terms related to the interactions between point charges of the electrostatic embedding
        ELSE IF (y==4 .OR. y==5) THEN
          READ(110+local,*) U(1)
        END IF
      CLOSE(110+local, status="delete")

      cmd = "rm " // TRIM(ADJUSTL(arq)) // ".gjf " // TRIM(ADJUSTL(arq)) // ".log"
      num = system(cmd)
    ELSE
      PRINT *, "Error running Gaussian in shell, on the forces calculus"
      CALL log_write("")
      CALL log_write(" ERROR: Error running Gaussian in shell, on the forces calculus.")
      CALL log_write("        The error happened on the execution of the following GJF file(s)")
      CALL log_write("        and we have as result the following LOG file(s):")
        IF (y==1) THEN
          WRITE(arq,*) list(1)
          cmd = "inp_" // TRIM(ADJUSTL(name_out)) // "_m" // TRIM(ADJUSTL(arq))
        ELSE IF (y==2) THEN
          WRITE(arq,*) list(1)
          cmd = "inp_" // TRIM(ADJUSTL(name_out)) // "_d" // TRIM(ADJUSTL(arq))
          WRITE(arq,*) list(2)
          cmd = TRIM(ADJUSTL(cmd)) // "-" // TRIM(ADJUSTL(arq))
        ELSE IF (y==3) THEN
          WRITE(arq,*) list(1)
          cmd = "inp_" // TRIM(ADJUSTL(name_out)) // "_t" // TRIM(ADJUSTL(arq))
          WRITE(arq,*) list(2)
          cmd = TRIM(ADJUSTL(cmd)) // "-" // TRIM(ADJUSTL(arq))
          WRITE(arq,*) list(3)
          cmd = TRIM(ADJUSTL(cmd)) // "-" // TRIM(ADJUSTL(arq))
        ELSE IF (y==4) THEN
          cmd = "inp_" // TRIM(ADJUSTL(name_out)) // "_full"
        ELSE IF (y==5) THEN
          cmd = "inp_" // TRIM(ADJUSTL(name_out)) // "_full2"
        END IF
        CALL log_write("          " // TRIM(ADJUSTL(cmd)) // ".gjf")
        INQUIRE(FILE=(TRIM(ADJUSTL(cmd)) // ".log"),EXIST=exists)
          IF(.NOT.exists) THEN
            CALL log_write("          " // TRIM(ADJUSTL(cmd)) // ".log --> Not found.")
          ELSE
            CALL log_write("          " // TRIM(ADJUSTL(cmd)) // ".log")
          END IF
      CALL log_write("        The possible GJF's and LOG's files were not deleted.")
      CALL log_close(1)
      STOP
    END IF
  END SUBROUTINE qm_prog_gaussian

END MODULE manbo_forces
