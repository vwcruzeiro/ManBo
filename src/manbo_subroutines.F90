MODULE manbo_subroutines

#ifdef USE_INTEL
  USE IFPORT
#endif /* USE_INTEL */
  USE manbo_variables
  ! Here we call manbo_variables module
  USE manbo_log
  ! Here we call manbo_variables module

  IMPLICIT NONE

  CONTAINS

  FUNCTION vector_module(vec)
  ! This function returns the module of a vector
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: vec
    DOUBLE PRECISION :: vector_module
    INTEGER :: i

    vector_module = 0.0
      DO i=1,SIZE(vec, DIM=1)
        vector_module = vector_module + vec(i)*vec(i)
      END DO
    vector_module = SQRT(vector_module)

    RETURN
  END FUNCTION vector_module

  FUNCTION cross_product(vec1,vec2)
  ! This function returns the cross product of two vectors
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(3), INTENT(in) :: vec1, vec2
    DOUBLE PRECISION, DIMENSION(3) :: cross_product

    cross_product(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
    cross_product(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
    cross_product(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)

    RETURN
  END FUNCTION cross_product

  FUNCTION sig(num)
  ! This function returns 1 if the number is positive and -1 if the number is negative
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(in) :: num
    INTEGER :: sig
    
    IF(num < 0) THEN
      sig = -1
    ELSE
      sig = 1
    END IF

    RETURN
  END FUNCTION sig

  FUNCTION valid_pair(m1,m2) RESULT(valid)
  ! This function returns true if a pair of two molecules (m1,m2) is within the cutoff
    IMPLICIT NONE
    INTEGER, INTENT(in) :: m1, m2
    INTEGER :: i, j, n1, n2, num
    LOGICAL :: valid
    
    valid = .FALSE.
    IF (group_monomers==1) THEN
      IF (m1<=n_mols .AND. m2<=n_mols) THEN
        IF (vector_module(mc(m1)%r - mc(m2)%r)<cutoff) valid = .TRUE.
      ELSE IF (m1<=n_mols .AND. m2>n_mols) THEN
        IF (vector_module(mc(m1)%r - mc_out(m2-n_mols)%r)<cutoff) valid = .TRUE.
      ELSE IF (m1>n_mols .AND. m2>n_mols) THEN
        IF (vector_module(mc_out(m1-n_mols)%r - mc_out(m2-n_mols)%r)<cutoff) valid = .TRUE.
      END IF
    ELSE
      num = 0
      DO i=1,3
        IF (m1<=n_mols) THEN
          n1 = orig_mols(m1)%m(i)
        ELSE
          n1 = orig_mols_out(m1-n_mols)%m(i)        
        END IF
        IF (n1==0) EXIT
        DO j=1,3
          IF (m2<=n_mols) THEN
            n2 = orig_mols(m2)%m(j)
          ELSE
            n2 = orig_mols_out(m2-n_mols)%m(j)        
          END IF
          IF (n2==0) EXIT
          IF (n1<=n_mols_orig .AND. n2<=n_mols_orig) THEN
            IF (vector_module(mc_orig(n1)%r - mc_orig(n2)%r)<cutoff) THEN
              num = num+1
              IF (num==2) EXIT
            END IF
          ELSE IF (n1<=n_mols_orig .AND. n2>n_mols_orig) THEN
            IF (vector_module(mc_orig(n1)%r - mc_out_orig(n2-n_mols_orig)%r)<cutoff) THEN
              num = num+1
              IF (num==2) EXIT
            END IF
          ELSE IF (n1>n_mols_orig .AND. n2>n_mols_orig) THEN
            IF (vector_module(mc_out_orig(n1-n_mols_orig)%r - mc_out_orig(n2-n_mols_orig)%r)<cutoff) THEN
              num = num+1
              IF (num==2) EXIT
            END IF
          END IF
        END DO
        IF (num==2) THEN
          valid = .TRUE.
          EXIT
        END IF
      END DO
    END IF

    RETURN
  END FUNCTION valid_pair
  
  SUBROUTINE delsubstr(str,substr)
  ! Deletes all occurrences of substring 'substr' from string 'str' and
  ! shifts characters left to fill holes.
    IMPLICIT NONE
    CHARACTER(LEN=*):: str,substr
    INTEGER :: lensubstr, ipos

    lensubstr=len_trim(substr)
    DO
      ipos=index(str,substr)
      IF(ipos == 0) EXIT
      IF(ipos == 1) THEN
        str=str(lensubstr+1:)
      ELSE
        str=str(:ipos-1)//str(ipos+lensubstr:)
      END IF
    END DO  
   RETURN
  END SUBROUTINE delsubstr

  SUBROUTINE calculate_center_of_mass(x)
  ! This subroutine calculates the center of mass of the molecules inside or outside the box
    INTEGER, INTENT(in) :: x
    ! x == 1 for molecules inside the box,
    ! x == 2 for molecules outside the box, 
    ! x == 3 for original molecules inside the box, and
    ! x == 4 for original molecules outside the box
    INTEGER :: i, j

    IF(x==1) THEN
      DO i=1,n_mols
        mc(i)%r = 0.0
        mc(i)%mass = 0.0
      END DO
      DO i=1,n_atoms
        mc(data_manbo(i)%mol)%r = mc(data_manbo(i)%mol)%r + (data_manbo(i)%mass)*(data_manbo(i)%r)
        mc(data_manbo(i)%mol)%mass = mc(data_manbo(i)%mol)%mass + data_manbo(i)%mass
      END DO
      DO i=1,n_mols
        mc(i)%r = mc(i)%r/mc(i)%mass
      END DO
    ELSE IF (x==2) THEN
      DO i=1,(n_mols_out - n_mols)
        mc_out(i)%r = 0.0
        mc_out(i)%mass = 0.0
      END DO
      DO i=1,n_atoms_out
        j = data_manbo_out(i)%mol - n_mols
        mc_out(j)%r = mc_out(j)%r + (data_manbo_out(i)%mass)*(data_manbo_out(i)%r)
        mc_out(j)%mass = mc_out(j)%mass + data_manbo_out(i)%mass
      END DO
      DO i=1,(n_mols_out - n_mols)
        mc_out(i)%r = mc_out(i)%r/mc_out(i)%mass
      END DO
    ELSE IF(x==3) THEN
      DO i=1,n_mols_orig
        mc_orig(i)%r = 0.0
        mc_orig(i)%mass = 0.0
      END DO
      DO i=1,n_atoms
        mc_orig(data_manbo(i)%mol_orig)%r = mc_orig(data_manbo(i)%mol_orig)%r + (data_manbo(i)%mass)*(data_manbo(i)%r)
        mc_orig(data_manbo(i)%mol_orig)%mass = mc_orig(data_manbo(i)%mol_orig)%mass + data_manbo(i)%mass
      END DO
      DO i=1,n_mols_orig
        mc_orig(i)%r = mc_orig(i)%r/mc_orig(i)%mass
      END DO
    ELSE IF(x==4) THEN
      DO i=1,(n_mols_out_orig - n_mols_orig)
        mc_out_orig(i)%r = 0.0
        mc_out_orig(i)%mass = 0.0
      END DO
      DO i=1,n_atoms_out
        j = data_manbo_out(i)%mol_orig - n_mols_orig
        mc_out_orig(j)%r = mc_out_orig(j)%r + (data_manbo_out(i)%mass)*(data_manbo_out(i)%r)
        mc_out_orig(j)%mass = mc_out_orig(j)%mass + data_manbo_out(i)%mass
      END DO
      DO i=1,(n_mols_out_orig - n_mols_orig)
        mc_out_orig(i)%r = mc_out_orig(i)%r/mc_out_orig(i)%mass
      END DO
    END IF
  END SUBROUTINE calculate_center_of_mass

  SUBROUTINE find_pair(pairs_list,val,m1o,m2o,m1,m2)
  ! This subroutine finds a valid pair (m1,m2) inside the box equivalent to (m1o,m2o) outside the box
    INTEGER, INTENT(out) :: m1,m2
    INTEGER, INTENT(in) :: val,m1o,m2o
    INTEGER, DIMENSION(val,2), INTENT(in) :: pairs_list
    DOUBLE PRECISION :: d1, d2
    INTEGER :: i, j
  
    d1 = vector_module(mc_out(m1o - n_mols)%r - mc_out(m2o - n_mols)%r)
    j = 1
    m1 = orig_repli_mol(m1o - n_mols)
    DO i=1,val
      IF (pairs_list(i,1) == m1) THEN
        IF (pairs_list(i,2) <= n_mols) THEN
          d2 = vector_module(mc(m1)%r - mc(pairs_list(i,2))%r)
        ELSE
          d2 = vector_module(mc(m1)%r - mc_out(pairs_list(i,2)-n_mols)%r)
        END IF
        IF (ABS(d1-d2)<1E-5) THEN
          m2 = pairs_list(i,2)
          j = 0
          EXIT
        END IF
      END IF
    END DO
    IF (j==1) THEN
      m1 = orig_repli_mol(m2o - n_mols)
      DO i=1,val
        IF (pairs_list(i,1) == m1) THEN
          IF (pairs_list(i,2) <= n_mols) THEN
            d2 = vector_module(mc(m1)%r - mc(pairs_list(i,2))%r)
          ELSE
            d2 = vector_module(mc(m1)%r - mc_out(pairs_list(i,2)-n_mols)%r)
          END IF
          IF (ABS(d1-d2)<1E-5) THEN
            m2 = pairs_list(i,2)
            EXIT
          END IF
        END IF
      END DO
    END IF
  END SUBROUTINE find_pair

  SUBROUTINE pair_num(pairs_list,val,m1,m2,local)
  ! This subroutine finds the numbering of a given pair (m1,m2)
    INTEGER, INTENT(in) :: val,m1,m2
    INTEGER, INTENT(out) :: local
    INTEGER, DIMENSION(val,2), INTENT(in) :: pairs_list
    INTEGER :: i
  
    DO i=1,val
      IF (pairs_list(i,1) == m1 .AND. pairs_list(i,2) == m2) THEN
        local = i
        EXIT
      END IF
    END DO
  END SUBROUTINE pair_num

  SUBROUTINE apply_boundary_conditions
  ! This subroutine applies periodics boundary conditions
    IMPLICIT NONE
    INTEGER :: i, j
    DOUBLE PRECISION :: lx2, ly2, lz2
    
    ! Variables to store the values of half of the box lengths in each dimension
    lx2 = box_dim(1)/2
    ly2 = box_dim(2)/2
    lz2 = box_dim(3)/2

      DO i=1,n_mols
        IF(ABS(mc(i)%r(1)) > lx2) THEN
          DO j=1,n_atoms
            IF(data_manbo(j)%mol==i) THEN
              data_manbo(j)%r(1) = data_manbo(j)%r(1) - box_dim(1)*sig(data_manbo(j)%r(1))
            END IF
          END DO
          mc(i)%r(1) = mc(i)%r(1) - box_dim(1)*sig(mc(i)%r(1))
        END IF
        IF(ABS(mc(i)%r(2)) > ly2) THEN
          DO j=1,n_atoms
            IF(data_manbo(j)%mol==i) THEN
              data_manbo(j)%r(2) = data_manbo(j)%r(2) - box_dim(2)*sig(data_manbo(j)%r(2))
            END IF
          END DO
          mc(i)%r(2) = mc(i)%r(2) - box_dim(2)*sig(mc(i)%r(2))
        END IF
        IF(ABS(mc(i)%r(3)) > lz2) THEN
          DO j=1,n_atoms
            IF(data_manbo(j)%mol==i) THEN
              data_manbo(j)%r(3) = data_manbo(j)%r(3) - box_dim(3)*sig(data_manbo(j)%r(3))
            END IF
          END DO
          mc(i)%r(3) = mc(i)%r(3) - box_dim(3)*sig(mc(i)%r(3))
        END IF
      END DO
  END SUBROUTINE apply_boundary_conditions

  SUBROUTINE reorient_box
  ! This subroutine reorients the positions of the atoms inside the box and gets the box dimensions.
  ! The system is translated in such a way the origin of the axis is at the center of mass of the system
  ! and then rotations are performed to place the farthest atom on the edge of the box
    IMPLICIT NONE
    INTEGER :: i
    DOUBLE PRECISION, DIMENSION(3) :: vec, r_old
    DOUBLE PRECISION :: tmass, theta, alp, num, val, val_old
    
    ! Computing the center of mass of the system
    vec = 0.0
    tmass = 0.0
    DO i=1,n_atoms
      vec = vec + data_manbo(i)%mass*data_manbo(i)%r
      tmass = tmass + data_manbo(i)%mass
    END DO
    vec = vec/tmass
    
    ! Translating the atoms so that the origin of the axis stays at the center of mass
    DO i=1,n_atoms
      data_manbo(i)%r = data_manbo(i)%r - vec
    END DO
    
    ! Getting the coordinates of the farthest atom
    vec = 0.0
    DO i=1,n_atoms
      IF (vector_module(data_manbo(i)%r)>vector_module(vec)) vec = data_manbo(i)%r
    END DO
    ! Normalizing the vector
    vec = vec/vector_module(vec)
    
    ! Rotating such that the farthest atom is on the edge of the box
    IF (vec(1)<0.0) THEN
      num = -1.0/SQRT(3.0)
    ELSE
      num = 1.0/SQRT(3.0)
    END IF
    val_old = vec(1)-num
    DO i=1,500
      theta = 2.0*3.14159265359*FLOAT(i)/500.0
      val = vec(1)*COS(theta)+vec(2)*SIN(theta)-num
      IF (val/val_old<0.0) THEN
        theta = 2.0*3.14159265359*(FLOAT(2*i-1)/2.0)/500.0
        EXIT
      END IF
    END DO
    IF (vec(2)<0.0) THEN
      num = -1.0/SQRT(3.0)
    ELSE
      num = 1.0/SQRT(3.0)
    END IF
    val_old = -vec(1)*SIN(theta)+vec(2)*COS(theta)-num
    DO i=1,500
      alp = 2.0*3.14159265359*FLOAT(i)/500.0
      val = -vec(1)*COS(alp)*SIN(theta)+vec(2)*COS(alp)*COS(theta)+vec(3)*SIN(alp)-num
      IF (val/val_old<0.0) THEN
        alp = 2.0*3.14159265359*(FLOAT(2*i-1)/2.0)/500.0
        EXIT
      END IF
    END DO
    DO i=1,n_atoms
      r_old(1) = data_manbo(i)%r(1)
      r_old(2) = data_manbo(i)%r(2)
      r_old(3) = data_manbo(i)%r(3)
      data_manbo(i)%r(1) =  r_old(1)*COS(theta)          + r_old(2)*SIN(theta)
      data_manbo(i)%r(2) = -r_old(1)*COS(alp)*SIN(theta) + r_old(2)*COS(alp)*COS(theta) + r_old(3)*SIN(alp)
      data_manbo(i)%r(3) =  r_old(1)*SIN(alp)*SIN(theta) - r_old(2)*SIN(alp)*COS(theta) + r_old(3)*COS(alp)
      r_old(1) = data_manbo(i)%v(1)
      r_old(2) = data_manbo(i)%v(2)
      r_old(3) = data_manbo(i)%v(3)
      data_manbo(i)%v(1) =  r_old(1)*COS(theta)          + r_old(2)*SIN(theta)
      data_manbo(i)%v(2) = -r_old(1)*COS(alp)*SIN(theta) + r_old(2)*COS(alp)*COS(theta) + r_old(3)*SIN(alp)
      data_manbo(i)%v(3) =  r_old(1)*SIN(alp)*SIN(theta) - r_old(2)*SIN(alp)*COS(theta) + r_old(3)*COS(alp)
    END DO
    
    ! Getting the box dimensions
    box_dim(1) = MAXVAL(data_manbo%r(1)) - MINVAL(data_manbo%r(1))
    box_dim(2) = MAXVAL(data_manbo%r(2)) - MINVAL(data_manbo%r(2))
    box_dim(3) = MAXVAL(data_manbo%r(3)) - MINVAL(data_manbo%r(3))
  END SUBROUTINE reorient_box

  SUBROUTINE propagate_positions
  ! This subroutine updates the positions of the atoms inside the box
    IMPLICIT NONE
    INTEGER :: i

      DO i=1,n_atoms
        data_manbo(i)%r = data_manbo(i)%r + delta_t*data_manbo(i)%v + &
                          0.5*delta_t*delta_t*aua_to_manbo*data_manbo(i)%f/data_manbo(i)%mass        
      END DO
      
      CALL calculate_center_of_mass(1)
      
      IF(apply_pbc) THEN
        CALL apply_boundary_conditions
      END IF
  END SUBROUTINE propagate_positions

  SUBROUTINE propagate_velocities
  ! This subroutine updates the velocities of the atoms inside the box
    IMPLICIT NONE
    INTEGER :: i
  
      DO i=1,n_atoms
        data_manbo(i)%v = data_manbo(i)%v + 0.5*delta_t*aua_to_manbo*data_manbo(i)%f/data_manbo(i)%mass
      END DO
  END SUBROUTINE propagate_velocities

  SUBROUTINE save_data_manbo(x)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: x
    ! x == 1 for step 0, and x != 1 for other steps
    INTEGER :: i
    CHARACTER(LEN=50) :: arq

    IF (x==1) THEN
      arq = TRIM(ADJUSTL(name_out)) // "_positions.xyz"
      OPEN(UNIT=15, FILE=arq)
      arq = TRIM(ADJUSTL(name_out)) // "_velocities.dat"
      OPEN(UNIT=18, FILE=arq)
      arq = TRIM(ADJUSTL(name_out)) // "_forces.dat"
      OPEN(UNIT=17, FILE=arq)
    ELSE
      arq = TRIM(ADJUSTL(name_out)) // "_positions.xyz"
      OPEN(UNIT=15, FILE=arq, POSITION="append")
      arq = TRIM(ADJUSTL(name_out)) // "_velocities.dat"
      OPEN(UNIT=18, FILE=arq, POSITION="append")
      arq = TRIM(ADJUSTL(name_out)) // "_forces.dat"
      OPEN(UNIT=17, FILE=arq, POSITION="append")
    END IF
    WRITE(15,*) n_atoms
    WRITE(15,'("Positions of atoms inside and outside of the box - Step ",i8)') md_step
    WRITE(18,*) n_atoms
    WRITE(18,'("Velocities of atoms inside the box - Step ",i8)') md_step
    WRITE(17,*) n_atoms
    WRITE(17,'("Forces of atoms inside the box - Step ",i8)') md_step
      DO i=1,n_atoms
        WRITE(15,'(2x,a2,4x,3(f14.10,4x))') atom_symbol(data_manbo(i)%an), data_manbo(i)%r(1),&
                                            data_manbo(i)%r(2), data_manbo(i)%r(3)
        WRITE(18,'(2x,3(f14.10,4x))') data_manbo(i)%v(1), data_manbo(i)%v(2),&
                                            data_manbo(i)%v(3)
        WRITE(17,'(2x,3(f14.10,4x))') data_manbo(i)%f(1), data_manbo(i)%f(2), data_manbo(i)%f(3)
      END DO
    CLOSE(17)
    CLOSE(18)
    CLOSE(15)
  END SUBROUTINE save_data_manbo
  
  SUBROUTINE save_rst_file()
    IMPLICIT NONE
    INTEGER :: i
    CHARACTER(LEN=50) :: arq
    CHARACTER(LEN=200) :: line, text
    LOGICAL :: declared

    arq = TRIM(ADJUSTL(name_out)) // ".rst.in"
    OPEN(UNIT=15, FILE=arq)
    
    OPEN(UNIT=13,FILE=name_in,STATUS='OLD')
      line = ""
      declared = .FALSE.
      DO WHILE (line/="END_OF_INPUT_VARIABLES")
        READ(13, '(a200)') text
        READ(text, *) line
        IF (line=="box_dim=") declared = .TRUE.
        IF (line=="n_steps=") THEN
          WRITE(15,*) "    n_steps=           ", (n_steps-md_step)
        ELSE IF (line=="END_OF_INPUT_VARIABLES" .AND. .NOT. declared) THEN
          WRITE(15,*) "    box_dim=           ", box_dim(1), box_dim(2), box_dim(3)
          WRITE(15,*) text
        ELSE
          WRITE(15,*) text
        END IF
      END DO
      READ(13, '(a200)') text
      WRITE(15,*) text
      READ(13, '(a200)') text
      WRITE(15,*) text
      DO i=1,n_atoms
        READ(13, '(a200)') text
        WRITE(15,'(3x,a2,4x,3(f14.10,3x),2x,i5)') atom_symbol(data_manbo(i)%an), data_manbo(i)%r(1),&
                                                  data_manbo(i)%r(2), data_manbo(i)%r(3),&
                                                  data_manbo(i)%mol_orig
      END DO
      READ(13, '(a200)') text
      WRITE(15,*) text
      DO i=1,n_mols_orig
        READ(13, '(a200)') text
        WRITE(15,'(3x,i4,5x,i3,10x,i2)') i, char_mul_orig(i)%q_mol, char_mul_orig(i)%mul
      END DO
      IF (use_embedding) THEN
        READ(13, '(a200)') text
        WRITE(15,*) text
        DO i=1,n_atoms
          READ(13, '(a200)') text
          WRITE(15,'(3x,f10.6)') data_manbo(i)%q
        END DO
      END IF
      READ(13, '(a200)') text
      WRITE(15,*) text
      DO i=1,n_atoms
          WRITE(15,'(3x,3(f14.10,3x))') data_manbo(i)%v(1), data_manbo(i)%v(2), data_manbo(i)%v(3)
      END DO
    CLOSE(13)
    
    CLOSE(15)
  END SUBROUTINE save_rst_file

  SUBROUTINE print_properties(x)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: x
    ! x == 1 for step 0, and x != 1 for other steps
    INTEGER :: i
    CHARACTER(LEN=50) :: arq
    DOUBLE PRECISION :: E_kin
    DOUBLE PRECISION, DIMENSION(3) :: p_sys, L_sys
    ! p_sys is the total linear momentum of the system, given in atomic unit (1.99285175E-24 kg.m/s)
    ! L_sys is the total angular momentum of the system, given in atomic unit (1.054571726E−34 J·s)

    arq = TRIM(ADJUSTL(name_out)) // "_properties.dat"
    IF (x==1) THEN
      OPEN(UNIT=15, FILE=arq)
    ELSE
      OPEN(UNIT=15, FILE=arq, POSITION="append")
    END IF
    WRITE(15,*) "-------------------------------------------------------------------------------------"
    WRITE(15,'(" Step ",i8," of ",i8)') md_step, n_steps
    WRITE(15,*) "-------------------------------------------------------------------------------------"
    temperature = 0.0
      DO i=1,n_atoms
        temperature = temperature + (data_manbo(i)%mass)*((data_manbo(i)%v(1))**2 + (data_manbo(i)%v(2))**2)
        temperature = temperature + (data_manbo(i)%mass)*((data_manbo(i)%v(3))**2)
      END DO
    temperature = convert_energy*temperature/n_atoms
    temperature = temperature/(3*const_boltzmann)
    WRITE(15,'(a,f12.4,a)') " System Temperature: ", temperature, " K"
      IF(apply_therm) THEN
        WRITE(15,'(a,f12.6)') " Temperature/Target Temperature: ", temperature/target_therm_temp
      END IF
    p_sys = 0.0
      DO i=1,n_atoms
        p_sys = p_sys + aup_to_manbo*(data_manbo(i)%mass)*(data_manbo(i)%v)
      END DO
    WRITE(15,*) " "
    WRITE(15,'(a,3(f12.6,2x),a)') " System's total linear momentum:  ", p_sys(1), p_sys(2),&
                                  p_sys(3), ' a.u. ("reduced Plank constant"/bohr)'
    WRITE(15,'(a,f12.6,2x,a)') " System's total linear momentum modulus:  ", vector_module(p_sys),&
                               ' a.u. ("reduced Plank constant"/bohr)'

    L_sys = 0.0
      DO i=1,n_atoms
        L_sys = L_sys + aul_to_manbo*(data_manbo(i)%mass)*cross_product(data_manbo(i)%r,data_manbo(i)%v)
      END DO
    WRITE(15,*) " "
    WRITE(15,'(a,3(f12.6,2x),a)') " System's total angular momentum:  ", L_sys(1), L_sys(2),&
                                  L_sys(3), ' a.u. ("reduced Plank constant"/bohr)'
    WRITE(15,'(a,f12.6,2x,a)') " System's total angular momentum modulus:  ", vector_module(L_sys),&
                               ' a.u. ("reduced Plank constant"/bohr)'
    WRITE(15,*) " "
    WRITE(15,'(a,f16.8,2x,a)') " One-Body Term Energy (E1):  ", E1, " hartrees"
    WRITE(15,'(a,f16.8,2x,a)') " Two-Body Term Energy (E2):  ", E2, " hartrees"
    WRITE(15,'(a,f16.8,2x,a)') " Tree-Body Term Energy (E3): ", E3, " hartrees"
    WRITE(15,'(a,f16.8,2x,a)') " Total energy: ", (E1+E2+E3), " hartrees"
    E_kin = 0.0
      DO i=1,n_atoms
        E_kin = E_kin + 0.5*(data_manbo(i)%mass)*((data_manbo(i)%v(1))**2)
        E_kin = E_kin + 0.5*(data_manbo(i)%mass)*((data_manbo(i)%v(2))**2)
        E_kin = E_kin + 0.5*(data_manbo(i)%mass)*((data_manbo(i)%v(3))**2)
      END DO
    E_kin = convert_energy*E_kin*convert_jou_har
    WRITE(15,'(a,f16.8,2x,a)') " Nuclei kinetic energy (E_kin): ", E_kin, " hartrees"
    WRITE(15,'(a,f16.8,2x,a)') " Total + kinetic energy (E1 + E2 + E3 + E_kin): ", (E1 + E2 + E3 + E_kin), " hartrees"
    WRITE(15,'(a,f16.6,2x,a)') " Total energy per molecule: ", convert_mol_energy*(E1 + E2 + E3 + E_kin)/n_mols, " kcal/mol"
    WRITE(15,'(a,f16.6,2x,a)') " E2 + E3 per molecule: ", convert_mol_energy*(E2 + E3)/n_mols, " kcal/mol"
    CLOSE(15)
  END SUBROUTINE print_properties

  SUBROUTINE apply_thermalization()
    IMPLICIT NONE
    INTEGER :: i
    CHARACTER(LEN=60) :: line

    ! Calculating the system's temperature
    temperature = 0.0
      DO i=1,n_atoms
        temperature = temperature + (data_manbo(i)%mass)*((data_manbo(i)%v(1))**2 + (data_manbo(i)%v(2))**2)
        temperature = temperature + (data_manbo(i)%mass)*((data_manbo(i)%v(3))**2)
      END DO
    temperature = convert_energy*temperature/n_atoms
    temperature = temperature/(3*const_boltzmann)

    CALL log_write("")
    WRITE(line,'(f8.4)') temperature
    CALL log_write("  System's Temperature: " // TRIM(ADJUSTL(line)) // " K")
    WRITE(line,'(f8.4)') (temperature/target_therm_temp)
    CALL log_write("  Temperature/Target Temperature: " // TRIM(ADJUSTL(line)))
    ! Here we wrote some informations in ManBo log

    DO i=1,n_atoms
      data_manbo(i)%v = SQRT(1 + (delta_t/t_therm)*(target_therm_temp/temperature - 1))*data_manbo(i)%v
    END DO
    ! Here we applied the thermalization on velocities
    
  END SUBROUTINE apply_thermalization

  SUBROUTINE generate_rand_vel()
  ! This subroutine generate random velocities at a certain temperature
    IMPLICIT NONE
    INTEGER :: i
    DOUBLE PRECISION :: rv_temp
    DOUBLE PRECISION, DIMENSION(3) :: rv_p, rv_L
    CHARACTER(LEN=60) :: line

    CALL SRAND(rand_vel_seed)

    DO i=1,n_atoms
      data_manbo(i)%v(1) = RAND()
      data_manbo(i)%v(2) = RAND()
      data_manbo(i)%v(3) = RAND()
    END DO

    rv_p = 0.0
    rv_L = 0.0
      DO i=1,n_atoms
        rv_p = rv_p + (data_manbo(i)%mass)*(data_manbo(i)%v)
        rv_L = rv_L + (data_manbo(i)%mass)*cross_product(data_manbo(i)%r,data_manbo(i)%v)
      END DO

    DO i=1,n_atoms
      data_manbo(i)%v = data_manbo(i)%v - (1/(n_atoms*data_manbo(i)%mass))*rv_p
    END DO

    rv_temp = 0.0
      DO i=1,n_atoms
        rv_temp = rv_temp + (data_manbo(i)%mass)*((data_manbo(i)%v(1))**2 + (data_manbo(i)%v(2))**2)
        rv_temp = rv_temp + (data_manbo(i)%mass)*((data_manbo(i)%v(3))**2)
      END DO
    rv_temp = convert_energy*rv_temp/n_atoms
    rv_temp = rv_temp/(3*const_boltzmann)

    DO i=1,n_atoms
      data_manbo(i)%v = SQRT(rand_vel_temp/rv_temp)*data_manbo(i)%v
    END DO

  END SUBROUTINE

  SUBROUTINE group_mons()
  ! This subroutine groups close molecules into single monomers
    IMPLICIT NONE
    INTEGER :: i, j, k, num, cnt, n, maxn
    INTEGER, DIMENSION(2) :: pair
    INTEGER, DIMENSION(3) :: trimer
    INTEGER, DIMENSION(n_mols_orig) :: complete
    DOUBLE PRECISION, DIMENSION(n_mols_orig-1,n_mols_orig) :: distp
    DOUBLE PRECISION, DIMENSION(n_mols_orig-2,n_mols_orig-1,n_mols_orig) :: distt
    DOUBLE PRECISION :: dist_max
    CHARACTER(LEN=60) :: line
    
    complete = 1
    ! Here 1 means no molecules has been added to a new group

    ! Reallocating important variables
    IF (.not. ALLOCATED(mc_orig)) ALLOCATE(mc_orig(n_mols_orig))
    
    ! Computing distances
    n_mols = n_mols_orig
    CALL calculate_center_of_mass(3) ! Centers of mass of the original molecules inside the box
    IF (group_monomers == 2 .AND. n_mols_orig > 1) THEN
      IF (mod(n_mols_orig,2) == 0) THEN
        n_mols = n_mols_orig/2
      ELSE
        n_mols = INT(n_mols_orig/2.0)+1
      END IF

      distp = 0.0
      maxn = 1
      DO i=1,n_mols_orig-1
         DO j=i+1,n_mols_orig
           distp(i,j) = vector_module(mc_orig(i)%r - mc_orig(j)%r)
           maxn = maxn + 1
         END DO
      END DO

      dist_max = MAXVAL(distp)
    ELSE IF (group_monomers == 3 .AND. n_mols_orig > 2) THEN
      IF (mod(n_mols_orig,3) == 0) THEN
        n_mols = n_mols_orig/3
      ELSE
        n_mols = INT(FLOAT(n_mols_orig)/3.0)+1
      END IF

      distt = 0.0
      maxn = 1
      DO i=1,n_mols_orig-2
         DO j=i+1,n_mols_orig-1
           DO k=j+1,n_mols_orig
             distt(i,j,k) = vector_module(mc_orig(i)%r - mc_orig(j)%r) + vector_module(mc_orig(i)%r - mc_orig(k)%r)
             distt(i,j,k) = distt(i,j,k) + vector_module(mc_orig(j)%r - mc_orig(k)%r)
             maxn = maxn + 1
           END DO
         END DO
      END DO

      dist_max = MAXVAL(distt)    
    END IF
    
    ! Reallocating important variables
    IF (md_step == 0) THEN
      IF (.not. ALLOCATED(orig_mols)) ALLOCATE(orig_mols(n_mols))
      DEALLOCATE(mc, char_mul)
      ALLOCATE(mc(n_mols), char_mul(n_mols))
      IF (ALLOCATED(orig_repli_mol)) THEN
        num = NINT(FLOAT(SIZE(orig_repli_mol)*n_mols)/FLOAT(n_mols_orig))
        DEALLOCATE(orig_repli_mol)
        ALLOCATE(orig_repli_mol(num))
      END IF
    END IF
    
    ! Setting all values of orig_mols to zero. If a group has only one or two molecules, then the vacant positions will be identified by the zeros.
    orig_mols%m(1) = 0
    orig_mols%m(2) = 0
    orig_mols%m(3) = 0

    ! Grouping molecules
    cnt = 1
    n = 1
    IF (group_monomers == 2 .AND. n_mols_orig > 1) THEN
      DO WHILE (cnt <= n_mols)
        pair = MINLOC(distp,(distp>0.0))
        IF (complete(pair(1)) == 1 .AND. complete(pair(2)) == 1) THEN
          DO num=1,n_atoms
            IF (data_manbo(num)%mol_orig == pair(1) .OR. data_manbo(num)%mol_orig == pair(2)) THEN
              data_manbo(num)%mol = cnt
            END IF
          END DO
          char_mul(cnt)%q_mol = char_mul_orig(pair(1))%q_mol + char_mul_orig(pair(2))%q_mol
          char_mul(cnt)%mul = char_mul_orig(pair(1))%mul-1 + char_mul_orig(pair(2))%mul-1 + 1
          complete(pair(1)) = 0
          complete(pair(2)) = 0
          orig_mols(cnt)%m(1) = pair(1)
          orig_mols(cnt)%m(2) = pair(2)
          cnt = cnt+1
        ELSE IF (cnt == n_mols .AND. mod(n_mols_orig,2) > 0) THEN
          DO i=1,n_mols_orig
            IF (complete(i) == 1) THEN
              DO num=1,n_atoms
                IF (data_manbo(num)%mol_orig == i) THEN
                  data_manbo(num)%mol = cnt
                END IF
              END DO
              char_mul(cnt) =  char_mul_orig(i)
              orig_mols(cnt)%m(1) = i
              EXIT
            END IF
          END DO
          cnt = cnt+1
        END IF
        distp(pair(1),pair(2)) = dist_max + 1.0
        n = n + 1
        IF(n == maxn+1) THEN
          PRINT *, "Error on grouping monomers. ManBo cannot run."
          CALL log_write("ERROR: Error on grouping monomers on manbo_subroutines.F90")
          CALL log_close(1)
          STOP
        END IF
      END DO

      CALL log_write("  Grouped each 2 closer molecules into single monomers.")
      CALL log_write("")
    ELSE IF (group_monomers == 3 .AND. n_mols_orig > 2) THEN
      DO WHILE (cnt <= n_mols)
        trimer = MINLOC(distt,(distt>0.0))
        IF (complete(trimer(1)) == 1 .AND. complete(trimer(2)) == 1 .AND. complete(trimer(3)) == 1) THEN
          DO num=1,n_atoms
            IF (data_manbo(num)%mol_orig == trimer(1) .OR. data_manbo(num)%mol_orig == trimer(2)&
                .OR. data_manbo(num)%mol_orig == trimer(3)) THEN
              data_manbo(num)%mol = cnt
            END IF
          END DO
          char_mul(cnt)%q_mol = char_mul_orig(trimer(1))%q_mol + char_mul_orig(trimer(2))%q_mol + char_mul_orig(trimer(3))%q_mol
          char_mul(cnt)%mul = char_mul_orig(trimer(1))%mul-1 + char_mul_orig(trimer(2))%mul-1 + char_mul_orig(trimer(3))%mul-1 + 1
          complete(trimer(1)) = 0
          complete(trimer(2)) = 0
          complete(trimer(3)) = 0
          orig_mols(cnt)%m(1) = trimer(1)
          orig_mols(cnt)%m(2) = trimer(2)
          orig_mols(cnt)%m(3) = trimer(3)
          cnt = cnt+1
        ELSE IF (cnt == n_mols .AND. mod(n_mols_orig,3) > 0) THEN
          char_mul(cnt)%q_mol = 0
          char_mul(cnt)%mul = 0
          j = 1
          DO i=1,n_mols_orig
            IF (complete(i) == 1) THEN
              DO num=1,n_atoms
                IF (data_manbo(num)%mol_orig == i) THEN
                  data_manbo(num)%mol = cnt
                END IF
              END DO
              char_mul(cnt)%q_mol = char_mul(cnt)%q_mol + char_mul_orig(i)%q_mol
              char_mul(cnt)%mul = char_mul(cnt)%mul + char_mul_orig(i)%mul-1
              orig_mols(cnt)%m(j) = i
              j = j + 1
            END IF
          END DO
          char_mul(cnt)%mul = char_mul(cnt)%mul + 1
          cnt = cnt+1
        END IF
        distt(trimer(1),trimer(2),trimer(3)) = dist_max + 1.0
        n = n + 1
        IF(n == maxn+1) THEN
          PRINT *, "Error on grouping monomers. ManBo cannot run."
          CALL log_write("ERROR: Error on grouping monomers on manbo_subroutines.F90")
          CALL log_close(1)
          STOP
        END IF
      END DO

      CALL log_write("  Grouped each 3 closer molecules into single monomers.")
      CALL log_write("")
    END IF
    
  END SUBROUTINE

END MODULE manbo_subroutines
