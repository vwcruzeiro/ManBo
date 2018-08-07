MODULE manbo_log

  USE manbo_variables
  ! Here we call manbo_variables module

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE log_start()
  ! This subroutine starts the log file of ManBo
    IMPLICIT NONE
    CHARACTER(LEN=50) :: arq

    arq = TRIM(ADJUSTL(name_out)) // ".log"

    OPEN(UNIT=100, FILE=arq)
    WRITE(100,*) "=============================================================================="
    WRITE(100,*) "                                   ManBo                                      "
    WRITE(100,*) "      Ab Initio Molecular Dynamics Program that uses Many Body Expansion      "
    WRITE(100,*) "=============================================================================="
    WRITE(100,*) "     Created by Vinicius Wilian Dias Cruzeiro and Herbert de Castro Georg     "
    WRITE(100,*) ""
    CALL SYSTEM_CLOCK(init_dt)
    CALL log_time("| ManBo has started at")
    WRITE(100,*) ""
  END SUBROUTINE log_start

  SUBROUTINE log_time(text)
  ! This subroutine write a string on the log file of ManBo including the time in the end of the string
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(in) :: text
    CHARACTER(LEN=10) :: month
    CHARACTER(LEN=30) :: line
    INTEGER, DIMENSION(8) :: time

    CALL date_and_time(VALUES=time)
      IF(time(2)==1) THEN
        month = "January"
      ELSE IF(time(2)==2) THEN
        month = "February"
      ELSE IF(time(2)==3) THEN
        month = "March"
      ELSE IF(time(2)==4) THEN
        month = "April"
      ELSE IF(time(2)==5) THEN
        month = "May"
      ELSE IF(time(2)==6) THEN
        month = "June"
      ELSE IF(time(2)==7) THEN
        month = "July"
      ELSE IF(time(2)==8) THEN
        month = "August"
      ELSE IF(time(2)==9) THEN
        month = "September"
      ELSE IF(time(2)==10) THEN
        month = "October"
      ELSE IF(time(2)==11) THEN
        month = "November"
      ELSE IF(time(2)==12) THEN
        month = "December"
      END IF

    WRITE(line, '(i2.2,":",i2.2,":",i2.2," of ",a," ",i2.2,", ",i4.4)') time(5), time(6), time(7),&
                         TRIM(ADJUSTL(month)), time(3), time(1)
    WRITE(100,'(a," ",a)') text, line
  END SUBROUTINE log_time

  SUBROUTINE log_write(text)
  ! This subroutine write a string in the log file of ManBo
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(in) :: text

    WRITE(100,'(a)') text
  END SUBROUTINE log_write

  SUBROUTINE log_close(stat)
  ! This subroutine closes the log file of ManBo
    IMPLICIT NONE
    INTEGER, INTENT(in) :: stat
    INTEGER :: final_dt
    INTEGER, DIMENSION(5) :: dif_dt

    WRITE(100,*) ""
    
    CALL log_time("| The execution ended at")
    CALL SYSTEM_CLOCK(final_dt)
    CALL dt_diference(init_dt,final_dt,dif_dt)
    WRITE(100,'("| The execution lasted ",i4," days, ",i2," hours, ",i2," minutes, and ",f6.3," seconds.")') dif_dt(1),&
                 dif_dt(2), dif_dt(3), (REAL(dif_dt(4)) + REAL(dif_dt(5))/1000)
      IF(stat==0) THEN
        WRITE(100,'(a)') "| ManBo was successfully executed."
      ELSE
        WRITE(100,'(a)') "| Manbo had an error and stopped running."
      END IF
    WRITE(100,*) ""
    CLOSE(100)
  END SUBROUTINE log_close

  SUBROUTINE dt_diference(idt,fdt,ddt)
  ! This subroutine calculates the diferente between two dates
    IMPLICIT NONE
    INTEGER, INTENT(in) :: idt, fdt
    INTEGER, DIMENSION(5), INTENT(out) :: ddt
    INTEGER :: ifdt

    ddt = 0

    ifdt = fdt - idt
    ddt(1) = INT(REAL(ifdt)/(3600*24*1000))
    ! ddt(1) is the number of days
    ifdt = ifdt - 3600*24*1000*ddt(1)
    ddt(2) = INT(REAL(ifdt)/(3600*1000))
    ! ddt(2) is the number of hours
    ifdt = ifdt - 3600*1000*ddt(2)
    ddt(3) = INT(REAL(ifdt)/(60*1000))
    ! ddt(3) is the number of minutes
    ifdt = ifdt - 60*1000*ddt(3)
    ddt(4) = INT(REAL(ifdt)/(1000))
    ! ddt(4) is the number of seconds
    ifdt = ifdt - 1000*ddt(4)
    ddt(5) = ifdt
    ! ddt(5) is the number of miliseconds
  END SUBROUTINE dt_diference
END MODULE manbo_log
