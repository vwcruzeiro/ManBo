!
! This module sets labels and masses (in atomic units) for atomic numbers  
! Obs.: The position of the element in the vector corresponds to the atomic number
!

MODULE manbo_periodic_table

  IMPLICIT NONE

  CHARACTER(len=2), DIMENSION(105), PARAMETER :: atom_symbol = &
  ! 1st Row
  (/ 'H ',                                                                                                 'He',&
  ! 2nd Row
     'Li', 'Be',                                                             'B ', 'C ', 'N ', 'O ', 'F ', 'Ne',&
  ! 3rd Row
     'Na', 'Mg',                                                             'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar',&
  ! 4th Row
     'K ', 'Ca', 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',&
  ! 5th Row
     'Rb', 'Sr', 'Y ', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I ', 'Xe',&
  ! 6th Row
     'Cs', 'Ba',&
  ! Lanthanides
                 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',&
  ! Continued 6th Row
                       'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',&
  ! 7th Row
     'Fr', 'Ra',&
  ! Actinides
                 'Ac', 'Th', 'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr',&
  ! Aditional types
                       'LP', 'Xx' /)


  DOUBLE PRECISION, DIMENSION(105), PARAMETER :: atom_mass = &
  ! 1st Row
  !  'H ' ,                                                                                                                 'He' 
  (/1.0079,                                                                                                                4.0026,&
  ! 2nd Row
  !  'Li' , 'Be' ,                                                                       'B ' , 'C ' , 'N ' , 'O ' , 'F ' , 'Ne' 
    6.9410,9.0122,                                                                      10.811,12.011,14.007,15.999,18.998,20.180,&
  ! 3rd Row
  !  'Na' , 'Mg' ,                                                                       'Al' , 'Si' , 'P ' , 'S ' , 'Cl' , 'Ar' 
    22.990,24.305,                                                                      26.982,28.086,30.974,32.066,35.453,39.948,&
  ! 4th Row
  !  'K ' , 'Ca' , 'Sc' , 'Ti' , 'V ' , 'Cr' , 'Mn' , 'Fe' , 'Co' , 'Ni' , 'Cu' , 'Zn' , 'Ga' , 'Ge' , 'As' , 'Se' , 'Br' , 'Kr' 
    39.098,40.078,44.956,47.867,50.942,51.996,54.938,55.845,58.933,58.693,63.546,65.390,69.723,72.610,74.922,78.960,79.904,83.800,&
  ! 5th Row
  !  'Rb' , 'Sr' , 'Y ' , 'Zr' , 'Nb' , 'Mo' , 'Tc' , 'Ru' , 'Rh' , 'Pd' , 'Ag' , 'Cd' , 'In' , 'Sn' , 'Sb' , 'Te' , 'I ' , 'Xe' 
    85.468,87.620,88.906,91.224,92.906,95.940,98.906,101.07,102.91,106.42,107.87,112.41,114.82,118.71,121.76,127.60,126.90,131.29,&
  ! 6th Row
  !  'Cs' , 'Ba' 
    132.91,137.33,&
  ! Lanthanides
  !                'La' , 'Ce' , 'Pr' , 'Nd' , 'Pm' , 'Sm' , 'Eu' , 'Gd' , 'Tb' , 'Dy' , 'Ho' , 'Er' , 'Tm' , 'Yb' , 'Lu' 
                  138.91,140.12,140.91,144.24,146.92,150.36,151.96,157.25,158.93,162.50,164.93,167.26,168.93,173.04,174.97,&
  ! Continued 6th Row
  !                       'Hf' , 'Ta' , 'W ' , 'Re' , 'Os' , 'Ir' , 'Pt' , 'Au' , 'Hg' , 'Tl' , 'Pb' , 'Bi' , 'Po' , 'At' , 'Rn' 
                         178.49,180.95,183.84,186.21,190.23,192.22,195.08,196.97,200.59,204.38,207.20,208.98,209.98,209.99,222.02,&
  ! 7th Row
  !  'Fr' , 'Ra' 
    223.02,226.03,&
  ! Actinides
  !                'Ac' , 'Th' , 'Pa' , 'U ' , 'Np' , 'Pu' , 'Am' , 'Cm' , 'Bk' , 'Cf' , 'Es' , 'Fm' , 'Md' , 'No' , 'Lr' 
                  227.03,232.04,231.04,238.03,237.05,239.05,241.06,244.06,249.08,252.08,252.08,257.10,258.10,259.10,262.11,&
  ! Aditional types
  !                       'LP' , 'Xx' 
                          0.00 , 0.00 /)

END MODULE manbo_periodic_table
