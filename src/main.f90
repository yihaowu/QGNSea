program main 
implicit none
include 'rbfconstypes.h'

character*512 :: setupfile
type(setup) :: setupinfo
external scf_rowlands, scf_moy_ramillien

read(*,'(A)')setupfile
write(*,*)setupfile
call read_setup(setupfile, setupinfo)
call region_setup(setupinfo)

call V2X(setupinfo)


end program
