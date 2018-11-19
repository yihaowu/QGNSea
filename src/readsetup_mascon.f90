subroutine read_setup(setupfile, settings)
    implicit none
    include 'rbfconstypes.h'
    character*512, intent(in) :: setupfile
    type(setup), intent(out) :: settings
    character*512 :: obssetupfile     ! setup file for the observations
    character*512 :: synthsetupfile   ! setup file for the synthesis data
    character*512 :: priorisetupfile    ! setup file for the priori mascons
    character*512 :: masconsblockfile ! setup file for mascons or other RBFs 
    character*512 :: outputdir  ! the folder name for the output files
    character*512 :: loadlovefile	! setup file for load love numbers up to degree 200
    character*512 :: oceanfile		! setup file for marking land and ocean
    character*512 :: recoveredEWH	! file name of the results
    character*512 :: recoveredERR	! file name of the results error
    
    integer*4 :: typeRBFs  !  type of the radial basis functions (e.g., Mascons/Point mass, Poisson wavelets, etc)
    integer*4 :: typeconstmat! type of the constraint matrices  (e.g., from Ramillien, rowlands, etc)
    logical :: synth ! whether you want to evaluate your solutions
    logical :: const_solution ! whether you want to add constraints
    
    integer*4 :: n_obs               ! number of the observations  (I suppose only one observation group used in Mascon code, e.g., disturbing potential difference)
    integer*4 :: n_synthobsgroups  ! number of the sythesis groups (there may be several sythesis groups)
    integer*4 :: n_synthobs        ! number of the sythesis data
    
    logical :: write_consmat,write_desmat,write_normat ! whether you want to output the constraint/design/normal matrices into disk space

    real*8 :: regpar_alpha ! the scale factor of the constraint matrices
    real*8 :: collength ! the corelation length used in the constraint matrices
    integer :: order   !maximum degree of loading love numbers
    
    ! settings for regionfile
    character*512 :: regionfile  ! the setup file for determining boundary of the target region
    type(region) :: regions(1)   !  (I assume there is only one target region)

    ! settings for the blocks 
     real*8 :: blocksize   !the width of blocks, in degrees
     real*8 :: border     !extending distance of interested area
     integer :: blocknum    !number of mascon blocks


    integer :: status
     
    namelist /masconsetup/ obssetupfile,synthsetupfile,priorisetupfile,masconsblockfile,loadlovefile
    namelist /masconsetup/ oceanfile,outputdir,regionfile,recoveredEWH,recoveredERR
    namelist /masconsetup/ synth,const_solution,write_consmat,write_desmat,write_normat
    namelist /masconsetup/ typeRBFs,typeconstmat,n_obs,n_synthobsgroups,n_synthobs,order,blocknum
    namelist /masconsetup/ regpar_alpha,collength,blocksize,border
    
        write(*,'(A)') 'reading setup...'

        ! read setup
	open(unit = 1,file = trim(setupfile),iostat = status, action = 'read')
	read(unit = 1,nml = masconsetup)
	rewind(1)
	close(1)

        settings%obssetupfile = obssetupfile
        settings%synthsetupfile = synthsetupfile
        settings%priorisetupfile = priorisetupfile
        settings%masconsblockfile = masconsblockfile
	settings%loadlovefile = loadlovefile
	settings%oceanfile = oceanfile
        settings%outputdir = outputdir
        settings%regionfile = regionfile
        settings%recoveredEWH = recoveredEWH
        settings%recoveredERR = recoveredERR
        
        settings%synth = synth
        settings%const_solution = const_solution
        settings%write_consmat = write_consmat
        settings%write_desmat = write_desmat
        settings%write_normat = write_normat
        
        settings%typeRBFs = typeRBFs
        settings%typeconstmat = typeconstmat
        settings%n_obs = n_obs
        settings%n_synthobsgroups = n_synthobsgroups
        settings%n_synthobs = n_synthobs
        settings%order = order
        settings%blocknum = blocknum
        
        settings%regpar_alpha = regpar_alpha
        settings%collength = collength
        settings%blocksize = blocksize
        settings%border = border        
end subroutine read_setup

subroutine region_setup(settings)
    implicit none
    include 'masconstypes.h'
    
    type(setup),intent(inout) :: settings
    real*8 :: lon1, lon2, lat1, lat2
	integer :: status

    open(1, file=trim(settings%regionfile), iostat = status, action = 'read')
        read(1,*) lon1, lon2, lat1, lat2
        settings%regions(1)%lon1 = lon1
        settings%regions(1)%lon2 = lon2
        settings%regions(1)%lat1 = lat1
        settings%regions(1)%lat2 = lat2
    close(1)
end subroutine region_setup


