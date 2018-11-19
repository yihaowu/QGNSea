type region
    real*8 :: lon1, lon2, lat1, lat2  
end type region

type setup
    character*512 :: obssetupfile     
    character*512 :: synthsetupfile   
    character*512 :: priorisetupfile  
    character*512 :: masconsblockfile 
    character*512 :: outputdir		
    character*512 :: recoveredEWH	
    character*512 :: recoveredERR	
    character*512 :: loadlovefile	
    character*512 :: oceanfile		
 
    integer*4 :: typeRBFs  
    integer*4 :: typeconstmat 
    logical :: synth 
    logical :: const_solution 
   
    integer*4 :: n_obs        
    integer*4 :: n_synthobsgroups 
    integer*4 :: n_synthobs       
    
    logical :: write_consmat,write_desmat,write_normat

    real*8 :: regpar_alpha 
    real*8 :: collength 
    integer*4 :: order  
    

    character*512 :: regionfile  
    type(region) :: regions(1)   


    real*8 :: blocksize  
    real*8 :: border     
    integer*4:: blocknum 
end type setup

