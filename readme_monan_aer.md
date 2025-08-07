
! -- Changes in aerosol MONAN
a) Registry.xml : 
  scalars and tendencies: aer01, aer02
  number of aerosols, gases and total tracers
               <package name="aerosols_monan" description="module for aerosols treatment in MPAS-MONAN"/>
               <var name="aer01" array_group="aerosols_monan" units="kg kg^{-1}"
                    description="aerosol 1 mixing ratio"
                    packages="aerosols_monan"/>

b) module atm_time_integration
   
   integer, pointer :: index_aer01,index_aer02

           if (index_aer01 > 0) then
               scalars_driving(index_aer01,1:nVertLevels,1:nCells+1) = mpas_atm_get_bdy_state( clock, block, nVertLevels, nCells, 'aer01', rk_timestep(rk_step) )
            end if






     !--- aerosol section
      call mpas_pool_get_dimension(state, 'index_aer01', index_aer01)
      call mpas_pool_get_dimension(state, 'index_aer02', index_aer02)



c)  module mpas_atmphys_driver

d) Makefile of physics

e) namelist.atmosphere 

      &aerosols_monan
        config_aer_monan = 1

f)



