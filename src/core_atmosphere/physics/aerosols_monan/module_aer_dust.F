! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module module_aerosol_dust
 use mpas_kind_types
 use mpas_pool_routines
 use mpas_timer, only : mpas_timer_start, mpas_timer_stop

 use mpas_atmphys_constants

 implicit none
 private
 public:: driver_aerosol_dust

 contains

!=================================================================================================================
 subroutine driver_aerosol_dust(itimestep,configs,mesh,diag_physics,sfc_input,its,ite,state,time_lev)
!=================================================================================================================
 implicit none
!input and inout arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: sfc_input
 type(mpas_pool_type),intent(in):: state

 integer,intent(in):: its,ite
 integer,intent(in):: itimestep
 integer,intent(in):: time_lev

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics

!local pointers:
 logical,pointer:: config_do_restart
 real(kind=RKIND),dimension(:),pointer:: areaCell
 integer,pointer:: nCellsSolve,nVertLevels
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars

 integer,pointer:: num_aer_monan,num_gen_aer_monan,num_dust_aer_monan
 integer,pointer:: num_gas_monan, num_tracer_monan
 integer,pointer:: config_aer_monan, config_dust_aer_monan 

 integer,pointer,dimension(:) :: aer_ind_monan

 integer,pointer:: index_dst01,index_dst02,index_dst03,index_dst04
 real(kind=RKIND),dimension(:,:)  ,pointer:: dst01,dst02,dst03,dst04

 real(kind=RKIND),dimension(:,:),pointer:: aer_src,aer_snk
 real(kind=RKIND),dimension(:,:),pointer:: zz, zgrid,rho_zz
 real(kind=RKIND),dimension(:,:),pointer:: sh2o,smois
 real(kind=RKIND),dimension(:),pointer:: xland,lai,vegfra,xice
 real(kind=RKIND),dimension(:),pointer:: u10,v10,ust
 integer,dimension(:),pointer:: isltyp,ivgtyp
 real(kind=RKIND),dimension(:),pointer:: rmol ! = 1./L
!local variables:
 integer:: initflag,iCell,k, n_aer
 real(kind=RKIND):: dx, flux
 real,dimension(its:ite):: dens ! density at the surface, kg/m^3   

!CCPP-compliant flags:
 character(len=StrKIND):: errmsg
 integer:: errflg


!initialization of CCPP-compliant flags:
 errmsg = ' '
 errflg = 0
 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)
 call mpas_pool_get_array(mesh,'areaCell',areaCell)
 call mpas_pool_get_array(mesh,'zz'   ,zz   ) ! d(zeta)/dz, vertical metric term, unitless
 call mpas_pool_get_array(mesh,'zgrid',zgrid) ! geometric height of layer interfaces, m MSL
 call mpas_pool_get_array(sfc_input,'xland' ,xland)
 call mpas_pool_get_array(sfc_input,'ivgtyp',ivgtyp)
 call mpas_pool_get_array(sfc_input,'isltyp',isltyp)
 call mpas_pool_get_array(sfc_input,'lai',lai)
 call mpas_pool_get_array(sfc_input,'vegfra',vegfra)
 call mpas_pool_get_array(sfc_input,'xice',xice)
 call mpas_pool_get_array(sfc_input,'sh2o',sh2o)
 call mpas_pool_get_array(sfc_input,'smois',smois)
 call mpas_pool_get_array(diag_physics,'u10',u10)
 call mpas_pool_get_array(diag_physics,'v10',v10)
 call mpas_pool_get_array(diag_physics,'ust',ust)

 call mpas_pool_get_config(configs,'config_do_restart',config_do_restart)
 call mpas_pool_get_config(configs,'config_aer_monan' ,config_aer_monan )
 call mpas_pool_get_config(configs,'config_dust_aer_monan' ,config_dust_aer_monan )

 call mpas_pool_get_dimension(mesh,'num_aer_monan'     , num_aer_monan     )
 call mpas_pool_get_dimension(mesh,'num_dust_aer_monan', num_dust_aer_monan)
 call mpas_pool_get_dimension(mesh,'num_gen_aer_monan' , num_gen_aer_monan )

 call mpas_pool_get_dimension(mesh,'num_gas_monan'   , num_gas_monan   )
 call mpas_pool_get_dimension(mesh,'num_tracer_monan', num_tracer_monan)

 call mpas_pool_get_array(diag_physics,'aer_src_monan', aer_src        )
 call mpas_pool_get_array(diag_physics,'aer_snk_monan', aer_snk        )
 call mpas_pool_get_array(diag_physics,'rmol'    ,rmol    )
 call mpas_pool_get_array(state,'rho_zz' ,rho_zz ,time_lev)



 !copy all MPAS arrays to rectanguler grid:
! call sfclayer_from_MPAS(configs,mesh,diag_physics,sfc_input,its,ite, &
!                         state,time_lev)

 dx = sqrt(maxval(areaCell))
 
 call mpas_pool_get_array(state,'scalars',scalars,time_lev)

 initflag = 1
 if(config_do_restart .or. itimestep > 1) initflag = 0


 
 call mpas_pool_get_array(diag_physics,'aer_ind_monan',aer_ind_monan)

 !ner_20250814:clay fraction and porosity tables
 !1-Sand;2-Loamy Sand;3-Sandy loam;4-Silt loam;5-Silt;6-Loam;7-Sandy clay loam;8-Silty clay loam; 9-Clay loam;10-Sandy clay;11 - Silty clay;12- Clay; 13- Organic material;14-water
 !ref:https://agupubs.onlinelibrary.wiley.com/doi/10.1002/2015JD024558
 !clayfrac=[3,6,10,13,5,18,27,34,34,42,47,58,0,0]
 !watsat=[0.339,0.421,0.434,0.476,0.476,0.439,0.404,0.464,0.465,0.406,0.468,0.468,0.439,1.0] !porosity
 INTEGER,PARAMETER,DIMENSION(14) ::  clayfrac=(/3,6,10,13,5,18,27,34,34,42,47,58,0,0/)!clay
 REAL,PARAMETER,DIMENSION(14) ::  clayfrac=(/0.339,0.421,0.434,0.476,0.476,0.439,0.404,0.464,0.465,0.406,0.468,0.468,0.439,1.0/)!porosity

 !rockfrc=! Rock drag partition

 clayfrc(isltyp(iCell)) = 0.0

 do iCell = 1, nCellsSolve

   
   
   a=1 #Tuning constant for threshold gravimetric soil moisture/a is a tunable constant typically around 0.5–2 (a = 1 was adopted in Kok et al., 2014b) 
   gwc_thr=0.01*a*(0.0017*clayfrc+0.0014*(clayfrc*0.01)**2); #clayfraction [0,1]
   k=1
   dens(iCell) = zz(k,iCell) * rho_zz(k,iCell) ! density at the surface, kg/m^3
   k = 1 ; n_aer=1
   print*,'vegty=',ivgtyp(iCell)
   if (ivgtyp(iCell) == 17) then
      flux = -0.001 / 86400.0  ! 1 g/m^2/day destruction over water
   else
      flux =  0.002 / 86400.0  ! 2 g/m^2/day production over land
      if(ivgtyp(iCell) == 2) flux = 4.*flux
   endif
   aer_src(n_aer,iCell) = flux/  (zgrid(k+1,iCell) - zgrid(k,iCell)) / zz(k,iCell) ! unit: kg/m^3/s
 enddo
 
 !---- dust aerosols section
     call mpas_pool_get_dimension(state,'index_dst01',index_dst01); dst01 => scalars(index_dst01,:,:)
     call mpas_pool_get_dimension(state,'index_dst02',index_dst02); dst02 => scalars(index_dst02,:,:)
     call mpas_pool_get_dimension(state,'index_dst03',index_dst03); dst03 => scalars(index_dst03,:,:)
     call mpas_pool_get_dimension(state,'index_dst04',index_dst04); dst04 => scalars(index_dst04,:,:)
     !print*,'index dust', index_dst01,index_dst04,maxval(dst01),minval(dst01),maxval(dst04),minval(dst04)
 


!copy local arrays to MPAS grid:
! call sfclayer_to_MPAS(configs,sfc_input,diag_physics,its,ite)

!call mpas_log_write('--- end subroutine driver_aerosol_monan.')

 end subroutine driver_aerosol_dust

!=================================================================================================================
 end module module_aerosol_dust
!=================================================================================================================
