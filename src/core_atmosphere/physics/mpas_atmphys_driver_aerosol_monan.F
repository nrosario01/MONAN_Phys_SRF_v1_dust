! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_driver_aer_monan
 use mpas_kind_types
 use mpas_pool_routines
 use mpas_timer, only : mpas_timer_start, mpas_timer_stop

 use mpas_atmphys_constants
 use mpas_atmphys_vars

 implicit none
 private
 public:: init_aerosol_monan,       &
          allocate_aerosol_monan,   &
          deallocate_aerosol_monan, &
          driver_aerosol_monan



! * coupling gust front speed from downdrafts with surface winds 
!   Saulo R. Freitas (saulo.freitas@inpe.br) / 2025-01-19
 
 contains


!=================================================================================================================
 subroutine allocate_aerosol_monan(configs)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs

!local pointers:
 integer, pointer :: config_aer_monan
 integer, pointer:: num_aer_monan, num_gas_monan, num_tracer_monan

!-----------------------------------------------------------------------------------------------------------------

 !call mpas_pool_get_config(configs, 'config_aer_monan', config_aer_monan)
 !call mpas_pool_get_dimension(mesh, 'num_aer_monan', num_aer_monan)
 !call mpas_pool_get_dimension(mesh, 'num_gas_monan', num_gas_monan)
 !call mpas_pool_get_dimension(mesh, 'num_tracer_monan', num_tracer_monan)
 

 !if(.not.allocated(aer_snk_p) ) allocate(aer_snk_p(num_aer_monan,ims:ime,jms:jme))
 !if(.not.allocated(aer_src_p) ) allocate(aer_src_p(num_aer_monan,ims:ime,jms:jme))
 !if(.not.allocated(aer_p)     ) allocate(aer_p    (num_aer_monan,ims:ime,kms:kme,jms:jme))


 end subroutine allocate_aerosol_monan

!=================================================================================================================
 subroutine deallocate_aerosol_monan(configs)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs

!local pointers:

 if(allocated(aer_snk_p) ) deallocate(aer_snk_p)
 if(allocated(aer_src_p) ) deallocate(aer_src_p)
 if(allocated(aer_p)     ) deallocate(aer_p    )


 end subroutine deallocate_aerosol_monan

!=================================================================================================================
 subroutine aerosol_from_MPAS(configs,mesh,diag_physics,sfc_input,its,ite, &
                               state,time_lev)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: sfc_input
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(in):: state

 integer,intent(in):: its,ite
 integer,intent(in):: time_lev

!local variables:
 integer:: i,j,k

!local pointers:
 logical,pointer:: config_frac_seaice
 character(len=StrKIND),pointer:: sfclayer_scheme

 real(kind=RKIND),pointer:: len_disp
 real(kind=RKIND),dimension(:),pointer:: meshDensity
 real(kind=RKIND),dimension(:),pointer:: skintemp,sst,xice,xland
 real(kind=RKIND),dimension(:),pointer:: hpbl,mavail
 real(kind=RKIND),dimension(:),pointer:: br,cpm,chs,chs2,cqs2,flhc,flqc,gz1oz0,hfx,qfx,  &
                                         qgh,qsfc,lh,mol,psim,psih,regime,rmol,ust,ustm, &
                                         wspd,znt,zol

!local pointers specific to monin_obukhov:
 real(kind=RKIND),dimension(:),pointer:: fh,fm

!local pointers specific to mynn:
 real(kind=RKIND),dimension(:),pointer:: ch,qcg,snowh

!----------------------------------------------------------------------------------------------------------------- 
!-srf: local variables and pointers specific to GF cumulus scheme: 
 real,parameter :: mx_buoy1 =   250.5 !J/kg
 real,parameter :: mx_buoy2 = 20004.0 !J/kg 
 real,parameter :: adjust_scale = 1.0
 character(len=StrKIND),pointer:: convection_scheme
 integer,pointer:: index_buoyx 
 integer,pointer:: config_gf_gustf ! parameter to couple the gustfront speed with the 
                                    ! surface layer fluxes    
 integer,pointer:: config_gf_cporg
 real(kind=RKIND),dimension(:),pointer:: rmfxdncu,v2d2
 real(kind=RKIND),dimension(:),pointer:: u_gustfront,v_gustfront
 real(kind=RKIND),dimension(:,:),pointer:: buoyx
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars
 real :: wspd_local, factor1, factor2, gspd_local
!-----------------------------------------------------------------------------------------------------------------

!input variables:
 call mpas_pool_get_config(configs,'config_aerosol_monan_scheme',sfclayer_scheme   )
 call mpas_pool_get_config(configs,'config_len_disp'       ,len_disp          )
 call mpas_pool_get_array(mesh,'meshDensity',meshDensity)

 call mpas_pool_get_array(diag_physics,'hpbl'    ,hpbl    )
 call mpas_pool_get_array(sfc_input   ,'xland'   ,xland   )

!inout variable

 do j = jts,jte
 do i = its,ite
    !input variables:
    dx_p(i,j)     = len_disp / meshDensity(i)**0.25
    hpbl_p(i,j)   = hpbl(i)
    xland_p(i,j)  = xland(i)

    !output variables:
    q2_p(i,j)     = 0._RKIND
    t2m_p(i,j)    = 0._RKIND
    th2m_p(i,j)   = 0._RKIND
    u10_p(i,j)    = 0._RKIND
    v10_p(i,j)    = 0._RKIND
    hfx_spr_p(i,j)= 0._RKIND
    lh_spr_p(i,j) = 0._RKIND
    sgsgustcu_p(i,j) = 0.0_RKIND

    !output variables (optional):
    cd_p(i,j)     = 0._RKIND
    cda_p(i,j)    = 0._RKIND
    ck_p(i,j)     = 0._RKIND
    cka_p(i,j)    = 0._RKIND
    ustm_p(i,j)   = ustm(i)
 enddo
 enddo

 end subroutine aerosol_from_MPAS

!=================================================================================================================
 subroutine aerosol_to_MPAS(configs,sfc_input,diag_physics,its,ite)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: sfc_input
 integer,intent(in):: its,ite

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics

!local variables:
 integer:: i,j

!local pointers:
 logical,pointer:: config_frac_seaice
 character(len=StrKIND),pointer:: sfclayer_scheme

 real(kind=RKIND),dimension(:),pointer:: br,cpm,chs,chs2,cqs2,flhc,flqc,gz1oz0,hfx,qfx,  &
                                         qgh,qsfc,lh,mol,psim,psih,regime,rmol,ust,wspd, &
                                         znt,zol
 real(kind=RKIND),dimension(:),pointer:: q2,t2m,th2m,u10,v10
 real(kind=RKIND),dimension(:),pointer:: cd,cda,ck,cka,ustm
 real(kind=RKIND),dimension(:),pointer:: xice

!local pointers specific to monin_obukhov:
 real(kind=RKIND),dimension(:),pointer:: fh,fm

!local pointers specific to mynn:
 real(kind=RKIND),dimension(:),pointer:: ch,qcg

 real(kind=RKIND),dimension(:),pointer:: lh_spr,hfx_spr

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_aerosol_monan_scheme',sfclayer_scheme   )

!inout variables:
 call mpas_pool_get_array(diag_physics,'br'    ,br    )
 call mpas_pool_get_array(diag_physics,'cpm'   ,cpm   )


 do j = jts,jte
 do i = its,ite
    !inout variables:
    br(i)     = br_p(i,j)
    cpm(i)    = cpm_p(i,j)
  
 enddo
 enddo

 
 end subroutine aerosol_to_MPAS

!=================================================================================================================
 subroutine init_aerosol_monan(configs,diag_physics,state,mesh,time_lev)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: state
 type(mpas_pool_type),intent(in):: mesh
 integer,intent(in):: time_lev

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics

!local variables and pointers:
 !logical, parameter:: allowed_to_read = .false. !actually not used in subroutine sfclayinit.
 !character(len=StrKIND),pointer:: sfclayer_scheme
 integer,pointer:: nCellsSolve,nVertLevels
 integer,pointer:: num_aer_monan,num_gen_aer_monan,num_dust_aer_monan
 integer,pointer:: num_gas_monan, num_tracer_monan
 integer,pointer:: config_aer_monan, config_dust_aer_monan 
 integer,pointer:: index_aer01,index_aer02
 integer,pointer:: index_dst01,index_dst02,index_dst03,index_dst04
 integer,pointer,dimension(:) :: aer_ind_monan
 
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars
 real(kind=RKIND),dimension(:,:)  ,pointer:: dst01,dst02,dst03,dst04
 real(kind=RKIND),dimension(:,:)  ,pointer:: aer01,aer02

!CCPP-compliant flags:
 character(len=StrKIND):: errmsg
 integer:: errflg,ind

!-----------------------------------------------------------------------------------------------------------------

!initialization of CCPP-compliant flags:
 errmsg = ' '
 errflg = 0
 call mpas_pool_get_dimension(mesh,'nCellsSolve'       ,nCellsSolve        )
 call mpas_pool_get_dimension(mesh,'nVertLevels'       ,nVertLevels        )
 call mpas_pool_get_dimension(mesh,'num_aer_monan'     , num_aer_monan     )
 call mpas_pool_get_dimension(mesh,'num_dust_aer_monan', num_dust_aer_monan)
 call mpas_pool_get_dimension(mesh,'num_gen_aer_monan' , num_gen_aer_monan )
 call mpas_pool_get_dimension(mesh,'num_gas_monan'     , num_gas_monan     )
 call mpas_pool_get_dimension(mesh,'num_tracer_monan'  , num_tracer_monan  )


 call mpas_pool_get_config(configs,'config_dust_aer_monan' ,config_dust_aer_monan )

 call mpas_pool_get_dimension(state,'index_aer01',index_aer01)
 call mpas_pool_get_dimension(state,'index_aer02',index_aer02)
 call mpas_pool_get_dimension(state,'index_dst01',index_dst01)
 call mpas_pool_get_dimension(state,'index_dst02',index_dst02)
 call mpas_pool_get_dimension(state,'index_dst03',index_dst03)
 call mpas_pool_get_dimension(state,'index_dst04',index_dst04)

 call mpas_pool_get_array(diag_physics,'aer_ind_monan',aer_ind_monan)

 ind = 0
 ind=ind+1 ; aer_ind_monan(ind) = index_aer01
 ind=ind+1 ; aer_ind_monan(ind) = index_aer02
 if(config_dust_aer_monan == 1) then 
     ind=ind+1 ; aer_ind_monan(ind) = index_dst01
     ind=ind+1 ; aer_ind_monan(ind) = index_dst02
     ind=ind+1 ; aer_ind_monan(ind) = index_dst03
     ind=ind+1 ; aer_ind_monan(ind) = index_dst04
 endif
 if(ind .ne. num_aer_monan) then 
            call mpas_log_write('check number of aerosols')

            stop " check number of aerosols"
 !print*,'ind=', aer_ind_monan
 endif 

 !dst01 => scalars(index_dst01,:,:)

 end subroutine init_aerosol_monan

!=================================================================================================================
 subroutine driver_aerosol_monan(itimestep,configs,mesh,diag_physics,sfc_input,its,ite, &
                                 state,time_lev)
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

 integer,pointer:: index_aer01,index_aer02
 real(kind=RKIND),dimension(:,:)  ,pointer:: aer01,aer02
 integer,pointer:: index_dst01,index_dst02,index_dst03,index_dst04
 real(kind=RKIND),dimension(:,:)  ,pointer:: dst01,dst02,dst03,dst04

 real(kind=RKIND),dimension(:,:)  ,pointer:: aer_src,aer_snk
 real(kind=RKIND),dimension(:,:)  ,pointer:: zz, zgrid
 real(kind=RKIND),dimension(:)    ,pointer:: xland
 integer,dimension(:),pointer:: isltyp,ivgtyp
!local variables:
 integer:: initflag,iCell,k, n_aer
 real(kind=RKIND):: dx, flux

!CCPP-compliant flags:
 character(len=StrKIND):: errmsg
 integer:: errflg

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write('')
!call mpas_log_write('--- enter subroutine driver_aerosol_monan:')

!initialization of CCPP-compliant flags:
 errmsg = ' '
 errflg = 0
 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)
 call mpas_pool_get_array(mesh,'areaCell',areaCell)
 call mpas_pool_get_array(mesh,'zz'   ,zz   ) ! d(zeta)/dz, vertical metric term, unitless
 call mpas_pool_get_array(mesh,'zgrid',zgrid) ! geometric height of layer interfaces, m MSL
 call mpas_pool_get_array(sfc_input,'xland' ,xland  )
 call mpas_pool_get_array(sfc_input,'ivgtyp',ivgtyp )
 call mpas_pool_get_array(sfc_input,'isltyp',isltyp )

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
 
!copy all MPAS arrays to rectanguler grid:
! call sfclayer_from_MPAS(configs,mesh,diag_physics,sfc_input,its,ite, &
!                         state,time_lev)

 dx = sqrt(maxval(areaCell))
 
 call mpas_pool_get_array(state,'scalars',scalars,time_lev)

 initflag = 1
 if(config_do_restart .or. itimestep > 1) initflag = 0

 !print*," inside driver_aerosol_monan", config_aer_monan, dx
 
 !---- generic aerosols section 
 call mpas_pool_get_dimension(state,'index_aer01',index_aer01); aer01 => scalars(index_aer01,:,:)
 call mpas_pool_get_dimension(state,'index_aer02',index_aer02); aer02 => scalars(index_aer02,:,:)
 
 call mpas_pool_get_array(diag_physics,'aer_ind_monan',aer_ind_monan)

 do iCell = 1, nCellsSolve
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
 
 if(config_dust_aer_monan == 1) then 
 !---- dust aerosols section
     call mpas_pool_get_dimension(state,'index_dst01',index_dst01); dst01 => scalars(index_dst01,:,:)
     call mpas_pool_get_dimension(state,'index_dst02',index_dst02); dst02 => scalars(index_dst02,:,:)
     call mpas_pool_get_dimension(state,'index_dst03',index_dst03); dst03 => scalars(index_dst03,:,:)
     call mpas_pool_get_dimension(state,'index_dst04',index_dst04); dst04 => scalars(index_dst04,:,:)
     !print*,'index dust', index_dst01,index_dst04,maxval(dst01),minval(dst01),maxval(dst04),minval(dst04)
 endif


!copy local arrays to MPAS grid:
! call sfclayer_to_MPAS(configs,sfc_input,diag_physics,its,ite)

!call mpas_log_write('--- end subroutine driver_aerosol_monan.')

 end subroutine driver_aerosol_monan

!=================================================================================================================
 end module mpas_atmphys_driver_aer_monan
!=================================================================================================================
