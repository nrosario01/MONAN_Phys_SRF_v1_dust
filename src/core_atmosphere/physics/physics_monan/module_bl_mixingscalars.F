module module_pbl_scalars
   use mpas_kind_types,only: kind_phys => RKIND
   use mpas_pool_routines
   !use mpas_atmphys_vars
   implicit none

   private
   public :: driver_pbl_scalars
   contains
   !
   !================================================================================================================= 
   ! A routine to calculate the dry mixing in PBL for scalars not included in the MYNN routine.
   ! This is heavily based on "mynn_mix_chem" routine in the file 'bl_mynn_subroutines.F'
   ! Saulo R. Freitas - 19 Sep 2024 (send comments to saulo.freitas@inpe.br)
   ! SRF: Included treatment for aerosols @ 07-Aug-2025 
   !================================================================================================================= 
   subroutine driver_pbl_scalars (itimestep,diag_physics,tend_physics,configs,state,mesh,time_lev &
                                 ,its,ite,kts,kte,dt_pbl)
   !=================================================================================================================     
      implicit none
      !--- input arguments
      type(mpas_pool_type),intent(in) :: configs
      type(mpas_pool_type),intent(in) :: diag_physics
      type(mpas_pool_type),intent(in) :: tend_physics
      type(mpas_pool_type),intent(in) :: state
      type(mpas_pool_type),intent(in) :: mesh
      integer,intent(in) :: itimestep,time_lev
      integer,intent(in) :: its,ite,kts,kte
      real, intent(in) :: dt_pbl
      
      !--- local vars
      integer, parameter :: max_nscalars = 100
      integer :: i,j,k,init_flag,nscalars,isc
      integer, dimension(max_nscalars) :: index_scalar
      real(kind=kind_phys):: dzk
      real(kind=kind_phys),dimension(kts:kte):: dz,rho,dfh
      
      !--- local pointers:
      logical,pointer:: config_do_restart
      integer,pointer:: num_scalars, index_cnvcf, index_buoyx,index_qv

      integer,pointer:: config_aer_monan
      integer,pointer:: num_gen_aer_monan,num_dust_aer_monan
      integer,pointer:: num_aer_monan => null()
      integer,pointer:: num_gas_monan, num_tracer_monan
      integer,pointer,dimension(:) :: aer_ind_monan
      real(kind=RKIND),dimension(:,:,:),pointer:: raerblten

      !--- kzh ="vertical diffusion coefficient of potential temperature" 
      !--- units="m^{2} s^{-1}"
      real(kind=RKIND),dimension(:,:),pointer :: kzh
      real(kind=RKIND),dimension(:,:),pointer :: rho_zz
      real(kind=RKIND),dimension(:,:),pointer :: zz
      real(kind=RKIND),dimension(:,:),pointer :: zgrid
      real(kind=RKIND),dimension(:,:),pointer :: qv

      !--- scalars from MPAS
      real(kind=RKIND),dimension(:,:,:),pointer:: scalars

      !--- local array for scalars
      real(kind=RKIND),dimension(:,:), allocatable :: local_scalars

      !--- CCPP-compliant flags:
      character(len=StrKIND):: errmsg
      integer:: errflg, n_aer
      real :: dt_inv

      !----initialization of CCPP-compliant flags:
      errmsg = ' '
      errflg = 0
      dt_inv = 1.0_RKIND/dt_pbl

      call mpas_pool_get_config(configs,'config_do_restart'  ,config_do_restart  )
      init_flag = 1
      if(config_do_restart .or. itimestep > 1) init_flag = 0

      !--- determine the number and which scalars will be mixed up
      nscalars = 0

      call mpas_pool_get_dimension(state,'index_buoyx',index_buoyx)
      if(index_buoyx > 0) then 
         nscalars = nscalars + 1
         index_scalar(nscalars) = index_buoyx
      endif
      call mpas_pool_get_dimension(state,'index_cnvcf',index_cnvcf)
      if(index_cnvcf > 0) then 
         nscalars = nscalars + 1
         index_scalar(nscalars) = index_cnvcf
      endif
      
      !---- monan aerosol section 
      call mpas_pool_get_config(configs, 'config_aer_monan', config_aer_monan)
      if(config_aer_monan == 1) then
          call mpas_pool_get_array(diag_physics,'aer_ind_monan', aer_ind_monan)
          call mpas_pool_get_dimension(mesh, 'num_aer_monan', num_aer_monan)
          do n_aer = 1, num_aer_monan
                nscalars = nscalars + 1
                index_scalar(nscalars) = aer_ind_monan(n_aer)
                !print*,'naer', n_aer, num_aer_monan, index_scalar(nscalars),nscalars
          enddo
          !--- tendencies 
          call mpas_pool_get_array(tend_physics,'raerblten',raerblten)
      endif
      !
      !--- if no scalar to be mixed, make a "U-turn" 
      if(nscalars == 0) return 

      !--- allocate array to store scalars that will be mixed 
      if(.not.allocated(local_scalars)) allocate(local_scalars(kts:kte,nscalars))

      !--- get the scalars
      call mpas_pool_get_array(state,'scalars',scalars,time_lev)

      !--- get vertical grid length and the air density
      call mpas_pool_get_array(state,'rho_zz' ,rho_zz ,time_lev)
      call mpas_pool_get_array(mesh, 'zz'     ,zz     )
      call mpas_pool_get_array(mesh, 'zgrid'  ,zgrid  )

      call mpas_pool_get_dimension(state,'index_qv',index_qv)
      qv => scalars(index_qv,:,:)

      !--- get vertical diffusivity array
      call mpas_pool_get_array(diag_physics,'kzh'  ,kzh  )

      do i = its,ite

          !--- copy scalars for the local array
          do isc = 1,nscalars
              local_scalars(kts:kte,isc) = scalars(index_scalar(isc),kts:kte,i)
          enddo 
          
          !zz_p(i,k,j)  = zz(k,i)
          !rho_p(i,k,j) = zz(k,i) * rho_zz(k,i)
          !rho_p(i,k,j) = rho_p(i,k,j)*(1._RKIND + qv_p(i,k,j))
          !dz_p(i,k,j)   = zgrid(k+1,i)-zgrid(k,i)

          do k = kts,kte
             dz(k)  =  zgrid(k+1,i)-zgrid(k,i)
             rho(k) =  zz(k,i) * rho_zz(k,i)*(1._RKIND + qv(k,i))
          enddo

          !---retrieve 'dfh' from "kzh" 
          dfh(kts) = 0.
          do k = kts+1,kte
              dzk    = 0.5*(dz(k)+dz(k-1))
              dfh(k) = kzh(k,i)/dzk
          enddo

          !---call to subroutine to apply PBL mixing
          !--- original routine
          !call mix_scalars    (kts,kte,dt_pbl,dz,rho,dfh,nscalars,local_scalars)
          !
          !--- optimized routine
          call mix_scalars_opt(kts,kte,dt_pbl,dz,rho,dfh,nscalars,local_scalars)

          !--- copy back the updated scalars
          do isc = 1, nscalars - num_aer_monan
              scalars(index_scalar(isc),kts:kte,i) = local_scalars(kts:kte,isc)  
          enddo
          
          !--- aerosol section - producing the tendencies instead of updating the scalar mixing ratio
          !--- (a tendency could be computed here by saving the previous values of the scalars)
          n_aer = 0
          do isc = nscalars - num_aer_monan + 1,nscalars
              n_aer = n_aer + 1
              !if(i==its) print*,'check index',index_scalar(isc),aer_ind_monan(n_aer)
              raerblten(n_aer,kts:kte,i) = dt_inv*(local_scalars(kts:kte,isc) - scalars(index_scalar(isc),kts:kte,i))
          enddo 

      enddo ! loop i

      !do n_aer = 1, num_aer_monan
      !         print*,'naer',   aer_ind_monan(n_aer), maxval(scalars( aer_ind_monan(n_aer),:,:)),minval(scalars( aer_ind_monan(n_aer),:,:))
      !enddo

      if(allocated(local_scalars)) deallocate(local_scalars)

   end subroutine driver_pbl_scalars
!-----------------------------------------------------------------------------

   subroutine mix_scalars(kts,kte,delt,dz,rho,dfh,nscalars,local_scalars)
     implicit none
     integer, intent(in) :: kts,kte,nscalars
     real(kind_phys), intent(in) :: delt
     real(kind_phys), dimension(:), intent(in) :: dfh,dz,rho
     real(kind_phys), dimension(:,:), intent(inout) :: local_scalars

!--- local vars
     integer                               :: k,isc
     real(kind_phys)                       :: dztop
     real(kind_phys), dimension(kts:kte)   :: dtz
     real(kind_phys), dimension(kts:kte)   :: a,b,c,d,x
     real(kind_phys), dimension(kts:kte)   :: rhoinv
     real(kind_phys), dimension(kts:kte+1) :: rhoz,khdz

     dztop = 0.5*(dz(kte)+dz(kte-1))

     do k = kts,kte
        dtz(k)=delt/dz(k)
     enddo
     
     !---  Prepare "constants" for diffusion equation.
     !---  khdz = rho*Kh/dz = rho*dfh
     rhoz(kts)  =rho(kts)
     rhoinv(kts)=1./rho(kts)
     khdz(kts)  =rhoz(kts)*dfh(kts)

     do k = kts+1,kte
        rhoz(k)  =(rho(k)*dz(k-1) + rho(k-1)*dz(k))/(dz(k-1)+dz(k))
        rhoz(k)  =   max(rhoz(k),1e-4)
        rhoinv(k)=1./max( rho(k),1e-4)
        !dzk      = 0.5*(dz(k)+dz(k-1))
        khdz(k)  = rhoz(k)*dfh(k)
     enddo
     rhoz(kte+1) = rhoz(kte)
     khdz(kte+1) = rhoz(kte+1)*dfh(kte)

     do isc = 1,nscalars
        k=kts

        a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
        b(k)=1.+dtz(k)*(khdz(k+1)+khdz(k))*rhoinv(k)
        c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)           
        d(k)=  local_scalars(k,isc) 

        do k=kts+1,kte-1
           a(k) =   -dtz(k)*khdz(k)*rhoinv(k)     
           b(k) = 1.+dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k)
           c(k) =   -dtz(k)*khdz(k+1)*rhoinv(k) 
           d(k) = local_scalars(k,isc) 
        enddo
        ! prescribed value at top
        a(kte) = 0.
        b(kte) = 1.
        c(kte) = 0.
        d(kte) = local_scalars(kte,isc)
        
        call tridiag3(kte,a,b,c,d,x)
        
        do k=kts,kte
           local_scalars(k,isc) = max(1.e-24,x(k))
        enddo
     enddo

   end subroutine mix_scalars
!-----------------------------------------------------------------------------
   subroutine mix_scalars_opt(kts,kte,delt,dz,rho,dfh,nscalars,local_scalars)
     implicit none
     integer, intent(in) :: kts,kte,nscalars
     real(kind_phys), intent(in) :: delt
     real(kind_phys), dimension(:), intent(in) :: dfh,dz,rho
     real(kind_phys), dimension(:,:), intent(inout) :: local_scalars

!--- local vars
     integer                               :: k,isc
     real(kind_phys)                       :: dztop
     real(kind_phys), dimension(kts:kte)   :: dtz
     real(kind_phys), dimension(kts:kte)   :: a,b,c,d,x
     real(kind_phys), dimension(kts:kte)   :: rhoinv
     real(kind_phys), dimension(kts:kte+1) :: rhoz,khdz

     dztop = 0.5*(dz(kte)+dz(kte-1))

     do k = kts,kte
        dtz(k)=delt/dz(k)
     enddo
     
     !---  Prepare "constants" for diffusion equation.
     !---  khdz = rho*Kh/dz = rho*dfh
     rhoz(kts)  =rho(kts)
     rhoinv(kts)=1./rho(kts)
     khdz(kts)  =rhoz(kts)*dfh(kts)

     do k = kts+1,kte
        rhoz(k)  =(rho(k)*dz(k-1) + rho(k-1)*dz(k))/(dz(k-1)+dz(k))
        rhoz(k)  =   max(rhoz(k),1e-4)
        rhoinv(k)=1./max( rho(k),1e-4)
        khdz(k)  = rhoz(k)*dfh(k)
     enddo
     rhoz(kte+1) = rhoz(kte)
     khdz(kte+1) = rhoz(kte+1)*dfh(kte)

     k=kts
     a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
     b(k)=1.+dtz(k)*(khdz(k+1)+khdz(k))*rhoinv(k)
     c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)           

     do k=kts+1,kte-1
           a(k) =   -dtz(k)*khdz(k)*rhoinv(k)     
           b(k) = 1.+dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k)
           c(k) =   -dtz(k)*khdz(k+1)*rhoinv(k) 
     enddo
     ! prescribed value at top
     a(kte) = 0.
     b(kte) = 1.
     c(kte) = 0.
        
     do isc = 1,nscalars
        d(kts:kte) = local_scalars(kts:kte,isc)

        call tridiag3(kte,a,b,c,d,x)
        
        do k=kts,kte
           local_scalars(k,isc) = max(1.e-24,x(k))
        enddo
     enddo

   end subroutine mix_scalars_opt
!-----------------------------------------------------------------------------
!>\ingroup gsd_mynn_edmf
   subroutine tridiag3(kte,a,b_in,c,d,x)
       !ccccccccccccccccccccccccccccccc                                                                   
       ! Aim: Inversion and resolution of a tridiagonal matrix                                            
       !          A X = D                                                                                 
       ! Input:                                                                                           
       !  a(*) lower diagonal (Ai,i-1)                                                                  
       !  b(*) principal diagonal (Ai,i)                                                                
       !  c(*) upper diagonal (Ai,i+1)                                                                  
       !  d                                                                                               
       ! Output                                                                                           
       !  x     results                                                                                   
       !ccccccccccccccccccccccccccccccc                                                                   
        implicit none
        integer,intent(in)   :: kte
        integer, parameter   :: kts=1
        real(kind_phys), dimension(:), intent(in   ) :: a,b_in,c
        real(kind_phys), dimension(:), intent(inout) :: d
        real(kind_phys), dimension(:), intent(out  ) :: x
        integer :: ind
        real :: b(kts:kte)

        b = b_in
        do ind=kte-1,kts,-1
         d(ind)=d(ind)-c(ind)*d(ind+1)/b(ind+1)
         b(ind)=b(ind)-c(ind)*a(ind+1)/b(ind+1)
        enddo

        do ind=kts+1,kte
         d(ind)=d(ind)-a(ind)*d(ind-1)/b(ind-1)
        enddo

        do ind=kts,kte
         x(ind)=d(ind)/b(ind)
        enddo

   end subroutine tridiag3
!-----------------------------------------------------------------------------
end module module_pbl_scalars
