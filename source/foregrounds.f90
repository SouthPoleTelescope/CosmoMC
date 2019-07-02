!2011/05/26: changed templates to Dls

module foregrounds
  use cmbtypes

  implicit none

  private 

  public :: foreground_params, utindex, InitForegroundData, cl_foreground, &
       GetForegroundParamsFromArray, GetAlphaPrior, GetCirrusFactor, &
       GetRadioAmpUncPrior, GetRadioAmpPrior, InitRadioAmpPrior, &
       GetRadioClAmpUncPrior, GetRadioClAmpPrior, InitRadioClAmpPrior, &
       cl_cib_foreground,cosmo_scale_ksz,cosmo_scale_tsz,index_tsz,index_ksz, &
       index_dg_po,index_dg_cl,index_cirrus,index_rg_po, &
       HaveForegroundsBeenInitialized,setForegroundsUninitialized,Bnu,&
       read_dl_template,read_cl_template,nForegroundParams,&
       getForegroundPriorLnL,InitFGModel,printForegrounds,&
       OpenReadBinaryFile,OpenWriteBinaryFile,OpenReadBinaryStreamFile,calFactorsToCIB,MaxNFreq
  
!set in settings.f90  

  integer, parameter :: MaxNFreq = 6
  integer, parameter :: NDecorrel = (MaxNFreq-1)*(MaxNFreq)/2
  integer, parameter :: nForegroundParams=39+NDecorrel
  integer, dimension(NDecorrel,NDecorrel) :: dmatrix_index
  logical :: cosmological_scaling_ksz
  logical :: cosmological_scaling_tsz,ApplyCirrusPrior90,ApplyCirrusPrior150,ApplyCirrusPrior220
  logical :: single_clustered_freq_scaling, combine_spire_for_tszcib_correlation 
  logical :: only1HaloTszCib, CIB_decorrelation_matrix
  logical :: tSZ_CIB_logFreq
  logical :: ShangModelCorrelationShape,applyCIBCalToCirrus,calFactorsToCIB
  logical :: use_decorrelation_matrix_form, use_sigma
  integer, parameter :: ntsz=1,nksz=2,ncirrus=1
  integer, parameter:: index_tsz=1, index_ksz = index_tsz+ntsz, &
       index_dg_po=index_ksz+nksz, &
       index_dg_cl=index_dg_po+2,&
       index_cirrus=index_dg_cl+4,&
       index_rg_po = index_cirrus + ncirrus
  type foreground_params
     !all czero's are D_3000 for that signal
     real(mcp) czero_tsz
     real(mcp) czero_ksz
     real(mcp) czero_ksz2
     real(mcp) czero_dg_po
     real(mcp) czero_dg_po_spire
     real(mcp) czero_dg_cl
     real(mcp) czero_dg_cl2
     real(mcp) czero_dg_cl_spire
     real(mcp) czero_dg_cl2_spire
     real(mcp) czero_cirrus
     real(mcp) czero_rg_po
     real(mcp) czero_rg_cl

     real(mcp) T_dg_po
     real(mcp) beta_dg_po
     real(mcp) sigmasq_dg_po
     real(mcp) T_spire
     real(mcp) beta_spire
     real(mcp) sigmasq_spire

     real(mcp) T_dg_cl
     real(mcp) beta_dg_cl
     real(mcp) sigmasq_dg_cl

     real(mcp) T_dg_cl2
     real(mcp) beta_dg_cl2
     real(mcp) sigmasq_dg_cl2

     real(mcp) alpha_rg
     real(mcp) sigmasq_rg

     real(mcp) T_cirrus
     real(mcp) beta_cirrus
     
     real(mcp) tsz_dg_cor_const
     real(mcp) tsz_dg_cor_lin_1e4
     real(mcp) tsz_dg_cor_quad_1e7

     real(mcp) tsz_spire_cor_const
     real(mcp) tsz_spire_cor_lin_1e4

     real(mcp) tsz_rg_cor
     real(mcp) dg_cl_ell_power
     real(mcp) cib_upturn_100ghz

     real(mcp) ksz_slope
     real(mcp) tsz_cib_slope

     real(mcp) decorrel_slope
     real(mcp) decorrel_matrix(NDecorrel)


  end type foreground_params

  real(mcp), parameter :: d3000 = 3000*3001/(2*pi)

  ! Kinetic and thermal SZ effect templates
  real(mcp), dimension(2:lmax) :: ksz_templ, tsz_templ, ksz2_templ
  real(mcp), dimension(2:lmax) :: clust_dg_templ,clust2_dg_templ
  real(mcp), dimension(2:lmax) :: cirrus_templ
  real(mcp), dimension(2:lmax) :: clust_rg_templ
  real(mcp), dimension(2:lmax) :: l_divide_3000

  real(mcp),dimension(2:lmax) :: cls_diff
  
  logical :: SuccessfulInitialization 
  
  real(mcp) DelAlphaPrior,CirrusFactorPrior,AddAlphaPoisson,AddAlphaSpire
  real(mcp) radio_amp,radio_unc,radio_cl_amp,radio_cl_unc
  real(mcp) HFIPoissonScale
contains
  subroutine printForegrounds(fgs)
    type(foreground_params),intent(in):: fgs
    
     print*,'czero_tsz',fgs%czero_tsz
     print*,'czero_ksz',fgs%czero_ksz

     print*,'czero_ksz2',fgs%czero_ksz2
     print*,'czero_dg_po',fgs%czero_dg_po
     print*,'czero_dg_po_spire',fgs%czero_dg_po_spire
     print*,'czero_dg_cl',fgs%czero_dg_cl
     print*,'czero_dg_cl2',fgs%czero_dg_cl2
     print*,'czero_dg_cl_spire',fgs%czero_dg_cl_spire
     print*,'czero_dg_cl2_spire',fgs%czero_dg_cl2_spire
     print*,'czero_cirrus',fgs%czero_cirrus
     print*,'czero_rg_po',fgs%czero_rg_po
     print*,'czero_rg_cl',fgs%czero_rg_cl

     print*,'T_dg_po',fgs%T_dg_po
     print*,'beta_dg_po',fgs%beta_dg_po
     print*,'sigmasq_dg_po',fgs%sigmasq_dg_po
     print*,'T_spire',fgs%T_spire
     print*,'beta_spire',fgs%beta_spire
     print*,'sigmasq_spire',fgs%sigmasq_spire
     print*,'T_dg_cl',fgs%T_dg_cl
     print*,'beta_dg_cl',fgs%beta_dg_cl
     print*,'sigmasq_dg_cl',fgs%sigmasq_dg_cl
     print*,'T_dg_cl2',fgs%T_dg_cl2
     print*,'beta_dg_cl2',fgs%beta_dg_cl2
     print*,'sigmasq_dg_cl2',fgs%sigmasq_dg_cl2
     print*,'alpha_rg',fgs%alpha_rg
     print*,'sigmasq_rg',fgs%sigmasq_rg
     print*,'T_cirrus',fgs%T_cirrus
     print*,'beta_cirrus',fgs%beta_cirrus

     print*,'tsz_dg_cor_const',fgs%tsz_dg_cor_const
     print*,'tsz_dg_cor_lin_1e4',fgs%tsz_dg_cor_lin_1e4
     print*,'tsz_dg_cor_quad_1e7',fgs%tsz_dg_cor_quad_1e7
     print*,'tsz_spire_cor_const',fgs%tsz_spire_cor_const
     print*,'tsz_spire_cor_lin_1e4',fgs%tsz_spire_cor_lin_1e4
     print*,'tsz_rg_cor',fgs%tsz_rg_cor
     print*,'dg_cl_ell_power',fgs%dg_cl_ell_power
     print*,'cib_upturn_100ghz',fgs%cib_upturn_100ghz
     print*,'ksz slope', fgs%ksz_slope
     print*,'tsz cib slope', fgs%tsz_cib_slope

     print*,'cib decorrelation slope',fgs%decorrel_slope
     print*,'cib decor matrix:',fgs%decorrel_matrix
  end subroutine printForegrounds

  subroutine setForegroundsUninitialized()
    SuccessfulInitialization=.false.
  end subroutine setForegroundsUninitialized
  function HaveForegroundsBeenInitialized()
    logical HaveForegroundsBeenInitialized
    HaveForegroundsBeenInitialized = SuccessfulInitialization
  end function HaveForegroundsBeenInitialized

  !
  ! Returns in cl_arr the matrix
  !       (TT, TE, TB)
  !       (ET, EE, EB)
  !       (BT, BE, BB)
  !
  ! corresponding to the point source cross/auto power spectrums
  ! at the given frequency and at one single l-value
  !
  ! the logical array which turns on/off different componenets. the order is
  !   (dutsy clustered, dusty poisson, radio poisson, ksz, tsz)
  !
  ! the array bytype returns the TT components in that same order
  !
  function cl_foreground(params,ifr,jfr,nfr,eff_fr,norm_fr,cib_cals,component_spectra)
    type(foreground_params) :: params
    integer :: i,ifr,jfr,nfr
    real(mcp), dimension(5,nfr) :: eff_fr
    real(mcp), dimension(nfr) :: cib_cals
    real(mcp), dimension(5) :: norm_fr
    real(mcp),dimension(2:lmax) :: cl_foreground
    real(mcp) :: fri,frj,norm,fr0,frqdep
    real(mcp),dimension(2:lmax) :: cl_dg_po, cl_dg_cl, cl_rg_po, cl_k_sz, &
         cl_t_sz, cl2_cirrus
    real(mcp),dimension(2:lmax) :: cl_tsz_rg_cor, cl_tsz_dgcl_cor
    real(mcp),dimension(2:lmax) :: cl_dg, cl_rg,cl_spire
    real(mcp),dimension(2:lmax) :: cltmpi,cltmpj
    real(mcp),dimension(2:lmax,7), optional, intent(out) :: component_spectra
    real(mcp) :: cib_cal_i,cib_cal_j

    cib_cal_i = cib_cals(ifr)
    cib_cal_j = cib_cals(jfr)

    cl_dg_cl(:) = (cib_cal_i*cib_cal_j)*cl_dusty_clustered(params,eff_fr(1,ifr),eff_fr(1,jfr),norm_fr(1),ifr,jfr)
    cl_dg_po(:) = (cib_cal_i*cib_cal_j)*cl_dusty_poisson(params,eff_fr(2,ifr),eff_fr(2,jfr),norm_fr(2),ifr,jfr)
    cl_rg_po(:) = cl_radio(params,eff_fr(3,ifr),eff_fr(3,jfr),norm_fr(3))
    cl_k_sz(:) = cl_ksz(params)
    cl_t_sz(:) = cl_tsz(params,eff_fr(5,ifr),eff_fr(5,jfr),norm_fr(5))
    cl_spire(:) = (cib_cal_i*cib_cal_j)*cl_dusty_spire(params,eff_fr(2,ifr),eff_fr(2,jfr))
    !TSZ-Dusty correlation

     if (combine_spire_for_tszcib_correlation) then 
        if (only1HaloTszCib) then
           cltmpi = cl_dusty_clustered(params,eff_fr(1,ifr),eff_fr(1,ifr),norm_fr(1),ifr,jfr,.true.) + &
                cl_dusty_clustered1_spire(params,eff_fr(2,ifr),eff_fr(2,ifr))
           cltmpj = cl_dusty_clustered(params,eff_fr(1,jfr),eff_fr(1,jfr),norm_fr(1),ifr,jfr,.true.) + &
                cl_dusty_clustered1_spire(params,eff_fr(2,jfr),eff_fr(2,jfr))
        else
           cltmpi = cl_dusty_clustered(params,eff_fr(1,ifr),eff_fr(1,ifr),norm_fr(1),ifr,jfr) + &
                cl_dusty_poisson(params,eff_fr(2,ifr),eff_fr(2,ifr),norm_fr(2),ifr,jfr) + &
                cl_dusty_spire(params,eff_fr(2,ifr),eff_fr(2,ifr))
           cltmpj = cl_dusty_clustered(params,eff_fr(1,jfr),eff_fr(1,jfr),norm_fr(1),ifr,jfr) + &
                cl_dusty_poisson(params,eff_fr(2,jfr),eff_fr(2,jfr),norm_fr(2),ifr,jfr) + &
                cl_dusty_spire(params,eff_fr(2,jfr),eff_fr(2,jfr))
        end if
        cltmpi = cltmpi*(cib_cal_i*cib_cal_i)
        cltmpj = cltmpj*(cib_cal_j*cib_cal_j)

       cl_tsz_dgcl_cor(:) = tsz_dgcl_cor(params) * &
            ( tsz_cib(params,eff_fr(2,jfr)) * sqrt( cl_tsz(params,eff_fr(5,ifr),eff_fr(5,ifr),norm_fr(5)) * cltmpj) + &
            tsz_cib(params,eff_fr(2,ifr)) * sqrt( cl_tsz(params,eff_fr(5,jfr),eff_fr(5,jfr),norm_fr(5)) * cltmpi))
    else
       if (only1HaloTszCib) then
          cltmpi = cl_dusty_clustered(params,eff_fr(1,ifr),eff_fr(1,ifr),norm_fr(1),ifr,jfr,.true.) 
          cltmpj = cl_dusty_clustered(params,eff_fr(1,jfr),eff_fr(1,jfr),norm_fr(1),ifr,jfr,.true.) 
       else
          cltmpi = cl_dusty_clustered(params,eff_fr(1,ifr),eff_fr(1,ifr),norm_fr(1),ifr,jfr) + &
               cl_dusty_poisson(params,eff_fr(2,ifr),eff_fr(2,ifr),norm_fr(2),ifr,jfr) 
          cltmpj = cl_dusty_clustered(params,eff_fr(1,jfr),eff_fr(1,jfr),norm_fr(1),ifr,jfr) + &
               cl_dusty_poisson(params,eff_fr(2,jfr),eff_fr(2,jfr),norm_fr(2),ifr,jfr) 
       endif
       cltmpi = cltmpi*(cib_cal_i*cib_cal_i)
       cltmpj = cltmpj*(cib_cal_j*cib_cal_j)
       cl_tsz_dgcl_cor(:) = tsz_dgcl_cor(params) * &
            ( tsz_cib(params,eff_fr(2,jfr)) * sqrt( cl_tsz(params,eff_fr(5,ifr),eff_fr(5,ifr),norm_fr(5)) * cltmpj) + &
            tsz_cib(params,eff_fr(2,ifr)) * sqrt( cl_tsz(params,eff_fr(5,jfr),eff_fr(5,jfr),norm_fr(5)) * cltmpi))
       
       !will be zero so skip if one not positive
       if (eff_fr(2,ifr) > 400 .or. eff_fr(2,jfr) > 400.0) then
          if (only1HaloTszCib) then
             cltmpi = cl_dusty_clustered1_spire(params,eff_fr(2,ifr),eff_fr(2,ifr))
             cltmpj = cl_dusty_clustered1_spire(params,eff_fr(2,jfr),eff_fr(2,jfr))
          else
             cltmpi = cl_dusty_spire(params,eff_fr(2,ifr),eff_fr(2,ifr))
             cltmpj = cl_dusty_spire(params,eff_fr(2,jfr),eff_fr(2,jfr))
          endif
          cltmpi = cltmpi*(cib_cal_i*cib_cal_i)
          cltmpj = cltmpj*(cib_cal_j*cib_cal_j)
          cl_tsz_dgcl_cor(:) = cl_tsz_dgcl_cor(:) + tsz_dgcl_cor(params) * &
               ( tsz_cib_spire(params,eff_fr(2,jfr)) * sqrt( cl_tsz(params,eff_fr(5,ifr),eff_fr(5,ifr),norm_fr(5)) * cltmpj) + &
               tsz_cib_spire(params,eff_fr(2,ifr)) * sqrt( cl_tsz(params,eff_fr(5,jfr),eff_fr(5,jfr),norm_fr(5)) * cltmpi))
       endif
    endif


    !TSZ-Radio correlation
    if (params%tsz_rg_cor .eq. 0) then
       cl_tsz_rg_cor(:) = 0.0
    else
       cl_tsz_rg_cor(:) = params%tsz_rg_cor * tsz_rg_cor() * ( &
            sqrt(cl_tsz(params,eff_fr(5,ifr),eff_fr(5,ifr),norm_fr(5)) * &
            cl_radio(params,eff_fr(1,jfr),eff_fr(1,jfr),norm_fr(1))) &
            + sqrt(cl_tsz(params,eff_fr(5,jfr),eff_fr(5,jfr),norm_fr(5)) * &
            cl_radio(params,eff_fr(1,ifr),eff_fr(1,ifr),norm_fr(1))) )
    endif
    
    norm = 1. /(d3000*cirrus_templ(3000))
    cl2_cirrus = cl_cirrus(params,eff_fr(1,ifr),eff_fr(1,jfr))

    if (applyCIBCalToCirrus) &
         cl2_cirrus = cl2_cirrus(:) * (cib_cal_i*cib_cal_j)
    
    cl_dg = cl_dg_po + cl_dg_cl
    cl_rg = cl_rg_po
    
    if (present(component_spectra)) then 

       component_spectra(:,1) = cl_dg_po
       component_spectra(:,2) = cl_dg_cl
       component_spectra(:,3) = cl_k_sz
       component_spectra(:,4) = cl_t_sz
       component_spectra(:,5) = cl_rg
       component_spectra(:,6) = cl_tsz_dgcl_cor
       component_spectra(:,7) = cl2_cirrus + cl_tsz_rg_cor + cl_spire
    endif

    cl_foreground = cl_dg + cl_rg + cl_t_sz + cl_k_sz + cl2_cirrus + cl_tsz_dgcl_cor + cl_tsz_rg_cor + cl_spire
  end function cl_foreground

  function tsz_cib(params,freq)
    real(mcp) tsz_cib
    type(foreground_params) :: params
    real(mcp) :: freq
    if (tSZ_CIB_logFreq) then 
              tsz_cib =  params%tsz_dg_cor_const + params%tsz_dg_cor_lin_1e4*1d-4 * log(freq) + params%tsz_dg_cor_quad_1e7*1d-7 * log(freq)**2
    else
       tsz_cib =  params%tsz_dg_cor_const + params%tsz_dg_cor_lin_1e4*1d-4 * freq + params%tsz_dg_cor_quad_1e7*1d-7 * freq * freq
    endif
  end function tsz_cib

  function tsz_cib_spire(params,freq)
    real(mcp) :: tsz_cib_spire
    type(foreground_params) :: params
    real(mcp) :: freq
    
    tsz_cib_spire= params%tsz_spire_cor_const + params%tsz_spire_cor_lin_1e4*1d-4 * freq 
  end function tsz_cib_spire

  function cl_cib_foreground(params,eff_fr,norm_fr,ifr,jfr)
    type(foreground_params) :: params
    integer,intent(in) :: ifr,jfr
    real(mcp) :: eff_fr
    real(mcp) :: norm_fr
    real(mcp),dimension(2:lmax) :: cl_cib_foreground
    real(mcp),dimension(2:lmax) :: cl_dg_po, cl_dg_cl,cl_dg_spire
    
    cl_dg_cl = cl_dusty_clustered(params,eff_fr,eff_fr,norm_fr,ifr,jfr)
    cl_dg_po = cl_dusty_poisson(params,eff_fr,eff_fr,norm_fr,ifr,jfr)
    cl_dg_spire = cl_dusty_spire(params,eff_fr,eff_fr)
    
    cl_cib_foreground = cl_dg_po + cl_dg_cl + cl_dg_spire

  end function cl_cib_foreground
  

  function cl_radio(params,fri,frj,fr0)
    real(mcp),dimension(2:lmax)  :: cl_radio
    type(foreground_params) :: params
    real(mcp) :: fri,frj,fr0
   
    cl_radio(:) = params%czero_rg_po/d3000/dBdT(fri,fr0)/dBdT(frj,fr0)*&
         (fri/fr0*frj/fr0)**(params%alpha_rg + &
         log(fri/fr0*frj/fr0)/2 * params%sigmasq_rg)

    if (params%czero_rg_cl .gt. 0) then 
       cl_radio(2:lmax) = cl_radio(2:lmax) + &
            params%czero_rg_cl/d3000/dBdT(fri,fr0)/&
            dBdT(frj,fr0)*&
            (fri/fr0*frj/fr0)**(params%alpha_rg)*&
            clust_rg_templ(:)
    endif
  end function cl_radio

  function cirrus_power3000(params,fri,frj)
    real(mcp)  :: cirrus_power3000
    type(foreground_params) :: params
    real(mcp) :: fri, frj, fr0
    real(mcp) :: frqdep
    
    fr0=220.0
    frqdep = ((fri*frj)/(fr0*fr0))**(params%beta_cirrus)
    frqdep = frqdep *&
         Bnu(fri,fr0,params%T_cirrus)*Bnu(frj,fr0,params%T_cirrus)
    
    frqdep = frqdep /( dBdT(fri,fr0) * dBdT(frj,fr0) )
    cirrus_power3000 = params%czero_cirrus*frqdep
  end function cirrus_power3000

  function cl_cirrus(params,fri,frj)
    real(mcp),dimension(2:lmax)  :: cl_cirrus
    type(foreground_params) :: params
    real(mcp) :: fri, frj, fr0
    real(mcp) :: power

    power = cirrus_power3000(params,fri,frj)
    
    cl_cirrus(:)=(power/(d3000*cirrus_templ(3000))) * cirrus_templ(:)

  end function cl_cirrus

  function cl_dusty_clustered1_spire(params,fri,frj)
    real(mcp),dimension(2:lmax)  :: cl_dusty_clustered1_spire
    type(foreground_params) :: params
    real(mcp) :: fri, frj,  fr0spire
    real(mcp) :: frqdep
    real(mcp) :: ff
    real(mcp) :: effalpha,effsigmasq,effalpha_2

    fr0spire=600.


    if (fri > 400 .and. frj > 400.0) then

       effalpha = params%T_spire + AddAlphaSpire*params%T_dg_po
       effsigmasq = params%sigmasq_spire + AddAlphaSpire*params%sigmasq_dg_po
       effalpha_2 = params%beta_spire + AddAlphaSpire*params%beta_dg_po
       frqdep = ((fri*frj)/(fr0spire*fr0spire))**( &
            effalpha_2 + log(fri/fr0spire*frj/fr0spire)/2 * effsigmasq )
       frqdep = frqdep *&
            Bnu(fri,fr0spire,effalpha)*Bnu(frj,fr0spire,effalpha)

       frqdep = frqdep * 1.0e6 /d3000/dBdT(fri,fr0spire)/dBdT(frj,fr0spire) 

       if (params%dg_cl_ell_power /= 0) then
          cl_dusty_clustered1_spire(:) = cl_dusty_clustered1_spire(:) + &
               frqdep * params%czero_dg_cl_spire * &
               clust_dg_templ(:)/(clust_dg_templ(3000)) * (l_divide_3000(:))**params%dg_cl_ell_power
       else
          cl_dusty_clustered1_spire(:) = cl_dusty_clustered1_spire(:) + &
               frqdep * params%czero_dg_cl_spire * &
               clust_dg_templ(:)/(clust_dg_templ(3000)) 
       endif

    else
       cl_dusty_clustered1_spire(:)=0
    endif

  end function cl_dusty_clustered1_spire

  function cl_dusty_spire(params,fri,frj)
    real(mcp),dimension(2:lmax)  :: cl_dusty_spire
    type(foreground_params) :: params
    real(mcp) :: fri, frj,  fr0spire
    real(mcp) :: frqdep
    real(mcp) :: ff
    real(mcp) :: effalpha,effsigmasq,effalpha_2
 
    fr0spire=600.
    cl_dusty_spire(:)=0

    if (fri > 400 .and. frj > 400.0) then
       
       effalpha = params%T_spire + AddAlphaSpire*params%T_dg_po
       effsigmasq = params%sigmasq_spire + AddAlphaSpire*params%sigmasq_dg_po
       effalpha_2 = params%beta_spire + AddAlphaSpire*params%beta_dg_po
       frqdep = ((fri*frj)/(fr0spire*fr0spire))**( &
            effalpha_2 + log(fri/fr0spire*frj/fr0spire)/2 * effsigmasq )
       frqdep = frqdep *&
            Bnu(fri,fr0spire,effalpha)*Bnu(frj,fr0spire,effalpha)
       
       frqdep = frqdep * 1.0e6 /d3000/dBdT(fri,fr0spire)/dBdT(frj,fr0spire) 
       
       cl_dusty_spire(:) =  frqdep * params%czero_dg_po_spire
       
       if (params%dg_cl_ell_power /= 0) then
          cl_dusty_spire(:) = cl_dusty_spire(:) + &
               frqdep * params%czero_dg_cl_spire * &
               clust_dg_templ(:)/(clust_dg_templ(3000)) * (l_divide_3000(:))**params%dg_cl_ell_power
       else
          cl_dusty_spire(:) = cl_dusty_spire(:) + &
               frqdep * params%czero_dg_cl_spire * &
               clust_dg_templ(:)/(clust_dg_templ(3000)) 
       endif

       if (params%czero_dg_cl2_spire /= 0) then
          cl_dusty_spire(:) = cl_dusty_spire(:)+ &
               clust2_dg_templ(:)/(clust2_dg_templ(3000))*params%czero_dg_cl2_spire*frqdep 
       endif
       
    endif
  end function cl_dusty_spire

  function cl_dusty_poisson(params,fri,frj,fr0,ifr,jfr)
    real(mcp),dimension(2:lmax)  :: cl_dusty_poisson
    type(foreground_params) :: params
    integer, intent(in) :: ifr,jfr
    real(mcp) :: fri, frj, fr0, fr0spire
    real(mcp) :: frqdep
    real(mcp) :: ff,ff1,ff2,decor
    real(mcp) :: effalpha,effsigmasq,effalpha_2

    !basic powerlaw
    frqdep = ((fri*frj)/(fr0*fr0))**( params%beta_dg_po)

    if (use_decorrelation_matrix_form) then
       if (use_sigma) then 
          ff1=(fri/fr0)
          !NB: we get a factor of 4x in the exponent because ff1 isn't squared
          ! we also get a factor of 1/2x because we will want the sqrt.
          !this is why ff has the /2 since there we do square, and don't want the sqrt
          ff1=ff1**(params%sigmasq_dg_po  * log(ff1))
          ff2=(frj/fr0)
          ff2=ff2**(params%sigmasq_dg_po  * log(ff2))
          
          ff = ((fri*frj)/(fr0*fr0))
          ff = ff ** (params%sigmasq_dg_po /2* log(ff2))
          decor=ff/(ff1*ff2)
       else
          decor = params%decorrel_matrix(dmatrix_index(ifr,jfr))
       endif
    else
       !not <1 like decor in the other case...this boosts the autospectra
       decor =  ((fri*frj)/(fr0*fr0))**( &
             log(fri/fr0*frj/fr0)/2 * params%sigmasq_dg_po )

    endif
    frqdep = frqdep*decor
       
    frqdep = frqdep *&
         Bnu(fri,fr0,params%T_dg_po)*Bnu(frj,fr0,params%T_dg_po)
    
    if (fri < 100.) &
         frqdep = frqdep * params%cib_upturn_100ghz
    if (frj < 100.) &
         frqdep = frqdep * params%cib_upturn_100ghz
        
    cl_dusty_poisson(:) = params%czero_dg_po/d3000/dBdT(fri,fr0)/dBdT(frj,fr0)*&
         frqdep


  end function cl_dusty_poisson

  function cl_dusty_clustered(params,fri,frj,fr0,ifr,jfr,only1halo)
    integer, intent(in) :: ifr,jfr
    real(mcp),dimension(2:lmax)  :: cl_dusty_clustered,sloped_decor
    type(foreground_params) :: params
    real(mcp) :: fri,frj,fr0,frqdep,frqdep0,effalpha_cl, effalpha_cl_2,frqdeplin,effsigmasq_cl
    real(mcp) :: effalpha,effsigmasq,effalpha_2,decor,ff,ff1,ff2
    logical, intent(in), optional :: only1halo
    logical :: Want2Halo, hasslope
    hasslope=.false.

    if (present(only1halo)) then
        Want2Halo = .not. only1halo
    else
        Want2Halo = .true.
    end if

    effalpha_cl = params%T_dg_cl + AddAlphaPoisson*params%T_dg_po
    effsigmasq_cl = params%sigmasq_dg_cl + AddAlphaPoisson*params%sigmasq_dg_po
    effalpha_cl_2 = params%beta_dg_cl + AddAlphaPoisson*params%beta_dg_po
    frqdep = ((fri*frj)/(fr0*fr0))**(effalpha_cl_2)

    frqdep = frqdep *&
         Bnu(fri,fr0,effalpha_cl)*Bnu(frj,fr0,effalpha_cl)
    
    frqdep = frqdep / dBdT(fri,fr0) / dBdT(frj,fr0)
    
    if (fri < 100.) &
         frqdep = frqdep * params%cib_upturn_100ghz
    if (frj < 100.) &
         frqdep = frqdep * params%cib_upturn_100ghz

    if (use_decorrelation_matrix_form) then
       decor = 1_mcp
       if (ifr /= jfr) then
          if (use_sigma) then 
             ff1=(fri/fr0)
             !NB: we get a factor of 4x in the exponent because ff1 isn't squared
             ! we also get a factor of 1/2x because we will want the sqrt.
             !this is why ff has the /2 since there we do square, and don't want the sqrt
             ff1=ff1**(params%sigmasq_dg_po  * log(ff1))
             ff2=(frj/fr0)
             ff2=ff2**(params%sigmasq_dg_po  * log(ff2))
             
             ff = ((fri*frj)/(fr0*fr0))
             ff = ff ** (params%sigmasq_dg_po /2* log(ff2))
             decor=ff/(ff1*ff2)
          else
             decor = params%decorrel_matrix(dmatrix_index(ifr,jfr))
          endif
          if (params%decorrel_slope /= 0 ) then
             sloped_decor(:) = (l_divide_3000(:)-1_mcp)*params%decorrel_slope + decor
             
             !limiting range
             where(sloped_decor<0_mcp) &
                  sloped_decor = 0_mcp
             where(sloped_decor>1_mcp) &
                  sloped_decor = 1_mcp
             
             !don't need this since will be set elsewhere
             hasslope=.true.
             decor=1_mcp
          endif
       endif
    else
       !not <1 like decor in the other case...this boosts the autospectra
       decor =  ((fri*frj)/(fr0*fr0))**( &
            log(fri/fr0*frj/fr0)/2 * effsigmasq_cl )
    endif
    frqdep = frqdep*decor


    cl_dusty_clustered(:) = clust_dg_templ(:) * (params%czero_dg_cl*frqdep/(d3000*clust_dg_templ(3000)))
    if (params%dg_cl_ell_power /= 0) &
         cl_dusty_clustered(:) =  cl_dusty_clustered(:) * (l_divide_3000(:))**params%dg_cl_ell_power
    if(hasslope) &
         cl_dusty_clustered(:) =  cl_dusty_clustered(:) * sloped_decor(:)
    
    if (Want2Halo .and. params%czero_dg_cl2 /= 0) then
       if (single_clustered_freq_scaling ) then 
          if (hasslope) then
             cl_dusty_clustered(:) = cl_dusty_clustered(:)+ clust2_dg_templ(:)/(d3000*clust2_dg_templ(3000))*params%czero_dg_cl2*frqdep * sloped_decor(:)
          else
             cl_dusty_clustered(:) = cl_dusty_clustered(:)+ clust2_dg_templ(:)/(d3000*clust2_dg_templ(3000))*params%czero_dg_cl2*frqdep 
          endif
       else
          effalpha_cl = params%T_dg_cl2 + AddAlphaPoisson*params%T_dg_po
          effsigmasq_cl = params%sigmasq_dg_cl2 + AddAlphaPoisson*params%sigmasq_dg_po
          effalpha_cl_2 = params%beta_dg_cl2 + AddAlphaPoisson*params%beta_dg_po
          frqdep = ((fri*frj)/(fr0*fr0))**(effalpha_cl_2)
          frqdep = frqdep *&
               Bnu(fri,fr0,effalpha_cl)*Bnu(frj,fr0,effalpha_cl)
          
          frqdep = frqdep / dBdT(fri,fr0) / dBdT(frj,fr0)
          
          if (fri < 100.) &
               frqdep = frqdep * params%cib_upturn_100ghz
          if (frj < 100.) &
               frqdep = frqdep * params%cib_upturn_100ghz

          !if the below isn't true, reuse values for decor, sloped_decor
          if (use_decorrelation_matrix_form .and. use_sigma .and. (ifr /= jfr)) then
             ff1=(fri/fr0)
             !NB: we get a factor of 4x in the exponent because ff1 isn't squared
             ! we also get a factor of 1/2x because we will want the sqrt.
             !this is why ff has the /2 since there we do square, and don't want the sqrt
             ff1=ff1**(params%sigmasq_dg_po  * log(ff1))
             ff2=(frj/fr0)
             ff2=ff2**(params%sigmasq_dg_po  * log(ff2))
             
             ff = ((fri*frj)/(fr0*fr0))
             ff = ff ** (params%sigmasq_dg_po /2* log(ff2))
             decor=ff/(ff1*ff2)
             
             if (params%decorrel_slope /= 0 ) then
                sloped_decor(:) = (l_divide_3000(:)-1_mcp)*params%decorrel_slope + decor
                
                !limiting range
                where(sloped_decor<0_mcp) &
                     sloped_decor = 0_mcp
                where(sloped_decor>1_mcp) &
                     sloped_decor = 1_mcp
                
                !don't need this since will be set elsewhere
                hasslope=.true.
                decor=1_mcp
             endif
          endif
          frqdep = frqdep*decor
          
          if (hasslope) then
             cl_dusty_clustered(:) = cl_dusty_clustered(:)+ clust2_dg_templ(:)/(d3000*clust2_dg_templ(3000))*params%czero_dg_cl2*frqdep * sloped_decor(:)
          else
             cl_dusty_clustered(:) = cl_dusty_clustered(:)+ clust2_dg_templ(:)/(d3000*clust2_dg_templ(3000))*params%czero_dg_cl2*frqdep 
          endif
       endif
       
    endif
  end function cl_dusty_clustered


  function cl_tsz(params,fri,frj,fr0)
    real(mcp),dimension(2:lmax) :: cl_tsz
    type(foreground_params) :: params
    real(mcp) :: fri,frj,fr0

    cl_tsz(:) = (params%czero_tsz /(d3000*tsz_templ(3000))  * tszFreqDep(fri,fr0) * tszFreqDep(frj,fr0) ) * tsz_templ(:)

  end function cl_tsz


!based on Laurie's email and templates!
!apply cosmological scaling
  function cosmo_scale_ksz(H0,sigma8,omegab,omegam,ns,tau)
    real(mcp) :: H0,sigma8,omegab,omegam,ns,tau
    real(mcp) :: cosmo_scale_ksz
    if (cosmological_scaling_ksz) then
       cosmo_scale_ksz = ((H0/71.0)**1.7 ) &
            * ( (sigma8/.8)**4.7 ) &
            * ( (omegab/.044)**2.1 ) &
            * ( (omegam/.264)**(-0.44) ) &
            * ( (ns/.96)**(-0.19) ) !&
    else
       cosmo_scale_ksz = 1.0
    endif
  end function cosmo_scale_ksz

!based on Laurie's email and templates!
!apply cosmological scaling
  function cosmo_scale_tsz(H0,sigma8,omegab)
    real(mcp) :: H0,sigma8,omegab
    real(mcp) :: cosmo_scale_tsz
    if (cosmological_scaling_tsz) then
       cosmo_scale_tsz = ((H0/71.0)**1.73 ) &
            * ( (sigma8/.8)**8.34 ) &
            * ( (omegab/.044)**2.81 ) 
    else
       cosmo_scale_tsz = 1.0
    endif
  end function cosmo_scale_tsz


  function cl_ksz(params)
    real(mcp),dimension(2:lmax) :: cl_ksz,slope
    type(foreground_params) :: params

    cl_ksz(:) = (params%czero_ksz  / (d3000*ksz_templ(3000))) *ksz_templ(:)
    if (params%czero_ksz2 /= 0) &
         cl_ksz(:) = cl_ksz(:) + (params%czero_ksz2  / (d3000*ksz2_templ(3000))) &
         * ksz2_templ(:)
    if (params%ksz_slope /= 0) then
       slope(:) = l_divide_3000(:)*params%ksz_slope
       slope(:) = slope(:) + (1.0 - slope(3000))
       where(slope<0_mcp) &
            slope = 0_mcp
       cl_ksz(:) = cl_ksz(:) * slope(:)
    endif
  end function cl_ksz

  function flat_tsz_cor()
    logical flat_tsz_cor
    flat_tsz_cor = (.not. ShangModelCorrelationShape)
  end function flat_tsz_cor
  
  ! The correlation fraction between tSZ and dusty sources, normalized to 1 at high ell
  function tsz_dgcl_cor(params)
    real(mcp),dimension(2:lmax) :: tsz_dgcl_cor
    type(foreground_params) :: params
    
    !motivated by shaw analysis of Sehgal sims
!    tsz_dgcl_cor = max(0., (.3 - .2 * exp(-(l-500.)/1000.))/.3)
    !motivated by simplicity
    
    
    if (ShangModelCorrelationShape) then
       tsz_dgcl_cor(:) = -0.0703 * (l_divide_3000(:)*l_divide_3000(:)) + &
            0.612 * l_divide_3000(:) + &
            0.458

    else
       tsz_dgcl_cor(:) = 1.0
       if (params%tsz_cib_slope /= 0) then
          tsz_dgcl_cor(:) = tsz_dgcl_cor(:) + (l_divide_3000(:)-1_mcp)*params%tsz_cib_slope
       endif
    endif
    
  end function tsz_dgcl_cor

  ! The correlation fraction between tSZ and dusty sources, normalized to 1 at high ell
  function tsz_rg_cor()

    real(mcp) :: tsz_rg_cor

    tsz_rg_cor = 1.0

  end function tsz_rg_cor

  !
  ! nu,nu0 in GHz
  !
  ! dBdT is proportional to derivative of planck function
  ! but is normalized so its equal to 1 at nu0
  !
  function dBdT(nu,nu0)

    real(mcp) x, x0, dBdT, dBdT0, nu, nu0

    x0 = nu0/56.78
    dBdT0 = x0**4 * exp(x0) / (exp(x0)-1)**2

    x = nu/56.78
    dBdT = x**4 * exp(x) / (exp(x)-1)**2 / dbdT0

  end function dBdT
  
  !proportional to the Planck function normalized to 1 at nu0
  function Bnu(nu,nu0,T)
    real(mcp) Bnu, nu,nu0,T
    !h/k
    !4.799237 × 10-11 s K
    !expect GHz
    ! so 4.799237e-2 K/GHz
    real(mcp), parameter :: hk = 4.799237e-2
    
    Bnu = (nu/nu0)**3
    Bnu = Bnu * (exp( hk*nu0/T)-1d0) / (exp( hk*nu/T)-1d0) 
    
  end function Bnu

  

  
  !
  ! nu, nu0 in GHz
  ! Gives the tsz frequency dependence normalized so its 1 at nu0
  !
  function tszFreqDep(nu,nu0)

    real(mcp) :: tszFreqDep, tszFreqDep0, nu, nu0, x, x0

    x = nu / 56.78
    x0 = nu0 / 56.78

    tszFreqDep0 = x0*(exp(x0)+1)/(exp(x0)-1) - 4
    tszFreqDep = x*(exp(x)+1)/(exp(x)-1) - 4
    tszFreqDep = tszFreqDep/tszFreqDep0

  end function tszFreqDep
  
!read template in Dl form from file
!since templates are in Cl's internally, convert to cl 
  function read_dl_template(filename)
    character*(*), intent(in) :: filename
    real(mcp), dimension(2:lmax) :: read_dl_template
    real(mcp) :: realtmp
    integer :: ll 
    read_dl_template(:)=0.0
    if (filename/='') then
       call OpenTxtFile(filename, tmp_file_unit)
       do
          read(tmp_file_unit,*,end=2) ll, realtmp
          if (ll>=2 .and. ll<=lmax) &
               read_dl_template(ll) = realtmp*(twopi/((ll+1.0)*ll))
       end do
2      Close(tmp_file_unit)
    end if
  end function read_dl_template

!read template in Cl form from file
  function read_cl_template(filename)
    character*(*), intent(in) :: filename
    real(mcp), dimension(2:lmax) :: read_cl_template
    real(mcp) :: realtmp
    integer :: ll
    read_cl_template(:)=0.0
    if (filename/='') then
       call OpenTxtFile(filename, tmp_file_unit)
       do
          read(tmp_file_unit,*,end=2) ll, realtmp
          if (ll>=2 .and. ll<=lmax) &
               read_cl_template(ll) = realtmp
       end do
2      Close(tmp_file_unit)
    end if
  end function read_cl_template


  !
  ! Loads the clustered, ksz, and tsz templates which are expected to be in
  ! foreground_folder with names cluster_*.dat, ksz.dat, and tsz.dat
  !
  subroutine InitForegroundData(fClustered,fClustered2,fKSZ,fKSZ2,fTSZ,&
       del_alpha,relative_alpha_cluster, relative_alpha_spire)
    integer :: l,dum,i
    character(len=120) :: file
    character*(*) :: fClustered,fClustered2, fKSZ, fKSZ2, fTSZ
    real(mcp) :: del_alpha
    logical :: relative_alpha_cluster,relative_alpha_spire
    SuccessfulInitialization=.true.

    
    ! clustered template
    clust_dg_templ = read_dl_template(fClustered)
    
    clust2_dg_templ = read_dl_template(fClustered2)


    ! KSZ template
    ksz_templ = read_dl_template(fKSZ)
    ! 2nd KSZ template
    ksz2_templ = read_dl_template(fKSZ2)

    ! TSZ template
    tsz_templ = read_dl_template(fTSZ)

    if (MPIRank == 0) then
       print*,'1halo (or all clustered) DG template:',trim(fClustered)
       print*,'2halo (2nd) DG template:',trim(fClustered2)
       print*,'kSZ template:',trim(fKSZ)
       print*,'ksz template2:',trim(fKSZ2)
       print*,'tSZ template:',trim(fTSZ)
    endif
    !prior on separation between clsutered and Poisson terms
    ! ignored if <= 0
    DelAlphaPrior = del_alpha
    AddAlphaPoisson = 0
    if (relative_alpha_cluster) then
       AddAlphaPoisson = 1
    end if

    AddAlphaSpire = 0
    if (relative_alpha_spire) then
       AddAlphaSpire = 1
    end if

    do i=2,lmax
       l_divide_3000(i) = real(i,mcp)/3000_mcp
    enddo

    cirrus_templ(200:lmax) = l_divide_3000(200:lmax)**(-3.2)
    cirrus_templ(2:199)=0.0
    
    clust_rg_templ = l_divide_3000**(-1.2)
    clust_rg_templ(2:49)=0.0

  end subroutine InitForegroundData
  
  subroutine InitRadioAmpPrior(iradio_amp,iradio_unc)
    real(mcp)::iradio_amp, iradio_unc
    radio_amp=iradio_amp
    radio_unc=iradio_unc
  end subroutine InitRadioAmpPrior
  
  function GetRadioAmpPrior()
    real(mcp) :: GetRadioAmpPrior
    GetRadioAmpPrior=radio_amp
  end function GetRadioAmpPrior
  function GetRadioAmpUncPrior()
    real(mcp) :: GetRadioAmpUncPrior
    GetRadioAmpUncPrior=radio_unc
  end function GetRadioAmpUncPrior

  subroutine InitRadioClAmpPrior(iradio_amp,iradio_unc)
    real(mcp)::iradio_amp, iradio_unc
    radio_cl_amp=iradio_amp
    radio_cl_unc=iradio_unc
  end subroutine InitRadioClAmpPrior
  
  function GetRadioClAmpPrior()
    real(mcp) :: GetRadioClAmpPrior
    GetRadioClAmpPrior=radio_cl_amp
  end function GetRadioClAmpPrior
  function GetRadioClAmpUncPrior()
    real(mcp) :: GetRadioClAmpUncPrior
    GetRadioClAmpUncPrior=radio_cl_unc
  end function GetRadioClAmpUncPrior
  
  !
  ! returns alpha prior on Poisson-cluster index
  !
  function GetAlphaPrior()
    real(mcp) :: GetAlphaPrior
    GetAlphaPrior = DelAlphaPrior
  end function GetAlphaPrior

  !
  ! returns cirrus prior pre-factor
  ! <0 is no prior.
  function GetCirrusFactor()
    real(mcp) :: GetCirrusFactor
    GetCirrusFactor = CirrusFactorPrior
  end function GetCirrusFactor


  !
  ! Index of the (i,j) entries of the upper triangular part of an n-by-n matrix
  !
  recursive function utindex(i,j,n)

    integer :: i, j, n, utindex

    if (i <= j) then
       utindex = (i-1)*n+j - i*(i-1)/2
    else
       utindex = utindex(j,i,n)
    end if

  end function utindex



  ! Read parameters from an array of reals ordered in the same order as the type
  ! Basically just cast the array as the type and the parameters line up
  function GetForegroundParamsFromArray(array)

    real(mcp), dimension(:) :: array
    type(foreground_params) :: GetForegroundParamsFromArray

    GetForegroundParamsFromArray = transfer(array,GetForegroundParamsFromArray)

  end function GetForegroundParamsFromArray

  !calculate external foreground prior LnL
  function  getForegroundPriorLnL(foregrounds)
    real(mcp) ::  getForegroundPriorLnL
    type(foreground_params) :: foregrounds
    real(mcp) :: del_alpha
    integer i,j
    real(mcp) :: radio_amp,radio_unc,fradio,radio_cl_amp,radio_cl_unc
    real(mcp) :: cirrus90, cirrus150,cirrus220
    real(mcp) :: prior90, prior150, prior220
    real(mcp) :: cirrus_factor
    real(mcp) :: freq
    integer anyhigh

    getForegroundPriorLnL=0
    
    del_alpha = GetAlphaPrior()
    ! dusty index prior
    if ( del_alpha > 0) then
       getForegroundPriorLnL = getForegroundPriorLnL + (foregrounds%T_dg_po-foregrounds%T_dg_cl)**2/(2*(del_alpha)**2)
!       print*,'dusty index prior:',getForegroundPriorLnL
    end if
    
    if (use_decorrelation_matrix_form) then 
       anyhigh = 0
       if (any(foregrounds%decorrel_matrix < 0) .or. any(foregrounds%decorrel_matrix > 1)) &
            anyhigh=1
       do i=1,MaxNFreq-2 
          do j=i+1,MaxNFreq-1
             if (foregrounds%decorrel_matrix(dmatrix_index(i,j)) < foregrounds%decorrel_matrix(dmatrix_index(i,j+1))) &
                  anyhigh=1
          enddo
       enddo
       if (anyhigh /= 0) &
            getForegroundPriorLnL = getForegroundPriorLnL +  1e6
 !      print*,'post decor prior:',getForegroundPriorLnL
    endif
    
    !cirrus prior
    !updated 2011/05/28 to RKs 08+09 values
!9/22 updated with adhoc 23h value
!also decided to leave out 90 GHz prior since effectively 0 and not well matched to SED shape.
    !9/23 updated to make all three optional. 150/220 will still be default if keyword's aren't in ini file
    cirrus_factor = 1./3.
    if (foregrounds%czero_cirrus .ne. 0) then
       if (ApplyCirrusPrior90) then 
          freq=97.9
          cirrus90=cirrus_power3000(foregrounds,freq,freq)
          
          prior90=0.16*cirrus_factor
          getForegroundPriorLnL = getForegroundPriorLnL + &
               (log(cirrus90 / prior90))**2 / (2*(0.43)**2)
          !     print*,'post cirrus90 prior:',getForegroundPriorLnL       
       endif
       if (ApplyCirrusPrior150) then 
          freq=153.44
          cirrus150=cirrus_power3000(foregrounds,freq,freq)
          
          prior150=0.21*cirrus_factor
          getForegroundPriorLnL = getForegroundPriorLnL + &
               (log(cirrus150/prior150))**2 / (2*(0.564)**2)
          !print*,'post cirrus150 prior:',getForegroundPriorLnL,cirrus150,prior150
       endif
       if (ApplyCirrusPrior220) then 
          freq=219.67
          cirrus220=cirrus_power3000(foregrounds,freq,freq)
          
          prior220=2.19*cirrus_factor
          getForegroundPriorLnL = getForegroundPriorLnL + &
               (log(cirrus220/prior220))**2 / (2*(0.33)**2)
         ! print*,'post cirrus220 prior:',getForegroundPriorLnL,cirrus220,prior220
       endif
    endif
       
    radio_amp = GetRadioAmpPrior()
    radio_unc = GetRadioAmpUncPrior()
    if (radio_amp >= 0 .and. radio_unc > 0) then
       fradio =(foregrounds%czero_rg_po - radio_amp) / radio_unc
       getForegroundPriorLnL = getForegroundPriorLnL + fradio**2 / 2
       if (foregrounds%czero_rg_po < 0) then
          getForegroundPriorLnL = getForegroundPriorLnL + 1e6
       endif
       !print*,'post rg poisson prior:',getForegroundPriorLnL       
    endif
    
    radio_cl_amp = GetRadioClAmpPrior()
    radio_cl_unc = GetRadioClAmpUncPrior()
    if (radio_cl_amp >= 0 .and. radio_cl_unc > 0) then
       fradio =(foregrounds%czero_rg_cl - radio_cl_amp) / radio_cl_unc
       getForegroundPriorLnL = getForegroundPriorLnL + fradio**2 / 2
       if (foregrounds%czero_rg_cl < 0) then
          getForegroundPriorLnL = getForegroundPriorLnL + 1e6
       endif
       !print*,'post rg clus prior:',getForegroundPriorLnL       
    endif
    
  end function getForegroundPriorLnL
  
  
  subroutine InitLensingTemplate(LensingTemplate)
    character*(*),intent(in) :: LensingTemplate
    
    cls_diff = read_dl_template(LensingTemplate)
    
  end subroutine InitLensingTemplate


  subroutine InitFGModel()
    use IniFile
    character(LEN=Ini_max_string_len) SPTtSZTemplate, SPTkSZTemplate, &
         SPTkSZ2Template, SPTClus2Template, &
         SPTClusTemplate, LensingTemplate
    logical relative_alpha_cluster, relative_alpha_spire
    integer i,j,k

!         ShangModelCorrelationShape
!    logical ApplyCirrusPrior90,ApplyCirrusPrior150,ApplyCirrusPrior220
!    logical cosmological_scaling_tsz, cosmological_scaling_ksz
    real(mcp)   radio_amp, radio_unc, del_alpha

    if (SuccessfulInitialization) then 
       return
       !skip - already called successfully
    endif

    k=0
    do j=2,MaxNFreq 
       do i=1,j-1
          k=k+1
          dmatrix_index(i,j)=k
          dmatrix_index(j,i)=k
       enddo
    enddo
    
    LensingTemplate = Ini_Read_String('add_lensing_template', .false.)
    call InitLensingTemplate(LensingTemplate)
        
    SPTtSZTemplate = Ini_Read_String('spt_dataset_tSZ', .false.)
    SPTkSZTemplate = Ini_Read_String('spt_dataset_kSZ', .false.)
    SPTkSZ2Template = Ini_Read_String('spt_dataset_kSZ2', .false.)
    SPTClusTemplate = Ini_Read_String('spt_dataset_clustered', .false.)
    SPTClus2Template = Ini_Read_String('spt_dataset_clustered2', .false.)
        
    del_alpha = Ini_read_Real('spt_prior_clusterpoisson',-1.0)
    radio_amp = Ini_read_Real('radio_ampl_mean',-1.0)
    radio_unc = Ini_read_Real('radio_ampl_unc',-1.0)
    call InitRadioAmpPrior(radio_amp,radio_unc)
    if (MPIRank == 0) &
         print*,'priors, dAlpha,radio amp, radio sigma',&
         del_alpha,radio_amp,radio_unc
    
    applyCIBCalToCirrus = Ini_Read_Logical('apply_cib_cal_to_cirrus',.true.)
    calFactorsToCIB     = Ini_Read_Logical('cal_factors_for_cib',.false.)
    
    CIB_decorrelation_matrix = Ini_Read_Logical('use_cib_decorrelation_matrix',.false.)

    tSZ_CIB_logFreq = Ini_Read_Logical('tsz_cib_logfreq',.true.)
    only1HaloTszCib = Ini_Read_Logical('only_1halo_tsz_cib',.false.)
    combine_spire_for_tszcib_correlation = Ini_Read_Logical(&
         'combine_spire_for_tszcib_correlation',.false.)
    ShangModelCorrelationShape = Ini_Read_Logical(&
         'shang_model_correlation_shape',.false.)
    use_decorrelation_matrix_form = Ini_Read_Logical(&
         'use_decorrelation_matrix_form',.false.)
    use_sigma = Ini_Read_Logical('use_sigma',.false.)
    if (MPIRank == 0) &
         print*,'tSZ-CIB assumptions:',only1HaloTszCib, &
         combine_spire_for_tszcib_correlation,ShangModelCorrelationShape
    
    single_clustered_freq_scaling = Ini_Read_Logical(&
         'single_clustered_freq_scaling',.true.)
    relative_alpha_cluster = Ini_Read_Logical('relative_alpha_cluster',.true.)
    relative_alpha_spire = Ini_Read_Logical('relative_alpha_spire',.true.)
    if (MPIRank == 0) &
         print*,'CIB freq dep. assumptions:',&
         single_clustered_freq_scaling ,relative_alpha_cluster,&
         relative_alpha_spire
    
    cosmological_scaling_ksz = Ini_Read_Logical(&
         'cosmological_scaling_ksz',.false.)
    cosmological_scaling_tsz = Ini_Read_Logical(&
         'cosmological_scaling_tsz',.false.)
    if (MPIRank == 0) &
         print*,'Cosmo scaling ksz, tsz:',&
         cosmological_scaling_ksz,cosmological_scaling_tsz
    
    ApplyCirrusPrior90 = Ini_Read_Logical('apply_prior_cirrus_90ghz',.false.)
    ApplyCirrusPrior150 = Ini_Read_Logical('apply_prior_cirrus_150ghz',.true.)
    ApplyCirrusPrior220 = Ini_Read_Logical('apply_prior_cirrus_220ghz',.true.)
    if (MPIRank == 0) &
         print*,'Cirrus priors used:',&
         ApplyCirrusPrior90,ApplyCirrusPrior150,ApplyCirrusPrior220    
    
    if (SPTtSZTemplate/='' .AND. SPTkSZTemplate/='' .AND. SPTkSZ2Template/='' &
         .AND. SPTClusTemplate/='') then
       call InitForegroundData(SPTClusTemplate, SPTClus2Template, &
            SPTkSZTemplate,SPTkSZ2Template,SPTtSZTemplate,del_alpha, &
            relative_alpha_cluster, relative_alpha_spire)
    else
       print*,'called initFGModel w/o all args'
       print*,SPTtSZTemplate,SPTkSZTemplate,SPTkSZ2Template,SPTClusTemplate
       call   setForegroundsUninitialized()
    end if
  end subroutine InitFGModel

  subroutine OpenReadBinaryFile(aname,aunit,record_length)
    character(LEN=*), intent(IN) :: aname
    integer, intent(in) :: aunit
    integer*8,intent(in) :: record_length
    open(unit=aunit,file=aname,form='unformatted',access='direct',recl=record_length,  err=500)
    return
    
500 call MpiStop('File not found: '//trim(aname))
  end subroutine openReadBinaryFile

  subroutine OpenReadBinaryStreamFile(aname,aunit)
    character(LEN=*), intent(IN) :: aname
    integer, intent(in) :: aunit
    open(unit=aunit,file=aname,form='unformatted',access='stream', err=500)
    return
    
500 call MpiStop('File not found: '//trim(aname))
  end subroutine OpenReadBinaryStreamFile

  subroutine OpenWriteBinaryFile(aname,aunit,record_length)
    character(LEN=*), intent(IN) :: aname
    integer, intent(in) :: aunit
    integer*8,intent(in) :: record_length
    open(unit=aunit,file=aname,form='unformatted',status='replace',access='direct',recl=record_length, err=500)
    return
    
500 call MpiStop('File not able to be written to: '//trim(aname))
  end subroutine OpenWriteBinaryFile



end module
