module spt

  use settings
  use MatrixUtils
  use cmbtypes
  use foregrounds


  implicit none

  ! The bandpowers, in order 90x90, 90x150, 90x220, 150x150, 150x220, 220x220
  integer :: nfreq, nband, nall,maxnbin
  integer,dimension(:),allocatable :: nbins,offsets
  ! cov == The cholesky factored bandpower covariance  
  double precision, dimension(:,:),allocatable :: cov,windows
  integer :: spt_windows_lmin,spt_windows_lmax
  double precision, dimension(:), allocatable :: cl_to_dl_conversion,dl_th,spec
  real(mcp),dimension(:,:), allocatable :: spt_eff_fr
  real(mcp),dimension(:), allocatable :: spt_prefactor
  real(mcp),dimension(5) ::  spt_norm_fr
  integer, dimension(:,:),allocatable :: indices

  logical :: SuccessfulSPTInitialization
  logical :: planckish 
  logical :: CallFGPrior
  logical :: printDlSPT 
  logical :: printDlSPTComponents 

contains
  subroutine setSPTUninitialized()
    SuccessfulSPTInitialization=.false.
  end subroutine setSPTUninitialized


  ! Load the bandpowers, covariance, and window functions
  subroutine InitSptData(iPlanckish)
    use IniFile
    logical,optional, intent(in) :: iPlanckish
    character(LEN=Ini_max_string_len) :: desc_file,bp_file, cov_file, window_folder
    logical binaryCov,binaryWindows
    integer i,j,k,dum,lwin,nwin
    real(mcp) rtmp
    double precision, allocatable, dimension(:) :: locwin
    integer*8 :: delta,offset
    integer*4 :: efflmin,efflmax,j0,j1
    printDlSPT = Ini_Read_Logical('print_cmb_TT',.false.)
    if (printDlSPT .and. MPIRank /= 0) then
       call MPIStop('Warning - print_cmb_TT is not thread-safe, quitting...')
    endif
    printDlSPTComponents = Ini_Read_Logical('print_foreground_terms',.false.)
    if (printDlSPTComponents .and. MPIRank /= 0) then
       call MPIStop('Warning - print_cmb_TT is not thread-safe, quitting...')
    endif

    desc_file = Ini_Read_String('spt_dataset_description',.false.)
    bp_file = Ini_Read_String('spt_dataset_bandpowers',.false.)

    !default to binary since faster i/o
    binaryCov=.True.
    binaryWindows=.True.
    cov_file = Ini_Read_String('spt_dataset_binary_covariance',.false.)
    if (cov_file == '') then 
       cov_file = Ini_Read_String('spt_dataset_covariance',.false.)
       binaryCov=.False.
    endif
    window_folder = Ini_Read_String('spt_dataset_binary_windir',.false.)
    if (window_folder == '') then 
       window_folder = Ini_Read_String('spt_dataset_windir',.false.)
       binaryWindows=.False.
    endif

    if (bp_file=='' .or. desc_file=='' .or. window_folder=='' .or. cov_file=='') then
       print*,'Missing required spt key: received: ',desc_file,bp_file,cov_file,window_folder
       stop
    endif
    !get info
    call OpenTxtFile(desc_file, tmp_file_unit)
    read (tmp_file_unit,*) nall,nfreq
    if (nfreq > MaxNFreq) &
         call MpiStop('spt initialized with more than allowed Nfrequencies')

    nband = (nfreq)*(nfreq+1)/2
    allocate(nbins(nband))
    do i=1,nband
       read (tmp_file_unit,*) j
       nbins(i)=j
    end do
    maxnbin=maxval(nbins(:))
    do i=1,5
       read (tmp_file_unit,*) rtmp
       spt_norm_fr(i) = rtmp
    enddo
    planckish=.false.
    if (present(iPlanckish)) then
       planckish=iPlanckish
    endif
    if (planckish) then 
       spt_norm_fr(5) = 143.
       print*,'using 143 as tSZ center freq to match planck'
    endif
    read (tmp_file_unit,*) spt_windows_lmin, spt_windows_lmax



    if (lmax .lt. spt_windows_lmax) spt_windows_lmax=lmax
    if (spt_windows_lmin .lt. 2) then
       print*, "Unexpected lmin for windows!"
       stop
    endif
    nwin = spt_windows_lmax-spt_windows_lmin+1
    allocate(cl_to_dl_conversion(spt_windows_lmin:spt_windows_lmax),dl_th(spt_windows_lmin:spt_windows_lmax))
    do j=spt_windows_lmin,spt_windows_lmax
       cl_to_dl_conversion(j) = (1.0 + j)*(j / TWOPI)
    enddo
    allocate(windows(spt_windows_lmin:spt_windows_lmax,1:nall),spt_eff_fr(5,nfreq), &
         spt_prefactor(nfreq),indices(2,nband),offsets(nband),&
         spec(1:nall),cov(1:nall,1:nall))
    do j=1,nfreq
       do i=1,5
          read (tmp_file_unit,*) rtmp
          spt_eff_fr(i,j) = rtmp
       enddo
    enddo
    do j=1,nfreq
       read (tmp_file_unit,*) spt_prefactor(j)
    enddo
    close(tmp_file_unit)

    ! Read in bandpowers
    call OpenTxtFile(bp_file, tmp_file_unit)
    do i=1,nall
       read (tmp_file_unit,*) dum,spec(i)
    end do
    close(tmp_file_unit)
    
    
    ! Read in covariance
    if (binaryCov) then 
       call OpenReadBinaryFile(cov_file,tmp_file_unit,nall*8_8)
       do i=1,nall
          read(tmp_file_unit,rec=i)cov(:,i)
       enddo
    else
       call OpenTxtFile(cov_file, tmp_file_unit)
       do i=1,nall
          do j=1,nall
             read (tmp_file_unit,*) cov(j,i)
          end do
       end do
    endif
    close (tmp_file_unit)
       


    call Matrix_CholeskyDouble(cov)
    
    ! Read in windows
    if (binaryWindows) then
       
       call OpenReadBinaryStreamFile(trim(window_folder)//trim('window_binary'),tmp_file_unit)
       read(tmp_file_unit,pos=1)efflmin,efflmax
       allocate(locwin(efflmin:efflmax))
       if ((efflmax .lt. spt_windows_lmin) .or. (efflmin .gt. spt_windows_lmax)) &
            call MpiStop('unallowed l-ranges for binary window functions')
       j0=efflmin
       if (spt_windows_lmin > j0) j0=spt_windows_lmin
       j1=efflmax
       if (spt_windows_lmax < j1) j1=spt_windows_lmax
       if (j1 < j0) &
            call MpiStop('unallowed l-ranges for binary window functions - no allowed ells')
       delta=(efflmax-efflmin+1)*8_8
       offset=2 * 4_8+1
       do i=1,nall
          read(tmp_file_unit,pos=((i-1)*delta + offset)) locwin
          windows(j0:j1,i)=locwin(j0:j1)
       end do
       close(tmp_file_unit)
       deallocate(locwin)
    else
       do i=1,nall
          call OpenTxtFile(trim(window_folder)//trim(numcat('window_',i)),tmp_file_unit)
          do j=spt_windows_lmin,spt_windows_lmax
             read (tmp_file_unit,*) dum, windows(j,i)
          end do
          close(tmp_file_unit)
       end do
    endif
    CallFGPrior=.false.
    if  (HaveForegroundsBeenInitialized() == .false.) then
       !this is first dataset to use foregrounds
       !so include Prior
       CallFGPrior = .true.
    endif
    i=0
    do j=1,nfreq
       do k=j,nfreq 
          i=i+1
          indices(1,i)=j
          indices(2,i)=k
       end do
    end do
    offsets(1)=1
    do i=2,nband
       offsets(i)=offsets(i-1)+nbins(i-1)
    enddo
       
    SuccessfulSPTInitialization = .true.
  end subroutine InitSptData
  

  function SptLnLike(cmbcls, freq_params,tszfac,kszfac)
    use omp_lib
    real(mcp), dimension(:,:), intent(in) :: cmbcls
    real(mcp), dimension(:) :: freq_params
    real(mcp) :: SptLnLike
    real(mcp), intent(in) :: tszfac, kszfac
    
    real(mcp), dimension(5) :: dum
    double precision, dimension(1:nall) :: deltacb
    double precision,dimension(1:maxnbin)::tmpcb
    real(mcp), dimension(nfreq) :: a_calib,cib_cals
    real(mcp) :: cmbcl
    integer l,b,i,j,k,ii,jj,kk,thisoffset,thisnbin,fid
    type(foreground_params) :: foregrounds
    real(mcp) :: norm
    real*4, dimension(2) :: arr
    real*4, dimension(10) :: comp_arr
    
    real(mcp),dimension(2:lmax) :: cl_fgs
    real(mcp), dimension(2:lmax,7) :: component_spectra


    if (HaveForegroundsBeenInitialized() .eq. .false.) then
       write(*,*)'trying to call SPT likelihood w/o initializing foregrounds'
       stop
    endif
    
    foregrounds = GetForegroundParamsFromArray(freq_params(1:nForegroundParams))
    foregrounds.czero_ksz = foregrounds.czero_ksz*kszfac
    foregrounds.czero_tsz = foregrounds.czero_tsz*tszfac

    if (planckish) &
         foregrounds%tsz_dg_cor_const = -1*foregrounds%tsz_dg_cor_const
    if (calFactorsToCIB) then
       a_calib(:)=1
       cib_cals(:) = freq_params(nForegroundParams+1:nForegroundParams+nfreq)
    else
       cib_cals(:)=1
       a_calib(:) = freq_params(nForegroundParams+1:nForegroundParams+nfreq)
    endif
    deltacb=0
    if (printDlSPT) then
       call printForegrounds(foregrounds)
       print*,'a_calib spt:',a_calib
    endif

    !$OMP PARALLEL DO  DEFAULT(NONE), &
    !$OMP  SHARED(indices,foregrounds,spt_eff_fr,spt_norm_fr,cl_to_dl_conversion,nbins,nfreq,cmbcls,spt_windows_lmax,spt_windows_lmin,windows,spt_prefactor,a_calib,deltacb,offsets,spec,nband,printDlSPT,printDlSPTComponents,cib_cals), &
    !$OMP  private(i,j,k,cl_fgs,dl_th,tmpcb,thisnbin,l,thisoffset,fid,arr,component_spectra,comp_arr), &
    !$OMP SCHEDULE(STATIC)
    do i=1,nband
       j=indices(1,i)
       k=indices(2,i)
       thisoffset=offsets(i)
       thisnbin=nbins(i)
       tmpcb(:)=0
       
       !first get theory spectra
       if (printDlSPTComponents) then
          cl_fgs(2:lmax) = cl_foreground(foregrounds,j,k,nfreq,spt_eff_fr, &
               spt_norm_fr,cib_cals,component_spectra)
       else
          cl_fgs(2:lmax) = cl_foreground(foregrounds,j,k,nfreq,spt_eff_fr, &
               spt_norm_fr,cib_cals)
       endif
       !add CMB
       cl_fgs(2:lmax)=cl_fgs(2:lmax)+cmbcls(2:lmax,1)
       
       !go to Dls
       dl_th(:) = cl_fgs(spt_windows_lmin:spt_windows_lmax)*cl_to_dl_conversion(:)
       if (printDlSPT) then
          fid=33+omp_get_thread_num()
          call OpenWriteBinaryFile(trim(numcat('suxp_spt_',i)),fid,4_8 * 2)
          do l=spt_windows_lmin,spt_windows_lmax
             arr(1)=l
             arr(2)=dl_th(l)
             write(fid,rec=l-spt_windows_lmin+1) arr(1:2)
          enddo
          close(fid)
       endif
       if (printDlSPTComponents) then
          fid=33+omp_get_thread_num()
          call OpenWriteBinaryFile(trim(numcat('suxp_spt_components_',i)),fid,4_8 * 10)
          do l=spt_windows_lmin,spt_windows_lmax
             comp_arr(1)=l
             comp_arr(2)=dl_th(l)
             comp_arr(3) = cmbcls(l,1) * cl_to_dl_conversion(l-spt_windows_lmin+1)
             comp_arr(4) = component_spectra(l,1) * cl_to_dl_conversion(l-spt_windows_lmin+1)
             comp_arr(5) = component_spectra(l,2) * cl_to_dl_conversion(l-spt_windows_lmin+1)
             comp_arr(6) = component_spectra(l,3) * cl_to_dl_conversion(l-spt_windows_lmin+1)
             comp_arr(7) = component_spectra(l,4) * cl_to_dl_conversion(l-spt_windows_lmin+1)
             comp_arr(8) = component_spectra(l,5) * cl_to_dl_conversion(l-spt_windows_lmin+1)
             comp_arr(9) = component_spectra(l,6) * cl_to_dl_conversion(l-spt_windows_lmin+1)
             comp_arr(10) = component_spectra(l,7) * cl_to_dl_conversion(l-spt_windows_lmin+1)
             write(fid,rec=l-spt_windows_lmin+1) comp_arr(1:10)
          enddo
          close(fid)
       endif
       

       !now bin with window functions
       call dgemv('T',spt_windows_lmax-spt_windows_lmin+1,thisnbin,1.0d0,&
            windows(:,thisoffset:thisoffset+thisnbin-1),spt_windows_lmax-spt_windows_lmin+1,&
            dl_th,1,0.0d0,tmpcb,1)
       
       !apply prefactors
       tmpcb(:) = tmpcb(:) * spt_prefactor(k)*spt_prefactor(j)*a_calib(j)*a_calib(k)
       if (printDlSPT) then
          open(fid,file=trim(numcat('est_bandpowers_spt_',i)))
          write(fid,*)'# ',j,k
          do l=1,thisnbin
             write(fid,*)l,tmpcb(l),spec(thisoffset+l-1)
          enddo
          close(fid)
       endif
       
       
       !get delta Dl's for covariance calculation
       deltacb(thisoffset:thisoffset+thisnbin-1) = tmpcb(1:thisnbin) - spec(thisoffset:thisoffset+thisnbin-1)
    enddo

    SptLnLike = Matrix_GaussianHalfChisqCholDouble(cov,deltacb) 
      if (CallFGPrior) then
       SptLnLike = SptLnLike + getForegroundPriorLnL(foregrounds)
    endif


    if (printDlSPT) then
       print *, 'Spt Likes',SptLnLike
       print*,'SPT chisq:',Matrix_GaussianHalfChisqCholDouble(cov,deltacb)*2
       print *,'Spt Chisq w prior',(SptLnLike)*2
    endif

  end function SptLnLike
  

end module spt
