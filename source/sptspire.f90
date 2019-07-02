module sptspire

  use settings
  use MatrixUtils
  use cmbtypes
  use foregrounds


  implicit none

  ! The bandpowers, in order 90x90, 90x150, 90x220, 150x150, 150x220, 220x220
  integer :: nfreq, nband,  nall,maxnbin
  integer,dimension(:),allocatable :: nbins,offsets

  ! cov == The cholesky factored bandpower covariance  
  double precision, dimension(:,:),allocatable :: cov,windows
  integer :: sptspire_windows_lmin,sptspire_windows_lmax
  double precision, dimension(:), allocatable :: cl_to_dl_conversion,dl_th,spec
  real(mcp),dimension(:,:), allocatable :: sptspire_eff_fr
  real(mcp),dimension(:), allocatable :: sptspire_prefactor,sptspire_ukcmb_to_jysr
  integer,  dimension(:), allocatable :: sptspire_cldl,sptspire_useband
  integer, dimension(:,:),allocatable :: indices
  real(mcp),dimension(5) ::  sptspire_norm_fr

  logical :: SuccessfulSPTSPIREInitialization
  logical :: CallFGPrior
  logical :: printDlCMB
  logical :: planckish
  integer :: n_swapped_cals
  logical :: ApplySpireCalPrior
  logical :: DropSptSpire, DropSptOnly, DropSpireOnly
contains
  subroutine setSPTSPIREUninitialized()
    SuccessfulSPTSPIREInitialization=.false.
  end subroutine setSPTSPIREUninitialized


  ! Load the bandpowers, covariance, and window functions
  subroutine InitSptspireData(iPlanckish)
    use IniFile
    logical,optional, intent(in) :: iPlanckish
    character(LEN=Ini_max_string_len) :: desc_file,bp_file, cov_file, window_folder
    integer i,j,k,dum,lwin,l,nwin
    real(mcp) rtmp,ll,w
    logical binaryCov,binaryWindows
    double precision, dimension(2) :: arr

    DropSptSpire = Ini_Read_Logical('skip_sptxspire',.false.)
    DropSptOnly = Ini_Read_Logical('skip_sptonly',.false.)
    DropSpireOnly = Ini_Read_Logical('skip_spireonly',.false.)
    if (DropSptSpire) then
       if (DropSptOnly .or. DropSpireOnly) then
          print*,'sptspire.f90: too many drops requested.'
          call mpistop()
       endif
    else
       if (DropSptOnly .and. DropSpireOnly) then
          print*,'sptspire.f90: too many drops requested.'
          call mpistop()
       endif
    endif
    
    n_swapped_cals = Ini_Read_Int('pretend_first_N_spt_cals_for_spire',0)
    ApplySpireCalPrior = n_swapped_cals > 0
    if (n_swapped_cals < 0) n_swapped_cals = -1*n_swapped_cals
    if (n_swapped_cals /= 0 .and. n_swapped_cals /= 1 .and. n_swapped_cals /= 3) then
       print*,'Bad value for pretend_first_N_spt_cals_for_spire',n_swapped_cals
       call mpistop()
    endif
    if (n_swapped_cals > 0) then 
       if (nfreq /= 6) then 
          print*,'can only handle swapped cal factors when have 6 freqs'
          call mpistop()
       endif
       print*,'Using swapped spire cals: ',n_swapped_cals
       print*,'Applying cal prior of 7%:',ApplySpireCalPrior
    endif

    desc_file = Ini_Read_String('sptspire_dataset_description',.false.)
    bp_file = Ini_Read_String('sptspire_dataset_bandpowers',.false.)

    !default to binary since faster i/o
    binaryCov=.True.
    binaryWindows=.True.
    cov_file = Ini_Read_String('sptspire_dataset_binary_covariance',.false.)
    if (cov_file == '') then 
       cov_file = Ini_Read_String('sptspire_dataset_covariance',.false.)
       binaryCov=.False.
    endif
    window_folder = Ini_Read_String('sptspire_dataset_binary_windir',.false.)
    if (window_folder == '') then 
       window_folder = Ini_Read_String('sptspire_dataset_windir',.false.)
       binaryWindows=.False.
    endif
    
    printDlCMB = Ini_Read_Logical('print_cmb_TT',.false.)
    if (printDlCMB .and. MPIRank /= 0) then
       call MPIStop('Warning - print_cmb_TT is not thread-safe, quitting...')
    endif

    if (bp_file=='' .or. desc_file=='' .or. window_folder=='' .or. cov_file=='') then
       print*,'Missing required sptspire key: received: ',trim(desc_file),trim(bp_file),trim(cov_file),trim(window_folder)
       stop
    endif
    
    !get info
    call OpenTxtFile(desc_file, tmp_file_unit)
    read (tmp_file_unit,*) nall,nfreq
    if (nfreq > MaxNFreq) &
         call MpiStop('sptspire initialized with more than allowed Nfrequencies')
    nband = (nfreq)*(nfreq+1)/2
    allocate(nbins(nband),indices(2,nband),offsets(nband))
    do i=1,nband
       read (tmp_file_unit,*) j
       nbins(i)=j
    end do
    maxnbin=maxval(nbins(:))
    do i=1,5
       read (tmp_file_unit,*) rtmp
       sptspire_norm_fr(i) = rtmp
    enddo
    planckish=.false.
    if (present(iPlanckish)) then
       planckish=iPlanckish
    endif
    if (planckish) then
       sptspire_norm_fr(5) = 143.
       print*,'using 143 as tSZ center freq to match planck'
    endif

    read (tmp_file_unit,*) sptspire_windows_lmin, sptspire_windows_lmax


    if (lmax .lt. sptspire_windows_lmax) sptspire_windows_lmax=lmax
    if (sptspire_windows_lmin .lt. 2) then
       print*, "Unexpected lmin for windows!"
       stop
    endif
    allocate(cl_to_dl_conversion(sptspire_windows_lmin:sptspire_windows_lmax),dl_th(sptspire_windows_lmin:sptspire_windows_lmax))
    do j=sptspire_windows_lmin,sptspire_windows_lmax
       cl_to_dl_conversion(j) = (1.0 + j)*(j / TWOPI)
    enddo
    allocate(windows(sptspire_windows_lmin:sptspire_windows_lmax,1:nall),sptspire_eff_fr(5,nfreq), &
         sptspire_prefactor(nband),sptspire_cldl(nband),&
         sptspire_ukcmb_to_jysr(nfreq),&
         spec(1:nall),cov(1:nall,1:nall))
    do j=1,nfreq
       do i=1,5
          read (tmp_file_unit,*) rtmp
          sptspire_eff_fr(i,j) = rtmp
       enddo
    enddo
    do j=1,nband
       read (tmp_file_unit,*) &
            sptspire_prefactor(j), sptspire_cldl(j)
    enddo
    do j=1,nfreq
       read (tmp_file_unit,*) sptspire_ukcmb_to_jysr(j)
    enddo
    ! Read in bandpowers
    call OpenTxtFile(bp_file, tmp_file_unit)
    do i=1,nall
       read (tmp_file_unit,*) dum,spec(i)
    end do
    close(tmp_file_unit)
    
    
    ! Read in covariance
    if (binaryCov) then 
       call OpenReadBinaryFile(cov_file,tmp_file_unit,nall*8_8)
       !print*,nall,nall*8_8
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

    !now we have everything needed to implement the dropping of certain bands.
    !assuming 6 freqs as per earlier cut
    if (DropSptOnly) then 
       sptspire_prefactor(1:3) = -1 
       sptspire_prefactor(7:8) = -1 
       sptspire_prefactor(12)  = -1 
    endif
    if (DropSpireOnly) then 
       sptspire_prefactor(16:21) = -1 
    endif
    if (DropSptSpire) then 
       sptspire_prefactor(4:6) = -1 
       sptspire_prefactor(9:11) = -1
       sptspire_prefactor(13:15)  = -1
    endif    
    do i=1,21
       if (sptspire_prefactor(i) < 0) then 
          cov(offsets(i):offsets(i)+nbins(i)-1,:) = 0
          cov(:,offsets(i):offsets(i)+nbins(i)-1) = 0
          do j=offsets(i),offsets(i)+nbins(i)-1
             cov(j,j)=1
          enddo
       endif
    enddo

    ! back to normal business
    call Matrix_CholeskyDouble(cov)
    
    ! Read in windows
    windows(:,:)=0

    if (binaryWindows) then
       do i=1,nall
          call OpenReadBinaryFile(trim(window_folder)//trim(numcat('window_',i)),tmp_file_unit,2*8_8)
          j=1
          do 
             read(tmp_file_unit,rec=j, err=111) arr
             j=j+1
             l=idnint(arr(1))
             if (dabs(l-arr(1)) > 1e-4) then
                stop 'non-integer l in window file'
             endif
             if (l>=sptspire_windows_lmin .and. l<=sptspire_windows_lmax) then
                windows(l,i)=arr(2)
             endif
          enddo
111       close(tmp_file_unit)
       end do
    else
       do i=1,nall
          call OpenTxtFile(trim(window_folder)//trim(numcat('window_',i)),tmp_file_unit)
          do 
             read(tmp_file_unit,*, end=1) ll, w
             l=nint(ll)
             if (abs(l-ll) > 1e-4) stop 'non-integer l in window file'
             if (l>=sptspire_windows_lmin .and. l<=sptspire_windows_lmax) then
                windows(l,i)=w
             endif
          enddo
1         close(tmp_file_unit)
       end do
    endif
    
    CallFGPrior=.false.
    if  (HaveForegroundsBeenInitialized() == .false.) then
       !this is first dataset to use foregrounds                                
       !so include Prior                                                        
       CallFGPrior = .true.
    endif

       
    SuccessfulSPTSPIREInitialization = .true.
  end subroutine InitSptspireData
  

  function SptspireLnLike(cmbcls, freq_params,tszfac,kszfac)
    use omp_lib
    real(mcp), dimension(:,:), intent(in) :: cmbcls
    real(mcp), dimension(:) :: freq_params
    real(mcp) :: SptspireLnLike
    real(mcp), intent(in) :: tszfac, kszfac
    
    real(mcp), dimension(5) :: dum
    double precision, dimension(1:nall) :: deltacb
    double precision,dimension(1:maxnbin)::tmpcb
    real(mcp), dimension(nfreq) :: a_calib,cib_cals
    
    integer l,b,i,j,k,ii,jj,kk,thisoffset,thisnbin,fid
    type(foreground_params) :: foregrounds
    real(mcp) :: norm
    real(mcp),dimension(3) :: spirescalings
    real(mcp),dimension(2:lmax) :: cl_fgs
    real*4, dimension(2) :: arr

    if (HaveForegroundsBeenInitialized() .eq. .false.) then
       write(*,*)'trying to call SPTSPIRE likelihood w/o initializing foregrounds'
       stop
    endif


    
    foregrounds = GetForegroundParamsFromArray(freq_params(1:nForegroundParams))
    foregrounds.czero_ksz = foregrounds.czero_ksz*kszfac
    foregrounds.czero_tsz = foregrounds.czero_tsz*tszfac
    if (planckish) &
         foregrounds%tsz_dg_cor_const = -1*foregrounds%tsz_dg_cor_const
    !these are extra calibration factors that you can allow
    if (calFactorsToCIB) then
       a_calib(:)=1
       cib_cals(:) = freq_params(nForegroundParams+1:nForegroundParams+nfreq)
    else
       cib_cals(:)=1
       a_calib(:) = freq_params(nForegroundParams+1:nForegroundParams+nfreq)
    endif
    if (n_swapped_cals == 1) then 
       cib_cals(1:3)=1
       a_calib(1:3)=1
       a_calib(4:nfreq) = a_calib(4:nfreq) * freq_params(nForegroundParams+1)
    endif
    if (n_swapped_cals == 3) then 
       cib_cals(1:3)=1
       a_calib(1:3)=1
       a_calib(4:nfreq) = a_calib(4:nfreq) * freq_params(nForegroundParams+1:nForegroundParams+3)
    endif
    
    if (printDlCMB) then
       print*,'mcp',mcp
       call OpenWriteBinaryFile('suxp_cl_cmb',33,4_8 * 2)
       do i=2,lmax
          arr(1)=i
          arr(2)=cmbcls(i,1)
          write(33,rec=i-1)arr(1:2)
       enddo
          
!       open(33,file='suxp_cl_cmb')
!       do i=2,lmax
!          write(33,*)i,cmbcls(i,1)
!       enddo
       close(33)
       call    printForegrounds(foregrounds)
       print*,'cals:',a_calib(:)
    endif
    

!    i=0
!    thisoffset=1
!    do j=1,nfreq
 !      do k=j,nfreq
 !         i=i+1

    !initialize to 0. needed because some may not be filled in with drop logicals
    deltacb(:)=0

!lmax is a parameter and doesn't need to be listed
    !$OMP PARALLEL DO  DEFAULT(NONE), &
    !$OMP  SHARED(indices,foregrounds,sptspire_eff_fr,sptspire_norm_fr,cl_to_dl_conversion,nbins,nfreq), &
    !$OMP  SHARED(cmbcls,sptspire_windows_lmax,sptspire_windows_lmin,windows,sptspire_prefactor,a_calib,cib_cals), &
    !$OMP  SHARED( sptspire_cldl,printDlCMB,spec,   deltacb,offsets,sptspire_ukcmb_to_jysr,nband), &
    !$OMP  private(i,j,k,cl_fgs,dl_th,tmpcb,thisnbin,l,thisoffset,fid,arr), &
    !$OMP SCHEDULE(STATIC)
    do i=1,nband
       j=indices(1,i)
       k=indices(2,i)

       !drop logicals implemented by setting the prefactor to be negative
       ! skip the dropped bandpowers
       if (sptspire_prefactor(i) < 0) cycle


       thisoffset=offsets(i)
       
       thisnbin=nbins(i)

       if (thisnbin <= 0) cycle !move to next loop
       tmpcb(:)=0
       
       !first get theory spectra
       cl_fgs(2:lmax) = cl_foreground(foregrounds,j,k,nfreq,sptspire_eff_fr,sptspire_norm_fr,cib_cals)
       
       !add CMB
       cl_fgs(2:lmax)=cl_fgs(2:lmax)+cmbcls(2:lmax,1)
       
       !go to Dls
       if (sptspire_cldl(i) .eq. 0) then 
          dl_th(:) = cl_fgs(sptspire_windows_lmin:sptspire_windows_lmax)*&
               sptspire_ukcmb_to_jysr(k)*sptspire_ukcmb_to_jysr(j)
       else             
          dl_th(:) = cl_fgs(sptspire_windows_lmin:sptspire_windows_lmax)*&
               cl_to_dl_conversion(:)
       endif
       
       if (printDlCMB) then
          fid=33+omp_get_thread_num()

          call OpenWriteBinaryFile(trim(numcat('suxp_sptspire_',i)),fid,4_8 * 2)
          do l=sptspire_windows_lmin,sptspire_windows_lmax
             arr(1)=l
             arr(2)=dl_th(l)
             write(fid,rec=l-sptspire_windows_lmin+1) arr(1:2)
          enddo
!          open(fid,file=trim(numcat('suxp_sptspire_',i)))
!          write(fid,*)'# ',j,k
!          do l=sptspire_windows_lmin,sptspire_windows_lmax
!             write(fid,*)l,dl_th(l)
!          enddo
          close(fid)
       endif
       
       
       !now bin with window functions
       call dgemv('T',sptspire_windows_lmax-sptspire_windows_lmin+1,thisnbin,1.0d0,&
            windows(:,thisoffset:thisoffset+thisnbin-1),&
            sptspire_windows_lmax-sptspire_windows_lmin+1,&
            dl_th,1,0.0d0,tmpcb,1)
       
       !          if (i == 1) then 
       !             print*,'l=2100,2300',dl_th(2098),dl_th(2298)
!             print*,'bin1/2',tmpcb(1),tmpcb(2)
!             open(33,file='suxp')
!             do l=sptspire_windows_lmin,sptspire_windows_lmax
!                write(33,'(2E)')l,dl_th(l)
!             enddo
!             close(33)
!          endif
          
          !apply prefactors
       tmpcb(:) = tmpcb(:) * sptspire_prefactor(i) *a_calib(j)*a_calib(k)
          
       if (printDlCMB) then
          open(fid,file=trim(numcat('est_bandpowers_sptspire_',i)))
          write(fid,*)'# ',j,k
          do l=1,thisnbin
             write(fid,*)l,tmpcb(l)
          enddo
          close(fid)
       endif
       
       !get delta Dl's for covariance calculation
       deltacb(thisoffset:thisoffset+thisnbin-1) = tmpcb(1:thisnbin) - spec(thisoffset:thisoffset+thisnbin-1)
       
!       thisoffset=thisoffset+thisnbin
    enddo

!    print*,'deltacb stats',nall,thisoffset-1
!    print*,'deltas:',deltacb(nall-4:nall)
!    print*,'b',tmpcb(thisnbin-4:thisnbin)
!    print*,'c',spec(nall-5:nall)
!    print*,'deltaspt:',deltacb(1:5)
!    print*,'b',tmpcb(1:5)
!    print*,'c',spec(1:5)
!    print*,'factorspire:',a_calib(6),sptspire_prefactor(21)
!    print*,'factorspt:',a_calib(1),sptspire_prefactor(1)
!    print*,sptspire_windows_lmin,sptspire_windows_lmax

    SptspireLnLike = Matrix_GaussianLogLikeCholDouble(cov,deltacb) 
!    print*,'SPTSPIRELnL pre fg:',SptspireLnLike
    if (CallFGPrior) then 
       SptspireLnLike = SptspireLnLike + getForegroundPriorLnL(foregrounds)
 !      print*,'call FG prior'
    endif
    if (ApplySpireCalPrior) then 
       SptSpireLnLike  = SptSpireLnLike + sum(((freq_params(nForegroundParams+1:nForegroundParams+n_swapped_cals)-1)/0.07)**2)
    endif


!    print*,'SPTSPIRELnL:',SptspireLnLike
    if (printDlCMB) then
       print *, 'Sptspire Likes',SptspireLnLike,SptspireLnLike-getForegroundPriorLnL(foregrounds)
       print *,Matrix_GaussianLogLikeCholDouble(cov,deltacb*0.0)
       print *,'Sptspire Chisq',(SptspireLnLike-Matrix_GaussianLogLikeCholDouble(cov,deltacb*0.0))*2
    endif
  !  stop
  end function SptspireLnLike
  
end module sptspire
