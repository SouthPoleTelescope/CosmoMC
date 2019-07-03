  !Likelihood code used in 
  !SPTpol+SZ l=2000-11000 power spectrum
  !For questions, please contact Christian Reichardt
  module CMB_SPT_hiell_2019
    use Likelihood_Cosmology
    use CosmoTheory
    use CosmologyTypes
    use FileUtils
    use MatrixUtils
    use foregrounds

    implicit none

    ! The bandpowers, in order 90x90, 90x150, 90x220, 150x150, 150x220, 220x220
    integer :: nfreq, nband, nall,maxnbin
    integer,dimension(:),allocatable :: nbins,offsets
    ! cov == The cholesky factored bandpower covariance  
    double precision, dimension(:,:),allocatable :: cov, windows, beam_err, cov_w_beam
    integer :: spt_windows_lmin,spt_windows_lmax
    double precision, dimension(:), allocatable :: cl_to_dl_conversion,dl_th,spec,ells
    double precision, dimension(:,:), allocatable :: spt_eff_fr
    double precision, dimension(:), allocatable :: spt_prefactor
    double precision, dimension(5) ::  spt_norm_fr
    integer, dimension(:,:),allocatable :: indices
    
    logical :: SuccessfulSPTInitialization
    logical :: normalizeSZ_143GHz
    logical :: CallFGPrior
    
    logical :: SuccessfulSPTHiellInitialization
    logical :: printDlSPT, printDlSPTComponents
    logical :: plankish
    logical :: CallFGPrioer
    
    logical :: binaryCov, binaryWindows, binaryBeamErr
    
    !real(mcp) :: meanBeam
    real(mcp) :: meanTcal, sigmaTcal
    real(mcp) :: meanPcal, sigmaPcal
    real(mcp) :: meankappa, sigmakappa
    real(mcp) :: meanAlphaEE, sigmaAlphaEE
    real(mcp) :: meanAlphaTE, sigmaAlphaTE


    Type, extends(TCMBLikelihood) :: TSPTHiEllLike
  contains
    procedure :: ReadIni => SPT_HiEll_ReadIni
    procedure :: InitSPTHiEllData
    procedure :: LogLike => SPTHiEllLnLike
 end Type TSPTHiEllLike

contains

  subroutine setSPTUninitialized()
    SuccessfulSPTHiellInitialization=.false.
  end subroutine setSPTUninitialized


 subroutine SPT_HiEll_ReadIni(this, Ini)
   use IniFile
   use IniObjects
   implicit none
   class(TSPTHiEllLike) :: this
   class(TSettingIni) :: Ini
   character (LEN=Ini_max_string_len) :: desc_file
   character (LEN=Ini_max_string_len) :: bp_file, param_file
   character (LEN=Ini_max_string_len) :: cov_file,beamerr_file
   character (LEN=Ini_max_string_len) :: window_folder

   call InitFGModel(Ini)
   
   param_file = Ini%Read_String_Default( 'spt_hiell_params_file','')
   call this%loadParamNames(param_file)

   desc_file = Ini%Read_String_Default('spt_hiell_description','')
   bp_file = Ini%Read_String_Default('spt_hiell_bandpowers','')
   
   
   !default to binary since faster i/o
   binaryCov=.True.
   binaryWindows=.True.
   binaryBeamErr=.True.
   cov_file =  Ini%Read_String_Default('spt_hiell_binary_covariance','')
   if (cov_file == '') then 
      cov_file = Ini%Read_String_Default('spt_hiell_covariance','')
      binaryCov=.False.
   endif
   window_folder = Ini%Read_String_Default('spt_hiell_binary_windir','')
   if (window_folder == '') then 
      window_folder = Ini%Read_String_Default('spt_hiell_windir','')
      binaryWindows=.False.
   endif
   beamerr_file = Ini%Read_String_Default('spt_hiell_binary_beamerr','')
   if (beamerr_file == '') then 
      beamerr_file = Ini%Read_String_Default('spt_hiell_beamerr','')
      binaryBeamErr=.False.
   endif
      
   normalizeSZ_143GHz = Ini%Read_Logical('normalizeSZ_143ghz',.false.)
   !do we want extra debug prints
   printDlSPT = Ini%Read_Logical('print_spectrum',.false.)
   printDlSPTComponents = Ini%Read_Logical('print_spectrum_components',.false.)
   if ((printDlSPT .or. printDlSPTComponents) .and. MPIRank /= 0) then
      call MPIStop('Warning - print_spectrum/print_spectrum_components is not MPI thread-safe, quitting...')
   endif
   
   if (bp_file=='' .or. desc_file=='' .or. window_folder=='' .or. cov_file=='' .or. beamerr_file == '') then
      print*,'Missing required spt hiell key: received: ',bp_file,desc_file,window_folder,cov_file,beamerr_file
      stop
   endif
   
   call this%InitSPTHiEllData(desc_file, bp_file, cov_file, beamerr_file, window_folder)
 end subroutine SPT_HiEll_ReadIni
 
 subroutine InitSPTHiEllData(this, desc_file, bp_file, cov_file, beamerr_file, window_folder)
   use IniFile
   implicit none
   class(TSPTHiEllLike) :: this
   character(LEN=Ini_max_string_len) :: desc_file, bp_file, cov_file,beamerr_file
   character(LEN=Ini_max_string_len) :: window_folder
   integer, parameter :: tmp_file_unit=82
   integer i,j,k,l,dum,lwin
   integer*4 :: neff
   integer*8 :: offset,delta
   integer*4 :: efflmin,efflmax
   real*8 :: arr(2)
   real*8 rtmp
   Type(TTextFile) :: F
   integer*4 :: errcode
   logical wexist

   !Obtain necessary info from the desc_file pertaining
   !to which freqs we're using, ell range, and number of bins per spectrum.
   call F%Open(desc_file)
   read(F%unit,*) nall,nfreq !number of bandpowers, number of freqs
   if (nfreq > MaxNFreq) &
        call MpiStop('spt initialized with more than allowed Nfrequencies')
   nband = (nfreq)*(nfreq+1)/2
   allocate(nbins(nband))
   
   do i=1,nband
      read (F%unit,*) j
      nbins(i)=j
   end do
   if (nall .ne. sum(nbins)) call MpiStop('mismatched number of bandpowers')
   maxnbin=maxval(nbins(:))
   do i=1,5
      read (F%unit,*) rtmp
      spt_norm_fr(i) = rtmp
   enddo
   if (normalizeSZ_143GHz) then 
      spt_norm_fr(5) = 143.
      print*,'using 143 as tSZ center freq'
   endif
   read (F%unit,*) spt_windows_lmin, spt_windows_lmax !Min and Max ell in window file
   do j=1,nfreq
       do i=1,5
          read (tmp_file_unit,*) rtmp
          spt_eff_fr(i,j) = rtmp
       enddo
    enddo
    do j=1,nfreq
       read (tmp_file_unit,*) spt_prefactor(j)
    enddo
   call F%Close()

   if (feedback > 1) then 
      print *, 'nall: ', nall
      print *, 'nfreq: ', nfreq
      print *, 'spt_windows_lmin: ', spt_windows_lmin
      print *, 'spt_windows_lmax: ', spt_windows_lmax
      print *, 'window_folder: ', trim(window_folder)
   endif


   allocate(this%cl_lmax(CL_T,CL_T), source=0)
   this%cl_lmax(CL_T,CL_T) = spt_windows_lmax

   if (spt_windows_lmin < 2 .or. spt_windows_lmin >= spt_windows_lmax) then
      call mpistop('Invalid lranges for sptpol')
   end if

   !ells vector is 2 ell longer in order to do cl derivatives.
   !As a result, so is cl_to_dl_conversion
   allocate( ells(spt_windows_lmin:spt_windows_lmax), &
             cl_to_dl_conversion(spt_windows_lmin:spt_windows_lmax) )

   allocate(windows(spt_windows_lmin:spt_windows_lmax,nall), &
        spec(nall))

   allocate(cov(nall,nall), beam_err(nall,nall),cov_w_beam(nall,nall))

   !Define an array with the l*(l+1)/2pi factor to convert to Dl from Cl.
   do j=spt_windows_lmin,spt_windows_lmax
      ells(j) = j
   enddo
   cl_to_dl_conversion(:) = (ells*(ells+1d0))/TWOPI

   ! Read in bandpowers
   !Should be 90x90, 90x`150, 90x220,150x150, 150x220, 220x220 in that order

   call F%Open(bp_file)
   do i=1,nall
      read (F%unit,*) dum,spec(i)
   end do
   call F%close()
   
   
   
   ! Read in covariance
   if (binaryCov) then 
      call OpenReadBinaryFile(cov_file,tmp_file_unit,nall*8_8)
      do i=1,nall
         read(tmp_file_unit,rec=i)cov(:,i)
      enddo
      close (tmp_file_unit)
   else
      call F%open(cov_file)
      do i=1,nall
         do j=1,nall
            read (F%unit,*) cov(j,i)
         end do
      end do
      call F%close()
   endif
   if (feedback > 1) print *, 'First entry of covariance matrix: ', cov(1,1)
   if (binaryBeamErr) then 
      call OpenReadBinaryFile(beamerr_file,tmp_file_unit,nall*8_8)
      do i=1,nall
         read(tmp_file_unit,rec=i)beam_err(:,i)
      enddo
      close (tmp_file_unit)
   else
      call F%open(beamerr_file)
      do i=1,nall
         do j=1,nall
            read (F%unit,*) beam_err(j,i)
         end do
      end do
      call F%close()
   endif
   if (feedback > 1) print *, 'First entry of beam correlation matrix: ', beam_err(1,1)
   
   
   
   ! Read in windows
   if (binaryWindows) then
      do i=1,nall
         inquire(FILE=trim(window_folder)//trim(numcat('window_',i)),EXIST=wexist)
         if (.not. wexist) then
            print*,'SPTpol, missing window file:', trim(window_folder)//trim(numcat('window_',i))
            call mpistop()
         endif
         call OpenReadBinaryFile(trim(window_folder)//trim(numcat('window_',i)),tmp_file_unit,2*8_8)
         j=1
         do 
            read(tmp_file_unit,rec=j, err=111) arr
            j=j+1
            l=idnint(arr(1))
            if (dabs(l-arr(1)) > 1e-4) then
               stop 'non-integer l in window file'
            endif
            if (l>=spt_windows_lmin .and. l<=spt_windows_lmax) then
               windows(l,i)=arr(2)
            endif

         enddo
111      close(tmp_file_unit)
         if (feedback > 4) then
            print*,'Spt-hiell wfile: ',trim(window_folder)//trim(numcat('window_',i))
            print*,'SPT-hiell window:',i,j,sum(windows(:,i))
         endif
      end do
   else
      do i=1,nall
         call F%Open(trim(window_folder)//trim(numcat('window_',i)))
         do j=spt_windows_lmin,spt_windows_lmax
            read (F%unit,*) dum, windows(j,i)
         end do
         call F%Close()
      end do
   end if


   SuccessfulSPTHiellInitialization = .true.

   if (feedback > 1) then
      print *, 'Successfully initialized SPT_HIELL data...'
   endif

 end subroutine InitSPTHiEllData



 function SPTHiEllLnLike(this, CMB, Theory, DataParams) 
   use MpiUtils
   implicit none
   
   class(TSPTHiEllLike) :: this
   Class(CMBParams) :: CMB
   Class(TCosmoTheoryPredictions), target :: Theory
   double precision :: DataParams(:) 
   double precision, dimension(spt_windows_lmax) :: dl_cmb
   double precision, dimension(spt_windows_lmax,7) :: component_spectra
   double precision :: PriorLnLike
   double precision :: dum
   double precision :: SPTHiEllLnLike
   double precision, parameter :: d3000 = 3000*3001/TWOPI
   double precision, parameter :: beta = 0.0012309
   double precision, parameter :: dipole_cosine = -0.4033
   double precision, dimension(1:nall) :: deltacb,cbs
   double precision, dimension(1:maxnbin) :: tmpcb
   double precision, dimension(1) :: junk, detcov
   double precision, dimension(2) :: PoissonLevels 
   double precision, dimension(2) :: ADust
   double precision, dimension(2) :: alphaDust
   double precision, dimension(3) ::  CalFactors !90, 150, 220
   double precision, dimension(10) ::  comp_arr
   type(foreground_params) :: foregroundParams
   integer :: i,j,k, l,kk, thisoffset,thisnbin
   double precision :: norm
   integer fid
   real*4, dimension(2) :: arr
   real*4, dimension(7) :: arr7
   integer*4 :: errcode
   integer, dimension(1)::ivec


   double precision, dimension(spt_windows_lmin:spt_windows_lmax) :: dl_fgs
   integer, parameter :: iFG = 4

   ! get CMB spectrum
   call Theory%ClArray(dl_cmb(:),CL_T,CL_T)      

   if (HaveForegroundsBeenInitialized() .eq. .false.) then
      write(*,*)'trying to call SPT likelihood w/o initializing foregrounds'
      call mpistop()
   endif
   if (iFG .gt. 100) call mpistop() !haven't done it yet
   CalFactors = DataParams(1:3)
   foregroundParams = GetForegroundParamsFromArray(DataParams(iFG:iFG+nForegroundParams))
   !add this to use negative for correlation

   if (printDlSPT) then
      call printForegrounds(foregroundParams)
   endif


   !$OMP PARALLEL DO  DEFAULT(NONE), &
   !$OMP  SHARED(cbs,indices,foregroundParams,spt_eff_fr,spt_norm_fr,cl_to_dl_conversion,nbins,nfreq,dl_cmb,spt_windows_lmax,spt_windows_lmin,windows,spt_prefactor,CalFactors,deltacb,offsets,spec,nband,printDlSPT,printDlSPTComponents), &
   !$OMP  private(i,j,k,dl_fgs,dl_th,tmpcb,thisnbin,l,thisoffset,fid,arr,component_spectra,comp_arr), &
   !$OMP SCHEDULE(STATIC)
   do i=1,nband
      j=indices(1,i)
      k=indices(2,i)
      thisoffset=offsets(i)
      thisnbin=nbins(i)
      tmpcb(:)=0

      !first get theory spectra
      if (printDlSPTComponents) then
         dl_fgs(spt_windows_lmin:spt_windows_lmax) = dl_foreground(foregroundParams,j,k,nfreq,spt_eff_fr, &
              spt_norm_fr,component_spectra) 
      else
         dl_fgs(spt_windows_lmin:spt_windows_lmax) = dl_foreground(foregroundParams,j,k,nfreq,spt_eff_fr, &
              spt_norm_fr) 
      endif
      !add CMB
      dl_fgs(spt_windows_lmin:spt_windows_lmax)=dl_fgs(spt_windows_lmin:spt_windows_lmax)+dl_cmb(spt_windows_lmin:spt_windows_lmax)
      


      if (printDlSPT) then
         fid=33+GetMpiRank()
         call OpenWriteBinaryFile(trim(numcat('suxp_spt_',i)),fid,4_8 * 2)
         do l=spt_windows_lmin,spt_windows_lmax
            arr(1)=l
            arr(2)=dl_fgs(l)
            write(fid,rec=l-spt_windows_lmin+1) arr(1:2)
         enddo
         close(fid)
      endif
      if (printDlSPTComponents) then
         fid=33+GetMpiRank()
         call OpenWriteBinaryFile(trim(numcat('suxp_spt_components_',i)),fid,4_8 * 10)
         do l=spt_windows_lmin,spt_windows_lmax
            comp_arr(1)=l
            comp_arr(2)=dl_fgs(l)
            comp_arr(3) = dl_cmb(l)
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
      tmpcb = tmpcb * spt_prefactor(k)*spt_prefactor(j)*CalFactors(j)*CalFactors(k)
      if (printDlSPT) then
         open(fid,file=trim(numcat('est_bandpowers_spt_',i)))
         write(fid,*)'# ',j,k
         do l=1,thisnbin
            write(fid,*)l,tmpcb(l),spec(thisoffset+l-1)
         enddo
         close(fid)
      endif

      cbs(thisoffset:thisoffset+thisnbin-1) = tmpcb(1:thisnbin) 
   enddo

   deltacb = cbs - spec

   do i=1,nall
      do j=1,nall
         cov_w_beam(i,j) = beam_err(i,j)*cbs(i)*cbs(j)
      enddo
   enddo

   cov_w_beam = cov + cov_w_beam

   SPTHiEllLnLike =  Matrix_GaussianLogLikeDouble(cov_w_beam, deltacb)
   
   if (CallFGPrior) then
      SPTHiEllLnLike = SPTHiEllLnLike + getForegroundPriorLnL(foregroundParams)
   endif


   if (feedback > 1)  then
      print *, 'SPTHiEllLnLike lnlike = ', SPTHiEllLnLike
      detcov = Matrix_GaussianLogLikeDouble(cov_w_beam, deltacb*0)
      print *, 'SPTHiEllLike chisq (after priors) = ', 2*(SPTHiEllLnLike-detcov)
   endif
 end function SPTHiEllLnLike

end module CMB_SPT_hiell_2019
