  !Likelihood code used in Henning et al 2017
  !SPTpol 500d TE and EE power spectrum
  !For questions, please contact Christian Reichardt or Jason Henning
  module CMB_SPTpol_TEEE_2017
    use Likelihood_Cosmology
    use CosmoTheory
    use CosmologyTypes
    use FileUtils
    use MatrixUtils

    implicit none

    ! The bandpowers, in order:
    !                 150x150TE, 150x150EE
    ! cov == The cholesky factored bandpower covariance 
    integer :: nfreq, nband, nbin, nall, bands_per_freq
    double precision, dimension(:,:), allocatable :: spec, cov, beam_err
    double precision, dimension(:,:), allocatable :: windows, dl_derivative, raw_spectra
    integer :: spt_windows_lmin, spt_windows_lmax,n_beam_terms
    integer, parameter :: nmc = 100000 !was 200
    !logical, parameter :: DoBeamGrid = .false.
    integer, parameter :: N_BEAM_EXPECTED = 2
    double precision, dimension(:), allocatable :: cl_to_dl_conversion,ells,rawspec_factor,deriv_factor

    logical :: SuccessfulSPTpolEEInitialization
    logical :: do_wimp
    logical :: correct_aberration
    logical :: SPTpol_Tcal_prior
    logical :: SPTpol_Pcal_prior
    logical :: SPTpol_kappa_prior
    logical :: SPTpol_alphaEE_prior
    logical :: SPTpol_alphaTE_prior
    logical :: printDlSPT
    logical :: SPTpol_EEonly
    logical :: SPTpol_TEonly

    !real(mcp) :: meanBeam
    real(mcp) :: meanTcal, sigmaTcal
    real(mcp) :: meanPcal, sigmaPcal
    real(mcp) :: meankappa, sigmakappa
    real(mcp) :: meanAlphaEE, sigmaAlphaEE
    real(mcp) :: meanAlphaTE, sigmaAlphaTE


    Type, extends(TCMBLikelihood) :: TSPTpolEELike
  contains
    procedure :: ReadIni => SPTpol_TEEE_ReadIni
    procedure :: InitSPTpolData
    procedure :: LogLike => SPTpolEELnLike
 end Type TSPTPolEELike

contains

 subroutine SPTpol_TEEE_ReadIni(this, Ini)
   use IniFile
   use IniObjects
   implicit none
   class(TSPTpolEELike) :: this
   class(TSettingIni) :: Ini
   character (LEN=Ini_max_string_len) :: desc_file
   character (LEN=Ini_max_string_len) :: bp_file, param_file
   character (LEN=Ini_max_string_len) :: cov_file,beamerr_file
   character (LEN=Ini_max_string_len) :: window_folder

   SPTpol_EEonly = Ini%Read_Logical('sptpol_EEonly', .false.)
   SPTpol_TEonly = Ini%Read_Logical('sptpol_TEonly', .false.)
   
   correct_aberration = Ini%Read_Logical('correct_aberration', .false.)
   SPTpol_Tcal_prior = Ini%Read_Logical('sptpol_tcal_prior', .false.)
   meanTcal = Ini%Read_Real('sptpol_meanTcal', 1.0)
   sigmaTcal = Ini%Read_Real('sptpol_sigmaTcal', 0.005)
   sigmaTcal = log(1+sigmaTcal)

   SPTpol_Pcal_prior = Ini%Read_Logical('sptpol_pcal_prior', .false.)
   meanPcal = Ini%Read_Real('sptpol_meanPcal', 1.0)
   sigmaPcal = Ini%Read_Real('sptpol_sigmaPcal', 0.02)
   sigmaPcal = log(1+sigmaTcal)

   SPTpol_kappa_prior = Ini%Read_Logical('sptpol_kappa_prior', .false.)
   meankappa = Ini%Read_Real('sptpol_meankappa', 0.0)
   sigmakappa = Ini%Read_Real('sptpol_sigmakappa', 0.001)

   SPTpol_alphaEE_prior = Ini%Read_Logical('sptpol_alphaEE_prior', .false.)
   meanAlphaEE = Ini%Read_Real('sptpol_meanAlphaEE', -2.42)
   sigmaAlphaEE = Ini%Read_Real('sptpol_sigmaAlphaEE', 0.02)

   SPTpol_alphaTE_prior = Ini%Read_Logical('sptpol_alphaTE_prior', .false.)
   meanAlphaTE = Ini%Read_Real('sptpol_meanAlphaTE', -2.42)
   sigmaAlphaTE = Ini%Read_Real('sptpol_sigmaAlphaTE', 0.02)

   if (feedback > 1) then
      print *, 'SPTpol priors:'
      print *, 'Use Tcal prior: ', SPTpol_Tcal_prior
      print *, 'meanTcal: ', meanTcal
      print *, 'sigmaTcal: ', sigmaTcal

      print *, 'Use Pcal prior: ', SPTpol_Pcal_prior
      print *, 'meanPcal: ', meanPcal
      print *, 'sigmaPcal: ', sigmaPcal

      print *, 'Use kappa prior: ', SPTpol_kappa_prior
      print *, 'meankappa: ', meankappa
      print *, 'sigmakappa: ', sigmakappa

      print *, 'Use alphaEE prior: ', SPTpol_alphaEE_prior
      print *, 'meanAlphaEE: ', meanAlphaEE
      print *, 'sigmaAlphaEE: ', sigmaAlphaEE

      print *, 'Use alphaTE prior: ', SPTpol_alphaTE_prior
      print *, 'meanAlphaTE: ', meanAlphaTE
      print *, 'sigmaAlphaTE: ', sigmaAlphaTE
   endif

   param_file = Ini%Read_String_Default( 'sptpol_TEEE_params_file','')
   call this%loadParamNames(param_file)


   !I could simplify this by moving all the relevant info from
   !desc_file into the ini file.
   desc_file = Ini%Read_String_Default('sptpol_TEEE_desc_file', '')

   !do we want extra debug prints
   printDlSPT = Ini%Read_Logical('print_spectrum',.false.)
   if (printDlSPT .and. MPIRank /= 0) then
      call MPIStop('Warning - print_spectrum is not MPI thread-safe, quitting...')
   endif

   !Read in the bandpower and cov files, and determine where the windows are.
   bp_file = Ini%Read_String_Default('sptpol_TEEE_bp_file','')
   cov_file = Ini%Read_String_Default('sptpol_TEEE_cov_file','')
   window_folder = Ini%Read_String_Default('sptpol_TEEE_window_dir','')
   !also beam profiles
   beamerr_file = Ini%Read_String_Default('sptpol_TEEE_beam_file','')
   if (bp_file=='' .or. desc_file=='' .or. window_folder=='' .or. cov_file=='' .or. beamerr_file == '') then
      print*,'Missing required sptpol key: received: ',bp_file,desc_file,window_folder,cov_file,beamerr_file
      stop
   endif

   call this%InitSPTpolData(desc_file, bp_file, cov_file, beamerr_file, window_folder)
 end subroutine SPTpol_TEEE_ReadIni

 subroutine InitSPTpolData(this, desc_file, bp_file, cov_file, beamerr_file, window_folder)
   use IniFile
   implicit none
   class(TSPTpolEELike) :: this
   character(LEN=Ini_max_string_len) :: desc_file, bp_file, cov_file,beamerr_file
   character(LEN=Ini_max_string_len) :: window_folder
   integer, parameter :: tmp_file_unit=82
   integer i,j,k,l,dum,lwin
   integer*4 :: neff
   integer*8 :: offset,delta
   integer*4 :: efflmin,efflmax
   real*8 :: arr(2)
   Type(TTextFile) :: F
   integer*4 :: errcode
   logical wexist

   !Obtain necessary info from the desc_file pertaining
   !to which freqs we're using, ell range, and number of bins per spectrum.
   call F%Open(desc_file)
   read(F%unit,*) nbin,nfreq !Number of bandpowers per spectrum, Number of freqs            
   read (F%unit,*) spt_windows_lmin, spt_windows_lmax !Min and Max ell in window file
   call F%Close()

   if (feedback > 1) then 
      print *, 'nbin: ', nbin
      print *, 'nfreq: ', nfreq
      print *, 'spt_windows_lmin: ', spt_windows_lmin
      print *, 'spt_windows_lmax: ', spt_windows_lmax
      print *, 'window_folder: ', trim(window_folder)
   endif
   if (nfreq .ne. 1) call mpistop('Sorry, current code wont work for multiple freqs')
   bands_per_freq = 3! Should be three for SPTpol (TT,TE,EE, although mostly ignore TT).
   if (bands_per_freq .gt. 3) call mpistop('dont have bb implemented')

   !nband = bands_per_freq*(nfreq)*(nfreq+1)/2 !With cross-freq spectra.
   nband = bands_per_freq*(nfreq)
   nall = nbin*(nfreq * (bands_per_freq - 1)) !Cov doesn't know about TT.
   if (feedback > 1) then
      print *, 'nband: ', nband
      print *, 'nall: ', nall
      print *, 'bands_per_freq: ', bands_per_freq
   endif

   allocate(this%cl_lmax(CL_E,CL_E), source=0)
   this%cl_lmax(CL_T,CL_T) = spt_windows_lmax+1
   this%cl_lmax(CL_E,CL_T) = spt_windows_lmax+1
   this%cl_lmax(CL_E,CL_E) = spt_windows_lmax+1

   if (spt_windows_lmin < 2 .or. spt_windows_lmin >= spt_windows_lmax) then
      call mpistop('Invalid lranges for sptpol')
   end if

   !ells vector is 2 ell longer in order to do cl derivatives.
   !As a result, so is cl_to_dl_conversion
   allocate(ells(spt_windows_lmin-1:spt_windows_lmax+1), &
        cl_to_dl_conversion(spt_windows_lmin-1:spt_windows_lmax+1), &
        rawspec_factor(spt_windows_lmin-1:spt_windows_lmax+1), &
        deriv_factor(spt_windows_lmin:spt_windows_lmax))

   allocate(windows(spt_windows_lmin:spt_windows_lmax,nall), &
        spec(nbin,nband))

   allocate(cov(1:nall,1:nall))

   !Define an array with the l*(l+1)/2pi factor to convert to Dl from Cl.
   do j=spt_windows_lmin-1,spt_windows_lmax+1
      ells(j) = j
   enddo
   cl_to_dl_conversion(:) = (ells*(ells+1d0))/TWOPI
   rawspec_factor = ells**3 / cl_to_dl_conversion
   deriv_factor = 0.5d0 / ells(spt_windows_lmin:spt_windows_lmax)**2

   ! Read in bandpowers
   !Should be TE, EE, TT, in that order from SPTpol analysis.

   call F%Open(bp_file)
   do i=1,nband
      do j=1,nbin
         read (F%unit,*) dum,spec(j,i)
      end do
   end do
   call F%close()


   
   ! Read in covariance
   !Should be TE, EE
   call OpenReadBinaryFile(cov_file,tmp_file_unit,nall*8_8)
   do i=1,nall
      read(tmp_file_unit,rec=i)cov(:,i)
      !        print*,i,cov(1,i)
!      if (feedback > 3) print*,i,cov(i,i)
   enddo
   close (tmp_file_unit)

   !Check if we want to fit TE or EE by itself.  If so,
   !blow up the relevant blocks of the covariance by a large number.

   print *, 'TE only: ', sptpol_TEonly
   print *, 'EE only: ', sptpol_EEonly
   
   if (sptpol_EEonly .or. sptpol_TEonly) then
      print *, 'Exploding off-diagonal cov blocks...'
      !Explode the off-diagonal blocks.
      cov(1:nbin,nbin+1:2*nbin) = cov(1:nbin,nbin+1:2*nbin)*1d12
      cov(nbin+1:2*nbin,1:nbin) = cov(nbin+1:2*nbin,1:nbin)*1d12
      
      !Explode TE auto-block if we only want EE.
      if (sptpol_EEonly) then
         print *, 'Exploding TE auto-cov block...'
         cov(1:nbin,1:nbin) = cov(1:nbin,1:nbin)*1d24
      endif

      !Explode EE auto-block if we only want TE.
      if (sptpol_TEonly) then
         print *, 'Exploding EE auto-cov block...'
         cov(nbin+1:2*nbin,nbin+1:2*nbin) = cov(nbin+1:2*nbin,nbin+1:2*nbin)*1d24
      endif
   
   endif


   print *, 'First entry of covariance matrix: ', cov(1,1)
   
   !Cholesky decompose covariance.
   ! removed since moved to the 'normal' loglike that redoes this
   call Matrix_CholeskyDouble(cov)

   ! Read in windows
   !Should be TE, EE

   !do i=1,nall
   !   inquire(FILE=trim(window_folder)//trim(numcat('window_',i)),EXIST=wexist)
   !   if (.not. wexist) then
   !      print*,'SPTpol, missing window file:', trim(window_folder)//trim(numcat('window_',i))
   !      call mpistop()
   !   endif
   !   call OpenReadBinaryFile(trim(window_folder)//trim(numcat('window_',i)),tmp_file_unit,2*8_8)
   !   j=1
   !   do 
   !      read(tmp_file_unit,rec=j, err=111) arr
   !      j=j+1
   !      l=idnint(arr(1))
   !      if (dabs(l-arr(1)) > 1e-4) then
   !         stop 'non-integer l in window file'
   !      endif
   !      if (l>=spt_windows_lmin .and. l<=spt_windows_lmax) then
   !         windows(l,i)=arr(2)
   !      endif
   !      
   !   enddo
!111   close(tmp_file_unit)
   !   if (feedback > 4) then
   !      print*,'Sptpol wfile: ',trim(window_folder)//trim(numcat('window_',i))
   !      print*,'SPTpol window:',i,j,sum(windows(:,i))
   !   endif
   !end do

   do i=1,nall
      call OpenTxtFile(trim(window_folder)//trim(numcat('window_',i)),tmp_file_unit)
      do j=spt_windows_lmin,spt_windows_lmax
         read (tmp_file_unit,*) dum, windows(j,i)
      end do
   end do


   !get beam error term
   !call OpenReadBinaryStreamFile(trim(beamerr_file),tmp_file_unit)
   !read(tmp_file_unit,pos=1)neff,n_beam_terms

   !print *, 'Beam error file: ', beamerr_file 
   neff = nall
   n_beam_terms = N_BEAM_EXPECTED
   
   print *, 'neff: ', nall
   print *, 'n_beam_terms: ', N_BEAM_EXPECTED
   
   !if (neff .ne. nbin*(bands_per_freq-1)) call mpistop('SPTpol: mismatched beam error file claimed Nbins')
   !if (n_beam_terms .lt. 1 .or. n_beam_terms .gt. nbin) call mpistop('SPTpol: invalid beam error file claimed Nbeam_terms')
   !if (n_beam_terms .ne. N_BEAM_EXPECTED) call mpistop('SPTpol: expected 7 beam error terms.')
   
   allocate(beam_err(neff,n_beam_terms))

   !delta=(neff)*8_8
   !offset=2 * 4_8+1
   !do i=1,(n_beam_terms)
   !   read(tmp_file_unit,pos=((i-1)*delta + offset))beam_err(:,i) 
   !   print*,'beam term',i,beam_err(:,i)
   !enddo
   !close(tmp_file_unit)

   call F%Open(beamerr_file)
   do i=1, N_BEAM_EXPECTED
      do j=1,neff
         read (F%unit,*) dum, beam_err(j,i)
      end do
   end do
   call F%close()

   !print*,'beam term', 1, beam_err(:,1)
   !print*,'beam term', 2, beam_err(:,2)

   SuccessfulSPTpolEEInitialization = .true.

   if (feedback > 1) then
      print *, 'Successfully initialized SPTPOL_TEEE data...'
   endif

 end subroutine InitSPTpolData



 function SPTpolEELnLike(this, CMB, Theory, DataParams) 
   implicit none
   class(TSPTpolEELike) :: this
   Class(CMBParams) :: CMB
   Class(TCosmoTheoryPredictions), target :: Theory
   real(mcp) :: DataParams(:) 
   real(mcp), dimension(spt_windows_lmax+1,2) :: dls
   real(mcp) :: PriorLnLike
   real(mcp) :: dum
   real(mcp) :: SPTpolEELnLike
   double precision, parameter :: d3000 = 3000*3001/TWOPI
   double precision, parameter :: beta = 0.0012309
   double precision, parameter :: dipole_cosine = -0.4033
   double precision, dimension(1:nall) :: deltacb,tmp2cb,BeamFac
   double precision, dimension(1:nbin) :: tmpcb
   double precision, dimension(1) :: junk, detcov
   real(mcp), dimension(2) :: PoissonLevels 
   real(mcp), dimension(2) :: ADust
   real(mcp), dimension(2) :: alphaDust
   real(mcp), dimension(3) ::  CalFactors !TT, TE, EE
   integer :: i,j,k, l,kk
   real(mcp) :: norm
   integer fid
   real*4, dimension(2) :: arr
   real*4, dimension(7) :: arr7
   integer*4 :: errcode
   real(mcp), dimension(nmc) :: lnl,beamlnl
   integer, dimension(1)::ivec
   real(mcp) :: minlnl,loclnl

   double precision, dimension(spt_windows_lmin:spt_windows_lmax) :: dl_fgs
   double precision, dimension(spt_windows_lmin:spt_windows_lmax,2) :: cl_derivative
   double precision, dimension(spt_windows_lmin:spt_windows_lmax,2) :: aberration
   double precision, dimension(spt_windows_lmin-1:spt_windows_lmax+1,2) :: raw_spectra
   !DataParams for SPTpol likelihood are: [kappa, CzeroTE, Czero_EE, Adust_TE, alpha_TE, Adust_EE, alpha_EE,MapTcal, MapPcal, BeamFac]
   integer, parameter :: iKappa = 1, iPS_TE=2, iPS_EE=3, iDust_TE=4,iDustAlpha_TE=5, iDust_EE=6, iDustAlpha_EE=7, iMapTcal=8, iMapPcal=9,iBeam=10
   double precision, dimension(N_BEAM_EXPECTED) :: BeamFactors

   !first we need the cl arrays
   !actually this function returns Dls. huh
   !so convert to Cls. We could eliminate this step by changing the later lines, but I was lazy
!!!!!!
   !Looks like we're defining the order of the cl arrays to be TE, EE here...
!!!!!!
   call Theory%ClArray(dls(:,1),CL_T,CL_E)      
   call Theory%ClArray(dls(:,2),CL_E,CL_E)


   !First calculate Cl derivatives for this position in parameter space.
   !raw_spectra = l^3 C_l
   !CAMB theory curves are in Cl not Dl!!! So we don't need to do the dl_to_cl_conversion.
   !zero it; will be set if correct_aberrations is true
   aberration(:,:)=0
   do i=1,2
      !raw_spectra(:,i) = ells(spt_windows_lmin-1:spt_windows_lmax+1)**3 * cls(spt_windows_lmin-1:spt_windows_lmax+1, i)
      raw_spectra(:,i) = rawspec_factor(:) *  dls(spt_windows_lmin-1:spt_windows_lmax+1, i)

      !The derivative at point n is roughly (raw_spectra[n+1] - raw_spectra[n-1])/2.  Then we devide by l^2 for the proper scaling for the
      !kappa parameter as described in Manzotti, et al. 2014, equation (32).
      cl_derivative(:,i) = deriv_factor(:) * (raw_spectra(spt_windows_lmin+1:spt_windows_lmax+1,i) - raw_spectra(spt_windows_lmin-1:spt_windows_lmax-1,i))

      if (correct_aberration) then
         if (feedback > 1 .and. i .eq. 1) then
            print *,'SPTpol: Correcting for aberration...'
         end if
         !Also get derives of the Dls for use with aberration corrections.
         aberration(:,i) = (dls(spt_windows_lmin+1:spt_windows_lmax+1,i) - dls(spt_windows_lmin-1:spt_windows_lmax-1,i))/2.
         aberration(:,i) = (-1*beta*dipole_cosine)* ells(spt_windows_lmin:spt_windows_lmax)*aberration(:,i)
      ENDIF
   enddo

   !DataParams for SPTpol likelihood are: [MapTcal, MapPcal, Czero_EE, CzeroTE, kappa, Adust_TE, alpha_TE, Adust_EE, alpha_EE]

   PoissonLevels(:) = DataParams(iPS_TE:iPS_EE)/d3000 !TE/EE Poisson
   ADust(1) = DataParams(iDust_TE) !TE dust amplitude, in Dl, NOT Cl yet...
   alphaDust(1) = DataParams(iDustAlpha_TE) ! TE dust spectral index (for Dl).
   ADust(2) = DataParams(iDust_EE) !EE dust amplitude, in Dl, NOT Cl yet...
   alphaDust(2) = DataParams(iDustAlpha_EE) ! EE dust spectral index (for Dl).

   !TT, TE, EE
   do i=1,3 
      CalFactors(i) = (DataParams(iMapTcal)*DataParams(iMapTcal)) * DataParams(iMapPcal)**(i-1)
   enddo
   BeamFactors = DataParams(iBeam:iBeam+N_BEAM_EXPECTED-1)

   tmpcb(:)=0
   do k=1,bands_per_freq - 1 !Don't care about TT. Just use it for leakage, which is turned off by default.

      !First get model foreground spectrum (in Cl).
      !Note all the foregrounds are recorded in Dl at l=3000, so we 
      !divide by d3000 to get to a normalized Cl spectrum.

      !Start with Poisson power and subtract the kappa parameter for super sample lensing.
      dl_fgs(:) = (PoissonLevels(k) - DataParams(iKappa)*cl_derivative(:, k) ) * cl_to_dl_conversion(spt_windows_lmin:spt_windows_lmax)

      !Now add model CMB.
      dl_fgs(:) = dl_fgs(:) + dls(spt_windows_lmin:spt_windows_lmax, k)
      !Do we want to correct for aberration?
      dl_fgs(:) = dl_fgs(:) + aberration(:,k)

      ! add dust foreground model (defined in Dl)
      dl_fgs(:) = dl_fgs(:) + Adust(k)*(ells(spt_windows_lmin:spt_windows_lmax)/80.0d0)**(alphaDust(k)+2.0d0)

      if (printDlSPT) then

         fid=33
         call OpenWriteBinaryFile(trim(numcat('like_tests/suxp_spt_',k)),fid,4_8 * 2)
         do l=spt_windows_lmin,spt_windows_lmax
            arr(1)=l
            arr(2)=dl_fgs(l)
            write(fid,rec=l-spt_windows_lmin+1) arr(1:2)
         enddo
         close(fid)

         call OpenWriteBinaryFile(trim(numcat('like_tests/suxp_spt_components_',k)),fid,4_8 * 7)
         do l=spt_windows_lmin,spt_windows_lmax
            arr7(1)=l
            arr7(2)=dl_fgs(l)
            arr7(3)=dls(l,k)
            arr7(4)=PoissonLevels(k)*cl_to_dl_conversion(l)
            arr7(5)=DataParams(iKappa)*cl_derivative(l, k) *cl_to_dl_conversion(l)
            arr7(6)= Adust(k)*(l/80d0)**(alphaDust(k)+2.0d0)
            arr7(7)=aberration(l,k)
            write(fid,rec=l-spt_windows_lmin+1) arr7(1:7)
         enddo
         close(fid)
      endif

      !Now bin into bandpowers with the window functions.
      call dgemv('T',spt_windows_lmax-spt_windows_lmin+1,nbin,1.0d0,&
           windows(:,1+nbin*(k-1):nbin+nbin*(k-1)),spt_windows_lmax-spt_windows_lmin+1,&
           dl_fgs,1,0d0,tmpcb,1)

      !print *, 'tmpcb(:) before CalFactors: ', tmpcb(:)
      
      !scale theory spectrum by calibration:
      tmpcb = tmpcb(:) / CalFactors(k+1)

      !print *, 'tmpcb after CalFactors: ', tmpcb
      
      if (printDlSPT) then
         fid=33
         open(fid,file=trim(numcat('like_tests/est_bandpowers_sptpol_',k)))
         do l=1,nbin
            write(fid,*)l,tmpcb(l),spec(l,k)
         enddo
         close(fid)
      endif
      tmp2cb(1+(k-1)*nbin:nbin+(k-1)*nbin) = tmpcb(:)
   end do
   
   BeamFac = 1d0
   do l=1,N_BEAM_EXPECTED      
      BeamFac = BeamFac(:) * (1 + beam_err(:,l) * BeamFactors(l))
   enddo

   !print *, 'BeamFac: ', BeamFac
   
   deltacb(:) = tmp2cb(:) * BeamFac(:)

   !print *, 'deltacb: ', deltacb(:)
   
   do k=1,bands_per_freq-1
      deltacb(1+(k-1)*nbin:nbin+(k-1)*nbin) = deltacb(1+(k-1)*nbin:nbin+(k-1)*nbin) - spec(:,k)
   enddo

   SPTpolEELnLike =  Matrix_GaussianLogLikeCholDouble(cov,deltacb)
   
   if (feedback > 1) then
      detcov = Matrix_GaussianLogLikeCholDouble(cov,deltacb*0)
      !print *, 'tmpcb after CalFactors: ', tmpcb
      !print *, 'EE Bandpowers: ', spec(:,2)
      !print *, 'deltacb: ', deltacb
      !print *, 'SPTpolLnLike lnlike = ', SPTpolEELnLike
      print *, 'SPTpolLnLike log(det(cov))/2 = ', detcov
      !print *, 'SPTpolLnLike chisq (before priors) = ', 2*(SPTpolEELnLike-detcov)
   endif
   
   !beam prior
   PriorLnLike =  0.5d0 * sum(BeamFactors**2)
   
   !            integer, parameter :: iKappa = 1, iPS_TE=2, iPS_EE=3, iDust_TE=4,iDustAlpha_TE=5, iDust_EE=6, iDustAlpha_EE=7, iMapTcal=8, iMapPcal=9, iBeam=10
   !Add Gaussian prior for temperature calibration.
   if (SPTpol_Tcal_prior) then
      PriorLnLike = PriorLnLike + 0.5d0*(log(DataParams(iMapTcal) / meanTcal)/sigmaTcal)**2
   endif

   !Add Gaussian prior for polarization efficiency.
   if (SPTpol_Pcal_prior) then
      PriorLnLike = PriorLnLike +   0.5d0*(log(DataParams(iMapPcal) / meanPcal)/sigmaPcal)**2
   endif

   !Add Gaussian prior for kappa.
   if (SPTpol_kappa_prior) then
      PriorLnLike = PriorLnLike + 0.5d0*((DataParams(iKappa) - meankappa) /sigmakappa)**2
   endif

   !Add Gaussian prior for alphaTE.
   if (SPTpol_alphaTE_prior) then
      PriorLnLike = PriorLnLike + 0.5d0*((DataParams(iDustAlpha_TE) - meanAlphaTE) /sigmaAlphaTE)**2
   endif

   !Add Gaussian prior for alphaEE.
   if (SPTpol_alphaEE_prior) then
      PriorLnLike = PriorLnLike + 0.5d0*((DataParams(iDustAlpha_EE) - meanAlphaEE) /sigmaAlphaEE)**2
   endif

   if (feedback > 1) print *, 'Prior contribution to sptpol:',PriorLnLike

   SPTpolEELnLike = SPTpolEELnLike + PriorLnLike

   if (feedback > 1)      print *, 'SPTpolLnLike lnlike = ', SPTpolEELnLike
   if (feedback > 1)      print *, 'SPTpolLnLike chisq (after priors) = ', 2*(SPTpolEELnLike-detcov)
 end function SPTpolEELnLike



 subroutine OpenTxtFile(aname, aunit)

   character(LEN=*), intent(IN) :: aname
   integer, intent(in) :: aunit     
   open(unit=aunit,file=aname,form='formatted',status='old', action='read', err=500)
   return                                                                        
500 call MpiStop('File not found: '//trim(aname))                               
 end subroutine OpenTxtFile

 subroutine OpenReadBinaryStreamFile(aname,aunit)
   character(LEN=*), intent(IN) :: aname
   integer, intent(in) :: aunit
   open(unit=aunit,file=aname,form='unformatted',access='stream', err=500)
   ! be aware the data is in LE: ,convert='LITTLE_ENDIAN')
   return

500 call MpiStop('File not found: '//trim(aname))
 end subroutine OpenReadBinaryStreamFile

 subroutine OpenReadBinaryFile(aname,aunit,record_length)
   character(LEN=*), intent(IN) :: aname
   integer, intent(in) :: aunit
   integer*8,intent(in) :: record_length
   open(unit=aunit,file=aname,form='unformatted',access='direct',recl=record_length,  err=500)
   return

500 call MpiStop('File not found: '//trim(aname))
 end subroutine openReadBinaryFile

 subroutine OpenWriteBinaryFile(aname,aunit,record_length)
   character(LEN=*), intent(IN) :: aname
   integer, intent(in) :: aunit
   integer*8,intent(in) :: record_length
   open(unit=aunit,file=aname,form='unformatted',status='replace',access='direct',recl=record_length, err=500)
   return

500 call MpiStop('File not able to be written to: '//trim(aname))
 end subroutine OpenWriteBinaryFile

 subroutine Matrix_CholeskyDouble(M, err)
   !Note upper triangular is not zeroed
   double precision, intent(inout):: M(:,:)
   integer n, info
   integer, optional :: err

   n=Size(M,DIM=1)
   if (Size(M,DIM=2)/=n) call MpiStop('Matrix_CholeskyDouble: non-square matrix')

   call dpotrf ('L', n, M, n, info)
   if (present(err)) then
      err = info
   else
      if (info/=0) &
           call MpiStop('Matrix_CholeskyDouble: not positive definite '//trim(IntToStr(info)))
   end if
 end subroutine Matrix_CholeskyDouble
 !    
 ! Returns -Log Likelihood for Gaussian: (d^T Cov^{-1} d + log|Cov|)/2     
 ! Cov is already the cholesky factorization of the covariance matrix   
 !                       
 function Matrix_GaussianLogLikeCholDouble(Cov, d) result(LogLike)
   double precision, intent(inout):: Cov(:,:)
   double precision, intent(in):: d(:)
   double precision, allocatable :: tmp(:)
   double precision :: LogLike
   integer info,i,n

   call Matrix_Start('GaussianLogLikeCholDouble')
   n = size(COV,DIM=1)
   if (Size(COV,DIM=2)/=n) call MpiStop('Matrix_GaussianLogLikeCholDouble: non-square matrix')
   if (Size(d)/=n) call MpiStop('Matrix_GaussianLogLikeCholDouble: covariance and d different size')

   LogLike = 0
   !Log Det term:                                  
   do i=1, n
      !Don't divide det term by 2 since we're dealing with the Cholesky-factored cov,
      !not the whole matrix.
      LogLike = LogLike  + log(Cov(i,i))
   end do


   !Solve for Cov^{-1}d [could use faster symmetric method]            
   allocate(tmp(n))
   tmp = d
   call DPOTRS('L', n, 1, Cov, n, tmp, n, INFO )
   if (INFO/=0) call MpiStop('Matrix_GaussianLogLikeCholDouble: error in solving for cov^{-1}d')

   !Add together             \

   LogLike = LogLike + dot_product(tmp,d)/2._dm
   deallocate(tmp)

   call Matrix_End('GaussianLogLikeCholDouble')

 end function Matrix_GaussianLogLikeCholDouble



end module CMB_SPTpol_TEEE_2017
