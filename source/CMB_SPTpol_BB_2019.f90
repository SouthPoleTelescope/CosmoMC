 !Likelihood code used in Sayre et al 2019
  !SPTpol 500d B power spectrum
  !For questions, please contact Christian Reichardt or JT Sayre
module CMB_SPTpol_BB_2019
    use Likelihood_Cosmology
    use CosmoTheory
    use CosmologyTypes
    use FileUtils
    use MatrixUtils

    implicit none

    ! The bandpowers, in order:
    !                 150x150BB, 95x150 BB, 95 x 95 BB
    ! cov == The  bandpower covariance 
    integer :: nfreq, nband, nbin, nall, bands_per_freq
    double precision, dimension(:,:), allocatable :: spec, cov, beam_err
    double precision, dimension(:,:), allocatable :: windows
    integer :: spt_windows_lmin, spt_windows_lmax,n_beam_terms
    double precision, dimension(:), allocatable :: eff_dust_freqs,eff_sync_freqslcenter
    double precision, dimension(:,:), allocatable :: effFreqs, effSyncFreqs

    !logical, parameter :: DoBeamGrid = .false.
    integer, parameter :: N_BEAM_EXPECTED = 7
    double precision, dimension(:), allocatable :: cl_to_dl_conversion,ells, Dls_poisson, Dls_galdust,Dls_tensor
       
    double precision :: blind_abb_offset
    double precision, dimension(:), allocatable :: blind_r_offsets

    logical :: SuccessfulSPTpolBBInitialization
    logical :: SPTpol_cal_prior
    logical :: SPTpol_Add_prior
    logical :: printDlSPT
!    logical :: SPTpol_150ghzonly
!    logical :: SPTpol_95ghzonly
    logical :: SPTpol_Drop150x150ghz, SPTpol_Drop90x150ghz, SPTpol_Drop90x90ghz
    logical :: BlindR
    logical :: BlindAbb
    logical :: IndepPoisson

    !real(mcp) :: meanBeam
    !real(mcp) :: meanBcal90, sigmaBcal90
    !real(mcp) :: meanBcal150, sigmaBcal150
    real(mcp) :: meanAdd, sigmaAdd
    real(mcp), dimension(2,2) :: InvCalCov
    

    Type, extends(TCMBLikelihood) :: TSPTpolBBLike
  contains
    procedure :: ReadIni => SPTpol_BB_ReadIni
    procedure :: InitSPTpolBBData
    procedure :: LogLike => SPTpolBBLnLike
 end Type TSPTPolBBLike

contains

 subroutine SPTpol_BB_ReadIni(this, Ini)
   use IniFile
   use IniObjects
   implicit none
   class(TSPTpolBBLike) :: this
   class(TSettingIni) :: Ini
   character (LEN=Ini_max_string_len) :: desc_file
   character (LEN=Ini_max_string_len) :: bp_file, param_file
   character (LEN=Ini_max_string_len) :: cov_file,beamerr_file
   character (LEN=Ini_max_string_len) :: window_file
   character (LEN=Ini_max_string_len) :: blind_r_file
   character (LEN=Ini_max_string_len) :: blind_abb_file
   character (LEN=Ini_max_string_len) :: r_template_file
   

   BlindR = Ini%Read_Logical('sptpol_blind_r', .false.)
   BlindAbb = Ini%Read_Logical('sptpol_blind_abb', .false.)
   IndepPoisson = Ini%Read_Logical('sptpol_indep_poisson', .true.)
   
!   SPTpol_150ghzonly = Ini%Read_Logical('sptpol_150ghz_only', .false.)
!!   SPTpol_95ghzonly = Ini%Read_Logical('sptpol_95ghz_only', .false.)
!   if ( SPTpol_95ghzonly .and.  SPTpol_150ghzonly ) then
!      call mpistop('Error: SPTpol cant be 95ghz-only and 150ghz-only')
!   endif
   SPTpol_Drop150x150ghz = Ini%Read_Logical('sptpol_drop_150x150ghz', .false.)
   SPTpol_Drop90x150ghz = Ini%Read_Logical('sptpol_drop_90x150ghz', .false.)
   SPTpol_Drop90x90ghz = Ini%Read_Logical('sptpol_drop_90x90ghz', .false.)
   if (SPTpol_Drop150x150ghz .and. SPTpol_Drop90x150ghz .and. SPTpol_Drop90x90ghz) then
      call mpistop('Error: SPTpol has no bandpowers left after drop flags')
   endif


   SPTpol_cal_prior = Ini%Read_Logical('sptpol_cal_prior', .false.)
   InvCalCov(1,1) = Ini%Read_Real('sptpol_invCal_90x90', 0.0004)
   InvCalCov(1,2) = Ini%Read_Real('sptpol_invCal_90x150', 0.)
   InvCalCov(2,1) = InvCalCov(1,2) 
   InvCalCov(2,2) = Ini%Read_Real('sptpol_invCal_150x150', 0.0004)



   SPTpol_Add_prior = Ini%Read_Logical('sptpol_Add_prior', .false.)
   meanAdd = Ini%Read_Real('sptpol_meanAdd', 0.0132)
   sigmaAdd = Ini%Read_Real('sptpol_sigmaAdd', 0.0055)


   if (feedback > 1) then
      print *, 'SPTpol priors:'
      print *, 'Use cal prior: ', SPTpol_cal_prior
      print *, 'this is inverse of Pol. cal uncertainty matrix:'
      print *, InvCalCov(1,1),InvCalCov(1,2),InvCalCov(2,2)

      print *, 'Use Add prior: ', SPTpol_Add_prior
      print *, 'this is an amplitude at 150 GHz for the Galactic terms'
      print *, 'meanAdd: ', meanAdd
      print *, 'sigmaAdd: ', sigmaAdd
   endif

   param_file = Ini%Read_String_Default( 'sptpol_BB_params_file','')
   call this%loadParamNames(param_file)


   !I could simplify this by moving all the relevant info from
   !desc_file into the ini file.
   desc_file = Ini%Read_String_Default('sptpol_BB_desc_file', '')

   !do we want extra debug prints
   printDlSPT = Ini%Read_Logical('print_spectrum',.false.)
   if (printDlSPT .and. MPIRank /= 0) then
      call MPIStop('Warning - print_spectrum is not MPI thread-safe, quitting...')
   endif

   !Read in the bandpower and cov files, and determine where the windows are.
   bp_file = Ini%Read_String_Default('sptpol_BB_bp_file','')
   cov_file = Ini%Read_String_Default('sptpol_BB_cov_file','')
   window_file = Ini%Read_String_Default('sptpol_BB_window_file','')
   blind_r_file = Ini%Read_String_Default('sptpol_blind_r_file','')
   blind_abb_file = Ini%Read_String_Default('sptpol_blind_abb_file','')
   r_template_file = Ini%Read_String_Default('r_template_file','')
   if (BlindR .and. (blind_r_file == '')) &
        call mpistop('Stop: SPTpol: Told to blind r, but file doesnt exist')
   if (BlindAbb .and. (blind_abb_file == '')) &
        call mpistop('Stop: SPTpol: Told to blind Abb, but file doesnt exist')
   !also beam profiles
   beamerr_file = Ini%Read_String_Default('sptpol_BB_beam_file','')
   if (bp_file=='' .or. desc_file=='' .or. window_file=='' .or. cov_file=='' .or. beamerr_file == '') then
      print*,'Missing required sptpol key: received: ',bp_file,desc_file,window_file,cov_file,beamerr_file
      stop
   endif

   call this%InitSPTpolBBData(desc_file, bp_file, cov_file, beamerr_file, window_file,blind_r_file,blind_abb_file,r_template_file)
 end subroutine SPTpol_BB_ReadIni

 subroutine InitSPTpolBBData(this, desc_file, bp_file, cov_file, beamerr_file, window_file,blind_r_file,blind_abb_file,r_template_file)
   use IniFile
   implicit none
   class(TSPTpolBBLike) :: this
   character(LEN=Ini_max_string_len) :: desc_file, bp_file, cov_file,beamerr_file
   character(LEN=Ini_max_string_len) :: window_file,blind_r_file,blind_abb_file
   character(LEN=Ini_max_string_len) :: str,r_template_file
   integer, parameter :: tmp_file_unit=82
   integer i,j,k,l,dum,lwin,idx
   integer*4 :: neff,i0,i1
   integer*8 :: offset,delta
   integer*4 :: efflmin,efflmax
   real*8 :: arr(2),tmp
   Type(TTextFile) :: F
   integer*4 :: errcode
   logical wexist
   real(mcp) :: lcen,bb90,bb90x150,bb150
   integer ::   llo,lhi
   real(mcp) :: tt,ee,bb,te
   integer :: ltmp
       


   !Obtain necessary info from the desc_file pertaining
   !to which freqs we're using, ell range, and number of bins per spectrum.
   call F%Open(desc_file)
   read(F%unit,*) nbin,nfreq !Number of bandpowers per spectrum, Number of freqs            
   read (F%unit,*) spt_windows_lmin, spt_windows_lmax !Min and Max ell in window file
   allocate(eff_dust_freqs(nfreq))
   do i=1,nfreq
      read (F%unit,*) eff_dust_freqs(i) !Eff freq of sptpol's bands for dust, monotonically decreasing
   enddo
      allocate(eff_sync_freqs(nfreq))
   do i=1,nfreq
      read (F%unit,*) eff_sync_freqs(i) !Eff freq of sptpol's bands for radio sources, monotonically decreasing
   enddo
   call F%Close()
   
   
   if (feedback > 1) then 
      print *, 'nbin: ', nbin
      print *, 'nfreq: ', nfreq
      print *, 'spt_windows_lmin: ', spt_windows_lmin
      print *, 'spt_windows_lmax: ', spt_windows_lmax
      print *, '95  GHz band eff freq for dust: ',eff_dust_freqs(2)
      print *, '150 GHz band eff freq for dust: ',eff_dust_freqs(1)
      print *, 'window_file: ', trim(window_file)
   endif
   if (nfreq .ne. 2) call mpistop('Sorry, current code assumes 95 and 150 GHz only')

   nband = (nfreq)*(nfreq+1)/2 !With cross-freq spectra.
   nall = nbin*nband
   if (feedback > 1) then
      print *, 'nband, nbin: ', nband, nbin
      print *, 'nall: ', nall
   endif

   allocate(this%cl_lmax(CL_B,CL_B), source=0)
   this%cl_lmax(CL_B,CL_B) = spt_windows_lmax+1

   if (spt_windows_lmin < 2 .or. spt_windows_lmin >= spt_windows_lmax) then
      call mpistop('Invalid lranges for sptpol')
   end if

   !ells vector is 2 ell longer in order to do cl derivatives.
   !As a result, so is cl_to_dl_conversion
   allocate(ells(spt_windows_lmin:spt_windows_lmax), &
        cl_to_dl_conversion(spt_windows_lmin:spt_windows_lmax),&
        Dls_galdust(spt_windows_lmin:spt_windows_lmax),&
        Dls_poisson(spt_windows_lmin:spt_windows_lmax),&
        Dls_tensor(spt_windows_lmin:spt_windows_lmax)&
        )

   allocate(windows(spt_windows_lmin:spt_windows_lmax,nall), &
        spec(nbin,nband),lcenter(nbin))

   allocate(cov(1:nall,1:nall))

   
   !Populate effFreqs Array
   allocate(effFreqs(2,nband),effSyncFreqs(2,nband))
   k=0
   do i=1,nfreq
      do j=i,nfreq
         k=k+1
         effFreqs(1,k) = eff_dust_freqs(i)
         effFreqs(2,k) = eff_dust_freqs(j)
         effSyncFreqs(1,k) = eff_sync_freqs(i)
         effSyncFreqs(2,k) = eff_sync_freqs(j)
      enddo
   enddo
   deallocate(eff_dust_freqs,eff_sync_freqs) ! no longer needed

   !Define an array with the l*(l+1)/2pi factor to convert to Dl from Cl.
   do j=spt_windows_lmin,spt_windows_lmax
      ells(j) = j
   enddo
   cl_to_dl_conversion(:) = (ells*(ells+1d0))/TWOPI
   Dls_poisson(:) =  (ells*(ells+1d0))/(3000.0*3001.0)  ! normalized at ell=3000
   Dls_galdust(:) = ((ells+1d0)/81d0) * (80d0/ells)**1.42 !normalized at ell=80
   
   ! Read in bandpowers
   !Should be 150x150, 95x150, 95x95 in that order
   if (feedback > 1) print*,'reading SPTpol BB bandpowers'
   call F%Open(bp_file)
   k=1
   do  while (k .le. nbin)
      read (F%unit,'(a)',end=111) str
      str = trim(adjustl(str)) !drop leading/trailing whitespace
      if (feedback > 2) print*,'readin',trim(str)
      idx = index(str,'#')+index(str,'!')
      if (idx .ne. 1) then
         if (feedback > 2) print*,'parsing str',k
         read(str,*) lcen,llo,lhi,bb90,bb90x150,bb150
         spec(k,1)=bb150
         spec(k,2)=bb90x150
         spec(k,3)=bb90
         lcenter(k)=floor(lcen)
         if (feedback > 2) print*,'parsing str',k,'got',spec(k,:)
         k=k+1
      endif
   enddo
111 continue
   call F%close()
   if (feedback > 1) print*,'reading SPTpol BB covariance'
   ! Read in covariance
   !Should be TE, EE
   call OpenReadBinaryFile(cov_file,tmp_file_unit,nall*8_8)
   do i=1,nall
      read(tmp_file_unit,rec=i)cov(:,i)
      !        print*,i,cov(1,i)
      if (feedback > 3) print*,i,'cov(i,i)',cov(i,i)
      if (feedback > 3) print*,i,'cov(1,i)',cov(1,i)
   enddo
   close (tmp_file_unit)

   !Check if we want to fit TE or EE by itself.  If so,
   !blow up the relevant blocks of the covariance by a large number.

   if (SPTpol_Drop150x150ghz)then 
      !first zero
      do i=1,nbin
         ltmp=i
         tmp=cov(ltmp,ltmp)
         cov(ltmp,:)=0
         cov(:,ltmp)=0
         cov(ltmp,ltmp)=1d12*tmp
      enddo
      if (feedback > 1) print*,'dropping 150x150',cov
   endif
   if (SPTpol_Drop90x150ghz)then 
      !first zero
      do i=1,nbin
         ltmp=i+nbin
         tmp=cov(ltmp,ltmp)
         cov(ltmp,:)=0
         cov(:,ltmp)=0
         cov(ltmp,ltmp)=1d12*tmp
      enddo
      if (feedback > 1) print*,'dropping 90x150',cov

   endif
   if (SPTpol_Drop90x90ghz)then 
      !first zero
      do i=1,nbin
         ltmp=i+2*nbin
         tmp=cov(ltmp,ltmp)
         cov(ltmp,:)=0
         cov(:,ltmp)=0
         cov(ltmp,ltmp)=1d12*tmp
      enddo
      if (feedback > 1) print*,'dropping 90x90',cov

   endif

   if (feedback > 1) then
      print*,'Dropping 150x150',SPTpol_Drop150x150ghz
      print*,'Dropping 90x150',SPTpol_Drop90x150ghz
      print*,'Dropping 90x90',SPTpol_Drop90x90ghz
      print *, 'First entry of covariance matrix: ', cov(1,1)
   endif

   call Matrix_CholeskyDouble(cov)

   if (feedback > 1) print*,'reading SPTpol BB bpwfs'
   call OpenReadBinaryStreamFile(trim(window_file),tmp_file_unit)
   read(tmp_file_unit,pos=1)i0,i1
   if (i0 .ne. spt_windows_lmin .or. i1 .ne. spt_windows_lmax) then
      print*,'mismatched sptpol BB window ranges, quitting'
      call mpistop()
   end if
   read(tmp_file_unit,pos=1+2*4_8)windows
   close(tmp_file_unit)
   
   if(feedback > 2) print*,'sptpol bpwfs element 501',windows(501,:)

   if (feedback > 1) print*,'reading SPTpol BB beam errors'
   !get beam error term
   call OpenReadBinaryStreamFile(trim(beamerr_file),tmp_file_unit)
   read(tmp_file_unit,pos=1)neff,n_beam_terms

   print*, 'neff: ', neff
   print*, 'n_beam_terms: ', n_beam_terms
   
   if (neff .ne. nall) call mpistop('SPTpol: mismatched beam error file claimed Nbandpowers')
   if (n_beam_terms .lt. 1 .or. n_beam_terms .gt. nbin) call mpistop('SPTpol: invalid beam error file claimed Nbeam_terms')
   if (n_beam_terms .ne. N_BEAM_EXPECTED) call mpistop('SPTpol: expected  a different number of beam error terms.')
   
   allocate(beam_err(neff,n_beam_terms))

   delta =(neff)*8_8
   offset=2 * 4_8+1
   do i=1,(n_beam_terms)
      read(tmp_file_unit,pos=((i-1)*delta + offset))beam_err(:,i)
      if(feedback > 2) print*,'beam term',i,beam_err(1:5,i)
   enddo
   close(tmp_file_unit)

   allocate(blind_r_offsets(nall))
   blind_r_offsets(:)=0
   if (BlindR) then
      call OpenReadBinaryFile(trim(blind_r_file),tmp_file_unit,nall*8_8)

      read(tmp_file_unit,rec=1)blind_r_offsets(:)
      close(tmp_file_unit)
      
      if (feedback > 1) &
           print*,'SPTpol: blinding r'
   endif


   blind_abb_offset=0d0
   if (BlindAbb) then
      call OpenReadBinaryFile(trim(blind_abb_file),tmp_file_unit,8_8)
      read(tmp_file_unit,rec=1)blind_abb_offset
      close(tmp_file_unit)
      print*,'SPTpol: blinding Abb'
   endif

   Dls_tensor(:)=0
   if (r_template_file .ne. '') then
      call F%Open(r_template_file)
      do  
         read (F%unit,'(a)',end=112) str
         str = trim(adjustl(str)) !drop leading/trailing whitespace
         if (feedback > 2) print*,'readin',trim(str)
         idx = index(str,'#')+index(str,'!')
         if (idx .ne. 1) then
            if (feedback > 2) print*,'parsing str',k
            read(str,*) ltmp,tt,ee,bb,te
            if (ltmp .le. spt_windows_lmax .or. ltmp .ge. spt_windows_lmin) &
                 Dls_tensor(ltmp) = bb
         endif
      enddo
112   continue
   endif

   if (feedback > 3) then
      print*,'SPTpol R template at 3:',Dls_tensor(3)
      print*,'SPTpol R template at 503:',Dls_tensor(503)
      print*,'SPTpol R template at 2100:',Dls_tensor(2100)
      print*,'SPTpol R template at 3000:',Dls_tensor(3000)
   endif
        
   

   SuccessfulSPTpolBBInitialization = .true.

   if (feedback > 1) then
      print *, 'Successfully initialized SPTPOL_BB data...'
   endif

 end subroutine InitSPTpolBBData



 function SPTpolBBLnLike(this, CMB, Theory, DataParams) 
   implicit none
   class(TSPTpolBBLike) :: this
   Class(CMBParams) :: CMB
   Class(TCosmoTheoryPredictions), target :: Theory
   real(mcp) :: DataParams(:) 
   real(mcp), dimension(spt_windows_lmax+1) :: dls
   real(mcp) :: PriorLnLike
   real(mcp) :: dum
   real(mcp) :: SPTpolBBLnLike
   double precision, parameter :: d3000 = 3000*3001/TWOPI
   double precision, dimension(1:nall) :: deltacb,tmp2cb,BeamFac
   double precision, dimension(1:nbin) :: tmpcb
   double precision, dimension(1) :: junk, detcov
   real(mcp), dimension(3) :: PoissonLevels
   real(mcp) :: ADust,chisq
   real(mcp),dimension(3) :: CalFactors !BB
   integer :: i,j,k, l
   real(mcp) :: norm
   integer fid
   real*4, dimension(2) :: arr
   real*4, dimension(5) :: arr5
   integer*4 :: errcode
   real(mcp) :: lnl
   real(mcp) :: y1, y2

   double precision, dimension(spt_windows_lmin:spt_windows_lmax) :: dl_fgs,Dls_dust150ghz

   !DataParams for SPTpol likelihood are: [Abb, Add, Poisson150, Poisson 90x150, Poisson90, alphaRG, MapBcal150, MapBcal90, Beam Factors]
   integer, parameter :: iAbb = 1, iR=2,iConstbb=3, iAdd = iAbb+3, &
        iPoisson90=iAdd+3, iPoisson90x150=iAdd+2,iPoisson150=iAdd+1, &
        iAlphaRG=iAdd+4
        iMapBcal90=iPoisson90+3,iMapBcal150=iPoisson90+2, &
        iBeam=iMapBcal90+1
   
   
   double precision, dimension(N_BEAM_EXPECTED) :: BeamFactors

   !first we need the Dl arrays
   call Theory%ClArray(dls(:),CL_B,CL_B)      
   if (feedback > 3) print *, 'CMB (initial) at lcenter:',dls(lcenter(:))

   !scale BB to new value:
   if (DataParams(iAbb) .ne. 1) &
        dls(:) = dls(:)*(DataParams(iAbb)+blind_abb_offset)
   if (DataParams(iAbb) == 0) &
        dls(:)=0

   !add constant if desired
   if (feedback > 3) print *, 'CMB (pre const) at lcenter:',dls(lcenter(:))
   dls(:) = dls(:) + DataParams(iConstbb)
   
   !print*,'BB at 2000',dls(2000)
   !add template for r (faster way to do r-varying chains while fixing LCDM)
   if (feedback > 3) print *, 'CMB (no R) at lcenter:',dls(lcenter(:))
   dls(spt_windows_lmin:spt_windows_lmax) = dls(spt_windows_lmin:spt_windows_lmax) + DataParams(iR)*dls_tensor
   


   
   
   !print*,'input params:',DataParams

   if (feedback > 3) print *, 'CMB at lcenter:',dls(lcenter(:))

   if (IndepPoisson) then 
      PoissonLevels(:) = DataParams(iPoisson150:iPoisson90)
   else
      PoissonLevels(1) = DataParams(iPoisson150) 
      PoissonLevels(2) = DataParams(iPoisson150) * PowerLawFreqScaling(DataParams(iAlphaRG),effSyncFreqs(1,1),effSyncFreqs(:,2))
      PoissonLevels(3) = DataParams(iPoisson150) * PowerLawFreqScaling(DataParams(iAlphaRG),effSyncFreqs(1,1),effSyncFreqs(:,3))
   end if
   
   ADust = DataParams(iAdd) !150GHz dust amplitude, in Dl
   Dls_dust150ghz = ADust * Dls_galdust

   CalFactors(1) = (DataParams(iMapBcal150)*DataParams(iMapBcal150))
   CalFactors(2) = (DataParams(iMapBcal90)*DataParams(iMapBcal150))
   CalFactors(3) = (DataParams(iMapBcal90)*DataParams(iMapBcal90))
   BeamFactors = DataParams(iBeam:iBeam+N_BEAM_EXPECTED-1)

   tmpcb(:)=0
   do k=1,nband !150x150, 95x150, 95x95ghz

      !First get model foreground spectrum (in Dl).
      !Note all the foregrounds are recorded in Dl at l=3000, so we 
      !divide by d3000 to get to a normalized Cl spectrum.

      !Start with Poisson power 
      dl_fgs(:) = PoissonLevels(k) * Dls_poisson(:)
      !print*,k,'Poisson at 1000',dl_fgs(1000)
      
      !Add galactic dust
      dl_fgs(:) = dl_fgs(:) + Dls_dust150ghz * dustFreqScalingFrom150GHz(effFreqs(:,k))
      !print*,k,'Gal at 1000',dl_fgs(1000) - (PoissonLevels(k) * Dls_poisson(1000))

      !Now add model CMB.
      dl_fgs(:) = dl_fgs(:) + dls(spt_windows_lmin:spt_windows_lmax)



      !print*,'at 2000:',k,dl_fgs(2000)

      if (printDlSPT) then
         fid=33
         call OpenWriteBinaryFile(trim(numcat('like_tests/suxp_spt_',k)),fid,4_8 * 2)
         do l=spt_windows_lmin,spt_windows_lmax
            arr(1)=l
            arr(2)=dl_fgs(l)
            write(fid,rec=l-spt_windows_lmin+1) arr(1:2)
         enddo
         close(fid)

         call OpenWriteBinaryFile(trim(numcat('like_tests/suxp_spt_components_',k)),fid,4_8 * 5)
         do l=spt_windows_lmin,spt_windows_lmax
            arr5(1)=l
            arr5(2)=dl_fgs(l)
            arr5(3)=dls(l)
            arr5(4)=PoissonLevels(k)*Dls_poisson(l)
            arr5(5)=  Dls_dust150ghz(l) * dustFreqScalingFrom150GHz(effFreqs(:,k))
            write(fid,rec=l-spt_windows_lmin+1) arr5(1:5)
         enddo
         close(fid)
      endif

      !Now bin into bandpowers with the window functions.
      call dgemv('T',spt_windows_lmax-spt_windows_lmin+1,nbin,1.0d0,&
           windows(:,1+nbin*(k-1):nbin+nbin*(k-1)),spt_windows_lmax-spt_windows_lmin+1,&
           dl_fgs,1,0d0,tmpcb,1)
      
      !print*,'binned bandpowers',k,tmpcb
      
      !scale theory spectrum by calibration:
      tmpcb = tmpcb(:) / CalFactors(k)
      !print*,'calib., binned bandpowers',k,tmpcb
      !print*,'data',spec(:,k)
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
   do i=1,N_BEAM_EXPECTED
      BeamFac = BeamFac(:) * (1 + beam_err(:,i) * BeamFactors(i))
   enddo

   !print *, 'BeamFac: ', BeamFac
   
   deltacb(:) = tmp2cb(:) * BeamFac(:)

   !add r blinding if needed
   !blind r offset mimics unknown r value, always positive.
   !thus we subtract it from the model
   !blinded model r will be >= true best-fit r.
   if ( (CMB%InitPower(amp_ratio_index) .ne. 0) .or. (DataParams(iR) .ne. 0) ) &
        deltacb(:) = deltacb - blind_r_offsets

   if (feedback > 3) print *, 'modelcb: ', deltacb(:)
   
   do k=1,nband
      deltacb(1+(k-1)*nbin:nbin+(k-1)*nbin) = deltacb(1+(k-1)*nbin:nbin+(k-1)*nbin) - spec(:,k)
   enddo
   !print *, 'deltacb: ', deltacb(:)

   SPTpolBBLnLike =  Matrix_GaussianLogLikeCholDouble(cov,deltacb)
   
   if (feedback > 1) then
      detcov = Matrix_GaussianLogLikeCholDouble(cov,deltacb*0)
      !print *, 'tmpcb after CalFactors: ', tmpcb
      !print *, 'EE Bandpowers: ', spec(:,2)
      print *, 'deltacb: ', deltacb
      chisq=0
      do i=1,nall
         print *,'sigma for ',i,deltacb(i)/sqrt(cov(i,i))
         chisq=chisq+(deltacb(i)/sqrt(cov(i,i)))**2
      enddo
      print *, 'manual chisq:',chisq
      print *, 'SPTpolLnLike lnlike = ', SPTpolBBLnLike
      print *, 'SPTpolLnLike lnlike try2= ',Matrix_GaussianLogLikeCholDouble(cov,deltacb)
      print *, 'SPTpolLnLike lnlike try3= ',Matrix_GaussianLogLikeCholDouble(cov,deltacb)
      print *, 'SPTpolLnLike lnlike try4= ',Matrix_GaussianLogLikeCholDouble(cov,deltacb)
      print *, 'SPTpolLnLike BB log(det(cov))/2 = ', detcov
      print *, 'SPTpolLnLike chisq (before priors) = ', 2*(SPTpolBBLnLike-detcov)
      print *, 'spec 150x150',spec(:,1)
      print *, 'spec 90x150',spec(:,2)
      print *, 'spec 90x90',spec(:,3)
   endif
   
   !beam prior
   PriorLnLike =  0.5d0 * sum(BeamFactors**2)
   
   !            integer, parameter :: iKappa = 1, iPS_TE=2, iPS_EE=3, iDust_TE=4,iDustAlpha_TE=5, iDust_EE=6, iDustAlpha_EE=7, iMapTcal=8, iMapPcal=9, iBeam=10
   !Add Gaussian prior for temperature calibration.

   if (SPTpol_cal_prior) then
      y1 = log(DataParams(iMapBcal90))
      y2 = log(DataParams(iMapBcal150))
      !hardwired 2x2 matrix multiply:
      ! 0.5 yt C^-1 y
      PriorLnLike = PriorLnLike + 0.5d0*( InvCalCov(1,1) * y1*y1 + &
           2* InvCalCov(1,2) * y1*y2 + &
            InvCalCov(2,2) * y2*y2 )
   endif


   !Add Gaussian prior for Gal dust amplitude at 150 GHz and ell=80
   if (SPTpol_Add_prior) then
      PriorLnLike = PriorLnLike + 0.5d0*((DataParams(iAdd) - meanAdd) /sigmaAdd)**2
   endif

   if (feedback > 1) print *, 'Prior contribution to sptpol:',PriorLnLike

   SPTpolBBLnLike = SPTpolBBLnLike + PriorLnLike

   if (feedback > 1)      print *, 'SPTpolLnLike lnlike = ', SPTpolBBLnLike
   if (feedback > 1)      print *, 'SPTpolLnLike chisq (after priors) = ', 2*(SPTpolBBLnLike-detcov)
 end function SPTpolBBLnLike



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  The following are utility functions used internally
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
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
    !4.799237 x 10-11 s K
    !expect GHz
    ! so 4.799237e-2 K/GHz
    real(mcp), parameter :: hk = 4.799237e-2
    
    Bnu = (nu/nu0)**3
    Bnu = Bnu * (exp( hk*nu0/T)-1d0) / (exp( hk*nu/T)-1d0) 
    
  end function Bnu

  !dust scaling by frequency from 150x150Ghz
  function dustFreqScalingFrom150GHz(effFreqs)
    real(mcp), dimension(2), intent(in) :: effFreqs
    real(mcp) dustFreqScalingFrom150GHz
    real(mcp), parameter :: beta = 1.59, Tdust=19.6 !Kelvin

    dustFreqScalingFrom150GHz = ((effFreqs(1)*effFreqs(2))/(150d0*150d0))**beta
    dustFreqScalingFrom150GHz = dustFreqScalingFrom150GHz * &
         Bnu(effFreqs(1),150d0,Tdust) * Bnu(effFreqs(2),150d0,Tdust)
    dustFreqScalingFrom150GHz = dustFreqScalingFrom150GHz / &
         dBdT(effFreqs(1),150d0)/ dBdT(effFreqs(2),150d0)
    
  end function dustFreqScalingFrom150GHz

  function PowerLawFreqScaling(alpha,sigmasq,fr0,freqs)
    real(mcp), dimension(2), intent(in) :: freqs !output freq pair
    real(mcp), intent(in) :: alpha ! spectral index
    real(mcp), intent(in) :: sigmasq !variance in spectral index across source population
    real(mcp), intent(in) :: fr0 ! normalization frequency, ie comparing predicted power in fr0 x fr0 spectrum to fri x frj spectrum.
    real(mcp) :: PowerLawFreqScaling
    real(mcp) :: fri,frj
    fri=freqs(1)
    frj=freqs(2)
    PowerLawFreqScaling = 1/dBdT(fri,fr0)/dBdT(frj,fr0)*&
         (fri/fr0*frj/fr0)**(alpha + &
         log(fri/fr0*frj/fr0)/2 * sigmasq)


  end function PowerLawFreqScaling

  

end module CMB_SPTpol_BB_2019
