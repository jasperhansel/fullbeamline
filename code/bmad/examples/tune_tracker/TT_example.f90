!+
! Program tune_tracker
!
! Driver program for tune_tracker_mod.  This module creates tune trackers based on the settings in tune_tracker.in
! Tracking is done with damping on and fluctuations off.
! BPM data at element 1 is stored and FFT analysis is done to determine horizontal and vertical tune.
!
! The following log files are opened for each tune tracker:
!   bpm_meas_<id>.out  -- turn number, x at bpm, xdot at bpm
!   kck_meas_<id>.out  -- turn number, x at kicker, xdot at kicker
!   vco_stat_<id>.out  -- turn number, vco frequency in fractional tune units
!   mod_stat_<id>.out  -- tuen number, modulator amplitude
!
! The FFT is stored in the following two files.
! The units of the x-axis of the FFT output files is fractional tune.
!   xfft.out    -- FFT of horizontal data at element 1 over last half of turns
!   yfft.out    -- FFT of vertical data at element 1 over last half of turns
!-
PROGRAM tune_tracker_driver
  USE nr
  USE bmad
  USE mode3_mod
  USE tune_tracker_mod
  IMPLICIT NONE

  TYPE(lat_struct) ring
  TYPE(coord_struct), ALLOCATABLE :: orb(:)
  TYPE(coord_struct), ALLOCATABLE :: co(:)
  TYPE(tt_param_struct) :: tt_params(max_tt) !max_tt set on tune_tracker module

  CHARACTER(20) inputfile
  CHARACTER(200) lat_file,bpm_out_file,vco_out_file,kck_out_file,mod_out_file
  CHARACTER(200) xfft_out_file, yfft_out_file

  INTEGER nTTs      ! number of tune trackers
  INTEGER i,j       ! loop counter
  INTEGER nturns    ! number of turns
  INTEGER nfft
  INTEGER n_ele_track
  INTEGER id(max_tt)

  REAL(rp) Tring  ! period of ring
  REAL(rp) bpmdata
  REAL(rp) sinphi
  REAL(rp) aFracTune, bFracTune, zFracTune, kick
  REAL(rp) Deltaphi
  REAL(rp) TTw
  REAL(rp), ALLOCATABLE :: xdata(:), ydata(:)
  COMPLEX(rp), ALLOCATABLE :: xfft(:), yfft(:)
  CHARACTER(1) i_str

  LOGICAL error

  NAMELIST /TT_in/ lat_file, nturns, nTTs, tt_params

  inputfile = "tune_tracker.in"
  OPEN(UNIT=2,FILE=inputfile)
  READ(2,NML=TT_in)
  CLOSE(2)

  ! Parse the ring and set it up for tracking
  CALL bmad_parser(lat_file,ring)
  n_ele_track = ring%n_ele_track
  CALL set_on_off(rfcavity$, ring, on$)
  bmad_com%radiation_damping_on = .true.
  bmad_com%radiation_fluctuations_on = .false.
  CALL twiss3_at_start(ring,error)
  CALL twiss3_propagate_all(ring)
  CALL closed_orbit_calc(ring,orb,6)

  ! Save Closed Orbit
  ALLOCATE(co(0:n_ele_track))
  co = orb

  ! Check if the element for each tune tracker's kicker has kick attributes.
  DO i=1,nTTs
    IF( ring%ele(tt_params(i)%kck_loc)%key == hkicker$ ) THEN
      IF( tt_params(i)%orientation /= 'h' ) THEN
        WRITE(*,'(A,I2,A)') "ERROR: Orientation of tune tracker ", i," does not match kicker orientation"
        WRITE(*,'(A)') "TERMINATING"
        STOP
      ENDIF
      IF( .not. attribute_free(tt_params(i)%kck_loc,'KICK',ring,.false.) ) THEN
        WRITE(*,'(A,I2,A)') "ERROR: kick attribute for TT ", i, " kicker is not free."
        WRITE(*,'(A)') "TERMINATING"
        STOP
      ENDIF
    ELSEIF( ring%ele(tt_params(i)%kck_loc)%key == vkicker$ ) THEN
      IF( tt_params(i)%orientation /= 'v' ) THEN
        WRITE(*,'(A,I2,A)') "ERROR: Orientation of tune tracker ", i," does not match kicker orientation"
        WRITE(*,'(A)') "TERMINATING"
        STOP
      ENDIF
      IF( .not. attribute_free(tt_params(i)%kck_loc,'KICK',ring,.false.) ) THEN
        WRITE(*,'(A,I2,A)') "ERROR: kick attribute for TT ", i, " kicker is not free."
        WRITE(*,'(A)') "TERMINATING"
        STOP
      ENDIF
    ELSE
      IF( .not. has_hkick_attributes(ring%ele(tt_params(i)%kck_loc)%key) ) THEN
        WRITE(*,'(A,I2,A)') "ERROR: Element for TT ", i, " kicker does not have kick attributes."
        WRITE(*,'(A)') "TERMINATING"
        STOP
      ENDIF
      IF( tt_params(i)%orientation == 'h' ) THEN
        IF( .not. attribute_free(tt_params(i)%kck_loc,'HKICK',ring,.false.) ) THEN
          WRITE(*,'(A,I2,A)') "ERROR: hkick attribute for TT ", i, " kicker is not free."
          WRITE(*,'(A)') "TERMINATING"
          STOP
        ENDIF
      ELSEIF( tt_params(i)%orientation == 'v' ) THEN
        IF( .not. attribute_free(tt_params(i)%kck_loc,'VKICK',ring,.false.) ) THEN
          WRITE(*,'(A,I2,A)') "ERROR: vkick attribute for TT ", i, " kicker is not free."
          WRITE(*,'(A)') "TERMINATING"
          STOP
        ENDIF
      ELSEIF( tt_params(i)%orientation == 'z' ) THEN
        IF( .not. attribute_free(tt_params(i)%kck_loc,'PHI0',ring,.false.) ) THEN
          WRITE(*,'(A,I2,A)') "ERROR: phi0 attribute for TT ", i, " kicker is not free."
          WRITE(*,'(A)') "TERMINATING"
          STOP
        ENDIF
      ENDIF
    ENDIF
  ENDDO

  !calculate ring period
  Tring = ring%ele(n_ele_track)%s / c_light

  aFracTune = MOD(ring%ele(n_ele_track)%mode3%a%phi,2.0_rp*pi)/2.0_rp / pi
  bFracTune = MOD(ring%ele(n_ele_track)%mode3%b%phi,2.0_rp*pi)/2.0_rp / pi
  zFracTune = MOD(ring%ele(n_ele_track)%mode3%c%phi,2.0_rp*pi)/2.0_rp / pi
  WRITE(*,'(A50,F10.7)') "Calculated (twiss_and_track) Fractional Tune (a): ", aFracTune
  WRITE(*,'(A50,F10.7)') "Calculated (twiss_and_track) Fractional Tune (b): ", bFracTune
  WRITE(*,'(A50,F10.7)') "Calculated (twiss_and_track) Fractional Tune (z): ", zFracTune

  !initialize each tune tracker
  DO i=1,nTTs
    !calculate initial frequency of modulator
    WRITE(*,'(A14,I3,A21,F8.5,A22)') &
         "Tune tracker #",i, " VCO base frequency: ", tt_params(i)%modTfrac0, " oscillations per turn"
    tt_params(i)%modw0 = 2.0_rp*pi/(Tring/tt_params(i)%modTfrac0)   ! 2pi/(initial period of modulator)

    ! Determine phase between BPM and kicker
    IF(tt_params(i)%orientation == 'h') THEN
      !The following calculation works because the PLL leads the BPM data by pi/2
      Deltaphi = ring%ele(n_ele_track)%mode3%a%phi - &
                 ring%ele(tt_params(i)%bpm_loc)%mode3%a%phi + ring%ele(tt_params(i)%kck_loc)%mode3%a%phi
      tt_params(i)%Onum = 1  !element of coord_struct for BPM to observe
      !Normalize kick amplitude by beta at kicker
      tt_params(i)%kickAmplitude = tt_params(i)%kickAmplitude / SQRT(ring%ele(tt_params(i)%kck_loc)%mode3%a%beta)
    ELSEIF(tt_params(i)%orientation == 'v') THEN
      !The following calculation works because the PLL leads the BPM data by pi/2
      Deltaphi = ring%ele(n_ele_track)%mode3%b%phi - &
                 ring%ele(tt_params(i)%bpm_loc)%mode3%b%phi + ring%ele(tt_params(i)%kck_loc)%mode3%b%phi
      tt_params(i)%Onum = 3  !element of coord_struct for BPM to observe
      !Normalize kick amplitude by beta at kicker
      tt_params(i)%kickAmplitude = tt_params(i)%kickAmplitude / SQRT(ring%ele(tt_params(i)%kck_loc)%mode3%b%beta)
    ELSEIF(tt_params(i)%orientation == 'z') THEN
      Deltaphi = -pi/2.0_rp + ring%ele(n_ele_track)%mode3%c%phi - &
                 ring%ele(tt_params(i)%bpm_loc)%mode3%c%phi + ring%ele(tt_params(i)%kck_loc)%mode3%c%phi
      tt_params(i)%Onum = 1  !element of coord_struct for BPM to observe
    ENDIF
    tt_params(i)%phi_to_kicker = Deltaphi
    tt_params(i)%phi_to_kicker = MOD(tt_params(i)%phi_to_kicker,2.0*pi) 

    tt_params(i)%Dt = Tring
    tt_params(i)%offset = orb(tt_params(i)%bpm_loc)%vec(tt_params(i)%Onum)

    id(i) = init_dTT(tt_params(i), orb(0))
  ENDDO

  ! Allocate data arrays for FFT
  ALLOCATE(xdata(nturns))
  ALLOCATE(ydata(nturns))
  ALLOCATE(xfft(nturns/2))
  ALLOCATE(yfft(nturns/2))

  ! Open output files for each kicker
  DO i=1,nTTs
    WRITE(i_str,'(I1)') i
    bpm_out_file = "bpm_meas_"//i_str//".out"
    kck_out_file = "kck_meas_"//i_str//".out"
    vco_out_file = "vco_stat_"//i_str//".out"
    mod_out_file = "mod_stat_"//i_str//".out"
    OPEN(UNIT=(i*100+20),FILE=bpm_out_file)
    OPEN(UNIT=(i*100+21),FILE=kck_out_file)
    OPEN(UNIT=(i*100+22),FILE=vco_out_file)
    OPEN(UNIT=(i*100+23),FILE=mod_out_file)

    write(i*100+20,'(a1,a,i6,"   ",a)') "#", "Measurements at bpm ", tt_params(i)%bpm_loc, ring%ele(tt_params(i)%bpm_loc)%name
    write(i*100+20,'(a1,a9,2a14)') "#", "turn", tt_params(i)%orientation, tt_params(i)%orientation//"'"

    write(i*100+21,'(a1,a,i6,"   ",a)') "#", "Measurements at kicker ", tt_params(i)%kck_loc, ring%ele(tt_params(i)%kck_loc)%name
    write(i*100+21,'(a1,a9,2a14)') "#", "turn", tt_params(i)%orientation, tt_params(i)%orientation//"'"

    write(i*100+22,'(a1,a9,a)') "#", "turn", "   VCO Frequency (fractional, unit)"

    write(i*100+23,'(a1,a9,a)') "#", "turn", "   Kicker modulator sinphi"
  ENDDO

  ! Print information about the location of each TT's kicker and BPM
  DO i=1,nTTs
    WRITE(*,'(A,I2,A,A,A,A)') "Kicker #", i, " at element: ", ring%ele(tt_params(i)%kck_loc)%name, &
                                  " Type: ", key_name(ring%ele(tt_params(i)%kck_loc)%key)
    WRITE(*,'(A,I2,A,A,A,A)') "   BPM #", i, " at element: ", ring%ele(tt_params(i)%bpm_loc)%name, &
                                  " Type: ", key_name(ring%ele(tt_params(i)%bpm_loc)%key)
  ENDDO


  ! Main Loop
  sinphi = 0.0_rp
  DO i=1, nturns
    CALL track_all(ring,orb)
    orb(0) = orb(n_ele_track)

    DO j=1,nTTs
      bpmdata = orb(tt_params(j)%bpm_loc)%vec(tt_params(j)%Onum)
      WRITE((j*100+23),'(I10,3ES14.6)') i, sinphi
      sinphi = TT_update(bpmdata,id(j))
      kick = sinphi*tt_params(j)%kickAmplitude
      IF( ring%ele(tt_params(j)%kck_loc)%key == hkicker$ ) THEN
        ring%ele(tt_params(j)%kck_loc)%value(kick$) = kick
      ELSEIF( ring%ele(tt_params(j)%kck_loc)%key == vkicker$ ) THEN
        ring%ele(tt_params(j)%kck_loc)%value(kick$) = kick
      ELSE
        IF(tt_params(j)%orientation == 'h') THEN
          ring%ele(tt_params(j)%kck_loc)%value(hkick$) = kick
        ELSEIF(tt_params(j)%orientation == 'v') THEN
          ring%ele(tt_params(j)%kck_loc)%value(vkick$) = kick
        ELSEIF(tt_params(j)%orientation == 'z') THEN
          ring%ele(tt_params(j)%kck_loc)%value(phi0$) = kick
        ENDIF
      ENDIF

      WRITE((j*100+20),'(I10,2ES14.6)') i, orb(tt_params(j)%bpm_loc)%vec(tt_params(j)%Onum), &
                                           orb(tt_params(j)%bpm_loc)%vec(tt_params(j)%Onum+1)
      WRITE((j*100+21),'(I10,2ES14.6)') i, orb(tt_params(j)%kck_loc)%vec(tt_params(j)%Onum), &
                 orb(tt_params(j)%kck_loc)%vec(tt_params(j)%Onum+1)-co(tt_params(j)%kck_loc)%vec(tt_params(j)%Onum+1)
      TTw = get_dTT('wf',id(j))
      WRITE((j*100+22),'(I10,ES14.6)') i, TTw*Tring/2.0_rp/pi
    ENDDO

    ! FFT data gathered at element 1
    xdata(i) = orb(1)%vec(1) - co(1)%vec(1)
    ydata(i) = orb(1)%vec(3) - co(1)%vec(3)

    IF( MOD(i,1000) == 0 ) WRITE(*,*) "Progress: ", i, "/", nturns
  ENDDO

  DO i=1,nTTs
    CLOSE(i*100+20)
    CLOSE(i*100+21)
    CLOSE(i*100+22)
    CLOSE(i*100+23)
    WRITE(*,'(A,I2,A,F10.7,A)') "Tune tracker ",i," VCO final frequency: ", &
                            get_dTT('wf',id(i))*Tring/2.0_rp/pi, " oscillations per turn"
    CALL dest_dTT(id(i),orb(n_ele_track))
  ENDDO

  ! FFT Calculation
  nfft = nturns / 2
  DO i=1,nfft
    xfft(i) = CMPLX(xdata(nturns-nfft+i),0.)
    yfft(i) = CMPLX(ydata(nturns-nfft+i),0.)
  ENDDO
  CALL four1_dp(xfft,1)
  CALL four1_dp(yfft,1)
  xfft_out_file = "xfft.out"
  yfft_out_file = "yfft.out"
  OPEN(UNIT=24,FILE=xfft_out_file)
  OPEN(UNIT=25,FILE=yfft_out_file)
  DO i=1,nfft
    WRITE(24,'(F11.7,ES14.5)') (i-1)/REAL(nfft),ABS(xfft(i))
    WRITE(25,'(F11.7,ES14.5)') (i-1)/REAL(nfft),ABS(yfft(i))
  ENDDO
  CLOSE(24)
  CLOSE(25)
  DEALLOCATE(xdata)
  DEALLOCATE(ydata)
  DEALLOCATE(xfft)
  DEALLOCATE(yfft)

END PROGRAM tune_tracker_driver





