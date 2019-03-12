
program unfold
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE io_files,  ONLY : prefix, tmp_dir, diropn
  USE mp_global, ONLY : npool, mp_startup,  intra_image_comm
  USE wvfct,     ONLY : nbnd, npwx
  USE klist,     ONLY : xk, nks, ngk, igk_k
  USE io_global, ONLY : ionode, ionode_id, stdout
  USE mp,        ONLY : mp_bcast, mp_barrier
  USE mp_world,  ONLY : world_comm
  USE wavefunctions,        ONLY : evc
  USE io_files,             ONLY : nwordwfc, iunwfc
  USE gvect, ONLY : ngm, g 
  USE noncollin_module, ONLY : npol, noncolin
  USE environment,ONLY : environment_start, environment_end
  USE fft_base,  only : dffts
  USE scatter_mod,  only : gather_grid
  USE fft_interfaces, ONLY : invfft

  !
  IMPLICIT NONE
  CHARACTER (len=256) :: outdir
  CHARACTER(LEN=256), external :: trimcheck
  character(len=256) :: filename
  INTEGER            :: npw, iunitout,ios,ik,i,ibnd, ig, is
  LOGICAL            :: exst
  INTEGER :: first_k, last_k, first_band, last_band

  NAMELIST / inputpp / outdir, prefix, first_k, last_k, first_band, last_band

  !
  CALL mp_startup ( )
  CALL environment_start ( 'UNFOLD' )

  prefix = 'pwscf'
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'

  IF ( npool > 1 ) CALL errore('bands','pools not implemented',npool)
  !
  IF ( ionode )  THEN
     !
     ! set defaults:
     first_k = 0
     last_k = 0
     first_band = 0
     last_band = 0
     !
     CALL input_from_file ( )
     ! 
     READ (5, inputpp, err = 200, iostat = ios)
200  CALL errore ('UNFOLD', 'reading inputpp namelist', ABS (ios) )
     !
     tmp_dir = trimcheck (outdir)
     ! 
  END IF

  !
  ! ... Broadcast variables
  !
  CALL mp_bcast( tmp_dir, ionode_id, world_comm )
  CALL mp_bcast( prefix, ionode_id, world_comm )
  CALL mp_bcast( first_k, ionode_id, world_comm )
  CALL mp_bcast( last_k, ionode_id, world_comm )
  CALL mp_bcast( first_band, ionode_id, world_comm )
  CALL mp_bcast( last_band, ionode_id, world_comm )

  !
  !   Now allocate space for pwscf variables, read and check them.
  !
  CALL read_file
  call openfil_pp

  exst=.false.

  if (first_k <= 0) first_k = 1 
  if (last_k <= 0) last_k = nks
  if (first_band <= 0) first_band = 1 
  if (last_band <= 0) last_band = nbnd

  CALL init_us_1

  ! IF (ionode) CALL diropn (iuwfcr, filename, lrwfcr, exst)

  DO ik = first_k, last_k
     
     npw = ngk(ik)
     CALL davcio (evc, 2*nwordwfc, iunwfc, ik, - 1)

     do ibnd = first_band, last_band

     enddo
  enddo

  CALL environment_end ( 'UNFOLD' )

  CALL stop_pp
  STOP

end program unfold
