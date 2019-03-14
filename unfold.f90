
program unfold
  !-----------------------------------------------------------------------
  !
  USE constants, only: ANGSTROM_AU, eps12
  USE kinds, ONLY : DP
  USE io_files,  ONLY : prefix, tmp_dir, diropn
  USE mp_global, ONLY : npool, mp_startup,  intra_image_comm
  USE wvfct,     ONLY : nbnd, npwx, et
  USE klist,     ONLY : xk, nks, ngk, igk_k
  USE io_global, ONLY : ionode, ionode_id, stdout
  USE mp,        ONLY : mp_bcast, mp_barrier
  USE mp_world,  ONLY : world_comm, mpime, nproc
  USE wavefunctions,        ONLY : evc
  USE io_files,             ONLY : nwordwfc, iunwfc
  USE gvect, ONLY : ngm, g 
  USE noncollin_module, ONLY : npol, noncolin
  USE environment,ONLY : environment_start, environment_end
  USE fft_base,  only : dffts
  USE scatter_mod,  only : gather_grid
  USE fft_interfaces, ONLY : invfft
  USE gvect, ONLY: g
  USE cell_base, ONLY: at, alat 
  USE matrix_inversion, only: invmat

  !
  IMPLICIT NONE
  CHARACTER (len=256) :: outdir
  CHARACTER(LEN=256), external :: trimcheck
  character(len=256) :: filename
  INTEGER            :: npw, iunitout,ios,ik,i,ibnd, ig, is, ipw
  LOGICAL            :: exst
  INTEGER :: first_k, last_k, first_band, last_band, nk_sub, nbnd_sub

  REAL(DP) :: at_puc(3,3), SC_inv(3,3), SC(3,3)
  REAL(DP), ALLOCATABLE :: gvec(:,:), wnk(:,:), enk(:,:), filter(:), evc2(:,:)

  NAMELIST / inputpp / outdir, prefix, first_k, last_k, first_band, last_band, SC

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
  call mp_bcast( SC, ionode_id, world_comm )

  !
  !   Now allocate space for pwscf variables, read and check them.
  !
  CALL read_file
  call openfil_pp
  CALL init_us_1

  call invmat(3, SC, SC_inv)

  at_puc = transpose(matmul(SC_inv, at))

  if (ionode) then
    print '(5x,a)', "Supercell lattice vectors (Ang.)"
    print '(6x,a,3F12.8)', 't1 = ', at(:,1)*alat/ANGSTROM_AU
    print '(6x,a,3F12.8)', 't2 = ', at(:,2)*alat/ANGSTROM_AU
    print '(6x,a,3F12.8)', 't3 = ', at(:,3)*alat/ANGSTROM_AU
    print *, ''
    print '(5x,a)', "Unitcell lattice vectors (Ang.)"
    print '(6x,a,3F12.8)', 'a1 = ', at_puc(:,1)*alat/ANGSTROM_AU
    print '(6x,a,3F12.8)', 'a2 = ', at_puc(:,2)*alat/ANGSTROM_AU
    print '(6x,a,3F12.8)', 'a3 = ', at_puc(:,3)*alat/ANGSTROM_AU
    print *, ''
    print '(5x,a)', "Supercell Matrix S, [t1 t2 t3] = S*[a1 a2 a3]^T"
    print '(6x,3F12.0)', SC(1,:)
    print '(6x,3F12.0)', SC(2,:)
    print '(6x,3F12.0)', SC(3,:)
    print *, ''
  endif

  if (first_k <= 0) first_k = 1 
  if (last_k <= 0) last_k = nks
  if (first_band <= 0) first_band = 1 
  if (last_band <= 0) last_band = nbnd

  nk_sub = last_k - first_k + 1
  nbnd_sub = last_band - first_band + 1

  allocate(gvec(3,npwx),wnk(nbnd_sub, nk_sub),enk(nbnd_sub, nk_sub),evc2(npwx,nbnd))
  allocate(filter(npwx))
  ! IF (ionode) CALL diropn (iuwfcr, filename, lrwfcr, exst)

  DO ik = first_k, last_k

    npw = ngk(ik)
    CALL davcio (evc, 2*nwordwfc, iunwfc, ik, - 1)

    gvec(1:3,1:npw) = g(1:3, igk_k(1:npw,ik))
    call cryst_to_cart(npw, gvec, at_puc, -1)
    do ipw = 1, npw
      if (abs(gvec(1,ipw)-int(gvec(1,ipw))).lt.eps12.and.&
          abs(gvec(2,ipw)-int(gvec(2,ipw))).lt.eps12.and.&
          abs(gvec(3,ipw)-int(gvec(3,ipw))).lt.eps12) then
          filter(ipw) = 1.d0
      else
          filter(ipw) = 0.d0
      endif
    enddo

    evc2 = dble(dconjg(evc)*evc)

    wnk(:,ik-first_k+1) = 0.d0
    do ibnd = first_band, last_band
      enk(ibnd-first_band+1, ik-first_k+1) = et(ibnd,ik)
      wnk(ibnd-first_band+1, ik-first_k+1) = sum(evc2(1:npw,ibnd)*filter(1:npw))
    enddo
  enddo

  if (ionode) then
    open(11, file="wnk.dat", form="formatted", status="unknown")
    write(11, "(a,2I6)") "# nbnd, nk = ", nbnd_sub, nk_sub
    open(12, file="enk.dat", form="formatted", status="unknown")
    write(12, "(a,2I6)") "# nbnd, nk = ", nbnd_sub, nk_sub
    do ibnd = 1, nbnd_sub
      do ik = 1, nk_sub
        write(11, "(F12.6)", advance='no') wnk(ibnd, ik)
        write(12, "(F12.6)", advance='no') enk(ibnd, ik)
      enddo
      write(11, *)
      write(12, *)
    enddo
    close(11)
    close(12)
  endif

  CALL environment_end ( 'UNFOLD' )

  CALL stop_pp
  STOP

end program unfold
