
program unfold
  !-----------------------------------------------------------------------
  !
  USE constants, only: ANGSTROM_AU, eps12
  USE kinds, ONLY : DP
  USE mp_global, ONLY : mp_startup, mp_global_end
  USE mp_world,  ONLY : nproc, mpime, world_comm
  USE mp_pools,  ONLY : npool, nproc_pool, me_pool, my_pool_id, inter_pool_comm, intra_pool_comm
  USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
  USE io_files,  ONLY : prefix, tmp_dir, diropn, nwordwfc, iunwfc
  USE io_global, ONLY : ionode, ionode_id, stdout
  USE environment,ONLY : environment_start, environment_end

  USE ions_base, ONLY : nat, nsp, ityp, tau
  USE wvfct,     ONLY : nbnd, npwx, et
  USE klist,     ONLY : xk, nks, ngk, igk_k
  USE wavefunctions, ONLY : evc
  USE gvect, ONLY : ngm, g 
  USE noncollin_module, ONLY : npol, noncolin
  USE cell_base, ONLY: at, alat 
  USE matrix_inversion, only: invmat
  USE uspp_param, ONLY : upf

  !
  implicit none
  character (len=256) :: outdir
  character(len=256), external :: trimcheck
  character(len=256) :: filename
  integer :: npw, iunitout,ios,ik,i,ibnd, ig, is, ipw, ipol, i1, i2
  logical :: exst, any_uspp
  integer :: first_k, last_k, first_band, last_band, nk_sub, nbnd_sub

  real(dp) :: at_puc(3,3), SC_inv(3,3), SC(3,3)

  real(dp), allocatable :: wnk(:,:), enk(:,:), gvec(:,:), filter(:)

  NAMELIST / inputpp / outdir, prefix, first_k, last_k, first_band, last_band, SC

  !
  CALL mp_startup ( start_images = .false. )
  CALL environment_start ( 'UNFOLD' )

  prefix = 'pwscf'
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
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
    READ (5, inputpp, iostat = ios)
    CALL errore ('UNFOLD', 'reading inputpp namelist', ABS (ios) )
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
  CALL mp_bcast( SC, ionode_id, world_comm )

  call read_file
  call openfil_pp 

  if (npool>1) then
    call errore('unfold', 'npool>1 not allowed.', 1)
  endif
  
  any_uspp = any(upf(1:nsp)%tvanp)

  if (any_uspp) then
    call errore('unfold', 'uspp not supported.', 1)
  endif

  call invmat(3, SC, SC_inv)
  at_puc = matmul(at, SC_inv)

  if (first_k <= 0) first_k = 1 
  if (last_k <= 0) last_k = nks
  if (first_band <= 0) first_band = 1 
  if (last_band <= 0) last_band = nbnd

  nk_sub = last_k - first_k + 1
  nbnd_sub = last_band - first_band + 1

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
    print '(5x,a)', "Supercell Matrix S, [t1 t2 t3]^T = S*[a1 a2 a3]^T"
    print '(6x,3F12.0)', SC(1,:)
    print '(6x,3F12.0)', SC(2,:)
    print '(6x,3F12.0)', SC(3,:)
    print *, ''
    print '(5x,a,2I6)', 'first_k, last_k = ', first_k, last_k
    print '(5x,a,2I6)', 'first_band, last_band = ', first_band, last_band
    print *, ''
    print *, 'noncolin = ', noncolin
    print *, 'npol = ', npol
  endif

  allocate(wnk(nbnd_sub,nk_sub),enk(nbnd_sub,nk_sub))
  allocate(gvec(3,npwx))
  allocate(filter(npwx))

  if (ionode) then
    print '(5x,a)', 'Calculating wnk ...'
  endif

  do ik = first_k, last_k
    if (ionode) then
      print '(5x,a,I6)', 'k-point # ', ik
    endif

    call davcio (evc, 2*nwordwfc, iunwfc, ik, -1)
    npw = ngk(ik)

    ! G-vectors in cartesian coordinate, used in |k+G>
    gvec(1:3,1:npw) = g(1:3,igk_k(1:npw,ik))

    ! transform to cryst. coordinate of the primitive unit cell
    call cryst_to_cart(npw, gvec, at_puc, -1)
    filter = 0.d0
    do ipw = 1, npw
      if (abs(gvec(1,ipw)-int(gvec(1,ipw))).lt.eps12.and.&
          abs(gvec(2,ipw)-int(gvec(2,ipw))).lt.eps12.and.&
          abs(gvec(3,ipw)-int(gvec(3,ipw))).lt.eps12) then
          filter(ipw) = 1.d0
      else
          filter(ipw) = 0.d0
      endif
    enddo

    wnk(:,ik-first_k+1) = 0.d0
    do ibnd = first_band, last_band
      enk(ibnd-first_band+1, ik-first_k+1) = et(ibnd,ik)

      do ipol = 1, npol

        i1 = (ipol-1)*npwx+1
        i2 = i1+npw-1

        wnk(ibnd-first_band+1, ik-first_k+1) = wnk(ibnd-first_band+1, ik-first_k+1) + sum(filter(1:npw) * dble(dconjg(evc(i1:i2,ibnd))*evc(i1:i2,ibnd)))
      enddo

    enddo
  enddo

  call mp_sum (wnk, intra_pool_comm)

  if (ionode) then
    print '(5x,a)', 'Writing wnk.dat ...'
    open(11, file="wnk.dat", form="formatted", status="unknown")
    write(11, "(a,2I6)") "# nbnd, nk = ", nbnd_sub, nk_sub
    open(12, file="enk.dat", form="formatted", status="unknown")
    write(12, "(a,2I6)") "# nbnd, nk = ", nbnd_sub, nk_sub
    do ibnd = 1, nbnd_sub
      do ik = 1, nk_sub
        write(11, "(F24.12,1x)", advance='no') wnk(ibnd, ik)
        write(12, "(F12.6,1x)", advance='no') enk(ibnd, ik)
      enddo
      write(11, *)
      write(12, *)
    enddo
    close(11)
    close(12)
  endif

  CALL environment_end ( 'UNFOLD' )

  CALL stop_pp ( )

  ! this is needed because openfil is called above
  CALL close_files ( .false. )
  stop

end program unfold

