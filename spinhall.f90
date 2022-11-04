program bidebug2
    c*************************************************
    c     Read the hopping.1 file and evaluate the Kubo 
    c     formula for the spin Hall conductivity in Bismuth bilayer.
    c     The way the spin Hall conductivity is implemented here
    c     is only correct if there is no spin-flip-SOI.
    c     Otherwise it is an approximation.
    c     Orbitals 4,5,6 and 10,11,12 must be the down-orbitals.
    c     Frank Freimuth
    c
    c    Compilation on 64-bit cluster: mpiifort -CB -r8 spinhall.F -L /tmp_mnt/local/intel/mkl/current/lib/em64t -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -o spinhall.x
    c
    c     Compilation on cluster: mpif90 -r8 bidebug2.F -L /tmp_mnt/local/intel/mkl/current/lib/32 -lmkl -lmkl_lapack -lguide -lpthread -o bidebug2.x
    c    mpiifort -CB bidebug2.F -L /tmp_mnt/local/intel/mkl/current/lib/em64t -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -o bidebug2.x
    c
    c     Compilation on Juropa: mpif90 -O3 -r8 bidebug2.F -lmkl_lapack -lmkl -lz -lguide -lpthread -o bidebug2.x 
    c
    c     Compilation on Jugene: /bgsys/drivers/ppcfloor/comm/bin/mpixlf95_r -O2 -qfixed=72 -qrealsize=8 -qnosave -qdpc=e -I/bgsys/drivers/ppcfloor/comm/include -qarch=450 -qtune=450 bidebug2.F -L/opt/ibmcmp/xlmass/bg/4.4/bglib -lmassv -L/bgsys/drivers/ppcfloor/lib -L/bgsys/local/lib -L/usr/local/bg_soft/lapack/lib -llapack -L/bgsys/local/lib -lesslbg -o bidebug2.x
    c*************************************************
          implicit none
    
          complex,allocatable:: hops   (:,:,:)
          complex,allocatable:: rsnabla(:,:,:,:)
          complex,allocatable:: rspauli(:,:,:,:)
          complex,allocatable:: pauli(:,:,:)
          complex,allocatable:: paulifft(:,:,:)
          complex,allocatable:: paulifft2(:,:,:)
          real               :: rdum,idum
          integer            :: ix,iy,iz,band1,band2,h,num_wann
          integer            :: ik1,ik2,ik3
          real               :: twopi
          real               :: phas
          complex            :: fac,fac2
          complex,allocatable:: ham(:,:)
          real               :: vl,vu
          integer            :: ne,j
          real               :: abstol
          real,allocatable   :: eigvals(:)
          complex,allocatable:: eigvecs(:,:)
          integer            :: info
          complex,allocatable:: work(:)
          integer            :: lwork
          integer,allocatable:: iwork(:)
          real,allocatable   :: rwork(:)
          integer,allocatable:: ifail(:)
          real               :: kpoints(3)
          real               :: scale
          integer            :: maxhopx2,maxhopy2,maxhopz2,dire
          real,allocatable   :: fermienergy(:)
          real,allocatable   :: deviation(:)
          integer            :: grid,i1,i2,i3,i4,orb
          real,allocatable   :: conductivity(:),conductivity2(:)
          real,allocatable   :: conductivity13(:),conductivity23(:)
          real,allocatable   :: conductivity_ahe(:)
          real,allocatable   :: conductivity13_ahe(:),conductivity23_ahe(:)
          real,allocatable   :: conductivity_fsur(:)
          real,allocatable   :: conductivity13_fsur(:)
          real,allocatable   :: conductivity23_fsur(:)
          integer            :: ierr,isize,irank,kp1,kp2,kp3
          integer            :: ix1,ix2,ix3,num_occ
          complex,allocatable:: momentum(:,:,:)
          complex,allocatable:: momentum2(:,:,:)
    
          complex,allocatable:: spinmomentum(:,:,:)
          complex,allocatable:: spinmomentum2(:,:,:)
    
          complex            :: berry
    
          integer,allocatable:: sortarray(:)
          integer            :: n1,n2,n3,n4,dir
          real,allocatable   :: occupation(:),occupation2(:)
          integer,allocatable   :: nrpts(:)
          real               :: occupation_number
          integer            :: step,i,ii,num_steps
          real               :: fermi_min,fermi_max
          logical            :: l_tb,l_nabla,l_bfield
          real               :: bfield(3)
          integer            :: rvecnum,num_lines,length
          integer,allocatable:: irvec(:,:)
          real,allocatable   :: kpts(:,:)
          real               :: kder,amat(3,3)
          real               :: volume,cross(3)
          real               :: bohrincm,condq
          integer            :: maxdim,num_kpts
          real               :: minenerg,maxenerg,gamma	
          logical            :: l_bandstruc,l_fermisurf
          real,allocatable   :: magnetic(:,:)
    
          INCLUDE 'mpif.h'
          integer stt(MPI_STATUS_SIZE)
    
          CALL MPI_INIT(ierr)
    
          CALL MPI_COMM_RANK (MPI_COMM_WORLD,irank,ierr)
          CALL MPI_COMM_SIZE (MPI_COMM_WORLD,isize,ierr)
    
    
          abstol=2.0*tiny(abstol)
          twopi=2*3.141592654
    
          bohrincm=0.529*1.e-8
          condq=38.74*1.e-6
    
          if(irank.eq.0)then
             open(300,file='ahe_inp')
             read(300,*)amat(1,:)      
             read(300,*)amat(2,:)
             read(300,*)amat(3,:)
             read(300,*)fermi_min,fermi_max,num_steps
             read(300,*)grid
             read(300,*)maxdim
             read(300,*)l_tb,l_nabla
             read(300,*)l_bfield,bfield
             read(300,*)occupation_number
             read(300,*)l_fermisurf,gamma
             close(300)
          endif   
    
          call mpi_bcast(amat,9,MPI_DOUBLE_PRECISION,
         &               0,mpi_comm_world,ierr)
          call mpi_bcast(fermi_min,1,MPI_DOUBLE_PRECISION,
         &               0,mpi_comm_world,ierr)
          call mpi_bcast(occupation_number,1,MPI_DOUBLE_PRECISION,
         &               0,mpi_comm_world,ierr)
          call mpi_bcast(fermi_max,1,MPI_DOUBLE_PRECISION,
         &               0,mpi_comm_world,ierr)
          call mpi_bcast(num_steps,1,MPI_INTEGER,
         &               0,mpi_comm_world,ierr)
          call mpi_bcast(maxdim,1,MPI_INTEGER,
         &               0,mpi_comm_world,ierr)
          call mpi_bcast(grid,1,MPI_INTEGER,
         &               0,mpi_comm_world,ierr)
          call mpi_bcast(l_tb,1,MPI_LOGICAL,
         &               0,mpi_comm_world,ierr)
          call mpi_bcast(l_nabla,1,MPI_LOGICAL,
         &               0,mpi_comm_world,ierr)
          call mpi_bcast(l_bfield,1,MPI_LOGICAL,
         &               0,mpi_comm_world,ierr)
          call mpi_bcast(l_fermisurf,1,MPI_LOGICAL,
         &               0,mpi_comm_world,ierr)
          call mpi_bcast(bfield,3,MPI_DOUBLE_PRECISION,
         &               0,mpi_comm_world,ierr)
          call mpi_bcast(gamma,1,MPI_DOUBLE_PRECISION,
         &               0,mpi_comm_world,ierr)
    
          if(irank.eq.0)then
             open(200,file='hopping.1')
             num_lines=0
             num_wann=0
             do
                read(200,fmt=*,end=311)ix,iy,iz,band1,band2,rdum,idum
                num_lines=num_lines+1
                num_wann=max(num_wann,band1)
             enddo   
     311     continue
             rvecnum=num_lines/(num_wann*num_wann)
             write(*,*)"num_lines=",num_lines
             write(*,*)"rvecnum=",rvecnum
             allocate( hops(1:num_wann,1:num_wann,rvecnum) )
             allocate( irvec(3,rvecnum) )
             hops=0.0
             rewind(200)
             num_lines=0
             do 
                read(200,fmt=*,end=300)ix,iy,iz,band1,band2,rdum,idum
                num_lines=num_lines+1
                rvecnum=(num_lines-1)/(num_wann*num_wann)+1
                irvec(1,rvecnum)=ix
                irvec(2,rvecnum)=iy
                irvec(3,rvecnum)=iz
                hops( band1,band2,rvecnum )=cmplx(rdum,idum)
             enddo
     300     continue
             close(200)
    
            allocate(rspauli(1:num_wann, 1:num_wann, 3, rvecnum))
            open(400,file='./rspauli.1')
            num_lines=0
            Do
               read(400, fmt=*,end=500) ix,iy,iz,band1,band2,dir,rdum,idum
               num_lines=num_lines+1
               rvecnum=(num_lines-1)/(num_wann*num_wann*3)+1
               rspauli(band1, band2, dir, rvecnum)=cmplx(rdum,idum)
            End Do
     500    continue
            close(400)
            write(*,*) 'rvecnum',rvecnum
            allocate(nrpts(rvecnum))
            open(14,file='nrpts_inp')
            do j=1,rvecnum/15
                    read(14,'(15I5)') (nrpts(15*(j-1)+i) ,i=1,15)
            enddo
            read(14,'(<mod(rvecnum,15)>I5)') (nrpts(15*(rvecnum/15)+i),
         &                                   i=1,mod(rvecnum,15))
            close(14)
          write(*,*) nrpts
          write(*,*) 'size',size(nrpts)
    
          endif
    
          if(irank.eq.0)then
             write(*,*)"rvecnum=",rvecnum
             write(*,*)"num_wann=",num_wann
          endif
    
          call mpi_bcast(num_wann,1,MPI_INTEGER,
         &               0,mpi_comm_world,ierr)
    
          call mpi_bcast(rvecnum,1,MPI_INTEGER,
         &               0,mpi_comm_world,ierr)
    
    
          if(.not.allocated(nrpts))allocate(nrpts(rvecnum))
          call mpi_bcast(nrpts,rvecnum,MPI_INTEGER,
         &               0,mpi_comm_world,ierr)
    
          if(.not.allocated(irvec))allocate(irvec(3,rvecnum))
          length=3*rvecnum
          call mpi_bcast(irvec,length,MPI_INTEGER,
         &               0,mpi_comm_world,ierr)
    
          if(.not.allocated(hops))allocate(hops(num_wann,num_wann,rvecnum))
          length=num_wann*num_wann*rvecnum
          call mpi_bcast(hops,length,MPI_DOUBLE_COMPLEX,
         &               0,mpi_comm_world,ierr)
    
          if(l_bfield)then
             allocate(magnetic(3,num_steps))
             allocate(paulifft(num_wann,num_wann,3))
             allocate(paulifft2(num_wann,num_wann,3))
             if(.not.allocated(rspauli))
         &      allocate(rspauli(num_wann,num_wann,3,rvecnum))
             length=num_wann*num_wann*rvecnum*3
             call mpi_bcast(rspauli,length,MPI_DOUBLE_COMPLEX,
         &               0,mpi_comm_world,ierr)
          endif   
    
    
          cross(1)=amat(1,2)*amat(2,3)-amat(1,3)*amat(2,2)
          cross(2)=amat(1,3)*amat(2,1)-amat(1,1)*amat(2,3)
          cross(3)=amat(1,1)*amat(2,2)-amat(1,2)*amat(2,1)
    
          volume=cross(1)*amat(3,1)+cross(2)*amat(3,2)+cross(3)*amat(3,3)
          if(maxdim.eq.2)then
             volume=volume/amat(3,3)/(0.529e-8)
          endif   
    
          allocate( sortarray(num_steps) )
          allocate( deviation(num_steps) )
          allocate( occupation(num_steps) )
          allocate( occupation2(num_steps) )
          allocate( conductivity(num_steps) )
          allocate( conductivity2(num_steps) )
          allocate( conductivity13(num_steps) )
          allocate( conductivity23(num_steps) )
          allocate( conductivity_ahe(num_steps) )
          allocate( conductivity13_ahe(num_steps) )
          allocate( conductivity23_ahe(num_steps) )
          allocate( fermienergy(num_steps) )
    
          magnetic=0.0
          occupation=0.0
          conductivity=0.0
          conductivity13=0.0
          conductivity23=0.0
          conductivity_ahe=0.0
          conductivity13_ahe=0.0
          conductivity23_ahe=0.0
    
    
          allocate(        ham(num_wann,num_wann))
          allocate(    pauli(num_wann,num_wann,3))
          allocate(momentum (num_wann,num_wann,3))
          allocate(momentum2(num_wann,num_wann,3))
    
          allocate(spinmomentum (num_wann,num_wann,3))
          allocate(spinmomentum2(num_wann,num_wann,3))
    
          allocate(eigvals(num_wann))
          allocate(eigvecs(num_wann,num_wann))
          print*,"num_wann=",num_wann
          lwork=12.0*num_wann
          allocate( work(lwork) )
          allocate( rwork(17*num_wann) )
          allocate( iwork(15*num_wann) )
          allocate( ifail(15*num_wann) )
    
          if(l_bfield)then
            do ii=1,rvecnum
             do i=1,num_wann
              do j=1,num_wann
                hops(j,i,ii)=hops(j,i,ii)+bfield(1)*rspauli(j,i,1,ii)*0.5
         &                   +bfield(2)*rspauli(j,i,2,ii)*0.5
         &                   +bfield(3)*rspauli(j,i,3,ii)*0.5
              enddo !j
             enddo !i
            enddo
          endif
    
    
          open(13,file="Bulk_Berry")
          ik1=0
          do kp1=0,grid-1
           kpoints(1)=-0.5+real(kp1)/real(grid)
           do kp2=0,grid-1
            kpoints(2)=-0.5+real(kp2)/real(grid)
            do kp3=0,(grid-1)*(maxdim-2)
             ik1=ik1+1
             if(mod(ik1-1,isize).ne.irank)cycle
            ! kpoints(1)=-0.5+real(kp1)/real(grid)
            ! kpoints(2)=-0.5+real(kp2)/real(grid)
             kpoints(3)=-0.5+real(kp3)/real(grid)
           
            !kpoints(1)=-0.00d0*pi+(kp1-1)*2.0d0*pi/grid
            !kpoints(2)=-0.00d0*pi+(kp2-1)*2.0d0*pi/grid
            !kpoints(3)=-0.00d0*pi+(kp3-1)*2.0d0*pi/grid
    
            ! write(*,'(3i5,A,3f16.8)') kp1,kp2,kp3,"kpoint=",kpoints
             ham=0.0  
             pauli=cmplx(0.d0,0.d0)
    
             do ii=1,rvecnum
                ix=irvec(1,ii)
                iy=irvec(2,ii)
                iz=irvec(3,ii)
                phas=     iy*kpoints(2)
                phas=phas+iz*kpoints(3)
                phas=phas+ix*kpoints(1)
                phas=phas*twopi
                fac=cmplx(cos(phas),sin(phas))
                do i=1,num_wann
                 do j=1,num_wann
                  ham(j,i)=ham(j,i)+fac*hops(j,i,ii)/nrpts(ii)
                  do dir=1,3
                     pauli(j,i,dir)=pauli(j,i,dir)+fac*rspauli(j,i,dir,ii)
                  end do
                 enddo !j
                enddo !i
             enddo !ii
    
             !pauli(:,:,:)=cmplx(0.d0,0.d0)
             !Do i=1,num_wann        
             !   If ( i <= num_wann/2 ) then
             !      pauli(i,i,3)=cmplx(1.d0,0.d0)
             !   Else
             !      pauli(i,i,3)=cmplx(-1.d0,0.d0)
             !   End If
             !End Do 
    !        Call print_rmatrix("Pauli 3rd:", num_wann,num_wann,
    !    &                      DBLE(pauli(:,:,3)))
             !write(*,'(A)') "Hamiltonian:"
             !write(*,'(2f16.8)') ham
             call zheevx('V','A','U',num_wann,ham,
         &        num_wann,
         &        vl,vu,1,num_wann,abstol,ne,eigvals,eigvecs,num_wann,
         &        work,lwork,rwork,iwork,ifail,info)
             if(info.ne.0)stop 'zheevx'
             !write(*,'(A)') "eigenvects"
             !write(*,'(2f16.8)') eigvecs
    
              momentum=0.0
              do ii=1,rvecnum
                 ix=irvec(1,ii)
                 iy=irvec(2,ii)
                 iz=irvec(3,ii)
                 phas=     iy*kpoints(2)
                 phas=phas+iz*kpoints(3)
                 phas=phas+ix*kpoints(1)
                 phas=phas*twopi
                 !write(*,'(A, i5)') 'rvecno.=', ii
                 fac=cmplx(-sin(phas),cos(phas))
                 !write(*,'(A,2f16.8)') "phase=", fac
                 do dir=1,3
                   kder=amat(dir,1)*ix+amat(dir,2)*iy+amat(dir,3)*iz
                   !write(*,'(i3, A, f16.8)') dir, 'k-derivative:', kder
                   fac2=fac*kder                  
                   !write(*,*) fac2
                   do n2=1,num_wann
                     do n1=1,num_wann
                       momentum(n1,n2,dir)=
         &              momentum(n1,n2,dir)+
         &              fac2*hops(n1,n2,ii)
                     enddo
                   enddo
                 enddo
              enddo               !ii
             !write(*,'(A)') "J_cx"
             !write(*,'(2f16.8)') momentum(:,:,1)
             !write(*,'(A)') "J_cy"
             !write(*,'(2f16.8)') momentum(:,:,2)
    
              spinmomentum=0.0
              IF (.true.) then
                Do dir=1,3
                         spinmomentum(:,:,dir)=
         &                (MATMUL(pauli(:,:,3),momentum(:,:,dir))
         &                +MATMUL(momentum(:,:,dir),pauli(:,:,3)))/2.d0
                End Do
              Else
                do ii=1,rvecnum
                   ix=irvec(1,ii)
                   iy=irvec(2,ii)
                   iz=irvec(3,ii)
                   phas=     iy*kpoints(2)
                   phas=phas+iz*kpoints(3)
                   phas=phas+ix*kpoints(1)
                   phas=phas*twopi
                   fac=cmplx(-sin(phas),cos(phas))
                   do dir=1,3
                     kder=amat(dir,1)*ix+amat(dir,2)*iy+amat(dir,3)*iz
                     fac2=fac*kder                  
                     do n2=1,num_wann
                       do n1=1,num_wann                   
                         spinmomentum(n1,n2,dir)=
         &                spinmomentum(n1,n2,dir)+
         &                fac2*hops(n1,n2,ii)
                       enddo
                     enddo
                   enddo
                enddo               !ii
              End If
             !write(*,'(A)') "spinJ_cx"
             !write(*,'(2f16.8)') spinmomentum(:,:,1)
             !write(*,'(A)') "spinJ_cy"
             !write(*,'(2f16.8)') spinmomentum(:,:,2)
    
              momentum2=0.0
              do dir=1,3
                 do n2=1,num_wann
                  do n4=1,num_wann
                   do n1=1,num_wann
                    do n3=1,num_wann
                  momentum2(n4,n2,dir)=momentum2(n4,n2,dir)+
         +         momentum(n3,n1,dir)*
         *           eigvecs(n1,n2)*conjg(eigvecs(n3,n4))
                    enddo !n4
                   enddo !n3
                  enddo !n2
                 enddo !n1
              enddo  !dir
             !write(*,'(A)') "J_cx_avrg"
             !write(*,'(2f16.8)') momentum2(:,:,1)
             !write(*,'(A)') "J_cy_avrg"
             !write(*,'(2f16.8)') momentum2(:,:,2)
    
    
              spinmomentum2=0.0
              do dir=1,3
                 do n2=1,num_wann
                  do n4=1,num_wann
                   do n1=1,num_wann
                    do n3=1,num_wann
                  spinmomentum2(n4,n2,dir)=spinmomentum2(n4,n2,dir)+
         +         spinmomentum(n3,n1,dir)*
         *           eigvecs(n1,n2)*conjg(eigvecs(n3,n4))
                    enddo !n4
                   enddo !n3
                  enddo !n2
                 enddo !n1
              enddo  !dir 
    
              do step=1,num_steps
    
                fermienergy(step) = fermi_min +
         &        (fermi_max-fermi_min)*real(step)/real(num_steps)
    
    
                iy=0
                do ix=1,num_wann
                   if(eigvals(ix).le.fermienergy(step))then
                      iy=iy+1
                   endif
                enddo
    
                num_occ=iy
     
                occupation(step)=occupation(step)+num_occ
                berry=cmplx(0.d0,0.d0)
                do ik2=num_occ+1,num_wann
                 do orb=1,num_occ
                   if (abs(eigvals(orb)-eigvals(ik2)) .gt. 0.00001d0) then
                  conductivity(step)=conductivity(step)-
         +            aimag(momentum2(orb,ik2,1)*
         &          conjg(spinmomentum2(orb,ik2,2)))/
         /          (eigvals(orb)-eigvals(ik2))**2
                  conductivity13(step)=conductivity13(step)-
         +            aimag(momentum2(orb,ik2,1)*
         &          conjg(spinmomentum2(orb,ik2,3)))/
         /           (eigvals(orb)-eigvals(ik2))**2
                  conductivity23(step)=conductivity23(step)-
         +            aimag(momentum2(orb,ik2,2)*
         &          conjg(spinmomentum2(orb,ik2,3)))/
         /             (eigvals(orb)-eigvals(ik2))**2
                    ! if (step .eq. 240) then    ! set step number 200
                    !    write(13,'(2i5,4f16.8)') kp1,kp2, kpoints(1),
         &          !       kpoints(2),conductivity(step)
                    ! endif
                   end if
                 enddo !ik2
                enddo 
                ! write(*,'(3i5,1X,2f16.8)') kp1,kp2,kp3,berry
               ! if (step .eq. 240) then    ! set step number 200
               ! write(13,'(4f16.8)') kpoints(1),kpoints(2),berry
               ! endif             
     
                do ik2=num_occ+1,num_wann
                 do orb=1,num_occ
                  conductivity_ahe(step)=conductivity_ahe(step)-
         +            aimag(momentum2(orb,ik2,1)*
         &          conjg(momentum2(orb,ik2,2)))/
         /          (eigvals(orb)-eigvals(ik2))**2
                  conductivity13_ahe(step)=conductivity13_ahe(step)-
         +            aimag(momentum2(orb,ik2,1)*
         &          conjg(momentum2(orb,ik2,3)))/
         /           (eigvals(orb)-eigvals(ik2))**2
    
                  conductivity23_ahe(step)=conductivity23_ahe(step)-
         +            aimag(momentum2(orb,ik2,2)*
         &          conjg(momentum2(orb,ik2,3)))/
         /             (eigvals(orb)-eigvals(ik2))**2
    
                 enddo !ik2
                enddo 
    
    
              enddo !step
    
    
    
            enddo !kp3
           enddo !kp2
          enddo !kp1
    
             conductivity_ahe  =conductivity_ahe/
         /        grid**maxdim/volume/bohrincm*twopi*condq*2.0
             conductivity13_ahe=conductivity13_ahe/
         /        grid**maxdim/volume/bohrincm*twopi*condq*2.0
             conductivity23_ahe=conductivity23_ahe/
         /        grid**maxdim/volume/bohrincm*twopi*condq*2.0
    
    
             conductivity  =conductivity/
         /        grid**maxdim/volume/bohrincm*twopi*condq*2.0
             conductivity13=conductivity13/
         /        grid**maxdim/volume/bohrincm*twopi*condq*2.0
             conductivity23=conductivity23/
         /        grid**maxdim/volume/bohrincm*twopi*condq*2.0
    
    
    
          occupation=occupation/grid**maxdim
          magnetic=magnetic/grid**maxdim	
    
          call MPI_REDUCE(
         &        occupation,occupation2,num_steps,
         &        MPI_DOUBLE_PRECISION,MPI_SUM,0,
         &        mpi_comm_world,ierr)      
          occupation=occupation2
    
    
          call MPI_REDUCE(
         &        conductivity,conductivity2,num_steps,
         &        MPI_DOUBLE_PRECISION,MPI_SUM,0,
         &        mpi_comm_world,ierr)      
          conductivity=conductivity2
    
          call MPI_REDUCE(
         &        conductivity13,conductivity2,num_steps,
         &        MPI_DOUBLE_PRECISION,MPI_SUM,0,
         &        mpi_comm_world,ierr)      
          conductivity13=conductivity2
    
          call MPI_REDUCE(
         &        conductivity23,conductivity2,num_steps,
         &        MPI_DOUBLE_PRECISION,MPI_SUM,0,
         &        mpi_comm_world,ierr)      
          conductivity23=conductivity2
    
            call MPI_REDUCE(
         &        conductivity_ahe,conductivity2,num_steps,
         &        MPI_DOUBLE_PRECISION,MPI_SUM,0,
         &        mpi_comm_world,ierr)      
            conductivity_ahe=conductivity2
    
            call MPI_REDUCE(
         &        conductivity13_ahe,conductivity2,num_steps,
         &        MPI_DOUBLE_PRECISION,MPI_SUM,0,
         &        mpi_comm_world,ierr)      
            conductivity13_ahe=conductivity2
    
            call MPI_REDUCE(
         &        conductivity23_ahe,conductivity2,num_steps,
         &        MPI_DOUBLE_PRECISION,MPI_SUM,0,
         &        mpi_comm_world,ierr)      
            conductivity23_ahe=conductivity2
    
    
    
    
          if(irank.eq.0)then
                open(123,file='output_ahe_condquant',recl=10000)
                do step=1,num_steps
                   write(123,*)"fermienergy=",fermienergy(step)
     
                      write(123,*)"occupation=",occupation(step)
     
                   write(123,*)"conductivity=",conductivity_ahe(step)/
         &                          (0.5*77.48e-6)
                   write(123,*)"conductivity13=",conductivity13_ahe(step)/
         &                          (0.5*77.48e-6)
                   write(123,*)"conductivity23=",conductivity23_ahe(step)/
         &                          (0.5*77.48e-6)
                   write(123,*)"************************************"
                enddo !step
                close(123)
    
                open(123,file='output_she_condquant',recl=10000)
                do step=1,num_steps
                   write(123,*)"fermienergy=",fermienergy(step)
     
                      write(123,*)"occupation=",occupation(step)
     
                   write(123,*)"conductivity=",conductivity(step)/
         &                          (0.5*77.48e-6)
                   write(123,*)"conductivity13=",conductivity13(step)/
         &                          (0.5*77.48e-6)
                   write(123,*)"conductivity23=",conductivity23(step)/
         &                          (0.5*77.48e-6)
                   write(123,*)"************************************"
                enddo !step
                close(123)
    
          endif
          call mpi_barrier(mpi_comm_world,ierr)
          call MPI_Finalize(ierr)
    
    
    
    
          end program bidebug2
    
    !!*****************************************************************************************
          Subroutine PRINT_MATRIX(message, M, N, A)
    !!
    !! print out complex matrix
    !!-----------------------------------------------------------------------------
    !!
          character*(*)   message
          integer         M, N, I, J
          complex*16      A(M, N)
    
          write(*,*)
          write(*,*) message
          Do I = 1, M
             write(*,9998) ( A(I, J), J=1, N)
          End Do
    !
     9998 FORMAT( 12(:,1X,'(',F6.2,',',F6.2,')'))
          Return
          End Subroutine PRINT_MATRIX
    
    !!*****************************************************************************************
          Subroutine PRINT_RMATRIX(message, M, N, A)
    !!
    !! print out M*N real matrix 
    !!---------------------------------------------------------------------------------
    !!
          Character*(*) message
          Integer       M, N, I, J
          Real*8        A(M,N)
    
          write(*,*)
          write(*,*) message
          Do I= 1, M
             Write(*, 9998) ( A(I,J), J=1,N)
          End Do
    
     9998 FORMAT ( 18(:,1X, F6.3))
    
          End Subroutine PRINT_RMATRIX
    
    