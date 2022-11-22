PROGRAM PLOT_PWFN

  ! =================================================================!
  ! A program to plot plane wave data of electrons and positron from !
  ! pwfn.data to VASP CHGCAR format, so that electric and positronic ! 
  ! densities can be inspected with VMD.                             !
  ! =================================================================!
  
  !=============== TODO =================================================================!
  !       1) argparse-routine should regognize whether we want to plot                   !
  !          electrons, positrons or both. Both is the case now, and it is expensive.    ! 
  !       2) Currently this works only for 64 electron cell with 8 k-points              !
  !          and 5 bands, from which only 4 bands are read in. The number of bands       !
  !          should instead be deduced from the number of electrons.                     !
  !       3) argparse- routine could plot only certain k-points and bands for electrons, !
  !          just because it would be cool.                                              !
  ! =====================================================================================!

  implicit none  

  integer :: Nx,Ny,Nz,Ne,Nk,Nb,Ng,ik,ib,ie,ix,iy,iz
  real, allocatable :: Ggrid(:,:),Kvec(:,:),DEN(:,:,:),POSDEN(:,:,:)
  complex,allocatable :: CP(:),CE(:,:,:),elwf(:,:,:,:,:),poswf(:,:,:)
  real :: Rvec(3,3)
  
  ! Read the G-vector grid, plane-wave coefficients and the simulation cell parameters.
  call read_pwfn

  ! Read in the real-space grid parameters
  if(iargc()==0)then
     write(*,*) "Writing electron density with grid  30*30*30"
     Nx=30
     Ny=30
     Nz=30
  else
     call argparse 
  endif 

  allocate(elwf(Nk,Nb,Nx,Ny,Nz),poswf(Nx,Ny,Nz))
  
  ! Fourier transform the single-electron wave functions
  do ik=1,Nk
     do ib=1,4
        call fourier(CE(ik,ib,:),Kvec(:,ik),elwf(ik,ib,:,:,:))
     enddo
  enddo
  !ik=1
  !ib=1
  ! Fourier transform positron wave function
  !call fourier(CP,Kvec(:,1),poswf)

  allocate(DEN(Nx,Ny,Nz),POSDEN(Nx,Ny,Nz))

  DEN=0
  POSDEN=0

  write(*,*) "Computing densitites..."
  
  ! Loop over real space grid, calculate density (= square of the wave function)
  do ix=1,Nx ! 1st lattice vector direction
     write(*,*) "ix: ", ix
     do iy=1,Ny ! 2nd lattice vector direction
        do iz=1,Nz ! 3rd lattice vector direction
           do ik=1,Nk ! k-points
              do ib=1,Nb ! Bands
                 ! Electron density
                DEN(ix,iy,iz)=DEN(ix,iy,iz)+ABS(elwf(ik,ib,ix,iy,iz))**2
              enddo
           enddo
           ! Positron density (only in gamma point).
           !POSDEN(ix,iy,iz)=POSDEN(ix,iy,iz)+poswf(ix,iy,iz) ! NOW THIS IS WF, NOT DENSITY
        enddo
     enddo
  enddo
  
  ! Write densities into CHGCAR-format
  !call write_CHGCAR("CHGCAR_el",DEN)
  call write_CHGCAR("CHGCAR",DEN)

  call finish
  
contains

  SUBROUTINE finish
    ! Deallocates arrays.
    
    deallocate(CE,CP,elwf,poswf,Ggrid,Kvec,POSDEN,DEN)
    
  END SUBROUTINE finish
    
  SUBROUTINE write_CHGCAR(name,density)

    ! Write CHGCAR-file based on input
    ! input parameters: 
    !        - name:    Name of the CHGCAR-file (9 characters!)
    !        - density: The density to be plotted, array of dimension NX x Ny x Nz

    character(len=6), intent(in) :: name
    real,intent(IN) :: density(Nx,Ny,Nz)
    integer :: i,j,k
    real :: a2au

    a2au = 1.889725989
    
    ! Open output file
    open(unit=30, FILE=name, status="REPLACE")

    ! Header
    write(30,*) "Sicell"

    ! Lattice constant, cell parameters
    write(30,*) 5.43000000000000
    write(30,*) (Rvec(i,1)/(a2au*5.43),i=1,3)
    write(30,*) (Rvec(i,2)/(a2au*5.43),i=1,3)
    write(30,*) (Rvec(i,3)/(a2au*5.43),i=1,3)

    ! Atoms
    write(30,*) "Si"
    write(30,*) 2
    write(30,*) "Direct"
    write(30,*) 0.000000,0.000000,0.000000
    write(30,*) 0.250000,0.250000,0.250000

    write(30,*)

    ! Grid size
    write(30,*) Nx,Ny,Nz

    ! Density
    write(30,'(5F20.13)') (((REAL(density(i,j,k)),i=1,Nx),j=1,Ny),k=1,Nz)

    close(30)
    
  END SUBROUTINE write_CHGCAR
  
  SUBROUTINE fourier(coeff,Kvector, wf)

    ! Fourier transform of wave functions in a plane wave basis into real space.

    complex, intent(in) :: coeff(Ng)
    real, intent(in) :: Kvector(3,3)
    complex, intent(inout) :: wf(Nx,Ny,Nz)
    integer :: i,ii,j,k,l,count
    complex :: im
    real :: r1(3),r2(3),r3(3)
    
    im=(0,1)
    count=1
    
    ! Loop throuught the real space grid
    do i=1,Nx
       write(*,*) "k: ", ik, "/", Nk, "b: ", ib, "/", Nb, " step:" , count, "/", Nx*Ny*Nz
       do j=1,Ny
          do k=1,Nz
             count=count+1
             ii=ii+1

             r1 = (i-1)*Rvec(:,1)
             r2 = (j-1)*Rvec(:,2)
             r3 = (k-1)*Rvec(:,3)

             ! Fourier sum
             do l=1,Ng
                wf(i,j,k) = wf(i,j,k)+coeff(l)&
                     *EXP(im*DOT_PRODUCT(Ggrid(:,l)+Kvector(:,ik),r1)/Nx)&
                     *EXP(im*DOT_PRODUCT(Ggrid(:,l)+Kvector(:,ik),r2)/Ny)&
                     *EXP(im*DOT_PRODUCT(Ggrid(:,l)+Kvector(:,ik),r3)/Nz)

             enddo
          enddo
       enddo
    enddo
    
  END SUBROUTINE fourier
  
  SUBROUTINE argparse
    
    ! Parses the command line arguments
    ! TODO: add options for specializing which particles to plot,
    !       e.g. single electrons, only positrons, all electrons etc.

    integer :: i,ii,jj
    character(len=3) particle,N

    if(iargc()/=3)then
       write(*,*) "Dimensions must be 3, currently. Give three grid parameters."
       CALL EXIT(0)
    else
       do ii=1,iargc()
          call getarg(ii,N)
          read(N,*)jj
          if(ii==1)Nx=jj
          if(ii==2)Ny=jj
          if(ii==3)Nz=jj
       enddo
       write(*,*) "Grid dimensions:"
       write(*,*) Nx,Ny,Nz
    endif

  END SUBROUTINE argparse
  
  SUBROUTINE read_pwfn
  
    ! Reads in the pwfn.data file and saves the necessary data.
    
  integer :: Natoms,i,j,k,iost,dummy,dummy2
  real :: n_e
  
  open(unit=20, file="pwfn.data", status="old", action="read", form="formatted", iostat=iost)

  call skip(20,29)
  ! Number of electrons
  read(20,*) n_e
  write(Ne,*) n_e
  
  call skip(20,4)
  ! Number of atoms
  read(20,*) Natoms

  call skip(20,1+Natoms+1)
  ! Lattice vectors
  do i=1,3
     read(20,*) (Rvec(j,i),j=1,3)
  enddo

  call skip(20,4)
  ! Number of G-vectors
  read(20,*) Ng
  allocate(Ggrid(3,Ng),CP(Ng))
  
  call skip(20,1)
  ! G-vectors
  do i=1,Ng
     read(20,*) (Ggrid(j,i),j=1,3)
  enddo

  call skip(20,4)
  ! Number of k-points
  read(20,*) Nk
  Ne=Ne*Nk/2
  allocate(Kvec(3,Nk))
  
  do i=1,Nk
     call skip(20,1)
     ! Nb: number of bands
     ! Kvec: k-vector coordinates
     read(20,*) dummy, Nb, dummy2, (Kvec(j,i),j=1,3)
     if(i==1)allocate(CE(Nk,Nb,Ng))
     do j=1,Nb
        call skip(20,3)
        ! Electron coefficients
        do k=1,Ng
           read(20,*) CE(i,j,k)
        enddo ! bands
     enddo ! k-points
  enddo
  call skip(20,5)
  ! Positron coefficients
  do i=1,Ng
     read(20,*) CP(i)
  enddo

end SUBROUTINE read_pwfn


  SUBROUTINE skip(iunit,nskip)
    !==================================================!
    ! Skip lines in io file (From CASINO distribution) !
    !==================================================!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: iunit,nskip
    INTEGER iskip
    do iskip=1,nskip
       read(iunit,fmt=*,err=1,end=1)
    enddo
    return
1   write(*, *)"SKIP','Error reading file with the skip routine."
  END SUBROUTINE skip

END PROGRAM PLOT_PWFN
