PROGRAM parse

  IMPLICIT NONE

  integer :: i,j,k,l,NG,iost,ii,jj
  ! The real space grid
  integer :: Nx, Ny, Nz
  ! g-vectors, reciprocal cell and real grid 
  real, dimension(:,:), allocatable :: ggrid, recgrid, realgrid
  ! The positron density fourier coefficients
  real, dimension(:,:), allocatable :: poswf
  ! Positron density in real space
  complex, dimension(:, :, :), allocatable :: posden
  ! Imaginary number
  complex :: im
  ! Real space point
  real, dimension(3) :: r1,r2,r3
  real :: realnumber
  ! Conversion factor
  real :: a2au, au2a
  ! Pi
  real :: pi
  ! Dummy for command line arguments
  character(len=3) particle,N

  ! The real grid
  Nx=1
  Ny=1
  Nz=1
  
  if(iargc()==0)then
     write(*,*) "Writing electron density with grid  30*30*30"
  else
     call getarg(1,particle)
     if(particle/="1" .and. particle/="2".and. particle/="3")then
        write(*,*)"No approprate particle named. Quitting"
        call EXIT(0)
     endif
        write(*,*) "Writing spin "//particle//"-density."
     if(iargc()>1)then
        if(iargc()/=4)then
           write(*,*) "Dimensions must be 3, currently. Give three grid parameters."
           CALL EXIT(0)
        else
           do ii=2,iargc()
              call getarg(ii,N)
              read(N,*)jj
              if(ii==2)Nx=jj
              if(ii==3)Ny=jj
              if(ii==4)Nz=jj
           enddo
           write(*,*) "Grid dimensions:"
           write(*,*) Nx,Ny,Nz
        endif
     endif
  endif

  
  pi = 3.141592653589793238
  
  a2au = 1.889725989
  au2a = 0.529177208
  
  im = (0,1)

  allocate(posden(Nx, Ny, Nz))

  allocate(realgrid(3,3))
  
  posden=0
  
  ! Open the file including CASINO expectation value 
  open(unit=20, file="expval.data", status="old", action="read", form="formatted", iostat=iost)

  call skip(20,18)
  ! Read in the simulation cell parameters
  read(20,*) ((realgrid(i,j), i=1,3), j=1,3)

  call skip(20,1)

  ! Read how primitive cell is folded into a supercell
  read(20,*) (r1(i),i=1,3)
  realgrid(:,1) = r1(1)*realgrid(:,1)
  realgrid(:,2) = r1(2)*realgrid(:,2)
  realgrid(:,3) = r1(3)*realgrid(:,3)
  
  ! Number of G vectors
  call skip(20,11)
  read(20,*) NG

  allocate(poswf(2,NG))
  allocate(ggrid(3,NG))
  allocate(recgrid(3,3))
  
  call skip(20,1)

  ! Read reciprocal cell
  read(20,*) ((recgrid(i,j), i=1,3), j=1,3)

  call skip(20,1)

  ! Read reciprocal lattice vectors SHOULD THIS BE DONE IN A LOOP?
  do j=1,NG
     read(20,*) (ggrid(i,j), i=1,3)
  enddo

  ! For positron:
  if(trim(particle)=="3")then
     call skip(20,15+NG+7+NG+7)
  else if(trim(particle)=="1")then! For electron, spin 1:
     call skip(20,15)
  else if(trim(particle)=="2")then! For electron, spin 2:
     call skip(20,15+NG+7)
  endif
  ! Read positron density Fourier coefficients
  do j=1,NG
     read(20,*) (poswf(i,j), i=1,2)
  enddo

  ! Close density file
  close(20)
  
  ! Loop over the real grid, at each point, calculate the fourier sum of the density
  ii=0
  do i=1,Nx
     write(*,*) ii
     do j=1,Ny
        do k=1,Nz

           ii=ii+1

           r1 = (i-1)*realgrid(:,1)
           r2 = (j-1)*realgrid(:,2)
           r3 = (k-1)*realgrid(:,3)
           
           ! Fourier sum
           do l=1,NG
              posden(i,j,k) = posden(i,j,k)+(poswf(1,l)+im*poswf(2,l))&
                   *EXP(-im*DOT_PRODUCT(ggrid(:,l),r1)/Nx)&
                   *EXP(-im*DOT_PRODUCT(ggrid(:,l),r2)/Ny)&
                   *EXP(-im*DOT_PRODUCT(ggrid(:,l),r3)/Nz)
              
           enddo
        enddo
     enddo
  enddo

  ! Open output file
  open(unit=30, FILE="CHGCAR", status="REPLACE")

  ! Header
  write(30,*) "Sicell"

  ! Lattice constant, cell parameters
  write(30,*) 5.43000000000000
  write(30,*) (realgrid(i,1)/(a2au*5.43),i=1,3)
  write(30,*) (realgrid(i,2)/(a2au*5.43),i=1,3)
  write(30,*) (realgrid(i,3)/(a2au*5.43),i=1,3)

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
  write(30,'(5F20.13)') (((REAL(posden(i,j,k)),i=1,Nx),j=1,Ny),k=1,Nz)

  close(30)
  deallocate(posden,ggrid,realgrid,recgrid)

  contains
  
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
  
  
END PROGRAM parse
