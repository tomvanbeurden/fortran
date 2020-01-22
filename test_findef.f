C----|--1---------2---------3---------4---------5---------6---------7-|
     
      PROGRAM TEST

c     cm(1)= bulk modulus
c     cm(2)= shear modulus
c     cm(3)=Young’s modulus along tt−axis
c     cm(4)=Young’s modulus along ll−axis
c     cm(5)=Young’s modulus along ww−axis
c     cm(6)=Young’s modulus along tl−axis
c     cm(7)=Young’s modulus along lw−axis
c     cm(8)=Young’s modulus along tw−axis
c     cm(9 )= Initial Yield stress tt−axis
c     cm(10)= Initial Yield stress tl−axis
c     cm(11)= Initial Yield stress tw−axis
c     cm(12)= Plateau stress tt−axis
c     cm(13)= Plateau stress tl−axis
c     cm(14)= Plateau stress tw−axis
c     cm(15)= Non−local influence radius
c     cm(16)= Element size along tt−axis
c     cm(17)= Densification strain
c     cm(18)= Hardening modulus
c     cm(19)= Number of utilized history variables
c     cm(20)= Softening shape parameter
c     cm(21)= Poisson’s ratio on tl−axis
c     cm(22)= Poisson’s ratio on lw−axis
c     cm(23)= Poisson’s ratio on tw−axis


      


      IMPLICIT NONE
c      double precision :: cm(23) = 0.0d0
      real :: cm(23) = 0.0d0
      real :: eps(6) = 0.0d0
      real :: sig(6) = 0.0d0
      real :: epsp(3) = 0.0d0
      real :: hsv(29) = 0.0d0
      real :: dt1 = 0.1d0
      integer :: nnpcrv(17) = -1
      integer :: imax = 0
      character*5 :: etype = "solid" 
      logical :: failel = .false.
      logical :: reject = .false.
      integer :: idele = -1
      integer :: num_hv
      integer :: elsiz = -1
      real :: pi
      real :: alpha
      real :: u = 0.0d0, u_x = 0.0d0, u_y = 0.0d0
      real :: u_list (100000)
      real :: sigma11_list(100000)

      integer :: i=1
      print *, "Hello from main" 

      cm(1)= 110000000.0d0
      cm(2)= 120000000.0d0
      cm(3)= 150000000.0d0
      cm(4)= 4000.0d0
      cm(5)= 4000.0d0
      cm(6)= 75000000.0d0
      cm(7)= 45000000.0d0
      cm(8)= 13000.0d0
      cm(9 )= -500000.0d0
      cm(10)= -350000.0d0
      cm(11)= -260000.0d0
      cm(12)= -300000.0d0
      cm(13)= -140000.0d0
      cm(14)= -100000.0d0
      cm(15)= 0.3d0
      cm(16)= 0.1d0
      cm(17)= 0.8d0
      cm(18)= 500000000.0d0
      cm(19)= 20.0d0
      cm(20)= 0.2d0
      cm(21)= 0.33d0
      cm(22)= 0.000088d0
      cm(23)= 0.33d0
      
      pi = 4.d0*datan(1.d0)
      alpha = pi/180*90
      num_hv = int(cm(19))! Number of history variables in use

      hsv(12) = 1.0d0
      hsv(16) = 1.0d0
      hsv(20) = 1.0d0



c      print "(2f12.2)", cm(1), cm(10)
      

      open (10, file='output_file.txt', status='unknown')


      do while(u_x .gt. -0.010d0)
        hsv(num_hv+1) = 1.0+u_x
        hsv(num_hv+2) = u_y
        hsv(num_hv+5) = 1.0d0
        hsv(num_hv+9) = 1.0d0
c
        call umat43(cm, eps, sig, epsp, hsv, dt1, etype,
     1   failel, elsiz, idele, reject)
c     
        write(10,*)"", u_x, sig(1), sig(2), sig(3), sig(4), sig(5), 
     1   sig(6)

        u = u+0.000002
        u_x = -sin(alpha)*u
        u_y = cos(alpha)*u
        sigma11_list(i) = sig(1)
        u_list(i) = u_x

c        
        
        i = i+1

      end do

 
      
      close(10)


      END
