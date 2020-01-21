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
      real :: cm(23) = 0.0d0
      real :: eps(6) = 0.0d0
      real :: sig(6) = 0.0d0
      real :: epsp(3) = 0.0d0
      real :: hsv(20) = 0.0d0
      real :: dt1 = 0.1d0
      integer :: nnpcrv(17) = -1
      integer :: imax = 0
      character :: etype = "solid"
      logical :: failel = .false.
      logical :: reject = .false.
      integer :: idele = -1
      integer :: num_hv = 11
      integer :: elsiz = -1;

	  print *, "Hello from main" 

      cm(1)= 110000000d0
      cm(2)= 120000000d0
      cm(3)= 130000000d0
      cm(4)= 140000000d0
      cm(5)= 150000000d0
      cm(6)= 160000000d0
      cm(7)= 170000000d0
      cm(8)= 180000000d0
      cm(9 )= 0.5000000d0
      cm(10)= 0.4000000d0
      cm(11)= 0.35000000d0
      cm(12)= 0.3000000d0
      cm(13)= 0.25000000d0
      cm(14)= 0.20000000d0
      cm(15)= 0.3d0
      cm(16)= 0.1d0
      cm(17)= 0.8d0
      cm(18)= 500.0000000d0
      cm(19)= 11.0d0
      cm(20)= 0.2d0
      cm(21)= 0.33d0
      cm(22)= 0.000088d0
      cm(23)= 0.33d0

      

c      print "(2f12.2)", cm(1), cm(10)


      call umat43(cm, eps, sig, epsp, hsv, dt1, etype,
     1 failel, elsiz, idele, reject)





      END
