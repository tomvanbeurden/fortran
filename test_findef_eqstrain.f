C----|--1---------2---------3---------4---------5---------6---------7-|
c *******************************************************************
c | Edit by Tom van Beurden, 29−01−2020                       |
c | Eindhoven University of Technology & BMW Group            |
c *******************************************************************    
C     Test file to run umat file externally.


      PROGRAM TEST

      IMPLICIT NONE

      double precision :: cm(23) = 0.0d0 !Initiate material parameters
c      double precision :: eps(6) = 0.0d0 !Initiate strain 
      double precision :: sig(6) = 0.0d0 !Initiate stress
      double precision :: epsp(3) = 0.0d0 !Initiate plastic strain
      double precision :: hsv(30) = 0.0d0 !Initiate history array
      double precision :: dt1 = 0.1d0     !Define timestep
      integer :: nnpcrv(17) = -1  ! ??-> LSDyna parameter, req. input, unused
      character*5 :: etype = "solid" !Element type
      logical :: failel = .false.  ! ??-> LSDyna parameter, req. input, unused
      logical :: reject = .false.  ! ??-> LSDyna parameter, req. input, unused
      integer :: idele = -1        ! ??-> LSDyna parameter, req. input, unused
      integer :: num_hv       !number of history variables in use
      integer :: elsiz = -1   !unused parameter
      double precision :: pi
      double precision :: alpha	     !Loading angle
      double precision :: alpha_deg = 90.0d0 !Loading angle in degrees
      double precision :: u = 0.0d0, u_x = 0.0d0, u_y = 0.0d0 !define displacement components

 
c     Define parameters:
      cm(1)= 110000000.0d0  ! bulk
      cm(2)= 120000000.0d0  ! shear
      cm(3)= 150.0d6        ! Ett
      cm(4)= 0.040d6!Ell
      cm(5)= 0.040d6!Eww
      cm(6)= 75.0d6!Gtl
      cm(7)= 0.013d6!Glw
      cm(8)= 65.0d6!Gtw
      cm(9 )= -0.50d6
      cm(10)= -0.35d6
      cm(11)= -0.26d6
      cm(12)= -0.30d6
      cm(13)= -0.14d6
      cm(14)= -0.10d6
      cm(15)= 0.3d0
      cm(16)= 0.1d0
      cm(17)= 0.8d0
      cm(18)= 500.0d6
      cm(19)= 21.0d0
      cm(20)= 0.2d0
      cm(21)= 0.33d0
      cm(22)= 8.8d-5
      cm(23)= 0.33d0
      

      pi = 4.d0*datan(1.d0)
      alpha = pi/180.0d0*alpha_deg
      num_hv = int(cm(19))! Number of history variables in use

c     For my material model I use the deformation gradient of the prev. timestep. For this reason the old deformationg gradient is stored in the history array, and the diagonal components need to be initiated to 1.0.
      hsv(4) = 1.0d0
      hsv(8) = 1.0d0
      hsv(12) = 1.0d0

c     I also keep track of the plastic deformation gradient, initiating to unity:
      hsv(13) = 1.0d0
      hsv(17) = 1.0d0
      hsv(21) = 1.0d0

c     Open file to write output
      open (10, file='output_file.txt', status='unknown')

c     Loop until the displacement in x direction reaches a certain value. In the loop the new deformation gradient is prescribed. The new deformation gradient is in this case stored in the history array from hsv(num_hv+1) to hsv(num_hv+9). I apply compression under an angle wrt the 11 direction, so the F11 component and F21 component change. 
      do while(u_x .gt. -0.85d0)
        hsv(num_hv+1) = 1.0d0+u_x
        hsv(num_hv+3) = u_y
        hsv(num_hv+5) = 1.0d0
        hsv(num_hv+9) = 1.0d0
c
c       call material model:
        call umat43(cm, sig, epsp, hsv, dt1, etype,
     1   failel, elsiz, idele, reject)
c     
c       write output
        write(10,*)"", u_x, sig(1), sig(2), sig(3), sig(4), sig(5), 
     1   sig(6)

c       Displacement increment
        u = u+0.0002d0
        u_x = -sin(alpha)*u
        u_y = cos(alpha)*u

c        


        

      end do

 
      
      close(10)


      END
