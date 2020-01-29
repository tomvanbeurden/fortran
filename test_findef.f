C----|--1---------2---------3---------4---------5---------6---------7-|
C     Tom van Beurden (27-01-2019)
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
      double precision :: alpha_deg = 45.5d0 !Loading angle in degrees
      double precision :: u = 0.0d0, u_x = 0.0d0, u_y = 0.0d0 !define displacement components

 
c     Define parameters:
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
      cm(19)= 21.0d0
      cm(20)= 0.2d0
      cm(21)= 0.33d0
      cm(22)= 0.000088d0
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
      do while(u_x .gt. -0.83d0)
        hsv(num_hv+1) = 1.0d0+u_x
        hsv(num_hv+2) = u_y
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
