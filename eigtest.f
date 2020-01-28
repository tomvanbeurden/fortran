C----|--1---------2---------3---------4---------5---------6---------7-|
C     Tom van Beurden (28-01-2019)
C     Test file to run eigval solver
C
C     To compile  : gfortran eigtest.f eigcomplex.f
C     To run      : ./a.out




     
      PROGRAM EIGTEST


      IMPLICIT NONE
      double precision :: Fp_r(3,3)
      double precision :: Fp_c(3,3)
      integer :: nm = 3, n = 3, matz = 1
      double precision:: wr(3), zr(3,3), zr1(3), zr2(3),zr3(3)
      double precision:: wi(3), zi(3,3)
      integer ierr
      double precision fv1(3), fv2(3), fv3(3)
      double precision zr1_norm,zr2_norm,zr3_norm
      wr = 0.0d0
      wi = 0.0d0
      zr = 0.0d0
      zi = 0.0d0

      
      Fp_c = 0.0d0
c      
      Fp_r(1,1) = 0.89d0
      Fp_r(2,1) = 0.1d0
      Fp_r(3,1) = 0.0d0
      Fp_r(1,2) = -0.1d0
      Fp_r(2,2) = 1.0d0
      Fp_r(3,2) = 0.0d0
      Fp_r(1,3) = 0.0d0
      Fp_r(2,3) = 0.0d0
      Fp_r(3,3) = 1.0d0
      
      print*,"Fp1", Fp_r(1,1), Fp_r(1,2), Fp_r(1,3)
      print*,"Fp2", Fp_r(2,1), Fp_r(2,2), Fp_r(2,3)
      print*,"Fp3", Fp_r(3,1), Fp_r(3,2), Fp_r(3,3)

c      print*,"Fp1_c", Fp_c(1,1), Fp_c(1,2), Fp_c(1,3)
c      print*,"Fp2_c", Fp_c(2,1), Fp_c(2,2), Fp_c(2,3)
c      print*,"Fp3_c", Fp_c(3,1), Fp_c(3,2), Fp_c(3,3)
      
      call cg(nm,n,Fp_r,Fp_c,wr,wi,matz,zr,zi,fv1,fv2,fv3,ierr)

      print*,"eigval", wr(1), wr(2), wr(3)
      print*,"eigval_i", wi(1), wi(2), wi(3)
      
 
      print*,"eigvec1_r", zr(1,1), zr(1,2), zr(1,3)
      print*,"eigvec2_r", zr(2,1), zr(2,2), zr(2,3)
      print*,"eigvec3_r", zr(3,1), zr(3,2), zr(3,3)

      print*,"eigvec1_c", zi(1,1), zi(1,2), zi(1,3)
      print*,"eigvec2_c", zi(2,1), zi(2,2), zi(2,3)
      print*,"eigvec3_c", zi(3,1), zi(3,2), zi(3,3)

      END
