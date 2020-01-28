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
      double precision norm_1, norm_2, norm_3

      complex*16 :: compmat(3,3), compvec(3)
      complex*16 :: eigvec_inv(3,3), multi(3,3)
      complex*16 :: eigval_diag_exp(3,3)
      logical :: flag = .false.
      wr = 0.0d0
      wi = 0.0d0
      zr = 0.0d0
      zi = 0.0d0

      
      Fp_c = 0.0d0
c      
c      Fp_r(1,1) = 1.0d0
c      Fp_r(2,1) = 5.0d0
c      Fp_r(3,1) = 8.0d0
c      Fp_r(1,2) = -2.0d0
c      Fp_r(2,2) = 7.50d0
c      Fp_r(3,2) = -2.5d0
c      Fp_r(1,3) = 1.6d0
c      Fp_r(2,3) = -4.8d0
c      Fp_r(3,3) = 3.4d0

      Fp_r = 0.0d0
      Fp_r(1,1) = 0.89d0
      Fp_r(2,1) = 0.1d0
      Fp_r(1,2) = -0.1d0
      Fp_r(2,2) = 1.0d0
      Fp_r(3,2) = 5.0d0
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
      
      
c      print*,"compmat", compmat(1,1)
      compmat(1,1) = cmplx(zr(1,1),zi(1,1))
c      print*,"compmat", compmat(1,1)
      compmat = cmplx(zr,zi)
      compvec = cmplx(wr,wi)


      print*,"eigvec", compvec(1),compvec(2),compvec(3)
      print*,"exp" , exp(compvec(1))

      print*,"eigvec_real", real(compvec)
      print*,"eigcomp", compmat(1,1), compmat(1,2), compmat(1,3)
      print*,"eigcomp", compmat(2,1), compmat(2,2), compmat(2,3)
      print*,"eigcomp", compmat(3,1), compmat(3,2), compmat(3,3)
c    
      norm_1 = norm2([norm2(real(compmat(:,1))),
     1  norm2(aimag(compmat(:,1)))])
      norm_2 = norm2([norm2(real(compmat(:,2))),
     1  norm2(aimag(compmat(:,2)))])
      norm_3 = norm2([norm2(real(compmat(:,3))),
     1  norm2(aimag(compmat(:,3)))])     


c      compmat(:,1) = compmat(:,1)/norm_1
c      compmat(:,2) = compmat(:,2)/norm_2
c      compmat(:,3) = compmat(:,3)/norm_3

      print*, "norm", norm_1, norm_2, norm_3
      call M33INV_comp (compmat, eigvec_inv,  flag)

      print*,"eigveci", eigvec_inv(1,1),eigvec_inv(1,2),eigvec_inv(1,3)
      print*,"eigveci", eigvec_inv(2,1),eigvec_inv(2,2),eigvec_inv(2,3)
      print*,"eigveci", eigvec_inv(3,1),eigvec_inv(3,2),eigvec_inv(3,3)

c      print*,"eigcomp", compmat(1,1), compmat(1,2), compmat(1,3)
c      print*,"eigcomp", compmat(2,1), compmat(2,2), compmat(2,3)
c      print*,"eigcomp", compmat(3,1), compmat(3,2), compmat(3,3)

      eigval_diag_exp = 0.0d0
      eigval_diag_exp(1,1) = exp(compvec(1))
      eigval_diag_exp(2,2) = exp(compvec(2))
      eigval_diag_exp(3,3) = exp(compvec(3))


      
c      print*,"eigval_diag_exp", eigval_diag_exp(1,1), 
c     1 eigval_diag_exp(1,2), eigval_diag_exp(1,3)
c      print*,"eigval_diag_exp", eigval_diag_exp(2,1), 
c     1 eigval_diag_exp(2,2), eigval_diag_exp(2,3)
c      print*,"eigval_diag_exp", eigval_diag_exp(3,1), 
c     1 eigval_diag_exp(3,2), eigval_diag_exp(3,3)
      multi = matmul(compmat,matmul(eigval_diag_exp,eigvec_inv))

      print*,"multi", multi(1,1),multi(1,2),multi(1,3)
      print*,"multi", multi(2,1),multi(2,2),multi(2,3)
      print*,"multi", multi(3,1),multi(3,2),multi(3,3)

      END













