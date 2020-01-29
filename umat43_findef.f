C----|--1---------2---------3---------4---------5---------6---------7-|
C     Tom van Beurden (27-01-2019)

c_________________________________________________________________________  
c 
c                        MAIN CODE
c_________________________________________________________________________


      subroutine umat43(cm, sig, epsp, hsv, dt1, etype,
     1 failel, elsiz, idele, reject)
c
c*********************************************************************
c| Livermore Software Technology Corporation (LSTC)             |
c| −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−− |
c| Copyright 1987−2008 Livermore Software Tech. Corp            |
c| Alll rights reserved                                         |
c| −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−− |
c| Edit by Dennis van Iersel, 19−09−2018                        |
c| Eindhoven University of Technology & BMW Group               |
c| Based on implementation by Popp 2007                         |
c*********************************************************************
c
c     Phenomenological metallic honeycomb material model
c
c     Variables
c
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
c         .
c         .
c     cm(lmc)= last material constant
c
c
c     Note :
c     − eps is in strain Voigt notation, meaning that
c     e.g. eps(4) is the sum of tensorial xy and yx strain increments.
c     − eps is computed as true (logarithmic) strain
c
c     sig(1)= local x stress
c     sig(2)= local y stress
c     sig(3)= local z stress
c     sig(4)= local xy stress
c     sig(5)= local yz stress
c     sig(6)= local zx stress
c     hsv(1)= local plastic strain TT−axis
c     hsv(2)= non−local plastic strain TT−axis
c     hsv(3)= yield law value
c     hsv(4)=F11  (old deformation gradient)
c     hsv(5)=F21  (old deformation gradient)
c     hsv(6)=F31  (old deformation gradient)
c     hsv(7)=F12  (old deformation gradient)
c     hsv(8)=F22  (old deformation gradient)
c     hsv(9)=F32  (old deformation gradient)
c     hsv(10)=F13 (old deformation gradient)
c     hsv(11)=F23 (old deformation gradient)
c     hsv(12)=F33 (old deformation gradient)
c     hsv(13)=Fp11 (plastic deformation gradient)
c     hsv(14)=Fp21 (plastic deformation gradient)
c     hsv(15)=Fp31 (plastic deformation gradient)
c     hsv(16)=Fp12 (plastic deformation gradient)
c     hsv(17)=Fp22 (plastic deformation gradient)
c     hsv(18)=Fp32 (plastic deformation gradient)
c     hsv(19)=Fp13 (plastic deformation gradient)
c     hsv(20)=Fp23 (plastic deformation gradient)
c     hsv(21)=Fp33 (plastic deformation gradient)
c         .
c     hsv(nhv)=nhvth history variable
c     hsv(nhv+1)=F11 (deformation gradient)
c     hsv(nhv+2)=F21 (deformation gradient)
c     hsv(nhv+3)=F31 (deformation gradient)
c     hsv(nhv+4)=F12 (deformation gradient)
c     hsv(nhv+5)=F22 (deformation gradient)
c     hsv(nhv+6)=F32 (deformation gradient)
c     hsv(nhv+7)=F13 (deformation gradient)
c     hsv(nhv+8)=F23 (deformation gradient)
c     hsv(nhv+9)=F33 (deformation gradient)
c
c     dt1=current time step size
c     capa=reduction factor for transverse shear
c     etype:
c        eq. ”solid” for solid elements
c        eq. ”sph” for smoothed particle hydronamics
c        eq. ”sld2d” for shell forms 13 (2D solids - plane strain)
c        eq. ”sldax” for shell forms 14  and 15 (2D solids)
c        eq. ”shlt” for shell forms 25, 26, and 27 (shell with thickness stretch)
c        eq. ”shell” for all other shell elements plus thick shell forms 1 and 2
c        eq. ”tshel” for thcik shell forms forms 3 and 5
c        eq. ”hbeam” for beam element forms 1, 11, 14
c        eq. ”tbeam ” for beam element form 3 (truss)
c        eq. ”dbeam” for beam element form 6 (discrete)
c        eq. ”beam ” for all other beam elements (currently not used)
c
c     time=current problem time.
c     temp=current temperature
c
c     cma=additional memory for material data defined by LMCA at 
c       6th field of 2nd crad of *DATA_USER_DEFINED
c
c     All transformations into the element local system are 
c     performed prior to entering this subroutine. Transformations
c     back tot he global system are perforemd after exiting this
c     routine.
c
c     All history variables are initialized to zero in the input
c     phase. Initialization of history variables to nonzeroe values
c     may be done during the first call to this subroutine for each
c     element.
c
c     Energy calculations for the dyna3d energy balance are done
c     outside of this subroutine
c
c
c
c
c
      double precision cm(23), hsv(30)
      integer nnpcrv
      double precision dt1
      integer imax
      character*5 etype
      logical failel, reject
      integer*8 idele, num_hv, elsiz

      double precision tol, m, dlambda
      double precision Ett, Ell, Eww, Gtl, Glw, Gtw, ytt, ytl, ytw
      double precision nutl, nult, nutw, nuwt, nulw, nuwl, Delta
      double precision yld0tt, yld0tl, yld0tw, yldcrtt, yldcrtl,yldcrtw
      double precision flow11, flow12, flow13
      double precision epsp_tt, ftrial, epsd, hd, g
      double precision sig_trial(6), sig(6)
      double precision sig_old(6),var_nonloc,var_loc,dnl,lr,le,eps_star
      double precision C4(6,6)
      double precision F_old(3,3), F_new(3,3), F_dot(3,3), F_mid(3,3)
      double precision F_midinv(3,3), L_mid(3,3), D_mid(3,3),W_mid(3,3)
      double precision Fp_old(3,3), Fp_new(3,3), Dp(3,3), deltaLp(3,3)
      double precision Dp_v(6),zeros(3,3)

c     for eigproblem Fp
      double precision wr(3), wi(3), zr(3,3),zi(3,3)
      double precision fv1(3),fv2(3),fv3(3)
      integer ierr
      complex*16 :: eigval(3), eigval_diag_exp(3,3), eigvec(3,3)
      complex*16 :: eigvec_inv(3,3)
      double precision dFp(3,3), dFp_i(3,3), Fp_warning

      logical :: SING_flag_M33inv = .false.


c
c     Initialize/set parameters:−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
      num_hv = int(cm(19))! Number of history variables in use
      toll = 0.000001! Tolerance for yieldsurface
      imax = 5000! Max. number of iterations for plastic multiplier
      m = 1.0d0 ! Parameter in yield surface
      dlambda= 0.0d0 ! Plastic multiplier
      flow11 = 0.0d0
      flow12 = 0.0d0
      flow13 = 0.0d0
      sig_trial = 0.0d0

c
c      if(ncycle.eq.1) then
c        call usermsg(’mat43’)
c      end if
c
c     Assign material parameters/read history variables:−−−−−−−
c
c     Moduli of elasticity:
      Ett=cm(3)
      Ell=cm(4)
      Eww=cm(5)
      Gtl=cm(6)
      Glw=cm(7)
      Gtw=cm(8)
c
c     Poisson’s ratios:
      nutl= cm(21)
      nulw= cm(22)
      nutw= cm(23)
      nult = nutl*Ell/Ett
      nuwl = nulw*Eww/Ell
      nuwt = nutw*Eww/Ett
c
c     Yield stresses:
      yld0tt=cm(9) !Initial yield stress
      yld0tl=cm(10)
      yld0tw=cm(11)
      yldcrtt=cm(12) !Plateau stress
      yldcrtl=cm(13)
      yldcrtw=cm(14)
      ytt=yld0tt
      ytl=yld0tl
      ytw=yld0tw
c
c     Yield stress function shape
      hd=cm(18) ! Hardening modulus during densification 
      g=cm(20) ! Softening shape parameter
      epsd = cm(17) ! Densification strain
c
c     Stored plastic strain:
      epsp_tt=hsv(1)

c     Stress of previous timestep:
      sig_old = sig

c     Deformation gradient of previous timestep
      F_old(1,1)= hsv(4) ! F11 (old deformation gradient)
      F_old(2,1)= hsv(5) ! F21 (old deformation gradient)
      F_old(3,1)= hsv(6) ! F31 (old deformation gradient)
      F_old(1,2)= hsv(7) ! F12 (old deformation gradient)
      F_old(2,2)= hsv(8) ! F22 (old deformation gradient)
      F_old(3,2)= hsv(9) ! F32 (old deformation gradient)
      F_old(1,3)= hsv(10)! F13 (old deformation gradient)
      F_old(2,3)= hsv(11)! F23 (old deformation gradient)
      F_old(3,3)= hsv(12)! F33 (old deformation gradient)

c     Plastic deformation gradient of previous timestep
      Fp_old(1,1)= hsv(13)! Fp11 (old plastic deformation gradient)
      Fp_old(2,1)= hsv(14)! Fp21 (old plastic deformation gradient)
      Fp_old(3,1)= hsv(15)! Fp31 (old plastic deformation gradient)
      Fp_old(1,2)= hsv(16)! Fp12 (old plastic deformation gradient)
      Fp_old(2,2)= hsv(17)! Fp22 (old plastic deformation gradient)
      Fp_old(3,2)= hsv(18)! Fp32 (old plastic deformation gradient)
      Fp_old(1,3)= hsv(19)! Fp13 (old plastic deformation gradient)
      Fp_old(2,3)= hsv(20)! Fp23 (old plastic deformation gradient)
      Fp_old(3,3)= hsv(21)! Fp33 (old plastic deformation gradient)

c     Deformation gradient of current timestep
      F_new(1,1)= hsv(num_hv+1)! F11 (old deformation gradient)
      F_new(2,1)= hsv(num_hv+2)! F21 (old deformation gradient)
      F_new(3,1)= hsv(num_hv+3)! F31 (old deformation gradient)
      F_new(1,2)= hsv(num_hv+4)! F12 (old deformation gradient)
      F_new(2,2)= hsv(num_hv+5)! F22 (old deformation gradient)
      F_new(3,2)= hsv(num_hv+6)! F32 (old deformation gradient)
      F_new(1,3)= hsv(num_hv+7)! F13 (old deformation gradient)
      F_new(2,3)= hsv(num_hv+8)! F23 (old deformation gradient)
      F_new(3,3)= hsv(num_hv+9)! F33 (old deformation gradient)
      

      F_mid = (F_new + F_old)*0.5d0 !deformation gradient at time t+1/2dt
      F_dot = (F_new - F_old)*1.0d0/dt1 !time derrivative of deformation gradient
      
c      print*,"hsv", hsv(1),hsv(2),hsv(3),hsv(4),hsv(5),hsv(6),hsv(7)
c      print*,"Fmid",F_mid(1,1),F_mid(2,2),F_mid(3,3),F_mid(1,2)
c      print*,"Fold",F_old(1,1),F_old(2,2),F_old(3,3),F_old(1,2)
c      print*,"Fnew",F_new(1,1),F_new(2,2),F_new(3,3),F_new(1,2)
      call M33inv (F_mid, F_midinv, SING_flag_M33inv) 

      if ( SING_flag_M33inv) then
       print*, "ERROR_M33inv:_Singular_matrix!"
      endif

      L_mid = matmul(F_dot,F_midinv) !velocity gradient at time t+1/2dt
      D_mid = 0.5d0*(L_mid+transpose(L_mid)) !Rate of deformation at time t+1/2dt
      W_mid = 0.5d0*(L_mid-transpose(L_mid)) !Spin tensor at time t+1/2dt

c
c     Compute material properties: −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
c
c     Compute stiffnesstensor orthotropic material:
      Delta = 1.0d0-nutl*nult-nulw*nuwl-nuwt*nutw-2.0d0*nuwl*nult*nutw
      C4= 0.0d0
      C4(1,1)= (1.0d0-nulw*nuwl)/Delta*Ett
      C4(2,2)= (1.0d0-nutw*nuwt)/Delta*Ell
      C4(3,3)= (1.0d0-nutl*nult)/Delta*Eww
      C4(1,2)= (nutl+nutw*nuwl)/Delta*Ell
      C4(2,3)= (nulw+nult*nutw)/Delta*Eww
      C4(1,3)= (nutw+nutl*nulw)/Delta*Eww
      C4(2,1)= C4(1,2)
      C4(3,2)= C4(2,3)
      C4(3,1)= C4(1,3)
      C4(4,4)= 2.0d0*Gtl
      C4(5,5)= 2.0d0*Glw
      C4(6,6)= 2.0d0*Gtw
c
c     Non−local correction initial yield stress
      lr = cm(15)
      le = cm(16)
      eps_star = epsd*(.5d0-(.5 d0*le-(le**3.0d0)/(12.0d0*lr**2.0d0)+
     1   (le**5.0d0)/(160.0d0*lr**4.0d0))/(16.0d0*lr/15.0d0))
      dnl = min(var_nonloc/(eps_star), 0.99d0)
      yld0tt = yld0tt-dnl*(yld0tt-yldcrtt)
      yld0tl = yld0tl-dnl*(yld0tl-yldcrtl)
      yld0tw = yld0tw-dnl*(yld0tw-yldcrtw)
c
      if (yld0tt .ge. 0.0d0)
     1 print *, "WARNING:_compressive_strength_should_be_negative!"

c
c     Compute linear strain tensor: −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
c
c     Check elementtype: 
      if (etype .eq. 'solid') then
        call get_stressincrement(sig_trial, sig_old, F_old, F_new, 
     1  F_mid, D_mid, C4, dlambda, flow11, flow12, flow13, dt1, Dp_v)
c
c        Check for compressive stresses
        if(sig(1) .le. 0.0d0) then
c
c         Compute yield stresses
          call compute_yield_stress(epsp_tt, ytt, yld0tt, yldcrtt, 
     1     epsd, hd, g)
          call compute_yield_stress(epsp_tt, ytl, yld0tl, yldcrtl, 
     2     epsd, hd, g)
          call compute_yield_stress(epsp_tt, ytw, yld0tw, yldcrtw, 
     3     epsd, hd, g)
c
c
c         Check if yieldsurface boundary is exceeded
          call compute_yield_coupled(ftrial,sig_trial(1),
     1     sig_trial(4), sig_trial(6), m, ytt, ytl, ytw)

          if (ftrial .gt. tol) then
c
c           Compute plastic flow direction and step size
            call get_flowvector(sig_old,flow11,flow12,flow13)
            call get_plasticmultiplier(dlambda, sig_old, 
     1        F_old,F_new, F_mid, D_mid, C4, flow11, flow12, flow13,   
     2        dt1, ytt, ytl, ytw, m)
c
c
            call get_stressincrement(sig, sig_old, F_old, F_new, F_mid,
     1       D_mid, C4, dlambda, flow11, flow12, flow13, dt1, Dp_v)

            call voigt2full(Dp,Dp_v)
c	    Dp = 0.0d0
c	    Dp(1,1) = dlambda*flow11
c	    Dp(1,2) = 0.5d0*dlambda*flow12
c	    Dp(1,3) = 0.5d0*dlambda*flow13
c	    Dp(2,1) = Dp(1,2)
c           Dp(3,1) = Dp(1,3)
            deltaLp = dt1*(Dp+W_mid)
            zeros = 0.0d0


c           Complex eigensolver (open source tool EISPACK)
c           call cg(nm,n,Fp_r,Fp_c,wr,wi,matz,zr,zi,fv1,fv2,fv3,ierr)
            call cg(3,3,deltaLp,zeros,wr, wi, 1,zr,zi,
     1       fv1,fv2,fv3,ierr)

            eigvec = cmplx(zr,zi)
            eigval = cmplx(wr,wi)
            eigval_diag_exp = 0.0d0
            eigval_diag_exp(1,1) = exp(eigval(1))
            eigval_diag_exp(2,2) = exp(eigval(2))
            eigval_diag_exp(3,3) = exp(eigval(3))

            call M33INV_comp (eigvec, eigvec_inv,  SING_flag_M33inv)

            dFp_i = aimag(matmul(eigvec,
     1        matmul(eigval_diag_exp,eigvec_inv)))
            
            Fp_warning=abs(dFp_i(1,1))+abs(dFp_i(1,2))+abs(dFp_i(1,3))
     1       +abs(dFp_i(2,1))+abs(dFp_i(2,2))+abs(dFp_i(2,3))
     1       +abs(dFp_i(3,1))+abs(dFp_i(3,2))+abs(dFp_i(3,3))

            if (Fp_warning .gt. 0.000001)
     1       print*,"warning: complex Fp", Fp_warning


            dFp = real(matmul(eigvec,
     1        matmul(eigval_diag_exp,eigvec_inv)))
            Fp_new = matmul(dFp,Fp_old)
            
          else
c           Update stresses for elastic compression:
            sig = sig_trial
            Fp_new = Fp_old
          endif
        else
c         Update stresses for tensile case (linearelastic):
          sig = sig_trial
          Fp_new = Fp_old
        endif
c
c       Store relevant paremeters: −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
c
c       Store plastic strains in history variables:
c        hsv(1)= epsp_tt
        hsv(1)= 1-F_new(1,1)

c
c       Store non-local variables
        hsv(2)= var_nonloc
        var_loc=epsp_tt
        hsv(3)= ftrial
c
c
c       Store old deformation gradient
        hsv(4)= hsv(num_hv+1) !F11_old (deformation gradient)
        hsv(5)= hsv(num_hv+2) !F21_old (deformation gradient)
        hsv(6)= hsv(num_hv+3) !F31_old (deformation gradient)
        hsv(7)= hsv(num_hv+4) !F12_old (deformation gradient)
        hsv(8)= hsv(num_hv+5) !F22_old (deformation gradient)
        hsv(9)= hsv(num_hv+6) !F32_old (deformation gradient)
        hsv(10)= hsv(num_hv+7) !F13_old (deformation gradient)
        hsv(11)= hsv(num_hv+8) !F23_old (deformation gradient)
        hsv(12)= hsv(num_hv+9) !F33_old (deformation gradient)

c       Store plastic deformation gradient
        hsv(13)= Fp_new(1,1) !Fp11_old (new pl deformation gradient)
        hsv(14)= Fp_new(2,1) !Fp21_old (new pl deformation gradient)
        hsv(15)= Fp_new(3,1) !Fp31_old (new pl deformation gradient)
        hsv(16)= Fp_new(1,2) !Fp12_old (new pl deformation gradient)
        hsv(17)= Fp_new(2,2) !Fp22_old (new pl deformation gradient)
        hsv(18)= Fp_new(3,2) !Fp32_old (new pl deformation gradient)
        hsv(19)= Fp_new(1,3) !Fp13_old (new pl deformation gradient)
        hsv(20)= Fp_new(2,3) !Fp23_old (new pl deformation gradient)
        hsv(21)= Fp_new(3,3) !Fp33_old (new pl deformation gradient)
c
        
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c EDIT: commented out this part below:
c     Material model only available for solids:
c      else
c        cerdat(1)= etype
c        call lsmg(3,MSG_SOL+1151,ioall, ierdat,rerdat,cerdat,0)
c      endif
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      endif
c
c
      return
      end


c_________________________________________________________________________
c 
c                      STRESS INCREMENT
c_________________________________________________________________________

      subroutine get_stressincrement(sig_new_v, sig_old_v, F_old,F_new,
     1  F_mid, D, C4,dlambda, flow11, flow12, flow13, dt1, Dp_v)
      
      double precision sig_new_v(6), sig_old_v(6)
      double precision dlambda, flow11, flow12, flow13
      double precision F_old(3,3), F_new(3,3), F_mid(3,3), D(3,3)
      double precision C4(6,6), dt1

      double precision sig_old(3,3), sig_old_bar(3,3)
      double precision F_oldinv(3,3), F_midinv(3,3)
      double precision Deinv(3,3), Dp(3,3), De(3,3),De_bar(3,3),Dp_v(6) 
      double precision De_bar_v(6), CdoubleD_v(6), CdoubleD(3,3)
      double precision sig_new_bar(3,3), sig_new(3,3)
      logical :: SING_flag_M33inv = .false.

      call voigt2full(sig_old,sig_old_v)
c      sig_old(1,1) = sig_old_v(1)
c      sig_old(2,2) = sig_old_v(2)
c      sig_old(3,3) = sig_old_v(3)
c      sig_old(1,2) = sig_old_v(4)
c      sig_old(2,3) = sig_old_v(5)
c      sig_old(1,3) = sig_old_v(6)
c      sig_old(2,1) = sig_old(1,2)
c      sig_old(3,1) = sig_old(1,3)
c      sig_old(3,2) = sig_old(2,3)
c
      call M33inv(F_old, F_oldinv, SING_flag_M33inv)
c      if ( SING_flag_M33inv) then
c       print*, "ERROR_M33inv:_Singular_matrix!"
c      endif
c
      call M33inv(F_mid, F_midinv, SING_flag_M33inv)
c      if ( SING_flag_M33inv) then
c       print*, "ERROR_M33inv:_Singular_matrix!"
c      endif
c
      sig_old_bar = matmul(matmul(F_oldinv,sig_old),
     1 transpose(F_oldinv)) 
c     
      Dp_v = 0.0d0
      Dp_v(1) = dlambda*flow11
      Dp_v(4) = 0.5d0*dlambda*flow12
      Dp_v(6) = 0.5d0*dlambda*flow13
      
      call voigt2full(Dp, Dp_v)
c
      De = D - Dp 
c
      call M33inv(De, Deinv, SING_flag_M33inv)
      
c      if ( SING_flag_M33inv) then
c       print*, "ERROR_M33inv:_Singular_matrix!"
c      endif
      De_bar = matmul(matmul(F_midinv,De),
     1 transpose(F_midinv))
c
      call full2voigt(De_bar,De_bar_v)
c      De_bar_v(1) = De_bar(1,1)
c      De_bar_v(2) = De_bar(2,2)
c      De_bar_v(3) = De_bar(3,3)
c      De_bar_v(4) = De_bar(1,2)
c      De_bar_v(5) = De_bar(2,3)
c      De_bar_v(6) = De_bar(1,3)
c
      CdoubleD_v = matmul(C4,De_bar_v)
      call voigt2full(CdoubleD,CdoubleD_v)
c      CdoubleD(1,1) = CdoubleD_v(1)
c      CdoubleD(2,2) = CdoubleD_v(2)
c      CdoubleD(3,3) = CdoubleD_v(3)
c      CdoubleD(1,2) = CdoubleD_v(4)
c      CdoubleD(2,3) = CdoubleD_v(5)
c      CdoubleD(1,3) = CdoubleD_v(6)
c      CdoubleD(2,1) = CdoubleD(1,2)
c      CdoubleD(3,1) = CdoubleD(3,1)
c      CdoubleD(3,2) = CdoubleD(2,3)
c
      sig_new_bar = sig_old_bar + dt1*CdoubleD
      sig_new = matmul(matmul(F_new,sig_new_bar),transpose(F_new))
c
      call full2voigt(sig_new,sig_new_v)
c      sig_new_v(1) = sig_new(1,1)
c      sig_new_v(2) = sig_new(2,2)
c      sig_new_v(3) = sig_new(3,3)
c      sig_new_v(4) = sig_new(1,2)
c      sig_new_v(5) = sig_new(2,3)
c      sig_new_v(6) = sig_new(1,3)
c
      end



c_________________________________________________________________________
c 
c                      YIELD STRESS COMUPTATION
c_________________________________________________________________________

      subroutine compute_yield_stress(epsp, yld, yld0, yldcr, 
     1 epsd, hd, g)
c *******************************************************************
c | Edit by Dennis van Iersel, 27−09−2018                     |
c | Eindhoven University of Technology & BMW Group            |
c | Based on paper by Mohr & Dojojo 2003                      |
c *******************************************************************
c
c     Routine computes the current yield stress according to
c     accumulated plastic strain
c
      implicit none
      double precision epsp, yld, yld0, yldcr, epsd, g, hd
      double precision s_avg, sf, sc, lambdap, lambdac

c
c     Softening/platuea stress region:
      if(epsp .ge. 0.0d0 .and. epsp .le. epsd) then
        s_avg =(1.0d0-g)*yldcr+g*yld0
        yld = yldcr+(yld0-yldcr)*((epsd-epsp)/(epsd))
     1   **((yld0-s_avg)/(s_avg-yldcr))
c     Densification region:
      elseif(epsp .gt. epsd) then
        yld = yldcr + dsign(1.0d0,yldcr)*hd*(epsp-epsd)**2.0d0
      endif

c!!!!!!!!!!!!!!!!!!!!!!!
c CHANGED dsign(yldcr,1.0d0) to sign(1.0,yldcr)
c Because of double error and "sign(A,B) returns value of A with sign of B"
c!!!!!!!!!!!!!!!!!!!!!!!      
      return
      end


c_________________________________________________________________________
c 
c                      YIELD LAW COMUPTATION
c_________________________________________________________________________

      subroutine compute_yield_coupled(f,sig11,sig12,sig13,m,
     1 ytt,ytl,ytw)

c ******************************************************************
c | Eit by Dennis van Iersel, 26-09-2018                      |
c | Eindhoven University of Technology & BMW Group            |
c | Based on paper by Mohr & Dojojo 2004                      |
c ******************************************************************
c

      implicit none
      double precision m, f , sig11, sig12, sig13, ytt, ytl, ytw
c

      f =(sig11/ytt)+((sig12/ytl)**2.0d0+(sig13/ytw)**2.0d0)**(m*0.5d0)
     1  -1.0d0
c
      return
      end




c_________________________________________________________________________
c 
c                       PLASTIC FLOW VECTOR
c_________________________________________________________________________

      subroutine get_flowvector(sig, flow11, flow12, flow13 )
c ******************************************************************
c | Edit by Dennis van Iersel, 26−09−2018                |
c | Eindhoven University of Technology & BMW Group       |
c | Based on implementation by Popp 2007                 |
c | Based on paper by Mohr & Dojojo 2004                 |
c ******************************************************************
      implicit none
      double precision sig(6)
      double precision flow11, flow12, flow13
      double precision S11, S12, S13, rJ1, rJ2, rJ3
      double precision solu, soluInv, eval1, eval2, eval3, evec1
      double precision evec2, evec3, evecNorm, evecNormInv, tol
c
      tol = 0.000001d0
c
c     Reduced stress tensor
      S11 = sig(1)
      S12 = sig(4)
      S13 = sig(6)
c     S22, S23 & S33 are zero because of the reduced tensor
c
C     Stress tensor invariants( simplified for zero components):
      rJ1=S11
      rJ2=-S12**2-S13**2
c     rJ3 =0.0 d0
c
c     Principal stresses (eigenvalues) of reduced stress tensor:
      eval1 =0.0d0
      eval2 =.5d0*(rJ1+sqrt(rJ1**2.0d0-4.d0*rJ2))
      eval3 =.5d0*(rJ1-sqrt(rJ1**2.0d0-4.d0*rJ2))
c
c     Identificaiton of principle compressive stress:
      solu =0.0d0
      if(eval1.lt.eval3) then
        if(eval1 .lt. eval2) then
          solu=eval1
        else
          solu=eval2
        endif
      else
        if(eval2 .lt. eval3) then
          solu=eval2
        else
          solu=eval3
        endif
      endif
c
      if(solu .eq. 0.0d0) solu = tol
        soluInv = 1.0d0/solu
c
c     Determination of corresponding eigenvector:
      evec1 = 1.0d0
      evec2 = S12*soluInv
      evec3 = S13*soluInv
      evecNorm = sqrt(evec1**2.0d0+evec2**2.0d0+evec3**2.0d0)
      evecNormInv = 1.0d0/evecNorm
      evec1 = evec1*evecNormInv
      evec2 = evec2*evecNormInv
      evec3 = evec3*evecNormInv
c
c     Compute plastic flowvector:
      flow11 = -evec1
      flow12 = -evec2
      flow13 = -evec3

c!!!!!!!!!!!!!!!!!!!!!!!
c EDIT1: sign(1.d0,evec1) to sign(1.0,evec1)
c EDIT2: deleted sign alltogether 
c!!!!!!!!!!!!!!!!!!!!!!!  
c
      return
      end




c_________________________________________________________________________
c 
c                       PLASTIC MULTIPLIER
c_________________________________________________________________________

      subroutine get_plasticmultiplier(dlambda, sig_old_v, 
     1  F_old,F_new, F_mid, D, C4, flow11, flow12, flow13, dt1, ytt, 
     2  ytl, ytw, m)

c ******************************************************************
c | Edit by Dennis van Iersel, 26−09−2018                     |
c | Eindhoven University of Technology & BMW Group            |
c | Based on implementation by Popp 2007                      |
c | Based on paper by Mohr & Dojojo 2004                      |
c ******************************************************************
      implicit none
      double precision dlambda, sig_old_v(6)
      double precision F_old(3,3), F_new(3,3), F_mid(3,3), D(3,3)
      double precision C4(6,6), Dp_v(6)
      double precision flow11, flow12, flow13, dt, ytt, ytl, ytw, dt1, m

      double precision sig_test(6)
      double precision xr, fr, xi, ximin, fi, fimin, tol
      integer iter, imax

c     Initiate parameters:
      fr = 100.0d0
      iter = 0
      xi = 0.00001d0
      ximin = 0.0d0
      imax = 5000
      tol = 0.00001d0
c
c     Secant method for determination of plastic multiplier:
      do while(abs(fr) .gt. tol .and. iter .lt. imax )
        iter = iter +1
c
        call get_stressincrement(sig_test, sig_old_v, F_old,F_new, 
     1  F_mid, D, C4, ximin, flow11, flow12, flow13, dt1, Dp_v) 
c
        call compute_yield_coupled(fimin, sig_test(1), sig_test(4), 
     1   sig_test(6),m, ytt, ytl, ytw)
c
        call get_stressincrement(sig_test, sig_old_v, F_old,F_new, 
     1  F_mid, D, C4, xi, flow11, flow12, flow13, dt1, Dp_v) 
c
        call compute_yield_coupled(fi, sig_test(1), sig_test(4), 
     1   sig_test(6),m, ytt, ytl, ytw)
c
c       Root update with secant method:
        xr = 0.0d0
        if(abs(fimin-fi) .gt. 0.0d0)
     1    xr=xi-fi*(ximin-xi)/(fimin-fi)
c
c
        if(abs(fimin-fi) .eq. 0.0d0)
     1    print*, "ERROR_ Use_double_precision!"
c
        call get_stressincrement(sig_test, sig_old_v, F_old,F_new, 
     1  F_mid, D, C4, xr, flow11, flow12, flow13, dt1, Dp_v) 

c
        call compute_yield_coupled(fr, sig_test(1), sig_test(4), 
     1   sig_test(6),m, ytt, ytl, ytw)
c
c       Update for next secant step:
        ximin=xi
        xi=xr
      end do
c
c     Error messages :
      if(iter .eq. imax)
     1 print*, "ERROR:Max._number_of_iterations_in_secant!" 
      if(abs(fr) .gt. tol)
     1 print* , "ERROR:Secant_method_not_converged!"
c
c     Return plastic multiplier:
      dlambda=xr
c
      return
      end

c_________________________________________________________________________
c 
c                       INVERSE 3x3
c_________________________________________________________________________

!***********************************************************************************************************************************
!  M33INV  -  Compute the inverse of a 3x3 matrix.
!
!  A       = input 3x3 matrix to be inverted
!  AINV    = output 3x3 inverse of matrix A
!  OK_FLAG = (output) .TRUE. if the input matrix could be inverted, and .FALSE. if the input matrix is singular.
!  Open source subroutine retreived from http://web.hku.hk/~gdli/UsefulFiles/matrix/m33inv_f90.txt
!***********************************************************************************************************************************

      SUBROUTINE M33INV (A, AINV, SING_FLAG)

      IMPLICIT NONE

      DOUBLE PRECISION  A(3,3)
      DOUBLE PRECISION  AINV(3,3)
      LOGICAL SING_FLAG

      DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-10
      DOUBLE PRECISION :: DET
      DOUBLE PRECISION, DIMENSION(3,3) :: COFACTOR

      DET =   A(1,1)*A(2,2)*A(3,3)  
     1       - A(1,1)*A(2,3)*A(3,2)  
     2       - A(1,2)*A(2,1)*A(3,3)  
     3       + A(1,2)*A(2,3)*A(3,1)  
     4       + A(1,3)*A(2,1)*A(3,2)  
     5       - A(1,3)*A(2,2)*A(3,1)

      IF (ABS(DET) .LE. EPS) THEN
         AINV = 0.0D0
         SING_FLAG = .TRUE.
         RETURN
      END IF

      COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
      COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
      COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
      COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

      AINV = TRANSPOSE(COFACTOR) / DET

      SING_FLAG = .FALSE.

      RETURN

      END SUBROUTINE M33INV


c_________________________________________________________________________
c 
c                       COMPLEX INVERSE 3x3
c_________________________________________________________________________

!***********************************************************************************************************************************
!  M33INV  -  Compute the inverse of a 3x3 matrix.
!
!  A       = input 3x3 matrix to be inverted
!  AINV    = output 3x3 inverse of matrix A
!  OK_FLAG = (output) .TRUE. if the input matrix could be inverted, and .FALSE. if the input matrix is singular.
!  Open source subroutine retreived from http://web.hku.hk/~gdli/UsefulFiles/matrix/m33inv_f90.txt
!***********************************************************************************************************************************

      SUBROUTINE M33INV_comp (A, AINV, SING_FLAG)

      IMPLICIT NONE

      complex*16  A(3,3)
      complex*16  AINV(3,3)
      LOGICAL SING_FLAG

      DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-10
      complex*16 :: DET
      complex*16, DIMENSION(3,3) :: COFACTOR

      DET =   A(1,1)*A(2,2)*A(3,3)  
     1       - A(1,1)*A(2,3)*A(3,2)  
     2       - A(1,2)*A(2,1)*A(3,3)  
     3       + A(1,2)*A(2,3)*A(3,1)  
     4       + A(1,3)*A(2,1)*A(3,2)  
     5       - A(1,3)*A(2,2)*A(3,1)

      IF (ABS(DET) .LE. EPS) THEN
         AINV = 0.0D0
         SING_FLAG = .TRUE.
         RETURN
      END IF

      COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
      COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
      COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
      COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

      AINV = TRANSPOSE(COFACTOR) / DET

      SING_FLAG = .FALSE.

      RETURN

      END SUBROUTINE M33INV_comp


c_________________________________________________________________________
c 
c                       voigt2full
c_________________________________________________________________________
      SUBROUTINE voigt2full (x_full, x_voigt)

      IMPLICIT NONE
      double precision x_full(3,3), x_voigt(6)

      x_full(1,1) = x_voigt(1)
      x_full(2,2) = x_voigt(2)
      x_full(3,3) = x_voigt(3)
      x_full(1,2) = x_voigt(4)
      x_full(2,3) = x_voigt(5)
      x_full(1,3) = x_voigt(6)
      x_full(2,1) = x_full(1,2)
      x_full(3,2) = x_full(2,3)
      x_full(3,1) = x_full(1,3)

      RETURN
      END

c_________________________________________________________________________
c 
c                       full2voigt
c_________________________________________________________________________
      SUBROUTINE full2voigt (x_full, x_voigt)

      IMPLICIT NONE
      double precision x_full(3,3), x_voigt(6)

      x_voigt(1) = x_full(1,1)
      x_voigt(2) = x_full(2,2)
      x_voigt(3) = x_full(3,3)
      x_voigt(4) = x_full(1,2)
      x_voigt(5) = x_full(2,3)
      x_voigt(6) = x_full(1,3)

      RETURN
      END
