C----|--1---------2---------3---------4---------5---------6---------7-|


c_________________________________________________________________________	
c	
c                        MAIN CODE
c_________________________________________________________________________


      subroutine umat43(cm, eps, sig, epsp, hsv, dt1, etype,
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
c     eps(1)= local x strain INCREMENT
c     eps(2)= local y strain INCREMENT
c     eps(3)= local z strain INCREMENT
c     eps(4)= local xy strain INCREMENT
c     eps(5)= local yz strain INCREMENT
c     eps(6)= local zx strain INCREMENT
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
c     hsv(1)= plastic strain TT−axis
c     hsv(2)= plastic strain TL−axis
c     hsv(3)= plastic strain TW−axis
c     hsv(4)=non−local plastic strain TT−axis
c     hsv(5)= local linear x strain previous time step
c     hsv(6)= local linear y strain previous time step
c     hsv(7)= local linear z strain previous time step
c     hsv(8)= local linear xy strain previous time step
c     hsv(9)= local linear yz strain previous time step
c     hsv(10)= local linear xz strain previous time step
c     hsv(11)= yield law value
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
      dimension cm(*), eps(*), hsv(*)
      integer nnpcrv
      real dt1
      integer imax
      character*5 etype
      logical failel, reject
      integer*8 idele, num_hv, elsiz

      real tol, m, dlambda
      real Ett, Ell, Eww, Gtl, Glw, Gtw, ytt, ytl, ytw
      real nutl, nult, nutw, nuwt, nulw, nuwl, Delta
      real yld0tt, yld0tl, yld0tw, yldcrtt, yldcrtl, yldcrtw
      real flow11, flow12, flow13
      real epsp_tt, epsp_tl, epsp_tw, ftrial, epsd, hd, g
      real eps_lin(6), eps_old(6), deps(6), sig_trial(6), sig(6)
      real var_nonloc, var_loc, dnl, lr, le, eps_star

      double precision c(3,3), u(3,3), V(3,3), D(3), C4(6,6)
c
c     Initialize/set parameters:−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
      num_hv = int(cm(19))! Number of history variables in use
      toll = 0.000001! Tolerance for yieldsurface
      imax = 100! Max. number of iterations for plastic multiplier
      m = 1.0d0 ! Parameter in yield surface
      dlambda= 0.0d0 ! Plastic multiplier
      flow11 = 0.0d0
      flow12 = 0.0d0
      flow13 = 0.0d0
      sig_trial = 0.0d0
      deps =0.0 d0
c     Initiate linear strain tensor parameters
      c = 0.0d0
      u = 0.0d0
      V = 0.0d0
      D = 0.0d0
c
c      if(ncycle.eq.1) then
c        call usermsg(’mat43’)
c      end if
c
c     Assign material parameters/read history variables:−−−−−−−------------
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
      epsp_tl=hsv(2)
      epsp_tw=hsv(3)
c
c     Strains of previous timestep:
      eps_old(1)= hsv(5)
      eps_old(2)= hsv(6)
      eps_old(3)= hsv(7)
      eps_old(4)= hsv(8)
      eps_old(5)= hsv(9)
      eps_old(6)= hsv(10)
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
      C4(4,4)= Gtl
      C4(5,5)= Glw
      C4(6,6)= Gtw
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
      if (etype .eq. "solid") then
c
c       Compute right cauchy−green deformation tensor(symmetric)
c       [ c1=c11 , c2=c22 , c3=c33 , c4=c12=c21 , c5=c23=c32 , c6=c13=c31 ]
c
        c(1,1)= hsv(num_hv+1)*hsv(num_hv+1)+hsv(num_hv+2)*hsv(num_hv+2)
     1   +hsv(num_hv+3)*hsv(num_hv+3)
        c(2,2)= hsv(num_hv+4)*hsv(num_hv+4)+hsv(num_hv+5)*hsv(num_hv+5)
     1   +hsv(num_hv+6)*hsv(num_hv+6)
        c(3,3)= hsv(num_hv+7)*hsv(num_hv+7)+hsv(num_hv+8)*hsv(num_hv+8)
     1   +hsv(num_hv+9)*hsv(num_hv+9)
        c(1,2)= hsv(num_hv+1)*hsv(num_hv+4)+hsv(num_hv+2)*hsv(num_hv+5)
     1   +hsv(num_hv+3)*hsv(num_hv+6)
        c(2,3)= hsv(num_hv+4)*hsv(num_hv+7)+hsv(num_hv+5)*hsv(num_hv+8)
     1   +hsv(num_hv+6)*hsv(num_hv+9)
        c(1,3)= hsv(num_hv+1)*hsv(num_hv+7)+hsv(num_hv+2)*hsv(num_hv+8)
     1   +hsv(num_hv+3)*hsv(num_hv+9)
        c(2,1)= c(1,2)
        c(3,2)= c(2,3)
        c(3,1)= c(1,3)
c
c       Determine eigenvalues/vectors of Right deformation tensor:
        call DSYEVJ3(c, V, D) ! downloaded open source tool
c
        D = sqrt(D)
c
c       Compute right stretch tensor:
        u(1,1) = D(1)
        u(2,2) = D(2)
        u(3,3) = D(3)
c
        u = matmul(matmul(V,u),transpose(V))
c
c       Compute linear (Biot) strain tensor
        eps_lin(1)=u(1,1)-1
        eps_lin(2)=u(2,2)-1
        eps_lin(3)=u(3,3)-1
        eps_lin(4)=u(1,2)
        eps_lin(5)=u(2,3)
        eps_lin(6)=u(1,3)
c
c       Compute stresses (return mapping algorithm): −−−−−−−−−−−−−−−−−−−
c
c       Compute cauchy stress:
c       (pay attention to material axis: 1=t,2=l,3=w,4=tl,5=lw,6=tw)
        deps = eps_lin-eps_old
        sig_trail= sig+matmul(C4,deps)
c
c       Check for compressive stresses
        if(sig(1) .le. 0.0d0) then
c
c         Compute yyield stresses
          call compute_yield_stress(epsp_tt, ytt, yld0tt, yldcrtt, epsd, hd, g)
          call compute_yield_stress(epsp_tt, ytl, yld0tl, yldcrtl, epsd, hd, g)
          call compute_yield_stress(epsp_tt, ytw, yld0tw, yldcrtw, epsd, hd, g)






      endif
      print *, "Hello from subroutine"
c
c
      print "(3f12.2)", cm(1), cm(10), C4(2,3)
c
c
      print *, "Goodbye from subroutine"
c
      END









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
      real m, f , sig11, sig12, sig13, ytt, ytl, ytw
c
      f =(sig11/ytt)+((sig12/ytl)**2.0d0+(sig13/ytw)**2.0d0)**(m*0.5d0)
     1 −1.0d0
c
      return
      end