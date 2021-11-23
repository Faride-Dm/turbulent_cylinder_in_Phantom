!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: setup
!
!  DESCRIPTION:
! this module does setup
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id: 65419773682803ee88a3ea112e9bdca4117e5bfb $
!
!  RUNTIME PARAMETERS:
!    ang_Bomega       -- Angle (degrees) between B and rotation axis
!    angvel           -- angular velocity in rad/s
!    cs_cylinder        -- sound speed in cylinder in code units
!    density_contrast -- density contrast in code units
!    dist_unit        -- distance unit (e.g. au)
!    dusttogas        -- dust-to-gas ratio
!    form_binary      -- the intent is to form a central binary
!    mass_unit        -- mass unit (e.g. solarm)
!    masstoflux       -- mass-to-magnetic flux ratio in units of critical value
!    np               -- actual number of particles in cylinder
!    pmass_dusttogas  -- dust-to-gas particle mass ratio
!    r_cylinder         -- radius of cylinder in code units
!    rho_pert_amp     -- amplitude of density perturbation
!    totmass_cylinder   -- mass of cylinder in code units
!
!  DEPENDENCIES: boundary, centreofmass, dim, eos, infile_utils, io,
!    kernel, options, part, physcon, prompting, ptmass, setup_params,
!    spherical, timestep, unifcyl           , units
!+
!--------------------------------------------------------------------------
module setup
 use part,    only:mhd
 use dim,     only:use_dust,maxvxyzu
 use options, only:calc_erot
 implicit none
 public :: setpart

 private
 !--private module variables
 real :: xmini(3), xmaxi(3)
 real :: density_contrast,totmass_cylinder,r_cylinder,cs_cylinder,h_cylinder
 real :: angvel, masstoflux,dusttogas,pmass_dusttogas,ang_Bomega,Bzero
 real :: rho_pert_amp
 real(kind=8)                 :: udist,umass
 integer                      :: np
 logical                      :: binary
 character(len=20)            :: dist_unit,mass_unit
 character(len= 1), parameter :: labelx(3) = (/'x','y','z'/)
 !smhr_begin
 logical                      :: perturbation
 !smhr_end
 
 !fdm_begin
 !logical                      :: turbulence
 !fdm_end

contains

!----------------------------------------------------------------
!+
!  setup for a cylinder-in-a-box
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use physcon,      only:pi,solarm,hours,years,au
 !fdm_begin
 use velfield,     only:set_velfield_from_cubes
 use io,           only:fatal,master
 !fdm_end
 use setup_params, only:rhozero,npart_total,rmax,ihavesetupB
 use io,           only:master
 !smhr_begin
 !  use unifdis,      only:set_unifdis
 use unifcyl,      only:set_unifcyl
 ! use spherical,    only:set_unifdis_sphereN
 use cylindrical,  only:set_unifdis_cylinderN
 !smhr_end
 use boundary,     only:set_boundary,xmin,xmax,ymin,ymax,zmin,zmax,dxbound,dybound,dzbound
 use prompting,    only:prompt
 use units,        only:set_units,select_unit,utime,unit_density,unit_Bfield
 !fdm_begin
 !  use eos,          only:polyk2,ieos
 use setvfield,    only:normalise_vfield
 use datafiles,    only:find_phantom_datafile
 use eos,          only:polyk2,ieos,dens_cylinder,dens_medium,dens_contrast
 !fdm_end
 use part,         only:Bxyz,Bextx,Bexty,Bextz,igas,idust,set_particle_type
 use timestep,     only:dtmax,tmax,dtmax_dratio,dtmax_min
 use ptmass,       only:icreate_sinks,r_crit,h_acc,h_soft_sinksink,rho_crit_cgs,h_soft_sinkgas
 use centreofmass, only:reset_centreofmass
 use options,      only:nfulldump,rhofinal_cgs
 use kernel,       only:hfact_default
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real(kind=8)       :: h_acc_in
 real               :: totmass,vol_box,psep,psep_box
 !smhr_begin
!  real               :: vol_cylinder,dens_cylinder,dens_medium,cs_medium,angvel_code,przero
 real               :: vol_cylinder,cs_medium,angvel_code,przero
 real               :: znew,zold,diff,lambda,kwave
 !smhr_end
 real               :: totmass_box,t_ff,r2,area,rmasstoflux_crit
 real               :: rxy2,rxyz2,phi,dphi,lzbox,lybox,lxbox,hi,epotgrav
 !fdm_begin
 !hi,epotgrav are added to real 
 character(len=20), parameter :: filevx = 'cube_v1.dat'
 character(len=20), parameter :: filevy = 'cube_v2.dat'
 character(len=20), parameter :: filevz = 'cube_v3.dat'
 character(len=120)           :: filex,filey,filez,filein,fileset
 !fdm_end
 integer            :: i,nx,np_in,npartcylinder,npmax,ierr
 logical            :: iexist
 logical            :: make_sinks = .true.
 character(len=100) :: filename
 character(len=40)  :: fmt
 character(len=10)  :: string,h_acc_char

 npmax = size(xyzh(1,:))
 filename=trim(fileprefix)//'.setup'
 print "(/,1x,63('-'),1(/,a),/,1x,63('-'),/)",&
   '  cylinder-in-box setup: Almost Archimedes'' greatest achievement.'
 inquire(file=filename,exist=iexist)
 if (iexist) then
    call read_setupfile(filename,ierr)
    np_in = np
    if (ierr /= 0) then
       if (id==master) call write_setupfile(filename)
       stop
    endif
 elseif (id==master) then
    print "(a,/)",trim(filename)//' not found: using interactive setup'
    dist_unit = 'pc'
    mass_unit = 'solarm'
    ierr = 1
    do while (ierr /= 0)
       call prompt('Enter mass unit (e.g. solarm,jupiterm,earthm)',mass_unit)
       call select_unit(mass_unit,umass,ierr)
       if (ierr /= 0) print "(a)",' ERROR: mass unit not recognised'
    enddo
    ierr = 1
    do while (ierr /= 0)
       call prompt('Enter distance unit (e.g. au,pc,kpc,0.1pc)',dist_unit)
       call select_unit(dist_unit,udist,ierr)
       if (ierr /= 0) print "(a)",' ERROR: length unit not recognised'
    enddo
    !
    ! units
    !
    call set_units(dist=udist,mass=umass,G=1.d0)
    !
    ! prompt user for settings
    !
    npmax = int(2.0/3.0*size(xyzh(1,:))) ! approx max number allowed in cylinder given size(xyzh(1,:))
    if (npmax < 300000) then
       np = npmax
    else if (npmax < 1000000) then
       np = 300000
    else
       np = 1000000
    endif
    call prompt('Enter the approximate number of particles in the cylinder',np,0,npmax)
    np_in      = np
    r_cylinder = 0.1
    call prompt('Enter the radius of the cylinder in units of '//dist_unit,r_cylinder,0.01)
    h_cylinder = 1.5
    call prompt('Enter the length of the cylinder in units of '//dist_unit,h_cylinder,1.0)
    lzbox      = 2.0
    call prompt('Enter the z-length of the box in units of '//dist_unit,lzbox,h_cylinder)
    lxbox      = 0.4
    call prompt('Enter the xy-length of the box in units of '//dist_unit,lxbox,4.*r_cylinder)
    lybox = lxbox
!     is_box   = .true.
!     call prompt('Enter the box size in units of cylinder radius: ',lbox,1.)
!     do i=1,3
!        ! note that these values will be saved to the .setup file rather than lbox so user can convert
!        ! box to a rectangle if desired
!        xmini(i) = -0.5*(lbox*r_cylinder)
!        xmaxi(i) = -xmini(i)
!     enddo
    
    xmini(1) = -lxbox
    xmaxi(1) = lxbox
    xmini(2) = -lybox
    xmaxi(2) = lybox
    xmini(3) = -lzbox
    xmaxi(3) = lzbox

    density_contrast = 30.0
    call prompt('Enter density contrast between cylinder and box ',density_contrast,1.)
    dens_contrast = density_contrast

    totmass_cylinder = 50.0
    call prompt('Enter total mass in cylinder in units of '//mass_unit,totmass_cylinder,0.)

    binary = .false.
!     call prompt('Do you intend to form a binary system?',binary)

    if (binary) then
       cs_cylinder = 0.1623
    else
       cs_cylinder = 3.4771871
    endif
    write(string,"(es10.3)") udist/utime
    call prompt('Enter sound speed in cylinder in units of '//trim(adjustl(string))//' cm/s',cs_cylinder,0.)

    if (binary) then
       angvel = 1.006d-12
    else
       angvel = 1.77d-13
    endif
    angvel = 0.0
!     call prompt('Enter angular rotation speed in rad/s ',angvel,0.)

    if (mhd) then
!        masstoflux =   5.0
       ang_Bomega = 180.0
!        call prompt('Enter mass-to-flux ratio in units of critical value ',masstoflux,0.)
       Bzero = 113.39
       call prompt('Enter magnitude of the magnetic field', Bzero, 0.,1133.9)
       call prompt('Enter the angle (degrees) between B and the rotation axis? ',ang_Bomega)
    endif

 !fdm_begin
 
 !   if (mhd) then
 !   	Bz_0 = 1.4142e-5
 !   	print *, ''
 !   	print *, 'Setting up MHD turbulence: (with uniform intial magnetic field in z-direction)'
 !   	print *, ''
 !   	if (id==master) call prompt('Enter initial magnetic field strength ',Bz_0)
 !   endif
 
 !fdm_end
    
    if (use_dust) then
       dusttogas = 0.01
       call prompt('Enter dust-to-gas ratio ',dusttogas,0.)
       pmass_dusttogas = dusttogas*10.0
       call prompt('Enter dust-to-gas particle mass ratio ',pmass_dusttogas,0.)
    endif

    perturbation = .true.
    call prompt('Do you intend to perturb the filament?',perturbation)
    
    if ((binary) .or. (perturbation)) then
       rho_pert_amp = 0.1
       call prompt('Enter the amplitute of the density perturbation ',rho_pert_amp,0.0,0.99)
    endif
    !
 ! ask about sink particle details; these will not be saved to the .setup file since they exist in the .in file
    call prompt('Do you wish to dynamically create sink particles? ',make_sinks)
    if (make_sinks) then
       call prompt('Enter the critical density for creation of sink particles ',rho_crit_cgs, 1.e-18)
       if (binary) then
          h_acc_char  = '3.35au'
       else
!           h_acc_char  = '1.0d-2'
          h_acc_char  = '108au'
       endif
       call prompt('Enter the accretion radius of the sink (with units; e.g. au,pc,kpc,0.1pc) ',h_acc_char)
       call select_unit(h_acc_char,h_acc_in,ierr)
       h_acc = h_acc_in
!        print*, h_acc_in,h_acc, h_acc_char
       if (ierr==0 ) h_acc = h_acc/udist
       r_crit        = 2.0*h_acc
!fdm_begin
       h_soft_sinksink = 0.00 ! see Federrath et al. (2010), ApJ, 713, 269
       h_soft_sinkgas = 0.00 ! see Federrath et al. (2010), ApJ, 713, 269
!fdm_end
       icreate_sinks = 1
       if (binary) h_soft_sinksink = 0.4*h_acc
    else
       icreate_sinks = 0
    endif
    call prompt('Enter ieos ',ieos)
    
! fdm_begin
    
    !turbulence = .true.
    !Mach = 1.
    !call prompt('Do you intend to use magnetic turbulence?',turbulence)
    !call prompt('Enter the magnitude of the mach number',Mach,0.,3.)
    
! fdm_end
    
    !
    ! write default input file
    !
    call write_setupfile(filename)
    print "(a)",'>>> rerun phantomsetup using the options set in '//trim(filename)//' <<<'
 else
    stop
 endif
 !
 ! units
 !
 call set_units(dist=udist,mass=umass,G=1.d0)
 !
 ! boundaries
 !
 call set_boundary(xmini(1),xmaxi(1),xmini(2),xmaxi(2),xmini(3),xmaxi(3))
 !
 ! general parameters
 !
 time        = 0.
 hfact       = hfact_default
 if (maxvxyzu >=4 ) then
    gamma    = 5./3.
 else
    gamma    = 1.
 endif
 rmax        = r_cylinder
 angvel_code = angvel*utime
 vol_box     = dxbound*dybound*dzbound
!  vol_cylinder  = 4./3.*pi*r_cylinder**3
 vol_cylinder  = pi*r_cylinder*r_cylinder*h_cylinder
 rhozero     = totmass_cylinder / vol_cylinder
 dens_cylinder = rhozero
 dens_medium = dens_cylinder/density_contrast
 cs_medium   = cs_cylinder*sqrt(density_contrast)
 totmass_box = (vol_box - vol_cylinder)*dens_medium
 totmass     = totmass_box + totmass_cylinder
 t_ff        = sqrt(3.*pi/(32.*dens_cylinder)) ! need to be checked
 !
 ! magnetic field
 !
 rmasstoflux_crit = 2./3.*0.53*sqrt(5./pi)
 if (mhd) then
!     area = pi*r_cylinder**2
!     if (masstoflux > tiny(masstoflux)) then
!fdm_begin
       !Bzero = totmass_cylinder/(area*masstoflux*rmasstoflux_crit)
!        call prompt('Enter magnitude of the magnetic field', Bzero, 0.0,50.0)
!fdm_end
!     else
!        Bzero = 0.
!     endif
    ihavesetupB = .true.
 else
    Bzero = 0.
 endif
 Bextx  = 0.
 Bexty  = 0.
 Bextz  = Bzero
 przero = cs_cylinder**2*dens_cylinder

 print "(a,i10)",' Input npart_cylinder = ',np
 print "(1x,50('-'))"
 print "(a)",'  Quantity         (code units)  (physical units)'
 print "(1x,50('-'))"
 fmt = "((a,1pg10.3,3x,1pg10.3),a)"
 print fmt,' Total mass       : ',totmass,totmass*umass,' g'
 print fmt,' Mass in cylinder   : ',totmass_cylinder,totmass_cylinder*umass,' g'
 print fmt,' Radius of cylinder : ',r_cylinder,r_cylinder*udist,' cm'
 print fmt,' Length of cylinder : ',h_cylinder,h_cylinder*udist,' cm'
 print fmt,' Density cylinder   : ',dens_cylinder,dens_cylinder*unit_density,' g/cm^3'
 print fmt,' Density medium   : ',dens_medium,dens_medium*unit_density,' g/cm^3'
 print fmt,' cs in cylinder     : ',cs_cylinder,cs_cylinder*udist/utime,' cm/s'
 print fmt,' cs in medium     : ',cs_medium,cs_medium*udist/utime,' cm/s'
 print fmt,' Free fall time   : ',t_ff,t_ff*utime/years,' yrs'
 print fmt,' Angular velocity : ',angvel_code,angvel,' rad/s'
 print fmt,' Omega*t_ff       : ',angvel_code*t_ff
 if (mhd) then
    print fmt,' B field (z)      : ',Bzero,Bzero*unit_Bfield*1.d6,' micro-G'
    print fmt,' Alfven speed     : ',Bzero/sqrt(dens_cylinder),Bzero/sqrt(dens_cylinder)*udist/utime,' cm/s'
    if (Bzero > 0.) then
       print fmt,' plasma beta      : ',przero/(0.5*Bzero*Bzero)
!        print fmt,' mass-to-flux     : ',totmass_cylinder/(area*Bzero)/rmasstoflux_crit
    endif
 endif
 if (use_dust) then
    print fmt,' dust-to-gas ratio: ',dusttogas,' '
    print fmt,' dust-to-gas particle mass ratio: ',pmass_dusttogas,' '
 endif
 print "(1x,50('-'))"
 !
 ! setup particles in the cylinder; use this routine to get N_cylinder as close to np as possible
 !
 call set_unifdis_cylinderN('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,psep,&
                    hfact,npart,np,xyzh,r_cylinder,h_cylinder,vol_cylinder,npart_total)
 print "(a,es10.3)",' Particle separation in cylinder = ',psep
 npartcylinder = npart
 if (np_in/=npartcylinder) np = npartcylinder
 !
 ! setup surrounding low density medium
 !
 psep_box = psep*(density_contrast)**(1./3.)  ! calculate psep in box
 call set_unifcyl('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,psep_box, &
                   hfact,npart,xyzh,rcylmin=r_cylinder,hcyl=h_cylinder,&
                   nptot=npart_total,cyl_surround_medium=.true.)
 print "(a,es10.3)",' Particle separation in low density medium = ',psep_box
 print "(a,i10,a)",' added ',npart-npartcylinder,' particles in low-density medium'
 print*, ""
 !
 ! set particle properties
 !
 npartoftype(:)    = 0
 npartoftype(igas) = npart
 massoftype(igas)  = totmass/npart_total
 do i = 1,npartoftype(igas)
    call set_particle_type(i,igas)
 enddo
 !
 ! reset to centre of mass
 !
 call reset_centreofmass(npart,xyzh,vxyzu)
 !
 ! Set dust
 if (use_dust) then
    ! particle separation in dust cylinder & sdjust for close-packed lattice
    psep = (vol_cylinder/pmass_dusttogas)**(1./3.)/real(nx)
    psep = psep*sqrt(2.)**(1./3.)
    call set_unifdis_sphereN('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,psep,&
                    hfact,npart,np,xyzh,r_cylinder,vol_cylinder,npart_total)
    npartoftype(idust) = npart_total - npartoftype(igas)
    massoftype(idust)  = totmass_cylinder*dusttogas/npartoftype(idust)
    !
    do i = npartoftype(igas)+1,npart
       call set_particle_type(i,idust)
    enddo
    !
    print "(a,4(i10,1x))", ' particle numbers: (gas_total, gas_cylinder, dust, total): ' &
                        , npartoftype(igas),npartcylinder,npartoftype(idust),npart
    print "(a,2es10.3)"  , ' particle masses: (gas,dust): ',massoftype(igas),massoftype(idust)
 else
    print "(a,3(i10,1x))", ' particle numbers: (cylinder, low-density medium, total): ' &
                        , npartcylinder, npart-npartcylinder,npart
    print "(a,es10.3)",' particle mass = ',massoftype(igas)
 endif
 !
 ! temperature set to give a pressure equilibrium
 !
 polyk  = cs_cylinder**2
 polyk2 = cs_medium**2
 !
 !--Stretching the spatial distribution to perturb the density profile (for binary==.true. only)
 !
 if (binary) then
    do i = 1,npart
       rxy2  = xyzh(1,i)*xyzh(1,i) + xyzh(2,i)*xyzh(2,i)
       rxyz2 = rxy2 + xyzh(3,i)*xyzh(3,i)
       if (rxyz2 <= r_cylinder**2) then
          phi       = atan(xyzh(2,i)/xyzh(1,i))
          if (xyzh(1,i) < 0.0) phi = phi + pi
          dphi      = 0.5*rho_pert_amp*sin(2.0*phi)
          phi       = phi - dphi
          xyzh(1,i) = sqrt(rxy2)*cos(phi)
          xyzh(2,i) = sqrt(rxy2)*sin(phi)
       endif
    enddo
 endif
 if (perturbation) then
    lambda = 0.75 ! in parsec
    kwave = 2.*pi/lambda
    do i = 1, npart
       rxy2  = xyzh(1,i)*xyzh(1,i) + xyzh(2,i)*xyzh(2,i)
       if ((xyzh(3,i) > -h_cylinder/2. .and. xyzh(3,i) < h_cylinder/2.) .and. &
       &  (rxy2 < r_cylinder**2)) then
          diff = 0.
          znew = xyzh(3,i)
          do
            zold = znew
            znew = xyzh(3,i) + rho_pert_amp*(sin(kwave*znew))/kwave
            diff = abs((znew - zold)/lambda)
            if (diff < 1.0e-12) exit
          enddo
          xyzh(3,i) = znew
       endif
    enddo
 endif 
 !
 ! velocity field corresponding to uniform rotation
 ! fdm_begin
 ! lines 474-494 are changed to comments and --Set velocities (from pre-made velocity cubes) is added.  
 !
 !do i=1,npart
 !   r2 = dot_product(xyzh(1:2,i),xyzh(1:2,i))
 !   if (r2 < r_cylinder**2) then
 !      vxyzu(1,i) = -angvel_code*xyzh(2,i)
 !      vxyzu(2,i) =  angvel_code*xyzh(1,i)
 !      vxyzu(3,i) = 0.
 !      if (maxvxyzu >= 4) vxyzu(4,i) = 1.5*polyk
 !   else
 !      vxyzu(1:3,i) = 0.
 !      if (maxvxyzu >= 4) vxyzu(4,i) = 1.5*polyk2
 !   endif
 !   if (mhd) then
 !      Bxyz(:,i) = 0.
 !      Bxyz(1,i) = Bzero*sin(ang_Bomega*pi/180.0)
 !      Bxyz(3,i) = Bzero*cos(ang_Bomega*pi/180.0)
 !   endif
 !enddo
 
 ! fdm_begin
 
 hi            = hfact*(massoftype(1)/rhozero)**(1./3.)
 ! hi and Mach is addeed to: subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,hi,time,fileprefix,Mach)
 epotgrav      = 3./5.*totmass**2/rmax      ! Gravitational potential energy

 
 !--Set velocities (from pre-made velocity cubes)
 filex = find_phantom_datafile(filevx,'velfield')
 filey = find_phantom_datafile(filevy,'velfield')
 filez = find_phantom_datafile(filevz,'velfield')
 
 call set_velfield_from_cubes(xyzh,vxyzu,npartoftype(igas),filex,filey,filez,&
                              1.,rmax,.false.,ierr)
 if (ierr /= 0) call fatal('setup','error setting up velocity field')

 !--Normalise the energy
 call normalise_vfield(npart,vxyzu,ierr,ke=epotgrav)
 if (ierr /= 0) call fatal('setup','error normalising velocity field')


 ! fdm_end
 
 !
 ! set default runtime parameters if .in file does not exist
 !
 filename=trim(fileprefix)//'.in'
 inquire(file=filename,exist=iexist)
 dtmax = t_ff/20.  ! Since this variable can change, always reset it if running phantomsetup
 if (.not. iexist) then
    if (binary) then
       tmax      = 13.33
    else
       tmax      = 10.
    endif
!     ieos         = 1
    nfulldump    = 1
    calc_erot    = .false.
    dtmax_dratio = 1.258
!     rho_crit_cgs  = 1.e-19
    if (make_sinks) then
       dtmax_min = dtmax/8.0
    else
       dtmax_min = 0.0
       rhofinal_cgs = 0.15
    endif
 endif

end subroutine setpart

!----------------------------------------------------------------
!+
!  write parameters to setup file
!+
!----------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only: write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter           :: iunit = 20
 integer                      :: i

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for cylinder-in-box setup routines'
 write(iunit,"(/,a)") '# units'
 call write_inopt(dist_unit,'dist_unit','distance unit (e.g. au)',iunit)
 call write_inopt(mass_unit,'mass_unit','mass unit (e.g. solarm)',iunit)
 write(iunit,"(/,a)") '# resolution'
 call write_inopt(np,'np','actual number of particles in cylinder',iunit)
 write(iunit,"(/,a)") '# options for box'
 do i=1,3
    call write_inopt(xmini(i),labelx(i)//'min',labelx(i)//' min',iunit)
    call write_inopt(xmaxi(i),labelx(i)//'max',labelx(i)//' max',iunit)
 enddo
 write(iunit,"(/,a)") '# intended result'
 call write_inopt(binary,'form_binary','the intent is to form a central binary',iunit)
 call write_inopt(perturbation,'perturbation','the intent is to perturb cylinder density',iunit)
 write(iunit,"(/,a)") '# options for cylinder'
 call write_inopt(r_cylinder,'r_cylinder','radius of cylinder in code units',iunit)
 call write_inopt(density_contrast,'density_contrast','density contrast in code units',iunit)
 call write_inopt(totmass_cylinder,'totmass_cylinder','mass of cylinder in code units',iunit)
 call write_inopt(cs_cylinder,'cs_cylinder','sound speed in cylinder in code units',iunit)
 call write_inopt(angvel,'angvel','angular velocity in rad/s',iunit)
 if (mhd) then
!     call write_inopt(masstoflux,'masstoflux','mass-to-magnetic flux ratio in units of critical value',iunit)
    call write_inopt(Bzero,'Bzero','magnitude of the magnetic field in code units',iunit)
    call write_inopt(ang_Bomega,'ang_Bomega','Angle (degrees) between B and rotation axis',iunit)
 endif
 if (use_dust) then
    call write_inopt(dusttogas,'dusttogas','dust-to-gas ratio',iunit)
    call write_inopt(pmass_dusttogas,'pmass_dusttogas','dust-to-gas particle mass ratio',iunit)
 endif
 if (perturbation) then
    call write_inopt(rho_pert_amp,'rho_pert_amp','amplitude of density perturbation',iunit)
 endif
 close(iunit)

end subroutine write_setupfile

!----------------------------------------------------------------
!+
!  Read parameters from setup file
!+
!----------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils, only: open_db_from_file,inopts,read_inopt,close_db
 use io,           only: error
 use units,        only: select_unit
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter            :: iunit = 21
 integer                       :: i,nerr
 type(inopts), allocatable     :: db(:)

 print "(a)",' reading setup options from '//trim(filename)
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(mass_unit,'mass_unit',db,ierr)
 call read_inopt(dist_unit,'dist_unit',db,ierr)
 call read_inopt(binary,'form_binary',db,ierr)
 call read_inopt(np,'np',db,ierr)
 do i=1,3
    call read_inopt(xmini(i),labelx(i)//'min',db,ierr)
    call read_inopt(xmaxi(i),labelx(i)//'max',db,ierr)
 enddo
 call read_inopt(r_cylinder,'r_cylinder',db,ierr)
 call read_inopt(density_contrast,'density_contrast',db,ierr)
 call read_inopt(totmass_cylinder,'totmass_cylinder',db,ierr)
 call read_inopt(cs_cylinder,'cs_cylinder',db,ierr)
 call read_inopt(angvel,'angvel',db,ierr)
 if (mhd) then
!     call read_inopt(masstoflux,'masstoflux',db,ierr)
    call read_inopt(Bzero,'Bzero',db,ierr)
    call read_inopt(ang_Bomega,'ang_Bomega',db,ierr)
 endif
 if (use_dust) then
    call read_inopt(dusttogas,'dusttogas',db,ierr)
    call read_inopt(pmass_dusttogas,'pmass_dusttogas',db,ierr)
 endif
 if (binary) then
    call read_inopt(rho_pert_amp,'rho_pert_amp',db,ierr)
 endif
 if (perturbation) then
    call read_inopt(rho_pert_amp,'rho_pert_amp',db,ierr)
 endif
 call close_db(db)
 !
 ! parse units
 call select_unit(mass_unit,umass,nerr)
 if (nerr /= 0) then
    call error('setup_cylinderinbox','mass unit not recognised')
    ierr = ierr + 1
 endif
 call select_unit(dist_unit,udist,nerr)
 if (nerr /= 0) then
    call error('setup_cylinderinbox','length unit not recognised')
    ierr = ierr + 1
 endif
 !
 if (ierr > 0) then
    print "(1x,a,i2,a)",'Setup_cylinderinbox: ',nerr,' error(s) during read of setup file.  Re-writing.'
 endif
 
 ! fdm_begin
 
 ! setup preferred values of .in file
 !filename= trim(fileprefix)//'.in'
 !inquire(file=filename,exist=iexist)
 !if (.not. iexist) then
 !   tmax         = 1.00   ! run for 20 turbulent crossing times
 !   dtmax        = 0.0025
 !   nfulldump    = 5      ! output 4 full dumps per crossing time
 !   beta         = 4      ! legacy from Price & Federrath (2010), haven't checked recently if still required
 !endif
 !npart = 0
 !npart_total = 0
 
 ! fdm_end 

end subroutine read_setupfile
!----------------------------------------------------------------
end module setup

