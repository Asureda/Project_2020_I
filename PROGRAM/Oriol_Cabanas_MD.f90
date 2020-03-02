PROGRAM DYMANICS
    use
    !----------------------------------------------------------------------------!
    !              MOLECULAR DYNAMICS SIMULATION                                 !
    !----------------------------------------------------------------------------!
    !Program to simulate de time integration of a system on N number of particles!
    !and a density. The program computes in reduced units and real units the     !
    !kinetik energy, potencial energy, total energy, temperature, momentum,      !
    !pressure and the readial distribution function.                             !
    !It can be added a contact with a termal bath                                !
    !                                                                            !
    !      Oriol Cabanas Jan, 2020 (Last mod. Jan 2020)                          !
    !----------------------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER :: n_particles,M,i,j,k,n,n_radial,n_gr_meas,n_meas,n_meas_gr,seed
    REAL*8 :: density,L,a,sigma,epsilon,mass,T_ini,T_therm,T_therm_prov,mom(3),mo,momentum,kin,pot,temp_instant,pressure,dx_radial
    REAL*8 :: ta,tb,h,t,t_prod,k_B,time_re,energy_re,temp_re,dist_re,press_re,n_mols,total_mass,n_avog, rho_re
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: r,v,F
    REAL*8, DIMENSION(:),ALLOCATABLE :: g_r
    REAL*8, EXTERNAL :: PBC1, NO_PBC
    REAL*8 :: KINETIK
    COMMON/PARAMETERS/n_particles,M,density,L,a
    OPEN(20,FILE='data.log')
    !-------------------------------------------------------------
    !         PARAMITTERS FOR THE SIMULATION
    !-------------------------------------------------------------
    !n_particles=256  !number of particles
    !density=2d-1   !reduced units
    ta=0d0          !md time
    !tb=5d0          !md time
    !h=3d-4            !md time
    !t_prod=0.0       !md time
    !n_meas=1       !steps between every mesure of observables
    !----------
    !T_ini=10d0         !reduced units
    !T_therm=10d0       !reduced units
    !sigma=3.4         !Angstroms
    !epsilon=0.998       !kJ/mol
    !mass=40d0          !g/mol
    !k_B=0.008314462   !kJ/mol.K
    !n_avog=6.022d23   !part/mol
    !----------------------------------------------------------------
    !Some calculations for the radial distribution function
    !dx_radial=0.01   !Discretization of the radial values
    !n_radial=int(L/dx_radial)  !Number of radial values
    n_gr_meas=0       !number of measures (allways must be zero)
    !--------------------------------------------------------------
    !dimensional factors
    ! temp_re=epsilon/k_B   !Kelvin
    ! energy_re=epsilon     !kJ/mol
    ! dist_re=sigma         !Angstroms
    ! time_re=1d-13*sqrt(mass*sigma**2d0/epsilon)  !Seconds
    ! press_re=1d27*epsilon/(n_avog*sigma**3d0)          !Pascals
    ! n_mols=n_particles/n_avog
    ! total_mass=n_mols*mass
    ! rho_re=mass*1d24/(sigma**3d0*n_avog)
    !--------------------------------------------------------------
    !L=((n_particles*1d0)/density)**(1d0/3d0)    !Computing the magnitude of the system
    !M=nint(((n_particles*1d0)/4d0)**(1d0/3d0))   !Number of cells
    !a=L/(M*1d0)                                  !dimension of the cells
    WRITE(20,*)'GEOMETRIC COMPUTED PARAMETTERS'
    WRITE(20,*)'N=',n_particles,'L=',L,'M=',m,'a=',a
    !--------------------------------------------------------------
    !              DIMENSIONING VARIABLES
    !--------------------------------------------------------------
    ! ALLOCATE(r(n_particles,3),v(n_particles,3),F(n_particles,3),g_r(0:n_radial+1))
    ! seed=1996
    CALL SRAND(seed)
    !--------------------------------------------------------------
    !             INITIAL CONDITIONS
    !---------------------------------------------------------------
    CALL FCC_INITIALIZE(r)       !Initial positions
    CALL UNIFORM_VELO(v,T_ini)         !Initial velocities
    !Checking initial momentum
    DO i=1,3
        mo=0d0
        DO j=1,n_particles
            mo=mo+v(j,i)
        END DO
        MOM(I)=Mo
    END DO
    WRITE(20,*) 'px=',mom(1),'py=',mom(2),'pz=',mom(3)
    WRITE(20,*) 'initial temperature=',2d0*KINETIK(v)/(3d0*n_particles)
    WRITE(20,*)'reduced dimensions'
    WRITE(20,*)'time//energy//dist//temperature//pressure//rho_re'
    WRITE(20,*) time_re,energy_re,dist_re,temp_re,press_re,rho_re
    WRITE(20,*)'density//sigma/epsilon/mass/initial momentum'
    WRITE(20,*) density,sigma,epsilon,mass,(mom(1)**2d0+mom(2)**2d0+mom(3)**2d0)**0.5
    WRITE(20,*)'density in reduced'
    WRITE(20,*) density
    WRITE(20,*)'density in real units'
    WRITE(20,*) density*rho_re
    !---------------------------------------------------------
    !        MELTING ITERATIONS AT HIGH TEMPERATURE
    !----------------------------------------------------------
    !T_therm_prov=100d0  !Temperature of the initial melting
    CALL VELO_RESCALING(v,T_therm_prov)
    DO i=1,100
        !VELOCITY VERLET STEP
        CALL VELO_VERLET(r,v,f,h,pot,pressure)
        !THERMAL NORMALIZATION
        CALL ANDERSEN(v,T_therm_prov,h,n_particles)
    END DO
    print*,'final melting'
    !----------------------------------------------------------
    !         STARTING VELOCITY VERLET INTEGRATION
    !----------------------------------------------------------
    CALL VELO_RESCALING(v,T_ini)
    n=int(tb-ta)/h  !number of time inegration steps
    print*,'time iterations',n
    t=ta
    pressure=0.0
    g_r=0d0
    n_meas_gr=-10!setps between dirtribution function measures
    print*,'time iterations',n,'n meas',n_meas,'n dist',n_meas_gr
    !---------------------------------------------------------
    !        CREATING OUTPUT FILES
    !---------------------------------------------------------
    OPEN(11,FILE='termodynamics_reduced.dat')
    OPEN(12,FILE='termodynamics_real.dat')
    OPEN(13,FILE='distriv_funct.dat')
    OPEN(14,FILE='positions.xyz')
    print*,'initial'
    DO i=1,n
        t=ta+i*h !iteration time
        !VELOCITY VERLET STEP
        CALL VELO_VERLET(r,v,f,h,pot,pressure)
        !THERMAL NORMALISATION
        IF(T_therm.gt.0d0) THEN
            !CALL VELO_RESCALING(v,T_therm)
            CALL ANDERSEN(v,T_therm,h,n_particles)
        END IF
        !production cycles
        IF((mod(i,n_meas).eq.0).and.(n_meas.gt.0)) THEN
            kin=KINETIK(v)
            DO k=1,3
                mo=0d0
                DO j=1,n_particles
                    mo=mo+v(j,k)
                END DO
                MOM(k)=Mo
            END DO
            
            momentum=(mom(1)**2d0+mom(2)**2d0+mom(3)**2d0)**0.5
            temp_instant=2d0*kin/(3d0*n_particles)
            pressure=(density*temp_instant+pressure/(3d0*L**3d0))
            WRITE(11,*)t,kin,pot,(kin+pot),temp_instant,pressure,momentum
            WRITE(12,*)t*time_re,kin*energy_re,pot*energy_re,(kin+pot)*energy_re,temp_instant*temp_re,pressure*press_re
            !print*,'meas'
        END IF
        IF((mod(i,n_meas_gr).eq.0).and.(n_meas_gr.gt.0)) THEN
            !MEASURIN THE RADIAL DISTRIBUTION
            CALL RAD_DIST_INTER(r,g_r,dx_radial,n_radial,n_particles,L)
            n_gr_meas=n_gr_meas+1 !Adding one measure
        END IF
        IF(mod(i,1000).eq.0)THEN
            print*,'iteration',i,'of',n !Som prints to see if the program is still working well
        END IF
    END DO
    n_gr_meas=n_gr_meas+1
    CALL RAD_DIST_INTER(r,g_r,dx_radial,n_radial,n_particles,L) !Measuting the final distribution
    !----------------------------------------------------------------------------------------------
    !    FINAL prints
    !---------------------------------------------------------------------------------------------
    !COMPUTING THE ACTUAL DISTRIBUTION FUNCTION
    CALL RAD_DIST_FINAL(g_r,dx_radial,n_radial,n_gr_meas,density)
    !radial distribution function
    DO i=1,n_radial
        WRITE(13,*)dx_radial*i,g_r(i)
    END DO
    !wRITING THE FINAL POSITION DISTRIBUTION
    WRITE(14,*)n_particles
    WRITE(14,*)
    DO i=1,n_particles
        WRITE(14,*)'C',r(i,:)
    END DO
    print*,'program end'
END PROGRAM
SUBROUTINE VELO_VERLET(r,v,f,h,pot,pressure)
    INTEGER i,n_particles,M
    REAL*8 h,r0(n_particles,3),r(n_particles,3),v0(n_particles,3),v(n_particles,3)
    REAL*8 F(n_particles,3),f0(n_particles,3),pot,kin,cutoff,density,L,a,PBC2,pressure
    REAL*8, EXTERNAL :: PBC1
    COMMON/PARAMETERS/n_particles,M,density,L,a
    r0=r
    v0=v
    f0=f
    kin=0d0
    cutoff=0.99*L*5d-1
    !print*,'in verlet'
    CALL INTERACTION_CUTOFF(r,f0,pot,PBC1,cutoff,pressure)
    DO i=1,n_particles
        r(i,:)=r0(i,:)+v0(i,:)*h+5d-1*f0(i,:)*h*h
        r(i,1)=PBC2(r(i,1),L)
        r(i,2)=PBC2(r(i,2),L)
        r(i,3)=PBC2(r(i,3),L)
    END DO
    CALL INTERACTION_CUTOFF(r,f,pot,PBC1,cutoff,pressure)
    DO i=1,n_particles
        v(i,:)=v(i,:)+5d-1*(f0(i,:)+f(i,:))*h
    END DO
    !print*,'out verlet'
    RETURN
END SUBROUTINE
SUBROUTINE INTERACTION_CUTOFF(positions,F,E,PBC,cutoff,pressure)
!COMPUTING THE TOTAL INTERACTION ENERGY GIVEN AN EXTENAL FUNCTION FOR POTENCIAL AND BOUNDARY CONDITIONS
    IMPLICIT NONE
    INTEGER :: n_particles,M,i,j
    REAL*8 :: cutoff,density,L,a,PBC,E,pot,pressure
    REAL*8 :: dx,dy,dz,d,ff
    REAL*8, DIMENSION(n_particles,3) :: positions, F
    COMMON/PARAMETERS/n_particles,M,density,L,a
    !print*,'in inter'
    F=0d0
    E=0d0
    !print*,'hola-1'
    pressure=0.0
    !print*,'hola'
    DO i=1,n_particles
        DO j=i+1,n_particles
            dx=PBC(positions(i,1)-positions(j,1),L)
            dy=PBC(positions(i,2)-positions(j,2),L)
            dz=PBC(positions(i,3)-positions(j,3),L)
            d=(dx**2d0+dy**2d0+dz**2d0)**0.5
            CALL L_J(d,ff,pot,cutoff)
            F(i,1)=F(i,1)+ff*dx
            F(i,2)=F(i,2)+ff*dy
            F(i,3)=F(i,3)+ff*dz
            F(j,1)=F(j,1)-ff*dx
            F(j,2)=F(j,2)-ff*dy
            F(j,3)=F(j,3)-ff*dz
            E=E+pot
            !print*,'hola1'
            pressure=pressure+(ff*dx**2d0+ff*dy**2d0+ff*dz**2d0)
            !print*,'hola 2'
        END DO
    END DO
    !print*,'en inter'
    RETURN
END SUBROUTINE
SUBROUTINE L_J(d,F,pot,cutoff)
    IMPLICIT NONE
    REAL*8 d,cutoff,F,pot
    F=(48d0/d**14d0)-(24d0/d**8d0)
    pot=0d0
    IF (d<cutoff) THEN
        pot=4d0*((1d0/d**12d0)-(1d0/d**6d0))
    END IF
    RETURN
END SUBROUTINE
FUNCTION NO_PBC(x,L)
!FUNCTION THAT RETURNS THE DISTANCE WITH NO PERIODIC BOUNDARY CONDITIONS
    IMPLICIT NONE
    REAL*8 x,L,NO_PBC
    NO_PBC=x
    RETURN
END FUNCTION
FUNCTION PBC1(x,L)
!FUNCTION THAT RETURNS THE DISTANCE GIVEN PBC AND L
    IMPLICIT NONE
    REAL*8 x,L,PBC1
    PBC1=x-int(2d0*x/L)*L
    RETURN
END FUNCTION
FUNCTION PBC2(x,L)
!FUNCTION THAT RETURNS THE DISTANCE GIVEN PBC AND L
    IMPLICIT NONE
    REAL*8 x,L,PBC2
    IF(x.lt.0)THEN
        x=L+mod(x,L)
    ELSE IF(x.gt.L) THEN
        x=mod(x,L)
    END IF
    PBC2=x
    RETURN
END FUNCTION
SUBROUTINE FCC_INITIALIZE(positions)
    IMPLICIT NONE
    INTEGER :: n_particles,M,n,i,j,k
    REAL*8 :: density,L,a
    REAL*8 :: positions(n_particles,3)
    COMMON/PARAMETERS/n_particles,M,density,L,a
    n=1
    !print*,n
    DO i=0,M-1
        DO j=0,M-1
            DO k=0,M-1
                !print*,size(positions)
                positions(n,:)=a*(/i,j,k/)
                !print*,positions(n,:)
                positions(n+1,:)=positions(n,:)+a*(/0.5,0.5,0.0/)
                positions(n+2,:)=positions(n,:)+a*(/0.5,0.0,0.5/)
                positions(n+3,:)=positions(n,:)+a*(/0.0,0.5,0.5/)
                n=n+4
            END DO
        END DO
    END DO
    PRINT*, 'particles positioned', n-1, 'of a total imput', n_particles
    RETURN
END SUBROUTINE FCC_INITIALIZE
SUBROUTINE UNIFORM_VELO(velocity,T)
    IMPLICIT NONE
    REAL*8 :: density,L,a
    INTEGER :: n_particles,M,i,j,seed
    REAL*8 :: velocity(n_particles,3),vi,vtot,T
    COMMON/PARAMETERS/n_particles,M,density,L,a
    seed=13
    CALL SRAND(seed)
    DO i=1,3
        vtot=0
        DO j=1,n_particles-1
            vi=2*RAND()-1
            velocity(j,i)=vi
            vtot=vtot+vi
        END DO
        velocity(n_particles,i)=-vtot
    END DO
    !Resacling the velocities to the temperature
    CALL VELO_RESCALING(velocity,T)
    RETURN
END SUBROUTINE UNIFORM_VELO
SUBROUTINE VELO_RESCALING(velocity,T)
    IMPLICIT NONE
    INTEGER n_particles,M
    REAL*8 :: density,L,a
    REAL*8 velocity(n_particles,3),T,KINETIK,alpha
    COMMON/PARAMETERS/n_particles,M,density,L,a
    alpha=sqrt(3d0*n_particles*T/(2d0*KINETIK(velocity)))
    velocity=alpha*velocity
END SUBROUTINE VELO_RESCALING
FUNCTION KINETIK(velocity)
    IMPLICIT NONE
    INTEGER n_particles,M,i
    REAL*8 :: density,L,a
    REAL*8 :: velocity(n_particles,3),KINETIK
    COMMON/PARAMETERS/n_particles,M,density,L,a
    KINETIK=0d0
    DO i=1,n_particles
        KINETIK=KINETIK+5d-1*(velocity(i,1)**2d0+velocity(i,2)**2d0+velocity(i,3)**2d0)
    END DO
    RETURN
END FUNCTION KINETIK
SUBROUTINE RAD_DIST_INTER(r,vec,dx_radial,n_radial,n_particles,L)
    IMPLICIT NONE
    INTEGER n_radial,n_particles,i,coef,j
    REAL*8 dx_radial,dist,vec(0:n_radial+1), r(n_particles,3),dx,dy,dz,PBC1,L
    DO i=1,n_particles
        DO j=1,n_particles
            IF (i.ne.j) THEN
                dx=PBC1(r(j,1)-r(i,1),L)
                dy=PBC1(r(j,2)-r(i,2),L)
                dz=PBC1(r(j,3)-r(i,3),L)
                dist=(dx**2d0+dy**2d0+dz**2d0)**0.5
                coef=int(0.5+dist/dx_radial)
                IF ((coef.gt.0).and.(coef.le.n_radial)) THEN
                    vec(coef)=vec(coef)+1.0
                END IF
            END IF
        END DO
    END DO
    RETURN
END
SUBROUTINE RAD_DIST_FINAL(vec,dx_radial,n_radial,n_gr_meas,rho)
    IMPLICIT NONE
    INTEGER n_radial,i,n_gr_meas
    REAL*8 dx_radial,rho,vec(0:n_radial+1),result(0:n_radial+1),aux
    vec=vec/(1d0*n_gr_meas)
    DO i=2,n_radial
        aux=(rho*4d0*3.1415*((((i)*dx_radial)**3d0)-(((i-1)*dx_radial)**3d0)))/3d0
        result(i)=vec(i)/aux
    END DO
    vec=result 
    RETURN
END
SUBROUTINE ANDERSEN(v,temp,dt,n_particles)
    IMPLICIT NONE
    INTEGER n_particles,i
    REAL*8 temp,dt,nu,sigma,n1,n2,n3,n4,n5,n6
    REAL*8, DIMENSION(n_particles,3) :: v
    nu=0.1/dt
    sigma=sqrt(temp)
    DO i=1,n_particles
        IF(RAND().lt.nu*dt)THEN
            n1=RAND();n2=RAND()
            n3=RAND();n4=RAND()
            n5=RAND();n6=RAND()
            v(n_particles,:)=(/sigma*sqrt(-2d0*log10(n1))*cos(2d0*3.1415*n2),sigma*sqrt(-2d0*log10(n1))*sin(2d0*3.1415*n2),&
                &sigma*sqrt(-2d0*log10(n3))*cos(2d0*3.1415*n4)/)
        END IF
    END DO
    RETURN
END

