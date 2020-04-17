!GRUP I: Àlex, Oriol, Laia, Sílvia i Elena

PROGRAM SEQUENTIAL_MD
  !Cridem els mòduls que necessitarem:
  use READ_DATA
  use ALLOCATE_VARS
  use Inicialitzar
  use Distribucio_Uniforme_vel
  use Interaction_Cutoff_Modul
  use Verlet_Algorithm
  use Andersen_modul
  use Distribucio_Radial
  use Reescala_velocitats
  use SAMPLE
  
  !Declarem les variables següents:
  IMPLICIT NONE
  
  !Cridem la seed per generar els números aleatoris:
  call srand(seed)
  
  !Inicilitzem el sistema de N partícules a temps = 0
  call read_all_data()
  call other_global_vars()
  call INITIALIZE_VARS()
  
  !Definim posicions inicials (definides amb una estructura FCC) i velocitats:
  call FCC_Initialize(r)
  call Uniform_velocity(v)
  
  !Fem un reescalat de les velocitats un cop generades:
  call VELO_RESCALING_MOD(v,T_therm_prov)

  !Fem el melting, per tal de trencar l'estructura FCC creada inicialment, i aconseguir un fluid
  ! amb l'objectiu d'analitzar les seves característiques termodinàmiques i cinètiques.
   DO i=1,n_melting
    call velo_verlet(r,v,F)
     call andersen(v,T_therm_prov)
   END DO

  !Ara, comencem a guardar els resultats en els diferents arxius:
  !Un cop acabat el melting, ja podem printejar la dinàmica del sistema
  !Es guardarà aquella informació que es trobi al worker0 (master):
  open(51,file='thermodynamics_reduced.dat')
  open(52,file='thermodynamics_real.dat')
  open(53,file='distrib_funct.dat')
  open(54,file='positions.xyz')

  !Apliquem l'algorisme de Verlet, que serveix per integrar la quantitat de moviment
  ! amb el termòstat d'Andersen que ens permet aconseguir treballar a una temperatura
  ! constant a partir de les velocitats:
  DO i=1,n_verlet
    t=t_a+i*h
    call VELO_VERLET(r,v,F)
    if(is_thermostat.eqv..true.)THEN
      call andersen(v,T_therm)
    end if
  
  !Cridem el mòdul que ens permet escriure per pantalla els valors de:
  ! la distribució radial, l'energia cinètica, potencial i total, la pressió
  ! i la temperatura instantània; tant en unitats reals com reduïdes:
  CALL SAMPLES()
  end do
  CALL gdr()
  print*,'PROGRAM END'
END PROGRAM SEQUENTIAL_MD
