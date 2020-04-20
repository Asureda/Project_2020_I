# Molecular Dynamics Symulations

Computational progect to develop sequenctial and paralelized programs to symulate molecular dynamics. 

## Fisrst steps 🚀
Information to install and execute the programs


### Pre-requisites 📋

Sequential compilers:

```
ifort (Default)
gfortran (Must configure Makefile)
```

Parallel compilers:
```
intel openmpi (Default)
```

### Instalation 🔧

The programs are ready-to-use. The user have to download the repository in a local computer folder or computing cluster, configure the compiller and flags options in the Makefile

Sequential program
```
Makefile:  configure compiller and flags variables (ifort by default)

```
Paralel program (compututing cluster)
```
Makefile:  configure the compiler anf flags variables (mpifort by default)
"run_sub.sh" (1): Check the execution order ( mpirun by default)
"run_sub.sh" (2): Configure the submit options ( BSC by default)
```

## Execution ⚙️

Sequential program
```
(1) Configure the symulation parameters (INPUT folder)
(2) Execute the "run.sh" script.
(3) Collect the results in the OUTPUT folder
    The results folder name is the date-time whem the task was executed

```
Paralel program (compututing cluster)
```
(1) Configure the symulation parameters (INPUT folder)
(2) Execute the "run_sub.sh" script.
(3) Collect the results in the OUTPUT folder
    The results folder name is the date-time whem the task was executed
```
### Program-check 🔩

In the OUTPUT foler is provided a run_check subfolder with input configuration parameters and graphs. 
Put the same paremeters in the INPUT files, run the program and compare the graphs. 
They should be similar except for a random factor

### About the different programs ⌨️

Brief description of the main characteristics

MDS
```
blabla
```
MDP1
```
blabla
```
MDP2
```
blabla
```
MDP3
```
blabla
```
MDP4
```
blabla
```

## Despliegue 📦

_Agrega notas adicionales sobre como hacer deploy_

## Construido con 🛠️

_Menciona las herramientas que utilizaste para crear tu proyecto_

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - El framework web usado
* [Maven](https://maven.apache.org/) - Manejador de dependencias
* [ROME](https://rometools.github.io/rome/) - Usado para generar RSS

## Contribuyendo 🖇️

Por favor lee el [CONTRIBUTING.md](https://gist.github.com/villanuevand/xxxxxx) para detalles de nuestro código de conducta, y el proceso para enviarnos pull requests.

## Wiki 📖

Puedes encontrar mucho más de cómo utilizar este proyecto en nuestra [Wiki](https://github.com/tu/proyecto/wiki)

## Versionado 📌

Usamos [SemVer](http://semver.org/) para el versionado. Para todas las versiones disponibles, mira los [tags en este repositorio](https://github.com/tu/proyecto/tags).

## Autores ✒️

_Menciona a todos aquellos que ayudaron a levantar el proyecto desde sus inicios_

* **Andrés Villanueva** - *Trabajo Inicial* - [villanuevand](https://github.com/villanuevand)
* **Fulanito Detal** - *Documentación* - [fulanitodetal](#fulanito-de-tal)

También puedes mirar la lista de todos los [contribuyentes](https://github.com/your/project/contributors) quíenes han participado en este proyecto. 

## Licencia 📄

Este proyecto está bajo la Licencia (Tu Licencia) - mira el archivo [LICENSE.md](LICENSE.md) para detalles

## Expresiones de Gratitud 🎁

* Comenta a otros sobre este proyecto 📢
* Invita una cerveza 🍺 o un café ☕ a alguien del equipo. 
* Da las gracias públicamente 🤓.
* etc.



---
⌨️ con ❤️ por [Villanuevand](https://github.com/Villanuevand) 😊
