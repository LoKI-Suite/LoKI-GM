<p align="center">----------------------------- README FILE FOR THE LoKI-GM SIMULATION TOOL -----------------------------<br>
<align="center">-----------------------------  10 steps to get acquainted with the tool  -----------------------------<br>
<align="center">(updated for version LoKI-GM_26.07)</p>

1. What's LoKI-GM ?   
   The LisbOn KInetics Global Model (LoKI-GM) is a simulation tool, intended to model non-equilibrium low-temperature plasmas (LTPs), produced from different gases or gas mixtures under a wide range of working conditions.
   Overall, LoKI-GM provides the combined chemical and transport description of plasma charged/neutral species, both in volume and surface phases, for user-defined working conditions: mixture compositions, pressure, dimensions and excitation features (glow (DC/HF), post-discharge, pulsed).

   LoKI-GM comprises two modules, that can run self-consistently coupled or as standalone tools:
   - a Boltzmann solver (LoKI-B) for the electron Boltzmann equation;
   - a Chemical solver (LoKI-C) for the global kinetic model(s) of pure gases or gas mixtures.
   
   The relevant documentation can be found [here](https://github.com/LoKI-Suite/LoKI-GM_RC/tree/master/Documentation).

2. What's the programming language of LoKI-GM ?   
   LoKI-GM is developed with flexible and upgradable object-oriented programming under MATLAB, leveraging its matrix-based architecture.   
   The tool adopts an ontology that privileges the separation between tool and data.

3. What are the input data of LoKI-GM ?  
   The general input to LoKI-GM specifies the working conditions, including the applied reduced electric field and frequency, gas pressure and temperature, electron density (or a related quantity such as discharge current or power density), and, in the case of gas mixtures, the fractions of the various components. 
   
   For LoKI-B, the input data additionally comprise the population distributions of the electronic, vibrational, and rotational states of the atomic/molecular species considered (determined self-consistently when coupled with LoKI-C), as well as the relevant sets of electron-scattering cross sections obtained from the open-access database LXCat (http://www.lxcat.net/). 
   
   For LoKI-C, the required inputs further include the geometrical dimensions of the plasma reactor, gas inflow/outflow conditions, wall/external temperatures relevant to the gas/plasma thermal model, the surface site density when a microkinetic mesoscopic surface model is employed, and the rate coefficients for all reactions considered in the kinetic scheme, including the electron rate coeffcients evaluated in LoKI-B.

4. What are the output results of LoKI-GM ?   
   As output, LoKI-B provides the isotropic and anisotropic components of the electron distribution function (the former usually termed the electron energy distribution function, EEDF), the corresponding electron macroscopic parameters (transport parameters and rate coefficients), and the electron power absorbed from the electric field and transferred to the various collisional channels. The latter parameters can be calculated using either the distribution function obtained from the solution to the EBE or some other form prescribed by the user, e.g. a generalized Maxwellian EEDF.   
   LoKI-C yields the densities of states and the reaction rates associated with the processes included in the kinetic scheme.
   
5. How to find your way in the code ?   
   After pulling the files in the repository, the LoKI-GM folder contains   

   A) Subfolder "Documentation", containing the user guide, the reference manual and relevant papers   
   PLEASE READ THEM BEFORE USING THE CODE !!!   

   B) Subfolder "Code" containing   
   &ensp;(a) Several '\*.m' files corresponding to the MATLAB code, of which 'loki.m' is the main file of LoKI-GM.   
   &ensp;(b) Two files to launch the GUI-based input for LoKI-B: the MATLAB file 'InputGUILokib.m', which provides a user-friendly graphical user interface for LoKI-B setup and run; and the batch file 'InputGUILokib.bat', which allows launching the GUI-based input from a cmd window (or by double-clicking the corresponding '\*.bat' file).   
   &ensp;(c) A subfolder "Input", containing the input files required for the simulations, organised as follows   
   &ensp;&ensp;i. Default configuration files in text format 'default_lokibc_setup.in', 'default_lokic_pulse_setup.in', 'default_lokib_setup.in' and 'default_lokib_pulse_setup.in', respectively, for LoKI-GM (the first two files, for steady-state simulations and quasi-stationary simulations, respectively) and for LoKI-B as standalone tool (the latter two, for electron kinetics stationary simulations and electron kinetics time-dependent simulations, respectively).   
   &ensp;&ensp;ii. Default configuration files in JSON format 'default_lokib_setup.json' and 'default_lokib_pulse_setup.json', for LoKI-B as standalone tool (for electron kinetics stationary simulations and electron kinetics time-dependent simulations, respectively).   
   &ensp;&ensp;iii. A subfolder "Databases" with '\*.txt' files, containing different properties (masses, energies of states, atomic/molecular constants, ...) for the gases used in the simulations.   
   &ensp;&ensp;iv. Several subfolders "Argon", "CO", ... "Nitrogen", "Oxygen" with '\*.txt' files, containing:  
   &ensp;&ensp;&ensp; - the data needed to perform simulations in different gases or gaseous mixtures: electron-scattering cross sections usually obtained from the open-access website LXCat (https://lxcat.net/), kinetic schemes, etc;   
   &ensp;&ensp;&ensp; - setup input files to run swarm simulations;  
   &ensp;&ensp;&ensp; - setup input files to run chemistry simulations.  
   &ensp;&ensp; Subfolders "Argon", "CO", "CO2", Nitrogen", "Oxygen" contain also '\*.json' files with electron-scattering cross sections (obtained from the demo version of LXCat 3.0, https://demo.lxcat.net/) and input setups for swarm simulations.  
   &ensp;(d) A subfolder "PropertyFunctions", with several '\*.m' auxiliary functions for calculating some predefined distribution of states (Boltzmann, Treanor, ...), the statistical weights of states due to their degeneracy, the energy of states according to some models, transport parameters and thermal conductivities of different species, etc.  
   &ensp;(e) A subfolder "RateCoeffFunctions", organized into one-level subfolders with self-explanatory names, with several '\*.m' auxiliary functions for calculating the values of the different rate coefficients used in LoKI-C.   
   &ensp;(f) A subfolder "OtherAuxFunctions", with several '\*.m' auxiliary functions, e.g. to calculate some working conditions.   
   &ensp;(g) A subfolder "Utilities", with several scripts to help process or parse data, etc.   
   &ensp;(h) A subfolder "Output" (eventually), where LoKI-GM will write the output files resulting from the simulations.

6. How to run LoKI-B ?   
   The minimum requirements to run the code is a computer with an installation of MATLAB (oldest recommended version R2024b; we cannot ensure that all the features of LoKI-GM will work properly under different versions).

   LoKI-GM runs upon calling the MATLAB function 'loki(setupFile)'.
   The end user interacts with the code by specifying a particular "setup" for the simulation.    
   This setup is sent to the 'loki()' function through the required input argument 'setupFile'.   
   The setup files should be located in [repository folder]/LoKI-GM_26.07/Code/Input/ with  
   &ensp;&ensp;- '.in' extension (this is just a recommendation in order to keep the input folder organised; the '.in' setup files are just plain text files);  
   &ensp;&ensp;- '.json' extension (in this case it is important to use this extension, for the parsing to work properly).

   The distribution of LoKI-GM includes some default configuration files for the benefit of the user. By using one of these files, it is very easy to make a first run of the code, just following the sequence of steps below.
   
   For example, using the 'default_lokibc_setup.in' setup file:   
   A) Open MATLAB   

   B) Navigate to the "Code" folder of your local copy of the repository: 
   ```bash  
   >> cd [repository folder]/LoKI-GM_26.07/Code/
   ```   
   C) Execute the following command in the MATLAB command window:
   ```bash
   >> loki('default_lokibc_setup.in')
   ```   
   D) The graphical user interface (GUI) should appear showing the solution(s) for the default setup file.


   To run LoKI-B only, users have different alternatives:   
   C') Using the MATLAB command window with one of the following instructions:
   ```bash
   >> loki('default_lokib_setup.in')
   ``` 
   ```bash
   >> lokibcl('default_lokib_setup.in')
   ``` 
   C'') Through a GUI-based input by using one of the following options:
   - from the MATLAB command window:
   ```bash
   >> InputGUILokib
   ``` 
   - from a cmd window (or by double-clicking the corresponding '\*.bat' file)
   ```bash
   >> InputGUILokib.bat
   ``` 

7. What's LoKI-GM distribution license ?   
   LoKI-GM is an open-source tool, licensed under the GNU general public license.  
   LoKI-GM is freely available for users to perform electron kinetics calculations, and for expert researchers who are invited to continue testing the tool and/or to contribute for its development and improvement.

8. How to contact the developers ?   
   You are much welcome to report any problem or bug you may find using the code, or to send any feedback with suggestions for further developments.

   After downloading LoKI-GM, and especially if you intend to interact with us, you are invited to send a short message   
   to: loki@tecnico.ulisboa.pt   
   with subject: <b>LoKI-GM</b>   
   just giving your <b>name</b> and <b>affiliation</b>.   

9. How to reference the code ?   
   LoKI-GM is the result of the efforts of the Portuguese group N-Plasmas Reactive: Modeling and Engineering (N-PRiME) of Instituto de Plasmas e Fusão Nuclear with Instituto Superior Técnico, and the Departamento de Física of Facultad de Ciencias with Universidad de Córdoba. These groups decided to share the outcome of its research on code development with the members of the Low-Temperature Plasmas community.

   When using the code in your work, please give proper credits to the main developers, by adding the following citations:     

   For LoKI-B   
   [] Tejero A et al "The LisbOn KInetics Boltzmann solver" 2019 Plasma Sources Sci. Technol. 28 043001 (https://doi.org/10.1088/1361-6595/ab0537)       
   [] Tejero A et al "On the quasi-stationary approach to solve the electron Boltzmann equation in pulsed plasmas" 2021 Plasma Sources Sci. Technol. 30 065008 (https://doi.org/10.1088/1361-6595/abf858)        
   [available as open-access papers]

   For LoKI-GM   
   [] Alves L L et al "LoKI-GM: a global model tool for plasma chemistry studies" 2026 Plasma Sources Sci. Technol. (in preparation)

10. Acknowledgments   
This work was funded by Portuguese FCT - Fundação para a Ciência e a Tecnologia, initially under project PTDC/FISPLA/1243/2014 (KIT-PLASMEBA), and currently under projects [PSI.COM](https://doi.org/10.54499/2022.04128.PTDC), [UID/50010/2025](https://doi.org/10.54499/UID/50010/2025), [UID/PRR/50010/2025](https://doi.org/10.54499/UID/PRR/50010/2025), [UID/PRR2/50010/2025](https://doi.org/10.54499/UID/PRR2/50010/2025), and [LA/P/0061/2020](https://doi.org/10.54499/LA/P/0061/2020). 
