# Orbit-Insertion

Script that simulates the flight of a 3-stage rocket for orbit insertion.

* ```voo_insercao_orbita.m```: Main function. Declares structural and payload masses, propulsive parameters, aerodynamic and ambiental parameters, Earth model patrameters, and initial conditions. Uses the Rocket Equation to compute the required Delta V. Simulates and draws the trajectory.
* ```aerodinamica_N_estagios.m```: Computes the drag force for each stage of a rocket of up to 3 stages.
* ```atm_padrao.m```: Computes the standard atmosphere data for geometrical altitudes from 0 to 2000km.
* ```desenha_mapa_trajetoria.m```: Draws a projection map along the provided trajectory.
* ```det_orb.m```: Determines orbital parameters from position and velocity observations.
* ```dinamica_foguete.m```: Computes the translation dynamics of a rocket considering the PCPF reference frame.
* ```earth_sphere.m```: Generate an earth-sized sphere.
* ```grav_axissimetrico.m```: Computes force of gravity for an axissimetrical body.
* ```long_ECEF2ECL.m```: Computes the celestial longitude for a given planet-fixed longitude
* ```modelo_aerodinamico.m```: Calculates the drag coefficient.
* ```propulsao_N_estagios.m```: Computes the propulsive parameters as a function of time.
* ```RvelPolar2RvelRet.m```: Converts velocity from the LVLH coordinate system to the ECI or ECEF coordinate systems.
* ```Vrel2Vine.m```: Converts velocity, elevation and azimuth values from relative velocity to inertial velocity.

![untitled](https://github.com/rewenila/orbit-insertion/assets/32466162/db315587-0cb0-42ab-9faa-50660d0246ab)
![untitled2](https://github.com/rewenila/orbit-insertion/assets/32466162/88b5584e-2b1f-41b9-985f-cc75dd9de528)

