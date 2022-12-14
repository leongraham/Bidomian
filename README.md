# Bidomian
2D Bidomain model for simulating the electrical activity of heart tissue. FitzHugh-Nagumo membrane and Beeler-Reuter membrane models are implemented.
Python code using escript from the University of Queensland for solving the PDEs. Extracellular and transmembrane potential are calculated. Data output are saved as VTK (.vtu) files, which can be visualised in Paraview or VisIt.

Setup:

#####Conductivity parameters##### 

Values (mS/cm) used from Roberts et  al. (1979):

sigma_int_l =  2.8

sigma_int_t =  0.26

sigma_ext_l =  2.2

sigma_ext_t =  1.3

#####Domain#####

1.0x1.0 cm tissue region with 100x100 nodes.

#####Stimulation#####

100 micro amps stimulation in centre of domain for a duration of 1 ms.

#####Membrane capacitance and surface-to-volume ratio#####

Cm=1.0 uF/cm^2

beta=2000 cm^-1

#####Time-step size#####

For FitzHugh-Nagumo membrane, typically dt=0.01 ms

For Beeler-Reuter membrane model, dt=0.001 ms for stability.

#####Fibre angle#####

This can be set with the parameter theta (degrees).

