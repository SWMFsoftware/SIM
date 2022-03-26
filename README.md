Copyright (C) 2002 Regents of the University of Michigan,
portions used with permission.
For more information, see http://csem.engin.umich.edu/tools/swmf

# SIM - Simple Ionosphere Model

SIM is an ionosphere electrodynamics solver being designed to use 
with non-dipolar or multipolar planetary magnetic fields. It works
stand-alone, and we are working to couple it with BATSRUS. 

- It solves for the ionospheric potential over the whole sphere 
instead of two separate hemispheres. This also means that linear 
system created for the Krylov solver is larger. As SIM is a serial
code, it is usually slower than Ridley_serial. 

- SIM allows for inhomogeneous Pedersen and Hall conductances.





