Notes on c-goldstein version 1:

What does goldstein stand for? ans. Global Ocean-Linear Drag Salt and
Temperature Equation Integrator.

General notes on running the code

Compilation and running

The command: 
make goldstein
produces an executable file called goldstein. To run the program from the
input file goin, use the command
goldstein < goin
You may need to edit the FLAGS= line in the Makefile to specify appropriate
command line options to the compiler. It is strongly recommended to always
specify an option such as -r8 or autodblpad which upgrades all reals to
real*8. 
You must create a directory called results as the executable file 
writes several output
files to the location ../results/lout.* where lout is a 3 character name for
the run. The values of * and when they appear for given values of npstp etc
are as follows:

Output

At multiples of iwstp the subroutine outm writes the following files:

lout.n     these are 'restart' files which contain ocean variables T,S
           u and v, atm. variables T and Q and ice variables h, A and T.
           and can subsequently
           be used as input files in continuing runs, n=1,2,3,...,0, the 11th
           will overwrite the first etc
           Output is in ascii and can be plotted or even read by eye line by 
           line but accuracy is sufficient to allow exact restarts with error
           at only machine error level.
lout.psi.n if there is topographic variation or oscillatory forcing the
           barotropic streamfunction is written 

At multiples of itstp the subroutine diag2 writes:

lout.t     averaged values of ocean T in various regions in 'time series' form
           for the whole run.
lout.s     equivalent averaged salinities 
lout.osw   or similar, time series of data with,in addition, one or more spatial 
           dimensions, see mains.F
lout.airt  averaged values of air temperature
lout.q	   averaged values of specific humidity

The program stops after nsteps steps and main program 'mains' writes:

lout.opsi  meridional overturning streamfunction
lout.zpsi  zonal      overturning streamfunction
lout.psi   barotropic streamfunction
lout.rho   density
lout.flux  surface heat flux
lout.cost  a 2d array showing the frequency of convection during the run

plus the following dimensionalized data:

2d surface fields of:

lout.swflux     absorbed incoming solar radiation
lout.sens       sensible heat flux into ocean
lout.lato       latent heat flux into ocean
lout.fluxo	net surface heat flux over ocean
lout.lwflux	re-radiated longwave heat flux
lout.plwflux	outgoing planetary longwave radiation 
lout.relh	relative humidity
lout.lata	latent heat flux into atmosphere

lout.fluxo	net surface heat flux over ocean
lout.fluxa	net heat flux into atmosphere
lout.fluxt	total surface heat flux from atmosphere to ocean/sea ice 	

lout.pptn       precipitation rate
lout.evap       evaporation rate
lout.runoff     river runoff
lout.pminse     precipitation minus evaporation plus runoff
lout.fwflux     net freshwater flux into ocean (P-E+R+freeze/melt)

lout.fluxi	net surface heat flux over sea ice

At multiples of npstp the subroutine diag is called which writes some 
information to unit 6, normally the screen. diag is called again at the
next step so you can spot any oscillatory instability easily. The second
time, ie one step after multiples of npstp, mains calculates an average
rate of change of all dynamic variables (if npstp .ge. 2).

Dimensions / units 

The ocean module uses dimensionless variables although temp and salinity
can be considered to be in Celcius and psu resp. The atmospheric and 
sea-ice transport routines are also dimensionless but the coupling
routine surflux works exclusively with dimensional variables in SI units.
The scales for most dimensionless quantities
are calculated and printed by the routine gseto using some reasonable
choices for scale values. The forcing amplitudes and functions are set
in gseto.f and gseta.f which do all the 'set-up' for ocean and 
atmosphere/sea ice respectively. Use the scales in gset to non-
dimensionalize any forcings you put in, and to interpret the results.

Plotting

Selected fields can be plotted using the UNIRAS routines in ../plotting,
routine names self-explanatory.

Other notes

Model subroutines additional to those in the ocean-only version:

gseto.f (in place of gset.f): initializes ocean
gseta.f: initializes atmosphere

tstepo.f (in place of tste2.f): ocean timestep (updates 3-d fields of
				ocean temperature and salinity)

surflux.f: computes heat and freshwater fluxes at
	   atmosphere-ocean/seaice/land interface
tstepa.f: atmospheric timestep (updates 2-d fields of surface  air
	  temperature and specific humidity)

thermo_seaice.f: updates sea ice thickness, fractional area due to
		 sea ice thermodynamics (every atmospheric timestep)
transp_seaice.f: updates sea ice thickness, fractional area due to
		 free drift with surface ocean current (every ocean ts)

Input data (in /code):

read in gseto.f:
taux_annav ... zonal component of wind stress
tauy_annav ... meridional component of wind stress
etopo5.dat ... ETOPO5 modified topography
worbob2.k1 ... land-sea mask

read in gseta.f:
albedo.dat ... climatological albedo
u_ncep.dat ... zonal component of surface wind 
v_ncep.dat ... meridional component of surface wind
diffatmt_DIF.dat OR diffatmt_ADVDIF.dat ... atmos. diffusivity of heat
diffatmq.dat ... atmos. diffusivity of moisture

References to ocean-only versions of the code:

Edwards, N.R. 1996. Unsteady similarity solutions and oscillating
ocean gyres. J. Marine Res., Vol. 54, pp. 793-826.

Edwards, N.R., Willmott, A.J. and Killworth, P.D. 1998. On the role of
topography and wind stress on the stability of the thermohaline circulation.
J. Physical Oceanography, Vol. 28, pp. 756-778.

Edwards, N.R. and Shepherd, J.G. 2001. Multiple thermohaline states due
to variable diffusivity in a hierarchy
of simple models. Ocean Modelling, Vol. 3, pp. 67-94.

Reference to global, ocean-only version:

Edwards, N.R. and Shepherd, J.G. 2002. Bifurcations of the thermohaline
circulation in a simplified three-dimensional model of
the world ocean and the effects of interbasin connectivity.
Climate Dynamics, Vol. 19, pp. 31-42.

The EMBM and sea ice code is similar to that outlined in:

Weaver et al. (2001): The UVic Earth System Climate Model: Model Description,
Climatology, and Applications to Past, Present and Future Climates.
Atmosphere-Ocean 39(4), 361-428.

