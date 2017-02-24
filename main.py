from particle import *
from kernels import W_4, gW_4, lW_v
from IISPH import IISPH_Algorithm
import particleInit


sysConsts = { "viscosity"       : 0.0 ,
              "rho0"            : 1000,
              "interactionlen"  : 0.52,
              "diameter"        : 0.25,
              "gravity"         : 0,
              "dt"              : 0.0001,
              "domain"          : (10,10),
              "walls"           : (0.5,0.5),
	        }

pVars = {     "mass"            : 0.0625 * 1000,
              "rho"             : 1000,
			  "pressure"        : 0,
			  "pressure_est"    : 0}

######### Initializes the "particles" #########
blockOfFluid = particleInit.ParticleInitData(
	systemConstants   =  sysConsts,
	particleVariables =  pVars
)	 				

simpleWall = particleInit.BoundaryInitData(
	systemConstants   =  sysConsts
)

IISPH = IISPH_Algorithm(W_4,gW_4,lW_v,omega = 0.5)

psys = ParticleSystem(
	interactionAlgo      = IISPH,
    particleInitData     = blockOfFluid,
    boundaryInitData     = simpleWall)
print("Starting Simulation.")
psys.run()      
