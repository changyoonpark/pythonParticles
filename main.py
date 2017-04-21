from particle import *
from kernels import W_4, gW_4, lW_v
from IISPH import IISPH_Algorithm
import particleInit


sysConsts = { "viscosity"       : 0.01,
              "rho0"            : 1000,
              "interactionlen"  : 0.55,
              "diameter"        : 0.25,
              "gravity"         : 10,
              "dt"              : 0.01,
              "domain"          : (24,12),
              "walls"           : (0.5,0.5),
              "densityDeviation": [[0,1]],
              "gamma"           : 0.5,
	        }

pVars = {     
        "rho"             : 1000,
			  "pressure"        : 0,
			  "pressure_est"    : 0}

######### Initializes the "particles" #########
blockOfFluid = particleInit.ParticleInitData(
	systemConstants   =  sysConsts,
	particleVariables =  pVars
)	 				

# diskOfFluid = particleInit.ParticleDiskInitData(
# 	systemConstants = sysConsts,
# 	particleVariables = pVars
# )

simpleWall = particleInit.BoundaryInitData(
	systemConstants   =  sysConsts
)

IISPH = IISPH_Algorithm(W_4,gW_4,lW_v,omega = 0.4)

psys = ParticleSystem(
	  interactionAlgo      = IISPH,
    particleInitData     = blockOfFluid,
    boundaryInitData     = simpleWall)
print("Starting Simulation.")
psys.run()      

