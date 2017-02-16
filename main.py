from particle import *
from kernels import W_4, gW_4, lW_v
from IISPH import IISPH_Algorithm
import particleInit



######### Initializes the "particles" #########
blockOfFluid = particleInit.ParticleInitData(

	systemConstants   =  {"viscosity"       : 0.3,
	                      "rho0"            : 1000,
	                      "interactionlen"  : 0.51,
	                      # "diameter"        : 0.5128,
	                      # "diameter"        : 0.5003,
	                      "diameter"        : 0.25,
	                      "gravity"         : 10,
	                      "dt"              : 0.005,
	                      "domain"          : (30,30),
	                      "walls"           : (0.5,0.5),
	                     },

	particleVariables =  {"rho"             : 0,
	 					  "pressure"        : 10,
	 					  "pressure_est"    : 0,
	 					  "temperature"     : 0,
	 					  "mass"            : 0.0625 * 1000,
	 					  "f_adv"           : Vec2(0,0),
	 					  "v_adv"           : Vec2(0,0),
	 					  "p_est"           : Vec2(0,0),
	 					  "v_est"           : Vec2(0,0),
	 					  "rho_est"         : 0,
	 					  "rho_adv"         : 0,
	 					  "d_ii"            : Vec2(0,0),
	 					  "sum_p_j_dij"     : Vec2(0,0),
	 					  "fp_dt**2/m"      : Vec2(0,0)}
)	 				

IISPH = IISPH_Algorithm(W_4,gW_4,lW_v,omega = 0.5)


psys = ParticleSystem(
	interactionAlgo      = IISPH,
    particleInitData     = blockOfFluid)
print("Starting Simulation.")
psys.run()
