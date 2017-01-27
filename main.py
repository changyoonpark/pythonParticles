from particle import *
import particleInit
import collisionModel



######### Initializes the "particles" #########
blockOfFluid = particleInit.ParticleInitData(

	systemConstants   =  {"viscosity"       : 0.1,
	                      "rho0"     : 1000,
	                      "kappa"           : 10,
	                      "gamma"           : 3,
	                      "interactionlen"  : 0.5,
	                      "dt"              : 0.015},

	particleVariables =  {"rho"         : 1000,
	 					  "pressure"        : 0,
	 					  "temperature"     : 0,
	 					  "mass"            : 0.1,
	 					  "f_adv"           : Vec2(0,0),
	 					  "v_adv"           : Vec2(0,0)}

)	 				


operationFuncs = [collisionModel.stupidSpring]

psys = ParticleSystem(
	operationFuncs       = operationFuncs,
    particleInitData     = blockOfFluid)

psys.run()
