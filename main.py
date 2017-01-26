from particle import *
import particleInit
import collisionModel


interactionLength = 0.5


######### Initializes the "particles" #########
blockOfFluid = particleInit.ParticleInitData(

	systemConstants   =  {"viscosity"     : 0.1,
	                      "restdensity"   : 1000,
	                      "kappa"         : 10},

	particleVariables =  {"density"       : 1000,
	 					  "pressure"      : 0,
	 					  "temperature"   : 0,
	 					  "mass"          : 0.1},

    interactionLength = interactionLength)

######### Iterate through the particles and do these #########
beforeCollisionFunctions = []
afterCollisionFunctions  = [] 

psys = ParticleSystem(
	collisionModel       = collisionModel.stupidSpring,
    beforeCollisionFuncs = beforeCollisionFunctions,
    afterCollisionFuncs  = afterCollisionFunctions,
    particleInitData     = blockOfFluid,
    dt                   = 0.015)

psys.run()
