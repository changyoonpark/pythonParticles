from particle import *
import particleInit
import collisionModel

stupidBlock = particleInit.ParticleInitData(viscosity = 0.1,
	 				  					    particleMass = 0.1,
	 				  		                interactionLength = 0.5)

psys = ParticleSystem(collisionModel = collisionModel.stupidSpring,
	 				  particleInitData = stupidBlock,
	 				  dt = 0.015)

psys.run()
