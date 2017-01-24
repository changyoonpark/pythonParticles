from particle import *
from collisionModel import stupidSpring
psys = ParticleSystem(collisionModel = stupidSpring)
psys.run()
