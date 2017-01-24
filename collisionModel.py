from helpers import *

def stupidSpring(pairData,interactionLength):
	dist = pairData.dist
	relvel = pairData.relvel
	reldir = pairData.reldir
	particlei = pairData.particlei
	particlei.fext = particlei.fext - 300 * pow(dist - interactionLength,1) * reldir
	particlei.fext = particlei.fext - 0.5 * relvel
