from helpers import *

def stupidSpring(systemConstants,pairsData,particle):
	
	for neighborParticle in particle.neighborList:
		pairData = pairsData.get((particle,neighborParticle))
		if pairData is not None :
			il = systemConstants["interactionlen"]
			particle.fext = particle.fext - 300 * \
				pow(pairData.dist - il,1) * pairData.reldir
			particle.fext = particle.fext - 0.5 * pairData.relvel
