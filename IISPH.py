from kernels import W_4, gradW_4


# def operation(systemConstants,pairsData,particle):
def predictAdvection(systemConstants,pairsData,particle):
#   calculate density
	particle.particleVariables["rho"] = 0
	for neighbor in particle.neighborList:
		pairData = pairsData[(particle,neighbor)]
		particle.particleVariables["rho"] += particle.particleVariables["mass"] * W_4(pairData,systemConstants["interactionlen"])
		
#	calculate advection force
		


	particle.particleVariables["v_adv"] = 
	particle.vel + systemConstants["dt"] * particle.particleVariables["f_adv"] / particle.particleVariables["mass"]
	