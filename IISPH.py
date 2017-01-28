from kernels import W_4, gW_4, lW_v


# def operation(systemConstants,pairsData,particle):
def predictAdvection(systemConstants,pairsData,particle):

#   calculate density
	particle.particleVariables["rho"] = 0
	for neighbor in particle.neighborList:
		pairData = pairsData[(particle,neighbor)]
		particle.particleVariables["rho"] += 
			particle.particleVariables["mass"] * W_4(pairData,systemConstants["interactionlen"])

#	calculate advection force (viscosity and gravity is applied.)
	for neighbor in particle.neighborList:
		pairData = pairsData[(particle,neighbor)]
		particle.particleVariables["f_adv"] += \
			systemConstants["viscosity"] * \
			lW_v(pairData,systemConstants["interactionlen"]) * \
			neighbor.particleVariables["mass"] / \
			neighbor.particleVariables["rho"] * \
			pairData.relvel

#   consider gravity.
	particle.particleVariables["f_adv"] -= \
		particle.particleVariables["mass"] * systemConstants["gravity"]

#   calculate itermediate velocity
	particle.particleVariables["v_adv"] += \
		particle.vel + \
		particle.particleVariables["dt"] * particle.particleVariables["f_adv"] / particle.particleVariables["mass"]	

#   calculate dii
	for neighbor in particle.neighborList:
		pairData = pairsData[(particle,neighbor)]
		particle.particleVariables["d_ii"] -= \
			systemConstants["dt"]**2 * \
			neighbor.particleVariables["mass"] / (particle.particleVariables["rho"]**2) * \
			gW_4(pairData,systemConstants["interactionlen"])

#   calculate advection density
	particle.particleVariables["rho_adv "] += particle.particleVariables["rho"]
	for neighbor in particle.neighborList:
		pairData = pairsData[(particle,neighbor)]
		v_adv_ij = particle.particleVariables["v_adv"] - neighbor.particleVariables["v_adv"]
		particle.particleVariables["rho_adv"] += \
			systemConstants["dt"] * neighbor.particleVariables["mass"] * \
			v_adv_ij.dot(gW_4(pairData,systemConstants["interactionlen"]))

#   initialize pressure estimate.
	particle.particleVariables["pressure"] = 0.5 * particle.particleVariables["pressure"] 









	