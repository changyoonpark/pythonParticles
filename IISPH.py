from kernels import W_4, gW_4, lW_v


def calculateDensity(systemConstants,pairsData,particleSystem):
	for particle in particleSystem:
	#   calculate density
		particle.particleVariables["rho"] = 0
		for neighbor in particle.neighborList:
			pairData = pairsData[(particle,neighbor)]
			particle.particleVariables["rho"] += 
				particle.particleVariables["mass"] * W_4(pairData,systemConstants["interactionlen"])

def calculate_dii_dij(systemConstants,pairsData,particleSystem):
	for particle in particleSystem:
	#   calculate dii
		particle.particleVariables["d_ii"] = Vec2(0,0)
		for neighbor in particle.neighborList:
			pairData = pairsData[(particle,neighbor)]
			particle.particleVariables["d_ii"] -= \
				neighbor.particleVariables["mass"] / (particle.particleVariables["rho"]**2) * \
				gW_4(pairData,systemConstants["interactionlen"])

		particle.particleVariables["d_ii"] = particle.particleVariables["d_ii"] * systemConstants["dt"]**2

	#   calculate the sum of dij's
		particle.particleVariables["d_ij"] = Vec2(0,0)
		for neigbor in particle.neighborList:
			pairData = pairsData[(particle,neighbor)]
			particle.particleVariables["d_ij"] -= \
				neighbor.particleVariables["mass"] / (neighbor.particleVariables["rho"]**2) * \
				gW_4(pairData,systemConstants["interactionlen"])

		particle.particleVariables["d_ij"] = particle.particleVariables["d_ij"] * systemConstants["dt"]**2

def calcaulateAdvectionForce(systemConstants,pairsData,particleSystem):
	for particle in particleSystem:
	#	calculate advection force (viscosity and gravity is applied.)
		particle.particleVariables["f_adv"] = Vec2(0,0)
		for neighbor in particle.neighborList:
			pairData = pairsData[(particle,neighbor)]
			particle.particleVariables["f_adv"] += \
				systemConstants["viscosity"] * \
				lW_v(pairData,systemConstants["interactionlen"]) * \
				neighbor.particleVariables["mass"] / \
				neighbor.particleVariables["rho"] * \
				neighbor.relvel

	#   consider gravity.
		particle.particleVariables["f_adv"] -= \
			particle.particleVariables["mass"] * systemConstants["gravity"]

def calcaulteAdvectionVel(systemConstants,pairsData,particleSystem):
#   calculate itermediate velocity
	for particle in particleSystem:
		particle.particleVariables["v_adv"] = \
			particle.vel + \
			particle.particleVariables["dt"] * particle.particleVariables["f_adv"] / particle.particleVariables["mass"]	

def calculateAdvectionDensity(systemConstants,pairsData,particleSystem):
	for particle in particleSystem:
	#   calculate advection density
		particle.particleVariables["rho_adv"] = particle.particleVariables["rho"]
		for neighbor in particle.neighborList:
			pairData = pairsData[(particle,neighbor)]
			v_adv_ij = particle.particleVariables["v_adv"] - neighbor.particleVariables["v_adv"]
			particle.particleVariables["rho_adv"] += \
				systemConstants["dt"] * neighbor.particleVariables["mass"] * \
				v_adv_ij.dot(gW_4(pairData,systemConstants["interactionlen"]))

	#   initialize pressure estimate.
		particle.particleVariables["pressure"] = 0.5 * particle.particleVariables["pressure"] 

def calculate_aii(systemConstants,pairsData,particleSystem):
	for particle in particleSystem:	
	# 	calculate a_ii for relaxed jacobi
		particle.particleVariables["a_ii"] = 0
		for neighbor in particle.neighborList:
			pairData = pairsData[(particle,neighbor)]
			particle.particleVariables["a_ii"] += \
				neighbor.particleVariables["mass"] * (particle.particleVariables["d_ii"] - particle.particleVariables["d_ij"]).dot(gW_4(pairData,systemConstants["interactionlen"]))


def calcAverageDensity(particleSystem):
	rhoAve = 0
	for particle in particleSystem:
		rhoAve += particle.particleVariables["rho"]
	return rhoAve / len(particleSystem)

def pressureSolver(systemConstants,pairsData,particleSystem):
# 	iteration counter
	l = 0
	while ( abs((calcAverageDensity(particleSystem) - systemConstants["rho0"])/systemConstants["rho0"]) > 0.02 ):

		for particle in particleSystem:
			












