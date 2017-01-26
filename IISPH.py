from kernels import W_4, gradW_4


def EOS(parameters):
	k = parameters["kappa"]
	rho0 = parameters["density"]
	rho = paremters