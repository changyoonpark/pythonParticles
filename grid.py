	# 	inline double W(Vector3D xp) const {
	# 		return N((xp.x - x.x) / consts.h) * N((xp.y - x.y) / consts.h) * N((xp.z - x.z) / consts.h);
	# 	}

	# 	inline Vector3D gW(Vector3D xp) const {
	# 		double xeval = (xp.x - x.x) / consts.h;
	# 		double yeval = (xp.y - x.y) / consts.h;
	# 		double zeval = (xp.z - x.z) / consts.h;
	# 		return Vector3D(  (1. / consts.h) * Nx(xeval) * N (yeval) * N (zeval),
	# 						  (1. / consts.h) * N (xeval) * Nx(yeval) * N (zeval),
	# 						  (1. / consts.h) * N (xeval) * N (yeval) * Nx(zeval));
	# 	}

	# private:


	# 	inline double N(double _x) const {
	# 		double abx = std::abs(_x);
	# 		if (abx < 1. && abx >= 0.)      return 0.5 * abx * abx * abx - abx * abx + 2./3.;
	# 		else if (abx >= 1. && abx < 2.) return - (1. / 6.) * abx * abx * abx + abx * abx - 2. * abx + 4./3.;
	# 		else return 0.;			
	# 	}

	# 	inline double Nx(double _x) const {
	# 		double abx = std::abs(_x);
	# 		double sign = (_x >= 0 ? 1.0 : -1.0);
	# 		if (abx < 1. && abx >= 0.)      return sign * (3./2. * abx * abx - 2. * abx);
	# 		else if (abx >= 1. && abx < 2.) return sign * (- 0.5 * abx * abx + 2. * abx - 2.);
	# 		else return 0.;			
	# 	}
from helpers import *

class Node:

	def __init__(self,nodePos,h):
		self.nodePos = nodePos		
		self.h = h
		self.quantities = dict()


	def W(self,particlePos):
		q = (particlePos - self.nodePos) / self.h
		return N( q.x ) * N( q.y )

	def gW(self,particlePos):
		q = (particlePos - self.nodePos) / self.h
		return Vec2(Nx(q.x) * N (q.x),
					N (q.y) * Nx(q.y))

	def N(self,x):
		abx = abs(x)
		if (abx < 1 and abx >= 0):
			return 0.5 * abx * abx * abx - abx * abx + 2.0/3.0

		elif abx >= 1 and abx < 2 :
			return - (1. / 6.) * abx * abx * abx + abx * abx - 2. * abx + 4./3.

		else :
			return 0

	def Nx(self,x):
		if x == 0:
			return 0

		abx = abs(x)
		if (abx < 1 and abx >= 0):
			return x * (3./2. * abx - 2.)
		elif (abx >= 1 and abx < 2):
			return (abx / x) * (- 0.5 * abx * abx + 2. * abx - 2.);
		else :
			return 0



class Grid:

	def __init__(self, domainExtent, h):
		self.nx = int(domainExtent[0] / h) + 1
		self.ny = int(domainExtent[1] / h) + 1
		self.domainExtent = domainExtent
		self.h = h
		self.nodes = dict()

		for x in range(0,self.nx):
			for y in range(0,self.ny):
				nodePos = Vec2(x * self.h, y * self.h);
				self.nodes[self.hashFunction(nodePos)] = Node(nodePos, h);

	def hashFunction(pos):
		return (int(pos.x / self.h),int(pos.y / self.h))

	
