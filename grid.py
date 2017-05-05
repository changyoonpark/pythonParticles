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

	def W(self,particlePos, smoothing):
		q = (particlePos - self.nodePos) / (self.h * (smoothing / 2.0))
		return self.N( q.x ) * self.N( q.y )

	def gW(self,particlePos, smoothing):
		q = (particlePos - self.nodePos) / (self.h * (smoothing / 2.0))
		return Vec2(self.Nx(q.x) * self.N (q.x),
					self.N (q.y) * self.Nx(q.y)) / self.h

	# def W(self,particlePos):
	# 	q = (particlePos - self.nodePos) / (2. * self.h)
	# 	return self.N( q.x ) * self.N( q.y )

	# def gW(self,particlePos):
	# 	q = (particlePos - self.nodePos) / (2. * self.h)
	# 	return Vec2(self.Nx(q.x) * self.N (q.x),
	# 				self.N (q.y) * self.Nx(q.y)) / self.h

	# def nablW(self,particlePos):
	# 	q = (particlePos - self.nodePos) / self.h
	# 	return  

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

	# def Nxx(self,x):
	# 	if x == 0:
	# 		return 0

	# 	abx = abs(x)
	# 	if (abx < 1 and abx >= 0):
	# 		return x * (3./2. * abx - 2.)
	# 	elif (abx >= 1 and abx < 2):
	# 		return (abx / x) * (- 0.5 * abx * abx + 2. * abx - 2.);
	# 	else :
	# 		return 0


class Grid:

	def __init__(self, particleSet, h, smoothing):
		# self.nx = int(domainExtent[0] / h) + 1
		# self.ny = int(domainExtent[1] / h) + 1
		# self.domainExtent = domainExtent
		self.h = h
		self.nodes = dict()

		for p in particleSet:
			(m, n) = self.hashFunction(p.pos)
			for i in range(1-smoothing,1+smoothing):
				for j in range(1-smoothing,1+smoothing):
					nodePos = Vec2((m+i) * self.h, (n+j) * self.h)
					if self.nodes.get(self.hashFunction(nodePos)) is None :
						self.nodes[self.hashFunction(nodePos)] = Node(nodePos, h);
						# print("adding node at {}".format(self.hashFunction(nodePos)))

	def sampleScalarGradientFromGrid(self,pos,quantity,smoothing):
		(m, n) = self.hashFunction(pos)
		q = Vec2(0,0)
		for i in range(1-smoothing,1+smoothing):
			for j in range(1-smoothing,1+smoothing):
				nodePos = Vec2((m+i) * self.h, (n+j) * self.h)
				node = self.nodes.get(self.hashFunction(nodePos))
				if node is None :
					continue
				q += node.quantities[quantity] * node.gW(pos, smoothing)

		return q


	def sampleScalarFromGrid(self,pos,quantity, smoothing):
		(m, n) = self.hashFunction(pos)
		q = 0
		for i in range(1-smoothing,1+smoothing):
			for j in range(1-smoothing,1+smoothing):
				nodePos = Vec2((m+i) * self.h, (n+j) * self.h)
				node = self.nodes.get(self.hashFunction(nodePos))
				if node is None :
					continue
				q += node.quantities[quantity] * node.W(pos, smoothing)

		return q
	
	# def sampleVecFromGrid(self,pos,quantity):
	# 	(m, n) = self.hashFunction(pos)
	# 	q = Vec2(0,0)
	# 	for i in range(-1, 3):
	# 		for j in range(-1, 3):
	# 			nodePos = Vec2((m+i) * self.h, (n+j) * self.h)
	# 			node = self.nodes.get(self.hashFunction(nodePos))
	# 			if node is None :
	# 				continue
	# 			q += node.quantities[quantity] * node.W(pos)

	# 	return q

	def hashFunction(self,pos):
		return (int( (pos.x + 1E-10) / self.h ),int( (pos.y + 1E-10) / self.h ))

	
