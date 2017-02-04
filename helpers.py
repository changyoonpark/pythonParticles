from math import sqrt

class ParticlePair:
	def __init__ (self,pi,pj,relpos,relvel,dist,reldir):
		self.particlei = pi
		self.particlej = pj
		self.relpos = relpos
		self.relvel = relvel
		self.dist = dist
		self.reldir = reldir

	def __str__ (self):
		st = ""
		st += "Particle i : {}, Particle j : {}\n".format(self.particlei.pID,self.particlej.pID)
		st += "Rel Vel : {}\n".format(self.relvel)		
		st += "Rel Pos : {}\n".format(self.relpos)
		st += "Rel Dir : {}\n".format(self.reldir)		
		return st

	def swapIdx(self):
		return (self.particlej,self.particlei)

class Vec2:

	def __init__ (self,x,y,tup = None):
		if tup is not None :
			self.x = tup[0]
			self.y = tup[1]
		else:
			self.x = x
			self.y = y

	def length(self):
		return sqrt(self.x*self.x+self.y*self.y)

	def dot(self,other):
		return self.x * other.x + self.y * other.y

	def dir(self):
		return self / self.length()
		# return Vec2(self.x / self.length(), self.y / self.length())

	def __mul__ (self,other):
		return Vec2(self.x * other, self.y * other)

	def __rmul__ (self,other):
		return Vec2(self.x * other, self.y * other)

	def __truediv__(self,other):
		return Vec2(self.x / other, self.y / other)

	def __add__ (self,other):
		return Vec2(self.x + other.x, self.y + other.y)

	def __sub__ (self,other):
		return Vec2(self.x - other.x, self.y - other.y)
		
	def __neg__ (self):
		return Vec2(-self.x,-self.y)

	def __str__ (self):
		return "Vector : ({},{})".format(self.x, self.y)
