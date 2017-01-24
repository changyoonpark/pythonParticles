from multiprocessing import Process, Queue

class foo:
	def __init__ (self):
		self.counter = 0
		self.q = Queue()

	def inc(self,*by):
		for el in by:
			self.counter += el
		self.q.put(self.counter)
		# print(self.counter)

	def multithreadadd(self,a):
		thread1 = Process(target = self.inc, args=a[0])
		thread1.start()
		thread1.join()
		thread2 = Process(target = self.inc, args=a[1])
		thread2.start()
		thread2.join()
		foo = self.q.get()
		print(self.q.get())



instance = foo()
instance.multithreadadd([[2,3],[5,10]])
# print(instance.counter)

