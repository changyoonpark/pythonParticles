from multiprocessing.dummy import Pool as ThreadPool 


def foo (num1, num2):
	print(num1 + num2)


pool = ThreadPool(4)
data = [(1,2),(2,3),(5,6)]
pool.starmap(foo,data)
