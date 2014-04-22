import numpy

def exp1(x,n, debug=False):
	"""Approximate the exponential function with a Taylor Series"""
	if n<0:
		print "***Invalid input -- n must be non-negative integer"
		return numpy.nan
	term=1.0;
	sum=term;
	for i in range(1,n+1):
		term=term*x/i
		if debug:
			print "j = %i, term = %.15f" %(i,term)
			print "partial_sum updated from %.15f to %.15f" %(sum,sum+term)
		sum=sum+term
	return sum

def sin1(x,n,debug=False):
	"""Approximate the sine function with a Maclaurin Series"""
	if n<0:
		print "***Invalid input -- n must be non-negative integer"
		return numpy.nan	
	term=1.0
	sum=0.0;
	for i in range(1,n+1):
		term=term*x/i
		if i%2==0:
			term_fixed=0.0 #Set even terms to zero
		else:
			term_fixed=term
			if i%4==3:
				term_fixed=term_fixed*-1
		if debug:
			print "j = %i, term = %.15f" %(i,term_fixed)
			print "partial_sum updated from %.15f to %.15f" %(sum,sum+term_fixed)
		sum=sum+term_fixed
	return sum
