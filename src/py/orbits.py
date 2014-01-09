""" Decompose a permutation into cycles """
""" http://codegolf.stackexchange.com/questions/1668/decompose-a-permutation-into-cycles """

input();d=dict((i,int(x))for i,x in enumerate(raw_input().split()))
while d:
 x=list(d)[0]
 while x in d:print x,;x=d.pop(x)
 print
 
def permutationOrbits(List):
	d = dict((i,int(x)) for i,x in enumerate(List))
	out = []
	while d:
		x = list(d)[0]
		orbit = []
		while x in d:
			orbit += [x],
			x = d.pop(x)
		out += [CAT(orbit)+orbit[0]]
	return out
		
if __name__ == "__main__":
	permutationOrbits([2, 3, 4, 5, 6, 7, 0, 1])
	permutationOrbits([3,9,8,4,10,7,2,11,6,0,1,5])