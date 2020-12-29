def ppmErr(m,m0):
   return (m-m0)/m0*1e6

def getMassTolRange(m,ppm):
   dM = ppm*m/1e6
   return [m-dM,m+dM]

def similarityScore_gauss(deltaVal,tol):
   return exp(-0.5 * (abs(deltaVal)/tol)^2)

def similarityScore_laplace(deltaVal,tol):
   return exp(-1 * (abs(deltaVal)/tol))

