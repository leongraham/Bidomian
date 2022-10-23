#bidomain

from esys.escript import *
from esys.escript.linearPDEs import LinearPDE, SolverOptions
from esys.escript.pdetools import Locator
from esys.finley import Rectangle
from esys.finley import ReadMesh

#Crank-Nicolson: alpha=0.5, Semi-Implicit: alpha=1.0
alpha=1.0

#bidomain conductivity values
sigma_int_l = alpha*2.8
sigma_int_t = alpha*0.26
sigma_ext_l = alpha*2.2
sigma_ext_t = alpha*1.3

theta=math.pi*0.0/180

# current stimulus
xc=[0.5,0.5]
xr=[0.0,0.0]
r=0.02
qc=-100
stim_duration=1.0
off=0

# membrane model parameters
FH_Vt=5.0
FH_Vp=100.0
C1 = 0.5
C2 = 0.02
C3 = 0.015
C4 = 0.005

# membrane capacitance and surface to volume ratio 
Cm=1.0
beta=2000

# time-step
dt=0.01
tend=120.0
t=0
i=0
output_freq1=1.0
output_freq2=0.1

#solver options
TOL=1.0e-5
MAXIT=800

# tissue domain
mesh = Rectangle(l0=1.0,l1=1.0,n0=100, n1=100)
x = mesh.getX()

# boundary conditions
b_c = Vector(0.0, ContinuousFunction(mesh))
#b_c[1]=whereZero(length(x-xr)-r)
b_c[1] = whereZero(x[1])
b_c[0] = whereZero(x[1])

Y = Vector(0.0,Function(mesh))
D = Tensor(0.0,Function(mesh))
A = Tensor4(0.0, Function(mesh))
X = Tensor(0.0, Function(mesh))

#translate from local tensor to global tensor
#Kglobal = A*Klocal*A^T
#where A is the rotation matrix
k00_i=sigma_int_l*cos(theta)*cos(theta) + sigma_int_t*sin(theta)*sin(theta)
k01_i=(sigma_int_l-sigma_int_t)*cos(theta)*sin(theta)
k10_i=(sigma_int_l-sigma_int_t)*cos(theta)*sin(theta)
k11_i=sigma_int_l*sin(theta)*sin(theta) + sigma_int_t*cos(theta)*cos(theta)

k00_e=sigma_ext_l*cos(theta)*cos(theta) + sigma_ext_t*sin(theta)*sin(theta)
k01_e=(sigma_ext_l-sigma_ext_t)*cos(theta)*sin(theta)
k10_e=(sigma_ext_l-sigma_ext_t)*cos(theta)*sin(theta)
k11_e=sigma_ext_l*sin(theta)*sin(theta) + sigma_ext_t*cos(theta)*cos(theta)

#fill A tensor

A[0,0,0,0] = k00_i 
A[0,0,0,1] = k01_i
A[0,1,0,0] = k10_i
A[0,1,0,1] = k11_i 

A[0,0,1,0] = k00_i
A[0,0,1,1] = k01_i
A[0,1,1,0] = k10_i
A[0,1,1,1] = k11_i

A[1,0,0,0] = k00_i
A[1,0,0,1] = k01_i
A[1,1,0,0] = k10_i
A[1,1,0,1] = k11_i 

A[1,0,1,0] = k00_i + k00_e
A[1,0,1,1] = k01_i + k01_e
A[1,1,1,0] = k10_i + k10_e
A[1,1,1,1] = k11_i + k11_e 

stim = Scalar(0.0, ContinuousFunction(mesh))
iion = Scalar(0.0, ContinuousFunction(mesh))
w = Scalar(0.0, ContinuousFunction(mesh))
Vm = Scalar(0.0, ContinuousFunction(mesh))
V = Vector(0.0, ContinuousFunction(mesh))
Vm = V[0]
Ve = V[1]

#electrogram output
xp=[[0.1,0.1],[0.2,0.1],[0.3,0.1],[0.4,0.1],[0.5,0.1],[0.6,0.1],[0.7,0.1],[0.8,0.1],[0.9,0.1],[0.1,0.2],[0.2,0.2],[0.3,0.2],[0.4,0.2],[0.5,0.2],[0.6,0.2],[0.7,0.2],[0.8,0.2],[0.9,0.2],[0.1,0.3],[0.2,0.3],[0.3,0.3],[0.4,0.3],[0.5,0.3],[0.6,0.3],[0.7,0.3],[0.8,0.3],[0.9,0.3],[0.1,0.4],[0.2,0.4],[0.3,0.4],[0.4,0.4],[0.5,0.4],[0.6,0.4],[0.7,0.4],[0.8,0.4],[0.9,0.4],[0.1,0.5],[0.2,0.5],[0.3,0.5],[0.4,0.5],[0.6,0.5],[0.7,0.5],[0.8,0.5],[0.9,0.5],[0.1,0.6],[0.2,0.6],[0.3,0.6],[0.4,0.6],[0.5,0.6],[0.6,0.6],[0.7,0.6],[0.8,0.6],[0.9,0.6],[0.1,0.7],[0.2,0.7],[0.3,0.7],[0.4,0.7],[0.5,0.7],[0.6,0.7],[0.7,0.7],[0.8,0.7],[0.9,0.7],[0.1,0.8],[0.2,0.8],[0.3,0.8],[0.4,0.8],[0.5,0.8],[0.6,0.8],[0.7,0.8],[0.8,0.8],[0.9,0.8],[0.1,0.9],[0.2,0.9],[0.3,0.9],[0.4,0.9],[0.5,0.9],[0.6,0.9],[0.7,0.9],[0.8,0.9],[0.9,0.9]]
L=Locator(mesh,xp)
ts, Ve_p,=[], []
Ve_pc=L.getValue(Ve)
Ve_p_data=FileWriter('Ve.txt')

# output at t=0
saveVTK("V.%2.5i.vtu"%i,extracellular_voltage=Ve,transmembrane_voltage=Vm)
ts.append(t); Ve_p.append(Ve_pc)
i+=1

#... open PDE ...
mypde=LinearPDE(mesh)
mypde.setSymmetryOn()
mypde.setValue(A=A)
#mypde.getSolverOptions().setVerbosityOn()
mypde.getSolverOptions().setPreconditioner(SolverOptions.JACOBI)
mypde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
mypde.getSolverOptions().setTolerance(TOL)
mypde.getSolverOptions().setIterMax(MAXIT)

# current source
stim=stim+qc*whereNegative(length(x-xc)-r)

# iteration
while t<=tend:
      t+=dt
      print "t=",t
     
      #FitzHugh-Nagumo membrane model
      iion = (C1*Vm*(Vm/FH_Vt-1.0)*(1.0-Vm/FH_Vp) - C2*Vm*w)
      w += dt*C3*(Vm - C4*w)
     
      #bidomain matrix equation form
      #-grad dot A grad[Vm Ve] + D[Vm Ve] = y

      D[0,0] = beta*Cm/dt

      Y[0] = (beta*Cm*Vm/dt) + beta*iion - beta*stim
      Y[1] = beta*stim

      #gVm=grad(Vm)
      #gVe=grad(Ve)

      #X[0,0] = -k00_i*gVm[0] - k00_i*gVe[0]
      #X[0,1] = -k11_i*gVm[1] - k11_i*gVe[1]
      #X[1,0] = -k00_i*gVm[0] - (k00_i + k00_e)*gVe[0]
      #X[1,1] = -k11_i*gVm[1] - (k11_i + k11_e)*gVe[1]

      mypde.setValue(D=D,Y=Y,q=b_c)
      V=mypde.getSolution()
      Vm = V[0]
      Ve = V[1]

      # find potential at point source
      Ve_pc=L.getValue(Ve)
      #print "Ve at point=",Ve_pc
      ts.append(t); Ve_p.append(Ve_pc)

      #turn off stimulus
      if ((t>=stim_duration)and(off==0)):
         stim=stim-qc*whereNegative(length(x-xc)-r)
         off=1

      if (i%(output_freq1/dt)==0):
         saveVTK("V.%2.5i.vtu"%i,extracellular_voltage=V[1],transmembrane_voltage=V[0])
      i+=1

for i in xrange(len(ts)):
    if (i%(output_freq2/dt)==0):
       Ve_p_data.write("%f "%(ts[i]))
       #print ts[i]
       for j in xrange(len(Ve_p[0])):
           Ve_p_data.write("%f "%(Ve_p[i][j]))
           #print Ve_p[i][j]
       Ve_p_data.write("\n")
Ve_p_data.close()
