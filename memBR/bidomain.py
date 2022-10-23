#bidomain
 
from esys.escript import *
from esys.escript.linearPDEs import LinearPDE, SolverOptions
from esys.escript.pdetools import Locator
from esys.finley import Rectangle
from esys.finley import ReadMesh
 
#bidomain conductivity values
sigma_int_l =  2.8
sigma_int_t =  0.26
sigma_ext_l =  2.2
sigma_ext_t =  1.3
 
theta=math.pi*0.0/180
 
# current stimulus
xci=[0.5,0.5]
xce=[0.5,0.5]
xr=[0.0,0.0]
r=0.02 #0.02
intra_stim=100 #100
extra_stim=-100
stim_duration=1.0
off=0

# tissue domain
mesh = Rectangle(l0=1.0,l1=1.0,n0=100, n1=100)
x = mesh.getX()
 
#intial membrane potential
Vm = Scalar(-84.0, ContinuousFunction(mesh))

Isi = Scalar(0.0, ContinuousFunction(mesh))
Ise = Scalar(0.0, ContinuousFunction(mesh))
iion = Scalar(0.0, ContinuousFunction(mesh))
w = Scalar(0.0, ContinuousFunction(mesh))
V = Vector(0.0, ContinuousFunction(mesh))
#Vm = V[0]
Ve = V[1]

#membrane parameters
Gna = 4.0
Gnc = 0.003
Ena = 50.0
Gs = 0.09

#initial conditions for gating variables
ax1 = (0.0005*exp(0.083*(Vm + 50.0)))/(exp(0.057*(Vm + 50.0)) + 1.0)
bx1 = (0.0013*exp(-0.06*(Vm + 20.0)))/(exp(-0.04*(Vm + 20.0)) + 1.0)
am = (-1.0*(Vm + 47.0))/(exp(-0.1*(Vm + 47.0)) - 1.0)
bm = (40.0*exp(-0.056*(Vm + 72.0)))
ah = (0.126*exp(-0.25*(Vm + 77.0)))
bh = (1.7)/(exp(-0.082*(Vm + 22.5)) + 1.0)
aj = (0.055*exp(-0.25*(Vm + 78.0)))/(exp(-0.2*(Vm + 78.0)) + 1.0)
bj = (0.3)/(exp(-0.1*(Vm + 32.0)) + 1.0)
ad = (0.095*exp(-0.01*(Vm -5.0)))/(exp(-0.072*(Vm - 5.0)) + 1.0)
bd = (0.07*exp(-0.017*(Vm + 44.0)))/(exp(0.05*(Vm + 44.0)) + 1.0)
af = (0.012*exp(-0.008*(Vm + 28.0)))/(exp(0.15*(Vm + 28.0)) + 1.0)
bf = (0.0065*exp(-0.02*(Vm + 30.0)))/(exp(-0.2*(Vm + 30.0)) + 1.0)

#steady state
m = am/(am + bm)
h = ah/(ah + bh)
j = aj/(aj + bj)
d = ad/(ad + bd)
f = af/(af + bf)
x1 = ax1/(ax1 + bx1)

#initial resting intracellular Ca concentration
cai = 2e-4
 
# membrane capacitance and surface to volume ratio
Cm=1.0
beta=2000
 
# time-step
dt=0.01
tend=300.0
t=0
i=0
output_freq1=1.0
output_freq2=0.1
 
#solver options
TOL=1.0e-8
MAXIT=800
 
# boundary conditions
b_c = Vector(0.0, ContinuousFunction(mesh))
b_c[1]=whereZero(length(x-xr)-r)
b_c[0]=whereZero(length(x-xr)-r)
#b_c[1] = whereZero(x[1])
#b_c[0] = whereZero(x[1])
 
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
mypde.setValue(A=A,r=[-84.0,0.0],q=b_c)
#mypde.getSolverOptions().setVerbosityOn()
mypde.getSolverOptions().setPreconditioner(SolverOptions.JACOBI)
mypde.getSolverOptions().setSolverMethod(SolverOptions.PCG)
mypde.getSolverOptions().setTolerance(TOL)
mypde.getSolverOptions().setIterMax(MAXIT)
 
# current source
Isi=Isi+intra_stim*whereNegative(length(x-xci)-r)
Ise=Ise+extra_stim*whereNegative(length(x-xce)-r)
 
# iteration
while t<=tend:
      t+=dt
      print "t=",t
 
      #Beeler-Reuter membrane model
      ax1 = (0.0005*exp(0.083*(Vm + 50.0)))/(exp(0.057*(Vm + 50.0)) + 1.0)
      bx1 = (0.0013*exp(-0.06*(Vm + 20.0)))/(exp(-0.04*(Vm + 20.0)) + 1.0)
      am = (-1.0*(Vm + 47.0))/(exp(-0.1*(Vm + 47.0)) - 1.0)
      bm = (40.0*exp(-0.056*(Vm + 72.0)))
      ah = (0.126*exp(-0.25*(Vm + 77.0)))
      bh = (1.7)/(exp(-0.082*(Vm + 22.5)) + 1.0)
      aj = (0.055*exp(-0.25*(Vm + 78.0)))/(exp(-0.2*(Vm + 78.0)) + 1.0)
      bj = (0.3)/(exp(-0.1*(Vm + 32.0)) + 1.0)
      ad = (0.095*exp(-0.01*(Vm -5.0)))/(exp(-0.072*(Vm - 5.0)) + 1.0)
      bd = (0.07*exp(-0.017*(Vm + 44.0)))/(exp(0.05*(Vm + 44.0)) + 1.0)
      af = (0.012*exp(-0.008*(Vm + 28.0)))/(exp(0.15*(Vm + 28.0)) + 1.0)
      bf = (0.0065*exp(-0.02*(Vm + 30.0)))/(exp(-0.2*(Vm + 30.0)) + 1.0)

      Es = -82.3 - 13.0287*log(cai)
      Is = Gs*d*f*(Vm - Es)
      Ik1 = 0.35*(4.0*(exp(0.04*(Vm + 85.0)) - 1.0)/(exp(0.08*(Vm + 53.0)) + exp(0.04*(Vm + 53.0))) + 0.2*(Vm + 23.0)/(1.0 -exp(-0.04*(Vm + 23.0))))
      Ix1 = x1*0.8*(exp(0.04*(Vm + 77.0)) - 1.0)/exp(0.04*(Vm + 35.0))
      Ina = (Gna*m*m*m*h*j + Gnc)*(Vm - Ena)

      iion = -(1/Cm)*(Ik1 + Ix1 + Ina + Is)

      m = m + dt*(am*(1.0-m) - bm*m)
      h = h + dt*(ah*(1.0-h) - bh*h)
      j = j + dt*(aj*(1.0-j) - bj*j)
      d = d + dt*(ad*(1.0-d) - bd*d)
      f = f + dt*(af*(1.0-f) - bf*f)
      x1 = x1 + dt*(ax1*(1.0-x1) - bx1*x1)

      #calcium uptake
      cai = cai + dt*((-10e-7)*Is + 0.07*((10e-7) - cai))
 
      #bidomain matrix equation form
      #-grad dot A grad[Vm Ve] + D[Vm Ve] = y
 
      D[0,0] = beta*Cm/dt
 
      Y[0] = (beta*Cm*Vm/dt) + beta*iion + beta*Isi #intra stim
      Y[1] = beta*Ise #extra stim
 
      mypde.setValue(D=D,Y=Y)
      V=mypde.getSolution()
      Vm = V[0]
      Ve = V[1]
 
      #total=0
      #num_points=Ve.getNumberOfDataPoints()
      #for k in range(1, num_points):
      #    total=total+Ve.getTupleForDataPoint(k)[0]
      #Ve_average=total/num_points
 
      #Ve_ref=Ve-Ve_average
 
      # find potential at point source
      Ve_pc=L.getValue(Ve)
      ts.append(t); Ve_p.append(Ve_pc)
 
      #turn off stimulus
      if ((t>=stim_duration)and(off==0)):
         Isi=Isi-intra_stim*whereNegative(length(x-xci)-r)
         Ise=Ise-extra_stim*whereNegative(length(x-xce)-r)
         off=1
 
      if (i%(output_freq1/dt)==0):
         saveVTK("V.%2.5i.vtu"%i,extracellular_voltage=V[1],transmembrane_voltage=V[0])
      i+=1
 
for i in xrange(len(ts)):
    if (i%(output_freq2/dt)==0):
       Ve_p_data.write("%f "%(ts[i]))
       for j in xrange(len(Ve_p[0])):
           Ve_p_data.write("%f "%(Ve_p[i][j]))
       Ve_p_data.write("\n")
Ve_p_data.close()
