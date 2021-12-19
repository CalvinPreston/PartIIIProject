# Plot for apparent magnitude with redshift from low redshift (z=1)

import numpy as np
import magnitude_constants as cn
from scipy import integrate as intg
gamma = 0

class ExpHistory:
    """
    System of eqns:
    du0/dt = - ( gamma * u0 ) - ( 3 * u0 * H )
    du1/dt = - ( 3 * u1 * H )
    du2/dt = ( u0 * gamma ) - ( 4 * u2 * H )
    du3/dt = ( u3 * H )

    with

    u0 = rhoc, u1 = rhom, u2 = rhor, u3 = a

    H = ( ((8*np.pi*G)/3) * (u0 + u1 + u2 + cn.RHOL) )**0.5
    """
    def __init__(self, gamma):
        self.gamma = gamma

    def __call__(self, u, t):
        u0, u1, u2, u3 = u      # As DM density, Matter density, Radiation density, Scale factor
        gamma = self.gamma      # Avoiding need for self. at each calling of the variable
        if gamma == 0:
            return np.array([ - (3 * u0 * (cn.PreFac* (u0 + u1 + u2 + cn.RHOL))**(0.5)), (-3 * u1 *(cn.PreFac * (u0 + u1 + u2 + cn.RHOL))**(0.5)), - ( 4 * u2 * (cn.PreFac * (u0 + u1 + u2 + cn.RHOL))**(0.5)), u3*(cn.PreFac * (u0 + u1 + u2 + cn.RHOL))**(0.5)])
        else:
            return np.array([(-gamma * u0) - (3 * u0 * (( ((8*np.pi*cn.G)/3) * (u0 + u1 + u2 + cn.RHOL)))**(0.5)), (-3 * u1 *( ((8*np.pi*cn.G)/3) * (u0 + u1 + u2 + cn.RHOL))**(0.5)), (u0 * gamma) - ( 4 * u2 * ( ((8*np.pi*cn.G)/3) * (u0 + u1 + u2 + cn.RHOL))**(0.5)), u3*( ((8*np.pi*cn.G)/3) * (u0 + u1 + u2 + cn.RHOL))**(0.5)])


if __name__ == "__main__":
    from RK4Test import RungeKutta4
    from matplotlib import pyplot as plt
 
    n = int(9999)

    f = ExpHistory(gamma)
    ininital_conditions = [cn.RHOC * (cn.zstart+1)**3, cn.RHOB * (cn.zstart+1)**3 , 0 , 1 ] # Defining initial conditions

    time_points = np.linspace(cn.TStart, cn.TEnd, n + 1)

    for solver_class in [RungeKutta4]:
        solver = solver_class(f)
        solver.set_initial_conditions(ininital_conditions)
        u, t = solver.solve(time_points)
        z1 = 11/u[:, 3]-1
        H1 = ((((8*np.pi*cn.G)/3) * (u[:,0] + u[:,1] + u[:,2] + cn.RHOL)))**0.5


y1 = 1/((( ((8*np.pi*cn.G)/3) * (u[:,0] + u[:,1] + u[:,2] + cn.RHOL) ))**0.5)

rev1 = z1[::-1]
rev2 = y1[::-1]

mb=n*[0]

for i in range (1,n):
    mb[i] = cn.MB + 25 + 5 * np.log10( (cn.c) * (1+rev1[i]) * (intg.cumtrapz(rev2, rev1 )[i])/cn.Mpc )
mb.append(0)
plt.plot(rev1 , mb,color="red" )


gamma = cn.H0

if __name__ == "__main__":
    from RK4Test import RungeKutta4
    from matplotlib import pyplot as plt
 
    f = ExpHistory(gamma)
    ininital_conditions = [cn.RHOC * (cn.zstart+1)**3, cn.RHOB * (cn.zstart+1)**3 , 0 , 1 ] # Defining initial conditions

    time_points = np.linspace(cn.TStart, cn.TEnd, n + 1)

    for solver_class in [RungeKutta4]:
        solver = solver_class(f)
        solver.set_initial_conditions(ininital_conditions)
        u, t = solver.solve(time_points)
        z2 = 11/u[:, 3]-1
        H2 = (( ((8*np.pi*cn.G)/3) * (u[:,0] + u[:,1] + u[:,2] + cn.RHOL) ))**0.5

y2 = 1/((( ((8*np.pi*cn.G)/3) * (u[:,0] + u[:,1] + u[:,2] + cn.RHOL) ))**0.5)

rev1 = z2[::-1]
rev2 = y2[::-1]

mb=n*[0]

for i in range (1,n):
    mb[i] = cn.MB + 25 + 5 * np.log10( (cn.c) * (1+rev1[i]) * (intg.cumtrapz(rev2, rev1 )[i])/cn.Mpc )
mb.append(0)
plt.plot(rev1 , mb, color="orange")



gamma = 2*cn.H0

if __name__ == "__main__":
    from RK4Test import RungeKutta4
    from matplotlib import pyplot as plt
 
    f = ExpHistory(gamma)
    ininital_conditions = [cn.RHOC * (cn.zstart+1)**3, cn.RHOB * (cn.zstart+1)**3 , 0 , 1 ] # Defining initial conditions

    time_points = np.linspace(cn.TStart, cn.TEnd, n + 1)

    for solver_class in [RungeKutta4]:
        solver = solver_class(f)
        solver.set_initial_conditions(ininital_conditions)
        u, t = solver.solve(time_points)
        z3 = 11/u[:, 3]-1
        H3 = (( ((8*np.pi*cn.G)/3) * (u[:,0] + u[:,1] + u[:,2] + cn.RHOL) ))**0.5

y3 = 1/((( ((8*np.pi*cn.G)/3) * (u[:,0] + u[:,1] + u[:,2] + cn.RHOL) ))**0.5)

rev1 = z3[::-1]
rev2 = y3[::-1]

mb=n*[0]

for i in range (1,n):
    mb[i] = cn.MB + 25 + 5 * np.log10( (cn.c) * (1+rev1[i]) * (intg.cumtrapz(rev2, rev1 )[i])/cn.Mpc )
mb.append(0)
plt.plot(rev1 , mb , color="green" )

import pandas as pd
data = pd.read_csv('test.txt',sep='\s+',header=None)
data = pd.DataFrame(data)

P = data[0]
Q = data[3]

plt.scatter(P, Q)

plt.legend([r'$\Gamma = 0$',r'$\Gamma = H_0$',r'$\Gamma = 2H_0$',"1a supernovae"], title = "Value of dark matter decay constant")
plt.xlabel(r'$z$',fontsize = 16)
plt.ylabel(r'$ m_{B}(z)$', fontsize = 16)
plt.title("Apparent magnitude of 1a supernovae at low redshift with \n fits from varying late-time cosmologies")




plt.xlim([0.01, 1])
plt.ylim([15, 25])


plt.grid()
plt.show()
