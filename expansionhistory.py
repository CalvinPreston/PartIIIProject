import numpy as np
import constants as cn

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

    H = ( ((8*np.pi*G)/3) * (u0 + u1 + u2 + RHOL) )**0.5
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
 
    n = int(1e4)

    f = ExpHistory(gamma)
    ininital_conditions = [cn.RHOC * (cn.zstart+1)**3, cn.RHOB * (cn.zstart+1)**3 , 0 , 1 ] # Defining initial conditions

    time_points = np.linspace(cn.TStart, cn.TEnd, n + 1)

    plt.figure()

    for solver_class in [RungeKutta4]:
        solver = solver_class(f)
        solver.set_initial_conditions(ininital_conditions)
        u, t = solver.solve(time_points)
        z1 = 11/u[:, 3]-1
        H1 = cn.Mpc * (( ((8*np.pi*cn.G)/3) * (u[:,0] + u[:,1] + u[:,2] + cn.RHOL) ))**0.5
        plt.xscale("log")
        plt.xlabel(r'$z$')
        plt.ylabel(r'$ \frac{H(z)}{(1+z)^{\frac{3}{2}}} $') 
        plt.title("Expansion history as a function of redshift")


gamma = cn.H0

if __name__ == "__main__":
    from RK4Test import RungeKutta4
    from matplotlib import pyplot as plt
 
    f = ExpHistory(gamma)
    ininital_conditions = [cn.RHOC * (cn.zstart+1)**3, cn.RHOB * (cn.zstart+1)**3 , 0 , 1 ] # Defining initial conditions

    time_points = np.linspace(cn.TStart, cn.TEnd, n + 1)

    plt.figure()

    for solver_class in [RungeKutta4]:
        solver = solver_class(f)
        solver.set_initial_conditions(ininital_conditions)
        u, t = solver.solve(time_points)
        z2 = 11/u[:, 3]-1
        H2 = cn.Mpc * (( ((8*np.pi*cn.G)/3) * (u[:,0] + u[:,1] + u[:,2] + cn.RHOL) ))**0.5

gamma = 2*cn.H0

if __name__ == "__main__":
    from RK4Test import RungeKutta4
    from matplotlib import pyplot as plt
 
    f = ExpHistory(gamma)
    ininital_conditions = [cn.RHOC * (cn.zstart+1)**3, cn.RHOB * (cn.zstart+1)**3 , 0 , 1 ] # Defining initial conditions

    time_points = np.linspace(cn.TStart, cn.TEnd, n + 1)

    plt.figure()

    for solver_class in [RungeKutta4]:
        solver = solver_class(f)
        solver.set_initial_conditions(ininital_conditions)
        u, t = solver.solve(time_points)
        z3 = 11/u[:, 3]-1
        H3 = cn.Mpc * (( ((8*np.pi*cn.G)/3) * (u[:,0] + u[:,1] + u[:,2] + cn.RHOL) ))**0.5


plt.plot(z1, H1 / ( (1+z1)**1.5) )
plt.plot(z2, H2 / ( (1+z2)**1.5) )
plt.plot(z3, H3 / ( (1+z3)**1.5) )
plt.xscale("log")
plt.xlabel(r'$z$')
plt.ylabel(r'$ \frac{H(z)}{(1+z)^{\frac{3}{2}}} [km s^{-1} Mpc^{-1}]$', fontsize = 16)
plt.title("Expansion history as a function of redshift")
plt.legend(["Gamma = 0", "Gamma = $H_0$", "Gamma = $2H_0$ "])
plt.show()
