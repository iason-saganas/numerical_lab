#!/usr/bin/env python

# Author(s): Parth Nayak, Michael Walther

import matplotlib
matplotlib.use('qtagg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import ode
import warnings
import argparse

warnings.filterwarnings('ignore')
parser = argparse.ArgumentParser()

tstart  = 0
astart   = np.array([1., 1.])
parser.add_argument("--show_dots", help = "Whether to show the previously selected points on the parameter space, True if yes, False if no, default: False")
parser.add_argument("--markersize", help = "If show_dots is True, what is the marker size for the dots to use, markersize ~ (dot's sidelength in pt)^2, default: 8")
args     = parser.parse_args() 
if args.show_dots is None: args.show_dots = False 
if args.markersize is None: args.markersize = 8

def rkdp(dydx,y0,a,b,params=None,tol=1e-6):
    """
    Integrate the system of ODEs given by y'(x,y)=dydx(x,y)
    with initial condition y(a)=y0 from a to b using the
    4th/5th-order Dormand-Prince Runge-Kutta method with
    automatic stepsize adjustment controlled by tolerance tol,
    and passing the parameters list to the dydx function.
    """
    xo = [a]
    nbad = [0]
    nfun = [0]
    def func(x,y,params):
        nfun[0] += 1
        if x<xo[0]:
            nbad[0] += 1
        xo[0] = x
        return dydx(x,y,*params)
    xa = []
    ya = []
    ngood = [0]
    def goodstep(x,y):
        ngood[0] += 1
        xa.append(x)
        ya.append(np.copy(y))
    solver = ode(func)
    solver.set_integrator('dopri5',rtol=tol,atol=1e-15,nsteps=10000,max_step=0.02,verbosity=0)
    solver.set_solout(goodstep)
    solver.set_initial_value(y0,t=a)
    solver.set_f_params(params)
    solver.integrate(b)
    return np.array(xa),np.array(ya),nfun[0],ngood[0],nbad[0],ngood[0]+nbad[0]


def dydx_friedmann(x,y,omega_m,omega_l):
    return (1 - omega_m - omega_l + omega_m/y + y**2*omega_l)**0.5

def dydx_friedmann_full(x,y,omega_m,omega_l):
    return np.array([y[1], omega_l*y[0] - omega_m/(2.*y[0]**2)])

class Cursor:
    """
    A cross hair cursor.
    """
    def __init__(self, ax):
        self.ax = ax
        self.horizontal_line = ax.axhline(color='r', lw=0.8, ls='-')
        self.vertical_line = ax.axvline(color='r', lw=0.8, ls='-')
        # text location in axes coordinates
        self.text = ax.text(0.62, 0.9, '', transform=ax.transAxes)

    def set_cross_hair_visible(self, visible):
        need_redraw = self.horizontal_line.get_visible() != visible
        self.horizontal_line.set_visible(visible)
        self.vertical_line.set_visible(visible)
        self.text.set_visible(visible)
        return need_redraw

    def on_mouse_move(self, event):
        if not event.inaxes:
            need_redraw = self.set_cross_hair_visible(False)
            if need_redraw:
                self.ax.figure.canvas.draw()
        else:
            self.set_cross_hair_visible(True)
            x, y = event.xdata, event.ydata
            # update the line positions
            self.horizontal_line.set_ydata(y)
            self.vertical_line.set_xdata(x)
            self.text.set_text('$\Omega_M$ = %.2f, $\Omega_\Lambda$ = %.2f' % (x, y))
            self.ax.figure.canvas.draw()

    
class Onclick:
    def __init__(self, omega_m, omega_l, ax_evo, ax_param):
        self.omega_m = omega_m
        self.omega_l = omega_l
        self.ax = ax_evo
        self.ax_param = ax_param
        self.count = 0

    def onclick(self, event):
        self.omega_m = np.round(event.xdata, 3)
        self.omega_l = np.round(event.ydata, 3)
        print(f"Omega_M = {self.omega_m:.2f}, Omega_L = {self.omega_l:.2f}")
        self.ax.clear()
        if (self.omega_m < 1.5)&(np.abs(self.omega_l)<0.05):
            tend = 100
        else:
            tend = 15
        t_pos,a_pos,nf,ng,nb,nt = rkdp(dydx_friedmann_full, astart, tstart,  tend, params=(self.omega_m, self.omega_l))
        t_neg,a_neg,nf,ng,nb,nt = rkdp(dydx_friedmann_full, astart, tstart, -tend, params=(self.omega_m, self.omega_l))
        self.ax.plot(t_pos, a_pos[:, 0], color = 'black', label = '$\Omega_M = $%.2f, $\Omega_\Lambda = $%.2f'%(self.omega_m, self.omega_l))
        self.ax.plot(t_neg, a_neg[:, 0], color = 'black')
        hline = self.ax.axhline(color='k', lw=0.8, ls='--')
        hline.set_ydata(1.)
        vline = self.ax.axvline(color='k', lw=0.8, ls='--')
        vline.set_xdata(0.)
        if args.show_dots:
            self.ax_param.scatter(self.omega_m,self.omega_l,color=f'C{self.count}', s = float(args.markersize))
        if np.max(a_pos[:,0]) > 10.:
            self.ax.set_ylim(0., 10.)
            self.ax.set_xlim(np.min(t_neg[a_neg[:,0]<10])-0.15, t_pos[a_pos[:,0]>10][0]+0.15)
        else:
            self.ax.set_ylim(0., np.max(a_pos[:,0]+0.05))
        self.ax.set_xlabel('$H_0 t$')
        self.ax.set_ylabel('$a(t)$')
        self.ax.legend()
        self.ax.grid(markevery = 0.1, alpha = 0.2)
        self.ax.figure.canvas.draw()
        #self.ax.clear()
        self.count+=1
        #self.ax.draw()


boundary1  = np.loadtxt('boundary1.txt')
boundary2  = np.loadtxt('boundary2.txt')

omega_m_b1 = boundary1[:,0]
omega_l_b1 = boundary1[:,1]
omega_m_b2 = boundary2[:,0]
omega_l_b2 = boundary2[:,1]


fig, (ax_paramspace, ax_evolution) = plt.subplots(1, 2, figsize = (15, 10))
fig.canvas.setWindowTitle("Omega_M vs Omega_L")

ax_paramspace.plot(omega_m_b1,omega_l_b1, ':k', alpha = 0.5)
ax_paramspace.plot(omega_m_b2,omega_l_b2, ':k', alpha = 0.5)
ax_paramspace.fill_between(omega_m_b1, omega_l_b1, -1, facecolor = 'grey', alpha = 0.2)
ax_paramspace.fill_between(omega_m_b2, omega_l_b2, 3, facecolor = 'grey', alpha = 0.2)

ax_paramspace.text(0.3, 2.5, 'No big bang', rotation = 40, rotation_mode = 'anchor')
ax_paramspace.text(0.7, 2.15, 'loitering', rotation = 49, rotation_mode = 'anchor')
ax_paramspace.text(1.6, 0.1, 'eternal expansion', rotation = 5, rotation_mode = 'anchor')
ax_paramspace.text(1.66, -0.11, 'final recollapse', rotation = 4, rotation_mode = 'anchor')

ax_paramspace.text(1.1, -0.5, 'open', rotation = 2, rotation_mode = 'anchor', color="red")
ax_paramspace.text(1.6, -0.5, 'closed', rotation = 2, rotation_mode = 'anchor', color="red")

ax_paramspace.text(2.46, 1.11, 'Decelerating', rotation = 7, rotation_mode = 'anchor', color="blue")
ax_paramspace.text(2.46, 1.5, 'Accelerating', rotation = 7, rotation_mode = 'anchor', color="blue")


ax_paramspace.set_xlim(0,3)
ax_paramspace.set_ylim(-1,3)
ax_paramspace.set_xlabel('$\Omega_M$')
ax_paramspace.set_ylabel('$\Omega_\Lambda$')
ax_paramspace.grid(markevery = 0.1, alpha = 0.2)

def linear(x, m, t):
    return m*x +t

x_continuous = np.linspace(0, 3, 500)
ax_paramspace.plot(x_continuous, linear(x_continuous, t=0, m=1/2), color="blue", ls="-", lw=2, label="Line ")
ax_paramspace.plot(x_continuous, linear(x_continuous, t=1, m=-1), color="red", ls="-", lw=2, label="Line ")

ax_paramspace.errorbar(0.33, 0.6699, yerr=[[0.019], [0.021]], xerr=0.02, fmt="o", markersize=3, color="black")
ax_paramspace.text(0.3, 0.8, 'Solution for our universe', rotation = 0, rotation_mode = 'anchor', color="black")

cursor_paramspace = Cursor(ax_paramspace)
ax_paramspace.get_figure().canvas.mpl_connect('motion_notify_event', cursor_paramspace.on_mouse_move)

click_paramspace = Onclick(0.3, 0.7, ax_evolution, ax_paramspace)
ax_paramspace.get_figure().canvas.mpl_connect('button_press_event', click_paramspace.onclick)

plt.subplots_adjust(left = 0.045, right = 0.98, bottom = 0.055, top = 0.99, wspace = 0.15)
# plt.savefig("omega l vs omega m plot", dpi=300)
plt.show()