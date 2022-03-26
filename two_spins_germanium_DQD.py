import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

sq2 = np.sqrt(2)
ueVtoMHz = (1.6e-19*1e-6)/(1e6*6.62e-34)
GHztoueV = (1e9*6.62e-34)/(1.6e-19*1e-6)


def solve_states_ST(es_in, tc_in, B):

    Sum_g = 0.25 + 0.15 # q1 GHz and q2 GHz 
    Del_g = 0.25 - 0.15 # q1 GHz and q2 GHz
    alpha = 1
    e_anti = 0
    #B in  Tesla (?)
    #tc_in in GHz
    # es: AWG pulse, needs converted by the lever arm 
    tc = tc_in*GHztoueV
    tso = 0.1*GHztoueV
    sEz = 0.5* Sum_g *  57.88 * B
    dEz = 0.5* Del_g *  57.88 * B
    
    if type(es_in) is float:
        es = [es_in]
    else:
        es = es_in
    evs = np.zeros((5,len(es)))
    evecs = np.zeros((5,5,len(es)))          
    for i in range(len(es)):
        e = (es[i]-e_anti)*alpha*1e3
        # { S20, S, up-up, T0, down-down  }
        #  S = (up-down  - down-up) / sqrt2
        H = np.array([[    -e,  sq2*tc, tso,   0, tso],
                      [sq2*tc,       0,   0, dEz,   0],
                      [   tso,       0, sEz,   0,   0],
                      [     0,     dEz,   0,   0,   0],
                      [   tso,       0,   0,   0,-sEz],                    
        ])
        ev, evec = np.linalg.eigh(H)
        evs[:,i] = ev
        evecs[:,:,i] = evec
    return evs/GHztoueV, evecs


def gs_es(es_in, tc, B):
    evs, _ = solve_states_ST(es_in, tc, B)
    return evs[1,:]-evs[0,:]
    

#from scipy.optimize import fsolve
#e_so = fsolve(gs_es,[20], args=(0.9, 50))[0]

#es = np.linspace(0,21,2101)
#evs = solve_states_ST(es, 0.9, 50 )[0]
#plt.figure()
#for i in range(5):
#    plt.plot(es, evs[i,:]  )

#%%

from matplotlib.widgets import Slider, Button



fig, ax = plt.subplots(figsize=(6,6))
ax.set_xlabel('detuning'+' (meV)')
ax.set_ylabel('Energies (GHz)')
ax.set_title('Energy levels of two spins \n in a germainum double quantum dot ')


plt.subplots_adjust(left=0.2, bottom=0.35, top=0.9)
x = np.linspace(-5, 5, 1001)

B0 = 0.5
tc0 = 2  # unit GHz

evs = solve_states_ST(x, tc0, B0 )[0]
ls = [ ax.plot(x, evs[i,:], lw=1)[0] for i in range(5)]
#l, = plt.plot(x, mycurve(x, x00, alpha0, tc0), lw=1, c='r')
ax.set_ylim(-5, 5)
ax.margins(x=0)

axcolor = 'lightgoldenrodyellow'

ax_B = plt.axes([0.35, 0.15, 0.5, 0.03], facecolor=axcolor)
ax_tc = plt.axes([0.35, 0.20, 0.5, 0.03], facecolor=axcolor)

s_B = Slider(ax_B, 'B (mT)', 0.0, 1, valinit=B0)
s_tc = Slider(ax_tc, r'$t_c$'+' (GHz)', 0.0, 10.0, valinit=tc0)


def update(val):

    B = s_B.val
    tc = s_tc.val
    
    evs = solve_states_ST(x, tc, B )[0]
    for i, l in enumerate(ls):
        l.set_ydata( evs[i,:]   )
    fig.canvas.draw_idle()


s_B.on_changed(update)
s_tc.on_changed(update)

resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')


def reset(event):

    s_B.reset()
    s_tc.reset()
button.on_clicked(reset)




plt.show()






