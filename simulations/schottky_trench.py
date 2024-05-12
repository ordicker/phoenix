import matplotlib.pyplot as plt
import numpy as np
import meep as mp
from meep.materials import Si, Pt, SiO2


###############################################################################
#                                  Parameters                                 #
###############################################################################
pitch = 5  # pixel pitch in um
depth = 10  # pixel depth in um
d_absorber = 0.5  # absorber_depth in um
d_extention = d_absorber + 2
dpml = 1
fsrc = 0.25  # pulse center frequency = 4um
df = 0.10    # pulse width (in frequency)
nfreq = 100  # number of frequencies at which to compute flux
epsilon = 2
d_trench = 5
w_pt = 0.001
RESOLUTION = 100


###############################################################################
#                                    Setup                                    #
###############################################################################
x_size = pitch+2*dpml
cell = mp.Vector3(x_size, depth+2*(d_extention+dpml), 0)

pml_layers = [mp.PML(dpml)]

geometry = [
    # silicon
    mp.Block(
        mp.Vector3(x_size, depth, mp.inf),
        center=mp.Vector3(0, 0, 0),
        material=mp.Medium(epsilon=11.73)
    ),
    # PtSi interface
    mp.Block(
        mp.Vector3(x_size, 0.1, mp.inf),
        center=mp.Vector3(0, depth/2, 0),
        material=mp.Medium(
            epsilon=epsilon, D_conductivity=2*np.pi*fsrc*5/epsilon)
        # material=mp.Medium(epsilon=-250, D_conductivity=2*np.pi*fsrc*100/250)
    ),
#############################################
#                 trenches                  #
#############################################
    # PtSi
    mp.Block(
        mp.Vector3(2*w_pt+1, w_pt+d_trench, mp.inf),
        center=mp.Vector3(1, depth/2-d_trench/2, 0),
        material=mp.Medium(
            epsilon=epsilon, D_conductivity=2*np.pi*fsrc*5/epsilon)
        # material=mp.Medium(epsilon=-250, D_conductivity=2*np.pi*fsrc*100/250)
    ),
    # SiO2
    mp.Block(
        mp.Vector3(1, d_trench, mp.inf),
        center=mp.Vector3(1, depth/2-d_trench/2, 0),
        #material=SiO2
        material=mp.Medium(epsilon=3.9)
    ),

    # PtSi
    mp.Block(
        mp.Vector3(2*w_pt+1, w_pt+d_trench, mp.inf),
        center=mp.Vector3(-1, depth/2-d_trench/2, 0),
        material=mp.Medium(
            epsilon=epsilon, D_conductivity=2*np.pi*fsrc*5/epsilon)
        # material=mp.Medium(epsilon=-250, D_conductivity=2*np.pi*fsrc*100/250)
    ),
    # SiO2
    mp.Block(
        mp.Vector3(1, d_trench, mp.inf),
        center=mp.Vector3(-1, depth/2-d_trench/2, 0),
        #material=SiO2
        material=mp.Medium(epsilon=3.9)
    )
    
]


sources = [mp.Source(mp.GaussianSource(fsrc, fwidth=df),
                     component=mp.Ez,
                     center=mp.Vector3(0, -(depth+d_extention)/2, 0),
                     size=mp.Vector3(x_size, 0, 0))]


sim = mp.Simulation(
    cell_size=cell,
    boundary_layers=pml_layers,
    geometry=geometry,
    sources=sources,
    resolution=RESOLUTION,
)


# Flux

# reflected flux
refl_fr = mp.FluxRegion(center=mp.Vector3(
    0, -(depth/2+d_extention)+0.5, 0), size=mp.Vector3(x_size, 0, 0))
refl = sim.add_flux(fsrc, df, nfreq, refl_fr)

# transmitted flux
tran_fr = mp.FluxRegion(center=mp.Vector3(
    0, depth/2+d_extention, 0), size=mp.Vector3(x_size, 0, 0))
tran = sim.add_flux(fsrc, df, nfreq, tran_fr)

flux_freqs = mp.get_flux_freqs(refl)

###############################################################################
#                                     Run                                     #
###############################################################################


pt = mp.Vector3(0, depth/2+d_extention, 0)

sim.run(until_after_sources=mp.stop_when_fields_decayed(100, mp.Ez, pt, 1e-4))
# sim.run(until=200)

# get fluxes
refl_flux = mp.get_fluxes(refl)
tran_flux = mp.get_fluxes(tran)
flux_freqs = mp.get_flux_freqs(refl)

###############################################################################
#                                     Plot                                    #
###############################################################################

wl = []
Rs = []
Ts = []
for i in range(nfreq):
    wl = np.append(wl, 1/flux_freqs[i])
    Rs = np.append(Rs, -refl_flux[i])
    Ts = np.append(Ts, tran_flux[i])

plt.figure()
sim.plot2D()
plt.show()
plt.figure()
sim.plot2D(fields=mp.Ez)
plt.show()

if mp.am_master():
    plt.figure()
    plt.plot(wl, Rs, 'bo-', label='reflectance')
    plt.plot(wl, Ts, 'ro-', label='transmittance')
    plt.plot(wl, 1-Rs-Ts, 'go-', label='loss')
    # plt.axis([3.0, 5.0, 0, 1])
    plt.xlabel("wavelength (μm)")
    plt.legend(loc="upper right")
    plt.show()
