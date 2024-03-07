import matplotlib.pyplot as plt
import numpy as np
import meep as mp

###########################################################
#                       Parameters                        #
###########################################################
DX = 5   # simulation X size in um
DY = 10  # simulation Y size in um
dpml = 2
fsrc = 0.25  # pulse center frequency = 4um
df = 0.15    # pulse width (in frequency)
nfreq = 100  # number of frequencies at which to compute flux

RESOLUTION = 10

#######################################################
#                     Setup                           #
#######################################################
x_size = DX+2*dpml
cell = mp.Vector3(x_size, DY+2*dpml, 0)

pml_layers = [mp.PML(dpml)]

geometry = [
    # background n=1.0
    mp.Block(
        mp.Vector3(x_size, DY, mp.inf),
        center=mp.Vector3(0, 0, 0),
        material=mp.Medium(index=1.0)
    ),
    # n index = 2.3
    mp.Block(
        mp.Vector3(x_size, 3.0, mp.inf),
        center=mp.Vector3(0, 0, 0),
        material=mp.Medium(index=2.2)
    )
]


sources = [mp.Source(mp.GaussianSource(fsrc, fwidth=df),
                     component=mp.Ez,
                     center=mp.Vector3(0, -DY/2, 0),
                     size=mp.Vector3(x_size, 0, 0))]


# Run without structure
geometry_without = [
    # background n=1.0
    mp.Block(
        mp.Vector3(x_size, DY, mp.inf),
        center=mp.Vector3(0, 0, 0),
        material=mp.Medium(index=1.0)
    )
]
sim = mp.Simulation(
    cell_size=cell,
    boundary_layers=pml_layers,
    geometry=geometry_without,
    sources=sources,
    resolution=RESOLUTION,
)
# Flux

# reflected flux
refl_fr = mp.FluxRegion(center=mp.Vector3(
    0, -(DY/2)+0.5, 0), size=mp.Vector3(x_size, 0, 0))
refl_without = sim.add_flux(fsrc, df, nfreq, refl_fr)

# transmitted flux
tran_fr = mp.FluxRegion(center=mp.Vector3(
    0, DY/2-0.5, 0), size=mp.Vector3(x_size, 0, 0))
tran_without = sim.add_flux(fsrc, df, nfreq, tran_fr)

pt = mp.Vector3(0, DY/2, 0)

#sim.run(until=2000)
sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, pt, 1e-3))

# for normalization run, save flux fields data for reflection plane
norm_refl_data = sim.get_flux_data(refl_without)

# save incident power for transmission plane
norm_tran_flux = mp.get_fluxes(tran_without)

#################################################
#                  reset                        #
#################################################
sim.reset_meep()

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
    0, -(DY/2)+0.5, 0), size=mp.Vector3(x_size, 0, 0))
refl = sim.add_flux(fsrc, df, nfreq, refl_fr)

# transmitted flux
tran_fr = mp.FluxRegion(center=mp.Vector3(
    0, DY/2-0.5, 0), size=mp.Vector3(x_size, 0, 0))
tran = sim.add_flux(fsrc, df, nfreq, tran_fr)

flux_freqs = mp.get_flux_freqs(refl)

sim.load_minus_flux_data(refl, norm_refl_data)

#####################################################
#                         Run                       #
#####################################################


pt = mp.Vector3(0, DY/2, 0)

#sim.run(until=2000)
sim.run(until_after_sources=mp.stop_when_fields_decayed(100, mp.Ez, pt, 1e-4))


# get fluxes
refl_flux = mp.get_fluxes(refl)
tran_flux = mp.get_fluxes(tran)
flux_freqs = mp.get_flux_freqs(refl)

###################################################
#                       Plot                      #
###################################################

wl = []
Rs = []
Ts = []
for i in range(nfreq):
    wl = np.append(wl, 1/flux_freqs[i])
    Rs = np.append(Rs, -refl_flux[i]/norm_tran_flux[i])
    Ts = np.append(Ts, tran_flux[i]/norm_tran_flux[i])

if mp.am_master():
    plt.figure()
    plt.plot(wl, Rs, 'bo-', label='reflectance')
    plt.plot(wl, Ts, 'ro-', label='transmittance')
    plt.plot(wl, 1-Rs-Ts, 'go-', label='loss')
    # plt.axis([3.0, 5.0, 0, 1])
    plt.xlabel("wavelength (μm)")
    plt.legend(loc="upper right")
    plt.show()

plt.figure()
sim.plot2D()
plt.show()
plt.figure()
sim.plot2D(fields=mp.Ez)
plt.show()
