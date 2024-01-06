import matplotlib.pyplot as plt
import numpy as np
import meep as mp

cell = mp.Vector3(12, 24, 0)

geometry = [
    mp.Block(
        mp.Vector3(5, 10, mp.inf),
        center=mp.Vector3(0, 5, 0),
        material=mp.Medium(epsilon=3.4, D_conductivity=2*np.pi*0.42*0.5/3.4),
    )
]

#sources = [
#    mp.Source(
#        mp.ContinuousSource(frequency=0.15), component=mp.Ez, center=mp.Vector3(0, -5)
#    )
#]


fsrc = 1.0 # frequency of eigenmode or constant-amplitude source
bnum = 1    # band number of eigenmode
rot_angle = np.radians(90)
k_point = mp.Vector3(fsrc).rotate(mp.Vector3(z=1), rot_angle)


sources = [mp.EigenModeSource(src=mp.ContinuousSource(fsrc),
                              center=mp.Vector3(0,-8,0),
                              size=mp.Vector3(x=12),
                              direction=mp.AUTOMATIC if rot_angle == 0 else mp.NO_DIRECTION,
                              eig_kpoint=k_point,
                              eig_band=1,
                              eig_parity=mp.EVEN_Y+mp.ODD_Z if rot_angle == 0 else mp.ODD_Z,
                              eig_match_freq=True)]




pml_layers = [mp.PML(2.0)]

RESOLUTION = 10

sim = mp.Simulation(
    cell_size=cell,
    boundary_layers=pml_layers,
    geometry=geometry,
    sources=sources,
    resolution=RESOLUTION,
)

sim.run(until=200)

eps_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric)
plt.figure()
plt.imshow(eps_data.transpose(), interpolation="spline36", cmap="binary")
plt.axis("off")
plt.show()

ez_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Ez)
plt.figure()
plt.imshow(eps_data.transpose(), interpolation="spline36", cmap="binary")
plt.imshow(ez_data.transpose(), interpolation="spline36", cmap="RdBu", alpha=0.9)
plt.axis("off")
plt.show()
