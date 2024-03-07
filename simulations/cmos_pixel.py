import math
import argparse
import meep as mp
from meep.materials import W, SiO2, cSi

dair = 1.0                         # air gap thickness
dmcl = 1.7                         # micro lens thickness
dsub = 3.0                         # substrate thickness
dpml = 1.0                         # PML thickness

sz = dpml+dair+dmcl+dsub+dpml

lambda_min = 0.7                   # minimum source wavelength
lambda_max = 1.0                   # maximum source wavelength
fmin = 1/lambda_max
fmax = 1/lambda_min
fcen = 0.5*(fmin+fmax)
df = fmax-fmin

def main(args):
    sxy = args.Np*args.a
    cell_size = mp.Vector3(sxy,sxy,sz)

    boundary_layers = [mp.PML(dpml, direction=mp.Z, side=mp.High),
                       mp.Absorber(dpml, direction=mp.Z, side=mp.Low)]

    geometry = []

    if args.substrate:
        geometry = [mp.Sphere(material=SiO2, radius=dmcl, center=mp.Vector3(z=0.5*sz-dpml-dair-dmcl)),
                    mp.Block(material=cSi, size=mp.Vector3(mp.inf,mp.inf,dsub+dpml),
                             center=mp.Vector3(z=-0.5*sz+0.5*(dsub+dpml))),
                    mp.Block(material=W, size=mp.Vector3(mp.inf, args.wire_w, args.wire_h),
                             center=mp.Vector3(0,-0.5*sxy+0.5*args.wire_w,-0.5*sz+dpml+dsub+0.5*args.wire_h)),
                    mp.Block(material=W, size=mp.Vector3(mp.inf, args.wire_w, args.wire_h),
                             center=mp.Vector3(0,+0.5*sxy-0.5*args.wire_w,-0.5*sz+dpml+dsub+0.5*args.wire_h)),
                    mp.Block(material=W, size=mp.Vector3(args.wire_w, mp.inf, args.wire_h),
                             center=mp.Vector3(-0.5*sxy+0.5*args.wire_w,0,-0.5*sz+dpml+dsub+0.5*args.wire_h)),
                    mp.Block(material=W, size=mp.Vector3(args.wire_w, mp.inf, args.wire_h),
                             center=mp.Vector3(+0.5*sxy-0.5*args.wire_w,0,-0.5*sz+dpml+dsub+0.5*args.wire_h))]

    if args.substrate and args.texture:
        for nx in range(args.Np):
            for ny in range(args.Np):
                cx = -0.5*sxy+(nx+0.5)*args.a
                cy = -0.5*sxy+(ny+0.5)*args.a
                geometry.append(mp.Cone(material=SiO2,
                                        radius=0,
                                        radius2=args.cone_r,
                                        height=args.cone_h,
                                        center=mp.Vector3(cx,cy,0.5*sz-dpml-dair-dmcl-0.5*args.cone_h)))

    if args.substrate:
        geometry.append(mp.Block(material=SiO2,
                                 size=mp.Vector3(mp.inf, args.trench_w, args.trench_h),
                                 center=mp.Vector3(0, -0.5*sxy+0.5*args.trench_w, 0.5*sz-dpml-dair-dmcl-0.5*args.trench_h)))
        geometry.append(mp.Block(material=SiO2,
                                 size=mp.Vector3(mp.inf, args.trench_w, args.trench_h),
                                 center=mp.Vector3(0, +0.5*sxy-0.5*args.trench_w, 0.5*sz-dpml-dair-dmcl-0.5*args.trench_h)))
        geometry.append(mp.Block(material=SiO2,
                                 size=mp.Vector3(args.trench_w, mp.inf, args.trench_h),
                                 center=mp.Vector3(-0.5*sxy+0.5*args.trench_w, 0, 0.5*sz-dpml-dair-dmcl-0.5*args.trench_h)))
        geometry.append(mp.Block(material=SiO2,
                                 size=mp.Vector3(args.trench_w, mp.inf, args.trench_h),
                                 center=mp.Vector3(+0.5*sxy-0.5*args.trench_w, 0, 0.5*sz-dpml-dair-dmcl-0.5*args.trench_h)))

    sources = [mp.Source(mp.GaussianSource(fcen, fwidth=df),
                         component=mp.Ex,
                         center=mp.Vector3(z=0.5*sz-dpml-0.5*dair),
                         size=mp.Vector3(sxy,sxy))]

    sim = mp.Simulation(resolution=args.res,
                        cell_size=cell_size,
                        boundary_layers=boundary_layers,
                        geometry=geometry,
                        dimensions=3,
                        split_chunks_evenly=False,
                        k_point=mp.Vector3(),
                        sources=sources)

    nfreq = 50
    refl = sim.add_flux(fcen, df, nfreq, mp.FluxRegion(center=mp.Vector3(0,0,0.5*sz-dpml),size=mp.Vector3(sxy,sxy,0)))
    trans_grid = sim.add_flux(fcen, df, nfreq,
                              mp.FluxRegion(center=mp.Vector3(0,0,-0.5*sz+dpml+dsub+args.wire_h),size=mp.Vector3(sxy,sxy,0)))
    trans_sub_top = sim.add_flux(fcen, df, nfreq,
                                 mp.FluxRegion(center=mp.Vector3(0,0,-0.5*sz+dpml+dsub),size=mp.Vector3(sxy,sxy,0)))
    trans_sub_bot = sim.add_flux(fcen, df, nfreq,
                                 mp.FluxRegion(center=mp.Vector3(0,0,-0.5*sz+dpml),size=mp.Vector3(sxy,sxy,0)))

    sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ex, mp.Vector3(z=-0.5*sz+dpml+0.5*dsub), 1e-9))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-res', type=float, default=30, help='reslution (default: 30 pixels/um)')
    parser.add_argument('-substrate', action='store_true', default=False, help='add the substrate? (default: False)')
    parser.add_argument('-texture', action='store_true', default=False, help='add the texture? (default: False)')
    parser.add_argument('-a', type=float, default=0.4, help='lattice periodicity (default: 0.4 um)')
    parser.add_argument('-cone_r', type=float, default=0.2, help='cone radius (default: 0.2 um)')
    parser.add_argument('-cone_h', type=float, default=0.28247, help='cone height (default: 0.28247 um)')
    parser.add_argument('-wire_w', type=float, default=0.1, help='metal-grid wire width (default: 0.1 um)')
    parser.add_argument('-wire_h', type=float, default=0.2, help='metal-grid wire height (default: 0.2 um)')
    parser.add_argument('-trench_w', type=float, default=0.1, help='trench width (default: 0.1 um)')
    parser.add_argument('-trench_h', type=float, default=2.0, help='trench height (default: 2.0 um)')
    parser.add_argument('-Np', type=int, default=3, help='number of periods in supercell (default: 3)')
    args = parser.parse_args()
    main(args)
  
