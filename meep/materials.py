from meep.materials import Si, Pt
import numpy as np
import matplotlib.pyplot as plt

wvl_min = 3.0 # units of μm
wvl_max = 5.0 # units of μm
nwvls = 21
wvls = np.linspace(wvl_min, wvl_max, nwvls)

Si_epsilon = np.array([Si.epsilon(1/w)[0][0] for w in wvls])
 
plt.subplot(1,2,1)
plt.plot(wvls,np.real(Si_epsilon),'bo-')
plt.xlabel('wavelength (μm)')
plt.ylabel('real(ε)')
 
plt.subplot(1,2,2)
plt.plot(wvls,np.imag(Si_epsilon),'ro-')
plt.xlabel('wavelength (μm)')
plt.ylabel('imag(ε)')
 
plt.suptitle('Si from Meep materials library')
plt.subplots_adjust(wspace=0.4)
plt.show()

Pt_epsilon = np.array([Pt.epsilon(1/w)[0][0] for w in wvls])

plt.subplot(1,2,1)
plt.plot(wvls,np.real(Pt_epsilon),'bo-')
plt.xlabel('wavelength (μm)')
plt.ylabel('real(ε)')

plt.subplot(1,2,2)
plt.plot(wvls,np.imag(Pt_epsilon),'ro-')
plt.xlabel('wavelength (μm)')
plt.ylabel('imag(ε)')

plt.suptitle('Pt from Meep materials library')
plt.subplots_adjust(wspace=0.4)
plt.show()

