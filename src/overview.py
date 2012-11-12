import numpy as np
import matplotlib.pyplot as pl

t, x, y, z, vx, vy, vz, ekin, epot, d, d2 = np.loadtxt("data.txt").transpose()

pl.figure(figsize=(16, 12))

pl.subplot(221)
pl.title("Schwerpunkt")
pl.plot(x,y, 'o-')
pl.xlabel("x")
pl.ylabel("y")
pl.legend()
pl.grid()

pl.subplot(222)
pl.title("Gesamtimpuls")
pl.plot(vx, vy, 'o-')
pl.xlabel("x")
pl.ylabel("y")
pl.legend()
pl.grid()

pl.subplot(223)
pl.title("Energie")
pl.plot(t, ekin, label="Ekin")
pl.plot(t, epot, label="Epot")
pl.plot(t, (ekin+epot)*0.5, 'o-', label="0.5*Eges")
pl.xlabel("Zeit")
pl.ylabel("Energie")
pl.legend(loc="upper left")
pl.grid()

pl.subplot(224)
pl.title("Diffusion")
pl.plot(t, d, label="d")
pl.plot(t, d2, label="d2")
pl.loglog()
pl.xlabel("Zeit")
pl.ylabel("Diffusion")
pl.legend(loc="upper left")
pl.grid()
pl.show()
