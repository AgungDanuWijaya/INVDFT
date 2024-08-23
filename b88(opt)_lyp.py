import mysql.connector
from pyscf import gto, dft,scf,cc
import  math
import numpy as np
import numpy
from pyscf.geomopt.geometric_solver import (optimize)
def annn(rho01,gama):
    tau1=gama**0.5
    x = tau1 / (rho01 + 10E-20) ** (4.0 / 3.0)
    b = 0.0035939430400305887
    b88_g = -1.5 * (3.0 / 4.0 / math.pi) ** (1.0 / 3.0) - b * (x ** 2) / (1.0 + 6.0 * b * x * np.arcsinh(x))

    exc1 = rho01 ** (4.0 / 3.0 ) * b88_g
    return exc1


def eval_xc_gga(xc_code, rho, spin, relativity=0, deriv=2, verbose=None,omega=None):
    rho1 = rho[0]
    rho2 = rho[1]

    rho01, dx1, dy1, dz1 = rho1[:4]
    rho02, dx2, dy2, dz2 = rho2[:4]
    rho01=rho01+1E-20
    rho02 = rho02 + 1E-20
    w1 = rho01 / (rho01 + rho02)
    w2 = rho02 / (rho01 + rho02)
    gamma1 = dx1 ** 2 + dy1 ** 2 + dz1 ** 2
    gamma2 = dx2 ** 2 + dy2 ** 2 + dz2 ** 2



    exc1=[0]*len(rho01)
    exc2 = [0] * len(rho01)
    vgamma_1=[0] * len(rho01)
    vgamma_3 = [0] * len(rho01)

    vrho1=[0] * len(rho02)
    vrho2 = [0] * len(rho02)
    vgamma_2 = [0] * len(rho02)





    for inm in range(len(rho02)) :
        delta1 = (gamma1[inm] ) * 1.0E-25
        delta2 = (gamma2[inm] ) * 1.0E-25
        delta1_ = (rho01[inm]) * 1.0E-9
        delta2_ = (rho02[inm]) * 1.0E-9



        ex1 = annn(rho01[inm], gamma1[inm])
        ex1_ = annn(rho01[inm] + delta1_, gamma1[inm])
        ex1__ = annn(rho01[inm] , gamma1[inm]+delta1)
        vrho1[inm] = (ex1_ - ex1) / (delta1_+1e-250)
        exc1[inm] = (ex1/rho01[inm])*w1[inm]
        vgamma_1[inm] = (ex1__ - ex1)  / (delta1+1e-250)


        ex2 = annn(rho02[inm], gamma2[inm])
        ex2_ = annn(rho02[inm] + delta2_, gamma2[inm])
        ex2__ = annn(rho02[inm] , gamma2[inm]+delta2)
        vrho2[inm] = (ex2_ - ex2) / (delta2_+1e-250)
        exc2[inm] = (ex2/rho02[inm])*w2[inm]
        vgamma_2[inm] = (ex2__ - ex2) / (delta2+1e-250)




    pbe_xc = dft.libxc.eval_xc(',LYP', rho, spin, relativity, deriv,
                               verbose)

    vgamma_ = np.array([vgamma_1, vgamma_3, vgamma_2])
    vgamma_=vgamma_+np.array(pbe_xc[1][1]).T

    vgamma = np.transpose(vgamma_)

    vrho_ = np.array([vrho1, vrho2])
    vrho_ = vrho_ + numpy.array(pbe_xc[1][0]).T

    vrho = np.transpose(vrho_)

    exc1 = np.array([exc1])
    exc2 = np.array([exc2])



    blyp = dft.libxc.eval_xc('BLYP', rho, spin, relativity, deriv,
                             verbose)
    al=1
    bel=0
    exc = np.transpose(al*(exc1 + exc2 + pbe_xc[0])+bel*blyp[0])
    vrho=al*vrho+bel*blyp[1][0]
    vgamma=al*vgamma+bel*blyp[1][1]

    vxc = (vrho, vgamma, None, None)
    fxc = None  # 2nd order functional derivative
    kxc = None  # 3rd order functional derivative
    return exc, vxc, fxc, kxc

mol = gto.Mole()
mol.verbose = 4
mol.atom = """8	0.0000	0.0000	0.1173
1	0.0000	0.7572	-0.4692
1	0.0000	-0.7572	-0.4692"""
mol.charge = 0
mol.spin = 0
mol.basis = "6-31G"
mol.build()
mf = scf.UHF(mol)
mol = optimize(mf, maxsteps=100)
mfl = dft.UKS(mol)
mfl = mfl.define_xc_(eval_xc_gga, xctype='GGA')
mfl.level_shift = 0.5
mfl.grids.level = 0
mfl.conv_tol = 0.000000001
ae =  mfl.kernel()
print(ae)