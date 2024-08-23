import mysql.connector
from pyscf import gto, dft,scf,cc
import  math
import numpy as np
import numpy
from pyscf.geomopt.geometric_solver import (optimize)
mydb = mysql.connector.connect(
      host="127.0.0.1",
      user="mabok_janda",
      password="yut28092018DAM^",
      database="SIAKAD_MIPA_I",
    auth_plugin='mysql_native_password'
    )
mycursor = mydb.cursor()


def sigmoid(x):
    return x

def calcsig(h):
    return [sigmoid(x) for x in h]

def calc(input, w, p):
    hasil = [0] * len(w)
    for i in range(len(w)):
        hasil[i]=w[i] * np.power(np.abs(input), p[i])
        #print(hasil[i],w[i],input,p[i])

    return hasil

def f_ANN(input, c_node, w):
    node = c_node[1:]
    for k in range(len(node)):
        re = [0] * node[k]
        for i in range(len(input)):
            re = numpy.array(re)+numpy.array( calc(input[i], w[0][k][i], w[1][k][i]))
            #print(re)
        input = calcsig(re)
    return np.sum(np.array(input))


def annn(rho,gama) :
    c_node = [2, 4, 4, 1]

    we = [0] * (len(c_node) - 1)
    for i in range(len(c_node) - 1):
        re_i = [0] * c_node[i]
        for j in range(c_node[i]):
            re_j = [0] * c_node[i + 1]
            re_i[j] = re_j
        we[i] = re_i
    we1 = [0] * (len(c_node) - 1)
    for i in range(len(c_node) - 1):
        re_i1 = [0] * c_node[i]
        for j in range(c_node[i]):
            re_j1 = [0] * c_node[i + 1]
            re_i1[j] = re_j1
        we1[i] = re_i1
    param = [0, 0]
    w = [we, we1, param]

    w[0][0][0][0] = -0.5231459874049235
    w[0][0][0][1] = -0.030540649489540914
    w[0][0][0][2] = -0.09429639186978088
    w[0][0][0][3] = 0.10972701401348686
    w[0][0][1][0] = -0.9761554669131317
    w[0][0][1][1] = 0.06859492447986695
    w[0][0][1][2] = -0.5397096625582933
    w[0][0][1][3] = 0.16447375803143013
    w[0][1][0][0] = -0.5868412991249858
    w[0][1][0][1] = -0.24573924316478646
    w[0][1][0][2] = 1.5779277653547352
    w[0][1][0][3] = -0.13705751788915957
    w[0][1][1][0] = -0.23826207119704734
    w[0][1][1][1] = -0.9417828384414029
    w[0][1][1][2] = 0.4244697519081852
    w[0][1][1][3] = -0.26149405924438107
    w[0][1][2][0] = -0.5533676429423282
    w[0][1][2][1] = -0.7846519040746585
    w[0][1][2][2] = -0.1666278176984269
    w[0][1][2][3] = -0.5539037779813359
    w[0][1][3][0] = 0.3695554049737332
    w[0][1][3][1] = 0.8343855839826453
    w[0][1][3][2] = 0.893308919468841
    w[0][1][3][3] = 0.17563969705577887
    w[0][2][0][0] = -0.2577297792566702
    w[0][2][1][0] = 0.16936591984251423
    w[0][2][2][0] = -0.33640454102400696
    w[0][2][3][0] = 0.08913161477370546
    w[1][0][0][0] = 1.226274836345858
    w[1][0][0][1] = 1.203573283118754
    w[1][0][0][2] = 1.3809659484225714
    w[1][0][0][3] = 1.6759323050531303
    w[1][0][1][0] = 1.14280786176743
    w[1][0][1][1] = 1.053257219130628
    w[1][0][1][2] = 1.00336325688442
    w[1][0][1][3] = 1.1943411315577002
    w[1][1][0][0] = 0.7869354113931101
    w[1][1][0][1] = 0.6479703761752541
    w[1][1][0][2] = 1.1911171305606758
    w[1][1][0][3] = 0.6859457595092219
    w[1][1][1][0] = 1.570396970477044
    w[1][1][1][1] = 0.9920292762948137
    w[1][1][1][2] = 1.9382624009708072
    w[1][1][1][3] = 0.7314935020709852
    w[1][1][2][0] = 0.9307430701776422
    w[1][1][2][1] = 1.1730847661861914
    w[1][1][2][2] = 1.0278487424837677
    w[1][1][2][3] = 1.1461911007746244
    w[1][1][3][0] = 0.8170016901998531
    w[1][1][3][1] = 1.1065434713977282
    w[1][1][3][2] = 1.2265963844092267
    w[1][1][3][3] = 1.456937186086773
    w[1][2][0][0] = 1.0726197765814958
    w[1][2][1][0] = 1.2045121148939764
    w[1][2][2][0] = 0.8325032587971499
    w[1][2][3][0] = 1.0415099478614964
    w[2][0] = 1.0451240146757117
    w[2][1] = 0.1087130794657298


    input = [rho,rho ** w[2][0] * gama ** w[2][1]]
    if (rho < 1E-11):
        input = [0, 0]
    if (gama < 1E-11):
        input = [0, 0]
    if (gama < 1E-11 and rho < 1E-11):
        input = [0, 0]


    a = f_ANN(input, c_node, w)
    return a


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
    al=0.7
    bel=0.3
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