from ngsolve import *
from numpy.polynomial.legendre import leggauss as gauss_points
import matplotlib.pyplot as plt

def gen_mesh(h):
    from netgen.geom2d import unit_square
    return Mesh(unit_square.GenerateMesh(maxh=h))

def Solve(mesh, order, rhs):
    V = H1(mesh, order=order, dirichlet='left|right')
    u,v = V.TnT()
    a = BilinearForm(V)
    a += SymbolicBFI(grad(u)*grad(v))
    a.Assemble()
    f = LinearForm(V)
    f += SymbolicLFI(rhs*v)
    f.Assemble()
    sol = GridFunction(V)
    inv = a.mat.Inverse(V.FreeDofs(), inverse='sparsecholesky')
    sol.vec.data = inv * f.vec
    return sol, V, a, f

# 5.a: directly, by hand
def CalcFlux_5a(mesh, fun):
    ileft = 0
    iright = 0
    # [-1,1] intrule -> [0,1] intrule
    pw = [(0.5*(1+p), 0.5*w) for p,w in zip(*gauss_points(5))]
    for p, w in pw:
        ileft = ileft - w * fun(mesh(0,p,0))
        iright = iright + w * fun(mesh(1,p,0))
    F = {'left':ileft, 'right':iright}
    return F

# 5.b: directly, use "Integrate"
def CalcFlux_5b(mesh, fun):
    # crashes with less recent NGSolve versions!!
    L = Integrate(fun, mesh, BND, region_wise=True)
    return { bcname : L[bcnum] for bcnum, bcname in enumerate(mesh.GetBoundaries()) if bcname in ['left', 'right']}
# 5.b  workaround
def CalcFlux_5b_alt(mesh, fun):
    return { bcn : Integrate(fun, mesh, BND, definedon=mesh.Boundaries(bcn)) for bcn in ['left', 'right'] }
    
# 6.a: indirectly, use "Integrate"
def CalcFlux_6a(mesh, sol, w, rhs):
    F = {}
    for bcn, wf in [('left',1-x), ('right',x)]:
        w.Set(wf)
        F[bcn] = Integrate(grad(sol)*grad(w) - rhs*w, mesh, VOL)
    return F

# 6.b: indirectly, use matrix-vector operations
def CalcFlux_6b(mesh, sol, w, rhs, a, f):
    F = {}
    v = sol.vec.CreateVector()
    for bcn, wf in [('left',1-x), ('right',x)]:
        w.Set(wf)
        v.data = a.mat * sol.vec - f.vec
        F[bcn] = InnerProduct(v, w.vec)
    return F
        
do_calc = False                      
save_file = 'ue2_ex45.data'          
hs = [0.5**k for k in range(2,9)]    
ps = list(range(1,4))                
rhs = x
fluxes = { n : {} for n in ['5a','5b','6a','6b'] }
ex_fluxl = -1/6
ex_fluxr = -1/3
ex_flux = ex_fluxr + ex_fluxl

def save_results():
    import pickle
    pickle.dump( (hs, ps, fluxes), open(save_file, 'wb'))

def load_results():
    import pickle
    return pickle.load(open(save_file, 'rb'))

def do_calculations():
    for h in hs:
        mesh = gen_mesh(h)
        for p in ps:
            sol, V, a, f = Solve(mesh, p, rhs)
            dudxn = GridFunction(V)
            dudxn.Set(IfPos(x-0.5, 1, -1) * grad(sol)[0])
            w = GridFunction(V)
            fluxes['5a'][(h,p)] = CalcFlux_5a(mesh, dudxn)
            fluxes['5b'][(h,p)] = CalcFlux_5b_alt(mesh, dudxn)
            fluxes['6a'][(h,p)] = CalcFlux_6a(mesh, sol, w, rhs)
            fluxes['6b'][(h,p)] = CalcFlux_6b(mesh, sol, w, rhs, a, f)
            
            
def do_plot():
    plt.figure()
    plt.suptitle('Exercise 2, problems 5 & 6 ')
    colors = ['b','r','g', 'c', 'm', 'y', 'k']
    spc = 1
    flux_err = lambda var, h, p : abs(fluxes[var][(h,p)]['left'] + fluxes[var][(h,p)]['right'] - ex_flux)
    for var in fluxes:
        plt.subplot(2,2,spc)
        plt.title('variant '+var)
        spc = spc+1
        c = 0
        for p in ps:
            plt.semilogy(hs, [flux_err(var,h,p) for h in hs], label='p = '+str(p), color=colors[c])
            c = (c+1)%len(colors)
        plt.legend()
            
    plt.show()
    
if do_calc:
    do_calculations()
    save_results()
else:
    hs, ps, fluxes = load_results()
do_plot()
