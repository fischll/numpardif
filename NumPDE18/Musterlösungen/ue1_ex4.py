from netgen.geom2d import *
from ngsolve import *

# construct the geometry for the square
def SqGeom():
    coords = [(0,0), (2,0), (2,2), (0,2)]
    geom = SplineGeometry()
    for co in coords:
        geom.AddPoint(*co)
    for i in range(len(coords)):
        geom.Append(["line", i, (i+1)%len(coords)], leftdomain=1, rightdomain=0, bc='outer')
    geom.SetMaterial(1, 'inner')
    return geom

# construct the geometry for the L-shaped domain
def LGeom():
    coords = [(0,0), (2,0), (2,1), (1,1), (1,2), (0,2)]
    geom = SplineGeometry()
    for co in coords:
        geom.AddPoint(*co)
    for i in range(len(coords)):
        geom.Append(["line", i, (i+1)%len(coords)], leftdomain=1, rightdomain=0, bc='outer')
    geom.SetMaterial(1, 'inner')
    return geom

# Parameters
do_calc = True                          # if True, do calculations // if False, load results from file
save_file = 'ue1_ex4.data'              # where to save results to / load them from
hs = [0.5**k for k in range(1,7)]       # mesh size
ps = list(range(1, 8))                 # orders
err_types = {'l2' : lambda a,grada,b,gradb: (a-b)*(a-b) , \
             'h1' : lambda a,grada,b,gradb: (a-b)*(a-b) + (grada-gradb)*(grada-gradb)}
geoms = { 'square' : SqGeom(), 'L-shape' : LGeom() }  # geometries

# RHS from Ex. 4
rhs = {'square' : 1 , 'L-shape' : 1}    # analytic sol in unit_square
ref_sols  = {}
ref_grads  = {}

# # Polynomial exact solution for the square
# rhs = {'square' : x*(2-x)+y*(2-y) , 'L-shape' : 1}
# ref_sols  = {'square': 0.5*x*(2-x)*y*(2-y)}
# ref_grads = {'square': CoefficientFunction(((1-x)*y*(2-y), (1-y)*x*(1-x)))}

# # Non-polynomial exact solution for square ang L-shape
# from math import pi
# rhs = {'square' : 2*pi**2*sin(pi*x)*sin(pi*y), \
#        'L-shape' : 8*pi**2*sin(2*pi*x)*sin(2*pi*y) }
# ref_sols = {'square' : sin(pi*x)*sin(pi*y), \
#             'L-shape' : sin(2*pi*x)*sin(2*pi*y) }
# ref_grads = { 'square' : CoefficientFunction( (pi*cos(pi*x)*sin(pi*y), pi*sin(pi*x)*cos(pi*y)) ), \
#               'L-shape' : CoefficientFunction( (2*pi*cos(2*pi*x)*sin(2*pi*y), 2*pi*sin(2*pi*x)*cos(2*pi*y)) ) }



# given a geometry, a mesh size h and an order for the Finite Element Space,
# solve the laplace equation and return the solution
def Solve(geom, maxh, order, rhs):
    mesh = Mesh(geom.GenerateMesh(maxh=maxh))
    V = H1(mesh, order=order, dirichlet='outer')
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
    return sol

sols = {name : {} for name in geoms }
errors = { err_type : {name : {} for name in geoms } for err_type in err_types}

def do_calculations():
    with TaskManager():
        # compute solutions
        for name, G in geoms.items():
            sols[name] = {}
            for h in hs:
                for p in ps:
                    sols[name][(h,p)] = Solve(geom=G, maxh=h, order=p, rhs=rhs[name])
        # compute errors
        for name, G in geoms.items():
            # use a fine mesh for integration!
            mesh = Mesh(G.GenerateMesh(maxh=hs[-1]))
            # as reference solution, use the one we calculated on the
            # finest mesh with the highest order FEM space
            ref_sol  = ref_sols[name]  if name in ref_sols  else sols[name][(hs[-1], ps[-1])]
            ref_grad = ref_grads[name] if name in ref_grads else grad(ref_sol)
            for h in hs:
                for p in ps:
                    for err_type, fun in err_types.items():
                        s = sols[name][(h,p)]
                        error_cf = fun(s, grad(s), ref_sol, ref_grad)
                        errors[err_type][name][(h,p)] = sqrt(Integrate(error_cf, mesh))

def save_results():
    import pickle
    pickle.dump( (hs, ps, sols, errors), open(save_file, 'wb'))

def load_results():
    import pickle
    return pickle.load(open(save_file, 'rb'))
    
def do_plot():
    import matplotlib.pyplot as plt
    from math import log
    colors = ['b','r','g', 'c', 'm', 'y', 'k']
    c = 0 # position in color array
    n_errs = len(errors)
    lims = {'h1' : [1e-6, 1e0], 'l2' : [1e-9, 1e-1]}
    rates = {'h1' : {'square'  : {'func' : lambda h,p : (h)**(p),         'label' : 'O(h^p)'},
                     'L-shape' : {'func' : lambda h,p : (h)**(0.5*p),     'label' : 'O(h^(p/2))'} },
             'l2' : {'square' :  {'func' : lambda h,p : (h)**(p+1),       'label' : 'O(h^(p+1))'},
                     'L-shape' : {'func' : lambda h,p : (h)**(0.5*(p+1)), 'label' : 'O(h^((p+1)/2))'} } }
    for err_type in errors:
        n_geoms = len(errors[err_type])
        plt.figure()
        plt.suptitle(err_type + '-norm')
        sp_cnt = 1 # position in subplots .. 1-based
        for geom in errors[err_type]:
            # top row: p vs. error
            plt.subplot(2, n_geoms, sp_cnt)
            plt.title(geom + ' - p vs. error')
            #plt.ylim(lims[err_type])
            plt.plot([], [],'--', label=rates[err_type][geom]['label'])
            c = 0
            ylims = [1e5,-1e5]
            for h in hs:
                use_ps = ps[:-1] if h is hs[-1] else ps
                use_errs = [errors[err_type][geom][(h,p)] for p in use_ps]
                ylims = [min(use_errs+[ylims[0]]), max(use_errs+[ylims[1]])]
                plt.semilogy(use_ps, use_errs, label='h = 1/2^'+str(-int(log(h,2))), color=colors[c], marker='d')
                rate_func = rates[err_type][geom]['func']
                plt.semilogy(ps, [0.8 * use_errs[0] / rate_func(h, ps[0]) * rate_func(h,p) for p in ps], '--', color=colors[c])
                c = (c+1) % len(colors)
            plt.ylim([0.5*ylims[0], 2*ylims[1]])
            plt.legend()
            sp_cnt = sp_cnt + n_geoms
            # bottom row: p vs. error
            plt.subplot(2, n_geoms, sp_cnt)
            # plt.ylim(lims[err_type])
            plt.plot([], [],'--', label=rates[err_type][geom]['label'])
            c = 0
            for p in ps[slice(None,None,2)]:
                use_hs = hs[:-1] if p is ps[-1] else hs
                use_errs = [errors[err_type][geom][(h,p)] for h in use_hs]
                plt.loglog(use_hs, use_errs, label='p = '+str(p), color=colors[c], marker='d')
                rate_func = rates[err_type][geom]['func']
                plt.loglog(hs, [0.8 * use_errs[0] / rate_func(hs[0], p) * rate_func(h,p) for h in hs], '--', color=colors[c])
                c = (c+1) % len(colors)
            plt.title(geom + ' - h vs. error')
            plt.legend()
            sp_cnt = 1 + ( (sp_cnt + n_geoms ) % (2*n_geoms) )
    plt.show()

if do_calc:
    do_calculations()
    save_results()
else:
    hs, ps, sols, errors = load_results()
do_plot()
