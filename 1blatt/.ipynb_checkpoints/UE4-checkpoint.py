
# coding: utf-8

# In[1]:


import netgen.gui
from ngsolve import *
get_ipython().run_line_magic('gui', 'tk')


# In[2]:


def L2Norm(gfu, mesh):
    return sqrt(Integrate(gfu*gfu, mesh))
def H1Norm(gfu, mesh):
    return sqrt(Integrate(grad(gfu)*grad(gfu), mesh))


# In[3]:


from netgen.geom2d import *
geom = SplineGeometry()
pnts = [(0,0), (2,0), (2,1), (1,1), (1,2), (0,2)]
nums = [geom.AppendPoint(*p) for p in pnts]
lines = [(nums[0], nums[1], 10, 1, 0, "b"),
         (nums[1], nums[2], 10, 1, 0, "rb"),
         (nums[2], nums[3], 10, 1, 0, "ul"),
         (nums[3], nums[4], 10, 1, 0, "ru"),
         (nums[4], nums[5], 10, 1, 0, "ul"),
         (nums[5], nums[0], 10, 1, 0, "l")]

for p0, p1, bn, ml, mr, name in lines:
    geom.Append( [ "line", p0, p1 ],  bc=bn, leftdomain=ml, rightdomain=mr)
    #geom.Append( [ "line", p0, p1 ], bc=bn, leftdomain=ml, rightdomain=mr)

mesh = Mesh(geom.GenerateMesh(maxh=0.2))
Draw(mesh)


# In[4]:


#def FEmSpace
fes = H1(mesh, order=3, dirichlet =".*")
a = BilinearForm(fes)
f = LinearForm(fes)
gfu = GridFunction(fes)


# In[5]:


fes.ndof
print(mesh.GetBoundaries())


# In[6]:


#def test and trial-functions
u = fes.TrialFunction()
v = fes.TestFunction()
a += SymbolicBFI (grad(u)*grad(v))
funcf = 1
f += SymbolicLFI (funcf*v)


# In[7]:


a.Assemble()
f.Assemble()


# In[8]:


gfu.vec.data = a.mat.Inverse(freedofs = fes.FreeDofs()) * f.vec
Draw (gfu, mesh, "temperature")
Draw (-grad(gfu), mesh, "grad_temp")


# In[10]:


error_L2 = []
error_H1 = []
funcf = 1
for i in range(1,9):
    fes1 = H1(mesh, order=i, dirichlet =".*")
    a1 = BilinearForm(fes1)
    f1 = LinearForm(fes1)
    gfu1 = GridFunction(fes1)
    u1 = fes1.TrialFunction()
    v1 = fes1.TestFunction()
    a1 += SymbolicBFI (grad(u1)*grad(v1))
    f1 += SymbolicLFI (funcf*v1)
    a1.Assemble()
    f1.Assemble()
    gfu1.vec.data = a1.mat.Inverse(freedofs = fes1.FreeDofs()) * f1.vec
    
    fes2 = H1(mesh, order=i+2, dirichlet =".*")
    a2 = BilinearForm(fes2)
    f2 = LinearForm(fes2)
    gfu2 = GridFunction(fes2)
    u2 = fes2.TrialFunction()
    v2 = fes2.TestFunction()
    a2 += SymbolicBFI (grad(u2)*grad(v2))
    f2 += SymbolicLFI (funcf*v2)
    a2.Assemble()
    f2.Assemble()
    gfu2.vec.data = a2.mat.Inverse(freedofs = fes2.FreeDofs()) * f2.vec
    
    error_L2.append(L2Norm(gfu1-gfu2, mesh))
    error_H1.append(L2Norm(grad(gfu1)-grad(gfu2), mesh))


# In[ ]:


print(sqrt(Integrate((gfu1-gfu2)*(gfu1-gfu2),mesh)))

