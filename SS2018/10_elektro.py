from ngsolve import *
#from mygeom import MakeMesh

from netgen.geom2d import SplineGeometry





def TopRectangle(geom):

    # Make top subdomain (material 1)  \Omega_1:
    #          +-------------------------------------------+ (2,2)
    #          | 2        	                            1  |
    #          |       Subdomain 1                         |
    #          |                                           |
    #          |                                           |
    #          |                                           |
    #          |                                           |
    #          |                                           |
    #          |                                           |
    #          |                                           |
    #          |                                           |
    #          | 3                                      0  |
    #          +-------------------------------------------+
    #(-2,-0.35)|         Subdomain 2                       | (2,-0.35)
    #          |                                           |

    
    # First, make a list of all points:
    #             point 0,       point 1,       point 2,       point 3
    top_pnts = [ (+2.00,-0.35), (+2.00,+2.00), (-2.00,+2.00), (-2.00,-0.35)]

    # Add them to given geometry and collect the vertex numbers
    top_nums = [geom.AppendPoint(*p) for p in top_pnts]

    # Make a list of line segments in this format:
    #           beginning point number,
    #           |            end point number,
    #           |            |             boundary condition number,
    #           |            |             |     domain on left side,
    #           |            |             |     |   domain on right side
    #           |            |             |     |   |
    #           V            V             V     V   V
    lines  = [ (top_nums[0], top_nums[1],  10,   1,  0),
               (top_nums[1], top_nums[2],  10,   1,  0),
               (top_nums[2], top_nums[3],  10,   1,  0),
               (top_nums[3], top_nums[0],  10,   1,  2) ]

    for p0,p1,bn,ml,mr  in  lines:
        geom.Append([ "line", p0, p1 ],
                    bc=bn, leftdomain=ml, rightdomain=mr )

    return (geom, top_nums)
#---------------------------------------------------------------

def PunchCircle(geom):
    #
    #          +-------------------------------------------+
    #          |         	                               |
    #          |       Subdomain 1                         |
    #          |                                           |
    #          |                                           |
    #          |                                           |
    #          |                                           |
    #          |                3   2   1                  |
    #          |                   ---      Disk of radius |
    #          |                4(  o  )0   0.25 centered  |
    #          |                   ---      at o=origin    |
    #          |                5   6   7                  |
    #          +-------------------------------------------+

    disc_pnts = [(+0.25, +0.00),   # point 0
                 (+0.25, +0.25),   # point 1
                 (+0.00, +0.25),   # point 2
                 (-0.25, +0.25),   # point 3
                 (-0.25, +0.00),   # point 4
                 (-0.25, -0.25),   # point 5
                 (+0.00, -0.25),   # point 6
                 (+0.25, -0.25) ]  # point 7

    disc_nums = [geom.AppendPoint(*p) for p in disc_pnts]

    # Make SPLINE curve in this tuple format:
    #   (beginning control point, middle control point, final control point,
    #       boundary condition number, domain on left side,  domain on right side)
    #
    curves = [ (disc_nums[0], disc_nums[1], disc_nums[2],  1,0,1),
               (disc_nums[2], disc_nums[3], disc_nums[4],  1,0,1),
               (disc_nums[4], disc_nums[5], disc_nums[6],  1,0,1),
               (disc_nums[6], disc_nums[7], disc_nums[0],  1,0,1)  ]

    for p0,p1,p2,bc,left,right in curves:
        geom.Append( ["spline3", p0,p1,p2],
                     bc=bc, leftdomain=left, rightdomain=right)
    return (geom, disc_nums)
#---------------------------------------------------------------

# Input:
#   "topn" contains the global numbers of the top subdomain
#   "geom" the global geometry to which objects are added

def AddBottomRectangle(geom,topn):   
    
    #          +-------------------------------------------+ (2,2)
    #          | top[2]     	                    top[1] |
    #          |       Subdomain 1                         |
    #          |                                           |
    #          |                                           |
    #          |                                           |
    #          |                                           |
    #          |                                           |
    #          |                                           |
    #          |                                           |
    #          |                                           |
    #          | top[3]                             top[0] |
    #          +-------------------------------------------+ (2,-0.35)
    #          |                                           |
    #          |                                           |
    #          |                                           |
    #          |                                           |
    #          |                                           |
    #          |                                           |
    #          |bot[0]     Subdomain 2              bot[1] |
    #          +-------------------------------------------+ (2,-2)
    #     (-2,-2)

    
    bot_pnts = [ (-2.00, -2.00),  (+2.00, -2.00) ]
    botn = [geom.AppendPoint(*p) for p in bot_pnts]

    # Format:  (p0,     p1,        bn,   ml, mr)
    lines  = [ (topn[3], botn[0],    10,   2,   0),
               (botn[0], botn[1],    10,   2,   0),
               (botn[1], topn[0],    10,   2,   0) ]

    for p0,p1,bn,ml,mr  in  lines:
        geom.Append([ "line", p0, p1 ],
                    bc=bn, leftdomain=ml, rightdomain=mr )

    return (geom, botn)
#---------------------------------------------------------------

def PunchPlate(geom):

    #          | top[3]                             top[0] |
    #          +-------------------------------------------+ 
    #          |            +----------------+ (0.5,-0.5)  |
    #          |            +----------------+             |
    #          |   (-0.5,-0.6)                (0.5,-0.6)   |
    #          |                                           |
    #          |                                           |
    #          |                                           |
    #          | bot[0]     Subdomain 2             bot[1] |
    #          +-------------------------------------------+ 
    #     (-2,-2)

    plate_pnts = [ (+0.50, -0.60),  (+0.50, -0.50),
                   (-0.50, -0.50),  (-0.50, -0.60) ]
               
    plate_nums = [geom.AppendPoint(*p) for p in plate_pnts]

    # Format:  (    p0,          p1,          bn,   ml, mr)
    lines  = [ (plate_nums[0], plate_nums[1],  2,    0, 2),
               (plate_nums[1], plate_nums[2],  2,    0, 2),
               (plate_nums[2], plate_nums[3],  2,    0, 2),
               (plate_nums[3], plate_nums[0],  2,    0, 2) ]

    for p0,p1,bn,ml,mr  in  lines:
        geom.Append([ "line", p0, p1 ],
                    bc=bn, leftdomain=ml, rightdomain=mr )
    return geom
#---------------------------------------------------------------

def MakeMesh() :

    geometry = SplineGeometry()
    geometry, top = TopRectangle(geometry)
    geometry, dsc = PunchCircle(geometry)
    geometry, bot = AddBottomRectangle(geometry,top)
    geometry = PunchPlate(geometry)

    return  Mesh( geometry.GenerateMesh(maxh=1) )



"""
#========================================
def TopRectangle(geom):
    top_pnts = [(+2.00,-0.35), (+2.00,+2.00), (-2.00,+2.00), (-2.00,-0.35)]
    top_nums = [geom.AppendPoint(*p) for p in top_pnts]

    lines =  [ (top_nums[0], top_nums[1],   10,    1,  0),
               (top_nums[1], top_nums[2],   10,    1,  0),
               (top_nums[2], top_nums[3],   10,    1,  0),
               (top_nums[3], top_nums[0],   10,    1,  2) ]    
    
    for p0,p1,bn,ml,mr  in  lines:                                          
        geom.Append( [ "line", p0, p1 ],  bc=bn, leftdomain=ml, rightdomain=mr)
    return geom,(top_nums[3], top_nums[0])

#========================================
def BottomRectangle(geom,bottompoints):
    bottom_pnts = [(+2.00,-2.00), (-2.00,-2.00)]
    bottom_nums = [geom.AppendPoint(*p) for p in bottom_pnts]

    lines =  [ (bottompoints[1], bottom_nums[0],   10,    0,  2),
               (bottom_nums[0], bottom_nums[1],   10,    0,  2),
               (bottom_nums[1], bottompoints[0],   10,    0,  2)]
             #  (bottom_nums[3], bottom_nums[0],   10,    0,  2) ]    
    
    for p0,p1,bn,ml,mr  in  lines:                                          
        geom.Append( [ "line", p0, p1 ],  bc=bn, leftdomain=ml, rightdomain=mr)
    return geom
#========================================
def Circle(geom,r):
    
    #              point 0   point 1   point 2   point 3  point 4    point 5   point 6   point 7
    disc_pnts = [ ( r,  0), ( r,  r), ( 0,  r), (-r,  r), (-r,  0), (-r, -r), ( 0, -r), ( r, -r) ]
    
    disc_pnums = [geom.AppendPoint(*p) for p in disc_pnts]
    
    curves = [ (disc_pnums[0], disc_pnums[1], disc_pnums[2],  1,0,1),
               (disc_pnums[2], disc_pnums[3], disc_pnums[4],  1,0,1),
               (disc_pnums[4], disc_pnums[5], disc_pnums[6],  1,0,1),
               (disc_pnums[6], disc_pnums[7], disc_pnums[0],  1,0,1)  ]

    for p0,p1,p2,bc,left,right in curves:
        geom.Append( ["spline3", p0,p1,p2],bc=bc, leftdomain=left, rightdomain=right)
    return geom
#========================================
def Stick(geom):
    stick_pnts = [(-0.5,-0.5),(+0.5,-0.5),(0.5,-0.6),(-0.5,-0.6)]
    stick_nums = [geom.AppendPoint(*p) for p in stick_pnts]

    lines =  [ (stick_nums[0], stick_nums[1],   10,    2,  0),
               (stick_nums[1], stick_nums[2],   10,    2,  0),
               (stick_nums[2], stick_nums[3],   10,    2,  0),
               (stick_nums[3], stick_nums[0],   10,    2,  0) ]    
    for p0,p1,bn,ml,mr in lines:
        geom.Append(["line", p0,p1], bc=bn, leftdomain=ml,rightdomain=mr)
    return geom
#========================================
def MakeMesh():
    geometry = SplineGeometry()
    geometry, bottompoints = TopRectangle(geometry)
    geometry = Circle(geometry, 0.25)
    geometry = BottomRectangle(geometry,bottompoints)
    geometry = Stick(geometry)

    return Mesh(geometry.GenerateMesh())
#========================================
"""



mesh = MakeMesh()

V = H1(mesh, order=4, dirichlet=[1,2])
"""
u = GridFunction(V)
Draw(u)
for i in range(V.ndof):
    u.vec[:] = 0
    u.vec[i] = 1
    Redraw()
    input()
    
"""
epsl = CoefficientFunction([1, 10])#DomainConstantCF([1, 10])

u,v = V.TrialFunction(), V.TestFunction() #new

a = BilinearForm(V, symmetric=True)
#a += Laplace(epsl)
a += SymbolicBFI(epsl*grad(u)*grad(v)) #new

a.Assemble()

extendone = CoefficientFunction([1,0])#DomainConstantCF([1,0])
u1 = GridFunction(V, name="extension")
u1.Set(extendone, BND)
Draw(u1)

u=GridFunction(V, name="solution")
u.vec.data = u1.vec
a.Assemble()
f = u.vec.CreateVector()
f.data = a.mat*u1.vec
u.vec.data -= a.mat.Inverse(V.FreeDofs(),inverse="sparsecholesky")*f
Draw(u)
Draw(-epsl*grad(u),mesh,"E")
