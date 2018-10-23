from ngsolve import *



from netgen.geom2d import SplineGeometry

#========================================
def TopRectangle(geom):
    top_pnts = [(+2.00,-0.35), (+2.00,+2.00), (-2.00,+2.00), (-2.00,-0.35)]
    top_nums = [geom.AppendPoint(*p) for p in top_pnts]

    lines =  [ (top_nums[0], top_nums[1],   10,    1,  0),
               (top_nums[1], top_nums[2],   10,    1,  0),
               (top_nums[2], top_nums[3],   10,    1,  0),
               (top_nums[3], top_nums[0],   10,    1,  0) ]    
    
    for p0,p1,bn,ml,mr  in  lines:                                          
        geom.Append( [ "line", p0, p1 ],  bc=bn, leftdomain=ml, rightdomain=mr)
    return geom

#========================================
def BottomRectangle(geom):
    top_pnts = [(+2.00,-0.35), (+2.00,-2.00), (-2.00,-2.00), (-2.00,-0.35)]
    top_nums = [geom.AppendPoint(*p) for p in top_pnts]

    lines =  [ (top_nums[0], top_nums[1],   10,    0,  2),
               (top_nums[1], top_nums[2],   10,    0,  2),
               (top_nums[2], top_nums[3],   10,    0,  2),
               (top_nums[3], top_nums[0],   10,    0,  2) ]    
    
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
def MakeMesh():
    
    geometry = SplineGeometry()
    geometry = TopRectangle(geometry)
    geometry = Circle(geometry, 0.25)
    geometry = BottomRectangle(geometry)

    return Mesh(geometry.GenerateMesh())
#========================================
