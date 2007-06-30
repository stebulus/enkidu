"""Enkidu: generating PostScript/LaTeX figures.

Classes:
    vec         a vector/point (we don't distinguish)
    line        a line
    circle      a circle
    affine      an affine transformation
    figure      a diagram

Angles are always in degrees; 0 is east, 90 is north.
"""
from __future__ import division
from os.path import splitext, extsep
from math import atan2, cos, pi, sin, sqrt

deg = pi/180  # so that cos(ang*deg) means "cosine of ang degrees", etc.

def det(a, b, c, d):
    """The 2x2 determinant, that is, a*d - b*c."""
    return a*d - b*c

def withmin(func, items):
    """The item in items with the minimum value of func."""
    items = iter(items)
    minitem = items.next()
    minval = func(minitem)
    for it in items:
        curval = func(it)
        if minval > curval:
            minval = curval
            minitem = it
    return minitem

def withmax(func, items):
    """The item in items with the maximum value of func."""
    return withmin(lambda x: -func(x), items)

def xcoord(pt):
    """The x coordinate of pt.  (Usually written 'pt.x'.)"""
    return pt.x

def ycoord(pt):
    """The y coordinate of pt.  (Usually written 'pt.y'.)"""
    return pt.y

def leftmost(pts):
    """The leftmost of the points in pts."""
    return withmin(xcoord, pts)

def rightmost(pts):
    """The rightmost of the points in pts."""
    return withmax(xcoord, pts)

def bottommost(pts):
    """The bottommost of the points in pts."""
    return withmin(ycoord, pts)

def topmost(pts):
    """The topmost of the points in pts."""
    return withmax(ycoord, pts)

def outermost(pts):
    """The outermost two of the (assumed collinear) points in pts."""
    pt1 = withmax(distfrom(pts[0]), pts)
    pt2 = withmax(distfrom(pt1), pts)
    return pt1, pt2

class vec(object):
    """A vector or point (we don't distinguish) in the plane.

    Given vecs u and v, the expressions u+v, -v, u-v, 4*v, v*4,
    and u == v have the obvious meanings.  v/2.0 also has the
    obvious meaning; so does v/2, if __future__.division is in
    effect.  (Integer division as in v//2 is not supported.)

    Other properties of and operations on vecs:
        v.x         the x coordinate of v
        v.y         the y coordinate of v
        dot(u, v)   the dot product of u and v
        v.norm      the norm of v
        v.normsq    the norm of v, squared
        v.angle     the angle from the positive x semi-axis to v
                        (if v == vec.zero, then v.angle = 0.0)
        unit(v)     the unit vector in the direction of v
        rot90(v)    the image of v under clockwise rotation by 90 degrees

    Ways to make vecs:
        vec(x, y)           Cartesian coordinates
        vec.cartes(x, y)    synonym of vec(x, y)
        vec.polar(r, ang)   polar coordinates
        vec.zero            the zero vector
        foot(pt, lin)       the foot of the perpendicular from pt to lin
        circ.pt(ang)        the point on circle circ at angle ang
        meet(lin1, lin2)    the intersection of lines lin1 and lin2
        intersect(a, b)     list of points where objects a and b meet

    vecs are value objects, meaning that their value is fixed at
    construction time.  You *can* change v.x and v.y, but you're
    on your own if you do.
    """
    def __init__(self, x, y):
        """The point with Cartesian coordinates x and y."""
        self.x = x
        self.y = y
    @classmethod
    def cartes(cls, x, y):
        """The point with Cartesian coordinates x and y."""
        return cls(x, y)
    @classmethod
    def polar(cls, r, ang):
        """The point with polar coordinates (r,ang)."""
        ang = ang*deg
        return cls(r*cos(ang), r*sin(ang))

    def __str__(self):
        return '(%s, %s)' % (self.x, self.y)
    def __repr__(self):
        return 'vec(%r, %r)' % (self.x, self.y)

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y
    def __hash__(self):
        return hash((self.x, self.y))

    def __add__(self, other):
        return vec(self.x + other.x, self.y + other.y)
    def __neg__(self):
        return vec(-self.x, -self.y)
    def __sub__(self, other):
        return vec(self.x - other.x, self.y - other.y)
    def __mul__(self, other):
        return vec(self.x*other, self.y*other);
    __rmul__ = __mul__
    def __truediv__(self, other):
        return vec(self.x/other, self.y/other);

    @property
    def normsq(self):
        """The square of the norm of this vec."""
        return dot(self, self)
    @property
    def norm(self):
        """The norm of this vec."""
        return sqrt(self.normsq)
    @property
    def angle(self):
        """The angle this vec forms with the positive x semi-axis."""
        return atan2(self.y, self.x)/deg

    def _affine(self, aff):
        return vec(aff.a*self.x + aff.b*self.y + aff.c,
                   aff.d*self.x + aff.e*self.y + aff.f)

vec.zero = vec(0,0)

def dot(u, v):
    """The dot product of vecs u and v."""
    return u.x*v.x + u.y*v.y

def unit(v):
    """The unit vector in the direction of the vec v.

    Raises GeometryError if v == vec.zero.
    """
    try:
        return v/v.norm
    except ZeroDivisionError:
        raise GeometryError('zero vector has no direction')

def rot90(v):
    """The image of vec v under clockwise rotation by 90 degrees."""
    return vec(-v.y, v.x)

def distfrom(pt):
    """A function returning the distance to pt."""
    return lambda x: (pt - x).norm

def bbox(*pts):
    """The bounding box of the given points."""
    lft = min([pt.x for pt in pts])
    bot = min([pt.y for pt in pts])
    rt  = max([pt.x for pt in pts])
    top = max([pt.y for pt in pts])
    return lft, bot, rt, top

def inflatebox(factor, lft, bot, rt, top):
    """A box concentric with (lft,bot,rt,top) but larger by factor."""
    midx = (rt + lft)/2
    halfwidth = factor*(rt - lft)/2
    midy = (top + bot)/2
    halfheight = factor*(top - bot)/2
    return midx - halfwidth, midy - halfheight, \
           midx + halfwidth, midy + halfheight

class line(object):
    """A line.

    Properties of and operations on lines:
        lin.a, lin.b, lin.c     coefficients of the normal equation ax + by = c
        lin.normal              a vec normal to lin
        lin.dir                 a vec parallel to lin
        meet(lin1, lin2)        the intersection of lines lin1 and lin2
        lin.side(pt)            which side of lin the point pt is on
        lin.sameside(pt1, pt2)  whether pt1 and pt2 are on the same side of lin
        lin.pt(t)               the point on line with parameter t

    lines are value objects, meaning that their value is fixed at
    construction time.  You *can* change lin.a, lin.b, and lin.c,
    but you're on your own if you do.

    Ways to make lines:
        line(a, b, c)           the line ax + by = c
        line.ptnorm(pt, v)      the line through pt with normal v
        line.ptdir(pt, v)       the line through pt in direction v
        line.ptslope(pt, m)     the line through pt with slope m
        line.ptangle(pt, ang)   the line through pt at angle ang
        join(pt1, pt2)          the line through pt1 and pt2
        parallel(pt, lin)       the line through pt and parallel to lin
        perp(pt, lin)           the line through pt and perpendicular to lin
        perpbisect(pt1, pt2)    the perpendicular bisector of segment pt1 pt2
        angbisect(p, q, r)      the angle bisector of angle p q r
    """
    def __init__(self, a, b, c):
        """The line with equation a*x + b*y = c.

        Raises GeometryError if a == 0 and b == 0.
        """
        self.a = a
        self.b = b
        self.c = c
        self._tuple = a, b, c
        if a == 0 and b == 0:
            raise GeometryError('%s has no direction' % self)
    @classmethod
    def ptnorm(cls, pt, v):
        """The line through pt with normal vector v.

        Raises GeometryError if v == vec.zero.
        """
        return line(v.x, v.y, dot(v, pt))
    @classmethod
    def ptdir(cls, pt, v):
        """The line through pt with direction vector v.

        Raises GeometryError if v == vec.zero.
        """
        return cls.ptnorm(pt, rot90(v))
    @classmethod
    def ptslope(cls, pt, m):
        """The line through pt with slope m."""
        return cls.ptdir(pt, vec(1, m))
    @classmethod
    def ptangle(cls, pt, ang):
        """The line through pt at angle ang."""
        return cls.ptdir(pt, vec(cos(ang*deg), sin(ang*deg)))

    def __str__(self):
        return '<line %sx + %sy = %s>' % self._tuple
    def __repr__(self):
        return 'line(%r, %r, %r)' % self._tuple

    @property
    def normal(self):
        """A normal vector for this line."""
        return vec(self.a, self.b)
    @property
    def dir(self):
        """A direction vector for this line."""
        return rot90(self.normal)

    def side(self, pt):
        """The side of this line that pt is on.

        Returns 0 if pt is on this line, 1 if pt is on one side, and -1
        if pt is on the other.  (Which side is which is not specified.)

        Note that, due to the inexactness of floating point arithmetic,
        this method does not give a really reliable way of testing
        whether a point lies on a line.
        """
        t = self.a*pt.x + self.b*pt.y - self.c
        if t < 0:
            return -1
        elif t > 0:
            return 1
        else:
            return 0

    def sameside(self, pt1, pt2):
        """Whether pt1 and pt2 are on the same side of this line."""
        return self.side(pt1)*self.side(pt2) > 0

    def pt(self, t):
        """The point at position t on this line.

        The parameterization used here is not specified, except that
        it is an affine map from R to R^2.  For example,
            lin.pt(.5) == (lin.pt(0) + lin.pt(1))/2
        (up to the exactness of the underlying arithmetic).
        """
        return foot(vec.zero, self) + t*self.dir

    def _affine(self, aff):
        # fixme test this
        aff = aff.invert()
        return line(aff.a*self.a + aff.d*self.b,
                    aff.b*self.a + aff.e*self.b,
                    self.c - aff.c*self.a - aff.f*self.b)

def join(pt1, pt2):
    """The line through (the terminal points of) vecs pt1 and pt2.

    Raises GeometryError if pt1 == pt2.
    """
    return line.ptdir(pt1, pt2-pt1)

def meet(lin1, lin2):
    """The intersection of lines lin1 and lin2.

    Raises GeometryError if lin1 and lin2 are parallel.
    """
    d = det(lin1.a, lin1.b, lin2.a, lin2.b)
    try:
        # Cramer's rule
        return vec(det(lin1.c, lin1.b, lin2.c, lin2.b)/d,
                   det(lin1.a, lin1.c, lin2.a, lin2.c)/d)
    except ZeroDivisionError:
        raise GeometryError('%s and %s are parallel' % (lin1, lin2))

def parallel(pt, lin):
    """The line through pt and parallel to lin."""
    return line.ptnorm(pt, lin.normal)

def perp(pt, lin):
    """The line through pt and perpendicular to lin."""
    return line.ptdir(pt, lin.normal)

def foot(pt, lin):
    """The foot of the perpendicular from pt to lin."""
    return meet(lin, perp(pt, lin))

def perpbisect(pt1, pt2):
    """The perpendicular bisector of the segment from pt1 to pt2.

    Raises GeometryError if pt1 == pt2.
    """
    return line.ptnorm((pt1+pt2)/2, pt1-pt2)

def angbisect(p, q, r):
    """The (internal) angle bisector of angle p q r.

    Raises GeometryError if the three points are not distinct.
    """
    return perpbisect(q + unit(p-q), q + unit(r-q))

class circle(object):
    """A circle.

    Properties of and operations on circles:
        circ.o          the centre of circ
        circ.r          the radius of circ
        circ.pt(ang)    the point on circ at angle ang

    Ways to make circles:
        circle(o, r)            the circle with centre o and radius r
        circle.ondiam(a, b)     the circle on segment a b as a diameter
        circle.through(a, b, c) the circle through points a, b, c
    """
    def __init__(self, o, r):
        """The circle with centre o and radius r.

        Raises GeometryError if r < 0.
        """
        self.o = o
        self.r = r
        if r < 0:
            raise GeometryError('radius %s is negative' % r)
    @classmethod
    def ondiam(cls, pt1, pt2):
        """The circle with pt1 and pt2 as endpoints of a diameter."""
        return cls((pt1+pt2)/2, (pt1-pt2).norm/2)
    @classmethod
    def through(cls, a, b, c):
        """The circle through points a, b, c.

        Raises GeometryError if a, b, c are collinear.
        """
        o = meet(perpbisect(a, b), perpbisect(b, c))
        return cls(o, (o-a).norm)

    def pt(self, ang):
        """The point on this circle at angle ang."""
        return self.o + vec.polar(self.r, ang)

# fixme incircle

class inversion(object):
    """An inversion of the plane.

    Properties of and operations on inversions:
        inv.circle              the circle of inversion
        inv(pt)                 the inverse image of a point
        inv.linecircle(circ)    the inverse image of a line
                                    not through the centre of inversion
        inv.circleline(circ)    the inverse image of a circle
                                    through the centre of inversion
        inv.circlecircle(circ)  the inverse image of a circle
                                    not through the centre of inversion

    Ways to make inversions:
        inversion(circ)     the inversion in the given circle
    """
    def __init__(self, circ):
        """The inversion in the given circle."""
        self.circle = circ

    def __eq__(self, other):
        return self.circle == other.circle
    def __hash__(self):
        return hash(self.circle)

    def __call__(self, pt):
        """The image of pt under this inversion.

        Raises ZeroDivisionError if pt == self.circle.o.
        """
        v = pt - self.circle.o
        return self.circle.o + unit(v)*(self.circle.r**2 / v.norm)

    def linecircle(self, lin):
        """The inverse image of a line not through the centre of inversion.

        If lin actually passes through the centre of inversion,
        you're on your own.
        """
        return circle.ondiam(self.circle.o, self(foot(self.circle.o, lin)))

    def circleline(self, circ):
        """The inverse image of a circle through the centre of inversion.

        If circ doesn't actually pass through the centre of inversion,
        you're on your own.
        """
        v = circ.o - self.circle.o
        return line.ptnorm(self(circ.o + v), v)

    def circlecircle(self, circ):
        """The inverse image of a circle not through the centre of inversion.

        If circ actually passes through the centre of inversion,
        you're on your own.
        """
        v = unit(circ.o - self.circle.o)*circ.r
        return circle.ondiam(self(circ.o + v), self(circ.o - v))

def _intersect_lineline(lin1, lin2):
    try:
        return [meet(lin1, lin2)]
    except GeometryError:
        return []

def _intersect_circleline(circ, lin):
    xpart = foot(circ.o, lin) - circ.o
    try:
        ylen = sqrt(circ.r**2 - xpart.norm**2)
    except ValueError:  # sqrt of negative value
        return []
    if ylen == 0:
        return [circ.o + xpart]
    else:
        ypart = ylen*unit(lin.dir)
        return [circ.o + xpart + ypart, circ.o + xpart - ypart]

def _intersect_circlecircle(circ1, circ2):
    odiff = circ2.o - circ1.o
    d = odiff.norm
    try:
        xlen = (d + (circ1.r**2 - circ2.r**2)/d)/2
    except ZeroDivisionError:  # d == 0; circles are concentric
        return []
    try:
        ylen = sqrt(circ1.r**2 - xlen**2)
    except ValueError:  # sqrt of negative; circles do not meet
        return []
    xdir = unit(odiff)  # odiff != vec.zero because d != 0
    xpart = xlen*xdir
    if ylen == 0:
        return [circ1.o + xpart]
    else:
        ypart = ylen*rot90(xdir)
        return [circ1.o + xpart + ypart, circ1.o + xpart - ypart]

_intersectfuncs = {
        (line, line): _intersect_lineline,
        (circle, line): _intersect_circleline,
        (line, circle): lambda lin, circ: _intersect_circleline(circ, lin),
        (circle, circle): _intersect_circlecircle,
    }

def intersect(obj1, obj2):
    """The list of points where obj1 and obj2 intersect.

    The order of the points is not defined.  If obj1 and obj2 do not
    intersect, an empty list is returned.

    Raises TypeError if I haven't implemented the intersection
    calculation for objects of the given types.  (As of now, obj1
    and obj2 must be lines and/or circles.)
    """
    t1 = type(obj1)
    t2 = type(obj2)
    try:
        func = _intersectfuncs[t1,t2]
    except KeyError:
        raise TypeError("can't compute intersection of %r (%r) and %r (%r)"
                        % (obj1, t1, obj2, t2))
    return func(obj1, obj2)

class affine(object):
    """An affine transformation of the plane.

    The map with coefficients a,b,c,d,e,f is the map that sends (x,y)
    to (a*x + b*y + c, d*x + e*y + f).

    Properties of and operations on affines:
        aff.[a-f]       the coefficients of aff
        aff1*aff2       the composition of aff1 and aff2
        aff(x)          the image of x under aff
                        (x must implement _affine, as e.g. vec and line do)

    Ways to make affines:
        affine(a,b,c,d,e,f)     directly from the coefficients
        affine.identity         the identity transformation
        translate(dx, dy)       the translation by (dx,dy)
        translate(v)            the translation by v
        scale(x, y)             scaling horizontally by x and vertically by y
        rotate(ang, pt)         the rotation about pt by ang degrees
        reflect(lin)            the reflection in lin
        aff.invert()            the inverse of aff
    """
    # fixme: construct affine given a triangle and its image
    # fixme: shears
    def __init__(self, a, b, c, d, e, f):
        """The affine map with the specified coefficients."""
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.e = e
        self.f = f
        self._tuple = a, b, c, d, e, f

    def __str__(self):
        return '<aff %s %s %s ; %s %s %s>' % self._tuple
    def __repr__(self):
        return 'affine(%r, %r, %r, %r, %r, %r)' % self._tuple

    def __eq__(self, other):
        return self._tuple == other._tuple
    def __hash__(self):
        return hash(self._tuple)

    def __call__(self, obj):
        """Compute the image of obj under this affine map.

        This method delegates to obj._affine().
        """
        return obj._affine(self)

    def __mul__(self, other):
        """The composition of this affine map with another."""
        return affine(self.a*other.a + self.b*other.d,
                      self.a*other.b + self.b*other.e,
                      self.a*other.c + self.b*other.f + self.c,
                      self.d*other.a + self.e*other.d,
                      self.d*other.b + self.e*other.e,
                      self.d*other.c + self.e*other.f + self.f)

    def invert(self):
        """The inverse of this affine map.

        Raises ZeroDivisionError if this map is uninvertible.
        """
        d = det(self.a, self.b, self.d, self.e)
        return affine(self.e/d, -self.b/d,
                        det(self.b, self.c, self.e, self.f)/d,
                      -self.d/d, self.a/d,
                        -det(self.a, self.c, self.d, self.f)/d)

affine.identity = affine(1, 0, 0, 0, 1, 0)

def translate(*args):
    """A translation.

    Called as translate(v), the translation by that vector.
    Called as translate(dx, dy), the translation by those x and y amounts.
    """
    if len(args) == 2:
        dx, dy = args
    elif len(args) == 1:
        dx, dy = args[0].x, args[0].y
    else:
        raise TypeError('translate() takes 1 or 2 arguments (%d given)'
                        % len(args))
    return affine(1, 0, dx, 0, 1, dy)

def scale(x, y):
    """Scaling horizontally by x and vertically by y."""
    return affine(x, 0, 0, 0, y, 0)

def rotate(ang, pt=vec.zero):
    """The rotation about pt by ang degrees."""
    ang = ang*deg
    c = cos(ang)
    s = sin(ang)
    rot = affine(c, -s, 0, s, c, 0)
    t = translate(pt)
    return t*rot*t.invert()

def reflect(lin):
    """The reflection in lin."""
    t = translate(lin.pt(0))
    r = rotate(lin.dir.angle)
    return t*r*scale(1, -1)*r.invert()*t.invert()

class GeometryError(Exception):
    """Raised when a necessary geometrical condition is not satisfied."""
    pass

class figure(object):  # fixme inline docs, implementation
    """A figure."""
    _PSPROC = """\
.4 setlinewidth
/dotradius 1.2 def
/*drawingmatrix* matrix currentmatrix def
/drawing {  % proc drawing -
    1 dict begin
        /*savedmatrix* matrix currentmatrix def
        *drawingmatrix* setmatrix
        exec
        *savedmatrix* setmatrix
    end
} def
/drawcoords {  % x y drawcoords x' y'
    transform *drawingmatrix* itransform
} def
/stroke0 /stroke load def
/stroke { { stroke0 } drawing } def
/box {
    4 copy 4 2 roll 8 -1 roll 5 -1 roll exch 5 1 roll 8 1 roll
    moveto lineto lineto lineto closepath
} def
/overbox {
    { fullwidth neg fullheight neg fullwidth 2 mul fullheight 2 mul box }
    drawing
} def
/dot {
    moveto { currentpoint dotradius 0 rmoveto dotradius 0 360 arc } drawing
} def
/clipdot { 2 copy dot stroke dot overbox eoclip newpath } def
/arrowd {  % d x1 y1 x2 y2
    4 2 roll drawcoords 4 2 roll drawcoords
    {
        3 index 3 index moveto
        % d x1 y1 x2 y2
        4 copy 3 -1 roll sub 3 1 roll exch sub atan
        % d x1 y1 x2 y2 ang
        3 1 roll translate rotate
        % d x1 y1
        pop pop
        % d
        currentlinewidth sub 0 translate
        0 0 lineto -4 2 moveto 0 0 lineto -4 -2 lineto stroke
    } drawing
} def
/arrow {  % x1 y1 x2 y2
    0 5 1 roll arrowd
} def
/arrowtodot {  % x1 y1 x2 y2
    dotradius neg 5 1 roll arrowd
} def
"""

    def __init__(self, width, height, margin, lft, bot, rt, top,
                 distort=False):
        """A figure of the given dimensions.
        
        The drawing region is width x height PostScript points; in
        the figure coordinate system it extends from (lft,bot) to
        (rt,top).  The generated EPS file includes an extra margin
        points on every side.
        """
        self._pscode = []
        self._texcode = []
        self.affine = affine.identity

        fullwidth = width + 2*margin
        fullheight = height + 2*margin
        cwidth = rt-lft
        cheight = top-bot
        if distort:
            self.xscale = width/cwidth
            self.yscale = height/cheight
        elif abs(cwidth/cheight) < abs(width/height):
            self.xscale = self.yscale = width/cwidth
        else:
            self.xscale = self.yscale = height/cheight
        bl = vec(lft, bot)
        tr = vec(rt, top)
        self.lftedge = line.ptdir(bl, vec(0, 1))
        self.rtedge = line.ptdir(tr, vec(0, 1))
        self.botedge = line.ptdir(bl, vec(1, 0))
        self.topedge = line.ptdir(tr, vec(1, 0))

        self.ps('%!PS-Adobe-3.0 EPSF-3.0')
        self.ps('%%%%BoundingBox: 0 0 %d %d' % (fullwidth, fullheight))
        self.ps(self._PSPROC)
        self.ps('/fullwidth %f def' % fullwidth)
        self.ps('/fullheight %f def' % fullheight)
        self.ps('/width %f def' % width)
        self.ps('/height %f def' % height)

        self.ps('gsave')
        self.translate(margin, margin)
        self.scale(self.xscale, self.yscale)
        self.translate(-lft, -bot)

        self.tex('\\setlength{\\unitlength}{1bp}%')
        self.tex('\\begin{picture}(%d,%d)' % (fullwidth, fullheight))

    def emit(self, ps, tex, fakeext=None):
        """Write the figure to the given files."""
        for code in self._pscode:
            print >>ps, code
        print >>ps, 'grestore'
        print >>ps, '%%EOF'

        if fakeext is None:
            includedname = ps.name
        else:
            includedname = splitext(ps.name)[0] + extsep + fakeext
        for code in self._texcode:
            print >>tex, code
        print >>tex, '\\put(0,0){\includegraphics*{%s}}' % includedname
        print >>tex, '\\end{picture}'

    def ps(self, code):
        """Add PostScript code."""
        self._pscode.append(code)

    def tex(self, code):
        """Add TeX code."""
        self._texcode.append(code)

    def clipdot(self, pt):
        """Add a clipped dot at the given point.

        A "clipped dot" is a little circle whose interior is clipped
        out, so that future drawing will not leave marks in it.
        """
        self.ps('%f %f clipdot' % (pt.x, pt.y))

    def polyline(self, *pts, **kwargs):
        """Draw a sequence of line segments.

        The path is closed if close=True is specified.
        """
        if len(pts) < 2:
            return
        self.ps('%f %f moveto' % (pts[0].x, pts[0].y))
        for pt in pts[1:]:
            self.ps('%f %f lineto' % (pt.x, pt.y))
        if kwargs.get('close', False):
            self.ps('closepath')
        self.ps('stroke')

    def polygon(self, *pts):
        """Draw a closed sequence of line segments."""
        return self.polyline(*pts, **{'close': True})

    def arrow(self, u, v):
        """Draw an arrow from u to v."""
        self.ps('%f %f %f %f arrow' % (u.x, u.y, v.x, v.y))

    def arrowtodot(self, u, v):
        """Draw an arrow from u to v which looks okay if there's a dot at v."""
        self.ps('%f %f %f %f arrowtodot' % (u.x, u.y, v.x, v.y))

    def linelimits(self, lin):
        """The points where lin intersects the edges of the drawing region.

        Returns an empty list if lin does not intersect the drawing
        region.  The order of the points in the list is not specified.
        """
        cat = 0
        try:
            lftpt = meet(lin, self.lftedge)
            rtpt = meet(lin, self.rtedge)
        except GeometryError:
            cat = 1
        try:
            toppt = meet(lin, self.topedge)
            botpt = meet(lin, self.botedge)
        except GeometryError:
            cat = 2
        if cat == 0:
            pts = [(lftpt,0), (rtpt,0), (toppt,1), (botpt,1)]
            pts.sort(key=lambda (a,b): dot(lin.dir, a))
            if pts[0][1] == pts[1][1]:
                return []
            else:
                return [pts[1][0], pts[2][0]]
        elif cat == 1:
            return [toppt, botpt]
        elif cat == 2:
            return [lftpt, rtpt]

    def line(self, lin):
        """Draw the given line."""
        self.polyline(*self.linelimits(lin))
        self.ps('stroke')

    def circle(self, circ, startang=0, endang=360):
        """Draw an arc of the given circle."""
        self.ps('%f %f %f %f %f arc stroke'
                % (circ.o.x, circ.o.y, circ.r, startang, endang))

    def curve(self, a, b, c, d):
        """Draw a cubic Bezier spline with the given control points."""
        self.ps('%f %f moveto %f %f %f %f %f %f curveto stroke'
                % (a.x, a.y, b.x, b.y, c.x, c.y, d.x, d.y))

    def qurve(self, a, b, c):
        """Draw a quadratic Bezier spline with the given control points."""
        self.curve(a, (a + 2*b)/3, (2*b + c)/3, c)

    def dashed(self):
        """Make future lines dashed."""
        self.ps('[1 2] .5 setdash')

    def undashed(self):
        """Make future lines solid."""
        self.ps('[] 0 setdash')

    _NORMAL_OFFSET = {
        'l': vec.polar(2, 0),
        'bl': vec.polar(2, 45),
        'b': vec.polar(2, 90),
        'br': vec.polar(2, 135),
        'r': vec.polar(2, 180),
        'tr': vec.polar(2, 225),
        't': vec.polar(2, 270),
        'tl': vec.polar(2, 315),
    }

    def label(self, text, align, pt, offset=vec.zero):
        """Place text at position pt.

        The align argument is the same as the argument to \\makebox
        in the LaTeX picture environment, or None for centering.

        Adjustments to the label's position can be given via the
        offset argument, which is a vec in drawing coordinates (i.e.,
        its x and y components are in PostScript points).
        """
        if align is None:
            baseoffset = vec.zero
        else:
            baseoffset = self._NORMAL_OFFSET[align]
        pos = self.affine(pt) + baseoffset + offset
        self.tex('\\put(%f,%f){\\makebox(0,0)[%s]{%s}}'
                 % (pos.x, pos.y, align, text))

    def translate(self, dx, dy):
        """Translate the figure coordinate system by (dx,dy)."""
        self.affine = self.affine * translate(dx,dy)
        self.ps('%f %f translate' % (dx,dy))

    def scale(self, x, y):
        """Scale the figure coordinate system by (x,y)."""
        self.affine = self.affine * scale(x,y)
        self.ps('%f %f scale' % (x,y))

    def rtangmark(self, a, b, c, size=4):
        """Add a right angle mark to the angle a-b-c."""
        # fixme does this work in distorted coords?
        leg1 = size*unit(b-a)/self.xscale
        leg2 = size*unit(c-b)/self.xscale
        self.ps('%f %f moveto %f %f 2 copy neg exch neg exch rmoveto '
            '%f %f rlineto rlineto stroke'
            % (b.x, b.y, leg1.x, leg1.y, leg2.x, leg2.y))

    def angmark(self, a, b, c, size=4):
        """Add a circular angle mark to the angle a-b-c."""
        # fixme does this work in distorted coords?
        self.circle(circle(b, size/self.xscale), (a-b).angle, (c-b).angle)

    def segmark(self, pt1, pt2, n, size=4, gap=2, angle=60, prop=1/2):
        """Add tick marks to the middle of the segment pt1-pt2."""
        # fixme move bulk of PS code into figure._PSPROC
        pt = (1-prop)*pt1 + prop*pt2
        v = vec.polar(size/2, angle)
        self.ps('%f %f moveto' % (pt.x, pt.y))
        self.ps('{')
        self.ps('currentpoint translate')
        self.ps('%f rotate' % (pt2 - pt1).angle)
        self.ps('%f currentlinewidth add %d mul 2 div neg 0 rmoveto' % (gap, n))
        self.ps('%f %f' % (v.x, v.y))
        self.ps('1 1 %d {' % n)
        self.ps('   pop')
        self.ps('   2 copy neg exch neg exch rmoveto')
        self.ps('   2 copy 2 mul exch 2 mul exch rlineto')
        self.ps('   2 copy neg exch neg exch rmoveto')
        self.ps('   %f currentlinewidth add 0 rmoveto' % gap)
        self.ps('} for')
        self.ps('pop pop stroke')
        self.ps('} drawing')

# class figure(object):
# 
# 
# /vadd {  % x1 y1 x2 y2 vadd x y
#     exch 4 -1 roll add 3 1 roll add
# } def
# /vmul {  % x y c vmul x' y'
#     dup 4 -1 roll mul 3 1 roll mul
# } def
# /vsub {  % x1 y1 x2 y2 vsub x y
#     -1 vmul vadd
# } def
# /lincomb {  % x1 y1 x2 y2 a b lincomb x' y'
#     exch 4 1 roll vmul 5 2 roll vmul vadd
# } def
# /vdot {  % x1 y1 x2 y2 vdot c
#     exch 4 -1 roll mul 3 1 roll mul add
# } def
# /det {  % a b c d det ad-bc
#     4 -1 roll mul 3 1 roll mul sub
# } def
# 
# /line {  % a b c line lin
#     [ 4 1 roll ]
# } def
# /linenormal {  % lin linenormal a b
#     aload pop pop
# } def
# /lineconst {  % lin lineconst c
#     aload pop 3 1 roll pop pop
# } def
# /pointnormal {  % x1 y1 a b pointnormal lin
#     4 2 roll 3 index 3 index vdot line
# } def
# /pointslope {  % x1 y1 m pointslope lin
#     -1 pointnormal
# } def
# /join {  % x1 y1 x2 y2 join lin
#     3 index 3 index vsub rot90 pointnormal
# } def
# /intersect {  % lin1 lin2 intersect x y
#     6 dict begin
#         /lin2 exch def
#         /lin1 exch def
#         /a [ lin1 linenormal pop lin2 linenormal pop ] cvx def
#         /b [ lin1 linenormal exch pop lin2 linenormal exch pop ] cvx def
#         /c [ lin1 lineconst lin2 lineconst ] cvx def
#         /detab a b det def
#         c b det detab div
#         a c det detab div
#     end
# } def
# /q2c {  % x1 y1 x2 y2 q2c x1' y1' x2' y2' x3' y3'
#     % p1 = (1/3) q0 + (2/3) q1
#     % p2 = (1/3) q2 + (2/3) q1
#     3 index 3 index currentpoint 2 3 div 1 3 div lincomb
#     6 2 roll 4 2 roll
#     3 index 3 index 2 3 div 1 3 div lincomb
#     4 2 roll
# } def
# /qurveto {  % x1 y1 x2 y2 qurveto -
#     q2c curveto
# } def
# /approxcurve {  % f df a b n approxcurve -
#     2 dict begin
#         /n exch def
#         /b exch def
#         /a exch def
#         /df exch def
#         /f exch def
#         /dx b a sub n div def
#         /tangent { dup f 1 index df pointslope } def
#         a dup f moveto
#         1 1 n {
#             dx mul a add /x exch def
#             /pt1 [ x dup f ] cvx def
#             /control [ x tangent x dx sub tangent intersect ] cvx def
#             control pt1 qurveto
#         } for
#     end
# } def
# 
# /dashed {  % proc dashed
#     gsave [1 2] .5 setdash exec grestore
# } def
# """

_USAGE = """\
Usage:
    python diagram.py base [ext]
Run Enkidu program diagram.py, generating base.eps and base.tex.
If ext is specified, refer to base.eps as "base.ext" in the base.tex
file.  (This "fake extension" feature is useful with PDFLaTeX; see
the documentation.)"""

def main(fig):
    import sys
    if len(sys.argv) not in [2,3]:
        print >>sys.stderr, _USAGE
        sys.exit(2)
    basename = sys.argv[1]
    try:
        fakeext = sys.argv[2]
    except IndexError:
        fakeext = None
    ps = open(basename + '.eps', 'w')
    try:
        tex = open(basename + '.tex', 'w')
        try:
            fig.emit(ps, tex, fakeext)
        finally:
            tex.close()
    finally:
        ps.close()
