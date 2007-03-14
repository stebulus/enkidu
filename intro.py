#+
from __future__ import division
from enkidu import *

circ = circle(vec.zero, 1)
a = circ.pt(180)
b = circ.pt(0)
p = circ.pt(110)
n = foot(p, join(a, b))

lft, bot, rt, top = bbox(a, b, circ.pt(90))
margin = 3
width = 158 - 2*margin
height = (top - bot)*width/(rt - lft)
f = figure(width, height, margin, lft, bot, rt, top)

for pt in [a,b,n,p]:
    f.clipdot(pt)
f.polyline(a, b)
f.polyline(p, n)
f.circle(circ, 0, 180)
f.circle(circle.ondiam(a, n), 0, 180)
f.circle(circle.ondiam(n, b), 0, 180)
f.circle(circle.ondiam(p, n), 0, 360)
f.label('$A$', 't', a)
f.label('$B$', 't', b)
f.label('$N$', 't', n)
f.label('$P$', 'b', p)

if __name__ == '__main__':
    main(f)
#-
