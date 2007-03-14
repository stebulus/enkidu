from __future__ import division
from enkidu import *

#+
circ = circle(vec.zero, 1)
a = circ.pt(90)
b = circ.pt(-90)
p = circ.o - 2*vec(circ.r, 0)
q = withmax(distfrom(a), intersect(circ, join(p, a)))
p2 = meet(join(p, circ.o), join(q, b))

#-
lft, bot, rt, top \
    = inflatebox(1.2, *bbox(p, *[circ.pt(90*ang) for ang in range(4)]))
margin = 3
width = 158 - 2*margin
height = (top - bot)*width/(rt - lft)
f = figure(width, height, margin, lft, bot, rt, top)

#+
for pt in [a, b, p, q, p2]:
    f.clipdot(pt)
f.circle(circ)
f.polyline(p, a, b, q)
f.line(join(p, circ.o))
f.label('$A$', 'b', a)
f.label('$B$', 't', b)
f.label('$P$', 'b', p)
f.label('$Q$', 'br', q)
f.label("$P'$", 'bl', p2)
#-

if __name__ == '__main__':
    main(f)
