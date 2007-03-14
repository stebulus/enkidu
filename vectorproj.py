from __future__ import division
from enkidu import *

#+
u = vec(2, 1)
v = vec(0, 1.5)
projuv = v*dot(u, v)/v.normsq
projvu = u*dot(u, v)/u.normsq
vecs = [u, v, projuv, projvu]

#-
lft, bot, rt, top = bbox(vec.zero, *vecs)
margin = 12
width = 132 - 2*margin
height = (top - bot)*width/(rt - lft)
f = figure(width, height, margin, lft, bot, rt, top)

#+
for pt in vecs + [vec.zero]:
    f.clipdot(pt)
for pt in vecs:
    f.arrowtodot(vec.zero, pt)
f.dashed()
f.polyline(u, projuv)
f.polyline(v, projvu)
f.undashed()
f.label('$0$', 't', vec.zero)
f.label('$u$', 'l', u)
f.label('$v$', 'b', v)
f.label(r'$\frac{\langle u,v\rangle}{\|v\|^2}v$', 'r', projuv)
f.label(r'$\frac{\langle u,v\rangle}{\|u\|^2}u$', 'tl', projvu)
#-

if __name__ == '__main__':
    main(f)
