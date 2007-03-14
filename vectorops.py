from __future__ import division
from enkidu import *

#+
u = vec(2, 1)
v = vec(-1, 1)
vecs = [u, v, u+v, -u, (v-u)/2, 2*v]

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
f.polyline(-u, v, u+v, u)
f.undashed()
f.label('\\texttt{vec.zero}', 'tl', vec.zero)
f.label('\\texttt{u}', 'l', u)
f.label('\\texttt{v}', 'b', v)
f.label('\\texttt{2*v}', 'br', 2*v)
f.label('\\texttt{u+v}', 'bl', u+v)
f.label('\\texttt{-u}', 'tr', -u)
f.label('\\texttt{(v-u)/2}', 'br', (v-u)/2)
#-

if __name__ == '__main__':
    main(f)
