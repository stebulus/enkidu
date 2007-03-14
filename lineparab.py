from __future__ import division
from enkidu import *

#+
focus = vec(0, 1)
n = vec(1, 0)
q = vec(1/2, 2/9)
d = line.ptslope(n, 0)
ell = join(q, n)

qfoot = foot(q, d)
qcirc = circle(q, (qfoot - q).norm)
r = withmin(distfrom(n), intersect(qcirc, join(n, focus)))
p = meet(ell, parallel(focus, join(q, r)))
pfoot = foot(p, d)

#-
lft, bot, rt, top = inflatebox(1.2, *bbox(focus, n, p, q, r, qfoot, pfoot))
margin = 3
width = 134 - 2*margin
height = (top - bot)*width/(rt - lft)
f = figure(width, height, margin, lft, bot, rt, top)

#+
for pt in [focus, n, p, pfoot, q, qfoot, r]:
    f.clipdot(pt)
f.line(d)
f.line(ell)
f.polyline(n, focus, p, pfoot)
f.polyline(r, q, qfoot)
f.circle(qcirc)
f.rtangmark(q, qfoot, n)
f.rtangmark(p, pfoot, n)
f.label('$P$', 'b', p)
f.label("\\rlap{$P'$}\\phantom{$P$}", 't', pfoot)
f.label('$Q$', 'b', q)
f.label("\\rlap{$Q'$}\\phantom{$Q$}", 't', qfoot)
f.label('$N$', 't', n)
f.label('$F$', 'b', focus)
f.label('$R$', 'l', r, offset=vec(0,3))
f.label('$d$', 'r', leftmost(f.linelimits(d)))
f.label('$\\ell$', 'br', leftmost(f.linelimits(ell)))
#-

if __name__ == '__main__':
    main(f)
