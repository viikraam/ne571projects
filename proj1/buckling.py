from cyl2d import secant
from critNoRefl import LWRMaterial as lm

def obj(bsquare):
    return (lm.nSigF[0] / (lm.sigA[0] + lm.diffC[0]*bsquare) +
            lm.gtrans[1][0] / (lm.sigA[0] + lm.diffC[0]*bsquare) *
            lm.nSigF[1] /
            (lm.sigA[1] + lm.diffC[1]*bsquare) ) -1.0

# secant method loop
res = 100.
x0 = .0001
x1 =.0000000001
while abs(res) > 1e-9:
    y0 = obj(x0)
    y1 = obj(x1)
    print(y0, y1)
    try:
        newx = secant(x0, x1, y0, y1)
    except:
        break
    res = obj(newx)
    x0 = x1
    x1 = newx

print("b^2 should be {}".format(newx))

