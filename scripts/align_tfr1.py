import re, sys, json
def seq(path):
    txt = open(path).read()
    block = txt.split("\nSQ   ",1)[1].split("\n",1)[1]
    s = "".join(block.split("//")[0].split())
    return s
h = seq(sys.argv[1]); m = seq(sys.argv[2])
print("human", len(h), "mouse", len(m))
# simple NW, match 1 / mismatch -1 / gap -2 (linear)
n, k = len(h), len(m)
GAP = -2
prev = [GAP*j for j in range(k+1)]
ptr = []
for i in range(1, n+1):
    cur = [GAP*i] + [0]*k
    row = bytearray(k+1)
    hi = h[i-1]
    for j in range(1, k+1):
        d = prev[j-1] + (1 if hi == m[j-1] else -1)
        u = prev[j] + GAP
        l = cur[j-1] + GAP
        if d >= u and d >= l: cur[j] = d; row[j] = 0
        elif u >= l: cur[j] = u; row[j] = 1
        else: cur[j] = l; row[j] = 2
    ptr.append(row); prev = cur
# traceback
i, j = n, k; pairs = []
while i > 0 and j > 0:
    t = ptr[i-1][j]
    if t == 0: pairs.append((i, j)); i -= 1; j -= 1
    elif t == 1: pairs.append((i, None)); i -= 1
    else: pairs.append((None, j)); j -= 1
while i > 0: pairs.append((i, None)); i -= 1
while j > 0: pairs.append((None, j)); j -= 1
pairs.reverse()
ident = sum(1 for a,b in pairs if a and b and h[a-1]==m[b-1])
aln = len(pairs)
print("aligned cols %d, identities %d, overall identity %.1f%%" % (aln, ident, 100*ident/min(n,k)))
# ectodomain 89..760 human
eco = [(a,b) for a,b in pairs if a and 89 <= a <= 760]
ei = sum(1 for a,b in eco if b and h[a-1]==m[b-1])
print("ectodomain 89-760: %d human residues, %d identical to mouse = %.1f%%" % (len(eco), ei, 100*ei/len(eco)))
asy = [(a,b) for a,b in pairs if a and 121 <= a <= 760]
ai = sum(1 for a,b in asy if b and h[a-1]==m[b-1])
print("assay construct 121-760: %d residues, %.1f%% identical" % (len(asy), 100*ai/len(asy)))
# divergence hotspots: 20-residue windows with lowest identity in 121-760
win = 20
rows = []
for s in range(121, 761-win):
    w = [(a,b) for a,b in asy if s <= a < s+win]
    if not w: continue
    idw = sum(1 for a,b in w if b and h[a-1]==m[b-1])
    rows.append((idw/len(w), s, s+win-1))
rows.sort()
print("\nLEAST conserved 20-aa windows (avoid for a mouse-cross-reactive epitope):")
for r in rows[:8]: print("  %d-%d  %.0f%% id" % (r[1], r[2], 100*r[0]))
print("\nMOST conserved 20-aa windows (candidate epitope cores):")
for r in rows[-8:]: print("  %d-%d  %.0f%% id" % (r[1], r[2], 100*r[0]))
# human histidines in the assay construct (native His = pre-existing pH-sensitive positions)
his = [a for a,b in asy if h[a-1]=='H']
print("\nhuman His in 121-760: %d -> %s" % (len(his), his))
json.dump({"pairs_len":aln}, open("/dev/null","w"))
