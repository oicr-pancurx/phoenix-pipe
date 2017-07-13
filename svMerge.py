#!/.mounts/labs/steinlab/public/gavin/local/bin/python
import sys, glob
from operator import attrgetter
import operator
from collections import Counter
import re
from operator import attrgetter
import math

class Interval(object):
    __slots__ = ['start', 'stop', 'value']
    def __init__(self, start, stop, value):
        self.start = start
        self.stop  = stop
        self.value = value

    def __lt__(self, rhs):
        return self.start < rhs.start

    def __repr__(self):
        return 'Interval({}, {}, {})'.format(self.start, self.stop, self.value)

class IntervalTree(object):
    __slots__ = ['intervals', 'left', 'right', 'center']

    def __init__(self, ivals, depth = 16, minbucket = 64, lftext = 0, rgtext=0, maxbucket=512):
        self.intervals = []
        depth -= 1
        self.center = -1
        self.left   = None
        self.right  = None

        if depth == 0 or (len(ivals) < minbucket and len(ivals) < maxbucket):
            self.intervals = tuple(ivals)
        else:
            if lftext == 0 and rgtext == 0:
                ivals.sort()

            lftp, rgtp, centp = 0, 0, 0

            if lftext or rgtext:
                lftp = lftext
                rgtp = rgtext
            else:
                lftp = ivals[0].start
                rgtp = max(ivals, key=attrgetter('stop'))

            centp = ivals[len(ivals) / 2].start
            self.center=  centp

            lfts = []
            rgts = []

            for i in ivals:
                if i.stop < centp:
                    lfts.append(i)
                elif i.start > centp:
                    rgts.append(i)
                else:
                    self.intervals.append(i)

            if lfts:
                self.left = IntervalTree(lfts, depth, minbucket, lftp, centp, maxbucket)
            
            if rgts:
                self.right = IntervalTree(rgts, depth, minbucket, centp, rgtp, maxbucket)

    def overlaps(self, start, stop, debug = False):
        if debug:
            if self.intervals:
                print "First interval: ",self.intervals[0]
        if self.intervals and self.intervals[0].start <= stop:
            for i in self.intervals:
                if debug:
                    print "  Checking: ",i
                if i.stop >= start and i.start <= stop:
                    yield i
        if debug:
            print "Left: ",start," < ",self.center
        if self.left and start <= self.center:
            for i in self.left.overlaps(start, stop, debug):
                yield i

        if debug:
            print "Right: ",start," < ",self.center
        if self.right and stop >= self.center:
            for i in self.right.overlaps(start, stop, debug):
                yield i

    def contained(self, start, stop):
        if self.intervals and self.intervals[0].start <= stop:
            for i in self.intervals:
                if i.start >= start and i.stop <= stop:
                    yield i

        if self.left and start <= self.center:
            for i in self.left.contained(start, stop):
                yield i


        if self.right and stop >= self.center:
            for i in self.right.contained(start, stop):
                yield i

mapping = {
    'ITX':'INV',
    'DEL':'DEL',
    'CTX':'TRA',
    'INS':'DUP',
    'INV':'INV'
}

class SV(object):
    def __init__(self, index, ref1, pos1, ref2, pos2, tid, tool, **extra):
        self.ref1 = ref1
        self.pos1 = pos1
        self.ref2 = ref2
        self.pos2 = pos2
        self.index = index
        self.extra = extra
        self.tool = tool
        self.tid  = tid

    def __getitem__(self, k):
        return self.extra[k]

    def __str__(self):
        return "    Position 1 = {} {:,} Position 2 = {} {:,} Type = {} Extra: {}".format(self.ref1, self.pos1, self.ref2, self.pos2, self.tid, ", ".join("{} = {}".format(k, v ) for k, v in self.extra.iteritems()))


def parse_crest(fname):
    print "Parsing crest"
    print "  ",fname
    with open(fname) as fp:
        svs = []
        #h = fp.readline().rstrip().split("\t")
        for l in fp:
            chr1, pos1, strand1, c1, chr2, pos2, strand2, c2, tid, sz, tumor, normal, tmp = l.rstrip().split("\t")
            pos1, pos2 = int(pos1), int(pos2)
            svs.append(SV(len(svs) + 1, chr1, pos1, chr2, pos2, mapping[tid], 'crest', soft_clip1=c1, soft_clip2=c2, size=sz, normal=normal, tumor=tumor, strand1 = strand1, strand2=strand2))
    return svs

def parse_delly(prefix):
    svs = []
    header = []
    print "Parsing Delly"
    for f in glob.glob(prefix + '*.filtered.txt'):
        print "  ",f
        with open(f) as fp:
            h = fp.readline().rstrip().split("\t")
            skip = set((0, 1, 9, 10))
            header = [h[i] for i in xrange(0, len(h)) if i not in skip]
            cols = [i for i in xrange(0, len(h)) if i not in skip]
            for l in fp:
                t = l.rstrip().split("\t")
                chr1, pos1, chr2, pos2 = t[0], int(t[1]), t[9], int(t[10])
                extra = {h[i]:t[i] for i in cols}
                svs.append(SV(len(svs), chr1, pos1, chr2, pos2, t[20], 'delly', **extra))
    return svs, header

def merge_internal_svs(svs, window, get_best):
    trees = {}
    found = False
    for i, s in enumerate(svs):
        chr1, pos1, chr2, pos2 = s.ref1, s.pos1, s.ref2, s.pos2
        if chr1 == chr2 and min(pos1, pos2) + window > max(pos1, pos2) - window:
            dist = int(math.floor((max(pos1, pos2) - min(pos1, pos2) + 1) / 2.0))
            trees.setdefault(chr1, []).append(Interval(min(pos1, pos2) - window, min(pos1, pos2) + dist, i))
            trees.setdefault(chr2, []).append(Interval(max(pos1, pos2) - dist, max(pos1, pos2) + window, i))
        else:
            trees.setdefault(chr1, []).append(Interval(pos1 - window, pos1 + window, i))
            trees.setdefault(chr2, []).append(Interval(pos2 - window, pos2 + window, i))
    trees = {k:IntervalTree(sorted(v)) for k, v in trees.iteritems()}
    resolved = []
    merged = set()
    for i, s in enumerate(svs):
        if i in merged:
            continue

        overlaps1 = set(iv.value for iv in trees[s.ref1].overlaps(s.pos1, s.pos1) if iv.value != i and svs[iv.value].tid == s.tid)
        overlaps2 = set(iv.value for iv in trees[s.ref2].overlaps(s.pos2, s.pos2) if iv.value != i and svs[iv.value].tid == s.tid)
        shared = overlaps1 & overlaps2
        if len(shared) > 0:
            #print "    ",i, overlaps1, overlaps2, shared
            #print "    ", s

            #print "  Overlaps1: ",overlaps1
            #print "  Overlaps2: ",overlaps2
            #print "  Shared:    ",shared
            #print "    ",s
            #for k, v in s.extra.iteritems():
            #    print "     ",k, v
            #for si in shared:
            #    print "    ",svs[si]
            #    for k, v in svs[si].extra.iteritems():
            #        print "     ",k, v
            #print "  ",[(si, int(svs[si]['soft_clip1']) + int(svs[si]['soft_clip2'])) for si in shared]
            #print ""
            best = get_best(shared, svs)
            #best = 
            #print "Best: ",best
            found = True
            resolved.append(svs[best[0]])
            for ii in shared:
                merged.add(ii)
        else:
            resolved.append(s)
            #return True
    print "  Count before = {} count after merging = {}".format(len(svs), len(resolved))
    for i, s in enumerate(resolved):
        s.index = i
    return resolved

def merge_svs(crest, delly, window):
    trees = {}
    for i, s in enumerate(delly):
        chr1, pos1, chr2, pos2 = s.ref1, s.pos1, s.ref2, s.pos2
        if chr1 == chr2 and min(pos1, pos2) + window > max(pos1, pos2) - window:
            print "Overlap Problem Didn't expect this!!!"
            sys.exit()
        #print "Delly: ",s
        #print "  {}: {:,} - {:,}".format(chr1, pos1-window, pos1+window)
        #print "  {}: {:,} - {:,}".format(chr2, pos2-window, pos2+window)
        trees.setdefault(chr1, []).append(Interval(pos1 - window, pos1 + window, i))
        trees.setdefault(chr2, []).append(Interval(pos2 - window, pos2 + window, i))

    trees = {k:IntervalTree(sorted(v)) for k, v in trees.iteritems()}
    overlapping = set()
    dmap, cmap, skips = {}, {}, set()
    for s in crest:
        if s.ref1 not in trees or s.ref2 not in trees:
            continue
        overlaps1 = set(iv.value for iv in trees[s.ref1].overlaps(s.pos1, s.pos1) if iv.value != i and delly[iv.value].tid == s.tid)
        overlaps2 = set(iv.value for iv in trees[s.ref2].overlaps(s.pos2, s.pos2) if iv.value != i and delly[iv.value].tid == s.tid)
        shared = overlaps1 & overlaps2
        #print "Overlap"
        #print "  ",s
        #for iv in trees[s.ref1].overlaps(s.pos1, s.pos1):
        #    print "    1: ",delly[iv.value]
        #for iv in trees[s.ref2].overlaps(s.pos2, s.pos2):
        #    print "    2: ",delly[iv.value]
        #for o in overlaps1:
        #    print "    1: ",delly[o]
        #for o in overlaps2:
        #    print "    2: ",delly[o]

        if len(shared) == 1:
            index, = shared
            dmap[index] = s.index
            cmap[s.index] = index
            #print "  ",delly[index]
        elif len(shared) > 1:
            skips.add(s.index)

            #print "More than one overlap ERROR!!"
            #sys.exit()


    print "Found {} overlapping SVS, Crest Total = {}, Delly Total = {}".format(len(dmap), len(crest), len(delly))

    return dmap, cmap, skips

crest_svs  = parse_crest(sys.argv[1])
delly_svs, delly_header = parse_delly(sys.argv[2])


best_crest = lambda shared, svs: max(((si, int(svs[si]['soft_clip1']) + int(svs[si]['soft_clip2'])) for si in shared), key=lambda k: k[1])
#best_delly = lambda si, svs: max(((si, int(svs[si]['soft_clip1']) + int(svs[si]['soft_clip2'])) for si in shared), key=lambda k: k[1])
best_delly = lambda shared, svs: sys.exit()

WINDOW = 25
print "Merging Delly for Internal Overlaps"
#delly_svs = merge_internal_svs(delly_svs, WINDOW, best_delly)

#sys.exit()

print "Merging Crest for Internal Overlaps"
crest_svs = merge_internal_svs(crest_svs, WINDOW, best_crest)

print "Merging Crest and Delly"
dmap, cmap, skips = merge_svs(crest_svs, delly_svs, WINDOW)

print "Writing Output"
output = []

header = ['chrom1', 'pos1', 'chrom2', 'pos2', 'type']
header += ['delly_{}'.format(x) for x in delly_header]
header += ['crest_{}'.format(x) for x in ('soft_clip1', 'strand1', 'soft_clip2', 'strand2', 'size', 'tumor_pairs','normal_pairs')]

fout = open(sys.argv[3], 'w')
fout.write("\t".join(header))
fout.write("\n")

#print header
output = []
for c in crest_svs:
    row = None
    if c.index in skips:
        continue
    if c.index in cmap:
        d = delly_svs[cmap[c.index]]
        row = [d.ref1, d.pos1, d.ref2, d.pos2, d.tid]
        row += [d[k] for k in delly_header]
        row += [c['soft_clip1'], c['strand1'], c['soft_clip2'], c['strand2'], c['size'], c['tumor'], c['normal']]
    else:
        row = [c.ref1, c.pos1, c.ref2, c.pos2, c.tid]
        row += ['NA' for k in delly_header]
        row += [c['soft_clip1'], c['strand1'], c['soft_clip2'], c['strand2'], c['size'], c['tumor'], c['normal']]
    output.append(row)

for d in delly_svs:
    if d.index in dmap:
        continue
    row = [d.ref1, d.pos1, d.ref2, d.pos2, d.tid]
    row += [d[k] for k in delly_header]
    row += ['NA','NA','NA','NA','NA', 'NA', 'NA']
    output.append(row)

output.sort()
for row in output:
    fout.write("\t".join(map(str, row)))
    fout.write("\n")
