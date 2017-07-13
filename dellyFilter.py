import vcf, sys, glob
from operator import attrgetter
import operator
from collections import Counter

import re

chr_re = re.compile('^chr[0-9XY]+$')

def find_k_repeats(seq):
    counts = []
    for k in xrange(2, 5):
        for j in xrange(0, k):
            lb = seq[j:k + j]
            #print "First: ",lb
            lc = 1
            longest = (0, '')
            for i in xrange(k + j, len(seq), k):
                b = seq[i:i + k]
                if len(b) != k:
                    break
                #print "  lb = ",lb, "b = ",b, "i = ",i, "j = ",j, "lc = ",lc
                if b == lb:
                    lc += 1
                else:
                    if lc > longest[0]:
                        longest = (lc, b, k)
                    lb = b
                    lc = 1
            if lc > longest[0]:
                longest = (lc, b, k)
            counts.append(longest)
    return max(counts, key=lambda k:k[0])

def check_repeats(seq):
    lb = seq[0]
    bc = 1
    longest = 1
    for i in xrange(1, len(seq)):
        #print "  Index = {} last_base = {} base = {} longest = {} bc = {}".format(i, lb, seq[i], longest, bc)
        if lb == seq[i]:
            #print "  equal"
            bc += 1
        else:
            #print "  different"
            if bc > longest:
                longest = bc
            bc = 1
            lb = seq[i]
    if bc > longest:
        longest = bc
    longest_run = longest
    if longest_run >= 10:
        return False
    #print seq
    #print "  Longest: ",longest
    #print "  Base:    ",lb
    #print "  Tandem:  ",rep
    rep = find_k_repeats(seq)

    if rep[0] * rep[2] >= 15:
        return False

    return True

#ftree = {}

#with open(sys.argv[1]) as fp:
#   for l in fp:
#        t = l.rstrip().split("\t")
#        ref, lft, rgt, feat = t[1], int(t[2]), int(t[3]), t[7]
#        ftree.setdefault(ref, []).append(Interval(lft, rgt, feat))
#for k, v in ftree.iteritems():
#    ftree[k] = IntervalTree(v)

def print_header(fin):
    tabs = ['chrom','pos','id','ref','alt','qual','filter']
    for k in info_keys:
        tabs.append(k)
    for s in fin.samples:
        for k in fmt_keys:
            tabs.append('{}_{}'.format(s, k))
    return "\t".join(tabs)

def print_tabs(r):
    tabs = []
    tabs.extend([r.CHROM, r.POS, r.ID, r.REF, r.ALT, r.QUAL, r.FILTER])
    for k in info_keys:
        tabs.append(r.INFO.get(k, ''))

    for s in r.samples:
        for k in fmt_keys:
            tabs.append(s[k])
    return "\t".join(map(str, tabs))

for fname in glob.glob('*.vcf'):
    if 'filtered' in fname:
        continue
    fin = vcf.Reader(open(fname))
    foname = fname.replace('.vcf', '.filtered.vcf')
    foname2 = fname.replace('.vcf', '.filtered.txt')
    fout = vcf.Writer(open(foname, 'w'), fin)
    fout2 = open(foname2, 'w')
    infos = fin.infos
    info_keys = fin.infos.keys()
    fmt_keys  = fin.formats.keys()
    formats = fin.formats
    count, total, reps = 0, 0, 0
    rep_fail, ref_fail, pri_fail, sz_fail = 0, 0, 0, 0
    pe_fail = 0
    filter_fail = 0
    non_standard = 0
    prat_fail = 0
    small_fail = 0
    mapq_fail = 0
    rv_fail = 0
    fout2.write(print_header(fin))
    fout2.write("\n")

    check1, check2, check3, check4 = 0, 0, 0, 0

    for r in fin:
        #print "Chrom = {} Pos = {:,} Id = {} Ref = {} Alt = {} Qual = {} Filter = {}".format(r.CHROM, r.POS, r.ID, r.REF, r.ALT, r.QUAL, r.FILTER)
        #for k, v in r.INFO.iteritems():
        #    print "  {:<10} = {} [{}]".format(k, v, infos[k].desc)
        ref = None
        pri = None
        info = r.INFO
        pe = info['PE']
        sr = info['SR'] if 'SR' in info else 0
        pri, ref = r.samples[0], r.samples[1]
        #for s in r.samples:
        #    primary = s.sample[-1] == 'P'
        #    if primary:
        #        pri = s
        #    else:
        #        ref = s

        #Check SVLEN
        svlen = info['SVLEN'] if 'SVLEN' in info else -1
        chr1, chr2 = r.CHROM, info['CHR2']
        lft = r.POS
        rgt = info['END']
        total += 1
        #Skip small SVs

        if not chr_re.match(chr1) or not chr_re.match(chr2):
            #print "Skip ",chr1,chr2
            non_standard += 1
            continue

        if 'DEL' in str(r.ALT) and 0 < svlen < 2000:
            sz_fail += 1
            #print "  Skip"
            continue

        if info['MAPQ'] < 20:
            mapq_fail += 1
            continue

        #Check the reference
        rDV, rDR = ref['DV'], ref['DR']
        pDV, pDR = pri['DV'], pri['DR']
        #if pe < pDV:
        #    pe_fail += 1
        #    continue

        #if pe < 5:
        #    continue

        """
        if 0 < svlen < 500:
            if pri['RV'] < 1:
                small_fail += 1
                continue
        """

        prat = (1.0 * pDV / pe)
        srat = (1.0 * sr / pri['RV']) if pri['RV'] > 0 else 0
        if prat < 0.75 or prat > 1.25:
            prat_fail += 1
            continue

        if pri['RV'] <= 1 or ref['RV'] != 0:
            rv_fail += 1
            continue

        passed = False
        if pDV >= 25 and rDV <= 1:
            passed = True
            check1 += 1
        elif 10 <= pDV < 25 and rDV == 0:
            passed = True
            check2 += 1
        elif 4 <= pDV < 10 and rDV == 0:
            passed = True
            check3 += 1
        #elif (pDV < 4 and pe >= 10 and srat >= 0.5):
        #    passed = True
        #    check4 += 1

        if passed:
            if 'IMPRECISE' not in info and not check_repeats(info['CONSENSUS']):
                #print "  ** Tandem Fail"
                rep_fail += 1
                continue
            #has_centromere = False
            #if chr1 == chr2:
            #    if chr1 in ftree:
            #        type1 = set(x.value for x in ftree[chr1].overlaps(lft, lft))
            #        type2 = set(x.value for x in ftree[chr2].overlaps(rgt, rgt))
            #        #print chr1, lft, chr2, rgt, type1, type2
            #        if type1 and type2:
            #            has_centromere = True
            #else:
            #    if chr2 in ftree and chr1 in ftree:
            #        type1 = set(x.value for x in ftree[chr1].overlaps(lft, lft))
            #        type2 = set(x.value for x in ftree[chr2].overlaps(rgt, rgt))
            #        #print chr1, lft, chr2, rgt, type1, type2
            #        if type1 and type2:
            #            has_centromere = True
            #if has_centromere:
            #    #print "  ** Repeat Region Fail"
            #    reps += 1
            #    continue
            fout.write_record(r)
            fout2.write(print_tabs(r))
            fout2.write("\n")
            count += 1
        else:
            filter_fail += 1

        #print "--------------" * 5

    if total > 0:
        print "{}".format(fname)
        print "  {:<15} {:>6} [{:.2f}%]".format("Passed", count, 100.0 * count / total)
    if count > 0:
        print "    {:<15} {:>6} [{:.2f}%]".format("Check 1", check1, 100.0 * check1 / count)
        print "    {:<15} {:>6} [{:.2f}%]".format("Check 2", check2, 100.0 * check2 / count)
        print "    {:<15} {:>6} [{:.2f}%]".format("Check 3", check3, 100.0 * check3 / count)
    if total > 0:
        print "  {:<15} {:>6} [{:.2f}%]".format("Size Fail", sz_fail, 100.0 * sz_fail / total)
        print "  {:<15} {:>6} [{:.2f}%]".format("Non standard", non_standard, 100.0 * non_standard / total)
        print "  {:<15} {:>6} [{:.2f}%]".format("RV Fail", rv_fail, 100.0 * rv_fail / total)
        print "  {:<15} {:>6} [{:.2f}%]".format("Filter Fail", filter_fail, 100.0 * filter_fail / total)
        print "  {:<15} {:>6} [{:.2f}%]".format("PE Fail", pe_fail, 100.0 * pe_fail / total)
        print "  {:<15} {:>6} [{:.2f}%]".format("MapQ Fail", mapq_fail, 100.0 * mapq_fail / total)
        print "  {:<15} {:>6} [{:.2f}%]".format("DVp / PE", prat_fail, 100.0 * prat_fail / total)
        print "  {:<15} {:>6} [{:.2f}%]".format("Repeat Fail", reps, 100.0 * reps / total)
        print "  {:<15} {:>6} [{:.2f}%]".format("Tandem Fail", rep_fail, 100.0 * rep_fail / total)

