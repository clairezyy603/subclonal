#!/usr/bin/env python3
"""
Build a synthetic test case with known answers and check that `subclonal freq`
recovers them.

Two lineages differing at planted positions are mixed at known ratios, reads
are simulated from the mixture, and the recovered frequencies are compared
against what was put in. This is how the field-offset bug in the original
parser was caught: it silently returned zero differences on real data, and
only a dataset with a known answer exposed it.

Needs: samtools, bwa. No external read simulator — reads are generated here.
"""

import os
import random
import subprocess
import sys

random.seed(11)

OUT = "testdata"
GLEN = 60000
READ = 150
DEPTH = 60

# position -> frequency of the alt allele in each sample
PLANTED = {
    1500:  {"S1": 0.00, "S2": 1.00, "S3": 0.50},
    7200:  {"S1": 1.00, "S2": 1.00, "S3": 0.00},
    15050: {"S1": 0.30, "S2": 0.00, "S3": 0.00},
    15070: {"S1": 0.30, "S2": 0.00, "S3": 0.00},   # adjacent, same carriers
    15090: {"S1": 0.30, "S2": 0.00, "S3": 0.00},   # -> one event, not three
    28400: {"S1": 0.00, "S2": 0.60, "S3": 0.60},
    41000: {"S1": 1.00, "S2": 0.00, "S3": 1.00},
}
SAMPLES = ["S1", "S2", "S3"]


def run(cmd, **kw):
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL, **kw)


def main():
    if not all(os.system("which %s > /dev/null 2>&1" % t) == 0
               for t in ("bwa", "samtools")):
        sys.exit("needs bwa and samtools on PATH")
    os.makedirs(OUT, exist_ok=True)

    ref = [random.choice("ACGT") for _ in range(GLEN)]
    alt = {}
    for p in PLANTED:
        others = [b for b in "ACGT" if b != ref[p]]
        alt[p] = random.choice(others)

    with open(OUT + "/ref.fna", "w") as fh:
        fh.write(">chr\n")
        s = "".join(ref)
        for i in range(0, len(s), 70):
            fh.write(s[i:i+70] + "\n")
    run(["bwa", "index", OUT + "/ref.fna"])

    manifest = open(OUT + "/samples.tsv", "w")
    manifest.write("sample\tbam\n")

    for s in SAMPLES:
        reads = []
        n = GLEN * DEPTH // READ
        for k in range(n):
            start = random.randint(0, GLEN - READ - 1)
            seq = ref[start:start + READ].copy()
            for p, freqs in PLANTED.items():
                if start <= p < start + READ and random.random() < freqs[s]:
                    seq[p - start] = alt[p]
            reads.append("@r%d\n%s\n+\n%s\n" % (k, "".join(seq), "I" * READ))
        fq = "%s/%s.fq" % (OUT, s)
        open(fq, "w").write("".join(reads))
        bam = "%s/%s.bam" % (OUT, s)
        p1 = subprocess.Popen(["bwa", "mem", "-t", "2", OUT + "/ref.fna", fq],
                              stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p2 = subprocess.Popen(["samtools", "sort", "-o", bam],
                              stdin=p1.stdout, stdout=subprocess.DEVNULL,
                              stderr=subprocess.DEVNULL)
        p1.stdout.close(); p2.communicate()
        run(["samtools", "index", bam])
        manifest.write("%s\t%s\n" % (s, bam))
        sys.stderr.write("  built %s\n" % s)
    manifest.close()

    run([sys.executable, "subclonal.py", "freq",
         "--reference", OUT + "/ref.fna",
         "--manifest", OUT + "/samples.tsv",
         "--min-depth", "10", "-o", OUT + "/freq.tsv"])

    import csv
    rows = {int(r["pos"]): r for r in
            csv.DictReader(open(OUT + "/freq.tsv"), delimiter="\t")}

    print("\n%-8s %-6s %s" % ("pos", "", "  ".join("%-16s" % s for s in SAMPLES)))
    print("-" * 70)
    ok = True
    for p in sorted(PLANTED):
        # a site where every sample has the same frequency is not variable
        # between samples and is correctly dropped
        informative = len(set(PLANTED[p].values())) > 1
        r = rows.get(p + 1)   # pileup is 1-based
        if r is None:
            status = "dropped" if not informative else "** MISSED **"
            if informative:
                ok = False
            print("%-8d %-6s %s" % (p, "", status))
            continue
        cells = []
        for s in SAMPLES:
            got = float(r["%s_freq" % s])
            want = PLANTED[p][s]
            good = abs(got - want) <= 0.08
            if not good:
                ok = False
            cells.append("%.2f want %.2f %s" % (got, want, "ok" if good else "X"))
        print("%-8d %-6s %s" % (p, "", "  ".join("%-16s" % c for c in cells)))

    run([sys.executable, "subclonal.py", "events",
         "--freq", OUT + "/freq.tsv", "--max-gap", "200",
         "-o", OUT + "/events.tsv"])
    ev = list(csv.DictReader(open(OUT + "/events.tsv"), delimiter="\t"))
    merged = [e for e in ev if int(e["n_sites"]) == 3]
    print("\nevents: %d sites -> %d events" % (len(rows), len(ev)))
    if merged:
        e = merged[0]
        print("  the three adjacent sites at 15050/15070/15090 merged into one "
              "event spanning %s bp" % e["span_bp"])
    else:
        print("  ** the adjacent sites did not merge **")
        ok = False

    print("\n%s" % ("PASS" if ok else "FAIL"))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
