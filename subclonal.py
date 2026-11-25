#!/usr/bin/env python3
"""
subclonal — read allele frequencies out of bacterial pileups without calling genotypes.

A genotype caller answers "which base is this sample?". In a bacterial
population sampled as colonies, that is the wrong question: a site at 60% A
and 40% G is not a noisy A, it is two lineages. Collapsing it to one call
discards the signal.

This reports the frequency at every variable site, per sample, and leaves the
interpretation to you.

Subcommands:
  freq     per-sample allele frequencies at variable sites
  events   merge adjacent sites that share a sample set into single events
"""

import argparse
import csv
import os
import shutil
import subprocess
import sys


# --------------------------------------------------------------------------
# freq
# --------------------------------------------------------------------------

def need(tool):
    if shutil.which(tool) is None:
        sys.exit("error: %s not found on PATH" % tool)


def parse_pileup_bases(field, ref):
    """Count A/C/G/T in one mpileup base string. Indels are skipped, not called."""
    counts = {"A": 0, "C": 0, "G": 0, "T": 0}
    i = 0
    n = len(field)
    while i < n:
        c = field[i]
        if c == "^":          # read start, next char is mapping quality
            i += 2
            continue
        if c == "$":          # read end
            i += 1
            continue
        if c in "+-":         # indel: +4ACGT — length then that many bases
            j = i + 1
            num = ""
            while j < n and field[j].isdigit():
                num += field[j]
                j += 1
            i = j + (int(num) if num else 0)
            continue
        if c in ".,":
            counts[ref] += 1
        elif c.upper() in counts:
            counts[c.upper()] += 1
        i += 1
    return counts


def pileup(reference, bams, names, regions, mapq, baseq, threads):
    cmd = ["samtools", "mpileup", "-f", reference, "-q", str(mapq),
           "-Q", str(baseq), "-a"]
    if regions:
        cmd += ["-r", regions]
    cmd += bams
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.DEVNULL, text=True, bufsize=1)
    for line in proc.stdout:
        p = line.rstrip("\n").split("\t")
        if len(p) < 4:
            continue
        contig, pos, ref = p[0], int(p[1]), p[2].upper()
        if ref not in "ACGT":
            continue
        per = []
        for k in range(len(bams)):
            base = 4 + k * 3
            if base + 1 >= len(p):
                per.append((0, {"A": 0, "C": 0, "G": 0, "T": 0}))
                continue
            depth = int(p[base - 1]) if p[base - 1].isdigit() else 0
            per.append((depth, parse_pileup_bases(p[base], ref)))
        yield contig, pos, ref, per
    proc.wait()


def cmd_freq(args):
    need("samtools")
    bams, names = [], []
    if args.manifest:
        with open(args.manifest) as fh:
            for rec in csv.DictReader(fh, delimiter="\t"):
                if not rec.get("sample") or not rec.get("bam"):
                    sys.exit("manifest needs columns: sample, bam")
                names.append(rec["sample"])
                bams.append(rec["bam"])
    else:
        for b in args.bam:
            bams.append(b)
            names.append(os.path.basename(b).rsplit(".", 1)[0])
    for b in bams:
        if not os.path.exists(b):
            sys.exit("error: no such file: %s" % b)
    sys.stderr.write("%d samples, reference %s\n" % (len(bams), args.reference))

    out = sys.stdout if args.out == "-" else open(args.out, "w", newline="")
    w = csv.writer(out, delimiter="\t")
    header = ["contig", "pos", "ref", "alt", "n_samples_variant"]
    for s in names:
        header += ["%s_freq" % s, "%s_depth" % s]
    w.writerow(header)

    kept = seen = 0
    for contig, pos, ref, per in pileup(args.reference, bams, names,
                                        args.region, args.mapq, args.baseq,
                                        args.threads):
        seen += 1
        # pooled non-reference counts decide the alt allele
        pooled = {"A": 0, "C": 0, "G": 0, "T": 0}
        for depth, c in per:
            for b in "ACGT":
                pooled[b] += c[b]
        alts = [(pooled[b], b) for b in "ACGT" if b != ref]
        alts.sort(reverse=True)
        if not alts or alts[0][0] == 0:
            continue
        alt = alts[0][1]

        row, nvar, usable = [], 0, 0
        for depth, c in per:
            tot = sum(c.values())
            if tot < args.min_depth:
                row += ["NA", tot]
                continue
            usable += 1
            f = c[alt] / tot
            row += ["%.3f" % f, tot]
            if f >= args.min_freq:
                nvar += 1
        if usable < args.min_samples:
            continue
        # a site where every sample agrees carries no information between samples
        if nvar == 0 or nvar == usable:
            if not args.keep_invariant:
                continue
        w.writerow([contig, pos, ref, alt, nvar] + row)
        kept += 1
        if kept % 500 == 0:
            sys.stderr.write("  %d variable sites\r" % kept)
            out.flush()

    sys.stderr.write("\n%d positions examined, %d variable sites written\n"
                     % (seen, kept))
    if out is not sys.stdout:
        out.close()
        sys.stderr.write("wrote %s\n" % args.out)


# --------------------------------------------------------------------------
# events
# --------------------------------------------------------------------------

def cmd_events(args):
    """Adjacent sites carried by the same samples are one event, not several."""
    rows = list(csv.DictReader(open(args.freq), delimiter="\t"))
    if not rows:
        sys.exit("error: %s is empty" % args.freq)
    samples = [c[:-5] for c in rows[0] if c.endswith("_freq")]
    sys.stderr.write("%d sites, %d samples\n" % (len(rows), len(samples)))

    def carriers(r):
        out = []
        for s in samples:
            v = r.get("%s_freq" % s, "NA")
            if v != "NA" and float(v) >= args.min_freq:
                out.append(s)
        return tuple(out)

    keyed = []
    for r in rows:
        keyed.append((r["contig"], int(r["pos"]), carriers(r), r))
    keyed.sort(key=lambda x: (x[0], x[1]))

    events, cur = [], None
    for contig, pos, key, r in keyed:
        if (cur and cur["contig"] == contig and cur["key"] == key
                and pos - cur["end"] <= args.max_gap):
            cur["end"] = pos
            cur["n_sites"] += 1
        else:
            if cur:
                events.append(cur)
            cur = {"contig": contig, "start": pos, "end": pos,
                   "key": key, "n_sites": 1}
    if cur:
        events.append(cur)

    out = sys.stdout if args.out == "-" else open(args.out, "w", newline="")
    w = csv.writer(out, delimiter="\t")
    w.writerow(["contig", "start", "end", "span_bp", "n_sites",
                "n_carriers", "carriers"])
    for e in events:
        w.writerow([e["contig"], e["start"], e["end"],
                    e["end"] - e["start"] + 1, e["n_sites"],
                    len(e["key"]), ",".join(e["key"])])
    sys.stderr.write("%d sites merged into %d events (max gap %d bp)\n"
                     % (len(rows), len(events), args.max_gap))
    merged = [e for e in events if e["n_sites"] > 1]
    if merged:
        sys.stderr.write("  %d events span more than one site; largest spans %d sites\n"
                         % (len(merged), max(e["n_sites"] for e in merged)))
    if out is not sys.stdout:
        out.close()
        sys.stderr.write("wrote %s\n" % args.out)


# --------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(
        description="Allele frequencies from bacterial pileups, without genotype calls.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
  subclonal.py freq --reference ref.fna --manifest samples.tsv -o freq.tsv
  subclonal.py events --freq freq.tsv --max-gap 200 -o events.tsv
""")
    sub = p.add_subparsers(dest="cmd", required=True)

    f = sub.add_parser("freq", help="per-sample allele frequencies")
    f.add_argument("--reference", required=True)
    f.add_argument("--manifest", help="TSV with columns: sample, bam")
    f.add_argument("--bam", nargs="*", default=[], help="BAMs (alternative to --manifest)")
    f.add_argument("--region", help="contig:start-end, restrict to one region")
    f.add_argument("-o", "--out", default="-")
    f.add_argument("--mapq", type=int, default=20,
                   help="minimum mapping quality (default 20) — raise it where "
                        "paralogues are a concern")
    f.add_argument("--baseq", type=int, default=20)
    f.add_argument("--min-depth", type=int, default=10,
                   help="below this a sample is NA at that site (default 10)")
    f.add_argument("--min-freq", type=float, default=0.10,
                   help="frequency at which a sample counts as carrying the "
                        "alt allele (default 0.10)")
    f.add_argument("--min-samples", type=int, default=2,
                   help="skip sites covered in fewer samples than this")
    f.add_argument("--keep-invariant", action="store_true",
                   help="also write sites where every sample agrees")
    f.add_argument("-t", "--threads", type=int, default=1)
    f.set_defaults(func=cmd_freq)

    e = sub.add_parser("events", help="merge adjacent sites with identical carriers")
    e.add_argument("--freq", required=True, help="output of `freq`")
    e.add_argument("--max-gap", type=int, default=200,
                   help="sites this close with the same carriers merge (default 200)")
    e.add_argument("--min-freq", type=float, default=0.10)
    e.add_argument("-o", "--out", default="-")
    e.set_defaults(func=cmd_events)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
