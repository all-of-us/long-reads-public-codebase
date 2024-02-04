"""
Given a VCF, fill in the <DEL>, <INV>, <DUP> ALT sequence

DEL and INV are easy. But I gotta check the DUPs' coordinates. 
They'll have a span, but can I just insert that sequence?
Also, do I have to worry about dup-to-ins? 
"""

import sys
import pysam
import truvari

REF_CLEAN = True # Set to false if you're working with the right reference
MAX_SV = 100_000_000 # Filter things smaller than this

RC = str.maketrans("ATCG", "TAGC")
def do_rc(s):
    """
    Reverse complement a sequence
    """
    return s.translate(RC)[::-1]

def resolve(entry, ref):
    """
    """
    if entry.start > ref.get_reference_length(entry.chrom):
        return None
    if entry.alts[0] in ['<CNV>', '<INS>']:
        return None

    seq = ref.fetch(entry.chrom, entry.start, entry.stop)
    if entry.alts[0] == '<DEL>':
        entry.ref = seq
        entry.alts = [seq[0]]
    elif entry.alts[0] == '<INV>':
        entry.ref = seq
        entry.alts = [do_rc(seq)]
    elif entry.alts[0] == '<DUP>':
        entry.info['SVTYPE'] = 'INS'
        entry.ref = seq[0]
        entry.alts = [seq]
        entry.stop = entry.start + 1
    entry.qual = 1

    return entry

if __name__ == '__main__':
    default_quals =  {"pbsv": 3,
                      "sniffles": 2,
                      "pav": 4}

    vcf = pysam.VariantFile(sys.argv[1])
    ref = pysam.FastaFile(sys.argv[2])
    n_header = vcf.header.copy()
    d_qual = 1
    for key, val in default_quals.items():
        if key in sys.argv[1]:
            d_qual = default_quals[key]
    if REF_CLEAN:
        for ctg in vcf.header.contigs.keys():
            if ctg not in ref.references:
                n_header.contigs[ctg].remove_header()
    
    out = pysam.VariantFile("/dev/stdout", 'w', header=n_header)
    seen = set()
    for entry in vcf:
        key = truvari.entry_to_hash(entry)
        if key in seen:
            continue
        seen.add(key)

        if REF_CLEAN and entry.chrom not in ref.references:
            continue

        if truvari.entry_size(entry) >= MAX_SV:
            continue

        entry.qual = d_qual
        if entry.alts[0].startswith("<"):
            entry = resolve(entry, ref)

        if entry is None or set(entry.alts[0]) == {'N'}:
            continue
        if entry.info['SVTYPE'] != 'INV':
            entry.info['SVLEN'] = abs(len(entry.ref) - len(entry.alts[0]))
        else:
            entry.info['SVLEN'] = len(entry.ref)
        # No more blank genotypes
        n_gt = tuple([_ if _ is not None else 0 for _ in entry.samples[0]['GT']])
        # Preserve phasing informatino
        is_phased = entry.samples[0].phased
        entry.samples[0]['GT'] = n_gt
        entry.samples[0].phased = is_phased

        entry.translate(n_header)
        try:
            out.write(entry)
        except Exception:
            sys.stderr.write(f"{entry}\n{type(entry)}\n")
