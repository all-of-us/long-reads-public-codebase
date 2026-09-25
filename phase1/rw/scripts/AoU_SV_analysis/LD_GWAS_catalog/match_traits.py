#!/usr/bin/env python3
"""Select GWAS traits that mention a disease or disorder term.

A trait matches if, after lower-casing and removing punctuation, it contains a term as
whole words (the term's last word may be plural). Traits containing "/" match if any
part does. README.md section 4, step 3.

Inputs:
  --ref   disease or disorder terms, one per line
  --gwas  GWAS trait names, one per line
Output:
  --out   the matching trait lines

Usage:
  python match_traits.py --ref disease_terms.txt --gwas gwas_traits.txt --out matched_traits.txt
"""
import argparse
import re
import unicodedata


def normalize(text: str) -> str:
    text = unicodedata.normalize("NFKD", text.lower().strip())
    return re.sub(r"[^\w\s]", "", text)

def load_trait_patterns(traits_file: str):
    regex_patterns = []

    with open(traits_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            trait = normalize(line)
            words = trait.split()
            if not words:
                continue
            # Allow the last word of the term to be plural.
            words[-1] = rf"{words[-1]}s?"
            pattern_str = rf"\b{' '.join(words)}\b"
            pattern = re.compile(pattern_str, re.IGNORECASE)
            regex_patterns.append(pattern)
    return regex_patterns


def match_traits(gwas_file: str, regex_patterns):
    matched_lines = []
    with open(gwas_file, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            sub_traits = [normalize(part) for part in line.split("/")]
            for sub_trait in sub_traits:
                if any(p.search(sub_trait) for p in regex_patterns):
                    matched_lines.append(line)
                    break
    return matched_lines


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref", "-r", required=True, help="Reference file containing disease or disorder terms.")
    parser.add_argument("--gwas", "-g", required=True, help="Identified GWAS traits")
    parser.add_argument("--out", "-o", required=True, help="Output file for matched GWAS traits.")
    return parser.parse_args()


def main():
    args = parse_args()
    regex_patterns = load_trait_patterns(args.ref)
    matched_lines = match_traits(args.gwas, regex_patterns)
    with open(args.out, "w") as f:
        f.write("\n".join(matched_lines))

    print(f"Total matched lines: {len(matched_lines)}")
    print(f"Results written to: {args.out}")


if __name__ == "__main__":
    main()
