"""Fetch influenza A(H3N2) HA coding sequences from GenBank for Dataset B."""
import argparse
import csv
import re
import sys
import time
from collections import Counter
from datetime import date
from pathlib import Path

from Bio import Entrez, SeqIO

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "public_genomic"
MIN_INTERVAL = 0.4  # NCBI policy: no more than 3 requests per second without an API key
_last_call = 0.0


def throttled(fn, **kw):
    global _last_call
    wait = MIN_INTERVAL - (time.monotonic() - _last_call)
    if wait > 0:
        time.sleep(wait)
    _last_call = time.monotonic()
    handle = fn(**kw)
    try:
        return handle.read() if fn is Entrez.efetch else Entrez.read(handle)
    finally:
        handle.close()


def season_window(season):
    # Northern-hemisphere season "2017-2018" = 1 Oct 2017 .. 30 Apr 2018
    y1, y2 = (int(x) for x in season.split("-"))
    if y2 != y1 + 1:
        sys.exit(f"season must look like 2017-2018, got {season}")
    return f"{y1}-10", f"{y2}-04"


def parse_month(collection_date):
    # GenBank dates come as 2018, Jan-2018, 12-Jan-2018 or ISO 2018-01-12; month precision is required
    s = collection_date.strip()
    m = re.fullmatch(r"(\d{4})-(\d{2})(?:-\d{2})?", s)
    if m:
        return f"{m[1]}-{m[2]}"
    months = "Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec".split()
    m = re.fullmatch(r"(?:\d{1,2}-)?([A-Z][a-z]{2})-(\d{4})", s)
    if m and m[1] in months:
        return f"{m[2]}-{months.index(m[1]) + 1:02d}"
    return None


def ha_cds(record):
    for f in record.features:
        if f.type != "CDS":
            continue
        gene = f.qualifiers.get("gene", [""])[0].upper()
        product = f.qualifiers.get("product", [""])[0].lower()
        if gene == "HA" or "hemagglutinin" in product:
            return f
    return None


def describe(record):
    src = next((f for f in record.features if f.type == "source"), None)
    q = src.qualifiers if src else {}
    loc = q.get("geo_loc_name", q.get("country", [""]))[0]  # INSDC renamed /country to /geo_loc_name
    cdate = q.get("collection_date", [""])[0]
    return {
        "accession": record.id,
        "strain": q.get("strain", [""])[0],
        "serotype": q.get("serotype", [""])[0],
        "host": q.get("host", [""])[0],
        "collection_date": cdate,
        "month": parse_month(cdate) if cdate else None,
        "country": loc.split(":")[0].strip(),
    }


def fetch_genbank(ids, batch):
    text = []
    for i in range(0, len(ids), batch):
        chunk = ids[i:i + batch]
        text.append(throttled(Entrez.efetch, db="nucleotide", id=",".join(chunk), rettype="gb", retmode="text"))
        print(f"  fetched {min(i + batch, len(ids))}/{len(ids)}", file=sys.stderr)
    return "".join(text)


def parse_records(gb_text, raw_path):
    raw_path.parent.mkdir(parents=True, exist_ok=True)
    raw_path.write_text(gb_text, encoding="utf-8")
    with raw_path.open(encoding="utf-8") as fh:
        return list(SeqIO.parse(fh, "genbank"))


def candidates(records, lo, hi, log):
    keep, reasons = [], Counter()
    for rec in records:
        info = describe(rec)
        f = ha_cds(rec)
        if info["serotype"].upper() != "H3N2":
            reasons["serotype not H3N2"] += 1
        elif not info["host"].lower().startswith(("homo sapiens", "human")):
            # serotype alone admits swine, canine and avian H3N2, whose HA lineages are far from seasonal human H3N2
            reasons["host not human"] += 1
        elif f is None:
            reasons["no HA CDS feature"] += 1
        elif "<" in str(f.location) or ">" in str(f.location):
            reasons["CDS marked partial"] += 1
        elif info["month"] is None:
            reasons["collection date lacks month"] += 1
        elif not (lo <= info["month"] <= hi):
            reasons["collected outside season window"] += 1
        else:
            info["seq"] = str(f.extract(rec.seq)).upper()
            keep.append(info)
    log.append(f"- records parsed: {len(records)}")
    for r, c in sorted(reasons.items()):
        log.append(f"- excluded, {r}: {c}")
    log.append(f"- in-season complete HA CDS candidates: {len(keep)}")
    return keep


def spread(items, n):
    # evenly spaced picks over the accession-sorted list, so one large submission batch cannot dominate
    items = sorted(items, key=lambda x: x["accession"])
    if len(items) <= n:
        return items
    return [items[(i * len(items)) // n] for i in range(n)]


def write_outputs(rows, log, log_name="selection_log.md"):
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "accessions.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["accession", "strain", "year", "country", "length", "collection_date"])
        for r in rows:
            w.writerow([r["accession"], r["strain"], r["month"][:4], r["country"], len(r["seq"]), r["collection_date"]])
    with (OUT / "h3n2_ha_cds.fasta").open("w", encoding="utf-8") as fh:
        for r in rows:
            fh.write(f">{r['accession']} {r['strain']} | {r['collection_date']} | {r['country']}\n{r['seq']}\n")
    (OUT / log_name).write_text("\n".join(log) + "\n", encoding="utf-8")


def cmd_select(args):
    lo, hi = season_window(args.season)
    y1, y2 = args.season.split("-")
    # publication-date window is wider than the season because sequences are deposited after collection
    query = (f'"Influenza A virus"[Organism] AND H3N2[All Fields] AND hemagglutinin[Title] '
             f'AND "complete cds"[Title] AND ("{y1}/10/01"[PDAT] : "{int(y2) + 1}/12/31"[PDAT])')
    log = [f"# Dataset B selection log", "",
           f"- run date: {date.today().isoformat()}", f"- mode: select", f"- season window (collection month): {lo} .. {hi}",
           f"- Entrez query: `{query}`"]
    # esearch returns at most 10000 ids per call, so page until every hit is collected
    res = throttled(Entrez.esearch, db="nucleotide", term=query, retmax=10000, idtype="acc")
    all_ids = list(res["IdList"])
    while len(all_ids) < int(res["Count"]):
        page = throttled(Entrez.esearch, db="nucleotide", term=query, retmax=10000, retstart=len(all_ids), idtype="acc")
        if not page["IdList"]:
            break
        all_ids += page["IdList"]
    all_ids = sorted(set(all_ids))
    if not all_ids:
        sys.exit("esearch returned no records; stopping")
    # sort locally and thin evenly, so the candidate pool does not depend on NCBI's default result order
    ids = [all_ids[(i * len(all_ids)) // args.max_candidates] for i in range(min(args.max_candidates, len(all_ids)))]
    log.append(f"- esearch total hits: {res['Count']}; accessions returned: {len(all_ids)}; "
               f"examined {len(ids)} evenly spaced in accession order")
    records = parse_records(fetch_genbank(ids, args.batch), OUT / "raw" / "candidates.gb")
    keep = candidates(records, lo, hi, log)
    if not keep:
        sys.exit("no candidates survived filtering; stopping")
    lengths = Counter(len(k["seq"]) for k in keep)
    modal = lengths.most_common(1)[0][0]
    log.append(f"- CDS length distribution among candidates: {dict(sorted(lengths.items()))}")
    keep = [k for k in keep if len(k["seq"]) == modal]
    log.append(f"- kept modal length {modal} nt (no alignment needed): {len(keep)}")
    rows = spread(keep, args.n)
    log.append(f"- selected by even spacing over accession order: {len(rows)}")
    write_outputs(rows, log)
    print("\n".join(log))


def cmd_refetch(args):
    with (OUT / "accessions.csv").open(encoding="utf-8") as fh:
        wanted = [r["accession"] for r in csv.DictReader(fh)]
    records = parse_records(fetch_genbank(wanted, args.batch), OUT / "raw" / "refetch.gb")
    by_id = {r.id: r for r in records}
    missing = [a for a in wanted if a not in by_id]
    if missing:
        sys.exit(f"accessions not returned by NCBI: {missing}")
    rows = []
    for a in wanted:
        info = describe(by_id[a])
        info["seq"] = str(ha_cds(by_id[a]).extract(by_id[a].seq)).upper()
        rows.append(info)
    write_outputs(rows, [f"# Dataset B refetch log", "", f"- run date: {date.today().isoformat()}",
                         f"- mode: refetch of {len(wanted)} accessions from accessions.csv"], "refetch_log.md")
    print(f"refetched {len(rows)} records")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--email", required=True, help="contact email sent to NCBI with each request")
    p.add_argument("--batch", type=int, default=100)
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("select", help="search GenBank and choose sequences")
    s.add_argument("--season", default="2017-2018")
    s.add_argument("--n", type=int, default=36)
    s.add_argument("--max-candidates", type=int, default=600)
    sub.add_parser("refetch", help="re-download exactly the accessions in accessions.csv")
    args = p.parse_args()
    Entrez.email = args.email
    Entrez.tool = "genomath_dm_course_project"
    {"select": cmd_select, "refetch": cmd_refetch}[args.cmd](args)


if __name__ == "__main__":
    main()
