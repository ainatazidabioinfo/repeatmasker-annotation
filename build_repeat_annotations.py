#!/usr/bin/env python3
# ============================================================
#  build_repeat_annotations.py
#
#  Construit, à partir d'un seul GTF RepeatMasker, des fichiers
#  d'annotation cohérents pour deux outils de quantification
#  des éléments transposables (ET) en single-cell :
#    - un GTF pour Stellarscope ;
#    - un CSV et un BED pour MATES.
#
#  Problème résolu : dans le GTF RepeatMasker, transcript_id
#  désigne une sous-famille (ex. AluY) et non un locus. Des
#  milliers de loci partagent donc le même identifiant, ce qui
#  empêche de comparer les sorties des deux outils locus par
#  locus. Le script attribue à chaque locus un identifiant
#  unique : NomET:chrom:start:end:strand.
#
#  Garantie : les trois fichiers sont écrits à partir de la MÊME
#  table filtrée. Les deux outils travaillent donc exactement sur
#  les mêmes loci, avec les mêmes identifiants.
#
#  Auteure : Aïnata Zida
# ============================================================

import argparse
import os
import sys

import pandas as pd


# ------------------------------------------------------------
# Lecture du GTF source
# ------------------------------------------------------------

def parse_gtf(gtf_path):
    """Lit le GTF RepeatMasker et renvoie une ligne par locus.

    Coordonnées conservées telles quelles : 1-based, fin incluse
    (convention GTF).
    """
    if not os.path.exists(gtf_path):
        sys.exit(f"[ERREUR] GTF source introuvable : {gtf_path}")

    rows = []
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, start, end = fields[0], int(fields[3]), int(fields[4])
            score, strand, attrs = fields[5], fields[6], fields[8]

            values = {"gene_id": "", "transcript_id": "",
                      "repClass": "", "repFamily": ""}
            for part in attrs.split(";"):
                part = part.strip()
                if not part or '"' not in part:
                    continue
                key = part.split()[0]
                if key in values:
                    values[key] = part.split('"')[1]

            rows.append((chrom, start, end, score, strand,
                         values["gene_id"], values["transcript_id"],
                         values["repClass"], values["repFamily"]))

    df = pd.DataFrame(rows, columns=[
        "chrom", "start", "end", "score", "strand",
        "gene_id", "transcript_id", "rep_class", "rep_family",
    ])
    return df


# ------------------------------------------------------------
# Filtrage commun et identifiant unique
# ------------------------------------------------------------

def build_locus_table(df):
    """Filtre les loci une seule fois et crée l'identifiant unique.

    Renvoie (table filtrée, dictionnaire des loci retirés par motif).
    """
    removed = {}
    n0 = len(df)

    # 1. Annotations incomplètes (MATES a besoin de la famille)
    complete = (df[["chrom", "gene_id", "rep_family"]] != "").all(axis=1)
    removed["annotation incomplète"] = int((~complete).sum())
    df = df[complete]

    # 2. Intervalles ou brins invalides
    valid = (df["start"] >= 1) & (df["end"] >= df["start"]) \
        & df["strand"].isin(["+", "-"])
    removed["intervalle ou brin invalide"] = int((~valid).sum())
    df = df[valid]

    # 3. Identifiant unique, puis suppression des doublons exacts
    df = df.copy()
    df["locus_id"] = (df["gene_id"] + ":" + df["chrom"] + ":"
                      + df["start"].astype(str) + ":"
                      + df["end"].astype(str) + ":" + df["strand"])
    n_before = len(df)
    df = df.drop_duplicates(subset="locus_id").reset_index(drop=True)
    removed["locus dupliqué"] = n_before - len(df)

    assert sum(removed.values()) == n0 - len(df)
    return df, removed


# ------------------------------------------------------------
# Écriture des fichiers
# ------------------------------------------------------------

def _ensure_dir(path):
    """Crée le dossier parent si besoin (rien à faire s'il n'y en a pas)."""
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)


def write_stellarscope_gtf(loci, output_path, source="RepeatMasker"):
    """GTF pour Stellarscope.

    L'attribut `locus` est obligatoire : c'est la valeur par défaut
    de --attribute dans `stellarscope assign` (hérité de Telescope).
    """
    _ensure_dir(output_path)
    attrs = (
        'gene_id "' + loci["gene_id"] + '"; '
        'transcript_id "' + loci["locus_id"] + '"; '
        'locus "' + loci["locus_id"] + '"; '
        'repClass "' + loci["rep_class"] + '"; '
        'repFamily "' + loci["rep_family"] + '";'
    )
    # Écriture ligne à ligne : passer par to_csv(quoting=QUOTE_NONE,
    # escapechar=...) peut échapper les guillemets (\"AluY\") selon la
    # version de pandas/Python, ce qui rend le GTF invalide.
    lines = (loci["chrom"] + "\t" + source + "\texon\t"
             + loci["start"].astype(str) + "\t" + loci["end"].astype(str)
             + "\t" + loci["score"] + "\t" + loci["strand"] + "\t.\t" + attrs)
    with open(output_path, "w") as f:
        f.write("\n".join(lines) + "\n")


def build_mates_reference(loci):
    """Référence attendue par MATES (0.1.8), huit colonnes :
    chromosome, start, end, TE_Name, index, strand, TE_Fam, length.

    Conversion des coordonnées :
      GTF  : 1-based, fin incluse
      BED  : 0-based, fin exclue  -> start - 1, end inchangé
    """
    mates = pd.DataFrame({
        "chromosome": loci["chrom"],
        "start": loci["start"] - 1,
        "end": loci["end"],
        "TE_Name": loci["gene_id"],      # sous-famille : L1P5, AluY...
        "index": loci.index.astype(int),  # identifiant numérique exigé par MATES
        "strand": loci["strand"],
        "TE_Fam": loci["rep_family"],     # famille : L1, Alu, ERVL...
    })
    mates["length"] = mates["end"] - mates["start"]
    return mates


def write_mates(loci, csv_path, bed_path):
    """CSV (étapes de traitement MATES) et BED (split_count_10X_data),
    avec les mêmes huit colonnes."""
    mates = build_mates_reference(loci)
    _ensure_dir(csv_path)
    mates.to_csv(csv_path, sep=",", index=False, header=False)
    _ensure_dir(bed_path)
    mates.to_csv(bed_path, sep="\t", index=False, header=False)
    return mates


# ------------------------------------------------------------
# Point d'entrée
# ------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Harmonise un GTF RepeatMasker pour Stellarscope et MATES "
                    "(identifiant unique par locus, mêmes loci partout).")
    parser.add_argument("gtf_source", help="GTF RepeatMasker d'entrée")
    parser.add_argument("out_gtf", help="GTF de sortie pour Stellarscope")
    parser.add_argument("out_csv", help="CSV de sortie pour MATES")
    parser.add_argument("out_bed", help="BED de sortie pour MATES")
    parser.add_argument("--source", default="RepeatMasker",
                        help="valeur de la colonne 2 du GTF de sortie")
    args = parser.parse_args(argv)

    raw = parse_gtf(args.gtf_source)
    if raw.empty:
        sys.exit("[ERREUR] Aucun locus lu dans le GTF source")

    n_ids = raw["transcript_id"].nunique()
    print(f"[INFO] {len(raw)} loci lus, {n_ids} transcript_id distincts "
          f"({100 * (1 - n_ids / len(raw)):.0f} % de doublons d'identifiant)")

    loci, removed = build_locus_table(raw)
    for reason, n in removed.items():
        if n:
            print(f"[INFO] {n} loci retirés : {reason}")
    if loci.empty:
        sys.exit("[ERREUR] Aucun locus valide après filtrage")

    write_stellarscope_gtf(loci, args.out_gtf, source=args.source)
    write_mates(loci, args.out_csv, args.out_bed)
    print(f"[OK] {len(loci)} loci écrits dans les trois fichiers")


if __name__ == "__main__":
    main()
