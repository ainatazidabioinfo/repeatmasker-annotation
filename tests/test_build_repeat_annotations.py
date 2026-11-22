"""Tests du script build_repeat_annotations.py (lancer : pytest)."""

import os
import sys

import pandas as pd

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, ".."))

import build_repeat_annotations as bra  # noqa: E402

GTF = os.path.join(HERE, "data", "mini_rmsk.gtf")


def run(tmp_path):
    out_gtf = tmp_path / "out" / "stellarscope.gtf"
    out_csv = tmp_path / "out" / "mates.csv"
    out_bed = tmp_path / "out" / "mates.bed"
    bra.main([GTF, str(out_gtf), str(out_csv), str(out_bed)])
    gtf = pd.read_csv(out_gtf, sep="\t", header=None)
    csv = pd.read_csv(out_csv, header=None)
    bed = pd.read_csv(out_bed, sep="\t", header=None)
    return gtf, csv, bed


def test_lecture_et_doublons_d_identifiant():
    raw = bra.parse_gtf(GTF)
    assert len(raw) == 9
    assert raw["transcript_id"].nunique() == 4


def test_filtrage(tmp_path):
    loci, removed = bra.build_locus_table(bra.parse_gtf(GTF))
    assert removed == {"annotation incomplète": 1,
                       "intervalle ou brin invalide": 1,
                       "locus dupliqué": 1}
    assert len(loci) == 6
    assert loci["locus_id"].is_unique


def test_memes_loci_dans_les_trois_fichiers(tmp_path):
    gtf, csv, bed = run(tmp_path)
    assert len(gtf) == len(csv) == len(bed) == 6
    # Les fichiers MATES sont identiques au séparateur près
    assert csv.equals(bed)
    # Chaque locus du GTF se retrouve dans MATES après conversion
    gtf_keys = set(zip(gtf[0], gtf[3] - 1, gtf[4], gtf[6]))
    mates_keys = set(zip(csv[0], csv[1], csv[2], csv[5]))
    assert gtf_keys == mates_keys


def test_conversion_des_coordonnees(tmp_path):
    gtf, csv, _ = run(tmp_path)
    # Premier locus : 1001-1300 en GTF (1-based, fin incluse)
    assert (gtf.iloc[0, 3], gtf.iloc[0, 4]) == (1001, 1300)
    # devient 1000-1300 en BED (0-based, fin exclue), longueur 300
    assert (csv.iloc[0, 1], csv.iloc[0, 2], csv.iloc[0, 7]) == (1000, 1300, 300)


def test_attribut_locus_pour_stellarscope(tmp_path):
    gtf, _, _ = run(tmp_path)
    attrs = gtf.iloc[0, 8]
    assert 'locus "AluY:chr1:1001:1300:+"' in attrs
    assert 'transcript_id "AluY:chr1:1001:1300:+"' in attrs
