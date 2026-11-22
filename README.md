# repeatmasker-annotation

Construire, à partir d'une seule annotation RepeatMasker, des fichiers de référence **cohérents** pour deux outils de quantification locus-spécifique des éléments transposables (ET) en single-cell : [Stellarscope](https://github.com/nixonlab/stellarscope) et [MATES](https://github.com/mcgilldinglab/MATES).

## Le problème

Pour comparer deux outils de quantification, il faut qu'ils travaillent sur **exactement les mêmes loci, avec les mêmes identifiants**. Sinon, une différence de résultat peut venir de l'annotation, et non des algorithmes.

Or, dans le GTF RepeatMasker, `transcript_id` désigne une **sous-famille** (par exemple `AluY`), pas un locus. Sur l'annotation hg38 utilisée lors de mon stage, on compte **5 683 690 loci pour seulement 515 262 `transcript_id` distincts**, soit 91 % de doublons. Impossible, dans ces conditions, de relier un locus d'un outil au même locus dans l'autre.

Les deux outils attendent en plus des formats différents :

| Outil | Format attendu | Coordonnées |
|---|---|---|
| Stellarscope | GTF avec un attribut `locus` (valeur par défaut de `--attribute`) | 1-based, fin incluse |
| MATES | CSV et BED à huit colonnes | 0-based, fin exclue |

## La solution

1. Chaque locus reçoit un identifiant unique : `NomET:chrom:start:end:strand` (par exemple `AluY:chr1:1001:1300:+`).
2. Le filtrage est fait **une seule fois** : annotations incomplètes, intervalles ou brins invalides, loci dupliqués. Chaque retrait est compté et affiché, aucun locus ne disparaît sans trace.
3. Les trois fichiers sont écrits à partir de **la même table filtrée**, avec conversion des coordonnées pour MATES (`start - 1`).

```
GTF RepeatMasker ──► lecture ──► filtrage + locus_id unique ──┬──► GTF Stellarscope
                                                               ├──► CSV MATES
                                                               └──► BED MATES
```

## Installation

```bash
conda env create -f environment.yml
conda activate repeat-harmonizer
```

## Utilisation

```bash
python build_repeat_annotations.py <gtf_source> <out_gtf> <out_csv> <out_bed>
```

Exemple sur le jeu de test fourni :

```bash
python build_repeat_annotations.py tests/data/mini_rmsk.gtf \
    examples/stellarscope_TE.gtf examples/mates_TE.csv examples/mates_TE.bed
```

```
[INFO] 9 loci lus, 4 transcript_id distincts (56 % de doublons d'identifiant)
[INFO] 1 loci retirés : annotation incomplète
[INFO] 1 loci retirés : intervalle ou brin invalide
[INFO] 1 loci retirés : locus dupliqué
[OK] 6 loci écrits dans les trois fichiers
```

Les sorties obtenues sont dans [`examples/`](examples/). Le script a quatre arguments de sortie explicites pour s'intégrer simplement dans une règle Snakemake.

## Tests

```bash
pytest
```

Les tests vérifient notamment que les trois fichiers contiennent les mêmes loci, que chaque identifiant est unique et que la conversion des coordonnées est correcte. Le jeu de test (`tests/data/mini_rmsk.gtf`) est **synthétique** : il reproduit les cas rencontrés dans une vraie annotation (identifiants partagés, doublon, annotation incomplète, brin invalide).

## Contexte

Ce script a été écrit pendant mon stage de Master 2 en bioinformatique, consacré à la comparaison de méthodes de quantification des éléments transposables en scRNA-seq. Cette version publique a été nettoyée et testée. Elle ne contient aucune donnée ni information propre au laboratoire.

## Auteure

Aïnata Zida

## Licence

MIT

