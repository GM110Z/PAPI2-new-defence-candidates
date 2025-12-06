**operon-finder.py** : Uses the _operon.tsv file from https://github.com/GCA-VH-lab/FlaGs2 to only extract queries that are in multi-gene operons
Runs as: python operon-finder.py flags_file_operon.tsv --max-gap 30 --require-same-strand --out result.tsv

**jsontomatrix.py**: Uses the session.json file from clinker to produce a presence absence matrix of genes that are shared or unique across pairwise genome comparisons
Runs as: python jsontomatrix.py session.json outfile.csv
