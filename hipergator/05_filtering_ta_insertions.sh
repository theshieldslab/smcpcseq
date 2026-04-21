#!/usr/bin/bash

scriptPWD=${PWD}

cd /blue/robert.shields/arwalker/UA159

### Create annotation table
awk -F '\t' '
$3 == "CDS" {
    locus_tag = name = annotation = ".";
    split($9, attrs, ";");
    for (i in attrs) {
        if (attrs[i] ~ /locus_tag=/) {
            split(attrs[i], x, "="); locus_tag = x[2];
        } else if (attrs[i] ~ /^gene=/) {
            split(attrs[i], x, "="); name = x[2];
        } else if (attrs[i] ~ /^product=/ || attrs[i] ~ /^Note=/) {
            split(attrs[i], x, "="); annotation = x[2];
        }
    }
    # fallback: if no gene= was found, use locus_tag
    if (name == ".") {
        name = locus_tag;
    }
    start = $4;
    end = $5;
    strand = $7;
    gene_length = end - start + 1;
    print locus_tag "\t" name "\t" start "\t" end "\t" gene_length "\t" strand "\t" annotation;
}' /blue/robert.shields/anasolanomorales/references/UA159.gff > UA159_cds_annotations.tsv

cd ${scriptPWD}
