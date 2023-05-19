#!/bin/bash
cp /Users/alexwang/0data/0mango/accessory/accessory_chrom.ID ./
cp /Users/alexwang/0data/0mango/accessory/accessory.genome.all.fasta ./
cp /Users/alexwang/0data/0mango/accessory/accessory.all.gtf ./
python -m jcvi.formats.gff bed --type=gene --key=gene_id accessory.all.gtf -o accessory.bed

## extract cds
#gffread accessory.all.gtf -g accessory.genome.all.fasta -x accessory.cds.fasta
#awk '/>/{split($0,a,"=");print ">"a[2];next}{print}' accessory.cds.fasta > accessory.cds
#rm accessory.cds.fasta

## extract pep
gffread accessory.all.gtf -g accessory.genome.all.fasta -y accessory_protein.fasta
awk '/>/{split($0,a,"=");print ">"a[2];next}{print}' accessory_protein.fasta > accessory.pep
rm accessory_protein.fasta

## get anchors
python -m jcvi.compara.catalog ortholog --dbtype=prot accessory accessory --cscore=.70
python -m jcvi.compara.synteny screen --minspan=10 --simple accessory.accessory.anchors accessory.accessory.anchors.new
