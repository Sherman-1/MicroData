#!/bin/bash

#PBS -N FOLDSEEK
#PBS -q bim
#PBS -l host=node04
#PBS -l ncpus=60
#PBS -l mem=400Gb
#PBS -l walltime=1000:00:00
#PBS -koed


foldseek="/store/EQUIPES/BIM/MEMBERS/ambre.baumann/Outils/foldseek/bin/foldseek"

cd /datas/SIMON/AFDB_20_to_100

if [ ! -d "globular" ]; then mkdir globular; fi
if [ ! -d "moltens" ]; then mkdir moltens; fi

rm -rf globular/*
rm -rf moltens/*

cp /store/EQUIPES/BIM/MEMBERS/simon.herman/MicroData/SCOPe_database/final_output/pdb_files/abcde/*.ent globular/
cp /store/EQUIPES/BIM/MEMBERS/simon.herman/MicroData/SCOPe_database/final_output/pdb_files/g/*.ent moltens/

if [ ! -d "database" ]; then mkdir database; fi
if [ ! -d "database/tmp" ]; then mkdir database/tmp; fi


# Targets database
$foldseek createdb pdbs_microproteins/ database/afdb_20_to_100_fs_db
$foldseek createindex database/afdb_20_to_100_fs_db database/tmp

# Globulars
$foldseek createdb globular/ database/globular_db
$foldseek createindex database/globular_db database/tmp

# Moltens
$foldseek createdb moltens/ database/moltens_db
$foldseek createindex database/moltens_db database/tmp


# Globular vs afdb_20_to_100
$foldseek easy-search database/globular_db database/afdb_20_to_100_fs_db Globular_vs_afbd_micro.tsv database/tmp -e 0.1 --num-iterations 2 --exhaustive-search 1 --format-output query,target,alntmscore,qtmscore,ttmscore,evalue,pident,qlen,qstart,qend,tlen,tstart,tend,qcov,tcov,prob,lddt,lddtfull --format-mode 4 --threads 60

# Moltens vs afdb_20_to_100
$foldseek easy-search database/moltens_db database/afdb_20_to_100_fs_db Molten_vs_afbd_micro.tsv database/tmp -e 0.1 --num-iterations 2 --exhaustive-search 1 --format-output query,target,alntmscore,qtmscore,ttmscore,evalue,pident,qlen,qstart,qend,tlen,tstart,tend,qcov,tcov,prob,lddt,lddtfull --format-mode 4 --threads 60

