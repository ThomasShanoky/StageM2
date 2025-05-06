# kmerator report
*date: 2025-04-16 14:19*  
*login: tlouvet*

**kmerator version:** 2.3.1

**Command:**

```
kmerator \
  --selection NPM1 DNMT3A FLT3 TET2 NRAS TP53 RUNX1 IDH2 ASXL1 WT1 KRAS IDH1 PTPN11 SRSF2 CEBPA KIT NF1 STAG2 GATA2 EZH2 BCOR JAK2 SMC1A RAD21 SF3B1 CBL \
  --datadir /scratch/indexes/kmerator/ \
  --genome /scratch/indexes/jellyfish/human/GRCh38_with_MT_canonical.jf \
  --specie homo_sapiens \
  --kmer-length 31 \
  --release 113 \
  --output Scripts/ExpressionKmers/GenesKmers \
  --thread 1 \
  --tmpdir /tmp/kmerator_no9zruxf \
  --assembly GRCh38
```

**Working directory:** `/scratch/users/tlouvet`

**Specie:** `homo_sapiens`

**Assembly:** `GRCh38`

**Transcriptome release:** `113`

**Genes/transcripts succesfully done (26)**

- NPM1: NPM1:ENST00000296930 - kmers/contigs: 161/9 (level: gene)
- DNMT3A: DNMT3A:ENST00000321117 - kmers/contigs: 9337/8 (level: gene)
- FLT3: FLT3:ENST00000241453 - kmers/contigs: 3791/2 (level: gene)
- TET2: TET2:ENST00000380013 - kmers/contigs: 9536/3 (level: gene)
- NRAS: NRAS:ENST00000369535 - kmers/contigs: 3838/21 (level: gene)
- TP53: TP53:ENST00000269305 - kmers/contigs: 2379/11 (level: gene)
- RUNX1: RUNX1:ENST00000675419 - kmers/contigs: 5908/4 (level: gene)
- IDH2: IDH2:ENST00000330062 - kmers/contigs: 2372/7 (level: gene)
- ASXL1: ASXL1:ENST00000375687 - kmers/contigs: 7022/1 (level: gene)
- WT1: WT1:ENST00000452863 - kmers/contigs: 2978/2 (level: gene)
- KRAS: KRAS:ENST00000311936 - kmers/contigs: 4503/41 (level: gene)
- IDH1: IDH1:ENST00000345146 - kmers/contigs: 1821/21 (level: gene)
- PTPN11: PTPN11:ENST00000351677 - kmers/contigs: 3614/97 (level: gene)
- SRSF2: SRSF2:ENST00000359995 - kmers/contigs: 1795/6 (level: gene)
- CEBPA: CEBPA:ENST00000498907 - kmers/contigs: 2571/1 (level: gene)
- KIT: KIT:ENST00000288135 - kmers/contigs: 5110/3 (level: gene)
- NF1: NF1:ENST00000358273 - kmers/contigs: 11557/45 (level: gene)
- STAG2: STAG2:ENST00000371145 - kmers/contigs: 6149/1 (level: gene)
- GATA2: GATA2:ENST00000341105 - kmers/contigs: 3346/3 (level: gene)
- EZH2: EZH2:ENST00000320356 - kmers/contigs: 2471/7 (level: gene)
- BCOR: BCOR:ENST00000378444 - kmers/contigs: 6865/2 (level: gene)
- JAK2: JAK2:ENST00000381652 - kmers/contigs: 6966/2 (level: gene)
- SMC1A: SMC1A:ENST00000322213 - kmers/contigs: 9662/4 (level: gene)
- RAD21: RAD21:ENST00000297338 - kmers/contigs: 3041/34 (level: gene)
- SF3B1: SF3B1:ENST00000335508 - kmers/contigs: 6221/6 (level: gene)
- CBL: CBL:ENST00000264033 - kmers/contigs: 10904/10 (level: gene)
