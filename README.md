# Introduction
This repository contains three historical versions of curated maize TE libraries derived from the Maize TE Consortium (MTEC). I combined the three together and further clean the combined library with the following commands and curations. If you are looking for a comprehensive and high-quality maize TE library, look no further, this is the one (usually named like "maizeTE02052020" in the root directory).

## Files
- `maizeTE10102014` was download from the MTEC official website (http://maizesequence.org). The website is gone, but I managed to get a mirror of the first page. Please refer to the file [history/Maizedatabase_mirror.pdf](https://github.com/oushujun/MTEC/blob/master/history/Maizedatabase_mirror.pdf) for more information about the MTEC project.
- `TE_12-Feb-2015_15-35.fa` was shared by [Nicolas Blavet](https://github.com/blavetn) from `https://github.com/mcstitzer/maize_v4_TE_annotation/issues/9`.
- `Wessler-Bennetzen_2.fasta` was used to annotate the initial B73 genome ([Schnable et al. 2009](https://science.sciencemag.org/content/326/5956/1112.full)), which was shared by [Kapeel Chougule](https://scholar.google.com/citations?user=-ZMiorYAAAAJ&hl=en). I believe this is an earlier version of the MTEC library.
- `nonTE.repeat.fa` contains 5 non-TE repeats (knob180, knob TR-1, rDNA spacer, subtelomere 4-12-1, and CentC) in maize, which was shared by [Jianing Liu](https://www.genetics.uga.edu/directory/jianing-liu).


## Combine the three MTEC libraries + nonTE repeats
### 1. Reformat sequence IDs
`for i in history/Wessler-Bennetzen_2.fasta history/maizeTE10102014 history/TE_12-Feb-2015_15-35.fa; do perl -nle 's/\s+$//g; $_=(split)[0]; s/\-/_/g; print $_' $i > $i.mod; done`

### 2. Combine sequences with unique IDs
`perl bin/output_by_list.pl 1 <(cat history/*.mod) 1 <(cat history/*.mod|grep \>|sort -u) -FA > history/maizeTE11212019.ori`

### 3. Split the library into consensus and others
`perl bin/output_by_list.pl 1 history/maizeTE11212019.ori 1 <(grep consen history/maizeTE11212019.ori) -FA > history/maizeTE11212019.ori.consensus`

`perl bin/output_by_list.pl 1 history/maizeTE11212019.ori 1 <(grep consen history/maizeTE11212019.ori) -FA -ex > history/maizeTE11212019.ori.others`

### 4. Remove TEs in others that are represented by consensus TEs
`RepeatMasker -pa 36 -div 40 -lib history/maizeTE11212019.ori.consensus -cutoff 225 history/maizeTE11212019.ori.others`

`perl bin/make_masked.pl -rmout history/maizeTE11212019.ori.others.out -genome history/maizeTE11212019.ori.others -maxdiv 20 -minscore 200 -minlen 80 -t 30`

`perl bin/cleanup_tandem.pl -nc 1000 -nr 0.5 -minlen 80 -cleanN 1 -cleanT 1 -trf 0 -f history/maizeTE11212019.ori.others.new.masked > history/maizeTE11212019.ori.others.new.masked.cln`

`cat history/maizeTE11212019.ori.consensus history/maizeTE11212019.ori.others.new.masked.cln > history/maizeTE11212019.ori2`

### 5. Remove redundant sequences
`perl bin/cleanup_nested.pl -in history/maizeTE11212019.ori2 -cov 0.98 -minlen 80 -miniden 95 -iter 2 -t 36`

### 6. Remove nonTE repeats and tandem repeats
`RepeatMasker -pa 36 -div 40 -no_is -norna -nolow -lib history/nonTE.repeat.fa -cutoff 225 history/maizeTE11212019.ori2.cln`

`perl bin/cleanup_tandem.pl -nc 1000 -nr 0.5 -minlen 80 -cleanN 1 -cleanT 1 -trf 1 -f history/maizeTE11212019.ori2.cln.masked > history/maizeTE11212019.ori2.cln2`

`cat history/nonTE.repeat.fa history/maizeTE11212019.ori2.cln2 > history/maizeTE11212019.ori3`


## Improve the combined library
### 1. Reclassify unknown TEs
`python2 TEsorter.py history/maizeTE11212019.ori3 -p 36`

### 2. Find misclassified entries
The file `history/maizeTE11212019.ori3.rexdb.cls.tsv` contains new classifications of the library. Most of them are consistent with the old classification. What really improved are the LTR/unknown classification.

`perl -nle '($info, $cla)=(split)[0,2]; my $oldcla=$1 if $info=~/^([A-Z]+)_/; $cla=~s/EnSpm_CACTA/DTC/; $cla=~s/hAT/DTA/; $cla=~s/PIF_Harbinger/DTH/; $cla=~s/MuDR_Mutator/DTM/; $cla=~s/Tc1_Mariner/DTT/; $cla=~s/Gypsy/RLG/; $cla=~s/Copia/RLC/; print "$oldcla\t$cla\t$info" if $cla ne $oldcla' history/maizeTE11212019.ori3.rexdb.cls.tsv |less`

### 3. Some LTRs appear to have the same name but different classifications (RLG/RLC/RLX)
`grep RL history/maizeTE11212019.ori3|perl -nle 's/RL._//; print $_'|sort|uniq -c |sort -k1,1|tac|less
blastn -query list.fa -subject list.fa -outfmt=6 > list.fa.out`

### 4. Manually check misclassified sequences (`list.fa`).
These entries are put in the `history/removal.list` and removed:

|Seq_ID|Removal reason|
| ----------- | ----------- |
|RLC_chr3_D_28761151|rDNA-contained|
|DTM_Zm08959_AC199876_1|LTRcoding-contained|
|DTM_Zm22805IC_AC207689_1|LINE-contained|
|RIX_nugimu_AC203843_0|Duplicted_with_RIL_nugimu_AC203843_0|
|RLX_fageri_AC204875_8470|misclassified_as_LINE|
|DTA_ZM00171_consensus|misclassified_as_CACTA|
|DTA_ZM00205_consensus|misclassified_as_CACTA|
|DTA_ZM00284_consensus|misclassified_as_CACTA|
|RLX_teki_AC202867-7492|rDNA-contained|
|RLG_ajajog_AC191578_3186|A_RLG_nested_in_RLC_ajajog_AC191578_3186|
|RLC_iwim_AC203300_7761|misclassified_RLG_duplicated|
|RLC_kupu_AC216069_13264|misclassified_RLG_duplicated|
|RLX_pute_AC197188_5467|duplicated_RLC_pute_AC197188_5467|
|RLX_votaed_AC215881_13209|duplicated_RLC_votaed_AC215881_13209|
|RLC_votaed_AC215881_13209|5-6_LTR_nested_together|
|RLX_bobeg_AC193485_3670|5_LTR_nested_together|

`perl bin/output_by_list.pl 1 history/maizeTE11212019.ori3 1 history/removal.list -FA -ex > history/maizeTE11212019.ori3.cln`

### 5. PPP_PPO_AC185414 is changed to DTH_PPO_AC185414 manually

### 6. Update LTR classifications
`perl -nle '($info, $cla)=(split)[0,2]; my $oldcla=$1 if $info=~/^([A-Z]+)_/; $cla=~s/EnSpm_CACTA/DTC/; $cla=~s/hAT/DTA/; $cla=~s/PIF_Harbinger/DTH/; $cla=~s/MuDR_Mutator/DTM/; $cla=~s/Tc1_Mariner/DTT/; $cla=~s/Gypsy/RLG/; $cla=~s/Copia/RLC/; next unless /LTR/; my $info_new=$info; $info_new=~s/$oldcla/$cla/; print "$info|$info_new" if $cla ne $oldcla' history/maizeTE11212019.ori3.rexdb.cls.tsv > history/maizeTE11212019.ori3.rexdb.cls.tsv.LTR`

`for i in `cat history/maizeTE11212019.ori3.rexdb.cls.tsv.LTR`; do perl -i -slane 'my ($old, $new)=(split /\|/, $info); s/$old/$new/; print $_' -- -info=$i history/maizeTE11212019.ori3.cln; done`

### 7. Convert sequence names to RepeatMasker format
`perl -nle 'my $id=(split)[0]; $id=~s/RLC_(.*)/$1#LTR\/Copia/; $id=~s/RLG_(.*)/$1#LTR\/Gypsy/; $id=~s/RLX_(.*)/$1#LTR\/unknown/; $id=~s/DHH_(.*)/$1#DNA\/Helitron/; $id=~s/DTA_(.*)/$1#DNA\/DTA/; $id=~s/DTC_(.*)/$1#DNA\/DTC/; $id=~s/DTH_(.*)/$1#DNA\/DTH/; $id=~s/DTM_(.*)/$1#DNA\/DTM/; $id=~s/DTT_(.*)/$1#DNA\/DTT/; $id=~s/(RIT_.*)/$1#LINE\/RTE/; $id=~s/(RIL_.*)/$1#LINE\/L1/; $id=~s/(RIX_.*)/$1#LINE\/unknown/; $id=~s/(ZM_CACTA_noncoding.*)/$1#MITE\/DTC/; $id=~s/(ZM_Stowaway.*)/$1#DNA\/DTT/; $id=~s/(ZM_Tourist.*)/$1#DNA\/DTH/; $id=~s/(ZM_hAT_noncoding.*)/$1#MITE\/DTA/; $id=~s/(RST_.*)/$1#SINE\/tRNA/; print $id' history/maizeTE11212019.ori3.cln > history/maizeTE11222019.ori`

### 8. Rename short TIR (<= 600bp) to MITE
`perl bin/rename_MITE.pl history/maizeTE11222019.ori > history/maizeTE11222019.ori.rename`


## Remove gene sequences
### 1. Mask gene cds
`RepeatMasker -pa 36 -q -no_is -norna -nolow -div 40 -lib history/Zea_mays.AGPv4.cds.all.noTE.fa.mod.cln -cutoff 500 history/maizeTE11222019.ori.rename`

`perl bin/cleanup_tandem.pl -misschar n -Nscreen 1 -nc 1000 -nr 0.3 -minlen 80 -maxlen 5000000 -cleanN 1 -cleanT 1 -trf 0 -f history/maizeTE11222019.ori.rename.masked > history/maizeTE11222019.ori.rename.nogene`

### 2. Finalize
Manually add the sequence CL569186.1#subtelomere/4-12-1 back to `history/maizeTE11222019.ori.rename.nogene`

`cp history/maizeTE11222019.ori.rename.nogene maizeTE11222019`


## Updates
01/30/2020

Added four CRM sequences (CRM1-4) contributed by Na Wang from [Gernot and Presting (2008)](https://link.springer.com/article/10.1007%2Fs00438-007-0302-5).

02/03/2020

Added 3-letter names before all consensus seq IDs. IDs like "ZM00034_consensus" were inherited from the 2014 version MTEC, so kept it unchanged.

`perl -nle 's/>(.*)#(.*)\/(.*)/>$3_$1#$2\/$3/; print $_' maizeTE01302020 > maizeTE02032020`

02/05/2020

Fix namings

`perl -nle 's/>(.*)#(.*)\/(.*)/>$3_$1#$2\/$3/ if /consensus/; s/>(.*)/>$1#LTR\/CRM/ if /CRM/; print $_' history/maizeTE01302020 > maizeTE02052020`


