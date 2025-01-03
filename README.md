# quarTeT: Telomere-to-telomere Toolkit
quarTeT is a collection of tools for T2T genome assembly and basic analysis in automatic workflow.

Task include:
- [AssemblyMapper](#AssemblyMapper): reference-guided genome assembly
- [GapFiller](#GapFiller): long-reads based gap filling
- [TeloExplorer](#TeloExplorer): telomere identification
- [CentroMiner](#CentroMiner): centromere candidate prediction

## Version Change log
1.2.5
- Add new '--groupcontig' option for AssemblyMapper. Adding this option will output a folder containing contigs grouped by reference sequence (will group unassigned contigs into one).

1.2.4
- Add new '--teclade' and '--teminrepeattimes' option for AssemblyMapper to control the behavior of built-in TeloExplorer.
- Add new '-a' option for GapFiller to select unimap as aligner. (Also fix the bug that default aligner not set in GapFiller after v1.2.3, thanks to a927050047, [PR #46](https://github.com/aaranyue/quarTeT/pull/46)) 

1.2.3
- Add new '--extract-ref-flanks' option for AssemblyMapper, which allow chimeric contig output for gap filling. (see [issue #42](https://github.com/aaranyue/quarTeT/issues/42) for detail)
- Support Unimap as aligner. As a optimized version of minimap2, it still use the '--minimapoption'.
- Fix a bug that monopolizer contig length standard is too high.
- Deprecate the web server for maintenance issue (sorry!).

1.2.2
- Add new '--keep' option for AssemblyMapper, which add all unplaced contigs in the draft genome.
- Add a new output AGP file for GapFiller to describe the modified chromosome structure.
- Fix another bug that cause small number of N representing unknown bases are identified as gap.

1.2.1
- Fix a bug in CentroMiner that ploting halted when optional gene/TE annotation file is not given. 

1.2.0
- CentroMiner is refactored. It can receive gene annotation as well now. The output has only 2 folders now: TandemRepeat (fasta, gff3 of TR) and Candidate (candidate info and TR, TE, and gene content line chart, require ggplot2). TE data is removed from candidate info. It is recommended to check the line chart to decide which candidate you vote.
- '--noplot' option is added to each module. With this option, any ploting will be skiped. If you have problem in graphical issue, try this option.
- CloserScore in GapFiller detail is changed to CloserIdentity.
- Due to the major changes, the web server version of quarTeT will not be updated to the newest version for now.

1.1.8
- Gapfiller will throw a warning instead of error when flanking sequence contains gap.

1.1.7
- Support RepeatMasker's TE annotation format for CentroMiner module.

1.1.6
- add new option 'maximum TR length' (-r) for CentroMiner to avoid trf stuck. (Thanks to atotickov, [PR #22](https://github.com/aaranyue/quarTeT/pull/22))
- disable unstable 'join' mode for GapFiller by default. Use '--enablejoin' option to enable this mode. '--fillonly' option is removed.

1.1.5
- SVG output is moved to work dir instead of tmp dir. Intermediate file for figure drawing is saved to tmp dir instead of auto-remove.
- Fix a bug that running multiple quarTeT in one folder may cause error due to intermediate file overwrite.
- Fix a bug in AssemblyMapper that with option '--nofilter', contigs shorter than 50000 bp are still marked as too short and count in discarded length.
- Fix a bug that error in R figure drawing is not reported.

1.1.4
- Fix a bug in AssemblyMapper that large dict tmp file not write properly.
- Reduce more peak memory.
- Add a memory insufficient error report.

1.1.3
- Reduce peak memory.
- Add more error report.
- Fix some error without exit.

1.1.2
- Fix a bug that AssemblyMapper cannot overwrite existing telomere checking result.
- Fix a bug that small number of N repesenting unknown bases are identified as gap.

1.1.1
- Fix a bug that CentroMiner stuck after v1.0.4

1.1.0 
- AssemblyMapper: new option '--nofilter'. With this option, input contigs will not be split or discard even if have gaps or too short.
- GapFiller: support join, but this is not as reliable as fill. you can use option '--fillonly' and '--joinonly' to disable one of them.
- TeloExplorer: now compatible with latest tidk version 0.2.31.
- fix a bug that error report added in v1.0.4 didn't include stderr.

1.0.4 
- Add more report when called programs are failed.

1.0.3
- Fix a bug that when figure drawing is failed, there are no warning raised.

1.0.2 
- Fix a bug in TeloExplorer that when more than one possible telomere-like repeats are found, it will be considered as no telomere-like repeat found.

1.0.1 
- Fix a bug in CentroMiner that when no centromere-like region is found on a chromosome, genome overview plotting will unexceptly exit.

1.0.0 
- Initial release

## Getting Started
### Use quarTeT on Web
~~quarTeT can be easily accessed on [our web server](http://www.atcgn.com:8080/quarTeT/home.html).~~ (Currently deprecated.)
### Use quarTeT on local
quarTeT command-line program is availble for Linux.
#### Dependencies
- [Python3](https://www.python.org/) (>3.6, tested on 3.7.4 and 3.9.12) 
- [Minimap2](https://github.com/lh3/minimap2) (tested on 2.24-r1122 and 2.24-r1155-dirty) 
- (Optional)[Unimap](https://github.com/lh3/unimap) (tested on 0.1-r41)
- [MUMmer4](https://github.com/mummer4/mummer) (tested on 4.0.0rc1) 
- [trf](https://github.com/Benson-Genomics-Lab/TRF) (tested on 4.09) 
- [CD-hit](https://github.com/weizhongli/cdhit) (tested on 4.6 and 4.8.1) 
- [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (tested on 2.8.1 and 2.11.0) 
- [tidk](https://github.com/tolkit/telomeric-identifier) (tested on 0.2.1 and 0.2.31) 
- [gnuplot](https://github.com/gnuplot/gnuplot) (tested on 4.6 patchlevel 2 and 6) 
- [R](https://www.R-project.org/) (>3.5.0, tested on 3.6.0 and 4.2.2) 
    - RIdeogram (tested on 0.2.2)
    - ggplot2 (tested on 3.3.6 and 3.4.4)

All these dependencies can be easily install via conda:

`conda create -n quartet --channel conda-forge --channel bioconda Python Minimap2 MUMmer4 trf CD-hit BLAST tidk R R-RIdeogram R-ggplot2 gnuplot`

(Recently we discover that using conda to install R will result in blank PNG. However, SVG is correctly generated.)

#### Installation
quarTeT do not require installation.

Just clone this repository, and run `python3 {path}/quartet.py`

## Usage
```
quarTeT: Telomere-to-telomere Toolkit

Usage: python3 quartet.py <module> <parameters>

Modules:
AssemblyMapper  | am    Assemble draft genome.
GapFiller       | gf    Fill gaps in draft genome.
TeloExplorer    | te    Identify telomeres.
CentroMiner     | cm    Identify centromere candidates.

Use <module> -h for module usage.
```

### AssemblyMapper
AssemblyMapper is a reference-guided assemble tool.

A phased contig-level assembly and a close-related reference genome are required as input, both in fasta format.

Note that contigs should be phased. 

It's recommended to obtain such an assembly using [hifiasm](https://github.com/chhylp123/hifiasm).

you can convert `{prefix}.bp.hap1.p_ctg.gfa` and `{prefix}.bp.hap2.p_ctg.gfa` generated by hifiasm to FASTA format as input, respectively.  

```
Usage: python3 quartet.py AssemblyMapper <parameters>
  -h, --help            show this help message and exit
  -r REFERENCE_GENOME   (*Required) Reference genome file, FASTA format.
  -q CONTIGS            (*Required) Phased contigs file, FASTA format.
  -c MIN_CONTIG_LENGTH  Contigs shorter than INT (bp) will be removed, default: 50000
  -l MIN_ALIGNMENT_LENGTH
                        The min alignment length to be select (bp), default: 10000
  -i MIN_ALIGNMENT_IDENTITY
                        The min alignment identity to be select (%), default: 90
  -p PREFIX             The prefix used on generated files, default: quarTeT
  -t THREADS            Use number of threads, default: 1
  -a {minimap2,unimap,mummer}
                        Specify alignment program (support minimap2, unimap and mummer), default: minimap2
  --nofilter            Use original sequence input, no filtering.
  --keep                Keep the unplaced contigs in draft genome
  --groupcontig         Add an folder output of contigs grouped by destination.
  --extract-ref-flanks CHIMERA
                        Add an output of chimera contig containing reference flanks of x bp (check issue#42 for detail), default: 0 (off)
  --plot                Plot a colinearity graph for draft genome to reference alignments. (will cost more time)
  --noplot              Skip all ploting.
  --overwrite           Overwrite existing alignment file instead of reuse.
  --minimapoption MINIMAPOPTION
                        Pass additional parameters to minimap2/unimap program, default: -x asm5
  --nucmeroption NUCMEROPTION
                        Pass additional parameters to nucmer program.
  --deltafilteroption DELTAFILTEROPTION
                        Pass additional parameters to delta-filter program.
  --teclade {plant,animal,other}
                        Specify clade of this genome for telomere search. Plant will search TTTAGGG, animal will search TTAGGG, other will use tidk explore's suggestion, default: other
  --teminrepeattimes TE_MIN_REPEAT_TIMES
                        The min repeat times to considered as telomere, default: 100
```
Output files should be as follow:
```
{prefix}.draftgenome.fasta        | The pseudo-chromosome-level assembly, fasta format.
{prefix}.draftgenome.agp          | The structure of this assembly, AGP format.
{prefix}.draftgenome.stat         | The statistic of this assembly, including total size and each chromosome's size, GC content, gap count and locations.
{prefix}.draftgenome.png          | The figure draws relative length of chromosomes and gap locations for assembly.
{prefix}.contig.mapinfo           | The statistic of input contigs, including total mapped and discarded size, and each contig's destination.
{prefix}.contig_map_ref.png       | The alignment colinearity graph between contigs and reference genome.
{prefix}.draftgenome_map_ref.png  | The alignment colinearity graph between this assembly genome and reference genome. Only available with --plot.
```

### GapFiller
GapFiller is a long-reads based gapfilling tool.

A gap-tied genome and corresponding long-reads are required as input, both in fasta format.

If possible, using long-reads assembled and polished contigs instead of reads may improve the quality.

GapFiller support two mode: `fill` and `join`. `fill` means find a segment that can be placed in the gap and link them together. `join` means find an evidence that the gap edge is overlaped and can be directly merged into one. For now, `join` mode is unstable and disabled by default. If you want to enable it, be careful and check the result manually.

```
Usage: python3 quartet.py GapFiller <parameters>
  -h, --help            show this help message and exit
  -d DRAFT_GENOME       (*Required) Draft genome file to be filled, FASTA format.
  -g GAPCLOSER_CONTIG [GAPCLOSER_CONTIG ...]
                        (*Required) All contigs files (accept multiple file) used to fill gaps, FASTA format.
  -f FLANKING_LEN       The flanking seq length of gap used to anchor (bp), default: 5000
  -l MIN_ALIGNMENT_LENGTH
                        The min alignment length to be select (bp), default: 1000
  -i MIN_ALIGNMENT_IDENTITY
                        The min alignment identity to be select (%), default: 40
  -m MAX_FILLING_LEN    The max sequence length acceptable to fill any gaps, default: 1000000
  -a {minimap2,unimap}  Specify alignment program (support minimap2 and unimap), default: minimap2
  -p PREFIX             The prefix used on generated files, default: quarTeT
  -t THREADS            Use number of threads, default: 1
  --enablejoin          Enable join mode to close the gaps. (Unstable)
  --joinonly            Use only join mode without fill, should be used with --enablejoin.
  --overwrite           Overwrite existing alignment file instead of reuse.
  --minimapoption MINIMAPOPTION
                        Pass additional parameters to minimap2 program, default: -x asm5
  --noplot              Skip all ploting.
```
Output files should be as follow:
```
{prefix}.genome.filled.fasta        | The gap-filled genome, fasta format.
{prefix}.genome.filled.modified.agp | The modified chromosome structure, AGP format.
{prefix}.genome.filled.detail       | Detailed information for each gap, including gap closed and remains, total filled size and closer's ID, range, etc.
{prefix}.genome.filled.stat         | The statistic of filled genome, including total size and each chromosome's size, GC content, gap count and locations.
{prefix}.genome.filled.png          | The figure draws relative length of chromosomes and gap locations for assembly.
```

### TeloExplorer
TeloExplorer is a telomere identification tool.

A genome file in fasta format is required as input.
```
Usage: python3 quartet.py TeloExplorer <parameters>
  -h, --help            show this help message and exit
  -i GENOME             (*Required) Genome file to be identified, FASTA format.
  -c {plant,animal,other}
                        Specify clade of this genome. Plant will search TTTAGGG, animal will search TTAGGG, other will use tidk explore's suggestion, default: other
  -m MIN_REPEAT_TIMES   The min repeat times to be reported, default: 100
  -p PREFIX             The prefix used on generated files, default: quarTeT
  --noplot              Skip all ploting.
```
Output files should be as follow:
```
{prefix}.telo.info  | The statistic of telomere, including monomer, repeat times on both end of each chromosome.
{prefix}.telo.png   | The figure draws telomere location, alongside relative length of chromosomes and gap locations for assembly.
```

### CentroMiner
CentroMiner is a centromere prediction tool.

A genome file in fasta format is required as input.

Optionally, an additional input of TE annotation (or just LTR annotation) and gene annotation in gff3 format can improve the performance.

It's recommended to obtain TE annotation using [EDTA](https://github.com/oushujun/EDTA).

`{genome file}.mod.EDTA.TEanno.gff3` generated by EDTA can directly feed CentroMiner, unless you have sequence ID longer than 15 characters.

Note that the sequence ID in first column should be consistent with in genome. Some tools may change sequence ID if ID is too long.

The sequence ontology in the third column should include `LTR` or `long_terminal_repeat` (EDTA default), or have `LTR` in the eighth column (RepeatMasker default) to be recognized.
```
Usage: python3 quartet.py CentroMiner <parameters>
  -h, --help            show this help message and exit
  -i GENOME_FASTA       (*Required) Genome file, FASTA format.
  --TE TE               TE annotation file, gff3 format.
  --gene GENE           gene annotation file, gff3 format.
  -n MIN_PERIOD         Min period to be consider as centromere repeat monomer. Default: 100
  -m MAX_PERIOD         Max period to be consider as centromere repeat monomer. Default: 200
  -s CLUSTER_IDENTITY   Min identity between TR monomers to be clustered (Cannot be smaller than 0.8). Default: 0.8
  -d CLUSTER_MAX_DELTA  Max period delta for TR monomers in a cluster. Default: 10
  -e EVALUE             E-value threholds in blast. Default: 0.00001
  -g MAX_GAP            Max allowed gap size between two tandem repeats to be considered as in one tandem repeat region. Default: 50000
  -l MIN_LENGTH         Min size of tandem repeat region to be selected as candidate. Default: 100000
  -t THREADS            Limit number of using threads, default: 1
  -p PREFIX             Prefix used by generated files. Default: quarTeT
  --trf [TRF_PARAMETER ...]
                        Change TRF parameters: <match> <mismatch> <delta> <PM> <PI> <minscore> Default: 2 7 7 80 10 50
  -r MAX_TR_LENGTH      Maximum TR length (in millions) expected for trf. Default: 3
  --overwrite           Overwrite existing trf dat file instead of reuse.
  --noplot              Skip all ploting.
```
Output files should be as follow:
```
Candidates/             | The folder of all centromere candidates data on each chromosome. Check the pdf line chart and vote a best candidate.
TandemRepeat/           | The folder of all tandem repeat monomers identified by trf and cluster result on each chromosome.
```

**Note that manual selection is highly recommended.**

In some species, best candidate may unexpectedly fall in near telomere region or other tandem repeat region.

You can combine the gene, TR and TE line chart with Hi-C contact heatmap, pairwise colinearity etc. to better vote your candidates. Mostly the centromere region is riched in TR, surrounded by rich TE, and low in gene content, but there are also exceptions. 

## Citation
Yunzhi Lin, Chen Ye, Xingzhu Li, Qinyao Chen, Ying Wu, Feng Zhang, Rui Pan, Sijia Zhang, Shuxia Chen, Xu Wang, Shuo Cao, Yingzhen Wang, Yi Yue, Yongsheng Liu, Junyang Yue. quarTeT: a telomere-to-telomere toolkit for gap-free genome assembly and centromeric repeat identification. Horticulture Research 2023;10:uhad127, https://doi.org/10.1093/hr/uhad127


