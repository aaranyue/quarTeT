#!/usr/bin/env python3

import subprocess
import os
import sys
from multiprocessing.dummy import Pool
import argparse
import itertools
import quartet_util

def centroMiner(args):
    genomefile, tegfffile, genegfffile, minperiod, maxperiod, e, maxgap, minlength, prefix, threads, overwrite, noplot, match, mismatch, delta, PctMatch, PctIndel, minscore, identity, periodmaxdelta, wordlength, max_TR_length = args
    
    # split genome into chr
    subprocess.run(f'mkdir tmp', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    subprocess.run(f'mkdir tmp/splitchr', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    subprocess.run(f'mkdir tmp/trfdat', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    subprocess.run(f'mkdir tmp/blast', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    subprocess.run(f'mkdir tmp/Rdata', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    subprocess.run(f'mkdir TandemRepeat', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    subprocess.run(f'mkdir Candidates', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

    genomedict = quartet_util.readFastaAsDict(genomefile)
    chrlendict = {}
    print('[Info] Spliting genome to chromosome...')
    for Chr in genomedict:
        chrlendict[Chr] = len(genomedict[Chr])
        splitchrfastafile = f'tmp/splitchr/{prefix}.{Chr}.fasta'
        with open(splitchrfastafile, 'w') as Chrfasta:
            Chrfasta.write(f'>{Chr}\n{genomedict[Chr]}\n')
    genomedictkey = genomedict.keys()
    del genomedict

    def centroblaster(Chr):
        # run TRF and get all patterns
        datfile = f'tmp/trfdat/{prefix}.{Chr}.fasta.{match}.{mismatch}.{delta}.{PctMatch}.{PctIndel}.{minscore}.{maxperiod}.dat'
        splitchrfastafile = f'tmp/splitchr/{prefix}.{Chr}.fasta'
        if not os.path.exists(datfile) or overwrite == True:
            subprocess.run(f'trf {splitchrfastafile} {match} {mismatch} {delta} {PctMatch} {PctIndel} {minscore} {maxperiod} -l {max_TR_length} -d -h', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        subprocess.run(f'mv -t tmp/trfdat -f {prefix}.{Chr}.fasta.{match}.{mismatch}.{delta}.{PctMatch}.{PctIndel}.{minscore}.{maxperiod}.dat', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        with open(datfile, 'r') as trfResult:
            linelist = []
            for line in trfResult:
                if line.startswith('Sequence:'):
                    seqid = line.split(' ')[1].strip()
                if len(line.split(' ')) == 15:
                    start, end, period, copynum, consensusSize, pctmatch, pctindel, score, perc_A, perc_C, perc_G, perc_T, entropy, repeat_pattern, repeat_seq = line.split(' ')
                    if int(period) > int(minperiod):
                        linelist.append([seqid, start, end, period, copynum, consensusSize, pctmatch, pctindel, score, perc_A, perc_C, perc_G, perc_T, entropy, repeat_pattern, repeat_seq])
        trfastafile = f'TandemRepeat/{prefix}.{Chr}.tr.fasta'
        with open(trfastafile, 'w') as pattern:
                count = 1
                for [seqid, start, end, period, copynum, consensusSize, pctmatch, pctindel, score, perc_A, perc_C, perc_G, perc_T, entropy, repeat_pattern, repeat_seq] in linelist:
                    trfid = f'{Chr}@TR_{str(count).zfill(5)}'
                    pattern.write(f'>{trfid}\n{repeat_pattern}\n')
                    count += 1

        # blast chr with patterns
        clusterfastafile = f'TandemRepeat/{prefix}.{Chr}.tr.cluster.fasta'
        subprocess.run(f'cd-hit-est -i {trfastafile} -o {clusterfastafile} -c {identity} -n {wordlength} -S {periodmaxdelta} -M 0', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        blastdb = f'tmp/splitchr/{prefix}.{Chr}'
        subprocess.run(f'makeblastdb -dbtype nucl -in {splitchrfastafile} -out {blastdb}', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        blastresultfile = f'tmp/blast/{prefix}.{Chr}.tr.blast'
        subprocess.run(f'blastn -db {blastdb} -query {clusterfastafile} -out {blastresultfile} -outfmt 7 -evalue {e}', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

        # build TR annotation
        blastgff3file = f'TandemRepeat/{prefix}.{Chr}.tr.blast.gff3'
        with open(blastresultfile, 'r') as c, open(blastgff3file, 'w') as w:
            trintervals = []
            for line in c:
                if line.startswith('#'):
                    continue
                query, subject, iden, length, mismatches, gap, qstart, qend, sstart, send, evalue, bitscore = line.split('\t')
                ref = subject
                source = 'BLAST'
                typo = 'TR'
                if send > sstart:
                    start = sstart
                    end = send
                    strand = '+'
                else:
                    start = send
                    end = sstart
                    strand = '-'
                score = bitscore.strip()
                phase = '.'
                attributes = f'ID={query};length={length};identity={iden};evalue={evalue};'
                w.write(f'{ref}\t{source}\t{typo}\t{start}\t{end}\t{score}\t{strand}\t{phase}\t{attributes}\n')
                trintervals.append([int(start), int(end)])

        # read TE and gene annotation
        if tegfffile != None:
            with open(tegfffile, 'r') as te:
                teintervals = []
                for line in te:
                    if line.startswith('#') or len(line.split()) < 9:
                        continue
                    if line.split()[0] == Chr and ('LTR' in line.split()[2].upper() or 'LTR' in line.split()[8].upper() or 'long_terminal_repeat' in line.split()[2]):
                        teintervals.append([int(line.split()[3]), int(line.split()[4])])
                        
        if genegfffile != None:
            with open(genegfffile, 'r') as ge:
                geintervals = []
                for line in ge:
                    if line.startswith('#') or len(line.split()) < 9:
                        continue
                    if line.split()[0] == Chr and 'gene' in line.split()[2]:
                        geintervals.append([int(line.split()[3]), int(line.split()[4])])
                
        
        # analysis blast result and build region info
        with open(blastresultfile, 'r') as blast:
            linelist = []
            for line in blast:
                if line.startswith('#'):
                    continue
                linelist.append(line.split('\t'))

        linelist.sort(key=lambda x:(min(int(x[8]), int(x[9]))))
        continuousregion = []
        prevstart = 1
        prevend = 1
        TRrange = []
        subTRrange = {}
        subTRlength = {}
        for [query, subject, ident, alignmentlength, mismatches, gapopens, querystart, queryend, subjectstart, subjectend, evalue, bitscore] in linelist:
            thisstart = min(int(subjectstart), int(subjectend))
            if thisstart > prevend + int(maxgap):
                # collect region and restart
                prevlength = float(prevend - prevstart)
                if prevlength > minlength:
                    TRlength = len(set(itertools.chain(*[list(x) for x in TRrange])))
                    for tr, ran in subTRrange.items():
                        subTRlength[tr] = len(set(itertools.chain(*[list(x) for x in ran])))
                    continuousregion.append([prevstart, prevend, TRlength, subTRlength])
                prevstart = thisstart
                TRrange = []
                subTRrange = {}
                subTRlength = {}
            prevend = max(int(subjectstart), int(subjectend))
            TRrange.append(range(min(int(subjectstart), int(subjectend)), max(int(subjectstart), int(subjectend)) + 1))
            if query not in subTRrange:
                subTRrange[query] = []
            subTRrange[query].append(range(min(int(subjectstart), int(subjectend)), max(int(subjectstart), int(subjectend)) + 1))
        # collect final region
        prevlength = float(prevend - prevstart)
        if prevlength > minlength:
            TRlength = len(set(itertools.chain(*[list(x) for x in TRrange])))
            for tr, ran in subTRrange.items():
                subTRlength[tr] = len(set(itertools.chain(*[list(x) for x in ran])))
            continuousregion.append([prevstart, prevend, TRlength, subTRlength])

        # sort and write candidates
        TRdict = quartet_util.readFastaAsDict(trfastafile)
        continuousregion.sort(key=lambda x:x[2]/(x[1]-x[0]), reverse=True)
        candidatefile = f'Candidates/{prefix}.{Chr}.candidate'
        with open(candidatefile, 'w') as out:
            out.write(f'# Chr\tstart\tend\tlength\tTRlength\tTRcoverage\n')
            out.write(f'#\tsubTR\tperiod\tsubTRlength\tsubTRcoverage\tpattern\n')
            for [start, end, TRlength, subTRlength] in continuousregion:
                length = end-start+1
                TRcoverage = round(TRlength/length * 100, 2)
                out.write(f'{Chr}\t{start}\t{end}\t{length}\t{TRlength}\t{TRcoverage}%\n')
                for [tr, subtrlength] in sorted(subTRlength.items(), key=lambda x:x[1], reverse=True):
                    subTRcoverage = round(subtrlength/length * 100, 2)
                    out.write(f'\t{tr}\t{len(TRdict[tr])}\t{subtrlength}\t{subTRcoverage}%\t{TRdict[tr]}\n')
        del TRdict
        
        # draw TR, TE, and gene line chart
        if noplot != True:
            tsvfile = f'tmp/Rdata/{prefix}.{Chr}.tsv'
            pdffile = f'Candidates/{prefix}.{Chr}.pdf'
            with open(tsvfile, 'w') as tsv:
                for i in range(0, chrlendict[Chr], 50000):
                    l = i - 50000
                    r = i + 50000
                    v1 = quartet_util.calculate_cover_length(trintervals, l, r)
                    tsv.write(f'{i}\t{v1}\tTR\n')
                    if tegfffile != None:
                        v2 = quartet_util.calculate_cover_length(teintervals, l, r)
                        tsv.write(f'{i}\t{v2}\tTE\n')
                    if genegfffile != None:
                        v3 = quartet_util.calculate_cover_length(geintervals, l, r)
                        tsv.write(f'{i}\t{v3}\tgene\n')
            rscript = f'library(ggplot2);options(scipen=999);data<-read.table("{tsvfile}");colnames(data)<-c("site","value","type");data$site<-as.numeric(data$site);data$value<-as.numeric(data$value);pdf("{pdffile}");ggplot(data,aes(x=site,y=value,group=type,color=type,shape=type))+geom_line()+labs(x="Position",y="Length (bp)")+theme(axis.text.x=element_text(angle=90,hjust=0.5))+scale_x_continuous(breaks=seq(0,max(data$site),1000000),minor_breaks=seq(0,max(data$site),500000))+facet_wrap(~type,scales="free_y",dir="v")'
            subprocess.run(f"echo '{rscript}' | Rscript -", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

        print(f'[Info] {Chr} done.')
        
    # multithread
    def print_error(value):
        print("error: ", value)
    print('[Info] Processing each chromosome...')
    p = Pool(min(len(genomedictkey), int(threads)))
    for Chr in genomedictkey:
        p.apply_async(centroblaster, (Chr,), error_callback=print_error)
    p.close()
    p.join()
    
    print(f'[Output] Tandem repeats data write to folder: TandemRepeat')
    print(f'[Output] Centromere candidates data write to folder: Candidates')

### RUN ###
if __name__ == '__main__':
    # Argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='genome_fasta',required=True, help='(*Required) Genome file, FASTA format.')
    parser.add_argument('--TE', dest='TE', default=None, help='TE annotation file, gff3 format.')
    parser.add_argument('--gene', dest='gene', default=None, help='gene annotation file, gff3 format.')
    parser.add_argument('-n', dest='min_period', type=int, default=100, help='Min period to be consider as centromere repeat monomer. Default: 100')
    parser.add_argument('-m', dest='max_period', type=int, default=200, help='Max period to be consider as centromere repeat monomer. Default: 200')
    parser.add_argument('-s', dest='cluster_identity', type=float, default=0.8, help='Min identity between TR monomers to be clustered (Cannot be smaller than 0.8). Default: 0.8')
    parser.add_argument('-d', dest='cluster_max_delta', type=float, default=10, help='Max period delta for TR monomers in a cluster. Default: 10')
    parser.add_argument('-e', dest='evalue', type=float, default=0.00001, help='E-value threholds in blast. Default: 0.00001')
    parser.add_argument('-g', dest='max_gap', type=int, default=50000, help='Max allowed gap size between two tandem repeats to be considered as in one tandem repeat region. Default: 50000')
    parser.add_argument('-l', dest='min_length', type=int, default=100000, help='Min size of tandem repeat region to be selected as candidate. Default: 100000')
    parser.add_argument('-t', dest='threads', default='1', help='Limit number of using threads, default: 1')
    parser.add_argument('-p', dest='prefix', type=str, default='quarTeT', help='Prefix used by generated files. Default: quarTeT')
    parser.add_argument('--trf', dest='trf_parameter', nargs='*', default=[2,7,7,80,10,50], help='Change TRF parameters: <match> <mismatch> <delta> <PM> <PI> <minscore> Default: 2 7 7 80 10 50')
    parser.add_argument('-r', dest='max_TR_length', type=int, default=3, help='Maximum TR length (in millions) expected for trf. Default: 3')
    parser.add_argument('--overwrite', dest='overwrite', action='store_true', default=False, help='Overwrite existing trf dat file instead of reuse.')
    parser.add_argument('--noplot', dest='noplot', action='store_true', default=False, help='Skip all ploting.')

    # parse input paramater
    genomefile = quartet_util.decompress(parser.parse_args().genome_fasta)

    tegfffile = quartet_util.decompress(parser.parse_args().TE)
    genegfffile = quartet_util.decompress(parser.parse_args().gene)

    minperiod = int(parser.parse_args().min_period)
    maxperiod = int(parser.parse_args().max_period)
    e = float(parser.parse_args().evalue)
    maxgap = int(parser.parse_args().max_gap)
    minlength = int(parser.parse_args().min_length)
    max_TR_length = int(parser.parse_args().max_TR_length)
    
    prefix = parser.parse_args().prefix
    threads = parser.parse_args().threads
    overwrite = parser.parse_args().overwrite
    noplot = parser.parse_args().noplot

    if len(parser.parse_args().trf_parameter) != 6:
        print('[Error] TRF parameter should be <match> <mismatch> <delta> <PM> <PI> <minscore>.')
        sys.exit(0)
    match = int(parser.parse_args().trf_parameter[0])
    mismatch = int(parser.parse_args().trf_parameter[1])
    delta = int(parser.parse_args().trf_parameter[2])
    PctMatch = int(parser.parse_args().trf_parameter[3])
    PctIndel = int(parser.parse_args().trf_parameter[4])
    minscore = int(parser.parse_args().trf_parameter[5])

    identity = float(parser.parse_args().cluster_identity)
    periodmaxdelta = int(parser.parse_args().cluster_max_delta)
    if identity < 0.8 or identity > 1:
        print('[Error] Cluster identity should be set in 0.8 ~ 1.')
        sys.exit(0)
    elif identity < 0.85:
        wordlength = 5
    elif identity < 0.88:
        wordlength = 6
    elif identity < 0.9:
        wordlength = 8
    else:
        wordlength = 10

    # check prerequisites
    quartet_util.check_prerequisite(['trf', 'cd-hit-est', 'makeblastdb', 'blastn'])

    # run
    args = [genomefile, tegfffile, genegfffile, minperiod, maxperiod, e, maxgap, minlength, prefix, threads, overwrite, noplot,
            match, mismatch, delta, PctMatch, PctIndel, minscore, identity, periodmaxdelta, wordlength, max_TR_length]
    print(f'[Info] Paramater: genomefile={genomefile}, tegfffile={tegfffile}, genegfffile={genegfffile}, minperiod={minperiod}, maxperiod={maxperiod}, e={e}, maxgap={maxgap}, minlength={minlength}, prefix={prefix}, threads={threads}, overwrite={overwrite}, noplot={noplot}, match={match}, mismatch={mismatch}, delta={delta}, PctMatch={PctMatch}, PctIndel={PctIndel}, minscore={minscore}, identity={identity}, periodmaxdelta={periodmaxdelta}, wordlength={wordlength}, max_TR_length={max_TR_length}')
    quartet_util.run(centroMiner, args)
