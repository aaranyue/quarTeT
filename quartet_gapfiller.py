#!/usr/bin/env python3

import sys
import argparse
import subprocess
import os
import re
import quartet_util

### MAIN PROGRAM ###
def GapFiller(args):
    draftgenomefile, gapclosercontigfilelist, flanking, minalignmentlength2, minalignmentidentity2, maxfillinglen, prefix, threads, minimapoption, overwrite, fillonly, joinonly = args

    # get gap's flanking seq
    print('[Info] Getting gaps flanking sequence...')
    draftgenomedict = quartet_util.readFastaAsDict(draftgenomefile)
    flankingdict = {}
    gapdict = {}
    for sid, seq in draftgenomedict.items():
        if 'N' in seq:
            i = 1
            gapsitelist = [r.span() for r in re.finditer(r'N+', seq)]
            for gapsite in gapsitelist:
                start = max(gapsite[0] - flanking, 0)
                end = min(gapsite[1] + flanking, len(seq))
                flankingdict[f'{sid}.{i}.L'] = seq[start:gapsite[0]]
                flankingdict[f'{sid}.{i}.R'] = seq[gapsite[1]:end]
                gapdict[f'{sid}.{i}'] = seq[gapsite[0]:gapsite[1]]
                i += 1
    if flankingdict == {}:
        print('[Error] Input genome does not have gap.')
    subprocess.run(f'mkdir tmp', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    flankingfastafile = f'tmp/{prefix}.gap.flanking.fasta'
    with open(flankingfastafile, 'w') as f:
        for sid, seq in flankingdict.items():
            if 'N' in seq:
                print('[Error] Flanking sequence contains gap. Recommend to lower -f parameter or check your file.')
                sys.exit(0)
            f.write(f'>{sid}\n{seq}\n')
    
    # process gapfilling file(s)
    gapcloserdict = {}
    for gapfillfile in gapclosercontigfilelist:
        gapfiller = os.path.basename(gapfillfile)
        print(f'[Info] gapfilling with {gapfiller}...')
        gapfillfasta = quartet_util.readFastaAsDict(gapfillfile)
        
        # If gapfiller itself has gap, split it
        gapfilldict = {}
        needsplit = False
        for sid, seq in gapfillfasta.items():
            if 'N' in seq:
                needsplit = True
                i = 1
                for tig in re.split(r'N+', seq):
                    gapfilldict[f'{sid}_tig{i}'] = tig
                    i += 1
            else:
                gapfilldict[sid] = seq
        if needsplit:
            print(f'[Info] Gaps found in {gapfiller}. Split it to contigs.')
            gapfillfile = f'tmp/{gapfiller}.splitcontig.fasta'
            with open(gapfillfile, 'w') as c:
                for tigid, seq in gapfilldict.items():
                    c.write(f'>{tigid}\n{seq}\n')
            gapfillfasta = gapfilldict

        pafgapfillfile = quartet_util.minimap(gapfillfile, flankingfastafile, prefix, f'flank_map_{gapfiller}', minimapoption, False, overwrite)

        allalignment = []
        with open(pafgapfillfile, 'r') as paf:
            for line in paf:
                if len(line.split()) < 11:
                    continue
                qryid, qrylen, qrystart, qryend, strand, refid, reflen, refstart, refend, match, alignlen = line.split()[:11]
                if int(alignlen) >= minalignmentlength2 and int(match) / int(alignlen) >= minalignmentidentity2:
                    allalignment.append({'qryid': qryid, 'qrylen': int(qrylen), 'qrystart': int(qrystart) + 1, 'qryend': int(qryend), 'strand': strand, 
                                         'refid': refid, 'reflen': int(reflen), 'refstart': int(refstart) + 1, 'refend': int(refend), 
                                         'match': int(match), 'alignlen': int(alignlen), 'identity': int(match) / int(alignlen),
                                         'gapid': '.'.join(qryid.split('.')[:-1]), 'LR': qryid.split('.')[-1]})
        if allalignment == []:
            print(f'[Info] No alignment for this file.')
            continue
        
        # process each gap
        print(f'[Info] Analysising Alignments...')
        for gapid in gapdict:
            Leftanchor = [aln for aln in allalignment if aln['gapid'] == gapid and aln['LR'] == 'L']
            Rightanchor = [aln for aln in allalignment if aln['gapid'] == gapid and aln['LR'] == 'R']
            if Leftanchor == [] or Rightanchor == []:
                continue
            # score each pair
            anchorpair = [(Laln, Raln) for Laln in Leftanchor for Raln in Rightanchor if Laln['refid'] == Raln['refid'] and Laln['strand'] == Raln['strand']]
            if anchorpair == []:
                continue
            bestscore = gapcloserdict[gapid]['score'] if gapid in gapcloserdict else 0
            for Laln, Raln in anchorpair:
                if Laln['strand'] == '-':
                    Laln, Raln = Raln, Laln
                if Laln['refend'] < Raln['refstart'] and joinonly != True:
                    score = Laln['identity'] + Raln['identity']
                    if score > bestscore:
                        fillstart = Laln['refend'] + Laln['qrylen'] - Laln['qryend'] + 1
                        fillend = Raln['refstart'] - Raln['qrystart']
                        fillseq = gapfillfasta[Laln['refid']][fillstart-1:fillend]
                        if len(fillseq) > maxfillinglen or fillseq == '':
                            continue
                        gapcloserdict[gapid] = {'sid': f'{gapfillfile}@{Laln["refid"]}', 'range': f'{fillstart}-{fillend}',
                                                'seq': fillseq if Laln['strand'] == '+' else quartet_util.reversedseq(fillseq), 'strand': Laln['strand'], 
                                                'score': score}
                elif Laln['refend'] >= Raln['refstart'] and fillonly != True:
                    if bestscore != 0:
                        continue
                    else:
                        gapcloserdict[gapid] = {'sid': f'{gapfillfile}@{Laln["refid"]}', 'range': 'join',
                                                'seq': '', 'strand': '', 
                                                'score': -1}
    
    print(f'[Info] All files processed.')
    if gapcloserdict == {}:
        print(f'[Error] No gap can be closed.')
        sys.exit(0)

    # make detail
    print(f'[Info] Generating filled fasta file and statistics...')
    filldetailfile = f'{prefix}.genome.filled.detail'
    with open(filldetailfile, 'w') as de:
        totalfilledlen = sum([len(gapcloserdict[x]['seq']) for x in gapcloserdict])
        de.write(f'# Gap Closed: {len(gapcloserdict)}\n# Total Filled length: {totalfilledlen}\n# Gap Remains: {len(gapdict)-len(gapcloserdict)}\n')
        de.write('# Seqid\tGap_identifier\tStatus\tCloserTigID\tCloserRange\tCloserLength\tCloserStrand\tCloserScore\n')
        for gapid in gapdict:
            if gapid in gapcloserdict:
                info = gapcloserdict[gapid]
                de.write(f'{".".join(gapid.split(".")[:-1])}\t{gapid.split(".")[-1]}\tClosed\t{info["sid"]}\t{info["range"]}\t{len(info["seq"])}\t{info["strand"]}\t{info["score"]}\n')
            else:
                de.write(f'{".".join(gapid.split(".")[:-1])}\t{gapid.split(".")[-1]}\tNot_closed\n')
    print(f'[Output] Filling detail write to: {filldetailfile}')
    
    # make fasta
    filledfastafile = f'{prefix}.genome.filled.fasta'
    with open(filledfastafile, 'w') as w:
        chrfastadict = {}
        for sid, seq in draftgenomedict.items():
            if sid not in ['.'.join(gapid.split('.')[:-1]) for gapid in gapcloserdict]:
                w.write(f'>{sid}\n{seq}\n')
                chrfastadict[sid] = seq
            else:
                seqlist = re.split(r'N+', seq)
                for i in range(len(seqlist) - 1):
                    gapid = f'{sid}.{i+1}'
                    if gapid in gapcloserdict:
                        seqlist.insert(2*i+1, gapcloserdict[gapid]['seq'])
                    else:
                        seqlist.insert(2*i+1, gapdict[gapid])
                seq = ''.join(seqlist)
                w.write(f'>{sid}\n{seq}\n')
                chrfastadict[f'{sid}'] = seq
    print(f'[Output] Filled genome fasta file write to: {filledfastafile}')
                
    # make stat
    filledstatfile = f'{prefix}.genome.filled.stat'
    with open(filledstatfile, 'w') as info:
        totallen = 0
        gc = 0
        for chrs in chrfastadict:
            totallen += len(chrfastadict[chrs])
            gc += chrfastadict[chrs].count('G') + chrfastadict[chrs].count('C')
        gccontent = gc / totallen
        info.write(f'# Total Size: {totallen}\n')
        info.write(f'# GC content: {gccontent}\n')
        info.write(f'# ChrID\tLength\tGapcount\tGaplocus\n')
        for sid, seq in chrfastadict.items():
            if 'N' in seq:
                count = 0
                locus = []
                gapsitelist = [r.span() for r in re.finditer(r'N+', seq)]
                for gapsite in gapsitelist:
                    start = gapsite[0] + 1
                    end = gapsite[1]
                    count += 1
                    locus.append(f'{start}-{end}')
                locus = "\t".join(locus)
                info.write(f'{sid}\t{len(seq)}\t{count}\t{locus}\n')
            else:
                info.write(f'{sid}\t{len(seq)}\t0\n')
    print(f'[Output] Filled genome stat write to: {filledstatfile}')
    
    # make plot
    agpfile = f'tmp/{prefix}.genome.agp'
    quartet_util.agpGap(filledfastafile, agpfile)
    quartet_util.drawgenome(agpfile, f'{prefix}.genome.filled')

### RUN ###
if __name__ == '__main__':
    # Argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', dest='draft_genome',required=True, help='(*Required) Draft genome file to be filled, FASTA format.')
    parser.add_argument('-g', dest='gapcloser_contig', nargs='+', required=True, help='(*Required) All contigs files (accept multiple file) used to fill gaps, FASTA format.')
    parser.add_argument('-f', dest='flanking_len', type=int, default=5000, help='The flanking seq length of gap used to anchor (bp), default: 5000')
    parser.add_argument('-l', dest='min_alignment_length', type=int, default=1000, help='The min alignment length to be select (bp), default: 1000')
    parser.add_argument('-i', dest='min_alignment_identity', type=float, default=40, help='The min alignment identity to be select (%%), default: 40')
    parser.add_argument('-m', dest='max_filling_len', type=int, default=1000000, help='The max sequence length acceptable to fill any gaps, default: 1000000')
    parser.add_argument('-p', dest='prefix', default='quarTeT', help='The prefix used on generated files, default: quarTeT')
    parser.add_argument('-t', dest='threads', default='1', help='Use number of threads, default: 1')
    parser.add_argument('--fillonly', dest='fillonly', action='store_true', default=False, help='Only fill the gaps without join.')
    parser.add_argument('--joinonly', dest='joinonly', action='store_true', default=False, help='Only join the gaps without fill.')
    parser.add_argument('--overwrite', dest='overwrite', action='store_true', default=False, help='Overwrite existing alignment file instead of reuse.')
    parser.add_argument('--minimapoption', dest='minimapoption', default='-x asm5', help='Pass additional parameters to minimap2 program, default: -x asm5')

    # parse input paramater
    draftgenomefile = quartet_util.decompress(parser.parse_args().draft_genome)
    gapclosercontigfilelist = [quartet_util.decompress(x) for x in parser.parse_args().gapcloser_contig]
    flanking = int(parser.parse_args().flanking_len)
    minalignmentlength2 = int(parser.parse_args().min_alignment_length)
    if minalignmentlength2 > flanking:
        print('[Error] min_alignment_length should be less than flanking_len.')
    minalignmentidentity2 = float(parser.parse_args().min_alignment_identity) / 100
    if minalignmentidentity2 < 0 or minalignmentidentity2 > 1:
        print('[Error] min_alignment_identity should be within 0~100.')
    maxfillinglen = int(parser.parse_args().max_filling_len)
    prefix = parser.parse_args().prefix
    threads = parser.parse_args().threads
    minimapoption = parser.parse_args().minimapoption + f' -t {threads}'
    overwrite = parser.parse_args().overwrite
    fillonly = parser.parse_args().fillonly
    joinonly = parser.parse_args().joinonly

    # check prerequisites
    quartet_util.check_prerequisite(['minimap2', 'Rscript'])

    # run
    args = [draftgenomefile, gapclosercontigfilelist, flanking, minalignmentlength2, minalignmentidentity2, 
            maxfillinglen, prefix, threads, minimapoption, overwrite, fillonly, joinonly]
    quartet_util.run(GapFiller, args)
