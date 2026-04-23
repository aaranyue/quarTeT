#!/usr/bin/env python3
# Last modified: V1.2.4

import sys
import argparse
import subprocess
import os
import re
import quartet_util

def evaluate_suffix_prefix_overlap(left_seq, right_seq, min_identity, min_overlap, max_check_len):
    left_check = left_seq[-max_check_len:] if len(left_seq) > max_check_len else left_seq
    right_check = right_seq[:max_check_len] if len(right_seq) > max_check_len else right_seq
    max_overlap = min(len(left_check), len(right_check))
    if max_overlap == 0:
        return {'status': 'no_overlap', 'length': 0, 'identity': 0}
    min_overlap = max(1, min(min_overlap, max_overlap))
    best_any = None
    best_long = None
    best_valid = None
    for overlap_len in range(max_overlap, 0, -1):
        left_overlap = left_check[-overlap_len:]
        right_overlap = right_check[:overlap_len]
        match_count = sum(1 for left_base, right_base in zip(left_overlap, right_overlap) if left_base == right_base)
        overlap_identity = match_count / overlap_len
        current = {'length': overlap_len, 'identity': overlap_identity, 'match_count': match_count}
        if best_any is None or current['match_count'] > best_any['match_count'] or (current['match_count'] == best_any['match_count'] and current['length'] > best_any['length']):
            best_any = current
        if overlap_len >= min_overlap:
            if best_long is None or current['identity'] > best_long['identity'] or (current['identity'] == best_long['identity'] and current['length'] > best_long['length']):
                best_long = current
            if overlap_identity >= min_identity:
                best_valid = current
                break
    if best_valid is not None:
        return {'status': 'pass', 'length': best_valid['length'], 'identity': best_valid['identity']}
    if best_any is None or best_any['match_count'] == 0:
        return {'status': 'no_overlap', 'length': 0, 'identity': 0}
    if best_long is None:
        return {'status': 'overlap_too_short', 'length': best_any['length'], 'identity': best_any['identity']}
    return {'status': 'overlap_identity_too_low', 'length': best_long['length'], 'identity': best_long['identity']}

def update_joinattempt(joinattemptdict, gapid, status, info):
    if gapid not in joinattemptdict:
        joinattemptdict[gapid] = {'status': status, **info}
        return
    if status == 'Joined':
        joinattemptdict[gapid] = {'status': status, **info}
        return
    if joinattemptdict[gapid]['status'] == 'Joined':
        return
    prevscore = joinattemptdict[gapid].get('score', -1)
    newscore = info.get('score', -1)
    if newscore >= prevscore:
        joinattemptdict[gapid] = {'status': status, **info}

### MAIN PROGRAM ###
def GapFiller(args):
    draftgenomefile, gapclosercontigfilelist, flanking, minalignmentlength2, minalignmentidentity2, maxfillinglen, prefix, threads, minimapoption, overwrite, enablejoin, joinonly, joinminoverlap, joinminidentity, noplot, aligner = args

    # get gap's flanking seq
    print('[Info] Getting gaps flanking sequence...')
    draftgenomedict = quartet_util.readFastaAsDict(draftgenomefile)
    flankingdict = {}
    gapdict = {}
    gapcomponentdict = {}
    for sid, seq in draftgenomedict.items():
        if 'N'*100 in seq:
            i = 1
            gapsitelist = [r.span() for r in re.finditer(r'N+', seq)]
            prev_gap_end = 0
            for gapsite in gapsitelist:
                start = max(gapsite[0] - flanking, 0)
                end = min(gapsite[1] + flanking, len(seq))
                leftseq = seq[start:gapsite[0]]
                rightseq = seq[gapsite[1]:end]
                if 'N'*100 in leftseq or 'N'*100 in rightseq:
                    print(f'[Warning] Flanking sequence of gap {sid}.{i} contains another gap. This indicates two gaps are too close and a very small contig is placed in between.')
                else:
                    flankingdict[f'{sid}.{i}.L'] = seq[start:gapsite[0]]
                    flankingdict[f'{sid}.{i}.R'] = seq[gapsite[1]:end]
                gapdict[f'{sid}.{i}'] = seq[gapsite[0]:gapsite[1]]
                left_component_start = prev_gap_end + 1
                left_component_end = gapsite[0]
                right_component_start = gapsite[1] + 1
                if i < len(gapsitelist):
                    right_component_end = gapsitelist[i][0]
                else:
                    right_component_end = len(seq)
                gapcomponentdict[f'{sid}.{i}'] = {
                    'left_seq': seq[left_component_start-1:left_component_end],
                    'right_seq': seq[right_component_start-1:right_component_end],
                    'left_range': (left_component_start, left_component_end),
                    'right_range': (right_component_start, right_component_end),
                }
                prev_gap_end = gapsite[1]
                i += 1
    if flankingdict == {}:
        print('[Error] Input genome does not have valid gap.')
        sys.exit(0)
    subprocess.run(f'mkdir tmp', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    flankingfastafile = f'tmp/{prefix}.gap.flanking.fasta'
    with open(flankingfastafile, 'w') as f:
        for sid, seq in flankingdict.items():
            f.write(f'>{sid}\n{seq}\n')
    del draftgenomedict
    
    # process gapfilling file(s)
    gapcloserdict = {}
    joinattemptdict = {}
    for gapfillfile in gapclosercontigfilelist:
        gapfiller = os.path.basename(gapfillfile)
        print(f'[Info] gapfilling with {gapfiller}...')
        gapfillfasta = quartet_util.readFastaAsDict(gapfillfile)
        
        # If gapfiller itself has gap, split it
        gapfilldict = {}
        needsplit = False
        for sid, seq in gapfillfasta.items():
            if 'N'*100 in seq:
                needsplit = True
                i = 1
                for tig in re.split(r'N{100,}', seq):
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

        # reduce memory 
        with open(f'tmp/{prefix}.gapfillfasta.fasta', 'w') as tmpgapfillfasta:
            for sid, seq in gapfillfasta.items():
                tmpgapfillfasta.write(f'>{sid}\n{seq}\n')
        del gapfillfasta
        del gapfilldict

        pafgapfillfile = quartet_util.minimap(gapfillfile, flankingfastafile, prefix, f'flank_map_{gapfiller}', minimapoption, False, overwrite, aligner)

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
        gapfillfasta = quartet_util.readFastaAsDict(f'tmp/{prefix}.gapfillfasta.fasta')
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
                    score = (Laln['identity'] + Raln['identity']) / 2
                    if score > bestscore:
                        fillstart = Laln['refend'] + Laln['qrylen'] - Laln['qryend'] + 1
                        fillend = Raln['refstart'] - Raln['qrystart']
                        fillseq = gapfillfasta[Laln['refid']][fillstart-1:fillend]
                        if len(fillseq) > maxfillinglen or fillseq == '':
                            continue
                        gapcloserdict[gapid] = {'sid': f'{gapfillfile}@{Laln["refid"]}', 'range': f'{fillstart}-{fillend}',
                                                'seq': fillseq if Laln['strand'] == '+' else quartet_util.reversedseq(fillseq), 'strand': Laln['strand'], 
                                                'score': score}
                        bestscore = score
                elif Laln['refend'] >= Raln['refstart'] and enablejoin == True:
                    joinscore = (Laln['identity'] + Raln['identity']) / 2
                    if bestscore != 0:
                        continue
                    componentinfo = gapcomponentdict.get(gapid)
                    if componentinfo == None:
                        continue
                    min_overlap = min(joinminoverlap, len(componentinfo['left_seq']), len(componentinfo['right_seq']))
                    overlapinfo = evaluate_suffix_prefix_overlap(
                        componentinfo['left_seq'],
                        componentinfo['right_seq'],
                        joinminidentity,
                        min_overlap,
                        flanking,
                    )
                    if overlapinfo['status'] == 'no_overlap':
                        update_joinattempt(joinattemptdict, gapid, 'Failed', {
                            'sid': f'{gapfillfile}@{Laln["refid"]}',
                            'strand': Laln['strand'],
                            'score': joinscore,
                            'reason': 'no_overlap',
                            'left_range': componentinfo['left_range'],
                            'right_range': componentinfo['right_range'],
                            'join_overlap': overlapinfo['length'],
                            'join_identity': overlapinfo['identity'],
                        })
                        continue
                    if overlapinfo['status'] == 'overlap_too_short':
                        update_joinattempt(joinattemptdict, gapid, 'Failed', {
                            'sid': f'{gapfillfile}@{Laln["refid"]}',
                            'strand': Laln['strand'],
                            'score': joinscore,
                            'reason': 'overlap_too_short',
                            'left_range': componentinfo['left_range'],
                            'right_range': componentinfo['right_range'],
                            'join_overlap': overlapinfo['length'],
                            'join_identity': overlapinfo['identity'],
                        })
                        continue
                    if overlapinfo['status'] == 'overlap_identity_too_low':
                        update_joinattempt(joinattemptdict, gapid, 'Failed', {
                            'sid': f'{gapfillfile}@{Laln["refid"]}',
                            'strand': Laln['strand'],
                            'score': joinscore,
                            'reason': 'overlap_identity_too_low',
                            'left_range': componentinfo['left_range'],
                            'right_range': componentinfo['right_range'],
                            'join_overlap': overlapinfo['length'],
                            'join_identity': overlapinfo['identity'],
                        })
                        continue
                    if overlapinfo['length'] >= len(componentinfo['right_seq']):
                        update_joinattempt(joinattemptdict, gapid, 'Failed', {
                            'sid': f'{gapfillfile}@{Laln["refid"]}',
                            'strand': Laln['strand'],
                            'score': joinscore,
                            'reason': 'overlap_too_long',
                            'left_range': componentinfo['left_range'],
                            'right_range': componentinfo['right_range'],
                            'join_overlap': overlapinfo['length'],
                            'join_identity': overlapinfo['identity'],
                        })
                        continue
                    else:
                        gapcloserdict[gapid] = {'sid': f'{gapfillfile}@{Laln["refid"]}', 'range': 'join',
                                                'seq': '', 'strand': Laln['strand'], 
                                                'score': joinscore,
                                                'join_overlap': overlapinfo['length'],
                                                'join_identity': overlapinfo['identity']}
                        update_joinattempt(joinattemptdict, gapid, 'Joined', {
                            'sid': f'{gapfillfile}@{Laln["refid"]}',
                            'strand': Laln['strand'],
                            'score': joinscore,
                            'left_range': componentinfo['left_range'],
                            'right_range': componentinfo['right_range'],
                            'join_overlap': overlapinfo['length'],
                            'join_identity': overlapinfo['identity'],
                        })
                        bestscore = gapcloserdict[gapid]['score']
    
    print(f'[Info] All files processed.')
    if gapcloserdict == {}:
        print(f'[Error] No gap can be closed.')
        sys.exit(0)

    # make detail
    print(f'[Info] Generating filled fasta file and statistics...')
    filldetailfile = f'{prefix}.genome.filled.detail'
    joindetailfile = f'{prefix}.genome.join.detail'
    with open(filldetailfile, 'w') as de:
        totalfilledlen = sum([len(gapcloserdict[x]['seq']) for x in gapcloserdict])
        de.write(f'# Gap Closed: {len(gapcloserdict)}\n# Total Filled length: {totalfilledlen}\n# Gap Remains: {len(gapdict)-len(gapcloserdict)}\n')
        de.write('# Seqid\tGap_identifier\tStatus\tCloserTigID\tCloserRange\tCloserLength\tCloserStrand\tCloserIdentity\n')
        for gapid in gapdict:
            if gapid in gapcloserdict:
                info = gapcloserdict[gapid]
                de.write(f'{".".join(gapid.split(".")[:-1])}\t{gapid.split(".")[-1]}\tClosed\t{info["sid"]}\t{info["range"]}\t{len(info["seq"])}\t{info["strand"]}\t{info["score"]}\n')
            else:
                de.write(f'{".".join(gapid.split(".")[:-1])}\t{gapid.split(".")[-1]}\tNot_closed\n')
    print(f'[Output] Filling detail write to: {filldetailfile}')
    joincount = len([gapid for gapid in gapcloserdict if gapcloserdict[gapid]['range'] == 'join'])
    joinfailcount = len([gapid for gapid in joinattemptdict if joinattemptdict[gapid]['status'] == 'Failed'])
    joinfailreasoncount = {
        'no_overlap': 0,
        'overlap_too_short': 0,
        'overlap_identity_too_low': 0,
        'overlap_too_long': 0,
    }
    for gapid in joinattemptdict:
        info = joinattemptdict[gapid]
        if info['status'] == 'Failed' and info['reason'] in joinfailreasoncount:
            joinfailreasoncount[info['reason']] += 1
    with open(joindetailfile, 'w') as jd:
        jd.write(f'# Gap Joined: {joincount}\n')
        jd.write(f'# Gap Join Failed: {joinfailcount}\n')
        jd.write(f'# Failed_no_overlap: {joinfailreasoncount["no_overlap"]}\n')
        jd.write(f'# Failed_overlap_too_short: {joinfailreasoncount["overlap_too_short"]}\n')
        jd.write(f'# Failed_overlap_identity_too_low: {joinfailreasoncount["overlap_identity_too_low"]}\n')
        jd.write(f'# Failed_overlap_too_long: {joinfailreasoncount["overlap_too_long"]}\n')
        jd.write('# Seqid\tGap_identifier\tStatus\tCloserTigID\tCloserStrand\tCloserIdentity\tLeftComponentRange\tRightComponentRange\tRightTrimmedRange\tJoinOverlapLength\tJoinOverlapIdentity\tReason\n')
        for gapid in gapdict:
            if gapid in joinattemptdict:
                info = joinattemptdict[gapid]
                left_range = info.get('left_range', ('NA', 'NA'))
                right_range = info.get('right_range', ('NA', 'NA'))
                if info['status'] == 'Joined':
                    right_keep_start = right_range[0] + info['join_overlap']
                    right_trimmed_range = f'{right_keep_start}-{right_range[1]}'
                    reason = 'NA'
                else:
                    right_trimmed_range = 'NA'
                    reason = info.get('reason', 'unknown')
                jd.write(
                    f'{".".join(gapid.split(".")[:-1])}\t{gapid.split(".")[-1]}\t{info["status"]}\t{info.get("sid", "NA")}\t{info.get("strand", "NA")}\t{info.get("score", "NA")}\t'
                    f'{left_range[0]}-{left_range[1]}\t'
                    f'{right_range[0]}-{right_range[1]}\t'
                    f'{right_trimmed_range}\t{info.get("join_overlap", "NA")}\t{info.get("join_identity", "NA")}\t{reason}\n'
                )
    print(f'[Output] Join detail write to: {joindetailfile}')
    
    # make fasta
    filledfastafile = f'{prefix}.genome.filled.fasta'
    filledragpfile = f'{prefix}.genome.filled.modified.agp'
    draftgenomedict = quartet_util.readFastaAsDict(draftgenomefile)
    closedsidset = {'.'.join(gapid.split('.')[:-1]) for gapid in gapcloserdict}
    with open(filledfastafile, 'w') as w, open(filledragpfile, 'w') as ragp:
        chrfastadict = {}
        for sid, seq in draftgenomedict.items():
            if sid not in closedsidset:
                w.write(f'>{sid}\n{seq}\n')
                chrfastadict[sid] = seq
            else:
                matches = list(re.finditer(r'N{100,}', seq))
                segments = []
                start = 0
                for match in matches:
                    end = match.start()
                    segments.append((seq[start:end], start + 1, end))
                    start = match.end()
                segments.append((seq[start:], start + 1, len(seq)))
                seqlist = []
                if segments[0][0] != '':
                    seqlist.append((segments[0][0], 'P', sid, segments[0][1], segments[0][2], '+'))
                for i in range(len(matches)):
                    gapid = f'{sid}.{i+1}'
                    nextseq, nextstart, nextend = segments[i+1]
                    if gapid in gapcloserdict:
                        if gapcloserdict[gapid]['range'] != 'join':
                            seqlist.append((gapcloserdict[gapid]['seq'], 'W', gapcloserdict[gapid]['sid'], gapcloserdict[gapid]['range'].split('-')[0], gapcloserdict[gapid]['range'].split('-')[1], gapcloserdict[gapid]['strand']))
                            if nextseq != '':
                                seqlist.append((nextseq, 'P', sid, nextstart, nextend, '+'))
                        else:
                            trimlen = gapcloserdict[gapid]['join_overlap']
                            trimmedseq = nextseq[trimlen:]
                            trimmedstart = nextstart + trimlen
                            if trimmedseq != '':
                                seqlist.append((trimmedseq, 'P', sid, trimmedstart, nextend, '+'))
                    else:
                        seqlist.append((gapdict[gapid], 'U' if len(gapdict[gapid]) == 100 else 'N', len(gapdict[gapid]), 'scaffold', 'yes', 'unspecified'))
                        if nextseq != '':
                            seqlist.append((nextseq, 'P', sid, nextstart, nextend, '+'))
                newseq = ''
                partnum = 1
                for subseq, agp5, agp6, agp7, agp8, agp9 in seqlist:
                    if subseq == '':
                        continue
                    ragp.write(f'{sid}\t{len(newseq)+1}\t{len(newseq)+len(subseq)}\t{partnum}\t{agp5}\t{agp6}\t{agp7}\t{agp8}\t{agp9}\n')
                    newseq += subseq
                    partnum += 1
                w.write(f'>{sid}\n{newseq}\n')
                chrfastadict[f'{sid}'] = newseq
    print(f'[Output] Filled genome fasta file write to: {filledfastafile}')
    print(f'[Output] Modified chromosome agp file write to: {filledragpfile}')
                
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
            if 'N'*100 in seq:
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
    if noplot != True:
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
    parser.add_argument('-a', dest='aligner', choices=['minimap2', 'unimap'], default='minimap2', help='Specify alignment program (support minimap2 and unimap), default: minimap2')
    parser.add_argument('-p', dest='prefix', default='quarTeT', help='The prefix used on generated files, default: quarTeT')
    parser.add_argument('-t', dest='threads', default='1', help='Use number of threads, default: 1')
    parser.add_argument('--enablejoin', dest='enablejoin', action='store_true', default=False, help='Enable join mode to close the gaps. (Unstable)')
    parser.add_argument('--joinonly', dest='joinonly', action='store_true', default=False, help='Use only join mode without fill, should be used with --enablejoin.')
    parser.add_argument('--join-min-overlap', dest='join_min_overlap', type=int, default=50, help='The min overlap length required for join on draft genome components (bp), default: 50')
    parser.add_argument('--join-min-identity', dest='join_min_identity', type=float, default=30, help='The min overlap identity required for join on draft genome components (%%), default: 30')
    parser.add_argument('--overwrite', dest='overwrite', action='store_true', default=False, help='Overwrite existing alignment file instead of reuse.')
    parser.add_argument('--minimapoption', dest='minimapoption', default='-x asm5', help='Pass additional parameters to minimap2 program, default: -x asm5')
    parser.add_argument('--noplot', dest='noplot', action='store_true', default=False, help='Skip all ploting.')

    # parse input paramater
    draftgenomefile = quartet_util.decompress(parser.parse_args().draft_genome)
    gapclosercontigfilelist = [quartet_util.decompress(x) for x in parser.parse_args().gapcloser_contig]
    flanking = int(parser.parse_args().flanking_len)
    minalignmentlength2 = int(parser.parse_args().min_alignment_length)
    if minalignmentlength2 > flanking:
        print('[Error] min_alignment_length should be less than flanking_len.')
        sys.exit(0)
    minalignmentidentity2 = float(parser.parse_args().min_alignment_identity) / 100
    if minalignmentidentity2 < 0 or minalignmentidentity2 > 1:
        print('[Error] min_alignment_identity should be within 0~100.')
        sys.exit(0)
    maxfillinglen = int(parser.parse_args().max_filling_len)
    prefix = parser.parse_args().prefix
    threads = parser.parse_args().threads
    minimapoption = parser.parse_args().minimapoption + f' -t {threads}'
    overwrite = parser.parse_args().overwrite
    aligner = parser.parse_args().aligner
    enablejoin = parser.parse_args().enablejoin
    joinonly = parser.parse_args().joinonly
    joinminoverlap = int(parser.parse_args().join_min_overlap)
    if joinminoverlap < 1:
        print('[Error] join_min_overlap should be at least 1.')
        sys.exit(0)
    joinminidentity = float(parser.parse_args().join_min_identity) / 100
    if joinminidentity < 0 or joinminidentity > 1:
        print('[Error] join_min_identity should be within 0~100.')
        sys.exit(0)
    noplot = parser.parse_args().noplot

    # check prerequisites
    quartet_util.check_prerequisite([aligner])

    # run
    args = [draftgenomefile, gapclosercontigfilelist, flanking, minalignmentlength2, minalignmentidentity2, 
            maxfillinglen, prefix, threads, minimapoption, overwrite, enablejoin, joinonly, joinminoverlap, joinminidentity, noplot, aligner]
    print(f'[Info] Paramater: draftgenomefile={draftgenomefile}, gapclosercontigfilelist={gapclosercontigfilelist}, flanking={flanking}, minalignmentlength2={minalignmentlength2}, minalignmentidentity2={minalignmentidentity2}, maxfillinglen={maxfillinglen}, aligner={aligner}, prefix={prefix}, threads={threads}, minimapoption={minimapoption}, overwrite={overwrite}, enablejoin={enablejoin}, joinonly={joinonly}, joinminoverlap={joinminoverlap}, joinminidentity={joinminidentity}, noplot={noplot}')
    quartet_util.run(GapFiller, args)
