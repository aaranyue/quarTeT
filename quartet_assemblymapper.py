#!/usr/bin/env python3
# Last modified: V1.2.5

import argparse
import subprocess
import re
import os
import sys
import quartet_util
                     
### MAIN PROGRAM ###
def AssemblyMapper(args):
    refgenomefile, qryfile, mincontiglength, minalignmentlength, minalignmentidentity, prefix, threads, aligner, nofilter, keep, groupcontig, chimera, plot, noplot, overwrite, nucmeroption, deltafilteroption, minimapoption, teclade, teminrepeattimes = args
    
    # split scaffolds to contigs and remove short contigs
    print('[Info] Filtering contigs input...')
    subprocess.run(f'mkdir tmp', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    inputdict = quartet_util.readFastaAsDict(qryfile)
    totaldict = inputdict.copy()
    contigsdict = {}
    if nofilter != True:
        for tigid, seq in inputdict.items():
            if 'N'*100 in seq:
                i = 1
                for tig in re.split(r'N{100,}', seq):
                    totaldict[f'{tigid}.{i}'] = tig
                    if len(tig) >= mincontiglength:
                        contigsdict[f'{tigid}.{i}'] = tig
                    i += 1
            elif len(seq) >= mincontiglength:
                contigsdict[tigid] = seq
        contigfile = f'tmp/{prefix}.contig.fasta'
        with open(contigfile, 'w') as c:
            for tigid, seq in contigsdict.items():
                c.write(f'>{tigid}\n{seq}\n')
    else:
        contigfile = qryfile
        contigsdict = inputdict

    # check telomere in contigs
    print('[Info] Checking telomere in contigs...')
    telofile = f'tmp/{prefix}.tig.telo.info'
    if not os.path.exists(telofile) or overwrite == True:
        subprocess.run(f'python3 {sys.path[0]}/quartet_teloexplorer.py -i {contigfile} -p {prefix}.tig --noplot -c {teclade} -m {teminrepeattimes}', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        subprocess.run(f'mv -t tmp/ -f {prefix}.tig.telo.info', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    monopolize = []
    forceleft = []
    forceright = []
    refdict = quartet_util.readFastaAsDict(refgenomefile)
    minchrlen = min([len(y) for x, y in refdict.items()])
    if os.path.exists(telofile):
        with open(telofile, 'r') as telo:
            for line in telo:
                if line.startswith('#'):
                    continue
                tigid, tiglen, status = line.split()[0:3]
                if status == 'both' and int(tiglen) >= minchrlen / 2:
                    monopolize.append(tigid)
                elif status == 'left':
                    forceleft.append(tigid)
                elif status == 'right':
                    forceright.append(tigid)
    else:
        print('[Warning] Cannot identify telomeres in contigs.')

    # reduce memory   
    with open(f'tmp/{prefix}.totaldict.fasta', 'w') as tmptotaldict:
        for sid, seq in totaldict.items():
            tmptotaldict.write(f'>{sid}\n{seq}\n')
    del totaldict
    with open(f'tmp/{prefix}.contigsdict.fasta', 'w') as tmpcontigsdict:
        for sid, seq in contigsdict.items():
            tmpcontigsdict.write(f'>{sid}\n{seq}\n')
    del contigsdict
    refdictkey = refdict.keys()
    del refdict
    del inputdict

    # get all alignments
    allAlignment = {}
    doplot = not noplot
    if aligner == 'mummer':
        coordsfile = quartet_util.mummer(refgenomefile, contigfile, prefix, 'contig_map_ref', nucmeroption, deltafilteroption, doplot, overwrite)
        
        with open(coordsfile, 'r') as f:
            for line in f:
                if len(line.split()) != 9:
                    continue
                refstart, refend, qrystart, qryend, reflen, qrylen, identity, refid, qryid = line.split()
                if f'{refid}#{qryid}' not in allAlignment:
                    alignment = {'refid': refid, 'qryid': qryid, 'weight': 0, 'sumposition': 0, 'sumpositive': 0, 'sumnegative': 0, 'score': 0}
                else:
                    alignment = allAlignment[f'{refid}#{qryid}']
                alignment['weight'] += 1
                alignment['sumposition'] += int(refstart)
                alignment['score'] += float(identity) * float(qrylen)
                if qrystart < qryend:
                    alignment['sumpositive'] += int(qrylen)
                else:
                    alignment['sumnegative'] += int(qrylen)
                allAlignment[f'{refid}#{qryid}'] = alignment

    elif aligner == 'minimap2' or aligner == 'unimap':
        paffile = quartet_util.minimap(refgenomefile, contigfile, prefix, 'contig_map_ref', minimapoption, doplot, overwrite, aligner)

        with open(paffile, 'r') as f:
            for line in f:
                qryid, qrylen, qrystart, qryend, strand, refid, reflen, refstart, refend, match, alignlen = line.split()[:11]
                if float(alignlen) < minalignmentlength:
                    continue
                if float(match) / float(alignlen) < minalignmentidentity:
                    continue
                if f'{refid}#{qryid}' not in allAlignment:
                    alignment = {'refid': refid, 'qryid': qryid, 'weight': 0, 'sumposition': 0, 'sumpositive': 0, 'sumnegative': 0, 'score': 0, 'refstart': 0, 'refend': 0}
                else:
                    alignment = allAlignment[f'{refid}#{qryid}']
                alignment['weight'] += int(alignlen)
                alignment['sumposition'] += (int(refstart) + 1) * int(alignlen)
                alignment['score'] += float(match)
                if strand == '+':
                    alignment['sumpositive'] += int(alignlen)
                else:
                    alignment['sumnegative'] += int(alignlen)
                if alignment['refstart'] < int(refstart):
                    alignment['refstart'] = int(refstart)
                if alignment['refend'] < int(refend):
                    alignment['refend'] = int(refend)
                allAlignment[f'{refid}#{qryid}'] = alignment
    if allAlignment == {}:
        print(f'[Error] All alignments are filtered. Recommended to adjust filter arguments.')
        sys.exit(0)

    # analysis each alignment and create map
    print('[Info] Analysising aligments...')
    map = {}
    for key, alignment in allAlignment.items():
        # decide monopolizer contig first
        if alignment['qryid'] in monopolize:
            strand = '+' if alignment['sumpositive'] >= alignment['sumnegative'] else '-'
            if alignment['qryid'] not in map:
                map[alignment['qryid']] = {'refid': alignment['refid'], 'score': alignment['score'], 'position': 0, 'strand': strand, 'refstart': alignment['refstart'], 'refend': alignment['refend']}
            else:
                # solve double mapping by score
                if alignment['score'] > map[alignment['qryid']]['score']:
                    map[alignment['qryid']] = {'refid': alignment['refid'], 'score': alignment['score'], 'position': 0, 'strand': strand, 'refstart': alignment['refstart'], 'refend': alignment['refend']}
    monopolizedchr = [y['refid'] for x, y in map.items()]

    for key, alignment in allAlignment.items():
        if alignment['qryid'] not in monopolize and alignment['refid'] not in monopolizedchr:
            strand = '+' if alignment['sumpositive'] >= alignment['sumnegative'] else '-'
            if (strand == '+' and alignment['qryid'] in forceleft) or (strand == '-' and alignment['qryid'] in forceright):
                position = 0
            elif (strand == '+' and alignment['qryid'] in forceright) or (strand == '-' and alignment['qryid'] in forceleft):
                position = 999999999
            else:
                position = alignment['sumposition'] / alignment['weight']
            if alignment['qryid'] not in map:
                map[alignment['qryid']] = {'refid': alignment['refid'], 'score': alignment['score'], 'position': position, 'strand': strand, 'refstart': alignment['refstart'], 'refend': alignment['refend']}
            else:
                # solve double mapping by score
                if alignment['score'] > map[alignment['qryid']]['score']:
                    map[alignment['qryid']] = {'refid': alignment['refid'], 'score': alignment['score'], 'position': position, 'strand': strand, 'refstart': alignment['refstart'], 'refend': alignment['refend']}

                     
    # write all contigs' destination
    contiginfo = []
    mappedbases, discardedbases = 0, 0
    totaldict = quartet_util.readFastaAsDict(f'tmp/{prefix}.totaldict.fasta')
    contigsdict = quartet_util.readFastaAsDict(f'tmp/{prefix}.contigsdict.fasta')
    for tigid, seq in totaldict.items():
        tiglen = len(seq)
        if tigid in map:
            target = map[tigid]['refid']
            mappedbases += tiglen
        else:
            if tiglen < mincontiglength and nofilter != True:
                target = 'TooShort'
                discardedbases += tiglen
            elif tigid in contigsdict:
                target = 'NoAlignment'
                discardedbases += tiglen
            else:
                target = 'Split'
        contiginfo.append([tigid, tiglen, target])
    contiginfo.sort(key=lambda x: x[0])
    contigmapinfofile = f'{prefix}.contig.mapinfo'
    with open(contigmapinfofile, 'w') as w:
        w.write(f'# Total mapped: {mappedbases} bp\n')
        w.write(f'# Total discarded: {discardedbases} bp\n')
        w.write(f'# ContigID\tContigLength\tTargetID\n')
        for [tigid, tiglen, target] in contiginfo:
            w.write(f'{tigid}\t{tiglen}\t{target}\n')
    print(f'[Output] Mapping result for each contigs write to: {contigmapinfofile}')

    # group contig option
    if groupcontig == True:
        subprocess.run(f'rm -rf {prefix}.groupcontig', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        subprocess.run(f'mkdir {prefix}.groupcontig', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        for [tigid, tiglen, target] in contiginfo:
            file_path = f"{prefix}.groupcontig/{target}.tig.fasta"
            with open(file_path, 'a') as file:
                file.write(f">{tigid}\n{totaldict[tigid]}\n")
        print(f'[Output] grouped contigs write to: {prefix}.groupcontig/')

    # chimera option
    if chimera != 0:
        refdict = quartet_util.readFastaAsDict(refgenomefile)
        with open(f'{prefix}.chimera_contig.fasta', 'w') as nw:
            for tigid in map:
                refid = map[tigid]['refid']
                refstart = map[tigid]['refstart']
                refend = map[tigid]['refend']
                strand = map[tigid]['strand']
                tigseq = totaldict[tigid]
                left = max(refstart-chimera, 1)
                right = min(refend+chimera, len(refdict[refid]))
                if strand == '+':
                    seqout = refdict[refid][left-1:refstart] + tigseq + refdict[refid][refend:right+1]
                else:
                    seqout = quartet_util.reversedseq(refdict[refid][refend:right+1]) + tigseq + quartet_util.reversedseq(refdict[refid][left-1:refstart])
                nw.write(f'>{tigid}_chimera_{refid}_{left}~{refstart}+this+{refend}~{right}\n{seqout}\n')
        print(f'[Output] Chimeric contigs write to: {prefix}.chimera_contig.fasta')

    # check if all Chr have contigs mapped.
    mappedchrlist = []
    for qryid, mapinfo in map.items():
        if mapinfo['refid'] not in mappedchrlist:
            mappedchrlist.append(mapinfo['refid'])
    for id in refdictkey:
        if id not in mappedchrlist:
            print(f'[Warning] {id} failed to be mapped by any contigs.') 

    # write draft genome fasta and agp and stat
    print('[Info] Generating draft genome and statistics...')
    chrfastadict = {}
    agplist = []
    gaplocus = {}
    counts = {}
    qryorder = sorted(map.keys(), key=lambda x: map[x]['position'],)
    for qryid in qryorder:
        seq = contigsdict[qryid] if map[qryid]['strand'] == '+' else quartet_util.reversedseq(contigsdict[qryid])
        newid = f"{map[qryid]['refid']}"
        if newid not in chrfastadict:
            chrfastadict[newid] = seq
            gaplocus[newid] = []
            counts[newid] = 1
            agplist.append([newid, 1, len(seq), counts[newid], 'W', qryid, 1, len(seq), map[qryid]['strand']])
        else:
            gapstart = len(chrfastadict[newid]) + 1
            chrfastadict[newid] += 'N'*100 + seq
            gaplocus[newid].append([gapstart, gapstart + 99])
            counts[newid] += 1
            agplist.append([newid, gapstart, gapstart + 99, counts[newid], 'U', 100, 'scaffold', 'yes', 'align_genus'])
            counts[newid] += 1
            agplist.append([newid, gapstart + 100, gapstart + len(seq) + 99, counts[newid], 'W', qryid, 1, len(seq), map[qryid]['strand']])

    draftgenomefastafile = f'{prefix}.draftgenome.fasta'
    with open(draftgenomefastafile, 'w') as fa:
        for chrs, seq in sorted(chrfastadict.items(), key=lambda x: x[0]):
            fa.write(f'>{chrs}\n{seq}\n')
        if keep == True:
            for [tigid, tiglen, target] in contiginfo:
                if target == 'TooShort' or target == 'NoAlignment':
                    fa.write(f'>{tigid}\n{contigsdict[tigid]}\n')
    if os.path.getsize(draftgenomefastafile) == 0:
        print('[Error] No Chromosome can be assembly.')
        sys.exit(0)
    print(f'[Output] draft genome fasta file write to: {draftgenomefastafile}')
    
    draftgenomeagpfile = f'{prefix}.draftgenome.agp'
    with open(draftgenomeagpfile, 'w') as agp:
        agplist.sort(key=lambda x: x[0])
        for inf in agplist:
            agp.write("\t".join([str(x) for x in inf]) + '\n')
        if keep == True:
            for [tigid, tiglen, target] in contiginfo:
                if target == 'TooShort' or target == 'NoAlignment':
                    agp.write(f"{tigid}\t1\t{tiglen}\t1\tW\t{tigid}\t1\t{tiglen}\t+\n")
    print(f'[Output] draft genome agp file write to: {draftgenomeagpfile}')

    draftgenomestatfile = f'{prefix}.draftgenome.stat'
    with open(draftgenomestatfile, 'w') as info:
        totallen = 0
        gc = 0
        for chrs in chrfastadict:
            totallen += len(chrfastadict[chrs])
            gc += chrfastadict[chrs].count('G') + chrfastadict[chrs].count('C')
        gccontent = gc / totallen
        info.write(f'# Total Size: {totallen}\n')
        info.write(f'# GC content: {gccontent}\n')
        info.write(f'# AssemblyID\tLength\tGapcount\tGaplocus\n')
        for chrs in sorted(gaplocus.keys()):
            gapcount = len(gaplocus[chrs])
            tmpl = []
            for [start, end] in gaplocus[chrs]:
                tmpl.append(f'{start}-{end}')
            locus = '\t'.join(tmpl)
            info.write(f'{chrs}\t{len(chrfastadict[chrs])}\t{gapcount}\t{locus}\n')
    print(f'[Output] draft genome stat write to: {draftgenomestatfile}')
    
    if doplot == True:
        quartet_util.drawgenome(draftgenomeagpfile, f'{prefix}.draftgenome')

    # show colinearity in plot
    if plot == True:
        print('[Info] Aligning draft and ref to plot...')
        if aligner == 'mummer':
            quartet_util.mummer(refgenomefile, draftgenomefastafile, prefix, 'draftgenome_mummer_ref', nucmeroption, deltafilteroption, plot, overwrite)
        if aligner == 'minimap2' or aligner == 'unimap':
            quartet_util.minimap(refgenomefile, draftgenomefastafile, prefix, 'draftgenome_map_ref', minimapoption, plot, overwrite, aligner)

### RUN ###
if __name__ == '__main__':
    # Argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', dest='reference_genome',required=True, help='(*Required) Reference genome file, FASTA format.')
    parser.add_argument('-q', dest='contigs', required=True, help='(*Required) Phased contigs file, FASTA format.')
    parser.add_argument('-c', dest='min_contig_length', type=int, default=50000, help='Contigs shorter than INT (bp) will be removed, default: 50000')
    parser.add_argument('-l', dest='min_alignment_length', type=int, default=10000, help='The min alignment length to be select (bp), default: 10000')
    parser.add_argument('-i', dest='min_alignment_identity', type=float, default=90, help='The min alignment identity to be select (%%), default: 90')
    parser.add_argument('-p', dest='prefix', default='quarTeT', help='The prefix used on generated files, default: quarTeT')
    parser.add_argument('-t', dest='threads', default='1', help='Use number of threads, default: 1')
    parser.add_argument('-a', dest='aligner', choices=['minimap2', 'unimap', 'mummer'], default='minimap2', help='Specify alignment program (support minimap2, unimap and mummer), default: minimap2')
    parser.add_argument('--nofilter', dest='nofilter', action='store_true', default=False, help='Use original sequence input, no filtering.')
    parser.add_argument('--keep', dest='keep', action='store_true', default=False, help='Keep the unplaced contigs in draft genome')
    parser.add_argument('--groupcontig', dest='groupcontig', action='store_true', default=False, help='Add an folder output of contigs grouped by destination.')
    parser.add_argument('--extract-ref-flanks', dest='chimera', type=int, default=0, help='Add an output of chimera contig containing reference flanks of x bp (check issue#42 for detail), default: 0 (off)')
    parser.add_argument('--plot', dest='plot', action='store_true', default=False, help='Plot a colinearity graph for draft genome to reference alignments. (will cost more time)')
    parser.add_argument('--noplot', dest='noplot', action='store_true', default=False, help='Skip all ploting.')
    parser.add_argument('--overwrite', dest='overwrite', action='store_true', default=False, help='Overwrite existing alignment file instead of reuse.')
    parser.add_argument('--minimapoption', dest='minimapoption', default='-x asm5', help='Pass additional parameters to minimap2/unimap program, default: -x asm5')
    parser.add_argument('--nucmeroption', dest='nucmeroption', default='', help='Pass additional parameters to nucmer program.')
    parser.add_argument('--deltafilteroption', dest='deltafilteroption', default='', help='Pass additional parameters to delta-filter program.')
    parser.add_argument('--teclade', dest='te_clade', choices=['plant', 'animal', 'other'], default='other', help='Specify clade of this genome for telomere search. Plant will search TTTAGGG, animal will search TTAGGG, other will use tidk explore\'s suggestion, default: other')
    parser.add_argument('--teminrepeattimes', dest='te_min_repeat_times', type=int, default=100, help='The min repeat times to considered as telomere, default: 100')

    # parse input paramater
    refgenomefile = quartet_util.decompress(parser.parse_args().reference_genome)
    qryfile = quartet_util.decompress(parser.parse_args().contigs)
    mincontiglength = int(parser.parse_args().min_contig_length)
    minalignmentlength = int(parser.parse_args().min_alignment_length)
    minalignmentidentity = float(parser.parse_args().min_alignment_identity) / 100
    if minalignmentidentity < 0 or minalignmentidentity > 1:
        print('[Error] min_alignment_identity should be within 0~100.')
        sys.exit(0)
    prefix = parser.parse_args().prefix
    threads = parser.parse_args().threads
    aligner = parser.parse_args().aligner
    nofilter = parser.parse_args().nofilter
    keep = parser.parse_args().keep
    groupcontig = parser.parse_args().groupcontig
    chimera = parser.parse_args().chimera
    plot = parser.parse_args().plot
    noplot = parser.parse_args().noplot
    overwrite = parser.parse_args().overwrite
    teclade = parser.parse_args().te_clade
    teminrepeattimes = parser.parse_args().te_min_repeat_times
    quartet_util.check_prerequisite(['delta-filter', 'show-coords'])
    if aligner == 'mummer':
        quartet_util.check_prerequisite(['nucmer'])
        nucmeroption = parser.parse_args().nucmeroption + f' -t {threads}'
        deltafilteroption = parser.parse_args().deltafilteroption + f' -l {minalignmentlength} -i {minalignmentidentity}'
        minimapoption = ''
    if aligner == 'minimap2' or aligner == 'unimap':
        quartet_util.check_prerequisite([aligner])
        nucmeroption = ''
        deltafilteroption = ''
        minimapoption = parser.parse_args().minimapoption + f' -t {threads}'
    
    # run
    args = [refgenomefile, qryfile, mincontiglength, minalignmentlength, minalignmentidentity, 
            prefix, threads, aligner, nofilter, keep, groupcontig, chimera, plot, noplot, overwrite, nucmeroption, deltafilteroption, minimapoption, teclade, teminrepeattimes]
    print(f'[Info] Paramater: refgenomefile={refgenomefile}, qryfile={qryfile}, mincontiglength={mincontiglength}, minalignmentlength={minalignmentlength}, minalignmentidentity={minalignmentidentity}, prefix={prefix}, threads={threads}, aligner={aligner}, nofilter={nofilter}, keep={keep}, groupcontig={groupcontig}, chimera={chimera}, plot={plot}, noplot={noplot}, overwrite={overwrite}, nucmeroption={nucmeroption}, deltafilteroption={deltafilteroption}, minimapoption={minimapoption}, teclade={teclade}, teminrepeattimes={teminrepeattimes}')  
    quartet_util.run(AssemblyMapper, args)
