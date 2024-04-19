import time
import sys
import math
import subprocess
import os
import re
from collections import defaultdict

def run(func, args):
    starttime = time.time()
    func(args)
    endtime = time.time()
    print('[Info] Complete!')
    timecost = endtime - starttime
    timecostd = math.floor(timecost / (60*60*24))
    timecosth = math.floor(timecost % (60*60*24) / (60*60))
    timecostm = math.floor(timecost % (60*60*24) % (60*60) / 60)
    timecosts = math.floor(timecost % (60*60*24) % (60*60) % 60)
    print(f'[Info] Time Cost: {timecostd}d, {timecosth}h, {timecostm}m, {timecosts}s')
    
def check_prerequisite(prerequisitelist: list):
    # print('Checking prerequisites...')
    prerequisitenotfound = []
    for prerequisite in prerequisitelist:
        cmd = subprocess.run(f'which {prerequisite}', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        if cmd.stdout == b'':
            prerequisitenotfound.append(prerequisite)
        # else:
        #     print(f'{prerequisite} located at: {cmd.stdout.decode("utf-8").strip()}')
    if prerequisitenotfound != []:
        for prerequisite in prerequisitenotfound:
            print(f'[Error] prerequisite not found: {prerequisite}')
        print(f'[Error] Please make sure these software have been installed, exported to $PATH, and authorized executable.')
        sys.exit(0)
    # print('All prerequisites found.')

def decompress(file):
    if 'gzip compressed data' in subprocess.run(f'file {file}', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True).stdout.decode():
        subprocess.run(f'gzip -d {file}', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        file = '.'.join(file.split('.')[:-1])
    return file

def readFastaAsDict(fastafile):
    fastaDict = {}
    fil = open(fastafile, 'r')
    allline = fil.read()
    fil.close()
    eachidseq = allline.split('>')
    for idseq in eachidseq:
        if idseq != '':
            sidraw, seqraw = idseq.split('\n', 1)
            sid = sidraw.split()[0].strip()
            seq = seqraw.replace('\n', '').upper()
            fastaDict[sid] = seq
    return fastaDict

def reversedseq(seq: str):
    seq = seq[::-1]
    seq = seq.replace('A', 'E')
    seq = seq.replace('T', 'A')
    seq = seq.replace('E', 'T')
    seq = seq.replace('C', 'E')
    seq = seq.replace('G', 'C')
    seq = seq.replace('E', 'G')
    return seq

def calculate_cover_length(intervals, left, right):
    # sort and merge intervals
    intervals.sort(key=lambda x: x[0])
    merged = []
    for i in intervals:
        if not merged or merged[-1][1] < i[0]:
            merged.append(i)
        else:
            merged[-1][1] = max(merged[-1][1], i[1])
    # calculate
    cover_length = 0
    for i in merged:
        if i[1] <= left or i[0] >= right:
            continue
        else:
            cover_length += min(i[1], right) - max(i[0], left)
    return cover_length

def changeSuffix(filename, newsuffix):
    namelist = filename.split('.')
    namelist[-1] = newsuffix
    return '.'.join(namelist)

def mummer(reffasta, qryfasta, prefix, suffix, nucmeroption, deltafilteroption, plot, overwrite):
    print('[Info] Running MUMmer...')
    subprocess.run(f'mv -t ./ -f tmp/{prefix}.{suffix}.delta', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    if not os.path.exists(f'{prefix}.{suffix}.delta') or overwrite == True:
        runsub(f'nucmer {nucmeroption} -p {prefix}.{suffix} {reffasta} {qryfasta}', 'nucmer')
    runsub(f'delta-filter {deltafilteroption} -m -q {prefix}.{suffix}.delta > {prefix}.{suffix}.filter.delta', 'delta-filter')
    runsub(f'show-coords -T -H {prefix}.{suffix}.filter.delta > {prefix}.{suffix}.coords', 'show-coords')
    if os.path.getsize(f'{prefix}.{suffix}.coords') == 0:
        print(f'[Error] No alignment found or all alignments are filtered. Recommended to adjust filter arguments.')
        sys.exit(0)
    if plot == True:
        runsub(f'mummerplot -t png -p {prefix}.{suffix} {prefix}.{suffix}.filter.delta', 'mummerplot')
        with open(f'{prefix}.{suffix}.gp', 'r') as gp:
            newgp = ''
            for line in gp:
                line = line.replace('w lp ls', 'w line ls')
                newgp += line
        with open(f'{prefix}.{suffix}.gp', 'w') as gp:
            gp.write(newgp)
        runsub(f'gnuplot {prefix}.{suffix}.gp', 'gnuplot')
        print(f'[Output] Colinearity graph write to: {prefix}.{suffix}.png')
    subprocess.run(f'mkdir tmp', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    subprocess.run(f'mv -t tmp -f {prefix}.{suffix}.gp {prefix}.{suffix}.fplot {prefix}.{suffix}.rplot {prefix}.{suffix}.delta {prefix}.{suffix}.filter.delta {prefix}.{suffix}.coords', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    return f'tmp/{prefix}.{suffix}.coords'

def minimap(reffasta, qryfasta, prefix, suffix, minimapoption, plot, overwrite):
    print('[Info] Running minimap2...')
    subprocess.run(f'mv -t ./ -f tmp/{prefix}.{suffix}.paf', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    if not os.path.exists(f'{prefix}.{suffix}.paf') or overwrite == True:
        cmdr = subprocess.run(f'minimap2 {minimapoption} -c -o {prefix}.{suffix}.paf {reffasta} {qryfasta}', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        if '[morecore]' in cmdr.stderr.decode("utf-8") or cmdr.returncode < 0:
            print(f'[Error] Memory insufficient.')
            sys.exit(0)
        elif cmdr.returncode != 0:
            print(f'[Error] Unexcepted error occur in minimap2 as follow:')
            print(f'cmd: {cmdr.args}')
            print(f'returncode: {cmdr.returncode}')
            print('stdout:')
            print(cmdr.stdout.decode("utf-8"))
            print('stderr:')
            print(cmdr.stderr.decode("utf-8"))
            sys.exit(1)
    if os.path.getsize(f'{prefix}.{suffix}.paf') == 0:
        print(f'[Error] No alignment found.')
        sys.exit(0)
    if plot == True:
        # convert paf to delta
        class Alignment:
            def __init__(self, refstart, refend, qrystart, qryend, cigar, strand, mismatches):
                self.refstart, self.refend = int(refstart) + 1, int(refend)
                self.qrystart  = int(qrystart) + 1 if strand == '+' else int(qryend)
                self.qryend = int(qryend) if strand == '+' else int(qrystart) + 1
                self.strand = strand
                self.mismatches = mismatches
                self.parsed_cigar = [(int(i[0]), i[1]) for i in list(re.findall(re.compile(r'(\d+)([MIDNSHP=X])'), cigar))]
        alns, reflen, qrylen = {}, {}, {}
        with open(f'{prefix}.{suffix}.paf', 'r') as f:
            for line in f:
                fields = line.split('\t')
                refid, qryid = fields[5], fields[0]
                reflen[refid] = int(fields[6])
                qrylen[qryid] = int(fields[1])
                for i in fields[12:]:
                    if i.startswith("cg:Z:"):
                        cs = i[5:]
                    if i.startswith("NM:i:"):
                        nm = int(i[5:])
                x = Alignment(fields[7], fields[8], fields[2], fields[3], cs, fields[4], nm)
                if refid not in alns:
                    alns[refid] = defaultdict(list)
                alns[refid][qryid].append(x)
        with open(f'{prefix}.{suffix}.delta', 'w') as d:
            d.write(f'{reffasta} {qryfasta}\nNUCMER\n')
            for refid in alns.keys():
                for qryid in alns[refid].keys():
                    d.write(f'>{refid} {qryid} {reflen[refid]} {qrylen[qryid]}\n')
                    for aln in alns[refid][qryid]:
                        d.write(f'{aln.refstart} {aln.refend} {aln.qrystart} {aln.qryend} {aln.mismatches} {aln.mismatches} 0\n')
                        offsets = []
                        cigar = aln.parsed_cigar[1:-1] if aln.parsed_cigar[0][1] == 'S' or aln.parsed_cigar[0][1] == 'H' else aln.parsed_cigar[:-1]
                        counter = 1
                        for op in cigar:
                            if op[1] == "M":
                                counter += op[0]
                            elif op[1] == "D":
                                offsets.append(counter)
                                num_I = op[0]
                                for i in range(1, num_I):
                                    offsets.append(1)
                                counter = 1
                            elif op[1] == 'I':
                                offsets.append(-1 * counter)
                                num_I = op[0]
                                for i in range(1, num_I):
                                    offsets.append(-1)
                                counter = 1
                        offsets.append(0)
                        offsets = [str(a) for a in offsets]
                        d.write('\n'.join(offsets) + '\n')
        runsub(f'show-coords -T -H {prefix}.{suffix}.delta > {prefix}.{suffix}.coords', 'show-coords')
        runsub(f'mummerplot -t png -p {prefix}.{suffix} {prefix}.{suffix}.delta', 'mummerplot')
        with open(f'{prefix}.{suffix}.gp', 'r') as gp:
            newgp = ''
            for line in gp:
                line = line.replace('w lp ls', 'w line ls')
                newgp += line
        with open(f'{prefix}.{suffix}.gp', 'w') as gp:
            gp.write(newgp)
        runsub(f'gnuplot {prefix}.{suffix}.gp', 'gnuplot')
        print(f'[Output] Colinearity graph write to: {prefix}.{suffix}.png')
    subprocess.run(f'mkdir tmp', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    subprocess.run(f'mv -t tmp -f {prefix}.{suffix}.gp {prefix}.{suffix}.fplot {prefix}.{suffix}.rplot {prefix}.{suffix}.delta {prefix}.{suffix}.filter.delta {prefix}.{suffix}.paf {prefix}.{suffix}.coords', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    return f'tmp/{prefix}.{suffix}.paf'

def agpGap(infile, outagp):
    fastadict = readFastaAsDict(infile)
    with open(outagp, 'w') as w:
        for sid, seq in fastadict.items():
            gapsitelist = [r.span() for r in re.finditer(r'N+', seq)]
            if gapsitelist == []:
                w.write(f'{sid}\t1\t{len(seq)}\t1\tW\t.\t.\t.\t.\n')
            else:
                sitelist = [(1, gapsitelist[0][0]), (gapsitelist[0][0]+1, gapsitelist[0][1])]
                for i in range(len(gapsitelist)):
                    if i + 1 == len(gapsitelist):
                        sitelist += [(gapsitelist[i][1]+1, len(seq))]
                    else:
                        sitelist += [(gapsitelist[i][1]+1, gapsitelist[i+1][0]), (gapsitelist[i+1][0]+1, gapsitelist[i+1][1])]
                count = 1
                for site in sitelist:
                    ty = 'W' if count % 2 == 1 else 'U'
                    w.write(f'{sid}\t{site[0]}\t{site[1]}\t{count}\t{ty}\t.\t.\t.\t.\n')
                    count += 1

def drawgenome(agpfile, outprefix, centrofile=None, telofile=None):
    subprocess.run(f'mkdir tmp', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    if centrofile == None and os.path.exists(f'{outprefix}.best.candidate'):
        centrofile = f'{outprefix}.best.candidate'
    if telofile == None and os.path.exists(f'{outprefix}.telo.info'):
        telofile = f'{outprefix}.telo.info'  
    with open(agpfile, 'r') as agp:
        agpblock = defaultdict(list)
        chrnumber = 0
        for line in agp:
            if chrnumber == 0 or line.split()[0] != chrid:
                chrnumber += 1
            chrid, seqstart, seqend, num, attr, tigid, tigstart, tigend, strand = line.split()
            agpblock[chrid].append([chrnumber, seqstart, seqend, num, attr, tigid, tigstart, tigend, strand])
    with open(f'tmp/{outprefix}.chr.txt', 'w') as chrfile:
        if centrofile == None: 
            chrfile.write('Chr\tStart\tEnd\n')
            for chrid in agpblock:
                chrfile.write(f'{agpblock[chrid][-1][0]}\t1\t{agpblock[chrid][-1][2]}\n')
        else:
            centrodict = {}
            with open(centrofile, 'r') as c:
                for line in c:
                    if line.startswith('\t') or line.startswith('#'):
                        continue
                    chr, start, end = line.split()[0:3]
                    centrodict[chr] = [start, end]
            chrfile.write('Chr\tStart\tEnd\tCE_start\tCE_end\n')
            for chrid in agpblock:
                if chrid in centrodict:
                    chrfile.write(f'{agpblock[chrid][-1][0]}\t1\t{agpblock[chrid][-1][2]}\t{centrodict[chrid][0]}\t{centrodict[chrid][1]}\n')
                else:
                    chrfile.write(f'{agpblock[chrid][-1][0]}\t1\t{agpblock[chrid][-1][2]}\t0\t0\n')
    haslabel = False
    with open(f'tmp/{outprefix}.label.txt', 'w') as labelfile:   
        labelfile.write('Type\tShape\tChr\tStart\tEnd\tcolor\n')     
        for chrid in agpblock:
            for line in agpblock[chrid]:
                if line[4] != 'W':
                    labelfile.write(f'gap\tbox\t{line[0]}\t{line[1]}\t{line[2]}\tff7f00\n')
                    haslabel = True
        if telofile != None:
            with open(telofile, 'r') as te:
                for line in te:
                    if line.startswith('#'):
                        continue
                    chrid, chrlen, status, leftnum, leftdirection, rightnum, leftdirection = line.split('\t')
                    if int(leftnum) != 0:
                        labelfile.write(f'telomere\ttriangle\t{agpblock[chrid][-1][0]}\t1\t10000\t0000ff\n')
                        haslabel = True
                    if int(rightnum) != 0:
                        labelfile.write(f'telomere\ttriangle\t{agpblock[chrid][-1][0]}\t{int(chrlen)-10000}\t{chrlen}\t0000ff\n')
                        haslabel = True
    with open(f'tmp/{outprefix}.genomedrawer.r', 'w') as r:
        if haslabel == True:
            Rscript = f'''library(RIdeogram)
chr <- read.table("tmp/{outprefix}.chr.txt", sep = "\\t", header = T, stringsAsFactors = F)
label <- read.table("tmp/{outprefix}.label.txt", sep = "\\t", header = T, stringsAsFactors = F)
ideogram(karyotype = chr, label = label, label_type = "marker", output = "{outprefix}.svg")
convertSVG("{outprefix}.svg", file = "{outprefix}", device = "png")'''
            r.write(Rscript)
        else:
            RscriptNolabel = f'''library(RIdeogram)
chr <- read.table("tmp/{outprefix}.chr.txt", sep = "\\t", header = T, stringsAsFactors = F)
ideogram(karyotype = chr, output = "{outprefix}.svg")
convertSVG("{outprefix}.svg", file = "{outprefix}", device = "png")'''
            r.write(RscriptNolabel)
    cmdr = subprocess.run(f'Rscript tmp/{outprefix}.genomedrawer.r', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    if cmdr.returncode != 0:
        print(f'[Warning] Unexcepted error occur in Rscript figure drawing as follow:')
        print(f'cmd: {cmdr.args}')
        print(f'returncode: {cmdr.returncode}')
        print('stdout:')
        print(cmdr.stdout.decode("utf-8"))
        print('stderr:')
        print(cmdr.stderr.decode("utf-8"))
    else:
        print(f'[Output] Chromosome plot write to: {outprefix}.png')
    
def runsub(cmd, name):
    cmdr = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    if cmdr.returncode != 0:
        print(f'[Error] Unexcepted error occur in {name} as follow:')
        print(f'cmd: {cmdr.args}')
        print(f'returncode: {cmdr.returncode}')
        print('stdout:')
        print(cmdr.stdout.decode("utf-8"))
        print('stderr:')
        print(cmdr.stderr.decode("utf-8"))
        sys.exit(1)
