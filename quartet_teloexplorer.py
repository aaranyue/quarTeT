#!/usr/bin/env python3

import sys
import os
import math
import argparse
import subprocess
import collections
import quartet_util

### MAIN PROGRAM ###
def teloExplorer(args):
    genomefile, clade, minrepeattimes, prefix = args
    subprocess.run(f'mkdir tmp', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    
    # run tidk explore
    telorangedict = {'plant': '-l 7', 'animal': '-l 6', 'other': '-m 5 -x 12'}
    print('[Info] Running tidk explore...')
    subprocess.run(f'tidk explore -d tmp/ -f {genomefile} -o {prefix} -e tsv {telorangedict[clade]} -t {minrepeattimes}', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

    # check suggested telomere
    with open(f'tmp/{prefix}.txt', 'r') as r:
        con = r.read()
        if len(con.split()) == 3:
            print(f'[Error] No telomere-like repeats found.')
            sys.exit(0)
        for i in range(3, len(con.split())-1, 3):
            telorepeat = con.split()[i]
            if clade == 'plant' and telorepeat not in 'TTTAGGG'*2 and telorepeat not in 'CCCTAAA'*2:
                print(f'[Warning] Telomere repeat found is {telorepeat} instead of TTTAGGG.')
            elif clade == 'animal' and telorepeat not in 'TTAGGG'*2 and telorepeat not in 'CCCTAA'*2:
                print(f'[Warning] Telomere repeat found is {telorepeat} instead of TTAGGG.')
            else:
                print(f'[Info] Telomere repeat found is {telorepeat}.')
                break
        
    # run tidk search
    print('[Info] Running tidk search...')
    subprocess.run(f'tidk search -d tmp/ -f {genomefile} -o {prefix} -e csv -s {telorepeat} -w 10000', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

    # sort telomere info
    print('[Info] Analysising...')
    with open(f'tmp/{prefix}_telomeric_repeat_windows.csv', 'r') as l:
        telodict = {}
        block = collections.defaultdict(list)
        for line in l:
            if not line.startswith('id'):
                chrid, start, forwardrepeatnum, reverserepeatnum, repeatpattern = line.split(',')
                block[chrid].append([int(forwardrepeatnum), int(reverserepeatnum)])
        for chrid in block:
            sidesize = math.floor(len(block[chrid])/2) if len(block[chrid]) < 30 else 15
            lefttelo = block[chrid][0:sidesize]
            totalforwardrepeatnum = 0
            totalreverserepeatnum = 0
            for forwardrepeatnum, reverserepeatnum in lefttelo:
                totalforwardrepeatnum += forwardrepeatnum
                totalreverserepeatnum += reverserepeatnum
            if max(totalforwardrepeatnum, totalreverserepeatnum) >= minrepeattimes:
                if totalforwardrepeatnum >= totalreverserepeatnum:
                    telodict[f'{chrid}.L'] = [totalforwardrepeatnum, '+']
                else:
                    telodict[f'{chrid}.L'] = [totalreverserepeatnum, '-']
                
            righttelo = block[chrid][len(block[chrid])-sidesize:]
            totalforwardrepeatnum = 0
            totalreverserepeatnum = 0
            for forwardrepeatnum, reverserepeatnum in righttelo:
                totalforwardrepeatnum += forwardrepeatnum
                totalreverserepeatnum += reverserepeatnum
            if max(totalforwardrepeatnum, totalreverserepeatnum) >= minrepeattimes:
                if totalforwardrepeatnum >= totalreverserepeatnum:
                    telodict[f'{chrid}.R'] = [totalforwardrepeatnum, '+']
                else:
                    telodict[f'{chrid}.R'] = [totalreverserepeatnum, '-']
    
    teloinfofile = f'{prefix}.telo.info'
    fasta = quartet_util.readFastaAsDict(genomefile)
    with open(teloinfofile, 'w') as t:
        both, side, no = 0, 0, 0
        status = {}
        for chrid in block:
            if f'{chrid}.L' in telodict and f'{chrid}.R' in telodict:
                both += 1
                status[chrid] = 'both'
            elif f'{chrid}.L' in telodict and f'{chrid}.R' not in telodict:
                side += 1
                status[chrid] = 'left'
            elif f'{chrid}.R' in telodict and f'{chrid}.L' not in telodict:
                side += 1
                status[chrid] = 'right'
            else:
                no += 1
                status[chrid] = 'no'
        t.write(f'# Telomere repeat monomer: {telorepeat}\n')
        t.write(f'# Both telomere found: {both}\n')
        t.write(f'# Only one telomere found: {side}\n')
        t.write(f'# No telomere found: {no}\n')
        t.write('# Chrid\tchrlen\tstatus\tleftnum\tleftdirection\trightnum\trightdirection\n')
        for chrid in block:
            leftid = f'{chrid}.L'
            if leftid in telodict:
                leftinfo = f'{telodict[leftid][0]}\t{telodict[leftid][1]}'
            else:
                leftinfo = '0\t'
            rightid = f'{chrid}.R'
            if rightid in telodict:
                rightinfo = f'{telodict[rightid][0]}\t{telodict[rightid][1]}'
            else:
                rightinfo = '0\t'
            t.write(f'{chrid}\t{len(fasta[chrid])}\t{status[chrid]}\t{leftinfo}\t{rightinfo}\n')
    print(f'[Output] Telomere information write to: {teloinfofile}')
                
    agpfile = f'tmp/{prefix}.genome.agp'
    if not os.path.exists(agpfile):
        quartet_util.agpGap(genomefile, agpfile)
    quartet_util.drawgenome(agpfile, f'{prefix}.telo', telofile=teloinfofile)

### RUN ###
if __name__ == '__main__':
    # Argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='genome',required=True, help='(*Required) Genome file to be identified, FASTA format.')
    parser.add_argument('-c', dest='clade', choices=['plant', 'animal', 'other'], default='other', help='Specify clade of this genome. Plant will search TTTAGGG, animal will search TTAGGG, other will use tidk explore\'s suggestion, default: other')
    parser.add_argument('-m', dest='min_repeat_times', type=int, default=100, help='The min repeat times to be reported, default: 100')
    parser.add_argument('-p', dest='prefix', default='quarTeT', help='The prefix used on generated files, default: quarTeT')

    # parse input paramater
    genomefile = quartet_util.decompress(parser.parse_args().genome)
    clade = parser.parse_args().clade
    minrepeattimes = parser.parse_args().min_repeat_times
    prefix = parser.parse_args().prefix

    # check prerequisites
    quartet_util.check_prerequisite(['tidk', 'Rscript'])

    # run
    args = [genomefile, clade, minrepeattimes, prefix]
    quartet_util.run(teloExplorer, args)
    