#!/usr/bin/env python

'''
Usage: call_peak.py [options] <rampagedir>

Options:
    -h --help                      Show help message.
    -v --version                   Show version.
    -l LENGTH                      Feature length for F-seq. [default: 30]
'''

import sys
import os.path
import tempfile
import glob
import operator
from collections import Counter
from seqlib.path import which, check_dir
from seqlib.helper import run_command
from seqlib.ngs import check_bed

__author__ = 'Xiao-Ou Zhang <xiaoou.zhang@umassmed.edu>'
__version__ = '0.0.1'


def fseq(options):
    '''
    Call peaks using F-seq
    '''
    # parse options
    if not which('fseq'):
        sys.exit('Error: No F-seq installed!')
    folder = check_dir(options['<rampagedir>'])
    flength = options['-l']
    percent = [0.95, 0.9, 0.85, 0.8, 0.75, 0.7]
    with open(os.path.join(folder, 'total_counts.txt'), 'r') as f:
        total = int(f.read().rstrip())
    # run F-seq
    flist = {'+': 'rampage_plus_5end.bed', '-': 'rampage_minus_5end.bed'}
    all_peak_f = os.path.join(folder, 'rampage_peaks.txt')
    with open(all_peak_f, 'w') as out:
        for strand in flist:
            peak_f = run_fseq(folder, flist[strand], strand, flength,
                              percent)
            with open(peak_f, 'r') as f:
                for line in f:
                    if total:  # calculate RPM
                        out.write(line.rstrip() + '\t%d\n' % total)
                    else:
                        out.write(line)


def run_fseq(folder, bed, strand, flength, percent):
    prefix = os.path.splitext(bed)[0]
    # create bed files
    temp_dir = tempfile.mkdtemp()
    bed_f = os.path.join(folder, bed)
    # run fseq
    command = 'fseq -f 0 -l %s -of bed -o %s %s' % (flength, temp_dir, bed_f)
    run_command(command, 'Error in F-seq!')
    # cat fseq files
    peak_f = os.path.join(folder, prefix + '_fseq.bed')
    cat_files(temp_dir, peak_f)
    # resize peaks
    resized_peak_f = os.path.join(folder, prefix + '_peak.bed')
    resize_peak(peak_f, bed_f, resized_peak_f, strand, percent)
    return resized_peak_f


def cat_files(temp_dir, outf):
    with open(outf, 'w') as out:
        for fname in glob.iglob(os.path.join(temp_dir, 'chr*')):
            with open(fname) as f:
                out.write(f.read())


def resize_peak(peak, bed, resized_peak, strand, percent):
    bed_f = check_bed(bed)
    with open(peak, 'r') as f, open(resized_peak, 'w') as out:
        for line in f:
            chrom, start, end = line.split()[:3]
            start = int(start)
            end = int(end) + 1
            # count total tags
            total = 0
            sites = []
            for read in bed_f.fetch(chrom, start, end):
                sites.append(int(read.split()[1]))
                total += 1
            if total == 0:
                continue
            # fetch peak tag
            sites = Counter(sites)
            loc, height = sites.most_common()[0]
            # count peak region tags
            peak_sites = 0
            for i in range(loc - 2, loc + 2):
                peak_sites += sites.get(i, 0)
            info = []
            # resize cluster region
            for p in percent:
                required_sites = int(total * p)
                region = []
                for site, num in sorted(sites.items(),
                                        key=operator.itemgetter(1),
                                        reverse=True):
                    region.append(site)
                    region.sort()
                    new_start, new_end = region[0], region[-1] + 1
                    new_total = 0
                    for i in range(new_start, new_end):
                        new_total += sites[i]
                    if new_total >= required_sites:
                        break
                info.append('|'.join(str(x) for x in (new_start, new_end,
                                                      new_total)))
            # output result
            out_format = '%s\t%d\t%d\tpeak\t0\t%s\t%d\t%d\t%d\t%d\t%s\n'
            out.write(out_format % (chrom, start, end, strand, loc,
                                    height, peak_sites, total,
                                    '\t'.join(info)))


if __name__ == '__main__':
    from docopt import docopt
    fseq(docopt(__doc__, version=__version__))
