"""If a read fragment spans the boundary of a simulated SV then potentially one mate will be deleted or duplicated and the other will not.

In this state the bam cannot be converted back into 2 fastq files because there's now a mismatch between the mates.

This script is a bit of a hack in that it just fills in the missing mate for all unpaired mates, but it makes it possible to re generate a simulated fastq.
"""
from collections import defaultdict
import copy
import pysam
import sys


def main():
    inbam, out_bam = sys.argv[1:]
    # First pass, make set of synthetics where only one side is represented
    syn_set = set()
    inbamf = pysam.AlignmentFile(inbam, 'rb')
    for read in inbamf.fetch():
        if '-' not in read.qname:
            # not synthetic
            continue
        first = read.is_read1
        fid = (read.qname, first)
        assert fid not in syn_set, 'i dont get it'
        if (read.qname, not first) in syn_set:
            syn_set.remove((read.qname, not first))
        else:
            syn_set.add(fid)
    # Make dict of originals to synthetics
    qname2synthetics = defaultdict(set)
    for qname, first in syn_set:
        original_qname = qname.split('-')[0]
        qname2synthetics[original_qname].add((qname, first))
    # pass through bam and fill in the needed reads
    with pysam.AlignmentFile(out_bam, 'wb', template=inbamf) as outbamf:
        for read in inbamf.fetch():
            if read.qname not in qname2synthetics:
                outbamf.write(read)
                continue
            for alt_qname, first in qname2synthetics[read.qname]:
                if read.is_read1 != first:
                    new_read = copy.deepcopy(read)
                    new_read.qname = alt_qname
                    outbamf.write(new_read)
main()
