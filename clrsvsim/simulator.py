import copy
import logging
import os
import random
import re
import string

import itertools
import pyfasta
import pysam
from cigar import Cigar
from preconditions import preconditions

BASES = ['A', 'C', 'G', 'T']
PYSAM_SORT_MEM = '6G'
MIN_INSERT_LEN = 150
CIGAR_STR_PATTERN = re.compile('^(\d+[A-Z])+$')
CIGAR_OP_PATTERN = re.compile('(\d+)([A-Z])')
COMPLEMENTS = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
REF_BUF_LEN = 1000

logger = logging.getLogger(__name__)


def complement(base):
    return COMPLEMENTS.get(base, base)


@preconditions(lambda a, b: len(a) == len(b) == 2 and a[0] <= a[0] and b[0] <= b[1])
def overlap(a, b):
    # Overlap between two integer intervals
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def get_inverse_sequence(input_bam, chrom, start, end, ref_genome_fa=None):
    """
    Get the reverse complement of a region from an alignment file.
    If a base has multiple alleles in the alignment, one base will be selected according to the allele frequencies.
    If the region is not covered in the file, the corresponding bases are fetched from a reference genome, if provided.

    Args:
        input_bam: The alignment containing the region to invert.
        chrom: The chromosome in which the region is located.
        start: Start of inversion (inclusive).
        end: End of inversion (exclusive).
        ref_genome_fa: If provided, will be used to infer bases for locations with no coverage in `input_bam`.

    Returns:
        The reverse complement of the region.
    """

    # Item i in this sequence is a list of all bases observed in the alignment in position i.
    # E.g for the reads {ACC, AGT}:
    # inverse_sequence_bases[0] = ['A', 'A']
    # inverse_sequence_bases[1] = ['C', 'G']
    # inverse_sequence_bases[2] = ['C', 'T']
    inverse_sequence_bases = [[]] * (end - start)

    # Collect all bases found in actual alignments.
    with pysam.Samfile(input_bam, 'rb') as insam:
        for pileupcolumn in insam.pileup(chrom, start, end, truncate=True, stepper='all'):
            pos = pileupcolumn.pos - start
            inverse_sequence_bases[pos] = [
                pileupread.alignment.query_sequence[pileupread.query_position] for pileupread in pileupcolumn.pileups
                if pileupread.query_position is not None and not pileupread.is_del and not pileupread.is_refskip
            ]

    # Select a base for every position; for positions that were not covered in the alignment, backfill from the ref genome, if provided.
    for pos in range(len(inverse_sequence_bases)):
        if not inverse_sequence_bases[pos]:
            base = complement(ref_genome_fa[chrom][start + pos]) if ref_genome_fa and ref_genome_fa.get(chrom) else 'N'
            inverse_sequence_bases[pos].append(base)

    # Reverse and complement.
    inverse_sequence_bases.reverse()
    return "".join([complement(random.choice(bases)) for bases in inverse_sequence_bases])


def random_base():
    return random.choice(BASES)


def _sort_index(unsorted, output_bam):
    pysam.sort("-o", output_bam, "-m", PYSAM_SORT_MEM, unsorted)
    os.unlink(unsorted)
    pysam.index(output_bam)


@preconditions(lambda read, breakpoint: read and 0 <= breakpoint < read.rlen)
def make_split_read(read, breakpoint, clip_left, hard_clip_threshold=1.0, sequence=None):
    """
    Create a split read (a continuous soft-clip from one end of a read until a breakpoint).
    Modifies both the CIGAR string and the actual sequence of the read.

    For example, If the read sequence is `ACACACAC` with a CIGAR of 8M, the breakpoint in position 3, and the sequence provided is `GTGTGT`,
    then if the clipping is a to the left of the breakpoint the modified read will have a CIGAR of 3S5M and its sequence will be `GTGCACAC`,
    and if the clipping is to the right of the breakpoint, the modified read will have a CIGAR of 3M5S and its sequence will be `ACAGTGTG`.

    Args:
        read: The read to modify.
        breakpoint: The breakpoint of the read.
        clip_left: Whether to clip every base to the left of the breakpoint or to the right of it.
        hard_clip_threshold: By default bases are soft-clipped. If more than `hard_clip_threshold` of the read is clipped, hard-clip instead
        sequence: An optional sequence to use for overriding bases in the clipped region.

    Returns:
        The split read.

    """
    split_read = copy.deepcopy(read)
    split_read.qname = split_read.query_name = read.qname + '-' + 'split'

    # CIGAR clipping.
    if read.cigarstring:
        cigar = Cigar(read.cigarstring)
        split_read.cigarstring = str(cigar.mask_left(breakpoint) if clip_left else cigar.mask_right(read.rlen - breakpoint))

    # Convert to hard-clipping, if needed.
    clip_len = breakpoint if clip_left else read.rlen - breakpoint
    if float(clip_len) / read.rlen > hard_clip_threshold:
        soft_clipped_cigar = '{}S'.format(clip_len)
        hard_clipped_cigar = '{}H'.format(clip_len)
        cigar = split_read.cigarstring
        if clip_left and cigar.startswith(soft_clipped_cigar):
            cigar = cigar.replace(soft_clipped_cigar, hard_clipped_cigar, 1)
        elif not clip_left and split_read.cigarstring.endswith(soft_clipped_cigar):
            cigar = cigar[:-len(hard_clipped_cigar)] + hard_clipped_cigar
        split_read.cigarstring = cigar

    if clip_left:
        # adjust reference match.
        split_read.reference_start += breakpoint

    # Sequence replacement.
    if sequence:
        split_seq = list(split_read.seq)
        if clip_left:
            split_seq[:breakpoint] = sequence[-breakpoint:]
        else:
            split_seq[breakpoint:] = sequence[:read.rlen - breakpoint]
        split_read.seq = ''.join(split_seq)

    return split_read


def modify_read(read, snp_rate, insert_rate, del_rate):
    """
    Add random mismatches and insertions/deletions to a read.

    Args:
        read: The read to modify.
        snp_rate: The fraction of bases that will be changed.
        insert_rate: The fraction of bases that will be inserted.
        del_rate: The fraction of bases that will be deleted, if there was no insertion in that position.

    Returns: A tuple of (the modified read, the number of changes made).
    """
    new_read = copy.deepcopy(read)
    # Modify the read name so it's not considered a duplicate.
    rand_str = ''.join([random.choice(string.ascii_uppercase) for _ in range(5)])
    new_read.qname = new_read.query_name = read.qname + '-' + rand_str

    if snp_rate == insert_rate == del_rate == 0:
        return new_read, 0

    seq = list(read.seq)  # The sequence to modify.
    deleted = [False] * len(seq)  # True if the base in position i was deleted.
    modified_bases = set()  # All positions that were modified.

    for i in range(len(seq)):
        if random.random() < snp_rate:
            seq[i] = random.choice(list(set(BASES) - set(read.seq[i])))  # A random base different from the current one.
            modified_bases.add(i)
        if random.random() < del_rate:
            deleted[i] = True
            modified_bases.add(i)
        elif random.random() < insert_rate:
            seq.insert(i, random_base())
            deleted.insert(i, False)
            modified_bases.add(i)

    not_deleted = [seq[i] for i in range(len(seq)) if not deleted[i]]
    new_read.seq = ''.join(not_deleted)
    return new_read, len(modified_bases)


@preconditions(lambda chrom, start, end, ratio_change: chrom and start >= 0 and end >= 0 and ratio_change >= 0 and ratio_change != 1.0)
def modify_copy_number(input_bam, output_bam, chrom, start, end, ref_genome, ratio_change, snp_rate=0, indel_rate=0, split_read_ratio=1.0,
                       random_seed=None):
    """
    Increase or decrease the amount of reads in a region; write and index a modified BAM.

    Args:
        input_bam: The aligned reads to modify.
        output_bam: The modified bam to write.
        chrom: Chromosome to apply modification to.
        start: Start of modification.
        end: End of modification.
        ref_genome: The reference genome. Will be used to specify soft-clipped sequence.
        ratio_change: The fraction to increase or decrease reads. Ratios above 1.0 will result in an increase of reads in the region,
        by duplicating current reads in the region with minor modifications. Ratios below 1.0 will result in a reduction of reads in
        the region by randomly subsampling existing reads.
        snp_rate: The fraction of bases that will be modified if new reads are added.
        indel_rate: The fraction of bases that will be inserted or deleted if new reads are added.
        split_read_ratio: The fraction of reads spanning a breakpoint that will be converted to split reads.
        random_seed: The seed to use for random operations.

    Mates of reads are not modified, if they are outside the impacted region.
    Copy number changes for a region shorter than a read are not supported.

    TODO: Allow supplying a BED for multiple regions.

    Returns:
        A tuple of: (number of reads in the region, number of reads written in the region, number of split reads written)
    """
    fa = pyfasta.Fasta(ref_genome)
    if ratio_change > 1:
        right_clip_seq = str(fa[chrom][start:start + REF_BUF_LEN])
        left_clip_seq = str(fa[chrom][end - REF_BUF_LEN:end])
    elif ratio_change < 1:
        left_clip_seq = str(fa[chrom][start - REF_BUF_LEN:start])
        right_clip_seq = str(fa[chrom][end: end + REF_BUF_LEN])

    def write_copies(read, outsam, num_copies):
        """
        Write multiple copies of a read, with random modifications to each copy.

        Args:
            read: The read to write copies of.
            outsam: The pysam object to write the read to.
            num_copies: The number of times to write the read. If fractional, the read will be written 'num_copies' on average.
            E.g. if num_copies == 0.75, in 75% of calls to write_copies a read will be written,

        Returns:
            The number of copies written.
        """

        copies_written = 0
        orig_read = read

        while num_copies >= 1:
            # Write the original read first, then keep writing modified versions of it.
            outsam.write(read)
            read, _ = modify_read(orig_read, snp_rate, indel_rate / 2, indel_rate / 2)  # assume insertion/deletion equally likely.
            num_copies -= 1
            copies_written += 1

        # Probabilistically write the remaining fractional part of num_copies.
        if random.random() < num_copies:
            outsam.write(read)
            copies_written += 1

        return copies_written

    def write_split_read(read, outsam, num_copies):
        """
        Convert a read to a split read and write multiple copies of it, with random modifications to each copy.

        Args:
            read: The read to write copies of.
            outsam: The pysam object to write the read to.
            num_copies: The number of times to write the read. See `write_copies`.

        Returns:
            The number of copies written.
        """

        spans_left_breakp = read.reference_start < start < read.reference_end
        spans_right_breakp = read.reference_start < end < read.reference_end
        if num_copies < 1:
            if spans_left_breakp and spans_right_breakp:
                # pick one with more matching bp
                left_matching_bp = (start - read.reference_start)
                right_matching_bp = (read.reference_end - end)
                clip_left = left_matching_bp < right_matching_bp
            elif spans_left_breakp:
                clip_left = False
            elif spans_right_breakp:
                clip_left = True
            else:
                raise ValueError('Internal disagreement as to whether read shoud be split.')
            if clip_left:
                breakpoint = end - read.reference_start
            else:
                breakpoint = start - read.reference_start
        elif num_copies > 1:
            if spans_left_breakp:
                breakpoint = start - read.reference_start
                clip_left = True
            else:
                breakpoint = end - read.reference_start
                clip_left = False

        # If the breakpoint is beyond the read, just write the original one and bail.
        # This happens with reads that have significant gaps between matching blocks.
        if breakpoint >= read.rlen:
            outsam.write(read)
            return 1

        # Use the reverse 'alternate sequence' in case of left clipping, so that it always terminates in the same base.
        # To visualize this, in the following reads only the replaced portion is visible:
        #
        # Left clip:                                   # Right clip:
        #                                              #
        #                 breakpoint                   #                 breakpoint
        #                     |                        #                      |
        #                     v                        #                      v
        #             ACGTACGT----------               #           -----------ACGT
        #                TACGT----------------         #               -------ACGTAC
        #                   GT---------------------    #                 -----ACGTACGTA
        sequence_for_replacement = left_clip_seq if clip_left else right_clip_seq
        split_read = make_split_read(read, breakpoint, clip_left, sequence=sequence_for_replacement)

        # Write variants of the read.
        reads_written = 0
        if num_copies >= 1:
            # If this is a duplication, first write the 'original' split read, then generate modifications of it.
            outsam.write(read)
            reads_written = 1 + write_copies(split_read, outsam, num_copies - 1)
        else:
            # Assume heterozygous deletion - that is, a 50% chance of writing the original read rather than the split one.
            reads_written += write_copies(read, outsam, num_copies)
            reads_written += write_copies(split_read, outsam, 1 - num_copies)

        return reads_written

    # Main loop for `modify_copy_number`: iterate over all reads; modify those that are in the affected region.
    reads_in_region = modified_reads_in_region = split_reads = 0
    unsorted = output_bam + '.unsorted'
    with pysam.Samfile(input_bam, 'rb') as insam, pysam.Samfile(unsorted, mode='wb', template=insam) as outsam:
        for read in insam:
            try:
                ref_name = read.reference_name
            except ValueError:
                # Unaligned reads and other edge cases.
                outsam.write(read)
                continue

            if ref_name != chrom:
                # Don't modify reads that are not in the affected region.
                outsam.write(read)
            else:
                if read.reference_start >= start and read.reference_end <= end:
                    # Read is fully within the impacted region - change its copy number.
                    modified_reads_in_region += write_copies(read, outsam, ratio_change)
                    reads_in_region += 1
                elif read.reference_end > start and read.reference_start < end:
                    # Read spans the breakpoint - convert to split read, and change copy number.
                    if random.random() < split_read_ratio:
                        split_reads += write_split_read(read, outsam, ratio_change)
                    else:
                        modified_reads_in_region += write_copies(read, outsam, ratio_change)
                else:
                    # Don't modify reads that are not in the affected region.
                    outsam.write(read)

    _sort_index(unsorted, output_bam)

    return reads_in_region, modified_reads_in_region, split_reads


def modify_read_for_insertion(read, position, sequence, snp_rate, indel_rate, padding=0, max_clip_len=None):
    # Half of the reads will be pre-insertion, the rest post-insertion
    clip_left = random.choice([True, False])

    if clip_left:
        position -= padding / 2
    else:
        position += padding / 2

    breakpoint = position - read.reference_start

    if breakpoint >= read.rlen or breakpoint <= 0:
        return read, 0

    if max_clip_len and int(max_clip_len) > 0:
        clip_len = breakpoint if clip_left else read.rlen - breakpoint
        if clip_len > max_clip_len:
            return read, 0

    # The actual insertion.
    read = make_split_read(read, breakpoint, clip_left, sequence=sequence)

    # Add noise.
    return modify_read(read, snp_rate, indel_rate / 2, indel_rate / 2)


def unpack_cigar(cigarstring):
    """
    Converts a Cigar string into components that are 1-base long, e.g.
    unpack_cigar('3S2M1H') == ['1S', '1S', '1S', '1M', '1M', 1H']
    """

    if not cigarstring or not CIGAR_STR_PATTERN.match(cigarstring):
        raise ValueError("Invalid Cigar string")

    cigar_operations = CIGAR_OP_PATTERN.findall(cigarstring)
    return list(itertools.chain.from_iterable([int(length) * ['1' + operation] for length, operation in cigar_operations]))


def get_max_clip_len(read):
    """Returns the longest clipped portion of a read, to the left or to the right."""

    if not read.cigarstring:
        raise ValueError("Missing Cigar string")

    cigar_tuples = list(Cigar(read.cigarstring).items())
    clip_lengths = [cigar_tuples[i][0] for i in (0, -1) if cigar_tuples[i][1] not in Cigar.ref_consuming_ops]
    return max(clip_lengths) if clip_lengths else 0


def invert_read(read, start, end, sequence, snp_rate, indel_rate, max_clip_len=None):
    """
    Invert (a portion of) a read.

    Args:
        read: The read to modify.
        start: The start of the inversion.
        end: The end of the inversion.
        sequence: The full sequence that is inverted in the sample that the read belong to. This sequence should be provided in its
                  reverse-complement form (e.g. as returned by `get_inverse_sequence`).
        snp_rate: The fraction of bases that will be randomly modified in reads that were modified.
        indel_rate: The fraction of bases that will be randomly inserted or deleted in in reads that were modified.
        max_clip_len: If more than "max_clip_len" of the read would be clipped on either end, return None,
                      since this read would not have been captured.

    Returns:
        A duplicate read to the provided one, where any position covered by the inverted region is replaced with the inversion.
    """
    inv_len = end - start
    if start >= read.reference_end or end <= read.reference_start or inv_len < 2:
        return read, 0

    read_with_inversion = copy.deepcopy(read)
    read_with_inversion.qname = read_with_inversion.query_name = read.qname + '-' + 'inv'

    if read.reference_start <= start < end <= read.reference_end:
        # Read spans the entire inversion.
        left_breakpoint = start - read.reference_start
        right_breakpoint = left_breakpoint + inv_len
        read_with_inversion.seq = "{left}{inv}{right}".format(
            left=read.seq[:left_breakpoint],
            inv="".join(reversed(read.seq[left_breakpoint:right_breakpoint])),
            right=read.seq[right_breakpoint:])

        # Clipped bases in reads must start at a read boundary; choose the closest one.
        # TODO: add a supplemental/secondary read where the shorter region is matched, and the longer one clipped.
        cigar_tuples = unpack_cigar(read.cigarstring)
        if left_breakpoint < read.rlen - right_breakpoint:
            start_clip, end_clip = 0, right_breakpoint
        else:
            start_clip, end_clip = left_breakpoint, read.rlen
        for i in range(start_clip, end_clip):
            cigar_tuples[i] = '1S'

        read_with_inversion.cigarstring = str(Cigar("".join(cigar_tuples)).merge_like_ops())

    elif start <= read.reference_start < read.reference_end <= end:
        # Inversion spans the entire read.
        pos_in_inversion = read.reference_start - start
        inv_seq = sequence[pos_in_inversion:pos_in_inversion + read.rlen]
        read_with_inversion = make_split_read(read_with_inversion, 0, clip_left=False, sequence=inv_seq)

        # If a read was reversed, modify its strand.
        read_with_inversion.is_reverse = not read.is_reverse

    elif start > read.reference_start:
        # Inversion starts mid-read, continuing to the end of it (or past it).
        breakpoint = start - read.reference_start
        read_with_inversion = make_split_read(read_with_inversion, breakpoint, clip_left=False, sequence=sequence)

    elif end < read.reference_end:
        # Inversion starts before the read, continuing into it.
        breakpoint = end - read.reference_start
        read_with_inversion = make_split_read(read_with_inversion, breakpoint, clip_left=True, sequence=sequence)

    if max_clip_len and int(max_clip_len) < get_max_clip_len(read_with_inversion):
        return None, 0

    # Add noise.
    return modify_read(read_with_inversion, snp_rate, indel_rate / 2, indel_rate / 2)


def insert_sequence(input_bam, output_bam, chrom, position, fasta_file,
                    insertion_ratio=0.5, snp_rate=0, indel_rate=0, padding=0, max_clip_len=None, random_seed=None):
    """
    Insert a sequence into an alignment; write and index a modified BAM.
    Inserts must be at least `MIN_INSERT_LEN` long.

    Example:

            Insertion:                               AAAATTTTT

            Reference genome,                      GGGGGGCCCCCCC
            insertion position marked:                   ^

            Padding: 4 bases                           GGCC

            Padded insertion:                    GGCCAAAATTTTTGGCC

            Sample after insertion:          GGGGGGCCAAAATTTTTGGCCCCCCC

            Original reads                          GGGGGC
            (read len = 6)                           GGGGCC
                                                      GGGCCC
                                                       GGCCCC
                                                        GCCCCC

            Possible reads after insertion:
            - Same read len                     GGGCCA
            - No insertion-only reads,           GGCCAA
              as these would not be               GCCAAA
              picked up by the capture             CCAAAA
                                                    CAAAAT
                                                         TTTTTG
                                                          TTTTGG
                                                           TTTGGC
                                                            TTGGCC
                                                             TGGCCC

            After alignment to ref genome:        GGGGGGCCCCCCC
            - Asterisks mark soft clips
                                                  *****G
                                                   ****GG
                                                    ***GGC
                                                     **GGCC
                                                      *GGCCC

                                                     GGGCC*
                                                      GGCC**
                                                       GCC***
                                                        CC****
                                                         C*****

    Args:
        input_bam: The aligned reads to modify.
        output_bam: The modified bam to write.
        chrom: Chromosome to apply modification to.
        position: Where the insertion should occur.
        fasta_file: A fasta file containing the sequence to insert.
        insertion_ratio: The fraction of input reads that would have an insertion; e.g. 0.5 for heterozygous, 1.0 for homozygous.
        snp_rate: The fraction of bases that will be randomly modified in reads that had an insertion.
        indel_rate: The fraction of bases that will be randomly inserted or deleted in in reads that had an insertion.
        padding: The number of bases around the insertion position that will pad the insertion element.
        max_clip_len: If more than "max_clip_len" of the read would lie within the insertion, then don't apply it, return the original read.
        random_seed: The seed to use for random operations.

    Returns:
        A tuple of: (number of reads in the region, number of reads modified in the region, number of bases modified)
    """

    # Load sequence to insert.
    fa = pyfasta.Fasta(fasta_file)
    if not len(fa):
        raise ValueError("No sequences found in {}".format(fasta_file))
    else:
        if len(fa) > 1:
            logger.warn("Multiple regions in {}, using {}".format(fasta_file, fa.keys()[0]))
        sequence = str(fa.values()[0])

    if len(sequence) < MIN_INSERT_LEN:
        raise ValueError("Min insertion size is {}; given sequence is {}".format(MIN_INSERT_LEN, len(sequence)))

    random_seed = int(random_seed) if random_seed else 2016
    random.seed(random_seed)
    unsorted = output_bam + '.unsorted'
    reads_in_region = modified_reads_in_region = modified_bases = 0
    with pysam.Samfile(input_bam, 'rb') as insam, pysam.Samfile(unsorted, mode='wb', template=insam) as outsam:
        for read in insam:
            try:
                ref_name = read.reference_name
            except ValueError:
                # Unaligned reads and other edge cases.
                outsam.write(read)
                continue

            if ref_name != chrom or not (read.reference_start < position < read.reference_end):
                # Don't modify reads that are not in the affected region.
                outsam.write(read)
            else:
                reads_in_region += 1
                if random.random() < insertion_ratio:
                    read, bases_modified_in_read = modify_read_for_insertion(
                        read, position, sequence, snp_rate, indel_rate, padding, max_clip_len)
                    modified_reads_in_region += 1
                    modified_bases += bases_modified_in_read
                outsam.write(read)

    _sort_index(unsorted, output_bam)

    return reads_in_region, modified_reads_in_region, modified_bases


def invert_sequence(input_bam, output_bam, chrom, start, end, inversion_ratio=0.5, ref_genome=None, snp_rate=0, indel_rate=0,
                    max_clip_len=None, random_seed=None):
    """
    Inverts the sequence between [start, end); write and index a modified BAM.

    Example:
                   Orig genome                  Genome with inversion

                12345ACACAC67890                  12345GTGTGT67890
                     ^     ^           -->             ^     ^
                   start  end                        start  end

    Reads (len=5):
                12345                  -->        12345
                 2345A                 -->         2345G
                  345AC                -->          345GT
                   45ACA               -->           45GTG
                    5ACAC              -->            5GTGT
                     ACACA             -->             GTGTG
                      CACAC            -->              TGTGT
                       ACAC6           -->               GTGT6
                        CAC67          -->                TGT67
                         AC678         -->                 GT678
                          C6789        -->                  T6789
                           67890       -->                   67890

    Reads spanning breakpoints will be soft-clipped.

    Args:
        input_bam: The aligned reads to mod ify.
        output_bam: The modified bam to write.
        chrom: Chromosome to apply modification to.
        start: Start of modification.
        end: End of modification.
        inversion_ratio: The fraction of reads in the region that will be inverted.
        ref_genome: The reference genome. If supplied, will be used to get inversion sequences missing from the alignment.
        snp_rate: The fraction of bases that will be modified, if a read is changed.
        indel_rate: The fraction of bases that will be inserted or deleted, if a read is changed.
        random_seed: The seed to use for random operations.

    Returns:
        A tuple of: (number of reads in the region, number of reads modified in the region, number of bases modified)
    """

    random_seed = int(random_seed) if random_seed else 2016
    random.seed(random_seed)

    ref_genome_fa = pyfasta.Fasta(ref_genome) if ref_genome else None
    inverse_sequence = get_inverse_sequence(input_bam, chrom, start, end, ref_genome_fa)

    unsorted = output_bam + '.unsorted'
    reads_in_region = modified_reads_in_region = modified_bases = 0
    with pysam.Samfile(input_bam, 'rb') as insam, pysam.Samfile(unsorted, mode='wb', template=insam) as outsam:
        for read in insam:
            try:
                ref_name = read.reference_name
            except ValueError:
                # Unaligned reads and other edge cases.
                outsam.write(read)
                continue

            if ref_name == chrom and (start <= read.reference_start < end) or (start <= read.reference_end < end):
                reads_in_region += 1
                if random.random() < inversion_ratio:
                    read, bases_modified_in_read = invert_read(read, start, end, inverse_sequence, snp_rate, indel_rate, max_clip_len)
                    if read:
                        modified_reads_in_region += 1
                        modified_bases += bases_modified_in_read

            if read:
                outsam.write(read)

    _sort_index(unsorted, output_bam)

    return reads_in_region, modified_reads_in_region, modified_bases
