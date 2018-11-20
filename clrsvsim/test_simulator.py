import unittest
from os import path

import pysam
from cigar import Cigar
from mock import Mock
from pyfasta import Fasta

from simulator import (
    make_split_read,
    modify_read,
    modify_read_for_insertion,
    invert_read,
    unpack_cigar,
    get_max_clip_len,
    get_inverse_sequence,
    overlap
)

TEST_DATA_DIR = path.join(path.dirname(path.realpath(__file__)), 'test_data')


class SplitReadTest(unittest.TestCase):

    def test_make_split_read(self):
        read = Mock()
        read.seq = 'A' * 20
        read.qual = '*' * len(read.seq)
        read.rlen = len(read.seq)
        read.qname = read.query_name = 'name'
        read.reference_start = 100
        read.cigarstring = '20M'
        alternate_seq = 'C' * 10 + 'T' * 10

        split_read = make_split_read(read, 5, True, sequence=alternate_seq)
        self.assertEqual(split_read.seq, 'T' * 5 + 'A' * 15)
        self.assertEqual(split_read.cigarstring, '5S15M')

        split_read = make_split_read(read, 5, False, sequence=alternate_seq)
        self.assertEqual(split_read.seq, 'A' * 5 + 'C' * 10 + 'T' * 5)
        self.assertEqual(split_read.cigarstring, '5M15S')

        split_read = make_split_read(read, 5, False, hard_clip_threshold=0.1, sequence=alternate_seq)
        self.assertEqual(split_read.seq, 'A' * 5 + 'C' * 10 + 'T' * 5)
        self.assertEqual(split_read.cigarstring, '5M15H')

        split_read = make_split_read(read, 5, False, hard_clip_threshold=0.9, sequence=alternate_seq)
        self.assertEqual(split_read.seq, 'A' * 5 + 'C' * 10 + 'T' * 5)
        self.assertEqual(split_read.cigarstring, '5M15S')

    def test_make_split_read_bam_file(self):
        sorted_bam = path.join(TEST_DATA_DIR, 'sorted.bam')
        with pysam.Samfile(sorted_bam, 'rb') as samfile:
            for read in samfile:
                if not read.cigarstring:
                    continue

                for breakpoint in (10, 50, 100):
                    if breakpoint >= read.rlen:
                        continue

                    for is_left_split in (True, False):
                        split_read = make_split_read(read, breakpoint, is_left_split)
                        cigar_items = list(Cigar(split_read.cigarstring).items())
                        clipped_item = cigar_items[0] if is_left_split else cigar_items[-1]
                        min_clip_len = breakpoint if is_left_split else read.rlen - breakpoint  # Can be longer if adjacent to another clip.
                        self.assertGreaterEqual(clipped_item[0], min_clip_len)
                        self.assertIn(clipped_item[1], ('S', 'H'))  # Will be soft-clipped unless already hard-clipped.

    def test_modify_read(self):
        read = Mock()
        read.seq = 'AAAAA'
        read.qname = 'test'

        # SNPs
        modified, changes = modify_read(read, 1, 0, 0)
        self.assertEqual(changes, len(modified.seq))
        self.assertEqual(len(modified.seq), len(read.seq))
        self.assertTrue(all([read.seq[i] != modified.seq[i] for i in range(len(read.seq))]))

        # Insertions
        modified, changes = modify_read(read, 0, 1, 0)
        self.assertEqual(changes, len(read.seq))
        self.assertEqual(len(modified.seq), len(read.seq) * 2)

        # Deletions
        modified, changes = modify_read(read, 0, 0, 1)
        self.assertEqual(changes, len(read.seq))
        self.assertEqual(len(modified.seq), 0)

    def test_modify_read_for_insertion(self):
        read = Mock()
        read.seq = 'AAAAAA'
        read.qual = '*' * len(read.seq)
        read.qname = 'test'
        read.rlen = len(read.seq)
        read.reference_start = 100
        read.cigarstring = '{}M'.format(read.rlen)

        ins_position = 103
        ins_seq = 'CCCCCCCC'

        modified, changes = modify_read_for_insertion(read, ins_position, ins_seq, 0, 0)
        self.assertEqual(changes, 0)
        # Read can either be modified to be on the left or right of the insertion
        self.assertIn(modified.seq, ('AAACCC', 'CCCAAA'))
        self.assertIn(modified.cigarstring, ('3M3S', '3S3M'))

        # Test padding
        modified, _ = modify_read_for_insertion(read, ins_position, ins_seq, 0, 0, padding=2)
        self.assertIn(modified.seq, ('AAAACC', 'CCAAAA'))
        self.assertIn(modified.cigarstring, ('4M2S', '2S4M'))

        # Insertion positions beyond the read boundaries should not modify it
        for position in (0, 1000):
            modified, _ = modify_read_for_insertion(read, position, ins_seq, 0, 0, )
            self.assertEqual(read, modified)

        # Test limiting the maximum clip length
        modified, _ = modify_read_for_insertion(read, ins_position, ins_seq, 0, 0, max_clip_len=3)
        self.assertIn(modified.cigarstring, ('3M3S', '3S3M'))  # clip len below max - allowed
        modified, _ = modify_read_for_insertion(read, ins_position, ins_seq, 0, 0, max_clip_len=2)
        self.assertEqual(read, modified)  # clip len above max - insertion should not happen

    def test_unpack_cigar(self):
        for bad_cigar_string in (None, '', 'ok', '1', '1s2m', '1S2'):
            self.assertRaises(ValueError, unpack_cigar, bad_cigar_string)

        for cigar, unpacked in [
            ('1M', ['1M']),
            ('2M', ['1M', '1M']),
            ('2M1S', ['1M', '1M', '1S']),
            ('100S', ['1S'] * 100)
        ]:
            self.assertEqual(unpack_cigar(cigar), unpacked)

    def test_get_max_clip_len(self):
        read = Mock()
        read.cigarstring = None

        self.assertRaises(ValueError, get_max_clip_len, read)

        for cigar, max_len in [
            ('4M', 0),
            ('1S2M', 1),
            ('1S1M2S', 2),
            ('2S1M1S', 2),
            ('1M1S', 1)
        ]:
            read.cigarstring = cigar
            self.assertEqual(get_max_clip_len(read), max_len)

    def test_invert_read(self):
        read = Mock()
        read.seq = '123456'
        read.qual = '*' * len(read.seq)
        read.qname = 'test'
        read.rlen = len(read.seq)
        read.reference_start = 100
        read.reference_end = read.reference_start + read.rlen
        read.cigarstring = '{}M'.format(read.rlen)

        def assert_inversion(read, start, end, sequence, expected_seq, expected_cigar):
            inv, _ = invert_read(read, start, end, sequence, 0, 0)
            msg_prefix = 'invert({}, {}--{})'.format(read.seq, start, end)
            self.assertEqual(inv.seq, expected_seq, '{}: {} != {}'.format(msg_prefix, inv.seq, expected_seq))
            self.assertEqual(inv.cigarstring, expected_cigar, '{}: {} != {}'.format(msg_prefix, inv.cigarstring, expected_cigar))

        # Inversions that are fully within the read
        for start, end, expected_seq, expected_cigar in [
            # fully contained, away from borders
            (102, 104, '124356', '2M4S'),
            (101, 104, '143256', '4S2M'),
            # fully contained, touching borders
            (100, 104, '432156', '4S2M'),
            (102, 106, '126543', '2M4S'),
            # spanning exactly the read
            (100, 106, '654321', '6S'),
            # edge cases
            (102, 102, '123456', '6M'),
            (102, 103, '123456', '6M'),
        ]:
            assert_inversion(read, start, end, '', expected_seq, expected_cigar)

        # The sequence that's inverted in the entire genome; only a subset will appear in each read.
        sequence = '9876543210'

        for start, expected_seq, expected_cigar in [
            # inversion ends before the read
            (80, '123456', '6M'),
            (90, '123456', '6M'),

            # inversion starts before the read, and extends into it
            (91, '023456', '1S5M'),
            (92, '103456', '2S4M'),
            (93, '210456', '3S3M'),
            (94, '321056', '4S2M'),
            (95, '432106', '5S1M'),

            # read is fully contained in the inversion
            (96, '543210', '6S'),
            (97, '654321', '6S'),
            (98, '765432', '6S'),
            (99, '876543', '6S'),
            (100, '987654', '6S'),

            # inversion starts mid-read
            (101, '198765', '1M5S'),
            (102, '129876', '2M4S'),
            (103, '123987', '3M3S'),
            (104, '123498', '4M2S'),
            (105, '123459', '5M1S'),

            # inversion starts past the read
            (106, '123456', '6M'),
            (110, '123456', '6M'),
        ]:
            assert_inversion(read, start, start + len(sequence), sequence, expected_seq, expected_cigar)

    def test_get_inverse_sequence(self):
        #   Reads represented in this file (start position = 100):
        #
        #   ACGTACGTAC
        #   ACGTCCGTAC
        #    CGTCCGTACT
        #    CGTCCGAACT
        #     GTCCGAACTT
        #      TCCGAACTTC
        #       CCGAACTTAA
        #       CCGAACTTAA
        #       CCGAACTTAG
        #        CGAACTTAGC
        #
        bam = path.join(TEST_DATA_DIR, 'sv_sim.bam')
        self.assertEqual(get_inverse_sequence(bam, '1', 100, 102), 'GT')
        self.assertEqual(get_inverse_sequence(bam, '1', 108, 111), 'AGT')

        # Inversion of an area with no reads, no ref genome provided
        self.assertEqual(get_inverse_sequence(bam, '1', 98, 102), 'GTNN')
        self.assertEqual(get_inverse_sequence(bam, '1', 0, 100), 'N' * 100)

        # Inversion of an area with no reads, ref genome provided
        ref_genome_fa = Fasta(path.join(TEST_DATA_DIR, 'sv_sim.fa'))
        self.assertEqual(get_inverse_sequence(bam, '1', 0, 4, ref_genome_fa), 'AAAA')
        self.assertEqual(get_inverse_sequence(bam, '1', 98, 102, ref_genome_fa), 'GTAA')

    def test_overlap(self):
        self.assertEqual(overlap((0, 0), (0, 0)), 0)
        self.assertEqual(overlap((0, 1), (0, 1)), 1)
        self.assertEqual(overlap((0, 1), (1, 1)), 0)
        self.assertEqual(overlap((0, 1), (0, 2)), 1)
        self.assertEqual(overlap((0, 1), (1, 2)), 0)
        self.assertEqual(overlap((0, 2), (1, 2)), 1)
        self.assertEqual(overlap((0, 2), (1, 3)), 1)
        self.assertEqual(overlap((0, 2), (0, 3)), 2)
        self.assertEqual(overlap((0, 2), (2, 4)), 0)
        self.assertEqual(overlap((0, 3), (1, 2)), 1)
        self.assertEqual(overlap((0, 4), (1, 3)), 2)
        self.assertEqual(overlap((0, 4), (2, 4)), 2)
        self.assertEqual(overlap((0, 4), (0, 2)), 2)

    # TODO: add tests for:
        # inversion directly from BAM
        # inversion of an area that has no reads in the BAM
        # max clip len
