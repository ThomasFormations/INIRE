#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys

import logging
import textwrap
import random

from collections import defaultdict

from pysam import FastaFile
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def eprint(*args, **kwargs):
    print(*args,  file=sys.stderr, **kwargs)


class Mutator():
    """
    Class associated to a particular mutation experimen

    For a given genome an a list of annotated intervals, constructs a mutated
    genome sequence.

    Each interval is associated with a particular mutation, for example
    chr1    150 200 S
    specifies that the associated genome interval sequence will be shuffled.
    Three possible mutations are defined :
      - S : shuffle
      - I : inversion
      - M : mask

    Parameters
    ----------
    fasta_handle: pysam FastaFile handle
        provide functionnality to fetch sequences (A cache mechanism is
        implemented in order to prevent multiple pysam FastaFile fetch invocation)
    intervals: list of :obj:`BedInterval`
        the specified mutations (interval + mutation type), see BedInterval class
    maximumCached: in optional
        the maximum number of cached sequences
    Attributes
    ----------
    handle: the pysam handle
    references: list
        the list of chromosome names (unused)
    intervals: list
        the list of BedIntervals
    cachedSequences: dict
        the chromosomes sequences stored in a dictionnary
    chromosome_mutations: dict
        a dictionnary storing the number of mutations for each chromosome
    """
    def __init__(self, fasta_handle, intervals, maximumCached=1):
        self.handle = fasta_handle
        self.maximumCached = maximumCached
        self.references = fasta_handle.references
        self.intervals = intervals
        self.cachedSequences = {}
        self.chromosome_mutations = defaultdict(int)

    def flush(self):
        self.cachedSequences = {}

    @property
    def chromosomes(self):
        return self.references

    def fetch(self, chromosome):
        if chromosome not in self.cachedSequences:
            if len(self.cachedSequences) >= (self.maximumCached-1):
                self.cachedSequences = {}
            self.cachedSequences[chromosome] = self.handle.fetch(chromosome)
        return self.cachedSequences[chromosome]

    def modify(self, chrom, sequence):
        self.cachedSequences[chrom] = sequence

    def shuffle(self, inter):
        """"
        Interval will be shuffled
        """
        seq = self.fetch(inter.chrom)
        subseq = seq[inter.start:inter.end]
        shuffled = ''.join(random.sample(subseq, len(subseq)))
        seq = replace_substring(seq, shuffled, inter.start, inter.end)
        self.cachedSequences[inter.chrom] = seq

    def mask(self, inter):
        """"
        Interval will be masked
        """
        seq = self.fetch(inter.chrom)
        masked = 'N' * inter.len
        seq = replace_substring(seq, masked, inter.start, inter.end)
        self.cachedSequences[inter.chrom] = seq

    def invert(self, inter):
        """"
        Interval will be rerverse complemented
        """
        seq = self.fetch(inter.chrom)
        subseq = seq[inter.start:inter.end]
        inverted = str(Seq(subseq).reverse_complement())
        seq = replace_substring(seq, inverted, inter.start, inter.end)
        self.cachedSequences[inter.chrom] = seq

    def mutate(self):
        """
        Mutate the sequence for each interval according to the mutation type
        """
        for interval in self.intervals:
            self.chromosome_mutations[interval.chrom] += 1
            if interval.op == "shuffle":
                self.shuffle(interval)
            elif interval.op == "mask":
                self.mask(interval)
            elif interval.op == "inversion":
                self.invert(interval)
            else:
                self.chromosome_mutations[interval.chrom] -= 1
                raise ValueError("%s is not a valid operation" % interval.op)

    def intervals_complement(self, chrom):
        """
        Constructs the complement of the intervals for a given chromosome
        Equivalent to bedtools complement
        """
        chrom_intervals = [inter for inter in self.intervals if inter.chrom == chrom]
        chrom_len = len(self.fetch(chrom))
        intervals = []
        previous = 0
        for inter in sorted(chrom_intervals, key=lambda x: x.start):
            intervals.append([previous, inter.start])
            previous = inter.end
        intervals.append([previous, chrom_len])
        return intervals

    def get_concatenated_seq(self, intervals, seq):
        """
        Construct the concatenated sequence of the intervals specified by intervals
        For debugging purpose
        Parameters
        ----------
        intervals: list
            a list of 2-dimetional array [start, end]
        seq: str
            the complete sequence of the chromosome
        Returns
        -------
        str
            the concatenated the sequence
        """
        concat_seq = ""
        for int in intervals:
            concat_seq += seq[int[0]:int[1]]
        return concat_seq

    def check(self, chrom):
        """
        Check that the complement of the mutated intervals remains unchanged

        """
        comp_intervals = self.intervals_complement(chrom)
        mutated_seq = self.fetch(chrom)
        mutated_comp_seq = self.get_concatenated_seq(comp_intervals, mutated_seq)
        original_seq = self.handle.fetch(chrom)
        original_comp_seq = self.get_concatenated_seq(comp_intervals, original_seq)
        if mutated_comp_seq != original_comp_seq:
            raise ValueError("Mutations occur outside input mutations in chrom %s" % chrom)
        else:
            eprint("Valid mutated chrom %s" % chrom)

    def get_SeqRecords(self):
        """Returns the set of mutated chromosomes as biopython SeqRecords"""
        seq_records = []
        for chrom in self.chromosomes:
            self.check(chrom)
            seq = self.fetch(chrom)
            num = self.chromosome_mutations[chrom]
            seq_record = SeqRecord(Seq(seq).upper(), id=chrom,
                                   description="mutated chromosome %d mutations" % num)
            seq_records.append(seq_record)
        return seq_records


def replace_substring(seq, newstring, start, end):
    """Replaces in a string, a substring, specified by positions, with a given string
       Both string should have the same size
    """
    if end > len(seq):
        raise ValueError("end: %d outside given string" % end)
    if len(newstring) != end - start:
        raise ValueError("substring does not have the correct size")
    return seq[:start] + newstring + seq[end:]


class BedInterval():
    """
    Tiny class for storing a bed interval with an associated mutation
    """
    def __init__(self, chrom, start, end, strand, operation):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.op = operation

    @property
    def len(self):
        return self.end - self.start

    def __str__(self):
        return "%s\t%d\t%d\t%s\t%s" % (self.chrom, self.start, self.end, self.strand,
                                       self.op)


def read_bed(bedfile):
    """ Read a bed file and strores the associated BedInterval in a list"""
    intervals = []
    with open(bedfile, "r") as fin:
        for line in fin:
            fields = line.strip().split()
            bedint = BedInterval(fields[0], fields[1], fields[2],
                                 fields[5], fields[7])
            intervals.append(bedint)
    return intervals


def main(inputbed, genome, outfasta):
    intervals = read_bed(inputbed)
    mutator = Mutator(FastaFile(genome), intervals)

    mutator.mutate()
    seq_records = mutator.get_SeqRecords()
    SeqIO.write(seq_records, outfasta, "fasta")


def parse_arguments():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent('''\
                                     Mutate a genome fasta sequence according to the mutations specified in a bed file
                                     '''))
    parser.add_argument('--bed',
                        required=True, help='the input bed file file')
    parser.add_argument('--genome',
                        required=True, help='the genome fasta file')
    parser.add_argument('--output',
                        required=True, help='the output fasta file')

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_arguments()

    main(args.bed, args.genome, args.output)
