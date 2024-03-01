#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import re
import textwrap
import pickle
import numpy as np
from pyfaidx import Fasta

import orca_predict
from orca_utils import genomeplot
from selene_sdk.sequences import Genome


def parse_ucsc_genomic_region(ucsc_string):
    pattern = re.compile(r'(?P<chromosome>[\w\d]+):(?P<start>\d+)-(?P<end>\d+)')
    match = pattern.match(ucsc_string)

    if match:
        return {
            'chromosome': match.group('chromosome'),
            'start': int(match.group('start')),
            'end': int(match.group('end'))
        }
    else:
        return None


def main(region, fasta, output_prefix):
    """
    """
    orca_predict.load_resources(models=['32M', '256M'], use_cuda=True)

    parsed_region = parse_ucsc_genomic_region(region)
    chrom = parsed_region['chromosome']
    start = parsed_region['start']
    end = parsed_region['end']

    genome = Fasta(fasta)
    sequence = str(genome[chrom][:])

    mpos = 16_000_000
    midpoint = int(len(sequence) / 2)
    wpos = midpoint
    outputs_ref = orca_predict.genomepredict(Genome.sequence_to_encoding(sequence)[None,:,:],
                                         chrom, mpos=mpos, wpos=midpoint,
                                         use_cuda=True)

    model_labels = ["H1-ESC", "HFF"]
    genomeplot(
        outputs_ref,
        show_genes=False,
        show_tracks=False,
        show_coordinates=True,
        model_labels=model_labels,
        file=output_prefix + ".pdf",
        )

    output_pkl = "%s.pkl" % output_prefix
    file = open(output_pkl, 'wb')
    pickle.dump(outputs_ref, file)

    predict = outputs_ref
    resolutions = ["%dMb" % r for r in [32, 16, 8, 4, 2, 1]]
    # Hff is the second prediction
    hff_predictions = predict['predictions'][1]
    for pred, resol in zip(hff_predictions, resolutions):
        output = "%s_predictions_%s.txt" % (output_prefix, resol)
        header = "# Orca=predictions region=%s resol=%s mpos=%d wpos=%d" % (region, resol,
                                                                            mpos, wpos)
        np.savetxt(output, pred, delimiter='\t', header=header, comments='')
    hff_normmats = predict['normmats'][1]
    for pred, resol in zip(hff_normmats, resolutions):
        output = "%s_normmats_%s.txt" % (output_prefix, resol)
        header = "# Orca=normmats region=%s resol=%s mpos=%d wpos=%d" % (region, resol,
                                                                         mpos, wpos)
        np.savetxt(output, pred, delimiter='\t', header=header, comments='')


def parse_arguments():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent('''\
                                     Run glint alignment of query on subject region
                                     '''))
    parser.add_argument('--region',
                        required=True, help='the region (ucsc format)')
    parser.add_argument('--fasta',
                        required=True, help='fasta file')
    parser.add_argument('--outprefix',
                        required=True, help='the output prefix')

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_arguments()

    main(args.region, args.fasta, args.outprefix)
