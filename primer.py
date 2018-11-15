"""
Author Hanqing Liu

This file contain functions for using Primer3

Input:
- A bed file indicates targeted region/sites
- Genome fasta file corresponding to the bed file
- Config file for primer 3
"""

import pandas as pd
import pathlib


def _read_bed(file_path):
    bed_df = pd.read_table(file_path,
                           index_col='region_id', header=None, comment='#',
                           names=['seq_name', 'start', 'end', 'region_id'])
    return bed_df


def _read_fasta_fai(fai_path):
    # read samtools faidx format
    fai_df = pd.read_table(fai_path,
                           index_col=0, header=None,
                           names=['seq_name', 'length', 'start_at',
                                  'line_seq_length', 'line_total_length'])
    return fai_df


def _query_genome(fasta_path, fai_df,
                  seq_name, region_start, region_end,
                  left_expand, right_expand, primer_name):
    # check fai, get position
    if seq_name not in fai_df.index:
        raise KeyError(f'{seq_name} not in the faidx file of genome fasta {fasta_path}')
    if region_end < region_start:
        raise ValueError(
            f'Region end {region_end} < Region Start {region_start} at primer {primer_name}, check your input')

    # calculate length
    seq_start_pos = fai_df.loc[seq_name, 'start_at']
    # in case the region start is close to ref sequence start
    chrom_query_start = max(0, seq_start_pos - left_expand)
    if chrom_query_start == 0:
        real_left_expand = seq_start_pos
    else:
        real_left_expand = left_expand
    chrom_query_length = real_left_expand + (region_end - region_start) + right_expand

    # check some error
    if chrom_query_length > 999999:
        raise ValueError(f'At primer {primer_name}, do you really want to design primer with lenght > 999999?')
    if fasta_path.endswith('.gz'):
        raise NotImplementedError(
            'Genome fasta is gziped, query sequence from gzip file could be super slow, '
            'unzip it should be much faster.')

    # get sequence
    with open(fasta_path) as f:
        f.seek(chrom_query_start)
        sequence_context = ''
        for line in f:
            if line[0] == '>':  # read to next seq
                break
            elif len(sequence_context) >= chrom_query_length:  # read enough length
                break
            else:
                sequence_context += line.strip()
        query_sequence = sequence_context[:min(chrom_query_length, len(sequence_context))]
    if len(sequence_context) >= chrom_query_length:
        real_right_expand = right_expand
    else:
        real_right_expand = len(sequence_context) - real_left_expand - (region_end - region_start)
    if real_right_expand < 0:
        raise ValueError(f'At primer {primer_name}, Region end is outside reference sequence.')

    query_result = pd.Series({  # following primer3 tag names
        'SEQUENCE_NAME': primer_name,
        'SEQUENCE_TEMPLATE': query_sequence,
        'SEQUENCE_TARGET': f'{real_left_expand},{region_end - region_start}'
    })
    return query_result


def prepare_primer3(bed_path, fasta_path, primer3_config,
                    left_expand=None, right_expand=None, both_expand=100,
                    max_length=99999, drop_too_long=False):
    # all the number within primer3 config are not carefully checked, remain this to primer3.
    if left_expand is None:
        left_expand = both_expand
    if right_expand is None:
        right_expand = both_expand
    if left_expand is None or right_expand is None:
        raise ValueError('Specify left_expand & right_expand, or both_expand')

    bed_df = _read_bed(bed_path)
    if not pathlib.Path(fasta_path+'.fai').exists():
        raise FileNotFoundError(f'{fasta_path} do not have .fai index, use samtools faidx to index the file first')
    fai_df = _read_fasta_fai(fasta_path+'.fai')

    primer_records = []
    for primer_name, (seq_name, start, end) in bed_df.iterrows():
        region_length = start - end
        total_length = left_expand + region_length + right_expand
        if region_length > max_length:
            print(f'{primer_name} is dropped due to exceed {max_length} bp')
            continue
        elif total_length > max_length:
            if drop_too_long:
                print(f'{primer_name} is dropped due to exceed {max_length} bp (including expansion)')
                continue
            print(f'{primer_name} expansion is shorten due to exceed {max_length} bp')
            extra_length = max_length - region_length
            left_portion = left_expand / (left_expand + right_expand)
            left_expand = int(extra_length * left_portion)
            right_expand = extra_length - left_expand

        primer_series = _query_genome(fasta_path=fasta_path,
                                      fai_df=fai_df,
                                      seq_name=seq_name,
                                      region_start=start,
                                      region_end=end,
                                      left_expand=left_expand,
                                      right_expand=right_expand,
                                      primer_name=primer_name)
        primer_records.append(primer_series)
    primer_template_df = pd.DataFrame(primer_records).set_index('SEQUENCE_NAME')

    # TODO run primer3_core for each row, parse output
