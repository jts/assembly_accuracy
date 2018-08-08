#! /usr/bin/env python

import os
import pysam
import argparse
import vcf
import math
import sys
from string import maketrans

def qscore(p):
    return -10 * math.log10(1 - p)

def var2key(var):
    return var.CHROM + ":" + str(var.POS) + ":" + var.REF + ">" + str(var.ALT[0])

def key2pos(key):
    fields = key.split(":")
    return int(fields[1])

def load_calls(filename):
    calls = dict()
    vcf_reader = vcf.Reader(open(filename, 'r'))
    for record in vcf_reader:
        key = var2key(record)
        calls[key] = record
    return calls

def make_aligned_strings(read, reference_file):

    qseq = read.query_alignment_sequence
    rseq = reference_file.fetch(read.reference_name, read.reference_start, read.reference_end)

    sys.stderr.write("alignment length: " + str(read.reference_end - read.reference_start) + "\n")

    qseq_a_list = list()
    rseq_a_list = list()

    qcurr = 0
    rcurr = 0

    for cigarOp, cigarLength in read.cigar:

        qc = 0
        rc = 0

        # https://www.biostars.org/p/76119/
        if(cigarOp == 0): #match
            qc = cigarLength
            rc = cigarLength
        elif(cigarOp == 1): #insertions
            qc = cigarLength
            rc = 0
        elif(cigarOp == 2): #deletion
            qc = 0
            rc = cigarLength
        elif(cigarOp == 4): #soft clipping
            # ignore since we use query_aligned_sequence
            continue
        elif(cigarOp == 5): # hard clipping
            # ignore
            continue
        else:
            raise ValueError("Unhandled cigar op:" + str(cigarOp))

        # append characters from query/reference to the alignment
        qseq_a_list.append(qseq[qcurr:qcurr+qc])
        rseq_a_list.append(rseq[rcurr:rcurr+rc])

        mc = max(qc, rc)
        qpad = mc - qc
        rpad = mc - rc

        # pad so the lengths match
        qseq_a_list.append("-" * qpad)
        rseq_a_list.append("-" * rpad)

        qcurr += qc
        rcurr += rc

    qaligned = "".join(qseq_a_list).upper()
    raligned = "".join(rseq_a_list).upper()

    return qaligned, raligned

def print_alignment(read, query_aligned, ref_aligned):
    reference_position = read.reference_start
    print "Alignment for %s:%s" % (read.reference_name, reference_position)

    for i in range(0, len(query_aligned), 80):
        print "QUERY    " + query_aligned[i:i+80]
        print "REF      " + ref_aligned[i:i+80] + " " + str(reference_position)
        print ""
        reference_position += len(ref_aligned[i:i+80].replace('-', ''))

def rev_comp_aligned(s):
    trans = maketrans("ACGT-", "TGCA-")
    return s.translate(trans)[::-1]

def print_edits(fp, read, query_aligned, ref_aligned):
    reference_position = read.reference_start
    hard_clip_tag = 5
    query_pre_hard_clip = 0
    if read.cigartuples[0][0] == hard_clip_tag:
        query_pre_hard_clip += read.cigartuples[0][1]

    query_post_hard_clip = 0
    if read.cigartuples[-1][0] == hard_clip_tag:
        query_post_hard_clip += read.cigartuples[-1][1]

    # calculate sequence length, including clips    
    sequence_length = query_pre_hard_clip + len(read.seq) + query_post_hard_clip
    
    #
    query_start_position = read.query_alignment_start + query_pre_hard_clip
    query_end_position = read.query_alignment_end + query_pre_hard_clip

    #print "PREHC: %d POSTHC: %d QS: %d QE: %d QL: %d RS: %d RE: %d" % (query_pre_hard_clip, query_post_hard_clip, query_start_position, query_end_position, sequence_length, read.reference_start, read.reference_end)

    query_position = query_start_position

    # switch strands if rc
    if read.is_reverse:
        query_position = sequence_length - query_end_position
        query_aligned = rev_comp_aligned(query_aligned)
        ref_aligned = rev_comp_aligned(ref_aligned)
    #print "QP: %d (after strand)" % (query_position)

    #print_alignment(read, query_aligned, ref_aligned)

    i = 0
    n = len(ref_aligned)
    while i < n:
        if query_aligned[i] == ref_aligned[i]:
            reference_position += 1
            query_position += 1
            i += 1
        else:
            # new difference found, find next matching base
            j = i
            while query_aligned[j] != ref_aligned[j] and j < n:
                j += 1

            # Get difference strings
            q_sub = query_aligned[i:j]
            r_sub = ref_aligned[i:j]

            q_gaps = q_sub.count("-")
            r_gaps = r_sub.count("-")

            offset = 0
            if q_gaps > 0 or r_gaps > 0:
                # append a single base to the start
                q_sub = query_aligned[i-1] + q_sub.replace("-", "")
                r_sub = ref_aligned[i-1] + r_sub.replace("-", "")
                offset = 1
            
            fp.write("%s\t%d\t.\t%s\t%s\t.\tPASS\t.\n" % (read.query_name, query_position - offset + 1, q_sub.upper(), r_sub.upper()))
            reference_position += (j - i - r_gaps)
            query_position += (j - i - q_gaps)
            i = j

parser = argparse.ArgumentParser( description='Calculate the accuracy of a genome assembly by comparing to a reference')
parser.add_argument('--reference', type=str, required=True)
parser.add_argument('--assembly', type=str, required=True)
parser.add_argument('--variants', type=str, required=False)
parser.add_argument('--min-mapping-quality', type=int, required=False, default=0)
parser.add_argument('--print-alignment', action='store_true')
parser.add_argument('--write-edits', type=str, required=False)
args = parser.parse_args()

out_bam = "assembly_analysis.sorted.bam"
os.system("minimap2 -Y -a -x asm5 %s %s | samtools sort -T assembly_analysis.tmp -o %s -" % (args.reference, args.assembly, out_bam))
os.system("samtools index %s" % (out_bam))

# Calculate variants with htsbox
out_vcf = "assembly_analysis.vcf"
os.system("htsbox pileup -c -V0.0001 -S300 -q10 -Q1 -s1 -f %s %s > %s" % (args.reference, out_bam, out_vcf))

# Read variants VCF to ignore as errors
variants = dict()
if args.variants is not None:
    variants = load_calls(args.variants)

# Read differences with respect to the reference
candidate_errors = load_calls(out_vcf)

# Open reference file
reference_file = pysam.FastaFile(args.reference)

# counters
matches = 0
mismatches = 0
insertions = 0
deletions = 0
error_blocks = 0

edits_fp = None
if args.write_edits is not None:
    edits_fp = open(args.write_edits, "w")

# Calculate the number of matching bases from the alignment
samfile = pysam.AlignmentFile(out_bam)
for read in samfile:

    try:
        if read.mapq < args.min_mapping_quality:
            continue

        query_aligned, ref_aligned = make_aligned_strings(read, reference_file)
        reference_position = int(read.reference_start)

        if args.print_alignment:
            print_alignment(read, query_aligned, ref_aligned)

        if edits_fp is not None:
            print_edits(edits_fp, read, query_aligned, ref_aligned)

        for i in range(0, len(query_aligned)):
            if query_aligned[i] == ref_aligned[i]:
                matches += 1

    except ValueError as inst:
        pass

vcf_reader = vcf.Reader(out_vcf)
vcf_writer = vcf.Writer(open("assembly_analysis.errors.vcf", 'w'), vcf_reader)

for key in candidate_errors:
    if key in variants:
        continue # not an error
   
    error = candidate_errors[key]
    if error.is_snp:
        mismatches += len(error.REF)
    if error.is_indel:
        ref_len = len(error.REF)
        alt_len = len(error.ALT[0])
        if alt_len > ref_len:
            insertions += alt_len - 1
        elif ref_len > alt_len:
            deletions += ref_len - 1
        else:
            mismatches += alt_len # mnp?
    error_blocks += 1
    vcf_writer.write_record(error)

alignment_length = (matches + mismatches + insertions + deletions)
identity = float(matches) / alignment_length
error_block_rate = 1 - (float(error_blocks) / alignment_length)

print "%s percent identity: %.4f (Q%.1f) [matches=%d mismatches=%d insertions=%d deletions=%d] variant rate: %.4lf (Q%.1f)" % (args.assembly, 100 * identity, qscore(identity), matches, mismatches, insertions, deletions, 100*error_block_rate, qscore(error_block_rate))
