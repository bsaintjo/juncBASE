#!/lab/64/bin/python
# compareSampleSets_getASEvents.py
# Author: Angela Brooks
# Program Completion Date:
# Modification Date(s):
# Copyright (c) 2011, Angela Brooks. anbrooks@gmail.com
# All rights reserved.
"""Will do test to compare two sets of groups in the data.

   Input will be an output file from createAS_CountTables.py
"""

import sys
import optparse
import pdb
import os
import random

import numpy as np
import statistics
import scipy.stats as scistats
from statsmodels.stats import multitest
import pandas as pd
import plotnine as pn

# TODO trace all instances and delete
grdevices = None

from helperFunctions import updateDictOfLists

#############
# CONSTANTS #
#############
NA = "NA"
DEF_SIGN_CUTOFF = 0.05
DEF_THRESH = 10
DEF_DPSI_THRESH = 5.0

DEF_START_IDX = 11

PROP_NON_NA = 0.666

DEF_TEST = "Wilcoxon"
NUM_ITERATIONS = 10000

INFINITY = 100000000000000000000000000000000000000000000

#################
# END CONSTANTS #
#################


###########
# CLASSES #
###########
class OptionParser(optparse.OptionParser):
    """
    Adding a method for required arguments.
    Taken from:
    http://www.python.org/doc/2.3/lib/optparse-extending-examples.html
    """

    def check_required(self, opt):
        option = self.get_option(opt)

        # Assumes the option's 'default' is set to None!
        if getattr(self.values, option.dest) is None:
            print("%s option not supplied" % option)
            self.print_help()
            sys.exit(1)


###############
# END CLASSES #
###############


########
# MAIN #
########
def main():
    opt_parser = OptionParser()
    # Add Options. Required options should have default=None
    opt_parser.add_option("--in_prefix",
                          dest="in_prefix",
                          type="string",
                          help="""Prefix of output files created from
                                  createAS_CountTables. In createAS_CountTables
                                  this was the -o option""",
                          default=None)
    opt_parser.add_option("-i",
                          dest="generic_file",
                          type="string",
                          help="""Run statistical tests on a generic table.
                                  A generic file with any type of value can also
                                  be used. The first line should be a header
                                  that starts with # and contains sample names.""",
                          default=None)
    opt_parser.add_option(
        "--generic",
        dest="samp_start_idx",
        type="int",
        help="""The samp_start_idx gives the 0-based index of the
                                  column containing the sample value.""",
        default=None)
    #   opt_parser.add_option("--left_intron",
    #                         dest="left_input",
    #                         type="string",
    #                         help="""Resulting length-normalized file from createAS_CountTables.py, which
    #                                 contains the exclusion and inclusion counts
    #                                 for just the left side of an intron retention
    #                                 event.""",
    #                         default=None)
    #   opt_parser.add_option("--right_intron",
    #                         dest="right_input",
    #                         type="string",
    #                         help="""Resulting length-normalized file from createAS_CountTables.py, which
    #                                 contains the exclusion and inclusion counts
    #                                 for just the right side of an intron retention
    #                                 event.""",
    #                         default=None)
    opt_parser.add_option("--all_psi_output",
                          dest="all_psi_output",
                          type="string",
                          help="""Output file that will contain the PSI values
                                  for all events and samples. The last two
                                  columns will correspond to the raw-pvalue and
                                  corrected p-value. If a generic file is used,
                                  this will be the output file""",
                          default=None)
    opt_parser.add_option("--thresh",
                          dest="threshold",
                          type="float",
                          help="""Threshold for minimum abundance
                                  in an event. Default=%d""" % DEF_THRESH,
                          default=DEF_THRESH)
    opt_parser.add_option(
        "--mt_correction",
        dest="mt_method",
        type="string",
        help="""Multiple testing correction Method: "BH" - Benjamini & Hochberg,
                                  "bonferroni".  Must select these strings as
                                  the option""",
        default=None)
    opt_parser.add_option("--which_test",
                          dest="which_test",
                          type="string",
                          help="""Which test to use. Either "t-test" or
                                  "Wilcoxon". Default=%s""" % DEF_TEST,
                          default=DEF_TEST)
    opt_parser.add_option("--permutation",
                          dest="permutation",
                          action="store_true",
                          help="""Will do permutation tests to get empircal
                                  p-value""",
                          default=False)
    opt_parser.add_option("--samp2batch",
                          dest="samp2batch_file",
                          type="string",
                          help="""If doing a permutation test, will account for
                                  potential batch effects""",
                          default=None)
    opt_parser.add_option(
        "--delta_thresh",
        dest="delta_thresh",
        type="float",
        help="""Minimum PSI(or generic value) difference between the maximum
                                  and minimum values for a given event to be
                                  considered a change. This
                                  should probably be less than the delta
                                  threshold used to filter significantly
                                  associated events. Default=%s""" %
        DEF_DPSI_THRESH,
        default=DEF_DPSI_THRESH)
    opt_parser.add_option("--sample_set1",
                          dest="sample_set1",
                          type="string",
                          help="""Comma delimited list of samples in set 1
                                  or a file with a list of names, one per line. 
                                  Names must be in header columns of input
                                  files.""",
                          default=None)
    opt_parser.add_option("--sample_set2",
                          dest="sample_set2",
                          type="string",
                          help="""Comma delimited list of samples in set 2
                                  or a file with a list of names, one per line.
                                  Names must be in header columns of input
                                  files.""",
                          default=None)
    opt_parser.add_option("--as_only",
                          dest="as_only",
                          action="store_true",
                          help="""Will output the psi table just to get a sense
                                  of alternative splicing. It will not perform
                                  any statistical analyses.
                                  Names must be in header columns of input
                                  files.""",
                          default=False)
    opt_parser.add_option(
        "--html_dir",
        dest="html_dir",
        type="string",
        help="""Optional: location to put html output table and
                                 associated images""",
        default=None)
    opt_parser.add_option(
        "--html_out_sign_thresh",
        dest="sign_thresh",
        type="float",
        help="""Significance threshold of q-value for printed out   
                                 html_table. DEF=%.2f""" % DEF_SIGN_CUTOFF,
        default=DEF_SIGN_CUTOFF)
    opt_parser.add_option(
        "--pdf",
        dest="make_pdf",
        action="store_true",
        help="""Optional: Will create images as pdf instead of
                                 .png as the default.""",
        default=None)

    (options, args) = opt_parser.parse_args()

    # validate the command line arguments
    #    opt_parser.check_required("-i")
    opt_parser.check_required("--all_psi_output")
    opt_parser.check_required("--mt_correction")
    #   opt_parser.check_required("--sample_set1")
    #   opt_parser.check_required("--sample_set2")

    if options.in_prefix:
        prefix = options.in_prefix
        input_file = open("%s_AS_exclusion_inclusion_counts_lenNorm.txt" %
                          prefix)
        left_input_file_name = "%s_left_intron_counts_lenNorm.txt" % prefix
        right_input_file_name = "%s_right_intron_counts_lenNorm.txt" % prefix
    else:
        if not options.generic_file:
            print("Must include either --in_prefix or -i")
            opt_parser.print_help()
            sys.exit(1)

        input_file = open(options.generic_file)
        left_input_file_name = None
        right_input_file_name = None

    sum_thresh = options.threshold

    delta_thresh = options.delta_thresh

    permutation = options.permutation
    samp2batch = None
    if options.samp2batch_file:
        samp2batch = parseBatchFile(options.samp2batch_file)

    # TODO Disable HTML output for now
    html_out_dir = options.html_out_dir
    if html_out_dir:
        print("HTML output has been disabled at this time", file=sys.stderr)
        html_out_dir = False

    html_out_table_name = None
    if html_out_dir:
        html_out_dir = formatDir(html_out_dir)
        if not os.path.exists(html_out_dir):
            os.mkdir(html_out_dir)
        html_out_table_name = html_out_dir + "/index.html"
    sign_thresh = options.sign_thresh

    html_out = None
    if html_out_table_name:
        html_out = open(html_out_table_name, "w")
        initiateHTML_table(html_out)

    image_file_type = "png"
    if options.make_pdf:
        image_file_type = "pdf"

    as_only = options.as_only
    if not as_only:
        opt_parser.check_required("--sample_set1")
        opt_parser.check_required("--sample_set2")

    in_sample_set1 = options.sample_set1
    in_sample_set2 = options.sample_set2

    # JuncBASE table default
    samp_start_idx = 11
    isGeneric = False
    if options.samp_start_idx:
        samp_start_idx = options.samp_start_idx
        isGeneric = True

    if permutation and isGeneric:
        print("Permutation test is only for JuncBASE tables")
        opt_parser.print_help()
        sys.exit(1)

    left_input_file = None
    right_input_file = None
    if left_input_file_name is None:
        print(
            "Warning:) No intron retention file given as input.  Will not calculate IR events."
        )
    else:
        left_input_file = open(left_input_file_name)
        right_input_file = open(right_input_file_name)

    all_psi_output = open(options.all_psi_output, "w")

    method = options.mt_method
    if method != "BH" and method != "bonferroni":
        print("Wrong) method indicated.")
        opt_parser.print_help()
        sys.exit(1)

    method_map = {"BH": "fdr_bh", "bonferroni": "bonferroni"}
    method = method_map[method]

    which_test = options.which_test
    if which_test != "Wilcoxon" and which_test != "t-test":
        print("Wrong method indicated.")
        opt_parser.print_help()
        sys.exit(1)

    test_map = {"Wilcoxon": scistats.wilcoxon, "t-test": scistats.ttest_ind}
    which_test = test_map[which_test]

    idx2sample = {}

    # {event_type:(set1_medianPSI, set2medianPSI),]}
    event_type2PSI_vals_4_set = {}

    # {event:psi_vals_idx}
    event2PSI_val_idx = {}

    # {event_type:[pval]}
    event_type2pvals = {}

    # {event::pval_idx}
    event2idx = {}

    # {event:{col:psi}}
    event2col2psi = {}

    # {event:{col:sum_counts}}
    header = None
    total_samples = None
    for line in input_file:
        line = formatLine(line)

        if line.startswith("#"):
            header = line
            headerList = header.split("\t")
            if html_out:
                writeHTMLHeader(html_out, headerList)
            sampleList = headerList[samp_start_idx:]
            # Get sample idx
            for i in range(len(sampleList)):
                idx2sample[i] = sampleList[i]

            if as_only:
                # These are arbitrarily chosen and not really used
                in_sample_set1 = sampleList[0]
                in_sample_set2 = sampleList[1]

            # If there were no batches, all samples are in the same batch
            if permutation:
                if samp2batch is None:
                    samp2batch = {}.fromkeys(sampleList, 0)
#           for sample in sample_set1:
#               idx2sample[sampleList.index(sample)] = sample
#           for sample in sample_set2:
#               idx2sample[sampleList.index(sample)] = sample

            sample_set1 = getSamples(in_sample_set1)
            sample_set2 = getSamples(in_sample_set2)

            sample_set1_checked = checkSamples(sampleList, sample_set1)
            sample_set2_checked = checkSamples(sampleList, sample_set2)

            # The threshold for the number of samples that need to have expressed AS
            # events in order to consider testing
            samp_set_thresh1 = float(len(sample_set1_checked)) * PROP_NON_NA
            samp_set_thresh2 = float(len(sample_set2_checked)) * PROP_NON_NA

            if permutation:
                # batch2setLabels : {batch:{"idx":[indexes in batch],
                #                            "samp_set":[parallele list indicating which sample set it is in]}
                (batch2setLabels,
                 batch2len) = buildBatchDict(sampleList, samp2batch,
                                             sample_set1_checked,
                                             sample_set2_checked)

            continue

        line_list = line.split("\t")

        event = "\t".join(line_list[0:samp_start_idx])
        counts = line_list[samp_start_idx:]
        if permutation:
            total_counts = []

        if event in event2idx:
            print("Warning: Skipping duplicate event: %s" % event)
            continue

        if isGeneric:
            event_type = "generic"
        else:
            event_type = getEventType(event)

        if event_type not in event_type2pvals:
            event_type2pvals[event_type] = []
            event_type2PSI_vals_4_set[event_type] = []

        total_samples = len(counts)

        # Fill PSI dict
        min_psi = INFINITY
        max_psi = -INFINITY
        set1_psis = []
        set2_psis = []
        all_psis = []
        na_count = 0
        for i in range(total_samples):
            if isGeneric:
                # psi is actually a generic value that is in the table
                psi = counts[i]
            else:
                (psi, sum_ct) = getPSI_sample_sum(counts[i], sum_thresh)
            if psi != NA:
                psi_val = float(psi)
                all_psis.append(psi_val)
                if psi_val < min_psi:
                    min_psi = psi_val
                if psi_val > max_psi:
                    max_psi = psi_val
            else:
                all_psis.append(NA)
                na_count += 1

            if event in event2col2psi:
                event2col2psi[event][i] = psi
            else:
                event2col2psi[event] = {i: psi}

            if isGeneric:
                if psi < sum_thresh:
                    continue
            else:
                # Compare samples groups together in a wilcoxon rank sum test
                [col_excl, col_incl] = list(map(int, counts[i].split(";")))

                total_count = col_excl + col_incl
                if permutation:
                    total_counts.append(total_count)
                if total_count < sum_thresh:
                    continue
                # Both samples have to be non-zero
#               if belowThreshold(sum_thresh, col_excl, col_incl):
#                   continue

            if idx2sample[i] in sample_set1:
                if event2col2psi[event][i] != NA:
                    set1_psis.append(event2col2psi[event][i])
            elif idx2sample[i] in sample_set2:
                if event2col2psi[event][i] != NA:
                    set2_psis.append(event2col2psi[event][i])

        if as_only:
            if (float(total_samples - na_count) / total_samples) < PROP_NON_NA:
                continue
        else:
            if len(set1_psis) <= samp_set_thresh1 or len(
                    set2_psis) <= samp_set_thresh2:
                continue

        if (max_psi - min_psi) < delta_thresh:
            continue

        if as_only:
            cur_len = len(event_type2pvals[event_type])
            event_type2pvals[event_type].append(1.0)
            event2idx[event] = cur_len
            psi_vals_cur_len = len(event_type2PSI_vals_4_set[event_type])
            event_type2PSI_vals_4_set[event_type].append((0.0, 0.0))
            event2PSI_val_idx[event] = psi_vals_cur_len
            continue

        psi_vals_cur_len = len(event_type2PSI_vals_4_set[event_type])
        event_type2PSI_vals_4_set[event_type].append(
            (statistics.median(set1_psis), statistics.median(set2_psis)))

        event2PSI_val_idx[event] = psi_vals_cur_len

        # Calculate p-val for intron retention later
        if event_type == "intron_retention":
            continue

#        cur_len2 = len(event_type2col2pvals[event_type][j])

#           if event in event2pairs2idx:
#               event2pairs2idx[event][(0,j)] = cur_len
#           else:
#               event2pairs2idx[event] = {(0,j):cur_len}

#           if event in event2col2idx:
#               event2col2idx[event][j] = cur_len2
#           else:
#               event2col2idx[event] = {j:cur_len2}
#

        cur_len = len(event_type2pvals[event_type])

        try:
            if permutation:
                null_dist = get_null_dist(total_counts, all_psis, which_test,
                                          batch2setLabels, batch2len)

                this_stat = which_test(set1_psis, set2_psis)[0]

                raw_pval = get_emp_pval(null_dist, this_stat)
            else:
                raw_pval = which_test(set1_psis, set2_psis)[1]
        except:
            continue

        if np.isnan(raw_pval):
            continue

        event_type2pvals[event_type].append(raw_pval)
        event2idx[event] = cur_len

    # Now calculate intron retention
    if (not as_only) and (not isGeneric):
        if left_input_file:
            left_events2counts = getIntronLeftRightCounts(
                left_input_file, samp_start_idx)
            right_events2counts = getIntronLeftRightCounts(
                right_input_file, samp_start_idx)
        else:
            left_events2counts = {}
            right_events2counts = {}

        for event in left_events2counts:
            if event not in right_events2counts:
                continue

            # If the event is not in this dictionary, the sum of the left and
            # right counts did not pass the thresholds.
            if event not in event2PSI_val_idx:
                continue

            set1_psis_left = []
            set2_psis_left = []
            set1_psis_right = []
            set2_psis_right = []

            left_total_counts = []
            right_total_counts = []
            left_all_psis = []
            right_all_psis = []

            left_min_psi = 200
            left_max_psi = -1
            right_min_psi = 200
            right_max_psi = -1
            for j in range(total_samples):
                try:
                    [left_col_excl, left_col_incl
                     ] = list(map(int,
                                  left_events2counts[event][j].split(";")))
                    [right_col_excl, right_col_incl] = list(
                        map(int, right_events2counts[event][j].split(";")))
                except:
                    pdb.set_trace()

                left_total = left_col_excl + left_col_incl
                right_total = right_col_excl + right_col_incl
                left_total_counts.append(left_total)
                right_total_counts.append(right_total)
                # Both samples have to be non-zero
                #               if (belowThreshold(sum_thresh, left_col_excl, left_col_incl)
                #                                  or
                #                   belowThreshold(sum_thresh, right_col_excl, right_col_incl)):
                #                   continue

                (left_psi,
                 sum_ct) = getPSI_sample_sum(left_events2counts[event][j],
                                             sum_thresh)
                (right_psi,
                 sum_ct) = getPSI_sample_sum(right_events2counts[event][j],
                                             sum_thresh)

                if left_psi != NA:
                    left_psi_val = float(left_psi)
                    left_all_psis.append(left_psi_val)
                    if left_psi_val < left_min_psi:
                        left_min_psi = left_psi_val
                    if left_psi_val > left_max_psi:
                        left_max_psi = left_psi_val
                else:
                    left_all_psis.append(NA)

                if right_psi != NA:
                    right_psi_val = float(right_psi)
                    right_all_psis.append(right_psi_val)
                    if right_psi_val < right_min_psi:
                        right_min_psi = right_psi_val
                    if right_psi_val > right_max_psi:
                        right_max_psi = right_psi_val
                else:
                    right_all_psis.append(NA)

                if left_total < sum_thresh or right_total < sum_thresh:
                    continue

                if idx2sample[j] in sample_set1:
                    if left_psi != NA:
                        set1_psis_left.append(left_psi)
                    if right_psi != NA:
                        set1_psis_right.append(right_psi)
                elif idx2sample[j] in sample_set2:
                    if left_psi != NA:
                        set2_psis_left.append(left_psi)
                    if right_psi != NA:
                        set2_psis_right.append(right_psi)

            if len(set1_psis_left) <= samp_set_thresh1 or len(set1_psis_right) <= samp_set_thresh1\
                or len(set2_psis_left) <= samp_set_thresh2 or len(set2_psis_right) <= samp_set_thresh2:
                continue

            if (left_max_psi - left_min_psi) < delta_thresh:
                continue
            if (right_max_psi - right_min_psi) < delta_thresh:
                continue

            cur_len = len(event_type2pvals["intron_retention"])

            try:
                if permutation:
                    null_dist = get_null_dist(left_total_counts, left_all_psis,
                                              which_test, batch2setLabels,
                                              batch2len)
                    this_stat = which_test(set1_psis, set2_psis)[0]
                    left_pval = get_emp_pval(null_dist, this_stat)

                    null_dist = get_null_dist(right_total_counts,
                                              right_all_psis, which_test,
                                              batch2setLabels, batch2len)
                    this_stat = which_test(set1_psis, set2_psis)[0]
                    right_pval = get_emp_pval(null_dist, this_stat)
                else:
                    left_pval = which_test(set1_psis_left, set2_psis_left)[1]

                    right_pval = which_test(set1_psis_right,
                                            set2_psis_right)[1]
            except:
                continue

            if np.isnan(left_pval) or np.isnan(right_pval):
                continue
            else:
                combined_pval = (left_pval +
                                 right_pval) - left_pval * right_pval

            event_type2pvals["intron_retention"].append(combined_pval)
            event2idx[event] = cur_len

    # All pairs have been evaluated, so now do multiple testing correction on
    # everything
    event_type2adjusted_pvals = {}
    event_type2col2adjusted_pvals = {}

    # Used for printing boxplots
    data_counter = 0

    for event_type in event_type2pvals:
        if as_only:
            event_type2adjusted_pvals[event_type] = list(
                event_type2pvals[event_type])
        else:
            event_type2adjusted_pvals[event_type] = list(
                multitest.multipletests(event_type2pvals[event_type],
                                        method=method)[1])

    # Now go through all events and print out pvals
    all_psi_output.write(header)
    if as_only:
        all_psi_output.write("\n")
    else:
        all_psi_output.write(
            "\tset1_med\tset2_med\tdelta_val\traw_pval\tcorrected_pval\n")

    for event in event2idx:
        if isGeneric:
            event_type = "generic"
        else:
            event_type = getEventType(event)

        this_idx = event2idx[event]
        if this_idx == NA:
            psi_vals = []
            for i in range(total_samples):
                psi_vals.append(event2col2psi[event][i])

            outline = "%s\t%s\tNA\tNA\n" % (event, "\t".join(psi_vals))

            all_psi_output.write(outline)
            continue

        psi_vals = []
        for i in range(total_samples):
            psi_vals.append(event2col2psi[event][i])

        outline = "%s\t%s" % (event, "\t".join(psi_vals))

        if as_only:
            outline += "\n"
            all_psi_output.write(outline)
            continue

        # Add median PSI and delta PSI values
        this_psi_vals_idx = event2PSI_val_idx[event]
        outline += "\t%.2f\t%.2f\t%.2f" % (
            event_type2PSI_vals_4_set[event_type][this_psi_vals_idx][0],
            event_type2PSI_vals_4_set[event_type][this_psi_vals_idx][1],
            event_type2PSI_vals_4_set[event_type][this_psi_vals_idx][1] -
            event_type2PSI_vals_4_set[event_type][this_psi_vals_idx][0])

        outline += "\t%f\t%f\n" % (
            event_type2pvals[event_type][this_idx],
            event_type2adjusted_pvals[event_type][this_idx])

        all_psi_output.write(outline)

        if html_out:
            if event_type2adjusted_pvals[event_type][this_idx] < sign_thresh:
                data_counter = printDataToHTML(grdevices, html_out_dir,
                                               html_out, outline,
                                               samp_start_idx, idx2sample,
                                               sample_set1, sample_set2,
                                               data_counter, image_file_type)

    all_psi_output.close()

    sys.exit(0)


############
# END_MAIN #
############


#############
# FUNCTIONS #
#############
def belowThreshold(sum_thresh, col_excl, col_incl):

    if col_excl + col_incl < sum_thresh:
        return True

    return False


def buildBatchDict(sampleList, samp2batch, sample_set1_samps,
                   sample_set2_samps):
    """
    batch2setLabels : {batch:{"idx":[indexes in batch],
                              "samp_set":[parallele list indicating which sample set it is in]}
    batch2len : dictionary to hold the length of the batch to prevent
                recalculating this later
    """
    batch2setLabels = {}
    num_samples = len(sampleList)
    for i in range(num_samples):
        samp = sampleList[i]
        batch = samp2batch[samp]
        if samp in sample_set1_samps:
            samp_set = 0
        elif samp in sample_set2_samps:
            samp_set = 1
        else:
            continue

        if batch not in batch2setLabels:
            batch2setLabels[batch] = {"idx": [i], "samp_set": [samp_set]}
        else:
            batch2setLabels[batch]["idx"].append(i)
            batch2setLabels[batch]["samp_set"].append(samp_set)

    batch2len = {}
    for batch in batch2setLabels:
        batch2len[batch] = len(batch2setLabels[batch]["idx"])

    return batch2setLabels, batch2len


def checkSamples(sampleList, sample_set):
    checkedSamples = []

    sampleList_set = set(sampleList)

    for samp in sample_set:
        if samp not in sampleList:
            print("Warning: Sample in sample set not in data: %s" % samp)
            continue
        checkedSamples.append(samp)

    return checkedSamples


def formatDir(i_dir):
    i_dir = os.path.realpath(i_dir)
    if i_dir.endswith("/"):
        i_dir = i_dir.rstrip("/")
    return i_dir


def formatLine(line):
    line = line.replace("\r", "")
    line = line.replace("\n", "")
    return line


def get_emp_pval(null_dist, this_stat):
    """
    Two-tailed emp_pval
    """

    high_ctr = 0
    low_ctr = 0

    for stat in null_dist:
        if this_stat > stat:
            high_ctr += 1
        if this_stat < stat:
            low_ctr += 1

    p_val = None
    if high_ctr < low_ctr:
        p_val = 2 * (float(high_ctr) / NUM_ITERATIONS)
    elif low_ctr < high_ctr:
        p_val = 2 * (float(low_ctr) / NUM_ITERATIONS)
    elif high_ctr == 0 and low_ctr == 0:
        p_val = 1 / NUM_ITERATIONS
    elif high_ctr == low_ctr:
        p_val = 1.0

    return p_val


def get_null_dist(total_counts, all_psis, which_test, batch2setLabels,
                  batch2len):
    stats = []
    for i in range(NUM_ITERATIONS):
        this_idx2sample_set = {}

        for batch in batch2setLabels:
            samp_shuffle = list(batch2setLabels[batch]["samp_set"])
            random.shuffle(samp_shuffle)

            # Assign idx2sampset
            for j in range(batch2len[batch]):
                this_idx2sample_set[batch2setLabels[batch]["idx"]
                                    [j]] = samp_shuffle[j]

        # Sample labels have been shuffled
        this_set1_psis = []
        this_set2_psis = []
        for idx in this_idx2sample_set:
            if all_psis[idx] == NA:
                continue
            this_incl = scistats.binom.rvs(n=total_counts[idx],
                                           p=all_psis[idx] / 100,
                                           size=1)
            this_psi = this_incl / total_counts[idx]

            if this_idx2sample_set[idx] == 0:
                this_set1_psis.append(this_psi)
            else:
                this_set2_psis.append(this_psi)

        stats.append(which_test(this_set1_psis, this_set2_psis)[0])
    return stats


def getPSI_sample_sum(excl_incl_ct_str, sum_thresh):
    if excl_incl_ct_str == NA:
        return NA, NA

    excl_str, incl_str = excl_incl_ct_str.split(";")

    try:
        excl = float(excl_str)
        incl = float(incl_str)
    except:
        print("Warning: Bad counts for PSI")
        return (NA, NA)

    if (excl + incl) < sum_thresh:
        return (NA, NA)

    psi = (incl / (incl + excl)) * 100

    psi_str = "%.2f" % psi
    total_ct = "%d" % int(excl + incl)

    return (psi_str, total_ct)


def getEventType(event):
    return event.split("\t")[1]


def getIntronLeftRightCounts(file, samp_start_idx):

    intron_event2counts = {}

    for line in file:
        line = formatLine(line)

        if line.startswith("#"):
            continue

        line_list = line.split("\t")

        event = "\t".join(line_list[0:samp_start_idx])
        counts = line_list[samp_start_idx:]

        # If the reference is NA, then do not calculate
        if counts[0] == NA:
            continue

        intron_event2counts[event] = counts

    return intron_event2counts


def getSamples(sample_set):

    samples = []
    # Check if it is a file
    if os.path.exists(sample_set):
        sample_file = open(sample_set)

        for line in sample_file:
            line = formatLine(line)
            samples.append(line)

        sample_file.close()
    else:
        samples = sample_set.split(",")

    return samples


def initiateHTML_table(html_out):
    html_out.write("""<html><head>
                      <title>compareSampleSets Results</title>
                      </head>
                      <body>
                      <table border="1">\n""")


def makePlot(grdevices, plotName, samp_set1_vals, samp_set2_vals,
             image_file_type):

    samp_vector = ["set1" for i in range(len(samp_set1_vals))]
    samp_vector.extend(["set2" for i in range(len(samp_set2_vals))])

    data_vector = samp_set1_vals + samp_set2_vals

    dframe = pd.DataFrame(list(zip(samp_vector, data_vector)),
                          columns=["sample", "value"])

    gg = pn.ggplot(dframe, pn.aes(x="sample", y="value")) + pn.geom_jitter(
        position="jitter", width=0.2,
        height=0.01) + pn.coord_cartesian(ylim=(0, 100)) + pn.theme_bw()

    # TODO Just infer format from plotName
    gg.save(filename=plotName, format=image_file_type)


def parseBatchFile(batch_file_name):
    samp2batch = {}
    batch_file = open(batch_file_name)
    for line in batch_file:
        line = formatLine(line)
        if line.startswith("#"):
            continue
        if "Batch" in line:
            continue

        samp, batch = line.split("\t")

        samp2batch[samp] = batch

    batch_file.close()

    return samp2batch


def printDataToHTML(grdevices, html_dir, html_out, outline, samp_start_idx,
                    idx2sample, sample_set1, sample_set2, data_counter,
                    image_file_type):

    outline = outline.rstrip("\n")

    # Need to remove the last few values from adding it to the values
    html_outline = ""
    outline_list = outline.split("\t")
    for i in range(len(outline_list)):
        html_outline += "<td>%s</td>" % outline_list[i]

    outline_list = outline_list[:-5]

    samp_set1_vals = []
    samp_set2_vals = []

    for i in range(len(outline_list)):

        idx_shift = i - samp_start_idx
        if idx_shift in idx2sample:
            if idx2sample[idx_shift] in sample_set1:
                if outline_list[i] != NA:
                    samp_set1_vals.append(float(outline_list[i]))
            if idx2sample[idx_shift] in sample_set2:
                if outline_list[i] != NA:
                    samp_set2_vals.append(float(outline_list[i]))

    plotName = "%s/%d.%s" % (html_dir, data_counter, image_file_type)
    makePlot(grdevices, plotName, samp_set1_vals, samp_set2_vals,
             image_file_type)

    image_tag = "<tr><td><a href=\"%d.%s\">link</a></td>" % (data_counter,
                                                             image_file_type)
    html_outline = image_tag + html_outline
    html_outline += "</tr>\n"
    html_out.write(html_outline)

    return data_counter + 1


def writeHTMLHeader(html_out, headerList):
    header = "<tr>"
    header_elems = ["boxplot"] + headerList + [
        "set1_med", "set2_med", "delta_val", "raw_pval", "corrected_pval"
    ]
    for item in header_elems:
        header += "<th>%s</th>" % item

    header += "</tr>\n"
    html_out.write(header)


#################
# END FUNCTIONS #
#################
if __name__ == "__main__": main()
