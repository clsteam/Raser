#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" guess sequence quality code

"""

# Standard library
import itertools
import collections
import operator
import fileinput


class GuessEncode(object):
    """
    Guess the encoding of a stream of qual lines.
    """

    def __init__(self, stdin_or_file="-"):
        """
        :param stdin_or_file: STDIN stream or fastq file (support gz or bz2 format)
        ex:
            'samtools view <BAM file> | cut -f 11 |'
            'cat <FASTQ file> |'
        """
        self.N_MOST_COMMON_THRESH = 4
        self.RANGES = {
            'Sanger': (33, 73, "phred33"),
            'Illumina-1.8': (33, 74, "phred33"),
            'Solexa': (59, 104, "solexa64"),
            'Illumina-1.3': (64, 104, "phred64"),
            'Illumina-1.5': (66, 105, "phred64")
        }
        self.stdin_or_file = stdin_or_file
        self.gmin = 99
        self.gmax = 0

    def judge(self):
        """
        :return: phred33 phred64 solexa64
        """
        pipe_input = fileinput.input(files=self.stdin_or_file, openhook=fileinput.hook_compressed)
        result = self._pall(pipe_input)
        if result:
            pipe_input.close()
            return result
        else:
            # Run the second time
            result = self._pall(pipe_input, 10000)
            pipe_input.close()
            if result:
                return result
            else:
                # default
                return "phred33"

    def _pall(self, pipe_input, max_nreads=400):
        """
        :param max_nreads: How many reads are read to determine the format, default 100.
        :return: phred33 0(phred64 solexa64)
        """
        valid = []
        for i, line in enumerate(itertools.islice(pipe_input, 3, max_nreads, 4)):
            line = str(line.rstrip(), "utf-8") if isinstance(line, bytes) else line.rstrip()
            qual_val_counts = collections.Counter(ord(qual_char) for qual_char in line)
            min_base_qual = min(qual_val_counts.keys())
            max_base_qual = max(qual_val_counts.keys())
            if min_base_qual < self.gmax or max_base_qual > self.gmax:
                self.gmin, self.gmax = min(min_base_qual, self.gmin), max(max_base_qual, self.gmax)
                valid = self._get_encodings_in_range(self.gmin, self.gmax)
                valid = self._heuristic_filter(valid, qual_val_counts)
                if len(valid) == 0:
                    # no encodings for range
                    return 0
        phred_result = set([self.RANGES[x][2] for x in valid])
        if len(phred_result) == 1:
            result = phred_result.pop()
            if result == "solexa64":
                return 0
            else:
                return result
        else:
            # 1> no encodings for range
            # 2> its format can't be judged by the input reads, you can try to increase the number of sequences.
            return 0

    def _get_encodings_in_range(self, rmin, rmax):
        valid_encodings = []
        for encoding, (emin, emax, phred) in self.RANGES.items():
            if rmin >= emin and rmax <= emax:
                valid_encodings.append(encoding)
        return valid_encodings

    def _heuristic_filter(self, valid, qual_val_counts):
        """
            Apply heuristics to particular ASCII value scores
           to try to narrow-down the encoding, beyond min/max.
        """
        if 'Illumina-1.5' in valid:
            if qual_val_counts[64] > 0 or qual_val_counts[65] > 0:
                # 66: Phread+64 quality score of 2 'B'
                #     used by Illumina 1.5+ as QC indicator
                valid.remove('Illumina-1.5')
            elif 66 in map(operator.itemgetter(0), qual_val_counts.most_common(self.N_MOST_COMMON_THRESH)):
                # A large number of 'B' quality scores (value 2, ASCII 66)
                # were detected, which makes it likely that this encoding is
                # Illumina-1.5, which has been returned as the only option.
                valid = ['Illumina-1.5']
        return valid