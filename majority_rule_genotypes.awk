#!/usr/bin/env awk
# Author: Nathan S. Watson-Haigh (Bioinformatics Hub, University of Adelaide)
#
# Usage
#  awk -v a_from=6 -v a_to=8 -v b_from=9 -v b_to=11 -f majority_rule_genotypes.awk < t/test.tsv
#
# Summary
# This script can be used to summarise, by majority rule, the genotype calls for two groups of columns.
# One such use case, is summarising the calls for parent A and parent B, each contained in multiple columns.
#
# Expected input
# This script expects 4 awk variables to be set by the calling script:
#  a_from=<int>
#  a_to=<int>
#  b_from=<int>
#  b_to=<int>
# These variables define which column numbers (1-based) contain the genotype call information for the 2 groups of columns (e.g. parent A and parent B)
#
# Limitations:
#  The columns for each group must be consecutive, so they can be defined as a range of columns defined by the start and end column numbers.
#  The columns for group A must be before the columns for group B
#  The group B columns must come straight after the columns for group A

BEGIN{IFS="\t"; FS="\t";
  num_a_cols=a_to-a_from+1;
  num_b_cols=b_to-b_from+1;
}
NR<=1 {
  # munge header line
  for(i=1; i<NF; i++) {
    if(i<=a_from || i>=b_to) {
      printf "%s\t", $i
    }
  }
  print $NF
}
NR>1 {
  # Ensure the default majority rule values are set at the start of each non-header line
  a_maj="N"
  b_maj="N"
  # Ensure we start with empty arrays for tracking the number/type of calls seen for each group
  delete a_calls
  delete b_calls

  for(i=1; i<=NF; i++) {
    if(i>=a_from && i<=a_to) {
      # We are in one of the group A columns, store the call
      a_calls[$i]++

      if(i==a_to) {
        # got to last col of group A, so calculate and output the majority rule
        for(a_call in a_calls){
          #printf "*%s-%s*\t", a_call, a_calls[a_call]
          if((a_calls[a_call]/num_a_cols) > 0.5){
            a_maj=a_call
            break
          }
        }
        printf "%s\t", a_maj
        continue
      }
      continue
    } else if(i>=b_from && i<=b_to) {
      # We are in one of the group B columns, store the call
      b_calls[$i]++

      if(i==b_to) {
        # got to last col of group B, so calculate and output the majority rule
        for(b_call in b_calls){
          #printf "*%s-%s*\t", b_call, b_calls[b_call]
          if((b_calls[b_call]/num_b_cols) > 0.5){
            b_maj=b_call
            break
          }
        }
        printf "%s\t", b_maj
        continue
      }
      continue
    }

    # Print out all other columns as is
    if(i==NF) {
      printf "%s\n", $i
    } else {
      printf "%s\t", $i
    }
  }
}
