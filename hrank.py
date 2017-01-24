# Finding the next lexicographically larger permutation.
# 1.Find the largest index k such that a[k] < a[k + 1]. If no such index exists, the permutation is the last permutation.
# 2.Find the largest index l greater than k such that a[k] < a[l].
# 3.Swap the value of a[k] with that of a[l].
# 4.Reverse the sequence from a[k + 1] up to and including the final element a[n].