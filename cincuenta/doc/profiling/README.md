# Profiling records

This directory preserves archival performance-profiling artifacts for Cincuenta.
It is not a location for runnable job examples.

## Frontier nb11 Darshan profile

[`frontier-nb11-darshan-summary.txt`](frontier-nb11-darshan-summary.txt) is raw
`darshan-parser` output from a 16-process Frontier run of the nb11 case.  The
profile was collected with Darshan log version 3.41 (library version 3.4.7) for
job 5629047, which ran from 2026-10-06 18:56:48 to 19:05:15 UTC (506.4904 s).
The executable was invoked as:

```text
cincuenta -f input.ain -p 12 -l cincuenta.log
```

The profile records ROMIO settings `romio_no_indep_rw=true;cb_nodes=4` and is
kept as the unmodified parser summary so its detailed I/O counters remain
available for later comparison.

This was cincuenta mpi job cut down to last around 10 minutes.  My
intepretation of the output is that IO was almost insignificant order
of less than 1/100 of the the total time.
