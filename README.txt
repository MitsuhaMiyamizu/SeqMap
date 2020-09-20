Compile:
g++ -O3 -m64 -o seqmap match.cpp

Run:
./seqmap 2 probes.fa trans.fa output.txt [options]

For Eland output format: (the default output format is /eland:2)
./seqmap 2 probes.fa trans.fa output.txt /eland  (or /eland:n, where n=1,2 or 3)

Output format for all matches (use option /output_all_matches)
TranscriptId  TranscriptCoord  TranscriptSeq  ProbeId  ProbeSeq  MismatchNumber [Strand]

Mostly used options:
use "./seqmap 3 probes.fa trans.fa output.txt" to search for 3bp mismatches.
use option "/cut:1,25" to take the first 25mer of the probes for the mapping (default is to map the full-length probes).
use option "/allow_insdel:1" to allow 1bp insertion/deletion in the mapping (default is disallowed).
use option "/forward_strand" to map to forward strand only (default is to map to both strands).

Visit "http://biogibbs.stanford.edu/~jiangh/SeqMap/" for documents and updates.
Visit "http://biogibbs.stanford.edu/~jiangh/SeqMap/FAQ.html" for a quick start.
