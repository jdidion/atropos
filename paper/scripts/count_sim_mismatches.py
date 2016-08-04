import sys
import editdistance

if int(sys.argv[2]) == 1:
    adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG"
else:
    adapter = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"

with open(sys.argv[1],'rt') as f:
    while True:
        line = next(f)
        if line.startswith('#') or line.startswith('@'): continue
        break

    total_reads = 0
    num_reads_mismatch = 0
    total_edit_dist = 0
    num_adapters = 0
    total_ref_bp = 0
    total_adapter_bp = 0
    total_adapter_edit_dist = 0
    
    while True:
        if total_reads > 0:
            try:
                next(f)
            except:
                break
        r1 = next(f).rstrip()
        r2 = next(f).rstrip()
        total_reads += 1
        total_ref_bp += len(r1)
        
        if r1 == r2:
            continue
        
        if len(r1) != len(r2):
            num_adapters += 1
            l1 = len(r1)
            n = min(len(r2) - l1, len(adapter))
            total_adapter_bp += n
            total_adapter_edit_dist += editdistance.eval(r2[l1:(l1+n)], adapter[0:n])
            r2 = r2[0:l1]
        
        if r1 != r2:
            num_reads_mismatch += 1
            total_edit_dist += editdistance.eval(r1,r2)
        
        
print("{} / {} ({:0.03%}) reads with adapter".format(num_adapters, total_reads, num_adapters / total_reads))
print("{} / {} ({:0.03%}) reads with mismatch(es)".format(num_reads_mismatch, total_reads, num_reads_mismatch / total_reads))
print("{} / {} ({:0.03%}) ref bp mismatch".format(total_edit_dist, total_ref_bp, total_edit_dist / total_ref_bp))
print("{} / {} ({:0.03%}) adapter bp mismatch".format(total_adapter_edit_dist, total_adapter_bp, total_adapter_edit_dist / total_adapter_bp))
