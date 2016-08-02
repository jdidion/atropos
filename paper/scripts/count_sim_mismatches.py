import sys
import editdistance
f=open(sys.argv[1],'rt')
while True:
  line = next(f)
  if line.startswith('#') or line.startswith('@'): continue
  break
next(f)
next(f)
total = 0
num_adapters = 0
num_mismatch = 0
total_len = 0
num_snps = 0
for line in f:
  r1=next(f).rstrip()
  r2=next(f).rstrip()
  total += 1
  total_len += len(r1)
  if r1 == r2: continue
  if len(r1) != len(r2):
    num_adapters += 1
    r2 = r2[0:len(r1)]
  if r1 != r2:
    num_mismatch += 1
    num_snps += editdistance.eval(r1,r2)
print("{}/{} adapter".format(num_adapters, total))
print("{}/{} mismatch".format(num_mismatch, total))
print("{}/{} bp mismatch".format(num_snps, total_len))
