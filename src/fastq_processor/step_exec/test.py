seq = "ACG CTG TTA TCC CTA AAG T"
new_seq = ""
for s in reversed(seq.replace(" ", "")):
    new_seq += s
print(new_seq)