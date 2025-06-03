def parse_blast_results(file):
    hits = {}
    with open(file) as f:
        for line in f:
            query, subject, *rest = line.strip().split("\t")
            if query not in hits:
                hits[query] = subject
    return hits

til_to_AT = parse_blast_results("til_to_AT.txt")
AT_to_til = parse_blast_results("AT_to_til.txt")

reciprocal_hits = []
for til, AT in til_to_AT.items():
    if AT in AT_to_til and AT_to_til[AT] == til:
        reciprocal_hits.append((til, AT))

# Write reciprocal hits to a file
with open("recip_til_AT.txt", "w") as f:
    for til, AT in reciprocal_hits:
        f.write(f"{til}\t{AT}\n")

