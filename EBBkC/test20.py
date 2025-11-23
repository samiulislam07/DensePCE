edges = set()

# 3-clique: nodes 0,1,2
for i in range(3):
    for j in range(i+1, 3):
        edges.add((i, j))

# 4-clique: nodes 3,4,5,6
for i in range(3, 7):
    for j in range(i+1, 7):
        edges.add((i, j))

# 5-clique: nodes 7,8,9,10,11
for i in range(7, 12):
    for j in range(i+1, 12):
        edges.add((i, j))

# Add some extra edges to connect the cliques and fill up to 20 nodes
extra_edges = [
    (2, 3), (6, 7), (11, 12), (12, 13), (13, 14), (14, 15), (15, 16), (16, 17), (17, 18), (18, 19),
    (0, 19), (5, 15), (10, 18)
]
for e in extra_edges:
    edges.add(tuple(sorted(e)))

# Write to file
with open("test20.edges", "w") as f:
    for u, v in sorted(edges):
        f.write(f"{u} {v}\n")