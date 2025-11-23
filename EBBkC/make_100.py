import random

random.seed(42)
N = 100
communities = [
    list(range(0, 20)),    # Community 1: nodes 0-19
    list(range(20, 40)),   # Community 2: nodes 20-39
    list(range(40, 60)),   # Community 3: nodes 40-59
    list(range(60, 80)),   # Community 4: nodes 60-79
    list(range(80, 100)),  # Community 5: nodes 80-99
]

edges = set()

# Add all edges within each community (clique)
for comm in communities:
    for i in range(len(comm)):
        for j in range(i+1, len(comm)):
            edges.add((comm[i], comm[j]))

# Add some random inter-community edges
for _ in range(200):
    c1, c2 = random.sample(communities, 2)
    u = random.choice(c1)
    v = random.choice(c2)
    if u != v:
        edges.add(tuple(sorted((u, v))))

# Write to file
with open("facebook_100.edges", "w") as f:
    for u, v in sorted(edges):
        f.write(f"{u} {v}\n")
