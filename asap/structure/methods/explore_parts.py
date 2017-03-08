from itertools import product


def get_shared_edges(part1, part2):
    """
    Get the shared edges of the two parts.
    """
    edges1, edges2 = part1.edges, part2.edges
    if not edges1 or not edges2:
        return []

    shared_edges = []
    for e1, e2 in product(edges1, edges2):
        if e1.IsSame(e2):
            unique = True
            for ei in shared_edges:
                if ei.IsSame(e1):
                    unique = False
                    break
            if unique:
                shared_edges.append(e1)

    return shared_edges


def get_shared_nodes(part1, part2):
    """
    Get the shared nodes of the two parts.
    """
    nodes1 = part1.nodes
    nodes2 = part2.nodes
    return list(nodes1 & nodes2)
