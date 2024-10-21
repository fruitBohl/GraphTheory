import networkx as nx
import matplotlib.pyplot as plt


def _matrix_mult_mod_q(a: tuple[int], b: tuple[int], q: int) -> tuple[int]:
    a1, a2, a3, a4 = a
    b1, b2, b3, b4 = b

    return (
        (a1 * b1 + a2 * b3) % q,
        (a1 * b2 + a2 * b4) % q,
        (a3 * b1 + a4 * b3) % q,
        (a3 * b2 + a4 * b4) % q,
    )


def _reduce_mod_q(sol: tuple[int], q: int) -> tuple[int]:
    return tuple([x % q for x in sol])


def _construct_matrix(sol: tuple[int], x: int, y: int) -> tuple[int]:
    a, b, c, d = sol

    return tuple(
        [
            (a + b * x + d * y) % q,
            (-b * y + c + d * x) % q,
            (-b * y - c + d * x) % q,
            (a - b * x - d * y) % q,
        ]
    )


def _generate_GL2(q: int) -> list[tuple[int]]:
    matrices = []

    for a in range(q):
        for b in range(q):
            for c in range(q):
                for d in range(q):
                    if (a * d - b * c) % q == 0:
                        continue
                    else:
                        matrices.append((a, b, c, d))
    return matrices


def _generate_PGL_from_GL(
    gl: list[tuple[int]], q: int
) -> dict[tuple[int], set(tuple[int])]:
    scalar_matrices = [(k, 0, 0, k) for k in range(1, q)]
    pgl = {}

    for matrix in gl:
        left_coset = set()

        for scalar_matrix in scalar_matrices:
            out = _matrix_mult_mod_q(matrix, scalar_matrix, q)
            left_coset.add(out)

        if left_coset not in pgl.values():
            pgl[matrix] = left_coset

    return pgl


def _check_matrices_not_in_same_coset(
    matrices: list[tuple[int]], pgl: dict[tuple[int], set(tuple[int])]
) -> bool:
    for i in range(len(matrices)):
        for j in range(i + 1, len(matrices)):
            if matrices[i] in pgl and matrices[j] in pgl:
                if pgl[matrices[i]] == pgl[matrices[j]]:
                    print(
                        f"Matrices {matrices[i]} and {matrices[j]} are in the same coset."
                    )
                    return False

    print("All matrices are in different cosets")
    return True


def _compute_cayley_graph_edges(
    pgl: dict[tuple[int], set(tuple[int])],
    S: list[tuple[int]],
) -> list[tuple[tuple[int], tuple[int]]]:
    edge_reps = set()
    coset_representitives = list(pgl.keys())

    for s in S:
        for coset_rep_x in coset_representitives:
            maybe_coset_rep = _matrix_mult_mod_q(coset_rep_x, s, q)

            for coset_rep_y in coset_representitives:
                if (maybe_coset_rep in pgl[coset_rep_y]) or (
                    maybe_coset_rep == coset_rep_y
                ):
                    edge_reps.add(frozenset({coset_rep_x, coset_rep_y}))

    return edge_reps


def _get_coset_representitive(matrix: tuple[int]) -> tuple[int]:
    possible_rep = []
    for coset_rep in pgl2q.keys():
        if matrix in pgl2q[coset_rep]:
            possible_rep.append(coset_rep)

    if len(possible_rep) != 1:
        print(f"Matrix {matrix} has more than one coset representitive.")
        breakpoint()
    return possible_rep[0]


def _find_bipartite_sets(edges: list[tuple[tuple[int], tuple[int]]]) -> tuple[set, set]:
    left = set()
    right = set()
    visited = set()

    start_edge = tuple(edges.pop())
    from_vertex = start_edge[0]

    from_set_left = True

    left.add(from_vertex)
    visited.add(from_vertex)

    while len(left) + len(right) < len(edges):
        to_vertices = set()

        for edge in edges:
            edge_tup = tuple(edge)

            if from_vertex in edge_tup:
                if from_vertex == edge_tup[0]:
                    if from_set_left:
                        right.add(edge_tup[1])
                        to_vertices.add(edge_tup[1])
                    else:
                        left.add(edge_tup[1])
                        to_vertices.add(edge_tup[1])
                else:
                    if from_set_left:
                        right.add(edge_tup[0])
                        to_vertices.add(edge_tup[0])
                    else:
                        left.add(edge_tup[0])
                        to_vertices.add(edge_tup[0])

        from_vertex = None
        for to_vertex in to_vertices:
            if to_vertex not in visited:
                from_vertex = to_vertex
                break
        if not from_vertex:
            break
        visited.add(from_vertex)
        from_set_left = not from_set_left

    for edge in edges:
        edge_tup = tuple(edge)
        if edge_tup[1] in left:
            right.add(edge_tup[0])
        if edge_tup[1] in right:
            left.add(edge_tup[0])
        if edge_tup[0] in left:
            right.add(edge_tup[1])
        if edge_tup[0] in right:
            left.add(edge_tup[1])

    return left, right


if __name__ == "__main__":
    p = 5
    q = 7

    solutions = [
        (1, 2, 0, 0),
        (1, -2, 0, 0),
        (1, 0, 2, 0),
        (1, 0, -2, 0),
        (1, 0, 0, 2),
        (1, 0, 0, -2),
    ]

    solutions_mod_q = [_reduce_mod_q(sol, q) for sol in solutions]

    x = 2
    y = 3

    matrices = [_construct_matrix(sol, x, y) for sol in solutions_mod_q]

    gl2q = _generate_GL2(q)
    pgl2q = _generate_PGL_from_GL(gl2q, q)
    element_reps = list(pgl2q.keys())

    _check_matrices_not_in_same_coset = _check_matrices_not_in_same_coset(
        matrices, pgl2q
    )

    matrix_coset_reps = [_get_coset_representitive(matrix) for matrix in matrices]

    edges = _compute_cayley_graph_edges(pgl2q, matrix_coset_reps)

    left, right = _find_bipartite_sets(edges)

    # number each coset representitive
    element_rep_to_number = {element_reps[i]: i for i in range(len(element_reps))}

    G = nx.Graph()
    G.add_nodes_from([i for i in range(len(element_reps))])
    G.add_edges_from(
        [
            (
                element_rep_to_number[tuple(edge)[0]],
                element_rep_to_number[tuple(edge)[1]],
            )
            for edge in edges
        ]
    )

    left = [element_rep_to_number[el] for el in left]

    pos = nx.bipartite_layout(G, nodes=left)
    nx.draw_networkx_edges(G, pos)
    plt.show()
