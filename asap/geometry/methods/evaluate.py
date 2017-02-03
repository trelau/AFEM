from numpy import array, floor


def find_span(n, p, u, uk):
    """
    Determine the knot span index.

    :param int n: Number of control points - 1.
    :param int p: Degree.
    :param float u: Parameter.
    :param ndarray uk: Knot vector.

    :return: Knot span.
    :rtype: int

    *Reference:* Algorithm A2.1 from "The NURBS Book".
    """
    # Special case
    if u >= uk[n + 1]:
        return n
    if u <= uk[p]:
        return p

    # Do binary search
    low = p
    high = n + 1
    mid = int(floor((low + high) / 2.))
    while u < uk[mid] or u >= uk[mid + 1]:
        if u < uk[mid]:
            high = mid
        else:
            low = mid
        mid = int(floor((low + high) / 2.))
    return mid


def basis_funs(i, u, p, uk):
    """
    Compute the non-vanishing basis functions.

    :param int i: Knot span index.
    :param float u: Parameter.
    :param int p: Degree.
    :param ndarray uk: Knot vector.

    :return: Non-vanishing basis functions.
    :rtype: ndarray

    Reference: Algorithm A2.2 from "The NURBS Book"
    """
    bf = [0.0] * (p + 1)
    bf[0] = 1.0
    left = [0.0] * (p + 1)
    right = [0.0] * (p + 1)
    for j in range(1, p + 1):
        left[j] = u - uk[i + 1 - j]
        right[j] = uk[i + j] - u
        saved = 0.0
        for r in range(0, j):
            temp = bf[r] / (right[r + 1] + left[j - r])
            bf[r] = saved + right[r + 1] * temp
            saved = left[j - r] * temp
        bf[j] = saved
    return array(bf, dtype=float)
