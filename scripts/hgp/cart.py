"""libcint's Cartesian component order, shared by the generators."""


def cart_components(l):
    out = []
    for lx in range(l, -1, -1):
        for ly in range(l - lx, -1, -1):
            out.append((lx, ly, l - lx - ly))
    return out


def cart_index(lx, ly, lz):
    l = lx + ly + lz
    idx = 0
    for k in range(l, lx, -1):
        idx += l - k + 1
    return idx + (l - lx - ly)
