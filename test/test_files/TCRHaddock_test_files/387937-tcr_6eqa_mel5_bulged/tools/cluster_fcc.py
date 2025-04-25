#!/usr/bin/env python
import argparse
import sys
from time import time, ctime

"""
Asymmetric Taylor-Butina Disjoint Clustering Algorithm.

Authors:
           RODRIGUES Joao
           TRELLET Mikael
           MELQUIOND Adrien
"""


class Element(object):
    """
    Defines a 'clusterable' Element
    """
    __slots__ = ["name", "cluster", "neighbors"]

    def __init__(self, name) -> None:
        """

        Args:
            name:
        """
        self.name = name
        self.cluster = 0
        self.neighbors = set()

    def add_neighbor(self, neighbor) -> None:
        """
        Adds another element to the neighbor list.

        Args:
            neighbor:

        Returns:

        """
        self.neighbors.add(neighbor)

    def assign_cluster(self, clust_id) -> None:
        """Assigns the Element to Cluster. 0 if unclustered"""
        self.cluster = clust_id


class Cluster(object):
    """
    Defines a Cluster. A Cluster is created with a name and a center (Element class)
    """
    __slots__ = ["name", "center", "members"]

    def __init__(self, name, center) -> None:
        """

        Args:
            name:
            center:
        """
        self.name = name
        self.center = center
        self.members = []
        self.populate()

    def __len__(self) -> int:
        """

        Returns:

        """
        return len(self.members)+1  # +1 Center

    def populate(self) -> None:
        """
        Populates the Cluster member list through the neighbor list of its center.

        Returns:

        """
        name = self.name
        # Assign center
        ctr = self.center
        ctr.assign_cluster(name)

        mlist = self.members
        # Assign members
        ctr_nlist = (n for n in ctr.neighbors if not n.cluster)
        for e in ctr_nlist:
            mlist.append(e)
            e.assign_cluster(name)

    def add_member(self, element) -> None:
        """
        Adds one single element to the cluster.

        Args:
            element:

        Returns:

        """
        self.members.append(element)
        element.assign_cluster(self.name)


def read_matrix(path: str, cutoff, strictness) -> dict:
    """
    Reads in a four column matrix (1 2 0.123 0.456\n) and creates an dictionary of Elements.

    The strictness factor is a <float> that multiplies by the cutoff to produce
    a new cutoff for the second half of the matrix. Used to allow some variability
    while keeping very small interfaces from clustering with anything remotely similar.

    Args:
        path:
        cutoff:
        strictness:

    Returns:

    """
    cutoff = float(cutoff)
    partner_cutoff = float(cutoff) * float(strictness)

    elements = dict()

    f = open(path, "r")
    for line in f:
        ref, mobi, d_rm, d_mr = line.split()
        ref = int(ref)
        mobi = int(mobi)
        d_rm = float(d_rm)
        d_mr = float(d_mr)

        # Create or Retrieve Elements
        if ref not in elements:
            r = Element(ref)
            elements[ref] = r
        else:
            r = elements[ref]

        if mobi not in elements:
            m = Element(mobi)
            elements[mobi] = m
        else:
            m = elements[mobi]

        # Assign neighbors
        if d_rm >= cutoff and d_mr >= partner_cutoff:
            r.add_neighbor(m)
        if d_mr >= cutoff and d_rm >= partner_cutoff:
            m.add_neighbor(r)

    f.close()

    return elements


def remove_true_singletons(element_pool) -> tuple:
    """
    Removes from the pool elements without any neighbor.

    Args:
        element_pool:

    Returns:

    """
    ep = element_pool

    ts = set([e for e in ep if not ep[e].neighbors])

    # Remove ts from everybody's neighbor list
    ts_e = set(ep[e] for e in ts)
    for e in element_pool:
        ep[e].neighbors = ep[e].neighbors.difference(ts_e)

    # Remove ts from pool
    for e in ts:
        del ep[e]

    return ts, ep


def cluster_elements(element_pool, threshold) -> tuple:
    """
    Groups Elements within a given threshold together in the same cluster.
    Args:
        element_pool:
        threshold:

    Returns:

    """
    clusters = []
    threshold -= 1  # Account for center
    ep = element_pool
    cn = 1  # Cluster Number
    while 1:
        # Clusterable elements
        ce = [e for e in ep if not ep[e].cluster]
        if not ce:  # No more elements to cluster
            break

        # Select Cluster Center
        # Element with largest neighbor list
        ctr_nlist, ctr = sorted(
            [(len([se for se in ep[e].neighbors if not se.cluster]), e) for e in ce]
        )[-1]

        # Cluster until length of remaining elements lists are above threshold
        if ctr_nlist < threshold:
            break

        # Create Cluster
        c = Cluster(cn, ep[ctr])
        cn += 1
        clusters.append(c)

    return ep, clusters


def output_clusters(handle, clusters) -> None:
    """
    Outputs the cluster name, center, and members.

    Args:
        handle:
        clusters:

    Returns:

    """
    write = handle.write
    for c in clusters:
        write(f"Cluster {c.name} -> {c.center.name} ")
        for m in sorted(c.members, key=lambda k: k.name):
            write(f"{m.name} ")
        write("\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("matrix_file", type=str, help="Matrix file")
    parser.add_argument("threshold", type=float, help="Threshold")
    parser.add_argument(
        "-o", "--output", dest="output_handle", action="store", type=str, default=sys.stdout,
        help="Output File [STDOUT]"
    )
    parser.add_argument(
        "-c", "--cluster-size", dest="clus_size", action="store", type=int, default=4,
        help="Minimum number of elements in a cluster [Default: 4]"
    )
    parser.add_argument(
        "-s", "--strictness", dest="strictness", action="store", type=float, default=0.75,
        help="Multiplier for cutoff for M->R inclusion threshold. [0.75 = eff. cutoff of 0.5625]"
    )
    args = parser.parse_args()

    # Read Matrix
    sys.stderr.write(f"+ BEGIN: {ctime()}\n")
    t_init = time()

    try:
        pool = read_matrix(args.matrix_file, args.threshold, args.strictness)
    except IOError:
        sys.stderr.write(f"File not found: {args.matrix_file}\n")
        sys.exit(1)

    time_taken = int(time()-t_init)
    sys.stderr.write(f"+ Read {len(pool)}x{len(pool)} distance matrix in {time_taken} seconds\n")

    # Cluster
    element_pool, clusters = cluster_elements(pool, args.clus_size)

    # Output Clusters
    o = args.output_handle
    if isinstance(o, str):
        o_handle = open(o, "w")
    else:
        o_handle = o

    sys.stderr.write(f"+ Writing {len(clusters)} Clusters\n")
    output_clusters(o_handle, clusters)
    if isinstance(o, str):
        o_handle.close()

    total_elements = len(element_pool)
    clustered = sum([len(c) for c in clusters])

    # Calculate coverage
    clust_coverage = clustered*100/max(float(total_elements), 1.0)
    sys.stderr.write(f"+ Coverage {clust_coverage:3.2f}% ({clustered}/{total_elements})\n")
    t_elapsed = time()-t_init
    sys.stderr.write(f"+ END: {ctime()} [{t_elapsed:3.2f} seconds]\n")


if __name__ == "__main__":
    main()
