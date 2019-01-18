# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import json
import copy


class DataProcess4():
    clusterHierData = []

    class ClusterNode:
        """
        A tree node class for representing a cluster.
        Leaf nodes correspond to original observations, while non-leaf nodes
        correspond to non-singleton clusters.
        The to_tree function converts a matrix returned by the linkage
        function into an easy-to-use tree representation.
        See Also
        --------
        to_tree : for converting a linkage matrix ``Z`` into a tree object.
        """

        def __init__(self, id, left=None, right=None, dist=0, count=1):
            if id < 0:
                raise ValueError('The id must be non-negative.')
            if dist < 0:
                raise ValueError('The distance must be non-negative.')
            if (left is None and right is not None) or \
                    (left is not None and right is None):
                raise ValueError('Only full or proper binary trees are permitted.'
                                 '  This node has one child.')
            if count < 1:
                raise ValueError('A cluster must contain at least one original '
                                 'observation.')
            self.id = id
            self.left = left
            self.right = right
            self.dist = dist
            if self.left is None:
                self.count = count
            else:
                self.count = left.count + right.count

        def __lt__(self, node):
            if not isinstance(node, ClusterNode):
                raise ValueError("Can't compare ClusterNode "
                                 "to type {}".format(type(node)))
            return self.dist < node.dist

        def __gt__(self, node):
            if not isinstance(node, ClusterNode):
                raise ValueError("Can't compare ClusterNode "
                                 "to type {}".format(type(node)))
            return self.dist > node.dist

        def __eq__(self, node):
            if not isinstance(node, ClusterNode):
                raise ValueError("Can't compare ClusterNode "
                                 "to type {}".format(type(node)))
            return self.dist == node.dist

        def get_id(self):
            """
            The identifier of the target node.
            For ``0 <= i < n``, `i` corresponds to original observation i.
            For ``n <= i < 2n-1``, `i` corresponds to non-singleton cluster formed
            at iteration ``i-n``.
            Returns
            -------
            id : int
                The identifier of the target node.
            """
            return self.id

        def get_count(self):
            """
            The number of leaf nodes (original observations) belonging to
            the cluster node nd. If the target node is a leaf, 1 is
            returned.
            Returns
            -------
            get_count : int
                The number of leaf nodes below the target node.
            """
            return self.count

        def get_left(self):
            """
            Return a reference to the left child tree object.
            Returns
            -------
            left : ClusterNode
                The left child of the target node.  If the node is a leaf,
                None is returned.
            """
            return self.left

        def get_right(self):
            """
            Returns a reference to the right child tree object.
            Returns
            -------
            right : ClusterNode
                The left child of the target node.  If the node is a leaf,
                None is returned.
            """
            return self.right

        def is_leaf(self):
            """
            Returns True if the target node is a leaf.
            Returns
            -------
            leafness : bool
                True if the target node is a leaf node.
            """
            return self.left is None

        def pre_order(self, func=(lambda x: x.id)):
            """
            Performs pre-order traversal without recursive function calls.
            When a leaf node is first encountered, ``func`` is called with
            the leaf node as its argument, and its result is appended to
            the list.
            For example, the statement::
               ids = root.pre_order(lambda x: x.id)
            returns a list of the node ids corresponding to the leaf nodes
            of the tree as they appear from left to right.
            Parameters
            ----------
            func : function
                Applied to each leaf ClusterNode object in the pre-order traversal.
                Given the i'th leaf node in the pre-ordeR traversal ``n[i]``, the
                result of func(n[i]) is stored in L[i]. If not provided, the index
                of the original observation to which the node corresponds is used.
            Returns
            -------
            L : list
                The pre-order traversal.
            """

            # Do a preorder traversal, caching the result. To avoid having to do
            # recursion, we'll store the previous index we've visited in a vector.
            n = self.count

            curNode = [None] * (2 * n)
            lvisited = set()
            rvisited = set()
            curNode[0] = self
            k = 0
            preorder = []
            while k >= 0:
                nd = curNode[k]
                ndid = nd.id
                if nd.is_leaf():
                    preorder.append(func(nd))
                    k = k - 1
                else:
                    if ndid not in lvisited:
                        curNode[k + 1] = nd.left
                        lvisited.add(ndid)
                        k = k + 1
                    elif ndid not in rvisited:
                        curNode[k + 1] = nd.right
                        rvisited.add(ndid)
                        k = k + 1
                    # If we've visited the left and right of this non-leaf
                    # node already, go up in the tree.
                    else:
                        k = k - 1

            return preorder

    def to_tree(self, Z, rd=False):
        """
        Converts a hierarchical clustering encoded in the matrix ``Z`` (by
        linkage) into an easy-to-use tree object.
        The reference r to the root ClusterNode object is returned.
        Each ClusterNode object has a left, right, dist, id, and count
        attribute. The left and right attributes point to ClusterNode objects
        that were combined to generate the cluster. If both are None then
        the ClusterNode object is a leaf node, its count must be 1, and its
        distance is meaningless but set to 0.
        Note: This function is provided for the convenience of the library
        user. ClusterNodes are not used as input to any of the functions in this
        library.
        Parameters
        ----------
        Z : ndarray
            The linkage matrix in proper form (see the ``linkage``
            function documentation).
        rd : bool, optional
            When False, a reference to the root ClusterNode object is
            returned.  Otherwise, a tuple (r,d) is returned. ``r`` is a
            reference to the root node while ``d`` is a dictionary
            mapping cluster ids to ClusterNode references. If a cluster id is
            less than n, then it corresponds to a singleton cluster
            (leaf node). See ``linkage`` for more information on the
            assignment of cluster ids to clusters.
        Returns
        -------
        L : list
            The pre-order traversal.
        """

        Z = np.asarray(Z, order='c')


        # Number of original objects is equal to the number of rows minus 1.
        n = Z.shape[0] + 1

        # Create a list full of None's to store the node objects
        d = [None] * (n * 2 - 1)

        # Create the nodes corresponding to the n original objects.
        for i in xrange(0, n):
            d[i] = DataProcess4.ClusterNode(i)

        nd = None

        for i in xrange(0, n - 1):
            fi = int(Z[i, 0])
            fj = int(Z[i, 1])
            # print str(fi) + " " + str(fj)

            if fi > i + n:
                raise ValueError(('Corrupt matrix Z. Index to derivative cluster '
                                  'is used before it is formed. See row %d, '
                                  'column 0') % fi)
            if fj > i + n:
                raise ValueError(('Corrupt matrix Z. Index to derivative cluster '
                                  'is used before it is formed. See row %d, '
                                  'column 1') % fj)
            nd = DataProcess4.ClusterNode(i + n, d[fi], d[fj], Z[i, 2])
            # fw.write("["+str(d[fi].id) + "," + str(d[fj].id) + "," + str(n+i)+"],\n")
            #          ^ id   ^ left ^ right ^ dist
            # if Z[i, 3] != nd.count:
            #     raise ValueError(('Corrupt matrix Z. The count Z[%d,3] is '
            #                       'incorrect.') % i)
            d[n + i] = nd

        if rd:
            return (nd, d)
        else:
            return nd

    def condensed_index(self, n, i, j):
        """
        Calculate the condensed index of element (i, j) in an n x n condensed
        matrix.
        """
        if i < j:
            return n * i - (i * (i + 1) / 2) + (j - i - 1)
        elif i > j:
            return n * j - (j * (j + 1) / 2) + (i - j - 1)

    def pdist(self, Z, n):
        D = np.ndarray(n * (n - 1) / 2, dtype=np.double)
        for i in range(n):
            for j in range(n):
                if i < j:
                    D[n * i - (i * (i + 1) / 2) + (j - i - 1)] = Z[i][j]
        return D

    def linkage(self, dists, dists1, Z, n, method):
        """
        Perform hierarchy clustering.
        Parameters
        ----------
        dists1 : ndarray
            A condensed matrix stores the pairwise distances of the observations.
        dists : ndarray
        A condensed matrix stores the pairwise distances(in percentage) of the observations.
        Z : ndarray
            A (n - 1) x 4 matrix to store the result (i.e. the linkage matrix).
        n : int
            The number of observations.
        method : int
            The linkage method. 0: single 1: complete 2: average 3: centroid
            4: median 5: ward 6: weighted
        """
        # inter-cluster dists
        D = np.ndarray(n * (n - 1) / 2, dtype=np.int)
        D_num = np.ndarray(n * (n - 1) / 2, dtype=np.int)
        # map the indices to node ids
        id_map = np.ndarray(n, dtype=np.int)

        D[:] = dists
        D_num[:] = dists1
        for i in range(n):
            id_map[i] = i

        for k in range(n - 1):
            # find two closest clusters x, y (x < y)
            percent_current_max = 0
            number_current_max = 0

            # get the max value of interconnection between different nodes/clusters
            for i in range(n - 1):
                # if the node is dropped, skip this node
                if id_map[i] == -1:
                    continue

                # get the max value of interconnection between id_map[i] node and other nodes/clusters
                i_start = self.condensed_index(n, i, i + 1)
                for j in range(n - i - 1):
                    if D[i_start + j] >= percent_current_max:
                        percent_current_max = D[i_start + j]
                        number_current_max = D_num[i_start + j]
                        x = i
                        y = i + j + 1

            # get the real index for the two nodes with maximun value
            id_x = id_map[x]
            id_y = id_map[y]

            # update the index for two nodes
            id_map[x] = -1  # cluster x will be dropped
            id_map[y] = n + k  # cluster y will be replaced with the new cluster

            # record the new node
            Z[k, 0] = min(id_x, id_y)
            Z[k, 1] = max(id_y, id_x)

            Z[k, 2] = percent_current_max
            Z[k, 3] = number_current_max

            # update the distance matrix
            self.updateClusterGenes(x, y)

            for i in range(n):
                id_i = id_map[i]
                if id_i == -1 or id_i == n + k:
                    continue

                # ni = 1 if id_i < n else <int>Z[id_i - n, 3]

                #D[self.condensed_index(n, i, y)] = max(
                    #D[self.condensed_index(n, i, x)],
                    #D[self.condensed_index(n, i, y)])
                #D_num[self.condensed_index(n, i, y)] = max(
                    #D_num[self.condensed_index(n, i, x)],
                    #D_num[self.condensed_index(n, i, y)])

                # calcuate the real distance
                #self.updateClusterGenes(x, y)
                D[self.condensed_index(n, i, y)] = self.calRealDis(i,x,y)
                D_num[self.condensed_index(n, i, y)] = self.calRealDisNum(i,x,y)
                if i < x:
                    D[self.condensed_index(n, i, x)] = -1
                    D_num[self.condensed_index(n, i, x)] = -1
        for t in range(n - 1):
            level = []
            level.append(int(Z[t, 0]))
            level.append(int(Z[t, 1]))
            level.append(t + n)
            level.append(int(Z[t, 3]))
            level.append(int(Z[t, 2]))
            self.clusterHierData.append(level)

        return Z

    def calRealDis(self, i, x, y):
        """
        calculate the real distance between different clusters based on the number of different genes.
        Instad of using the max value of two clusters/node for approxiamte estimation.
        """
        iGenes = self.getGenes(i)
        yGenes = self.getGenes(y)

        totalGenes = self.getTotalGenes(iGenes, yGenes)

        numOfIntersectedGene = self.getNumOfIntersectedGenes(i, yGenes)
        percentageOfOverlappingGene = numOfIntersectedGene*100.0/len(totalGenes)


        #return numOfIntersectedGene
        return percentageOfOverlappingGene
    def calRealDisNum(self, i, x, y):
        """
        calculate the real distance between different clusters based on the number of different genes.
        Instad of using the max value of two clusters/node for approxiamte estimation.
        """
        xGenes = self.getGenes(x)
        yGenes = self.getGenes(y)
        totalGenes = self.getTotalGenes(xGenes, yGenes)
        numOfIntersectedGene = self.getNumOfIntersectedGenes(i, totalGenes)

        return numOfIntersectedGene
    def getGenes(self, index):

        return self.go_info[index]["genes"].split(";")

    def getTotalGenes(self, xGenes, yGenes):

        totalGene = []

        for gene in xGenes:
            totalGene.append(gene)

        for gene in yGenes:
            if gene not in totalGene:
                totalGene.append(gene)

        return totalGene

    def getNumOfIntersectedGenes(self, index_i, totalGene):
        numberOfCommonGene = 0
        genes = self.getGenes(index_i)

        for gene in genes:
            if gene in totalGene:
                numberOfCommonGene += 1
        return numberOfCommonGene

    def updateClusterGenes(self, x,y):
        xGenes = self.getGenes(x)
        yGenes = self.getGenes(y)
        for gene in xGenes:
            if gene not in yGenes:
                yGenes.append(gene)

        self.go_info[y]["genes"] = ";".join(yGenes)

    def getGODependency(self, GO_inf):

        def recuriveGetGOId(GO_id):
            if GO_id in GO_hier:
                GO_hier_list[GO_id] = GO_hier[GO_id]

                for i in GO_hier[GO_id]["p"]:
                    if not GO_hier_list.has_key(i):
                        recuriveGetGOId(i.encode('ascii', 'ignore'))

        remote_server = False;

        if (remote_server):
            root_dir = "/home/ubuntu"
        else:
            root_dir = ""

        with open(root_dir + 'js/GO.js', "r") as fr_GO:
            for GO in fr_GO:
                GO_hier = json.loads(str(GO))

        GO_hier_list = {}

        for gos in GO_inf:
            recuriveGetGOId(gos["GO_id"].encode('ascii', 'ignore'))

        return json.dumps(GO_hier_list)

    def reOrder(self, gen_anno_reord, go_inf):
        go_inf_tmp = []
        for i in gen_anno_reord:
            go_inf_tmp.append(go_inf[i])
        return go_inf_tmp

    def createMatrix(self, go_inf):

        # logger.debug(go_inf)

        size = len(go_inf)
        matrix = [[0] * size for _ in range(size)];
        percent_matrix = [[0] * size for _ in range(size)];

        for i in range(0, size):
            for j in range(0, size):
                if i == j:
                    matrix[i][j] = 0
                    percent_matrix[i][j] = 0
                else:
                    genes_i = go_inf[i]["genes"].split(";")
                    genes_j = go_inf[j]["genes"].split(";")
                    overlappingGenes = set(genes_i).intersection(genes_j)
                    totalGenes = self.getTotalGenes(genes_i, genes_j)
                    matrix[i][j] = len(overlappingGenes)
                    # get count
                    percent_matrix[i][j] = len(overlappingGenes) * 100.0 / len(totalGenes)
                    # get count

        return matrix, percent_matrix

    def createMatrixReord(self, go_inf):
        size = len(go_inf)

        matrix_count = [[0] * size for _ in range(size)];
        matrix_geneName = {}

        for i in range(0, size):
            for j in range(0, size):
                if i == j:
                    matrix_count[i][j] = 0
                    matrix_geneName[str(i) + "-" + str(j)] = ""
                else:

                    genes_i = go_inf[i]["genes"].split(";")
                    genes_j = go_inf[j]["genes"].split(";")
                    overlappingGenes = set(genes_i).intersection(genes_j)
                    totalGenes = self.getTotalGenes(genes_i, genes_j)

                    geneStr = ""
                    for k in overlappingGenes:
                        geneStr += k + ";"

                    # get gene name
                    matrix_geneName[str(i) + "-" + str(j)] = geneStr[:-1]

                    # get count
                    matrix_count[i][j] = len(overlappingGenes)*100.0/len(totalGenes)
        #print matrix_count
        return {"matrix_count": matrix_count}

    def dataProcess(self, go_inf):

        self.clusterHierData = []  ##clear clusterHierData

        self.go_info = copy.deepcopy(go_inf)

        size = len(go_inf)

        if size == 0:
            raise Exception("go_inf is empty")

        num_matrix, percent_matrix = self.createMatrix(go_inf)

        D = np.ndarray(size * (size - 1) / 2, dtype=np.int)
        # D = self.pdist(matrix,size)
        D = self.pdist(percent_matrix, size)
        D_num = self.pdist(num_matrix, size)
        # Z = np.zeros((size - 1, 3))
        Z = np.zeros((size - 1, 4))
        Z = self.linkage(D, D_num, Z, size, 0)

        nd = [None] * (size * 2 - 1)
        nd = self.to_tree(Z)

        go_index_reord = nd.pre_order()  # array_order

        go_inf_reOrder = self.reOrder(go_index_reord, go_inf)

        matrix_reOrder = self.createMatrixReord(go_inf_reOrder)  # create matrix and clusterHierData
        #print matrix_reOrder
        go_hier = self.getGODependency(go_inf_reOrder)

        #print self.clusterHierData
        return {"matrix": matrix_reOrder, "go_index_reord": go_index_reord, "clusterHierData": self.clusterHierData,
                "go_inf": go_inf_reOrder, "go_hier": go_hier}