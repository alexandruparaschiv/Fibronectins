""" Responsible for grouping the fibronectin monomers into """


from interaction_matrix import InteractionMatrixBuilder

class FibrilsNetwork:

    def __init__(self,radial_cutoff,nematic_cutoff,df,data,periodic_ckdt_tree, atom_to_monomer,fibronectins):

        matrix = InteractionMatrixBuilder(radial_cutoff,df,data,periodic_ckdt_tree)
        interaction_matrix = matrix.get_interactions(data,periodic_ckdt_tree,atom_to_monomer,df)
        self.clusters = InteractionMatrixBuilder.clustering(interaction_matrix)
        self.aggregates = list(filter(lambda x: len(x)>1, self.clusters))
        self.network = []
        self.nematic_cutoff = nematic_cutoff
        self.generate_fibril_network(fibronectins)

    def generate_fibril_network(self,fibronectins):

        for aggregate in self.aggregates:
            alignment_matrix = [[0]*len(aggregate) for _ in range(len(aggregate))]
            for i in range(len(aggregate)):
                for j in range(len(aggregate)):
                    alignment = InteractionMatrixBuilder.get_alignment(aggregate[i],aggregate[j],fibronectins)
                    if alignment > self.nematic_cutoff:
                         alignment_matrix[i][j],alignment_matrix[j][i] = 1,1
            alignment_matrix = InteractionMatrixBuilder.clustering(alignment_matrix)

            fibrils = []

            for item in alignment_matrix:
                fibril = []
                for fibril_indices in item:
                    fibril.append(aggregate[fibril_indices])
                fibrils.append(fibril[:])
            self.network.extend(fibrils)
