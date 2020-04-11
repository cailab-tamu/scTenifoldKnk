import numpy as np


class gene(object):

    def __init__(self, geneID, geneType, binID = -1):

        """
        geneType: 'MR' master regulator or 'T' target
        bindID is optional
        """

        self.ID = geneID
        self.Type = geneType
        self.binID = binID
        self.Conc = []
        self.Conc_S = []
        self.dConc = []
        self.k = [] #For dynamics simulation it stores k1 to k4 for Rung-Kutta method, list of size 4 * num_c_to_evolve
        self.k_S = [] #For dynamics simulation it stores k1 to k4 for Rung-Kutta method, list of size 4 * num_c_to_evolve
        self.simulatedSteps_ = 0
        self.converged_ = False
        self.converged_S_ = False
        self.ss_U_ = 0 #This is the steady state concentration of Unspliced mRNA
        self.ss_S_ = 0 #This is the steady state concentration of Spliced mRNA

    def append_Conc (self, currConc):
        if currConc < 0:
            currConc = 0
        self.Conc.append(currConc)

    def append_Conc_S (self, currConc):
        if currConc < 0:
            currConc = 0
        self.Conc_S.append(currConc)

    def append_dConc (self, currdConc):
        self.dConc.append(currdConc)

    def append_k (self, list_currK):
        self.k.append(list_currK)

    def append_k_S (self, list_currK):
        self.k_S.append(list_currK)

    def del_lastK_Conc(self, K):
        for k in range(K):
            self.Conc.pop(-1)

    def del_lastK_Conc_S(self, K):
        for k in range(K):
            self.Conc_S.pop(-1)

    def clear_Conc (self):
        """
        This method clears all the concentrations except the last one that may
        serve as intial condition for rest of the simulations
        """
        self.Conc = self.Conc[-1:]

    def clear_dConc (self):
        self.dConc = []

    def incrementStep (self):
        self.simulatedSteps_ += 1

    def setConverged (self):
        self.converged_ = True

    def setConverged_S (self):
        self.converged_S_ = True

    def set_scExpression(self, list_indices):
        """
        selects input indices from self.Conc and form sc Expression
        """
        self.scExpression = np.array(self.Conc)[list_indices]

    def set_ss_conc_U(self, u_ss):
        if u_ss < 0:
            u_ss = 0

        self.ss_U_ = u_ss

    def set_ss_conc_S(self, s_ss):
        if s_ss < 0:
            s_ss = 0

        self.ss_S_ = s_ss

    def clear_k(self):
        self.k = []

    def clear_k_S(self):
        self.k_S = []
