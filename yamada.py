import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

class JAKSTATModel:
    """JAK-STAT Pathway Model based on Yamada et al. 2003"""

    def __init__(self):
        # Model compartments and constants
        self.compartment_cytoplasm = 1.0
        self.compartment_nucleus = 1.0
        self.const_species_IFN = 10.0

        # Complete reaction rate constants
        self.rates = {
            # Receptor module
            'k1f': 0.1,   'k1b': 0.05,  # R-JAK binding
            'k2f': 0.02,  'k2b': 0.02,  # IFN binding
            'k3f': 0.04,  'k3b': 0.2,   # Dimerization
            'k4': 0.005,                # Activation
            
            # STAT1 module
            'k5f': 0.008, 'k5b': 0.8,   # STAT1 binding
            'k6': 0.4,                  # STAT1 activation
            'k7f': 0.005, 'k7b': 0.5,   # pSTAT1 binding
            'k8f': 0.02,  'k8b': 0.1,   # STAT1 dimerization
            
            # Phosphatase module
            'k9f': 0.001, 'k9b': 0.2,   # SHP2 binding
            'k10': 0.003,               # Dephosphorylation
            'k11f': 0.001, 'k11b': 0.2, # PPX binding
            'k12': 0.003,               # STAT1 dephosphorylation
            'k13f': 0.001, 'k13b': 0.2, # Complex formation
            'k14': 0.003,               # Complex dissociation
            
            # Nuclear module
            'k15': 0.005,               # Nuclear import
            'k16f': 0.02, 'k16b': 0.1,  # Nuclear reactions
            'k17f': 0.001, 'k17b': 0.2, # PPN binding
            'k18': 0.003,               # Nuclear dephosphorylation
            
            # Transcription/Translation
            'k23': 0.01,  'k23b': 400,  # Transcription
            'k24': 0.001,               # mRNA export
            'k25': 0.01,                # Translation
            'k26': 0.0005,              # mRNA degradation
            'k27': 0.0005,              # Protein degradation
            
            # Feedback module
            'k28f': 0.02, 'k28b': 0.1,  # SOCS1 binding
            'k29f': 0.001, 'k29b': 0.2,  # Complex formation
        }

        # Initialize state vector with all species
        self.x0 = np.zeros(33)
        self.x0[0:33] = [
            10.0,   # 0:  R (Receptor)
            10.0,   # 1:  JAK
            0.0,    # 2:  RJ complex
            0.0,    # 3:  IFNRJ complex
            0.0,    # 4:  IFNRJ2 dimer
            0.0,    # 5:  IFNRJ2* activated
            1000.0, # 6:  STAT1c cytoplasmic
            0.0,    # 7:  IFNRJ2*_STAT1c
            0.0,    # 8:  STAT1c* phosphorylated
            0.0,    # 9:  IFNRJ2*_STAT1c*
            0.0,    # 10: STAT1c*_STAT1c* dimer
            100.0,  # 11: SHP2
            0.0,    # 12: IFNRJ2*_SHP2
            50.0,   # 13: PPX
            0.0,    # 14: STAT1c*_PPX
            0.0,    # 15: STAT1c_STAT1c*
            0.0,    # 16: STAT1n*_STAT1n*
            0.0,    # 17: STAT1n*
            60.0,   # 18: PPN
            0.0,    # 19: STAT1n*_PPN
            0.0,    # 20: STAT1n
            0.0,    # 21: STAT1n_STAT1n*
            0.0,    # 22: mRNAn
            0.0,    # 23: mRNAc
            0.0,    # 24: SOCS1
            0.0,    # 25: IFNRJ2*_SOCS1
            0.0,    # 26: Complex
            0.0,    # 27: STAT1c*_STAT1c*_PPX
            0.0,    # 28: STAT1n*_STAT1n*_PPN
            0.0,    # 29: IFNRJ2*_SOCS1_STAT1c
            0.0,    # 30: IFNRJ2*_SHP2_STAT1c
            0.0,    # 31: IFNRJ2*_SHP2_SOCS1
            0.0     # 32: IFNR
        ]

    def calculate_reaction_rates(self, x):
        """Calculate all reaction rates"""
        v = np.zeros(46)
        
        # 1. Receptor activation module
        v[0] = self.compartment_cytoplasm * (self.rates['k1f'] * x[0] * x[1] - self.rates['k1b'] * x[2])
        v[1] = self.compartment_cytoplasm * (self.rates['k2f'] * self.const_species_IFN * x[2] - self.rates['k2b'] * x[3])
        v[2] = self.compartment_cytoplasm * (self.rates['k3f'] * x[3] * x[3] - self.rates['k3b'] * x[4])
        v[3] = self.compartment_cytoplasm * self.rates['k4'] * x[4]

        # 2. STAT1 activation module
        v[4] = self.compartment_cytoplasm * (self.rates['k5f'] * x[6] * x[5] - self.rates['k5b'] * x[7])
        v[5] = self.compartment_cytoplasm * self.rates['k6'] * x[7]
        v[6] = self.compartment_cytoplasm * (self.rates['k7f'] * x[5] * x[8] - self.rates['k7b'] * x[9])
        v[7] = self.compartment_cytoplasm * (self.rates['k8f'] * x[8] * x[8] - self.rates['k8b'] * x[10])

        # 3. Phosphatase regulation module
        v[8:15] = [
            self.compartment_cytoplasm * (self.rates['k9f'] * x[5] * x[11] - self.rates['k9b'] * x[12]),
            self.compartment_cytoplasm * self.rates['k10'] * x[12],
            self.compartment_cytoplasm * (self.rates['k11f'] * x[13] * x[8] - self.rates['k11b'] * x[14]),
            self.compartment_cytoplasm * self.rates['k12'] * x[14],
            self.compartment_cytoplasm * (self.rates['k13f'] * x[13] * x[10] - self.rates['k13b'] * x[27]),
            self.compartment_cytoplasm * self.rates['k14'] * x[27],
            self.compartment_cytoplasm * (self.rates['k7f'] * x[6] * x[8] - self.rates['k7b'] * x[15])
        ]

        # 4. Nuclear transport and reactions
        v[15:23] = [
            self.compartment_nucleus * self.rates['k15'] * x[10],
            self.compartment_nucleus * (self.rates['k16f'] * x[17] * x[17] - self.rates['k16b'] * x[16]),
            self.compartment_nucleus * (self.rates['k17f'] * x[18] * x[17] - self.rates['k17b'] * x[19]),
            self.compartment_nucleus * self.rates['k18'] * x[19],
            self.compartment_nucleus * (self.rates['k17f'] * x[18] * x[16] - self.rates['k17b'] * x[28]),
            self.compartment_nucleus * self.rates['k18'] * x[28],
            self.compartment_nucleus * (self.rates['k7f'] * x[20] * x[17] - self.rates['k7b'] * x[21]),
            self.compartment_cytoplasm * self.rates['k14'] * x[20]
        ]

        # 5. Gene expression module
        v[23] = self.compartment_nucleus * self.rates['k23'] * x[16] / (self.rates['k23b'] + x[16])
        v[24] = self.compartment_nucleus * self.rates['k24'] * x[22]
        v[25] = self.compartment_cytoplasm * self.rates['k25'] * x[23]
        v[26] = self.compartment_cytoplasm * self.rates['k26'] * x[23]
        v[27] = self.compartment_cytoplasm * self.rates['k27'] * x[24]

        # 6. Feedback regulation module
        v[28:46] = [
            self.compartment_cytoplasm * (self.rates['k28f'] * x[24] * x[5] - self.rates['k28b'] * x[25]),
            self.compartment_cytoplasm * (self.rates['k29f'] * x[25] * x[6] - self.rates['k29b'] * x[29]),
            self.compartment_cytoplasm * (self.rates['k9f'] * x[29] * x[11] - self.rates['k9b'] * x[26]),
            self.compartment_cytoplasm * self.rates['k10'] * x[26],
            self.compartment_cytoplasm * self.rates['k10'] * x[31],
            self.compartment_cytoplasm * (self.rates['k9f'] * x[25] * x[11] - self.rates['k9b'] * x[31]),
            self.compartment_cytoplasm * (self.rates['k28f'] * x[12] * x[24] - self.rates['k28b'] * x[31]),
            self.compartment_cytoplasm * (self.rates['k5f'] * x[12] * x[6] - self.rates['k5b'] * x[30]),
            self.compartment_cytoplasm * self.rates['k6'] * x[30],
            self.compartment_cytoplasm * self.rates['k6'] * x[29],
            self.compartment_cytoplasm * self.rates['k10'] * x[25],
            self.compartment_cytoplasm * self.rates['k10'] * x[25],
            self.compartment_cytoplasm * self.rates['k10'] * x[25],
            self.compartment_cytoplasm * self.rates['k10'] * x[25],
            self.compartment_cytoplasm * (self.rates['k28f'] * x[5] * x[24] - self.rates['k28b'] * x[25]),
            self.compartment_cytoplasm * (self.rates['k28f'] * x[12] * x[24] - self.rates['k28b'] * x[31]),
            self.compartment_cytoplasm * (self.rates['k2f'] * self.const_species_IFN * x[0] - self.rates['k2b'] * x[32]),
            self.compartment_cytoplasm * (self.rates['k1f'] * x[32] * x[1] - self.rates['k1b'] * x[3])
        ]
        
        return v

    def system_equations(self, t, x):
        """Complete system of ODEs"""
        v = self.calculate_reaction_rates(x)
        dxdt = np.zeros(33)
        
        # Receptor module ODEs
        dxdt[0:6] = [
            (1/self.compartment_cytoplasm) * (-v[0] - v[44]),  # R
            (1/self.compartment_cytoplasm) * (-v[0] - v[45]),  # JAK
            (1/self.compartment_cytoplasm) * (v[0] - v[1]),    # RJ
            (1/self.compartment_cytoplasm) * (v[1] - 2*v[2] + v[45]),  # IFNRJ
            (1/self.compartment_cytoplasm) * (v[2] - v[3]),    # IFNRJ2
            (1/self.compartment_cytoplasm) * (v[3] - v[4] - v[8] - v[28])  # IFNRJ2*
        ]

        # STAT1 module ODEs
        dxdt[6:11] = [
            (1/self.compartment_cytoplasm) * (-v[4] - v[14] + v[11] + v[13]),  # STAT1c
            (1/self.compartment_cytoplasm) * (v[4] - v[5] - v[35]),  # IFNRJ2*_STAT1c
            (1/self.compartment_cytoplasm) * (v[5] - 2*v[7] - v[10]),  # STAT1c*
            (1/self.compartment_cytoplasm) * (v[6] - v[7]),  # IFNRJ2*_STAT1c*
            (1/self.compartment_cytoplasm) * (v[7] - v[12] - v[15])  # STAT1c*_STAT1c*
        ]

        # Regulatory module ODEs
        dxdt[11:16] = [
            (1/self.compartment_cytoplasm) * (-v[8] - v[33] - v[34] + v[9] + v[31] + v[36]),  # SHP2
            (1/self.compartment_cytoplasm) * (v[8] - v[9]),  # IFNRJ2*_SHP2
            (1/self.compartment_cytoplasm) * (-v[10] - v[12] + v[11]),  # PPX
            (1/self.compartment_cytoplasm) * (v[10] - v[11]),  # STAT1c*_PPX
            (1/self.compartment_cytoplasm) * (v[14] - v[13])  # STAT1c_STAT1c*
        ]

        # Nuclear module ODEs
        dxdt[16:22] = [
            (1/self.compartment_nucleus) * (v[15] - v[19]),  # STAT1n*_STAT1n*
            (1/self.compartment_nucleus) * (2*v[16] - v[17] - v[21]),  # STAT1n*
            (1/self.compartment_nucleus) * (-v[17] - v[19]),  # PPN
            (1/self.compartment_nucleus) * (v[17] - v[18]),  # STAT1n*_PPN
            (1/self.compartment_nucleus) * (v[18] + v[20] - v[22]),  # STAT1n
            (1/self.compartment_nucleus) * (v[21] - v[20])  # STAT1n_STAT1n*
        ]

        # Gene expression module ODEs
        dxdt[22:25] = [
            (1/self.compartment_nucleus) * (v[23] - v[24]),  # mRNAn
            (1/self.compartment_cytoplasm) * (v[24] - v[25] - v[26]),  # mRNAc
            (1/self.compartment_cytoplasm) * (v[25] - v[27] - v[28] - v[42] - v[43])  # SOCS1
        ]

        # Feedback module ODEs
        dxdt[25:33] = [
            (1/self.compartment_cytoplasm) * (v[28] - v[29] - v[33] - v[40]),  # IFNRJ2*_SOCS1
            (1/self.compartment_cytoplasm) * (v[30] - v[31] - v[32]),  # IFNRJ2*_SHP2_SOCS1_STAT1c
            (1/self.compartment_cytoplasm) * (v[12] - v[13]),  # STAT1c*_STAT1c*_PPX
            (1/self.compartment_nucleus) * (v[19] - v[20]),  # STAT1n*_STAT1n*_PPN
            (1/self.compartment_cytoplasm) * (v[29] - v[30] - v[37]),  # IFNRJ2*_SOCS1_STAT1c
            (1/self.compartment_cytoplasm) * (v[35] - v[36]),  # IFNRJ2*_SHP2_STAT1c
            (1/self.compartment_cytoplasm) * (v[33] + v[34] - v[31] - v[32]),  # IFNRJ2*_SHP2_SOCS1
            (1/self.compartment_cytoplasm) * (v[44] - v[45])  # IFNR
        ]
        
        return dxdt

def main():
    """Main function to run JAK-STAT pathway simulation"""
    model = JAKSTATModel()

    # Simulation parameters
    t_span = (0, 100)
    t_eval = np.linspace(0, 100, 1000)

    # Solve system
    solution = solve_ivp(model.system_equations, t_span, model.x0, t_eval=t_eval, method='BDF', rtol=1e-3)
    plot_results(solution)

def plot_results(solution):
    """Plot simulation results for STAT1n-STAT1n*"""
    stat1n_stat1n_star_index = 21
    stat1n_stat1n_star_name = 'STAT1n_STAT1n*'

    plt.figure(figsize=(10, 6))
    plt.plot(solution.t, solution.y[stat1n_stat1n_star_index], label=stat1n_stat1n_star_name)
    
    plt.xlabel('Time (s)')
    plt.ylabel('Concentration')
    plt.title('STAT1n_STAT1n* Dynamics in JAK-STAT Pathway')
    plt.legend(loc='upper right')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('stat1n_stat1n_star_dynamics.png', dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    main()