using Catalyst

# Define the reaction system
@parameters k0 k1 k2 k4
@variables t
@species IFNγ(t) Ir(t) IIr(t) STAT1Uc(t) STAT1Dc(t)

# Define the reactions
rxs = [
    # IFN-γ degradation
    Reaction(k0, [IFNγ], nothing),
    
    # IFN-γ + Ir binding to form IIr
    Reaction(k1, [IFNγ, Ir], [IIr]),
    
    # IIr degradation
    Reaction(k2, [IIr], nothing),
    
    # IIr + STAT1Uc -> STAT1Dc
    Reaction(k4, [IIr, STAT1Uc], [STAT1Dc, IIr])
]

# Create the reaction system
@named rs = ReactionSystem(rxs, t)

# Convert to ODESystem
osys = convert(ODESystem, rs) 