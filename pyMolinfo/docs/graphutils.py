# Graph Utils
# ----------------
# import libs


def combine_chains1(chains):
    '''
    Combine multiple chains into a single chain, using the following rules:
    1. Cx-{n} in chain 1 is a placeholder for the chain 2
    2. {n}-Cn in chain 2 is a placeholder for the chain 1
    3. {n} is removed from the combined chain
    4. the index of the atoms in the second chain is reset and starts from the last atom in the first chain
    5. n is integer and varies from 1 to infinity, and exist in both chains

    Parameters
    ----------
    chains : list
        List of chains to be combined

    Returns
    -------
    list
        Combined chain

    Examples
    --------
    # Chain 1
    Chain1 = ["C1=C2","C2-{1}","C3=C4","C4-C5","C5=C6","C6-C1"]
    # Chain 2
    Chain2 = ["XX1-C2","{1}-C2","C2=C3","C3-XX5"]
    # chains
    Chains = [Chain1, Chain2]
    # combine chains
    >>> combine_chains(Chains)
    '''
    pass
