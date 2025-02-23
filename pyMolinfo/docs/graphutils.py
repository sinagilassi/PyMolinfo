# Graph Utils
# ----------------
# import libs


import re


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


Chain1 = ["C1=C2", "C2-{1}", "C3=C4", "C4-C5", "C5=C6", "C6-C1"]
Chain2 = ["XX1-C2", "{1}-C2", "C2=C3", "C3-XX4"]
Chains = [Chain1, Chain2]


def extract_highest_index(chain):
    """Extracts the highest numerical index from atom labels in the given chain."""
    indices = []
    for bond in chain:
        # Extract numbers from atoms
        atoms = re.findall(r'[A-Za-z]+(\d+)', bond)
        indices.extend(map(int, atoms))
    return max(indices) if indices else 0


# res
# res_1 = extract_highest_index(Chain1)
# res_2 = extract_highest_index(Chain2)
# print(res_1, res_2)


def update_indices(chain, offset, gate_map):
    """Update indices of atoms in the chain by adding an offset while keeping gates mapped."""
    updated_chain = []
    for bond in chain:
        atoms = bond.split('-')
        new_bond = []
        for atom in atoms:
            if atom in gate_map:
                # Replace gate {n} with correct atom
                new_bond.append(gate_map[atom])
            else:
                match = re.match(r'([A-Za-z]+)(\d+)', atom)
                if match:
                    prefix, num = match.groups()
                    new_atom = f"{prefix}{int(num) + offset}"  # Shift index
                    new_bond.append(new_atom)
                else:
                    new_bond.append(atom)
        updated_chain.append("-".join(new_bond))  # Ensure all bonds are intact
    return updated_chain


def combine_chains(chains):
    """
    Combine multiple chains into a single chain following the given rules.

    Parameters
    ----------
    chains : list of lists
        A list containing multiple chains, where each chain is a list of strings.

    Returns
    -------
    list
        A combined chain with placeholders resolved.
    """
    combined_chain = chains[0][:]  # Choose the first chain as the reference
    highest_index = extract_highest_index(combined_chain)  # Find highest index

    for chain in chains[1:]:
        gate_map = {}  # Store {n} connections to real atom numbers
        temp_chain = []

        for bond in chain:
            match = re.search(r'\{(\d+)\}', bond)  # Find {n} placeholders
            if match:
                n = match.group(1)
                for i, prev_bond in enumerate(combined_chain):
                    if f"-{{{n}}}" in prev_bond:  # If Cx-{n} exists in previous chain
                        real_atom = prev_bond.split(
                            '-')[0]  # Get the real atom
                        gate_map[f"{{{n}}}"] = real_atom  # Store mapping
                        combined_chain[i] = prev_bond.replace(
                            f"-{{{n}}}", f"-C{highest_index+1}")  # Replace gate
                    elif f"{{{n}}}-" in prev_bond:  # If {n}-Cx exists in previous chain
                        real_atom = prev_bond.split(
                            '-')[1]  # Get the real atom
                        gate_map[f"{{{n}}}"] = real_atom  # Store mapping
                        combined_chain[i] = prev_bond.replace(
                            f"{{{n}}}-", f"C{highest_index+1}-")  # Replace gate
            else:
                temp_chain.append(bond)

        # Update indices for the new chain and merge
        highest_index += 1  # Ensure new atoms start from max index + 1
        updated_chain = update_indices(temp_chain, highest_index, gate_map)
        combined_chain.extend(updated_chain)

    return combined_chain


# Example usage:
Chain1 = ["C1=C2", "C2-{1}", "C3=C4", "C4-C5", "C5=C6", "C6-C1"]
Chain2 = ["XX1-C2", "{1}-C2", "C2=C3", "C3-XX4"]

Chains = [Chain1, Chain2]
combined = combine_chains(Chains)
print(combined)
