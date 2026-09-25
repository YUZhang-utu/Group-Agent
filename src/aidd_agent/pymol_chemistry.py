"""Portable chemical-component metadata prepared on the Chat host, not by the LLM."""
import hashlib
from pathlib import Path


def component_metadata(path, expected=None):
    """Read deposited aromatic flags and bonds; never infer charges from residue names."""
    import gemmi
    path = Path(path)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if expected and digest != expected:
        raise ValueError('Structure provenance mismatch during chemical typing')
    if path.suffix.lower() not in {'.cif', '.mmcif'}:
        return dict(status='unavailable', reason='No mmCIF component dictionary', sha256=digest)
    block = gemmi.cif.read_file(str(path)).sole_block()
    components = {}
    for row in block.find(['_chem_comp_atom.comp_id', '_chem_comp_atom.atom_id',
                           '_chem_comp_atom.type_symbol', '_chem_comp_atom.pdbx_aromatic_flag']):
        comp, name, element, aromatic = [gemmi.cif.as_string(x) for x in row]
        components.setdefault(comp, dict(atoms={}, bonds=[]))['atoms'][name] = dict(
            element=element.title(), aromatic=aromatic == 'Y')
    for row in block.find(['_chem_comp_bond.comp_id', '_chem_comp_bond.atom_id_1',
                           '_chem_comp_bond.atom_id_2']):
        comp, left, right = [gemmi.cif.as_string(x) for x in row]
        if comp in components:
            components[comp]['bonds'].append([left, right])
    return dict(status='available' if components else 'unavailable',
                reason='' if components else 'No deposited component aromatic flags',
                sha256=digest, components=components)


def apply_components(atoms, edges, metadata, aromatic):
    """Map by residue identity and unique atom names, preserving live coordinates.

    Existing bonds are required: removing a bond or atom in the viewer must not
    silently recreate it from a source template. Conflicting elements are rejected.
    """
    from collections import defaultdict
    groups = defaultdict(lambda: defaultdict(list))
    for i, atom in atoms.items():
        groups[(atom.get('segi', ''), atom['chain'], atom['resi'], atom['resn'])][atom['name']].append(i)
    typed = set(); ambiguous = 0
    graph = {frozenset(edge) for edge in edges}
    components = metadata.get('components', {})
    for identity, names in groups.items():
        template = components.get(identity[-1], {})
        definitions = template.get('atoms', {})
        for name, indices in names.items():
            if len(indices) != 1:
                ambiguous += 1
                continue
            i = indices[0]; definition = definitions.get(name, {})
            if definition.get('element') != atoms[i]['elem']:
                continue
            required = [b if a == name else a for a, b in template.get('bonds', []) if name in (a, b)]
            # Only complete heavy-atom neighborhoods receive source aromatic typing.
            required = [n for n in required if definitions.get(n, {}).get('element') != 'H']
            if any(len(names.get(n, [])) != 1 or frozenset((i, names[n][0])) not in graph for n in required):
                continue
            typed.add(i)
            if definition.get('aromatic'):
                aromatic.add(i)
    return aromatic, dict(component_typed_atoms=len(typed), ambiguous_atom_names=ambiguous,
                         component_status=metadata.get('status', 'unavailable'),
                         component_reason=metadata.get('reason', 'No source component metadata'),
                         source_sha256=metadata.get('sha256'))
