from openff.toolkit.topology import Molecule
from openff.toolkit.utils.toolkits import ToolkitRegistry, AmberToolsToolkitWrapper
toolkit_precedence = [AmberToolsToolkitWrapper]
toolkit_registry = ToolkitRegistry(toolkit_precedence)
molecule = Molecule.from_file("ethanol.sdf")
molecule.assign_partial_charges("am1bcc", use_conformers=molecule.conformers, toolkit_registry=toolkit_registry)

