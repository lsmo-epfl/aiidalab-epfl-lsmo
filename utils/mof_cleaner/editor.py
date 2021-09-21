import ipywidgets as ipw
from aiidalab_widgets_base import StructureManagerWidget
from ase import Atoms
from traitlets import (
    Instance,
    Int,
    List,
    Unicode,
    Union,
    default,
    link,
    observe,
    validate,
)

from .MOF_cleaner import clean_MOF_free_solvents, find_atomic_overlap


class SolventOverlapCleaner(ipw.VBox):
    """Widget that allows for the basic structure editing."""

    structure = Instance(Atoms, allow_none=True)
    selection = List(Int)

    def __init__(self, title=""):
        self.title = title
        button_freesolvent = ipw.Button(description="Find free solvent")
        button_freesolvent.on_click(self.find_freesolvent)
        self.tolerance = ipw.FloatSlider(
            value=1.0,
            min=0.1,
            max=2.0,
            step=0.1,
            description="Tolerance",
            continuous_update=False,
            orientation="horizontal",
            readout=True,
            readout_format=".1f",
        )
        button_atomicoverlap = ipw.Button(description="Find atomic overlap")
        button_atomicoverlap.on_click(self.find_atomicoverlap)
        self.moleculesize_ratio = ipw.FloatText(
            description="Solvent size ratio to the biggest component",
            value=0.5,
            step=0.1,
            style={"description_width": "initial"},
            layout={"width": "400px"},
        )
        self.framework_min_size = ipw.FloatText(
            description="Minimum size for the framework components",
            value=20,
            step=1,
            style={"description_width": "initial"},
            layout={"width": "400px"},
        )
        self.solvent_max_size = ipw.FloatText(
            description="Maximum size of the solvents",
            value=100,
            step=5,
            style={"description_width": "initial"},
            layout={"width": "300px"},
        )

        self.check_for_metal = ipw.Checkbox(
            value=True,
            description="only MOFs: Check metal containing components",
            style={"description_width": "initial"},
            layout={"width": "initial"},
        )

        self.output = ipw.HTML("")
        super().__init__(
            children=[
                ipw.HBox(
                    [
                        ipw.VBox(
                            [
                                ipw.HTML("<b>Solvent finder</b>"),
                                self.moleculesize_ratio,
                                self.framework_min_size,
                                self.solvent_max_size,
                                self.check_for_metal,
                                button_freesolvent,
                            ]
                        ),
                        ipw.VBox(
                            [
                                ipw.HTML("<b>Atomic overlap finder</b>"),
                                self.tolerance,
                                button_atomicoverlap,
                            ],
                            layout={"margin": "0px 0px 0px 100px"},
                        ),
                    ]
                ),
                self.output,
            ]
        )

    def find_freesolvent(self, _=None):
        self.output.value = ""
        if self.structure is None:
            self.output.value = "No structure provided."
            return
        solvent_indices, solvent_atomtypes, log = clean_MOF_free_solvents(
            self.structure,
            moleculesize_ratio=self.moleculesize_ratio.value,
            framework_min_size=self.framework_min_size.value,
            solvent_max_size=self.solvent_max_size.value,
            check_metal=self.check_for_metal.value,
        )

        self.output.value = log
        self.selection = solvent_indices

    def find_atomicoverlap(self, _=None):
        self.output.value = ""
        if self.structure is None:
            self.output.value = "No structure provided."
            return
        overlap_indices, overlap_atomtypes, log = find_atomic_overlap(
            self.structure, tolerance=self.tolerance.value
        )
        self.output.value = log
        self.selection = overlap_indices
