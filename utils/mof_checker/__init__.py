import ipywidgets as ipw
import functools as fct
from traitlets import Instance, Int, List, observe
from ase import Atoms, Atom
from mofchecker import MOFChecker
import tempfile
from itertools import chain
from IPython.display import clear_output


ENABLED_CHECKS = [
    "no_atomic_overlaps",
    "no_undercoordinated_carbon",
    "no_overcoordinated_carbon",
    "no_overcoordinated_hydrogen",
    "no_overcoordinated_nitrogen",
    "no_undercoordinated_nitrogen",
    "no_undercoordinated_rare_earth",
    "no_floating_molecule",
    "no_false_terminal_oxo",
]

def flatten_list(lst):
    for element in lst:
        if type(element) is list:
            for element2 in flatten_list(element):
                yield element2
        else:
            yield element

class MofCheckerWidget(ipw.VBox):
    """Widget that allows to check MOF structure correctness."""
    structure = Instance(Atoms, allow_none=True)
    selection = List(Int)

    def __init__(self, title=''):
        self.title = title
        button_check = ipw.Button(description="Analyse structure")
        button_check.on_click(self.check_structure)
        self.checks = []
        self._output = ipw.Output()
        self.missing_h_button = ipw.Button(description="Add", layout={"width": "initial"})

        super().__init__(children=[button_check,
                                   self._output,
                                  ])
    
    def check_structure(self, _=None):
        
        if self.structure is None:
            return

        self.mfchk = MOFChecker.from_ase(self.structure, primitive=False)
        self.checks = []
        if self.structure:
            self.checks = [check for (check_name, check) in self.mfchk.checks.items() if check_name in ENABLED_CHECKS]
        
        issue_found = False
        
        with self._output:
            clear_output()
            for check in self.checks:
                if not check.is_ok:
                    issue_found = True
                    try:
                        selection = list(chain.from_iterable(check.flagged_indices))
                    except TypeError:
                        selection = check.flagged_indices
                    button = ipw.Button(description="Select", layout={"width": "initial"})
                    button.on_click(fct.partial(self._select_atoms, selection=selection))
                    text = ipw.HTML("Found " + check.name.lower())
                    display(ipw.HBox([text, button]))
                    
            if self.mfchk.undercoordinated_c_candidate_positions + self.mfchk.undercoordinated_n_candidate_positions:
                self.missing_h_button.disabled = False
                self.missing_h_button.on_click(self.add_missing_hydrogens)
                display(ipw.HBox([ipw.HTML("Missing hydrogens found"), self.missing_h_button]))
            
            if not issue_found:
                display(ipw.HTML("No issues found \u2705"))
            
    def _select_atoms(self, _=None, selection=None):
        self.selection = selection

    def add_missing_hydrogens(self, _=None):
        """Add atoms."""
        atoms = self.structure.copy()
        selection = self.selection

        h_positions = self.mfchk.undercoordinated_c_candidate_positions + self.mfchk.undercoordinated_n_candidate_positions
        for pos in flatten_list(h_positions):
            atoms += Atom('H', pos)

        self.structure = atoms
        self.selection = selection

    @observe('structure')
    def observe_structure(self, _=None):
        self.missing_h_button.disabled = True
        with self._output:
            clear_output()

class CheckMofStructure(ipw.VBox):
    structure = Instance(Atoms, allow_none=True)

    def __init__(self, **kwargs):
        self.text = ipw.HTML()
        self.check_button = ipw.Button(description="Final Structure Check", button_style='warning')
        self.check_button.on_click(self.check_structure)
        super().__init__([self.check_button, self.text], **kwargs)

    def check_structure(self, _=None):
        if self.structure is None:
            self.text.value = ''
        else:
            self.text.value = ""
            self.mfchk = MOFChecker.from_ase(self.structure, primitive=False)
            issue_found = False
            for (check_name, check) in self.mfchk.checks.items():
                if check_name in ENABLED_CHECKS:
                    if not check.is_ok:
                        issue_found = True
                    self.text.value += "\u2705 OK: " if check.is_ok else "\u274C Fail: "
                    self.text.value += f"{check.description} </br>"
            self.check_button.button_style = 'danger' if issue_found else 'success'

    @observe('structure')
    def observe_structure(self, _=None):
        self.check_button.button_style = 'warning'
        self.text.value = ""
