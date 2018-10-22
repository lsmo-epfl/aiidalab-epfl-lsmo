from aiida.orm import CalculationFactory, DataFactory
from aiida.orm import load_node
from aiida.orm.code import Code
from aiida.orm.data.base import Bool, Str, Float
from aiida.work import workfunction as wf
from aiida.work.run import run, submit
from aiida.work.workchain import WorkChain, ToContext, if_, while_, Outputs

# subworkflows
from workflows.charges import DdecChargesWorkChain
from aiida_cp2k.workflows import Cp2kRobustGeoOptWorkChain
from aiida_raspa.workflows import RaspaConvergeWorkChain
from aiida_zeopp.workflows import ZeoppBlockPocketsWorkChain

# data objects
ArrayData = DataFactory('array')
CifData = DataFactory('cif')
NetworkParameters = DataFactory('zeopp.parameters')
ParameterData = DataFactory('parameter')
RemoteData = DataFactory('remote')
StructureData = DataFactory('structure')

def multiply_unit_cell (cif, threshold):
    from math import cos, sin, sqrt, pi
    import numpy as np
    deg2rad=pi/180.
    threshold /= 2.0

    struct=cif.values.dictionary.itervalues().next()

    a = float(struct['_cell_length_a'])
    b = float(struct['_cell_length_b'])
    c = float(struct['_cell_length_c'])

    alpha = float(struct['_cell_angle_alpha'])*deg2rad
    beta  = float(struct['_cell_angle_beta'])*deg2rad
    gamma = float(struct['_cell_angle_gamma'])*deg2rad

    # this code makes sure that the structure is replicated enough times in x, y, z direction
    # in order to be compatible with the threshold value
    
    # first step is computing cell parameters according to  https://en.wikipedia.org/wiki/Fractional_coordinates
    v = sqrt(1-cos(alpha)**2-cos(beta)**2-cos(gamma)**2+2*cos(alpha)*cos(beta)*cos(gamma))
    cell=np.zeros((3,3))
    cell[0,:] = [a, 0, 0]
    cell[1,:] = [b*cos(gamma), b*sin(gamma),0]
    cell[2,:] = [c*cos(beta), c*(cos(alpha)-cos(beta)*cos(gamma))/(sin(gamma)),c*v/sin(gamma)]
    cell=np.array(cell)

    # diagonalizing the cell matrix
    diag = np.diag(cell)
    # and computing nx, ny and nz
    return tuple(int(i) for i in np.ceil(threshold/diag*2.))


@wf
def from_cif_to_structuredata(cif, threshold):
    """Helper function that converts CifData object into StructureData"""
    repeat = multiply_unit_cell(cif, threshold=threshold.value)
    structure = cif.get_ase().repeat(repeat)
    return StructureData(ase=structure).store()


class Isotherm(WorkChain):
    """Workchain that for a given matherial will compute an isotherm of a
    certain gaz adsorption."""
    @classmethod
    def define(cls, spec):
        super(Isotherm, cls).define(spec)
        
        # structure, adsorbant, pressures
        spec.input('structure', valid_type=CifData)
        spec.input("probe_molecule", valid_type=ParameterData)
        spec.input("pressures", valid_type=ArrayData)
        spec.input("min_cell_size", valid_type=Float)

        # cp2k
        spec.input('cp2k_code', valid_type=Code)
        spec.input("_cp2k_options", valid_type=dict, default=None, required=False)
        spec.input('cp2k_parent_folder', valid_type=RemoteData, default=None, required=False)

        # ddec
        spec.input('ddec_code', valid_type=Code)
        spec.input("_ddec_options", valid_type=dict, default=None, required=False)

        # zeopp
        spec.input('zeopp_code', valid_type=Code)
        spec.input("_zeopp_options", valid_type=dict, default=None, required=False)

        # raspa
        spec.input("raspa_code", valid_type=Code)
        spec.input("raspa_parameters", valid_type=ParameterData)
        spec.input("_raspa_options", valid_type=dict, default=None, required=False)

        # settings
        spec.input("_interactive", valid_type=bool, default=False, required=False)
        spec.input("_usecharges", valid_type=bool, default=False, required=False)

        # workflow
        spec.outline(
            cls.init,
            cls.run_geo_opt,
            cls.parse_geo_opt,
            if_(cls.should_use_charges)(
                cls.run_point_charges,
                cls.parse_point_charges,
            ),
            cls.run_geom_zeopp,
            cls.parse_geom_zeopp,
            cls.run_henry_raspa,
            while_(cls.should_run_loading_raspa)(
                cls.run_loading_raspa,
                cls.parse_loading_raspa,
            ),
            cls.return_results,
        )

        # TODO: once the workflow is ready, explicitely specify the outputs
        spec.dynamic_output()

    def init(self):
        """Initialize variables and the pressures we want to compute"""
        self.ctx.structure = self.inputs.structure
        self.ctx.pressures = self.inputs.pressures.get_array("pressures")
        self.ctx.current_p_index = 0
        self.ctx.result = []

        self.ctx.raspa_parameters = self.inputs.raspa_parameters.get_dict()

        if self.inputs._usecharges:
            self.ctx.raspa_parameters['GeneralSettings']['UseChargesFromCIFFile'] = "yes"
        self.ctx.restart_raspa_calc = None

    def run_geo_opt(self):
        """Optimize geometry."""
        threshold = self.inputs.min_cell_size
    # uncomment this for the test runs
#        params_dict = {
#                'MOTION':{
#                    'MD':{
#                        'STEPS': 5,
#                        },
#                    'GEO_OPT': {
#                        'MAX_ITER': 5,
#                    },
#                    'CELL_OPT': {
#                        'MAX_ITER': 5,
#                    },
#                },
#        }
        parameters = ParameterData(dict=params_dict)

        inputs = {
            'code'      : self.inputs.cp2k_code,
            'structure' : from_cif_to_structuredata(cif=self.ctx.structure, threshold=threshold),
            'parameters': parameters,
            '_options'  : self.inputs._cp2k_options,
            '_label'    : "Cp2kRobustGeoOptWorkChain",
        }

        # Create the calculation process and launch it
        running = submit(Cp2kRobustGeoOptWorkChain, **inputs)
        self.report("pk: {} | Running Cp2kRobustGeoOptWorkChain to optimize geometry".format(running.pid))
        return ToContext(geo_opt_calc=Outputs(running))
    
    def parse_geo_opt(self):
        """Extract optimized structure and put it into self.ctx.structure"""
        self.ctx.structure = self.ctx.geo_opt_calc['output_structure']

    def should_use_charges(self):
        """Whether it is needed to employ charges."""
        return self.inputs._usecharges

    def run_point_charges(self):
        """Compute the charge-density of a structure that can be later
        used for extracting ddec point charges."""

        inputs = {
            'structure'          : self.ctx.structure,
            'cp2k_code'          : self.inputs.cp2k_code,
            'cp2k_parent_folder' : self.ctx.geo_opt_calc['remote_folder'],
            '_cp2k_options'      : self.inputs._cp2k_options,
            'ddec_code'          : self.inputs.ddec_code,
            '_ddec_options'      : self.inputs._ddec_options,
            '_label'             : "DdecChargesWorkChain",
        }

        # Create the calculation process and launch it
        running = submit(DdecChargesWorkChain, **inputs)
        self.report("pk: {} | Running DdecChargesWorkChain to compute the point charges".format(running.pid))
        return ToContext(point_charges_calc=Outputs(running))

    def parse_point_charges(self):
        """Extract structure with charges and put it into self.ctx.structure"""
        self.ctx.structure = self.ctx.point_charges_calc['output_structure']

    def run_geom_zeopp(self):
        """This is the main function that will perform a raspa calculation for the current pressure."""
        # network parameters
        sigma = self.inputs.probe_molecule.dict.sigma

        # Create the input dictionary
        inputs = {
            'probe_radius' : Float(sigma),
            'structure'    : self.ctx.structure,
            'zeopp_code'   : self.inputs.zeopp_code,
            '_options'     : self.inputs._zeopp_options,
            '_label'       : "ZeoppBlockPocketsWorkChain",
        }

        # Create the calculation process and launch it
        running = submit(ZeoppBlockPocketsWorkChain, **inputs)
        self.report("pk: {} | Running geometry analysis with zeo++".format(running.pid))

        return ToContext(zeopp=Outputs(running))

    def parse_geom_zeopp(self):
        """Extract the pressure and loading average of the last completed raspa calculation."""
        self.ctx.raspa_parameters['GeneralSettings']['HeliumVoidFraction'] = \
        self.ctx.zeopp["output_parameters"].dict.POAV_Volume_fraction

    def should_run_loading_raspa(self):
        """We run another raspa calculation only if the current iteration is smaller than
        the total number of pressures we want to compute."""
        return self.ctx.current_p_index < len(self.ctx.pressures)

    def run_loading_raspa(self):
        """This function will run RaspaConvergeWorkChain for the current pressure"""
        pressure = self.ctx.pressures[self.ctx.current_p_index]
        self.ctx.raspa_parameters['GeneralSettings']['ExternalPressure'] = pressure

        parameters = ParameterData(dict=self.ctx.raspa_parameters).store()
        # Create the input dictionary
        inputs = {
            'code'       : self.inputs.raspa_code,
            'structure'  : self.ctx.structure,
            'parameters' : parameters,
            '_options'    : self.inputs._raspa_options,
            '_label'     : "run_loading_raspa",
        }
        try:
            inputs['block_component_0'] = self.ctx.zeopp['block']
        except:
            pass

        if self.ctx.restart_raspa_calc is not None:
            inputs['retrieved_parent_folder'] = self.ctx.restart_raspa_calc

        # Create the calculation process and launch it
        running = submit(RaspaConvergeWorkChain, **inputs)
        self.report("pk: {} | Running raspa for the pressure {} [bar]".format(running.pid, pressure/1e5))
        self.ctx.current_p_index += 1
        return ToContext(raspa_loading=Outputs(running))

    def parse_loading_raspa(self):
        """Extract the pressure and loading average of the last completed raspa calculation"""
        self.ctx.restart_raspa_calc = self.ctx.raspa_loading['retrieved_parent_folder']
        pressure = self.ctx.raspa_parameters['GeneralSettings']['ExternalPressure']/1e5
        loading_average = self.ctx.raspa_loading["component_0"].dict.loading_absolute_average
        self.ctx.result.append((pressure, loading_average))

    def run_henry_raspa(self):
        """This is the main function that will perform a raspaa calculation for the current pressure"""
        raspa_parameters = self.inputs.raspa_parameters.get_dict()
        raspa_parameters['GeneralSettings'].pop('ExternalPressure')
        for i, comp in enumerate(raspa_parameters['Component']):
            name = comp['MoleculeName']
            raspa_parameters['Component'][0] = {
                "MoleculeName"                     : name,
                "MoleculeDefinition"               : "TraPPE",
                "WidomProbability"                 : 1.0,
                "CreateNumberOfMolecules"          : 0,
            }

        parameters = ParameterData(dict=raspa_parameters).store()

        # Create the input dictionary
        inputs = {
            'code'       : self.inputs.raspa_code,
            'structure'  : self.ctx.structure,
            'parameters' : parameters,
            '_options'   : self.inputs._raspa_options,
            '_label'     : "RaspaConvergeWorkChain",
        }

        # Create the calculation process and launch it
        running = submit(RaspaConvergeWorkChain, **inputs)
        self.report("pk: {} | Running raspa for the Henry coefficients".format(running.pid))

        return ToContext(raspa_henry=Outputs(running))

    def return_results(self):
        """Attach the results of the raspa calculation and the initial structure to the outputs."""
        self.out("result", ParameterData(dict={"isotherm": self.ctx.result}).store())
        self.report("Workchain <{}> completed successfully".format(self.calc.pk))
        return

# EOF
