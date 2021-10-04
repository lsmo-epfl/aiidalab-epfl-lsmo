import ipywidgets as ipw


def get_start_widget(appbase, jupbase):
    # http://fontawesome.io/icons/
    template = """
    <table>
    <tr>
        <th style="text-align:center">Prepare the structure</th>
        <th style="width:70px" rowspan=2></th>
        <th style="text-align:center">Pore analysis</th>
        <th style="width:70px" rowspan=2></th>
        <th style="text-align:center">Isotherm calculations</th>
    <tr>

    <td valign="top"><ul>
    <li><a href="{appbase}/multistage_geo_opt_ddec.ipynb" target="_blank">Geometry Optimization and Charges</a>
    <li><a href="{appbase}/results.ipynb?process_label=Cp2kMultistageDdecWorkChain" target="_blank">Results</a>
    </ul></td>

    <td valign="top"><ul>
    <li><a href="{appbase}/pore_analysis.ipynb" target="_blank">Pore Analysis</a>
    <li><a href="{appbase}/results.ipynb?process_label=NetworkCalculation" target="_blank">Results</a>
    </ul></td>

    <td valign="top"><ul>
    <li><a href="{appbase}/compute_isotherm.ipynb" target="_blank">Compute one</a>
    <li><a href="{appbase}/analyse_isotherm_results.ipynb" target="_blank">Analyse the results</a>
    <li><a href="{appbase}/results.ipynb?process_label=IsothermWorkChain" target="_blank">Results</a>
    </ul></td>
    </tr></table>
"""

    html = template.format(appbase=appbase, jupbase=jupbase)
    return ipw.HTML(html)


# EOF
