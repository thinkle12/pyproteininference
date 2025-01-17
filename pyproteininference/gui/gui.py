# macOS packaging support
from multiprocessing import freeze_support  # noqa

freeze_support()  # noqa

import logging
import sys
from multiprocessing import Manager, Queue
from typing import Tuple

import webview
from nicegui import app, run, ui, native

import pyproteininference.pipeline
from pyproteininference.gui.configuration import Configuration

logger = logging.getLogger()

logging.basicConfig(
    stream=sys.stderr,
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)


class LogElementHandler(logging.Handler):
    """A logging handler that emits messages to a log element."""

    def __init__(self, element: ui.log, level: int = logging.NOTSET) -> None:
        self.element = element
        super().__init__(level)

    def emit(self, record: logging.LogRecord) -> None:
        try:
            msg = self.format(record)
            self.element.push(msg)
        except Exception:
            self.handleError(record)


def run_inference_analysis_async(q: Queue, config) -> str:
    """Run some heavy computation that updates the progress bar through the queue."""
    pipeline = pyproteininference.pipeline.ProteinInferencePipeline.create_from_gui_config(q, config)
    pipeline.execute()
    return "Protein Inference Analysis Completed"


@ui.page('/')
def index():

    config = Configuration()

    with ui.tabs().classes('w-full') as tabs:
        execute = ui.tab('Run')
        configure = ui.tab('Configure')
        output = ui.tab('Output')
        about = ui.tab('About')
        tabs.set_value('Run')

    with ui.tab_panels(tabs, value=execute).classes('w-full'):

        with ui.tab_panel(output):

            with ui.row().classes('w-full border'):
                current_step = ui.markdown(
                    '_To view output, first configure the correct parameters on the Configure tab '
                    'and then select the input and output files on the Run tab, '
                    'and then click Run._'
                )
                ui.space()
                progress_spinner = ui.spinner('dots', size='lg', color='red')
            progressbar_output = ui.linear_progress(value=0, show_value=False).props("instant-feedback")

            progressbar_output.visible = False
            progress_spinner.visible = False
            ui.separator()

        with ui.tab_panel(execute):
            with ui.grid(columns=3):

                async def pick_input_files() -> None:
                    result = await app.native.main_window.create_file_dialog(allow_multiple=True)
                    config.input_files = result

                ui.label('Input result files:')
                ui.label('').bind_text_from(
                    config, 'input_files', backward=lambda files: " ".join([file.split('/')[-1] for file in files])
                )
                ui.button('Choose result file(s)', on_click=pick_input_files, icon='folder')

                async def pick_fasta_file() -> None:
                    result = await app.native.main_window.create_file_dialog(allow_multiple=False)
                    config.fasta_file = result

                ui.label('Fasta Database:')
                ui.label('').bind_text_from(
                    config, 'fasta_file', backward=lambda file: "" if file == "" else file[0].split('/')[-1]
                )
                ui.button('Choose fasta file', on_click=pick_fasta_file, icon='file')

                async def pick_output_file() -> None:
                    result = await app.native.main_window.create_file_dialog(
                        dialog_type=webview.SAVE_DIALOG, save_filename="inference_output.txt"
                    )
                    config.output_file = result

                ui.label('Output file:')
                ui.label('').bind_text_from(
                    config, 'output_file', backward=lambda file: "" if file == "" else file.split('/')[-1]
                )
                ui.button('Set output file', on_click=pick_output_file, icon='save')

                async def start_computation():

                    if len(config.input_files) < 1:
                        ui.notify("Please choose input result file(s) for the analysis and try again.")
                    elif len(config.output_file) < 1:
                        ui.notify("Please choose the output filename for the analysis and try again.")
                    else:
                        tabs.set_value('Output')
                        progressbar_output.visible = True
                        progress_spinner.visible = True
                        current_step.set_content("Starting analysis...")
                        result = await run.cpu_bound(run_inference_analysis_async, queue, config)
                        ui.notify(result)
                        current_step.set_content("Analysis is complete.")
                    progressbar_output.visible = False
                    progress_spinner.visible = False

                def update_progressbar(progress_value: Tuple[int, str]):
                    progressbar_output.set_value(progress_value[0])
                    current_step.set_content(progress_value[1])

                # Create a queue to communicate with the heavy computation process
                queue = Manager().Queue()
                # Update the progress bar on the main process
                ui.timer(
                    0.5,
                    callback=lambda: update_progressbar(
                        queue.get() if not queue.empty() else (progressbar_output.value, current_step.content)
                    ),
                )

                ui.button('Run Analysis', on_click=start_computation)

        with ui.tab_panel(configure):
            ui.markdown('## General')
            with ui.grid(columns=3):
                ui.label('False Discovery Rate')
                ui.number(label="FDR", min=0, max=1, precision=3, step=0.01).bind_value(config, 'false_discovery_rate')
                ui.markdown('_Enter a value between 0 and 1.  For example, 0.01 = 1%_')

                ui.label('Use Protein Picker')
                ui.switch().bind_value(config, 'picker')
                ui.markdown(
                    'Protein Picker is an algorithm that treats target and decoy proteins as pairs when '
                    'analyzed as separate target / decoy searches (as opposed to using a joint, concatenated '
                    'target/decoy database for a single analysis. '
                    'If both the target and decoy proteins are identified from the searches when protein '
                    'picker is run, then the target and decoy scores are compared with one another. The '
                    'one with the better score is kept to continue on in the analysis while the one with '
                    'the worse score gets filtered out of the analysis.'
                )

                ui.label('Use FASTA for Protein Information')
                ui.switch().bind_value(config, 'use_alt_proteins')
                ui.markdown(
                    '_If selected and a fasta database is supplied, look up additional potential peptide to '
                    'protein mappings from the fasta file and incorporate them into the results.'
                )

                ui.label('XML Parser Method')
                ui.select(
                    {
                        "openms": "OpenMS",
                        "pyteomics": "Pyteomics",
                    },
                    label="XML Parser Method",
                ).bind_value(config, 'xml_input_parser')
                ui.markdown(
                    '_If selected and a fasta database is supplied, look up additional potential peptide to '
                    'protein mappings from the fasta file and incorporate them into the results.'
                )

            ui.separator()

            ui.markdown('## Inference')
            with ui.grid(columns=3):
                ui.label('Inference Type')
                ui.select(
                    {
                        "parsimony": "Parsimony",
                        "exclusion": "Exclusion",
                        "inclusion": "Inclusion",
                        "peptide_centric": "Peptide Centric",
                        "first_protein": "First Protein",
                    },
                    label="Inference Type",
                ).bind_value(config, 'inference_type')
                ui.markdown(
                    "_Choose which inference method to run. `Parsimony` will map peptides to the minimal set of "
                    "proteins. `Exclusion` will remove all non-distinguishing peptides from further analysis. "
                    "`Inclusion` will map all peptides to all possible proteins. `Peptide Centric` will create "
                    "protein groups based on proteins that contain the same set of peptides from the search. "
                    "Finally, `First Protein` simply chooses the first protein from the list (and is included "
                    "just for completeness)._"
                )

                ui.label('Grouping Method')
                ui.select(
                    {
                        "subset_peptides": "Subset Peptides",
                        "shared_peptides": "Shared Peptides",
                        "parsimonious_grouping": "Most Parsimonious Grouping",
                        "None": "None",
                    },
                    label="Grouping Method",
                ).bind_value(config, 'grouping_type')
                ui.markdown(
                    "_How to group proteins for a given 'Inference Type'. This parameter only effects grouped "
                    "proteins and has no impact on protein leads, except for 'Most Parsimonious Grouping', which "
                    "will report as a group lead the set of identical proteins which are explained by the peptides "
                    "contained within the group. This will only influence how proteins are grouped in 'Parsimony' "
                    "and 'Inclusion'._"
                )

            ui.separator()

            ui.markdown('## Data Filtering')
            with ui.grid(columns=3):
                ui.label('PSM PEP Filter')
                ui.number(label="PEP", min=0, max=1, precision=3, step=0.01).bind_value(config, 'pep_restriction')
                ui.markdown(
                    'Optional: Select a PSM PEP (Posterior Error Probability) value to filter on. '
                    'All PSMs with PEP values greater than this number will be filtered out from '
                    'protein inference analysis._'
                )

                ui.label('PSM q-value Filter')
                ui.number(label="q-value", min=0, max=1, precision=3, step=0.01).bind_value(
                    config, 'q_value_restriction'
                )
                ui.markdown(
                    'Optional: Select a PSM q-value to filter on. All PSMs with PEP values greater than this number '
                    'will be filtered out from protein inference analysis._'
                )

                ui.label('PSM Minimum Length Filter')
                ui.number(label="Peptide Length", min=0, max=10000, precision=0, step=1).bind_value(
                    config, 'peptide_length_restriction'
                )
                ui.markdown(
                    'Optional: Select a minimum residue length for peptides. All PSMs with lengths less than this '
                    'number will be filtered out from protein inference analysis._'
                )

                ui.label('Alternative Protein Limit')
                ui.number(label="Alternative Protein Limit", min=0, max=10000, precision=0, step=1).bind_value(
                    config, 'max_allowed_alternative_proteins'
                )
                ui.markdown(
                    'Optional: Set the maximum number of proteins that a single PSM can map to while performing '
                    'the inference analysis._'
                )

            ui.markdown('## Scoring')
            with ui.grid(columns=3):
                ui.label('Protein Scoring Method')
                ui.select(
                    {
                        "multiplicative_log": "Multiplicative Log",
                        "best_peptide_per_protein": "Best Peptide Per Protein",
                        "iterative_downweighted_log": "Iterative Downweighted Log",
                        "downweighted_multiplicative_log": "Downweighted Multiplicative Log",
                        "top_two_combined": "Top Two Combined",
                        "geometric_mean": "Geometric Mean",
                        "additive": "Additive",
                    },
                    label="Protein Scoring Method",
                ).bind_value(config, 'protein_score')
                ui.markdown(
                    "Select how protein scores are generated. The recommended scoring method is "
                    "'Multiplicative Log'. If you choose 'Additive', 'PSM Score' below must be set to "
                    "'Peptide Filter Score' and PSM score type must be set to 'Additive' as well."
                )

                ui.label('PSM Scoring Method')
                ui.select(
                    {
                        "posterior_error_prob": "Posterior Error Probability",
                        "q-value": "q-value",
                        "score": "Peptide Filter Score (e.g. Percolator Score)",
                        "custom": "Custom (Enter Below)",
                    },
                    label="PSM Scoring Method",
                ).bind_value(config, 'psm_score')
                ui.markdown(
                    "Select which PSM score is used when generating protein scores. The recommended PSM Score "
                    "method is 'Posterior Error Probability'. If you choose 'Additive' as the protein score, "
                    "please choose 'Peptide Filter Score' here."
                )

                ui.label('Custom PSM Score')
                ui.input(label="Custom PSM Score").bind_value(config, 'psm_score_custom')
                ui.markdown(
                    '_Enter the column name (for tabular input) or XML attribute (for idXML/mzID/pepXML) of the '
                    'custom score to be used, and set "PSM Scoring Method" above to "Custom".  If you are using '
                    'one of the pre-selected scores above, leave this blank._'
                )

                ui.label('PSM Score Type')
                ui.select(
                    {
                        "multiplicative": "Multiplicative",
                        "additive": "Additive",
                    },
                    label="PSM Score Type",
                ).bind_value(config, 'psm_score_type')
                ui.markdown(
                    "Select the type of score that corresponds with PSM Score above. Select 'Multiplicative' "
                    "unless your PSM score above is 'Percolator Score',"
                    " then choose 'Additive'."
                )
            ui.separator()

            ui.markdown('## Identifier Information')
            with ui.grid(columns=3):
                ui.label('Decoy Identifier')
                ui.input(label="Decoy Identifier").bind_value(config, 'decoy_symbol')
                ui.markdown(
                    '_Enter the prefix used to denote decoy hits within the dataset. '
                    'Often this is one of `DECOY_`, `REV_`, or `##`._'
                )

                ui.label('Isoform Delimiter')
                ui.input(label="Isoform Delimiter").bind_value(config, 'isoform_symbol')
                ui.markdown(
                    '_Enter the symbol used to distinguish isoform proteins and canonical proteins. '
                    'This is usually `-`._'
                )

                ui.label('Reviewed Identifier Prefix')
                ui.input(label="Reviewed Identifier Prefix").bind_value(config, 'reviewed_identifier_symbol')
                ui.markdown(
                    '_Enter the symbol used to distinguish reviewed proteins and unreviewed proteins. '
                    'This is usually `sp|`. If your database does not contain reviewed Proteins, '
                    'this can be left as `sp|`_'
                )

                ui.label('Split Identifiers in Fasta')
                ui.switch().bind_value(config, 'identifier_splitting')
                ui.markdown(
                    '_If the database contains identifiers in full Uniprot format starting with "sp|" or "tr|", '
                    'while the identifiers in the results omit these prefixes, enabling this option will split off '
                    'leading "sp|" or "tr|" from the database identifiers to ensure correct matching to the results. '
                    'This setting is usually not necessary._'
                )

            ui.separator()

            ui.markdown('## Peptide Digest Method')
            with ui.grid(columns=3):
                ui.label('Proteolytic Enzyme')
                ui.select(
                    {
                        "arg-c": "Arg C",
                        "asp-n": "Arg N",
                        "bnps-skatole": "BNPS Skatole",
                        "caspase 1": "Caspase 1",
                        "caspase 2": "Caspase 2",
                        "caspase 3": "Caspase 3",
                        "caspase 4": "Caspase 4",
                        "caspase 5": "Caspase 5",
                        "caspase 6": "Caspase 6",
                        "caspase 7": "Caspase 7",
                        "caspase 8": "Caspase 8",
                        "caspase 9": "Caspase 9",
                        "caspase 10": "Caspase 10",
                        "chymotrypsin high specificity": "Chymotrypsin High Specificity",
                        "chymotrypsin low specificity": "Chymotrypsin Low Specificity",
                        "clostripain": "Clostripain",
                        "cnbr": "CNBR",
                        "enterokinase": "Enterokinase",
                        "factor xa": "Factor XA",
                        "formic acid": "Formic Acid",
                        "glutamyl endopeptidase": "Glutamyl Endopeptidase",
                        "granzyme b": "Granzyme B",
                        "hydroxylamine": "Hydroxylamine",
                        "iodosobenzoic acid": "Iodosobenzoic Acid",
                        "lysc": "LsyC",
                        "ntcb": "NTCB",
                        "pepsin ph1.3": "Pepsin ph1.3",
                        "pepsin ph2.0": "Pepsin ph2.0",
                        "proline endopeptidase": "Proline Endopeptidase",
                        "proteinase k": "Proteinase K",
                        "staphylococcal peptidase i": "Staphylococcal Peptidase I",
                        "thermolysin": "Thermolysin",
                        "thrombin": "Thrombin",
                        "trypsin": "Trypsin",
                    },
                    label="Proteolytic Enzyme",
                ).bind_value(config, 'digest_type')
                ui.markdown(
                    "_Select the digest type that was used on your searches. Can be any expasy rule as "
                    "defined here: https://pyteomics.readthedocs.io/en/latest/_modules/pyteomics/parser.html. "
                    "This is most commonly 'Trypsin'._"
                )

                ui.label('Allowed Missed Cleavages')
                ui.number(label="Allowed Missed Cleavages", min=0, max=50, precision=0, step=1).bind_value(
                    config, 'missed_cleavages'
                )
                ui.markdown(
                    '_The number of missed cleavages allowed as entered into the search parameters for this dataset._'
                )

            ui.separator()

            ui.markdown('## Parsimony Method Options')
            with ui.grid(columns=3):
                ui.label('Shared Peptide Assignment')
                ui.select(
                    {
                        "all": "All",
                        "best": "Best",
                    },
                    label="Shared Peptide Assignment",
                ).bind_value(config, 'shared_peptides')
                ui.markdown(
                    "_How to assign shared peptides for parsimony. 'All' assigns shared peptides to all possible"
                    " proteins in the output. 'Best' assigns shared peptides to the best scoring protein"
                    " which is a 'winner take all' approach. This is specific to the 'Parsimony' Inference type._"
                )

            ui.separator()

            ui.markdown('## Peptide-Centric Method Options')
            with ui.grid(columns=3):
                ui.label('Maximum Protein Group Size')
                ui.number(label="Maximum Protein Group Size", min=0, max=10000, precision=0, step=1).bind_value(
                    config, 'max_identifiers'
                )
                ui.markdown(
                    "_The maximum number of proteins to constitute a 'group' in 'Peptide Centric' inference. "
                    "Ignored for other inference methods._"
                )

        with ui.tab_panel(about):
            ui.markdown(
                r"""
            # pyProteinInference
            Version 1.1.0
            
            Trent Hinkle and Corey Bakalarski
            
            &copy; 2024 Genentech, Inc.
            """
            )


ui.run(native=True, title="pyProteinInference", reload=False, port=native.find_open_port())
