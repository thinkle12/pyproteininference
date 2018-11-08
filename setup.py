from setuptools import setup
import glob

setup(
    setup_requires=['pbr>=1.8', 'setuptools>=17.1'],
    pbr=True,
    scripts=glob.glob("scripts/*.py"),
    name='protein_inference',
    version='1.0.0',
    packages=['digest', 'protein_inference','scripts'],
    url='',
    license='',
    author='hinklet',
    author_email='hinklet@gene.com',
    description='Python Package for running custom Protein Inference on output from Percolator',
    package_data = {
    'glpkinout/': ['glpkinout/glpkout_159260_Bioplex2_b10090.mod',
                   'glpkinout/glpkout_159260_Bioplex2_b10090.sol'],
    'data/': ['data/159260_Bioplex2_b10090_percolator_decoy_psm.txt',
                           'data/159260_Bioplex2_b10090_percolator_target_psm.txt'
                           'data/UniprotKBConcat1708_HUMAN.fasta'],
    'output/': ['output/all_idwl_159260_Bioplex2_b10090_q.csv',
                'output/csep_idwl_159260_Bioplex2_b10090_q.csv',
                'output/leads_idwl_159260_Bioplex2_b10090_q.csv',
                'output/qvalues_all_idwl_159260_Bioplex2_b10090_q.csv',
                'output/qvalues_csep_idwl_159260_Bioplex2_b10090_q.csv',
                'output/qvalues_leads_idwl_159260_Bioplex2_b10090_q.csv'],
    'plots/': ['plots/159260_Bioplex2_b10090_plot_idwl.pdf'],
    'scripts':['scripts/Command_Line_PI_runner.py',
                                'scripts/practice_digest.py',
                                'scripts/runner.py']}
)
