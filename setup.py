from setuptools import setup, find_packages
import pkg_resources

with open("README.md", "r") as readme_file:
    readme = readme_file.read()

REQUIREMENTS = ["numpy>=1.18.4", "pandas>=1.0.3", "changeo>=1.0.0", "anndata>=0.7.1", "scanpy>=1.4.6", "scrublet>=0.2.1", "python-Levenshtein>=0.12.0", "distance>=0.1.3", "joblib>=0.14.1", "scikit-learn>=0.23.0", "scipy>=1.4.1", "numba>=0.48.0", "seaborn>=0.10.1", "python-igraph>=0.8.2", "networkx>=2.4", "leidenalg>=0.8.0"]
setup(
    name="sc-dandelion",
    version="0.0.1",
    author="zktuong",
    author_email="kt16@sanger.ac.uk",
    description="sc-BCR analysis tool",
    long_description=readme,
    long_description_content_type="text/markdown",
    url="https://github.com/zktuong/dandelion/",
    packages=find_packages(),
    setup_requires=['setuptools_scm'],
    scripts = ['dandelion/utilities/_blastn.py', 'dandelion/utilities/_igblastn.py', 'dandelion/utilities/_makeblastdb.py', 'dandelion/utilities/_tigger-genotype.py'],
    install_requires=REQUIREMENTS,
    package_data={'dandelion': ['bin/blastn', 'bin/igblastn', 'bin/makeblastdb', 'bin/tigger-genotype.R', 'database/blast/human/*', 'database/blast/mouse/*', 'database/igblast/database/*', 'database/igblast/fasta/*', 'database/igblast/optional_file/*', 'database/igblast/internal_data/human/*','database/igblast/internal_data/human/*', 'database/igblast/internal_data/mouse/*', 'database/igblast/internal_data/rabbit/*', 'database/igblast/internal_data/rat/*', 'database/igblast/internal_data/rhesus_monkey/*','database/germlines/imgt/IMGT.yaml','database/germlines/imgt/human/constant/*.fasta', 'database/germlines/imgt/human/leader/*.fasta', 'database/germlines/imgt/human/vdj/*.fasta', 'database/germlines/imgt/human/vdj_aa/*.fasta', 'database/germlines/imgt/mouse/constant/*.fasta', 'database/germlines/imgt/mouse/leader/*.fasta', 'database/germlines/imgt/mouse/vdj/*.fasta', 'database/germlines/imgt/mouse/vdj_aa/*.fasta']},
    # data_files=[('bin', ['bin/blastn', 'bin/igblastn', 'bin/makeblastdb', 'bin/tigger-genotype.R']),
                # ('database', ['database/blast/human/*', 'database/blast/mouse/*', 'database/igblast/database/*', 'database/igblast/fasta/*', 'database/igblast/optional_file/*', 'database/igblast/internal_data/human/*','database/igblast/internal_data/human/*', 'database/igblast/internal_data/mouse/*', 'database/igblast/internal_data/rabbit/*', 'database/igblast/internal_data/rat/*', 'database/igblast/internal_data/rhesus_monkey/*','database/germlines/imgt/IMGT.yaml','database/germlines/imgt/human/constant/*.fasta', 'database/germlines/imgt/human/leader/*.fasta', 'database/germlines/imgt/human/vdj/*.fasta', 'database/germlines/imgt/human/vdj_aa/*.fasta', 'database/germlines/imgt/mouse/constant/*.fasta', 'database/germlines/imgt/mouse/leader/*.fasta', 'database/germlines/imgt/mouse/vdj/*.fasta', 'database/germlines/imgt/mouse/vdj_aa/*.fasta'])],
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)",
    ],
    zip_safe=False
)

