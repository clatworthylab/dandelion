from setuptools import setup, find_packages
import pkg_resources

with open("README.md", "r") as readme_file:
    readme = readme_file.read()
exec(open('dandelion/version.py').read())

requirements = ["numpy>=1.18.4", "pandas>=1.0.3", "changeo>=1.0.0", "anndata>=0.7.1", "scanpy>=1.4.6", "distance>=0.1.3", "joblib>=0.14.1", "scikit-learn>=0.23.0", "scikit-bio>=0.5.6", "scipy>=1.4.1", "numba>=0.48.0", "seaborn>=0.10.1", "networkx>=2.4", "leidenalg>=0.8.0", "plotnine>=0.6.0", "polyleven>=0.5", "h5py<3.0.0"]
setup(
    name="sc-dandelion",
    version=__version__,
    author="zktuong",
    author_email="kt16@sanger.ac.uk",
    description="sc-BCR analysis tool",
    long_description=readme,
    long_description_content_type="text/markdown",
    url="https://github.com/zktuong/dandelion/",
    packages=find_packages(),
    setup_requires=['setuptools_scm'],
    install_requires=requirements,
    package_data={'dandelion': ['bin/tigger-genotype.R']},
    data_files=[('bin', ['bin/tigger-genotype.R'])],
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)",
    ],
    zip_safe=False
)

