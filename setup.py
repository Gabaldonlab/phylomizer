from ez_setup import use_setuptools
use_setuptools()

from setuptools import setup, find_packages

import versioneer
versioneer.versionfile_source = 'source/_version.py'
versioneer.versionfile_build = 'source/_version.py'
versioneer.tag_prefix = '1.0.0'
versioneer.parentdir_prefix = 'phylomizer-'

setup(name="phylomizer",
      description="Phylogenetic tree reconstruction pipeline",
      long_description="""
        Phylogenetic tree reconstruction pipeline containing three independent
        modules which resemble the steps followed by a classical phylogenetist
        interested on building the phylogeny of a given gene family - The
        pipeline automatically performs all steps and control if something
        goes wrong - if that's the case then it can restart from the failure
        point rather than from the beginning of the whole process
      """,
      author="Salvador Capella-Gutierrez, Toni Gabaldon",
      author_email="[scapella, tgabaldon]@crg.es",
      license="GNU GPL 3.0",
      url="http://big.crg.cat/comparative_genomics",
      packages=find_packages(exclude=('others','bin','*.pyc')),
      include_package_data=True,
      install_requires=['Bio', 'subprocess', 'socket'],
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
)
