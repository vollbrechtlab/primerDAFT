from setuptools import setup, find_packages

setup(name='primerDAFT',
      version='0.2',
      description='A Python based PCR Primer Design and Filtering Package',
      url='https://github.com/vollbrechtlab/primerDAFT',
      author='Takao Shibamoto',
      author_email='takaos@iastate.edu',
      license='MIT',
      packages=find_packages(exclude=['venv','test','sphinx','docs','diagram']),
      scripts=['bin/asd.py'],
      install_requires=[
          'primer3-py',
          'biopython',
          'pandas',
          'pysam'
      ],
      zip_safe=False)
