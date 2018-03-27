from setuptools import setup

setup(name='primerDAFT',
      version='0.1',
      description='A Python based PCR Primer Design and Filtering Package',
      url='https://github.com/vollbrechtlab/primerDAFT',
      author='Takao Shibamoto',
      author_email='takaos@iastate.edu',
      license='MIT',
      packages=['primerDAFT'],
      scripts=['bin/asd.py'],
      install_requires=[
          'primer3-py',
          'biopython',
          'pandas',
          'pysam'
      ],
      zip_safe=False)
