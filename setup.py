from setuptools import setup

setup(name='primer blast dx',
      version='0.1',
      description='simple primer blast',
      url='https://github.com/takao42/primer-blast-dx',
      author='Takao',
      author_email='takaos@iastate.edu',
      license='MIT',
      packages=['primer_blast_dx'],
      install_requires=[
          'primer3-py',
          'biopython',
          'pandas',
          'pysam'
      ],
      zip_safe=False)
