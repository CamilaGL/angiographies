from setuptools import setup, find_namespace_packages

setup(name='angiographies',
      packages=find_namespace_packages(include=["angiographies", "angiographies.*"]),
      version='0.1',
      description='My framework for processing brain angiographies.',
      url='https://github.com/CamilaGL/angiographies',
      author='Camila Garcia',
      entry_points={
          'console_scripts': [
              'angiographies_vmtknetwork = angiographies.skeletonisation.vmtknetwork:main',
              
          ],
     }
    )