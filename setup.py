from setuptools import setup, find_namespace_packages

setup(name='angiographies',
      packages=find_namespace_packages(include=["angiographies", "angiographies.*"]),
      version='0.5.6',
      description='My framework for processing brain angiographies.',
      url='https://github.com/CamilaGL/angiographies',
      author='Camila Garcia',
      entry_points={
          'console_scripts': [
              'angiographies_vmtknetwork = angiographies.skeletonisation.vmtknetwork:main',
              'angiographies_networkediting = angiographies.skeletonisation.networkediting:main',
              'angiographies_orderedthinning = angiographies.skeletonisation.orderedthinning:main',
              'angiographies_thinning = angiographies.skeletonisation.orderedthinning2:main',
              'angiographies_skeletongraph = angiographies.skeletonisation.skeletongraph:main',
              'angiographies_polydatamerger = angiographies.skeletonisation.polydatamerger:main',
              'angiographies_nidusextractor = angiographies.skeletonisation.nidusextractor:main',
              'angiographies_convexhull = angiographies.skeletonisation.convexhull:main',
              'angiographies_pipeline = angiographies.skeletonisation.pipeline:main',
              'angiographies_overlap = angiographies.utils.overlapimg:main',
              'angiographies_evaluate = angiographies.utils.evaluate:main',
              
          ],
     }
    )