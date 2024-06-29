from setuptools import setup
setup(
   name='INS-Analysis',
   version='1.0',
   description='Python package for Inelastic Neutron Scattering data analysis',
   author='Jose Andres Cortes',
   author_email='jose.cortes@uta.edu',
   packages=['INS-Analysis'],  #same as name
   install_requires=['pandas', 'numpy', 'scipy'], #external packages as dependencies,
)