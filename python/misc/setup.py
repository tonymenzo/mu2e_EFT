from setuptools import setup, find_packages

with open ('mu2e/_version.py') as f:
    exec(f.read())

setup(name='mu2e',
      version=__version__,
      author='Wick Haxton, Ken McElvain, Tony Menzo, Evan Rule, Jure Zupan',
      author_email='menzoad@mail.uc.edu',
      url='',
      description='A python package for muon to electron conversion theory analysis',
      long_description=""" This package contains classes for Wilson coefficients
                           for muon to electron conversion interactions, 
                           in effective theories at various scales. It allows for  
                           renormalization-group running, and matching between the 
                           different effective theories.""",
      license='MIT',
      packages=find_packages(),
      install_requires=['numpy', 'scipy', 'setuptools']
    )
