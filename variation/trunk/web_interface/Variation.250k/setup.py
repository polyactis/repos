from setuptools import setup, find_packages
import os

version = '0.1'

setup(name='Variation.250k',
      version=version,
      description="Variation 250k web portal",
      long_description=open("README.txt").read() + "\n" +
                       open(os.path.join("docs", "HISTORY.txt")).read(),
      # Get more strings from http://www.python.org/pypi?%3Aaction=list_classifiers
      classifiers=[
        "Framework :: Plone",
        "Programming Language :: Python",
        "Topic :: Software Development :: Libraries :: Python Modules",
        ],
      keywords='SNP polymorphism 250k DB',
      author='Yu Huang',
      author_email='yuhuang@usc.edu',
      url='http://papaya.usc.edu/',
      license='GPL',
      packages=find_packages(exclude=['ez_setup']),
      namespace_packages=['Variation'],
      include_package_data=True,
      zip_safe=False,
      install_requires=[
          'setuptools',
          # -*- Extra requirements: -*-
      ],
      entry_points="""
      # -*- Entry points: -*-
      """,
      )
