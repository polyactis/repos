from setuptools import setup, find_packages
import sys, os

version = '1.0.09'

setup(name='optilux.codeexamples',
      version=version,
      description="Code examples accompanying the Optilux Cinemas case study",
      long_description="""\
""",
      # Get more strings from http://www.python.org/pypi?%3Aaction=list_classifiers
      classifiers=[
        "Framework :: Plone",
        "Framework :: Zope2",
        "Framework :: Zope3",
        "Programming Language :: Python",
        "Topic :: Software Development :: Libraries :: Python Modules",
        ],
      keywords='Zope Plone code examles',
      author='Martin Aspeli',
      author_email='optilude@gmx.net',
      url='http://packpub.com',
      license='GPL',
      packages=find_packages(exclude=['ez_setup']),
      namespace_packages=['optilux'],
      include_package_data=True,
      zip_safe=False,
      install_requires=[
          # -*- Extra requirements: -*-
      ],
      entry_points="""
      # -*- Entry points: -*-
      """,
      )
