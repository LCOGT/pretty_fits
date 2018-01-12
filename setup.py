from setuptools import setup

setup(
    name='pretty_fits',
    version='0.1.1',
    packages=['pretty_fits'],
    install_requires=[
        'astropy',
        'astroscrappy',
        'Pillow',
        'alipy'
    ],
    dependency_links=[
    'git+git://github.com/LCOGT/alipy.git#egg=alipy',
    'git+git://github.com/LCOGT/fits2image.git#egg=fits2image'
    ],
    entry_points='''
        [console_scripts]
        pretty_fits=pretty_fits:run
    ''',

classifiers=[

    'Development Status :: 3 - Alpha',

    'Intended Audience :: Developers',
    'Topic :: Software Development :: Build Tools',

    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.2',
    'Programming Language :: Python :: 3.3',
    'Programming Language :: Python :: 3.4',
],


)
