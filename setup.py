from setuptools import setup

setup(
    name="quarTeT",
    version="1.2.5",
    author="Echoring",
    author_email="linyunzhi20@gmail.com",
    description="A telomere-to-telomere toolkit for gap-free genome assembly and centromeric repeat identification",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/aaranyue/quarTeT",
    py_modules=['quartet', 'quartet_assemblymapper', 'quartet_centrominer', 'quartet_gapfiller', 'quartet_teloexplorer', 'quartet_util'],
    python_requires='>=3.6',
    entry_points={
        'console_scripts': [
            'quartet=quartet:main',
        ],
    },
)
