from setuptools import setup, find_packages


setup(
    name="more_dynamics",
    version="0.1.0",
    packages=find_packages(include=["more_dynamics", "more_dynamics.*"]),
    include_package_data=True,
    install_requires=[
        "casadi>=3.5.5",
    ],
)