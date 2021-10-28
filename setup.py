import setuptools
import os

base_dir = os.path.dirname(os.path.realpath(__file__))

with open(os.path.join(base_dir, "README.md"), encoding="utf-8") as in_f:
    long_description = in_f.read()

about = {}
with open(os.path.join(base_dir, "clip", "__about__.py"), encoding="utf-8") as in_f:
    exec(in_f.read(), about)

setuptools.setup(
    name=about["__title__"],
    version=about["__version__"],
    author=about["__author__"],
    author_email=about["__email__"],
    description=about["__summary__"],
    long_description=long_description,
    long_description_content_type="text/markdown",
    url=about["__url__"],
    package_dir={"": "."},
    packages=["clip"],
    python_requires=">=3.6",
    entry_points={"console_scripts": ["clippy = clip:main",],},
)
