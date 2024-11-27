<p align="center">
   <img src="src/homologyviz/images/logo.png" alt="HomologyViz" width="450">
</p>

<div align="center">

![PyPI Version](https://img.shields.io/pypi/v/homologyviz)
![GitHub License](https://img.shields.io/github/license/ivanmugu/homologyviz)
![Python](https://img.shields.io/badge/python-3.11-blue.svg)

</div>

---

# Make a graphical representation of BLASTn alignments
Homology Visualization (HomologyViz) uses GenBank files (.gb) to align the sequences and plot the genes. HomologyViz uses the information from the `CDS` features section to plot the genes. To customize the colors for plotting genes, you can add a `Color` tag in the `CDS` features with a color in hexadecimal. For example, add the tag `/Color="#00ff00"` to show a green gene. Or, you can edit the colors interactively in the plot.

HomologyViz is an easy-to-use option for people with little coding knowledge. The program uses Dash for an interactive and web-friendly experience. HomologyViz is a flexible app that allows you to export your graph in different formats and sizes for publication quality.

## Requirements

- [blastn](https://www.ncbi.nlm.nih.gov/books/NBK569861/) must be installed
  locally and in the path

HomologyViz has been tested in Chrome using macOS.

## Installation

First, create a virtual environment with `conda` or `venv`. Then, install
homologyviz using pip as follows:

```bash
pip install homologyviz
```

## Usage

To run the app type:

```bash
homologyviz
```

Output:
<p align="center">
   <img src="src/homologyviz/images/HomologyViz_app.png" alt="HomologyViz" width="600">
</p>


## Credits

Inspired by easyfig: Sullivan et al (2011) Bioinformatics 27(7):1009-1010

## License

BSD 3-Clause License

## Notes

I am developing HomologyViz in my free time, so if you find a bug, it may take me some time to fix it. However, I will fix the problems as soon as possible. Also, if you have any suggestions, let me know, and I will try to implement them.