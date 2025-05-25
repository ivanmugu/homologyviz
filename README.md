<p align="center">
   <img src="https://raw.githubusercontent.com/ivanmugu/homologyviz/refs/heads/main/src/homologyviz/assets/logo.png" alt="HomologyViz" width="450">
</p>

<div align="center">

[![PyPI Version](https://img.shields.io/pypi/v/homologyviz)](https://pypi.org/project/homologyviz/)
![GitHub License](https://img.shields.io/github/license/ivanmugu/homologyviz)
![Python](https://img.shields.io/badge/python-3.11-blue.svg)
[![Built with Dash + Plotly](https://img.shields.io/badge/Built%20with-Dash%20%2B%20Plotly-1f6feb?logo=plotly&logoColor=white)](https://dash.plotly.com/)

</div>

---

# 🧬 Create visual representations of BLASTn alignments

**HomologyViz** is a Python-based web app for visualizing pairwise **BLASTn alignments** between DNA sequences with gene annotations and customizable color maps. It allows you to explore homologous regions and gene structures interactively, with **publication-ready output**.

HomologyViz reads **GenBank files (.gb)**, uses the `CDS` features to plot genes as arrows, and performs local BLASTn alignments between sequences. Homology regions are colored by identity percentage using a flexible, user-selectable color scale.

You may optionally customize gene colors directly in the GenBank file by adding a `/Color` tag to a `CDS` feature. For example, `/Color="#00ff00"` will render that gene in green. You can also **edit gene colors interactively** using the web interface.

HomologyViz is ideal for researchers and students with **little or no coding experience**. It runs locally and requires no cloud deployment or server hosting.

---

## 🚀 Features

- Interactive plots of aligned DNA sequences
- Gene annotations rendered as directional arrows
- Homology regions colored by % identity using customizable color scales
- Click-to-select genes or regions and change their color
- Exportable visualizations (SVG, PNG, PDF, HTML)
- Built with [Dash](https://dash.plotly.com/) and [Plotly](https://plotly.com/python/)

---

## 🖼 Example

<p align="center">
   <img src="https://raw.githubusercontent.com/ivanmugu/homologyviz/refs/heads/main/src/homologyviz/images/HomologyViz_app.png" alt="HomologyViz" width="600">
</p>

---

## ✅ Requirements

- Python 3.11+
- [`blastn`](https://www.ncbi.nlm.nih.gov/books/NBK569861/) (part of [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)) must be installed and available in your system `PATH`.

> HomologyViz has been tested in **Google Chrome** on **macOS**, but should work across all modern browsers and platforms.

---

## 📦 Installation

You can install HomologyViz using `pip`:

```bash
pip install homologyviz
```

Or clone this repository:

```bash
git clone https://github.com/yourusername/homologyviz.git
cd homologyviz
pip install -e .
```

---

## 🧪 Usage

Once installed, launch the app by typing:

```bash
homologyviz
```

---

## 📚 Documentation

You can find the full documentation for **HomologyViz** at:

👉 [https://homologyviz.readthedocs.io](https://homologyviz.readthedocs.io/en/latest/)

---

## 👤 Author

Created by Iván Muñoz Gutiérrez, PhD<br>
University of California, Irvine<br>
School of Biological Sciences<br>
Github: [@ivanmugu](https://github.com/ivanmugu)

---

## 🙏 Acknowledgments

- Built on Dash, Plotly, and Biopython
- Inspired by easyfig: [Sullivan et al (2011) Bioinformatics 27(7):1009-1010](https://academic.oup.com/bioinformatics/article/27/7/1009/230508)

---

## 📄 License

BSD 3-Clause License

---

## 📝 Notes

I am developing HomologyViz in my free time, so if you find a bug, it may take me some time to fix it. However, I will fix the problems as soon as possible. Also, if you have any suggestions, let me know, and I will try to implement them.
