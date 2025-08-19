# Nardini Online

🧬 **Blazing fast protein analysis**

We provide a simple, accessible interface for running **NARDINI** - a powerful tool for analyzing Intrinsically Disordered Regions (IDRs) in proteins. Built with Modal and FastAPI, it offers both a user-friendly Jupyter notebook interface and a REST API for programmatic access.

## 🚀 Features

### ⚡ High Performance
- **Parallel Processing**: Processes multiple sequences simultaneously using Modal's distributed computing
- **Smart Caching**: Previously analyzed sequences are cached to avoid reprocessing
- **Scalable Infrastructure**: Runs on Modal's cloud platform for unlimited scalability

### 🎯 User-Friendly Interface
- **Jupyter Notebook**: Interactive notebook with step-by-step guidance
- **REST API**: Programmatic access for integration with other tools
- **Real-time Progress**: Monitor analysis progress with live updates
- **Easy Downloads**: Automated result retrieval and organization

### 📊 Comprehensive Results
- **TSV Data Files**: Detailed statistical results for each sequence
- **PNG Visualizations**: Plots and charts for data interpretation
- **Merged Reports**: Combined results across all analyzed sequences
- **Organized Output**: Structured file organization for easy access

## 🚀 Getting Started

**Ready to analyze your protein sequences?** 

Simply open our **[Google Colab Notebook](https://colab.research.google.com/drive/1jKDu1AfOHI-4BM1hHELIDFQCqbPz86M9)** - no installation required! The notebook provides a step-by-step walkthrough that will guide you through the entire process.

*For detailed setup and deployment instructions, see our [Installation Guide](GUIDE.md).*

## 📖 How to Use RunFasta

### 📔 Interactive Notebook (Recommended for Lab Members)

The easiest way to get started is through our Google Colab notebook:

1. **Open the Colab notebook** (link provided above)
2. **Upload your FASTA file** containing protein sequences
3. **Run the analysis** - the notebook will handle everything automatically
4. **Monitor progress** in real-time as sequences are processed
5. **Download your results** - get comprehensive analysis files and visualizations

The notebook is designed to be beginner-friendly with clear explanations at each step.

### 🔧 For Advanced Users

If you prefer working with APIs or integrating RunFasta into your own tools, we also provide a REST API interface. See our [Technical Documentation](GUIDE.md) for detailed API information.

## 📊 Simple Workflow

1. **Prepare your FASTA file** - Make sure your protein sequences are in FASTA format
2. **Open the Colab notebook** - Click the link above to get started
3. **Upload and analyze** - Follow the guided steps in the notebook
4. **Get your results** - Download comprehensive analysis files and visualizations
5. **Interpret findings** - Review the statistical analysis and plots to understand your protein patterns

## 🙏 Credits

**✨ Created by Tanuj Vasudeva and Ethan Caine, 2025 ✨**

### Acknowledgments

- **Dr. John Woolford** at Carnegie Mellon University for project support
- **Modal** for hosting infrastructure
- **Katherine Parry** for helpful advice
- **NARDINI Authors** for the original analysis tool

## 📚 References

Cohan, M. C., Shinn, M. K., Lalmansingh, J. M., & Pappu, R. V. (2021). Uncovering non-random binary patterns within sequences of intrinsically disordered proteins. *Journal of Molecular Biology*, 434(2), 167373.

## 📄 License

This project is licensed under the MIT License. We hope to advance scientific research by making our tools freely accessible to the global research community.

---

**Need help?** Check out the demo notebook for step-by-step instructions, or refer to the API documentation in [GUIDE.md](GUIDE.md) for programmatic usage.