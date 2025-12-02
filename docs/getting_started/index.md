# Getting Started with ``MatsubaraFunctions.jl``

Welcome! In this book, you will find minimal and advanced examples, showing how to use the MatsubaraFunctions.jl package for working with Matsubara frequency meshes and related functions in Julia.

You can either just read this documentation online, or you can clone the repository and run the code examples interactively on your own machine. To do so, follow these steps:

1. **Install Julia**: If you haven't already, download and install Julia from the [official website](https://julialang.org/downloads/).

2. **Clone the Repository**: Open your terminal and run the following command to clone the repository:
   ```bash
   gh repo clone dominikkiese/MatsubaraFunctions.jl
   ```
3. **Navigate to the Directory**: Change into the cloned repository's directory:
   ```bash
   cd MatsubaraFunctions.jl
   ```

4. **Install IJulia**: Install IJulia in your Julia distribution and install a named kernel that includes ``--project=@.`` so that the notebooks in this repo will activate the correct project automatically.

    From the Julia REPL, run this once on your machine:
    ```julia
        using Pkg
        Pkg.add("IJulia")
        using IJulia
        # install a kernel that launches Julia with --project=@.
        installkernel("Julia (project)", "--project=@.")
    ```

    That way, the default kernel installed by IJulia is set up to call Julia with ``--project=@.``; that causes a Project.toml in the notebook folder (or parent) to become the active project automatically.

6. **Launch and edit Jupyter Notebook**: You can now navigate to ``docs/getting_started/`` and edit and run the notebooks either from within your editor of choice (e.g., VS Code) or by launching Jupyter Notebook from the terminal:
   ```bash
   jupyter notebook
   ```
    This requires you to have jupyter installed. Make sure to select the kernel named "Julia (project)" when opening the notebooks.