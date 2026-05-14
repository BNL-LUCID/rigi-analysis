"""Command-line interface dispatcher for rigi_analysis.

This module provides the central entry point `rigi-analysis-run` that dynamically
loads and executes individual scripts from the `mutation`, `sv`, and `workflows` modules.
"""

import sys
import importlib


def main():
    """Execute the requested module's main function.

    Parses the script name from the command line arguments, dynamically searches
    for the corresponding module within the rigi_analysis subpackages, rewrites
    sys.argv to ensure argparse within the target script formats help messages
    correctly, and finally calls the module's `main()` function.

    Raises:
        SystemExit: If no script name is provided, if the script cannot be found,
            or if the target script does not define a `main()` function.
    """
    if len(sys.argv) < 2:
        print("Usage: rigi-analysis-run <script_name> [args...]")
        print("\nAvailable modules are dynamically resolved from mutation, sv, or workflows.")
        sys.exit(1)
        
    script_name = sys.argv[1]
    
    # Remove the script_name from sys.argv so argparse in the target script 
    # receives the correct arguments
    sys.argv.pop(1)
    
    # Update sys.argv[0] so that argparse prints the correct command usage
    sys.argv[0] = f"rigi-analysis-run {script_name}"
    
    # Try to import from mutation first, then sv, then workflows
    namespaces = [
        f"rigi_analysis.mutation.{script_name}",
        f"rigi_analysis.sv.{script_name}",
        f"rigi_analysis.workflows.{script_name}"
    ]
    
    module = None
    for ns in namespaces:
        try:
            module = importlib.import_module(ns)
            break
        except ImportError as e:
            # If the module itself wasn't found, continue checking other namespaces.
            # If it was found but failed to import something else, raise the error.
            if getattr(e, 'name', None) == ns:
                continue
            raise
            
    if module is None:
        print(f"Error: Could not find script '{script_name}' in rigi_analysis.")
        print("Searched in: rigi_analysis.mutation, rigi_analysis.sv, rigi_analysis.workflows")
        sys.exit(1)
        
    if not hasattr(module, 'main'):
        print(f"Error: Script '{script_name}' does not have a main() function.")
        sys.exit(1)
        
    # Execute the script's main function
    module.main()


if __name__ == "__main__":
    main()
