using Coverage
using Pkg
Pkg.test("MatsubaraFunctions"; coverage=true)
coverage = process_folder() # defaults to src/; alternatively, supply the folder name as argument
covered_lines, total_lines = get_summary(coverage)

clean_folder(".")