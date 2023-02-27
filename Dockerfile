FROM julia:1.9.0-beta4

RUN mkdir /home/outputs

RUN julia -e "using Pkg; Pkg.add(\"DataFrames\"); \
  Pkg.add(\"StatFiles\"); \
  Pkg.add(\"BenchmarkTools\"); \
  Pkg.add(\"StatsBase\"); \
  Pkg.add(\"Bootstrap\"); \
  Pkg.add(\"Random\"); \
  Pkg.add(\"CSV\"); \
  Pkg.add(\"Tidier\"); \
  Pkg.add(\"Inequality\");"


COPY data/ /home/data/
COPY benchmarks.jl /home/benchmarks.jl

CMD  julia -e "include(\"/home/benchmarks.jl\")"

# Windows: docker run -v ${PWD}/outputs/:/home/outputs/ r_image