# Instructions

Referencing https://raw.githubusercontent.com/rvaser/spoa/master/README.md

## Exact steps

> python scripts/gen_data.py > data.in  
>
> make && ./pdbs data.in 1 3 10 1 0 50 0.50 0.90 log.out query.in
>
> cd poa_util
>
> git clone --recursive https://github.com/rvaser/spoa spoa
>
> cd spoa
>
> mkdir build
> 
> cd build
>
> cmake -DCMAKE_BUILD_TYPE=Release ..
>
> make
>
> g++ poa.cpp -std=c++11 -Ispoa/include/ -Lspoa/build/lib/ -lspoa -o poa
>
> ./poa 0 5 -4 -8 ../query.poa
