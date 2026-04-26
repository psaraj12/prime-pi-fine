# prime-pi-fine

A high-performance prime-counting tool for computing **π(x)** — the number of primes less than or equal to `x`.

This tool uses a persistent fine-grained prefix table so that repeated prime-counting queries become very fast after an initial build.

## What it does

Instead of recomputing prime counts from scratch for every query, `prime-pi-fine` builds a cumulative checkpoint table:

```text
checkpoint[i] = number of primes up to that checkpoint
```

For a query `π(x)`, it does:

```text
π(x) = stored prefix count + small local count
```

The local count is limited to at most one checkpoint interval.

Current checkpoint interval:

```text
30,000 numbers
```

So each query is approximately:

```text
O(1) lookup + tiny segmented sieve
```

## Key idea

This tool trades initial setup cost for very fast repeated queries.

It is similar to building an index in a database:

```text
Build once
Query many times
```

## Features

- Computes exact `π(x)`
- Wheel-30 segmented sieve
- Fine-grained checkpoints every 30,000 numbers
- Persistent binary table file
- Fast repeated queries
- Parallel build using OpenMP
- Extend mode to grow an existing table
- Auto-aligns build/extend targets to safe checkpoint boundaries
- Correctness-safe block/checkpoint alignment

## Build

```bash
g++ -O3 -march=native -fopenmp -std=c++17 prime_pi_fine_autoalign.cpp -o prime_pi_fine
```

If your system does not support OpenMP, remove `-fopenmp`, but build performance will be slower.

## Usage

### Build a table

```bash
./prime_pi_fine build 1000000000000 12 pi.dat
```

Arguments:

```text
build <end> [threads] [table_file]
```

Example:

```bash
./prime_pi_fine build 1000000000000 12 pi.dat
```

This builds a table up to approximately `10^12`.

The program may round the target upward to the next safe checkpoint boundary.

Example:

```text
Requested end 1000000000000 is not checkpoint-aligned; using 1000000019999 instead.
```

### Query π(x)

```bash
./prime_pi_fine query 123456789012 pi.dat
```

Example output:

```text
π(123456789012) = ...
Query time: ... µs
```

### Extend an existing table

```bash
./prime_pi_fine extend 2000000000000 12 pi.dat
```

Arguments:

```text
extend <new_end> [threads] [table_file]
```

This extends the existing table instead of rebuilding from zero.

### Show table info

```bash
./prime_pi_fine info pi.dat
```

## Example workflow

```bash
# Compile
g++ -O3 -march=native -fopenmp -std=c++17 prime_pi_fine_autoalign.cpp -o prime_pi_fine

# Build table
./prime_pi_fine build 1000000000000 12 pi.dat

# Query values
./prime_pi_fine query 100 pi.dat
./prime_pi_fine query 1000000 pi.dat
./prime_pi_fine query 1000000000 pi.dat
./prime_pi_fine query 1000000000000 pi.dat

# View metadata
./prime_pi_fine info pi.dat
```

## Known reference values

You can use these to sanity-check correctness:

```text
π(10)             = 4
π(100)            = 25
π(1,000)          = 168
π(10,000)         = 1,229
π(100,000)        = 9,592
π(1,000,000)      = 78,498
π(1,000,000,000)  = 50,847,534
π(1,000,000,000,000) = 37,607,912,018
```

## Design overview

### 1. Wheel-30 representation

Only numbers coprime to 30 are represented:

```text
1, 7, 11, 13, 17, 19, 23, 29 mod 30
```

This avoids storing obvious composites divisible by 2, 3, or 5.

### 2. Segmented sieve

The build process sieves large blocks for cache-efficient prime counting.

### 3. Fine prefix checkpoints

The table stores cumulative prime counts every 30,000 numbers.

This makes queries fast because only the remainder after the nearest checkpoint must be counted.

### 4. Auto-alignment

For correctness, sieve blocks and checkpoint boundaries must align cleanly.

The program automatically rounds requested build/extend targets upward to a safe checkpoint end.

## Performance

Observed performance on a 12-thread machine:

```text
Build up to ~10^12: around 5 minutes
Query time: microseconds after table is built
```

Actual performance depends on CPU, memory bandwidth, compiler, operating system, and thread count.

## When this tool is useful

Good fit:

```text
Many repeated π(x) queries
Prime distribution experiments
Number theory exploration
Benchmarking
Precomputed prime-count lookup
```

Not ideal:

```text
One-off π(x) query only
Very small x where direct methods are enough
```

## Important notes

- This computes prime counts, not Goldbach decompositions.
- The table file can become large for very high ranges.
- Query speed depends on the table already covering `x`.
- Use the same program version for build, extend, and query to avoid format mismatches.
- Avoid extending older non-aligned table files; rebuild with this version.

## License

MIT License


## Author

Created by Santhiagu as part of high-performance number theory experiments.
